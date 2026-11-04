#!/usr/bin/env node
// Focused loopback check: axis display must not filter donor means, summaries, heatmap or CSV.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const options = {};
const allowed = ['url', 'output', 'source', 'runtime-kind', 'image', 'reference', 'chrome'];
for (let i = 2; i < process.argv.length; i += 2) {
  const name = process.argv[i]?.replace(/^--/, '');
  assert(allowed.includes(name) && process.argv[i + 1], 'Unexpected or incomplete argument');
  options[name] = process.argv[i + 1];
}
for (const name of ['url', 'output', 'source', 'runtime-kind']) assert(options[name], `--${name} required`);
const target = new URL(options.url);
assert(target.protocol === 'http:' && ['127.0.0.1', 'localhost', '[::1]'].includes(target.hostname) &&
  !target.username && !target.password && !target.search && !target.hash, 'Unadorned loopback HTTP URL required');
assert(/^[a-f0-9]{40}$/.test(options.source), 'Full source commit required');
assert(['source-overlay', 'native-image'].includes(options['runtime-kind']), 'Explicit runtime kind required');
if (options['runtime-kind'] === 'native-image') assert(/^sha256:[a-f0-9]{64}$/.test(options.image || ''), 'Image ID required');
if (!target.pathname.endsWith('/')) target.pathname += '/';
const reference = options.reference ? JSON.parse(await fs.readFile(options.reference, 'utf8')) : null;
const report = { status: 'running', startedAt: new Date().toISOString(), target: target.href,
  identity: { source: options.source, runtimeKind: options['runtime-kind'], image: options.image || null,
    binding: 'Caller supplied; verify against the companion runtime inspection.' }, checks: [], pageErrors: 0,
  externalRequestsBlocked: 0,
  limitations: ['Loopback checks do not establish public access, scientific typing accuracy, or image security.'] };
const save = async () => {
  await fs.mkdir(path.dirname(path.resolve(options.output)), { recursive: true });
  await fs.writeFile(options.output, JSON.stringify(report, null, 2) + '\n', { mode: 0o600 });
};
let browser, page;
const asArray = value => Array.isArray(value) ? value : [value];
async function snapshot() {
  return page.evaluate(() => {
    const plot = document.getElementById('spatial-enrichment_bars');
    const heat = document.getElementById('spatial-enrichment_heatmap')?.data?.find(t => t.type === 'heatmap');
    const traces = plot?.data || [];
    return { bars: traces.filter(t => t.type === 'bar').map(t => ({
      name: t.name, group: t.legendgroup, x: t.x, y: t.y, text: t.text, plus: t.error_y?.array, minus: t.error_y?.arrayminus })),
    points: traces.filter(t => t.type === 'scatter').map(t => ({ group: t.legendgroup, x: t.x, y: t.y, text: t.text, symbol: t.marker?.symbol })),
    xaxis: plot?.layout?.xaxis, legend: plot?.layout?.legend, yRange: plot?._fullLayout?.yaxis?.range,
    fittedRange: plot?.layout?.yaxis?.autorange === false ? plot.layout.yaxis.range : null,
    heat: heat && { x: heat.x, y: heat.y, z: heat.z } };
  });
}
async function csvRows() {
  return page.evaluate(async () => {
    const link = document.querySelector('#spatial-dl_enrichment_csv');
    if (!link?.href) throw new Error('CSV link missing');
    const response = await fetch(link.href);
    if (!response.ok) throw new Error('CSV request failed');
    return (await response.text()).split(/\r?\n/).filter(line => line && !line.startsWith('#')).join('\n');
  });
}
async function setRadio(name, value) {
  const selector = `input[name=${JSON.stringify(name)}][value=${JSON.stringify(value)}]`;
  const already = await page.$eval(selector, el => el.checked);
  if (already) return;
  const prior = await page.evaluate(() => window.__enrichmentBarsUpdates);
  await page.$eval(selector, el => el.click());
  await page.waitForFunction((name, value, prior) =>
    window.Shiny?.shinyapp?.$inputValues?.[name] === value && window.__enrichmentBarsUpdates > prior &&
    !document.documentElement.classList.contains('shiny-busy') &&
    !document.getElementById('spatial-enrichment_bars')?.classList.contains('recalculating'), {}, name, value, prior);
}
async function checkScenario(name) {
  await setRadio('spatial-enrich_axis_range', 'bulk');
  const bulk = await snapshot(), firstCsv = await csvRows();
  await setRadio('spatial-enrich_axis_range', 'all');
  const all = await snapshot(), secondCsv = await csvRows();
  assert(JSON.stringify(all.bars) === JSON.stringify(bulk.bars), 'Axis choice changed summary bars or intervals');
  assert(JSON.stringify(all.heat) === JSON.stringify(bulk.heat), 'Axis choice changed the donor heatmap');
  assert(secondCsv === firstCsv, 'Axis choice changed CSV data rows');
  assert.equal(all.xaxis.type, 'linear');
  assert(JSON.stringify(all.xaxis.ticktext) === JSON.stringify(all.heat.x), 'Phenotype labels differ between charts');
  assert.equal(all.legend?.groupclick, 'togglegroup');
  assert(all.points.length && all.points.every(t => t.symbol === 'circle'), 'Full range must show every donor point');
  const values = all.points.flatMap(t => asArray(t.y));
  assert.equal(values.length, all.heat.z.flat().filter(Number.isFinite).length, 'Donor points missing');
  assert(all.yRange[0] <= Math.min(...values) && all.yRange[1] >= Math.max(...values), 'Full range misses donor points');
  for (const trace of all.points) {
    const texts = asArray(trace.text), xs = asArray(trace.x);
    for (let i = 0; i < texts.length; i++) {
      const [phenotype, donor] = texts[i].split('<br>');
      const status = donor.match(/\(([^)]+)\)/)?.[1], bar = all.bars.find(t => t.name === status);
      const index = bar && asArray(bar.text).indexOf(phenotype);
      assert(bar && index >= 0 && Math.abs(xs[i] - asArray(bar.x)[index]) <= 0.071, 'Donor point is not aligned with its bar');
      assert.equal(trace.group, bar.group, 'Donor point must follow its status legend');
    }
  }
  const clipped = bulk.points.filter(t => ['triangle-up-open', 'triangle-down-open'].includes(t.symbol));
  const bulkCount = bulk.points.reduce((n, t) => n + asArray(t.y).length, 0);
  assert.equal(bulkCount, values.length, 'Fitted range dropped donor points');
  for (const trace of clipped) {
    assert(bulk.fittedRange && asArray(trace.text).every(t => t.includes('Outside the displayed axis range')),
      'Off-scale marker lacks its true-value hover');
  }
  report.checks.push({ name, passed: true, donorPoints: values.length,
    offScalePoints: clipped.reduce((n, t) => n + asArray(t.y).length, 0),
    summariesUnchanged: true, heatmapUnchanged: true, csvRowsUnchanged: true });
  return all;
}
try {
  browser = await puppeteer.launch({ executablePath: options.chrome || '/usr/bin/google-chrome', headless: true,
    args: ['--disable-dev-shm-usage', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'] });
  page = await browser.newPage(); page.setDefaultTimeout(120000);
  page.on('pageerror', () => { report.pageErrors++; });
  await page.setRequestInterception(true);
  page.on('request', request => {
    const requested = new URL(request.url());
    if (['http:', 'https:'].includes(requested.protocol) && requested.origin !== target.origin) {
      report.externalRequestsBlocked++; void request.abort('blockedbyclient');
    } else void request.continue();
  });
  await page.setViewport({ width: 1600, height: 1100, deviceScaleFactor: 1 });
  await page.goto(target.href, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
  await page.evaluate(() => {
    window.__enrichmentBarsUpdates = 0;
    jQuery(document).on('shiny:value.enrichmentAxis', event => {
      if (event.name === 'spatial-enrichment_bars') window.__enrichmentBarsUpdates++;
    });
  });
  await page.click('#tabs a[data-value="Spatial"]');
  await page.waitForFunction(() => {
    const bars = document.getElementById('spatial-enrichment_bars');
    const heat = document.getElementById('spatial-enrichment_heatmap');
    return bars?._fullLayout && heat?._fullLayout && bars.data?.some(t => t.type === 'bar') &&
      heat.data?.some(t => t.type === 'heatmap') && !document.documentElement.classList.contains('shiny-busy');
  });
  assert.equal(await page.$('#spatial-enrich_clip'), null, 'Obsolete object exclusion still present');
  const choices = await page.$$eval('input[name="spatial-enrich_axis_range"]', inputs => inputs.map(x => ({ value: x.value, checked: x.checked })));
  assert.deepEqual(choices, [{ value: 'bulk', checked: true }, { value: 'all', checked: false }]);
  const cohort = await checkScenario('peri_enrichment_cohort');
  if (reference) {
    const column = cohort.heat.x.indexOf('T Cell CD8a⁺');
    assert(column >= 0, 'Reference phenotype missing');
    let matched = 0;
    for (let i = 0; i < cohort.heat.y.length; i++) {
      const id = String(cohort.heat.y[i]).split(' ')[0];
      if (reference[id] === undefined) continue;
      assert(Math.abs(cohort.heat.z[i][column] - reference[id]) < 1e-6, 'Donor mean differs from unfiltered source');
      matched++;
    }
    assert.equal(matched, Object.keys(reference).length);
    report.checks.push({ name: 'unfiltered_reference_donor_means', passed: true, donors: matched });
  }
  await setRadio('spatial-enrich_baseline', 'own');
  await checkScenario('peri_enrichment_own_tissue');
  await setRadio('spatial-enrich_region', 'core');
  await checkScenario('core_proportion');
  await setRadio('spatial-enrich_region', 'peri_prop');
  await checkScenario('peri_proportion');
  await setRadio('spatial-enrich_region', 'peri');
  for (const width of [1200, 1600]) {
    await page.setViewport({ width, height: 1100, deviceScaleFactor: 1 });
    const controls = await page.evaluate(() => {
      const ids = ['spatial-enrich_region', 'spatial-enrich_baseline', 'spatial-enrich_stat', 'spatial-enrich_axis_range'];
      const card = document.getElementById(ids[0]).closest('.card').getBoundingClientRect();
      return ids.map(id => {
        const r = document.getElementById(id).getBoundingClientRect();
        return { x: r.x, y: r.y, right: r.right, bottom: r.bottom,
          fits: r.width > 0 && r.left >= card.left && r.right <= card.right + 1 };
      });
    });
    assert(controls.every(control => control.fits), 'Enrichment controls exceed the card');
    for (let i = 0; i < controls.length; i++) for (let j = i + 1; j < controls.length; j++) {
      const a = controls[i], b = controls[j];
      assert(a.right <= b.x + 1 || b.right <= a.x + 1 || a.bottom <= b.y + 1 || b.bottom <= a.y + 1,
        'Enrichment controls overlap');
    }
    report.checks.push({ name: `controls_fit_${width}`, passed: true });
  }
  assert.equal(report.pageErrors, 0);
  assert.equal(report.externalRequestsBlocked, 0);
  report.status = 'passed';
} catch (error) {
  report.status = 'failed'; report.failure = String(error); process.exitCode = 1;
} finally {
  report.finishedAt = new Date().toISOString(); await browser?.close(); await save();
}
console.log(JSON.stringify({ status: report.status, checks: report.checks, output: options.output, failure: report.failure || null }));

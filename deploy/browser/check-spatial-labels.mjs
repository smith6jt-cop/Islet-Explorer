#!/usr/bin/env node
// Check rendered donor labels, row headings and percentage ticks on a specified loopback candidate.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';
const options = {};
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  assert(['url', 'output', 'shots', 'chrome'].includes(key) && process.argv[i + 1], 'Invalid argument');
  options[key] = process.argv[i + 1];
}
for (const key of ['url', 'output', 'shots']) assert(options[key], `--${key} required`);
const url = new URL(options.url);
assert(url.protocol === 'http:' && ['localhost', '127.0.0.1'].includes(url.hostname) &&
  !url.username && !url.password && !url.search && !url.hash, 'Unadorned loopback URL required');
await fs.mkdir(path.dirname(path.resolve(options.output)), { recursive: true });
await fs.mkdir(options.shots, { recursive: true });
const report = { status: 'running', target: url.href, checks: [], page_errors: 0,
  limitations: ['Bind this report to the separately inspected candidate container and artifact manifests.'] };
const save = () => fs.writeFile(options.output, JSON.stringify(report, null, 2) + '\n', { mode: 0o600 });
const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
let browser, page;
async function idle(id) {
  await page.waitForFunction(id => document.getElementById(id)?._fullLayout &&
    !document.getElementById(id).classList.contains('recalculating') &&
    !document.documentElement.classList.contains('shiny-busy'), {}, id);
}
async function radio(name, value, output) {
  await page.$eval(`input[name="${name}"][value="${value}"]`, el => el.click());
  await pause(600); await idle(output);
  assert.equal(await page.evaluate(name => Shiny.shinyapp.$inputValues[name], name), value);
}
async function check(name, run) {
  const detail = await run(); report.checks.push({ name, passed: true, detail }); await save();
}
async function heat() {
  const result = await page.$eval('#spatial-enrichment_heatmap', el => {
    const trace = el.data.find(t => t.type === 'heatmap');
    return { donors: trace.y, labels: Array.from(el.querySelectorAll('.ytick text')).map(t => t.textContent),
      boxes: Array.from(el.querySelectorAll('.ytick text')).map(t => {
        const r = t.getBoundingClientRect(); return { top: r.top, bottom: r.bottom };
      }) };
  });
  assert(result.donors.length > 0);
  assert.deepEqual([...result.labels].sort(), [...result.donors].sort(), 'Every donor must have a label');
  const boxes = result.boxes.sort((a, b) => a.top - b.top);
  assert(boxes.slice(1).every((b, i) => b.top >= boxes[i].bottom), 'Donor labels overlap');
  return { donors: result.donors.length, labels: result.labels.length };
}
async function composition(norm) {
  // The dynamic plot container is rebound in a second Shiny flush. An idle
  // document can briefly still contain the previous density plot or norm.
  await page.waitForFunction(norm => {
    const el = document.getElementById('spatial-immune_distance_kde');
    const footer = el?.layout?.annotations?.find(a => a.text?.startsWith('Signed distance'))?.text || '';
    return el?.data?.some(t => t.type === 'bar') &&
      footer.endsWith(norm === 'prop' ? 'y = proportion of the cells at that distance' : 'y = cells');
  }, {}, norm);
  await idle('spatial-immune_distance_kde');
  const result = await page.$eval('#spatial-immune_distance_kde', el => {
    const rect = e => { const r = e.getBoundingClientRect(); return { top: r.top, bottom: r.bottom }; };
    const headings = Array.from(el.querySelectorAll('.annotation-text')).filter(e =>
      !['ND', 'Aab+', 'T1D'].includes(e.textContent) && !e.textContent.startsWith('Signed distance'));
    const ticks = {};
    for (const e of el.querySelectorAll('g[class]')) if (/^y\d*tick$/.test(e.getAttribute('class'))) {
      const key = e.getAttribute('class'); (ticks[key] ||= []).push(e.textContent);
    }
    const axes = Object.entries(el._fullLayout).filter(([k]) => /^yaxis\d*$/.test(k))
      .map(([k, a]) => ({ name: k, domain: a.domain, range: a.range, title: a.title.text }));
    return { headings: headings.map(e => ({ text: e.textContent, ...rect(e) })), ticks, axes,
      selected: document.querySelectorAll('input[name="spatial-dist_comp_items"]:checked').length };
  });
  assert.equal(result.headings.length, result.selected, 'One horizontal heading per selected row');
  const boxes = result.headings.sort((a, b) => a.top - b.top);
  assert(boxes.slice(1).every((b, i) => b.top > boxes[i].bottom), 'Row headings overlap');
  assert(result.axes.every(a => !a.title), 'Long rotated y titles must not recur');
  assert(Object.values(result.ticks).every(t => t.length === new Set(t).size), 'Repeated y tick labels');
  for (const a of result.axes) for (const b of result.axes)
    if (JSON.stringify(a.domain) === JSON.stringify(b.domain)) assert.deepEqual(a.range, b.range);
  return { rows: result.headings.length, labeled_axes: Object.keys(result.ticks).length,
    lymphatic_ticks: result.ticks[Object.keys(result.ticks)[boxes.findIndex(b => b.text.includes('Lymphatic'))]] };
}
async function shot(name, id) {
  await page.$eval('#' + id, el => {
    const scroller = document.getElementById('main-content');
    const bottom = document.querySelector('.tab-pane.active .tab-controls').getBoundingClientRect().bottom;
    scroller.scrollTop += el.getBoundingClientRect().top - bottom - 12;
  });
  await (await page.$('#' + id)).screenshot({ path: path.join(options.shots, name + '.png') });
}
try {
  browser = await puppeteer.launch({ executablePath: options.chrome || '/usr/bin/google-chrome', headless: true,
    args: ['--no-sandbox', '--disable-dev-shm-usage'] });
  page = await browser.newPage(); page.setDefaultTimeout(180000);
  page.on('pageerror', () => report.page_errors++);
  await page.setViewport({ width: 1600, height: 2200, deviceScaleFactor: 1 });
  await page.goto(url.href); await idle('plot-plt');
  await page.click('#tabs a[data-value="Spatial"]'); await idle('spatial-enrichment_heatmap');
  const donorChoices = await page.$eval('#spatial-donor', el => Object.keys(el.selectize.options).length);
  const first = await heat(); assert.equal(first.donors, donorChoices, 'All cohort donors must appear');
  for (const width of [1200, 1600]) {
    await page.setViewport({ width, height: 2200, deviceScaleFactor: 1 }); await pause(600);
    for (const region of ['peri', 'core', 'peri_prop']) {
      await radio('spatial-enrich_region', region, 'spatial-enrichment_heatmap');
      await check(`heatmap_${region}_${width}`, heat);
    }
    await shot(`heatmap-${width}`, 'spatial-enrichment_heatmap');
    await radio('spatial-dist_ptype', 'composition', 'spatial-immune_distance_kde');
    for (const norm of ['prop', 'count']) {
      await radio('spatial-dist_norm', norm, 'spatial-immune_distance_kde');
      await check(`composition_${norm}_${width}`, () => composition(norm));
      if (norm === 'prop') await shot(`composition-${width}`, 'spatial-immune_distance_kde');
    }
  }
  assert.equal(report.page_errors, 0); report.status = 'passed';
} catch (error) { report.status = 'failed'; report.error = String(error); process.exitCode = 1; }
finally { await browser?.close(); await save(); }
console.log(JSON.stringify({ status: report.status, checks: report.checks.length, error: report.error }));

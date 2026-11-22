#!/usr/bin/env node
// Sparse T-cell regression: entering a measure comparison restores individual object points.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import puppeteer from 'puppeteer-core';
const options = {};
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  assert(['url', 'output', 'reference', 'chrome', 'scope'].includes(key) && process.argv[i + 1]);
  options[key] = process.argv[i + 1];
}
assert(options.url && options.output);
const url = new URL(options.url);
assert(url.protocol === 'http:' && ['localhost', '127.0.0.1'].includes(url.hostname) &&
  !url.username && !url.password && !url.search && !url.hash);
const scope = options.scope || 'points';
assert(['points', 'combined'].includes(scope));
assert(scope !== 'combined' || options.reference, 'Combined scientific checks require an independent reference');
const referenceBytes = options.reference ? await fs.readFile(options.reference) : null;
const reference = referenceBytes ? JSON.parse(referenceBytes.toString('utf8')) : null;
const lookup = reference && new Map(reference.data.map(r => [r['Case ID'] + '|' + r.islet_key, r]));
const report = { status: 'running', scope, checks: [], page_errors: 0, target: url.href,
  reference_sha256: referenceBytes && crypto.createHash('sha256').update(referenceBytes).digest('hex'),
  reference_checked: Boolean(reference), limitations: ['Bind to the inspected app/runtime/data combination.'] };
const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
let browser, page;
function expectedX(key, col) {
  const r = lookup.get(key); assert(r, 'Unknown donor/object key');
  if (!col.startsWith('combined_prop_')) return r[col];
  const pair = reference.pairs.find(p => p.stem === col.slice('combined_prop_'.length));
  assert(pair, 'Combined pair missing from reference');
  const core = Number.isFinite(r.total_cells_core) ? r.total_cells_core : 0;
  const peri = Number.isFinite(r.total_cells_peri) ? r.total_cells_peri : 0;
  // Reference proportions are already in percent. Weight by each compartment's cells.
  return ((r[pair.core_col] || 0) * core + (r[pair.peri_col] || 0) * peri) / (core + peri);
}
function csvRows(text) {
  const input = text.split(/\r?\n/).filter(line => line && !line.startsWith('#')).join('\n');
  const rows = []; let row = [], field = '', quoted = false;
  for (let i = 0; i < input.length; i++) {
    const c = input[i];
    if (c === '"') {
      if (quoted && input[i + 1] === '"') { field += '"'; i++; } else quoted = !quoted;
    } else if (c === ',' && !quoted) { row.push(field); field = ''; }
    else if (c === '\n' && !quoted) { row.push(field); rows.push(row); row = []; field = ''; }
    else field += c;
  }
  row.push(field); rows.push(row);
  const header = rows.shift(); return rows.map(r => Object.fromEntries(header.map((h, i) => [h, r[i]])));
}
async function idle() {
  await pause(300);
  await page.waitForFunction(() => !document.documentElement.classList.contains('shiny-busy'));
}
async function select(id, value) {
  await page.waitForFunction((id, value) => document.getElementById(id)?.selectize?.options[value], {}, id, value);
  await page.$eval('#' + id, (el, value) => {
    if (!el.selectize.options[value]) throw new Error('Choice unavailable');
    el.selectize.setValue(value);
  }, value);
  await idle();
}
try {
  browser = await puppeteer.launch({ executablePath: options.chrome || '/usr/bin/google-chrome', headless: true,
    args: ['--no-sandbox', '--disable-dev-shm-usage'] });
  page = await browser.newPage(); page.setDefaultTimeout(180000);
  page.on('pageerror', () => report.page_errors++);
  await page.setViewport({ width: 1600, height: 1200 });
  await page.goto(url.href);
  await page.waitForFunction(() => document.getElementById('plot-dist')?.data?.length &&
    !document.documentElement.classList.contains('shiny-busy'));
  await select('plot-mode', 'Microenvironment');
  await select('plot-which', 'Capillary|islet_core');
  await page.$eval('input[name="plot-axis_range"][value="all"]', el => el.click()); await idle();
  for (const side of ['right', 'left']) {
    const output = side === 'right' ? 'plot-dist' : 'plot-plt';
    await select(`plot-${side}_type`, 'status');
    await page.$eval('#plot-show_points', el => { if (el.checked) el.click(); }); await idle();
    assert.equal(await page.$eval('#plot-show_points', el => el.checked), false);
    await select(`plot-${side}_type`, 'xy');
    await page.waitForFunction(output => document.getElementById('plot-show_points').checked &&
      Shiny.shinyapp.$inputValues['plot-show_points'] === true &&
      document.getElementById(output)?.data?.some(t => t.key?.length), {}, output);
    report.checks.push({ name: `${side}_comparison_restores_object_points`, passed: true });
    const columns = scope === 'combined' ? await page.$eval(`#plot-${side}_x_measure`, el =>
      Object.keys(el.selectize.options).filter(k => k.startsWith('combined_prop_'))) :
      ['prop_All T cells', 'prop_CD8a Tcell', 'peri_prop_CD4_Tcell', 'peri_prop_CD8a_Tcell'];
    assert(columns.length > 0, 'No comparison choices');
    for (const col of columns) {
      await select(`plot-${side}_x_measure`, col);
      const label = await page.$eval(`#plot-${side}_x_measure`, (el, col) =>
        el.selectize.options[col][el.selectize.settings.labelField], col);
      await page.waitForFunction((id, label) => document.getElementById(id)?.layout?.xaxis?.title?.text === label &&
        !document.getElementById(id).classList.contains('recalculating') &&
        !document.documentElement.classList.contains('shiny-busy'), {}, output, label);
      const data = await page.$eval('#' + output, el => el.data.filter(t => t.key?.length)
        .map(t => ({ key: t.key, x: t.x })));
      let points = 0, nonzero = 0;
      for (const trace of data) for (let i = 0; i < trace.x.length; i++) {
        points++; if (trace.x[i] > 0) nonzero++;
        if (lookup) assert(Math.abs(trace.x[i] - expectedX(trace.key[i], col)) < 1e-8, 'Source x differs');
      }
      assert(points > 0 && nonzero > 0, 'Sparse T-cell values must not collapse to zero');
      let csv_rows;
      if (scope === 'combined') {
        assert(label.includes('core+peri %'), 'Combined region/units missing from label');
        const csv = await page.$eval(side === 'right' ? '#plot-dl_distribution_csv' : '#plot-dl_summary',
          async el => { const response = await fetch(el.href); if (!response.ok) throw new Error('CSV download failed');
            return response.text(); });
        assert(csv.includes(label), 'CSV must identify its combined x measure');
        const rows = csvRows(csv); csv_rows = rows.length;
        assert.equal(csv_rows, points, 'CSV/object point counts differ');
        for (const r of rows) assert(Math.abs(Number(r.x_value) - expectedX(r['Case ID'] + '|' + r.islet_key, col)) < 1e-8,
          'Combined CSV value differs from independent source calculation');
      }
      report.checks.push({ name: `${side}_${col}`, passed: true, points, nonzero, csv_rows });
    }
    // An explicit hide remains available after entering comparison mode.
    await page.$eval('#plot-show_points', el => el.click());
    await page.waitForFunction(id => {
      const el = document.getElementById(id);
      return !el.data.some(t => t.key?.length) && el.layout.title.text.includes('individual points hidden');
    }, {}, output);
    const medians = await page.$eval('#' + output, el => el.data.find(t => t.marker?.symbol === 'diamond').x);
    report.checks.push({ name: `${side}_explicit_hide_is_labeled`, passed: true,
      donor_medians: medians.length, zero_medians: medians.filter(x => x === 0).length });
  }
  assert.equal(report.page_errors, 0); report.status = 'passed';
} catch (error) { report.status = 'failed'; report.error = String(error); process.exitCode = 1; }
finally {
  await browser?.close(); await fs.mkdir(path.dirname(path.resolve(options.output)), { recursive: true });
  await fs.writeFile(options.output, JSON.stringify(report, null, 2) + '\n', { mode: 0o600 });
}
console.log(JSON.stringify({ status: report.status, checks: report.checks.length, error: report.error }));

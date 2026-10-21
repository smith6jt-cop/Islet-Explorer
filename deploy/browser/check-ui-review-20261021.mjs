// Focused check of the 2026-10-21 UI review: tab order, palette labels, grouped control panels, the frozen header /
// tab bar / control panels, the Spatial UMAP colourings and the linked Spatial UMAPs. Use an isolated local worker:
//   node deploy/browser/check-ui-review-20261021.mjs --url http://127.0.0.1:3841/islet-explorer-test/ \
//     --output /tmp/ui-review-check.json [--shots <dir>]
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const arg = name => process.argv.includes(name) ? process.argv[process.argv.indexOf(name) + 1] : null;
const url = arg('--url') || 'http://127.0.0.1:3841/islet-explorer-test/';
const output = arg('--output') || '/tmp/ui-review-check.json';
const shots = arg('--shots');
const report = {checks: {}, errors: []};
const browser = await puppeteer.launch({executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-dev-shm-usage']});
const page = await browser.newPage();
await page.setViewport({width: 1680, height: 1050});
page.setDefaultTimeout(180000);
page.on('pageerror', error => report.errors.push(String(error)));
let failed = 0;
async function check(name, fn) {
  let result;
  try { result = await fn(); } catch (error) { result = {passed: false, error: String(error)}; }
  const passed = result === true || result?.passed === true;
  if (!passed) failed += 1;
  report.checks[name] = result;
  console.log(passed ? 'PASS' : 'FAIL', name, typeof result === 'object' ? JSON.stringify(result) : '');
}
async function shot(name) { if (shots) await page.screenshot({path: path.join(shots, name + '.png')}); }
const idle = () => page.waitForFunction(() => !document.documentElement.classList.contains('shiny-busy'));
async function plotted(id) {
  await page.waitForFunction(id => {
    const gd = document.getElementById(id);
    return gd && !gd.classList.contains('recalculating') && gd.data?.some(t => (t.x?.length || 0) > 1);
  }, {}, id);
}
async function tab(name) {
  await page.click(`#tabs a[data-value="${name}"]`);
  await page.waitForFunction(name => document.querySelector('#tabs li.active a')?.dataset.value === name, {}, name);
}
async function choose(id, value) {
  await page.$eval('#' + id, (el, value) => {
    if (el.selectize) el.selectize.setValue(value);
    else { el.value = value; el.dispatchEvent(new Event('change', {bubbles: true})); }
  }, value);
}

await page.goto(url, {waitUntil: 'networkidle2'});
await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
await plotted('plot-plt');
await idle();
await shot('plot_top');

await check('tab_order_viewer_last', () => page.evaluate(() => {
  const tabs = [...document.querySelectorAll('#tabs > li > a')].map(a => a.dataset.value);
  return {passed: JSON.stringify(tabs) === JSON.stringify(['Plot', 'Trajectory', 'Spatial', 'Viewer']), tabs};
}));
await check('palette_without_default_label', () => page.evaluate(() => {
  const el = document.getElementById('sidebar_donor_palette');
  const names = el.selectize ? Object.values(el.selectize.options).map(o => o.label) : [...el.options].map(o => o.text);
  const value = el.selectize ? el.selectize.getValue() : el.value;
  return {passed: !names.some(n => /default/i.test(n)) && value === 'Bright' && names[0] === 'Bright', names, value};
}));
await check('plot_data_controls_grid', () => page.evaluate(() => {
  const grid = document.querySelector('#tabs ~ .tab-content .tab-pane.active .control-panel[open] .control-grid');
  const groups = [...(grid?.children || [])];
  const tops = groups.map(g => Math.round(g.getBoundingClientRect().top));
  const widths = groups.map(g => Math.round(g.getBoundingClientRect().width));
  const titles = groups.map(g => g.querySelector('.control-group-title')?.textContent);
  const inputsFill = [...grid.querySelectorAll('.selectize-control, .irs')].every(el => {
    const g = el.closest('.control-group').getBoundingClientRect(), r = el.getBoundingClientRect();
    return r.width >= g.width - 2;
  });
  return {passed: groups.length === 4 && new Set(tops).size === 1 && Math.max(...widths) - Math.min(...widths) <= 2 &&
    inputsFill, titles, widths, height: Math.round(grid.getBoundingClientRect().height)};
}));
await check('logo_dark_and_light', () => page.evaluate(() => {
  const dark = document.querySelector('.brand-logo-dark'), light = document.querySelector('.brand-logo-light');
  return {passed: document.documentElement.dataset.theme === 'dark' && dark.offsetParent !== null &&
    light.offsetParent === null && /logo-header-dark/.test(dark.src)};
}));

// Scroll the main panel: header, tab bar and the tab's controls stay put, cards move under them.
async function frozen(label) {
  return page.evaluate(label => {
    const scroller = document.getElementById('main-content');
    const pane = document.querySelector('.tab-content > .tab-pane.active');
    const controls = pane.querySelector('.tab-controls');
    const parts = {header: document.querySelector('.app-header'), tabs: document.getElementById('tabs'), controls};
    const before = Object.fromEntries(Object.entries(parts).map(([k, el]) => [k, Math.round(el.getBoundingClientRect().top)]));
    const card = pane.querySelector('.card');
    const cardBefore = card.getBoundingClientRect().top;
    scroller.scrollTop = 700;
    const after = Object.fromEntries(Object.entries(parts).map(([k, el]) => [k, Math.round(el.getBoundingClientRect().top)]));
    const cardMoved = card.getBoundingClientRect().top < cardBefore - 300;
    const stacked = after.tabs >= after.header + parts.header.getBoundingClientRect().height - 1 &&
      after.controls >= after.tabs + parts.tabs.getBoundingClientRect().height - 1;
    // The frozen header is at the top of the scroll panel, and nothing scrolled under is drawn over it.
    const top = document.elementFromPoint(400, after.controls + 10);
    const result = {passed: scroller.scrollTop > 600 && cardMoved && stacked &&
      Math.abs(after.header - scroller.getBoundingClientRect().top) <= 1 && controls.contains(top),
      label, before, after, scrollTop: scroller.scrollTop};
    return result;
  }, label);
}
await check('plot_frozen_controls', () => frozen('Plot'));
await shot('plot_scrolled');
await page.evaluate(() => { document.getElementById('main-content').scrollTop = 0; });

await tab('Trajectory');
await plotted('traj-traj_scatter');
await idle();
await check('trajectory_frozen_controls', () => frozen('Trajectory'));
await page.evaluate(() => { document.getElementById('main-content').scrollTop = 0; });

await tab('Spatial');
await plotted('spatial-leiden_umap');
await plotted('spatial-spatial_umap_immune');
await idle();
await shot('spatial_top');
await check('spatial_frozen_controls', () => frozen('Spatial'));
await page.evaluate(() => {
  const card = document.getElementById('spatial-leiden_umap').closest('.card');
  const scroller = document.getElementById('main-content');
  scroller.scrollTop += card.getBoundingClientRect().top - 420;
});
await shot('spatial_umaps');

await check('spatial_umap_colour_choices', () => page.evaluate(() => {
  const el = document.getElementById('spatial-umap_color_left');
  const values = Object.keys(el.selectize.options);
  const want = ['leiden', 'donor_status', 'n_aab', 'donor', 'age', 'age_group', 'gender', 'BMI', 'HbA1c_percent',
    'C_peptide_ng_mL', 't1d_duration_years', 'immune_frac_peri', 'immune_frac_core', 'islet_diam_um', 'total_cells_core'];
  return {passed: want.every(v => values.includes(v)), values,
    left: el.selectize.getValue(), right: document.getElementById('spatial-umap_color_right').selectize.getValue()};
}));

const colourings = ['donor_status', 'donor', 'age', 'age_group', 'gender', 'BMI', 'HbA1c_percent', 'C_peptide_ng_mL',
  't1d_duration_years', 'n_aab', 'immune_frac_core', 'islet_diam_um', 'total_cells_core', 'leiden'];
const drawn = {};
for (const value of colourings) {
  const before = await page.$eval('#spatial-leiden_umap', gd => gd.layout?.title?.text);
  await choose('spatial-umap_color_left', value);
  await page.waitForFunction(old => {
    const gd = document.getElementById('spatial-leiden_umap');
    return gd.layout?.title?.text !== old && !gd.classList.contains('recalculating');
  }, {timeout: 30000}, before).catch(() => null);
  await idle();
  drawn[value] = await page.$eval('#spatial-leiden_umap', gd => {
    const pts = gd.data.find(t => t.meta?.trajectory_role === 'islet');
    const colours = pts?.marker?.color;
    const legend = gd.data.filter(t => t.showlegend && t.name).map(t => t.name);
    return {title: gd.layout.title.text, points: pts?.x?.length || 0, legend: legend.slice(0, 6), n_legend: legend.length,
      colourbar: Boolean(pts?.marker?.colorbar), grey: Array.isArray(colours) ?
        colours.filter(c => c === '#9d9d9d' || c === null).length : null};
  });
  if (['donor_status', 't1d_duration_years', 'donor'].includes(value)) await shot('spatial_umap_' + value);
}
await check('spatial_umap_colourings_draw', () => ({
  passed: Object.values(drawn).every(d => d.points > 100 && (d.colourbar || d.n_legend > 0)), drawn}));

// Linked UMAPs: hover rings the object in both plots; a click opens its segmentation map below.
const ids = ['spatial-leiden_umap', 'spatial-spatial_umap_immune'];
await page.waitForFunction(ids => ids.every(id => document.getElementById(id)?.querySelector('.trajectory-link-highlight')), {}, ids);
const key = await page.$eval('#spatial-leiden_umap', gd => {
  const t = gd.data.find(t => t.meta?.trajectory_role === 'islet');
  return t.key[Math.floor(t.key.length / 3)];
});
async function emit(id, event, key) {
  await page.evaluate((id, event, key) => {
    const gd = document.getElementById(id);
    const curveNumber = gd.data.findIndex(t => t.meta?.trajectory_role === 'islet' && t.key.includes(key));
    const data = gd.data[curveNumber];
    gd.emit(event, {points: [{curveNumber, pointNumber: data.key.indexOf(key), data}]});
  }, id, event, key);
}
await emit('spatial-leiden_umap', 'plotly_hover', key);
await check('spatial_hover_rings_both', () => page.waitForFunction((ids, key) => ids.every(id =>
  document.getElementById(id).querySelector('.trajectory-link-highlight')?.getAttribute('data-linked-key') === key),
  {timeout: 5000}, ids, key).then(() => ({passed: true, key})));
await emit('spatial-leiden_umap', 'plotly_unhover', key);
await emit('spatial-spatial_umap_immune', 'plotly_click', key);
await check('spatial_click_opens_segmentation', async () => {
  await page.waitForSelector('#spatial-segmentation_viewer_panel #islet_segmentation_view img', {timeout: 60000});
  await idle();
  return page.evaluate(key => {
    const title = document.querySelector('#spatial-segmentation_viewer_panel h3')?.textContent || '';
    const islet = key.split('|')[1];
    const pinned = ['spatial-leiden_umap', 'spatial-spatial_umap_immune'].every(id =>
      document.getElementById(id).querySelector('.trajectory-link-highlight')?.getAttribute('data-linked-key') === key);
    return {passed: title.startsWith(islet) && pinned, title, key};
  }, key);
});
await page.evaluate(() => {
  const panel = document.getElementById('spatial-segmentation_viewer_panel');
  document.getElementById('main-content').scrollTop += panel.getBoundingClientRect().top - 420;
});
await shot('spatial_segmentation');

// The selection is shared: the Plot tab rings the same object; Clear selection releases it everywhere.
await tab('Plot');
await plotted('plot-plt');
await check('selection_shared_with_plot_tab', () => page.waitForFunction(key =>
  document.getElementById('plot-plt').querySelector('.trajectory-link-highlight')?.getAttribute('data-linked-key') === key,
  {timeout: 20000}, key).then(() => ({passed: true})));
await tab('Spatial');
await page.click('#spatial-clear_link_selection');
await check('clear_selection_closes_panel', () => page.waitForFunction(ids =>
  !document.querySelector('#spatial-segmentation_viewer_panel .card') &&
  ids.every(id => !document.getElementById(id).querySelector('.trajectory-link-highlight')?.getAttribute('data-linked-key')),
  {timeout: 20000}, ids).then(() => ({passed: true})));

await tab('Viewer');
await check('viewer_frozen_controls', () => page.evaluate(() => {
  const c = document.querySelector('.tab-content > .tab-pane.active .tab-controls');
  return {passed: Boolean(c) && getComputedStyle(c).position === 'sticky'};
}));

report.failed = failed;
await fs.writeFile(output, JSON.stringify(report, null, 2));
await browser.close();
if (report.errors.length) console.log('page errors:', report.errors);
process.exit(failed || report.errors.length ? 1 : 0);

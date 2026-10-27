// Focused check of the 2026-10-21 UI review: tab order, palette labels, grouped control panels, the frozen header /
// tab bar / control panels, the Spatial UMAP colourings and the linked Spatial UMAPs. Use an isolated local worker:
//   node deploy/browser/check-ui-review-20261021.mjs --url http://127.0.0.1:3840/islet-explorer/ \
//     --output /tmp/ui-review-check.json [--shots <dir>]
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

// Display names of phenotype labels (app/shiny_app/R/00_globals.R phenotype_text(); owner, 2026-10-27): broad name
// first, every marker + and - a superscript. Values (data labels) are unchanged.
const PHENOTYPE_DISPLAY = {'Beta cell': 'Endocrine Beta', 'Alpha cell': 'Endocrine Alpha', 'Delta cell': 'Endocrine Delta',
  'CHGA+ hormone-negative': 'Endocrine CHGA+ hormone-negative', 'CD8a Tcell': 'T Cell CD8a+', 'CD4 Tcell': 'T Cell CD4+',
  'T cell': 'T Cell', 'B cell': 'B Cell', 'CD8 T cell': 'T Cell CD8+', 'CD4 T cell': 'T Cell CD4+',
  'All endocrine cells': 'Endocrine (all)', 'All T cells': 'T Cell (all)', 'All macrophages': 'Macrophage (all)'};
const phenotypeText = label => (PHENOTYPE_DISPLAY[label] ?? label)
  .replace(/^(CD68[+-]CD163[+-])Macrophage_(LGALS3[+-]CD4[+-])$/, 'Macrophage $1$2')
  .replace(/\+/g, '⁺').replace(/([A-Za-z0-9])-(?=$|[ A-Z])/g, '$1⁻');
const arg = name => process.argv.includes(name) ? process.argv[process.argv.indexOf(name) + 1] : null;
const url = arg('--url') || 'http://127.0.0.1:3840/islet-explorer/';
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
  // The tab list sits in the header row since 2026-10-26, so the tab content is found on its own.
  const grid = document.querySelector('.tab-content .tab-pane.active .control-panel[open] .control-grid');
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

// Plot types (2026-10-21, second round): every type draws in both cards, its CSV downloads, and the object types
// stay linked plots.
const types = ['diameter', 'status', 'donor_medians', 'ecdf', 'donor_curves', 'xy', 'clinical', 'density'];
const typeResults = {};
for (const [card, id] of [['left', 'plot-plt'], ['right', 'plot-dist']]) {
  for (const type of types) {
    const before = await page.$eval('#' + id, gd => gd.layout?.title?.text || '');
    await choose(`plot-${card}_type`, type);
    await page.waitForFunction((id, old) => {
      const gd = document.getElementById(id);
      return gd && (gd.layout?.title?.text || '') !== old && !gd.classList.contains('recalculating');
    }, {timeout: 60000}, id, before).catch(() => null);
    await idle();
    const drawn = await page.$eval('#' + id, gd => ({
      title: (gd.layout?.title?.text || '').replace(/<[^>]+>/g, ' ').slice(0, 140),
      traces: gd.data.filter(t => (t.x?.length || 0) > 1 || (t.y?.length || 0) > 1).length,
      keyed: gd.data.filter(t => t.meta?.trajectory_role === 'islet').length }));
    const heading = await page.$eval(`#plot-${card}_title`, el => el.textContent);
    const csv = await page.$eval(card === 'left' ? '#plot-dl_summary' : '#plot-dl_distribution_csv',
      async el => { const t = await (await fetch(el.href)).text(); const lines = t.split('\n').filter(Boolean);
        return {header: lines.filter(l => l.startsWith('#')).length, rows: lines.filter(l => !l.startsWith('#')).length - 1,
                first: lines.find(l => !l.startsWith('#'))?.slice(0, 100)}; });
    typeResults[`${card}:${type}`] = {...drawn, heading, csv};
    if (card === 'right') await shot(`plot_type_${type}`);
  }
}
await check('plot_types_draw_in_both_cards', () => {
  const linked = ['diameter', 'status', 'xy', 'clinical'];
  const bad = Object.entries(typeResults).filter(([k, r]) => {
    const type = k.split(':')[1];
    return !(r.traces > 0 && r.csv.rows > 0 && r.csv.header > 3 && r.heading.length > 3 &&
             (!linked.includes(type) || r.keyed > 0) && !/Loading|No data/.test(r.title));
  }).map(([k]) => k);
  return {passed: bad.length === 0, bad, results: typeResults};
});
await choose('plot-left_type', 'diameter');
await choose('plot-right_type', 'status');
await idle();

await check('composition_by_donor_card', async () => {
  await page.evaluate(() => {
    const el = document.getElementById('plot-composition_by_donor');
    document.getElementById('main-content').scrollTop += el.getBoundingClientRect().top - 500;
  });
  await plotted('plot-composition_by_donor');
  await shot('composition_by_donor');
  const fig = await page.$eval('#plot-composition_by_donor', gd => {
    const bars = gd.data.filter(t => t.type === 'bar');
    const donors = bars[0]?.x?.length || 0;
    const sums = Array.from({length: donors}, (_, i) => bars.reduce((a, t) => a + (t.y[i] || 0), 0));
    return {phenotypes: bars.map(t => t.name), donors, sums: sums.map(v => Math.round(v * 10) / 10)};
  });
  // Leave out (third round): the core view leaves out the four endocrine types by default, the band view Exocrine.
  const leftOut = await page.$$eval('input[name="plot-composition_leave_out"]:checked', els => els.map(e => e.value));
  const csv = await page.$eval('#plot-dl_composition_csv', async el => (await (await fetch(el.href)).text()).split('\n').length);
  await choose('plot-composition_region', 'peri');
  await page.waitForFunction(() => /peri-islet band/.test(document.getElementById('plot-composition_by_donor').layout.title.text));
  const periLeftOut = await page.$$eval('input[name="plot-composition_leave_out"]:checked', els => els.map(e => e.value));
  const endocrine = ['Beta cell', 'Alpha cell', 'Delta cell', 'CHGA+ hormone-negative'];
  // legend entries are display names; the Leave out values stay data labels
  return {passed: fig.donors >= 15 && fig.sums.every(v => Math.abs(v - 100) < 0.6) && fig.phenotypes.includes('No subtype') &&
    !fig.phenotypes.some(n => ['Endocrine', 'T cell', 'Immune', ...endocrine].map(phenotypeText).includes(n)) && csv > 100 &&
    endocrine.every(n => leftOut.includes(n)) && JSON.stringify(periLeftOut) === JSON.stringify(['Exocrine']),
    ...fig, csv_lines: csv, left_out: leftOut, peri_left_out: periLeftOut};
});
await page.evaluate(() => { document.getElementById('main-content').scrollTop = 0; });

// Third round: near-empty measures are not offered, and the default axis fits 99% of objects with the rest pinned.
await check('near_empty_measures_not_offered', () => page.evaluate(() => {
  const values = Object.keys(document.getElementById('plot-which').selectize.options);
  const gone = ['peri_prop_CHGAplus_hormone-negative', 'prop_CD4 Tcell', 'prop_B cell', 'prop_Exocrine'];
  return {passed: !gone.some(v => values.includes(v)) && values.includes('peri_prop_CD4_Tcell'), values};
}));
await check('axis_fits_the_bulk', async () => {
  await choose('plot-which', 'peri_prop_CD4_Tcell');
  await page.waitForFunction(() => /T Cell CD4/.test(document.getElementById('plot-plt').layout?.title?.text || ''));
  await idle();
  const read = () => page.$eval('#plot-plt', gd => ({range: gd.layout.yaxis.range, title: gd.layout.title.text,
    pinned: gd.data.filter(t => t.marker?.symbol === 'triangle-up-open' || (t.mode === 'markers' && t.marker?.symbol === 'triangle-up-open')).length}));
  const bulk = await read();
  await page.evaluate(() => { document.querySelectorAll('.tab-pane.active details.control-panel')[1].open = true; });
  await page.click('input[name="plot-axis_range"][value="all"]');
  await page.waitForFunction(() => !/off scale/.test(document.getElementById('plot-plt').layout?.title?.text || ''));
  await idle();
  const all = await read();
  await page.click('input[name="plot-axis_range"][value="bulk"]');
  await page.evaluate(() => { document.querySelectorAll('.tab-pane.active details.control-panel')[1].open = false; });
  await choose('plot-which', 'prop_Beta cell');
  await idle();
  return {passed: /off scale/.test(bulk.title) && bulk.range[1] < 20 && !/off scale/.test(all.title), bulk, all};
});

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
    // Since 2026-10-26 the tab bar is inside the header row, so the controls stack right under the header; the
    // older layout (tab bar under the header) is still accepted.
    const headerBottom = after.header + parts.header.getBoundingClientRect().height - 1;
    const stacked = parts.header.contains(parts.tabs) ? after.controls >= headerBottom :
      after.tabs >= headerBottom && after.controls >= after.tabs + parts.tabs.getBoundingClientRect().height - 1;
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

// The nPOD donor fields of v10 (2026-10-22) are among them.
const colourings = ['donor_status', 'donor', 'age', 'age_group', 'gender', 'BMI', 'HbA1c_percent', 'C_peptide_ng_mL',
  't1d_duration_years', 'n_aab', 'race_ethnicity', 'cause_of_death', 'recovery_type', 'ABO_group', 'ICU_days',
  'transit_time_min', 'pancreas_weight_g', 'RIN', 'peak_glucose_mg_dL', 'immune_frac_core', 'islet_diam_um',
  'total_cells_core', 'leiden'];
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

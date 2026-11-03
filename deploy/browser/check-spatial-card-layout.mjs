#!/usr/bin/env node
// Focused loopback Spatial UI acceptance. Runtime identities are caller-supplied;
// bind this result to the separately inspected image/container receipt.
import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const options = {};
const allowed = ['url', 'output', 'shots', 'source', 'runtime-kind', 'image', 'chrome', 'scope'];
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  assert(allowed.includes(key) && process.argv[i + 1], 'Unexpected or incomplete argument');
  options[key] = process.argv[i + 1];
}
for (const key of ['url', 'output', 'shots', 'source', 'runtime-kind']) assert(options[key], `--${key} required`);
const url = new URL(options.url);
assert(['127.0.0.1', 'localhost', '[::1]'].includes(url.hostname) && url.protocol === 'http:', 'Loopback HTTP required');
assert(!url.username && !url.password && !url.search && !url.hash, 'Unadorned preview URL required');
assert(/^[a-f0-9]{40}$/.test(options.source), 'Full source commit required');
assert(['source-overlay', 'native-image'].includes(options['runtime-kind']), 'Explicit runtime kind required');
if (options['runtime-kind'] === 'native-image') assert(/^sha256:[a-f0-9]{64}$/.test(options.image || ''), 'Native image digest required');
const scope = options.scope || 'all';
assert(['all', 'distance'].includes(scope), 'Scope must be all or distance');
const hash = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
const report = {
  schema: 'islet.spatial_card_layout/1', scope, status: 'running', started_utc: new Date().toISOString(),
  target: { url: url.href, source: options.source, runtime_kind: options['runtime-kind'], image: options.image || null,
    identity_basis: 'Caller supplied; verify against the companion runtime inspection.' },
  checks: [], cases: [], screenshots: [], page_errors: [], application_http_errors: [], external_requests_blocked: 0,
  limitations: ['No scientific retyping, public-access, security or reboot acceptance is established here.',
    'Captured screenshots require separate direct visual inspection for complete legend signs and axis labels.',
    'Download controls are checked for placement and href; downloaded file contents are not tested.'],
};
await fs.mkdir(path.dirname(path.resolve(options.output)), { recursive: true });
await fs.mkdir(path.resolve(options.shots), { recursive: true });
const save = () => fs.writeFile(options.output, JSON.stringify(report, null, 2) + '\n', { mode: 0o600 });
let browser, page;
async function check(name, run) {
  try { report.checks.push({ name, passed: true, detail: await run() }); }
  catch (error) { report.checks.push({ name, passed: false, error: error.stack || String(error) }); await save(); throw error; }
  await save();
}
async function idle() {
  await page.waitForFunction(() => !document.documentElement.classList.contains('shiny-busy') &&
    !document.getElementById('spatial-tissue_scatter')?.classList.contains('recalculating'));
}
async function tissueReady() {
  await page.waitForFunction(() => {
    const el = document.getElementById('spatial-tissue_scatter'), img = el?.querySelector('img');
    return img?.complete && img.naturalWidth > 0 && !el.classList.contains('shiny-output-error') &&
      Shiny.shinyapp.$values['spatial-tissue_scatter']?.coordmap?.panels?.length;
  });
  await idle();
  let prior = '', stable = 0;
  const started = Date.now(), history = [];
  for (let i = 0; i < 80 && stable < 3; i++) {
    await delay(250); await idle();
    const sample = await page.$eval('#spatial-tissue_scatter img', async el => {
      const rect = el.getBoundingClientRect();
      const sample = {src: el.src, naturalWidth: el.naturalWidth, naturalHeight: el.naturalHeight,
        width: rect.width, height: rect.height,
        domain: Shiny.shinyapp.$values['spatial-tissue_scatter']?.coordmap?.panels?.[0]?.domain,
        checkedMarkers: [...document.querySelectorAll('input[name="spatial-mac_markers_filter"]:checked')].map(x => x.value),
        markerInputs: Object.fromEntries(Object.entries(Shiny.shinyapp.$inputValues).filter(([key]) => /mac.*markers/.test(key))) };
      const canvas = document.createElement('canvas'); canvas.width = el.naturalWidth; canvas.height = el.naturalHeight;
      const context = canvas.getContext('2d'); context.drawImage(el, 0, 0);
      const bytes = context.getImageData(0, 0, canvas.width, canvas.height).data;
      const digest = await crypto.subtle.digest('SHA-256', bytes);
      sample.pixel_sha256 = [...new Uint8Array(digest)].map(x => x.toString(16).padStart(2, '0')).join('');
      canvas.width = 0; canvas.height = 0;
      return sample;
    });
    const {src, domain, checkedMarkers, markerInputs, pixel_sha256, ...dimensions} = sample;
    // Preserve the original strict criterion until the recorded image/dimension
    // evidence distinguishes genuine redraw feedback from a volatile image URL.
    const next = JSON.stringify({src, ...dimensions});
    history.push({elapsed_ms: Date.now() - started, source_sha256: hash(src), pixel_sha256,
      domain, checkedMarkers, markerInputs, ...dimensions});
    stable = next === prior ? stable + 1 : 0; prior = next;
  }
  (report.tissue_stability_attempts ||= []).push({phase: report.phase || null, stable, history});
  if (stable < 3) await save();
  assert(stable >= 3, 'Tissue image did not stabilize');
}
async function geometry() {
  return page.evaluate(() => {
    const el = document.getElementById('spatial-tissue_scatter'), img = el.querySelector('img');
    const map = Shiny.shinyapp.$values['spatial-tissue_scatter'].coordmap, panel = map.panels[0];
    const ir = img.getBoundingClientRect(), er = el.getBoundingClientRect();
    const sx = ir.width / map.dims.width, sy = ir.height / map.dims.height;
    return { panel: { x: ir.x + panel.range.left * sx, y: ir.y + panel.range.top * sy,
      width: (panel.range.right - panel.range.left) * sx, height: (panel.range.bottom - panel.range.top) * sy },
      output: { width: er.width, height: er.height }, domain: panel.domain,
      image: { width: ir.width, height: ir.height, naturalWidth: img.naturalWidth, naturalHeight: img.naturalHeight },
      stableOutputElement: el === window.__spatialReviewOutput,
      downloads: ['spatial-download_spatial', 'spatial-dl_tissue_scatter'].map(id => {
        const link = document.getElementById(id); return { id, below: link.getBoundingClientRect().top >= er.bottom,
          href: Boolean(link.getAttribute('href')) };
      }),
    };
  });
}
async function scrollTo(selector) {
  await page.$eval(selector, el => {
    const scroller = document.getElementById('main-content');
    const sticky = document.querySelector('.tab-pane.active .tab-controls');
    const offset = Math.max(document.querySelector('.app-header')?.getBoundingClientRect().bottom || 0,
      sticky?.getBoundingClientRect().bottom || 0) + 12;
    scroller.scrollTop += el.getBoundingClientRect().top - offset;
  });
}
async function screenshot(name, selector = null) {
  const destination = path.join(options.shots, `${name}.png`);
  if (selector) await (await page.$(selector)).screenshot({ path: destination });
  else await page.screenshot({ path: destination });
  report.screenshots.push({ name, path: path.resolve(destination), sha256: hash(await fs.readFile(destination)) });
}
async function captureTissue(scenario, width) {
  await page.setViewport({ width, height: 1300, deviceScaleFactor: 1 });
  await tissueReady(); await scrollTo('#spatial-tissue_scatter');
  const measurement = await geometry();
  assert(measurement.panel.width >= measurement.output.width - 150, 'Map does not fill available interior width');
  assert(Math.abs(measurement.image.width - measurement.output.width) <= 2, 'Raster does not span the card output');
  assert(measurement.stableOutputElement && measurement.downloads.every(x => x.below && x.href), 'Output or downloads changed');
  const raw = await page.$eval('#spatial-tissue_scatter img', img => img.src);
  assert(raw.startsWith('data:image/png;base64,'));
  const bytes = Buffer.from(raw.split(',')[1], 'base64');
  const destination = path.join(options.shots, `tissue-${scenario}-${width}-raw.png`);
  await fs.writeFile(destination, bytes);
  report.screenshots.push({ name: `tissue-${scenario}-${width}-raw`, path: path.resolve(destination), sha256: hash(bytes) });
  await screenshot(`tissue-${scenario}-${width}-card`, '#spatial-tissue_scatter');
  const expectedLegendColumns = Math.max(1, Math.min(4, Math.floor((measurement.output.width - 80) / 220)));
  assert(expectedLegendColumns >= 2, 'Desktop legend should use multiple columns');
  const detail = { scenario, viewport: width, expectedLegendColumns,
    legend_column_basis: 'Source layout calculation; confirm rendered columns in the captured PNG.', ...measurement };
  report.cases.push(detail); return detail;
}
async function dataPanel(open) {
  const selector = '.tab-pane.active .tab-controls details.control-panel';
  const handle = await page.$(selector);
  if (await handle.evaluate(el => el.open) !== open) await handle.evaluate(el => el.querySelector('summary').click());
  await page.waitForFunction((selector, open) => document.querySelector(selector)?.open === open, {}, selector, open);
}
async function plotSnapshot() {
  return page.evaluate(() => {
    const fields = trace => ({ name: trace.name, x: Array.from(trace.x || []), y: Array.from(trace.y || []),
      error: trace.error_y, dash: trace.line?.dash });
    const gradient = document.getElementById('spatial-immune_gradient')?.data || [];
    return { donor: document.getElementById('spatial-donor').value,
      gradientDonor: gradient.filter(t => /^Donor /.test(t.name || '')).map(fields),
      gradientCohort: gradient.filter(t => !/^Donor /.test(t.name || '')).map(fields),
      others: Object.fromEntries(['spatial-leiden_umap', 'spatial-spatial_umap_immune', 'spatial-cluster_composition',
        'spatial-insulitis_bars', 'spatial-beta_immune'].map(id => [id, (document.getElementById(id)?.data || []).map(fields)])),
      tissue: document.querySelector('#spatial-tissue_scatter img').src,
    };
  });
}
try {
  browser = await puppeteer.launch({ executablePath: options.chrome || '/usr/bin/google-chrome', headless: true,
    args: ['--disable-dev-shm-usage', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'] });
  report.chrome = await browser.version(); page = await browser.newPage(); page.setDefaultTimeout(180000);
  await page.setRequestInterception(true);
  page.on('request', request => {
    const target = new URL(request.url());
    if (['http:', 'https:'].includes(target.protocol) && target.origin !== url.origin) {
      report.external_requests_blocked++; void request.abort('blockedbyclient');
    } else void request.continue();
  });
  page.on('pageerror', error => report.page_errors.push(String(error)));
  page.on('response', response => {
    if (response.status() >= 400 && !new URL(response.url()).pathname.endsWith('/favicon.ico'))
      report.application_http_errors.push({ status: response.status(), url: response.url() });
  });
  await page.setViewport({ width: 1200, height: 1300, deviceScaleFactor: 1 });
  await page.goto(url.href, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
  await page.click('a[data-value="Spatial"]'); await tissueReady();
  await page.evaluate(() => { window.__spatialReviewOutput = document.getElementById('spatial-tissue_scatter'); });
  await check('visible_donor_outside_collapsed_data_controls', async () => {
    const state = await page.$eval('#spatial-donor', el => ({ outside: !el.closest('details'),
      heading: el.closest('.spatial-donor-selection').querySelector('h5').textContent,
      visible: el.closest('.shiny-input-container').getBoundingClientRect().height > 0,
      closed: [...document.querySelectorAll('.tab-pane.active .tab-controls details')].every(x => !x.open) }));
    assert(state.outside && state.visible && state.closed);
    assert.equal(state.heading, 'Donor For Tissue And Immune Gradient'); return state;
  });
  if (scope === 'all') for (const width of [1200, 1600]) await check(`default_tissue_${width}`, () => captureTissue('default', width));
  await check('four_marker_dense_selection', async () => {
    await dataPanel(true);
    for (const marker of ['CD68', 'CD163', 'LGALS3', 'CD4']) {
      report.phase = `select_marker_${marker}`; await save();
      const selector = `input[name="spatial-mac_markers_filter"][value="${marker}"]`;
      if (!await page.$eval(selector, el => el.checked)) await page.click(selector);
      await page.waitForFunction(selector => document.querySelector(selector)?.checked, {}, selector);
      await tissueReady();
    }
    await page.click('#spatial-pheno_filter .pheno-groups-overall .pheno-select-all');
    await tissueReady();
    const state = await page.$$eval('input[name="spatial-pheno_filter"]', inputs => inputs.map(x => ({ value: x.value, checked: x.checked })));
    report.dense_selection_debug = { phenotypes: state,
      markers: await page.$$eval('input[name="spatial-mac_markers_filter"]', boxes => boxes.map(x => ({value: x.value, checked: x.checked}))) };
    assert(state.length >= 20 && state.every(x => x.checked));
    assert.equal(state.filter(x => x.value.startsWith('Unclassified ')).length, 4);
    assert(!state.some(x => /^(Conflict(?:\s|_|$)|All unclassified cells$|All macrophages$|No subtype$)/i.test(x.value)));
    await dataPanel(false); return state;
  });
  if (scope === 'all') for (const width of [1200, 1600]) await check(`dense_tissue_${width}`, () => captureTissue('dense', width));
  if (scope === 'all') await check('real_brush_and_exact_reset', async () => {
    await scrollTo('#spatial-tissue_scatter'); await tissueReady();
    const before = await geometry(), panel = before.panel;
    await page.evaluate(() => {
      window.__spatialBrushEvents = [];
      jQuery(document).on('shiny:inputchanged.spatialCardReview', event => {
        if (event.name === 'spatial-scatter_brush' && event.value)
          window.__spatialBrushEvents.push(JSON.parse(JSON.stringify(event.value)));
      });
    });
    await page.mouse.move(panel.x + panel.width * .2, panel.y + panel.height * .25); await page.mouse.down();
    await page.mouse.move(panel.x + panel.width * .7, panel.y + panel.height * .65, { steps: 10 }); await page.mouse.up();
    await page.waitForFunction(domain => window.__spatialBrushEvents.length &&
      JSON.stringify(Shiny.shinyapp.$values['spatial-tissue_scatter'].coordmap.panels[0].domain) !== domain,
      {}, JSON.stringify(before.domain));
    await tissueReady(); const zoom = await geometry();
    assert(zoom.panel.width >= zoom.output.width - 150);
    await page.$eval('#spatial-scatter_reset_zoom', el => el.click());
    await page.waitForFunction(domain => JSON.stringify(Shiny.shinyapp.$values['spatial-tissue_scatter'].coordmap.panels[0].domain) === domain,
      {}, JSON.stringify(before.domain));
    await tissueReady(); const reset = await geometry();
    const events = await page.evaluate(() => window.__spatialBrushEvents);
    assert.deepEqual(reset.domain, before.domain); assert(reset.stableOutputElement);
    return { before: before.domain, zoom: zoom.domain, reset: reset.domain, events };
  });
  if (scope === 'all') await check('donor_changes_tissue_and_gradient_overlay_only', async () => {
    await page.waitForFunction(() => document.getElementById('spatial-immune_gradient')?.data?.some(t => /^Donor /.test(t.name || '')));
    await idle(); const before = await plotSnapshot();
    const next = await page.$eval('#spatial-donor', el => Object.keys(el.selectize?.options || {}).find(value => value !== el.value));
    assert(next, 'Second donor missing');
    await page.$eval('#main-content', el => { el.scrollTop = 0; });
    const handle = await page.evaluateHandle(() => document.getElementById('spatial-donor').closest('.shiny-input-container').querySelector('.selectize-input'));
    await handle.asElement().click();
    const option = `.selectize-dropdown .option[data-value=${JSON.stringify(next)}]`;
    await page.waitForSelector(option, { visible: true });
    await page.click(option); await handle.dispose();
    await page.waitForFunction(next => document.getElementById('spatial-donor').value === next &&
      document.getElementById('spatial-immune_gradient')?.data?.some(t => (t.name || '').startsWith(`Donor ${Number(next)} (`)), {}, next);
    await tissueReady(); const after = await plotSnapshot();
    assert.notEqual(before.tissue, after.tissue); assert.notDeepEqual(before.gradientDonor, after.gradientDonor);
    assert.deepEqual(before.gradientCohort, after.gradientCohort); assert.deepEqual(before.others, after.others);
    return { beforeDonor: before.donor, afterDonor: after.donor, changedOverlay: after.gradientDonor.map(t => t.name),
      unaffectedCards: Object.keys(after.others), cohortUnchanged: true, tissueChanged: true };
  });
  for (const width of [1200, 1600]) await check(`distance_controls_and_headings_${width}`, async () => {
    await page.setViewport({ width, height: 1300, deviceScaleFactor: 1 }); await idle();
    await page.waitForFunction(() => {
      const plot = document.getElementById('spatial-immune_distance_kde');
      return plot?._fullLayout && Math.abs(plot._fullLayout.width - plot.getBoundingClientRect().width) < 2;
    });
    await scrollTo('.distance-profile-controls');
    const detail = await page.evaluate(() => {
      const rect = el => { const r = el.getBoundingClientRect(); return { x: r.x, y: r.y, right: r.right, bottom: r.bottom, width: r.width, height: r.height }; };
      return { groups: [...document.querySelectorAll('.distance-profile-controls .control-group')].filter(el => el.getBoundingClientRect().width > 0)
        .map(el => ({ heading: el.querySelector('.control-group-title').textContent, ...rect(el), fits: el.scrollWidth <= el.clientWidth + 1 })),
        headings: [...document.querySelectorAll('.tab-pane.active .card h5')].map(el => el.textContent) };
    });
    assert(detail.groups.length >= 5 && detail.groups.every(g => g.fits));
    for (let i = 0; i < detail.groups.length; i++) for (let j = i + 1; j < detail.groups.length; j++) {
      const a = detail.groups[i], b = detail.groups[j];
      assert(a.right <= b.x + 1 || b.right <= a.x + 1 || a.bottom <= b.y + 1 || b.bottom <= a.y + 1, 'Control groups overlap');
    }
    assert(detail.headings.includes('Distance-To-Islet Profile'));
    assert(detail.headings.every(h => !/(?:^|[\s-])[a-z]/.test(h)), 'Card heading contains a lowercase word');
    await screenshot(`distance-controls-${width}`); return detail;
  });
  for (const width of [1200, 1600]) await check(`distance_macrophage_separation_grid_${width}`, async () => {
    await page.setViewport({ width, height: 1300, deviceScaleFactor: 1 });
    await scrollTo('.distance-profile-controls');
    if (!await page.$eval('input[name="spatial-dist_ptype"][value="composition"]', el => el.checked))
      await page.click('input[name="spatial-dist_ptype"][value="composition"]');
    await page.waitForSelector('#spatial-dist_comp_items .pheno-mac-markers', { visible: true }); await idle();
    await page.waitForFunction(() => {
      const plot = document.getElementById('spatial-immune_distance_kde');
      return plot?._fullLayout && Math.abs(plot._fullLayout.width - plot.getBoundingClientRect().width) < 2;
    });
    const detail = await page.$eval('#spatial-dist_comp_items .pheno-group-macrophage', el => {
      const master = el.querySelector('.pheno-mac-toggle') || el.querySelector('.pheno-group-boxes');
      const markers = el.querySelector('.pheno-mac-markers'), grid = markers.querySelector('.shiny-options-group');
      return { label: markers.querySelector('.control-label').textContent.trim(), markerCount: grid.querySelectorAll('input').length,
        below: markers.getBoundingClientRect().top >= master.getBoundingClientRect().bottom,
        display: getComputedStyle(grid).display, columns: getComputedStyle(grid).gridTemplateColumns,
        markerTokens: [...el.querySelectorAll('.pheno-mac-children label')].flatMap(label => {
          const tokens = [], walk = document.createTreeWalker(label, NodeFilter.SHOW_TEXT);
          while (walk.nextNode()) {
            for (const match of walk.currentNode.textContent.matchAll(/[A-Z][A-Za-z0-9]*[⁺⁻]/g)) {
              const range = document.createRange(); range.setStart(walk.currentNode, match.index);
              range.setEnd(walk.currentNode, match.index + match[0].length);
              tokens.push({ text: match[0], lineBoxes: range.getClientRects().length });
            }
          }
          return tokens;
        }) };
    });
    assert(detail.label === 'Separate by:' && detail.markerCount === 4 && detail.below && detail.display === 'grid');
    assert(detail.markerTokens.length >= 40 && detail.markerTokens.every(token => token.lineBoxes === 1), 'A marker token split across lines');
    await screenshot(`distance-composition-controls-${width}`); return { width, ...detail };
  });
  assert(!report.page_errors.length && !report.application_http_errors.length && !report.external_requests_blocked);
  report.status = 'passed';
} catch (error) {
  report.status = 'failed'; report.failure = error.stack || String(error);
} finally {
  report.finished_utc = new Date().toISOString(); await save(); await browser?.close();
}
console.log(JSON.stringify({ status: report.status, checks: report.checks.length, output: options.output, failure: report.failure || null }));
if (report.status !== 'passed') process.exitCode = 1;

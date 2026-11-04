#!/usr/bin/env node
// Focused donor-composition acceptance for an explicitly identified loopback preview.
// Source/image identity is supplied by the caller; this harness does not attest the running image itself.
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { createRequire } from 'node:module';
import puppeteer from 'puppeteer-core';

const require = createRequire(import.meta.url);
const options = {};
const allowed = ['url', 'output', 'chrome', 'shots', 'source-commit', 'runtime-kind', 'image-id'];
for (let index = 2; index < process.argv.length; index += 2) {
  const name = process.argv[index]?.replace(/^--/, '');
  if (!allowed.includes(name) || !process.argv[index + 1]) {
    throw new Error('Usage: check-composition-controls.mjs --url LOOPBACK_URL --output REPORT.json --source-commit COMMIT --runtime-kind native-image|source-overlay [--image-id SHA256] [--shots DIR] [--chrome PATH]');
  }
  options[name] = process.argv[index + 1];
}
for (const name of ['url', 'output', 'source-commit', 'runtime-kind']) {
  if (!options[name]) throw new Error(`--${name} is required`);
}
if (!/^[0-9a-f]{40}$/.test(options['source-commit'])) throw new Error('--source-commit must be a full commit');
if (!['native-image', 'source-overlay'].includes(options['runtime-kind'])) throw new Error('Invalid --runtime-kind');
if (options['image-id'] && !/^sha256:[0-9a-f]{64}$/.test(options['image-id'])) throw new Error('Invalid --image-id');
if (options['runtime-kind'] === 'native-image' && !options['image-id']) throw new Error('Native runs require --image-id');
const target = new URL(options.url);
if (!['127.0.0.1', 'localhost', '[::1]'].includes(target.hostname) || target.protocol !== 'http:' ||
    target.username || target.password || target.search || target.hash) {
  throw new Error('Only plain HTTP loopback application URLs without credentials or queries are allowed');
}
if (!target.pathname.endsWith('/')) target.pathname += '/';
const output = path.resolve(options.output);
const shots = path.resolve(options.shots || `${output.replace(/\.json$/, '')}-screenshots`);
const report = {
  schema: 'islet.composition_controls_browser/1',
  status: 'running', startedAt: new Date().toISOString(),
  target: target.href,
  identity: { sourceCommit: options['source-commit'], runtimeKind: options['runtime-kind'],
    imageId: options['image-id'] || null, binding: 'caller-supplied; verify with the companion runtime receipt' },
  harness: { path: fileURLToPath(import.meta.url),
    sha256: crypto.createHash('sha256').update(await fs.readFile(fileURLToPath(import.meta.url))).digest('hex') },
  versions: { puppeteer: require('puppeteer-core/package.json').version },
  checks: [], screenshots: [],
  browserErrors: { page: 0, console: 0 },
  httpErrors: { missingDefaultFavicon: 0, application: 0 }, externalRequestsBlocked: 0,
  limitations: ['This focused loopback check does not establish public access, image security, scientific typing accuracy, or reboot persistence.',
    'Screenshots are generated for direct visual review; geometry assertions alone are not a visual review.'],
};
let browser;
let page;
const plot = '#plot-composition_by_donor';
const include = '#plot-composition_include';
const panel = '.composition-include details.control-panel';
const markerBoxes = 'input[name="plot-mac_markers_include"]';
const ownBoxes = `${include} input[name="plot-composition_include"]`;
const macGroup = `${include} .pheno-group-macrophage`;

async function save() {
  await fs.mkdir(path.dirname(output), { recursive: true });
  await fs.writeFile(output, `${JSON.stringify(report, null, 2)}\n`, { mode: 0o600 });
}
async function check(name, run) {
  try {
    const detail = await run();
    const passed = detail === true || detail?.passed === true;
    report.checks.push({ name, passed, ...(typeof detail === 'object' ? detail : {}) });
    await save();
    console.log(`${passed ? 'PASS' : 'FAIL'} ${name}`);
    if (!passed) throw new Error(`Check failed: ${name}`);
  } catch (error) {
    if (report.checks.at(-1)?.name !== name) report.checks.push({ name, passed: false, error: String(error) });
    await save();
    throw error;
  }
}
async function settled() {
  await page.evaluate(() => { window.__compositionStable = null; });
  await page.waitForFunction(() => {
    const gd = document.getElementById('plot-composition_by_donor');
    if (document.documentElement.classList.contains('shiny-busy') || gd?.classList.contains('recalculating')) {
      window.__compositionStable = null;
      return false;
    }
    const signature = JSON.stringify({
      title: gd?.layout?.title?.text,
      bars: (gd?.data || []).filter(trace => trace.type === 'bar').map(trace => [trace.name, trace.marker?.color, trace.y]),
      palette: document.getElementById('composition_palette')?.value,
      drilldown: document.getElementById('drilldown_palette')?.value,
    });
    if (window.__compositionStable?.signature !== signature) {
      window.__compositionStable = { signature, since: performance.now() };
      return false;
    }
    return performance.now() - window.__compositionStable.since >= 300;
  }, { polling: 100 });
}
async function waitBars() {
  await page.waitForFunction(() => {
    const gd = document.getElementById('plot-composition_by_donor');
    return gd?._fullLayout && !gd.classList.contains('recalculating') &&
      gd.data?.some(trace => trace.type === 'bar' && Array.from(trace.y || []).some(Number.isFinite));
  });
  await settled();
}
async function waitEmpty() {
  await page.waitForFunction(() => {
    const gd = document.getElementById('plot-composition_by_donor');
    return gd?._fullLayout && !gd.classList.contains('recalculating') &&
      !gd.data?.some(trace => trace.type === 'bar') && /No phenotypes selected/.test(gd.layout?.title?.text || '');
  });
  await settled();
}
async function snapshot() {
  return page.evaluate(() => {
    const gd = document.getElementById('plot-composition_by_donor');
    const bars = (gd?.data || []).filter(trace => trace.type === 'bar');
    const sums = new Map();
    for (const trace of bars) Array.from(trace.x || []).forEach((donor, i) => {
      const value = Number(trace.y[i]);
      if (Number.isFinite(value)) sums.set(donor, (sums.get(donor) || 0) + value);
    });
    const selected = [...document.querySelectorAll('input[name="plot-composition_include"]:checked')].map(box => box.value);
    const markers = [...document.querySelectorAll('input[name="plot-mac_markers_include"]:checked')].map(box => box.value);
    return { barCount: bars.length, donorCount: sums.size, sums: [...sums.values()],
      names: bars.map(trace => trace.name), colors: bars.map(trace => trace.marker?.color),
      yTitle: gd?.layout?.yaxis?.title?.text || '', title: gd?.layout?.title?.text || '', selected, markers };
  });
}
function normalized(state) {
  return state.barCount > 0 && state.donorCount >= 15 && state.sums.every(sum => Math.abs(sum - 100) < 1e-6);
}
const affectedPhenotypes = ['Endothelium Blood Vessel', 'Endothelium Lymphatic', 'Stromal'];
function distinctPhenotypeColors(state) {
  const colors = state.colors.map(color => String(color || '').toLowerCase());
  return colors.length === state.barCount && colors.every(color => color && color !== '#8f8f8f') &&
    new Set(colors).size === colors.length && affectedPhenotypes.every(name => state.names.includes(name));
}
async function scrollCard() {
  await page.evaluate(() => {
    const card = document.getElementById('plot-composition_by_donor').closest('.card');
    const scroller = document.getElementById('main-content');
    const sticky = document.querySelector('.tab-pane[data-value="Plot"] .tab-controls');
    const offset = Math.max(document.querySelector('.app-header')?.getBoundingClientRect().bottom || 0,
      sticky?.getBoundingClientRect().bottom || 0) + 12;
    scroller.scrollTop += card.getBoundingClientRect().top - offset;
  });
}
async function setPanel(open) {
  const current = await page.$eval(panel, element => element.open);
  if (current !== open) {
    await scrollCard();
    await page.click(`${panel} > summary`);
    await page.waitForFunction((selector, wanted) => document.querySelector(selector)?.open === wanted, {}, panel, open);
  }
}
async function screenshot(name) {
  await fs.mkdir(shots, { recursive: true });
  const destination = path.join(shots, `${name}.png`);
  await page.screenshot({ path: destination });
  report.screenshots.push({ name, path: destination,
    sha256: crypto.createHash('sha256').update(await fs.readFile(destination)).digest('hex') });
}
async function selectPalette(value) {
  await scrollCard();
  const control = await page.evaluateHandle(() => document.getElementById('composition_palette')
    .closest('.shiny-input-container').querySelector('.selectize-input'));
  await control.asElement().click();
  const option = `.selectize-dropdown .option[data-value=${JSON.stringify(value)}]`;
  await page.waitForSelector(option, { visible: true });
  await page.click(option);
  await control.dispose();
}

try {
  await save();
  browser = await puppeteer.launch({ executablePath: options.chrome || '/usr/bin/google-chrome', headless: true,
    args: ['--disable-dev-shm-usage', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'] });
  report.versions.chrome = await browser.version();
  page = await browser.newPage();
  page.setDefaultTimeout(120000);
  await page.setViewport({ width: 1600, height: 1100, deviceScaleFactor: 1 });
  await page.setRequestInterception(true);
  page.on('request', request => {
    const url = new URL(request.url());
    if (['http:', 'https:'].includes(url.protocol) && url.origin !== target.origin) {
      report.externalRequestsBlocked += 1;
      void request.abort('blockedbyclient');
    } else void request.continue();
  });
  page.on('pageerror', () => { report.browserErrors.page += 1; });
  page.on('console', message => { if (message.type() === 'error') report.browserErrors.console += 1; });
  page.on('response', response => {
    if (response.status() < 400) return;
    const url = new URL(response.url());
    if (response.status() === 404 && url.pathname.endsWith('/favicon.ico')) report.httpErrors.missingDefaultFavicon += 1;
    else report.httpErrors.application += 1;
  });
  await page.goto(target.href, { waitUntil: 'domcontentloaded', timeout: 120000 });
  await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
  await page.evaluate(() => {
    window.__compositionInputEvents = [];
    window.jQuery(document).on('shiny:inputchanged.compositionReview', event => {
      if (['plot-composition_include', 'plot-mac_markers_include', 'composition_palette',
        'drilldown_palette', 'plot-linked_click'].includes(event.name)) {
        window.__compositionInputEvents.push({ name: event.name, at: performance.now() });
      }
    });
  });
  await check('first_render_with_include_panel_collapsed', async () => {
    await waitBars();
    const state = await snapshot();
    const closed = await page.$eval(panel, element => !element.open);
    return { passed: closed && normalized(state) && distinctPhenotypeColors(state), closed, ...state };
  });
  await check('include_groups_axis_and_no_subtitle', async () => {
    await setPanel(true);
    const state = await snapshot();
    const groups = await page.$$eval(`${include} .pheno-group-title`, elements => elements.map(el => el.textContent.trim()));
    const labels = await page.$$eval(ownBoxes, boxes => boxes.map(box => box.value));
    return { passed: groups.length >= 4 && groups.includes('Macrophage') && groups.includes('Other') &&
      state.yTitle === 'Proportion of selected phenotypes' && !/<br>|<sup>/.test(state.title) &&
      !labels.some(label => /^(?:Conflict(?:\s|_|$)|Unclassified )/i.test(label) ||
        ['All unclassified cells', 'All macrophages', 'No subtype'].includes(label)),
    groups, yTitle: state.yTitle, title: state.title, unclassifiedSubclasses: labels.filter(label => /^Unclassified /.test(label)).length };
  });
  await check('none_clears_plot_and_stays_empty_across_marker_rerender', async () => {
    await page.click(`${include} .pheno-groups-overall .pheno-select-none`);
    await waitEmpty();
    await page.click(`${markerBoxes}[value="CD68"]`);
    await page.waitForSelector(`${macGroup} .pheno-mac-children input`);
    await page.waitForFunction(() => document.querySelector('input[name="plot-mac_markers_include"][value="CD68"]')?.checked);
    await waitEmpty();
    await setPanel(false);
    await setPanel(true);
    await settled();
    const state = await snapshot();
    return { passed: state.selected.length === 0 && state.barCount === 0 && state.markers.join() === 'CD68', ...state };
  });
  await check('split_macrophage_master_selects_children_not_parent_or_markers', async () => {
    await page.click(`${macGroup} .pheno-group-toggle`);
    await waitBars();
    const enabled = await snapshot();
    const children = await page.$$eval(`${macGroup} input[name="plot-composition_include"]`, boxes => boxes.map(box => box.value));
    const masterName = await page.$eval(`${macGroup} .pheno-group-toggle`, box => box.getAttribute('name'));
    const selectedChildren = children.every(value => enabled.selected.includes(value));
    await page.focus(`${macGroup} .pheno-group-toggle`);
    await page.keyboard.press('Space');
    await waitEmpty();
    const keyboardOff = await snapshot();
    await page.keyboard.press('Space');
    await waitBars();
    const keyboardOn = await snapshot();
    await page.click(`${macGroup} .pheno-group-toggle`);
    await waitEmpty();
    const disabled = await snapshot();
    return { passed: normalized(enabled) && selectedChildren && enabled.selected.length === children.length &&
      children.length >= 2 && !children.includes('Macrophage') && !masterName &&
      enabled.markers.join() === 'CD68' && disabled.markers.join() === 'CD68' && disabled.selected.length === 0 &&
      keyboardOff.selected.length === 0 && keyboardOff.markers.join() === 'CD68' &&
      normalized(keyboardOn) && keyboardOn.selected.length === children.length && keyboardOn.markers.join() === 'CD68',
    childCount: children.length, masterHasName: Boolean(masterName), enabled, keyboardOff, keyboardOn, disabled };
  });
  await check('all_restores_disjoint_categories_and_complete_bars', async () => {
    await page.click(`${include} .pheno-groups-overall .pheno-select-all`);
    await waitBars();
    const state = await snapshot();
    const choiceCount = await page.$$eval(ownBoxes, boxes => boxes.length);
    return { passed: normalized(state) && distinctPhenotypeColors(state) && state.selected.length === choiceCount && state.markers.join() === 'CD68' &&
      !state.selected.includes('Macrophage') && !state.selected.includes('All unclassified cells') &&
      !state.selected.some(label => /^Conflict(?:\s|_|$)/i.test(label)), choiceCount, ...state };
  });
  await check('macrophage_marker_boxes_are_below_master_with_separate_by_label', async () => {
    return page.$eval(macGroup, group => {
      const parent = group.querySelector('.pheno-group-toggle').closest('.checkbox').getBoundingClientRect();
      const markers = [...group.querySelectorAll('input[name="plot-mac_markers_include"]')];
      const label = group.querySelector('.pheno-mac-markers .control-label');
      const bounds = markers.map(box => box.closest('label').getBoundingClientRect());
      return { passed: markers.length === 4 && label?.textContent.trim() === 'Separate by:' &&
        bounds.every(box => box.top >= parent.bottom - 1) && label.getBoundingClientRect().top >= parent.bottom - 1,
      markerCount: markers.length, label: label?.textContent.trim(),
      parentBottom: parent.bottom, markerTops: bounds.map(box => box.top) };
    });
  });
  await check('palette_changes_bar_colors_and_synchronizes_spatial', async () => {
    const before = await snapshot();
    const current = await page.$eval('#composition_palette', element => element.selectize?.getValue() || element.value);
    const next = current === 'Colorblind Safe' ? 'Maximum Distinction' : 'Colorblind Safe';
    await selectPalette(next);
    await page.waitForFunction(old => {
      const gd = document.getElementById('plot-composition_by_donor');
      return !gd?.classList.contains('recalculating') &&
        JSON.stringify((gd?.data || []).filter(t => t.type === 'bar').map(t => t.marker?.color)) !== old;
    }, {}, JSON.stringify(before.colors));
    await settled();
    const after = await snapshot();
    await page.click('#tabs a[data-value="Spatial"]');
    await page.waitForFunction(() => document.querySelector('#tabs li.active a')?.dataset.value === 'Spatial');
    await page.waitForFunction(value => {
      const element = document.getElementById('spatial_palette');
      return (element?.selectize?.getValue() || element?.value) === value;
    }, {}, next);
    const spatial = await page.$eval('#spatial_palette', element => element.selectize?.getValue() || element.value);
    const expectedColors = Object.fromEntries(affectedPhenotypes.map(name => [name, after.colors[after.names.indexOf(name)]]));
    await page.waitForFunction(expected => {
      const bars = (document.getElementById('spatial-cluster_composition')?.data || []).filter(trace => trace.type === 'bar');
      return Object.entries(expected).every(([name, color]) => bars.some(trace =>
        trace.name === name && trace.marker?.color?.toLowerCase() === color?.toLowerCase()));
    }, {}, expectedColors);
    await page.click('#tabs a[data-value="Plot"]');
    await page.waitForFunction(() => document.querySelector('#tabs li.active a')?.dataset.value === 'Plot');
    await waitBars();
    return { passed: JSON.stringify(before.colors) !== JSON.stringify(after.colors) &&
      JSON.stringify(before.sums) === JSON.stringify(after.sums) && spatial === next && normalized(after) &&
      distinctPhenotypeColors(before) && distinctPhenotypeColors(after),
    beforePalette: current, afterPalette: next, spatialPalette: spatial,
    beforeColors: before.colors, afterColors: after.colors, spatialColors: expectedColors, donorCount: after.donorCount };
  });
  await check('palette_persists_when_linked_object_mounts_drilldown', async () => {
    const before = await snapshot();
    const wanted = await page.$eval('#composition_palette', element => element.selectize?.getValue() || element.value);
    const absentInitially = await page.$('#drilldown_palette') === null;
    await page.evaluate(() => { document.getElementById('main-content').scrollTop = 0; });
    await page.waitForFunction(() => document.getElementById('plot-plt')?.data?.some(trace =>
      trace.meta?.trajectory_role === 'islet' && trace.x?.length));
    const point = await page.$eval('#plot-plt', gd => {
      const rect = gd.querySelector('.svg-container').getBoundingClientRect();
      const x = gd._fullLayout.xaxis;
      const y = gd._fullLayout.yaxis;
      for (const trace of gd.data.filter(item => item.meta?.trajectory_role === 'islet')) {
        for (let i = 0; i < trace.x.length; i += 1) {
          const px = x.l2p(trace.x[i]);
          const py = y.l2p(trace.y[i]);
          const screenX = rect.x + x._offset + px;
          const screenY = rect.y + y._offset + py;
          if (px > 30 && px < x._length - 30 && py > 30 && py < y._length - 30 &&
              screenY > 200 && screenY < innerHeight - 30 && gd.contains(document.elementFromPoint(screenX, screenY))) {
            return { x: screenX, y: screenY };
          }
        }
      }
      return null;
    });
    if (!point) throw new Error('No visible linked Plot point was available for a physical click');
    await page.mouse.move(point.x, point.y);
    await page.waitForFunction(() => document.getElementById('plot-plt')._hoverdata?.some(item =>
      item.data.meta?.trajectory_role === 'islet'));
    await page.mouse.click(point.x, point.y);
    await page.waitForSelector('#drilldown_palette');
    await page.waitForFunction(() => {
      const output = document.getElementById('islet_segmentation_view');
      const image = output?.querySelector('img');
      return image?.complete && image.naturalWidth > 0 && !output.classList.contains('recalculating');
    });
    await settled();
    const readPalettes = () => page.evaluate(() => Object.fromEntries(
      ['composition_palette', 'spatial_palette', 'drilldown_palette'].map(id => {
        const element = document.getElementById(id);
        return [id, element?.selectize?.getValue() || element?.value];
      })));
    const palettes = await readPalettes();
    const after = await snapshot();
    const physicalClicks = await page.evaluate(() => window.__compositionInputEvents
      .filter(event => event.name === 'plot-linked_click').length);
    const firstMountPassed = absentInitially && physicalClicks > 0 && Object.values(palettes).every(value => value === wanted) &&
      JSON.stringify(before.colors) === JSON.stringify(after.colors);
    const remounts = [];
    if (firstMountPassed) {
      for (const [tabName, namespace] of [['Spatial', 'spatial'], ['Plot', 'plot']]) {
        await page.click(`#tabs a[data-value="${tabName}"]`);
        await page.waitForSelector(`#${namespace}-segmentation_viewer_panel #drilldown_palette`);
        await page.waitForFunction(id => {
          const image = document.querySelector(`#${id}-segmentation_viewer_panel #islet_segmentation_view img`);
          return image?.complete && image.naturalWidth > 0;
        }, {}, namespace);
        await settled();
        remounts.push({ tab: tabName, palettes: await readPalettes(), barColors: (await snapshot()).colors });
      }
      await waitBars();
    }
    const final = await snapshot();
    return { passed: firstMountPassed && remounts.length === 2 &&
      remounts.every(item => Object.values(item.palettes).every(value => value === wanted) &&
        JSON.stringify(item.barColors) === JSON.stringify(before.colors)) &&
      JSON.stringify(before.colors) === JSON.stringify(final.colors), absentInitially, physicalClicks,
    expectedPalette: wanted, palettes, remounts, beforeColors: before.colors, afterColors: final.colors };
  });
  await check('clear_marker_preserves_selected_population_without_duplicate_parent', async () => {
    await setPanel(true);
    await page.click(`${markerBoxes}[value="CD68"]`);
    await page.waitForSelector(`${macGroup} input[name="plot-composition_include"][value="Macrophage"]`);
    await page.waitForFunction(() => !document.querySelector('input[name="plot-mac_markers_include"][value="CD68"]')?.checked);
    await waitBars();
    const state = await snapshot();
    return { passed: normalized(state) && state.markers.length === 0 && state.selected.includes('Macrophage') &&
      !state.selected.some(value => /^Macrophage .*[+-]/.test(value) || /^CD68[+-]/.test(value)), ...state };
  });
  await setPanel(false);
  for (const width of [1200, 1600]) {
    await check(`compact_card_geometry_${width}`, async () => {
      await page.setViewport({ width, height: 1100, deviceScaleFactor: 1 });
      await scrollCard();
      await page.waitForFunction(() => {
        const gd = document.getElementById('plot-composition_by_donor');
        return Math.abs((gd?._fullLayout?.width || 0) - (gd?.getBoundingClientRect().width || 1)) < 3;
      });
      await settled();
      const geometry = await page.$eval(plot, element => {
        const card = element.closest('.card');
        const c = card.getBoundingClientRect();
        const p = element.getBoundingClientRect();
        const toolbar = card.querySelector('.composition-toolbar');
        const t = toolbar.getBoundingClientRect();
        const controls = [...toolbar.querySelectorAll('.selectize-control')].map(el => {
          const r = el.getBoundingClientRect();
          const input = el.querySelector('.selectize-input')?.getBoundingClientRect();
          const item = el.querySelector('.selectize-input > .item')?.getBoundingClientRect();
          return { left: r.left, right: r.right, width: r.width,
            rightTextGutter: input && item ? input.right - item.right : null };
        });
        return { cardWidth: c.width, plotWidth: p.width, abovePlot: p.top - c.top,
          toolbarHeight: t.height, closed: !card.querySelector('details.control-panel').open,
          overflow: card.scrollWidth > card.clientWidth + 1,
          controls, controlsInsideCard: controls.length === 3 && controls.every(r =>
            r.left >= c.left && r.right <= c.right + 1 && r.width >= 60 &&
            (r.rightTextGutter === null || r.rightTextGutter >= 30)) };
      });
      await screenshot(`composition-${width}-collapsed`);
      return { passed: geometry.closed && !geometry.overflow && geometry.controlsInsideCard &&
        geometry.abovePlot <= 190 && geometry.plotWidth >= geometry.cardWidth * 0.9, width, ...geometry };
    });
  }
  await setPanel(true);
  await scrollCard();
  await screenshot('composition-1600-expanded');
  await check('real_control_events_and_no_application_errors', async () => {
    const events = await page.evaluate(() => window.__compositionInputEvents || []);
    const counts = Object.fromEntries(['plot-composition_include', 'plot-mac_markers_include', 'composition_palette']
      .map(name => [name, events.filter(event => event.name === name).length]));
    const visibleErrors = await page.$$eval('.shiny-output-error', nodes => nodes.filter(node => node.getClientRects().length).map(node => node.id));
    return { passed: Object.values(counts).every(count => count > 0) && visibleErrors.length === 0 &&
      report.browserErrors.page === 0 && report.browserErrors.console <= report.httpErrors.missingDefaultFavicon &&
      report.httpErrors.application === 0 && report.externalRequestsBlocked === 0,
    inputEventCounts: counts, visibleErrorIds: visibleErrors };
  });
  report.status = 'passed';
} catch (error) {
  report.status = 'failed';
  report.failure = String(error);
  process.exitCode = 1;
  if (page) {
    report.failureState = await snapshot().catch(() => null);
    report.failureInputs = await page.evaluate(() => ({
      included: window.Shiny?.shinyapp?.$inputValues?.['plot-composition_include'],
      markers: window.Shiny?.shinyapp?.$inputValues?.['plot-mac_markers_include'],
      masterChecked: document.querySelector('#plot-composition_include .pheno-group-toggle')?.checked,
      events: window.__compositionInputEvents || [],
    })).catch(() => null);
    await screenshot('composition-failure').catch(() => {});
  }
} finally {
  report.finishedAt = new Date().toISOString();
  if (browser) await browser.close();
  await save();
  console.log(JSON.stringify({ status: report.status, checks: report.checks.length,
    passed: report.checks.filter(item => item.passed).length, output }));
}

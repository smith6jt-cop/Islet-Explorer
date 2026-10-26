#!/usr/bin/env node
// Browser acceptance against an explicitly selected local/staged application.
// Evidence contains checks and rule counts, never page text, scientific rows,
// download bodies, screenshots, console messages, image names, or credentials.
import fs from 'node:fs/promises';
import path from 'node:path';
import { createRequire } from 'node:module';
import puppeteer from 'puppeteer-core';

const require = createRequire(import.meta.url);
// Expanded cohort: the original 15 plus these additions. Donor 78 (HDL078/HPAP-078)
// was removed on 2026-10-15 (scripts/excluded_donors.csv) and must not appear.
const ADDITIONS = ['6090', '6232', '6520', '6539', '6558'];
const EXCLUDED = ['78'];
const EXPANDED_DONORS = 15 + ADDITIONS.length;
const options = {};
for (let index = 2; index < process.argv.length; index += 2) {
  const key = process.argv[index];
  if (!['--url', '--output', '--chrome', '--expanded', '--auth-file'].includes(key) || !process.argv[index + 1]) {
    throw new Error('Usage: check-app.mjs --url APP_URL --output REPORT.json [--chrome CHROME_PATH] [--expanded true] [--auth-file JSON]');
  }
  options[key.slice(2)] = process.argv[index + 1];
}
if (!options.url || !options.output) throw new Error('--url and --output are required');
const appUrl = new URL(options.url);
if (!['http:', 'https:'].includes(appUrl.protocol) || appUrl.username || appUrl.password ||
    appUrl.search || appUrl.hash) throw new Error('Supply an HTTP(S) application URL without credentials or query parameters');
if (!appUrl.pathname.endsWith('/')) appUrl.pathname += '/';
// Through UFRC's proxy the whole site sits behind Basic Auth. The credential comes
// from the same owner-only 0600 JSON file as deploy/verify-reviewer-access.py and
// is sent only to the HTTPS origin; it never reaches the report.
let credentials = null;
if (options['auth-file']) {
  if (appUrl.protocol !== 'https:') throw new Error('--auth-file requires an https application URL');
  const authPath = path.resolve(options['auth-file']);
  const authStat = await fs.lstat(authPath);
  if (!authStat.isFile() || (authStat.mode & 0o777) !== 0o600) {
    throw new Error('--auth-file must be a regular file with mode 0600');
  }
  const parsed = JSON.parse(await fs.readFile(authPath, 'utf8'));
  if (!parsed || typeof parsed !== 'object' || Object.keys(parsed).sort().join(',') !== 'password,username' ||
      typeof parsed.username !== 'string' || typeof parsed.password !== 'string') {
    throw new Error('--auth-file must hold exactly {"username": ..., "password": ...}');
  }
  credentials = { username: parsed.username, password: parsed.password };
}

const report = {
  schemaVersion: 1,
  runState: 'running',
  startedAt: new Date().toISOString(),
  target: { origin: appUrl.origin, applicationPath: appUrl.pathname, basicAuth: Boolean(credentials) },
  versions: {
    puppeteer: require('puppeteer-core/package.json').version,
    axe: require('axe-core/package.json').version,
  },
  checks: [],
  accessibility: [],
  browserErrors: { page: 0, console: 0 },
  httpErrors: { missingDefaultFavicon: 0, optionalImageSidecar404: 0, application: 0, applicationPaths: [] },
  httpErrorClasses: {},
  externalRequestsBlocked: 0,
  limitations: [
    'Automated checks are a subset of accessibility review; visualizations and assistive-technology use still require manual review.',
    credentials
      ? 'This run went through the HTTPS origin with Basic Auth; it does not establish production-container security.'
      : 'This run does not establish external Basic Auth, UFRC proxy behavior, or production-container security.',
  ],
};
const axeSource = await fs.readFile(require.resolve('axe-core/axe.min.js'), 'utf8');
let browser;

async function saveReport() {
  await fs.mkdir(path.dirname(path.resolve(options.output)), { recursive: true });
  await fs.writeFile(options.output, `${JSON.stringify(report, null, 2)}\n`, { mode: 0o600 });
}

async function check(name, operation) {
  try {
    const detail = await operation();
    const passed = detail === true || (detail && detail.passed === true);
    report.checks.push({ name, passed, ...(typeof detail === 'object' ? detail : {}) });
    await saveReport();
    return passed;
  } catch (error) {
    report.checks.push({ name, passed: false, error: error.name === 'TimeoutError' ? 'timeout' : 'check_failed' });
    await saveReport();
    return false;
  }
}

async function instrument(page) {
  page.setDefaultTimeout(180000);
  if (credentials) await page.authenticate(credentials);
  await page.setViewport({ width: 1600, height: 1100, deviceScaleFactor: 1 });
  await page.setRequestInterception(true);
  page.on('request', request => {
    const url = new URL(request.url());
    if (['http:', 'https:'].includes(url.protocol) && url.origin !== appUrl.origin) {
      report.externalRequestsBlocked += 1;
      void request.abort('blockedbyclient');
    } else {
      void request.continue();
    }
  });
  page.on('pageerror', () => { report.browserErrors.page += 1; });
  page.on('console', message => { if (message.type() === 'error') report.browserErrors.console += 1; });
  page.on('response', response => {
    if (response.status() < 400) return;
    const url = new URL(response.url());
    if (response.status() === 404 && url.origin === appUrl.origin && url.pathname.endsWith('/favicon.ico')) {
      report.httpErrors.missingDefaultFavicon += 1;
    } else if (response.status() === 404 && url.origin === appUrl.origin &&
        url.pathname.startsWith(`${appUrl.pathname}images/`) && url.pathname.endsWith('.offsets.json')) {
      // Avivator tries an optional sidecar before scanning OME-TIFF IFDs itself.
      // The Viewer canvas/range checks below verify that its native fallback works.
      report.httpErrors.optionalImageSidecar404 += 1;
    } else {
      report.httpErrors.application += 1;
      const normalizedPath = url.pathname
        .replace(/[A-Za-z][A-Za-z0-9+.-]*_\d+(?=\.ome\.tiff?$)/gi, '{image}')
        .replace(/\/\d+(?=\/|$)/g, '/{id}');
      if (!report.httpErrors.applicationPaths.includes(normalizedPath)) report.httpErrors.applicationPaths.push(normalizedPath);
      const path = url.pathname.toLowerCase();
      const category = path.endsWith('.ome.tiff') || path.endsWith('.ome.tif') ? 'image' :
        path.includes('/avivator/') ? 'viewer' : path.includes('/favicon') ? 'favicon' : 'app';
      report.httpErrorClasses[`${response.status()}:${category}`] =
        (report.httpErrorClasses[`${response.status()}:${category}`] || 0) + 1;
    }
  });
}

async function accessibility(context, scope) {
  await context.evaluate(axeSource);
  const result = await context.evaluate(async () => {
    const audit = await window.axe.run(document, {
      runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'] },
      iframes: false,
      resultTypes: ['violations', 'incomplete'],
    });
    const summarize = results => results.map(item => ({
      rule: item.id, impact: item.impact, affectedNodes: item.nodes.length,
    }));
    return { violations: summarize(audit.violations), needsReview: summarize(audit.incomplete) };
  });
  report.accessibility.push({ scope, ...result });
  return { passed: !result.violations.some(item => ['serious', 'critical'].includes(item.impact)) };
}

async function rendered(page, selector, kind) {
  await page.waitForFunction((target, outputKind) => {
    const element = document.querySelector(target);
    if (!element || element.classList.contains('shiny-output-error') ||
        element.classList.contains('recalculating')) return false;
    if (outputKind === 'plotly') return Boolean(element._fullLayout && element.data?.length);
    const image = element.querySelector('img');
    return Boolean(image?.complete && image.naturalWidth > 0);
  }, {}, selector, kind);
  return true;
}

async function download(page, selector) {
  await page.waitForFunction(target => {
    const link = document.querySelector(target);
    return link && !link.classList.contains('disabled') && link.getAttribute('href')?.includes('/download/');
  }, {}, selector);
  return page.$eval(selector, async link => {
    const response = await fetch(link.href, { credentials: 'same-origin' });
    const text = await response.text();
    const csv = text.includes(',') && text.includes('\n') && !text.trimStart().startsWith('<');
    return {
      passed: response.ok && csv && /attachment/i.test(response.headers.get('content-disposition') || ''),
      status: response.status,
      bytes: new TextEncoder().encode(text).length,
    };
  });
}

try {
  await saveReport();
  browser = await puppeteer.launch({
    executablePath: options.chrome || '/usr/bin/google-chrome',
    headless: true,
    args: ['--disable-dev-shm-usage', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'],
  });
  report.versions.chrome = await browser.version();
  const page = await browser.newPage();
  await instrument(page);
  await page.goto(appUrl.href, { waitUntil: 'domcontentloaded', timeout: 180000 });
  await check('shiny_connected', async () => {
    await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
    return true;
  });
  await check('plot_population_short_vocabulary_and_hidden_help', async () => {
    await page.waitForSelector('#plot-which option');
    return page.evaluate(() => {
      const selector = document.querySelector('#plot-which');
      const labels = selector?.selectize
        ? Object.values(selector.selectize.options).map(option => option.label)
        : [...(selector?.options || [])].map(option => option.textContent.trim());
      const retired = ['Endocrine', 'Acinar', 'Endothelial', 'T cell', 'Immune',
        'Unknown', 'ECAD+', 'Structural', 'APCs'];
      // The QuPath hormone-positivity fractions were withdrawn on 2026-10-14 (not comparable across the
      // cohort), so Cell Populations is the default layer and carries its % / Counts / Density selector.
      const hormoneLabels = ['Insulin+ positivity fraction', 'Glucagon+ positivity fraction',
        'Somatostatin+ positivity fraction'];
      const metrics = [...document.querySelectorAll('input[name="plot-metric"]')].map(input => input.value);
      // Help sits on each card and control panel (card_help(), 2026-10-20), closed until asked for.
      const plotTab = document.querySelector('.tab-pane[data-value="Plot"]');
      const helpTitles = [...(plotTab?.querySelectorAll('details.card-help > summary') || [])]
        .map(summary => summary.title);
      // Since 2026-10-21 each plot card's title is its chosen plot type, so its help is named by side.
      const wanted = ['Data controls help', 'View controls help', 'Left plot help', 'Right plot help',
        'Composition by donor help'];
      return labels.length > 0 && !labels.some(label => hormoneLabels.includes(label)) &&
        !labels.some(label => retired.includes(label)) && metrics.length === 3 &&
        wanted.every(title => helpTitles.includes(title)) &&
        !plotTab.querySelector('details.card-help[open]') &&
        !document.querySelector('#plot-dynamic_selector .help-block, #plot-metric_selector .help-block, #plot-dynamic_selector .phenotype-hint');
    });
  });
  await check('card_help_opens_on_its_card', async () => {
    // The Data controls "?" sits outside the panel's <summary>: opening it and clicking its text must not close the
    // panel, the popover must stay inside the window, one popover is open at a time, and Escape closes it.
    const panelOpen = () => page.evaluate(() =>
      document.querySelector('.tab-pane[data-value="Plot"] .control-panel-wrap > details.control-panel').open);
    const before = await panelOpen();
    await page.click('.tab-pane[data-value="Plot"] summary[title="Data controls help"]');
    // An open panel is moved to <body> (www/card-help.js) and returns to its <details> when closed.
    // The toggle event is asynchronous: wait for the panel to arrive.
    await page.waitForSelector('body > .card-help-panel p', { visible: true, timeout: 10000 });
    await page.click('body > .card-help-panel p');
    const opened = await page.evaluate(() => {
      const help = document.querySelector('.tab-pane[data-value="Plot"] summary[title="Data controls help"]')
        .parentElement;
      const panel = document.querySelector('body > .card-help-panel');
      const box = panel.getBoundingClientRect();
      return {open: help.open, box: [box.left, box.top, box.right, box.bottom, window.innerWidth, window.innerHeight].map(Math.round),
        inside: box.left >= 0 && box.right <= window.innerWidth && box.top >= 0 &&
        box.bottom <= window.innerHeight + 1, rules: Boolean(panel.querySelector('#show_phenotype_rules_plot'))};
    });
    await page.click('.tab-pane[data-value="Plot"] summary[title="Left plot help"]');
    await page.waitForFunction(() => document.querySelectorAll('body > .card-help-panel').length === 1 &&
      document.querySelector('summary[title="Left plot help"]').parentElement.open, { timeout: 10000 });
    const single = await page.evaluate(() =>
      [...document.querySelectorAll('details.card-help[open] > summary')].map(summary => summary.title));
    await page.keyboard.press('Escape');
    const closed = await page.waitForFunction(() => !document.querySelector('details.card-help[open]') &&
      !document.querySelector('body > .card-help-panel'), { timeout: 10000 }).then(() => true, () => false);
    const after = await panelOpen();
    return {...opened, panel_kept: before === after, single, closed,
      passed: opened.open && opened.inside && opened.rules && before === after &&
        JSON.stringify(single) === JSON.stringify(['Left plot help']) && closed};
  });
  await check('reviewer_tabs', async () => page.evaluate(() => {
    const tabs = [...document.querySelectorAll('#tabs > li > a')].map(link => link.dataset.value);
    return JSON.stringify(tabs) === JSON.stringify(['Plot', 'Trajectory', 'Spatial', 'Viewer']);
  }));
  await check('document_and_accessibility_link', async () => page.evaluate(() =>
    document.title === 'Islet Explorer' && document.documentElement.lang === 'en' &&
    Boolean(document.querySelector('footer a[href="https://accessibility.ufl.edu/"]'))));
  await check('keyboard_skip_link', async () => {
    await page.focus('.skip-link');
    await page.keyboard.press('Enter');
    return page.evaluate(() => document.activeElement.id === 'main-content');
  });

  const tabs = [
    { name: 'Plot', plots: [['#plot-plt', 'plotly'], ['#plot-dist', 'plotly']], downloads: ['#plot-dl_summary'] },
    { name: 'Trajectory', plots: [['#traj-traj_scatter', 'plotly'], ['#traj-paga_layout', 'plotly']], downloads: ['#traj-dl_umap_csv', '#traj-dl_paga_csv'] },
    { name: 'Spatial', plots: [['#spatial-tissue_scatter', 'image']], downloads: ['#spatial-download_spatial'] },
    { name: 'Viewer', plots: [], downloads: [] },
  ];
  let viewerHref;
  for (let index = 0; index < tabs.length; index += 1) {
    const tab = tabs[index];
    if (index > 0) {
      const previous = tabs[index - 1].name;
      const keyboardWorked = await check(`keyboard_tab_${tab.name.toLowerCase()}`, async () => {
        await page.focus(`#tabs a[data-value="${previous}"]`);
        await page.keyboard.press('ArrowRight');
        await page.waitForFunction(name => document.querySelector('#tabs li.active a')?.dataset.value === name,
          { timeout: 3000 }, tab.name);
        return true;
      });
      if (!keyboardWorked) await page.click(`#tabs a[data-value="${tab.name}"]`);
    }
    // The PAGA layout is built only on request (2026-10-14).
    if (tab.name === 'Trajectory') {
      await page.$eval('#traj-paga_show_layout', box => { if (!box.checked) box.click(); });
    }
    for (const [selector, kind] of tab.plots) {
      await check(`render_${selector.slice(1)}`, () => rendered(page, selector, kind));
    }
    for (const selector of tab.downloads) {
      await check(`download_${selector.slice(1)}`, () => download(page, selector));
    }
    if (options.expanded === 'true' && tab.name === 'Spatial') {
      await check('expanded_donor_census_and_method', async () => page.evaluate((additions, excluded, expected) => {
        const select = document.querySelector('#spatial-donor');
        const values = select?.selectize ? Object.keys(select.selectize.options) : [...select.options].map(option => option.value);
        const bare = values.map(value => String(Number(value)));
        return { passed: values.length === expected && additions.every(id => values.includes(id)) &&
          !excluded.some(id => bare.includes(id)) &&
          /endocrine composition|Harmony/.test(document.querySelector('#spatial-embedding_method')?.textContent || ''), donors: values.length };
      }, ADDITIONS, EXCLUDED, EXPANDED_DONORS));
      await check('render_expanded_leiden_umap', () => rendered(page, '#spatial-leiden_umap', 'plotly'));
      // Owner, 2026-10-19: the two UMAPs are one size. Owner review, 2026-10-20: the Neighborhood cards were
      // replaced (insulitis prevalence, immune cells and beta-cell content, immune gradient, donor heatmap, cluster
      // membership); each must draw data.
      await check('spatial_umaps_one_size_and_review_cards_drawn', async () => {
        await rendered(page, '#spatial-spatial_umap_immune', 'plotly');
        const ids = ['insulitis_bars', 'beta_immune', 'immune_gradient', 'enrichment_bars', 'enrichment_heatmap',
          'cluster_composition'];
        for (const id of ids) await rendered(page, `#spatial-${id}`, 'plotly');
        return page.evaluate(ids => {
          const area = id => {
            const r = document.querySelector(`${id} .nsewdrag`)?.getBoundingClientRect();
            return r ? [Math.round(r.width), Math.round(r.height)] : null;
          };
          const traces = id => (document.getElementById(`spatial-${id}`)?.data || [])
            .filter(t => (t.x?.length || t.z?.length || 0) > 0).length;
          const result = { leiden_plot_area: area('#spatial-leiden_umap'), immune_plot_area: area('#spatial-spatial_umap_immune'),
            traces: Object.fromEntries(ids.map(id => [id, traces(id)])),
            old_cards_gone: !document.querySelector('#spatial-infiltration_bars, #spatial-distance_boxplot, #spatial_remove_outliers') };
          const [a, b] = [result.leiden_plot_area, result.immune_plot_area];
          result.passed = Boolean(a && b) && Math.abs(a[0] - b[0]) <= 2 && Math.abs(a[1] - b[1]) <= 2
            && Object.values(result.traces).every(n => n > 0) && result.old_cards_gone;
          return result;
        }, ids);
      });
      for (const donor of ADDITIONS) {
        await check(`expanded_donor_${donor}_tissue`, async () => {
          const previous = await page.$eval('#spatial-tissue_scatter img', image => image.src);
          await page.$eval('#spatial-donor', (element, value) => {
            if (element.selectize) element.selectize.setValue(value);
            else { element.value = value; element.dispatchEvent(new Event('change', { bubbles: true })); }
          }, donor);
          await page.waitForFunction(old => document.querySelector('#spatial-tissue_scatter img')?.src !== old, {}, previous);
          return rendered(page, '#spatial-tissue_scatter', 'image');
        });
      }
    }
    if (tab.name === 'Viewer') {
      if (options.expanded === 'true') {
        await check('expanded_image_inventory_and_ranges', async () => {
          await page.waitForSelector('#viewer-selected_image option');
          return page.evaluate(async (base, additions, excluded, expected) => {
            const select = document.querySelector('#viewer-selected_image');
            const names = select?.selectize ? Object.keys(select.selectize.options) : [...select.options].map(option => option.value);
            const donorOf = value => Number(value.match(/_(\d+)\.ome\.tiff?$/i)?.[1]);
            if (names.some(value => excluded.map(Number).includes(donorOf(value)))) return { passed: false, excluded_image: true };
            let checked = 0;
            for (const id of additions.map(Number)) {
              const name = names.find(value => donorOf(value) === id);
              if (!name) return { passed: false, checked };
              const response = await fetch(new URL(`images/${name}`, base), { headers: { Range: 'bytes=0-15' } });
              if (response.status !== 206) { await response.body?.cancel(); return { passed: false, checked }; }
              if ((await response.arrayBuffer()).byteLength !== 16) return { passed: false, checked };
              checked += 1;
            }
            return { passed: names.length === expected && checked === additions.length, images: names.length, checked };
          }, appUrl.href, ADDITIONS, EXCLUDED, EXPANDED_DONORS);
        });
      }
      await check('viewer_local_image_url', async () => {
        await page.waitForSelector('#viewer-avivator_viewer');
        // renderUI creates the frame before its first navigation completes.
        await page.waitForFunction(() => {
          const frame = document.querySelector('#viewer-avivator_viewer');
          try { return Boolean(new URL(frame?.contentWindow?.location.href).searchParams.get('image_url')); }
          catch { return false; }
        });
        const iframe = await page.$('#viewer-avivator_viewer');
        const frame = await iframe.contentFrame();
        viewerHref = frame.url();
        const image = new URL(viewerHref).searchParams.get('image_url');
        return Boolean(image && new URL(image, appUrl).origin === appUrl.origin);
      });
      if (viewerHref) {
        const iframe = await page.$('#viewer-avivator_viewer');
        const frame = await iframe.contentFrame();
        await check('viewer_range_16_bytes', async () => frame.evaluate(async () => {
          const image = new URL(window.location.href).searchParams.get('image_url');
          const response = await fetch(image, { headers: { Range: 'bytes=0-15' }, credentials: 'same-origin' });
          if (response.status !== 206 || !/^bytes 0-15\/[0-9]+$/.test(response.headers.get('content-range') || '')) {
            await response.body?.cancel();
            return { passed: false, status: response.status, bytes: 0 };
          }
          const body = await response.arrayBuffer();
          return { passed: response.status === 206 && /^bytes 0-15\/[0-9]+$/.test(response.headers.get('content-range') || '') && body.byteLength === 16,
            status: response.status, bytes: body.byteLength };
        }));
        await check('viewer_canvas_rendered', async () => {
          await frame.waitForSelector('canvas', { timeout: 120000 });
          return frame.evaluate(() => [...document.querySelectorAll('canvas')]
            .some(canvas => canvas.width > 0 && canvas.height > 0));
        });
        await check('viewer_control_semantics', async () => frame.evaluate(() => {
          const panels = [...document.querySelectorAll('[role="tabpanel"]')];
          const channels = [...document.querySelectorAll('[role="combobox"]')];
          return Boolean(document.querySelector('#root[role="main"]')) && panels.length > 0 &&
            panels.every(panel => {
              const tab = document.getElementById(panel.getAttribute('aria-labelledby'));
              return tab?.getAttribute('role') === 'tab' && tab.getAttribute('aria-controls') === panel.id;
            }) && channels.some(channel => channel.getBoundingClientRect().width > 0 &&
              channel.getAttribute('aria-label'));
        }));
        await check('accessibility_viewer_canvas', () => accessibility(frame, 'Viewer frame'));
      }
    }
    await check(`accessibility_${tab.name.toLowerCase()}`, () => accessibility(page, tab.name));
    await check(`session_${tab.name.toLowerCase()}`, () => page.evaluate(() =>
      window.Shiny?.shinyapp?.$socket?.readyState === 1 && !document.querySelector('#shiny-disconnected-overlay')));
  }

  if (options.expanded === 'true') {
    await page.click('#tabs a[data-value="Trajectory"]');
    // The trajectory fits every object with ten or more core cells, of any class: the class selector is hidden
    // on this tab and the status line names the cutoff.
    await check('trajectory_fits_ten_cell_objects', async () => {
      await rendered(page, '#traj-traj_scatter', 'plotly');
      return page.evaluate(() => {
        const status = document.querySelector('#traj-traj_status')?.textContent.replace(/,/g, '') || '';
        const m = status.match(/(\d+) total \| (\d+) fitted \(10 or more core cells/);
        const selector = document.querySelector('#object_classes');
        const hidden = !selector || selector.offsetParent === null;
        return { passed: Boolean(m) && Number(m[2]) > 0 && Number(m[2]) < Number(m[1]) && hidden,
                 total: m ? Number(m[1]) : null, fitted: m ? Number(m[2]) : null, selector_hidden: hidden };
      });
    });
    // Singles are counted, not analyzed (owner, 2026-10-18): no Single choice in the object selector, and the
    // Plot tab's per-donor card shows them.
    await page.click('#tabs a[data-value="Plot"]');
    await check('singles_counted_not_analysed', async () => {
      await rendered(page, '#plot-object_counts_plot', 'plotly');
      return page.evaluate(() => {
        const choices = [...document.querySelectorAll('#object_classes input')].map(el => el.value);
        const traces = document.querySelector('#plot-object_counts_plot')?.data?.length || 0;
        return { passed: !choices.includes('EndocrineSingle') && choices.length === 2 && traces > 0, choices, traces };
      });
    });
    // Owner, 2026-10-19: the class selector sits in the Data controls (Plot and Spatial copies); one set of point
    // controls draws points on both Plot cards; the Panc Floyd logo heads the page in the theme's copy.
    await check('plot_controls_selector_and_logo', async () => {
      await rendered(page, '#plot-plt', 'plotly');
      await rendered(page, '#plot-dist', 'plotly');
      return page.evaluate(() => {
        const panelOf = id => document.querySelector(id)?.closest('details.control-panel')
          ?.querySelector('summary')?.textContent.trim();
        const points = id => (document.querySelector(id)?.data || [])
          .some(t => String(t.mode || '').includes('markers') && (t.x?.length || 0) > 50);
        const logos = [...document.querySelectorAll('img.brand-logo')]
          .filter(img => img.offsetParent !== null && img.naturalWidth > 0);
        const result = {
          plot_selector_panel: panelOf('#object_classes'),
          spatial_selector_panel: panelOf('#spatial_object_classes'),
          scatter_points: points('#plot-plt'), distribution_points: points('#plot-dist'),
          separate_dist_point_controls: Boolean(document.querySelector('#plot-dist_pt_size')),
          visible_logos: logos.map(img => img.getAttribute('src')),
        };
        result.passed = result.plot_selector_panel === 'Data controls' && result.spatial_selector_panel === 'Data controls'
          && result.scatter_points && result.distribution_points && !result.separate_dist_point_controls
          && logos.length === 1;
        return result;
      });
    });
    // Owner, 2026-10-19: hovering an object in either Plot card rings it in both, as on the Trajectory tab.
    await check('plot_points_linked_on_hover', async () => {
      const ids = ['plot-plt', 'plot-dist'];
      const rings = () => page.evaluate(ids => ids.map(id =>
        document.getElementById(id).querySelector('.trajectory-link-highlight')?.getAttribute('data-linked-key') || null), ids);
      const key = await page.evaluate(ids => {
        const sets = ids.map(id => new Set((document.getElementById(id).data || [])
          .filter(t => t.meta?.trajectory_role === 'islet').flatMap(t => t.key)));
        return [...sets[0]].find(k => sets[1].has(k)) || null;
      }, ids);
      const result = { key, hovered: {} };
      // Returning to the tab redraws its plots: wait until Shiny has been idle for a second, so an unhover is not
      // sent to a widget that is about to be replaced.
      await page.waitForFunction(() => {
        const busy = document.documentElement.classList.contains('shiny-busy') ||
          document.querySelector('#plot-plt.recalculating, #plot-dist.recalculating');
        const now = Date.now();
        if (busy || !window.__plotIdleSince) window.__plotIdleSince = busy ? null : now;
        return !busy && window.__plotIdleSince && now - window.__plotIdleSince > 1000;
      }, { polling: 200, timeout: 60000 });
      for (const source of ids) {
        for (let attempt = 0; attempt < 20; attempt++) {
          await page.evaluate((id, key) => {
            const gd = document.getElementById(id);
            const curveNumber = gd.data.findIndex(t => t.meta?.trajectory_role === 'islet' && t.key.includes(key));
            const data = gd.data[curveNumber];
            gd.emit('plotly_hover', { points: [{ data, curveNumber, pointNumber: data.key.indexOf(key) }] });
          }, source, key);
          await new Promise(resolve => setTimeout(resolve, 300));
          result.hovered[source] = await rings();
          if (result.hovered[source].every(k => k === key)) break;
        }
        await page.evaluate(id => document.getElementById(id).emit('plotly_unhover', {}), source);
      }
      await new Promise(resolve => setTimeout(resolve, 300));
      result.after_unhover = await rings();
      result.passed = Boolean(key) && ids.every(source => result.hovered[source].every(k => k === key))
        && result.after_unhover.every(k => k === null);
      return result;
    });
  }

  const viewerPage = await browser.newPage();
  await instrument(viewerPage);
  for (const [name, image] of [['missing', null], ['external', 'https://example.invalid/not-a-real-image.ome.tiff']]) {
    await check(`viewer_${name}_url_rejected`, async () => {
      const url = new URL('avivator/index.html', appUrl);
      if (image) url.searchParams.set('image_url', image);
      await viewerPage.goto(url.href, { waitUntil: 'networkidle0' });
      return viewerPage.evaluate(() => Boolean(document.querySelector('#root[role="status"]')) &&
        !document.querySelector('canvas') &&
        !performance.getEntriesByType('resource').some(entry =>
          /\/assets\/index-[^/]+\.js$/.test(new URL(entry.name).pathname)));
    });
  }
  report.checks.push({ name: 'no_external_network_requests', passed: report.externalRequestsBlocked === 0,
    attempts: report.externalRequestsBlocked });
  report.checks.push({ name: 'no_uncaught_browser_errors', passed: report.browserErrors.page === 0,
    count: report.browserErrors.page });
  report.checks.push({ name: 'application_http_responses', passed: report.httpErrors.application === 0,
    failures: report.httpErrors.application });
} catch (error) {
  report.checks.push({ name: 'browser_run_completed', passed: false,
    error: error.name === 'TimeoutError' ? 'timeout' : 'browser_run_failed' });
} finally {
  if (browser) await browser.close();
  report.completedAt = new Date().toISOString();
  report.runState = 'completed';
  report.passed = report.checks.length > 0 && report.checks.every(item => item.passed);
  await saveReport();
  console.log(JSON.stringify({ passed: report.passed,
    failedChecks: report.checks.filter(item => !item.passed).map(item => item.name),
    accessibilityRules: [...new Set(report.accessibility.flatMap(item => item.violations.map(rule => rule.rule)))],
  }));
  process.exitCode = report.passed ? 0 : 1;
}

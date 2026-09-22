#!/usr/bin/env node
// Browser acceptance against an explicitly selected local/staged application.
// Evidence contains checks and rule counts, never page text, scientific rows,
// download bodies, screenshots, console messages, image names, or credentials.
import fs from 'node:fs/promises';
import path from 'node:path';
import { createRequire } from 'node:module';
import puppeteer from 'puppeteer-core';

const require = createRequire(import.meta.url);
const options = {};
for (let index = 2; index < process.argv.length; index += 2) {
  const key = process.argv[index];
  if (!['--url', '--output', '--chrome'].includes(key) || !process.argv[index + 1]) {
    throw new Error('Usage: check-app.mjs --url APP_URL --output REPORT.json [--chrome CHROME_PATH]');
  }
  options[key.slice(2)] = process.argv[index + 1];
}
if (!options.url || !options.output) throw new Error('--url and --output are required');
const appUrl = new URL(options.url);
if (!['http:', 'https:'].includes(appUrl.protocol) || appUrl.username || appUrl.password ||
    appUrl.search || appUrl.hash) throw new Error('Supply an HTTP(S) application URL without credentials or query parameters');
if (!appUrl.pathname.endsWith('/')) appUrl.pathname += '/';

const report = {
  schemaVersion: 1,
  runState: 'running',
  startedAt: new Date().toISOString(),
  target: { origin: appUrl.origin, applicationPath: appUrl.pathname },
  versions: {
    puppeteer: require('puppeteer-core/package.json').version,
    axe: require('axe-core/package.json').version,
  },
  checks: [],
  accessibility: [],
  browserErrors: { page: 0, console: 0 },
  httpErrors: { missingDefaultFavicon: 0, application: 0 },
  externalRequestsBlocked: 0,
  limitations: [
    'Automated checks are a subset of accessibility review; visualizations and assistive-technology use still require manual review.',
    'This run does not establish external Basic Auth, UFRC proxy behavior, or production-container security.',
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
    if (response.status() === 404 && url.origin === appUrl.origin && url.pathname === '/favicon.ico') {
      report.httpErrors.missingDefaultFavicon += 1;
    } else {
      report.httpErrors.application += 1;
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
  await check('reviewer_tabs', async () => page.evaluate(() => {
    const tabs = [...document.querySelectorAll('#tabs > li > a')].map(link => link.dataset.value);
    return JSON.stringify(tabs) === JSON.stringify(['Plot', 'Trajectory', 'Viewer', 'Spatial']);
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
    { name: 'Viewer', plots: [], downloads: [] },
    { name: 'Spatial', plots: [['#spatial-tissue_scatter', 'image']], downloads: ['#spatial-download_spatial'] },
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
    for (const [selector, kind] of tab.plots) {
      await check(`render_${selector.slice(1)}`, () => rendered(page, selector, kind));
    }
    for (const selector of tab.downloads) {
      await check(`download_${selector.slice(1)}`, () => download(page, selector));
    }
    if (tab.name === 'Viewer') {
      await check('viewer_local_image_url', async () => {
        await page.waitForSelector('#viewer-avivator_viewer');
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

#!/usr/bin/env node
// Responsiveness baseline for the reviewer app. Drives one fixed scenario against an
// explicitly selected application URL and records, per step, how long the browser waited
// and where the time went: the WebSocket round trip, the server's compute until its
// output message arrived, the bytes that crossed the socket, and the browser's own draw.
// Run it against several vantage points (workstation container, VM backend through a
// tunnel, the public proxy) and compare medians; see deploy/README.md § "Responsiveness
// baseline". The report holds timings and byte counts only: no page text, no scientific
// rows, no credentials.
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const options = { runs: '3', throttle: 'none', label: '' };
for (let index = 2; index < process.argv.length; index += 2) {
  const key = process.argv[index];
  if (!['--url', '--output', '--chrome', '--auth-file', '--runs', '--throttle', '--label'].includes(key) ||
      !process.argv[index + 1]) {
    throw new Error('Usage: bench-app.mjs --url APP_URL --output REPORT.json [--chrome PATH] [--auth-file JSON] ' +
      '[--runs N] [--throttle none|home|slow] [--label NAME]');
  }
  options[key.slice(2)] = process.argv[index + 1];
}
if (!options.url || !options.output) throw new Error('--url and --output are required');
const appUrl = new URL(options.url);
if (!['http:', 'https:'].includes(appUrl.protocol) || appUrl.username || appUrl.password ||
    appUrl.search || appUrl.hash) throw new Error('Supply an HTTP(S) application URL without credentials or query parameters');
if (!appUrl.pathname.endsWith('/')) appUrl.pathname += '/';
const runs = Number.parseInt(options.runs, 10);
if (!Number.isInteger(runs) || runs < 1 || runs > 20) throw new Error('--runs must be 1..20');
// Emulated link conditions are added on top of the real path (Chrome DevTools network
// emulation): "home" approximates residential broadband, "slow" a poor connection.
const THROTTLES = {
  none: null,
  home: { latency: 40, downloadThroughput: 25e6 / 8, uploadThroughput: 5e6 / 8 },
  slow: { latency: 120, downloadThroughput: 8e6 / 8, uploadThroughput: 2e6 / 8 },
};
if (!(options.throttle in THROTTLES)) throw new Error('--throttle must be none, home or slow');

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
  label: options.label || appUrl.host,
  startedAt: new Date().toISOString(),
  target: { origin: appUrl.origin, applicationPath: appUrl.pathname, basicAuth: Boolean(credentials) },
  throttle: options.throttle,
  runs: [],
  summary: {},
};

// Attached before any page script runs: records Shiny's lifecycle events against the
// page clock so a step can be split into round trip, server work and browser draw.
const PAGE_HOOK = `(() => {
  window.__bench = { events: [], values: [] };
  const wait = () => {
    if (!window.jQuery) { setTimeout(wait, 1); return; }
    const $ = window.jQuery;
    $(document).on('shiny:connected shiny:busy shiny:idle shiny:disconnected', e => {
      window.__bench.events.push({ type: e.type, t: performance.now() });
    });
    $(document).on('shiny:value shiny:error', e => {
      window.__bench.values.push({ type: e.type, name: e.name, t: performance.now() });
    });
  };
  wait();
})();`;

const now = page => page.evaluate(() => performance.now());
const eventsAfter = (page, t0) => page.evaluate(start => ({
  events: window.__bench.events.filter(e => e.t >= start),
  values: window.__bench.values.filter(v => v.t >= start),
}), t0);

function rendered(page, selector, kind, timeout = 180000) {
  return page.waitForFunction((target, outputKind) => {
    const element = document.querySelector(target);
    if (!element || element.classList.contains('shiny-output-error') ||
        element.classList.contains('recalculating')) return false;
    if (outputKind === 'plotly') return Boolean(element._fullLayout && element.data?.length);
    const image = element.querySelector('img');
    return Boolean(image?.complete && image.naturalWidth > 0);
  }, { timeout, polling: 50 }, selector, kind);
}

async function idleAfter(page, t0, timeout = 180000) {
  // A busy cycle that started after t0 and has since returned to idle, with no output
  // still marked as recalculating on the visible tab.
  await page.waitForFunction(start => {
    const events = window.__bench.events.filter(e => e.t >= start);
    const busy = events.find(e => e.type === 'shiny:busy');
    if (!busy) return false;
    const idle = events.find(e => e.type === 'shiny:idle' && e.t > busy.t);
    if (!idle) return false;
    const active = document.querySelector('.tab-pane.active');
    return !active || !active.querySelector('.recalculating');
  }, { timeout, polling: 50 }, t0);
}

class SocketMeter {
  constructor(client) {
    this.bytesIn = 0; this.framesIn = 0; this.bytesOut = 0; this.framesOut = 0; this.sockets = 0;
    // DevTools reports frame payloads after decompression, so the negotiated extension (permessage-deflate
    // or none) is read from the handshake response instead; it tells whether the socket is compressed.
    this.extensions = null;
    client.on('Network.webSocketCreated', () => { this.sockets += 1; });
    client.on('Network.webSocketHandshakeResponseReceived', ({ response }) => {
      const headers = Object.fromEntries(Object.entries(response.headers || {}).map(([k, v]) => [k.toLowerCase(), v]));
      this.extensions = headers['sec-websocket-extensions'] || 'none';
    });
    client.on('Network.webSocketFrameReceived', ({ response }) => {
      this.framesIn += 1; this.bytesIn += response.payloadData?.length || 0;
    });
    client.on('Network.webSocketFrameSent', ({ response }) => {
      this.framesOut += 1; this.bytesOut += response.payloadData?.length || 0;
    });
  }
  mark() { return { bytesIn: this.bytesIn, framesIn: this.framesIn, bytesOut: this.bytesOut, framesOut: this.framesOut }; }
  since(mark) {
    return { wsBytesIn: this.bytesIn - mark.bytesIn, wsFramesIn: this.framesIn - mark.framesIn,
      wsBytesOut: this.bytesOut - mark.bytesOut, wsFramesOut: this.framesOut - mark.framesOut };
  }
}

class HttpMeter {
  constructor(page) {
    this.entries = [];
    page.on('response', response => {
      const url = new URL(response.url());
      if (url.origin !== appUrl.origin) return;
      const headers = response.headers();
      const length = Number.parseInt(headers['content-length'] || '0', 10) || 0;
      this.entries.push({ path: url.pathname, status: response.status(), length,
        encoding: headers['content-encoding'] || '', t: Date.now() });
    });
  }
  mark() { return this.entries.length; }
  since(mark, predicate = () => true) {
    const slice = this.entries.slice(mark).filter(predicate);
    return { requests: slice.length, bytes: slice.reduce((sum, e) => sum + e.length, 0) };
  }
}

// One scenario step: run `action`, wait with `settle`, and describe the wait from the
// page's own event log and the socket meter.
async function step(context, name, action, settle, extra = {}) {
  const { page, socket } = context;
  const t0 = await now(page);
  const socketMark = socket.mark();
  const wall0 = Date.now();
  let ok = true; let error;
  try {
    await action();
    await settle(t0);
  } catch (caught) {
    ok = false; error = String(caught?.message || caught).split('\n')[0].slice(0, 200);
  }
  const t1 = await now(page);
  const { events, values } = await eventsAfter(page, t0);
  const busy = events.find(e => e.type === 'shiny:busy');
  const idle = busy ? events.find(e => e.type === 'shiny:idle' && e.t > busy.t) : null;
  const firstValue = values.find(v => v.type === 'shiny:value');
  const lastValue = [...values].reverse().find(v => v.type === 'shiny:value');
  const result = {
    name, ok, error,
    totalMs: Math.round(t1 - t0),
    wallMs: Date.now() - wall0,
    // time from the action until the server reported busy: one round trip plus queueing
    busyAfterMs: busy ? Math.round(busy.t - t0) : null,
    // time from busy until the first and the last output value arrived: server compute
    // plus the transfer of each message
    firstValueMs: busy && firstValue ? Math.round(firstValue.t - busy.t) : null,
    lastValueMs: busy && lastValue ? Math.round(lastValue.t - busy.t) : null,
    idleAfterMs: idle ? Math.round(idle.t - t0) : null,
    // what the browser still had to do after the last value arrived (plotly/image draw)
    drawAfterLastValueMs: lastValue ? Math.round(t1 - lastValue.t) : null,
    valuesReceived: values.filter(v => v.type === 'shiny:value').map(v => v.name),
    errorsReceived: values.filter(v => v.type === 'shiny:error').map(v => v.name),
    ...socket.since(socketMark),
    ...extra,
  };
  context.run.steps.push(result);
  return result;
}

async function setSelectize(page, selector, pick) {
  // pick: function(values, current) -> value to choose; returns the chosen value
  return page.evaluate((target, pickSource) => {
    const element = document.querySelector(target);
    if (!element) throw new Error(`${target} not found`);
    const choose = new Function('values', 'current', `return (${pickSource})(values, current);`);
    if (element.selectize) {
      const values = Object.keys(element.selectize.options);
      const value = choose(values, element.selectize.getValue());
      element.selectize.setValue(value);
      return value;
    }
    const values = [...element.options].map(option => option.value);
    const value = choose(values, element.value);
    element.value = value;
    element.dispatchEvent(new Event('change', { bubbles: true }));
    return value;
  }, selector, pick.toString());
}

async function navigationTiming(page) {
  return page.evaluate(() => {
    const nav = performance.getEntriesByType('navigation')[0];
    const resources = performance.getEntriesByType('resource');
    const transfer = entry => entry.transferSize || 0;
    const scripts = resources.filter(r => r.initiatorType === 'script' || r.name.endsWith('.js'));
    const styles = resources.filter(r => r.initiatorType === 'link' || r.name.endsWith('.css'));
    const largest = [...resources].sort((a, b) => transfer(b) - transfer(a)).slice(0, 8)
      .map(r => ({ path: new URL(r.name).pathname.replace(/\/\d+(?=\/|$)/g, '/{id}'),
        transferBytes: transfer(r), decodedBytes: r.decodedBodySize || 0, ms: Math.round(r.duration) }));
    return {
      ttfbMs: nav ? Math.round(nav.responseStart) : null,
      domContentLoadedMs: nav ? Math.round(nav.domContentLoadedEventEnd) : null,
      loadMs: nav ? Math.round(nav.loadEventEnd) : null,
      documentTransferBytes: nav ? transfer(nav) : null,
      documentDecodedBytes: nav ? nav.decodedBodySize : null,
      resources: resources.length,
      resourceTransferBytes: resources.reduce((sum, r) => sum + transfer(r), 0),
      resourceDecodedBytes: resources.reduce((sum, r) => sum + (r.decodedBodySize || 0), 0),
      scripts: scripts.length, scriptTransferBytes: scripts.reduce((sum, r) => sum + transfer(r), 0),
      styles: styles.length,
      cachedResources: resources.filter(r => r.transferSize === 0 && r.decodedBodySize > 0).length,
      largest,
    };
  });
}

async function socketTransport(page) {
  return page.evaluate(() => {
    // Shiny Server's client wraps the SockJS connection in several layers; look through
    // them (bounded) for the transport name and the connection URL.
    const seen = new Set(); const queue = [window.Shiny?.shinyapp?.$socket];
    let transport = null; let url = null;
    while (queue.length && seen.size < 300 && !(transport && url)) {
      const object = queue.shift();
      if (!object || typeof object !== 'object' || seen.has(object)) continue;
      seen.add(object);
      if (!transport && typeof object.transport === 'string') transport = object.transport;
      if (!url && typeof object.url === 'string' && /^(wss?|https?):/.test(object.url)) url = object.url;
      for (const key of Object.keys(object)) {
        const value = object[key];
        if (value && typeof value === 'object') queue.push(value);
      }
    }
    let scheme = null;
    try { scheme = url ? new URL(url).protocol : null; } catch { scheme = null; }
    return { transport, scheme };
  });
}

async function scenario(browser, runIndex, warm) {
  const context = warm ? null : await browser.createBrowserContext();
  const page = warm ? await browser.newPage() : await context.newPage();
  page.setDefaultTimeout(180000);
  if (credentials) await page.authenticate(credentials);
  await page.setViewport({ width: 1600, height: 1100, deviceScaleFactor: 1 });
  await page.evaluateOnNewDocument(PAGE_HOOK);
  const client = await page.createCDPSession();
  await client.send('Network.enable');
  if (THROTTLES[options.throttle]) {
    await client.send('Network.emulateNetworkConditions', { offline: false, ...THROTTLES[options.throttle] });
  }
  const socket = new SocketMeter(client);
  const http = new HttpMeter(page);
  const run = { index: runIndex, cache: warm ? 'warm' : 'cold', startedAt: new Date().toISOString(), steps: [] };
  const ctx = { page, socket, run };
  const errors = { page: 0, console: 0 };
  page.on('pageerror', () => { errors.page += 1; });
  page.on('console', message => { if (message.type() === 'error') errors.console += 1; });

  // 1. Page load until Shiny's socket is open, then until the Plot tab's two plots are drawn.
  const httpMark = http.mark();
  const t0 = await (async () => { const wall = Date.now(); return wall; })();
  await page.goto(appUrl.href, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
  const connectedAt = await page.evaluate(() => {
    const e = window.__bench.events.find(ev => ev.type === 'shiny:connected');
    return e ? e.t : performance.now();
  });
  run.load = { connectedMs: Math.round(connectedAt), wallToConnectedMs: Date.now() - t0 };
  const socketMark0 = socket.mark();
  await rendered(page, '#plot-plt', 'plotly');
  await rendered(page, '#plot-dist', 'plotly');
  const plotsAt = await now(page);
  const firstIdle = await page.evaluate(start => {
    const e = window.__bench.events.find(ev => ev.type === 'shiny:idle' && ev.t > start);
    return e ? e.t : null;
  }, connectedAt);
  run.load.plotTabDrawnMs = Math.round(plotsAt);
  run.load.firstIdleMs = firstIdle ? Math.round(firstIdle) : null;
  run.load.plotTabDrawnAfterConnectMs = Math.round(plotsAt - connectedAt);
  Object.assign(run.load, socket.since(socketMark0));
  run.load.wsSockets = socket.sockets;
  const initial = (await eventsAfter(page, 0)).values.filter(v => v.type === 'shiny:value');
  run.load.initialValues = initial.length;
  // Arrival of each output after the socket opened: which outputs the first screen waits for.
  run.load.valueTimeline = initial.map(v => ({ name: v.name, afterConnectMs: Math.round(v.t - connectedAt) }));
  run.load.navigation = await navigationTiming(page);
  run.load.http = http.since(httpMark);
  run.load.socket = { ...(await socketTransport(page) || {}), extensions: socket.extensions };
  if (warm) {
    await page.close();
    return { run, errors };
  }

  // 2. Change the plotted measure (both Plot cards re-render).
  await step(ctx, 'plot_change_measure', async () => {
    const chosen = await setSelectize(page, '#plot-which', (values, current) => values.find(v => v !== current) || current);
    run.steps.chosen = chosen;
  }, async start => {
    await idleAfter(page, start);
    await rendered(page, '#plot-plt', 'plotly');
    await rendered(page, '#plot-dist', 'plotly');
  });
  await step(ctx, 'plot_change_measure_again', async () => {
    await setSelectize(page, '#plot-which', (values, current) => values.filter(v => v !== current).slice(-1)[0] || current);
  }, async start => {
    await idleAfter(page, start);
    await rendered(page, '#plot-plt', 'plotly');
    await rendered(page, '#plot-dist', 'plotly');
  });

  // 3. Drill-down click on the first keyed point: loads the case geometry and the
  //    object's single cells, then draws the segmentation map and the composition bars.
  await step(ctx, 'drilldown_click', async () => {
    const sent = await page.evaluate(() => {
      const gd = document.querySelector('#plot-plt');
      const trace = (gd?.data || []).find(t => Array.isArray(t.key) && t.key.length);
      if (!trace) return false;
      const key = trace.key[Math.floor(trace.key.length / 2)];
      window.Shiny.setInputValue('plot-linked_click', { key, source: 'plot-plt', nonce: Date.now() }, { priority: 'event' });
      return true;
    });
    if (!sent) throw new Error('no keyed trace on #plot-plt');
  }, async start => {
    await idleAfter(page, start);
    await page.waitForFunction(() => {
      const ids = ['#islet_drilldown_view', '#islet_segmentation_view'];
      return ids.some(id => {
        const element = document.querySelector(id);
        const image = element?.querySelector('img');
        return image?.complete && image.naturalWidth > 0 && !element.classList.contains('recalculating');
      });
    }, { polling: 50 });
  });

  // 4. Trajectory tab: scatter, then the PAGA layout (Python bridge).
  await step(ctx, 'tab_trajectory', async () => page.click('#tabs a[data-value="Trajectory"]'),
    async () => { await rendered(page, '#traj-traj_scatter', 'plotly'); });
  await step(ctx, 'trajectory_paga_layout', async () => {
    await page.$eval('#traj-paga_show_layout', box => { if (!box.checked) box.click(); });
  }, async () => { await rendered(page, '#traj-paga_layout', 'plotly'); });

  // 5. Spatial tab: tissue scatter (ggplot PNG) and the two UMAPs; then another donor.
  await step(ctx, 'tab_spatial', async () => page.click('#tabs a[data-value="Spatial"]'), async () => {
    await rendered(page, '#spatial-tissue_scatter', 'image');
    await rendered(page, '#spatial-leiden_umap', 'plotly');
    await rendered(page, '#spatial-spatial_umap_immune', 'plotly');
  });
  await step(ctx, 'spatial_change_donor', async () => {
    const previous = await page.$eval('#spatial-tissue_scatter img', image => image.src);
    run.steps.previousTissue = previous;
    await setSelectize(page, '#spatial-donor', (values, current) => values.find(v => v !== current) || current);
  }, async start => {
    await page.waitForFunction(old => document.querySelector('#spatial-tissue_scatter img')?.src !== old, { polling: 50 },
      run.steps.previousTissue);
    await idleAfter(page, start);
    await rendered(page, '#spatial-tissue_scatter', 'image');
  });
  delete run.steps.previousTissue; delete run.steps.chosen;

  // 6. Viewer tab: the Avivator frame, its first tiles, and the Range traffic until quiet.
  const viewerHttpMark = http.mark();
  const viewer = await step(ctx, 'tab_viewer_first_tiles', async () => page.click('#tabs a[data-value="Viewer"]'), async () => {
    await page.waitForFunction(() => {
      const frame = document.querySelector('#viewer-avivator_viewer');
      try { return Boolean(new URL(frame?.contentWindow?.location.href).searchParams.get('image_url')); }
      catch { return false; }
    }, { polling: 50 });
    await page.waitForFunction(() => {
      const frame = document.querySelector('#viewer-avivator_viewer');
      return Boolean(frame?.contentDocument?.querySelector('canvas'));
    }, { polling: 50 });
    await page.waitForResponse(response => response.status() === 206, { timeout: 120000 });
  });
  const quietStart = Date.now();
  try { await page.waitForNetworkIdle({ idleTime: 2000, timeout: 90000 }); } catch { /* still streaming */ }
  viewer.rangeQuietAfterMs = Date.now() - quietStart - 2000;
  Object.assign(viewer, (() => {
    const ranges = http.since(viewerHttpMark, e => e.status === 206);
    return { rangeRequests: ranges.requests, rangeBytes: ranges.bytes };
  })());

  run.errors = errors;
  await page.close();
  await context.close();
  return { run, errors };
}

function median(values) {
  const sorted = values.filter(v => typeof v === 'number' && Number.isFinite(v)).sort((a, b) => a - b);
  if (!sorted.length) return null;
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : Math.round((sorted[mid - 1] + sorted[mid]) / 2);
}

function summarize() {
  const cold = report.runs.filter(r => r.cache === 'cold');
  const warm = report.runs.filter(r => r.cache === 'warm');
  const loadKeys = ['connectedMs', 'plotTabDrawnMs', 'plotTabDrawnAfterConnectMs', 'firstIdleMs', 'wsBytesIn'];
  const summary = { coldRuns: cold.length, warmRuns: warm.length, load: {}, warmLoad: {}, steps: {} };
  for (const key of loadKeys) summary.load[key] = median(cold.map(r => r.load?.[key]));
  for (const key of ['ttfbMs', 'domContentLoadedMs', 'loadMs', 'resources', 'resourceTransferBytes', 'resourceDecodedBytes']) {
    summary.load[key] = median(cold.map(r => r.load?.navigation?.[key]));
  }
  for (const key of ['connectedMs', 'plotTabDrawnMs']) summary.warmLoad[key] = median(warm.map(r => r.load?.[key]));
  summary.warmLoad.cachedResources = median(warm.map(r => r.load?.navigation?.cachedResources));
  summary.load.socket = cold[0]?.load?.socket || null;
  const names = [...new Set(cold.flatMap(r => r.steps.map(s => s.name)))];
  for (const name of names) {
    const steps = cold.map(r => r.steps.find(s => s.name === name)).filter(Boolean);
    summary.steps[name] = {
      ok: steps.filter(s => s.ok).length, of: steps.length,
      totalMs: median(steps.map(s => s.totalMs)),
      busyAfterMs: median(steps.map(s => s.busyAfterMs)),
      lastValueMs: median(steps.map(s => s.lastValueMs)),
      drawAfterLastValueMs: median(steps.map(s => s.drawAfterLastValueMs)),
      wsBytesIn: median(steps.map(s => s.wsBytesIn)),
      rangeBytes: median(steps.map(s => s.rangeBytes)),
      rangeRequests: median(steps.map(s => s.rangeRequests)),
    };
  }
  report.summary = summary;
}

async function save() {
  await fs.mkdir(path.dirname(path.resolve(options.output)), { recursive: true });
  await fs.writeFile(options.output, `${JSON.stringify(report, null, 2)}\n`, { mode: 0o600 });
}

const browser = await puppeteer.launch({
  executablePath: options.chrome || '/usr/bin/google-chrome',
  headless: true,
  args: ['--disable-dev-shm-usage', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'],
});
report.chrome = await browser.version();
try {
  for (let index = 0; index < runs; index += 1) {
    const { run } = await scenario(browser, index, false);
    report.runs.push(run);
    summarize();
    await save();
    console.error(`[${report.label}] cold run ${index + 1}/${runs}: connected ${run.load.connectedMs} ms, ` +
      `plots ${run.load.plotTabDrawnMs} ms, steps ${run.steps.map(s => `${s.name}=${s.ok ? s.totalMs : 'FAIL'}`).join(' ')}`);
  }
  // Warm visit: the default browser context keeps its cache between these navigations.
  for (let index = 0; index < Math.min(runs, 2); index += 1) {
    const { run } = await scenario(browser, index, true);
    report.runs.push(run);
    console.error(`[${report.label}] warm visit ${index + 1}: connected ${run.load.connectedMs} ms, plots ${run.load.plotTabDrawnMs} ms`);
  }
  summarize();
  report.finishedAt = new Date().toISOString();
  await save();
} finally {
  await browser.close();
}
console.log(JSON.stringify(report.summary, null, 2));

// Focused interaction regression check; use an isolated local Shiny worker.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import puppeteer from 'puppeteer-core';

const arg = name => process.argv[process.argv.indexOf(name) + 1];
const url = process.argv.includes('--url') ? arg('--url') : 'http://127.0.0.1:3841/';
const output = process.argv.includes('--output') ? arg('--output') : '/tmp/trajectory-links-check.json';
const report = {checks: [], errors: [], hoverMs: []};
const browser = await puppeteer.launch({executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-dev-shm-usage']});
const page = await browser.newPage();
await page.setViewport({width: 1680, height: 1200});
page.setDefaultTimeout(180000);
page.on('pageerror', error => report.errors.push(String(error)));
const ids = ['traj-traj_scatter', 'traj-traj_umap_donor', 'traj-traj_umap_feature'];
async function check(name, fn) { await fn(); report.checks.push(name); console.log('PASS', name); }
async function ready() {
  await page.waitForFunction(ids => ids.every(id => {
    const gd = document.getElementById(id);
    return gd && !gd.classList.contains('recalculating') &&
      gd.data?.some(t => t.meta?.trajectory_role === 'islet' && t.x.length) &&
      gd.querySelector('.trajectory-link-highlight');
  }) && !document.documentElement.classList.contains('shiny-busy'), {}, ids);
}
async function choose(id, value) {
  const current = await page.$eval('#' + id, el => el.selectize ? el.selectize.getValue() : el.value);
  if (current === value) { await ready(); return; }
  await page.evaluate(() => { window.previousLinkOverlay = document.querySelector('#traj-traj_scatter .trajectory-link-highlight'); });
  await page.$eval('#' + id, (el, value) => {
    if (el.selectize) el.selectize.setValue(value);
    else { el.value = value; el.dispatchEvent(new Event('change', {bubbles: true})); }
  }, value);
  await page.waitForFunction(() => document.querySelector('#traj-traj_scatter .trajectory-link-highlight') !== window.previousLinkOverlay);
  await ready();
}
async function highlighted(key) {
  await page.waitForFunction((ids, key) => ids.every(id => {
    const gd = document.getElementById(id), overlay = gd.querySelector('.trajectory-link-highlight');
    if (!key) return !overlay?.getAttribute('data-linked-key');
    const trace = gd.data.find(t => t.meta?.trajectory_role === 'islet' && t.key.includes(key));
    if (!trace) return !overlay?.getAttribute('data-linked-key');
    const i = trace.key.indexOf(key);
    return overlay?.getAttribute('data-linked-key') === key &&
      +overlay.getAttribute('data-point-x') === trace.x[i] && +overlay.getAttribute('data-point-y') === trace.y[i];
  }), {}, ids, key);
}
async function emit(id, event, key) {
  await page.evaluate((id, event, key) => {
    const gd = document.getElementById(id);
    const curveNumber = gd.data.findIndex(t => t.meta?.trajectory_role === 'islet' && t.key.includes(key));
    const data = gd.data[curveNumber];
    gd.emit(event, {points: [{curveNumber, pointNumber: data.key.indexOf(key), data}]});
  }, id, event, key);
}
async function ranges() {
  return page.evaluate(ids => ids.map(id => {
    const layout = document.getElementById(id)._fullLayout;
    return [layout.xaxis.range, layout.yaxis.range];
  }), ids);
}
async function ringAlignment(key) {
  const offsets = await page.evaluate((ids, key) => ids.map(id => {
    const gd = document.getElementById(id), curve = gd.data.findIndex(t => t.key?.includes(key));
    const i = gd.data[curve].key.indexOf(key);
    const trace = gd.querySelector('.scatterlayer .trace' + gd._fullData[curve].uid);
    const point = [...trace.querySelectorAll('path.point')].find(p => p.__data__.i === i);
    const p = point.getBoundingClientRect(), r = gd.querySelector('.trajectory-link-highlight circle').getBoundingClientRect();
    return {id, dx: Math.abs(p.x + p.width / 2 - r.x - r.width / 2),
      dy: Math.abs(p.y + p.height / 2 - r.y - r.height / 2)};
  }), ids, key);
  assert.ok(offsets.every(o => o.dx < 1 && o.dy < 1), JSON.stringify(offsets));
}
try {
  await page.goto(url, {waitUntil: 'domcontentloaded'});
  await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
  await page.click('a[data-value="Trajectory"]');
  await ready();
  const keys = await page.$eval('#' + ids[0], gd => gd.data.filter(t => t.meta?.trajectory_role === 'islet').flatMap(t => t.key));
  const key = keys.find(k => /^78\|Islet_\d{4,}$/.test(k));
  assert.ok(key, 'expanded-cohort short donor ID and long islet ID must be present');
  const other = keys.find(k => !k.startsWith('78|'));
  await page.evaluate(() => {
    window.linkInputs = [];
    window.jQuery(document).on('shiny:inputchanged.trajectoryCheck', event => { if (!event.isDefaultPrevented()) window.linkInputs.push(event.name); });
  });
  await check('both UMAPs and scatter link exact IDs in every direction without server hover traffic', async () => {
    const before = await ranges();
    for (const id of ids) {
      await page.evaluate(() => { window.linkInputs = []; });
      const elapsed = await page.evaluate(async (id, key, ids) => {
        const gd = document.getElementById(id);
        const curveNumber = gd.data.findIndex(t => t.meta?.trajectory_role === 'islet' && t.key.includes(key));
        const data = gd.data[curveNumber], pointNumber = data.key.indexOf(key);
        const start = performance.now();
        gd.emit('plotly_hover', {points: [{data, curveNumber, pointNumber}]});
        await new Promise((resolve, reject) => {
          function poll() {
            const matched = ids.every(id => {
              const graph = document.getElementById(id), trace = graph.data.find(t => t.meta?.trajectory_role === 'islet' && t.key.includes(key));
              const i = trace.key.indexOf(key);
              const overlay = graph.querySelector('.trajectory-link-highlight');
              return overlay?.getAttribute('data-linked-key') === key && +overlay.getAttribute('data-point-x') === trace.x[i] && +overlay.getAttribute('data-point-y') === trace.y[i];
            });
            if (matched) requestAnimationFrame(resolve);
            else if (performance.now() - start > 5000) reject(new Error('Highlight timeout'));
            else requestAnimationFrame(poll);
          }
          requestAnimationFrame(poll);
        });
        return performance.now() - start;
      }, id, key, ids);
      report.hoverMs.push(elapsed);
      assert.ok(elapsed < 100, 'linked hover should paint within 100 ms');
      assert.deepEqual(await page.evaluate(() => window.linkInputs.filter(n => !n.startsWith('.clientdata'))), []);
      await page.$eval('#' + id, gd => gd.emit('plotly_unhover', {}));
      await highlighted(null);
    }
    assert.deepEqual(await ranges(), before, 'hover must not change plot ranges');
  });
  await check('click pins and opens the exact islet; hover previews then restores it', async () => {
    await emit(ids[1], 'plotly_click', key);
    await highlighted(key);
    await page.waitForFunction(() => document.querySelector('#traj-segmentation_viewer_panel')?.innerText.includes('Islet_'));
    const text = await page.$eval('#traj-segmentation_viewer_panel', el => el.innerText);
    assert.ok(text.includes(key.split('|')[1]));
    await page.waitForSelector('#islet_segmentation_view img');
    await emit(ids[2], 'plotly_hover', other);
    await highlighted(other);
    await page.$eval('#' + ids[2], gd => gd.emit('plotly_unhover', {}));
    await highlighted(key);
    await ringAlignment(key);
  });
  await check('rings follow zoom and resize; tab switching restores the pin', async () => {
    const before = await ranges();
    await page.evaluate(async (id, key) => {
      const gd = document.getElementById(id), t = gd.data.find(t => t.key?.includes(key));
      const i = t.key.indexOf(key), xa = gd._fullLayout.xaxis, ya = gd._fullLayout.yaxis;
      const dx = Math.abs(xa.range[1] - xa.range[0]) / 4, dy = Math.abs(ya.range[1] - ya.range[0]) / 4;
      await Plotly.relayout(gd, {'xaxis.range': [t.x[i] - dx, t.x[i] + dx], 'yaxis.range': [t.y[i] - dy, t.y[i] + dy]});
    }, ids[1], key);
    await page.setViewport({width: 1500, height: 1100});
    await page.waitForFunction(ids => ids.every(id => {
      const gd = document.getElementById(id), ring = gd.querySelector('.trajectory-link-highlight');
      const circle = ring.querySelector('circle');
      return Math.abs(+circle.getAttribute('cx') - (gd._fullLayout.xaxis._offset + gd._fullLayout.xaxis.l2p(+ring.getAttribute('data-point-x')))) < 0.5 &&
        Math.abs(+circle.getAttribute('cy') - (gd._fullLayout.yaxis._offset + gd._fullLayout.yaxis.l2p(+ring.getAttribute('data-point-y')))) < 0.5;
    }), {}, ids);
    await page.evaluate(async (id, range) => {
      await Plotly.relayout(document.getElementById(id), {'xaxis.range': range[0], 'yaxis.range': range[1]});
    }, ids[1], before[1]);
    await page.setViewport({width: 1680, height: 1200});
    await page.click('a[data-value="Plot"]');
    await highlighted(null);
    await page.click('a[data-value="Trajectory"]');
    await ready(); await highlighted(key);
    await ringAlignment(key);
    const coveredDownloads = await page.evaluate(() => ['traj-dl_umap_donor', 'traj-dl_umap_feature'].filter(id => {
      const el = document.getElementById(id), r = el.getBoundingClientRect();
      return !el.contains(document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2));
    }));
    assert.deepEqual(coveredDownloads, [], 'UMAP canvases must not cover their download buttons');
    await page.screenshot({path: output.replace(/\.json$/, '.linked.png'), fullPage: true});
  });
  await check('a missing UMAP coordinate leaves the other plots correctly linked', async () => {
    await page.evaluate((id, key) => {
      const gd = document.getElementById(id), t = gd.data.find(t => t.key?.includes(key)), i = t.key.indexOf(key);
      window.missingLinkPoint = {t, i, y: t.y[i]};
      t.y[i] = null; gd.emit('plotly_afterplot', {});
    }, ids[1], key);
    await page.waitForFunction(id => !document.querySelector('#' + id + ' .trajectory-link-highlight').getAttribute('data-linked-key'), {}, ids[1]);
    await emit(ids[0], 'plotly_hover', key);
    await page.waitForFunction((id, key) => document.querySelector('#' + id + ' .trajectory-link-highlight').getAttribute('data-linked-key') === key, {}, ids[2], key);
    await page.evaluate(id => {
      const {t, i, y} = window.missingLinkPoint; t.y[i] = y;
      document.getElementById(id).emit('plotly_afterplot', {});
    }, ids[1]);
    await page.$eval('#' + ids[0], gd => gd.emit('plotly_unhover', {}));
    await highlighted(key);
  });
  await check('feature, scale, size and color changes preserve pin and point identity', async () => {
    await choose('traj-pt_scale', 'raw'); await highlighted(key);
    await choose('traj-traj_y_mode', 'marker');
    await page.waitForSelector('#traj-traj_feature');
    await choose('traj-traj_feature', 'INS'); await highlighted(key);
    await choose('traj-traj_color_by', 'donor_id'); await highlighted(key);
    await choose('traj-traj_color_by', 'leiden_0.3'); await highlighted(key);
    await choose('traj-traj_point_size', 'uniform'); await highlighted(key);
    assert.match(await page.$eval('#traj-traj_umap_feature_title', el => el.textContent), /^UMAP: INS\b/);
  });
  await check('UMAP PNG downloads retain publication dimensions', async () => {
    for (const id of ['traj-dl_umap_donor', 'traj-dl_umap_feature']) {
      const header = await page.$eval('#' + id, async el => Array.from(new Uint8Array(await (await fetch(el.href)).arrayBuffer()).slice(0, 24)));
      const bytes = Buffer.from(header);
      assert.equal(bytes.subarray(1, 4).toString(), 'PNG');
      assert.equal(bytes.readUInt32BE(16), 3000); assert.equal(bytes.readUInt32BE(20), 2100);
    }
  });
  await check('clear and close remove pin and details', async () => {
    await page.click('#traj-clear_link_selection'); await highlighted(null);
    await emit(ids[2], 'plotly_click', other); await highlighted(other);
    await page.waitForSelector('#traj-clear_segmentation');
    await page.click('#traj-clear_segmentation'); await highlighted(null);
  });
  await check('physical click opens details once and dragging does not select', async () => {
    const panel = await page.$('#' + ids[0]);
    await panel.evaluate(el => el.scrollIntoView({block: 'center'}));
    const target = await page.$eval('#' + ids[0], gd => {
      const t = gd.data.find(t => t.meta?.trajectory_role === 'islet'), rect = gd.querySelector('.svg-container').getBoundingClientRect();
      const x = gd._fullLayout.xaxis, y = gd._fullLayout.yaxis;
      const i = t.x.findIndex((v, i) => {
        const px = x.l2p(v), py = y.l2p(t.y[i]);
        return px > 20 && px < x._length - 20 && py > 20 && py < y._length - 20;
      });
      return {x: rect.x + x._offset + x.l2p(t.x[i]), y: rect.y + y._offset + y.l2p(t.y[i])};
    });
    await page.mouse.move(target.x, target.y);
    await page.waitForFunction(id => document.getElementById(id)._hoverdata?.some(p => p.data.meta?.trajectory_role === 'islet'), {}, ids[0]);
    const picked = await page.$eval('#' + ids[0], gd => {const p=gd._hoverdata.find(p=>p.data.meta?.trajectory_role==='islet');return p.data.key[p.pointNumber];});
    await page.evaluate(() => { window.linkInputs = []; });
    await page.mouse.click(target.x, target.y);
    await highlighted(picked);
    const clicks = await page.evaluate(() => window.linkInputs.filter(n => n === 'traj-linked_click').length);
    assert.equal(clicks, 1);
    await page.mouse.down(); await page.mouse.move(target.x + 40, target.y + 25, {steps: 5}); await page.mouse.up();
    assert.equal(await page.evaluate(() => window.linkInputs.filter(n => n === 'traj-linked_click').length), 1);
  });
  await check('range recalculation clears a filtered-out pin', async () => {
    // Pick the smallest visible islet, then exclude it from the fitted cohort.
    await page.click('#traj-clear_link_selection');
    const csvText = await page.$eval('#traj-dl_umap_csv', async el => (await fetch(el.href)).text());
    // The export uses ordinary CSV quoting; only keys and numeric counts matter.
    const rows = csvText.split('\n').filter(s => s && !s.startsWith('#')).map(s => s.split(',').map(x => x.replace(/^"|"$/g, '')));
    const header = rows.shift(), idx = name => header.indexOf(name);
    const sorted = rows.filter(r => Number.isFinite(+r[idx('total_cells')])).sort((a,b) => +a[idx('total_cells')] - +b[idx('total_cells')]);
    const row = sorted[0], selected = row[idx('case_id')] + '|' + row[idx('islet_key')];
    await emit(ids[1], 'plotly_click', selected); await highlighted(selected);
    await page.$eval('#traj-cell_range', (el, lower) => {
      window.jQuery(el).data('ionRangeSlider').update({from: lower}); window.jQuery(el).trigger('change');
    }, +row[idx('total_cells')] + 5);
    // Slider inputs are debounced; wait for R to receive the new range before
    // pressing Recalculate, as a user would after moving the slider.
    await page.waitForFunction(() => document.querySelector('#traj-traj_range_pending')?.innerText.includes('Range not applied'));
    await page.click('#traj-recalculate_trajectory'); await ready(); await highlighted(null);
  });
  await check('no visible output errors or browser exceptions', async () => {
    const errors = await page.$$eval('.shiny-output-error', elements => elements.filter(e=>e.getClientRects().length).map(e=>e.textContent));
    assert.deepEqual(errors, []); assert.deepEqual(report.errors, []);
  });
  await page.screenshot({path: output.replace(/\.json$/, '.png'), fullPage: true});
  report.passed = true;
} catch (error) {
  report.failure = error.stack || String(error); report.passed = false;
  await page.screenshot({path: output.replace(/\.json$/, '.failed.png'), fullPage: true}).catch(() => {});
  process.exitCode = 1;
} finally {
  await browser.close();
  await fs.writeFile(output, JSON.stringify(report, null, 2) + '\n', {mode: 0o600});
  console.log(JSON.stringify(report, null, 2));
}

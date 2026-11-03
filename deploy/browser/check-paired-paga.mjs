#!/usr/bin/env node
// Focused local PAGA acceptance. Inputs identify the caller's source/runtime;
// they are provenance supplied by the caller, not inferred from the web page.
import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import fs from 'node:fs/promises';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const options = {};
const allowed = ['url', 'output', 'shots', 'source', 'runtime-kind', 'image', 'chrome', 'expected-plotly'];
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  assert(allowed.includes(key) && process.argv[i + 1], 'Unexpected or incomplete argument');
  options[key] = process.argv[i + 1];
}
for (const key of ['url', 'output', 'shots', 'source', 'runtime-kind']) assert(options[key], `--${key} required`);
const url = new URL(options.url);
assert(['127.0.0.1', 'localhost'].includes(url.hostname) && url.protocol === 'http:', 'Loopback HTTP preview required');
assert(!url.username && !url.password && !url.search && !url.hash, 'Unadorned preview URL required');
assert(['source-overlay', 'native-image'].includes(options['runtime-kind']), 'Explicit runtime kind required');
const ids = ['traj-paga_layout', 'traj-paga_layout_right'];
const controls = ['traj-paga_color_by', 'traj-paga_color_right'];
const modes = ['pseudotime', 'donor_status', 'leiden', 'case_id', 'feature'];
const labels = ['Pseudotime', 'Donor Status', 'Leiden Cluster', 'Donor', 'Selected Feature'];
const hash = value => crypto.createHash('sha256').update(JSON.stringify(value)).digest('hex');
const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
const report = {
  schema: 'islet.paired_paga_ui/1', started_utc: new Date().toISOString(), status: 'running',
  target: {url: url.href, source: options.source, runtime_kind: options['runtime-kind'], image: options.image || null,
    identity_basis: 'Caller-supplied identities; bind to image/container inspection separately.'},
  checks: [], cases: [], screenshots: [], page_errors: [], application_http_errors: [],
  limitations: ['Focused browser checks do not rerun scientific fitting, source-data validation, security or deployment acceptance.',
    'Screenshot capture is evidence for separate visual review; this program does not claim human visual inspection.'],
};
await fs.mkdir(path.dirname(path.resolve(options.output)), {recursive: true});
await fs.mkdir(path.resolve(options.shots), {recursive: true});
const save = () => fs.writeFile(options.output, JSON.stringify(report, null, 2) + '\n', {mode: 0o600});
let browser, page;
async function check(name, operation) {
  try { const detail = await operation(); report.checks.push({name, passed: true, detail}); }
  catch (error) { report.checks.push({name, passed: false, error: error.stack || String(error)}); }
  await save();
}
async function ready() {
  await page.waitForFunction(ids => ids.every(id => {
    const gd = document.getElementById(id);
    return gd?._fullLayout && !gd.classList.contains('recalculating') &&
      !gd.classList.contains('shiny-output-error') &&
      gd.data?.some(t => t.meta?.trajectory_role === 'islet' && t.key?.length) &&
      gd.querySelector('.trajectory-link-highlight');
  }) && !document.documentElement.classList.contains('shiny-busy'), {timeout: 180000}, ids);
}
async function snapshot() {
  return page.evaluate(ids => ids.map(id => {
    const gd = document.getElementById(id), l = gd._fullLayout;
    const rect = node => { if (!node) return null; const r = node.getBoundingClientRect();
      return {left:r.left, right:r.right, top:r.top, bottom:r.bottom, width:r.width, height:r.height}; };
    const widget = rect(gd.querySelector('.svg-container'));
    const guides = [...gd.querySelectorAll('g.legend, g.colorbar')].map(g => ({
      kind: g.classList.contains('legend') ? 'legend' : 'colorbar', bounds: rect(g),
      title: rect(g.querySelector('g.cbtitle > text')), fill: rect(g.querySelector('.cbfill')),
    })).filter(g => g.bounds.width && g.bounds.height);
    const traces = gd.data.filter(t => t.meta?.trajectory_role === 'islet');
    const points = traces.flatMap(t => t.key.map((key, i) => ({key, x:t.x[i], y:t.y[i]})))
      .sort((a,b) => a.key.localeCompare(b.key));
    const pointElements = [...gd.querySelectorAll('.scatterlayer .point')].map(rect)
      .filter(r => r.width && r.height);
    const map = {left:widget.left + l.xaxis._offset, top:widget.top + l.yaxis._offset,
      width:l.xaxis._length, height:l.yaxis._length};
    map.right = map.left + map.width; map.bottom = map.top + map.height;
    return {id, widget, map, guides, x_title:rect(gd.querySelector('.g-xtitle')),
      y_title:rect(gd.querySelector('.g-ytitle')), plot_title:rect(gd.querySelector('.g-gtitle')),
      margin:l.margin, ranges:[l.xaxis.range,l.yaxis.range], points,
      point_bottom:Math.max(...pointElements.map(r => r.bottom)),
      colors:traces.map(t => ({color:t.marker.color, colorscale:t.marker.colorscale || null})),
      legend_labels:[...gd.querySelectorAll('g.legend text.legendtext')].map(e => e.textContent),
      numeric_colorbars:traces.filter(t => t.marker?.showscale).length,
    };
  }), ids);
}
async function settled() {
  await ready();
  let prev = null, same = 0, current;
  for (let attempt = 0; attempt < 35; attempt++) {
    current = await snapshot();
    const signature = JSON.stringify(current.map(s => [s.map, s.margin, s.guides, s.x_title]));
    same = signature === prev ? same + 1 : 0;
    if (same >= 2) return current;
    prev = signature; await pause(150);
  }
  throw new Error('Paired plot geometry did not settle');
}
async function choose(side, value) {
  const id = controls[side];
  const prior = await page.$eval('#' + id, e => e.value);
  if (prior === value) return {side, value, changed: false};
  const before = await snapshot();
  await page.$eval('#' + id + '-selectized', e => e.scrollIntoView({block:'center'}));
  await page.click('#' + id + '-selectized');
  await page.keyboard.down('Control');
  await page.keyboard.press('KeyA');
  await page.keyboard.up('Control');
  await page.keyboard.type(labels[modes.indexOf(value)]);
  const option = `#${id} ~ .selectize-control .selectize-dropdown [data-value="${value}"]`;
  await page.waitForSelector(option, {visible:true});
  // Click a visible option, using normal mouse input rather than selectize.setValue.
  const targets = await page.$$(option);
  for (const target of targets) {
    if (await target.evaluate(e => e.getClientRects().length > 0 && getComputedStyle(e).visibility !== 'hidden')) {
      await target.click(); break;
    }
  }
  await page.waitForFunction((id, value) => document.getElementById(id)?.value === value &&
    Shiny.shinyapp.$inputValues[id] === value, {}, id, value);
  await page.waitForFunction((id, colors) => {
    const gd = document.getElementById(id);
    return gd && !gd.classList.contains('recalculating') && JSON.stringify(gd.data.filter(t =>
      t.meta?.trajectory_role === 'islet').map(t => ({color:t.marker.color,
      colorscale:t.marker.colorscale || null}))) !== colors;
  }, {timeout:60000}, ids[side], JSON.stringify(before[side].colors));
  const after = await settled();
  assert.notEqual(hash(after[side].colors), hash(before[side].colors), 'Physical color choice must change trace colors');
  assert.equal(hash(after[1-side].colors), hash(before[1-side].colors), 'Other panel coloring must remain unchanged');
  assert.deepEqual(after[side].points, before[side].points, 'Changing color must preserve keyed coordinates');
  return {side, value, changed:true, selected_color_hash:hash(after[side].colors), other_color_hash:hash(after[1-side].colors)};
}
function verifyGeometry(plots) {
  assert(Math.abs(plots[0].map.width - plots[1].map.width) < 1, 'Paired map widths differ');
  assert(Math.abs(plots[0].map.height - plots[1].map.height) < 1, 'Paired map heights differ');
  for (const p of plots) {
    assert(p.map.width > 200 && p.map.height > 200, `${p.id}: map unexpectedly small`);
    assert(p.guides.length > 0, `${p.id}: missing legend/colorbar`);
    for (const guide of p.guides) {
      assert(guide.bounds.top >= p.map.bottom + 2, `${p.id}: guide overlaps map`);
      assert(guide.bounds.top >= p.point_bottom + 1, `${p.id}: guide overlaps points`);
      assert(guide.bounds.top >= p.x_title.bottom + 6, `${p.id}: guide overlaps axis title`);
      assert(guide.bounds.left >= p.widget.left - 1 && guide.bounds.right <= p.widget.right + 1,
        `${p.id}: guide clipped horizontally`);
      assert(guide.bounds.bottom <= p.widget.bottom + 1, `${p.id}: guide clipped below widget`);
      if (guide.title && guide.fill) assert(guide.fill.top - guide.title.bottom >= 5,
        `${p.id}: colorbar title overlaps gradient`);
    }
    assert(p.plot_title.left >= p.widget.left - 1 && p.plot_title.right <= p.widget.right + 1,
      `${p.id}: plot title clipped`);
  }
}
async function shot(name, selector = '.paga-comparison-card') {
  const filename = path.join(path.resolve(options.shots), name + '.png');
  await page.$eval(selector, e => e.scrollIntoView({block:'center'}));
  await pause(150);
  const element = await page.$(selector);
  await element.screenshot({path:filename});
  const bytes = await fs.readFile(filename);
  report.screenshots.push({name, path:filename, sha256:crypto.createHash('sha256').update(bytes).digest('hex')});
  return filename;
}
async function physicalLink(side) {
  await page.$eval('#' + ids[side], e => e.scrollIntoView({block:'center'}));
  const target = await page.$eval('#' + ids[side], gd => {
    const t = gd.data.find(t => t.meta?.trajectory_role === 'islet');
    const r = gd.querySelector('.svg-container').getBoundingClientRect(), l = gd._fullLayout;
    const candidates = t.key.map((key,i) => ({key, x:r.left+l.xaxis._offset+l.xaxis.l2p(t.x[i]),
      y:r.top+l.yaxis._offset+l.yaxis.l2p(t.y[i])})).filter(p =>
      p.x > r.left+l.xaxis._offset+25 && p.x < r.left+l.xaxis._offset+l.xaxis._length-25 &&
      p.y > r.top+l.yaxis._offset+25 && p.y < r.top+l.yaxis._offset+l.yaxis._length-25);
    return candidates[Math.floor(candidates.length / 2)];
  });
  assert(target, 'Interior PAGA point required for physical hover');
  await page.mouse.move(target.x, target.y, {steps:5});
  await page.waitForFunction(id => document.getElementById(id)?._hoverdata?.some(p =>
    p.data.meta?.trajectory_role === 'islet'), {timeout:10000}, ids[side]);
  const picked = await page.$eval('#' + ids[side], gd => {const p=gd._hoverdata.find(p=>p.data.meta?.trajectory_role==='islet');
    return p.data.key[p.pointNumber];});
  await page.waitForFunction((ids, key) => ids.every(id => document.querySelector('#' + id +
    ' .trajectory-link-highlight')?.getAttribute('data-linked-key') === key), {}, ids, picked);
  const before = await page.evaluate(() => window.pagaInputEvents.filter(e=>e.name==='traj-linked_click').length);
  await page.mouse.click(target.x, target.y);
  await page.waitForFunction((key, before) => {
    const clicks = window.pagaInputEvents.filter(e => e.name === 'traj-linked_click');
    return clicks.length > before && clicks.at(-1).value?.key === key;
  }, {}, picked, before);
  await page.waitForFunction(key => {
    const el = document.querySelector('#traj-segmentation_viewer_panel');
    const image = el?.querySelector('#islet_segmentation_view img');
    return el?.innerText.includes(key.split('|')[1]) && image?.complete && image.naturalWidth > 0;
  }, {timeout:180000}, picked);
  // Moving away removes hover, so the remaining highlight proves the server pin.
  await page.mouse.move(20, 20);
  await page.waitForFunction((ids, key) => ids.every(id => document.querySelector('#' + id +
    ' .trajectory-link-highlight')?.getAttribute('data-linked-key') === key), {}, ids, picked);
  const clicks = await page.evaluate(() => window.pagaInputEvents.filter(e=>e.name==='traj-linked_click'));
  assert.equal(clicks.length, before + 1, 'One physical click must emit exactly one selection');
  assert.equal(clicks.at(-1).value.source, ids[side]);
  await shot(`physical-selection-${side}`, '#traj-segmentation_viewer_panel');
  return {side, source:ids[side], key:picked, one_click_event:true, segmentation_image_loaded:true};
}

try {
  await save();
  browser = await puppeteer.launch({executablePath:options.chrome || '/usr/bin/google-chrome', headless:true,
    args:['--no-sandbox','--disable-dev-shm-usage','--enable-unsafe-swiftshader']});
  page = await browser.newPage(); page.setDefaultTimeout(180000);
  await page.setViewport({width:1600,height:1100,deviceScaleFactor:1});
  page.on('pageerror', e => report.page_errors.push(String(e)));
  page.on('response', response => {
    if (response.status() >= 400 && !response.url().endsWith('/favicon.ico'))
      report.application_http_errors.push({url:response.url(),status:response.status()});
  });
  await page.goto(url.href,{waitUntil:'domcontentloaded'});
  await page.waitForFunction(() => window.Shiny?.shinyapp?.$socket?.readyState === 1);
  await page.waitForFunction(() => Boolean(window.Plotly?.version));
  report.versions = {chrome:await browser.version(), node:process.version, plotly:await page.evaluate(()=>Plotly.version)};
  await check('pinned_plotly_version', async () => assert.equal(report.versions.plotly, options['expected-plotly'] || '2.11.1'));
  await page.click('a[data-value="Trajectory"]');
  await page.waitForSelector('#traj-paga_show_layout');
  await page.$eval('#traj-paga_show_layout', e => e.scrollIntoView({block:'center'}));
  await page.click('#traj-paga_show_layout');
  await settled();
  await page.evaluate(() => {
    window.pagaInputEvents = [];
    jQuery(document).on('shiny:inputchanged.pagaAcceptance', e => {
      if (!e.isDefaultPrevented() && /^traj-(paga_color|linked_click)/.test(e.name))
        window.pagaInputEvents.push({name:e.name,value:e.value});
    });
  });
  const baseline = (await snapshot())[0].points;
  assert(baseline.length > 0 && new Set(baseline.map(p=>p.key)).size === baseline.length);
  report.cohort = {count:baseline.length,keyed_coordinates_sha256:hash(baseline)};
  for (const width of [1200,1600]) {
    await page.setViewport({width,height:1100,deviceScaleFactor:1});
    await settled();
    for (let i=0;i<modes.length;i++) {
      const left=modes[i],right=modes[(i+1)%modes.length];
      const entry={width,left,right}; report.cases.push(entry);
      await check(`paired_${width}_${left}_${right}`, async () => {
        entry.interactions = [await choose(0,left),await choose(1,right)];
        const plots = await settled();
        entry.plots = plots.map(({points,colors,...p})=>({...p,points:points.length,colors_sha256:hash(colors)}));
        entry.screenshot = await shot(`paga-${width}-${left}-${right}`);
        for (const p of plots) assert.deepEqual(p.points,baseline,'Fitted cohort/coordinates changed');
        assert(plots[0].ranges.flat().every((v,i)=>Math.abs(v-plots[1].ranges.flat()[i])<1e-7),
          'Paired axis ranges differ');
        verifyGeometry(plots);
        return {points:baseline.length, equal_map_area:true, guides_clear:true, physical_selections:entry.interactions};
      });
      await save();
    }
  }
  for (const side of [0,1]) await check(`physical_hover_click_${side}`,()=>physicalLink(side));
  await check('csv_matches_fitted_coordinates_and_both_selections',async()=>{
    await page.waitForFunction(()=>document.querySelector('#traj-dl_paga_csv')?.href.includes('/download/'));
    const csv = await page.$eval('#traj-dl_paga_csv',async link=>{
      const r=await fetch(link.href);return {ok:r.ok,body:await r.text(),disposition:r.headers.get('content-disposition')};
    });
    assert(csv.ok && /attachment/i.test(csv.disposition));
    const current = await page.evaluate(controls=>controls.map(id=>document.getElementById(id).value),controls);
    assert(csv.body.includes(`Color by — left: ${current[0]} | right: ${current[1]}`));
    const lines=csv.body.split(/\r?\n/).filter(line=>line && !line.startsWith('#'));
    const split=line=>line.match(/(?:"(?:[^"]|"")*"|[^,])+/g).map(v=>v.replace(/^"|"$/g,'').replace(/""/g,'"'));
    const headers=split(lines.shift()), rows=lines.map(split).map(values=>Object.fromEntries(headers.map((h,i)=>[h,values[i]])));
    assert.equal(rows.length,baseline.length);
    const byKey=new Map(rows.map(r=>[`${Number(r.case_id)}|${r.islet_id.replace(/^\d+_/,'')}`,r]));
    for(const point of baseline){const row=byKey.get(point.key);assert(row);for(const [field,axis] of [['fa_1','x'],['fa_2','y']])
      assert.equal(Number(Number(row[field]).toPrecision(5)),point[axis]);}
    return {rows:rows.length,bytes:Buffer.byteLength(csv.body),sha256:hash(csv.body),left:current[0],right:current[1]};
  });
  await check('no_visible_errors',async()=>{
    const errors=await page.$$eval('.shiny-output-error',els=>els.filter(e=>e.getClientRects().length).map(e=>e.textContent));
    assert.deepEqual(errors,[]);assert.deepEqual(report.page_errors,[]);assert.deepEqual(report.application_http_errors,[]);
  });
  report.input_events = await page.evaluate(()=>window.pagaInputEvents);
  report.status=report.checks.every(c=>c.passed)?'passed':'failed';
} catch(error) {
  report.status='failed';report.failure=error.stack||String(error);process.exitCode=1;
} finally {
  if(report.status!=='passed'){
    process.exitCode=1;
    if(page)await shot('failure').catch(()=>{});
  }
  report.finished_utc=new Date().toISOString();await save();
  if(browser)await browser.close();
  console.log(JSON.stringify({status:report.status,checks:report.checks.length,failed:report.checks.filter(c=>!c.passed).map(c=>c.name),output:options.output}));
}

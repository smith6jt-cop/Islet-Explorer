#!/usr/bin/env node
// Focused title/metric acceptance using physical controls on a specified preview.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath} from 'node:url';
import puppeteer from 'puppeteer-core';

const options = {};
for (let i=2;i<process.argv.length;i+=2) {
  const key=process.argv[i].replace(/^--/,'');
  assert(['url','output','shots','source','runtime-kind','image','chrome'].includes(key) && process.argv[i+1], 'Invalid argument');
  options[key]=process.argv[i+1];
}
for (const key of ['url','output','shots','source','runtime-kind']) assert(options[key], `--${key} required`);
assert(/^[0-9a-f]{40}$/.test(options.source), 'Full source commit required');
assert(['native-image','source-overlay'].includes(options['runtime-kind']), 'Explicit runtime kind required');
if(options['runtime-kind']==='native-image')assert(/^sha256:[0-9a-f]{64}$/.test(options.image || ''), 'Native image SHA required');
const url=new URL(options.url);
assert(url.protocol==='http:' && ['127.0.0.1','localhost'].includes(url.hostname) &&
  !url.username && !url.password && !url.search && !url.hash,'Loopback HTTP URL required');
const hash=bytes=>crypto.createHash('sha256').update(bytes).digest('hex');
const sourceFile=fileURLToPath(import.meta.url);
const report={schema:'islet.plot_titles_ui/1',status:'running',started_utc:new Date().toISOString(),
  target:{url:url.href,source:options.source,runtime_kind:options['runtime-kind'],image:options.image || null,
    identity_basis:'Caller supplied; pair with the runtime inspection receipt.'},
  harness:{path:sourceFile,sha256:hash(await fs.readFile(sourceFile))},checks:[],screenshots:[],page_errors:[],
  limitations:['Focused title/control behavior and geometry checks; not scientific data or deployment acceptance.',
    'Captured screenshots require separate direct visual review.']};
await fs.mkdir(path.dirname(path.resolve(options.output)),{recursive:true});
await fs.mkdir(path.resolve(options.shots),{recursive:true});
const save=()=>fs.writeFile(options.output,JSON.stringify(report,null,2)+'\n',{mode:0o600});
const pause=ms=>new Promise(resolve=>setTimeout(resolve,ms));
let browser,page;
async function ready() {
  await page.waitForFunction(()=>['plot-plt','plot-dist'].every(id=>{
    const el=document.getElementById(id);return el?._fullLayout && el.data?.length &&
      !el.classList.contains('recalculating') && !el.classList.contains('shiny-output-error') &&
      Math.abs(el._fullLayout.width-el.getBoundingClientRect().width)<3;
  }) && !document.documentElement.classList.contains('shiny-busy'));
}
async function dataPanel(open) {
  await page.waitForSelector('#plot-mode');
  const detail=await page.$('#plot-mode');
  const panel=await detail.evaluateHandle(e=>e.closest('details'));
  if(await panel.evaluate(e=>e.open)!==open) {
    const summary=await panel.asElement().$('summary');await summary.click();
    await page.waitForFunction(open=>document.getElementById('plot-mode').closest('details').open===open,{},open);
  }
  await panel.dispose();await detail.dispose();
  await pause(200);
}
async function choose(id,value) {
  await page.waitForFunction(id=>{const e=document.getElementById(id);return e?.selectize && e.getClientRects().length===0 && e.selectize.$control[0].getClientRects().length>0;},{},id);
  if(await page.$eval('#'+id,e=>e.value)===value)return;
  const label=await page.$eval('#'+id,(e,v)=>e.selectize.options[v]?.[e.selectize.settings.labelField],value);
  assert(label,`Missing choice ${id}/${value}`);
  const input='#'+id+'-selectized';
  await page.$eval(input,e=>e.scrollIntoView({block:'center'}));
  await page.click(input);await page.keyboard.down('Control');await page.keyboard.press('KeyA');await page.keyboard.up('Control');
  await page.keyboard.type(label);
  const option=`#${id} ~ .selectize-control .selectize-dropdown [data-value=${JSON.stringify(value)}]`;
  await page.waitForSelector(option,{visible:true});await page.click(option);
  await page.waitForFunction((id,v)=>document.getElementById(id)?.value===v && Shiny.shinyapp.$inputValues[id]===v,{},id,value);
  await ready();
  if(id==='plot-mode')await page.waitForFunction(v=>['left','right'].every(side=>document.getElementById(`plot-${side}_title`)?.textContent.startsWith(v)),{},value);
}
async function metric(value) {
  const selector=`input[name="plot-metric"][value="${value}"]`;
  await page.waitForSelector(selector,{visible:true});
  await page.$eval(selector,e=>e.scrollIntoView({block:'center'}));await page.click(selector);
  await page.waitForFunction(v=>Shiny.shinyapp.$inputValues['plot-metric']===v,{},value);
  await ready();
}
async function snapshot() {
  return page.evaluate(()=>{
    const rect=e=>{const r=e.getBoundingClientRect();return {left:r.left,top:r.top,right:r.right,bottom:r.bottom,width:r.width,height:r.height};};
    const clean=t=>String(t||'').replace(/<[^>]*>/g,' ').replace(/\s+/g,' ').trim();
    const which=document.getElementById('plot-which');
    return {layer:document.getElementById('plot-mode').value,selection:which.value,
      selection_label:which.selectize.options[which.value]?.[which.selectize.settings.labelField],
      metric:document.querySelector('input[name="plot-metric"]:checked')?.value || null,
      cards:['left','right'].map((side,i)=>{
        const heading=document.getElementById(`plot-${side}_title`),card=heading.closest('.plot-type-card');
        const plot=document.getElementById(i?'plot-dist':'plot-plt');
        return {side,type:document.getElementById(`plot-${side}_type`).value,heading:heading.textContent,
          plot_title:clean(plot._fullLayout.title.text),axis_title:clean(plot._fullLayout.yaxis.title.text),
          heading_box:rect(heading.closest('h5')),card_box:rect(card),
          toolbar_boxes:[...card.querySelectorAll('.plot-type-toolbar .selectize-input,.plot-type-toolbar .btn')].map(rect),
          plot_title_box:rect(plot.querySelector('.g-gtitle')),widget_box:rect(plot.querySelector('.svg-container'))};
      })};
  });
}
async function settled() {
  await ready();let last='',same=0,s;
  for(let i=0;i<30;i++) {
    s=await snapshot();const now=JSON.stringify(s);
    same=now===last?same+1:0;if(same>=2)return s;last=now;await pause(120);
  }
  throw new Error('Plot titles did not settle');
}
function verify(state,expected) {
  assert.equal(state.layer,expected.layer);
  assert.equal(state.selection,expected.selection);
  if(expected.metric)assert.equal(state.metric,expected.metric);
  const typeNames={diameter:'Against Diameter',status:'By Donor Status'};
  const measure=expected.layer==='Microenvironment'?state.selection.split('|')[0]:state.selection_label;
  const suffix={Percentage:'proportion',Counts:'count',Density:'density'}[expected.metric];
  for(const c of state.cards) {
    assert.equal(c.heading,`${state.layer} — ${typeNames[c.type]}`);
    assert(c.plot_title.includes(measure),`${c.side}: selected measure absent from plot title`);
    if(suffix)assert(c.plot_title.includes(suffix),`${c.side}: incorrect metric title`);
    if(expected.metric==='Counts')assert(/Cells per islet/i.test(c.axis_title));
    if(expected.metric==='Density')assert(/Cells per mm/i.test(c.axis_title));
    if(expected.metric==='Percentage')assert(/%.*cells/i.test(c.axis_title));
    if(state.layer==='Microenvironment')assert(c.axis_title.includes(measure) && c.axis_title.includes('area'));
    const h=c.heading_box;
    assert(h.left>=c.card_box.left-1 && h.right<=c.card_box.right+1,`${c.side}: heading outside card`);
    for(const b of c.toolbar_boxes)assert(Math.min(h.right,b.right)-Math.max(h.left,b.left)<=1 ||
      Math.min(h.bottom,b.bottom)-Math.max(h.top,b.top)<=1,`${c.side}: heading overlaps controls`);
  }
  const removeType=t=>t.replace(/ vs Islet Size.*$| across donor groups.*$/,'');
  assert.equal(removeType(state.cards[0].plot_title),removeType(state.cards[1].plot_title),'Cards identify different measures');
}
async function screenshot(name) {
  await dataPanel(false);
  await page.evaluate(()=>{
    const card=document.querySelector('.plot-type-card'),main=document.getElementById('main-content');
    const sticky=document.querySelector('.tab-pane.active .tab-controls');
    main.scrollTop+=card.getBoundingClientRect().top-(sticky.getBoundingClientRect().bottom+12);
  });
  await pause(150);const filename=path.join(path.resolve(options.shots),name+'.png');
  await page.screenshot({path:filename});
  report.screenshots.push({path:filename,sha256:hash(await fs.readFile(filename)),name});
}
async function run(name,expected) {
  await dataPanel(false);
  await page.waitForFunction(expected=>{
    const clean=t=>String(t||'').replace(/<[^>]*>/g,' ').replace(/\s+/g,' ').trim();
    const which=document.getElementById('plot-which');
    if(!which || which.value!==expected.selection)return false;
    const measure=expected.layer==='Microenvironment'?expected.selection.split('|')[0]:expected.label;
    const metric={Percentage:'proportion',Counts:'count',Density:'density'}[expected.metric];
    return ['plot-plt','plot-dist'].every(id=>{const el=document.getElementById(id),title=clean(el?._fullLayout?.title?.text);
      return el && !el.classList.contains('recalculating') && title.includes(measure) && (!metric || title.includes(metric));});
  },{timeout:90000},expected);
  const state=await settled();
  const entry={name,expected,state};report.checks.push(entry);
  try {verify(state,expected);entry.passed=true;}
  catch(e){entry.passed=false;entry.error=e.stack||String(e);}
  await screenshot(name);await save();
}
try {
  await save();browser=await puppeteer.launch({executablePath:options.chrome||'/usr/bin/google-chrome',headless:true,
    args:['--no-sandbox','--disable-dev-shm-usage']});
  page=await browser.newPage();page.setDefaultTimeout(120000);
  page.on('pageerror',e=>report.page_errors.push(String(e)));
  await page.setViewport({width:1600,height:1100});
  await page.goto(url.href,{waitUntil:'domcontentloaded'});
  await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1);
  await page.click('a[data-value="Plot"]');await ready();
  await page.evaluate(()=>{window.plotTitleEvents=[];jQuery(document).on('shiny:inputchanged.titleCheck',e=>{
    if(!e.isDefaultPrevented() && /^plot-(mode|which|metric|left_type|right_type)$/.test(e.name))
      window.plotTitleEvents.push({name:e.name,value:e.value});});});
  report.versions={chrome:await browser.version(),plotly:await page.evaluate(()=>Plotly.version)};
  for(const width of [1200,1600]) {
    await page.setViewport({width,height:1100});
    await dataPanel(true);await choose('plot-mode','Cell Populations');
    const populations=await page.$eval('#plot-which',e=>Object.values(e.selectize.options)
      .filter(o=>/^prop_/.test(o.value)&& /Beta|Alpha/.test(o.value)).map(o=>({value:o.value,label:o[e.selectize.settings.labelField]})));
    assert(populations.length>=2,'Two core endocrine populations required');
    for(const population of populations.slice(0,2))for(const m of ['Percentage','Counts','Density']) {
      await dataPanel(true);await choose('plot-which',population.value);await metric(m);
      await run(`${width}-${population.value.replace(/[^a-zA-Z0-9]+/g,'_')}-${m}`,{layer:'Cell Populations',metric:m,selection:population.value,label:population.label});
    }
    await choose('plot-left_type','status');await choose('plot-right_type','diameter');
    await run(`${width}-cell-populations-swapped`,{layer:'Cell Populations',metric:'Density',selection:populations[1].value,label:populations[1].label});
    await dataPanel(true);await choose('plot-mode','Microenvironment');
    const structures=await page.$eval('#plot-which',e=>Object.values(e.selectize.options)
      .filter(o=>/\|islet_core$/.test(o.value)).slice(0,2).map(o=>o.value));
    assert(structures.length>=2,'Two structural measures required');
    for(const structure of structures) {
      await dataPanel(true);await choose('plot-which',structure);
      await run(`${width}-${structure.replace(/[^a-zA-Z0-9]+/g,'_')}`,{layer:'Microenvironment',selection:structure});
    }
    await choose('plot-left_type','diameter');await choose('plot-right_type','status');
    await run(`${width}-microenvironment-restored-types`,{layer:'Microenvironment',selection:structures[1]});
  }
  report.input_events=await page.evaluate(()=>window.plotTitleEvents);
  report.visible_errors=await page.$$eval('.shiny-output-error',els=>els.filter(e=>e.getClientRects().length).map(e=>e.textContent));
  report.status=report.checks.every(c=>c.passed)&&!report.page_errors.length&&!report.visible_errors.length?'passed':'failed';
}catch(error){report.status='failed';report.failure=error.stack||String(error);}
finally{
  report.finished_utc=new Date().toISOString();await save();
  if(browser)await browser.close();if(report.status!=='passed')process.exitCode=1;
  console.log(JSON.stringify({status:report.status,checks:report.checks.length,output:options.output}));
}

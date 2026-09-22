// Candidate-specific acceptance, including explicit scientific empty states.
import fs from 'node:fs/promises';
import assert from 'node:assert/strict';
import {createRequire} from 'node:module';
import puppeteer from 'puppeteer-core';
const require=createRequire(import.meta.url);
const args=process.argv.slice(2);
const value=k=>args.includes(k)?args[args.indexOf(k)+1]:undefined;
const url=value('--url'),output=value('--output');
assert(url && output,'--url and --output are required');
assert(['127.0.0.1','localhost'].includes(new URL(url).hostname),'Loopback preview required');
const browser=await puppeteer.launch({executablePath:value('--chrome')||'/usr/bin/google-chrome',headless:true,
  args:['--no-sandbox','--disable-dev-shm-usage','--enable-unsafe-swiftshader']});
const evidence={url,checks:{},errors:[],started:new Date().toISOString()};
const pause=ms=>new Promise(resolve=>setTimeout(resolve,ms));
try {
  const page=await browser.newPage();await page.setViewport({width:1600,height:1100});
  page.on('pageerror',e=>evidence.errors.push(e.message));
  await page.goto(url,{waitUntil:'networkidle2',timeout:120000});
  await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1,{timeout:60000});
  assert.equal(await page.title(),'Islet Explorer');
  evidence.checks.app_connected=true;
  evidence.app_ready_seconds=(Date.now()-Date.parse(evidence.started))/1000;
  await page.click('a[data-value="Trajectory"]');
  await page.waitForFunction(()=>document.querySelector('#traj-paga_layout')?.data?.some(t=>t.key?.length),{timeout:60000});
  assert.equal(await page.$('#traj-traj_y_mode'),null);
  assert.equal(await page.$('#traj-traj_feature'),null);
  const donors=await page.$eval('#traj-traj_donor',el=>Object.keys(el.selectize.options));
  assert.equal(donors.length,22);assert(donors.includes('78')&&donors.includes('6539'));
  evidence.checks.donors=21;evidence.checks.ratio_only_controls=true;
  const status=await page.$eval('#traj-traj_status',el=>el.textContent);
  assert(status.includes('islet identity review required'));evidence.checks.identity_review_visible=true;
  await page.waitForFunction(()=>document.querySelectorAll('#traj-multi_heatmap_markers input[type="checkbox"]').length===31,{timeout:30000});
  const panelMarkers=await page.$$eval('#traj-multi_heatmap_markers input[type="checkbox"]',els=>els.map(e=>({value:e.value,checked:e.checked})));
  assert(panelMarkers.every(m=>m.checked));
  for(const marker of ['INS','GCG','SST','CHGA'])assert(panelMarkers.some(m=>m.value===marker));
  evidence.checks.all_31_markers_selected=true;
  await page.waitForFunction(()=>document.querySelector('#traj-traj_scatter_panel')?.textContent.includes('Cell-type ratios unavailable'),{timeout:30000});
  assert.equal(await page.$('#traj-traj_scatter'),null);
  evidence.checks.unavailable_ratios_are_information=true;
  const csv=async id=>{
    await page.waitForFunction(id=>document.querySelector(id)?.href.includes('/download/'),{timeout:30000},id);
    return page.$eval(id,async el=>{const r=await fetch(el.href);if(!r.ok)throw new Error('Download HTTP '+r.status);return r.text();});
  };
  const parse=text=>{
    const rows=[];let row=[],field='',quoted=false;
    for(let i=0;i<text.length;i++){
      const c=text[i];
      if(c==='"'){if(quoted&&text[i+1]==='"'){field+='"';i++;}else quoted=!quoted;}
      else if(c===','&&!quoted){row.push(field);field='';}
      else if(c==='\n'&&!quoted){row.push(field.replace(/\r$/,''));if(row.some(Boolean))rows.push(row);row=[];field='';}
      else field+=c;
    }
    if(field||row.length){row.push(field);rows.push(row);}
    const fields=rows.shift()||[];return rows.map(values=>Object.fromEntries(fields.map((name,i)=>[name,values[i]])));
  };
  const before=parse(await csv('#traj-dl_paga_csv'));
  const reference=new Map(before.map(r=>[r.islet_id,[r.pseudotime_raw,r.pseudotime_percentile]]));
  const extent=await page.$eval('#traj-paga_layout',el=>el.layout.xaxis.range);
  await page.$eval('#traj-traj_donor',el=>el.selectize.setValue('78'));
  await page.waitForFunction(()=>document.querySelector('#traj-paga_layout')?.data?.filter(t=>t.key?.length).every(t=>t.key.every(k=>k.startsWith('78|'))),{timeout:60000});
  const selected=parse(await csv('#traj-dl_paga_csv'));
  assert(selected.length>0 && selected.every(r=>r.case_id==='78'));
  for(const row of selected)assert.deepEqual([row.pseudotime_raw,row.pseudotime_percentile],reference.get(row.islet_id));
  assert.deepEqual(await page.$eval('#traj-paga_layout',el=>el.layout.xaxis.range),extent);
  const ratios=parse(await csv('#traj-dl_trajectory_csv'));
  assert(ratios.length>0 && ratios.every(r=>r.case_id==='78'));
  assert(ratios.every(r=>r.denominator!==''||r.raw_ratio===''));
  const heatmap=parse(await csv('#traj-dl_multi_heatmap_csv'));
  assert(heatmap.length>0 && heatmap.every(r=>r.case_id==='78'));
  assert.deepEqual(Object.keys(heatmap[0]).filter(k=>k.startsWith('core__')),panelMarkers.map(m=>'core__'+m.value));
  evidence.checks.all_31_markers_exported=true;
  evidence.checks.donor_filter_all_downloads=true;evidence.checks.fixed_coordinates=true;
  evidence.checks.unsupported_ratios_missing=true;
  // A valid cell-count range can contain no islets for the selected donor.
  const upper=Math.max(...before.map(r=>Number(r.measured_cells)));
  assert(upper>Math.max(...selected.map(r=>Number(r.measured_cells))));
  const setRange=range=>page.evaluate(range=>{
    window.jQuery('#traj-cell_range').data('ionRangeSlider').update({from:range[0],to:range[1]});
    Shiny.setInputValue('traj-cell_range',range,{priority:'event'});
  },range);
  await setRange([upper,upper]);
  await page.waitForFunction(()=>document.querySelector('#traj-traj_scatter_panel')?.textContent.includes('No islets in this selection'),{timeout:30000});
  for (const id of ['#traj-dl_paga_csv','#traj-dl_trajectory_csv','#traj-dl_multi_heatmap_csv']) assert.equal(parse(await csv(id)).length,0);
  assert.deepEqual(await page.$eval('#traj-paga_layout',el=>el.layout.xaxis.range),extent);
  await page.waitForFunction(()=>document.querySelector('#traj-traj_multi_heatmap_panel')?.textContent.includes('No supported islets'),{timeout:30000});
  assert.deepEqual(await page.$$eval('.shiny-output-error',els=>els.filter(el=>el.getClientRects().length>0).map(el=>el.textContent)),[]);
  await setRange([10,upper]);
  await page.waitForFunction(()=>document.querySelector('#traj-paga_layout')?.data?.some(t=>t.key?.length),{timeout:30000});
  await page.click('#traj-traj_markers_clear');
  await page.waitForFunction(()=>document.querySelector('#traj-traj_multi_heatmap_panel')?.textContent.includes('Choose heatmap markers'),{timeout:30000});
  await page.click('#traj-traj_markers_all');
  await page.waitForSelector('#traj-traj_multi_heatmap img',{timeout:30000});
  assert.equal(await page.$$eval('#traj-multi_heatmap_markers input:checked',els=>els.length),31);
  // A marker change affects only heatmap rows and its download.
  await page.click('#traj-multi_heatmap_markers input[value="INS"]');
  await pause(350);
  let subset=parse(await csv('#traj-dl_multi_heatmap_csv'));
  assert(!('core__INS' in subset[0]) && Object.keys(subset[0]).filter(k=>k.startsWith('core__')).length===30);
  await page.click('#traj-traj_markers_all');await pause(350);
  const restored=parse(await csv('#traj-dl_paga_csv'));
  assert.equal(restored.length,selected.length);
  for(const row of restored)assert.deepEqual([row.pseudotime_raw,row.pseudotime_percentile],reference.get(row.islet_id));
  evidence.checks.empty_filter_and_heatmap_recovery=true;
  // Trigger the same Plotly event that real graph point clicks generate.
  const key=await page.$eval('#traj-paga_layout',el=>el.data.find(t=>t.key?.length).key[0]);
  await page.$eval('#traj-paga_layout',(el,key)=>el.emit('plotly_click',{points:[{data:el.data.find(t=>t.key?.includes(key)),pointNumber:0,curveNumber:el.data.findIndex(t=>t.key?.includes(key)),key}]}),key);
  await page.waitForSelector('#traj-open_image',{timeout:60000});
  const cells=parse(await csv('#traj-dl_islet_cells'));
  assert(cells.length>=10 && 'INS' in cells[0] && 'CD31' in cells[0]);
  evidence.checks.cell_measurements_download=true;
  await page.click('#traj-open_image');
  await page.waitForSelector('#viewer-avivator_viewer',{timeout:60000});
  await page.waitForFunction(()=>document.querySelector('#viewer-avivator_viewer')?.src.includes('0078'),{timeout:60000});
  const images=await page.$eval('#viewer-selected_image',el=>Object.keys(el.selectize.options));
  assert.equal(images.length,21);
  evidence.checks.viewer_donors=21;evidence.checks.click_through_hpap=true;
  const frame=await (await page.$('#viewer-avivator_viewer')).contentFrame();
  await frame.waitForSelector('canvas',{timeout:60000});evidence.checks.image_canvas=true;
  const range=await page.evaluate(async url=>{
    const r=await fetch(url+'local_images/ND_6539.ome.tiff',{headers:{Range:'bytes=0-15'}});
    if(r.status!==206){await r.body?.cancel();return {status:r.status};}
    return {status:r.status,bytes:(await r.arrayBuffer()).byteLength};
  },url);
  assert.deepEqual(range,{status:206,bytes:16});evidence.checks.corrected_image_range=true;
  await page.click('a[data-value="Trajectory"]');
  await page.$eval('#traj-traj_donor',el=>el.selectize.setValue('6539'));
  await page.waitForFunction(()=>!document.querySelector('#traj-open_image'),{timeout:30000});
  evidence.checks.cell_drilldown_respects_donor=true;
  assert.equal(await page.$('a[data-value="Research results"]'),null);
  assert.equal(await page.$('#research-download_results'),null);
  const researchStatus=await page.evaluate(async url=>{const r=await fetch(url+'research_results/H1.png');await r.body?.cancel();return r.status;},url);
  assert.equal(researchStatus,404);
  evidence.checks.research_ui_downloads_and_route_absent=true;
  const unhandled=await page.$$eval('.shiny-output-error',els=>els.filter(e=>e.getClientRects().length>0).map(e=>e.textContent));
  assert.deepEqual(unhandled,[]);assert.deepEqual(evidence.errors,[]);
  await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});
  evidence.accessibility=[];
  for (const tab of ['Plot','Trajectory','Viewer','Spatial']) {
    await page.click(`a[data-value="${tab}"]`);await pause(300);
    try {
      await page.waitForFunction(()=>!Array.from(document.querySelectorAll('.recalculating')).some(el=>el.getClientRects().length>0),{timeout:30000});
    } catch(error) {
      evidence.pending_outputs=await page.$$eval('.recalculating',els=>els.filter(el=>el.getClientRects().length>0).map(el=>({id:el.id,classes:el.className})));
      throw error;
    }
    const violations=await page.evaluate(async()=>{
      const r=await axe.run(document,{runOnly:{type:'tag',values:['wcag2a','wcag2aa']}});
      return r.violations.map(v=>({id:v.id,impact:v.impact,nodes:v.nodes.map(n=>({target:n.target,html:n.html,summary:n.failureSummary}))}));
    });
    evidence.accessibility.push({tab,violations});
  }
  assert(!evidence.accessibility.some(scope=>scope.violations.some(v=>['serious','critical'].includes(v.impact))),'Serious or critical accessibility violation');
  evidence.node=process.version;
  evidence.passed=true;evidence.finished=new Date().toISOString();
} catch(error) { evidence.passed=false;evidence.failure=String(error);process.exitCode=1; }
finally {evidence.browser=await browser.version();await browser.close();await fs.writeFile(output,JSON.stringify(evidence,null,2)+'\n',{mode:0o600});}

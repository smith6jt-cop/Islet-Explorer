// Acceptance for the corrected-boundary candidate; no stale trajectory is allowed.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {createRequire} from 'node:module';
import puppeteer from 'puppeteer-core';
const require=createRequire(import.meta.url), args=process.argv.slice(2);
const value=k=>args.includes(k)?args[args.indexOf(k)+1]:undefined;
let url=value('--url');
const output=value('--output'),root=value('--data-root');
assert(url&&output&&root,'--url, --output, --data-root required');
const targetUrl=new URL(url),allowedHost=value('--allowed-host');
assert(!args.includes('--allowed-host')||(allowedHost&&!allowedHost.startsWith('--')),'--allowed-host requires a hostname');
assert(['127.0.0.1','localhost'].includes(targetUrl.hostname)||targetUrl.hostname===allowedHost,
  'URL must use loopback or the exact --allowed-host hostname');
assert(!targetUrl.username&&!targetUrl.password&&!targetUrl.href.includes('?')&&!targetUrl.href.includes('#'),
  'URL must not include credentials, query parameters, or a fragment');
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
const browser=await puppeteer.launch({executablePath:value('--chrome')||'/usr/bin/google-chrome',headless:true,
  args:['--no-sandbox','--disable-dev-shm-usage','--enable-unsafe-swiftshader']});
const evidence={url,checks:{},errors:[],timings_seconds:{},started:new Date().toISOString()};
const pause=ms=>new Promise(resolve=>setTimeout(resolve,ms));
let page;
try{
  page=await browser.newPage();await page.setViewport({width:1600,height:1100});
  page.on('pageerror',e=>evidence.errors.push(e.message));
  await page.goto(url,{waitUntil:'networkidle2',timeout:120000});
  if (page.url() !== url) {
    assert.equal(new URL(page.url()).hostname,new URL(url).hostname,'Unexpected redirect host');
    evidence.resolved_url=page.url();url=page.url();
  }
  await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1,{timeout:60000});
  assert.equal(await page.title(),'Islet Explorer');
  evidence.timings_seconds.app_ready=(Date.now()-Date.parse(evidence.started))/1000;
  const waitIdle=()=>page.waitForFunction(()=>!Array.from(document.querySelectorAll('.recalculating')).some(el=>el.getClientRects().length>0),{timeout:60000});
  const csv=async id=>{
    await page.waitForFunction(id=>document.querySelector(id)?.href.includes('/download/'),{timeout:30000},id);
    return page.$eval(id,async el=>{const r=await fetch(el.href);if(!r.ok)throw new Error('Download HTTP '+r.status);return r.text();});
  };
  let t=Date.now();await page.click('a[data-value="Trajectory"]');
  await page.waitForFunction(()=>document.querySelector('#traj-traj_status')?.textContent.includes('awaiting a reviewed refit'),{timeout:60000});
  await page.waitForSelector('#traj-traj_multi_heatmap img',{timeout:60000});await waitIdle();
  evidence.timings_seconds.first_trajectory_heatmap=(Date.now()-t)/1000;
  assert.equal(await page.$('#traj-traj_y_mode'),null);
  assert.equal(await page.$('#traj-traj_feature'),null);
  assert.equal(await page.$('#traj-paga_layout'),null);
  assert.equal(await page.$('#traj-traj_scatter'),null);
  assert((await page.$eval('#traj-heatmap_axis_note',el=>el.textContent)).includes('section area'));
  const donors=await page.$eval('#traj-traj_donor',el=>Object.keys(el.selectize.options));
  assert.equal(donors.length,22);assert(donors.includes('78'));
  const markers=await page.$$eval('#traj-multi_heatmap_markers input[type="checkbox"]',els=>els.map(e=>({value:e.value,checked:e.checked})));
  assert.equal(markers.length,31);assert(markers.every(m=>m.checked));
  for(const marker of ['INS','GCG','SST','CHGA'])assert(markers.some(m=>m.value===marker));
  const all=parse(await csv('#traj-dl_paga_csv'));assert.equal(all.length,8568);
  assert(all.every(r=>r.pseudotime_raw===''&&r.pseudotime_percentile===''&&r.supported==='FALSE'));
  evidence.checks.all_21_donors_and_31_markers=true;evidence.checks.stale_trajectory_withdrawn=true;
  t=Date.now();await page.$eval('#traj-traj_donor',el=>el.selectize.setValue('78'));
  await page.waitForFunction(()=>document.querySelector('#traj-boundary_review')?.textContent.includes('1,431'),{timeout:30000});await waitIdle();
  evidence.timings_seconds.hpap_filter_heatmap=(Date.now()-t)/1000;
  const selected=parse(await csv('#traj-dl_paga_csv'));
  assert.equal(selected.length,1431);assert(selected.every(r=>r.case_id==='78'&&r.islet_key.startsWith('AdaptiveIslet_78_')));
  const ratios=parse(await csv('#traj-dl_trajectory_csv'));
  assert.equal(ratios.length,1431);assert(ratios.every(r=>r.case_id==='78'&&r.denominator===''&&r.raw_ratio===''&&r.plotted_value===''));
  const heatmap=parse(await csv('#traj-dl_multi_heatmap_csv'));
  assert.equal(heatmap.length,1431);assert(heatmap.every(r=>r.case_id==='78'));
  assert.deepEqual(Object.keys(heatmap[0]).filter(k=>k.startsWith('core__')),markers.map(m=>'core__'+m.value));
  evidence.checks.corrected_hpap_islets=1431;evidence.checks.donor_filters_all_downloads=true;evidence.checks.unvalidated_ratios_missing=true;
  const inspect=async row=>{
    const key=row.case_id+'|'+row.islet_key;
    const search='#traj-inspect_islet_id + .selectize-control .selectize-input input';
    await page.click(search);await page.type(search,row.islet_key,{delay:5});
    await page.waitForSelector(`.selectize-dropdown [data-value="${key}"]`,{visible:true,timeout:30000});
    await page.click(`.selectize-dropdown [data-value="${key}"]`);
    await page.click('#traj-inspect_islet');await page.waitForSelector('#traj-dl_islet_cells',{timeout:30000});
    await page.waitForFunction(name=>document.querySelector('#traj-segmentation_viewer_panel')?.textContent.includes(name),{timeout:30000},row.islet_key);
    await waitIdle();
  };
  const chosen=selected.find(r=>Number(r.measured_cells)>=20);
  t=Date.now();await inspect(chosen);await page.waitForSelector('#islet_segmentation_view img',{timeout:30000});
  evidence.timings_seconds.islet_outline_and_cells=(Date.now()-t)/1000;
  const cells=parse(await csv('#traj-dl_islet_cells'));
  const sourceCells=parse(await fs.readFile(path.join(root,'cells',chosen.islet_id+'.csv'),'utf8'));
  assert.equal(cells.length,sourceCells.length);
  assert.deepEqual(cells.map(r=>r.object_id+'|'+r.cell_region).sort(),sourceCells.map(r=>r.object_id+'|'+r.cell_region).sort());
  assert.equal(cells.filter(r=>r.cell_region==='core'&&r.marker_measurements_complete==='TRUE').length,Number(chosen.measured_cells));
  t=Date.now();const outline=JSON.parse(await csv('#traj-dl_islet_boundary'));
  evidence.timings_seconds.single_islet_download=(Date.now()-t)/1000;
  assert.deepEqual(outline,JSON.parse(await fs.readFile(path.join(root,'islet_boundaries',chosen.islet_id+'.geojson'),'utf8')));
  assert.equal(outline.features.length,1);assert.equal(outline.features[0].properties.name,chosen.islet_key);
  evidence.checks.exact_outline_and_cell_downloads=true;evidence.inspected_islet=chosen.islet_id;
  await page.screenshot({path:output.replace(/\.json$/,'.png'),fullPage:true});
  await page.click('#traj-open_image');await page.waitForSelector('#viewer-avivator_viewer',{timeout:60000});
  await page.waitForFunction(()=>document.querySelector('#viewer-avivator_viewer')?.src.includes('0078'),{timeout:60000});
  const images=await page.$eval('#viewer-selected_image',el=>Object.keys(el.selectize.options));assert.equal(images.length,21);
  const frame=await (await page.$('#viewer-avivator_viewer')).contentFrame();await frame.waitForSelector('canvas',{timeout:60000});
  assert.equal(await page.$('#viewer-boundary_review'),null);
  evidence.checks.image_click_through_and_canvas=true;
  await page.$eval('#viewer-selected_image',el=>el.selectize.setValue('ND_6539.ome.tiff'));
  await page.waitForFunction(()=>document.querySelector('#viewer-avivator_viewer')?.src.includes('6539'),{timeout:30000});
  assert.equal(await page.$('#viewer-boundary_review'),null);
  evidence.checks.viewer_has_no_boundary_warning=true;
  evidence.checks.viewer_can_change_donor_after_click_through=true;
  const range=await page.evaluate(async url=>{const r=await fetch(url+'local_images/ND_6539.ome.tiff',{headers:{Range:'bytes=0-15'}});if(r.status!==206){await r.body?.cancel();return {status:r.status};}return {status:r.status,bytes:(await r.arrayBuffer()).byteLength};},url);
  assert.deepEqual(range,{status:206,bytes:16});evidence.checks.image_range_loading=true;
  await page.click('a[data-value="Trajectory"]');
  await page.$eval('#traj-traj_donor',el=>el.selectize.setValue('6539'));
  await page.waitForFunction(()=>document.querySelector('#traj-boundary_review')?.textContent.includes('227'),{timeout:30000});await waitIdle();
  assert.equal(await page.$('#traj-open_image'),null);
  for(const id of ['#traj-dl_paga_csv','#traj-dl_trajectory_csv','#traj-dl_multi_heatmap_csv']){
    const rows=parse(await csv(id));assert(rows.length>0&&rows.every(r=>r.case_id==='6539'));
  }
  evidence.checks.pending_donors_labeled_and_filtered=true;
  const upper=Math.max(...all.map(r=>Number(r.measured_cells)));
  await page.evaluate(upper=>{window.jQuery('#traj-cell_range').data('ionRangeSlider').update({from:upper,to:upper});Shiny.setInputValue('traj-cell_range',[upper,upper],{priority:'event'});},upper);
  await page.waitForFunction(()=>document.querySelector('#traj-traj_multi_heatmap_panel')?.textContent.includes('No supported islets'),{timeout:30000});
  for(const id of ['#traj-dl_paga_csv','#traj-dl_trajectory_csv','#traj-dl_multi_heatmap_csv'])assert.equal(parse(await csv(id)).length,0);
  await page.evaluate(upper=>{window.jQuery('#traj-cell_range').data('ionRangeSlider').update({from:0,to:upper});Shiny.setInputValue('traj-cell_range',[0,upper],{priority:'event'});},upper);
  await page.waitForSelector('#traj-traj_multi_heatmap img',{timeout:30000});
  await page.click('#traj-traj_markers_clear');
  await page.waitForFunction(()=>document.querySelector('#traj-traj_multi_heatmap_panel')?.textContent.includes('Choose heatmap markers'),{timeout:30000});
  await page.click('#traj-traj_markers_all');await page.waitForSelector('#traj-traj_multi_heatmap img',{timeout:30000});
  evidence.checks.empty_filter_and_marker_recovery=true;
  assert.equal(await page.$('a[data-value="Research results"]'),null);
  assert.equal(await page.evaluate(async url=>{const r=await fetch(url+'research_results/H1.png');await r.body?.cancel();return r.status;},url),404);
  evidence.checks.research_results_absent=true;
  await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});evidence.accessibility=[];
  for(const tab of ['Plot','Trajectory','Viewer','Spatial']){
    await page.click(`a[data-value="${tab}"]`);await pause(300);await waitIdle();
    const errors=await page.$$eval('.shiny-output-error',els=>els.filter(el=>el.getClientRects().length>0).map(el=>el.textContent));assert.deepEqual(errors,[]);
    const violations=await page.evaluate(async()=>{const r=await axe.run(document,{runOnly:{type:'tag',values:['wcag2a','wcag2aa']}});return r.violations.map(v=>({id:v.id,impact:v.impact,nodes:v.nodes.map(n=>n.target)}));});
    evidence.accessibility.push({tab,violations});
  }
  assert(!evidence.accessibility.some(s=>s.violations.some(v=>['serious','critical'].includes(v.impact))));
  assert.deepEqual(evidence.errors,[]);evidence.passed=true;
}catch(error){evidence.passed=false;evidence.failure=String(error);process.exitCode=1;
  if(page){evidence.visible_errors=await page.$$eval('.shiny-output-error',els=>els.map(e=>({id:e.id,text:e.textContent})));await page.screenshot({path:output.replace(/\.json$/,'.failed.png'),fullPage:true});}}
finally{evidence.node=process.version;evidence.browser=await browser.version();evidence.finished=new Date().toISOString();await browser.close();await fs.writeFile(output,JSON.stringify(evidence,null,2)+'\n',{mode:0o600});console.log(JSON.stringify(evidence,null,2));}

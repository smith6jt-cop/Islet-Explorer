// Run against an isolated schema-4 preview; synthetic evidence is identified.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import puppeteer from 'puppeteer-core';
const args=process.argv.slice(2), value=key=>args.includes(key)?args[args.indexOf(key)+1]:undefined;
const url=value('--url'),root=value('--data-root'),output=value('--output'),synthetic=args.includes('--synthetic');
assert(url&&root&&output,'--url, --data-root and --output are required');
const target=new URL(url);
assert(['127.0.0.1','localhost'].includes(target.hostname)||target.hostname===value('--allowed-host'),'Only loopback or explicit --allowed-host is permitted');
assert(!target.username&&!target.password&&!target.search&&!target.hash,'No URL credentials, query or fragment');
const parse=text=>{
  const rows=[];let row=[],field='',quoted=false;
  for(let i=0;i<text.length;i++){const c=text[i];if(c==='"'){if(quoted&&text[i+1]==='"'){field+='"';i++;}else quoted=!quoted;}
    else if(c===','&&!quoted){row.push(field);field='';}else if(c==='\n'&&!quoted){row.push(field.replace(/\r$/,''));if(row.some(Boolean))rows.push(row);row=[];field='';}else field+=c;}
  if(field||row.length){row.push(field);rows.push(row);}const fields=rows.shift()||[];
  return rows.map(values=>Object.fromEntries(fields.map((name,i)=>[name,values[i]])));
};
const read=async name=>parse(await fs.readFile(path.join(root,name),'utf8'));
const identity=JSON.parse(await fs.readFile(path.join(root,'candidate.json'),'utf8'));
assert.equal(String(identity.schema_version),'4');
assert.equal(Boolean(identity.synthetic),synthetic,'Synthetic fixture must be explicit in both CLI and candidate identity');
const catalog=JSON.parse(await fs.readFile(path.join(root,'typing/catalog.json'),'utf8'));
const names=entries=>entries.map(entry=>typeof entry==='string'?entry:entry.name);
const coordinates=await read('trajectory/coordinates.csv'),features=synthetic?coordinates:await read('trajectory/composition_features.csv'),counts=await read('trajectory/counts.csv'),denominators=await read('trajectory/endocrine_counts.csv');
const markerPanel=JSON.parse(await fs.readFile(path.join(root,'trajectory/marker_panel.json'),'utf8')).markers;
const embeddingManifest=synthetic?{}:JSON.parse(await fs.readFile(path.join(root,'analysis/islet_embedding_manifest.json'),'utf8'));
const classifierEligible=features.filter(row=>row.eligible.toLowerCase()==='true').length;
const embeddedIslets=Number(embeddingManifest.embedded_islets??embeddingManifest.eligible_islets??0);
const supportedPrimary=coordinates.filter(row=>row.pseudotime_raw!=='').length;
const supportedIds=new Set(coordinates.filter(row=>row.pseudotime_raw!=='').map(row=>row.islet_id));
const defaultBetaValues=counts.filter(row=>supportedIds.has(row.islet_id)&&row.compartment==='core'&&row.phenotype==='Beta cell'&&row.plotted_value!=='').map(row=>Number(row.plotted_value));
const savedRanks=coordinates.filter(row=>row.pseudotime_raw!=='').map(row=>Number(row.pseudotime_percentile)).sort((a,b)=>a-b);
assert.equal(savedRanks.length,supportedPrimary);assert(savedRanks.every(Number.isFinite));
for(let index=0;index<savedRanks.length;index++)assert(Math.abs(savedRanks[index]-(index+.5)/savedRanks.length)<1e-12,'Saved reference ranks must be the fixed evenly spaced midpoint grid');
const displayNumber=value=>Number(value).toLocaleString('en-US');
const browser=await puppeteer.launch({executablePath:value('--chrome')||'/usr/bin/google-chrome',headless:true,args:['--no-sandbox','--disable-dev-shm-usage']});
const evidence={url,synthetic,checks:{},errors:[],started:new Date().toISOString()};let page;
try{
  page=await browser.newPage();await page.setViewport({width:1600,height:1100});
  page.on('pageerror',error=>evidence.errors.push(error.message));
  await page.goto(url,{waitUntil:'networkidle2',timeout:120000});
  await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1,{timeout:60000});
  if(!synthetic){
    await page.click('a[data-value="Spatial"]');
    await page.waitForFunction(expected=>document.querySelector('#spatial-leiden_summary')?.textContent.includes(expected),{timeout:60000},displayNumber(embeddedIslets));
    await page.waitForSelector('#spatial-leiden_umap.js-plotly-plot',{timeout:60000});
    await page.waitForSelector('#spatial-spatial_umap_donor.js-plotly-plot',{timeout:60000});
    await page.waitForSelector('#spatial-cluster_composition.js-plotly-plot',{timeout:60000});
    const spatial=await page.evaluate(()=>{
      const umap=document.querySelector('#spatial-leiden_umap');
      const composition=document.querySelector('#spatial-cluster_composition');
      const resolution=document.querySelector('#spatial-leiden_res');
      return {
        points:(umap?.data||[]).reduce((n,trace)=>n+(String(trace.mode||'').includes('markers')?(trace.x?.length||0):0),0),
        resolutions:resolution?.selectize?Object.keys(resolution.selectize.options):[],
        composition:(composition?.data||[]).map(trace=>({name:trace.name,n:trace.x?.length||trace.y?.length||0})),
        summary:document.querySelector('#spatial-leiden_summary')?.textContent||''
      };
    });
    assert.equal(spatial.points,embeddedIslets);assert.equal(spatial.resolutions.length,4);
    assert(spatial.composition.length>0&&spatial.composition.every(trace=>trace.n>0));
    assert(spatial.summary.includes('21 donors'));
    assert(spatial.summary.includes(displayNumber(classifierEligible)));
    assert(spatial.summary.includes('does not contact the image edge'));
    assert(spatial.summary.includes('at least one canonical core cell'));
    evidence.checks.spatial_umap_leiden_and_cluster_composition=true;
  }
  await page.click('a[data-value="Trajectory"]');
  await page.waitForSelector('#traj-traj_multi_heatmap img',{timeout:60000});
  await page.waitForFunction(()=>document.querySelector('#traj-traj_status')?.textContent.includes('does not contact the image edge'),{timeout:60000});
  const trajectoryStatus=await page.$eval('#traj-traj_status',el=>el.textContent||'');
  assert(trajectoryStatus.includes('at least one canonical core cell'));
  assert(trajectoryStatus.includes(displayNumber(classifierEligible)));
  const idle=()=>page.waitForFunction(()=>!Array.from(document.querySelectorAll('.recalculating')).some(el=>el.getClientRects().length),{timeout:60000});
  const choose=async(id,value)=>{await page.$eval(id,(el,value)=>el.selectize?el.selectize.setValue(value):(el.value=value,el.dispatchEvent(new Event('change',{bubbles:true}))),value);await new Promise(resolve=>setTimeout(resolve,100));await idle();};
  const csv=async id=>{await page.waitForFunction(id=>document.querySelector(id)?.href.includes('/download/'),{timeout:30000},id);return parse(await page.$eval(id,async el=>(await fetch(el.href)).text()));};
  const markers=await page.$$eval('#traj-multi_heatmap_markers input[type="checkbox"]',els=>els.map(el=>({value:el.value,checked:el.checked})));
  assert.deepEqual(markers.map(m=>m.value),markerPanel);assert(markers.every(m=>m.checked));assert.equal(markers.length,31);
  const donors=await page.$eval('#traj-traj_donor',el=>Object.keys(el.selectize.options));
  assert.equal(donors.length,22);assert(donors.includes('78'));
  const choices=await page.$eval('#traj-traj_ratio_type',el=>Object.keys(el.selectize.options));
  assert.deepEqual(new Set(choices),new Set([...names(catalog.base_types),...names(catalog.states||[])]));
  const labels=await page.$eval('#traj-traj_ratio_type',el=>Object.values(el.selectize.options).map(option=>option.label));
  assert.deepEqual(new Set(labels),new Set(choices),'Catalog names, not group headings, must label selectable identities');
  assert.equal(await page.$('#traj-traj_y_mode'),null);assert.equal(await page.$('#traj-traj_feature'),null);
  const scatter=await page.$eval('#traj-traj_scatter',plot=>{
    const traces=(plot.data||[]).filter(trace=>String(trace.mode||'').includes('markers'));
    const y=traces.flatMap(trace=>Array.from(trace.y||[]).map(Number).filter(Number.isFinite));
    const x=traces.flatMap(trace=>Array.from(trace.x||[]).map(Number).filter(Number.isFinite)).sort((a,b)=>a-b);
    return {points:x.length,x,minimumY:Math.min(...y),maximumY:Math.max(...y),yRange:plot.layout?.yaxis?.range||[]};
  });
  if(!synthetic)assert.equal(scatter.points,supportedPrimary);else assert(scatter.points>0);
  assert(scatter.minimumY>=0);
  assert.equal(scatter.x.length,savedRanks.length);
  for(let index=0;index<scatter.x.length;index++)assert(Math.abs(scatter.x[index]-savedRanks[index])<1e-12,'Default scatter must use the saved fixed reference ranks');
  if(!synthetic){
    assert(defaultBetaValues.some(value=>value>0),'Canonical supported coordinates must include beta-containing islets');
    assert(scatter.maximumY>0,'The default Beta-cell/core view must display the nonzero canonical beta ratios');
    assert(Math.abs(scatter.maximumY-Math.max(...defaultBetaValues))<1e-10,'Default Beta-cell/core maximum differs from the canonical count table');
    assert(Number(scatter.yRange[1])>=scatter.maximumY,'The y-axis must include the observed nonzero beta ratios');
  }
  for(const id of ['#traj-traj_point_size','#traj-traj_point_size_slider','#traj-traj_alpha','#traj-traj_color_by','#traj-traj_donor_palette'])assert(await page.$(id));
  assert(await page.$('#traj-traj_strip_heatmap.js-plotly-plot'));
  const stripText=await page.$eval('#traj-traj_strip_heatmap',plot=>(plot.data||[]).flatMap(trace=>Array.from(trace.text||[])).join(' '));
  assert(!stripText.includes('Islets: 0'),'Evenly ranked full-cohort points must populate every strip bin');
  evidence.checks.fixed_rank_is_evenly_spaced_and_beta_ratios_are_nonzero=true;
  evidence.checks.dynamic_identities_states_and_all_markers=true;
  const baseline=await csv('#traj-dl_paga_csv');
  const secondary=baseline.filter(row=>row.supported.toLowerCase()==='false');
  if(synthetic){assert(secondary.length>0);assert(secondary.every(row=>row.pseudotime_raw===''&&row.graph_x!==''&&row.graph_y!==''));}
  evidence.checks.disconnected_graph_positions_do_not_create_pseudotime=true;
  const donor=baseline.find(row=>row.case_id==='78')?.case_id||baseline[0].case_id;
  await choose('#traj-traj_donor',donor);
  await page.waitForFunction(donor=>document.querySelector('#traj-traj_y_note')?.textContent.includes('selected islets'),{},donor);await idle();
  const selected=await csv('#traj-dl_paga_csv');assert(selected.length>0&&selected.every(row=>row.case_id===donor));
  for(const row of selected){const original=coordinates.find(r=>r.islet_id===row.islet_id);assert(Math.abs(Number(row.pseudotime_raw)-Number(original.pseudotime_raw))<1e-12);assert(Math.abs(Number(row.pseudotime_percentile)-Number(original.pseudotime_percentile))<1e-12);}
  const cellType=names(catalog.base_types).find(name=>name==='Beta cell')||names(catalog.base_types)[0];
  await choose('#traj-traj_ratio_type',cellType);await choose('#traj-traj_ratio_region','both');
  const ratios=await csv('#traj-dl_trajectory_csv');assert(ratios.every(row=>row.case_id===donor));
  for(const row of ratios){const e=denominators.find(d=>d.islet_id===row.islet_id),a=counts.find(d=>d.islet_id===row.islet_id&&d.phenotype===cellType&&d.compartment==='core'),b=counts.find(d=>d.islet_id===row.islet_id&&d.phenotype===cellType&&d.compartment==='peri');
    if(!e||!a||!b||!e.endocrine_count||Number(e.endocrine_count)===0||a.count===''||b.count===''){assert.equal(row.raw_ratio,'');assert.equal(row.plotted_value,'');}
    else{assert(Math.abs(Number(row.raw_ratio)-(Number(a.count)+Number(b.count))/Number(e.endocrine_count))<1e-10);assert(Math.abs(Number(row.plotted_value)-Math.log1p(Number(row.raw_ratio)))<1e-10);}}
  const expression=await csv('#traj-dl_multi_heatmap_csv');assert(expression.length>0&&expression.every(row=>row.case_id===donor&&row.measurement==='typed_expression'&&row.population==='Endocrine union'));
  assert.deepEqual(new Set(expression.map(row=>row.marker)),new Set(markerPanel));
  for(const row of expression){if(Number(row.n_cells)>0&&row.mean_log1p!=='')assert.equal(Number(row.value),Number(row.mean_log1p));else assert.equal(row.value,'');}
  assert(expression.some(row=>row.value!==''),'Typed expression must contain actual measurements');
  await choose('#traj-heatmap_measure','positive_ratio');
  const positive=await csv('#traj-dl_multi_heatmap_csv');assert(positive.length>0&&positive.every(row=>row.case_id===donor&&row.measurement==='positive_ratio'&&row.population==='All cells'));
  assert.deepEqual(new Set(positive.map(row=>row.marker)),new Set(markerPanel));
  for(const row of positive){if(row.positive_count!==''&&Number(row.denominator)>0)assert(Math.abs(Number(row.value)-Math.log1p(Number(row.positive_count)/Number(row.denominator)))<1e-10);else assert.equal(row.value,'');}
  await choose('#traj-heatmap_measure','typed_expression');
  if((catalog.states||[]).length){await choose('#traj-traj_ratio_type',names(catalog.states)[0]);assert((await csv('#traj-dl_trajectory_csv')).every(row=>row.phenotype===names(catalog.states)[0]));}
  evidence.checks.donor_filters_all_outputs_and_coordinates_stay_fixed=true;evidence.checks.both_heatmap_measurements_and_ratio_arithmetic=true;
  await page.click('#traj-traj_markers_clear');await page.waitForFunction(()=>document.querySelector('#traj-traj_multi_heatmap_panel')?.textContent.includes('Choose heatmap markers'),{timeout:30000});
  assert.equal((await csv('#traj-dl_multi_heatmap_csv')).length,0);
  await page.click('#traj-traj_markers_all');await page.waitForSelector('#traj-traj_multi_heatmap img',{timeout:30000});
  await choose('#traj-traj_donor','all');assert.equal((await csv('#traj-dl_paga_csv')).length,baseline.length);
  if(!synthetic){
    const target=baseline.find(row=>row.islet_id==='112_AdaptiveIslet_112_00007')||baseline[0];
    await page.evaluate(row=>window.Shiny.setInputValue('traj-traj_native_click',
      {key:`${row.case_id}|${row.islet_key}`,nonce:Date.now()},{priority:'event'}),target);
    await page.waitForFunction(()=>document.querySelector('#traj-dl_islet_cells')?.href.includes('/download/'),{timeout:60000});
    const cells=await csv('#traj-dl_islet_cells');assert(cells.length>0);
    for(const column of ['object_id','INS','GCG','SST','CHGA','base_type','classification_status'])assert(column in cells[0]);
    evidence.checks.schema4_islet_clickthrough_and_cell_download=true;
  }
  assert.equal(await page.$('a[data-value="Research results"]'),null);
  const errors=await page.$$eval('.shiny-output-error',els=>els.filter(el=>el.getClientRects().length).map(el=>el.textContent));assert.deepEqual(errors,[]);assert.deepEqual(evidence.errors,[]);
  evidence.checks.empty_marker_recovery_and_research_absent=true;
  await page.screenshot({path:output.replace(/\.json$/,'.png'),fullPage:true});evidence.passed=true;
}catch(error){evidence.passed=false;evidence.failure=String(error);process.exitCode=1;if(page){evidence.visible_errors=await page.$$eval('.shiny-output-error',els=>els.map(el=>({id:el.id,text:el.textContent})));await page.screenshot({path:output.replace(/\.json$/,'.failed.png'),fullPage:true});}}
finally{evidence.browser=await browser.version();evidence.finished=new Date().toISOString();await browser.close();await fs.writeFile(output,JSON.stringify(evidence,null,2)+'\n');console.log(JSON.stringify(evidence,null,2));}

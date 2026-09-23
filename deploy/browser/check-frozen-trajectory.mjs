import fs from 'node:fs/promises';
import puppeteer from 'puppeteer-core';
const value=k=>{const index=process.argv.indexOf(k);return index<0?undefined:process.argv[index+1];};
const mode=value('--mode')||'restored',url=value('--url')||'http://127.0.0.1:7777/',output=value('--output')||'/tmp/islet-trajectory-browser-audit.json';
const report={mode,url,started:new Date().toISOString(),checks:[],errors:[],observations:{}};
function parse(text){let data=text.split(/\r?\n/).filter(x=>!x.startsWith('#')).join('\n');const rows=[];let row=[],field='',quoted=false;for(let i=0;i<data.length;i++){let c=data[i];if(c==='"'){if(quoted&&data[i+1]==='"'){field+='"';i++;}else quoted=!quoted;}else if(c===','&&!quoted){row.push(field);field='';}else if((c==='\n'||c==='\r')&&!quoted){row.push(field);if(row.some(x=>x!==''))rows.push(row);row=[];field='';}else field+=c;}if(field||row.length){row.push(field);rows.push(row);}const head=rows.shift()||[];return rows.map(row=>Object.fromEntries(head.map((k,i)=>[k,row[i]??''])));}
const median=a=>{const s=a.filter(Number.isFinite).sort((a,b)=>a-b);return s.length%2?s[(s.length-1)/2]:(s[s.length/2-1]+s[s.length/2])/2;};
let browser,page;
const save=()=>fs.writeFile(output,JSON.stringify(report,null,2)+'\n',{mode:0o600});
async function check(name,fn){try{const result=await fn();report.checks.push({name,passed:result.passed,...result});}catch(e){report.checks.push({name,passed:false,error:String(e)});}await save();}
const byId=rows=>new Map(rows.map(r=>[r.combined_islet_id||r.islet_id,Number(r.pseudotime)]));
function invariance(a,b){let common=0,max=0;const map=byId(a);for(const r of b){const k=r.combined_islet_id||r.islet_id;if(map.has(k)){common++;max=Math.max(max,Math.abs(map.get(k)-Number(r.pseudotime)));}}return {passed:common>0&&max<1e-12,common,maxAbsoluteDifference:max};}
try{
 browser=await puppeteer.launch({executablePath:'/usr/bin/google-chrome',headless:true,args:['--no-sandbox','--disable-dev-shm-usage','--disable-gpu']});
 page=await browser.newPage();await page.setViewport({width:1680,height:1100});page.setDefaultTimeout(180000);page.on('pageerror',e=>report.errors.push(String(e)));
 await page.goto(url,{waitUntil:'domcontentloaded',timeout:120000});await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1);
 await page.click('a[data-value="Trajectory"]');
 const idle=async()=>{await new Promise(r=>setTimeout(r,700));await page.waitForFunction(()=>!document.documentElement.classList.contains('shiny-busy')&&!Array.from(document.querySelectorAll('.recalculating')).some(e=>e.getClientRects().length),{timeout:180000});};
 await idle();
 const snapshot=()=>page.evaluate(()=>({
 status:document.querySelector('#traj-traj_status')?.textContent||'',note:document.querySelector('#traj-pseudotime_note')?.textContent||'',yNote:document.querySelector('#traj-traj_y_note')?.textContent||'',graphNote:document.querySelector('#traj-traj_graph_note')?.textContent||document.querySelector('#traj-paga_note')?.textContent||'',
 errors:Array.from(document.querySelectorAll('.shiny-output-error')).filter(e=>e.getClientRects().length).map(e=>({id:e.id,text:e.textContent})),
 controls:Array.from(document.querySelectorAll('[id^="traj-"]')).filter(e=>/^(INPUT|SELECT)$/.test(e.tagName)).map(e=>({id:e.id,value:e.value})),
 plots:Array.from(document.querySelectorAll('[id^="traj-"]')).filter(e=>e.classList.contains('shiny-plot-output')||e.classList.contains('plotly')).map(e=>({id:e.id,plotlyTraces:e.data?.length||0,markerPoints:(e.data||[]).filter(t=>String(t.mode||'').includes('markers')).reduce((n,t)=>n+(t.x?.length||0),0),imageLoaded:Boolean(e.querySelector('img')?.naturalWidth)}))}));
 report.observations.initial=await snapshot();await page.screenshot({path:output.replace(/\.json$/,'.initial.png'),fullPage:true});await save();
 if(mode==='restored'){
 const choose=async(id,val)=>{await page.$eval(id,(el,v)=>{if(el.selectize)el.selectize.setValue(v);else{el.value=v;el.dispatchEvent(new Event('change',{bubbles:true}));}},val);await idle();};
 const csv=async(id,label)=>{await page.waitForFunction(id=>document.querySelector(id)?.href.includes('/download/'),{},id);const text=await page.$eval(id,async el=>(await fetch(el.href)).text());await fs.writeFile(output.replace(/\.json$/,'')+'.'+label+'.csv',text,{mode:0o600});return parse(text);};

 // Dynamic selectors can arrive after Shiny first reports idle. Select an
 // actual marker and require scientific outputs before inspecting the plots.
 await page.waitForFunction(()=>document.querySelector('#traj-traj_feature')?.selectize?.options?.INS);
 await choose('#traj-traj_feature','INS');
 await check('all_trajectory_plots_rendered',async()=>{
   const required=['traj-traj_scatter','traj-traj_heatmap','traj-traj_multi_heatmap','traj-traj_umap_donor','traj-traj_umap_feature','traj-paga_layout'];
   await page.waitForFunction(ids=>{
     const scatter=document.querySelector('#traj-traj_scatter');
     const hasPoints=(scatter?.data||[]).some(t=>String(t.mode||'').includes('markers')&&t.x?.length>0);
     return hasPoints&&ids.every(id=>{const el=document.getElementById(id),img=el?.querySelector('img');return el&&!el.classList.contains('recalculating')&&!el.classList.contains('shiny-output-error')&&(el.data?.length>0||(img?.complete&&img.naturalWidth>0));});
   },{},required);
   const s=await snapshot();
   return {passed:s.errors.length===0,plots:s.plots,visibleErrors:s.errors,status:s.status};
 });
 let rank=await csv('#traj-dl_umap_csv','rank-ins');
 await check('rank_scatter_uses_documented_bounded_jitter',async()=>{const scatter=await page.$eval('#traj-traj_scatter',p=>(p.data||[]).filter(t=>String(t.mode||'').includes('markers')).flatMap(t=>Array.from(t.x||[]).map((x,i)=>({x:Number(x),key:t.key?.[i]}))));const map=new Map(rank.map(r=>[r.case_id+'|'+r.islet_key,+r.pseudotime]));const lo=Math.min(...map.values()),hi=Math.max(...map.values()),bound=(hi-lo)*.012;const diffs=scatter.map(r=>Math.abs(r.x-map.get(r.key)));const exactEndpointCount=scatter.filter(r=>r.x===lo||r.x===hi).length;return {passed:scatter.length===rank.length&&diffs.every(v=>Number.isFinite(v)&&v<=bound+1e-12)&&exactEndpointCount===0,scatterPoints:scatter.length,csvRows:rank.length,jitterBound:bound,maxAbsoluteJitter:Math.max(...diffs),exactEndpointCount};});
 await choose('#traj-traj_feature','GCG');let other=await csv('#traj-dl_umap_csv','rank-gcg');await check('rank_axis_feature_switch_invariant',async()=>invariance(rank,other));
 await choose('#traj-pt_scale','raw');await choose('#traj-traj_feature','INS');const raw=await csv('#traj-dl_umap_csv','raw-ins');await choose('#traj-traj_feature','GCG');other=await csv('#traj-dl_umap_csv','raw-gcg');await check('raw_axis_feature_switch_invariant',async()=>invariance(raw,other));
 await choose('#traj-traj_y_mode','ratio');await page.waitForSelector('#traj-traj_ratio_type');await idle();const ratios=await csv('#traj-dl_umap_csv','ratio-raw');await check('ratio_mode_renders_and_keeps_axis',async()=>({...(invariance(raw,ratios)),rows:ratios.length,feature:ratios[0]?.feature_name,finiteValues:ratios.filter(r=>Number.isFinite(+r.feature_value)).length}));
 await choose('#traj-traj_y_mode','marker');await page.waitForSelector('#traj-traj_feature');await choose('#traj-traj_feature','INS');await choose('#traj-pt_scale','rank');
 await page.evaluate(()=>window.Shiny.setInputValue('traj-multi_heatmap_markers',['INS'],{priority:'event'}));await idle();
 const hm=await csv('#traj-dl_multi_heatmap_csv','heatmap-ins'),points=await csv('#traj-dl_umap_csv','heatmap-points');
 await check('heatmap_export_matches_selected_feature_bins',async()=>{const n=hm.length,lo=Math.min(...points.map(r=>+r.pseudotime)),hi=Math.max(...points.map(r=>+r.pseudotime)),bins=Array.from({length:n},()=>[]);for(const r of points){const i=Math.min(n-1,Math.floor((+r.pseudotime-lo)/(hi-lo)*n));bins[i].push(+r.feature_value);}const means=bins.map(b=>b.reduce((s,x)=>s+x,0)/b.length),mu=means.reduce((s,x)=>s+x,0)/n,sd=Math.sqrt(means.reduce((s,x)=>s+(x-mu)**2,0)/(n-1));const deltas=hm.map((r,i)=>Math.abs(+r.bin_mean-means[i]));const zd=hm.map((r,i)=>Math.abs(+r.z_clamped-Math.max(-2.5,Math.min(2.5,(means[i]-mu)/sd))));return {passed:Math.max(...deltas)<1e-9&&Math.max(...zd)<1e-9&&hm.every((r,i)=>+r.n_islets_in_bin===bins[i].length),bins:n,maxMeanDifference:Math.max(...deltas),maxZDifference:Math.max(...zd)};});
 await choose('#traj-pt_scale','raw');
 for(const [lower,upper,n,med] of [[10,1902,2310,.2525],[25,1902,908,.3825],[10,200,2181,.2386]]){
   await page.$eval('#traj-cell_range',(el,v)=>{const slider=window.jQuery(el).data('ionRangeSlider');slider.update({from:v[0],to:v[1]});window.jQuery(el).trigger('change');},[lower,upper]);await idle();
   const rows=await csv('#traj-dl_umap_csv',`range-${lower}-${upper}`);await check(`documented_range_${lower}_${upper}`,async()=>({passed:rows.length===n&&Math.abs(median(rows.map(r=>+r.pseudotime))-med)<=.00005,rows:rows.length,expectedRows:n,medianPseudotime:median(rows.map(r=>+r.pseudotime)),documentedMedian:med,status:(await snapshot()).status,note:(await snapshot()).note}));
 }
 report.observations.final=await snapshot();await page.screenshot({path:output.replace(/\.json$/,'.final.png'),fullPage:true});
 for(const id of ['traj-traj_multi_heatmap','traj-paga_layout']){
   const panel=await page.$('#'+id);
   await panel.evaluate(el=>el.scrollIntoView({block:'center'}));
   await panel.screenshot({path:output.replace(/\.json$/,'.'+id+'.png')});
 }
 }
}catch(e){report.failure=String(e);if(page){report.observations.failureText=await page.evaluate(()=>document.body.innerText);await page.screenshot({path:output.replace(/\.json$/,'.failed.png'),fullPage:true});}}
finally{report.finished=new Date().toISOString();report.passed=!report.failure&&report.errors.length===0&&report.checks.every(c=>c.passed);if(browser)await browser.close();await save();console.log(JSON.stringify(report,null,2));if(!report.passed)process.exitCode=1;}

import fs from 'node:fs/promises';
import assert from 'node:assert/strict';
import puppeteer from 'puppeteer-core';

const args = process.argv.slice(2);
const value = key => args.includes(key) ? args[args.indexOf(key) + 1] : undefined;
const url = value('--url');
const output = value('--output');
assert(url && output, '--url and --output are required');
const target = new URL(url);
assert(['127.0.0.1', 'localhost'].includes(target.hostname), 'Only loopback previews are permitted');

const browser = await puppeteer.launch({
  executablePath: value('--chrome') || '/usr/bin/google-chrome',
  headless: true,
  args: ['--no-sandbox', '--disable-dev-shm-usage']
});
const evidence = {url, started: new Date().toISOString(), errors: []};
let page;
try {
  page = await browser.newPage();
  await page.setViewport({width: 1600, height: 1100});
  page.on('pageerror', error => evidence.errors.push(error.message));
  const started = Date.now();
  await page.goto(url, {waitUntil: 'networkidle2', timeout: 120000});
  await page.waitForFunction(
    () => window.Shiny?.shinyapp?.$socket?.readyState === 1,
    {timeout: 60000}
  );
  await page.waitForFunction(() => {
    const selector = document.querySelector('#plot-which');
    const plot = document.querySelector('#plot-plt');
    return selector?.selectize?.getValue() === 'prop_Beta cell' &&
      Array.isArray(plot?.data) &&
      plot.data.some(trace => Array.from(trace.y || []).some(Number.isFinite));
  }, {timeout: 60000});

  const state = await page.evaluate(() => {
    const selector = document.querySelector('#plot-which');
    const plot = document.querySelector('#plot-plt');
    const values = Object.keys(selector.selectize.options);
    const finiteY = (plot.data || []).reduce(
      (n, trace) => n + Array.from(trace.y || []).filter(Number.isFinite).length,
      0
    );
    return {
      selected: selector.selectize.getValue(),
      values,
      finiteY,
      title: plot.layout?.title?.text || '',
      loadingPlaceholderVisible: Array.from(document.querySelectorAll('#plot-plt text'))
        .some(node => node.textContent.includes('Loading data')),
      visibleErrors: Array.from(document.querySelectorAll('.shiny-output-error'))
        .filter(node => node.getClientRects().length)
        .map(node => ({id: node.id, text: node.textContent}))
    };
  });

  for (const identity of ['Beta cell', 'Alpha cell', 'Delta cell', 'CD4 T cell',
                           'CD8 T cell', 'Macrophage-associated cell']) {
    assert(state.values.includes(`prop_${identity}`), `Missing core identity: ${identity}`);
    assert(state.values.includes(`peri_prop_${identity}`), `Missing peri identity: ${identity}`);
  }
  assert(!state.values.includes('Ins_any'), 'All-missing legacy hormone field is still selectable');
  assert.equal(state.selected, 'prop_Beta cell');
  assert(state.finiteY > 0, 'The default Plot contains no measured values');
  assert.equal(state.loadingPlaceholderVisible, false);
  assert.deepEqual(state.visibleErrors, []);
  assert.deepEqual(evidence.errors, []);

  // An identity retained for display but absent from this candidate must end
  // in a stable empty state, then recover when a measured identity is chosen.
  await page.$eval('#plot-which', element => element.selectize.setValue('prop_Acinar'));
  await page.waitForFunction(() =>
    document.querySelector('#plot-plt')?.textContent
      .includes('No measured islets match this population and the current filters.'),
  {timeout: 30000});
  assert.equal(await page.$('.shiny-output-error'), null);
  await page.$eval('#plot-which', element => element.selectize.setValue('prop_Beta cell'));
  await page.waitForFunction(() => {
    const plot = document.querySelector('#plot-plt');
    return Array.isArray(plot?.data) &&
      plot.data.some(trace => Array.from(trace.y || []).some(Number.isFinite));
  }, {timeout: 30000});
  evidence.unavailableIdentityEmptyStateAndRecovery = true;
  evidence.loadMilliseconds = Date.now() - started;
  evidence.state = state;
  evidence.passed = true;
  await page.screenshot({path: output.replace(/\.json$/, '.png'), fullPage: true});
} catch (error) {
  evidence.passed = false;
  evidence.failure = String(error);
  process.exitCode = 1;
  if (page) await page.screenshot({path: output.replace(/\.json$/, '.failed.png'), fullPage: true});
} finally {
  evidence.browser = await browser.version();
  evidence.finished = new Date().toISOString();
  await browser.close();
  await fs.writeFile(output, JSON.stringify(evidence, null, 2) + '\n');
  console.log(JSON.stringify(evidence, null, 2));
}

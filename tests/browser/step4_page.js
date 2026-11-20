// Drives the real Step 4 page in a browser against tests/browser/serve_step4.groovy. Needs Node with Playwright:
//   python tools/run_tests.py --browser
// Screenshots are written next to the URL file (build/web-step4), to look at after a change to the page.
const { chromium } = require('playwright');
const fs = require('fs');
const assert = require('assert');
(async () => {
  const out = process.argv[2], url = fs.readFileSync(out + '/url.txt', 'utf8').trim();
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const problems = [];
  page.on('console', (m) => { if (m.type() === 'error' && !/400 \(Bad Request\)/.test(m.text())) problems.push('console: ' + m.text()); });
  page.on('pageerror', (e) => problems.push('page error: ' + e.message));
  const text = (selector) => page.textContent(selector);
  await page.goto(url); await page.waitForFunction(() => window.calciumStep4.S.state);
  assert.strictEqual(await text('#recording'), 'Fixture recording');
  assert.strictEqual(await page.inputValue('#sampling_hz'), '2'); assert.strictEqual(await page.inputValue('#amplitude_floor'), '0.1');
  assert.ok(await page.isHidden('#kcl-row')); await page.check('#kcl'); assert.ok(await page.isVisible('#kcl-row')); await page.uncheck('#kcl');
  assert.deepStrictEqual(await page.$$eval('#view option', (o) => o.map((e) => e.textContent)), ['Colour picture', 'Movie: change from rest', 'Movie: as recorded']);
  await page.screenshot({ path: out + '/page-before.png' });

  // the plugin's own checks answer before anything runs
  await page.click('#analyze'); await page.waitForFunction(() => /Confirm/.test(document.getElementById('status').textContent));
  await page.check('#confirmed'); await page.fill('#stimulus', 'Capsaicin'); await page.fill('#stimulus_start', '3'); await page.fill('#stimulus_end', '14'); await page.fill('#response_baseline_seconds', '2.5');
  assert.ok(await page.$eval('#stimulus_start', (e) => e.classList.contains('changed')));
  await page.click('#analyze'); await page.waitForFunction(() => window.calciumStep4.S.results, null, { timeout: 15000 });
  await page.waitForFunction(() => /^Complete/.test(document.getElementById('status').textContent));
  const sent = fs.readFileSync(out + '/settings.txt', 'utf8');
  assert.ok(/sampling_hz=2\.0\n/.test(sent) && /amplitude_floor=0\.10\n/.test(sent) && /stimulus_start=3\n/.test(sent) && /stimulus=Capsaicin\n/.test(sent), 'untouched boxes send exact values, edited ones what was typed');

  // results first: tiles, groups, the selected cell's verdict, the movie beside the traces
  assert.ok((await text('#tile-responding')).includes('67 %') && (await text('#tile-responding')).includes('2 of 3 evaluable'));
  assert.strictEqual((await text('#tile-registered')).replace(/\D/g, ''), '6');
  assert.deepStrictEqual(await page.$$eval('#group option', (o) => o.map((e) => e.textContent)), ['All cells (6)', 'Responded (2)', 'Did not respond (1)', 'Excluded or unscorable (3)']);
  assert.ok((await text('#verdict')).includes('Cell001: Responded') && (await text('#verdict')).includes('largest rise 1.00 at 6 s (needed 0.10)'));
  assert.strictEqual(await page.inputValue('#view'), 'change');
  await page.waitForFunction(() => /minus|−/.test(document.getElementById('picture-note').textContent));
  assert.ok((await text('#legend')).includes('responded (2)') && (await text('#legend')).includes('did not respond (1)'));
  await page.screenshot({ path: out + '/page-review.png' });

  // time: the slider, Play and the picture move together
  await page.$eval('#time', (e) => { e.value = 12; e.dispatchEvent(new Event('input')); });
  assert.ok((await text('#time-label')).includes('5.5 s') && (await text('#time-label')).includes('frame 12 of 30'));
  await page.waitForFunction(() => /^frame 12 /.test(document.getElementById('picture-note').textContent));
  await page.click('#play'); await page.waitForFunction(() => window.calciumStep4.S.frame > 14); await page.click('#play');
  assert.strictEqual(await text('#play'), '▶ Play');

  // groups and cells: the list, a tile, the arrow keys and a click on the picture
  await page.selectOption('#group', '1'); assert.deepStrictEqual(await page.$$eval('#cell option', (o) => o.map((e) => e.value)), ['Cell001', 'Cell005']);
  await page.keyboard.press('ArrowRight'); assert.strictEqual(await page.inputValue('#cell'), 'Cell005');
  await page.click('#tile-excluded'); assert.deepStrictEqual(await page.$$eval('#cell option', (o) => o.map((e) => e.value)), ['Cell003', 'Cell004', 'Cell006']);
  assert.ok((await text('#verdict')).includes('Cell003: Excluded'));
  await page.click('#tile-registered');
  const spot = await page.evaluate(() => { const c = document.getElementById('picture'), b = c.getBoundingClientRect(), i = window.calciumStep4.S.info;
    const z = Math.min((c.clientWidth - 8) / i.width, (c.clientHeight - 8) / i.height); return { x: b.left + (c.clientWidth - i.width * z) / 2 + 70 * z, y: b.top + (c.clientHeight - i.height * z) / 2 + 60 * z }; });
  await page.mouse.click(spot.x, spot.y); assert.strictEqual(await page.inputValue('#cell'), 'Cell002');
  await page.click('#zoom-cell'); assert.ok(await page.evaluate(() => window.calciumStep4.Pic.zoom > 1)); await page.click('#zoom-all');

  // heatmap, tables and the folder button
  await page.click('.tab[data-tab=heatmap]');
  assert.deepStrictEqual(await page.evaluate(() => window.calciumStep4.Heat.names), ['Cell001', 'Cell005', 'Cell002', 'Cell003', 'Cell004', 'Cell006']);
  await page.screenshot({ path: out + '/page-heatmap.png' });
  await page.selectOption('#order', 'measured'); assert.strictEqual(await page.evaluate(() => window.calciumStep4.Heat.names[1]), 'Cell002');
  const heat = await (await page.$('#heatmap')).boundingBox(); await page.mouse.click(heat.x + 200, heat.y + 26 + 15 * 4 + 7);
  assert.strictEqual(await page.inputValue('#cell'), 'Cell005'); assert.strictEqual(await text('.tab.on'), 'Cells and traces');
  await page.click('.tab[data-tab=calls]'); assert.strictEqual(await page.$$eval('#table-calls tbody tr', (r) => r.length), 6);
  await page.click('#table-calls tbody tr:nth-child(2)'); assert.strictEqual(await page.evaluate(() => window.calciumStep4.S.cell), 'Cell002');
  await page.screenshot({ path: out + '/page-calls.png' });
  await page.click('#open-folder'); await page.waitForFunction(() => true); await page.waitForTimeout(300);
  assert.strictEqual(fs.readFileSync(out + '/opened.txt', 'utf8'), 'simpleca');
  await page.click('#toggle-settings'); assert.ok(await page.isHidden('aside'));
  assert.deepStrictEqual(problems, []);
  await browser.close();
  console.log('PASS: Step 4 page in a browser: settings checks, exact untouched values, tiles, groups, verdict, movie linked to time, Play, cell picking, heatmap order, tables and the folder button');
})().catch((e) => { console.error(e); process.exit(1); });

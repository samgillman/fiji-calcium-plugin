// Exercise the real Step 4 page and local server. This fixture has synthetic saved
// results; the adapter checks separately exercise the real R engine.
// Run: python tools/run_tests.py --browser. CI uploads these screenshots.
'use strict';
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');
const assert = require('assert/strict');

(async () => {
  const out = process.argv[2], url = fs.readFileSync(path.join(out, 'url.txt'), 'utf8').trim();
  let browser, page;
  try {
    browser = await chromium.launch();
    page = await browser.newPage({ viewport: { width: 1440, height: 900 }, acceptDownloads: true });
    page.setDefaultTimeout(10000);
    const problems = [];
    let injectedFrameFailure = false, refusals = false, aroundAborted = false;
    page.on('console', (message) => {
      if (message.type() !== 'error') return;
      if (injectedFrameFailure && /503|410/.test(message.text())) return;
      if (refusals && /status of 400/.test(message.text())) return;   // a review mark refused on purpose
      if (aroundAborted && /ERR_FAILED/.test(message.text())) return;   // an answer for "Around the cell" lost on purpose
      problems.push('console: ' + message.text());
    });
    page.on('pageerror', (error) => problems.push('page error: ' + error.message));
    const text = (selector) => page.textContent(selector);
    // The Cell box takes typing: a name and Enter opens that cell, an empty box is no cell.
    const openCell = async (name) => { await page.fill('#cell', name); await page.press('#cell', 'Enter'); };
    const cellsOffered = () => page.locator('#cell-names option').evaluateAll((options) => options.map((o) => o.value));
    const waitFrame = (frame) => page.waitForFunction((wanted) => {
      const app = window.calciumStep4;
      return app.S.frame === wanted && app.Pic.shownFrame === wanted && !app.Pic.loading && document.getElementById('movie-message').hidden;
    }, frame, { polling: 100 });
    const slider = async (frame) => {
      await page.locator('#time').evaluate((element, wanted) => {
        element.value = wanted; element.dispatchEvent(new Event('input', { bubbles: true }));
      }, frame);
      await waitFrame(frame);
    };
    const gray = () => page.evaluate(() => {
      const canvas = document.createElement('canvas'); canvas.width = 250; canvas.height = 150;
      const context = canvas.getContext('2d'); context.drawImage(window.calciumStep4.Pic.shown, 0, 0);
      return context.getImageData(40, 60, 1, 1).data[0];
    });
    const settledPaint = () => page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
    // The figure editor. A number is taken when Enter is pressed or its box is left, never while it is typed: typeBox types one whole
    // and presses Enter. boxSays: what a number box holds, whether it is marked, and what the line under it says.
    const typeBox = async (id, value) => { await page.fill('#figure-' + id, value); await page.press('#figure-' + id, 'Enter'); };
    const boxSays = (id) => page.evaluate((id) => { const box = document.getElementById('figure-' + id), note = document.getElementById('figure-' + id + '-note'); return [box.value, box.getAttribute('aria-invalid') === 'true', note.hidden ? '' : note.textContent]; }, id);
    // A figure is written in a font with a name (Arial), in which the painter draws lowered and raised digits (the ₀ of ΔF/F₀) from
    // the font's ordinary digits, smaller, as runs of their own. The two helpers that read a drawing back (drawing, heatDrawing) note
    // of each text where it was written, how wide and at what size, and asWritten joins such runs into the words they are.
    await page.addInitScript(() => {   // (with the page, also when it is loaded anew)
      window.textSeen = (ctx, text, x, y, at) => { const m = ctx.getTransform(), p = at(ctx, x, y), lw = ctx.measureText(text).width;
        return { text: String(text), x: p[0], y: p[1], w: lw * Math.hypot(m.a, m.b), ink: ctx.fillStyle, lx: x, ly: y, lw, size: parseFloat(/([\d.]+)px/.exec(ctx.font)[1]), align: ctx.textAlign, m: [m.a, m.b, m.c, m.d, m.e, m.f].join() }; };
      window.asWritten = (texts) => { const out = [];
        for (const t of texts) { const last = out[out.length - 1];
          if (last && t.m === last.m && Math.abs(t.lx - last.end) < 0.01 && t.size !== last.sizeNow) {
            const low = t.ly > last.ly, small = t.size < last.size;
            last.text += small ? Array.from(t.text, (c) => (low ? '₀₁₂₃₄₅₆₇₈₉₊₋' : '⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻')['0123456789+−'.indexOf(c)]).join('') : t.text; last.end = t.lx + t.lw; last.sizeNow = t.size; last.w += t.w;
          } else out.push(Object.assign({ end: t.lx + t.lw, sizeNow: t.size }, t)); }
        return out; };
    });
    // A saved picture against release 2026.11.22.2. The release's drawing of a saved picture is kept here word for word, but for
    // what it reads from the page: the picture's scale and the colours of the view. Returns how many colour values of the page's
    // own picture (viewCanvas: a picture, or with "labelled" a frame of a set; of the area Export would save now, which is the
    // whole picture unless "Area" says otherwise) differ from it; none, when the picture is as it was.
    const againstRelease = (frame, labelled) => page.evaluate(async ([frame, labelled]) => {
      const app = window.calciumStep4, S = app.S, Pic = app.Pic, ON_PICTURE = { responded: '#3ddbb5', not: '#ffa24a', other: '#aeb8c2' }, finite = (v) => typeof v === 'number' && isFinite(v);
      const VIEW_NAMES = { picture: 'colour overview', raw: 'as recorded', flat: 'background removed', green: 'GFP', magenta: 'response in magenta', change: 'change from baseline' };
      const calls = S.results ? ((S.results.callsByWindow || {})[S.window] || S.results.calls || {}) : {}, callOf = (name) => (calls[name] || {}).Status, kind = (call) => call === 'Responded' ? 'responded' : call === 'Did not respond' ? 'not' : 'other';
      const inGroup = (call, group) => group === 0 || (group === 1 && call === 'Responded') || (group === 2 && call === 'Did not respond') || (group === 3 && call === 'Excluded') || (group === 4 && call === 'Unscorable');
      const made = await app.viewCanvas(frame, labelled, app.savedArea()), f = frame || S.frame; let image = S.view === 'background' ? Pic.background : Pic.base;
      if (!(S.view === 'picture' || S.view === 'background')) image = await createImageBitmap(await (await fetch('api/frame?f=' + f + '&mode=' + S.view + '&v=' + encodeURIComponent(app.Buf.version || ''), { cache: 'no-store' })).blob());
      const iw = S.info.width, ih = S.info.height, scale = made.width / iw, font = getComputedStyle(document.body).fontFamily;
      const lines = !S.outlines ? [] : S.results ? [[ON_PICTURE.responded, 'responded', true], [ON_PICTURE.not, 'did not respond', true], [ON_PICTURE.other, 'excluded or unscorable', true]] : [['#ffffff', 'your cells', true]], key = labelled ? [] : app.pictureKey().filter((item) => !item[2]).concat(lines);
      const unit = Math.max(1, iw / 1388), size = Math.round(20 * unit), gap = Math.round(30 * unit), measure = document.createElement('canvas').getContext('2d'); measure.font = size + 'px ' + font;
      const places = []; let x = 14 * unit, row = 0; for (const item of key) { const width = size * 1.4 + measure.measureText(item[1]).width + size * 1.4; if (x + width > iw - 8 && x > 14 * unit) { x = 14 * unit; row++; } places.push([x, row]); x += width; }
      const keyHeight = labelled ? 0 : Math.round((row + 1) * gap + gap * 1.5), canvas = document.createElement('canvas'); canvas.width = Math.round(iw * scale); canvas.height = Math.round((ih + keyHeight) * scale);
      const ctx = canvas.getContext('2d'); ctx.scale(scale, scale); ctx.fillStyle = '#000'; ctx.fillRect(0, 0, iw, ih + keyHeight); ctx.drawImage(image, 0, 0, iw, ih);
      const known = !!S.results; ctx.lineJoin = 'round';
      if (S.outlines) S.paths.forEach((path, i) => {
        const call = callOf(S.info.names[i]); if (known && !inGroup(call, S.group, S.info.names[i])) return; const which = kind(call);
        ctx.setLineDash([]); ctx.strokeStyle = 'rgba(0,0,0,0.55)'; ctx.lineWidth = 3.2 * unit * S.outlineWidth; ctx.stroke(path);
        ctx.setLineDash(known && which === 'other' ? [5 * unit, 4 * unit] : []); ctx.strokeStyle = known ? ON_PICTURE[which] : '#ffffff'; ctx.lineWidth = 1.8 * unit * S.outlineWidth; ctx.stroke(path);
      });
      ctx.setLineDash([]); ctx.font = size + 'px ' + font; ctx.textBaseline = 'middle';
      key.forEach((item, i) => { const kx = places[i][0], ky = ih + gap * 0.9 + places[i][1] * gap;
        if (item[2]) { ctx.strokeStyle = item[0]; ctx.lineWidth = 2 * unit; ctx.strokeRect(kx, ky - size * 0.45, size * 0.9, size * 0.9); } else { ctx.fillStyle = item[0]; ctx.fillRect(kx, ky - size * 0.45, size * 0.9, size * 0.9); }
        ctx.fillStyle = '#fff'; ctx.fillText(item[1], kx + size * 1.4, ky); });
      if (labelled) {
        const time = app.frameTimes()[f - 1], text = 'Frame ' + f + (finite(time) ? ' · ' + time.toFixed(1) + ' s' : '') + ' · ' + (VIEW_NAMES[S.view] || S.view);
        ctx.font = '600 ' + Math.round(size * 1.1) + 'px ' + font; const tw = ctx.measureText(text).width;
        ctx.fillStyle = 'rgba(0,0,0,0.7)'; ctx.fillRect(10 * unit, 10 * unit, tw + 20 * unit, size * 1.9); ctx.fillStyle = '#fff'; ctx.fillText(text, 20 * unit, 10 * unit + size * 0.98);
      } else { ctx.fillStyle = '#AEB8C2'; ctx.font = Math.round(size * 0.8) + 'px ' + font; ctx.fillText(S.info.recording, 14 * unit, ih + keyHeight - gap * 0.55); }
      if (!Number.isInteger(scale) || made.width !== canvas.width || made.height !== canvas.height) return { size: [made.width, made.height], different: -1 };
      const a = made.getContext('2d').getImageData(0, 0, made.width, made.height).data, b = ctx.getImageData(0, 0, canvas.width, canvas.height).data; let different = 0, lit = 0;
      for (let i = 0; i < a.length; i++) { if (a[i] !== b[i]) different++; if (i % 4 < 3 && b[i] > 60) lit++; }
      return { size: [made.width, made.height], different, lit };
    }, [frame || null, !!labelled]);

    const screenshots = async (stage) => {
      for (const [name, width, height] of [['desktop', 1440, 900], ['compact', 1024, 768], ['mobile', 390, 768]]) {
        await page.setViewportSize({ width, height });
        await page.evaluate(() => window.scrollTo(0, 0)); await settledPaint();
        assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1), name + ' layout must not overflow horizontally');
        assert.ok(await page.isVisible('#transport'), name + ' keeps playback controls visible');
        if (stage === 'after' && width >= 1024) assert.ok(await page.evaluate(() => {
          const picture = document.getElementById('picture-card').getBoundingClientRect();
          const plot = document.getElementById('plot-card').getBoundingClientRect();
          return Math.abs(picture.top - plot.top) <= 2 && picture.right <= plot.left &&
            picture.height >= 350 && plot.height >= 350;
        }), name + ' keeps movie and traces readable, with page scrolling when needed');
        if (stage === 'after' && width >= 1024) assert.ok(await page.evaluate(() => {
          const note = document.getElementById('picture-note').getBoundingClientRect();
          const card = document.getElementById('picture-card').getBoundingClientRect();
          const fields = document.querySelector('.settings-scroll').getBoundingClientRect();
          const run = document.querySelector('.run-box').getBoundingClientRect();
          return note.bottom <= card.bottom && fields.bottom <= run.top;
        }), name + ' keeps the frame caption visible and run controls outside the settings scroll');
        await page.screenshot({ path: path.join(out, 'page-' + stage + '-' + name + '.png') });
      }
      await page.setViewportSize({ width: 1440, height: 900 }); await settledPaint();
    };
    const importProtocol = async (protocol, filename) => {
      const response = page.waitForResponse((r) => r.url().endsWith('/api/protocol') && r.request().method() === 'POST');
      await page.setInputFiles('#protocol-import', { name: filename, mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(protocol)) });
      assert.ok((await response).ok(), 'the protocol import must pass canonical validation');
      await page.waitForFunction((name) => document.getElementById('protocol-name').value === name && !document.getElementById('protocol-import').value, protocol.name);
    };

    await page.goto(url);
    await page.waitForFunction(() => window.calciumStep4 && window.calciumStep4.S.state);
    await waitFrame(1);
    assert.equal(await text('#recording'), 'Fixture recording');
    // The page is given its draft in the words it writes, so nothing reads as changed before anything is touched.
    assert.ok(await page.isHidden('#discard'), 'an untouched page offers nothing to discard');
    const newEstimator = await page.evaluate(() => window.calciumStep4.S.info.newEstimator);
    assert.equal(newEstimator, 'larger_of_baseline_sd_and_trace_noise');
    assert.equal(await page.evaluate(() => window.calciumStep4.collectDraft().protocol.sd_estimator), newEstimator, 'a new analysis takes its noise from the whole trace when that is larger');
    assert.equal(await page.inputValue('#view'), 'raw');
    assert.ok(await page.isEnabled('#play') && await page.isEnabled('#time'), 'movie controls work before an analysis');
    assert.ok(await page.isHidden('#result-tiles') && await page.isHidden('#plot-card'), 'preview has no empty result dashboard');
    assert.deepEqual(await cellsOffered(), ['Cell001', 'Cell002', 'Cell003', 'Cell004', 'Cell005', 'Cell006'], 'cells can be selected before analysis');
    assert.equal(await page.getAttribute('#cell', 'placeholder'), 'None', 'and none is a choice: an empty box');
    assert.equal(await page.inputValue('#cell'), '', 'the page opens with no cell selected'); assert.ok(await page.isDisabled('#zoom-cell'));
    assert.ok(await page.locator('#windows-editor details').evaluateAll((details) => details.every((d) => !d.open)), 'advanced response rules start collapsed');
    assert.ok(await page.evaluate(() => {
      const transport = document.getElementById('transport').getBoundingClientRect();
      const picture = document.getElementById('picture').getBoundingClientRect();
      return transport.top >= picture.bottom - 1 && document.getElementById('picture-card').contains(document.getElementById('transport'));
    }), 'the frame slider is directly under the recording, in its card');
    // Each control where it is looked for: the cell's call on the plot card, the quality review as a chip in the summary row,
    // one Export menu per panel, and summary figures that are not a second cell filter.
    assert.ok(await page.evaluate(() => document.getElementById('plot-card').contains(document.getElementById('verdict')) && document.getElementById('plot-card').contains(document.getElementById('cell-quality'))));
    assert.ok(await page.evaluate(() => document.getElementById('result-tiles').contains(document.getElementById('quality-review'))));
    assert.equal(await page.locator('#result-tiles button.tile, #result-tiles [data-group]').count(), 0);
    assert.deepEqual(await page.locator('details.export').evaluateAll((menus) => menus.map((m) => m.id + ':' + m.querySelectorAll('button').length)), ['export-picture:4', 'export-plot:3', 'export-heatmap:3', 'export-metrics:2', 'export-calls:2', 'export-data:2', 'export-counts:2']);
    // (the mark form's Save keeps a review mark, not a file)
    assert.equal(await page.locator('main button', { hasText: /^(Save|Download)/ }).evaluateAll((buttons) => buttons.filter((b) => !b.closest('details.export, .save-options, .mark-form')).length), 0, 'every save button is in an Export menu');
    const initial = await page.evaluate(() => window.calciumStep4.collectDraft());
    const originalProtocol = initial.protocol;
    assert.equal(originalProtocol.schema_version, 2);
    assert.equal(originalProtocol.windows.length, 2);
    await screenshots('before');
    // Before an analysis a saved picture is exactly as the release drew it: white outlines, and a key that says "your cells".
    { const same = await againstRelease(); assert.ok(same.different === 0 && same.lit > 1000 && same.size[0] === 1000, 'a picture saved before an analysis is the release\'s, pixel for pixel: ' + JSON.stringify(same)); }
    assert.ok(await page.isHidden('#call-legend'), 'and under the picture nothing is said about outlines yet, as before');
    // Save frames… has two ways: "Every Nth frame", as it opens and as before, and "Baseline and peak", one picture of two frames.
    // Before an analysis both of its frames are the frame on screen.
    await slider(7); await page.locator('#export-picture > summary').click(); await page.click('#save-frames');
    assert.deepEqual(await page.locator('#frames-options [data-frames]').evaluateAll((buttons) => buttons.map((b) => b.textContent + (b.getAttribute('aria-pressed') === 'true' && b.classList.contains('on') ? '*' : ''))), ['Every Nth frame*', 'Baseline and peak']);
    assert.ok(await page.isVisible('#frames-every') && await page.isHidden('#frames-baseline') && await page.isHidden('#frames-peak') && await page.isHidden('#frames-download') && await text('#frames-go') === 'Save frames', 'Every Nth frame is as it was');
    // "Every" takes a whole number of frames that the recording has, as typed: "2.5" is not rounded to 3, nor "1e1" read as 10.
    // What it cannot use stays in the box, marked, with the frames there are beside it, and nothing is saved.
    { const frames = await page.evaluate(() => window.calciumStep4.S.info.frames), everySays = async () => [await page.inputValue('#frames-every'), await page.getAttribute('#frames-every', 'aria-invalid'), await page.isVisible('#frames-every-range') && await text('#frames-every-range'), await page.isDisabled('#frames-go'), await text('#frames-count')];
      for (const typed of ['0', '2.5', '1e1', 'abc', '', String(frames + 1)]) { await page.fill('#frames-every', typed); assert.deepEqual(await everySays(), [typed, 'true', '1 to ' + frames, true, ''], 'every "' + typed + '" frames is kept as typed, marked, and not saved'); }
      await page.fill('#frames-every', '10'); const said = await everySays(); assert.deepEqual(said.slice(0, 4), ['10', null, false, false]); assert.match(said[4], /^\d+ pictures and one with all of them, PNG, 300 dpi\.$/); }
    await page.click('[data-frames=pair]');
    assert.ok(await page.isHidden('#frames-every') && await page.isVisible('#frames-download') && await text('#frames-go') === 'Save picture');
    assert.deepEqual([await page.inputValue('#frames-baseline'), await page.inputValue('#frames-peak')], ['7', '7'], 'before an analysis both boxes hold the frame on screen');
    await page.click('[data-frames=every]'); await page.click('#frames-cancel'); assert.ok(await page.isHidden('#frames-options')); await slider(1);


    // Stimulus timing in frames: one number, "after frame N". F₀, the baseline and the response window follow it unless set separately,
    // and the seconds sent to the engine are the recording's own frame times.
    const editor0 = page.locator('.window-editor[data-index="0"]'), window0 = () => page.evaluate(() => window.calciumStep4.collectDraft().protocol.windows[0]);
    const times = await page.evaluate(() => window.calciumStep4.frameTimes());
    assert.equal(await page.inputValue('#stim-frame-0'), '6');
    assert.equal(await editor0.locator('.derived[data-linked="1"]').count(), 0, 'a loaded protocol with windows of its own opens with those rows set separately');
    assert.deepEqual([await page.inputValue('#base-from-0'), await page.inputValue('#base-to-0'), await page.inputValue('#resp-from-0'), await page.inputValue('#resp-to-0')], ['2', '6', '7', '15']);
    assert.match(await text('#suggestion-0'), /A placeholder from the frames used in Fiji\. Set it to when the stimulus was added/);
    for (const cell of ['#f0-cell', '#base-cell-0', '#resp-cell-0']) await page.locator(cell + ' .relink').click();
    assert.equal(await editor0.locator('.derived[data-linked="1"]').count(), 3);
    let edited = await window0();
    assert.ok(edited.Stimulus_onset_s === times[6] && edited.Baseline_start_s === times[0] && edited.Baseline_end_s === times[6] && edited.Start_s === times[6] && edited.End_s === times[29], 'rows that follow the stimulus use exact frame times');
    assert.equal(await page.inputValue('#baseline_start'), '1'); assert.equal(await page.inputValue('#baseline_end'), '6');
    assert.equal(await page.locator('#base-cell-0 .link-toggle').textContent(), 'same as F₀'); assert.equal(await page.locator('#resp-cell-0 .link-toggle').textContent(), '7–30 (end)');
    // One number moves everything.
    await page.fill('#stim-frame-0', '4'); edited = await window0();
    assert.ok(edited.Stimulus_onset_s === times[4] && edited.Baseline_end_s === times[4] && edited.Start_s === times[4] && edited.End_s === times[29]);
    assert.equal(await page.inputValue('#baseline_end'), '4'); assert.equal(await page.locator('#side .field-error').count(), 0); assert.equal(await editor0.locator('[aria-invalid]').count(), 0);
    assert.ok(await page.isHidden('#suggestion-0'));
    assert.equal(await text('#frames-preview-0'), 'In seconds: baseline 0–1.5 · stimulus 2 · response 2–19.5');
    assert.match(await text('#confirm-text'), /Timing checked: recorded timestamps, stimulus after frame 4/);
    // "Use frame" sets it to the frame on screen.
    await page.evaluate(() => window.calciumStep4.setFrame(9));
    assert.equal(await editor0.locator('.from-frame').textContent(), 'Use frame 9');
    await editor0.locator('.from-frame').click(); edited = await window0();
    assert.equal(await page.inputValue('#stim-frame-0'), '9'); assert.ok(edited.Stimulus_onset_s === times[9] && edited.Baseline_end_s === times[9] && edited.Start_s === times[9]);
    assert.equal(await page.locator('#side .field-error').count(), 0);
    // A row set separately that is wrong: one message beside it, naming the frames, with a one-click repair.
    await page.locator('#base-cell-0 .link-toggle').click(); await page.fill('#base-to-0', '12');
    assert.equal(await page.locator('#side .field-error').count(), 1);
    assert.match(await editor0.locator('.field-error').textContent(), /Stimulus: baseline frames 10–12 come after the stimulus\.\s*Use 1–9/);
    assert.equal(await page.getAttribute('#base-to-0', 'aria-invalid'), 'true'); assert.ok(await page.isDisabled('#analyze'));
    await page.waitForFunction(() => JSON.stringify(window.calciumStep4.collectDraft()) === window.calciumStep4.S.draftAck);
    assert.doesNotMatch(await text('#validation'), /Stimulus: baseline/, 'the same problem is not repeated in the list');
    await editor0.locator('.field-error .repair').click();
    assert.equal(await page.locator('#side .field-error').count(), 0); assert.equal((await window0()).Baseline_end_s, times[9]);
    // F₀ frames after the stimulus: a warning with its repair.
    await page.locator('#f0-cell .link-toggle').click(); await page.fill('#baseline_end', '11');
    const under = (id) => page.evaluate((field) => { const note = document.getElementById(field).closest('.row').nextElementSibling; return note ? [note.className, note.textContent] : null; }, id);
    assert.deepEqual(await under('f0-cell'), ['field-warning', 'F₀ frames 10–11 come after the stimulus.Use 1–9'], 'a warning sits under the field it is about');
    assert.match((await under('base-from-1'))[1], /^Control: the baseline holds 1 frame but at least 3 are needed/, 'the control\'s warning is under the control\'s baseline');
    assert.equal(await text('#validation'), '2 warnings above', 'the Run box only says how many there are');
    await page.locator('#side .field-warning .repair', { hasText: 'Use 1–9' }).click();
    assert.equal(await page.inputValue('#baseline_end'), '9'); assert.equal((await under('f0-cell') || [''])[0] === 'field-warning', false); assert.equal(await text('#validation'), '1 warning above');
    // Too few baseline frames would leave every cell without a call: the run is stopped and the reason is given.
    await page.locator('#base-cell-0 .relink').click(); await page.fill('#stim-frame-0', '2');
    assert.match(await editor0.locator('.field-error').textContent(), /Stimulus: the baseline holds 2 frames but at least 3 are needed\. No cell would get a call\./);
    assert.ok(await page.isDisabled('#analyze'));
    // The windows are drawn under the frame slider.
    await page.fill('#stim-frame-0', '6');
    assert.ok(await page.locator('#time-windows .baseline').count() === 1 && await page.locator('#time-windows .scoring').count() === 1 && await page.locator('#time-windows .stimulus').count() === 1);
    await page.click('#discard');
    await page.waitForFunction((wanted) => JSON.stringify(window.calciumStep4.collectDraft()) === wanted, JSON.stringify(initial));
    assert.equal(await page.inputValue('#stim-frame-0'), '6'); assert.equal(await editor0.locator('.derived[data-linked="1"]').count(), 0, 'discard restores the values and which rows are set separately');

    // Clearing the time source must not reuse the initial session clock.
    const openAdvanced = () => page.evaluate(() => { document.getElementById('protocol-details').open = true; });   // timestamps, rule type, protocol file, diagnostics, log
    await openAdvanced();
    await page.click('#timestamps-clear'); await page.fill('#sampling_hz', '');
    await page.waitForFunction(() => window.calciumStep4.S.state.frameTimes === null);
    await page.waitForFunction(() => /20 frames per second/.test(document.getElementById('playback-note').textContent));
    assert.doesNotMatch(await text('#time-label'), /\d(?:\.\d+)? s/);
    for (const invalid of ['0', '-1', 'abc', '', '2abc']) {
      await page.fill('#sampling_hz', invalid);
      assert.ok(await page.isDisabled('#analyze'));
      assert.ok((await page.locator('#side .field-error').allTextContents()).some((s) => /Frame rate/.test(s)), 'the frame rate message is under the frame rate');
      assert.match(await text('#validation'), /^Fix \d+ problems? above/);
      assert.equal(await page.getAttribute('#sampling_hz', 'aria-invalid'), 'true');
    }
    await page.fill('#sampling_hz', '2');
    await page.waitForFunction(() => window.calciumStep4.S.state.frameTimes && window.calciumStep4.S.state.frameTimes[15] === 7.5);
    await page.setInputFiles('#timestamps-file', { name: initial.timestamps.name, mimeType: 'text/csv', buffer: Buffer.from(initial.timestamps.content) });
    await page.waitForFunction(() => window.calciumStep4.S.state.frameTimes && window.calciumStep4.S.state.frameTimes[15] === 12.5);

    // Reopening the same session must continue its revision counter, or new edits are ignored.
    const priorRevision = await page.evaluate(() => window.calciumStep4.S.state.draftRevision);
    assert.ok(priorRevision > 0);
    await page.reload(); await page.waitForFunction(() => window.calciumStep4 && window.calciumStep4.S.state); await waitFrame(1);
    assert.ok(await page.evaluate((previous) => window.calciumStep4.S.revision >= previous, priorRevision));
    await page.click('#optional-settings > summary'); await page.fill('#condition', 'Reopened session');
    await page.waitForFunction((previous) => window.calciumStep4.S.state.draftRevision > previous, priorRevision);
    assert.equal(await page.evaluate(async () => (await (await fetch('api/session')).json()).draft.values.condition), 'Reopened session');
    const savedRevision = await page.evaluate(() => window.calciumStep4.S.state.draftRevision);
    await page.fill('#condition', '');
    await page.waitForFunction((previous) => window.calciumStep4.S.state.draftRevision > previous, savedRevision);

    // An actual pointer drag on the range control changes the image before analysis.
    await page.locator('#time').scrollIntoViewIfNeeded();
    const range = await page.locator('#time').boundingBox();
    await page.mouse.move(range.x + 8, range.y + range.height / 2); await page.mouse.down();
    await page.mouse.move(range.x + 8 + (range.width - 16) * 11 / 29, range.y + range.height / 2, { steps: 8 }); await page.mouse.up();
    const dragged = Number(await page.inputValue('#time'));
    assert.ok(dragged >= 11 && dragged <= 13, 'dragging selects the intended movie region');
    await waitFrame(dragged); assert.ok(await gray() > 200, 'the responding square is bright in the selected frame');
    await slider(3); assert.ok(await gray() < 60, 'scrubbing back changes the actual image to the resting frame');
    assert.match(await text('#time-label'), /1\.00 s.*Frame 3 \/ 30/);

    // Playback: 10× by default, on a clock in recording time, from a buffer that makes scrubbing back free. (Waits for the movie
    // are given 30 s: on a busy computer the frames come slowly, and what is checked is that they come, not how fast; the 10×
    // timing below is the one check of speed.)
    assert.equal(await page.inputValue('#speed'), '10'); const movieWait = { timeout: 30000 };
    let frameRequests = 0; const countFrames = (request) => { if (request.url().includes('/api/frame?')) frameRequests++; }; page.on('request', countFrames);
    await slider(1); const started = Date.now(); await page.click('#play', movieWait);
    await page.waitForFunction(() => window.calciumStep4.S.timer === null, null, { timeout: 15000 });
    const took = (Date.now() - started) / 1000;
    assert.equal(await page.inputValue('#time'), '30', 'playback stops at the last frame');
    assert.ok(took > 1.1 && took < 4, '14.5 s of recording at 10× plays in about 1.45 s, not ' + took + ' s');
    await waitFrame(30); frameRequests = 0;
    for (const f of [24, 17, 9, 2]) await slider(f);
    assert.equal(frameRequests, 0, 'scrubbing back over frames already seen makes no requests');
    page.off('request', countFrames);
    await page.selectOption('#speed', '50', movieWait); await page.click('#play', movieWait);
    await page.waitForFunction(() => /playing every 4th frame/.test(document.getElementById('picture-note').textContent), null, movieWait);
    await page.click('#play', movieWait);
    assert.ok(!(await page.evaluate(() => window.calciumStep4.S.timer)));
    // A manual seek stops the timer.
    await slider(3);
    await page.selectOption('#speed', '0.5', movieWait); await page.click('#play', movieWait);
    await page.waitForFunction(() => window.calciumStep4.S.frame >= 5, null, movieWait);
    await slider(12);
    assert.equal(await page.evaluate(() => window.calciumStep4.S.timer), null);
    assert.equal(await text('#play'), '▶ Play');
    await page.waitForTimeout(250);
    assert.equal(await page.inputValue('#time'), '12', 'Play must not overwrite a manual seek');
    await page.click('#frame-back'); await waitFrame(11);
    await page.click('#frame-forward'); await waitFrame(12);
    await page.selectOption('#view', 'picture');
    assert.ok(await page.isDisabled('#play') && await page.isDisabled('#time'));
    assert.match(await text('#playback-note'), /Static colour overview.*does not move/);
    await page.selectOption('#view', 'raw'); await waitFrame(12);

    // A slow old response cannot replace the frame now requested by the slider.
    await slider(3);
    await page.evaluate(() => window.calciumStep4.resetFrames());   // frames already seen are buffered; these checks need fresh requests
    const delayedFrame = '**/api/frame?f=8&mode=raw&*';
    await page.route(delayedFrame, async (route) => {
      const response = await route.fetch();
      await new Promise((resolve) => setTimeout(resolve, 200)); await route.fulfill({ response });
    });
    await page.evaluate(() => {
      window.browserMovieNotes = [];
      window.browserMovieObserver = new MutationObserver(() => window.browserMovieNotes.push(document.getElementById('picture-note').textContent));
      window.browserMovieObserver.observe(document.getElementById('picture-note'), { childList: true, subtree: true });
    });
    const slowRequest = page.waitForRequest((request) => request.url().includes('/api/frame?f=8&mode=raw&'));
    await page.locator('#time').evaluate((element) => { element.value = 8; element.dispatchEvent(new Event('input')); });
    await slowRequest; await slider(20); await page.unroute(delayedFrame);
    assert.ok(await page.evaluate(() => {
      window.browserMovieObserver.disconnect(); return !window.browserMovieNotes.some((s) => /^Frame 8\b/.test(s));
    }), 'a stale frame must not be labeled as the current image');
    assert.ok(await gray() < 60);

    // A failed frame fetch is visible and the next manual seek recovers.
    await page.evaluate(() => window.calciumStep4.resetFrames());
    const failedFrame = '**/api/frame?f=9&mode=raw&*'; injectedFrameFailure = true;
    await page.route(failedFrame, (route) => route.fulfill({ status: 503, body: 'Synthetic unavailable frame' }));
    let askedForFailedFrame = 0; const countFailedFrame = (request) => { if (/\/api\/frame\?f=9&mode=raw&/.test(request.url())) askedForFailedFrame++; };
    page.on('request', countFailedFrame);
    await page.locator('#time').evaluate((element) => { element.value = 9; element.dispatchEvent(new Event('input')); });
    await page.waitForFunction(() => /could not load/.test(document.getElementById('movie-message').textContent), null, { polling: 100 });
    assert.equal(await page.evaluate(() => window.calciumStep4.S.timer), null);
    assert.match(await text('#picture-note'), /Last displayed frame: 20/);
    await page.waitForTimeout(400); page.off('request', countFailedFrame);
    assert.equal(askedForFailedFrame, 1, 'a frame that cannot load is not asked for again and again');
    // nothing is on its way when the route is removed: a request caught at that moment is never answered
    await page.waitForFunction(() => window.calciumStep4.Buf.inflight.size === 0, null, { polling: 100 });
    await page.unroute(failedFrame); await slider(10); injectedFrameFailure = false;
    assert.equal(await page.evaluate(() => window.calciumStep4.Pic.error), null);

    // The recording closed in Fiji: Fiji refuses every frame in plain words, the movie says them, and frames come back once it
    // is open again. A frame already on its way when the recording closes is let through first; then every frame is refused.
    const closedText = 'The recording is closed in Fiji. Open it again to go on.';
    await page.waitForFunction(() => window.calciumStep4.Buf.inflight.size === 0 && window.calciumStep4.Buf.decoding.size === 0, null, { polling: 100 });
    const anyFrame = '**/api/frame?*'; injectedFrameFailure = true;
    await page.route(anyFrame, (route) => route.fulfill({ status: 410, contentType: 'application/json', body: JSON.stringify({ error: closedText }) }));
    await page.evaluate(() => window.calciumStep4.resetFrames());
    const askedForClosed = page.waitForRequest((request) => /\/api\/frame\?f=11&mode=raw&/.test(request.url()));
    await page.locator('#time').evaluate((element) => { element.value = 11; element.dispatchEvent(new Event('input')); });
    await askedForClosed;
    await page.waitForFunction((said) => document.getElementById('movie-message').textContent === said, closedText, { polling: 100 });
    assert.equal(await page.evaluate(() => window.calciumStep4.S.frame), 11);
    assert.ok(!/java|Exception/.test(await text('#movie-message')));
    await page.waitForFunction(() => window.calciumStep4.Buf.inflight.size === 0, null, { polling: 100 });
    await page.unroute(anyFrame); await slider(10); injectedFrameFailure = false;
    assert.equal(await page.evaluate(() => window.calciumStep4.Buf.closed), '', 'a frame that loads again clears the message');
    assert.equal(await page.evaluate(() => window.calciumStep4.Pic.error), null);

    // Portable JSON retains precise untouched values and imported timing overrides.
    await openAdvanced();
    const precise = JSON.parse(JSON.stringify(originalProtocol));
    precise.name = 'Browser protocol round trip'; precise.validation_status = 'user_defined';
    for (const window of precise.windows) { window.Protocol_name = precise.name; window.Validation_status = precise.validation_status; }
    precise.windows[0].Floor_local_fraction = 0.12345678901234567;
    const first = precise.windows[0];
    precise.field_overrides = [{ Field: 'Fixture recording', Role: 'Stimulus', Stimulus: first.Stimulus,
      Baseline_start_s: first.Baseline_start_s, Baseline_end_s: first.Baseline_end_s,
      Stimulus_onset_s: first.Stimulus_onset_s, Start_s: first.Start_s, End_s: first.End_s }];
    await page.check('#confirmed'); await importProtocol(precise, 'precise-protocol.json');
    assert.ok(!(await page.isChecked('#confirmed')), 'an imported protocol with other timing needs timing confirmation');
    assert.match(await text('#status'), /Imported “Browser protocol round trip”\. Its timing replaced yours/, 'an import says what it did');
    const protocolDownload = page.waitForEvent('download'); await page.click('#protocol-export');
    const downloaded = await protocolDownload;
    assert.match(downloaded.suggestedFilename(), /^Fixture_recording_response-protocol_\d{8}-\d{6}\.json$/, 'exported protocol names say which recording and when');
    const portablePath = path.join(out, 'roundtrip-protocol.json'); await downloaded.saveAs(portablePath);
    const portable = JSON.parse(fs.readFileSync(portablePath, 'utf8'));
    assert.equal(portable.windows[0].Floor_local_fraction, precise.windows[0].Floor_local_fraction, 'display formatting must not round an untouched protocol');
    assert.deepEqual(portable.field_overrides, precise.field_overrides);
    await page.check('#confirmed'); await importProtocol(portable, 'saved-protocol.json');
    assert.ok(await page.isChecked('#confirmed'), 'importing the same timing again keeps the confirmation'); assert.match(await text('#status'), /The timing is unchanged/);
    assert.equal(await page.evaluate(() => window.calciumStep4.collectDraft().protocol.windows[0].Floor_local_fraction), precise.windows[0].Floor_local_fraction);
    await importProtocol(originalProtocol, 'original-protocol.json');

    // A version 2 rule made on the page from a legacy one is a protocol made here: it takes the noise of a new protocol,
    // and is exported with the engine version that goes with it. The legacy rule itself is exported naming no estimator.
    const legacyProtocol = { schema: 'simpleca.response-protocol', schema_version: 1, engine_version: 'legacy-1', name: 'Legacy consecutive-sample protocol', validation_status: 'legacy_unvalidated', signal_reference: 'fixed_reference_fractional_dff0',
      windows: [{ Role: 'Stimulus', Stimulus: 'Capsaicin', Start_s: 3, End_s: 7, Baseline_s: 2.5, Floor_local_fraction: 0.1, SD_multiplier: 3, Consecutive_samples: 2, Max_gap_s: null }], field_overrides: null };
    await importProtocol(legacyProtocol, 'legacy-protocol.json');
    assert.equal(await page.inputValue('#protocol-mode'), '1');
    const exportProtocol = async (name) => { const saving = page.waitForEvent('download'); await page.click('#protocol-export'); const file = path.join(out, name); await (await saving).saveAs(file); return JSON.parse(fs.readFileSync(file, 'utf8')); };
    const legacyExported = await exportProtocol('legacy-exported.json');
    assert.equal(legacyExported.schema_version, 1); assert.ok(!('sd_estimator' in legacyExported), 'a legacy file names no estimator');
    await page.selectOption('#protocol-mode', '2');
    const madeFromLegacy = await page.evaluate(() => window.calciumStep4.collectDraft().protocol);
    assert.equal(madeFromLegacy.schema_version, 2); assert.equal(madeFromLegacy.sd_estimator, newEstimator, 'the rule type change carries the estimator of a new protocol');
    const madeExported = await exportProtocol('made-from-legacy.json');
    assert.equal(madeExported.sd_estimator, newEstimator); assert.equal(madeExported.engine_version, '2.1.0');
    assert.deepEqual(madeExported.windows.map((w) => [w.SD_estimator, w.Engine_version]), [[newEstimator, '2.1.0']]);
    await importProtocol(originalProtocol, 'original-protocol.json');

    // Legacy is an explicit mode and disables elapsed-time-only fields.
    const beforeSwitch = await page.evaluate(() => window.calciumStep4.collectDraft().protocol);
    await page.selectOption('#protocol-mode', '1'); await page.selectOption('#protocol-mode', '2');
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.collectDraft().protocol), beforeSwitch, 'rule type round-trip preserves persistence, metadata and thresholds');
    await page.selectOption('#protocol-mode', '1');
    assert.ok(await page.isDisabled('#window-0-Duration_mode') && await page.isDisabled('#window-0-Min_duration_s'));
    assert.ok(await page.isEnabled('#window-0-Consecutive_samples'));
    assert.equal(await page.inputValue('#window-1-Stimulus'), 'KCl');
    // In the legacy rule the baseline always ends, and the response window always starts, at the stimulus.
    await page.fill('#stim-frame-0', '10');
    assert.equal(await page.inputValue('#window-0-Baseline_end_s'), '5'); assert.equal(await page.inputValue('#window-0-Start_s'), '5'); assert.equal(await page.inputValue('#window-0-Stimulus_onset_s'), '5');
    assert.ok(await page.isDisabled('#base-to-0') && await page.isDisabled('#resp-from-0'));
    await page.uncheck('#control'); await page.check('#control');
    assert.equal(await page.inputValue('#window-1-Duration_mode'), 'samples');
    assert.equal(await page.inputValue('#window-1-Start_s'), '');
    assert.ok(await page.isDisabled('#analyze'));
    // A control that was just ticked asks for its one number, in one message beside that box.
    assert.match(await page.locator('.window-editor[data-index="1"] .field-error').textContent(), /Control: enter the frame after which it was added\./);
    assert.equal(await page.locator('.window-editor[data-index="1"] .field-error').count(), 1); assert.doesNotMatch(await text('#validation'), /needs a number/);
    await importProtocol(originalProtocol, 'analysis-protocol.json');

    // The real HTTP submission carries v2 and the timestamp snapshot.
    await page.click('#analyze'); assert.match(await text('#status'), /Timing checked/);
    await page.check('#confirmed'); await page.click('#analyze');
    await page.waitForFunction(() => window.calciumStep4.S.results && !window.calciumStep4.S.state.busy, null, { timeout: 15000 });
    await page.waitForFunction(() => document.getElementById('review-stage').textContent === 'Completed analysis');
    const sent = JSON.parse(fs.readFileSync(path.join(out, 'sent.json'), 'utf8'));
    const submitted = JSON.parse(sent.protocol_json);
    assert.equal(submitted.schema_version, 2);
    // The page writes its fixed "2.0.0" and "sample_sd" into the draft; the run takes the estimator named at the top of it.
    assert.equal(sent.sd_estimator, newEstimator); assert.equal(submitted.sd_estimator, newEstimator); assert.equal(submitted.engine_version, '2.1.0');
    assert.deepEqual(submitted.windows.map((w) => w.SD_estimator), [newEstimator, newEstimator]);
    assert.deepEqual(submitted.windows.map((w) => w.Stimulus), ['Capsaicin', 'Ionomycin']);
    assert.equal(submitted.windows[0].Floor_local_fraction, 0.1);
    assert.equal(submitted.windows[1].Floor_local_fraction, 0.4);
    assert.match(sent.timestamps_sha256, /^[a-f0-9]{64}$/);
    assert.equal(sent.input_mode, 'background_corrected_fluorescence');
    assert.match(await text('#tile-responding'), /66\.7%.*2 of 3 evaluable/);
    assert.equal((await text('#tile-registered')).replace(/\D/g, ''), '6');
    // All cells lead the headline; the quiet cells follow on a line of their own; every figure comes from the tally worked out in Fiji.
    // A share of cells reads one way everywhere on the page: one decimal and "%" (2 of 3 is 66.7%, not 67%).
    assert.match(await page.locator('#tile-responding > b').textContent(), /^66\.7%$/);
    assert.match(await page.locator('#tile-responding .tile-line').last().textContent(), /^Quiet enough to detect a 30 % rise: \d+\.\d% \(\d+ of \d+\)$/);
    // The cell that went below zero stays counted until the switch leaves it out; then every figure follows, and comes back.
    const counted = await text('#tile-responding'), tally = await page.evaluate(() => window.calciumStep4.S.results.tally);
    assert.deepEqual(tally.belowZero, ['Cell005']); assert.equal(await text('#below-zero-text'), 'Leave out 1 cell below zero'); assert.ok(await page.isVisible('#below-zero-label') && !(await page.isChecked('#below-zero')));
    const allRow = tally.rows.find((row) => row.key === 'all'); assert.deepEqual(allRow.counted.slice(0, 2), [2, 3]);
    await page.check('#below-zero');
    assert.match(await text('#tile-responding'), new RegExp('Capsaicin: ' + allRow.without[0] + ' of ' + allRow.without[1] + ' evaluable cells · 1 below zero left out'));
    assert.equal(await page.locator('#tile-responding > b').textContent(), (100 * allRow.without[0] / allRow.without[1]).toFixed(1) + '%');
    await page.uncheck('#below-zero'); assert.equal(await text('#tile-responding'), counted);
    // The rule line says what is applied: samples in a row, also for a length typed in seconds.
    assert.match(await text('#rule-summary-0'), /for 2 samples in a row$/);
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4.samplesFor, at = (step) => Array.from({ length: 20 }, (_, i) => i * step); return [f(1, at(1)), f(1, at(0.999946296215057)), f(1, at(0.5)), f(1, at(0.19999999999999996)), f(0.4, at(0.5)), f(1, null), f(0, at(1))]; }), [2, 3, 3, 6, 2, null, null]);
    assert.equal(await page.inputValue('#view'), 'raw', 'analysis preserves the selected movie view');
    await openCell('Cell001'); assert.match(await text('#verdict'), /Cell001: Responded/);
    // The duration notice is set on every cell by the engine but belongs to the recording: it is shown once and flags no cell.
    assert.deepEqual(await page.locator('#group option').allTextContents(), ['All cells (6)', 'Responded (2)', 'Did not respond (1)', 'Excluded (2)', 'Unscorable (1)', 'Needs quality review (1)']);
    assert.equal(await text('#quality-summary'), '1 cell flagged'); assert.ok(await page.isHidden('#quality-list'), 'the chip does not push the data down');
    await openCell('Cell005'); assert.match(await text('#cell-quality'), /Low F₀/);
    await page.locator('#quality-review summary').click(); assert.match(await text('#quality-list'), /Applies to the whole recording: Minimum response duration is shorter than one frame interval/);
    assert.equal(await text('#quality-lead'), '1 cell flagged, including 1 responder.');
    assert.ok(await page.evaluate(() => { const tabs = document.querySelector('.tabs').getBoundingClientRect(), row = document.getElementById('result-tiles').getBoundingClientRect(); return tabs.top - row.bottom < 40; }), 'the open list lies over the page');
    // Every flagged cell is on a row of its own, with its call and what it is flagged for; under the rows, each flag once.
    assert.deepEqual(await page.locator('#table-flagged tbody tr').evaluateAll((rows) => rows.map((r) => [...r.cells].map((c) => c.textContent))), [['Cell005', 'Responded', 'Low F₀ · Below zero']]);
    assert.ok(await page.locator('#table-flagged tbody tr.on').count() === 1, 'the open cell is marked in the list');
    assert.deepEqual(await page.locator('#quality-list .flag-filter').allTextContents(), ['Low F₀: 1 cell, 1 responder', 'Below zero: 1 cell, 1 responder']);
    assert.match(await page.locator('#quality-list .flag-what').first().textContent(), /^Low F₀ relative to baseline noise/);
    assert.equal(await page.locator('#quality-filter').count(), 0, 'the Show list is the one place to show the flagged cells');
    await page.locator('#quality-list .flag-filter').first().click();
    assert.equal(await page.inputValue('#group'), '5'); assert.deepEqual(await cellsOffered(), ['Cell005']);
    await page.locator('#quality-list .flag-filter').first().click(); assert.equal(await page.evaluate(() => window.calciumStep4.S.qualityFlag), null);
    // A click on a flagged cell's row opens that cell: Show goes to the flagged cells, so the arrows step through them, and the list closes.
    await page.selectOption('#group', '0'); await openCell(''); await page.locator('#quality-review summary').click(); await page.locator('#quality-review summary').click();
    await page.locator('#table-flagged tbody tr').first().click();
    assert.deepEqual(await page.evaluate(() => { const a = window.calciumStep4; return { cell: a.S.cell, group: document.getElementById('group').value, open: document.getElementById('quality-review').open, box: document.getElementById('cell').value }; }), { cell: 'Cell005', group: '5', open: false, box: 'Cell005' });
    assert.equal(await text('.tab.on'), 'Cells and traces'); await page.locator('#quality-review summary').click();
    await page.selectOption('#group', '0');
    await page.locator('#quality-review summary').click();
    await openCell('Cell001');
    await screenshots('after');
    // Low windows and the connection banner must remain scrollable without hiding the caption.
    await page.setViewportSize({width:1440,height:739});
    await page.evaluate(() => document.getElementById('lost').classList.add('on'));
    await settledPaint();
    assert.ok(await page.evaluate(() => document.documentElement.scrollHeight > innerHeight));
    // Bring the caption to the middle of the window, wherever the rows above it put it: "if needed" left it where it was when it
    // sat exactly on the lower edge (seen on Linux once the headline gained a line), and the check below then failed by a fraction.
    await page.locator('#picture-note').evaluate((e) => e.scrollIntoView({ block: 'center' })); await settledPaint();
    assert.ok(await page.locator('#picture-note').evaluate((e) => {const a=e.getBoundingClientRect(),b=e.closest('.card').getBoundingClientRect();return a.top>=0&&a.bottom<=b.bottom&&a.bottom<=innerHeight;}), 'the caption can be scrolled fully into view and stays inside its card');
    await page.screenshot({path:path.join(out,'page-short-window.png')});
    await page.evaluate(() => document.getElementById('lost').classList.remove('on'));
    await page.setViewportSize({width:1440,height:900});


    // Control calls, thresholds and headline counts follow the selected response window.
    await page.selectOption('#response-window', 'Control::Ionomycin');
    assert.match(await text('#verdict'), /Cell001: Did not respond/);
    assert.ok(await page.evaluate(() => Math.abs(window.calciumStep4.referenceLevels().find((r) => r.name === 'required rise').value - 1.1) < 1e-12), 'local .4 threshold on baseline .5 belongs at global 1.1');
    assert.match(await text('#tile-responding'), /Ionomycin \(control\): 1 of 3 evaluable/);
    assert.ok(await page.isHidden('#below-zero-label'), 'the switch belongs to the primary stimulus');
    await page.selectOption('#group', '1');
    assert.deepEqual(await cellsOffered(), ['Cell002']);
    await page.selectOption('#response-window', 'Stimulus::Capsaicin');
    assert.match(await text('#tile-responding'), /Capsaicin: 2 of 3 evaluable/);
    assert.deepEqual(await cellsOffered(), ['Cell001', 'Cell005']);
    await page.selectOption('#group', '3');
    assert.deepEqual(await cellsOffered(), ['Cell003', 'Cell004']);
    await page.selectOption('#group', '4');
    assert.deepEqual(await cellsOffered(), ['Cell006']);
    await page.selectOption('#group', '0'); await openCell('Cell002');

    // ---------------------------------------------------------------- the time course and its legend
    // One style table and one painter for the page and for a saved figure; one function for the cells a plot uses; one place for the legend's words.
    const plotNow = () => page.evaluate(() => { const g = window.calciumStep4.Plot.geom, l = g.legend; return { L: g.L, T: g.T, W: g.W, H: g.H, lo: g.lo, hi: g.hi, legend: g.items.map((i) => i.text),
      rows: l ? l.layout.height / l.layout.rowHeight : 0, legendAt: l ? [l.x, l.y] : null }; });
    // What the painter puts on its canvas, in order: every text, filled box and stroked line (with the corners of the path it
    // strokes), and where its parts went. true: the saved figure at 96 pixels per inch; false: the page's own plot, painted again
    // (the layer under the pointer's own line and numbers).
    const drawing = (figure) => page.evaluate((figure) => {
      const f = window.calciumStep4, P = CanvasRenderingContext2D.prototype, was = {}, seen = { texts: [], boxes: [], lines: [] }, paths = new WeakMap(); for (const name of ['fillText', 'fillRect', 'stroke', 'beginPath', 'moveTo', 'lineTo', 'rect']) was[name] = P[name];
      const at = (ctx, x, y) => { const m = ctx.getTransform(); return [m.a * x + m.c * y + m.e, m.b * x + m.d * y + m.f]; }, wanted = (ctx) => figure || ctx.canvas === f.Plot.layer;
      const reach = (ctx, x, y) => { const p = at(ctx, x, y), b = paths.get(ctx) || [p[0], p[1], p[0], p[1]]; paths.set(ctx, [Math.min(b[0], p[0]), Math.min(b[1], p[1]), Math.max(b[2], p[0]), Math.max(b[3], p[1])]); };
      P.beginPath = function () { paths.delete(this); return was.beginPath.apply(this, arguments); };
      P.moveTo = function (x, y) { reach(this, x, y); return was.moveTo.apply(this, arguments); }; P.lineTo = function (x, y) { reach(this, x, y); return was.lineTo.apply(this, arguments); };
      P.rect = function (x, y, w, h) { reach(this, x, y); reach(this, x + w, y + h); return was.rect.apply(this, arguments); };
      P.fillText = function (text, x, y) { if (wanted(this)) seen.texts.push(window.textSeen(this, text, x, y, at)); return was.fillText.apply(this, arguments); };
      P.fillRect = function (x, y, w, h) { const p = at(this, x, y); if (wanted(this)) seen.boxes.push({ x: p[0], y: p[1], w, h, ink: this.fillStyle }); return was.fillRect.apply(this, arguments); };
      P.stroke = function () { const b = paths.get(this) || [0, 0, 0, 0]; if (wanted(this)) seen.lines.push({ ink: this.strokeStyle, width: this.lineWidth, alpha: this.globalAlpha, dashed: this.getLineDash().length > 0, x0: b[0], y0: b[1], x1: b[2], y1: b[3] }); return was.stroke.apply(this, arguments); };
      try { if (figure) f.figureCanvas('plot', 96); else window.dispatchEvent(new Event('resize')); } finally { Object.assign(P, was); }
      seen.texts = window.asWritten(seen.texts); const d = figure ? f.Fig.drawn : f.Plot.geom, l = d.legend, top = figure ? d.top : 0;
      return Object.assign(seen, { top, L: d.L, T: d.T, W: d.W, H: d.H, lo: d.lo, hi: d.hi, x0: d.x0, x1: d.x1, ink: figure ? f.Fig.plot.ink : '#3e4c59', legend: d.items.map((i) => i.text), items: d.items.map((i) => ({ text: i.text, kind: i.kind, colour: i.colour, dashed: !!i.dash })),
        places: l ? l.layout.places.map((place) => [place.x, place.y]) : [], wrapped: l ? l.layout.places.map((place) => place.lines) : [], legendAt: l ? [l.x, l.y] : null, legendWidth: l ? l.layout.width : 0, legendBox: l ? [l.x, top + l.y, l.layout.width, l.layout.height] : null });
    }, figure);
    const figureStyle = (changes) => page.evaluate((changes) => { const f = window.calciumStep4; Object.assign(f.Fig.plot, f.STYLE.plot, changes); }, changes || {});   // the built-in look, with these entries changed
    // No extra words on a plot. Outside its legend the painter writes its two axis titles, its tick labels and "stimulus", and nothing
    // else; in the legend, the entries and nothing else. And entries and marks go one to one: everything drawn in a colour of its own
    // (a line, dashed or not; a band; a bar) has one entry, and no entry is without its mark. The axes, the zero line, the casing
    // under a line and the stimulus line (named above it) are not marks; the tick at the start of the response band belongs to it.
    const inLegend = (seen, x, y) => !!seen.legendBox && x >= seen.legendBox[0] - 1 && x <= seen.legendBox[0] + seen.legendBox[2] + 1 && y >= seen.legendBox[1] - 1 && y <= seen.legendBox[1] + seen.legendBox[3] + 1;
    const wordsOf = (seen) => seen.texts.filter((t) => !inLegend(seen, t.x, t.y) && !/^−?\d+(\.\d+)?$/.test(t.text)).map((t) => t.text).sort();
    const marksOf = (seen) => { const marks = new Set();
      for (const s of seen.lines) if (!inLegend(seen, s.x0, s.y0) && ![seen.ink, '#9aa5b1', '#b45f06', '#ffffff'].includes(s.ink)) marks.add(s.alpha < 1 ? 'single cells' : s.x0 === s.x1 ? 'box #d5eedd' : s.ink + (s.dashed ? ' dashed' : ''));
      for (const b of seen.boxes) if (!inLegend(seen, b.x, b.y) && b.ink !== '#ffffff') marks.add('box ' + b.ink);
      return [...marks].sort(); };
    const entriesOf = (seen) => seen.items.map((item) => item.text === 'single cells' ? item.text : item.kind === 'line' ? item.colour + (item.dashed ? ' dashed' : '') : 'box ' + item.colour).sort();
    const clean = (seen, view) => {
      assert.deepEqual(wordsOf(seen), ['Time (s)', 'stimulus', 'ΔF/F₀'], view + ': outside the legend, only the axis titles, the tick labels and "stimulus" are written');
      assert.deepEqual(seen.texts.filter((t) => inLegend(seen, t.x, t.y)).map((t) => t.text), seen.wrapped.flat(), view + ': the legend holds its entries and nothing else');
      assert.deepEqual(seen.wrapped.map((lines) => lines.join(' ')), seen.legendBox ? seen.legend : [], view + ': every entry whole, on one line or more');
      assert.deepEqual(marksOf(seen), entriesOf(seen), view + ': legend entries and drawn marks, one to one');
      return seen; };
    const cleanBoth = async (view) => { clean(await drawing(false), view + ', on the page'); return clean(await drawing(true), view + ', in the figure'); };
    await openCell(''); await page.click('[data-plot=mean]');   // the average alone, as the page opens
    // The built-in look: 7 × 4 in, left and bottom axes, no grid, no title, black ink in a figure.
    assert.deepEqual(await page.evaluate(() => { const s = window.calciumStep4.STYLE.plot; return [s.width, s.height, s.frame, s.grid, s.ink, s.title, s.legend]; }), [7, 4, 'open', false, '#000000', null, 'top']);
    // The legend: each group once, in few words, with the mean line as its swatch and its number of cells; "baseline" and "response window" as two small entries.
    let now = await plotNow(), fig = await drawing(true), onPage = await drawing(false);
    assert.deepEqual(now.legend, ['all cells ± SEM (n = 4)', 'baseline', 'response window']);
    assert.equal(now.rows, 1); assert.equal(now.legendAt[0], now.L, 'the legend is one row that starts on the Y axis');
    assert.equal(await page.locator('#legend').count(), 0, 'the legend is drawn in the plot; there is no second one above it');
    assert.match(await page.getAttribute('#plot', 'aria-label'), /against time: all cells ± SEM \(n = 4\); baseline; response window$/, 'and a screen reader is told what it says');
    // The page and the saved figure are one drawing: the same legend, the same range, each in its own ink.
    assert.deepEqual(onPage.legend, fig.legend); assert.deepEqual([onPage.lo, onPage.hi], [fig.lo, fig.hi]);
    assert.deepEqual([fig.legendAt[0], new Set(fig.places.map((place) => place[1])).size], [fig.L, 1], 'in the figure, too, the legend is one row from the Y axis');
    const named = (seen, words) => seen.texts.filter((t) => words.includes(t.text));
    assert.deepEqual(named(onPage, ['Time (s)', 'ΔF/F₀', 'baseline']).map((t) => t.ink), ['#3e4c59', '#3e4c59', '#3e4c59'], 'the page writes in its own ink');
    assert.ok(fig.texts.length > 8 && fig.texts.filter((t) => t.text !== 'stimulus').every((t) => t.ink === '#000000'), 'a saved figure has black text');
    assert.ok(!fig.texts.some((t) => t.text === 'Fixture recording'), 'and no title until one is typed');
    assert.ok(!fig.lines.some((l) => l.ink === '#dee4ea'), 'no grid lines');
    clean(fig, 'the average, in the figure'); clean(onPage, 'the average, on the page');
    // Axes: the Y axis ends on a tick or half a tick step beyond the last one, so the numbers are round and little of the plot is
    // empty (the mean ± SEM runs from 0 to 1: with 5 % of room, half a step of 0.2 at either end); the time axis ends where the data end.
    const yTicks = (seen) => seen.texts.filter((t) => /^[-−]?[\d.]+$/.test(t.text) && t.x < seen.L).sort((a, b) => a.y - b.y);
    const beyond = (seen) => { const t = yTicks(seen), step = t[1].y - t[0].y; return [(t[0].y - seen.top - seen.T) / step, (seen.top + seen.T + seen.H - t[t.length - 1].y) / step]; };   // in tick steps, top and bottom
    assert.deepEqual([fig.lo, fig.hi], [-0.1, 1.1]); assert.deepEqual(yTicks(fig).map((t) => t.text), ['1.0', '0.8', '0.6', '0.4', '0.2', '0.0']);
    assert.ok(beyond(fig).every((steps) => Math.abs(steps - 0.5) < 0.03), 'half a step beyond the last tick at either end, with every number round: ' + beyond(fig));
    assert.deepEqual([fig.x0, fig.x1], await page.evaluate(() => { const t = window.calciumStep4.S.results.time; return [t[0], t[t.length - 1]]; }), 'the time axis ends where the data end');
    // Marks of time sit in a strip above the plot area, not on the traces: the stimulus label, and the response window as a thin bar.
    const stimulus = fig.texts.find((t) => t.text === 'stimulus'), bars = fig.boxes.filter((b) => b.ink === '#e9a84c');
    assert.ok(stimulus && stimulus.y < fig.top + fig.T, 'the stimulus label is above the plot area');
    assert.ok(bars.some((b) => b.w > 100) && bars.every((b) => b.y + b.h <= fig.top + fig.T - 2), 'the response window is a bar above the plot, clear of its top edge');
    assert.ok(fig.lines.some((l) => l.ink === '#b45f06' && l.dashed && l.width < 1.5), 'the stimulus is a thin dashed line');
    // The cells a plot uses, decided in one place: the cells shown, and the steady ones an automatic range is taken from.
    await page.check('#every');
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.plotCells()), { shown: ['Cell001', 'Cell002', 'Cell005', 'Cell006'], steady: ['Cell001', 'Cell002', 'Cell006'] });
    now = await plotNow(); fig = await cleanBoth('the average with all traces');
    assert.deepEqual(now.legend, ['all cells ± SEM (n = 4)', 'single cells', 'baseline', 'response window'], 'single cells get one entry; no group is listed twice');
    // (their swatch shows both ways they are drawn: a gray line, and under it the wider line of the responders)
    const swatches = fig.lines.filter((l) => inLegend(fig, l.x0, l.y0) && ['#8a94a0', '#177465'].includes(l.ink));
    assert.deepEqual(swatches.map((l) => l.ink), ['#8a94a0', '#177465']); assert.ok(swatches[1].y0 > swatches[0].y0 && swatches[1].width > swatches[0].width && swatches.every((l) => l.alpha === 1));
    // The cell that goes below zero (to −1.5 at the last frame) does not set the range. Its trace ends at the plot's edge, and nothing is written about it.
    assert.deepEqual([now.lo, now.hi], [-0.1, 1.1]); assert.ok(!fig.texts.some((t) => /continue/.test(t.text)) && !('leaving' in await page.evaluate(() => window.calciumStep4.Plot.geom)));
    // Single cells: the others first, neutral gray and faint; then the responders, a little wider and stronger; the mean last, on a white casing.
    const faint = fig.lines.filter((l) => l.alpha < 1), others = faint.filter((l) => l.ink === '#8a94a0'), green = faint.filter((l) => l.ink === '#177465'), order = (l) => fig.lines.indexOf(l);
    assert.deepEqual([others.length, green.length], [2, 2]); assert.ok(Math.max(...others.map(order)) < Math.min(...green.map(order)), 'responders are drawn after the others');
    assert.ok(green[0].width > others[0].width && green[0].alpha > others[0].alpha, 'and a little wider and stronger');
    assert.deepEqual([others[0].alpha, green[0].alpha].map((a) => +a.toFixed(4)), [0.75, 0.9], 'both faded by all four traces drawn (1.5 / √4), a responder half as much again: a lone responder among many is not drawn at full strength');
    const casing = fig.lines.findIndex((l) => l.ink === '#ffffff'), mean = fig.lines.findIndex((l) => l.ink === '#0072b2');
    assert.ok(casing > Math.max(...green.map(order)) && mean === casing + 1 && fig.lines[casing].width > fig.lines[mean].width, 'the mean is drawn last, on a white casing');
    // "Leave out … below zero" reaches the plot: those cells are not drawn, not averaged, and the file beside a figure says so.
    await page.check('#below-zero'); now = await plotNow();
    assert.equal(now.legend[0], 'all cells ± SEM (n = 3)');
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.Plot.drawn), ['Cell001', 'Cell002', 'Cell006']);
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.exportMetadata().cells_left_out), ['Cell005']);
    // and the numbers saved with it say so: a file named mean_trace, with mean_of_all_cells, always holds every cell
    const plotFile = () => page.evaluate(() => { const f = window.calciumStep4; return [f.plotKind(), f.plotNumbers().split('\n')[0]]; });
    assert.deepEqual(await plotFile(), ['mean_trace_without_below_zero', 'Time_s,mean_of_all_cells_without_below_zero,mean_of_all_cells_without_below_zero_SEM,mean_of_all_cells_without_below_zero_n_cells'], 'a mean without the cells below zero is named so, in the file and its columns');
    await page.uncheck('#below-zero'); assert.equal((await plotNow()).legend[0], 'all cells ± SEM (n = 4)');
    assert.deepEqual(await plotFile(), ['mean_trace', 'Time_s,mean_of_all_cells,mean_of_all_cells_SEM,mean_of_all_cells_n_cells'], 'with every cell, the names the release saved');
    // One cell has no SEM: no band, no "± SEM", and the swatch shows the band only where one is drawn.
    await page.click('[data-plot=calls]'); now = await plotNow(); fig = await cleanBoth('responders against non-responders with all traces');
    assert.deepEqual(now.legend, ['responded ± SEM (n = 2)', 'did not respond (n = 1)', 'single cells', 'baseline', 'response window']);
    assert.ok(fig.boxes.some((b) => b.ink === 'rgba(23, 116, 101, 0.25)') && !fig.boxes.some((b) => /^rgba\(200, 102, 14/.test(b.ink)), 'a band behind the swatch of the group of two, none for the group of one');
    await figureStyle({ sem: 'none' }); assert.equal((await drawing(true)).legend[0], 'responded (n = 2)', 'with the SEM not shown the legend does not promise one'); await figureStyle();
    await page.click('[data-plot=mean]');
    // The automatic Y range, on traces made for the purpose: the mean runs from −0.2 to 0.4; of the three steady cells one climbs to
    // 20, one stays at 0.1 and one at −0.5. Without single cells the range is the mean's. With them it takes in their bulk (down
    // to −0.5), but goes no further above the mean than one and a half times the mean's own height (0.4 + 0.9): the mean still
    // fills a quarter of the axis, and the cell that climbs to 20 ends at the plot's edge. A selected cell is shown whole.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, r = f.S.results, kept = { traces: r.traces, mean: r.mean, sem: r.sem, every: f.S.every, cell: f.S.cell, mode: f.S.plotMode }, n = r.time.length, ramp = (a, b) => r.time.map((t, i) => a + (b - a) * i / (n - 1));
      const range = (every, mode, cell) => { f.S.every = every; f.S.plotMode = mode; f.S.cell = cell; window.dispatchEvent(new Event('resize')); return [f.Plot.geom.lo, f.Plot.geom.hi]; };
      r.traces = Object.assign({}, r.traces, { Cell001: ramp(0.1, 0.1), Cell002: ramp(0, 20), Cell006: ramp(-0.5, -0.5) }); r.mean = ramp(-0.2, 0.4); r.sem = ramp(0, 0);
      try { return { mean: range(false, 'mean', null), traces: range(true, 'mean', null), cell: range(false, 'cell', 'Cell002') }; }
      finally { Object.assign(r, { traces: kept.traces, mean: kept.mean, sem: kept.sem }); f.S.every = kept.every; f.S.cell = kept.cell; f.S.plotMode = kept.mode; window.dispatchEvent(new Event('resize')); } }),
      { mean: [-0.25, 0.45], traces: [-0.75, 1.5], cell: [-2.5, 22.5] });
    // Legend rows are even: as few rows as the width allows, with the entries shared out between them (2 + 2, 3 + 2, 3 + 3), so that
    // none is left alone on the last row. Every row starts at the left edge, with one gap between entries. No entry is wider than
    // the room: a single word that would be is shortened, with an ellipsis for what is left out. (Text of 10 pixels a letter; a
    // swatch takes 21 and a gap 14.)
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, ctx = { measureText: (text) => ({ width: 10 * text.length }) }, look = { px: 12 }, entries = (lengths) => lengths.map((n) => ({ kind: 'line', text: 'x'.repeat(n) }));
      const rows = (lengths, width) => { const l = f.legendLayout(ctx, entries(lengths), width, look, false), out = []; for (const place of l.places) { const row = place.y / l.rowHeight; (out[row] = out[row] || []).push(place.x); } return { rows: out, width: l.width, height: l.height / l.rowHeight }; };
      const stacked = f.legendLayout(ctx, entries([8, 30]), 120, look, true);
      return [rows([8, 8, 8, 13], 340), rows([6, 6, 6, 6, 10], 380), rows([6, 6, 6, 6, 6, 6], 400), rows([3, 3, 3, 36], 400), rows([8, 8, 8], 340), rows([40, 2, 2], 100), f.legendLayout(ctx, entries([40, 2, 2]), 100, look, false).places.map((place) => place.text), rows([30, 30, 8, 6, 12], 340),
        stacked.places.map((place) => [place.x, place.y]), stacked.places.map((place) => place.text), stacked.width]; }),
      [{ rows: [[0, 115], [0, 115]], width: 266, height: 2 },                      // three would fit on the first row and leave the fourth alone: two and two
        { rows: [[0, 95, 190], [0, 95]], width: 271, height: 2 },                  // four would fit: three and two
        { rows: [[0, 95, 190], [0, 95, 190]], width: 271, height: 2 },             // four would fit: three and three
        { rows: [[0, 65, 130], [0]], width: 381, height: 2 },                       // the last is too wide to share a row: it stays alone, in two rows still
        { rows: [[0, 115, 230]], width: 331, height: 1 },                           // what fits on one row stays on one row
        { rows: [[0], [0, 55]], width: 96, height: 2 },                             // an entry wider than the room has a row to itself,
        ['xxxxxx…', 'xx', 'xx'],                                                    // shortened to the room (21 + 70 of 100)
        { rows: [[0], [0], [0], [0, 95]], width: 321, height: 4 },                  // two long entries and three short ones: the last row has two, not one
        [[0, 0], [0, 18]], ['xxxxxxxx', 'xxxxxxxx…'], 111]);                        // stacked: one under the other, none wider than the room (120)
    // An entry of several words wider than the room goes on two lines or more, broken between words and never inside "± SEM" or the
    // number of cells, and nothing of it is lost; its row is taller by a line (15) for each line more, and the next row starts under it.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, ctx = { measureText: (text) => ({ width: 10 * text.length }) }, look = { px: 12 }, entries = (texts) => texts.map((text) => ({ kind: 'line', text }));
      const laid = (texts, width, stacked) => { const l = f.legendLayout(ctx, entries(texts), width, look, stacked); return [l.places.map((place) => [place.x, place.y, place.lines]), l.width, l.height]; };
      return [laid(['all cells ± SEM (n = 4)', 'baseline'], 120, true), laid(['all cells ± SEM (n = 4)', 'xx', 'xx'], 200, false), laid(['did not respond ± SEM (n = 116)'], 400, false)]; }),
      [[[[0, 0, ['all cells', '± SEM', '(n = 4)']], [0, 48, ['baseline']]], 111, 66],     // three lines (18 + 15 + 15), then "baseline" under them
        [[[0, 0, ['all cells ± SEM', '(n = 4)']], [0, 33, ['xx']], [55, 33, ['xx']]], 171, 51],
        [[[0, 0, ['did not respond ± SEM (n = 116)']]], 331, 18]]);                    // what fits stays on one line
    // One cell. One threshold line, the required rise, with its entry; the response of a cell that responded is a band with its entry.
    // No name is written at a line or in the strip above the plot.
    await page.click('[data-plot=cell]'); await openCell('Cell001'); await page.uncheck('#every'); now = await plotNow(); fig = await cleanBoth('one cell');
    assert.deepEqual(now.legend, ['Cell001', 'all cells ± SEM (n = 4)', 'required rise', 'baseline', 'response window', 'response']);
    const across = (seen, ink) => seen.lines.filter((l) => l.dashed && l.ink === ink && !inLegend(seen, l.x0, l.y0) && l.y0 === l.y1), valueAt = (seen, y) => seen.hi - (y - seen.top - seen.T) / seen.H * (seen.hi - seen.lo);
    assert.deepEqual([across(fig, '#177465').length, across(fig, '#6a3d9a').length, across(fig, '#a46ea0').length, across(fig, '#7c8794').length], [1, 0, 0, 0], 'one threshold line: the required rise');
    assert.ok(Math.abs(valueAt(fig, across(fig, '#177465')[0].y0) - 0.1) < 0.005 && across(fig, '#177465')[0].x1 - across(fig, '#177465')[0].x0 >= fig.W - 1, 'drawn at its level, across the plot');
    assert.ok(fig.lines.filter((l) => l.dashed && l.ink === '#177465').length === 2, 'the swatch of a dashed line is dashed');   // the line and its swatch
    assert.equal(fig.W + fig.L > 7 * 96 - 20, true, 'the plot runs to the right margin: no room is kept for names beside it');
    await page.check('#every'); now = await plotNow(); await cleanBoth('one cell with all traces');
    assert.deepEqual(now.legend, ['Cell001', 'all cells ± SEM (n = 4)', 'single cells', 'required rise', 'baseline', 'response window', 'response'], 'with single cells drawn too, every entry stays in the legend');
    await page.uncheck('#every');
    // "Thresholds", off until ticked, adds the threshold that did not set the required rise, under its own name, and the baseline mean.
    // The call table says which one set it (Threshold_source). Two levels with one value are one line. The verdict above the plot says
    // what the cell needed and what set that. (The cells of issue 31: noise 12.2 % over a minimum rise of 10 %; noise 7.6 % under it.)
    assert.ok(!(await page.isChecked('#thresholds')) && await page.isEnabled('#thresholds') && await page.isVisible('#thresholds'), 'offered with one cell in ΔF/F₀, and off');
    for (const mode of ['mean', 'calls']) { await page.click('[data-plot=' + mode + ']'); assert.ok(await page.isHidden('#thresholds-label'), 'not shown where no thresholds are drawn (' + mode + '): no greyed box'); }
    await page.click('[data-plot=cell]'); assert.ok(await page.isVisible('#thresholds-label'));
    const withCall = (changes, read) => page.evaluate(async ([changes, read]) => { const f = window.calciumStep4, row = f.S.results.callsByWindow[f.S.window][f.S.cell], kept = Object.assign({}, row); Object.assign(row, changes);
      try { f.pick(f.S.cell); return read === 'verdict' ? document.getElementById('verdict').textContent : null; } finally { if (read) { for (const key of Object.keys(row)) delete row[key]; Object.assign(row, kept); f.pick(f.S.cell); } } }, [changes, read || null]);
    const byNoise = { Baseline_SD: '0.0407', Noise_threshold: '0.1221', Floor_threshold: '0.1', Threshold_increase: '0.1221', Threshold_source: 'Baseline noise' }, byFloor = { Baseline_SD: '0.0253', Noise_threshold: '0.0759', Floor_threshold: '0.1', Threshold_increase: '0.1', Threshold_source: 'Floor' };
    const original = await page.evaluate(() => { const f = window.calciumStep4; return Object.assign({}, f.S.results.callsByWindow[f.S.window][f.S.cell]); });
    assert.match(await withCall(byNoise, 'verdict'), /\(needed 12\.2%, set by its noise\)$/); assert.match(await withCall(byFloor, 'verdict'), /\(needed 10\.0%, the minimum rise\)$/);
    // For a responder the verdict says what the call rests on, from the saved call: the onset, how many frames the rise held and
    // the latency; the largest rise and its time follow. A largest rise outside the rise that held says so: a single higher frame
    // before it (cell 69 of Run2: largest rise at 92.0 s, held from 121.0 s) is "not held", one after its end "after the held rise".
    const cell69 = { Onset_s: '120.987414777279', Run_start_s: '120.987414777279', Run_end_s: '121.987', Run_samples: '2', Response_latency_s: '101.989390969276', Peak_s: '91.9904310703278', Peak_increase: '0.405565745971214', Threshold_increase: '0.303418646217159', Threshold_source: 'Baseline noise' };
    assert.equal(await withCall(cell69, 'verdict'), 'Cell001: Responded · rise held for 2 frames from 121.0 s, 102.0 s after the stimulus · largest rise 40.6% at 92.0 s, not held (needed 30.3%, set by its noise)');
    assert.match(await withCall(Object.assign({}, cell69, { Peak_s: '121.5' }), 'verdict'), /· largest rise 40\.6% at 121\.5 s \(needed/, 'a largest rise within the held rise needs no word');
    assert.match(await withCall(Object.assign({}, cell69, { Peak_s: '150' }), 'verdict'), /· largest rise 40\.6% at 150\.0 s, after the held rise \(needed/);
    assert.match(await withCall({}, 'verdict'), /^Cell001: Responded · rise held from 4\.5 s · largest rise/, 'a call saved without run length or latency says what it has');
    assert.equal(await withCall({ Status: 'Did not respond', Reason: 'No sustained threshold crossing', Peak_increase: '0.05' }, 'verdict'), 'Cell001: Did not respond · no sustained threshold crossing · largest rise 5.0% at 6.0 s (needed 10.0%, the minimum rise)', 'a cell that did not respond keeps the engine\'s reason');
    await page.check('#thresholds');
    for (const [call, entries, lines, view] of [
      [byNoise, ['required rise', 'minimum rise', 'baseline mean'], { '#177465': 0.1221, '#6a3d9a': 0.1, '#7c8794': 0 }, 'thresholds, the noise set the required rise'],
      [byFloor, ['required rise', 'noise threshold', 'baseline mean'], { '#177465': 0.1, '#a46ea0': 0.0759, '#7c8794': 0 }, 'thresholds, the minimum rise set the required rise'],
      [Object.assign({}, byFloor, { Noise_threshold: '0.1' }), ['required rise', 'baseline mean'], { '#177465': 0.1, '#7c8794': 0 }, 'thresholds, both the same'],
      [Object.assign({}, byNoise, { Noise_threshold: '0.12210000001' }), ['required rise', 'minimum rise', 'baseline mean'], { '#177465': 0.1221, '#6a3d9a': 0.1, '#7c8794': 0 }, 'thresholds, the two columns differing in their last digits: the table says which set it, not a comparison'],
      [Object.assign({}, byNoise, { Threshold_source: '' }), ['required rise', 'minimum rise', 'baseline mean'], { '#177465': 0.1221, '#6a3d9a': 0.1, '#7c8794': 0 }, 'thresholds, a table without Threshold_source']]) {
      await withCall(call); now = await plotNow(); fig = await cleanBoth(view);
      assert.deepEqual(now.legend, ['Cell001', 'all cells ± SEM (n = 4)'].concat(entries, ['baseline', 'response window', 'response']), view);
      for (const ink of ['#177465', '#6a3d9a', '#a46ea0', '#7c8794']) { const drawn = across(fig, ink);
        assert.equal(drawn.length, ink in lines ? 1 : 0, view + ': lines in ' + ink); if (drawn.length) assert.ok(Math.abs(valueAt(fig, drawn[0].y0) - lines[ink]) < 0.005, view + ': the line in ' + ink + ' is at its level'); }
    }
    assert.ok(!fig.texts.some((t) => /3 × SD|×/.test(t.text)), 'the noise threshold is not called "3 × SD"');
    assert.equal(await page.evaluate(() => window.calciumStep4.exportMetadata().thresholds_drawn), true, 'the file beside a figure says that the thresholds were drawn');
    // (the fixture's own call: no noise at all, so its noise threshold lies on the baseline mean and the two are one line)
    await page.evaluate((row) => { const f = window.calciumStep4, now = f.S.results.callsByWindow[f.S.window][f.S.cell]; for (const key of Object.keys(now)) delete now[key]; Object.assign(now, row); f.pick(f.S.cell); }, original);
    assert.deepEqual((await plotNow()).legend.slice(2, 4), ['required rise', 'noise threshold']); await cleanBoth('thresholds, no noise');
    // Typed Y limits that leave a threshold outside the plot: the line is not drawn, and the legend does not list it.
    await page.fill('#y-max', '0.05'); now = await plotNow(); await cleanBoth('thresholds, with a Y limit under the required rise');
    assert.ok(!now.legend.includes('required rise') && now.legend.includes('noise threshold')); await page.click('#y-auto');
    await page.uncheck('#thresholds'); assert.equal(await page.evaluate(() => window.calciumStep4.exportMetadata().thresholds_drawn), false);
    // The response band is never too thin to see: a run of one sample is still 3 pixels wide (2.25 pt), with a tick at the top of the
    // plot where it starts. A cell that did not respond has no band and no entry.
    await withCall({ Run_start_s: '6', Run_end_s: '6' }); fig = await cleanBoth('one cell, a run of one sample');
    const runBand = fig.boxes.filter((b) => b.ink === '#d5eedd' && !inLegend(fig, b.x, b.y)), runTick = fig.lines.filter((l) => l.ink === '#177465' && !l.dashed && l.x0 === l.x1 && !inLegend(fig, l.x0, l.y0));
    const xOf = (seen, t) => seen.L + (t - seen.x0) / (seen.x1 - seen.x0) * seen.W;
    assert.ok(runBand.length === 1 && Math.abs(runBand[0].w - 3) < 1e-6 && Math.abs(runBand[0].x - xOf(fig, 6)) < 1e-6 && Math.abs(runBand[0].h - fig.H) < 1e-6, 'the band of a one-sample run is 3 pixels wide, from where the run starts: ' + JSON.stringify(runBand));
    assert.ok(runTick.length === 1 && Math.abs(runTick[0].x0 - xOf(fig, 6)) < 1e-6 && Math.abs(runTick[0].y0 - fig.top - fig.T) < 1e-6 && runTick[0].y1 - runTick[0].y0 === 6, 'with a tick at the top of the plot at its start: ' + JSON.stringify(runTick));
    await page.evaluate((row) => { const f = window.calciumStep4, now = f.S.results.callsByWindow[f.S.window][f.S.cell]; for (const key of Object.keys(now)) delete now[key]; Object.assign(now, row); f.pick(f.S.cell); }, original);
    fig = await drawing(true); assert.ok(fig.boxes.some((b) => b.ink === '#d5eedd' && !inLegend(fig, b.x, b.y) && Math.abs(b.x - xOf(fig, 4.5)) < 1e-6 && Math.abs(b.x + b.w - xOf(fig, 8)) < 1e-6), 'a longer run keeps its own width, 4.5 to 8 s');
    await openCell('Cell002'); now = await plotNow(); fig = await cleanBoth('one cell that did not respond');
    assert.deepEqual(now.legend, ['Cell002', 'all cells ± SEM (n = 4)', 'required rise', 'baseline', 'response window']); assert.ok(!fig.boxes.some((b) => b.ink === '#d5eedd'));
    // A selected cell without a trace: said beside the verdict, not in the plot, and not in a saved figure.
    assert.ok(await page.isHidden('#cell-note')); await openCell('Cell003'); now = await plotNow(); fig = await cleanBoth('one cell without a trace');
    assert.deepEqual(now.legend, ['all cells ± SEM (n = 4)', 'baseline', 'response window']); assert.match(await text('#cell-note'), /^No ΔF\/F₀ trace/); assert.ok(await page.isVisible('#cell-note'));
    assert.ok(await page.evaluate(() => document.getElementById('cell-card').contains(document.getElementById('cell-note')) && document.getElementById('plot').clientHeight === window.calciumStep4.Plot.layer.height), 'the note is part of the cell panel, and the plot is drawn at the size it leaves');
    // The selected cell is always shown whole; negative ticks carry a true minus sign; a dotted line marks zero.
    await openCell('Cell005'); assert.ok(await page.isHidden('#cell-note')); fig = await drawing(true);
    assert.ok(fig.lo <= -1.5);
    assert.ok(yTicks(fig).some((t) => t.text.startsWith('−')) && !fig.texts.some((t) => /^-[\d.]+$/.test(t.text)), 'a true minus sign, not a hyphen');
    assert.ok(fig.lines.some((l) => l.ink === '#9aa5b1' && l.dashed), 'a dotted line at zero');
    // With "Thresholds", a threshold line on zero (the baseline mean of ΔF/F₀) is drawn alone: no dotted zero line mixed into it.
    await page.check('#thresholds'); fig = await drawing(true);
    assert.ok(across(fig, '#177465').length === 1 && fig.lines.some((l) => l.dashed && l.y0 === l.y1 && l.ink !== '#9aa5b1' && !inLegend(fig, l.x0, l.y0) && Math.abs(valueAt(fig, l.y0)) < 1e-3) && !fig.lines.some((l) => l.ink === '#9aa5b1' && !inLegend(fig, l.x0, l.y0)), 'with a threshold drawn on zero, no dotted zero line under it');
    await page.uncheck('#thresholds');
    // The thresholds are not in the stimulus line's brown.
    assert.ok(await page.evaluate(() => { const s = window.calciumStep4.STYLE.plot, rgb = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16)), apart = (a, b) => Math.hypot(...rgb(a).map((v, i) => v - rgb(b)[i]));
      return [s.colourLevelFloor, s.colourLevelNoise, s.colourLevelBaseline].every((c) => apart(c, s.colourStimulus) > 120); }), 'the minimum rise and the other thresholds are told apart from the stimulus line by colour');
    assert.ok(await page.evaluate(() => { const s = window.calciumStep4.STYLE.plot, rgb = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16)), apart = (a, b) => Math.hypot(...rgb(a).map((v, i) => v - rgb(b)[i]));
      return apart(s.colourLevelFloor, s.colourLevelBaseline) > 70 && apart(s.colourLevelFloor, s.colourMean) > 70; }), 'and the minimum rise, drawn with the baseline mean under Thresholds, from it and from the mean');
    // Placement in a figure: rows from the Y axis; too narrow, it wraps by whole entries; "Right of the plot" stacks the entries; hidden shows none.
    assert.equal(fig.legendAt[0], fig.L);
    await figureStyle({ width: 3 }); fig = clean(await drawing(true), 'a narrow figure');
    const rows = [...new Set(fig.places.map((place) => place[1]))];
    assert.ok(rows.length > 1 && rows.every((y) => Math.min(...fig.places.filter((place) => place[1] === y).map((place) => place[0])) === 0) && fig.legendAt[0] + fig.legendWidth <= 3 * 96, 'a narrow figure wraps the legend by whole entries: ' + JSON.stringify(fig.places));
    await figureStyle({ legend: 'right' }); fig = clean(await drawing(true), 'the legend right of the plot');
    assert.ok(fig.legendAt[0] >= fig.L + fig.W && fig.legendAt[0] + fig.legendWidth <= 7 * 96 && fig.places.every((place, i) => place[0] === 0 && (i === 0 || place[1] > fig.places[i - 1][1])), 'right of the plot, stacked, inside the figure');
    await figureStyle({ legend: 'none' }); fig = await drawing(true);
    assert.ok(fig.legendAt === null && !fig.texts.some((t) => t.text === 'response window'), 'hidden: no legend'); assert.deepEqual(wordsOf(fig), ['Time (s)', 'stimulus', 'ΔF/F₀'], 'and no name is written anywhere else instead');
    // What the figure editor can leave out goes from the legend with it: the response band, the required rise.
    await figureStyle({ run: false, levels: false }); fig = clean(await drawing(true), 'a figure without the response band and the required rise');
    assert.deepEqual(fig.legend, ['Cell005', 'all cells ± SEM (n = 4)', 'baseline', 'response window']);
    await figureStyle({ title: 'Typed title' }); assert.ok((await drawing(true)).texts.some((t) => t.text === 'Typed title'), 'a typed title is drawn');
    await figureStyle(); await openCell('Cell002');

    // Dragging begins on another cell's trace; pointer capture keeps it scrubbing.
    await page.locator('#plot').scrollIntoViewIfNeeded();
    await page.click('[data-plot=cell]'); await page.check('#every');   // the page opens on the average alone
    const points = await page.evaluate(() => {
      const { S, Plot } = window.calciumStep4, g = Plot.geom, box = document.getElementById('plot').getBoundingClientRect();
      const x = (t) => box.left + g.L + (t - g.x0) / (g.x1 - g.x0) * g.W;
      return { start: { x: x(S.results.time[12]), y: box.top + g.T + (g.hi - S.results.traces.Cell001[12]) / (g.hi - g.lo) * g.H },
        end: { x: x(S.results.time[21]), y: box.top + g.T + g.H / 2 } };
    });
    // The numbers on the plot follow the pointer; nothing is drawn on the traces while nobody is reading it.
    const plotAtRest = () => page.waitForFunction(() => { const P = window.calciumStep4.Plot; return !P.line && P.readout === null; }, null, { timeout: 6000 });
    const readout = () => page.evaluate(() => window.calciumStep4.Plot.readout), secondsAt = (sample) => page.evaluate((i) => window.calciumStep4.S.results.time[i].toFixed(1) + ' s', sample);
    const from = await secondsAt(12), to = await secondsAt(21);
    await page.mouse.move(5, 5); await plotAtRest();   // no line and no box while nothing moves and the pointer is elsewhere
    const woke = await page.evaluate(() => { const f = window.calciumStep4, was = f.S.frame; f.setFrame(was === 10 ? 11 : 10, true); const seen = { line: f.Plot.line, readout: f.Plot.readout, frame: f.S.frame }; f.setFrame(was, true); return seen; });
    assert.ok(woke.line && woke.readout.startsWith(await secondsAt(woke.frame - 1)), 'moving the frame from anywhere shows where the picture is, with that frame\'s numbers');
    await plotAtRest();   // and a moment after it stops, they go
    await page.mouse.move(points.start.x, points.start.y);
    assert.ok((await readout()).startsWith(from + ' · '), 'over the plot the numbers are those under the pointer'); assert.equal(await page.evaluate(() => window.calciumStep4.Plot.line), false, 'the pointer brings its own dashed line, not the frame\'s'); assert.equal(await page.evaluate(() => window.calciumStep4.Plot.at), 12);
    await page.mouse.move(points.end.x, points.end.y); assert.ok((await readout()).startsWith(to + ' · '), 'and move with it');
    assert.notEqual(await page.evaluate(() => window.calciumStep4.S.frame), 22, 'looking does not move the recording');
    await page.mouse.move(5, 5); assert.equal(await readout(), null, 'and go when the pointer leaves');
    await page.mouse.move(points.start.x, points.start.y); await page.mouse.down();
    assert.equal(await page.evaluate(() => window.calciumStep4.S.cell), 'Cell001');
    assert.ok(await page.evaluate(() => window.calciumStep4.Plot.scrubbing));
    assert.ok((await readout()).startsWith(from + ' · Cell001 '), 'while dragging the numbers are the frame\'s'); assert.ok(await page.evaluate(() => window.calciumStep4.Plot.line), 'and the solid line is at the frame');
    await page.mouse.move(points.end.x, points.end.y, { steps: 6 }); assert.ok((await readout()).startsWith(to + ' · Cell001 '));
    await page.mouse.up();
    await waitFrame(22);
    assert.equal(await page.evaluate(() => window.calciumStep4.Plot.scrubbing), false);
    assert.ok((await readout()).startsWith(to + ' · Cell001 '), 'after the drag the numbers stay under the pointer');
    await page.mouse.move(5, 5); await plotAtRest();   // a click elsewhere leaves nothing behind on the plot
    await page.mouse.move(points.end.x, points.end.y);
    assert.match(await text('#time-label'), /15\.50 s.*Frame 22 \/ 30/);

    // No cell selected is a state of its own (the page opens that way), and the plot follows the selection.
    const app = (what) => page.evaluate((key) => { const a = window.calciumStep4; return key === 'cell' ? a.S.cell : key === 'plot' ? a.S.plotMode : a.Pic.zoom; }, what);
    await openCell('');
    assert.equal(await app('cell'), null); assert.equal(await text('#verdict'), 'No cell selected.'); assert.ok(await page.isDisabled('#zoom-cell'));
    assert.equal(await app('plot'), 'mean', 'with no cell the plot shows the average');
    await openCell('Cell002');
    assert.equal(await app('plot'), 'cell', 'choosing a cell shows that cell over the average');
    // The Cell box: a number alone and Enter opens that cell; text that names no cell is put back to the open one.
    await openCell('1'); assert.equal(await app('cell'), 'Cell001'); assert.equal(await page.inputValue('#cell'), 'Cell001');
    await openCell('no such cell'); assert.equal(await app('cell'), 'Cell001'); assert.equal(await page.inputValue('#cell'), 'Cell001');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4.cellTyped; return [f('cell002'), f(' 2 ', true), f('2', false) === undefined, f('', true), f('77', true) === undefined]; }), ['Cell002', 'Cell002', true, null, true], 'a name in any case, a number on Enter only, nothing for an empty box, no match otherwise');
    await openCell('Cell002');
    await page.evaluate(() => document.activeElement && document.activeElement.blur()); await page.keyboard.press('Escape');
    assert.equal(await app('cell'), null, 'Escape leaves no cell selected'); assert.equal(await app('plot'), 'mean');
    await page.keyboard.press('ArrowRight'); assert.equal(await app('cell'), 'Cell001', 'the arrows start from the first cell');
    await page.click('[data-plot=calls]'); await openCell('Cell002'); assert.equal(await app('plot'), 'calls', 'a plot type chosen by hand stays');
    await page.click('[data-plot=cell]'); await openCell('Cell001');
    // The picture: a click on a cell opens it, a click on empty space shows no cell, a drag (pan) is not a click, and a double-click
    // on empty space shows the whole picture and keeps the cell. The Cell box: its ✕ empties it; Esc in it only undoes the typing.
    // (a spot well inside a cell other than the open one, and a spot with no outline within 8 pixels; neither under the controls that lie on the picture's corner)
    const spots = await page.evaluate(() => { const a = window.calciumStep4, canvas = document.getElementById('picture'), b = canvas.getBoundingClientRect(), bare = (x, y) => document.elementFromPoint(x, y) === canvas, at = (x, y) => (bare(x, y) ? a.cellAt({ clientX: x, clientY: y }) : NaN);
      const inside = (x, y, i) => [[0, 0], [5, 0], [-5, 0], [0, 5], [0, -5]].every(([dx, dy]) => at(x + dx, y + dy) === i);
      const clear = (x, y) => [-8, 0, 8].every((dx) => [-8, 0, 8].every((dy) => at(x + dx, y + dy) < 0));
      let cell = null, empty = null;
      for (let y = Math.round(b.top) + 12; y < b.bottom - 12 && !(cell && empty); y += 3) for (let x = Math.round(b.left) + 12; x < b.right - 12 && !(cell && empty); x += 3) {
        const i = at(x, y); if (i >= 0 && !cell && a.S.info.names[i] !== a.S.cell && inside(x, y, i)) cell = { x, y, name: a.S.info.names[i] }; else if (i < 0 && !empty && clear(x, y)) empty = { x, y }; }
      return { cell, empty }; });
    assert.ok(spots.cell && spots.empty, 'the picture has a cell and empty space to click');
    await page.mouse.click(spots.cell.x, spots.cell.y); assert.equal(await app('cell'), spots.cell.name, 'a click on a cell opens it');
    await page.mouse.click(spots.empty.x, spots.empty.y);
    assert.equal(await app('cell'), null, 'a click on empty space shows no cell'); assert.equal(await page.inputValue('#cell'), ''); assert.equal(await app('plot'), 'mean');
    await page.mouse.click(spots.cell.x, spots.cell.y);
    await page.mouse.move(spots.empty.x, spots.empty.y); await page.mouse.down(); await page.mouse.move(spots.empty.x + 30, spots.empty.y + 20, { steps: 5 }); await page.mouse.up();
    assert.equal(await app('cell'), spots.cell.name, 'dragging the picture keeps the cell');
    assert.ok(await app('zoom') === 1 && await page.isDisabled('#zoom-all') && await page.evaluate(([x, y]) => window.calciumStep4.cellAt({ clientX: x, clientY: y }) < 0, [spots.empty.x, spots.empty.y]), 'and at 1× there is nothing to move: the whole picture stays where it is');
    await page.mouse.move(spots.empty.x, spots.empty.y); await page.keyboard.down('Control'); await page.mouse.wheel(0, -300); await page.keyboard.up('Control'); await settledPaint();
    assert.ok(await app('zoom') > 1 && await page.evaluate(([x, y]) => window.calciumStep4.cellAt({ clientX: x, clientY: y }) < 0, [spots.empty.x, spots.empty.y]), 'zoomed in, still on empty space');
    await page.mouse.dblclick(spots.empty.x, spots.empty.y); await settledPaint();
    assert.equal(await app('zoom'), 1, 'a double-click shows the whole picture'); assert.equal(await app('cell'), spots.cell.name, 'and keeps the cell'); assert.equal(await app('plot'), 'cell', 'and the plot');
    assert.equal(await page.inputValue('#cell'), spots.cell.name, 'the box shows the cell');
    const boxLook = await page.evaluate(() => { const look = (e) => { const c = getComputedStyle(e); return [c.height, c.borderTopColor, c.borderRadius, c.fontSize, c.backgroundColor].join(' '); }; return [look(document.getElementById('cell')), look(document.getElementById('y-min'))]; });   // (the frames boxes, once the page's number boxes, are text boxes now, as the other boxes of the page)
    assert.equal(boxLook[0], boxLook[1], 'the Cell box looks like the other boxes');
    // (the ✕ shows on hover, left of the list's arrow; found by clicking from the text towards the arrow)
    const cellBox = await page.locator('#cell').boundingBox(); let crossed = false;
    for (let dx = 46; dx >= 14 && !crossed; dx -= 2) { await page.mouse.click(cellBox.x + cellBox.width - dx, cellBox.y + cellBox.height / 2); crossed = (await page.inputValue('#cell')) === ''; }
    assert.equal(await app('cell'), null, 'the box’s ✕ shows no cell'); assert.equal(await page.inputValue('#cell'), '');
    await openCell(spots.cell.name); await page.click('#cell'); await page.keyboard.type('Cel'); await page.keyboard.press('Escape'); await settledPaint();
    assert.equal(await app('cell'), spots.cell.name, 'Esc in the box keeps the cell'); assert.equal(await page.inputValue('#cell'), spots.cell.name, 'and puts its name back');
    assert.ok(await page.evaluate(() => document.activeElement !== document.getElementById('cell')), 'and leaves the box');
    await openCell('Cell001');
    // GFP: the recording in green, the green of "Response in magenta" without the magenta.
    await page.selectOption('#view', 'green'); await waitFrame(22);
    assert.match(await text('#view-key'), /cells, background removed/);
    const greenFrame = await page.evaluate(async () => {
      const bitmap = await createImageBitmap(await (await fetch('api/frame?f=15&mode=green')).blob()), canvas = document.createElement('canvas'); canvas.width = bitmap.width; canvas.height = bitmap.height;
      const ctx = canvas.getContext('2d'); ctx.drawImage(bitmap, 0, 0); const d = ctx.getImageData(0, 0, bitmap.width, bitmap.height).data, sum = [0, 0, 0];
      for (let i = 0; i < d.length; i += 4) { sum[0] += d[i]; sum[1] += d[i + 1]; sum[2] += d[i + 2]; } return sum; });
    assert.ok(greenFrame[1] > 0 && greenFrame[0] < 0.15 * greenFrame[1] && greenFrame[2] < 0.15 * greenFrame[1], 'a GFP frame is green only: ' + greenFrame);
    for (const mode of ['magenta', 'change']) assert.equal(await page.evaluate(async (m) => (await fetch('api/frame?f=15&mode=' + m)).status, mode), 200);
    await page.selectOption('#view', 'raw'); await waitFrame(22);
    // Outline thickness sits with "All outlines", not with the file options.
    assert.ok(await page.evaluate(() => document.querySelector('#picture-card .toolbar').contains(document.getElementById('outline-width')) && !document.getElementById('save-options').contains(document.getElementById('outline-width'))));
    // Scrolling over the picture scrolls the page; Ctrl (or Cmd, or a pinch) + scroll zooms.
    const pictureBox = await page.locator('#picture').boundingBox(), zoomBefore = await app('zoom');
    await page.mouse.move(pictureBox.x + pictureBox.width / 2, pictureBox.y + pictureBox.height / 2); await page.mouse.wheel(0, -300); await settledPaint();
    assert.equal(await app('zoom'), zoomBefore, 'scrolling over the picture does not zoom it');
    await page.keyboard.down('Control'); await page.mouse.wheel(0, -300); await page.keyboard.up('Control'); await settledPaint();
    assert.ok(await app('zoom') > zoomBefore, 'Ctrl + scroll zooms'); await page.click('#zoom-all'); await page.evaluate(() => window.scrollTo(0, 0));

    // ---------------------------------------------------------------- the picture: zooming and moving, zoom to region, the small map
    // One view (a zoom, and the point of the recording in the middle) and one way to set it. What these checks read of it: the zoom,
    // where the picture lies in its panel (x, y, and z screen pixels to one of its own), the middle and the size of the room it
    // is fitted into, the panel's place in the window, and what the controls and the cursor say.
    const viewNow = () => page.evaluate(() => { const a = window.calciumStep4, c = document.getElementById('picture'), b = c.getBoundingClientRect();
      return Object.assign(a.place(c.clientWidth, c.clientHeight), { zoom: a.Pic.zoom, w: c.clientWidth, h: c.clientHeight, left: b.left, top: b.top, iw: a.S.info.width, ih: a.S.info.height, level: document.getElementById('zoom-level').textContent,
        cell: a.S.cell, cursor: getComputedStyle(c).cursor, region: document.getElementById('zoom-region').getAttribute('aria-pressed'), band: a.Pic.band }); });
    const toScreen = (v, x, y) => [v.left + v.x + x * v.z, v.top + v.y + y * v.z], toRecording = (v, x, y) => [(x - v.left - v.x) / v.z, (y - v.top - v.y) / v.z];
    const near = (a, b, within) => a.length === b.length && a.every((value, i) => Math.abs(value - b[i]) <= (within === undefined ? 0.5 : within));
    const drag = async (from, to, during) => { await page.mouse.move(from[0], from[1]); await page.mouse.down(); await page.mouse.move(to[0], to[1], { steps: 4 }); const seen = during ? await during() : null; await page.mouse.up(); return seen; };
    let seen = await viewNow(); const plotAt = await page.locator('#plot').boundingBox();
    // The controls: one row on the picture's top left corner, small buttons of the page's own height and corner, each with a tooltip and a
    // name: zoom out, the zoom level, zoom in, then Whole picture and Zoom to region in words, as Zoom to cell is. The picture lies under
    // the row, never beneath it.
    const toolsNow = () => page.evaluate(() => { const stage = document.getElementById('picture-stage').getBoundingClientRect(), all = document.getElementById('view-tools').getBoundingClientRect(), parts = Array.from(document.querySelectorAll('#view-tools button, #view-tools .zoom-level'));
      return { box: [all.left, all.top, all.right, all.bottom], corner: all.left - stage.left <= 8 && all.top - stage.top <= 8 && document.getElementById('picture-stage').contains(document.getElementById('view-tools')), row: new Set(parts.map((e) => Math.round(e.getBoundingClientRect().top))).size === 1,
        order: parts.map((e) => Math.round(e.getBoundingClientRect().left)), level: getComputedStyle(document.getElementById('zoom-level')).backgroundColor,
        parts: parts.map((e) => { const b = e.getBoundingClientRect(); return [e.id, (b.width === 28 ? 'square' : b.width > 60 ? 'words' : b.width) + ' ' + b.height, e.title ? 'tooltip' : '', e.getAttribute('aria-label') || e.textContent].join(' · '); }) }; });
    const tools = await toolsNow();
    assert.ok(tools.corner && tools.row && tools.order.every((x, i) => !i || x > tools.order[i - 1]), 'the zoom and move controls lie in one row on the picture\'s top left corner: ' + JSON.stringify(tools));
    assert.deepEqual(tools.parts, ['zoom-out · square 28 · tooltip · Zoom out', 'zoom-level · square 28 · tooltip · 1×', 'zoom-in · square 28 · tooltip · Zoom in', 'zoom-all · words 28 · tooltip · Whole picture', 'zoom-region · words 28 · tooltip · Zoom to region']);
    { const corner = toScreen(seen, 0, 0); assert.ok(corner[1] >= tools.box[3] + 2, 'the whole picture lies under the controls, not beneath them: ' + corner + ' against ' + tools.box); }
    assert.deepEqual([seen.zoom, seen.level, seen.cursor, await page.isDisabled('#zoom-out'), await page.isDisabled('#zoom-all'), await page.isEnabled('#zoom-in')], [1, '1×', 'crosshair', true, true, true], 'at 1× there is nothing to zoom out of');
    // The buttons zoom in steps about the middle of the view, from 1× to 16×; Whole picture goes back.
    await page.click('#zoom-in'); seen = await viewNow();
    assert.deepEqual([seen.zoom, seen.level, seen.cursor], [1.5, '1.5×', 'grab']); assert.ok(near(toScreen(seen, seen.iw / 2, seen.ih / 2), [seen.left + seen.midX, seen.top + seen.midY]), 'a button zooms about the middle of the view');
    assert.ok(await page.isEnabled('#zoom-out') && await page.isEnabled('#zoom-all'));
    await page.click('#zoom-out'); assert.equal(await app('zoom'), 1);
    for (let i = 0; i < 8; i++) await page.click('#zoom-in');
    seen = await viewNow(); assert.deepEqual([seen.zoom, seen.level, await page.isDisabled('#zoom-in')], [16, '16×', true], 'eight steps reach 16×, the most');
    await page.click('#zoom-all'); assert.equal(await app('zoom'), 1, 'Whole picture shows all of it');
    // The keys, while the pointer is over the picture: + and − zoom about the pointer, 0 shows the whole picture. Elsewhere, or while typing, they do nothing.
    seen = await viewNow(); const pointer = toScreen(seen, 0.4 * seen.iw, seen.ih / 2).map(Math.round), beneath = toRecording(seen, pointer[0], pointer[1]);
    await page.mouse.move(pointer[0], pointer[1]); await page.keyboard.press('+'); seen = await viewNow();
    assert.equal(seen.zoom, 1.5, '+ zooms in'); assert.ok(near(toScreen(seen, beneath[0], beneath[1]), pointer), 'about the point under the pointer: ' + toScreen(seen, beneath[0], beneath[1]) + ' for ' + pointer);
    await page.keyboard.press('='); assert.equal(await app('zoom'), 2); await page.keyboard.press('-'); assert.equal(await app('zoom'), 1.5, '− zooms out'); await page.keyboard.press('0'); assert.equal(await app('zoom'), 1, '0 shows the whole picture');
    await page.focus('#cell'); await page.keyboard.press('+'); assert.equal(await app('zoom'), 1, 'not while typing in a box'); await page.keyboard.press('Escape');
    await page.mouse.move(plotAt.x + plotAt.width / 2, plotAt.y + plotAt.height / 2); await page.keyboard.press('+'); assert.equal(await app('zoom'), 1, 'not with the pointer elsewhere');
    // Ctrl + scroll zooms about the pointer. A pinch on a trackpad arrives as small wheel steps marked Ctrl with no Ctrl key held, which
    // give the change of scale as exp(−deltaY / 100): fingers spread to twice their distance double the picture, and the page itself is
    // not zoomed. The same steps with Ctrl really held zoom at the wheel's pace. Safari sends the scale itself.
    const wheelSteps = (deltaY, times) => page.evaluate(([x, y, d, n]) => { const c = document.getElementById('picture'); let taken = true; for (let i = 0; i < n; i++) taken = !c.dispatchEvent(new WheelEvent('wheel', { deltaY: d, ctrlKey: true, clientX: x, clientY: y, bubbles: true, cancelable: true })) && taken; return taken; }, [pointer[0], pointer[1], deltaY, times]);
    await page.mouse.move(pointer[0], pointer[1]);
    assert.ok(await wheelSteps(-100 * Math.log(2) / 4, 4), 'a pinch is taken by the picture'); seen = await viewNow();
    assert.ok(Math.abs(seen.zoom - 2) < 1e-9, 'a pinch to twice the distance doubles the picture: ' + seen.zoom); assert.ok(near(toScreen(seen, beneath[0], beneath[1]), pointer), 'about the point under the fingers');
    await page.keyboard.down('Control'); await wheelSteps(-100 * Math.log(2) / 4, 4); await page.keyboard.up('Control');
    assert.ok(Math.abs(await app('zoom') - 2 * Math.pow(1.0015, 100 * Math.log(2))) < 1e-9, 'with Ctrl held the same steps are a wheel\'s');
    await page.keyboard.down('Control'); await page.mouse.wheel(0, -200); await page.keyboard.up('Control'); await settledPaint(); seen = await viewNow();
    assert.ok(Math.abs(seen.zoom - 2 * Math.pow(1.0015, 100 * Math.log(2) + 200)) < 1e-6 && near(toScreen(seen, beneath[0], beneath[1]), pointer), 'Ctrl + scroll zooms about the pointer');
    const pinched = await page.evaluate(([x, y]) => { const c = document.getElementById('picture'), a = window.calciumStep4, was = a.Pic.zoom, send = (type, scale) => c.dispatchEvent(Object.assign(new Event(type, { bubbles: true, cancelable: true }), { scale, clientX: x, clientY: y }));
      send('gesturestart', 1); send('gesturechange', 1.2); send('gesturechange', 1.5); send('gestureend', 1.5); return a.Pic.zoom / was; }, pointer);
    assert.ok(Math.abs(pinched - 1.5) < 1e-9, 'a pinch sent as a gesture zooms by its scale: ' + pinched); assert.ok(near(toScreen(await viewNow(), beneath[0], beneath[1]), pointer));
    await page.click('#zoom-all'); await page.evaluate(() => window.scrollTo(0, 0));
    // Zoomed in, a drag moves the picture under a hand (closed while it moves), with a cell selected too, and keeps the cell. The
    // picture cannot be dragged out of view: no edge of it passes the middle. Against that limit a drag back takes hold at once.
    assert.equal(await app('cell'), 'Cell001'); await page.click('#zoom-in'); await page.click('#zoom-in'); seen = await viewNow();
    const grip = [Math.round(seen.left + seen.midX), Math.round(seen.top + seen.midY)], moving = await drag(grip, [grip[0] - 40, grip[1] - 25], viewNow), moved = await viewNow();
    assert.ok(near([moved.x - seen.x, moved.y - seen.y], [-40, -25], 1e-6), 'a drag moves the picture with the pointer: ' + [moved.x - seen.x, moved.y - seen.y]);
    assert.deepEqual([seen.cursor, moving.cursor, moved.cursor], ['grab', 'grabbing', 'grab'], 'a hand, closed while it drags'); assert.deepEqual([moved.cell, await app('plot'), moved.zoom], ['Cell001', 'cell', 2], 'the cell stays selected');
    await drag(grip, [grip[0] + 3000, grip[1] + 3000]); let far = await viewNow();
    assert.ok(near([far.x, far.y], [far.midX, far.midY], 1e-6), 'dragged as far as it goes, the picture\'s corner is in the middle of the view: ' + [far.x, far.y]);
    await drag(grip, [grip[0] - 3000, grip[1] - 3000]); far = await viewNow();
    assert.ok(near([far.x + far.iw * far.z, far.y + far.ih * far.z], [far.midX, far.midY], 1e-6), 'and the opposite corner the other way');
    await drag(grip, [grip[0] + 12, grip[1] + 9]); assert.ok(near([(await viewNow()).x - far.x, (await viewNow()).y - far.y], [12, 9], 1e-6), 'a drag back from the limit moves the picture at once');
    assert.equal(await app('cell'), 'Cell001'); await page.click('#zoom-all');
    // Zoom to region: with the button on, a drag draws a rectangle and the view zooms to it; the part of the recording under the
    // rectangle then fills the room beside the controls, to the pixel, and the button turns itself off.
    // (a rectangle on the whole picture as it lies now: two corners on whole screen pixels, and the same rectangle in the recording's pixels)
    const rectangle = async () => { const view = await viewNow(), a = toScreen(view, 0.3 * view.iw, 0.25 * view.ih).map(Math.round), b = toScreen(view, 0.55 * view.iw, 0.6 * view.ih).map(Math.round); return { a, b, part: toRecording(view, a[0], a[1]).concat(toRecording(view, b[0], b[1])) }; };
    const landed = (view, part) => { const z = Math.min((view.wide - 8) / (part[2] - part[0]), (view.high - 8) / (part[3] - part[1])), wide = (part[2] - part[0]) * z, high = (part[3] - part[1]) * z;
      return Math.abs(view.z - z) < 1e-9 && near(toScreen(view, part[0], part[1]).concat(toScreen(view, part[2], part[3])), [view.left + view.midX - wide / 2, view.top + view.midY - high / 2, view.left + view.midX + wide / 2, view.top + view.midY + high / 2], 1); };
    let drawn = await rectangle(); const cornerA = drawn.a, cornerB = drawn.b, part = drawn.part;
    await page.click('#zoom-region'); seen = await viewNow(); assert.deepEqual([seen.region, seen.cursor], ['true', 'crosshair'], 'the button stays down and the pointer is a cross');
    const whileDrawn = await drag(cornerA, cornerB, viewNow); seen = await viewNow();
    assert.ok(near(whileDrawn.band, [cornerA[0] - whileDrawn.left, cornerA[1] - whileDrawn.top, cornerB[0] - whileDrawn.left, cornerB[1] - whileDrawn.top], 1e-6) && whileDrawn.zoom === 1, 'while it is drawn the rectangle follows the pointer and nothing moves');
    assert.ok(landed(seen, part) && seen.zoom > 2, 'the view lands on the rectangle that was drawn: ' + JSON.stringify([seen.zoom, toScreen(seen, part[0], part[1]), toScreen(seen, part[2], part[3])]));
    assert.deepEqual([seen.region, seen.band, seen.cursor, seen.cell], ['false', null, 'grab', 'Cell001'], 'the button turns itself off, and the cell is kept'); const byButton = seen.zoom;
    // Esc cancels, before a rectangle and in the middle of one, and keeps the cell (Esc otherwise shows no cell).
    await page.click('#zoom-all'); await page.click('#zoom-region'); await page.keyboard.press('Escape'); seen = await viewNow(); assert.deepEqual([seen.region, seen.cell, seen.zoom], ['false', 'Cell001', 1], 'Esc leaves Zoom to region');
    await page.click('#zoom-region'); await drag(cornerA, cornerB, () => page.keyboard.press('Escape')); seen = await viewNow(); assert.deepEqual([seen.region, seen.band, seen.cell, seen.zoom], ['false', null, 'Cell001', 1], 'also while a rectangle is being drawn');
    // A click on a cell still selects it while the button is on, and the button stays on; a rectangle too small to be meant does nothing.
    const other = await page.evaluate(() => { const a = window.calciumStep4, i = a.S.info.names.indexOf('Cell002'), b = a.S.bounds[i]; return [(b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2]; });
    await page.click('#zoom-region'); seen = await viewNow(); await page.mouse.click(...toScreen(seen, other[0], other[1]).map(Math.round)); seen = await viewNow();
    assert.deepEqual([seen.cell, seen.region, seen.zoom], ['Cell002', 'true', 1], 'a click on a cell selects it');
    await drag(cornerA, [cornerA[0] + 5, cornerA[1] + 2]); seen = await viewNow(); assert.deepEqual([seen.zoom, seen.region, seen.cell], [1, 'true', 'Cell002'], 'a rectangle of a few pixels is not a region');
    await page.keyboard.press('Escape'); await openCell('Cell001');
    // Shift + drag does the same without the button; in every view, and while the movie plays.
    await page.keyboard.down('Shift'); await drag(cornerA, cornerB); await page.keyboard.up('Shift'); seen = await viewNow();
    assert.ok(landed(seen, part) && seen.zoom === byButton && seen.region === 'false', 'Shift + drag zooms to the rectangle'); await page.click('#zoom-all');
    for (const mode of ['flat', 'magenta', 'change', 'picture']) { await page.selectOption('#view', mode); if (mode !== 'picture') await waitFrame(22); await settledPaint(); drawn = await rectangle();   // (the key under the picture changes with the view, and with it the panel's height)
      await page.keyboard.down('Shift'); await drag(drawn.a, drawn.b); await page.keyboard.up('Shift'); seen = await viewNow(); assert.ok(landed(seen, drawn.part) && seen.zoom > 2 && await page.isVisible('#view-map'), 'zoom to region in the view ' + mode);
      await page.keyboard.press('0'); assert.equal(await app('zoom'), 1); }
    await page.selectOption('#view', 'raw'); await waitFrame(22); await settledPaint();
    await slider(1); await page.click('#play'); await page.waitForFunction(() => window.calciumStep4.S.frame >= 2);
    await page.click('#zoom-in'); drawn = await rectangle(); await page.click('#zoom-region'); await drag(drawn.a, drawn.b); seen = await viewNow(); const playing = await drag(grip, [grip[0] - 20, grip[1] - 10], viewNow);
    assert.ok(seen.zoom > 2 && seen.region === 'false' && playing.cursor === 'grabbing' && near([(await viewNow()).x - seen.x], [-20], 1e-6) && await page.evaluate(() => !!window.calciumStep4.S.timer), 'the buttons, a rectangle and a drag work while the movie plays, and it plays on');
    await page.click('#play'); await page.click('#zoom-all'); await slider(22);
    // Zoomed in, the frame on screen is wanted at its full size at once, not only when the frame next changes. (The fixture's frames
    // are small enough to be fetched whole at any zoom: for this check the page is told the recording is 16 times as wide.)
    const widthsAsked = [], noteWidth = (request) => { const asked = /api\/frame\?f=22&mode=raw&w=(\d+)&/.exec(request.url()); if (asked) widthsAsked.push(Number(asked[1])); }; page.on('request', noteWidth);
    await page.evaluate(() => { const a = window.calciumStep4; a.S.info.width *= 16; a.S.info.height *= 16; a.resetFrames(); }); await waitFrame(22);
    assert.ok(widthsAsked.length > 0 && widthsAsked.every((width) => width < 1000), 'at 1× a large recording is fetched reduced: ' + widthsAsked);
    await page.click('#zoom-in'); await page.waitForFunction(() => !window.calciumStep4.Pic.loading); page.off('request', noteWidth);
    assert.equal(widthsAsked[widthsAsked.length - 1], 4000, 'zooming in fetches the frame on screen at full size: ' + widthsAsked);
    await page.evaluate(() => { const a = window.calciumStep4; a.S.info.width /= 16; a.S.info.height /= 16; a.resetFrames(); }); await page.click('#zoom-all'); await waitFrame(22);
    // The small map: there while the picture is zoomed in, on the bottom right corner, in the picture's shape. Its frame is the part of
    // the recording that is in the panel. A press in it puts that point in the middle of the view, a drag takes the view along, and the cell stays.
    assert.ok(await page.isHidden('#view-map'), 'no map at 1×'); for (let i = 0; i < 3; i++) await page.click('#zoom-in');
    const mapNow = () => page.evaluate(() => { const a = window.calciumStep4, m = document.getElementById('view-map'), b = m.getBoundingClientRect(), stage = document.getElementById('picture-stage').getBoundingClientRect(), f = a.Pic.map, ctx = m.getContext('2d'), dpr = m.width / m.clientWidth;
      const light = (x, y) => { const d = ctx.getImageData(Math.floor(x * dpr), Math.floor(y * dpr), 1, 1).data; return d[0] + d[1] + d[2]; };
      return { box: [b.left, b.top, b.width, b.height], fromCorner: [stage.right - b.right, stage.bottom - b.bottom], frame: [f.x, f.y, f.width, f.height], edge: light(f.x + f.width / 2, f.y + 0.75), around: light(f.x + f.width / 2, Math.max(1, f.y / 2)), area: a.viewArea() }; });
    seen = await viewNow(); let map = await mapNow(); const perPixel = [map.box[2] / seen.iw, map.box[3] / seen.ih];   // the map's pixels to one of the recording, across and down (each side is a whole number of pixels)
    assert.ok(await page.isVisible('#view-map') && map.fromCorner.every((gap) => gap >= 4 && gap <= 12) && Math.abs(map.box[2] / map.box[3] - seen.iw / seen.ih) < 0.03, 'the map is on the bottom right corner, in the picture\'s shape: ' + JSON.stringify(map));
    { // Zoomed in, the controls are the same row as at 1×, and nothing of the picture lies under them: the strip stays the panel's own
      // colour, a click there finds no cell, and the zoom level looks as it does at 1×. The map is faint until the pointer is on it.
      const now = await toolsNow(), strip = await page.evaluate(() => { const c = document.getElementById('picture'), dpr = c.width / c.clientWidth, top = window.calciumStep4.place(c.clientWidth, c.clientHeight).roomTop, d = c.getContext('2d').getImageData(0, 0, c.width, Math.floor(top * dpr)).data;
        let other = 0; for (let i = 0; i < d.length; i += 4) if (d[i] !== 0x10 || d[i + 1] !== 0x16 || d[i + 2] !== 0x1c) other++; const b = c.getBoundingClientRect();
        return { other, top, cell: window.calciumStep4.cellAt({ clientX: b.left + c.clientWidth / 2, clientY: b.top + top - 2 }), faint: getComputedStyle(document.getElementById('view-map')).opacity }; });
      assert.deepEqual([now.parts.map((part) => part.replace(/ · [\d.]+×$/, '')), now.order, now.level], [tools.parts.map((part) => part.replace(/ · [\d.]+×$/, '')), tools.order, tools.level], 'zoomed in, the controls keep their place, order and look');
      assert.ok(seen.y < strip.top && strip.other === 0 && strip.cell === -1, 'the picture is not drawn under the controls: ' + JSON.stringify(strip));
      assert.equal(strip.faint, '0.6', 'the map is faint at rest');
      // and it never covers the selected cell: where the cell would lie under it, the map moves to the bottom left corner
      const sweep = await page.evaluate(() => { const a = window.calciumStep4, c = document.getElementById('picture'), m = document.getElementById('view-map'), i = a.S.info.names.indexOf(a.S.cell), b = a.S.bounds[i], was = [a.Pic.zoom, a.Pic.cx, a.Pic.cy], seen = { covered: 0, left: 0, views: 0 };
        for (let gx = 0; gx <= 10; gx++) for (let gy = 0; gy <= 10; gy++) { a.setView(3, gx / 10 * a.S.info.width, gy / 10 * a.S.info.height); const at = a.place(c.clientWidth, c.clientHeight), box = c.getBoundingClientRect(), r = m.getBoundingClientRect();
          const cell = [box.left + at.x + b.x0 * at.z, Math.max(box.top + at.roomTop, box.top + at.y + b.y0 * at.z), box.left + at.x + b.x1 * at.z, box.top + at.y + b.y1 * at.z]; seen.views++;
          if (m.classList.contains('left')) seen.left++; if (cell[2] > r.left && cell[0] < r.right && cell[3] > r.top && cell[1] < r.bottom) seen.covered++; }
        a.setView(...was); return seen; });
      assert.ok(sweep.covered === 0 && sweep.left > 0 && sweep.left < sweep.views, 'the map never lies on the selected cell: ' + JSON.stringify(sweep)); }
    const inPanel = [Math.max(0, -seen.x / seen.z), Math.max(0, (seen.roomTop - seen.y) / seen.z), Math.min(seen.iw, (seen.w - seen.x) / seen.z), Math.min(seen.ih, (seen.h - seen.y) / seen.z)];   // (in view: under the strip of the controls)
    assert.ok(near(map.frame, [inPanel[0] * perPixel[0], inPanel[1] * perPixel[1], (inPanel[2] - inPanel[0]) * perPixel[0], (inPanel[3] - inPanel[1]) * perPixel[1]], 1e-6), 'its frame is the part of the recording in the panel: ' + map.frame);
    assert.ok(near([map.area.x, map.area.y, map.area.width, map.area.height], [inPanel[0], inPanel[1], inPanel[2] - inPanel[0], inPanel[3] - inPanel[1]], 1e-9) && map.edge > 300 && map.around < 150, 'drawn as a light frame, with the field outside it dimmed: ' + [map.edge, map.around]);   // (light summed over red, green and blue; the frame is 1.5 pixels wide, so a pixel on it is at least half white)
    const press = [Math.round(map.box[0] + 0.3 * map.box[2]), Math.round(map.box[1] + 0.5 * map.box[3])], release = [Math.round(map.box[0] + 0.65 * map.box[2]), Math.round(map.box[1] + 0.45 * map.box[3])];
    const middle = (view) => toRecording(view, view.left + view.midX, view.top + view.midY), onMap = (point) => [(point[0] - map.box[0]) / perPixel[0], (point[1] - map.box[1]) / perPixel[1]];
    await page.mouse.click(press[0], press[1]); seen = await viewNow(); assert.ok(near(middle(seen), onMap(press), 1e-6) && seen.zoom === 3, 'a click in the map puts that point in the middle of the view: ' + middle(seen) + ' for ' + onMap(press));
    const stillHeld = await drag(press, release, viewNow); seen = await viewNow();
    assert.ok(near(middle(stillHeld), onMap(release), 1e-6) && near(middle(seen), onMap(release), 1e-6), 'a drag in the map takes the view along, before the button is let go: ' + middle(stillHeld) + ' for ' + onMap(release));
    map = await mapNow(); assert.ok(near([map.frame[0] + map.frame[2] / 2, map.frame[1] + map.frame[3] / 2], [release[0] - map.box[0], release[1] - map.box[1]], 2), 'and the frame follows');
    assert.deepEqual([seen.cell, seen.zoom], ['Cell001', 3]); await page.click('#zoom-all'); assert.ok(await page.isHidden('#view-map'));
    // The first three times the pointer comes onto the picture and stays a second, a line by the controls says how to zoom and move
    // with the pointer alone (Cmd on a Mac). It goes when the pointer leaves. Then no more, in this session or a later one: Fiji keeps the
    // count. A pointer that only crosses the picture sees the line but does not use up one of the three.
    const hintWords = 'Pinch or ' + (await page.evaluate(() => /Mac|iPhone|iPad/.test(navigator.platform)) ? 'Cmd' : 'Ctrl') + ' + scroll to zoom. Drag to move.', hinted = () => page.evaluate(() => window.calciumStep4.Pic.hinted);
    await page.mouse.move(5, 5); await page.evaluate(() => { window.calciumStep4.Pic.hinted = 0; });
    await page.mouse.move(grip[0], grip[1]); assert.ok(await page.isVisible('#view-hint')); await page.mouse.move(5, 5); await page.waitForTimeout(1100); assert.equal(await hinted(), 0, 'crossing the picture does not count');
    for (let time = 1; time <= 3; time++) { await page.mouse.move(grip[0], grip[1]);
      assert.ok(await page.isVisible('#view-hint') && await text('#view-hint') === hintWords, 'the hint is shown, time ' + time + ': ' + await text('#view-hint'));
      assert.ok(await page.evaluate(() => { const hint = document.getElementById('view-hint').getBoundingClientRect(), stage = document.getElementById('picture-stage').getBoundingClientRect(), by = document.getElementById('view-tools').getBoundingClientRect(); return hint.left >= stage.left && hint.right <= stage.right && hint.height < 40 && (Math.abs(hint.top - by.bottom) <= 8 || Math.abs(hint.left - by.right) <= 8); }), 'one line, by the controls');
      await page.waitForFunction((n) => window.calciumStep4.Pic.hinted === n, time); await page.mouse.move(5, 5); assert.ok(await page.isHidden('#view-hint'), 'and goes when the pointer leaves'); }
    await page.mouse.move(grip[0], grip[1]); assert.ok(await page.isHidden('#view-hint'), 'the fourth time it is not shown'); await page.mouse.move(5, 5);
    await page.waitForFunction(async () => (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs['step4.viewhint'] === '3');

    // ---------------------------------------------------------------- the picture: plain outlines and numbered cells
    // What the picture's painter puts on a canvas: every outline it strokes (colour, width, dashes), every text (with its place, height
    // and width in the canvas's own pixels) and every open square of a key. false: the page's own picture, drawn again. true: a saved
    // picture (viewCanvas).
    const cellsDrawn = (saved) => page.evaluate(async (saved) => {
      const f = window.calciumStep4, P = CanvasRenderingContext2D.prototype, was = { stroke: P.stroke, fillText: P.fillText, strokeRect: P.strokeRect }, shown = document.getElementById('picture'), seen = { strokes: [], texts: [], squares: [] };
      const scaled = (ctx) => { const m = ctx.getTransform(), dpr = ctx.canvas === shown ? shown.width / shown.clientWidth : 1; return { x: (x, y) => (m.a * x + m.e) / dpr, y: (x, y) => (m.d * y + m.f) / dpr, by: m.a / dpr }; };
      P.stroke = function (path) { if (path instanceof Path2D) seen.strokes.push({ on: this.canvas, ink: this.strokeStyle, width: this.lineWidth * scaled(this).by, dashed: this.getLineDash().length > 0 }); return was.stroke.apply(this, arguments); };
      P.fillText = function (text, x, y) { const s = scaled(this); seen.texts.push({ on: this.canvas, text: String(text), x: s.x(x, y), y: s.y(x, y), high: parseFloat(/([\d.]+)px/.exec(this.font)[1]) * s.by, wide: this.measureText(text).width * s.by }); return was.fillText.apply(this, arguments); };
      P.strokeRect = function () { seen.squares.push({ on: this.canvas, ink: this.strokeStyle, dashed: this.getLineDash().length > 0 }); return was.strokeRect.apply(this, arguments); };
      let canvas = shown; try { if (saved) canvas = await f.viewCanvas(); else window.dispatchEvent(new Event('resize')); } finally { Object.assign(P, was); }
      const mine = (list) => list.filter((entry) => entry.on === canvas).map((entry) => { delete entry.on; return entry; });
      return { strokes: mine(seen.strokes), texts: mine(seen.texts), squares: mine(seen.squares), size: [canvas.width, canvas.height], record: f.pictureRecord({}), chips: Array.from(document.querySelectorAll('#call-legend .chip')).map((chip) => chip.textContent + (chip.querySelector('i.dash') ? ' (dashed)' : '')) }; }, saved);
    const lines = (drawn) => drawn.strokes.filter((s) => s.ink !== 'rgba(0, 0, 0, 0.55)').map((s) => s.ink + (s.dashed ? ' dashed' : '')), numbers = (drawn) => drawn.texts.filter((t) => /^\d+$/.test(t.text));
    // The two controls sit with the outline thickness: "Outlines" (By call, as before, or Plain) and "Number cells", off until ticked.
    assert.deepEqual(await page.evaluate(() => { const bar = document.querySelector('#picture-card .toolbar'), colours = document.getElementById('outline-colours'), box = document.getElementById('number-cells');
      return [bar.contains(colours) && bar.contains(box), document.getElementById('outline-width').nextElementSibling.contains(colours), colours.closest('label').textContent.replace(/By call|Plain/g, '').trim(), Array.from(colours.options).map((o) => o.textContent + (o.selected ? '*' : '')), box.closest('label').textContent.trim(), box.checked]; }),
      [true, true, 'Outlines', ['By call*', 'Plain'], 'Number cells', false]);
    // As it opens: every outline in the colour of its call, excluded and unscorable cells dashed, the selected cell in blue with its
    // name, and no numbers. On screen and saved.
    let onScreen = await cellsDrawn(false), kept = await cellsDrawn(true);
    assert.deepEqual(lines(onScreen), ['#ffa24a', '#aeb8c2 dashed', '#aeb8c2 dashed', '#3ddbb5', '#aeb8c2 dashed', '#000000', '#5bc0ff']); assert.deepEqual(onScreen.texts.map((t) => t.text), ['Cell001']);
    assert.deepEqual(onScreen.chips, ['responded', 'did not respond', 'excluded or unscorable (dashed)']);
    assert.deepEqual(lines(kept), ['#3ddbb5', '#ffa24a', '#aeb8c2 dashed', '#aeb8c2 dashed', '#3ddbb5', '#aeb8c2 dashed']); assert.deepEqual(numbers(kept), []);
    assert.deepEqual(kept.squares, [{ ink: '#3ddbb5', dashed: false }, { ink: '#ffa24a', dashed: false }, { ink: '#aeb8c2', dashed: false }]);
    assert.deepEqual(kept.record, { size: '8x6', area: 'whole picture', outlines_shown: true, outline_colours: { responded: '#3ddbb5', did_not_respond: '#ffa24a', excluded_or_unscorable: '#aeb8c2' }, cells_numbered: false, format: 'png', dpi: 300, compression: 'lzw' }, 'the file beside a saved picture says how the cells were drawn');
    // And exactly as the release saved it, pixel for pixel: a picture, a frame of a set, a picture of one group with thick outlines, the still overview.
    { const same = []; same.push(await againstRelease()); same.push(await againstRelease(11, true));
      await page.selectOption('#group', '2'); await page.selectOption('#outline-width', '2.8'); same.push(await againstRelease()); await page.selectOption('#group', '0'); await page.selectOption('#outline-width', '1'); await openCell('Cell001');
      await page.selectOption('#view', 'picture'); same.push(await againstRelease()); await page.selectOption('#view', 'raw'); await waitFrame(22);
      assert.deepEqual(same.map((one) => one.different), [0, 0, 0, 0], 'with By call and no numbers, saved pictures and frames are the release\'s: ' + JSON.stringify(same)); assert.ok(same.every((one) => one.lit > 1000 && one.size[0] === 1000)); }
    // Plain: every outline white on its dark halo, all of one width, the dashes kept; the key reads "your cells", and names the dashes
    // while an excluded or unscorable cell is drawn. The selected cell is still shown as selected.
    await page.selectOption('#outline-colours', 'plain'); onScreen = await cellsDrawn(false); kept = await cellsDrawn(true);
    assert.deepEqual(lines(onScreen), ['#ffffff', '#ffffff dashed', '#ffffff dashed', '#ffffff', '#ffffff dashed', '#000000', '#5bc0ff']);
    assert.equal(new Set(onScreen.strokes.filter((s) => s.ink === '#ffffff').map((s) => s.width.toFixed(4))).size, 1, 'no outline is heavier than another: a responder is not marked');
    assert.deepEqual(onScreen.chips, ['your cells', 'excluded or unscorable (dashed)']);
    assert.deepEqual(lines(kept), ['#ffffff', '#ffffff', '#ffffff dashed', '#ffffff dashed', '#ffffff', '#ffffff dashed']); assert.deepEqual(kept.squares, [{ ink: '#ffffff', dashed: false }, { ink: '#ffffff', dashed: true }]);
    assert.ok(kept.texts.some((t) => t.text === 'your cells') && kept.texts.some((t) => t.text === 'excluded or unscorable') && !kept.texts.some((t) => /respond/.test(t.text)), 'the saved key reads "your cells": ' + kept.texts.map((t) => t.text));
    assert.deepEqual(kept.record.outline_colours, { your_cells: '#ffffff', excluded_or_unscorable: '#ffffff' }); assert.notEqual((await againstRelease()).different, 0, '(and the comparison with the release does tell a plain picture apart)');
    await page.selectOption('#group', '1'); onScreen = await cellsDrawn(false); assert.deepEqual([lines(onScreen), onScreen.chips], [['#ffffff', '#000000', '#5bc0ff'], ['your cells']], 'Show still chooses the cells; with no dashed outline drawn the key is "your cells" alone'); await page.selectOption('#group', '0');
    // Number cells: every outline that is drawn gets its number, the number the page knows the cell by (Cell004: 4), white on a dark
    // halo, in the middle of a cell large enough to show it (as every cell of the fixture is) and readable (9 to 13 pixels high on
    // screen). The selected cell keeps its name tag. Saved pictures carry the same numbers, for the selected cell too, which they do not mark.
    await page.check('#number-cells'); onScreen = await cellsDrawn(false); kept = await cellsDrawn(true); seen = await viewNow();
    assert.deepEqual(onScreen.texts.map((t) => t.text), ['2', '3', '4', '5', '6', 'Cell001']); assert.deepEqual(numbers(kept).map((t) => t.text), ['1', '2', '3', '4', '5', '6']);
    assert.deepEqual(await page.evaluate(() => { const a = window.calciumStep4; return a.S.info.names.map(a.cellNumber).concat(a.cellTyped('4', true)); }), ['1', '2', '3', '4', '5', '6', 'Cell004'], 'the number that opens the cell when typed into the Cell box');
    const boundsOf = await page.evaluate(() => window.calciumStep4.S.bounds);
    for (const label of numbers(onScreen)) { const b = boundsOf[Number(label.text) - 1], spot = [seen.x + (b.x0 + b.x1) / 2 * seen.z, seen.y + (b.y0 + b.y1) / 2 * seen.z];
      assert.ok(near([label.x, label.y], spot, 0.01) && label.high > 8.99 && label.high < 13.01, 'the number ' + label.text + ' is in the middle of its cell, at a readable size: ' + JSON.stringify(label)); }
    for (const label of numbers(kept)) { const b = boundsOf[Number(label.text) - 1], scale = kept.size[0] / seen.iw; assert.ok(near([label.x / scale, label.y / scale], [(b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2], 0.01) && label.high / scale > 9.99, 'in a saved picture too: ' + JSON.stringify(label)); }
    assert.deepEqual([kept.record.cells_numbered, kept.record.outlines_shown], [true, true]);
    // A cell too small to show its number with its outline around it has the number beside it, off the outline; numbers never lie on
    // one another. The same rule saved. (cell 5 made a few pixels across, a neighbour of cell 4)
    const shrunk = await page.evaluate(() => { const a = window.calciumStep4, was = a.S.bounds.map((b) => Object.assign({}, b)), c = a.S.centres[4]; a.S.bounds[4] = { x0: c[0] - 1.5, x1: c[0] + 1.5, y0: c[1] - 1.5, y1: c[1] + 1.5 }; return was; });
    onScreen = await cellsDrawn(false); kept = await cellsDrawn(true); await page.evaluate((was) => { window.calciumStep4.S.bounds = was; }, shrunk);
    for (const [drawn, scale, left, top, where] of [[onScreen, seen.z, seen.x, seen.y, 'on screen'], [kept, kept.size[0] / seen.iw, 0, 0, 'saved']]) {
      const five = numbers(drawn).find((t) => t.text === '5'), c = boundsOf[4], mid = [left + (c.x0 + c.x1) / 2 * scale, top + (c.y0 + c.y1) / 2 * scale], halo = 0.3 * five.high;
      assert.ok(Math.abs(five.x - mid[0]) >= (five.wide + halo) / 2 + 1.5 * scale - 0.01 || Math.abs(five.y - mid[1]) >= five.high / 2 + 1.5 * scale - 0.01, where + ': the number of a small cell lies beside it: ' + JSON.stringify([five, mid]));
      const all = numbers(drawn); for (let i = 0; i < all.length; i++) for (let j = i + 1; j < all.length; j++) { const a = all[i], b = all[j];
        assert.ok(Math.abs(a.x - b.x) >= (a.wide + b.wide) / 2 - 0.01 || Math.abs(a.y - b.y) >= a.high - 0.01, where + ': no number on another: ' + JSON.stringify([a, b])); } }
    // Numbers that would lie on one another are moved apart, and no number runs off the picture.
    const crowded = await page.evaluate(() => { const a = window.calciumStep4, was = a.S.centres; a.S.centres = was.map((c, i) => (i === 1 || i === 2 ? [125, 110] : i === 3 ? [a.S.info.width - 1, a.S.info.height - 1] : c)); return was; });   // two cells with one middle, clear of the others; and one in the picture's corner
    onScreen = await cellsDrawn(false); kept = await cellsDrawn(true); await page.evaluate((was) => { window.calciumStep4.S.centres = was; }, crowded);
    for (const [drawn, scale, where] of [[onScreen, seen.z, 'on screen'], [kept, kept.size[0] / seen.iw, 'saved']]) { const label = (text) => numbers(drawn).find((t) => t.text === text), a = label('2'), b = label('3'), c = label('4'), left = drawn === onScreen ? seen.x : 0, top = drawn === onScreen ? seen.y : 0;
      assert.ok(Math.abs(a.x - b.x) >= (a.wide + b.wide) / 2 - 0.01 || Math.abs(a.y - b.y) >= a.high - 0.01, where + ': two numbers at one place are moved apart: ' + JSON.stringify([a, b]));
      for (const one of [a, b, c]) assert.ok(one.x - one.wide / 2 >= left - 0.01 && one.y - one.high / 2 >= top - 0.01 && one.x + one.wide / 2 <= left + seen.iw * scale + 0.01 && one.y + one.high / 2 <= top + seen.ih * scale + 0.01, where + ': the number ' + one.text + ' stays inside the picture: ' + JSON.stringify(one)); }
    // The numbers go with the outlines: Show chooses them, and without outlines there are none (the two controls then wait for the outlines).
    await page.selectOption('#group', '1'); assert.deepEqual((await cellsDrawn(false)).texts.map((t) => t.text), ['5', 'Cell001']); assert.deepEqual(numbers(await cellsDrawn(true)).map((t) => t.text), ['1', '5']); await page.selectOption('#group', '0');
    await page.uncheck('#outlines'); onScreen = await cellsDrawn(false); kept = await cellsDrawn(true);
    assert.deepEqual([lines(onScreen), onScreen.texts.map((t) => t.text), lines(kept), numbers(kept), kept.squares], [['#000000', '#5bc0ff'], ['Cell001'], [], [], []]); assert.deepEqual([kept.record.outlines_shown, kept.record.outline_colours, kept.record.cells_numbered], [false, {}, false]);
    assert.ok(await page.isDisabled('#outline-colours') && await page.isDisabled('#number-cells') && await page.isHidden('#call-legend')); await page.check('#outlines'); assert.ok(await page.isEnabled('#outline-colours') && await page.isEnabled('#number-cells'));
    // Both choices are kept for the next recording, as the outline thickness is. A picture saved through Export carries them, with
    // its record; so do frames (checked with Save frames below, which runs with the choices as the page opens).
    await page.waitForFunction(async () => { const prefs = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs; return prefs['step4.plain'] === '1' && prefs['step4.numbers'] === '1'; });
    await page.locator('#export-picture > summary').click(); await page.click('#save-picture-folder');
    await page.waitForFunction(() => /^Saved \S+\.png/.test(document.getElementById('picture-saved').textContent), null, { polling: 100, timeout: 20000 });
    const plainName = (await text('#picture-saved')).match(/^Saved (\S+)\.png/)[1], findFile = (dir, name) => { for (const entry of fs.readdirSync(dir, { withFileTypes: true })) { const full = path.join(dir, entry.name); if (entry.isDirectory()) { const below = findFile(full, name); if (below) return below; } else if (entry.name === name) return full; } return null; };
    const plainRecord = JSON.parse(fs.readFileSync(findFile(out, plainName + '.metadata.json'), 'utf8')).picture, plainPng = fs.readFileSync(findFile(out, plainName + '.png'));
    assert.deepEqual([plainRecord.outlines_shown, plainRecord.outline_colours, plainRecord.cells_numbered], [true, { your_cells: '#ffffff', excluded_or_unscorable: '#ffffff' }, true], 'the record written beside the saved picture'); assert.equal(plainPng.readUInt32BE(16), 1000);
    await page.selectOption('#outline-colours', 'call'); await page.uncheck('#number-cells'); onScreen = await cellsDrawn(false);
    assert.deepEqual([lines(onScreen), onScreen.texts.map((t) => t.text)], [['#ffa24a', '#aeb8c2 dashed', '#aeb8c2 dashed', '#3ddbb5', '#aeb8c2 dashed', '#000000', '#5bc0ff'], ['Cell001']], 'back to By call and no numbers, the picture is as it opened');
    await page.waitForFunction(async () => { const prefs = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs; return prefs['step4.plain'] === '0' && prefs['step4.numbers'] === '0'; });
    // The plot's time axis follows the same rule, and its hint says so.
    const plotBox = await page.locator('#plot').boundingBox(), timeZoomed = () => page.evaluate(() => window.calciumStep4.Plot.x0 !== null);
    await page.mouse.move(plotBox.x + plotBox.width / 2, plotBox.y + plotBox.height / 2); await page.mouse.wheel(0, -300); await settledPaint();
    assert.equal(await timeZoomed(), false, 'scrolling over the plot does not zoom it');
    await page.keyboard.down('Control'); await page.mouse.wheel(0, -300); await page.keyboard.up('Control'); await settledPaint();
    assert.equal(await timeZoomed(), true, 'Ctrl + scroll zooms the time axis'); await page.mouse.dblclick(plotBox.x + plotBox.width / 2, plotBox.y + plotBox.height / 2);
    assert.equal(await timeZoomed(), false, 'double-click resets it'); await page.evaluate(() => window.scrollTo(0, 0));
    assert.match(await text('#plot-card .caption .help'), /Ctrl or ⌘ \+ scroll to zoom/);
    // The settings: a label never wraps. Beside its control at the panel's usual width, above it when the panel is narrow.
    const settingRows = () => page.evaluate(() => [...document.querySelectorAll('#side .row')].filter((r) => r.getBoundingClientRect().height).map((r) => { const label = r.querySelector('span'), a = label.getBoundingClientRect(), b = r.children[1].getBoundingClientRect(); return { label: label.textContent, wrapped: a.height > 30 || label.scrollWidth > label.clientWidth + 1, stacked: b.top >= a.bottom - 1 }; }));
    const sideBefore = await page.evaluate(() => { const was = { width: document.documentElement.style.getPropertyValue('--side-w'), open: [...document.querySelectorAll('#side details')].map((d) => d.open) }; for (const d of document.querySelectorAll('#side details')) d.open = true; return was; });
    for (const [width, stacked] of [[350, false], [300, true], [260, true], [460, false]]) {
      await page.evaluate((w) => document.documentElement.style.setProperty('--side-w', w + 'px'), width); const rows = await settingRows();
      assert.ok(rows.length >= 12, 'the settings rows are on the page'); assert.deepEqual(rows.filter((r) => r.wrapped).map((r) => r.label), [], 'no label wraps at ' + width);
      assert.deepEqual(rows.filter((r) => r.stacked !== stacked).map((r) => r.label), [], (stacked ? 'every label above its control at ' : 'every label beside its control at ') + width);
    }
    await page.evaluate((was) => { if (was.width) document.documentElement.style.setProperty('--side-w', was.width); else document.documentElement.style.removeProperty('--side-w'); [...document.querySelectorAll('#side details')].forEach((d, i) => { d.open = was.open[i]; }); }, sideBefore);
    // The frame rate is shown to four significant figures; the stored value keeps every digit.
    assert.deepEqual(await page.evaluate(() => [1.0000537066691941, 0.999195, 2, 29.97, 0.5, 1000.04, '', null].map(window.calciumStep4.rateText)), ['1.000', '0.9992', '2', '29.97', '0.5', '1000', '', '']);
    // The settings panel folds away from its own edge.
    // One shape and one size: every button, list and field is the control height with the page's one corner, every box has that corner,
    // and nothing is a pill or a circle. Tabs, menu entries and links written as buttons are the stated exceptions.
    const oddShapes = () => page.evaluate(() => {
      const root = getComputedStyle(document.documentElement), corner = root.getPropertyValue('--radius').trim(), high = parseFloat(root.getPropertyValue('--control-h')), odd = [];
      const name = (el) => el.id || (typeof el.className === 'string' && el.className) || el.tagName;
      for (const el of document.querySelectorAll('button, select, input:not([type=checkbox]):not([type=range]):not([type=file]):not([type=radio]), details.export > summary, .quality-chip > summary, .link-toggle, .file-button')) {
        const box = el.getBoundingClientRect(); if (!box.width || !box.height) continue;
        if (!el.matches('.tab, .export-menu button, .repair, .jump, .flag-filter, .save-note button') && Math.abs(box.height - high) > 0.6) odd.push('height ' + box.height + ': ' + name(el));
        if (!el.matches('.tab, .seg button, .export-menu button') && getComputedStyle(el).borderTopLeftRadius !== corner) odd.push('corner ' + getComputedStyle(el).borderTopLeftRadius + ': ' + name(el));
      }
      for (const el of document.querySelectorAll('.card, .tile, .scroller, .run-box, .export-menu, .quality-pop, .badge, .pill')) { if (el.getBoundingClientRect().width && getComputedStyle(el).borderTopLeftRadius !== corner) odd.push('corner ' + getComputedStyle(el).borderTopLeftRadius + ': ' + name(el)); }
      for (const el of document.querySelectorAll('body *')) { const r = getComputedStyle(el).borderTopLeftRadius; if (el.getBoundingClientRect().width && (r.endsWith('%') || parseFloat(r) > parseFloat(corner))) odd.push('rounder than the page corner (' + r + '): ' + name(el)); }
      return odd;
    });
    // The results: the sentence on a row of its own; under it the counts, the flagged list and the results folder on one line, level with each other.
    assert.deepEqual(await page.evaluate(() => { const top = (id) => Math.round(document.getElementById(id).getBoundingClientRect().top), bottom = (id) => Math.round(document.getElementById(id).getBoundingClientRect().bottom), row = ['tile-registered', 'tile-excluded', 'tile-unscorable', 'quality-summary', 'open-folder'].filter((id) => document.getElementById(id).getBoundingClientRect().width);
      return { under: row.every((id) => top(id) >= bottom('tile-responding')), tops: [...new Set(row.map(top))].length, bottoms: [...new Set(row.map(bottom))].length, shown: row.length >= 4 }; }), { under: true, tops: 1, bottoms: 1, shown: true });
    assert.deepEqual(await oddShapes(), [], 'controls and boxes share one height and one corner');
    await page.locator('#export-plot > summary').click(); assert.deepEqual(await oddShapes(), [], 'also with a menu open'); await page.keyboard.press('Escape');
    // Its arrow has a strip of its own beside the work area and a clear gap on both sides, so it never sits on or against the results, with the settings shown or hidden.
    assert.ok(await page.evaluate(() => document.getElementById('side-edge').contains(document.getElementById('toggle-settings')) && !document.querySelector('main').contains(document.getElementById('toggle-settings'))));
    const arrowClear = () => page.evaluate(() => { const arrow = document.getElementById('toggle-settings').getBoundingClientRect(), work = document.querySelector('main').getBoundingClientRect(); return arrow.width > 0 && arrow.left >= 2 && work.left - arrow.right >= 2; });
    assert.ok(await arrowClear(), 'arrow beside the work area, settings shown'); assert.ok(await page.isVisible('#resize-side'));
    assert.ok(await page.evaluate(() => document.getElementById('toggle-settings').getBoundingClientRect().left - document.getElementById('side').getBoundingClientRect().right >= 2), 'and clear of the settings');
    await page.click('#toggle-settings'); assert.ok(await page.isHidden('#side')); assert.equal(await page.getAttribute('#toggle-settings', 'aria-label'), 'Show settings');
    assert.ok(await arrowClear(), 'arrow beside the work area, settings hidden'); assert.ok(await page.isHidden('#resize-side'));
    assert.ok(await page.evaluate(() => document.querySelector('main').getBoundingClientRect().left <= 24), 'the work area keeps the page margin');
    await page.click('#toggle-settings'); assert.ok(await page.isVisible('#side')); assert.equal(await page.getAttribute('#toggle-settings', 'aria-label'), 'Hide settings');
    // Save frames also writes one picture with all the saved frames, built from the same labelled frames.
    await page.locator('#export-picture > summary').click(); await page.click('#save-frames'); await page.fill('#frames-every', '10');
    assert.match(await text('#frames-count'), /^3 pictures and one with all of them/);
    await page.click('#frames-go');
    await page.waitForFunction(() => /^Saved 3 frames and one combined picture in plots\//.test(document.getElementById('picture-saved').textContent), null, { polling: 100, timeout: 20000 });
    const framesFolder = (await text('#picture-saved')).match(/plots\/(frames_\S+)/)[1];
    const findFolder = (dir, name) => { for (const entry of fs.readdirSync(dir, { withFileTypes: true })) { if (!entry.isDirectory()) continue; const full = path.join(dir, entry.name); if (entry.name === name) return full; const below = findFolder(full, name); if (below) return below; } return null; };
    const framesOnDisk = findFolder(out, framesFolder); assert.ok(framesOnDisk, 'the frames folder is on disk: ' + framesFolder);
    assert.deepEqual(fs.readdirSync(framesOnDisk).filter((name) => /\.png$/.test(name)).sort(), ['all_frames_combined.png', 'frame_01.png', 'frame_11.png', 'frame_21.png']);
    const framesRecord = JSON.parse(fs.readFileSync(path.join(framesOnDisk, 'about_these_frames.metadata.json'), 'utf8')), layoutSaved = framesRecord.combined_picture;
    assert.ok(layoutSaved && layoutSaved.columns * layoutSaved.rows >= 3 && layoutSaved.file === 'all_frames_combined', 'how the combined picture is laid out is recorded');
    assert.match(framesFolder, /^frames_as_recorded_every_10_\d{8}-\d{6}$/); assert.equal(framesRecord.area, 'whole picture', 'frames of the whole picture, and their record says so');
    assert.deepEqual([framesRecord.outlines_shown, framesRecord.outline_colours, framesRecord.cells_numbered], [true, { responded: '#3ddbb5', did_not_respond: '#ffa24a', excluded_or_unscorable: '#aeb8c2' }, false], 'and how the cells were drawn in the frames');

    // ---------------------------------------------------------------- the picture: saving the area on screen, and a baseline and peak pair
    // "Area", in the picture's Export menu and in Save frames: "Whole picture", as it opens and as before, or "Area on screen". One
    // choice, shown in both places. At 1× the area on screen is the whole picture; zoomed in, both rows stand out, since the choice
    // then decides what is saved.
    const areaRows = () => page.evaluate(() => ['save-area', 'frames-area'].map((id) => { const box = document.getElementById(id), row = box.closest('label');
      return [row.textContent.replace(/Whole picture|Area on screen/g, '').trim(), Array.from(box.options).map((o) => o.textContent + (o.selected ? '*' : '')).join(' | '), !!box.closest(id === 'save-area' ? '#export-picture .export-menu' : '#frames-options'), row.classList.contains('matters')]; }));
    const chooseArea = async (value) => { await page.locator('#export-picture > summary').click(); await page.selectOption('#save-area', value); await page.keyboard.press('Escape'); };
    const areaNow = () => page.evaluate(() => { const a = window.calciumStep4; return { saved: a.savedArea(), seen: a.viewArea(), zoom: a.Pic.zoom }; }), frameOnScreen = await page.evaluate(() => window.calciumStep4.S.frame);
    assert.deepEqual(await areaRows(), [['Area', 'Whole picture* | Area on screen', true, false], ['Area', 'Whole picture* | Area on screen', true, false]]);
    assert.equal((await areaNow()).saved, null); await chooseArea('screen');
    assert.deepEqual([(await areaNow()).saved, await page.inputValue('#frames-area'), await app('cell')], [null, 'screen', 'Cell001'], 'at 1× the area on screen is the whole picture; and the list in Save frames follows');
    drawn = await rectangle(); await page.keyboard.down('Shift'); await drag(drawn.a, drawn.b); await page.keyboard.up('Shift'); let area = await areaNow();
    { const s = area.seen, edges = [Math.floor(s.x + 1e-6), Math.floor(s.y + 1e-6), Math.ceil(s.x + s.width - 1e-6), Math.ceil(s.y + s.height - 1e-6)];
      assert.deepEqual(area.saved, { x: edges[0], y: edges[1], width: edges[2] - edges[0], height: edges[3] - edges[1], zoom: area.zoom }, 'Area on screen: the part of the recording in the panel, to whole pixels, and the zoom');
      assert.ok(area.zoom > 2 && edges[0] <= drawn.part[0] && edges[1] <= drawn.part[1] && edges[2] >= drawn.part[2] && edges[3] >= drawn.part[3] && area.saved.width < 0.6 * seen.iw, 'it holds the rectangle the view was zoomed to: ' + JSON.stringify(area.saved)); }
    assert.deepEqual((await areaRows()).map((row) => row[3]), [true, true], 'zoomed in, the Area rows stand out');
    // With "Whole picture" a zoomed view changes nothing: saved pictures and frames are the release's, pixel for pixel. With "Area on screen" they are not.
    const zoomedArea = area.saved; assert.equal((await againstRelease()).different, -1, '(a picture of the area is another picture)'); await chooseArea('whole');
    { assert.equal((await areaNow()).saved, null, 'Whole picture: the whole field, however the picture is zoomed');
      const same = [await againstRelease(), await againstRelease(11, true)]; assert.deepEqual(same.map((one) => one.different), [0, 0], 'zoomed in with Whole picture, saved pictures and frames are the release\'s: ' + JSON.stringify(same)); }
    // What a saved picture's painter does, in the file's own pixels: the outlines it strokes, the texts it writes and the boxes it fills,
    // the part of the recording it draws (from: left, top, width, height in the recording; to: where in the file; smoothed or not),
    // and where it cuts. what: 'view' (viewCanvas: frame, labelled, area) or 'pair' (pairCanvas: baseline, peak, area).
    const savedDrawing = (what, ...given) => page.evaluate(async ([what, given]) => {
      const f = window.calciumStep4, P = CanvasRenderingContext2D.prototype, was = { stroke: P.stroke, fillText: P.fillText, fillRect: P.fillRect, drawImage: P.drawImage, rect: P.rect, clip: P.clip }, seen = { strokes: [], texts: [], boxes: [], images: [], cuts: [] };
      const at = (ctx, x, y) => { const m = ctx.getTransform(); return [m.a * x + m.e, m.d * y + m.f, m.a]; };
      P.stroke = function (path) { if (path instanceof Path2D) seen.strokes.push({ on: this.canvas, ink: this.strokeStyle, width: this.lineWidth * at(this, 0, 0)[2], dashed: this.getLineDash().length > 0 }); return was.stroke.apply(this, arguments); };
      P.fillText = function (text, x, y) { const [px, py, by] = at(this, x, y); seen.texts.push({ on: this.canvas, text: String(text), x: px, y: py, high: parseFloat(/([\d.]+)px/.exec(this.font)[1]) * by, wide: this.measureText(text).width * by, centred: this.textAlign === 'center' }); return was.fillText.apply(this, arguments); };
      P.fillRect = function (x, y, w, h) { const [px, py, by] = at(this, x, y); seen.boxes.push({ on: this.canvas, ink: this.fillStyle, x: px, y: py, wide: w * by, high: h * by }); return was.fillRect.apply(this, arguments); };
      P.drawImage = function (image, sx, sy, sw, sh, dx, dy, dw, dh) { if (arguments.length === 9) { const [px, py, by] = at(this, dx, dy); seen.images.push({ on: this.canvas, from: [sx, sy, sw, sh], to: [px, py, dw * by, dh * by], smoothed: this.imageSmoothingEnabled }); } return was.drawImage.apply(this, arguments); };
      P.rect = function (x, y, w, h) { const [px, py, by] = at(this, x, y); this.lastRect = [px, py, w * by, h * by]; return was.rect.apply(this, arguments); };
      P.clip = function () { seen.cuts.push({ on: this.canvas, box: this.lastRect }); return was.clip.apply(this, arguments); };
      let canvas; try { canvas = await (what === 'pair' ? f.pairCanvas : f.viewCanvas)(...given); } finally { Object.assign(P, was); }
      const mine = (list) => list.filter((entry) => entry.on === canvas).map((entry) => { delete entry.on; return entry; });
      return { strokes: mine(seen.strokes), texts: mine(seen.texts), boxes: mine(seen.boxes), images: mine(seen.images), cuts: mine(seen.cuts), size: [canvas.width, canvas.height] }; }, [what, given]);
    const outlinesOf = (made) => made.strokes.filter((s) => s.ink !== 'rgba(0, 0, 0, 0.55)'), barOf = (made) => { const words = made.texts.filter((t) => / µm$/.test(t.text)), bars = made.boxes.filter((b) => b.ink === '#ffffff' && Math.abs(b.high - 6 * 4) < 0.01 && b.wide > b.high); return { words, bars }; };   // (the bar is 6 pixels of a saved picture's pen thick; these pictures are 4 times the recording)
    // A picture of the area on screen is as wide as a picture of the whole field. The area is enlarged to fill it and cut cleanly at
    // its edges; the outlines, the numbers, the key and a frame's label are drawn at the saved picture's own scale: the same pens and
    // the same text heights as in a picture of the whole field, not stretched with the area. The selected cell is not marked, as in
    // every saved picture. A cell is numbered when its middle lies in the area, and every number stays inside the picture.
    await page.check('#number-cells'); const wholeMade = await savedDrawing('view', 12, false, null), areaMade = await savedDrawing('view', 12, false, zoomedArea), areaFrame = await savedDrawing('view', 12, true, zoomedArea), wholeFrame = await savedDrawing('view', 12, true, null);
    { const k = 1000 / zoomedArea.width, high = Math.round(1000 * zoomedArea.height / zoomedArea.width), one = areaMade.images[0];
      assert.deepEqual([areaMade.size[0], wholeMade.size[0], areaFrame.size], [1000, 1000, [1000, high]], 'as wide as a picture of the whole field, in the shape of the area');
      assert.ok(areaMade.images.length === 1 && near(one.from, [zoomedArea.x, zoomedArea.y, zoomedArea.width, zoomedArea.height], 1e-9) && near(one.to, [0, 0, 1000, high], 1e-6), 'drawn from the recording\'s own pixels of that area: ' + JSON.stringify(one));
      assert.equal(one.smoothed, false, 'an enlarged area keeps the recording\'s pixels as squares, whatever the zoom on screen'); assert.equal(wholeMade.images[0].smoothed, true, 'the whole field is drawn as it always was');
      assert.ok(areaMade.cuts.length === 1 && near(areaMade.cuts[0].box, [0, 0, 1000, high], 1e-6) && wholeMade.cuts.length === 0, 'outlines end cleanly at the edge of the area: ' + JSON.stringify(areaMade.cuts));
      assert.deepEqual(outlinesOf(areaMade).map((s) => [s.ink, s.dashed, +s.width.toFixed(4)]), outlinesOf(wholeMade).map((s) => [s.ink, s.dashed, +s.width.toFixed(4)]), 'the outlines keep their pens'); assert.ok(near([outlinesOf(areaMade)[0].width], [7.2], 1e-6));
      const inside = boundsOf.map((b, i) => [i + 1, (b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2]).filter(([, x, y]) => x >= zoomedArea.x && x <= zoomedArea.x + zoomedArea.width && y >= zoomedArea.y && y <= zoomedArea.y + zoomedArea.height);
      assert.ok(inside.length >= 2 && inside.length < 6, 'the area holds some of the cells: ' + inside.map((one) => one[0]));
      assert.deepEqual(numbers(areaMade).map((t) => t.text), inside.map((one) => String(one[0])), 'a cell whose middle is in the area gets its number'); assert.deepEqual(numbers(wholeMade).map((t) => t.text), ['1', '2', '3', '4', '5', '6']);
      for (const label of numbers(areaMade)) { const [, x, y] = inside.find((one) => String(one[0]) === label.text);
        assert.ok(near([label.x, label.y], [(x - zoomedArea.x) * k, (y - zoomedArea.y) * k], 0.01) && label.high >= 12 * 4 - 0.01 && label.high <= 22 * 4 + 0.01, 'in the middle of its cell, at a saved picture\'s size: ' + JSON.stringify(label));
        assert.ok(label.x - label.wide / 2 >= 0 && label.x + label.wide / 2 <= 1000 && label.y - label.high / 2 >= 0 && label.y + label.high / 2 <= high, 'inside the picture'); }
      const words = (made) => made.texts.filter((t) => !/^\d+$/.test(t.text) && !/ µm$/.test(t.text)).map((t) => [t.text, +t.x.toFixed(2), +t.high.toFixed(2)]);
      assert.deepEqual(words(areaMade), words(wholeMade), 'the key is the whole picture\'s, word for word and size for size'); assert.ok(words(areaMade).some((word) => word[0] === 'responded' && word[2] === 80));
      assert.deepEqual(words(areaFrame), words(wholeFrame), 'and so is a frame\'s label'); assert.deepEqual(words(areaFrame).map((word) => word[0]), ['Frame 12 · 5.5 s · as recorded']);
      assert.ok(!areaMade.strokes.some((s) => s.ink === '#5bc0ff') && !areaMade.texts.some((t) => t.text === 'Cell001'), 'the selected cell is not marked in a saved picture');
      // The scale bar of a picture or frame of an area: a round length, at most a fifth of the area's width, written over a bar of
      // exactly that length at the area's pixel size (the fixture's pixels are 0.8 µm). A picture has it at the right of its key,
      // under the picture, where it covers nothing; a frame of a set, which has no key, on its bottom right corner. A picture of
      // the whole field has none, as before.
      const micron = await page.evaluate(() => window.calciumStep4.S.info.micronsPerPixel), across = zoomedArea.width * micron, round = [1, 2, 5, 10, 20, 50, 100].filter((length) => length <= across / 5).pop();
      assert.equal(micron, 0.8); assert.deepEqual(await page.evaluate(() => [4.9, 5, 9.99, 10, 55, 100, 895, 1000, 2.5].map(window.calciumStep4.scaleBarLength)), [0.5, 1, 1, 2, 10, 20, 100, 200, 0.5], 'round lengths: 1, 2 or 5 times a power of ten, at most a fifth of the width');
      for (const made of [areaMade, areaFrame]) { const bar = barOf(made);
        assert.ok(bar.words.length === 1 && bar.bars.length === 1 && bar.words[0].text === round + ' µm', 'one scale bar that says its length: ' + JSON.stringify(bar.words));
        assert.ok(near([bar.bars[0].wide], [round / micron * k], 1e-6), 'as long as it says, at the area\'s pixel size: ' + bar.bars[0].wide + ' for ' + round + ' µm');
        const [word, line] = [bar.words[0], bar.bars[0]], right = Math.max(line.x + line.wide, word.x + word.wide / 2), others = made.texts.filter((t) => t !== word && !/^\d+$/.test(t.text));
        assert.ok(word.centred && near([word.x], [line.x + line.wide / 2], 1e-6) && word.y < line.y && right <= 1000 - 14 * 4 + 0.01 && right > 900, 'the length over the bar, at the right: ' + JSON.stringify([word, line]));
        if (made === areaMade) assert.ok(word.y - word.high / 2 >= high && line.y + line.high <= made.size[1] && others.every((t) => t.y - t.high / 2 >= line.y + line.high || t.y + t.high / 2 <= word.y - word.high / 2 || t.x + t.wide <= Math.min(line.x, word.x - word.wide / 2)), 'a picture has it in the strip under the picture, clear of the key: ' + JSON.stringify(others));
        else assert.ok(line.y + line.high <= high && word.y - word.high / 2 > high / 2, 'a frame has it on its bottom right corner'); }
      assert.deepEqual([barOf(wholeMade).words.length, barOf(wholeFrame).words.length], [0, 0], 'a picture of the whole field has no scale bar, as before'); }
    await page.uncheck('#number-cells');
    // The pixels of a saved area are the recording's: the same as in a picture of the whole field. Checked without outlines, clear of the
    // frame's label and of the scale bar, against a frame of the whole field saved at the recording's own size: each of its pixels
    // in the area is a square of 5 × 5 in the picture of the area, whether it was seen at 5× (where the page shows pixels as squares) or
    // at 2× (where the page smooths): the saved picture enlarges the recording 5 times either way, and an enlarged pixel is not smoothed.
    const saveAs = (size, dpi) => page.evaluate(([size, dpi]) => { for (const [id, value] of [['save-size', size], ['save-dpi', dpi]]) { const box = document.getElementById(id); box.value = value; box.dispatchEvent(new Event('change', { bubbles: true })); } }, [size, dpi]);
    await page.uncheck('#outlines'); await saveAs('screen', '150');
    const squares = await page.evaluate(async () => { const a = window.calciumStep4, part = { x: 20, y: 45, width: 50, height: 30 }, whole = await a.viewCanvas(12, true), result = { whole: [whole.width, whole.height] };
      const w = whole.getContext('2d').getImageData(0, 0, whole.width, whole.height).data;
      for (const zoom of [5, 2]) { const cut = await a.viewCanvas(12, true, Object.assign({ zoom }, part)), c = cut.getContext('2d').getImageData(0, 0, cut.width, cut.height).data; let different = 0, compared = 0, lit = 0;
        for (let y = 50; y < 92; y++) for (let x = 0; x < cut.width; x++) for (let channel = 0; channel < 3; channel++) { const mine = c[4 * (y * cut.width + x) + channel], theirs = w[4 * ((part.y + Math.floor(y / 5)) * whole.width + part.x + Math.floor(x / 5)) + channel]; compared++; if (mine !== theirs) different++; if (mine > 200) lit++; }
        result[zoom] = { size: [cut.width, cut.height], compared, different, lit }; }
      return result; });
    assert.deepEqual([squares.whole, squares[5].size], [[250, 150], [250, 150]]); assert.ok(squares[5].different === 0 && squares[5].lit > 5000 && squares[5].lit < squares[5].compared - 5000, 'a saved area is the same crop of the recording, pixel for pixel: ' + JSON.stringify(squares));
    assert.ok(squares[2].different === 0 && squares[2].lit > 5000, 'also an area seen at 2×: its pixels are squares with hard edges, not blended: ' + JSON.stringify(squares[2]));
    await page.check('#outlines'); await saveAs('8x6', '300');
    // Through the Export menu, zoomed in. With "Whole picture": the file is the whole field, its name does not say area, and its record
    // says "whole picture". With "Area on screen": the name says area, the picture has the area's shape, and the record gives the area in
    // pixels of the recording and in µm, the zoom, and the scale bar.
    const exported = async () => { await page.locator('#export-picture > summary').click(); await page.click('#save-picture-folder'); await page.waitForFunction((before) => { const now = document.getElementById('picture-saved').textContent; return /^Saved \S+\.png/.test(now) && now !== before; }, await text('#picture-saved'), { polling: 100, timeout: 20000 });
      const name = (await text('#picture-saved')).match(/^Saved (\S+)\.png/)[1], png = fs.readFileSync(findFile(out, name + '.png')); return { name, size: [png.readUInt32BE(16), png.readUInt32BE(20)], record: JSON.parse(fs.readFileSync(findFile(out, name + '.metadata.json'), 'utf8')) }; };
    let file = await exported(); assert.match(file.name, new RegExp('^Fixture_recording_picture_as_recorded_frame' + frameOnScreen + '_\\d{8}-\\d{6}$')); assert.deepEqual([file.size[0], file.record.picture.area, 'area_pixels' in file.record.picture], [1000, 'whole picture', false]);
    const wholeHeight = file.size[1]; await page.waitForTimeout(1000); await chooseArea('screen'); file = await exported();   // (a second later: the names carry the time to the second)
    assert.match(file.name, new RegExp('^Fixture_recording_picture_as_recorded_frame' + frameOnScreen + '_area_\\d{8}-\\d{6}$'), 'the file name says area'); assert.ok(file.size[0] === 1000 && file.size[1] !== wholeHeight && file.size[1] > Math.round(1000 * zoomedArea.height / zoomedArea.width), 'the area, with its key under it: ' + file.size);
    { const p = file.record.picture, micro = (v) => +(v * 0.8).toFixed(3);
      assert.deepEqual([p.area, p.area_pixels, p.area_um, p.zoom, p.scale_bar_um], ['area on screen', { left: zoomedArea.x, top: zoomedArea.y, width: zoomedArea.width, height: zoomedArea.height }, { left: micro(zoomedArea.x), top: micro(zoomedArea.y), width: micro(zoomedArea.width), height: micro(zoomedArea.height) }, +zoomedArea.zoom.toFixed(2), [1, 2, 5, 10, 20].filter((length) => length <= zoomedArea.width * 0.8 / 5).pop()], 'the record gives the area in pixels and in µm, the zoom and the scale bar: ' + JSON.stringify(p));
      assert.deepEqual([p.size, p.outlines_shown, p.format, p.dpi, file.record.frame, file.record.view], ['8x6', true, 'png', 300, frameOnScreen, 'raw']); }
    // Save frames with "Area on screen": every frame shows the area that was on screen when saving began, also when the view is moved
    // meanwhile (here, as the first frame is drawn). The folder's name says area and its record gives the area.
    await page.evaluate(() => { const a = window.calciumStep4, P = CanvasRenderingContext2D.prototype, was = P.drawImage; a.framesDrawn = [];
      P.drawImage = function (image, sx, sy, sw, sh) { if (arguments.length === 9) { a.framesDrawn.push([sx, sy, sw, sh]); if (a.framesDrawn.length === 1) document.getElementById('zoom-in').click(); } return was.apply(this, arguments); }; a.framesUnhook = () => { P.drawImage = was; }; });
    await page.locator('#export-picture > summary').click(); await page.click('#save-frames'); assert.equal(await page.inputValue('#frames-area'), 'screen'); await page.fill('#frames-every', '10'); await page.click('#frames-go');
    await page.waitForFunction(() => /^Saved 3 frames and one combined picture in plots\/frames_as_recorded_area_every_10_/.test(document.getElementById('picture-saved').textContent), null, { polling: 100, timeout: 20000 });
    { const drawnFrames = await page.evaluate(() => { const a = window.calciumStep4; a.framesUnhook(); return a.framesDrawn; }), later = await areaNow(), folder = findFolder(out, (await text('#picture-saved')).match(/plots\/(frames_\S+)/)[1]), about = JSON.parse(fs.readFileSync(path.join(folder, 'about_these_frames.metadata.json'), 'utf8'));
      assert.deepEqual(drawnFrames, [1, 2, 3].map(() => [zoomedArea.x, zoomedArea.y, zoomedArea.width, zoomedArea.height]), 'the three frames are drawn from one area: ' + JSON.stringify(drawnFrames));
      assert.ok(later.zoom > zoomedArea.zoom && later.saved.width < zoomedArea.width, 'although the view was moved while they were saved');
      assert.deepEqual([about.area, about.area_pixels, about.frames], ['area on screen', { left: zoomedArea.x, top: zoomedArea.y, width: zoomedArea.width, height: zoomedArea.height }, [1, 11, 21]]);
      const sizes = ['frame_01.png', 'frame_11.png', 'frame_21.png'].map((name) => { const png = fs.readFileSync(path.join(folder, name)); return png.readUInt32BE(16) + 'x' + png.readUInt32BE(20); });
      assert.deepEqual(sizes, [1, 2, 3].map(() => '1000x' + Math.round(1000 * zoomedArea.height / zoomedArea.width)), 'and are pictures of one size: ' + sizes); }
    await page.click('#zoom-all'); await chooseArea('whole');

    // Baseline and peak: one picture of two frames side by side. Its two boxes open filled in. Baseline: the last of the F₀ frames
    // (here 1–5). Peak: the selected cell's peak in the response window (Cell001: 6.0 s, frame 13; Cell005: 4.0 s, frame 9); with no
    // cell selected, the frame of the response window where the cells that responded are brightest on average; and where none
    // responded, the cells that got a call. The way last used is kept; a new page opens on "Every Nth frame".
    const openFrames = async (way) => { await page.locator('#export-picture > summary').click(); await page.click('#save-frames'); if (way) await page.click('[data-frames=' + way + ']'); };
    const pairBoxes = async () => [await page.inputValue('#frames-baseline'), await page.inputValue('#frames-peak'), await page.isEnabled('#frames-go') && await page.isEnabled('#frames-download')];
    const brightest = (statuses) => page.evaluate((statuses) => { const a = window.calciumStep4, r = a.S.results, calls = r.callsByWindow[a.S.window], row = Object.values(calls)[0], cells = Object.keys(calls).filter((name) => r.traces[name] && statuses.includes(calls[name].Status)); let frame = 0, top = -Infinity;
      r.time.forEach((t, i) => { if (t < Number(row.Start_s) || t > Number(row.End_s)) return; const values = cells.map((name) => r.traces[name][i]).filter((v) => typeof v === 'number' && isFinite(v)), mean = values.reduce((sum, v) => sum + v, 0) / values.length; if (values.length && mean > top) { top = mean; frame = i + 1; } }); return frame; }, statuses);
    await openFrames('pair'); assert.deepEqual(await pairBoxes(), ['5', '13', true], 'with Cell001 selected: the last F₀ frame, and the frame of its peak'); await page.click('#frames-cancel');
    await openCell('Cell005'); await openFrames(); assert.ok(await page.isVisible('#frames-baseline'), 'the way last used is kept'); assert.deepEqual(await pairBoxes(), ['5', '9', true], 'with Cell005 selected: its own peak');
    // A frame typed over a suggestion is kept while the suggestion stays the same, so it need not be typed again; a new suggestion replaces it.
    await page.fill('#frames-peak', '11'); await page.click('#frames-cancel'); await openFrames(); assert.deepEqual(await pairBoxes(), ['5', '11', true], 'a typed frame is kept'); await page.keyboard.press('Escape'); assert.ok(await page.isHidden('#frames-options') && await app('cell') === 'Cell005', 'Esc closes the panel and does nothing else');
    await openCell(''); assert.equal(await app('cell'), null); await openFrames();
    assert.equal(await brightest(['Responded']), 13); assert.deepEqual(await pairBoxes(), ['5', '13', true], 'with no cell selected: where the cells that responded are brightest on average'); await page.click('#frames-cancel');
    const callsWere = await page.evaluate(() => { const a = window.calciumStep4, calls = a.S.results.callsByWindow[a.S.window], were = {}; for (const name of Object.keys(calls)) were[name] = calls[name].Status; calls.Cell001.Status = 'Excluded'; calls.Cell005.Status = 'Did not respond'; return were; });
    await openFrames(); assert.equal(await brightest(['Did not respond']), 9); assert.deepEqual(await pairBoxes(), ['5', '9', true], 'where no cell responded: the cells that got a call'); await page.click('#frames-cancel');
    await page.evaluate((were) => { const a = window.calciumStep4, calls = a.S.results.callsByWindow[a.S.window]; for (const name of Object.keys(were)) calls[name].Status = were[name]; }, callsWere);
    // A box that does not hold a frame of the recording keeps what was typed and is marked; beside it, the frames there are; and
    // nothing can be saved until it does.
    await openCell('Cell001'); await openFrames(); assert.deepEqual(await pairBoxes(), ['5', '13', true]);
    for (const typed of ['0', '31', '2.5', 'abc', '', '1e1', '0x10']) { await page.fill('#frames-peak', typed);
      assert.deepEqual([await page.inputValue('#frames-peak'), await page.getAttribute('#frames-peak', 'aria-invalid'), await page.isVisible('#frames-peak-range') && await text('#frames-peak-range'), await page.isDisabled('#frames-go'), await page.isDisabled('#frames-download'), await page.getAttribute('#frames-baseline', 'aria-invalid'), await page.isHidden('#frames-baseline-range')],
        [typed, 'true', '1 to 30', true, true, null, true], 'a peak frame of "' + typed + '" is kept as typed, marked, and not saved'); }
    await page.fill('#frames-baseline', '40'); assert.deepEqual([await page.getAttribute('#frames-baseline', 'aria-invalid'), await text('#frames-baseline-range')], ['true', '1 to 30']);
    await page.click('#frames-cancel'); await openFrames(); assert.deepEqual(await pairBoxes(), ['5', '13', true], 'what was typed and is no frame is not kept: the panel opens on its suggestions again'); await page.fill('#frames-baseline', '40'); await page.fill('#frames-peak', 'abc');
    await page.fill('#frames-baseline', '5'); await page.fill('#frames-peak', '13'); assert.deepEqual([await pairBoxes(), await page.getAttribute('#frames-peak', 'aria-invalid'), await page.isHidden('#frames-peak-range')], [['5', '13', true], null, true]);
    assert.ok(await page.evaluate(() => { const panel = document.getElementById('frames-options').getBoundingClientRect(), button = document.querySelector('#export-picture > summary').getBoundingClientRect(), tops = new Set(Array.from(document.querySelectorAll('#frames-options .actions button')).filter((b) => !b.hidden).map((b) => Math.round(b.getBoundingClientRect().top)));
      return (panel.bottom <= button.top + 1 || panel.top >= button.bottom - 1) && Math.abs(panel.right - button.right) < 24 && panel.left >= 0 && tops.size === 1; }), 'the panel sits against the Export button, its three buttons in one row');
    const pairSays = (await text('#frames-count')).match(/^One picture, (\d+) × (\d+) px, PNG, 300 dpi\.$/); assert.ok(pairSays, 'the panel says what will be saved: ' + await text('#frames-count'));
    // The picture: "baseline" and "peak" over two panels, each with its frame and time; each panel drawn as a saved frame is, from the
    // whole of its own frame; under both the key, once, with one scale bar at its right. Clear of the label a frame saved alone has on
    // its corner, the panels are, pixel for pixel, the frames 5 and 13 as Save frames saves them.
    await page.check('#number-cells'); const pairMade = await savedDrawing('pair', 5, 13, null); await page.uncheck('#number-cells');
    { const [left, right] = pairMade.images, said = pairMade.texts.map((t) => t.text), bar = barOf(pairMade);
      assert.deepEqual(pairMade.size, [Number(pairSays[1]), Number(pairSays[2])], 'as large as the panel said'); assert.equal(pairMade.images.length, 2);
      assert.ok(near(left.from, [0, 0, 250, 150], 1e-9) && near(right.from, [0, 0, 250, 150], 1e-9) && near([left.to[2], left.to[3], right.to[2], right.to[3]], [1000, 600, 1000, 600], 1e-6) && left.to[0] === 0 && right.to[0] > 1000 && right.to[0] + 1000 === pairMade.size[0] && left.to[1] === right.to[1] && left.to[1] > 0, 'two panels of a saved frame\'s size, side by side under their names: ' + JSON.stringify([left.to, right.to]));
      assert.ok([left.to[0], left.to[1], right.to[0]].every(Number.isInteger) && left.smoothed && right.smoothed && pairMade.cuts.length === 2, 'each on whole pixels of the file, and cut at its own edges');
      const over = (word, panel) => { const at = pairMade.texts.find((t) => t.text === word); return !!at && at.x >= panel.to[0] && at.x < panel.to[0] + 500 && at.y < panel.to[1]; };
      assert.ok(over('baseline', left) && over('Frame 5 · 2.0 s', left) && over('peak', right) && over('Frame 13 · 6.0 s', right), 'what each panel is, with its frame and time, over it: ' + said);
      for (const word of ['baseline', 'peak', 'as recorded', 'responded', 'did not respond', 'excluded or unscorable', 'Fixture recording']) assert.equal(said.filter((text) => text === word).length, 1, 'written once: ' + word);
      assert.ok(pairMade.texts.filter((t) => ['as recorded', 'responded', 'Fixture recording'].includes(t.text)).every((t) => t.y > left.to[1] + 600), 'the key is under the pair');
      assert.ok(bar.words.length === 1 && bar.bars.length === 1 && bar.words[0].text === '20 µm' && near([bar.bars[0].wide], [20 / 0.8 * 4], 1e-6) && bar.bars[0].x > right.to[0] + 500 && bar.words[0].y - bar.words[0].high / 2 >= right.to[1] + 600, 'one scale bar, at the right of the key under the pair: a fifth of the 200 µm field is 40, so 20 µm');
      assert.deepEqual(numbers(pairMade).map((t) => t.text), ['1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6'], 'the cells are drawn on both panels as the page draws them, here with their numbers');
      const panels = await page.evaluate(async ([leftAt, rightAt]) => { const a = window.calciumStep4, pair = await a.pairCanvas(5, 13, null), alone = [await a.viewCanvas(5, true), await a.viewCanvas(13, true)], band = [200, 600];   // (the rows of a frame under the label on its corner)
        const rows = (canvas, x, y) => canvas.getContext('2d').getImageData(x, y + band[0], 1000, band[1] - band[0]).data, differ = (one, other) => { let n = 0; for (let i = 0; i < one.length; i++) if (one[i] !== other[i]) n++; return n; };
        const mine = [rows(pair, leftAt[0], leftAt[1]), rows(pair, rightAt[0], rightAt[1])], theirs = alone.map((canvas) => rows(canvas, 0, 0)); let lit = 0; for (let i = 0; i < mine[1].length; i += 4) if (mine[1][i] > 200) lit++;
        return { baseline: differ(mine[0], theirs[0]), peak: differ(mine[1], theirs[1]), swapped: differ(mine[0], theirs[1]), lit, alone: alone.map((canvas) => [canvas.width, canvas.height]) }; }, [left.to, right.to]);
      assert.deepEqual([panels.baseline, panels.peak, panels.alone], [0, 0, [[1000, 600], [1000, 600]]], 'the panels are frames 5 and 13 as saved alone, pixel for pixel: ' + JSON.stringify(panels));
      assert.ok(panels.swapped > 10000 && panels.lit > 4000, 'and they are the right way round: the cell is lit in the peak panel only (a square of 80 by 80 pixels): ' + JSON.stringify(panels)); }
    // With "Area on screen" both panels show the one area, and the bar is the area's.
    { const part = { x: 20, y: 45, width: 50, height: 30, zoom: 5 }, made = await savedDrawing('pair', 5, 13, part), bar = barOf(made);
      assert.ok(made.images.length === 2 && made.images.every((image) => near(image.from, [20, 45, 50, 30], 1e-9) && !image.smoothed && near([image.to[2], image.to[3]], [1000, 600], 1e-6)) && bar.words.length === 1 && bar.words[0].text === '5 µm' && near([bar.bars[0].wide], [5 / 0.8 * 20], 1e-6), 'a pair of the area: ' + JSON.stringify([made.images, bar])); }
    // Saved from the panel: into the plots folder, as a picture is, under a name that says baseline and peak with their frames. The
    // record names both frames with their times, in place of the frame on screen; the view; and the gray levels the view's frames
    // are shown between, which are the same for both panels.
    assert.ok(await page.isHidden('#frames-options'), '(a click elsewhere on the page, as on Number cells above, closes the panel)'); await openFrames(); assert.deepEqual(await pairBoxes(), ['5', '13', true]);
    await page.click('#frames-go'); await page.waitForFunction(() => /^Saved \S+baseline\S+\.png/.test(document.getElementById('picture-saved').textContent), null, { polling: 100, timeout: 20000 }); assert.ok(await page.isHidden('#frames-options'));
    { const name = (await text('#picture-saved')).match(/^Saved (\S+)\.png/)[1], png = fs.readFileSync(findFile(out, name + '.png')), record = JSON.parse(fs.readFileSync(findFile(out, name + '.metadata.json'), 'utf8')), range = await page.evaluate(async () => (await (await fetch('api/range?mode=raw')).json()).ranges);
      assert.match(name, /^Fixture_recording_picture_as_recorded_baseline_frame5_peak_frame13_\d{8}-\d{6}$/, 'the file name says baseline and peak with their frames'); assert.deepEqual([png.readUInt32BE(16), png.readUInt32BE(20)], pairMade.size);
      assert.deepEqual([record.baseline, record.peak, record.view, 'frame' in record, 'time_s' in record], [{ frame: 5, time_s: 2 }, { frame: 13, time_s: 6 }, 'raw', false, false], 'the record names both frames: ' + JSON.stringify(record.baseline) + JSON.stringify(record.peak));
      assert.deepEqual(record.brightness_range, range); assert.ok(range.length === 1 && range[0].shown === 'gray levels' && range[0].from === 100 && range[0].to > 100, 'and the gray levels the frames are shown between: ' + JSON.stringify(range));
      assert.deepEqual([record.picture.area, record.picture.scale_bar_um, record.picture.outlines_shown, record.picture.cells_numbered, record.picture.format, record.picture.dpi, record.cell], ['whole picture', 20, true, false, 'png', 300, 'Cell001']); }
    // "Download a copy" hands the same picture and its record to the browser.
    { const got = [], both = new Promise((resolve, reject) => { const late = setTimeout(() => { page.off('download', listener); reject(new Error('Expected the pair and its record as two downloads')); }, 10000), listener = (download) => { got.push(download.suggestedFilename()); if (got.length === 2) { clearTimeout(late); page.off('download', listener); resolve(); } }; page.on('download', listener); });
      await openFrames(); await page.click('#frames-download'); await both; got.sort();
      assert.ok(/^Fixture_recording_picture_as_recorded_baseline_frame5_peak_frame13_\d{8}-\d{6}\.metadata\.json$/.test(got[0]) && got[1] === got[0].replace(/\.metadata\.json$/, '.png'), 'the pair downloads as a picture with its record: ' + got); }
    // "Every Nth frame" is put back, as a new page has it; the picture is as these checks found it.
    await openFrames('every'); await page.click('#frames-cancel'); assert.deepEqual([await app('zoom'), await app('cell'), await app('plot'), (await areaNow()).saved, await page.inputValue('#save-area')], [1, 'Cell001', 'cell', null, 'whole']); assert.equal(await page.evaluate(() => window.calciumStep4.S.frame), frameOnScreen);

    // Actual time determines heatmap widths/clicks, and the acquisition gap stays blank.
    await page.click('.tab[data-tab=heatmap]');
    assert.equal(await page.inputValue('#palette'), 'signed', 'the heatmap opens in blue–white–red');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4; return [f.S.palette, f.S.range, f.S.order, f.S.rowNames, f.S.hideExcluded, f.heatmapData(f.STYLE.heatmap).colourAt(0).join()]; }), ['signed', 'robust', 'call', false, false, '255,255,255'], 'as its lists and boxes stand, with zero drawn white');
    await page.selectOption('#order', 'orig'); await settledPaint();
    const heat = await page.evaluate(() => {
      const app = window.calciumStep4, canvas = document.getElementById('heatmap'), box = canvas.getBoundingClientRect(), g = app.Heat.geom;
      const time = app.S.results.time, x = (t) => g.L + (t - time[0]) / (time[time.length - 1] - time[0]) * g.W;
      const y = g.T + (4.5 / g.rows) * g.H, dpr = canvas.width / canvas.clientWidth;
      const gapPixel = Array.from(canvas.getContext('2d').getImageData(Math.round(x(9.75) * dpr), Math.round(y * dpr), 1, 1).data);
      return { x: box.left + x(15.5), y: box.top + y, gapPixel, sampleBounds: app.samplingSupport().bounds.slice(14, 16) };
    });
    assert.deepEqual(heat.gapPixel.slice(0, 3), [229, 229, 229], 'unrecorded time must not stretch a neighboring sample');
    assert.ok(heat.sampleBounds[0][1] < heat.sampleBounds[1][0]);
    await page.screenshot({ path: path.join(out, 'page-heatmap.png') });
    await page.mouse.click(heat.x, heat.y); await waitFrame(22);
    assert.equal(await page.inputValue('#cell'), 'Cell005'); assert.equal(await text('.tab.on'), 'Cells and traces');
    await page.click('.tab[data-tab=heatmap]'); await page.selectOption('#palette', 'plasma'); await page.selectOption('#range', 'full');
    await page.screenshot({ path: path.join(out, 'page-heatmap-plasma.png') });
    await page.selectOption('#palette', 'signed'); await page.selectOption('#range', 'robust');
    await page.locator('#export-heatmap > summary').click(); await page.click('#edit-heatmap');
    await page.waitForFunction(() => document.getElementById('figure-preview').width > 300);
    assert.equal(await page.textContent('#figure-name'), 'Heatmap figure'); assert.equal(await page.inputValue('#figure-palette'), 'signed', 'the editor shows the colours chosen on the page');
    assert.equal(await page.inputValue('#figure-title'), '', 'the heatmap figure has no title until one is typed either');
    await page.selectOption('#figure-palette', 'plasma'); assert.equal(await page.inputValue('#palette'), 'plasma', 'colours chosen in the editor are the page\'s too');
    const heatBefore = await page.evaluate(() => document.getElementById('figure-preview').toDataURL());
    await typeBox('barMax', '9'); await typeBox('width', '4'); await typeBox('height', '3');
    await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was && /^1200 × 900 pixels/.test(document.getElementById('figure-pixels').textContent), heatBefore);
    assert.deepEqual(await page.evaluate(() => { const c = window.calciumStep4.figureCanvas('heatmap', 150); return [c.width, c.height]; }), [600, 450]);
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.figureKept('heatmap')), { width: 4, height: 3 }, 'the colour scale typed for this recording is not carried to the next');
    await page.selectOption('#figure-save-dpi', '600'); await typeBox('width', '20'); await typeBox('height', '20');
    await page.waitForFunction(() => /too large to save/.test(document.getElementById('figure-pixels').textContent)); assert.ok(await page.isDisabled('#figure-save') && await page.isDisabled('#figure-download'), 'a figure too large to draw cannot be saved');
    // Reset to defaults puts back every box the editor shows: also the colours, which are the page's, and the resolution, which is every saved picture's.
    await page.click('#figure-reset');
    assert.deepEqual([await page.inputValue('#figure-palette'), await page.inputValue('#figure-save-dpi'), await page.inputValue('#figure-width'), await page.inputValue('#figure-barMax')], ['signed', '300', '7', '']); await page.click('#figure-close');
    assert.ok(await page.isHidden('#figure-editor') && await page.inputValue('#palette') === 'signed' && await page.inputValue('#save-dpi') === '300');

    // ---------------------------------------------------------------- the heatmap
    // One painter for the page and for a saved figure; one function for the cells and their order; one for the colour range; one key.
    // What a heatmap drawing puts on its canvas: every text and filled box, with where its parts went.
    // true: the saved figure at 96 pixels per inch; false: the page's own heatmap, drawn again.
    const heatDrawing = (figure) => page.evaluate((figure) => {
      const f = window.calciumStep4, P = CanvasRenderingContext2D.prototype, was = {}, seen = { texts: [], boxes: [], lines: [] }, paths = new WeakMap(); for (const name of ['fillText', 'fillRect', 'stroke', 'beginPath', 'moveTo', 'lineTo']) was[name] = P[name];
      const at = (ctx, x, y) => { const m = ctx.getTransform(); return [m.a * x + m.c * y + m.e, m.b * x + m.d * y + m.f]; }, wanted = (ctx) => figure || ctx.canvas === document.getElementById('heatmap');
      const reach = (ctx, x, y) => { const p = at(ctx, x, y), b = paths.get(ctx) || [p[0], p[1], p[0], p[1]]; paths.set(ctx, [Math.min(b[0], p[0]), Math.min(b[1], p[1]), Math.max(b[2], p[0]), Math.max(b[3], p[1])]); };
      P.beginPath = function () { paths.delete(this); return was.beginPath.apply(this, arguments); };
      P.moveTo = function (x, y) { reach(this, x, y); return was.moveTo.apply(this, arguments); }; P.lineTo = function (x, y) { reach(this, x, y); return was.lineTo.apply(this, arguments); };
      P.fillText = function (text, x, y) { if (wanted(this)) seen.texts.push(window.textSeen(this, text, x, y, at)); return was.fillText.apply(this, arguments); };
      P.fillRect = function (x, y, w, h) { const p = at(this, x, y); if (wanted(this)) seen.boxes.push({ x: p[0], y: p[1], w, h, ink: this.fillStyle }); return was.fillRect.apply(this, arguments); };
      P.stroke = function () { const b = paths.get(this); if (b && wanted(this)) seen.lines.push({ x0: b[0], y0: b[1], x1: b[2], y1: b[3] }); return was.stroke.apply(this, arguments); };   // (the corners of each stroked path)
      try { if (figure) f.figureCanvas('heatmap', 96); else window.dispatchEvent(new Event('resize')); } finally { Object.assign(P, was); }
      seen.texts = window.asWritten(seen.texts); const d = figure ? f.Fig.drawn : f.Heat.geom;
      return Object.assign(seen, { top: figure ? d.top : 0, L: d.L, T: d.T, W: d.W, H: d.H, rows: d.names, named: d.named, groups: d.groups, low: d.low, high: d.high, below: d.below, above: d.above, x0: d.x0, x1: d.x1,
        key: d.items.map((i) => i.text), keyAt: d.key ? [d.key.x, d.key.y, d.key.layout.height / d.key.layout.rowHeight] : null, numbers: d.bar ? d.bar.numbers.map((n) => n.text) : null, barX: d.bar ? d.bar.x : null, times: d.times, strip: d.strip, flag: d.flag });
    }, figure);
    const heatStyle = (changes) => page.evaluate((changes) => { const f = window.calciumStep4; Object.assign(f.Fig.heatmap, f.STYLE.heatmap, changes); }, changes || {});   // the built-in look, with these entries changed
    const heatRows = () => page.evaluate(() => window.calciumStep4.Heat.geom.names);
    await heatStyle();
    // The built-in look: 7 × 5 in; the picture fills the width, with a narrow colour bar at the right as tall as the picture.
    assert.deepEqual(await page.evaluate(() => { const s = window.calciumStep4.STYLE.heatmap, c = window.calciumStep4.figureCanvas('heatmap', 100); return [s.width, s.height, c.width, c.height]; }), [7, 5, 700, 500]);
    // Export stays at the right end of the first row of the heatmap's controls; where the row is too narrow, the choices wrap, not Export.
    for (const width of [1440, 1024]) { await page.setViewportSize({ width, height: 900 }); await settledPaint();
      assert.ok(await page.evaluate(() => { const box = (sel) => document.querySelector(sel).getBoundingClientRect(), exp = box('#export-heatmap > summary');
        return Math.abs(exp.top - box('#order').top) <= 2 && exp.left >= box('#range').right && exp.left >= box('#hide-excluded').right; }), width + ' px: Export stays at the end of the first row of the heatmap\'s controls'); }
    await page.setViewportSize({ width: 1440, height: 900 }); await settledPaint();
    await page.selectOption('#order', 'call'); let heat1 = await heatDrawing(true), heatPage = await heatDrawing(false);
    assert.ok(heat1.W / (7 * 96) > 0.7 && 7 * 96 - (heat1.L + heat1.W) < 90, 'the picture fills the figure: ' + [heat1.L, heat1.W]);
    const barBoxes = heat1.boxes.filter((b) => Math.abs(b.x - heat1.barX) < 0.01 && b.h <= 1.5);
    assert.ok(heat1.barX > heat1.L + heat1.W && barBoxes.length >= heat1.H - 1 && barBoxes[0].w <= 12 && Math.abs(barBoxes[0].y - (heat1.top + heat1.T)) < 1, 'the colour bar is narrow, right of the picture and as tall as it');
    // Rows sorted by call: the largest rise on top of each group; a cell below zero last in its group; cells without a trace last.
    assert.deepEqual(heat1.rows, ['Cell001', 'Cell005', 'Cell002', 'Cell006', 'Cell003', 'Cell004']);
    assert.deepEqual(heatPage.rows, heat1.rows, 'the page and the figure have the same rows');
    { // One rule numbers the time axis of the time course and of the heatmap, on the page and saved: the same round times, no more than about six.
      const plotTimes = await page.evaluate(() => { const f = window.calciumStep4, g = f.Plot.geom; f.figureCanvas('plot', 96); return { x: [g.x0, g.x1], page: g.times, figure: f.Fig.drawn.times }; });
      assert.deepEqual(plotTimes.x, [heatPage.x0, heatPage.x1], '(the plot shows the whole recording)');
      assert.deepEqual([heatPage.times, heat1.times, plotTimes.figure], [plotTimes.page, plotTimes.page, plotTimes.page], 'the heatmap is numbered along time as the time course is: ' + JSON.stringify([heatPage.times, heat1.times, plotTimes]));
      assert.ok(plotTimes.page.length >= 3 && plotTimes.page.length <= 8, 'about six numbers: ' + plotTimes.page); }
    // (Cell005 is below zero; with the larger rise it would lead its group, and it still comes last in it)
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, row = f.S.results.callsByWindow[f.S.window].Cell005, was = row.Peak_increase; row.Peak_increase = '2';
      try { window.dispatchEvent(new Event('resize')); return f.Heat.geom.names.slice(0, 2); } finally { row.Peak_increase = was; window.dispatchEvent(new Event('resize')); } }), ['Cell001', 'Cell005'], 'a cell below zero goes last within its group');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, row = f.S.results.callsByWindow[f.S.window].Cell001, was = row.Peak_increase; row.Peak_increase = '0.5'; f.S.withoutBelowZero = false;
      const T = f.S.results.traces, kept = T.Cell005; T.Cell005 = kept.map((v) => Math.max(v, 0));   // no longer below zero: the larger rise is on top
      try { window.dispatchEvent(new Event('resize')); return f.Heat.geom.names.slice(0, 2); } finally { row.Peak_increase = was; T.Cell005 = kept; window.dispatchEvent(new Event('resize')); } }), ['Cell005', 'Cell001'], 'the largest rise in the response window is on top');
    // "Onset time": the responders in the order of their onsets (Cell005 at 4 s, Cell001 at 4.5 s), the rest as before.
    assert.deepEqual(await page.locator('#order option').evaluateAll((options) => options.map((o) => o.value)), ['call', 'onset', 'tpeak', 'amp', 'orig']);
    await page.selectOption('#order', 'onset'); assert.deepEqual(await heatRows(), ['Cell005', 'Cell001', 'Cell002', 'Cell006', 'Cell003', 'Cell004']);
    await page.selectOption('#order', 'amp'); assert.deepEqual(await heatRows(), ['Cell001', 'Cell006', 'Cell002', 'Cell005', 'Cell003', 'Cell004'], 'by peak amplitude, too, the cell below zero comes after the others');
    await page.selectOption('#order', 'call');
    // Groups: the rows of one call together. Nothing is written or drawn beside the rows but the axis title, "Cells": no brackets, no
    // counts, no row counter that reads like cell numbers. The strip shows the groups, and the key says how many cells each holds.
    assert.deepEqual(heat1.groups, [{ kind: 'responded', first: 0, n: 2 }, { kind: 'not', first: 2, n: 1 }, { kind: 'other', first: 3, n: 3 }]);
    // The strip: one solid bar for each group of rows, with no gaps between its rows; the mark of the cells below zero likewise; a group
    // of one cell is never shorter than stripLeast, centred on its row, so it can be found. (the boxes filled in the strip's and the mark's columns)
    const stripBars = (seen, at) => seen.boxes.filter((b) => Math.abs(b.x - (seen.L - at)) < 0.01).map((b) => ({ y: b.y - seen.top - seen.T, h: b.h, ink: b.ink }));
    for (const seen of [heat1, heatPage]) { const rowH = seen.H / seen.rows.length, bars = stripBars(seen, seen.strip).sort((a, b) => a.y - b.y);
      assert.deepEqual(bars.map((b) => [Math.round(b.y / rowH), Math.round(b.h / rowH)]), [[0, 2], [2, 1], [3, 3]], 'one bar for each group: ' + JSON.stringify(bars));
      assert.deepEqual(stripBars(seen, seen.flag).map((b) => [Math.round(b.y / rowH), Math.round(b.h / rowH)]), [[1, 1]], 'and one for the run of cells below zero'); }
    await heatStyle({ stripLeast: 60 }); { const seen = await heatDrawing(true), rowH = seen.H / seen.rows.length, bars = stripBars(seen, seen.strip), one = bars.find((b) => Math.abs(b.h - 60 * 96 / 72) < 0.01);
      assert.ok(rowH < 60 * 96 / 72 && one && Math.abs(one.y + one.h / 2 - 2.5 * rowH) < 0.01, 'a group of one row is drawn at least stripLeast high, centred on its row: ' + JSON.stringify([rowH, bars]));
      assert.equal(bars.indexOf(one), bars.length - 1, 'and on top of the longer runs beside it'); assert.ok(bars.every((b) => b.y >= -0.01 && b.y + b.h <= seen.H + 0.01), 'inside the picture\'s height'); }
    await heatStyle();
    const besideRows = (seen, t) => t.x < seen.L && t.y > seen.top + seen.T && t.y < seen.top + seen.T + seen.H, beside = (seen) => seen.texts.filter((t) => besideRows(seen, t)).sort((a, b) => a.y - b.y).map((t) => t.text);
    const besideLines = (seen) => seen.lines.filter((l) => l.x1 < seen.L && l.y1 > seen.top + seen.T && l.y0 < seen.top + seen.T + seen.H);
    for (const seen of [heat1, heatPage]) { assert.deepEqual(beside(seen), ['Cells'], 'left of the picture: the axis title and nothing else'); assert.deepEqual(besideLines(seen), [], 'and no line is drawn there: no brackets'); }
    assert.deepEqual(heat1.named, [], 'no cell names until they are asked for');
    // Every word of the heatmap: its key, "stimulus", the three titles, and numbers (the time axis, and the colour bar with its ≤ and ≥).
    const heatWords = (seen) => seen.texts.map((t) => t.text).filter((text) => !/^([≤≥] )?−?\d+(\.\d+)?$/.test(text) && !seen.key.includes(text) && !seen.named.includes(text)).sort();
    assert.deepEqual(heatWords(heat1), ['Cells', 'Time (s)', 'stimulus', 'ΔF/F₀']); assert.deepEqual(heatWords(heatPage), ['Cells', 'Time (s)', 'stimulus', 'ΔF/F₀']);
    await page.selectOption('#order', 'orig'); assert.deepEqual((await heatDrawing(true)).groups, [], 'in the original order the calls are mixed: no groups'); await page.selectOption('#order', 'call');
    // The key, in the figure's ink: the calls that are drawn, each with its number of cells; the mark of a cell below zero; the two bars above the picture. The stimulus is labelled above the picture.
    assert.deepEqual(heat1.key, ['responded (n = 2)', 'did not respond (n = 1)', 'excluded or unscorable (n = 3)', 'below zero', 'baseline', 'response window']); assert.deepEqual(heatPage.key, heat1.key);
    assert.ok(heat1.texts.filter((t) => heat1.key.includes(t.text)).every((t) => t.ink === '#000000') && heatPage.texts.filter((t) => heatPage.key.includes(t.text)).every((t) => t.ink === '#3e4c59'), 'black in a figure and in the page\'s ink on the page');
    const keyTexts = (seen) => seen.texts.filter((t) => seen.key.includes(t.text));
    assert.ok(heatPage.keyAt[2] === 1 && heatPage.keyAt[0] < heatPage.L && keyTexts(heatPage).every((t) => t.y < heatPage.T), 'the key is one row above the picture, from the left edge: ' + heatPage.keyAt);
    assert.ok(Number.isInteger(heat1.keyAt[2]) && new Set(keyTexts(heat1).map((t) => Math.round(t.y))).size === heat1.keyAt[2] && keyTexts(heat1).every((t) => t.y < heat1.top + heat1.T && t.x < 7 * 96), 'in a figure too narrow for one row it wraps by whole entries: ' + heat1.keyAt);
    assert.ok(heat1.texts.some((t) => t.text === 'stimulus' && t.y < heat1.top + heat1.T) && heat1.boxes.some((b) => b.ink === '#e9a84c' && b.w > 100 && b.y + b.h <= heat1.top + heat1.T), 'the stimulus label and the response window are above the picture');
    assert.equal(await page.locator('#heatmap-legend').count(), 0, 'the key is drawn in the heatmap; there is no second one above it');
    assert.match(await page.getAttribute('#heatmap', 'aria-label'), /one row for each of 6 cells.*Colours from ≤ 0\.0 to 1\.0\. responded \(n = 2\); did not respond \(n = 1\); excluded or unscorable \(n = 3\); below zero; baseline; response window$/, 'and a screen reader is told what it says');
    // The colour range comes from the cells that are neither excluded nor below zero, with either choice under Range.
    assert.deepEqual([heat1.low, heat1.high], [0, 1]); assert.deepEqual([heatPage.low, heatPage.high], [0, 1]);
    assert.deepEqual(await page.locator('#range option').allTextContents(), ['1st to 99th percentile', 'Smallest to largest']);
    await page.selectOption('#range', 'full'); assert.deepEqual(await page.evaluate(() => { const g = window.calciumStep4.Heat.geom; return [g.low, g.high, g.below]; }), [0, 1, true], 'the −1.5 of the cell below zero is not the smallest value the range is taken from');
    await page.selectOption('#range', 'robust');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, T = f.S.results.traces, kept = T.Cell005; T.Cell005 = kept.map(() => -3);   // off for the whole recording
      try { const d = f.heatmapData(f.STYLE.heatmap); return [d.low, d.high, d.below]; } finally { T.Cell005 = kept; } }), [0, 1, true], 'a cell that is below zero throughout does not set the range either');
    // Blue–white–red: white at zero, the full blue at the lower limit, the full red at the upper, each side on its own scale.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, d = f.heatmapData(Object.assign({}, f.STYLE.heatmap, { barMin: -0.5, barMax: 2 })); return [0, -0.5, 2, -9, 9].map((v) => d.colourAt(v).join()); }),
      ['255,255,255', '59,76,192', '180,4,38', '59,76,192', '180,4,38']);
    assert.ok(await page.evaluate(() => { const f = window.calciumStep4, d = f.heatmapData(Object.assign({}, f.STYLE.heatmap, { barMin: -0.5, barMax: 2 })), blue = d.colourAt(-0.25), red = d.colourAt(1);
      return Math.abs(blue[0] - 157) <= 2 && Math.abs(red[1] - 130) <= 2; }), 'half of the lower limit is half way to blue, and half of the upper limit half way to red');
    // The colour bar: both ends are always written, and an end says ≤ or ≥ whenever a value that is drawn lies beyond it.
    assert.deepEqual(heat1.numbers.slice(0, 2), ['1.0', '≤ 0.0']);
    const barEnds = (seen) => seen.numbers.slice(0, 2).map((text) => { const t = seen.texts.find((t) => t.text === text && t.x > seen.barX); return t ? Math.round(t.y - seen.top - seen.T) : null; });
    assert.deepEqual(barEnds(heat1), [0, Math.round(heat1.H)], 'the two numbers sit at the two ends of the bar');
    await heatStyle({ barMax: 0.5 }); assert.deepEqual((await heatDrawing(true)).numbers.slice(0, 2), ['≥ 0.5', '≤ 0.0'], '"Scale to" keeps the marks');
    await heatStyle({ barMin: -0.37, barMax: 0.83 }); heat1 = await heatDrawing(true);
    assert.deepEqual([heat1.low, heat1.high], [-0.37, 0.83]); assert.deepEqual(heat1.numbers.slice(0, 2), ['≥ 0.83', '≤ −0.37'], 'limits that are not round numbers are written in full, with their marks'); assert.deepEqual(barEnds(heat1), [0, Math.round(heat1.H)]);
    await heatStyle({ barMin: -2, barMax: 3 }); assert.deepEqual((await heatDrawing(true)).numbers.slice(0, 2), ['3', '−2'], 'no mark at an end that nothing lies beyond');
    // The numbers of the colour bar are a full line apart: a round number within a line of an end is not written, on the page and in a
    // figure alike (as Run1's "≥ 0.12" against "0.1" was). (the gaps between neighbouring numbers, in heights of their text)
    const barGaps = (seen) => { const bar = seen.texts.filter((t) => t.x > seen.barX && t.y >= seen.top + seen.T - 1 && t.y <= seen.top + seen.T + seen.H + 1).sort((a, b) => a.y - b.y);
      return bar.slice(1).map((t, i) => (t.y - bar[i].y) / (t.size * Math.abs(Number(t.m.split(',')[3])))); };
    for (const [style, where] of [[{ barMin: -0.6, barMax: 0.12, height: 3 }, 'a figure with ends close to round numbers'], [{ barMin: -0.6, barMax: 0.12 }, 'a 7 × 5 in figure'], [null, 'the page']]) {
      let seen; if (style) { await heatStyle(style); seen = await heatDrawing(true); } else { await heatStyle(); seen = await heatDrawing(false); }
      const gaps = barGaps(seen); assert.ok(gaps.length >= 1 && gaps.every((gap) => gap >= 2 - 1e-6), where + ': the colour bar\'s numbers are a line apart: ' + JSON.stringify([seen.numbers, gaps]));
      if (style) assert.deepEqual(seen.numbers.slice(0, 2), ['≥ 0.12', '≤ −0.6'], 'and both ends are written'); }
    await heatStyle({ barMin: 2, barMax: 1 }); heat1 = await heatDrawing(true); assert.deepEqual([heat1.low, heat1.high], [0, 1], 'a minimum above the maximum is not a range: the automatic one stays');
    await heatStyle(); await page.selectOption('#palette', 'plasma'); heat1 = await heatDrawing(true);
    assert.deepEqual([heat1.low, heat1.numbers[1]], [0, '≤ 0.0'], 'a palette that starts at zero says in the figure, too, that values below zero take its first colour'); await page.selectOption('#palette', 'signed');
    // Time-axis limits and tick spacing from the editor.
    await heatStyle({ xMin: 4, xMax: 8, xStep: 2 }); heat1 = await heatDrawing(true);
    assert.deepEqual([heat1.x0, heat1.x1], [4, 8]);
    assert.deepEqual(heat1.texts.filter((t) => t.y > heat1.top + heat1.T + heat1.H && /^\d+$/.test(t.text)).map((t) => t.text), ['4', '6', '8']);
    assert.ok(!heat1.key.includes('baseline') && heat1.key.includes('response window') && !heat1.texts.some((t) => t.text === 'stimulus'), 'the baseline and the stimulus, before 4 s, are not shown or listed: ' + heat1.key);
    assert.deepEqual([heat1.below, heat1.numbers[1]], [false, '0.0'], 'and the −1.5 at the last frame, outside the time shown, puts no ≤ on the bar');
    await heatStyle();
    // Cell names: the real names, every row while the rows are tall enough, every Nth when crowded; the figure carries them.
    await page.check('#row-names'); heatPage = await heatDrawing(false); heat1 = await heatDrawing(true);
    assert.deepEqual(heatPage.named, heatPage.rows); assert.deepEqual(heat1.named, heat1.rows, 'the saved figure carries the names');
    assert.deepEqual(beside(heat1).filter((text) => text !== 'Cells'), heat1.rows, 'written beside their rows, in the order of the rows, and nothing else with them');
    assert.equal(besideLines(heat1).length, 1, 'each with a tick to its row'); assert.deepEqual(heatWords(heat1), ['Cells', 'Time (s)', 'stimulus', 'ΔF/F₀']);
    assert.ok(await page.evaluate(() => { const c = document.getElementById('heatmap'), box = document.getElementById('heatmap-box'); return c.clientHeight <= box.clientHeight && box.scrollHeight <= box.clientHeight + 1; }), 'the heatmap stays in view with its time axis: nothing to scroll');
    await heatStyle({ height: 1.65 }); heat1 = await heatDrawing(true);   // rows thinner than the 12-pixel text
    const everyNth = Math.ceil(12 * 1.1 / (heat1.H / heat1.rows.length));
    assert.ok(everyNth > 1 && everyNth < heat1.rows.length, 'the rows of this figure are crowded: ' + heat1.H);
    assert.deepEqual(heat1.named, heat1.rows.filter((name, c) => c % everyNth === 0), 'crowded rows are named every Nth');
    await heatStyle();
    // Hide excluded cells: their rows go, and so does their word in the key.
    assert.ok(!(await page.isChecked('#hide-excluded')), 'off until it is ticked'); await page.check('#hide-excluded'); heat1 = await heatDrawing(true);
    assert.deepEqual(heat1.rows, ['Cell001', 'Cell005', 'Cell002', 'Cell006']); assert.deepEqual(heat1.groups[2], { kind: 'other', first: 3, n: 1 }); assert.equal(heat1.key[2], 'unscorable (n = 1)'); assert.ok(beside(heat1).includes('Cells'));
    // The file kept beside a saved heatmap says which cell each row is and which colour range was used.
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.pictureRecord({ figure: 'heatmap' }).heatmap), { rows: ['Cell001', 'Cell005', 'Cell002', 'Cell006'], colour_range: [0, 1], cell_names: true, excluded_hidden: true });
    // The editor has the same boxes, and they are the page's own: ticking there ticks here.
    await page.locator('#export-heatmap > summary').click(); await page.click('#edit-heatmap'); await page.waitForFunction(() => !window.calciumStep4.Fig.queued && document.getElementById('figure-preview').width > 300);
    assert.deepEqual(await page.evaluate(() => ['xMin', 'xMax', 'xStep', 'barMin', 'barMax', 'row-names', 'hide-excluded', 'order'].map((id) => { const box = document.getElementById('figure-' + id); return box ? box.type : null; })),
      ['text', 'text', 'text', 'text', 'text', 'checkbox', 'checkbox', 'select-one']);
    assert.ok(await page.isChecked('#figure-row-names') && await page.isChecked('#figure-hide-excluded'), 'they show what is ticked on the page');
    assert.equal(await page.inputValue('#figure-yTitle'), 'Cells'); assert.deepEqual([await page.getAttribute('#figure-barMin', 'placeholder'), await page.getAttribute('#figure-xMax', 'placeholder')], ['auto', 'auto']);
    await page.uncheck('#figure-hide-excluded'); await page.uncheck('#figure-row-names');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4; return [f.S.hideExcluded, f.S.rowNames, document.getElementById('hide-excluded').checked, document.getElementById('row-names').checked, f.Heat.geom.names.length]; }), [false, false, false, false, 6]);
    await typeBox('barMin', '-0.5'); await typeBox('barMax', '2'); await typeBox('xMin', '2'); await typeBox('xStep', '3');
    await page.waitForFunction(() => { const d = window.calciumStep4.Fig.drawn; return !window.calciumStep4.Fig.queued && d.kind === 'heatmap' && d.low === -0.5 && d.high === 2 && d.x0 === 2; });
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.figureKept('heatmap')), {}, 'the time limits and the colour limits typed for this recording are not carried to the next');
    await page.click('#figure-reset'); await page.click('#figure-close'); await heatStyle();
    // The tables: the saved columns named in words with their units, numbers rounded by their unit (time 2 decimals, ΔF/F₀ 3,
    // percentages 1 with "%", gray levels 2), one dash for a missing value, the saved name and value when pointing at them.
    // The saved files keep their names, values and every column.
    const tableNow = (id) => page.evaluate((id) => { const t = document.getElementById(id), cells = (row) => [...row.cells];
      return { head: cells(t.tHead.rows[0]).map((c) => c.textContent), titles: cells(t.tHead.rows[0]).map((c) => c.title), rows: [...t.tBodies[0].rows].map((r) => cells(r).map((c) => c.textContent)),
        tips: [...t.tBodies[0].rows].map((r) => cells(r).map((c) => c.title)), frozen: cells(t.tHead.rows[0]).map((c) => c.classList.contains('cell-column')) }; }, id);
    const at = (t, first, column, also) => { const r = t.rows.findIndex((row) => row[0] === first && (!also || row.includes(also))); return { text: t.rows[r][t.head.indexOf(column)], tip: t.tips[r][t.head.indexOf(column)] }; };
    assert.deepEqual(await page.evaluate(() => { const d = window.calciumStep4.decimal; return [d(-5.186e-17, 3), d(-0.0004, 3), d(-1.5, 2), d(0.999895989894867, 2), d(4.99948, 2), d(200 / 3, 1, true), d(NaN, 2)]; }),
      ['0.000', '0.000', '−1.50', '1.00', '5.00', '66.7%', '–'], 'fixed decimals, the true minus, no "−0.000" (a saved −5e−17), a dash for nothing');
    await page.click('.tab[data-tab=calls]'); assert.equal(await page.locator('#table-calls tbody tr').count(), 12);
    let calls = await tableNow('table-calls');
    assert.deepEqual(calls.head, ['Cell', 'Role', 'Stimulus', 'Call', 'Reason', 'Required rise', 'Required rise set by', 'Largest rise', 'Largest rise at, s', 'Onset, s', 'Baseline SD', 'Flags'], 'the short set, in words, with units');
    assert.deepEqual([at(calls, 'Cell001', 'Required rise', 'Stimulus').text, at(calls, 'Cell001', 'Largest rise', 'Stimulus').text, at(calls, 'Cell001', 'Largest rise at, s', 'Stimulus').text, at(calls, 'Cell001', 'Onset, s', 'Stimulus').text, at(calls, 'Cell001', 'Flags', 'Stimulus').text],
      ['10.0%', '100.0%', '6.00', '4.50', 'Duration under a frame'], 'rises in % as in the verdict, times with 2 decimals, flags in words');
    assert.equal(at(calls, 'Cell001', 'Largest rise', 'Stimulus').tip, '1', 'pointing at a number shows the saved value');
    assert.deepEqual(at(calls, 'Cell001', 'Required rise set by', 'Stimulus'), { text: 'Minimum rise', tip: 'Saved: Floor' }, 'what set the required rise, in the words of the verdict');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, data = { header: ['Cell', 'Threshold_source', 'Noise_source'],
      rows: [['Cell001', 'Baseline noise', 'Whole-trace noise'], ['Cell002', 'Baseline noise', 'Noise floor'], ['Cell003', 'Baseline noise', 'Baseline SD'], ['Cell004', 'Floor', 'Whole-trace noise'], ['Cell005', 'Baseline noise', '']] };
      const column = f.tableColumns('calls', data, true).find((c) => c.saved === 'Threshold_source');
      const plain = f.tableColumns('calls', { header: ['Cell', 'Threshold_source'], rows: [['Cell001', 'Baseline noise']] }, true).find((c) => c.saved === 'Threshold_source');
      return data.rows.map((row) => { const s = f.shownCell(column, row, () => 'Responded'); return [s.text, s.title || '']; }).concat([[f.shownCell(plain, ['Cell001', 'Baseline noise'], () => 'Responded').text, '']]); }),
      [['Whole-trace noise', 'Saved: Baseline noise · Noise_source: Whole-trace noise'], ['Noise floor', 'Saved: Baseline noise · Noise_source: Noise floor'], ['Baseline SD', 'Saved: Baseline noise · Noise_source: Baseline SD'],
        ['Minimum rise', 'Saved: Floor'], ['Baseline noise', ''], ['Baseline noise', '']], 'the noise that set the required rise is named as the Fiji window names it; a table without Noise_source says Baseline noise');
    assert.ok(await page.locator('#table-calls').evaluate((t) => { const even = t.tBodies[0].rows[1], style = (e) => getComputedStyle(e).backgroundColor;
      return style(even.cells[0]) === style(even) && style(even.cells[0]) !== style(t.tBodies[0].rows[0].cells[0]); }, undefined, { timeout: 10000 }), 'the column that stays put takes the stripe of its row');
    assert.deepEqual(at(calls, 'Cell002', 'Largest rise at, s', 'Stimulus'), { text: '–', tip: 'Saved: 3 · Did not respond' }, 'the timing of a cell that did not respond is a dash, its saved value in the tooltip');
    assert.equal(at(calls, 'Cell002', 'Largest rise', 'Stimulus').text, '0.0%', 'its largest rise, the number compared with the threshold, is shown');
    assert.equal(at(calls, 'Cell003', 'Required rise', 'Stimulus').text, '–', 'one dash for a value that was not saved');
    assert.ok(calls.head.every((h) => !/_/.test(h)) && calls.rows.every((row) => ['Responded', 'Did not respond', 'Excluded', 'Unscorable'].includes(row[3])), 'no code names on screen; one set of words for the call');
    assert.match(calls.titles[calls.head.indexOf('Baseline SD')], /response baseline.*Not the SD of the F₀ frames.*Saved as Baseline_SD, a fraction/, 'the tooltip knows the table: this Baseline SD is not the one of the metrics');
    assert.match(calls.titles[calls.head.indexOf('Largest rise at, s')], /comes before the onset/, 'the time of the largest rise and the onset are told apart');
    assert.ok(!calls.head.includes('Run start, s'), 'a column that is the same as another by definition is shown once');
    await page.check('#all-call-columns'); calls = await tableNow('table-calls');
    assert.ok(calls.head.length > 12 && calls.head.includes('Run start, s') && calls.head.includes('Minimum rise threshold'), 'every saved column, the twins too');
    assert.equal(at(calls, 'Cell001', 'Window cut short', 'Stimulus').text, 'no', 'TRUE and FALSE read yes and no');
    assert.match(calls.titles[calls.head.indexOf('Run start, s')], /the same as Onset_s by definition/);
    await page.locator('#table-calls').evaluate((e) => e.parentElement.scrollLeft = 800);
    assert.ok(await page.locator('#table-calls tbody .cell-column').first().evaluate((e) => {const a=e.getBoundingClientRect(),b=e.closest('.scroller').getBoundingClientRect();return a.left>=b.left-1&&a.left<b.left+5;}));
    await page.uncheck('#all-call-columns');
    // A click on a header sorts by it; cells with no value go last both ways.
    const sortedBy = async (id, column) => { const before = await tableNow(id), k = before.head.indexOf(column); await page.locator('#' + id + ' thead th').nth(k).click(); return (await tableNow(id)).rows.map((row) => row[k]); };
    assert.deepEqual(await sortedBy('table-calls', 'Onset, s'), ['4.00', '4.50'].concat(Array(10).fill('–')), 'rising, the cells with no onset last');
    assert.deepEqual(await sortedBy('table-calls', 'Onset, s'), ['4.50', '4.00'].concat(Array(10).fill('–')), 'falling, and those cells still last');
    await page.screenshot({ path: path.join(out, 'page-calls.png') });
    // The calls tab has its own Export: the saved table as it is, with the record of how the page was set beside it.
    await page.locator('#export-calls > summary').click(); await page.click('#calls-save');
    await page.waitForFunction(() => /^Saved in plots: /.test(document.getElementById('calls-saved').textContent));
    { const saved = (await text('#calls-saved')).match(/^Saved in plots: (\S+\.csv)/)[1], file = fs.readFileSync(path.join(out, 'plots', saved), 'utf8');
      assert.match(saved, /^Fixture_recording_responder_calls_\d{8}-\d{6}\.csv$/);
      assert.deepEqual(file.trim().split('\n')[0], (await page.evaluate(() => window.calciumStep4.S.results.tables.calls.header.join(','))), 'every saved column under its saved name');
      assert.ok(fs.existsSync(path.join(out, 'plots', saved.replace(/\.csv$/, '.metadata.json')))); }
    await page.click('.tab[data-tab=counts]');
    assert.deepEqual(await page.evaluate(() => { const t = document.getElementById('table-tally'); return { head: [...t.tHead.rows[0].cells].slice(0, 4).map((c) => c.textContent), first: [...t.tBodies[0].rows[0].cells].slice(0, 4).map((c) => c.textContent), rows: t.tBodies[0].rows.length }; }),
      { head: ['Cells', '% responding', 'Responded', 'Evaluable'], first: ['All cells', '66.7%', '2', '3'], rows: 2 }, 'the same tally on the Responder counts tab, its share as in the headline');
    // The counts have one row for each window: shown turned, a row for each column, so the 27 saved columns need no sideways scroll.
    let counts = await tableNow('table-counts');
    assert.deepEqual(counts.head, ['', 'Capsaicin', 'Ionomycin (control)']);
    assert.deepEqual(counts.rows.find((row) => row[0] === '% responding'), ['% responding', '66.7%', '33.3%'], 'the saved Percent_responding reads as the headline does');
    assert.deepEqual(counts.rows.find((row) => row[0] === 'Minimum rise'), ['Minimum rise', '10.0%', '10.0%']);
    assert.ok(counts.frozen[0] && !counts.rows.some((row) => row[0] === 'Recording'), 'the names stay in view; the short set leaves out the recording');
    await page.check('#all-count-columns'); counts = await tableNow('table-counts');
    assert.ok(counts.rows.some((row) => row[0] === 'Gap rule') && counts.rows.length === 25, 'every saved column but the stimulus, which heads the table, and the recording, the same in both');
    await page.uncheck('#all-count-columns');
    await page.locator('#export-counts > summary').click(); await page.click('#counts-save');
    await page.waitForFunction(() => /^Saved in plots: \S+responder_counts_\S+\.csv/.test(document.getElementById('counts-saved').textContent));
    await page.click('.tab[data-tab=metrics]');
    // The summary: each metric once, under a name of its own, with its unit and where it is measured; columns that are the same by
    // definition share a row (pointing at the name gives the saved columns); n has a column of its own; the timing and shape of a
    // response are given only for cells that responded; All cells leaves out the excluded cells.
    const summary = await tableNow('table-metrics');
    assert.deepEqual(summary.head, ['Metric', 'Unit', 'Measured', 'All cells', 'n', 'Responded', 'n', 'Did not respond', 'n'], 'no columns for how cells were chosen when they were all chosen one way');
    assert.ok(summary.frozen[0], 'the metric names stay in view'); assert.equal(new Set(summary.rows.map((row) => row[0])).size, summary.rows.length, 'no two rows share a name');
    const rowOf = (name) => summary.rows.find((row) => row[0] === name), tipOf = (name) => summary.tips[summary.rows.indexOf(rowOf(name))][0];
    assert.deepEqual(rowOf('Peak'), ['Peak', 'ΔF/F₀', 'whole recording', '0.750 ± 0.250', '4', '1.000 ± 0.000', '2', '0.000', '1'], 'mean ± SEM with the same decimals, n beside it');
    assert.equal(tipOf('Peak'), 'Saved as Peak_dFF0, Response_Amplitude', 'twins by definition share a row');
    assert.deepEqual(rowOf('Largest rise').slice(0, 5), ['Largest rise', '%', 'response window, of the local level', '66.7 ± 33.3%', '3'], 'from the calls when there are no window metrics');
    assert.deepEqual(rowOf('F₀').slice(3, 5), ['10.00 ± 0.00', '4'], 'gray levels with 2 decimals; the two excluded cells are not in All cells');
    assert.deepEqual(rowOf('Onset').slice(3), ['', '', '4.25 ± 0.25', '2', '', ''], 'the onset is the timing of a response: given once, under Responded'); assert.equal(tipOf('Onset'), 'Saved as Onset_s, Run_start_s');
    assert.ok(await page.locator('#table-metrics').evaluate((t) => [...t.tBodies[0].rows].every((r) => [...r.cells].slice(3).every((c) => c.classList.contains('num')))), 'a dash stands right like the numbers of its column');
    assert.equal(await page.locator('#metrics-table').evaluate((s) => getComputedStyle(s).maxWidth), 'none', 'the list of tables is not cut off ("Summary: mean ± SEM of each metri…")');
    for (const name of ['Time to peak', 'Time to 50% of peak', 'Rise time', 'FWHM']) { assert.match(rowOf(name)[2], /^whole recording/, name + ' is named for the whole recording'); assert.deepEqual([rowOf(name)[3], rowOf(name)[7]], ['', ''], name + ' is given once, not again under All cells or for cells that did not respond'); }
    assert.ok(!summary.rows.some((row) => row[0] === 'Window AUC'), 'a column this analysis did not save is not a row on screen');
    assert.match(await text('#metric-help'), /cells that responded/); assert.ok(await page.isHidden('#all-metric-columns-label'));
    const long = await page.evaluate(() => window.calciumStep4.metricSummary().long), groupsOf = (await page.evaluate(() => window.calciumStep4.S.results.chosenBlind)) ? 5 : 3;
    assert.deepEqual(long[0], ['Metric', 'Unit', 'Measured', 'Saved_as', 'Group', 'Mean', 'SD', 'SEM', 'N_cells'], 'the saved summary keeps its columns');
    assert.equal(long.length - 1, (await page.evaluate(() => window.calciumStep4.SUMMARY_ROWS.length)) * groupsOf, 'and has the same rows for every recording');
    assert.deepEqual(long.find((row) => row[0] === 'Window area under the curve' && row[4] === 'All cells').slice(1), ['ΔF/F₀ × s', 'response window, signed', '', 'All cells', '', '', '', 0], 'a metric not saved has no cells, not no row; one unit for both areas');
    { const rise = long.find((row) => row[0] === 'Peak rise' && row[4] === 'All cells');
      assert.deepEqual(rise.slice(0, 3).concat(rise.slice(5, 6)), ['Peak rise', 'fraction of baseline', 'response window', 2 / 3], 'in the saved unit and words, at full precision: a fraction stays a fraction'); }
    assert.deepEqual(long.find((row) => row[0] === 'Signal-to-noise' && row[4] === 'All cells').slice(0, 3), ['Signal-to-noise', 'peak ÷ baseline SD', 'whole recording'], 'the saved summary keeps its words');
    assert.ok(long.some((row) => row[0] === 'Raw baseline SD') && long.some((row) => row[0] === 'Area under the curve'), 'two names made unique');
    assert.equal(new Set(long.slice(1).map((row) => row[0] + '|' + row[4])).size, long.length - 1, 'Metric and Group are a key');
    assert.ok(await page.locator('#metrics-table option').evaluateAll((options) => options.every((o) => !/Legacy/.test(o.textContent))), 'the word Legacy is not shown');
    await page.locator('#export-metrics > summary').click(); await page.click('#metrics-save');
    await page.waitForFunction(() => /^Saved in plots: /.test(document.getElementById('metrics-saved').textContent));
    { const saved = (await text('#metrics-saved')).match(/^Saved in plots: (\S+\.csv)/)[1]; assert.match(saved, /_metrics_summary_/);
      assert.ok(fs.existsSync(path.join(out, 'plots', saved.replace(/\.csv$/, '.metadata.json'))), 'the summary is saved with the record of how the page was set'); }
    await page.selectOption('#metrics-table', 'groups');
    assert.match(await text('#table-metrics'), /No rows were saved/);
    await page.selectOption('#metrics-table', 'metrics'); assert.match(await text('#metric-help'), /signed/); assert.doesNotMatch(await text('#metric-help'), /Legacy/);
    let metrics = await tableNow('table-metrics');
    assert.deepEqual(metrics.head, ['Cell', 'Call', 'Normalization', 'Peak, ΔF/F₀', 'Time to peak, s', 'Rise time, s', 'FWHM, s', 'AUC, ΔF/F₀ × s', 'Signal-to-noise', 'Baseline SD, ΔF/F₀', 'F₀, gray levels'],
      'the call beside the cell; the saved Status under a heading of its own');
    assert.deepEqual([at(metrics, 'Cell001', 'Call').text, at(metrics, 'Cell001', 'Normalization').text, at(metrics, 'Cell001', 'Time to peak, s').text, at(metrics, 'Cell001', 'F₀, gray levels').text], ['Responded', 'Included', '6.00', '10.00']);
    assert.deepEqual(at(metrics, 'Cell002', 'Time to peak, s'), { text: '–', tip: 'Saved: 2.5 · Did not respond' }, 'whole-recording timing of a cell that did not respond: a dash');
    assert.ok(await page.isVisible('#all-metric-columns-label'));
    await page.check('#all-metric-columns'); metrics = await tableNow('table-metrics');
    assert.ok(metrics.head.includes('Peak above baseline, ΔF/F₀') && at(metrics, 'Cell005', 'FWHM cut off').text === 'yes', 'every saved column with the box ticked');
    await page.uncheck('#all-metric-columns');
    // Processed data: Frame and Time stay in view when the cells scroll sideways; headers say what the file calls them.
    await page.click('.tab[data-tab=data]');
    const data = await tableNow('table-data');
    assert.deepEqual([data.head.slice(0, 3), data.titles.slice(0, 3), data.rows[1].slice(0, 2), data.frozen.slice(0, 3)],
      [['Frame', 'Time, s', 'Cell001'], ['Frame. Saved as Frame.', 'Time. Saved as Time_s.', 'Cell001, ΔF/F₀. Saved as Cell001.'], ['2', '0.50'], [true, true, false]]);
    assert.ok(data.rows.every((row) => row.slice(2).every((v) => v === '–' || /^−?\d+\.\d{3}$/.test(v))), 'ΔF/F₀ with 3 decimals');
    assert.ok(await page.locator('#table-data').evaluate((t) => { const box = t.parentElement; box.style.maxWidth = '180px'; box.scrollLeft = 400;
      const s = box.getBoundingClientRect(), frame = t.tBodies[0].rows[0].cells[0].getBoundingClientRect(), time = t.tBodies[0].rows[0].cells[1].getBoundingClientRect(); const ok = box.scrollLeft > 0 && Math.abs(frame.left - s.left) < 3 && Math.abs(time.left - frame.right) < 2;
      box.style.maxWidth = ''; box.scrollLeft = 0; return ok; }), 'Frame and Time stay in place side by side');
    await page.click('.tab[data-tab=metrics]');
    await page.click('#log'); await page.waitForFunction(() => /Opened .* on the Fiji computer/.test(document.getElementById('open-feedback').textContent));

    const opened = page.waitForResponse((response) => response.url().endsWith('/api/open') && response.status() === 200);
    await page.click('#open-folder'); await opened;
    assert.equal(fs.readFileSync(path.join(out, 'opened.txt'), 'utf8'), 'browser-results');

    // A changed timestamp draft does not retime or recalculate completed results.
    await page.click('.tab[data-tab=review]');
    const applied = await page.evaluate(() => ({ fingerprint: window.calciumStep4.S.results.appliedFingerprint,
      time: window.calciumStep4.frameTimes(), calls: JSON.stringify(window.calciumStep4.S.results.callsByWindow) }));
    const staleState = await page.evaluate(() => window.calciumStep4.S.state);
    let releaseStalePoll, announceTrappedPoll, trapTimeout, trapped = false;
    const stalePollGate = new Promise((resolve) => { releaseStalePoll = resolve; });
    const pendingPoll = new Promise((resolve, reject) => {
      announceTrappedPoll = resolve;
      trapTimeout = setTimeout(() => reject(new Error('Expected a periodic state request to test stale responses')), 10000);
    });
    await page.route('**/api/state', async (route) => {
      if (trapped) { await route.continue(); return; }
      trapped = true; clearTimeout(trapTimeout); announceTrappedPoll(); await stalePollGate;
      await route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(staleState) });
    });
    await pendingPoll;
    const pendingClock = 'Frame,Time_s\n' + applied.time.map((t, i) => (i + 1) + ',' + (t + 100)).join('\n') + '\n';
    await page.setInputFiles('#timestamps-file', { name: 'pending-clock.csv', mimeType: 'text/csv', buffer: Buffer.from(pendingClock) });
    await page.waitForFunction(() => window.calciumStep4.S.state.settingsPending && document.getElementById('review-stage').textContent === 'Changes pending');
    releaseStalePoll(); await page.waitForFunction(() => !window.calciumStep4.S.polling); await page.unroute('**/api/state');
    assert.equal(await page.evaluate(() => window.calciumStep4.exportMetadata().settingsPending), true, 'an older poll response must not erase pending provenance');
    await slider(16);
    assert.match(await text('#time-label'), /12\.50 s.*Frame 16 \/ 30/);
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.frameTimes()), applied.time);
    assert.equal(await page.evaluate(() => JSON.stringify(window.calciumStep4.S.results.callsByWindow)), applied.calls);
    assert.equal(await page.evaluate(() => window.calciumStep4.S.results.appliedFingerprint), applied.fingerprint);
    assert.match(await text('#applied-note'), /last completed analysis/);

    // PNG and metadata are exported together from the applied result snapshot.
    const downloads = [];
    const downloadedPlots = new Promise((resolve, reject) => {
      const timeout = setTimeout(() => { page.off('download', listener); reject(new Error('Expected plot PNG and provenance JSON downloads')); }, 10000);
      const listener = (download) => { downloads.push(download); if (downloads.length === 2) { clearTimeout(timeout); page.off('download', listener); resolve(); } };
      page.on('download', listener);
    });
    assert.ok(await page.isHidden('#save-plot'), 'the save buttons are inside the Export menu');
    await page.locator('#export-plot > summary').click(); await page.click('#save-plot'); await downloadedPlots;
    assert.ok(!(await page.locator('#export-plot').evaluate((menu) => menu.open)), 'choosing an entry closes the menu');
    await page.locator('#export-picture > summary').click(); await page.locator('#export-picture .save-options-toggle').click();
    assert.ok(await page.isVisible('#save-options') && !(await page.locator('#export-picture').evaluate((menu) => menu.open)), 'Options opens from the picture\'s Export menu');
    assert.ok(await page.evaluate(() => { const panel = document.getElementById('save-options').getBoundingClientRect(), button = document.querySelector('#export-picture > summary').getBoundingClientRect(); return (panel.bottom <= button.top + 1 || panel.top >= button.bottom - 1) && Math.abs(panel.right - button.right) < 24; }), 'the options sit against the Export button');
    assert.equal(await page.locator('#export-plot .save-options-toggle, #export-heatmap .save-options-toggle').count(), 0, 'the plot and the heatmap have Edit figure in place of Options');
    assert.doesNotMatch(await page.textContent('#save-pixels'), /Plot/, 'Options speaks of pictures and frames only');
    await page.keyboard.press('Escape'); await page.locator('#plot').click({ position: { x: 5, y: 5 } }); assert.ok(await page.isHidden('#save-options'));
    for (const download of downloads) await download.saveAs(path.join(out, download.suggestedFilename()));
    const names = downloads.map((d) => d.suggestedFilename()), image = names.find((n) => n.endsWith('.png'));
    assert.match(image, /^Fixture_recording_traces_Cell\w+_\d{8}-\d{6}\.png$/, 'plot names carry the recording, the cell and the time');
    assert.ok(names.includes(image.replace(/\.png$/, '.metadata.json')));
    const png = fs.readFileSync(path.join(out, image));
    assert.deepEqual([png.readUInt32BE(16), png.readUInt32BE(20)], [2100, 1200], 'a plot figure nobody has edited is 7 × 4 in at 300 dpi, whatever the window size');
    assert.ok(await page.evaluate(() => { const g = window.calciumStep4.Plot.geom, c = document.getElementById('plot'); return g.L + g.W <= c.clientWidth && g.T + g.H <= c.clientHeight; }), 'saving a figure leaves the page\'s own plot as it was');
    const metadata = JSON.parse(fs.readFileSync(path.join(out, image.replace(/\.png$/, '.metadata.json')), 'utf8'));
    assert.ok(metadata.results_folder && metadata.time_range.every((t) => typeof t === 'number'), 'metadata names the results folder and the shown time range');
    assert.equal(metadata.appliedFingerprint, applied.fingerprint);
    assert.equal(metadata.settingsPending, true); assert.equal(metadata.time_s, 12.5);
    assert.equal(metadata.manifest.protocol_hash, 'synthetic-browser-fixture');
    assert.equal(metadata.picture.size, '7x4'); assert.equal(metadata.picture.figure.kind, 'plot'); assert.equal(metadata.picture.figure.title, null, 'the file beside a figure records how it was set');
    // The figure editor: its preview is the picture that gets saved; what is set there is kept, except what belongs to this recording.
    // (drawn, not only asked for: the stage still holds the figure that was open before, and the plot and the heatmap no longer share a size)
    const previewReady = () => page.waitForFunction(() => !window.calciumStep4.Fig.queued && document.getElementById('figure-preview').width > 300 && /pixels/.test(document.getElementById('figure-pixels').textContent));
    const previewNow = () => page.evaluate(() => document.getElementById('figure-preview').toDataURL());
    await page.locator('#export-plot > summary').click(); await page.click('#edit-plot'); await previewReady();
    assert.ok(await page.isVisible('#figure-editor') && !(await page.locator('#export-plot').evaluate((menu) => menu.open)), 'Edit figure opens the editor from the Export menu');
    assert.equal(await page.textContent('#figure-name'), 'Plot figure'); assert.match(await page.textContent('#figure-pixels'), /^2100 × 1200 pixels/);
    assert.ok(await page.evaluate(() => { const b = document.getElementById('figure-preview').getBoundingClientRect(), s = document.getElementById('figure-stage').getBoundingClientRect(); return Math.abs(b.width / b.height - 7 / 4) < 0.01 && b.left >= s.left && b.right <= s.right && b.top >= s.top && b.bottom <= s.bottom; }), 'the preview has the figure\'s shape and fits beside the boxes');
    assert.deepEqual(await oddShapes(), [], 'the editor\'s boxes have the page\'s height and corner');
    assert.deepEqual(await page.evaluate(() => { const rows = Array.from(document.querySelectorAll('#figure-controls .row')).filter((r) => r.offsetParent), edges = new Set(), units = new Set();
      for (const row of rows) { const box = row.querySelector('select, input:not([type=checkbox])'), unit = row.querySelector('.unit i'); if (!box) continue; const b = box.getBoundingClientRect(); edges.add(Math.round(b.width) + ' wide, right at ' + Math.round(b.right)); if (unit) units.add(unit.getBoundingClientRect().left > b.right); }
      return [edges.size, Array.from(units)]; }), [1, [true]], 'every box of the editor has one width and one right edge, the units outside it');
    assert.ok(await page.evaluate(() => Array.from(document.body.children).every((part) => part.id === 'figure-editor' || part.tagName === 'SCRIPT' || part.inert)), 'the page behind the editor is out of reach');
    assert.equal(await page.inputValue('#figure-title'), '', 'a figure has no title until one is typed'); assert.equal(await page.getAttribute('#figure-title', 'placeholder'), 'No title', 'an empty title box does not read as a title set to "none"'); assert.equal(await page.inputValue('#figure-yTitle'), 'ΔF/F₀'); assert.equal(await page.getAttribute('#figure-yMax', 'placeholder'), 'auto');
    assert.deepEqual(await page.evaluate(() => { const S = window.calciumStep4.S, has = (id) => !!document.getElementById('figure-' + id); return { mode: S.plotMode, every: S.every, cell: has('cellWidth'), levels: has('levels'), traces: has('traceWidth'), calls: has('colourResponded'), mean: has('colourMean') }; }),
      { mode: 'mean', every: true, cell: false, levels: false, traces: true, calls: true, mean: true }, 'boxes are offered for what the plot on screen shows: the mean and all traces, no single cell');
    assert.ok(await page.isHidden('#figure-save-compression'), 'TIFF compression shows with TIFF');
    let shown = await previewNow();
    // The legend can sit right of the plot, its entries stacked.
    assert.deepEqual(await page.locator('#figure-legend option').evaluateAll((options) => options.map((o) => o.value + ': ' + o.textContent)), ['top: Above the plot', 'right: Right of the plot', 'none: Not shown']);
    await page.selectOption('#figure-legend', 'right'); await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was, shown);
    assert.ok(await page.evaluate(() => { const d = window.calciumStep4.Fig.drawn; return d.kind === 'plot' && d.legend.x >= d.L + d.W; }), 'the preview shows the legend right of the plot');
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.figureKept('plot')), { legend: 'right' }, 'and the choice is kept');
    const atRight = await previewNow(); await page.selectOption('#figure-legend', 'top'); await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was, atRight);
    for (const [id, value] of [['width', '6'], ['height', '5'], ['title', 'Fixture title'], ['yMax', '3'], ['text', '10']]) await typeBox(id, value);
    await page.selectOption('#figure-font', 'times');
    await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was && /^1800 × 1500 pixels/.test(document.getElementById('figure-pixels').textContent), shown);
    assert.ok(await page.evaluate(() => { const b = document.getElementById('figure-preview').getBoundingClientRect(); return Math.abs(b.width / b.height - 6 / 5) < 0.01; }), 'the preview follows the size');
    for (const [id, state, expected] of [['colourMean', '#c0392b', 'a colour'], ['grid', true, 'a tick']]) {   // every kind of box redraws the preview (grid lines are off until ticked)
      shown = await previewNow();
      if (state === true) await page.check('#figure-' + id); else await page.locator('#figure-' + id).evaluate((box, value) => { box.value = value; box.dispatchEvent(new Event('input', { bubbles: true })); }, state);
      await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was, shown).catch(() => assert.fail('the preview did not change with ' + expected));
    }
    await typeBox('width', '0.2'); assert.deepEqual(await boxSays('width'), ['0.2', true, '1 to 20 in']); assert.equal(await page.evaluate(() => window.calciumStep4.Fig.plot.width), 6, 'a size that makes no figure is not taken');
    await page.locator('#figure-width').blur(); assert.deepEqual(await boxSays('width'), ['0.2', true, '1 to 20 in'], 'and the box keeps what was typed, marked, with the widths it takes under it'); await typeBox('width', '6');
    await typeBox('yMax', ''); assert.equal(await page.evaluate(() => window.calciumStep4.Fig.plot.yMax), null, 'an emptied axis box is automatic again'); await typeBox('yMax', '3');
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.figureKept('plot')), { width: 6, height: 5, font: 'times', text: 10, grid: true, colourMean: '#c0392b' }, 'the title and the axis numbers stay with this recording; the rest is kept');
    await page.waitForFunction(async () => { const kept = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs['step4.figure.plot']; return !!kept && JSON.parse(kept).font === 'times' && JSON.parse(kept).grid === true; });   // Fiji keeps it for the next recording
    const figureFiles = new Promise((resolve, reject) => { const got = [], timeout = setTimeout(() => { page.off('download', listener); reject(new Error('Expected the figure and how it was made')); }, 10000);
      const listener = (download) => { got.push(download); if (got.length === 2) { clearTimeout(timeout); page.off('download', listener); resolve(got); } }; page.on('download', listener); });
    shown = await previewNow(); await page.selectOption('#figure-sem', 'none'); await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was, shown);
    await page.click('#figure-download'); const figureNames = [];
    for (const download of await figureFiles) { figureNames.push(download.suggestedFilename()); await download.saveAs(path.join(out, 'figure_' + download.suggestedFilename())); }
    { // the file beside a figure describes the figure, not the page: its own Y and time range and the means it drew, without "± SEM" when no band is drawn
      const record = JSON.parse(fs.readFileSync(path.join(out, 'figure_' + figureNames.find((n) => n.endsWith('.metadata.json'))), 'utf8')), drew = await page.evaluate(() => { const d = window.calciumStep4.Fig.drawn; return { kind: d.kind, x: [d.x0, d.x1], y: [d.lo, d.hi], plotted: d.plotted }; });
      assert.equal(drew.kind, 'plot'); assert.deepEqual(record.y_axis, drew.y, 'the Y range drawn in the figure'); assert.equal(record.y_axis[1], 3, 'with the Y max typed in the editor');
      assert.deepEqual(record.time_range, drew.x, 'the time range drawn in the figure'); assert.deepEqual(record.plotted, drew.plotted, 'the means the figure drew');
      assert.ok(record.plotted.length && record.plotted.every((label) => !/SEM/.test(label)), 'no SEM band, so no "± SEM"'); assert.equal(record.picture.figure.sem, 'none');
      assert.notDeepEqual(record.y_axis, await page.evaluate(() => [window.calciumStep4.Plot.geom.lo, window.calciumStep4.Plot.geom.hi]), 'not the page\'s own Y range'); }
    await page.selectOption('#figure-sem', 'band');
    const figurePng = fs.readFileSync(path.join(out, 'figure_' + figureNames.find((n) => n.endsWith('.png'))));
    assert.deepEqual([figurePng.readUInt32BE(16), figurePng.readUInt32BE(20)], [1800, 1500], 'the saved figure has the size set in the editor');
    const figureMade = JSON.parse(fs.readFileSync(path.join(out, 'figure_' + figureNames.find((n) => n.endsWith('.metadata.json'))), 'utf8')).picture;
    assert.equal(figureMade.size, '6x5'); assert.equal(figureMade.dpi, 300); assert.deepEqual([figureMade.figure.title, figureMade.figure.yMax, figureMade.figure.grid, figureMade.figure.font], ['Fixture title', 3, true, 'times']);
    await page.screenshot({ path: path.join(out, 'page-figure-editor.png') });
    const frameThen = await page.evaluate(() => window.calciumStep4.S.frame);
    await page.keyboard.press('ArrowRight'); await page.keyboard.press('Escape'); assert.ok(await page.isHidden('#figure-editor'));
    assert.ok(await page.evaluate((frame) => window.calciumStep4.S.frame === frame && !window.calciumStep4.S.cell && Array.from(document.body.children).every((part) => !part.inert), frameThen), 'keys pressed in the editor do nothing on the page; Escape closes the editor');
    await page.locator('#export-plot > summary').click(); await page.click('#edit-plot'); await previewReady();
    assert.equal(await page.inputValue('#figure-title'), 'Fixture title', 'the editor reopens as it was left'); await page.click('#figure-reset');
    assert.equal(await page.inputValue('#figure-title'), ''); assert.equal(await page.inputValue('#figure-width'), '7'); assert.deepEqual(await page.evaluate(() => window.calciumStep4.figureKept('plot')), {});
    await page.waitForFunction(async () => !(await (await fetch('api/session', { cache: 'no-store' })).json()).prefs['step4.figure.plot']);
    await page.evaluate(() => { const f = window.calciumStep4; f.figureRestore('plot', JSON.stringify({ width: 500, font: 'comic', ink: 'red', grid: 'no', title: 'kept?', text: 11 })); });   // a damaged preference changes only what is valid
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.figureKept('plot')), { text: 11 }); await page.click('#figure-reset'); await page.click('#figure-close');
    // ---------------------------------------------------------------- the figure editor's rules: numbers, limits, text that does not fit, Reset
    const figureNow = (kind) => page.evaluate((kind) => Object.assign({}, window.calciumStep4.Fig[kind]), kind || 'plot');
    const drawnNow = async () => { await page.waitForFunction(() => !window.calciumStep4.Fig.queued); await settledPaint();
      return page.evaluate(() => { const d = window.calciumStep4.Fig.drawn; return { x: [d.x0, d.x1], y: [d.lo, d.hi], share: d.share, body: d.body, room: d.room }; }); };
    const asSaid = (v) => String(+v.toPrecision(6)).replace('-', '−');   // a number as the editor writes it in a message
    const openEditor = async (kind) => { await page.locator('#export-' + kind + ' > summary').click(); await page.click('#edit-' + kind); await previewReady(); };
    await openEditor('plot');
    // What the editor says about itself: Reset means the built-in look; format and resolution are every saved picture's; the fonts.
    assert.equal(await page.textContent('#figure-reset'), 'Reset to defaults', 'the button says what it resets to');
    assert.deepEqual(await page.evaluate(() => { const head = Array.from(document.querySelectorAll('#figure-controls h3')).find((h) => h.firstChild.textContent === 'File'), rows = [];
      for (let row = head.nextElementSibling; row && row.tagName !== 'H3'; row = row.nextElementSibling) rows.push(row.querySelector('span').textContent); return [head.querySelector('span').textContent, rows]; }),
      ['for every saved picture', ['Resolution', 'Format', 'TIFF compression']], 'resolution and format stand under a heading that says they are every saved picture\'s');
    assert.deepEqual(await page.locator('#figure-font option').evaluateAll((options) => options.map((o) => o.value + ': ' + o.textContent)), ['arial: Arial', 'helvetica: Helvetica', 'times: Times New Roman', 'page: As on the page'], 'Arial first, the page\'s own font last, no Courier');
    assert.equal(await page.inputValue('#figure-font'), 'arial', 'a figure is written in Arial until another font is chosen');
    assert.ok(await page.isHidden('#figure-earlier'), 'nothing is kept from an earlier session, and the editor does not say that something is');
    // A number is taken whole. Typed key by key, 40 is never the 4 it passes through: the box keeps what was typed, is marked, and
    // says under it which widths it takes; the figure keeps its width; nothing is kept for the next recording.
    const prefsSent = []; page.on('request', (request) => { if (/api\/prefs$/.test(request.url())) prefsSent.push(request.postData()); });
    const keys = async (id, typed) => { await page.click('#figure-' + id); await page.keyboard.press('ControlOrMeta+a'); await page.keyboard.type(typed); };
    await keys('width', '40');
    assert.deepEqual(await boxSays('width'), ['40', true, '1 to 20 in'], 'while it is typed the box says which widths it takes'); assert.equal((await figureNow()).width, 7);
    await page.keyboard.press('Enter'); await page.keyboard.press('Tab'); await settledPaint();
    assert.deepEqual(await boxSays('width'), ['40', true, '1 to 20 in'], 'after Enter and after leaving the box, what was typed is still there, marked'); assert.equal((await figureNow()).width, 7, 'and the figure keeps its width');
    assert.deepEqual(prefsSent, [], 'a number that was not taken is not kept either');
    await typeBox('height', '5'); assert.deepEqual([await boxSays('width'), await boxSays('height'), (await figureNow()).width, (await figureNow()).height], [['40', true, '1 to 20 in'], ['5', false, ''], 7, 5], 'it stays in its box, marked, while other boxes are used'); await typeBox('height', '4');
    assert.ok(await page.evaluate(() => { const look = (id) => getComputedStyle(document.getElementById(id)), marked = look('figure-width'), plain = look('figure-height'); return marked.borderTopColor !== plain.borderTopColor && marked.backgroundColor !== plain.backgroundColor; }), 'a marked box looks marked');
    assert.equal(await page.getAttribute('#figure-width', 'aria-describedby'), 'figure-width-note', 'and a screen reader is given the line under it');
    await keys('width', '12'); assert.equal((await figureNow()).width, 7, 'a number the figure can use is not taken while it is typed either'); await page.keyboard.press('Enter');
    assert.equal((await figureNow()).width, 12); assert.deepEqual(await boxSays('width'), ['12', false, '']);
    await typeBox('width', '7,5'); assert.deepEqual([(await figureNow()).width, await boxSays('width')], [7.5, ['7.5', false, '']], 'a decimal comma is a decimal point'); await typeBox('width', '7');
    for (const typed of ['0', '0.5', '25', '']) { await typeBox('width', typed); assert.deepEqual([await boxSays('width'), (await figureNow()).width], [[typed, true, '1 to 20 in'], 7], 'a width of "' + typed + '" is not taken'); }
    await typeBox('width', '7');
    // Every number box has limits and says them: text that is no number (also "1e1" and "0x10", which a program would read as 10 and
    // 16) is left in the box, marked, with what the box takes under it, and changes nothing in the figure.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, was = JSON.stringify(f.Fig.plot), told = {}, type = (box, text) => { box.value = text; box.dispatchEvent(new Event('input', { bubbles: true })); box.dispatchEvent(new Event('change', { bubbles: true })); };
      for (const note of document.querySelectorAll('#figure-controls .range-note')) { const box = document.getElementById(note.id.replace(/-note$/, '')), held = box.value;
        told[box.id.replace('figure-', '')] = ['abc', '1e1', '0x10'].map((text) => { type(box, text); return box.value === text && box.getAttribute('aria-invalid') === 'true' && !note.hidden ? note.textContent : 'taken'; }).filter((words, i, all) => all.indexOf(words) === i).join(' / '); type(box, held); }
      return Object.assign(told, { unchanged: JSON.stringify(f.Fig.plot) === was, boxes: document.querySelectorAll('#figure-controls input[inputmode=decimal]').length }); }),
      { width: '1 to 20 in', height: '1 to 20 in', text: '4 to 48 pt', titleText: '4 to 72 pt', xMin: 'Not a number', xMax: 'Not a number', xStep: 'Above 0', yMin: 'Not a number', yMax: 'Not a number', yStep: 'Above 0',
        axisWidth: '0.25 to 6 pt', meanWidth: '0.25 to 10 pt', traceWidth: '0.1 to 6 pt', traceOpacity: '1 to 100 %', responderWidth: '0.1 to 6 pt', responderOpacity: '1 to 100 %', unchanged: true, boxes: 16 });
    // Single cells: the line, colour and opacity of the traces of cells that did not respond, and the line and opacity of the
    // responders' traces drawn over them. Opacity is typed in percent; an empty box is "auto", fainter the more cells there are.
    {
      const singles = async () => { const seen = await drawing(true), inks = await page.evaluate(() => [window.calciumStep4.STYLE.plot.colourTrace, window.calciumStep4.STYLE.plot.colourResponded]), PT = 96 / 72;
        const of = (ink) => seen.lines.filter((l) => l.ink === ink && l.alpha < 1 && !inLegend(seen, l.x0, l.y0)).map((l) => [+l.alpha.toFixed(3), +(l.width / PT).toFixed(3)]);
        return { others: [...new Set(of(inks[0]).map(String))], responders: [...new Set(of(inks[1]).map(String))] }; };
      const auto = await singles(); assert.ok(auto.others.length === 1 && auto.responders.length === 1, 'the fixture draws both kinds of trace: ' + JSON.stringify(auto));
      assert.deepEqual([await page.getAttribute('#figure-traceOpacity', 'placeholder'), await page.inputValue('#figure-traceOpacity'), await page.inputValue('#figure-responderWidth')], ['auto', '', '0.9']);
      await typeBox('traceOpacity', '35'); await typeBox('responderWidth', '2'); await typeBox('responderOpacity', '80');
      assert.deepEqual([(await figureNow()).traceOpacity, (await figureNow()).responderWidth, (await figureNow()).responderOpacity, await boxSays('traceOpacity'), await singles()],
        [0.35, 2, 0.8, ['35', false, ''], { others: ['0.35,0.6'], responders: ['0.8,2'] }], 'the boxes reach the drawing: opacity in percent, lines in points');
      for (const typed of ['0', '0.5', '150', 'abc']) { await typeBox('responderOpacity', typed); assert.deepEqual([await boxSays('responderOpacity'), (await figureNow()).responderOpacity], [[typed, true, '1 to 100 %'], 0.8], 'an opacity of "' + typed + '" is not taken'); }
      await typeBox('traceOpacity', ''); await typeBox('responderOpacity', ''); await typeBox('responderWidth', '0.9');
      assert.deepEqual([(await figureNow()).traceOpacity, (await figureNow()).responderOpacity, await boxSays('traceOpacity'), await singles()], [null, null, ['', false, ''], auto], 'emptied, opacity is automatic again');
    }
    // Axis limits. A minimum at or above the maximum is not taken: the box just typed in is marked and says so, and the figure keeps
    // the limits it had. A number that waits in its box is taken as soon as the other end lets it.
    const allTime = (await drawnNow()).x; assert.deepEqual(allTime, [0, 19.5]);
    await typeBox('xMin', '15'); assert.deepEqual((await drawnNow()).x, [15, 19.5]);
    await typeBox('xMax', '5'); assert.deepEqual(await boxSays('xMax'), ['5', true, 'Maximum must be above minimum']); assert.deepEqual([(await figureNow()).xMax, (await drawnNow()).x], [null, [15, 19.5]], 'the figure keeps the limits it had');
    await typeBox('xMin', '2'); assert.deepEqual([await boxSays('xMin'), await boxSays('xMax'), (await drawnNow()).x], [['2', false, ''], ['5', false, ''], [2, 5]], 'with a minimum under it, the maximum that waited in its box is taken');
    await typeBox('xMin', '8'); assert.deepEqual([await boxSays('xMin'), (await figureNow()).xMin, (await drawnNow()).x], [['8', true, 'Minimum must be below maximum'], 2, [2, 5]]);
    await typeBox('xMax', ''); assert.deepEqual([await boxSays('xMin'), (await drawnNow()).x], [['8', false, ''], [8, 19.5]], 'an emptied maximum is the automatic one again, and the minimum that waited is under it');
    // (one limit alone is held against the automatic other end, and the message gives that end)
    await typeBox('xMin', '30'); assert.deepEqual([await boxSays('xMin'), (await drawnNow()).x], [['30', true, 'Must be below 19.5'], [8, 19.5]]);
    await typeBox('xMin', ''); await typeBox('xMax', '-3'); assert.deepEqual([await boxSays('xMax'), (await drawnNow()).x], [['-3', true, 'Must be above 0'], allTime]); await typeBox('xMax', '');
    const autoY = (await drawnNow()).y;
    await typeBox('yMin', '5'); assert.deepEqual([await boxSays('yMin'), (await drawnNow()).y], [['5', true, 'Must be below ' + asSaid(autoY[1])], autoY], 'a Y minimum above everything that is drawn is not taken: no empty plot');
    await typeBox('yMax', '9'); assert.deepEqual([await boxSays('yMin'), (await drawnNow()).y], [['5', false, ''], [5, 9]]);
    await typeBox('yMax', '1'); assert.deepEqual([await boxSays('yMax'), (await drawnNow()).y], [['1', true, 'Maximum must be above minimum'], [5, 9]], 'Y min 5 with Y max 1 keeps the range it had');
    // (a limit the figure holds and can no longer use, because the other end has changed under it, is marked and not drawn)
    await typeBox('yMax', ''); assert.deepEqual([await boxSays('yMin'), (await figureNow()).yMin, (await drawnNow()).y], [['5', true, 'Must be below ' + asSaid(autoY[1])], 5, autoY]); await typeBox('yMin', '');
    // Tick spacing: "every 0" is not taken, nor a spacing that would give more ticks than are drawn (60); the box says the least it takes.
    await typeBox('xStep', '0'); assert.deepEqual([await boxSays('xStep'), (await figureNow()).xStep], [['0', true, 'Above 0'], null]);
    await typeBox('xStep', '0.1'); assert.deepEqual([await boxSays('xStep'), (await figureNow()).xStep], [['0.1', true, 'At least 0.33'], null], '19.5 s in steps of 0.1 would be 195 ticks');
    await typeBox('xStep', '0.5'); assert.deepEqual([await boxSays('xStep'), (await figureNow()).xStep], [['0.5', false, ''], 0.5]);
    const xNumbers = (seen) => seen.texts.filter((t) => t.y > seen.top + seen.T + seen.H && /^−?[\d.]+$/.test(t.text)).map((t) => t.text);
    assert.equal(xNumbers(await drawing(true)).length, 40, 'a tick every 0.5 s');
    await typeBox('xMax', '40'); assert.deepEqual(await boxSays('xStep'), ['0.5', true, 'At least 0.67'], 'with the axis twice as long the spacing is too fine: its box says so');
    assert.deepEqual(xNumbers(await drawing(true)), ['0', '5', '10', '15', '20', '25', '30', '35', '40'], 'and the ticks are the automatic ones');
    await typeBox('xMax', ''); assert.deepEqual(await boxSays('xStep'), ['0.5', false, '']); await typeBox('xStep', '');
    // Text that does not fit. A size that would leave the plot no room is not taken (the plot must stay four heights of the text
    // each way). Text that leaves the plot under 40 % of the figure is taken, and the box says that it is too large for the figure.
    await typeBox('height', '3'); const roomy = await drawnNow(); assert.ok(roomy.share > 0.4 && roomy.room > 4); assert.deepEqual(await boxSays('text'), ['9', false, '']);
    await typeBox('text', '30'); assert.deepEqual([await boxSays('text'), (await figureNow()).text], [['30', true, 'Too large for this figure size'], 9], 'text size 30 on 7 × 3 in is not taken');
    assert.deepEqual(await drawnNow(), roomy, 'and the figure is as it was: the legend and the labels have not swallowed the plot');
    await typeBox('text', '16'); const tight = await drawnNow();
    assert.ok(tight.share < 0.4 && tight.room >= 4, 'at 16 pt the plot has under 40 % of the figure: ' + JSON.stringify(tight));
    assert.deepEqual([await boxSays('text'), (await figureNow()).text], [['16', false, 'Too large for this figure size'], 16], 'that is taken, and said under the box');
    assert.ok(await page.isEnabled('#figure-save') && !/no room/.test(await page.textContent('#figure-pixels')), 'and it can be saved');
    {   // no two numbers of the Y axis on one another: a low plot has fewer of them
      const seen = await drawing(true), ticksY = yTicks(seen); assert.ok(ticksY.length >= 2 && ticksY.every((t, i) => i === 0 || t.y - ticksY[i - 1].y >= t.size), 'the Y numbers of a low plot are clear of each other: ' + ticksY.map((t) => Math.round(t.y)));
    }
    await typeBox('text', '9'); await drawnNow(); assert.deepEqual(await boxSays('text'), ['9', false, ''], 'with the text it had, nothing is said');
    // (the title's size is the one that is too large when the plot would have its room without the title)
    await typeBox('height', '4'); await page.fill('#figure-title', 'A title'); await typeBox('titleText', '60'); const titled = await drawnNow();
    assert.ok(titled.share < 0.4 && titled.body >= 0.4, JSON.stringify(titled)); assert.deepEqual([await boxSays('titleText'), await boxSays('text')], [['60', false, 'Too large for this figure size'], ['9', false, '']]);
    await typeBox('titleText', '10'); await page.fill('#figure-title', ''); await page.press('#figure-title', 'Enter');
    // (and a figure too small for its text is not taken either)
    await typeBox('width', '1'); assert.deepEqual([await boxSays('width'), (await figureNow()).width], [['1', true, 'Too small for this text size'], 7]); await typeBox('width', '7');
    // A figure whose plot has no room by another way (settings kept for a plot with a shorter legend, say) is shown as it is, with
    // the reason, and is not saved: not from the editor, not from the Export menu.
    await page.click('#figure-close'); await figureStyle({ width: 1.3, height: 1.3 });
    await page.locator('#export-plot > summary').click(); await page.click('#save-plot-folder'); await page.waitForFunction(() => /Not saved/.test(document.getElementById('plot-saved').textContent));
    assert.equal(await page.textContent('#plot-saved'), 'Not saved: The text is too large for this figure size. Change it in Edit figure.');
    await page.locator('#export-plot > summary').click(); await page.click('#edit-plot'); await page.waitForFunction(() => !window.calciumStep4.Fig.queued && /^390 × 390 pixels/.test(document.getElementById('figure-pixels').textContent));   // (a preview too small for previewReady)
    assert.match(await page.textContent('#figure-pixels'), /^390 × 390 pixels: no room for the plot\. Choose smaller text or a larger size\.$/);
    assert.ok(await page.isDisabled('#figure-save') && await page.isDisabled('#figure-download')); assert.deepEqual(await boxSays('text'), ['9', false, 'Too large for this figure size']);
    await typeBox('text', '5'); assert.ok((await drawnNow()).room >= 4 && await page.isEnabled('#figure-save'), 'smaller text gives the plot its room back'); await page.click('#figure-reset'); await previewReady();
    // Small figures. A title too long for the figure is shortened, with an ellipsis; so are the axis titles and the stimulus label.
    // A legend entry too long for its row goes on two lines or more instead, whole, at the top and in a legend at the right, which
    // takes no more than a third of the width and stays clear of the figure's edge. Nothing runs out of the figure, the legend still
    // wraps by whole entries, and the plot starts under it.
    await page.click('#figure-close'); const longTitle = '2026-11-18_C1-Control_Run2_Capsaicin-0.5uM_300fr_src008';
    await figureStyle({ width: 2, title: longTitle }); fig = await drawing(true);
    const inFigure = (seen, wide) => seen.texts.filter((t) => t.align === 'left' && t.m.startsWith('1,0,0,1')).every((t) => t.x >= 0 && t.x + t.w <= wide * 96 + 0.5), titleOf = (seen) => seen.texts[seen.texts.length - 1];   // (text that starts at its x: the title, the legend)
    assert.ok(/^2026-11-18.+…$/.test(titleOf(fig).text) && longTitle.startsWith(titleOf(fig).text.slice(0, -1)) && titleOf(fig).x === fig.L, 'the title is shortened: ' + titleOf(fig).text);
    assert.ok(inFigure(fig, 2), 'title and legend end inside a figure 2 in wide');
    assert.deepEqual(fig.legend, ['all cells ± SEM (n = 4)', 'single cells', 'baseline', 'response window']); assert.equal(new Set(fig.places.map((place) => place[1])).size, 4, 'each entry on a row of its own');
    const whole = (seen) => seen.wrapped.every((lines, i) => lines.join(' ') === seen.legend[i]) && !seen.texts.some((t) => inLegend(seen, t.x, t.y) && /…/.test(t.text));
    assert.ok(fig.wrapped[0].length > 1 && fig.wrapped.slice(1).every((lines) => lines.length === 1) && whole(fig), 'the one entry too long for its row goes on two lines, whole: ' + JSON.stringify(fig.wrapped));
    assert.ok(fig.legendBox[1] + fig.legendBox[3] <= fig.top + fig.T, 'the plot starts under the taller legend');
    await figureStyle({ title: longTitle }); fig = await drawing(true); assert.equal(titleOf(fig).text, longTitle, 'a title that fits is written whole');
    await figureStyle({ width: 3, legend: 'right' }); fig = await drawing(true);
    assert.ok(fig.legendWidth <= 96 + 0.5 && inFigure(fig, 3) && fig.wrapped.some((lines) => lines.length > 1) && whole(fig) && fig.W > 96, 'a legend at the right takes at most a third of the width, its entries whole: ' + [fig.legendWidth, fig.W, JSON.stringify(fig.wrapped)]);
    assert.ok(3 * 96 - (fig.legendAt[0] + fig.legendWidth) >= 6 && fig.legendAt[0] >= fig.L + fig.W, 'and stays clear of the plot and of the figure\'s edge: ' + [fig.legendAt[0], fig.legendWidth, fig.L + fig.W]);
    await figureStyle({ width: 3, xTitle: 'x'.repeat(80), yTitle: 'y'.repeat(80), stimulusText: 's'.repeat(80) }); fig = await drawing(true);
    assert.deepEqual(['x', 'y', 's'].map((letter) => { const t = fig.texts.find((t) => t.text.startsWith(letter + letter)); return !!t && /…$/.test(t.text) && t.lw <= (letter === 'y' ? 4 : 3) * 96; }), [true, true, true], 'the axis titles and the stimulus label are shortened to the figure');
    await figureStyle();
    {   // the heatmap: its title, and the colour bar's, which makes room for itself up to a quarter of the width
      await heatStyle({ width: 3, title: longTitle, barTitle: 'b'.repeat(80) }); const heat = await heatDrawing(true), bar = heat.texts.find((t) => /^bb/.test(t.text));
      assert.ok(/…$/.test(heat.texts[heat.texts.length - 1].text) && /…$/.test(bar.text) && bar.w <= 72 + 0.5 && heat.W > 96 && inFigure(heat, 3), 'a long colour bar title does not take the picture\'s place: ' + [bar.w, heat.W]); await heatStyle();
    }
    // Reset to defaults puts every box of the editor back, also the ones that are every saved picture's, and nothing stays kept.
    await openEditor('plot');
    await page.selectOption('#figure-save-format', 'tiff'); await page.selectOption('#figure-save-dpi', '600'); await page.selectOption('#figure-save-compression', 'zip'); await page.check('#figure-grid'); await typeBox('width', '5'); await typeBox('yMax', '3'); await typeBox('text', 'abc');
    await page.waitForFunction(async () => { const kept = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs; return kept['step4.format'] === 'tiff' && kept['step4.dpi'] === '600' && kept['step4.compression'] === 'zip' && /"width":5/.test(kept['step4.figure.plot'] || ''); });
    await page.click('#figure-reset'); await previewReady();
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, v = (id) => document.getElementById(id).value; return [JSON.stringify(f.Fig.plot) === JSON.stringify(f.STYLE.plot), v('figure-save-format'), v('figure-save-dpi'), v('save-format'), v('save-dpi'), v('save-compression'), v('figure-text'), document.querySelectorAll('#figure-controls [aria-invalid]').length]; }),
      [true, 'png', '300', 'png', '300', 'lzw', '9', 0]);
    assert.match(await page.textContent('#figure-pixels'), /^2100 × 1200 pixels/);
    await page.waitForFunction(async () => { const kept = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs; return !kept['step4.figure.plot'] && kept['step4.format'] === 'png' && kept['step4.dpi'] === '300' && kept['step4.compression'] === 'lzw'; });
    await page.click('#figure-close'); await page.click('.tab[data-tab=heatmap]'); await openEditor('heatmap');
    await page.selectOption('#figure-palette', 'plasma'); await page.selectOption('#figure-range', 'full'); await page.selectOption('#figure-order', 'amp'); await page.check('#figure-row-names'); await page.check('#figure-hide-excluded'); await page.uncheck('#figure-bar'); await typeBox('barMax', '2');
    // (the colour scale has no automatic end to be held against: one limit moves the other; two that are no range are not taken)
    await typeBox('barMin', '3'); assert.deepEqual([await boxSays('barMin'), (await figureNow('heatmap')).barMin], [['3', true, 'Minimum must be below maximum'], null]);
    await typeBox('xMin', '50'); assert.deepEqual(await boxSays('xMin'), ['50', true, 'Must be below 19.5']); await typeBox('xStep', '0.2'); assert.deepEqual(await boxSays('xStep'), ['0.2', true, 'At least 0.33']);
    await page.click('#figure-reset'); await previewReady();
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, v = (id) => document.getElementById(id).value, on = (id) => document.getElementById(id).checked;
      return [JSON.stringify(f.Fig.heatmap) === JSON.stringify(f.STYLE.heatmap), v('palette'), v('range'), v('order'), on('row-names'), on('hide-excluded'), f.S.palette, f.S.range, f.S.order, f.S.rowNames, f.S.hideExcluded, v('figure-palette'), on('figure-row-names'), document.querySelectorAll('#figure-controls [aria-invalid]').length]; }),
      [true, 'signed', 'robust', 'call', false, false, 'signed', 'robust', 'call', false, false, 'signed', false, 0], 'the heatmap\'s colours, range, order and rows are back as the page opens with them');
    await page.click('#figure-close'); await page.click('.tab[data-tab=review]');
    // Settings kept from an earlier session: the editor says in one line that they are in use, for as long as one of them is. A font
    // that is no longer offered is not taken back.
    await page.evaluate(() => window.calciumStep4.figureRestore('plot', JSON.stringify({ width: 6, font: 'courier' })));
    await openEditor('plot'); assert.ok(await page.isVisible('#figure-earlier')); assert.equal(await page.textContent('#figure-earlier'), 'Settings kept from an earlier session are in use.');
    assert.deepEqual([await page.inputValue('#figure-width'), await page.inputValue('#figure-font')], ['6', 'arial']);
    assert.ok(await page.evaluate(() => { const line = document.getElementById('figure-earlier'); return line.getBoundingClientRect().height < 2 * parseFloat(getComputedStyle(line).fontSize); }), 'in one line');
    await typeBox('width', '7'); assert.ok(await page.isHidden('#figure-earlier'), 'no longer, once the kept width is replaced');
    await page.evaluate(() => window.calciumStep4.figureRestore('plot', JSON.stringify({ width: 6 }))); await page.click('#figure-close'); await openEditor('plot'); assert.ok(await page.isVisible('#figure-earlier'));
    await page.click('#figure-reset'); await previewReady(); assert.ok(await page.isHidden('#figure-earlier') && await page.inputValue('#figure-width') === '7', 'Reset to defaults leaves nothing of them'); await page.click('#figure-close');
    // The page's own Y axis boxes go by the same rule as the editor's limits: a limit that leaves no range, or text that is no
    // number, is not used and marks its box (no empty plot from 5 to 5.5).
    const yBoxes = () => page.evaluate(() => { const g = window.calciumStep4.Plot.geom; return ['y-min', 'y-max'].map((id) => { const box = document.getElementById(id); return box.getAttribute('aria-invalid') === 'true' ? box.title : ''; }).concat([[g.lo, g.hi]]); });
    const pageY = (await yBoxes())[2];
    await page.fill('#y-min', '5'); assert.deepEqual(await yBoxes(), ['Minimum must be below maximum', '', pageY]);
    await page.fill('#y-max', '9'); assert.deepEqual(await yBoxes(), ['', '', [5, 9]]);
    await page.fill('#y-max', '1'); assert.deepEqual(await yBoxes(), ['Minimum must be below maximum', 'Maximum must be above minimum', pageY]);
    await page.fill('#y-min', ''); await page.fill('#y-max', '1e1'); assert.deepEqual(await yBoxes(), ['', 'Not a number', pageY]);
    await page.click('#y-auto'); assert.deepEqual(await yBoxes(), ['', '', pageY]);
    for (const name of await cellsOffered()) {   // a cell that has a ΔF/F₀ trace
      await openCell(name); if (await page.evaluate(() => { const S = window.calciumStep4.S; return !!S.cell && !!S.results.traces[S.cell] && S.plotMode === 'cell'; })) break;
    }
    await page.locator('#export-plot > summary').click(); await page.click('#edit-plot'); await previewReady();
    assert.deepEqual(await page.evaluate(() => ['cellWidth', 'colourCell', 'levels', 'colourMean'].map((id) => !!document.getElementById('figure-' + id))), [true, true, true, true], 'with one cell on screen its line, colour and thresholds can be set');
    shown = await previewNow(); await page.uncheck('#figure-levels'); await page.waitForFunction((was) => document.getElementById('figure-preview').toDataURL() !== was, shown);
    // The marks of the cell's call, each with a box of its own name: the response band, the required rise, and "Thresholds", which is
    // the page's own box (ticking it here ticks it there). A box takes its entry out of the legend with its mark, and no other.
    assert.deepEqual(await page.evaluate(() => ['run', 'levels', 'thresholds'].map((id) => { const box = document.getElementById('figure-' + id); return box.closest('.row').querySelector('span').textContent + ': ' + box.checked; })), ['Response: true', 'Required rise: false', 'Thresholds: false']);
    shown = await previewNow(); await page.check('#figure-thresholds'); await page.waitForFunction((was) => !window.calciumStep4.Fig.queued && document.getElementById('figure-preview').toDataURL() !== was, shown);
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4; return [f.S.thresholds, document.getElementById('thresholds').checked, f.Fig.drawn.items.map((i) => i.text).filter((text) => /rise|threshold|^response$/.test(text))]; }), [true, true, ['noise threshold', 'response']]);
    await page.uncheck('#figure-thresholds'); assert.equal(await page.evaluate(() => window.calciumStep4.S.thresholds), false);
    await page.check('#figure-thresholds'); await page.click('#figure-reset');   // Reset to defaults unticks it, here and on the page, and ticks Required rise again
    assert.deepEqual(await page.evaluate(() => [window.calciumStep4.S.thresholds, document.getElementById('thresholds').checked, document.getElementById('figure-thresholds').checked, document.getElementById('figure-levels').checked]), [false, false, false, true]);
    await page.click('#figure-close'); await page.locator('#plot').click({ position: { x: 5, y: 5 } }); await page.keyboard.press('Escape');
    assert.equal(await page.evaluate(() => window.calciumStep4.S.cell), null);
    await page.screenshot({ path: path.join(out, 'page-pending.png') });
    // Restoring applied settings includes uploaded timestamp contents and the exact protocol.
    await page.click('#discard'); await page.waitForFunction(() => !window.calciumStep4.S.state.settingsPending);
    assert.equal(await page.evaluate(() => window.calciumStep4.exportMetadata().settingsPending), false);
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.collectDraft()), await page.evaluate(() => window.calciumStep4.S.results.appliedDraft));
    for (const [id, value] of [['baseline_start','0'],['baseline_end','5000'],['baseline_start','8']]) {
      await page.fill('#'+id,value);assert.ok(await page.isDisabled('#analyze'));await page.click('#discard');
    }
    await page.locator('#windows-editor details summary').first().click();
    await page.fill('#window-0-Floor_local_fraction','-1');assert.ok(await page.isDisabled('#analyze'));await page.click('#discard');
    // The stimulus onset can be taken from the frame on screen; Discard clears it and any old note.
    await page.evaluate(() => window.calciumStep4.setFrame(9));
    const afterScreen = await page.evaluate(() => window.calciumStep4.frameTimes()[window.calciumStep4.S.frame]);   // the first frame after the one on screen
    await page.locator('#windows-editor .from-frame').first().click();
    assert.equal(await page.inputValue('#stim-frame-0'), '9'); assert.equal(Number(await page.inputValue('#window-0-Stimulus_onset_s')), afterScreen);
    await page.waitForFunction(() => window.calciumStep4.S.state.settingsPending);
    await page.evaluate(() => { window.calciumStep4.S.note = 'an old message'; });
    await page.click('#discard'); await page.waitForFunction(() => !window.calciumStep4.S.state.settingsPending);
    assert.equal(await page.evaluate(() => window.calciumStep4.S.note), null);
    assert.match(await page.title(), /^6 cells · Step 4 · Fixture recording · /);
    // The page loaded anew with a format, a resolution and a figure width kept from before. The page's own lists show them from the
    // start, so the editor, which shows those lists, opens on the format that is in use; it says that kept settings are in use; and
    // Reset to defaults clears them all.
    await openEditor('plot'); await page.selectOption('#figure-save-format', 'tiff'); await page.selectOption('#figure-save-dpi', '600'); await typeBox('width', '6');
    await page.waitForFunction(async () => { const kept = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs; return kept['step4.format'] === 'tiff' && kept['step4.dpi'] === '600' && /"width":6/.test(kept['step4.figure.plot'] || ''); });
    await page.reload(); await page.waitForFunction(() => window.calciumStep4 && window.calciumStep4.S.results && !document.getElementById('edit-plot').disabled && !document.getElementById('plot-card').hidden);
    assert.deepEqual(await page.evaluate(() => ['save-format', 'save-dpi'].map((id) => document.getElementById(id).value)), ['tiff', '600'], 'the page\'s lists show what was kept before Options is ever opened');
    await openEditor('plot');
    assert.deepEqual([await page.inputValue('#figure-save-format'), await page.inputValue('#figure-save-dpi'), await page.isVisible('#figure-save-compression'), await page.inputValue('#figure-width'), await page.isVisible('#figure-earlier')], ['tiff', '600', true, '6', true]);
    assert.match(await page.textContent('#figure-pixels'), /^3600 × 2400 pixels/, 'and the size in pixels is the kept resolution\'s');
    await page.click('#figure-reset'); await previewReady();
    assert.deepEqual([await page.inputValue('#figure-save-format'), await page.inputValue('#figure-save-dpi'), await page.isVisible('#figure-save-compression'), await page.inputValue('#figure-width'), await page.isVisible('#figure-earlier')], ['png', '300', false, '7', false]);
    await page.waitForFunction(async () => { const kept = (await (await fetch('api/session', { cache: 'no-store' })).json()).prefs; return !kept['step4.figure.plot'] && kept['step4.format'] === 'png' && kept['step4.dpi'] === '300'; });
    await page.click('#figure-close');
    // ---------------------------------------------------------------- the checks on the calls, and the review marks
    // The checks run in Fiji once the engine is done. This fixture has no ring and core values beside it, so they flagged no cell
    // and say why once, in the list of flagged cells; no cell's panel says it.
    refusals = true;   // the refusals below are answered with 400, which the browser logs
    const sha = (file) => require('crypto').createHash('sha256').update(fs.readFileSync(file)).digest('hex');
    const uncheckedNote = 'Not checked for a passing object: the raw gray values and backgrounds are not beside the analysis.';
    await page.selectOption('#group', '0'); await openCell('Cell001');
    assert.deepEqual((await page.locator('#quality-list li').allTextContents()).filter((t) => t === uncheckedNote), [uncheckedNote], 'why the cells were not checked: once, in the flagged list');
    assert.ok(await page.isHidden('#cell-checks'), 'and not in a cell\'s panel, since no rings were measured');
    // When the checks' note is all there is to say, the chip says it, and the note leads the list. (noteAlone: the chip and the
    // list with no flag of the engine, no cell below zero and, when asked, no flag of the checks)
    const noteAlone = (withoutCheckFlags) => page.evaluate((withoutCheckFlags) => { const f = window.calciumStep4, S = f.S, kept = { quality: S.quality, results: S.results }, plain = (rows) => Object.fromEntries(Object.entries(rows).map(([c, row]) => [c, Object.assign({}, row, { QC_flags: '' })]));
      try { S.quality = {}; S.results = Object.assign({}, kept.results, { calls: plain(kept.results.calls), callsByWindow: Object.fromEntries(Object.entries(kept.results.callsByWindow).map(([k, rows]) => [k, plain(rows)])),
          traces: Object.fromEntries(Object.entries(kept.results.traces).map(([c, t]) => [c, t.map((v) => Math.max(v, 0))])), checks: withoutCheckFlags ? Object.assign({}, kept.results.checks, { cells: {} }) : kept.results.checks });
        f.refreshReview(); return [document.getElementById('quality-summary').textContent, document.getElementById('quality-lead').textContent, document.getElementById('quality-list').textContent];
      } finally { S.quality = kept.quality; S.results = kept.results; f.refreshReview(); } }, withoutCheckFlags);
    assert.deepEqual(await noteAlone(false), ['Passing objects not checked', uncheckedNote, ''], 'the note leads the list and is not said again under it');
    // An analysis made before the checks has none: nothing is said about them, and nothing fails.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, kept = f.S.results.checks; try { f.S.results.checks = null; f.refreshReview(); f.pick('Cell001'); return [document.getElementById('cell-checks').hidden, document.getElementById('quality-list').textContent.includes('passing object')]; } finally { f.S.results.checks = kept; f.refreshReview(); f.pick('Cell001'); } }), [true, false]);

    // Real flags: this analysis's checks table holds them (the fixture server puts checks.csv in its place after the real checks),
    // and the page reads them as Fiji reads a saved table. Cell001 responded and an object passed through it; Cell002 did not
    // respond, with a rise on its edge; Cell005 responded while the field moved, in short separate crossings.
    const checksHeader = fs.readFileSync(path.join(out, 'Fixture recording_call_checks_browser-results.csv'), 'utf8').split(/\r?\n/)[0];
    assert.equal(checksHeader, 'Cell,Call,Flag,Rule,Event_first_frame,Event_last_frame,Detail,Near_ring_sectors_usable,Area_pixels,Ring_measured,Field_movement_measured,Passing_object_check,Not_checked_because,Crossings,Longest_crossing_samples,QC_flags');
    const passDetail = 'seen outside at the left at frame 8, 1.4 x the cell\'s rise, and at the right at frame 13, 1.2 x', edgeDetail = 'core rises 0.20 x the whole outline, ring 1.6 x the rim';
    const movedDetail = 'the field moved 2.1 px between frames 3 and 5; otherwise possible passing object, frames 4 to 7: seen outside at the top at frame 3, 1.1 x the cell\'s rise, and at the bottom at frame 8, 1.0 x; 3 crossings, the longest 2 samples';
    const checksTable = (cell002) => [checksHeader, 'Cell001,Responded,"possible passing object, frames 9 to 12",passes through,9,12,"' + passDetail + '",8,400,yes,yes,checked,,2,6,POSSIBLE_PASSING_OBJECT', cell002,
      'Cell003,Excluded,,,,,,8,400,yes,yes,not checked,no call,,,', 'Cell004,Excluded,,,,,,8,400,yes,yes,not checked,no call,,,',
      'Cell005,Responded,"field moved, frames 3 to 5; short separate crossings",field moved; short separate crossings,4,7,"' + movedDetail + '",8,400,yes,yes,checked,,3,2,FIELD_MOVED; SHORT_SEPARATE_CROSSINGS', 'Cell006,Unscorable,,,,,,8,400,yes,yes,not checked,no call,,,'].join('\n') + '\n';
    fs.writeFileSync(path.join(out, 'checks.csv'), checksTable('Cell002,Did not respond,"rise on the left edge, not in the cell, frames 20 to 24",on the edge,20,24,"' + edgeDetail + '",7,400,yes,yes,checked,,0,0,RISE_ON_EDGE'));
    fs.writeFileSync(path.join(out, 'slow'), '');   // the engine and the checks take a while, to be seen running
    assert.equal(await text('#mark'), 'Mark as not a responder'); const runFrom = await page.evaluate(() => window.calciumStep4.S.version);
    await page.check('#confirmed'); await page.click('#analyze');
    // While the engine runs, Cancel stops it; while the checks run it is off and says why. No mark is set while an analysis runs.
    await page.waitForFunction(() => { const s = window.calciumStep4.S.state; return s.busy && s.cancellable; }, null, { polling: 50 });
    assert.ok(await page.isVisible('#cancel') && await page.isEnabled('#cancel'), 'Cancel stops the engine');
    assert.ok(await page.isHidden('#mark'), 'no mark button while an analysis runs');
    const post = (body) => page.evaluate((body) => fetch('api/review', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) }).then(async (r) => [r.status, (await r.json()).error || null]), body);
    assert.deepEqual(await post({ folder: 'browser-results', cell: 'Cell001', call: 'Did not respond', reason: 'Noise' }), [400, 'Wait for the analysis to finish.']);
    await page.waitForFunction(() => { const s = window.calciumStep4.S.state; return s.busy && !s.cancellable; }, null, { polling: 50 });
    assert.ok(await page.isVisible('#cancel') && await page.isDisabled('#cancel'), 'Cancel is off while the checks run');
    assert.equal(await page.getAttribute('#cancel', 'title'), 'The checks run to the end.');
    await page.waitForFunction((v) => { const S = window.calciumStep4.S; return !S.state.busy && S.version > v && !!S.results.checks && S.results.checks.cells.Cell001.flags.length === 1; }, runFrom, { timeout: 15000 });
    fs.unlinkSync(path.join(out, 'slow'));
    assert.ok(await page.isHidden('#cancel'), 'Cancel belongs to a running analysis');
    // The flagged list: each check flag in its text as given, beside the engine's own flags; a filter for each kind, by its rule.
    const flaggedRows = () => page.locator('#table-flagged tbody tr').evaluateAll((rows) => rows.map((r) => [...r.cells].map((c) => c.textContent)));
    assert.deepEqual(await flaggedRows(), [['Cell001', 'Responded', 'possible passing object, frames 9 to 12'], ['Cell002', 'Did not respond', 'rise on the left edge, not in the cell, frames 20 to 24'],
      ['Cell005', 'Responded', 'Low F₀ · Below zero · field moved, frames 3 to 5 · short separate crossings']], 'each check flag joins the list, in its own words');
    assert.deepEqual(await page.locator('#quality-list .flag-filter').allTextContents(), ['Possible passing object: 1 cell, 1 responder', 'Rise on the edge: 1 cell, 0 responders', 'Low F₀: 1 cell, 1 responder', 'Below zero: 1 cell, 1 responder',
      'Field moved: 1 cell, 1 responder', 'Short separate crossings: 1 cell, 1 responder'], 'a filter for each kind of check flag, named by its rule');
    assert.deepEqual([await text('#quality-summary'), await text('#quality-lead')], ['3 cells flagged', '3 cells flagged, including 2 responders.']);
    assert.ok(!(await text('#quality-list')).includes('Not checked'), 'with the rings and the field measured there is no note');
    await page.locator('#quality-review summary').click();
    // the list grows with its filters: at this size nothing in it is cut off, and nothing scrolls out of sight
    assert.deepEqual(await page.evaluate(() => [...document.querySelectorAll('.quality-pop, #quality-list')].map((e) => e.scrollHeight <= e.clientHeight + 1)), [true, true], 'every filter and note of the flagged list is shown whole');
    // The table of flagged cells is never cut: with more rows than this fixture has (12, as on a real recording with more than 7),
    // every row is inside the table's own box; in a short window the whole list scrolls instead.
    assert.deepEqual(await page.evaluate(() => { const box = document.getElementById('quality-cells'), body = document.getElementById('table-flagged').tBodies[0], first = body.rows[0];
      try { while (body.rows.length < 12) body.appendChild(first.cloneNode(true));
        const edge = box.getBoundingClientRect(); return [body.rows.length, box.scrollHeight <= box.clientHeight + 1, [...body.rows].every((row) => { const r = row.getBoundingClientRect(); return r.top >= edge.top - 0.5 && r.bottom <= edge.bottom + 0.5; })];
      } finally { window.calciumStep4.refreshReview(); } }), [12, true, true], 'twelve flagged cells: every row of the table is shown');
    await page.locator('#quality-list .flag-filter', { hasText: 'Possible passing object' }).click();
    assert.equal(await page.inputValue('#group'), '5'); assert.deepEqual(await cellsOffered(), ['Cell001']);
    await page.evaluate(() => { document.getElementById('quality-review').open = false; }); await page.selectOption('#group', '0');
    // The cell's panel: its check flags under its call; a flag goes to the first frame it names, one that names none to the
    // first frame of the event.
    const checksLine = () => page.evaluate(() => { const box = document.getElementById('cell-checks'); return { text: box.hidden ? '' : box.textContent, links: [...box.querySelectorAll('button.jump')].map((b) => [b.textContent, b.title]) }; });
    await openCell('Cell001');
    assert.deepEqual(await checksLine(), { text: 'possible passing object, frames 9 to 12', links: [['possible passing object, frames 9 to 12', 'Go to frame 9\n' + passDetail]] });
    assert.ok(await page.isHidden('#cell-quality'), 'the engine\'s flags keep a line of their own');
    await openCell('Cell005');
    assert.deepEqual(await checksLine(), { text: 'field moved, frames 3 to 5 · short separate crossings', links: [['field moved, frames 3 to 5', 'Go to frame 3\n' + movedDetail], ['short separate crossings', 'Go to frame 4\n' + movedDetail]] },
      'field moved goes to the first frame it names, short separate crossings, which names none, to the event\'s first frame');
    assert.match(await text('#cell-quality'), /^Low F₀ relative to baseline noise/); assert.ok(!(await text('#cell-quality')).includes('field moved'));
    await page.evaluate(() => window.calciumStep4.setFrame(1)); await page.locator('#cell-checks .jump').first().click(); await waitFrame(3);
    await openCell('Cell001'); await page.selectOption('#view', 'picture'); await page.evaluate(() => window.calciumStep4.setFrame(1));
    await page.locator('#cell-checks .jump').first().click(); await waitFrame(9);
    assert.deepEqual(await page.evaluate(() => [window.calciumStep4.S.view, document.getElementById('view').value]), ['raw', 'raw'], 'from the colour overview to the recording, at the event\'s first frame');
    await page.selectOption('#view', 'flat'); await page.evaluate(() => window.calciumStep4.setFrame(2)); await page.locator('#cell-checks .jump').first().click(); await waitFrame(9);
    assert.equal(await page.evaluate(() => window.calciumStep4.S.view), 'flat', 'a view of the recording stays as it is');
    await page.selectOption('#view', 'raw');
    // The frame is the first of the flag's frames as Fiji gives them; from an answer without them, the first frame its text names.
    // A flag whose frames are none goes to the event's first frame, and without that either it is text.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, cells = f.S.results.checks.cells, one = cells.Cell001, five = cells.Cell005;
      const kept = [one.firstFrame, JSON.stringify(one.flags), five.firstFrame, JSON.stringify(five.flags)], seen = [];
      const links = (name) => { f.pick(name); return [...document.getElementById('cell-checks').children].map((e) => e.tagName === 'BUTTON' ? e.title.split('\n')[0] : 'text: ' + e.textContent); };
      try {
        for (const flag of one.flags.concat(five.flags)) delete flag.frames;
        one.firstFrame = 11; seen.push(links('Cell001'));
        one.flags[0].frames = [10, 12]; seen.push(links('Cell001'));
        one.flags[0].frames = null; seen.push(links('Cell001'));
        one.firstFrame = null; seen.push(links('Cell001'));
        delete one.flags[0].frames; seen.push(links('Cell001'));
        five.firstFrame = null; seen.push(links('Cell005'));
        return seen;
      } finally { one.firstFrame = kept[0]; one.flags = JSON.parse(kept[1]); five.firstFrame = kept[2]; five.flags = JSON.parse(kept[3]); f.pick('Cell001'); } }),
      [['Go to frame 9'], ['Go to frame 10'], ['Go to frame 11'], ['text: possible passing object, frames 9 to 12'], ['Go to frame 9'], ['Go to frame 3', 'text: short separate crossings']],
      'from the text without frames, Fiji\'s frames first, the event\'s first frame for none, text without either');

    // The codes of the checks join QC_flags, after the engine's own: in the calls table on the rows of the primary stimulus, in
    // Cell quality on every cell's row, each joined as that table joins its own; and so in their Exports, whose record says so.
    // The saved tables are not changed.
    const savedCsv = (note) => { const saved = note.match(/^Saved in plots: (\S+\.csv)/)[1], lines = fs.readFileSync(path.join(out, 'plots', saved), 'utf8').trim().split(/\r?\n/).map((l) => l.split(','));
      return { lines, at: (column, ...first) => lines.find((l) => first.every((v, i) => l[lines[0].indexOf(i ? 'Role' : 'Cell')] === v))[lines[0].indexOf(column)], record: JSON.parse(fs.readFileSync(path.join(out, 'plots', saved.replace(/\.csv$/, '.metadata.json')), 'utf8')) }; };
    const codesNote = 'QC_flags includes the codes of the call checks, after the engine\'s codes.';
    await page.click('.tab[data-tab=calls]');
    { const shown = await tableNow('table-calls'), flags = (cell, role) => at(shown, cell, 'Flags', role);
      assert.deepEqual([flags('Cell001', 'Stimulus'), flags('Cell001', 'Control'), flags('Cell002', 'Stimulus').text, flags('Cell005', 'Stimulus').text, flags('Cell003', 'Stimulus').text],
        [{ text: 'Duration under a frame · Possible passing object', tip: 'DURATION_BELOW_SAMPLE_INTERVAL; POSSIBLE_PASSING_OBJECT' }, { text: 'Duration under a frame', tip: 'DURATION_BELOW_SAMPLE_INTERVAL' },
          'Duration under a frame · Rise on the edge', 'Duration under a frame · Field moved · Short separate crossings', 'Duration under a frame'], 'the Flags column: the codes of the checks after the engine\'s, on the primary stimulus only'); }
    await page.locator('#export-calls > summary').click(); await page.click('#calls-save'); await page.waitForFunction(() => /^Saved in plots: /.test(document.getElementById('calls-saved').textContent));
    { const saved = savedCsv(await text('#calls-saved'));
      assert.deepEqual([saved.at('QC_flags', 'Cell001', 'Stimulus'), saved.at('QC_flags', 'Cell001', 'Control'), saved.at('QC_flags', 'Cell005', 'Stimulus'), saved.at('QC_flags', 'Cell003', 'Stimulus'), saved.record.qc_flags],
        ['DURATION_BELOW_SAMPLE_INTERVAL; POSSIBLE_PASSING_OBJECT', 'DURATION_BELOW_SAMPLE_INTERVAL', 'DURATION_BELOW_SAMPLE_INTERVAL; FIELD_MOVED; SHORT_SEPARATE_CROSSINGS', 'DURATION_BELOW_SAMPLE_INTERVAL', codesNote], 'the calls Export holds the codes in QC_flags, and its record says so'); }
    assert.deepEqual(await page.evaluate(() => window.calciumStep4.S.results.tables.calls.rows.map((row) => row[window.calciumStep4.S.results.tables.calls.header.indexOf('QC_flags')])), Array(12).fill('DURATION_BELOW_SAMPLE_INTERVAL'), 'the saved calls table as it was');
    await page.click('.tab[data-tab=metrics]'); await page.selectOption('#metrics-table', 'quality');
    { const shown = await tableNow('table-metrics');
      assert.deepEqual([at(shown, 'Cell005', 'Flags'), at(shown, 'Cell001', 'Flags').text, at(shown, 'Cell003', 'Flags').text],
        [{ text: 'Low F₀ · Field moved · Short separate crossings', tip: 'reference_below_3_baseline_sd_review;FIELD_MOVED;SHORT_SEPARATE_CROSSINGS' }, 'Possible passing object', '–'], 'Cell quality: the codes after the engine\'s, joined as its own'); }
    await page.locator('#export-metrics > summary').click(); await page.click('#metrics-save'); await page.waitForFunction(() => /^Saved in plots: /.test(document.getElementById('metrics-saved').textContent));
    { const saved = savedCsv(await text('#metrics-saved'));
      assert.deepEqual([saved.at('QC_flags', 'Cell005'), saved.at('QC_flags', 'Cell001'), saved.at('QC_flags', 'Cell003'), saved.record.qc_flags],
        ['reference_below_3_baseline_sd_review;FIELD_MOVED;SHORT_SEPARATE_CROSSINGS', 'POSSIBLE_PASSING_OBJECT', '', codesNote], 'the Cell quality Export too'); }
    // a table without QC_flags gets no codes, and its record says nothing of them
    await page.selectOption('#metrics-table', 'metrics'); await page.locator('#export-metrics > summary').click(); await page.click('#metrics-save'); await page.waitForFunction(() => /^Saved in plots: \S+table_metrics/.test(document.getElementById('metrics-saved').textContent));
    assert.ok(!('qc_flags' in savedCsv(await text('#metrics-saved')).record), 'no code added, no word of codes in the record');
    await page.selectOption('#metrics-table', 'summary');
    // The window metrics take the codes on the primary stimulus's rows, the per-cell metrics on every row, each with its own
    // separator; no other table, no table without QC_flags, and nothing without checks.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, kept = f.S.results.checks;
      const w = f.withCodes('windowMetrics', { header: ['Cell', 'Role', 'Stimulus', 'QC_flags'], rows: [['Cell001', 'Stimulus', 'Capsaicin', 'A'], ['Cell001', 'Control', 'Ionomycin', ''], ['Cell005', 'Stimulus', 'Capsaicin', '']] });
      const m = f.withCodes('metrics', { header: ['Cell', 'QC_flags'], rows: [['Cell001', ''], ['Cell005', 'B'], ['Cell003', '']] });
      const other = f.withCodes('counts', { header: ['Cell', 'QC_flags'], rows: [['Cell001', '']] }).added, noColumn = f.withCodes('metrics', { header: ['Cell', 'Peak'], rows: [['Cell001', '1']] }).added;
      try { f.S.results.checks = null; const calls = f.S.results.tables.calls, none = f.withCodes('calls', calls);
        return [w.added, w.data.rows.map((r) => r[3]), m.added, m.data.rows.map((r) => r[1]), other, noColumn, none.added, none.data === calls, JSON.stringify(f.callsRows()) === JSON.stringify([calls.header].concat(calls.rows))];
      } finally { f.S.results.checks = kept; } }),
      [true, ['A; POSSIBLE_PASSING_OBJECT', '', 'FIELD_MOVED; SHORT_SEPARATE_CROSSINGS'], true, ['POSSIBLE_PASSING_OBJECT', 'B;FIELD_MOVED;SHORT_SEPARATE_CROSSINGS', ''], false, false, false, true, true], 'which tables take the codes, and how');
    await page.click('.tab[data-tab=review]');

    // ---------------------------------------------------------------- around the cell
    // What the passing-object check used around the selected cell, asked of Fiji cell by cell. This fixture has no ring files
    // beside it, so Fiji says why there is nothing to show: the button is off, and its tooltip says it.
    const aroundButton = () => page.evaluate(() => { const b = document.querySelector('[data-plot=around]'); return [b.disabled, b.title]; });
    const aroundTitle = 'The cell, its core and the ring around it, in gray levels';
    const AROUND_STYLE = ['zoneWidth', 'colourBefore', 'colourAfter', 'colourCore', 'colourAgain', 'otherWidth', 'colourOtherSectors', 'flagged', 'colourFlagged'];   // the settings of this plot alone
    await openCell('Cell001'); await page.waitForFunction(() => window.calciumStep4.Around.got.has('Cell001'));
    const realWhy = await page.evaluate(() => window.calciumStep4.Around.got.get('Cell001').why);
    assert.deepEqual(await aroundButton(), [true, realWhy], 'Fiji\'s sentence in the tooltip of the button, which is off'); assert.match(realWhy, /^[A-Z][^()]*\.$/, 'a sentence: ' + realWhy);
    // Fiji's answers for six cells: rule A on Cell001 (the sectors seen before and after, with their frames, one in the far
    // ring, and one sector of the ring not usable), rule B on Cell002 (rim, ring, core; a band that ends at the pause of the
    // fixture's clock), rule C on Cell005 (the sectors on the start and end sides without frames, and the core; a band of one
    // frame), no flag on Cell006 and Cell004 (the core and the eight sectors of the ring, no band) and a cell with nothing to show.
    const fixtureClock = Array.from({ length: 30 }, (_, i) => i / 2 + (i >= 15 ? 5 : 0)), sides = ['right', 'lower right', 'bottom', 'lower left', 'left', 'upper left', 'top', 'upper right'];
    const placeOf = (z) => z === 16 ? 'core' : (z < 8 ? 'ring, ' : z < 16 ? 'far ring, ' : 'rim, ') + sides[z < 16 ? z % 8 : z - 17];
    const aroundOf = (cell, flag, lines, others, band, unusable = []) => ({ folder: 'browser-results', cell, name: cell, why: '', time: fixtureClock, rise: fixtureClock.map((t, i) => i >= 8 && i <= 11 ? 40 : 0),
      zones: Array.from({ length: 25 }, (_, z) => ({ zone: z, key: 'zone_' + z, place: placeOf(z), usable: !unusable.includes(z), trace: unusable.includes(z) ? null : fixtureClock.map((t, i) => z - 8 + (i === 10 ? 30 : 0)) })),
      eventFirst: 9, eventLast: 12, rule: flag ? flag.rule : '', text: flag ? flag.text : '', flags: flag ? [flag] : [], codes: flag ? [flag.code] : [],
      lines: lines.map(([zone, part, frame]) => ({ zone, place: placeOf(zone), part, frame })), others, band });
    const passing = (rule, text) => ({ rule, text, code: rule === 'on the edge' ? 'RISE_ON_EDGE' : 'POSSIBLE_PASSING_OBJECT' });
    const faked = { Cell001: aroundOf('Cell001', passing('passes through', 'possible passing object, frames 9 to 12'), [[1, 'before', 8], [13, 'after', 13]], [0, 2, 3, 4, 5, 6, 7], [9, 12], [3]),
      Cell002: aroundOf('Cell002', passing('on the edge', 'rise on the lower right edge, not in the cell, frames 7 to 15'), [[18, 'rim', null], [1, 'ring', null], [16, 'core', null]], [0, 2, 3, 4, 5, 6, 7], [7, 15]),
      Cell005: aroundOf('Cell005', passing('crosses the cell', 'possible passing object, frames 5 to 5'), [[4, 'before', null], [0, 'after', null], [16, 'core', null]], [1, 2, 3, 5, 6, 7], [5, 5]),
      Cell006: aroundOf('Cell006', null, [[16, 'core', null]], [0, 1, 2, 3, 4, 5, 6, 7], null), Cell004: aroundOf('Cell004', null, [[16, 'core', null]], [0, 1, 2, 3, 4, 5, 6, 7], null),
      Cell003: { folder: 'browser-results', cell: 'Cell003', name: 'Cell003', why: 'Passing object not checked: no call.' } };
    // aroundHeld: an answer kept on its way until the test lets it go; aroundLost: a cell whose next answer is lost, as when Fiji
    // is not connected
    const aroundAsked = [], aroundHeld = {}, aroundLost = new Set();
    await page.route('**/api/around*', async (route) => { const cell = new URL(route.request().url()).searchParams.get('cell'); aroundAsked.push(cell);
      if (aroundLost.has(cell)) { aroundLost.delete(cell); return route.abort('failed'); }
      if (aroundHeld[cell]) await aroundHeld[cell];
      return route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(faked[cell] || { cell, why: 'There is no cell ' + cell + ' in this analysis.' }) }); });
    await page.evaluate(() => window.calciumStep4.loadResults()); await openCell('Cell001');   // the results read again: what was asked before is asked again
    await page.waitForFunction(() => { const a = window.calciumStep4.Around.got.get('Cell001'); return !!a && !a.why; });
    assert.deepEqual(await aroundButton(), [false, aroundTitle], 'with something to show the button is on, with its own tooltip');
    await page.click('[data-plot=around]');
    const aroundShown = () => page.evaluate(() => { const f = window.calciumStep4, a = f.Plot.around, g = f.Plot.geom;
      return { legend: g.items.map((i) => i.text), inks: g.items.map((i) => i.colour), lines: a.lines.map((l) => [l.place, l.part, l.frame, f.zoneInk(f.STYLE.plot, l)]), others: a.others.map((o) => o.zone), band: a.band, frames: a.frames,
        hidden: [document.getElementById('signal-label').hidden, document.getElementById('every-label').hidden], label: document.getElementById('plot').getAttribute('aria-label') }; });
    // The colours of the parts, the same in the Fiji window: before and rim, after and ring, the core; a colour met a second time
    // for one cell takes the next of four. The other ring sectors and the flagged frames in their own; none is a call's colour.
    const BEFORE = '#0072b2', AFTER = '#aa3377', CORE = '#009e73', AGAIN = ['#e69f00', '#56b4e9', '#cc79a7', '#882255'], OTHERS = '#aeb8c2', FLAGGED = '#e2d6f3', CALL_INKS = ['#177465', '#c8660e', '#8a94a0'];
    { const shown = await aroundShown();
      assert.deepEqual([shown.legend, shown.lines, shown.others, shown.band, shown.frames], [['Cell001', 'ring, lower right', 'far ring, upper left', 'other ring sectors', 'flagged frames'],
        [['ring, lower right', 'before', 8, BEFORE], ['far ring, upper left', 'after', 13, AFTER]], [0, 2, 4, 5, 6, 7], [3.75, 5.75], [9, 12]], 'rule A: the cell, the sector seen before and the one seen after by their places, the other usable ring sectors once, the flagged frames half a frame step to each side');
      assert.deepEqual(shown.inks, ['#1f2933', BEFORE, AFTER, OTHERS, FLAGGED], 'the legend: the cell, the zones of its flag, the other ring sectors, the band, each in its colour');
      assert.deepEqual(shown.hidden, [true, true], 'Signal and All traces belong to the other plots');
      assert.match(shown.label, /^Plot of ΔF \(gray levels\) against time: Cell001; ring, lower right; far ring, upper left; other ring sectors; flagged frames$/, 'a screen reader is told what it shows'); }
    { const inks = [BEFORE, AFTER, CORE, ...AGAIN, OTHERS, FLAGGED, '#1f2933'].concat(CALL_INKS); assert.equal(new Set(inks).size, inks.length, 'the zones, the other sectors, the band, the cell and the calls in colours of their own'); }
    // any mix of parts: a colour met again for the cell takes the next of AGAIN, in order, the last for any more
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, time = [0, 1, 2], trace = [0, 1, 0];
      const a = { time, rise: trace, zones: Array.from({ length: 25 }, (_, z) => ({ zone: z, place: 'zone ' + z, usable: true, trace })), others: [],
        lines: [[1, 'before'], [2, 'rim'], [3, 'after'], [16, 'core'], [4, 'ring'], [5, 'before'], [6, 'other'], [7, 'after']].map(([zone, part]) => ({ zone, part, frame: null })) };
      return f.aroundParts(a).lines.map((l) => f.zoneInk(f.STYLE.plot, l)); }), [BEFORE, AGAIN[0], AFTER, CORE, AGAIN[1], AGAIN[2], AGAIN[3], AGAIN[3]], 'parts met again');
    // Nothing is written on it but the axis titles, the tick labels, "stimulus" and the legend; its entries and its marks go one to one.
    const cleanAround = (seen, view) => {
      assert.deepEqual(wordsOf(seen), ['Time (s)', 'stimulus', 'ΔF (gray levels)'], view + ': outside the legend, only the axis titles, the tick labels and "stimulus" are written');
      assert.deepEqual(seen.texts.filter((t) => inLegend(seen, t.x, t.y)).map((t) => t.text), seen.wrapped.flat(), view + ': the legend holds its entries and nothing else');
      assert.deepEqual(marksOf(seen), entriesOf(seen), view + ': legend entries and drawn marks, one to one');
      assert.deepEqual(seen.lines.concat(seen.boxes).filter((m) => CALL_INKS.includes(String(m.ink).toLowerCase())), [], view + ': nothing in a colour of a call'); return seen; };
    { const onPage = cleanAround(await drawing(false), 'around the cell, rule A, on the page'), fig = cleanAround(await drawing(true), 'around the cell, rule A, in the figure');
      assert.ok(onPage.lo <= Math.min(0, 1 - 8) && onPage.hi >= 13 - 8 + 30, 'the Y axis holds every line: ' + [onPage.lo, onPage.hi]); assert.deepEqual(fig.legend, onPage.legend);
      assert.deepEqual([fig.x0, fig.x1], [0, 19.5], 'the whole recording'); }
    // the numbers at the frame of the picture, in gray levels with one decimal: the cell's rise and each zone of the flag by its place
    assert.equal(await page.evaluate(() => { const f = window.calciumStep4, was = f.S.frame; f.setFrame(11, true); const said = f.Plot.readout; f.setFrame(was, true); return said; }),
      '5.0 s · Cell001 40.0 · ring, lower right 23.0 · far ring, upper left 35.0', 'the numbers of the cell and the zones of its flag at the frame shown');
    // Export: the picture and its numbers under the name of this plot and cell; the record says what was drawn.
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, rows = f.plotNumbers().trim().split('\n'); return [f.plotKind(), rows[0], rows[11], rows.length]; }),
      ['around_the_cell', 'Time_s,Cell001,ring_lower_right,far_ring_upper_left,ring_right,ring_bottom,ring_left,ring_upper_left,ring_top,ring_upper_right', '5,40,23,35,22,24,26,27,28,29', 31], 'the values: time, the cell, each zone of its flag, each other ring sector, under their places');
    { const sent = page.waitForRequest((r) => r.url().includes('/api/export') && r.method() === 'POST'); await page.locator('#export-plot > summary').click(); await page.click('#save-plot-folder');
      const body = JSON.parse((await sent).postData()); await page.waitForFunction(() => /^Saved /.test(document.getElementById('plot-saved').textContent));
      assert.match(body.name, /^Fixture_recording_around_the_cell_Cell001_\d{8}-\d{6}$/, 'the file is named for this plot and cell');
      assert.deepEqual([body.metadata.signal, body.metadata.plotted, body.metadata.around_the_cell],['dF (gray levels)', ['Cell001', 'ring, lower right', 'far ring, upper left', 'other ring sectors', 'flagged frames'],
        { cell: 'Cell001', rule: 'passes through', flag: 'possible passing object, frames 9 to 12', units: 'gray levels', lines: [{ zone: 1, place: 'ring, lower right', part: 'before', frame: 8 }, { zone: 13, place: 'far ring, upper left', part: 'after', frame: 13 }],
          other_ring_sectors: [0, 2, 4, 5, 6, 7].map((zone) => ({ zone, place: placeOf(zone) })), flagged_frames: [9, 12], flagged_frames_s: [4, 5.5] }], 'the record: the signal in gray levels, what is plotted, and the lines, sectors and flagged frames, with the times of the first and last');
      assert.ok(AROUND_STYLE.every((key) => key in body.metadata.picture.figure), 'its figure settings include those of this plot'); }
    // The record beside a figure of another plot has none of this plot's settings: it is as it was before this plot.
    { await page.click('[data-plot=cell]'); const sent = page.waitForRequest((r) => r.url().includes('/api/export') && r.method() === 'POST'); await page.locator('#export-plot > summary').click(); await page.click('#save-plot-folder');
      const figure = JSON.parse((await sent).postData()).metadata.picture.figure, own = await page.evaluate(() => Object.keys(window.calciumStep4.STYLE.plot));
      assert.deepEqual(Object.keys(figure), ['kind'].concat(own.filter((key) => !AROUND_STYLE.includes(key))), 'One cell + average: every setting but those of Around the cell, in their order');
      await page.waitForFunction(() => /^Saved /.test(document.getElementById('plot-saved').textContent)); await page.click('[data-plot=around]'); }
    // the figure editor offers the boxes of this plot: the zone lines and the flagged frames, not the mean's or the baseline's
    await page.locator('#export-plot > summary').click(); await page.click('#edit-plot'); await previewReady();
    assert.deepEqual(await page.evaluate(() => ['zoneWidth', 'flagged', 'cellWidth', 'meanWidth', 'baseline', 'window', 'sem'].map((key) => !!document.getElementById('figure-' + key))), [true, true, true, false, false, false, false], 'the figure editor\'s boxes for this plot');
    await page.uncheck('#figure-flagged'); await page.waitForFunction(() => !window.calciumStep4.Fig.queued);
    assert.ok(!(await drawing(true)).legend.includes('flagged frames'), 'Flagged frames off: no band and no entry'); await page.click('#figure-reset'); await page.click('#figure-close');
    // The Y axis boxes are this plot's own, in gray levels: what is typed here is not the other plots', nor the other way round.
    await page.fill('#y-max', '80'); await page.waitForFunction(() => window.calciumStep4.Plot.geom.hi === 80);
    await page.click('[data-plot=cell]'); assert.deepEqual([await page.inputValue('#y-max'), await page.evaluate(() => window.calciumStep4.Plot.geom.hi !== 80)], ['', true], 'One cell + average keeps its own Y axis');
    await page.click('[data-plot=around]'); assert.deepEqual([await page.inputValue('#y-max'), await page.evaluate(() => window.calciumStep4.Plot.geom.hi)], ['80', 80], 'and Around the cell its own');
    await page.click('#y-auto'); await page.waitForFunction(() => window.calciumStep4.Plot.geom.hi !== 80);
    // The Y title and limits typed in the figure editor are kept for each kind of plot: those of the other plots come back after
    // Around the cell, and its own after them.
    const yTyped = () => page.evaluate(() => { const p = window.calciumStep4.Fig.plot; return [p.yTitle, p.yMax]; });
    const editPlot = async () => { await page.locator('#export-plot > summary').click(); await page.click('#edit-plot'); await previewReady(); };
    await page.click('[data-plot=cell]'); await editPlot(); await page.fill('#figure-yTitle', 'Signal title'); await typeBox('yMax', '3'); await page.click('#figure-close');
    await page.click('[data-plot=around]'); assert.deepEqual(await yTyped(), [null, null], 'Around the cell has none of them');
    await editPlot(); await page.fill('#figure-yTitle', 'Gray title'); await page.click('#figure-close');
    await page.click('[data-plot=cell]'); assert.deepEqual(await yTyped(), ['Signal title', 3], 'One cell + average gets its own back');
    await page.click('[data-plot=around]'); assert.deepEqual(await yTyped(), ['Gray title', null], 'and Around the cell its own');
    await page.evaluate(() => { window.calciumStep4.Fig.plot.yTitle = null; });   // its own title again, for the figures drawn below
    // Rule B: the rim sector, the ring sector outside it and the core, in that order, with the flagged frames.
    await openCell('Cell002'); await page.waitForFunction(() => !!window.calciumStep4.Plot.around && window.calciumStep4.Plot.around.lines.length === 3);
    { const shown = await aroundShown();
      assert.deepEqual([shown.legend, shown.lines, shown.band], [['Cell002', 'rim, lower right', 'ring, lower right', 'core', 'other ring sectors', 'flagged frames'],
        [['rim, lower right', 'rim', null, BEFORE], ['ring, lower right', 'ring', null, AFTER], ['core', 'core', null, CORE]], [2.75, 7.25]], 'rule B: rim, ring and core by their places; the band half the median frame step beyond its last frame, also at the pause'); }
    cleanAround(await drawing(false), 'around the cell, rule B, on the page');
    // Rule C, "crosses the cell": the ring sectors on the start and end sides, which have no frame, and the core; a band of
    // one frame is half a frame to each side of it.
    await openCell('Cell005'); await page.waitForFunction(() => !!window.calciumStep4.Plot.around && window.calciumStep4.Plot.around.lines.some((l) => l.place === 'ring, left'));
    { const shown = await aroundShown();
      assert.deepEqual([shown.legend, shown.lines, shown.band, shown.frames], [['Cell005', 'ring, left', 'ring, right', 'core', 'other ring sectors', 'flagged frames'],
        [['ring, left', 'before', null, BEFORE], ['ring, right', 'after', null, AFTER], ['core', 'core', null, CORE]], [1.75, 2.25], [5, 5]], 'rule C: the sectors on the start and end sides and the core'); }
    cleanAround(await drawing(true), 'around the cell, rule C, in the figure');
    // No flag: the cell, the core and the eight ring sectors in gray; no band. While Fiji's answer is on its way the button is
    // off, without a tooltip.
    let letGo; aroundHeld.Cell006 = new Promise((resolve) => { letGo = resolve; });
    await openCell('Cell006'); assert.deepEqual(await aroundButton(), [true, ''], 'an answer on its way: the button is off, without a tooltip'); letGo(); delete aroundHeld.Cell006;
    await page.waitForFunction(() => !!window.calciumStep4.Plot.around && window.calciumStep4.Plot.around.others.length === 8);
    assert.deepEqual(await aroundButton(), [false, aroundTitle], 'and on again with its answer');
    { const shown = await aroundShown();
      assert.deepEqual([shown.legend, shown.lines, shown.band, shown.inks], [['Cell006', 'core', 'other ring sectors'], [['core', 'core', null, CORE]], null, ['#1f2933', CORE, OTHERS]], 'no flag: the cell, the core and the ring, no band'); }
    cleanAround(await drawing(false), 'around the cell, no flag, on the page');
    // Nothing to show: the plot says why, the button is off with the same words, and the Export entries are off.
    await openCell('Cell003'); await page.waitForFunction(() => window.calciumStep4.Around.got.has('Cell003'));
    const emptyPlot = () => page.evaluate(() => [document.getElementById('plot-empty').hidden, document.getElementById('plot-empty').textContent, ...['save-plot', 'save-plot-folder', 'edit-plot'].map((id) => document.getElementById(id).disabled), window.calciumStep4.Plot.geom, document.getElementById('plot').getAttribute('aria-label')]);
    assert.deepEqual([await page.evaluate(() => window.calciumStep4.S.plotMode), await aroundButton(), await emptyPlot()], ['around', [true, 'Passing object not checked: no call.'], [false, 'Passing object not checked: no call.', true, true, true, null, 'Passing object not checked: no call.']], 'nothing to show: why, in the plot, in the tooltip and to a screen reader; nothing to export');
    assert.ok(await page.waitForFunction(() => getComputedStyle(document.querySelector('[data-plot=around]')).backgroundColor === 'rgb(183, 211, 227)', null, { timeout: 10000 }).then(() => true, () => false), 'chosen with nothing to show, the button looks off');   // after its 120 ms change of colour
    // The control window: the checks are made for the primary stimulus, so the button is off on it and says so.
    await openCell('Cell001'); await page.waitForFunction(() => !!window.calciumStep4.Plot.around);
    await page.selectOption('#response-window', 'Control::Ionomycin');
    assert.deepEqual([await aroundButton(), await emptyPlot()], [[true, 'Checks are made for the primary stimulus.'], [false, 'Checks are made for the primary stimulus.', true, true, true, null, 'Checks are made for the primary stimulus.']], 'on the control window the button is off and says why');
    await openCell('Cell002'); assert.deepEqual(await aroundButton(), [true, 'Checks are made for the primary stimulus.'], 'for every cell');
    assert.deepEqual([...new Set(aroundAsked)].sort(), ['Cell001', 'Cell002', 'Cell003', 'Cell005', 'Cell006'], 'each cell asked for once it is chosen, and kept');
    assert.equal(aroundAsked.length, 5, 'and asked once');
    await page.selectOption('#response-window', 'Stimulus::Capsaicin'); await page.waitForFunction(() => !!window.calciumStep4.Plot.around && window.calciumStep4.S.cell === 'Cell002');
    assert.deepEqual([(await aroundButton())[0], await page.isHidden('#plot-empty')], [false, true], 'back on the primary stimulus, the plot again');
    // An answer lost because Fiji was not connected is asked again once Fiji answers a poll, without choosing the cell again.
    { aroundAborted = true; aroundLost.add('Cell004'); const from = aroundAsked.length; await openCell('Cell004');
      await page.waitForFunction(() => { const a = window.calciumStep4.Around.got.get('Cell004'); return !!a && Array.isArray(a.time); }, null, { timeout: 10000 });
      assert.deepEqual([aroundAsked.slice(from), await aroundButton(), await page.isHidden('#plot-empty')], [['Cell004', 'Cell004'], [false, aroundTitle], true], 'lost once and asked again: the plot, with nothing to say about Fiji');
      aroundAborted = false; }
    // Without a cell the plot goes back to the average, and the button says why it is off; the average has its own Y title back.
    await openCell(''); assert.deepEqual([await page.evaluate(() => window.calciumStep4.S.plotMode), await aroundButton(), await page.evaluate(() => ['signal-label', 'every-label'].map((id) => document.getElementById(id).hidden)), await yTyped()],
      ['mean', [true, 'No cell selected.'], [false, false], ['Signal title', 3]], 'no cell: the average with its Signal, All traces and Y title, and the button off');
    await editPlot(); await page.click('#figure-reset'); await page.click('#figure-close');
    assert.deepEqual([await yTyped(), await page.evaluate(() => window.calciumStep4.Fig.yKept)], [[null, null], {}], 'Reset to defaults: no Y title kept for any plot');
    await page.unroute('**/api/around*'); await page.evaluate(() => window.calciumStep4.loadResults()); await page.selectOption('#group', '0');

    // The mark button says what it does; there is none for a cell with no call to review, nor on the control window.
    const markButton = async (name) => { await openCell(name); return (await page.isVisible('#mark')) ? text('#mark') : null; };
    assert.deepEqual([await markButton('Cell001'), await markButton('Cell002'), await markButton('Cell003'), await markButton('Cell006')], ['Mark as not a responder', 'Mark as a responder', null, null]);
    await page.selectOption('#response-window', 'Control::Ionomycin'); assert.equal(await markButton('Cell002'), null, 'no mark on the control window');
    // The checks are about the calls for the primary stimulus: on the control window no cell is flagged by them, as in the Fiji window.
    assert.deepEqual([await flaggedRows(), await page.locator('#quality-list .flag-filter').allTextContents(), (await page.locator('#group option').allTextContents())[5]],
      [[['Cell005', 'Did not respond', 'Low F₀ · Below zero']], ['Low F₀: 1 cell, 0 responders', 'Below zero: 1 cell, 0 responders'], 'Needs quality review (1)'], 'on the control window only the engine\'s flags are listed and counted');
    await openCell('Cell001'); assert.ok(await page.isHidden('#cell-checks'), 'and the cell panel has none');
    await page.selectOption('#response-window', 'Stimulus::Capsaicin');
    // A check's filter left on for the primary stimulus is let go on the control window, before Show counts its cells.
    await page.locator('#quality-review summary').click(); await page.locator('#quality-list .flag-filter', { hasText: 'Possible passing object' }).click();
    await page.evaluate(() => { document.getElementById('quality-review').open = false; });
    await page.selectOption('#response-window', 'Control::Ionomycin');
    assert.deepEqual([(await page.locator('#group option').allTextContents())[5], await cellsOffered()], ['Needs quality review (1)', ['Cell005']], 'a check filter left on the primary stimulus is let go on the control window');
    await page.selectOption('#response-window', 'Stimulus::Capsaicin'); await page.selectOption('#group', '0');
    // The form opens in place of the button. On a flagged responder its reason is the flag's; on a cell that did not respond none is chosen, and Save waits for one.
    await openCell('Cell002'); await page.click('#mark');
    assert.deepEqual(await page.evaluate(() => [document.getElementById('mark-reason').value, document.getElementById('mark-save').disabled]), ['', true], 'no reason chosen for a cell that did not respond, flag or no flag');
    await page.selectOption('#mark-reason', 'Bad outline'); assert.ok(await page.isEnabled('#mark-save')); await page.click('#mark-cancel');
    await openCell('Cell005'); await page.click('#mark'); assert.equal(await page.inputValue('#mark-reason'), 'Movement', 'the field moved: Movement'); await page.click('#mark-cancel');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, c = f.S.results.checks.cells.Cell005, kept = c.flags;
      try { const by = (rule, code) => { c.flags = [{ text: rule, rule, code }]; return f.startingReason('Cell005', 'Responded'); };
        return [by('short separate crossings', 'SHORT_SEPARATE_CROSSINGS'), by('on the edge', 'RISE_ON_EDGE'), by('passes through', 'POSSIBLE_PASSING_OBJECT'), by('field moved', 'FIELD_MOVED'), by('crosses the cell', 'POSSIBLE_PASSING_OBJECT'),
          f.startingReason('Cell005', 'Did not respond'), f.startingReason('Cell004', 'Responded')]; } finally { c.flags = kept; } }),
      ['Noise', 'Passing object', 'Passing object', 'Movement', 'Passing object', '', ''], 'the reason by the flag\'s code, whatever rule found it; none on a cell that did not respond or has no flag');
    // The mark button is on a line of its own, under the call line and the flags, where the form opens; the call line has the
    // card's whole width, so it is the same with the button and with the form.
    await openCell('Cell001');
    assert.deepEqual(await page.evaluate(() => { const mark = document.getElementById('mark'), card = document.getElementById('cell-card'), verdict = document.getElementById('verdict'), look = getComputedStyle(card);
      const id = (e) => e ? e.id : null;
      return [id(mark.previousElementSibling), id(mark.nextElementSibling), Math.round(verdict.getBoundingClientRect().width) === Math.round(card.clientWidth - parseFloat(look.paddingLeft) - parseFloat(look.paddingRight)),
        mark.getBoundingClientRect().top >= document.getElementById('cell-checks').getBoundingClientRect().bottom - 0.5]; }), ['cell-checks', 'mark-form', true, true], 'the mark button under the call line and the flags');
    const callLine = () => page.evaluate(() => { const box = document.getElementById('verdict').getBoundingClientRect(); return [box.width, box.height]; }), lineBefore = await callLine();
    // The form is one row, as high as the button's own: opening it moves neither the plot's buttons nor the plot. Also with the
    // card as wide as in the pictures of real recordings (506 px inside); a narrower card lets the form wrap.
    const plotButtonsAt = () => page.evaluate(() => [document.querySelector('[data-plot=mean]').getBoundingClientRect().top, document.getElementById('plot').getBoundingClientRect().top]);
    const formRow = (inner) => page.evaluate(async (inner) => { const card = document.getElementById('cell-card'), form = document.getElementById('mark-form'), look = getComputedStyle(card);
      if (inner) card.style.maxWidth = (inner + parseFloat(look.paddingLeft) + parseFloat(look.paddingRight)) + 'px';
      try { await new Promise((done) => requestAnimationFrame(() => requestAnimationFrame(done)));
        const tops = [...form.children].map((e) => Math.round(e.getBoundingClientRect().top)), control = parseFloat(getComputedStyle(document.documentElement).getPropertyValue('--control-h'));
        // each box inside its label, one after another without overlapping, the last inside the form
        const boxes = [...form.children].map((e) => e.getBoundingClientRect()), inLabels = [...form.querySelectorAll('label')].every((l) => { const box = l.getBoundingClientRect(), c = l.querySelector('select, input').getBoundingClientRect(); return c.left >= box.left - 0.5 && c.right <= box.right + 0.5; });
        const apart = boxes.every((b, i) => i === 0 || boxes[i - 1].right <= b.left + 0.5) && boxes[boxes.length - 1].right <= form.getBoundingClientRect().right + 0.5;
        return { inner: Math.round(card.clientWidth - parseFloat(look.paddingLeft) - parseFloat(look.paddingRight)), oneRow: tops.every((t) => Math.abs(t - tops[0]) <= 1), fits: form.getBoundingClientRect().height <= control + parseFloat(getComputedStyle(form).paddingTop) + 0.5 && inLabels && apart };
      } finally { card.style.maxWidth = ''; } }, inner);
    const buttonsBefore = await plotButtonsAt();
    await page.click('#mark'); assert.deepEqual(await callLine(), lineBefore, 'opening the form leaves the call line as it is');
    { const here = await formRow(null), pictured = await formRow(506);
      assert.deepEqual([here.oneRow, here.fits, pictured.inner, pictured.oneRow, pictured.fits], [true, true, 506, true, true], 'the form on one row, here (' + here.inner + ' px inside the card) and at the pictures\' width');
      assert.deepEqual(await plotButtonsAt(), buttonsBefore, 'opening the form moves neither the plot buttons nor the plot'); }
    assert.deepEqual(await page.evaluate(() => { const r = document.getElementById('mark-reason'), n = document.getElementById('mark-note');
      return { reason: r.value, options: [...r.options].filter((o) => !o.hidden).map((o) => o.textContent), placeholder: n.placeholder, max: n.maxLength, focus: document.activeElement.id, button: document.getElementById('mark').hidden,
        labels: [...document.querySelectorAll('#mark-form label')].map((l) => l.firstChild.textContent.trim()), buttons: [...document.querySelectorAll('#mark-form button')].map((b) => b.textContent) }; }),
      { reason: 'Passing object', options: ['Noise', 'Movement', 'Passing object', 'Bad outline', 'Other'], placeholder: 'Optional', max: 120, focus: 'mark-reason', button: true, labels: ['Reason', 'Note'], buttons: ['Save', 'Cancel'] }, 'the form, its reason taken from the flag');
    assert.deepEqual(await oddShapes(), [], 'the form\'s boxes and buttons have the page\'s height and corner');
    // Tab goes through the form; Esc on any part of it cancels and keeps the cell; another cell closes it.
    const focused = () => page.evaluate(() => document.activeElement.id);
    await page.keyboard.press('Tab'); assert.equal(await focused(), 'mark-note'); await page.keyboard.press('Tab'); assert.equal(await focused(), 'mark-save'); await page.keyboard.press('Tab'); assert.equal(await focused(), 'mark-cancel');
    await page.keyboard.press('Escape');
    assert.deepEqual(await page.evaluate(() => [window.calciumStep4.S.cell, document.getElementById('mark-form').hidden, document.activeElement.id]), ['Cell001', true, 'mark'], 'Esc closes the form and keeps the cell');
    await page.click('#mark'); await page.keyboard.press('Escape'); assert.deepEqual(await page.evaluate(() => [window.calciumStep4.S.cell, document.getElementById('mark-form').hidden]), ['Cell001', true], 'also from the reason');
    await page.click('#mark'); await openCell('Cell005'); assert.ok(await page.isHidden('#mark-form'), 'another cell closes the form');
    // What the page can tell is wrong it says in Fiji's words, and sends nothing; Fiji says the same.
    let reviewPosts = 0; const countPosts = (request) => { if (request.url().includes('/api/review')) reviewPosts++; }; page.on('request', countPosts);
    await openCell('Cell001'); await page.click('#mark'); await page.selectOption('#mark-reason', 'Other');
    assert.equal(await page.getAttribute('#mark-note', 'placeholder'), '', 'Other needs a note: its box does not say Optional');
    const saidFor = async (note) => { await page.fill('#mark-note', note); await page.press('#mark-note', 'Enter'); return text('#mark-said'); };
    assert.deepEqual([await saidFor('   '), await saidFor('=SUM(A1)'), await saidFor('a\tb'), await saidFor('\u2066b')], ['Other needs a short note.', 'The note cannot start with =, +, - or @.', 'The note must be one line.', 'The note must be plain text, without hidden characters.']);
    assert.equal(reviewPosts, 0, 'nothing is sent while the page can tell it would be refused'); page.off('request', countPosts);
    // A refusal that comes or goes changes the height of the cell card and so of the plot: the plot is drawn again at its new size, not stretched.
    assert.deepEqual(await page.evaluate(async () => { const frame = () => new Promise((done) => requestAnimationFrame(() => requestAnimationFrame(done)));
      const plot = document.getElementById('plot'), note = document.getElementById('mark-note'), dpr = window.devicePixelRatio || 1;
      note.value = 'x'; note.dispatchEvent(new Event('input', { bubbles: true })); await frame();   // the refusal goes
      const without = plot.clientHeight; note.value = '   '; document.getElementById('mark-form').requestSubmit(); await frame();   // and comes back
      return { moved: plot.clientHeight !== without, fits: Math.round(plot.clientHeight * dpr) === plot.height }; }), { moved: true, fits: true }, 'the plot follows the height of the cell card');
    const folder = 'browser-results';
    assert.deepEqual([await post({ folder, cell: 'Cell001', call: 'Did not respond', reason: 'Other', note: '  ' }), await post({ folder, cell: 'Cell001', call: 'Did not respond', reason: 'Noise', note: '=SUM(A1)' }),
      await post({ folder, cell: 'Cell001', call: 'Did not respond', reason: 'Noise', note: 'a\tb' }), await post({ folder, cell: 'Cell001', call: 'Did not respond', reason: 'Noise', note: 'x'.repeat(121) })],
      [[400, 'Other needs a short note.'], [400, 'The note cannot start with =, +, - or @.'], [400, 'The note must be one line.'], [400, 'Keep the note to 120 characters or fewer.']]);
    assert.deepEqual([await post({ folder, cell: 'Cell003', call: 'Responded', reason: 'Noise' }), await post({ folder: 'an-older-analysis', cell: 'Cell001', call: 'Did not respond', reason: 'Noise' }),
      await post({ folder, cell: 'Cell002', call: 'Did not respond', reason: 'Noise' }), await post({ folder, cell: 'Cell002', remove: true })],
      [[400, 'Cell003 has no call to review: it is excluded.'], [400, 'This mark is for another analysis. Reload the page.'], [400, 'Cell002\'s automatic call is already Did not respond. The reviewed call must be Responded.'], [400, 'Cell002 has no review mark.']]);
    // A refusal from Fiji is shown as given, in the form, which stays open.
    await page.evaluate(() => { window.calciumStep4.S.results.review.folder = 'an-older-analysis'; });
    await page.selectOption('#mark-reason', 'Noise'); assert.equal(await page.getAttribute('#mark-note', 'placeholder'), 'Optional', 'any other reason: the note is optional');
    await page.fill('#mark-note', ''); await page.click('#mark-save');
    await page.waitForFunction(() => document.getElementById('mark-said').textContent === 'This mark is for another analysis. Reload the page.');
    assert.ok(await page.isVisible('#mark-form'), 'a refusal leaves the form open');
    await page.evaluate(() => { window.calciumStep4.S.results.review.folder = 'browser-results'; });

    // Save: the call after review leads everything that counts, groups, colours or lists by call, at once, in place.
    const callsSha = sha(path.join(out, 'browser-results', 'responder_calls.csv')), reviewFile = path.join(out, 'Fixture recording_reviewed_calls_browser-results.csv');
    let resultsLoads = 0; const countLoads = (request) => { if (/\/api\/results(\?|$)/.test(request.url())) resultsLoads++; }; page.on('request', countLoads);
    // where the Cell box and the arrows are: the Show list keeps its width when its entries say "after review"
    const cellBoxAt = () => page.evaluate(() => ['cell', 'previous', 'next'].map((id) => document.getElementById(id).getBoundingClientRect().left)), boxBefore = await cellBoxAt();
    // Held down, Enter saves once: its repeats neither save again nor press Remove mark, which has the focus once the mark is saved.
    let reviewSent = 0; const countReview = (request) => { if (request.url().includes('/api/review')) reviewSent++; }; page.on('request', countReview);
    await page.selectOption('#mark-reason', 'Passing object'); await page.fill('#mark-note', 'bright spot crosses the cell');
    await page.focus('#mark-note'); await page.keyboard.down('Enter');
    await page.waitForFunction(() => !!window.calciumStep4.S.results.review.marks.Cell001 && document.activeElement.id === 'mark');
    await page.keyboard.down('Enter'); await page.keyboard.down('Enter'); await page.keyboard.up('Enter'); await page.waitForTimeout(300);
    assert.deepEqual([reviewSent, await page.evaluate(() => !!window.calciumStep4.S.results.review.marks.Cell001)], [1, true], 'a held Enter saves the mark once and keeps it'); page.off('request', countReview);
    assert.deepEqual(await cellBoxAt(), boxBefore, 'the Cell box and the arrows stay where they were');
    await page.waitForTimeout(1500);   // two polls
    assert.equal(resultsLoads, 0, 'a mark is applied in place: the results are not loaded again');
    assert.deepEqual(await page.evaluate(() => { const S = window.calciumStep4.S; return [S.cell, S.version === S.state.version, document.activeElement.id, document.getElementById('mark').textContent, document.getElementById('mark-form').hidden]; }),
      ['Cell001', true, 'mark', 'Remove mark', true], 'the cell stays selected; its button now removes the mark');
    // A poll sent before the answer can come back after it, with the older version: it must not read the results again.
    let stalePolls = 0;
    await page.route('**/api/state', async (route) => { const response = await route.fetch(), body = await response.json(); body.version -= 1; stalePolls++; await route.fulfill({ response, body: JSON.stringify(body) }); }, { times: 1 });
    await page.waitForFunction(() => window.calciumStep4.S.version > 0); await page.waitForTimeout(1500); await page.unroute('**/api/state');
    assert.deepEqual([stalePolls, resultsLoads], [1, 0], 'a poll answered with the version before the mark does not load the results again');
    assert.match(await text('#verdict'), /^Cell001: Did not respond after review · passing object: bright spot crosses the cell · automatic call Responded · rise held from 4\.5 s · largest rise 100\.0% at 6\.0 s \(needed 10\.0%, the minimum rise\)$/, 'the call after review first, then the reason and the automatic call, then what the automatic call rests on');
    assert.equal(await page.locator('#verdict .dot').evaluate((d) => getComputedStyle(d).backgroundColor), 'rgb(200, 102, 14)', 'its dot has the colour of the call after review');
    assert.match(await page.locator('#verdict span[title^="Marked "]').getAttribute('title'), /^Marked \d{4}-\d\d-\d\d \d\d:\d\d$/, 'when it was marked, in the tooltip');
    assert.equal(await page.locator('#tile-responding > b').textContent(), '33.3%');
    assert.equal(await page.locator('#tile-responding > span').first().textContent(), 'Capsaicin: 1 of 3 evaluable cells after review; 2 before review', 'the headline after review, with the count before review');
    assert.equal(await page.locator('#tile-responding .tile-line').last().textContent(), 'Quiet enough to detect a 30 % rise: 33.3% (1 of 3; 2 before review)');
    const reviewedAll = (await page.evaluate(() => window.calciumStep4.S.results.tally.rows.find((row) => row.key === 'all'))).reviewed;
    await page.check('#below-zero');
    assert.equal(await page.locator('#tile-responding > span').first().textContent(), 'Capsaicin: ' + reviewedAll.without[0] + ' of ' + reviewedAll.without[1] + ' evaluable cells after review; ' + allRow.without[0] + ' before review · 1 below zero left out');
    // After review on the Responder counts tab follows the switch as the headline does: the headline's figure is in the table
    const afterReviewCell = async () => (await tableNow('table-tally')).rows.find((row) => row[0] === 'All cells')[6];
    const headlineFigure = async () => (await page.locator('#tile-responding > b').textContent()) + ' (' + (await page.locator('#tile-responding > span').first().textContent()).match(/: (\d+ of \d+) evaluable/)[1] + ')';
    assert.deepEqual([await afterReviewCell(), await headlineFigure()], ['0.0% (0 of 2)', '0.0% (0 of 2)'], 'with the cells below zero left out, After review is the headline\'s figure');
    await page.uncheck('#below-zero');
    assert.deepEqual([await afterReviewCell(), await headlineFigure()], ['33.3% (1 of 3)', '33.3% (1 of 3)'], 'and with every cell counted, it is again');
    assert.deepEqual(await page.locator('#group option').allTextContents(), ['All cells (6)', 'Responded after review (1)', 'Did not respond after review (2)', 'Excluded (2)', 'Unscorable (1)', 'Needs quality review (3)', 'Marked (1)']);
    await page.selectOption('#group', '6'); assert.deepEqual(await cellsOffered(), ['Cell001']); await page.selectOption('#group', '2'); assert.deepEqual(await cellsOffered(), ['Cell001', 'Cell002']); await page.selectOption('#group', '0');
    assert.deepEqual(await page.evaluate(() => [0, 1, 4].map((i) => window.calciumStep4.outlineLook(i).pen)), ['not', 'not', 'responded'], 'the outlines take the call after review');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4; return f.heatmapData(f.STYLE.heatmap).calls.map((c) => c.text + ': ' + c.n); }), ['responded: 1', 'did not respond: 2', 'excluded or unscorable: 3'], 'and the heatmap\'s strip and key');
    await page.click('[data-plot=calls]'); assert.deepEqual((await plotNow()).legend.slice(0, 2), ['responded (n = 1)', 'did not respond ± SEM (n = 2)'], 'and the time course\'s groups and legend, in their own words');
    await page.click('[data-plot=cell]');
    assert.deepEqual(await flaggedRows().then((rows) => rows[0]), ['Cell001', 'Did not respond', 'possible passing object, frames 9 to 12']);
    assert.deepEqual(await page.locator('#table-flagged thead th').allTextContents(), ['Cell', 'Call after review', 'Flagged for'], 'the call in the list is named as in the metric tables');
    assert.equal(await page.locator('#table-flagged tbody tr').first().locator('td').nth(1).getAttribute('title'), 'Automatic call: Responded');
    assert.equal((await page.locator('#quality-list .flag-filter').allTextContents())[0], 'Possible passing object: 1 cell, 0 responders');
    const saidAfter = await page.evaluate(() => window.calciumStep4.S.said);
    assert.match(saidAfter, /after review/); assert.equal(await text('#status'), saidAfter, 'what Fiji said after the mark is the status line');
    assert.ok(await page.locator('#status').evaluate((e) => e.scrollHeight <= e.clientHeight + 1), 'and it is read in full');
    assert.deepEqual(await oddShapes(), [], 'Remove mark has the page\'s height and corner');
    await page.screenshot({ path: path.join(out, 'page-review-marked.png') });
    // What a saved file's record says, from what the page sends with it: "after review" only for what follows the marks.
    const sentRecord = async (menu, button, note) => { const sent = page.waitForRequest((r) => r.url().includes('/api/export') && r.method() === 'POST'); await page.locator(menu).click(); await page.click(button);
      const record = JSON.parse((await sent).postData()).metadata; await page.waitForFunction((id) => /^Saved /.test(document.getElementById(id).textContent), note); return record; };
    // A refused removal is said under its cell, and only there: not on the control window, and no longer once another cell is chosen.
    const refusedHere = 'This mark is for another analysis. Reload the page.';
    await page.evaluate(() => { window.calciumStep4.S.results.review.folder = 'an-older-analysis'; });
    await page.click('#mark'); await page.waitForFunction((said) => document.getElementById('mark-said').textContent === said, refusedHere);
    await page.evaluate(() => { window.calciumStep4.S.results.review.folder = 'browser-results'; });
    await page.selectOption('#response-window', 'Control::Ionomycin'); assert.ok(await page.isHidden('#mark-said'), 'a refusal is not shown on the control window');
    assert.deepEqual(await page.locator('#table-flagged thead th').allTextContents(), ['Cell', 'Call', 'Flagged for'], 'the control window\'s list has its own calls');
    { const record = await sentRecord('#export-plot > summary', '#save-plot-folder', 'plot-saved'); assert.deepEqual([record.displayed_window, 'calls' in record, 'review_file' in record], ['Control::Ionomycin', false, false], 'a figure of the control window does not follow the marks'); }
    await page.selectOption('#response-window', 'Stimulus::Capsaicin'); assert.equal(await text('#mark-said'), refusedHere, 'the refusal is back with its cell');
    { const record = await sentRecord('#export-plot > summary', '#save-plot-folder', 'plot-saved'); assert.deepEqual([record.displayed_window, record.calls], ['Stimulus::Capsaicin', 'after review'], 'one of the primary stimulus does'); }
    await openCell('Cell002'); await openCell('Cell001'); assert.ok(await page.isHidden('#mark-said'), 'another cell chosen, the refusal goes');
    // The tables: the summary and the per-cell tables split by the call after review; the tally with the share after review.
    await page.click('.tab[data-tab=metrics]'); await page.selectOption('#metrics-table', 'summary');
    assert.deepEqual((await tableNow('table-metrics')).rows.find((row) => row[0] === 'Peak'), ['Peak', 'ΔF/F₀', 'whole recording', '0.750 ± 0.250', '4', '1.000', '1', '0.500 ± 0.500', '2'], 'the summary by the call after review');
    assert.equal((await sentRecord('#export-metrics > summary', '#metrics-save', 'metrics-saved')).calls, 'after review', 'the summary saved says so');
    await page.selectOption('#metrics-table', 'metrics'); { const perCell = await tableNow('table-metrics'); assert.equal(perCell.head[1], 'Call after review'); assert.equal(at(perCell, 'Cell001', 'Call after review').text, 'Did not respond'); }
    assert.ok(!('calls' in await sentRecord('#export-metrics > summary', '#metrics-save', 'metrics-saved')), 'the saved per-cell table, as saved, does not');
    await page.selectOption('#metrics-table', 'summary');
    await page.click('.tab[data-tab=counts]');
    assert.deepEqual(await page.evaluate(() => { const line = document.getElementById('counts-before'); return [line.hidden, line.textContent, line.nextElementSibling.contains(document.getElementById('table-counts'))]; }), [false, 'Before review', true], 'the engine\'s counts, as saved, are said to be before review');
    assert.ok(!('calls' in await sentRecord('#export-counts > summary', '#counts-save', 'counts-saved')), 'and their record does not say after review');
    { const tally = await tableNow('table-tally');
      assert.deepEqual(tally.head, ['Cells', '% responding', 'Responded', 'Evaluable', 'Cells below zero', 'Without them', 'After review', 'Calls changed']);
      assert.deepEqual(tally.rows[0].slice(0, 4).concat(tally.rows[0].slice(6)), ['All cells', '66.7%', '2', '3', '33.3% (1 of 3)', '1'], 'before review as saved, after review beside it'); }
    await page.click('.tab[data-tab=calls]'); { const table = await tableNow('table-calls');
      assert.deepEqual(table.head.slice(-2), ['Reviewed call', 'Review reason']); assert.equal(table.head.length, 14);
      assert.deepEqual([at(table, 'Cell001', 'Call', 'Stimulus').text, at(table, 'Cell001', 'Reviewed call', 'Stimulus').text, at(table, 'Cell001', 'Review reason', 'Stimulus').text], ['Responded', 'Did not respond', 'Passing object: bright spot crosses the cell'], 'the automatic call stays; the reviewed call and the reason beside it');
      assert.deepEqual([at(table, 'Cell001', 'Reviewed call', 'Control').text, at(table, 'Cell002', 'Reviewed call', 'Stimulus').text, at(table, 'Cell002', 'Review reason', 'Stimulus').text], ['', '', ''], 'empty on the control and on cells without a mark'); }
    assert.equal(await text('#review-count'), '1 cell marked');
    const methods = '1 of 3 cells was reclassified after manual review: 1 from responder to non-responder (passing object).';
    assert.equal(await page.getAttribute('#copy-methods', 'title'), methods); assert.equal(await text('#copy-methods'), 'Copy Methods sentence');
    await page.evaluate(() => { window.copied = []; navigator.clipboard.writeText = async (t) => { window.copied.push(t); }; });
    await page.click('#copy-methods'); assert.deepEqual(await page.evaluate(() => window.copied), [methods]); assert.equal(await text('#calls-saved'), 'Copied');
    await page.locator('#export-calls > summary').click(); await page.click('#calls-save');
    await page.waitForFunction(() => /^Saved in plots: /.test(document.getElementById('calls-saved').textContent));
    { const saved = (await text('#calls-saved')).match(/^Saved in plots: (\S+\.csv)/)[1], lines = fs.readFileSync(path.join(out, 'plots', saved), 'utf8').trim().split('\n'), header = lines[0].split(',');
      assert.match(saved, /^Fixture_recording_responder_calls_\d{8}-\d{6}\.csv$/, 'the file\'s name as before');
      assert.deepEqual(header.slice(-3), ['Reviewed_call', 'Review_reason', 'Review_note']); assert.equal(header.slice(0, -3).join(','), await page.evaluate(() => window.calciumStep4.S.results.tables.calls.header.join(',')), 'the saved columns as saved, the review on the right');
      assert.ok(lines.find((l) => l.startsWith('Fixture,Cell001,Stimulus,')).endsWith(',Did not respond,Passing object,bright spot crosses the cell')); assert.ok(lines.find((l) => l.startsWith('Fixture,Cell001,Control,')).endsWith(',,,'));
      const record = JSON.parse(fs.readFileSync(path.join(out, 'plots', saved.replace(/\.csv$/, '.metadata.json')), 'utf8'));
      assert.deepEqual([record.calls, record.review_file, record.review_sha256], ['after review', 'Fixture recording_reviewed_calls_browser-results.csv', sha(reviewFile)], 'the record says the calls are after review, and which review'); }
    { const lines = fs.readFileSync(reviewFile, 'utf8').trim().split(/\r?\n/); assert.equal(lines.length, 2); assert.match(lines[1], /^Cell001,Responded,Did not respond,Passing object,bright spot crosses the cell,/); }
    assert.equal(sha(path.join(out, 'browser-results', 'responder_calls.csv')), callsSha, 'the engine\'s table is never written to');
    await page.click('#copy-methods'); assert.equal(await text('#calls-saved'), 'Copied');
    // A second mark, the other way: Other needs its note.
    await page.click('.tab[data-tab=review]'); await openCell('Cell002'); await page.click('#mark'); assert.ok(await page.isDisabled('#mark-save'));
    // Enter on the reason saves, by the page's own key handler: a browser that does not submit a form from a list must save too.
    await page.selectOption('#mark-reason', 'Other');
    await page.locator('#mark-reason').evaluate((s) => s.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true, cancelable: true })));
    assert.equal(await text('#mark-said'), 'Other needs a short note.', 'Enter on the reason saves, so the missing note is said');
    await page.selectOption('#mark-reason', 'Noise'); await page.selectOption('#mark-reason', 'Other');
    await page.press('#mark-reason', 'Enter'); assert.equal(await text('#mark-said'), 'Other needs a short note.');
    await page.fill('#mark-note', 'responds after the window'); await page.press('#mark-note', 'Enter');
    await page.waitForFunction(() => !!window.calciumStep4.S.results.review.marks.Cell002);
    assert.equal(await text('#calls-saved'), '', '"Copied" goes once the sentence it copied has changed');
    assert.match(await text('#verdict'), /^Cell002: Responded after review · other: responds after the window · automatic call Did not respond · no sustained threshold crossing · /);
    assert.deepEqual([await page.locator('#tile-responding > b').textContent(), await page.locator('#tile-responding > span').first().textContent(), await text('#review-count'), await page.getAttribute('#copy-methods', 'title')],
      ['66.7%', 'Capsaicin: 2 of 3 evaluable cells after review; 2 before review', '2 cells marked', '2 of 3 cells were reclassified after manual review: 1 from responder to non-responder (passing object) and 1 from non-responder to responder (other).']);
    // Remove mark: once the last is gone, everything is as it was before any mark; the Show list leaves Marked and keeps the cell.
    await page.click('#mark'); await page.waitForFunction(() => !window.calciumStep4.S.results.review.marks.Cell002);
    // A cell marked out of the group the Show list is on stays selected, and Show stays on its group.
    await page.selectOption('#group', '1'); await openCell('Cell005');
    // Held down on the button, Enter opens the form, and its repeats do not then save the reason the flag chose.
    await page.focus('#mark'); await page.keyboard.down('Enter'); await page.waitForFunction(() => document.activeElement.id === 'mark-reason');
    await page.keyboard.down('Enter'); await page.keyboard.down('Enter'); await page.keyboard.up('Enter'); await page.waitForTimeout(300);
    assert.deepEqual(await page.evaluate(() => [document.getElementById('mark-form').hidden, document.getElementById('mark-reason').value, !!window.calciumStep4.S.results.review.marks.Cell005]), [false, 'Movement', false], 'a held Enter opens the form and saves nothing');
    await page.selectOption('#mark-reason', 'Noise'); await page.click('#mark-save');
    await page.waitForFunction(() => !!window.calciumStep4.S.results.review.marks.Cell005);
    assert.deepEqual([await page.evaluate(() => window.calciumStep4.S.cell), await page.inputValue('#group'), await page.locator('#group option:checked').textContent(), await page.inputValue('#cell'), await cellsOffered()],
      ['Cell005', '1', 'Responded after review (0)', 'Cell005', ['Cell005']], 'the cell stays selected after Save, and Show on its group');
    await page.click('#mark'); await page.waitForFunction(() => !window.calciumStep4.S.results.review.marks.Cell005);
    await page.selectOption('#group', '6'); await openCell('Cell001'); await page.click('#mark');
    await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks).length === 0);
    assert.deepEqual(await page.locator('#group option').allTextContents(), ['All cells (6)', 'Responded (2)', 'Did not respond (1)', 'Excluded (2)', 'Unscorable (1)', 'Needs quality review (3)']);
    assert.deepEqual([await page.inputValue('#group'), await page.evaluate(() => window.calciumStep4.S.cell), await text('#mark')], ['0', 'Cell001', 'Mark as not a responder']);
    assert.deepEqual([await page.locator('#tile-responding > b').textContent(), await page.locator('#tile-responding > span').first().textContent()], ['66.7%', 'Capsaicin: 2 of 3 evaluable cells'], 'the headline as before any mark');
    assert.deepEqual([await text('#status'), await page.locator('#status').evaluate((e) => e.classList.contains('full'))], ['Done.', false], 'the status line as before any mark');
    assert.match(await text('#verdict'), /^Cell001: Responded · rise held/);
    assert.deepEqual([await cellBoxAt(), await page.locator('#table-flagged thead th').allTextContents()], [boxBefore, ['Cell', 'Call', 'Flagged for']], 'the Show list, the Cell box and the flagged list as before any mark');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, m = f.exportMetadata(true), t = document.getElementById('table-tally'), c = document.getElementById('table-calls');
      return [['calls', 'review_file', 'review_sha256'].some((k) => k in m), t.tHead.rows[0].cells.length, c.tHead.rows[0].cells.length, JSON.stringify(f.callsRows()) === JSON.stringify((({ data }) => [data.header].concat(data.rows))(f.withCodes('calls', f.S.results.tables.calls))),
        document.getElementById('copy-methods').hidden, document.getElementById('review-count').hidden, document.getElementById('counts-before').hidden]; }), [false, 6, 12, true, true, true, true], 'the record, the tally, the calls table, its export, the Methods line and the counts as before any mark');
    assert.ok(!fs.existsSync(reviewFile), 'no marks, no review file');
    assert.equal(resultsLoads, 0, 'no mark or removal loaded the results again'); page.off('request', countLoads);
    // A pass through the responders: under Show = Responded a mark keeps Show and the cell, which stays in the list in its place
    // until another cell is chosen; the group's count, mean, traces and picture leave it out, and the next arrow goes on to the next responder.
    await page.selectOption('#group', '1'); await openCell('Cell001'); await page.click('#mark'); await page.selectOption('#mark-reason', 'Noise'); await page.click('#mark-save');
    await page.waitForFunction(() => !!window.calciumStep4.S.results.review.marks.Cell001);
    assert.deepEqual([await page.inputValue('#group'), await page.locator('#group option:checked').textContent(), await page.evaluate(() => window.calciumStep4.S.cell), await cellsOffered()],
      ['1', 'Responded after review (1)', 'Cell001', ['Cell001', 'Cell005']], 'Show stays on Responded and the marked cell stays selected, listed in its place');
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4; return [f.plotCells().shown, f.outlineLook(0), !!f.outlineLook(4)]; }), [['Cell005'], null, true], 'the group\'s mean, traces and outlines are of the group alone');
    // A double-click on empty space of the picture shows the whole picture and keeps the cell open, still listed: Show stays on Responded.
    { await page.evaluate(() => window.calciumStep4.showTab('review'));
      const empty = await page.evaluate(() => { const r = document.getElementById('picture').getBoundingClientRect();
        for (let y = r.bottom - 20; y > r.top + 40; y -= 7) for (let x = r.left + 20; x < r.right - 20; x += 11) if (window.calciumStep4.cellAt({ clientX: x, clientY: y }) < 0) return [x, y]; return null; });
      assert.ok(empty, 'an empty point on the picture');
      await page.mouse.dblclick(empty[0], empty[1]);
      assert.deepEqual([await page.inputValue('#group'), await page.evaluate(() => window.calciumStep4.S.cell), await cellsOffered()], ['1', 'Cell001', ['Cell001', 'Cell005']], 'a double-click on empty space keeps the marked cell listed and Show on its group'); }
    await page.click('#next');
    assert.deepEqual([await page.evaluate(() => window.calciumStep4.S.cell), await page.inputValue('#group'), await cellsOffered()], ['Cell005', '1', ['Cell005']], 'the next arrow goes on to the next responder, and the marked cell leaves the list');
    await page.selectOption('#group', '6'); await openCell('Cell001'); await page.click('#mark');
    await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks).length === 0);
    assert.deepEqual([await page.inputValue('#group'), await page.evaluate(() => window.calciumStep4.S.cell)], ['0', 'Cell001'], 'Marked gives way to All cells once no mark is left');
    assert.ok(!fs.existsSync(reviewFile), 'no marks, no review file');

    // A later analysis of the recording: the marks of this one are carried over where the automatic call is as it was, and dropped
    // where it changed, as Fiji says in the status line until the next run. Its checks say that the field's movement was not saved.
    assert.deepEqual([await post({ folder, cell: 'Cell001', call: 'Did not respond', reason: 'Noise', note: '' }), await post({ folder, cell: 'Cell002', call: 'Responded', reason: 'Other', note: 'responds after the window' })].map((answer) => answer[0]), [200, 200]);
    // Marks made outside the page reach it with the results, and the Methods line with them.
    await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks || {}).length === 2);
    const methodsLine = () => page.evaluate(() => { const count = document.getElementById('review-count'), copy = document.getElementById('copy-methods');
      return [count.hidden ? null : count.textContent, copy.hidden ? null : copy.title === window.calciumStep4.S.results.review.methods]; });
    assert.deepEqual(await methodsLine(), ['2 cells marked', true], 'the results loaded with marks: the number marked and the sentence to copy');
    // A refusal shown when a new analysis comes is not shown on it.
    await openCell('Cell001'); await page.evaluate(() => { window.calciumStep4.S.results.review.folder = 'an-older-analysis'; });
    await page.click('#mark'); await page.waitForFunction((said) => document.getElementById('mark-said').textContent === said, refusedHere);
    await page.evaluate(() => { window.calciumStep4.S.results.review.folder = 'browser-results'; });
    fs.writeFileSync(path.join(out, 'rerun'), '');
    fs.writeFileSync(path.join(out, 'checks.csv'), checksTable('Cell002,Responded,,,,,,7,400,yes,no,not checked,no rise at the peak,0,0,').replace(/,yes,yes,/g, ',yes,no,'));
    // A poll that comes while Fiji is still finishing the run (its version moved on, not yet its finish message or tone) reads the
    // results; once the run is over the status line is the finish message all the same.
    let finishing = 0; const versionBefore = await page.evaluate(() => window.calciumStep4.S.version);
    await page.route('**/api/state', async (route) => { const response = await route.fetch(), body = await response.json();
      if (!finishing && Number(body.version) > versionBefore && !body.busy) { finishing++; Object.assign(body, { busy: true, tone: 'plain', message: 'Preparing analysis', stage: 'Checking the calls', percent: 100 }); }
      await route.fulfill({ response, body: JSON.stringify(body) }); });
    await page.check('#confirmed'); await page.click('#analyze');
    await page.waitForFunction(() => { const S = window.calciumStep4.S; return !S.state.busy && S.results && S.results.folder === 'browser-results-2'; }, null, { timeout: 15000 });
    await page.unroute('**/api/state'); assert.equal(finishing, 1, 'one poll came while the run was finishing');
    const carried = 'Review marks from browser-results: 1 carried over, 1 dropped. Cell002: call changed.';
    // A carried mark applies: the status line is Fiji's finish message, with the figures after review and the marks carried over
    // and dropped, as after a mark and as in the Fiji window; also when the results are read again, as after a reload.
    const finished = await page.evaluate(() => window.calciumStep4.S.state.message);
    assert.ok(finished.startsWith('Complete. Protocol, inputs and tables saved in browser-results-2.') && /after review/.test(finished) && finished.endsWith(carried), 'Fiji\'s finish message: ' + finished);
    assert.deepEqual([await text('#status'), await page.locator('#status').evaluate((e) => e.classList.contains('full'))], [finished, true], 'the finish message and the marks carried over and dropped, in the status line');
    assert.ok(await page.locator('#status').evaluate((e) => e.offsetParent !== null && e.scrollHeight <= e.clientHeight + 1), 'in full, without pointing at it');
    await page.evaluate(() => window.calciumStep4.loadResults());
    assert.equal(await text('#status'), await page.evaluate(() => window.calciumStep4.S.state.message), 'and when the results are read again, as after a reload');
    // Without a good finish of the analysis shown (a page opened on another analysis, a failed or an empty run) it is "Done." as before.
    const saidWhen = async (change) => {
      await page.route('**/api/state', async (route) => { const response = await route.fetch(), body = await response.json(); await route.fulfill({ response, body: JSON.stringify(Object.assign(body, change)) }); });
      try { await page.waitForFunction((c) => Object.keys(c).every((k) => window.calciumStep4.S.state[k] === c[k]), change); await page.evaluate(() => window.calciumStep4.loadResults());
        return [await page.evaluate(() => window.calciumStep4.S.said), change.tone ? null : await text('#status')];
      } finally { await page.unroute('**/api/state'); await page.waitForFunction(() => window.calciumStep4.S.state.tone === 'good' && window.calciumStep4.S.state.folder === 'browser-results-2'); await page.evaluate(() => window.calciumStep4.loadResults()); } };
    assert.deepEqual([await saidWhen({ folder: 'another-analysis' }), await saidWhen({ tone: 'bad', message: 'The analysis failed.' }), await saidWhen({ tone: 'warn' })],
      [[null, 'Done. ' + carried], [null, null], [null, null]], 'another analysis\'s finish, a failure or an empty run: the status line as before any mark');
    assert.equal(await text('#status'), finished, 'and the finish message again once it is this analysis\'s good finish');
    assert.deepEqual([await page.evaluate(() => window.calciumStep4.S.cell), await page.isHidden('#mark-said'), await text('#mark')], ['Cell001', true, 'Remove mark'], 'the cell stays selected, without the refusal made on the earlier analysis');
    assert.deepEqual(await methodsLine(), ['1 cell marked', true], 'the Methods line of the marks carried over');
    await openCell('Cell001'); assert.match(await text('#verdict'), /^Cell001: Did not respond after review · noise · from an earlier analysis · automatic call Responded · /);
    assert.equal(await page.locator('#verdict span', { hasText: 'from an earlier analysis' }).getAttribute('title'), 'browser-results', 'the analysis it came from, in the tooltip');
    await openCell('Cell002'); assert.match(await text('#verdict'), /^Cell002: Responded · /); assert.equal(await text('#cell-checks'), 'passing object not checked: no rise at the peak', 'an evaluable cell the checks could not look at says why');
    await openCell('Cell003'); assert.ok(await page.isHidden('#cell-checks'), 'a cell with no call to review does not');
    assert.deepEqual((await page.locator('#quality-list li').allTextContents()).filter((t) => t === 'The movement of the field was not saved, so a moving field is not told from an object.').length, 1, 'the note of the checks, once');
    assert.ok(!(await flaggedRows()).some((row) => row[0] === 'Cell002'), 'a cell that was not checked is not a flagged cell');
    // When that note is all there is to say and the rings were measured, the chip calls it a note of the checks: the cells were checked for a passing object.
    const movedNote = 'The movement of the field was not saved, so a moving field is not told from an object.';
    assert.deepEqual(await noteAlone(true), ['Checks note', movedNote, ''], 'the note leads the list and is not said again under it');
    // On the control window: the note, which is about the whole analysis, and no cell's line about the checks.
    await page.selectOption('#response-window', 'Control::Ionomycin');
    assert.ok((await page.locator('#quality-list li').allTextContents()).includes(movedNote), 'the control window has the note too');
    await openCell('Cell002'); assert.ok(await page.isHidden('#cell-checks'), 'but no line about a cell not checked');
    await page.selectOption('#response-window', 'Stimulus::Capsaicin');
    await openCell('Cell001'); await page.click('#mark'); await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks).length === 0);
    assert.ok((await text('#status')).endsWith(carried), 'what was carried over is still said once the carried mark is removed');
    await page.evaluate(() => window.calciumStep4.loadResults());
    assert.equal(await text('#status'), 'Done. ' + carried, 'and when the results are read again, as after a reload, until the next run');
    assert.deepEqual(await methodsLine(), [null, null], 'and the Methods line goes with the last mark');
    // A review file that does not match its analysis is not used: there is no mark button, and the page says why in Fiji's words,
    // as given. Fiji's sentence is the one it refuses a mark with while the file is there; the results it sends carry no marks then.
    const badReview = path.join(out, 'Fixture recording_reviewed_calls_browser-results-2.csv');
    fs.writeFileSync(badReview, 'not a review file\n');
    const [refusedStatus, blocked] = await post({ folder: 'browser-results-2', cell: 'Cell001', call: 'Did not respond', reason: 'Noise' });
    { const lead = 'The review file ' + path.basename(badReview) + ' '; assert.deepEqual([refusedStatus, blocked.slice(0, lead.length)], [400, lead], 'Fiji names the file'); }
    await page.route('**/api/results', async (route) => { const response = await route.fetch(), body = await response.json(); body.review = Object.assign({}, body.review, { marks: {}, problem: blocked }); await route.fulfill({ response, body: JSON.stringify(body) }); });
    await page.evaluate(() => window.calciumStep4.loadResults()); await openCell('Cell001');
    assert.deepEqual([await page.isHidden('#mark'), await text('#mark-said')], [true, blocked], 'no mark button, and what Fiji says about the file');
    await page.unroute('**/api/results'); fs.unlinkSync(badReview); await page.evaluate(() => window.calciumStep4.loadResults()); await openCell('Cell001');
    assert.deepEqual([await page.isVisible('#mark'), await page.isHidden('#mark-said')], [true, true], 'the file moved away, the button is back');

    // ---------------------------------------------------------------- automatic marks
    // The option, under Advanced, off as the page opens. It is a setting of the analysis: ticking it changes the settings, Fiji
    // keeps it with the others, and unticking it takes the change back.
    const option = () => page.evaluate(() => { const box = document.getElementById('automatic-marks'), label = box.closest('label');
      return { advanced: !!box.closest('#protocol-details'), text: label.textContent.trim(), title: label.title, checked: box.checked }; });
    assert.deepEqual(await option(), { advanced: true, text: 'Mark possible passing objects as not responders', checked: false,
      title: 'After each analysis, a responder flagged as a possible passing object or as a rise on the edge is marked Did not respond, with the reason Passing object. A mark can be removed.' }, 'the option under Advanced, off');
    await page.evaluate(() => { document.getElementById('protocol-details').open = true; });
    assert.deepEqual(await page.evaluate(() => { const box = document.getElementById('automatic-marks').getBoundingClientRect(), other = document.getElementById('confirmed').getBoundingClientRect(); return [box.width > 0, box.width === other.width, box.height === other.height]; }), [true, true, true], 'its tick box is as the page\'s other tick boxes');
    const keptValues = () => page.evaluate(() => fetch('api/session').then((r) => r.json()).then((s) => s.draft.values));
    const synced = () => page.waitForFunction(() => { const f = window.calciumStep4; return !!f.S.state && f.S.draftAck === JSON.stringify(f.collectDraft()); });
    assert.ok(!('automatic_marks' in await keptValues()), 'off: no key in the draft Fiji keeps');
    await page.check('#automatic-marks'); await synced(); await page.waitForFunction(() => window.calciumStep4.S.state.settingsPending === true);
    assert.deepEqual([(await keptValues()).automatic_marks, await page.isVisible('#applied-note'), await page.isVisible('#discard')], [true, true, true], 'on: Fiji keeps it, and the settings are not those of the analysis shown');
    await page.uncheck('#automatic-marks'); await synced(); await page.waitForFunction(() => window.calciumStep4.S.state.settingsPending === false);
    assert.deepEqual(['automatic_marks' in await keptValues(), await page.isHidden('#applied-note'), await page.isHidden('#discard')], [false, true, true], 'off again: the settings of the analysis shown');
    await page.check('#automatic-marks'); await synced();
    await page.reload(); await page.waitForFunction(() => window.calciumStep4 && window.calciumStep4.S.results && window.calciumStep4.S.state && !document.getElementById('plot-card').hidden);
    assert.equal(await page.isChecked('#automatic-marks'), true, 'after a reload, the option as Fiji keeps it');
    // A run with the option on: the responder flagged as a possible passing object, Cell001, is marked Did not respond
    // automatically; neither the cell flagged that did not respond (Cell002) nor the one whose field moved (Cell005) is. The
    // earlier analysis has no marks now, so none is carried over.
    for (const name of ['reviewed_calls', 'percent_responding']) { const file = path.join(out, 'Fixture recording_' + name + '_browser-results.csv'); if (fs.existsSync(file)) fs.unlinkSync(file); }
    const ranOn = await page.evaluate(() => window.calciumStep4.S.version);
    await page.check('#confirmed'); await page.click('#analyze');
    await page.waitForFunction((v) => { const S = window.calciumStep4.S; return !S.state.busy && S.version > v && !!S.results && S.results.folder === 'browser-results-2'; }, ranOn, { timeout: 15000 });
    assert.deepEqual(await page.evaluate(() => { const r = window.calciumStep4.S.results; return [r.appliedDraft.values.automatic_marks, Object.entries(r.review.marks).map(([cell, m]) => [cell, m.markedBy, m.reviewed, m.reason, m.note])]; }),
      [true, [['Cell001', 'automatic', 'Did not respond', 'Passing object', '']]], 'one automatic mark: the flagged responder');
    const finishedOn = await page.evaluate(() => window.calciumStep4.S.state.message);
    assert.ok(finishedOn.startsWith('Complete. Protocol, inputs and tables saved in browser-results-2.') && / after review /.test(finishedOn) && finishedOn.endsWith(' Marked automatically: Cell001.'), 'Fiji\'s finish message names the cell marked: ' + finishedOn);
    assert.equal(await text('#status'), finishedOn, 'and it is the status line');
    await openCell('Cell001');
    assert.match(await text('#verdict'), /^Cell001: Did not respond after review · passing object · marked automatically · automatic call Responded · rise held/, 'the call line says the mark was made automatically');
    assert.deepEqual([await text('#mark'), (await page.locator('#group option').allTextContents()).slice(-1)[0], await text('#review-count'), await page.getAttribute('#copy-methods', 'title')],
      ['Remove mark', 'Marked (1)', '1 cell marked', '1 of 3 cells was reclassified automatically from responder to non-responder (1 possible passing object).'], 'Remove mark, the Marked group, and the Methods sentence of an automatic mark');
    assert.deepEqual([await page.locator('#tile-responding > b').textContent(), await page.locator('#tile-responding > span').first().textContent()], ['66.7%', 'Capsaicin: 2 of 3 evaluable cells after review; 3 before review'], 'the headline follows it as any mark');
    await page.click('.tab[data-tab=calls]');
    { const table = await tableNow('table-calls'); assert.deepEqual([at(table, 'Cell001', 'Reviewed call', 'Stimulus').text, at(table, 'Cell001', 'Review reason', 'Stimulus').text], ['Did not respond', 'Passing object, automatic'], 'Review reason says it was made automatically'); }
    await page.locator('#export-calls > summary').click(); await page.click('#calls-save'); await page.waitForFunction(() => /^Saved in plots: /.test(document.getElementById('calls-saved').textContent));
    { const saved = savedCsv(await text('#calls-saved'));
      assert.deepEqual(saved.lines[0].slice(-4), ['Reviewed_call', 'Review_reason', 'Review_note', 'Marked_by'], 'the Export: Marked_by after Review_note');
      assert.deepEqual([saved.at('Marked_by', 'Cell001', 'Stimulus'), saved.at('Review_reason', 'Cell001', 'Stimulus'), saved.at('Marked_by', 'Cell002', 'Stimulus'), saved.at('Marked_by', 'Cell001', 'Control')], ['automatic', 'Passing object', '', ''], 'automatic on the marked row, nothing on the others'); }
    { const lines = fs.readFileSync(path.join(out, 'Fixture recording_reviewed_calls_browser-results-2.csv'), 'utf8').trim().split(/\r?\n/);
      assert.deepEqual([lines.length, lines[0].split(',').pop(), lines[1].split(',').pop()], [2, 'Marked_by', 'automatic'], 'the review file says who marked'); }
    // Remove mark works on it as on any other. What the run marked stays in the finish message until the next run.
    await page.click('.tab[data-tab=review]'); await openCell('Cell001'); await page.click('#mark');
    await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks).length === 0);
    assert.deepEqual([await methodsLine(), (await text('#status')).endsWith(' Marked automatically: Cell001.'), await text('#mark')], [[null, null], true, 'Mark as not a responder'], 'the mark removed, the finish message kept');
    assert.match(await text('#verdict'), /^Cell001: Responded · /);
    await page.evaluate(() => window.calciumStep4.loadResults());
    assert.equal(await text('#status'), await page.evaluate(() => window.calciumStep4.S.state.message), 'also when the results are read again');
    // With the option off a run marks nothing, and nothing of it is shown or saved: no Marked_by, "Done." as before.
    await page.evaluate(() => { document.getElementById('protocol-details').open = true; }); await page.uncheck('#automatic-marks'); await synced();
    const ranOff = await page.evaluate(() => window.calciumStep4.S.version);
    await page.check('#confirmed'); await page.click('#analyze');
    await page.waitForFunction((v) => { const S = window.calciumStep4.S; return !S.state.busy && S.version > v && !!S.results && S.results.folder === 'browser-results-2'; }, ranOff, { timeout: 15000 });
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4, r = f.S.results, coded = f.withCodes('calls', r.tables.calls).data;
      return ['automatic_marks' in r.appliedDraft.values, Object.keys(r.review.marks).length, JSON.stringify(f.callsRows()) === JSON.stringify([coded.header].concat(coded.rows)), / Marked automatically| automatic mark/.test(f.S.state.message)]; }),
      [false, 0, true, false], 'off: no mark, the calls Export without Marked_by, nothing said of automatic marks');
    assert.equal(await text('#status'), 'Done.');
    assert.ok(!fs.existsSync(path.join(out, 'Fixture recording_reviewed_calls_browser-results-2.csv')), 'and no review file');
    // With the option on and no responder to mark, "Done." too: Fiji's finish message is the status line only once a run has
    // marked a cell. A mark set by hand and removed again leaves "Done." as well.
    const runNow = async () => { const from = await page.evaluate(() => window.calciumStep4.S.version); await page.check('#confirmed'); await page.click('#analyze');
      await page.waitForFunction((v) => { const S = window.calciumStep4.S; return !S.state.busy && S.version > v && !!S.results && S.results.folder === 'browser-results-2' && !!S.state.message; }, from, { timeout: 15000 }); };
    fs.writeFileSync(path.join(out, 'checks.csv'), checksTable('Cell002,Did not respond,,,,,,7,400,yes,yes,checked,,0,0,').replace(/^Cell001,.*$/m, 'Cell001,Responded,,,9,12,,8,400,yes,yes,checked,,2,6,'));
    await page.evaluate(() => { document.getElementById('protocol-details').open = true; }); await page.check('#automatic-marks'); await synced(); await runNow();
    assert.deepEqual(await page.evaluate(() => { const r = window.calciumStep4.S.results; return [r.appliedDraft.values.automatic_marks, Object.keys(r.review.marks).length, /^Complete\. /.test(window.calciumStep4.S.state.message)]; }), [true, 0, true], 'the option on, no mark made');
    assert.equal(await text('#status'), 'Done.', 'no automatic mark: "Done.", not the finish message');
    await openCell('Cell001'); await page.click('#mark'); await page.selectOption('#mark-reason', 'Noise'); await page.click('#mark-save');
    await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks).length === 1); assert.ok((await text('#status')).startsWith('Complete. '), 'with a mark, the finish message');
    await page.click('#mark'); await page.waitForFunction(() => Object.keys(window.calciumStep4.S.results.review.marks).length === 0);
    assert.equal(await text('#status'), 'Done.', 'the mark removed: "Done." again');
    await page.evaluate(() => { document.getElementById('protocol-details').open = true; }); await page.uncheck('#automatic-marks'); await synced();
    // An analysis whose call checks table was written by 2026.11.24.1, without QC_flags: its codes come from its Rule column, in
    // the Flags column and in the Export's QC_flags, and Around the cell shows its cells.
    const oldTable = checksTable('Cell002,Did not respond,"rise on the left edge, not in the cell, frames 20 to 24",on the edge,20,24,"' + edgeDetail + '",7,400,yes,yes,checked,,0,0,RISE_ON_EDGE').trim().split('\n').map((line) => line.replace(/,[^,]*$/, '')).join('\n') + '\n';
    // (Fiji's answer for Around the cell is given here: this fixture has no ring files)
    await page.route('**/api/around*', (route) => new URL(route.request().url()).searchParams.get('cell') === 'Cell001' ? route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(Object.assign({}, faked.Cell001, { folder: 'browser-results-2' })) }) : route.continue());
    fs.writeFileSync(path.join(out, 'checks.csv'), oldTable); await runNow();
    assert.equal(fs.readFileSync(path.join(out, 'Fixture recording_call_checks_browser-results-2.csv'), 'utf8').split(/\r?\n/)[0].split(',').length, 15, 'the table as 2026.11.24.1 wrote it');
    assert.deepEqual(await page.evaluate(() => { const c = window.calciumStep4.S.results.checks.cells; return [c.Cell001.codes, c.Cell002.codes, c.Cell005.codes, c.Cell003.codes]; }),
      [['POSSIBLE_PASSING_OBJECT'], ['RISE_ON_EDGE'], ['FIELD_MOVED', 'SHORT_SEPARATE_CROSSINGS'], []], 'the codes from the Rule column');
    await page.click('.tab[data-tab=calls]');
    { const shown = await tableNow('table-calls');
      assert.deepEqual([at(shown, 'Cell001', 'Flags', 'Stimulus'), at(shown, 'Cell005', 'Flags', 'Stimulus').text], [{ text: 'Duration under a frame · Possible passing object', tip: 'DURATION_BELOW_SAMPLE_INTERVAL; POSSIBLE_PASSING_OBJECT' }, 'Duration under a frame · Field moved · Short separate crossings'], 'in the Flags column'); }
    await page.locator('#export-calls > summary').click(); await page.click('#calls-save'); await page.waitForFunction((before) => { const now = document.getElementById('calls-saved').textContent; return /^Saved in plots: /.test(now) && now !== before; }, await text('#calls-saved'));
    { const saved = savedCsv(await text('#calls-saved'));
      assert.deepEqual([saved.at('QC_flags', 'Cell001', 'Stimulus'), saved.at('QC_flags', 'Cell002', 'Stimulus'), saved.record.qc_flags], ['DURATION_BELOW_SAMPLE_INTERVAL; POSSIBLE_PASSING_OBJECT', 'DURATION_BELOW_SAMPLE_INTERVAL; RISE_ON_EDGE', codesNote], 'and in the Export\'s QC_flags'); }
    await page.click('.tab[data-tab=review]'); await openCell('Cell001'); await page.waitForFunction(() => { const a = window.calciumStep4.Around.got.get('Cell001'); return !!a && Array.isArray(a.time); });
    await page.click('[data-plot=around]'); await page.waitForFunction(() => !!window.calciumStep4.Plot.around);
    assert.deepEqual(await page.evaluate(() => [window.calciumStep4.Plot.geom.items.map((i) => i.text), document.getElementById('plot-empty').hidden]), [['Cell001', 'ring, lower right', 'far ring, upper left', 'other ring sectors', 'flagged frames'], true], 'Around the cell shows the cell');
    await page.unroute('**/api/around*'); await page.click('[data-plot=cell]');
    refusals = false;
    assert.deepEqual(problems, []);
    console.log('PASS: movie preview and pointer scrubbing, Play/manual seek, static view, stale/error frames, portable protocols/precision, primary/control results, picture zoom/move/region/map, plain outlines and numbered cells, saving an area and a baseline and peak pair, irregular-time heatmap, pending settings, exports, figure editor numbers, limits, fit and Reset, responsive layouts, the flags of the checks and their codes in the tables, around the cell, review marks and carried marks, automatic marks, a call checks table of 2026.11.24.1');
  } catch (error) {
    if (page) {
      // what the movie was doing, read before the screenshot (a full-page screenshot resizes the page and reloads the frame)
      const movie = await page.evaluate(() => { const app = window.calciumStep4; return !app ? 'page not loaded' : JSON.stringify({ frame: app.S.frame, shown: app.Pic.shownFrame, error: app.Pic.error,
        message: document.getElementById('movie-message').textContent, messageHidden: document.getElementById('movie-message').hidden, visibility: document.visibilityState, note: document.getElementById('picture-note').textContent, buffer: app.Buf && { version: app.Buf.version, failed: Array.from(app.Buf.failed.keys()), inflight: Array.from(app.Buf.inflight.keys()), decoding: Array.from(app.Buf.decoding.keys()) },
        serverVersion: app.S.state && app.S.state.frameVersion }); }).catch((problem) => 'unreadable: ' + problem.message);
      console.error('Movie state at failure: ' + movie);
      await page.screenshot({ path: path.join(out, 'page-failure.png'), fullPage: true }).catch(() => {});
      fs.writeFileSync(path.join(out, 'browser-failure.txt'), String(error.stack || error));
    }
    throw error;
  } finally {
    if (browser) await browser.close();
  }
})().catch((error) => { console.error(error); process.exitCode = 1; });

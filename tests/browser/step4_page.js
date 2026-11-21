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
    let injectedFrameFailure = false;
    page.on('console', (message) => {
      if (message.type() !== 'error') return;
      if (injectedFrameFailure && /503/.test(message.text())) return;
      problems.push('console: ' + message.text());
    });
    page.on('pageerror', (error) => problems.push('page error: ' + error.message));
    const text = (selector) => page.textContent(selector);
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
    assert.equal(await page.inputValue('#view'), 'raw');
    assert.ok(await page.isEnabled('#play') && await page.isEnabled('#time'), 'movie controls work before an analysis');
    assert.ok(await page.isHidden('#result-tiles') && await page.isHidden('#plot-card'), 'preview has no empty result dashboard');
    assert.deepEqual(await page.locator('#cell option').evaluateAll((o) => o.map((x) => x.textContent)), ['None', 'Cell001', 'Cell002', 'Cell003', 'Cell004', 'Cell005', 'Cell006'], 'cells can be selected before analysis, and none is a choice');
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
    assert.deepEqual(await page.locator('details.export').evaluateAll((menus) => menus.map((m) => m.id + ':' + m.querySelectorAll('button').length)), ['export-picture:4', 'export-plot:3', 'export-heatmap:3', 'export-metrics:2', 'export-data:2']);
    assert.equal(await page.locator('main button', { hasText: /^(Save|Download)/ }).evaluateAll((buttons) => buttons.filter((b) => !b.closest('details.export, .save-options')).length), 0, 'every save button is in an Export menu');
    const initial = await page.evaluate(() => window.calciumStep4.collectDraft());
    const originalProtocol = initial.protocol;
    assert.equal(originalProtocol.schema_version, 2);
    assert.equal(originalProtocol.windows.length, 2);
    await screenshots('before');

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

    // Playback: 10× by default, on a clock in recording time, from a buffer that makes scrubbing back free.
    assert.equal(await page.inputValue('#speed'), '10');
    let frameRequests = 0; const countFrames = (request) => { if (request.url().includes('/api/frame?')) frameRequests++; }; page.on('request', countFrames);
    await slider(1); const started = Date.now(); await page.click('#play');
    await page.waitForFunction(() => window.calciumStep4.S.timer === null, null, { timeout: 15000 });
    const took = (Date.now() - started) / 1000;
    assert.equal(await page.inputValue('#time'), '30', 'playback stops at the last frame');
    assert.ok(took > 1.1 && took < 4, '14.5 s of recording at 10× plays in about 1.45 s, not ' + took + ' s');
    await waitFrame(30); frameRequests = 0;
    for (const f of [24, 17, 9, 2]) await slider(f);
    assert.equal(frameRequests, 0, 'scrubbing back over frames already seen makes no requests');
    page.off('request', countFrames);
    await page.selectOption('#speed', '50'); await page.click('#play');
    await page.waitForFunction(() => /playing every 4th frame/.test(document.getElementById('picture-note').textContent));
    await page.click('#play');
    assert.ok(!(await page.evaluate(() => window.calciumStep4.S.timer)));
    // A manual seek stops the timer.
    await slider(3);
    await page.selectOption('#speed', '0.5'); await page.click('#play');
    await page.waitForFunction(() => window.calciumStep4.S.frame >= 5);
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
    assert.deepEqual(submitted.windows.map((w) => w.Stimulus), ['Capsaicin', 'Ionomycin']);
    assert.equal(submitted.windows[0].Floor_local_fraction, 0.1);
    assert.equal(submitted.windows[1].Floor_local_fraction, 0.4);
    assert.match(sent.timestamps_sha256, /^[a-f0-9]{64}$/);
    assert.equal(sent.input_mode, 'background_corrected_fluorescence');
    assert.match(await text('#tile-responding'), /67%.*2 of 3 evaluable/);
    assert.equal((await text('#tile-registered')).replace(/\D/g, ''), '6');
    // All cells lead the headline; the quiet cells follow on a line of their own; every figure comes from the tally worked out in Fiji.
    assert.match(await page.locator('#tile-responding > b').textContent(), /^67%$/);
    assert.match(await page.locator('#tile-responding .tile-line').last().textContent(), /^Quiet enough to detect a 30 % rise: \d+% \(\d+ of \d+\)$/);
    // The cell that went below zero stays counted until the switch leaves it out; then every figure follows, and comes back.
    const counted = await text('#tile-responding'), tally = await page.evaluate(() => window.calciumStep4.S.results.tally);
    assert.deepEqual(tally.belowZero, ['Cell005']); assert.equal(await text('#below-zero-text'), 'Leave out 1 cell below zero'); assert.ok(await page.isVisible('#below-zero-label') && !(await page.isChecked('#below-zero')));
    const allRow = tally.rows.find((row) => row.key === 'all'); assert.deepEqual(allRow.counted.slice(0, 2), [2, 3]);
    await page.check('#below-zero');
    assert.match(await text('#tile-responding'), new RegExp('Capsaicin: ' + allRow.without[0] + ' of ' + allRow.without[1] + ' evaluable cells · 1 below zero left out'));
    assert.equal(await page.locator('#tile-responding > b').textContent(), Math.round(100 * allRow.without[0] / allRow.without[1]) + '%');
    await page.uncheck('#below-zero'); assert.equal(await text('#tile-responding'), counted);
    // The rule line says what is applied: samples in a row, also for a length typed in seconds.
    assert.match(await text('#rule-summary-0'), /for 2 samples in a row$/);
    assert.deepEqual(await page.evaluate(() => { const f = window.calciumStep4.samplesFor, at = (step) => Array.from({ length: 20 }, (_, i) => i * step); return [f(1, at(1)), f(1, at(0.999946296215057)), f(1, at(0.5)), f(1, at(0.19999999999999996)), f(0.4, at(0.5)), f(1, null), f(0, at(1))]; }), [2, 3, 3, 6, 2, null, null]);
    assert.equal(await page.inputValue('#view'), 'raw', 'analysis preserves the selected movie view');
    await page.selectOption('#cell', 'Cell001'); assert.match(await text('#verdict'), /Cell001: Responded/);
    // The duration notice is set on every cell by the engine but belongs to the recording: it is shown once and flags no cell.
    assert.deepEqual(await page.locator('#group option').allTextContents(), ['All cells (6)', 'Responded (2)', 'Did not respond (1)', 'Excluded (2)', 'Unscorable (1)', 'Needs quality review (1)']);
    assert.equal(await text('#quality-summary'), '1 cell flagged'); assert.ok(await page.isHidden('#quality-list'), 'the chip does not push the data down');
    await page.selectOption('#cell', 'Cell005'); assert.match(await text('#cell-quality'), /Low F₀/);
    await page.locator('#quality-review summary').click(); assert.match(await text('#quality-list'), /Applies to the whole recording: Minimum response duration is shorter than one frame interval/);
    assert.equal(await text('#quality-lead'), '1 cell flagged, including 1 responder.');
    assert.ok(await page.evaluate(() => { const tabs = document.querySelector('.tabs').getBoundingClientRect(), row = document.getElementById('result-tiles').getBoundingClientRect(); return tabs.top - row.bottom < 40; }), 'the open list lies over the page');
    await page.locator('#quality-list .flag-filter').first().click();
    assert.equal(await page.inputValue('#group'), '5'); assert.deepEqual(await page.locator('#cell option').evaluateAll((o) => o.map((x) => x.value).filter(Boolean)), ['Cell005']);
    await page.locator('#quality-list .flag-filter').first().click(); assert.equal(await page.evaluate(() => window.calciumStep4.S.qualityFlag), null);
    await page.selectOption('#group', '0');
    await page.locator('#quality-review summary').click();
    await page.selectOption('#cell', 'Cell001');
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
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value).filter(Boolean)), ['Cell002']);
    await page.selectOption('#response-window', 'Stimulus::Capsaicin');
    assert.match(await text('#tile-responding'), /Capsaicin: 2 of 3 evaluable/);
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value).filter(Boolean)), ['Cell001', 'Cell005']);
    await page.selectOption('#group', '3');
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value).filter(Boolean)), ['Cell003', 'Cell004']);
    await page.selectOption('#group', '4');
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value).filter(Boolean)), ['Cell006']);
    await page.selectOption('#group', '0'); await page.selectOption('#cell', 'Cell002');

    // Dragging begins on another cell's trace; pointer capture keeps it scrubbing.
    await page.locator('#plot').scrollIntoViewIfNeeded();
    await page.click('[data-plot=cell]'); await page.check('#every');   // the page opens on the average alone
    const points = await page.evaluate(() => {
      const { S, Plot } = window.calciumStep4, g = Plot.geom, box = document.getElementById('plot').getBoundingClientRect();
      const x = (t) => box.left + g.L + (t - g.x0) / (g.x1 - g.x0) * g.W;
      return { start: { x: x(S.results.time[12]), y: box.top + g.T + (g.hi - S.results.traces.Cell001[12]) / (g.hi - g.lo) * g.H },
        end: { x: x(S.results.time[21]), y: box.top + g.T + g.H / 2 } };
    });
    await page.mouse.move(points.start.x, points.start.y); await page.mouse.down();
    assert.equal(await page.evaluate(() => window.calciumStep4.S.cell), 'Cell001');
    assert.ok(await page.evaluate(() => window.calciumStep4.Plot.scrubbing));
    await page.mouse.move(points.end.x, points.end.y, { steps: 6 }); await page.mouse.up();
    await waitFrame(22);
    assert.equal(await page.evaluate(() => window.calciumStep4.Plot.scrubbing), false);
    assert.match(await text('#time-label'), /15\.50 s.*Frame 22 \/ 30/);

    // No cell selected is a state of its own (the page opens that way), and the plot follows the selection.
    const app = (what) => page.evaluate((key) => { const a = window.calciumStep4; return key === 'cell' ? a.S.cell : key === 'plot' ? a.S.plotMode : a.Pic.zoom; }, what);
    await page.selectOption('#cell', '');
    assert.equal(await app('cell'), null); assert.equal(await text('#verdict'), 'No cell selected.'); assert.ok(await page.isDisabled('#zoom-cell'));
    assert.equal(await app('plot'), 'mean', 'with no cell the plot shows the average');
    await page.selectOption('#cell', 'Cell002');
    assert.equal(await app('plot'), 'cell', 'choosing a cell shows that cell over the average');
    await page.evaluate(() => document.activeElement && document.activeElement.blur()); await page.keyboard.press('Escape');
    assert.equal(await app('cell'), null, 'Escape leaves no cell selected'); assert.equal(await app('plot'), 'mean');
    await page.keyboard.press('ArrowRight'); assert.equal(await app('cell'), 'Cell001', 'the arrows start from the first cell');
    await page.click('[data-plot=calls]'); await page.selectOption('#cell', 'Cell002'); assert.equal(await app('plot'), 'calls', 'a plot type chosen by hand stays');
    await page.click('[data-plot=cell]'); await page.selectOption('#cell', 'Cell001');
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
    const findFolder = (dir) => { for (const entry of fs.readdirSync(dir, { withFileTypes: true })) { if (!entry.isDirectory()) continue; const full = path.join(dir, entry.name); if (entry.name === framesFolder) return full; const below = findFolder(full); if (below) return below; } return null; };
    const framesOnDisk = findFolder(out); assert.ok(framesOnDisk, 'the frames folder is on disk: ' + framesFolder);
    assert.deepEqual(fs.readdirSync(framesOnDisk).filter((name) => /\.png$/.test(name)).sort(), ['all_frames_combined.png', 'frame_01.png', 'frame_11.png', 'frame_21.png']);
    const layoutSaved = JSON.parse(fs.readFileSync(path.join(framesOnDisk, 'about_these_frames.metadata.json'), 'utf8')).combined_picture;
    assert.ok(layoutSaved && layoutSaved.columns * layoutSaved.rows >= 3 && layoutSaved.file === 'all_frames_combined', 'how the combined picture is laid out is recorded');

    // Actual time determines heatmap widths/clicks, and the acquisition gap stays blank.
    await page.click('.tab[data-tab=heatmap]'); await page.selectOption('#order', 'orig'); await settledPaint();
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
    await page.click('.tab[data-tab=calls]'); assert.equal(await page.locator('#table-calls tbody tr').count(), 12);
    const coreColumns = await page.locator('#table-calls thead th').count();
    await page.check('#all-call-columns'); assert.ok(await page.locator('#table-calls thead th').count() > coreColumns);
    await page.locator('#table-calls').evaluate((e) => e.parentElement.scrollLeft = 800);
    assert.ok(await page.locator('#table-calls tbody .cell-column').first().evaluate((e) => {const a=e.getBoundingClientRect(),b=e.closest('.scroller').getBoundingClientRect();return a.left>=b.left-1&&a.left<b.left+5;}));
    await page.uncheck('#all-call-columns');
    await page.screenshot({ path: path.join(out, 'page-calls.png') });
    await page.click('.tab[data-tab=counts]');
    assert.deepEqual(await page.evaluate(() => { const t = document.getElementById('table-tally'); return { head: [...t.tHead.rows[0].cells].slice(0, 4).map((c) => c.textContent), first: [...t.tBodies[0].rows[0].cells].slice(0, 4).map((c) => c.textContent), rows: t.tBodies[0].rows.length }; }),
      { head: ['Cells', '% responding', 'Responded', 'Evaluable'], first: ['All cells', '66.7', '2', '3'], rows: 2 }, 'the same tally on the Responder counts tab');
    await page.click('.tab[data-tab=metrics]');
    // The summary: each metric once, with its unit and where it is measured; saved columns with the same numbers share a row;
    // the timing and shape of a response are given only for cells that responded.
    const summary = await page.evaluate(() => { const t = document.getElementById('table-metrics'), saved = window.calciumStep4.S.results.tables.metrics;
      return { head: [...t.tHead.rows[0].cells].map((c) => c.textContent), rows: [...t.tBodies[0].rows].map((r) => [...r.cells].map((c) => c.textContent)), frozen: t.tHead.rows[0].cells[0].classList.contains('cell-column'), savedHeader: saved.header }; });
    assert.deepEqual(summary.head.slice(0, 3), ['Metric', 'Unit', 'Measured']); assert.equal(summary.head[summary.head.length - 1], 'Saved as'); assert.ok(summary.frozen, 'the metric names stay in view');
    const savedAs = (row) => row[row.length - 1].split(', '), rowOf = (column) => summary.rows.find((row) => savedAs(row).includes(column)), allSaved = summary.rows.flatMap(savedAs);
    assert.equal(new Set(allSaved).size, allSaved.length, 'no saved column is summarised twice'); assert.ok(allSaved.every((name) => summary.savedHeader.includes(name)));
    assert.equal(rowOf('Peak_dFF0'), rowOf('Response_Amplitude'), 'columns with the same numbers share a row'); assert.deepEqual(rowOf('Peak_dFF0').slice(0, 3), ['Peak', 'ΔF/F₀', 'whole recording']);
    assert.notEqual(rowOf('Time_to_Peak'), rowOf('Time_to_50_Peak')); assert.deepEqual(rowOf('Time_to_Peak').slice(0, 3), ['Peak time', 's', 'from recording start']);
    const notResponded = summary.head.indexOf('Did not respond'), responded = summary.head.indexOf('Responded');
    for (const column of ['Time_to_Peak', 'Time_to_50_Peak', 'Rise_Time', 'FWHM']) { assert.equal(rowOf(column)[notResponded], '–', column + ' is not given for cells that did not respond'); assert.match(rowOf(column)[responded], /n = \d/); }
    assert.match(rowOf('Peak_dFF0')[notResponded], /n = \d/, 'the peak is given for every cell'); assert.match(await text('#metric-help'), /cells that responded/);
    await page.selectOption('#metrics-table', 'groups');
    assert.match(await text('#table-metrics'), /No rows were saved/);
    await page.selectOption('#metrics-table', 'metrics'); assert.match(await text('#metric-help'), /signed net integral/);
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
    await page.locator('#export-plot > summary').click(); await page.locator('#export-plot .save-options-toggle').click();
    assert.ok(await page.isVisible('#save-options') && !(await page.locator('#export-plot').evaluate((menu) => menu.open)), 'Options opens from the Export menu');
    assert.ok(await page.evaluate(() => { const panel = document.getElementById('save-options').getBoundingClientRect(), button = document.querySelector('#export-plot > summary').getBoundingClientRect(); return panel.bottom <= button.top + 1 && Math.abs(panel.right - button.right) < 24; }), 'the options sit against the Export button');
    await page.keyboard.press('Escape'); await page.locator('#plot').click({ position: { x: 5, y: 5 } }); assert.ok(await page.isHidden('#save-options'));
    for (const download of downloads) await download.saveAs(path.join(out, download.suggestedFilename()));
    const names = downloads.map((d) => d.suggestedFilename()), image = names.find((n) => n.endsWith('.png'));
    assert.match(image, /^Fixture_recording_traces_Cell\w+_\d{8}-\d{6}\.png$/, 'plot names carry the recording, the cell and the time');
    assert.ok(names.includes(image.replace(/\.png$/, '.metadata.json')));
    const png = fs.readFileSync(path.join(out, image));
    assert.ok(png.readUInt32BE(16) >= 1600, 'exports are at least 1600 pixels wide, whatever the window size');
    assert.equal(await page.evaluate(() => window.calciumStep4.S.exportScale), null, 'the screen is redrawn at its own size afterwards');
    const metadata = JSON.parse(fs.readFileSync(path.join(out, image.replace(/\.png$/, '.metadata.json')), 'utf8'));
    assert.ok(metadata.results_folder && metadata.time_range.every((t) => typeof t === 'number'), 'metadata names the results folder and the shown time range');
    assert.equal(metadata.appliedFingerprint, applied.fingerprint);
    assert.equal(metadata.settingsPending, true); assert.equal(metadata.time_s, 12.5);
    assert.equal(metadata.manifest.protocol_hash, 'synthetic-browser-fixture');
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
    assert.match(await page.title(), /6 cells · SimpleCa/);
    assert.deepEqual(problems, []);
    console.log('PASS: movie preview and pointer scrubbing, Play/manual seek, static view, stale/error frames, portable protocols/precision, primary/control results, irregular-time heatmap, pending settings, exports and responsive layouts');
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

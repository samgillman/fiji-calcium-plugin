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
    assert.equal(await page.locator('#cell option').count(), 6, 'cells can be selected before analysis');
    assert.ok(await page.locator('#windows-editor details').evaluateAll((details) => details.every((d) => !d.open)), 'advanced response rules start collapsed');
    assert.ok(await page.evaluate(() => {
      const transport = document.getElementById('transport').getBoundingClientRect();
      const picture = document.getElementById('picture').getBoundingClientRect();
      return transport.bottom <= picture.top;
    }), 'the frame slider is above the recording');
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
    assert.match(await text('#validation'), /F₀ frames 10–11 come after the stimulus\./);
    await page.locator('#validation .repair', { hasText: 'Use 1–9' }).click();
    assert.equal(await page.inputValue('#baseline_end'), '9'); assert.doesNotMatch(await text('#validation'), /F₀ frame/);
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
      assert.match(await text('#validation'), /Frame rate/);
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
    assert.equal(await page.inputValue('#view'), 'raw', 'analysis preserves the selected movie view');
    await page.selectOption('#cell', 'Cell001'); assert.match(await text('#verdict'), /Cell001: Responded/);
    // The duration notice is set on every cell by the engine but belongs to the recording: it is shown once and flags no cell.
    assert.deepEqual(await page.locator('#group option').allTextContents(), ['All cells (6)', 'Responded (2)', 'Did not respond (1)', 'Excluded (2)', 'Unscorable (1)', 'Needs quality review (1)']);
    assert.match(await text('#quality-summary'), /1 cell flagged, including 1 responder$/);
    await page.selectOption('#cell', 'Cell005'); assert.match(await text('#cell-quality'), /Low F₀/);
    await page.locator('#quality-review summary').click(); assert.match(await text('#quality-list'), /Applies to the whole recording: Minimum response duration is shorter than one frame interval/);
    await page.locator('#quality-list .flag-filter').first().click();
    assert.equal(await page.inputValue('#group'), '5'); assert.deepEqual(await page.locator('#cell option').evaluateAll((o) => o.map((x) => x.value)), ['Cell005']);
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
    await page.locator('#picture-note').scrollIntoViewIfNeeded();
    assert.ok(await page.locator('#picture-note').evaluate((e) => {const a=e.getBoundingClientRect(),b=e.closest('.card').getBoundingClientRect();return a.bottom<=b.bottom&&a.bottom<=innerHeight;}));
    await page.screenshot({path:path.join(out,'page-short-window.png')});
    await page.evaluate(() => document.getElementById('lost').classList.remove('on'));
    await page.setViewportSize({width:1440,height:900});


    // Control calls, thresholds and headline counts follow the selected response window.
    await page.selectOption('#response-window', 'Control::Ionomycin');
    assert.match(await text('#verdict'), /Cell001: Did not respond/);
    assert.ok(await page.evaluate(() => Math.abs(window.calciumStep4.referenceLevels().find((r) => r.name === 'required rise').value - 1.1) < 1e-12), 'local .4 threshold on baseline .5 belongs at global 1.1');
    assert.match(await text('#tile-responding'), /Ionomycin \(control\): 1 of 3 evaluable/);
    await page.selectOption('#group', '1');
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value)), ['Cell002']);
    await page.click('#tile-responding');
    assert.equal(await page.inputValue('#response-window'), 'Stimulus::Capsaicin');
    assert.match(await text('#tile-responding'), /Capsaicin: 2 of 3 evaluable/);
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value)), ['Cell001', 'Cell005']);
    await page.click('#tile-excluded');
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value)), ['Cell003', 'Cell004']);
    assert.ok(await page.locator('#tile-excluded').evaluate((e) => e.classList.contains('on')));
    assert.ok(!(await page.locator('#tile-unscorable').evaluate((e) => e.classList.contains('on'))));
    await page.click('#tile-unscorable');
    assert.deepEqual(await page.locator('#cell option').evaluateAll((options) => options.map((o) => o.value)), ['Cell006']);
    await page.click('#tile-registered'); await page.selectOption('#cell', 'Cell002');

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
    await page.click('.tab[data-tab=metrics]'); await page.selectOption('#metrics-table', 'groups');
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
    await page.click('#save-plot'); await downloadedPlots;
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

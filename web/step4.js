// Step 4 page. Everything shown here is read from the tables the analysis saved, or from the recording's own frames.
// The page draws; it does not calculate results. Settings go back to the plugin, which runs the same engine as the Fiji window.
'use strict';
(function () {
  const $ = (id) => document.getElementById(id);
  const CALL = { responded: '#177465', not: '#C8660E', other: '#8A94A0' };
  const ON_PICTURE = { responded: '#3ddbb5', not: '#ffa24a', other: '#aeb8c2' };
  const GROUPS = ['All cells', 'Responded', 'Did not respond', 'Excluded', 'Unscorable', 'Needs quality review', 'Marked'];   // Marked: offered only while review marks apply
  /** The colour overview and the background areas are still pictures; every other view is a frame of the recording. */
  const still = () => S.view === 'picture' || S.view === 'background';
  const kind = (call) => call === 'Responded' ? 'responded' : call === 'Did not respond' ? 'not' : 'other';
  const inGroup = (call, group, name) => group === 0 || (group === 1 && call === 'Responded') || (group === 2 && call === 'Did not respond') ||
    (group === 3 && call === 'Excluded') || (group === 4 && call === 'Unscorable') || (group === 5 && qualityFor(name).length > 0 && (!S.qualityFlag || qualityFor(name).includes(S.qualityFlag))) ||
    (group === 6 && !!markOf(name));
  const S = { unlinked: new Set(), validationShown: null, imported: null, info: null, state: null, results: null, version: -1, group: 0, cell: null, listed: [], keep: null, frame: 1, view: 'raw', viewChosen: false,
    every: false, thresholds: false, plotMode: 'mean', signal: 'dff', outlines: true, outlineWidth: 1, outlineColours: 'call', numberCells: false, saveArea: 'whole', timer: null, lost: 0, note: null, polling: false, paths: [], boxes: [], draft: null, draftAck: '', revision: 0, draftTimer: null, draftQueue: Promise.resolve(), window: null, saving: false, initialDraft: null, protocolDrafts: {}, quality: {}, metrics: {}, closed: false,
    marking: null, reviewing: false, carried: null, said: null, finishing: false };   // the open mark form ({ cell }); a mark on its way to Fiji; the carried marks of this run; what Fiji said after the last mark; results read before the run was over (poll)
  const Plot = { dirty: true, layer: document.createElement('canvas'), x0: null, x1: null, y: {}, hover: null, at: null, readout: null, line: false, movedAt: -1e9, rest: null, geom: null, scrubbing: false, drawn: [], series: [] };
  const LINGER = 1500;   // ms the plot keeps showing where the picture is after the frame last moved
  const HINT_TIMES = 3, HINT_READ = 1000;   // how often the line that says how to zoom and move is shown, over all sessions; a showing counts once it has been up this long (ms)
  // The picture. Its view: zoom (1 is the whole picture) and the point of the recording in the middle (cx, cy; see setView). What the
  // pointer is doing: drag (a press that may become a move or a rectangle), moved, held, region ("Zoom to region" is on), band (the
  // rectangle being drawn), over and pointer (whether and where the pointer is on the picture), ctrl (the Ctrl key is held), pinch (the
  // zoom at which a pinch sent as a gesture began), hinted (how often the hint was shown). map: the frame as the small map last drew it.
  const Pic = { base: null, shown: null, zoom: 1, cx: null, cy: null, drag: null, moved: false, held: null, shownFrame: null, error: null, region: false, band: null, over: false, pointer: null, ctrl: false, pinch: null, hinted: 0, map: null };
  // Movie frames: compressed images are kept for the session (about 40 MB for 300 full-size frames, much less at display size),
  // decoded ones only near the playhead. Keys carry the server's frame version, which changes with the baseline behind
  // "Change from baseline" and the display range, so a stale image is never reused.
  const Buf = { version: null, blobs: new Map(), bitmaps: new Map(), inflight: new Map(), decoding: new Map(), failed: new Map(), bytes: 0, closed: '' };   // closed: what Fiji says when the recording was closed there
  const PARALLEL = 4, DECODED_AROUND = 30, MAX_BYTES = 90e6, MAX_FPS = 30;
  Object.defineProperty(Pic, 'loading', { get: () => Buf.inflight.size > 0 || Buf.decoding.size > 0 });
  // The playback clock runs in recording seconds; frames are shown when the playhead reaches their timestamps.
  const Clock = { t: 0, last: 0, step: 1, times: [], real: false, waitingSince: 0 };
  const Heat = { geom: null };
  const EMPTY_PLOT = $('plot-empty').textContent;   // what the plot says before an analysis

  async function api(path, body) {
    const options = body === undefined ? { cache: 'no-store' } : { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) };
    const response = await fetch('api/' + path, options);
    const data = await response.json().catch(() => ({}));
    if (!response.ok) throw new Error(data.error || 'The plugin answered ' + response.status);
    return data;
  }
  const finite = (v) => typeof v === 'number' && isFinite(v);
  const calls = () => S.results ? ((S.results.callsByWindow || {})[S.window] || S.results.calls || {}) : {};
  const callRow = (name) => calls()[name] || {};
  /** The review marks of the analysis shown (results.review.marks), by cell; null while it has none. */
  const marked = () => { const marks = S.results && S.results.review && S.results.review.marks; return marks && Object.keys(marks).length ? marks : null; };
  /** Whether the response window on screen is the primary stimulus, the one the review marks are for. */
  const primaryShown = () => { const w = primaryWindow(); return !w || w.key === S.window; };
  /** A cell's call as the page shows it in a response window (the one on screen when none is given): after review for the
   *  primary stimulus while review marks apply (results.review.after), the engine's own for a control window. Everything that
   *  counts, groups, colours or lists cells by call reads it. */
  function shownCall(name, key) {
    const r = S.results; if (!r) return undefined;
    if (key === undefined) key = S.window;
    const primary = primaryWindow(), after = marked() && (!primary || key === primary.key) ? (r.review.after || {})[name] : undefined;
    return after || engineCall(name, key);
  }
  /** The engine's own call of a cell in a response window, as saved. */
  const engineCall = (name, key) => ((((S.results.callsByWindow || {})[key] || S.results.calls || {})[name]) || {}).Status;
  /** The engine's call of a cell for the primary stimulus: the automatic call a review mark is made against. */
  const automaticCall = (name) => engineCall(name, (primaryWindow() || {}).key);
  const callOf = (name) => shownCall(name);
  /** The call of a cell for the primary stimulus, after review while marks apply: what the per-cell tables and the summary are split by. */
  const primaryCall = (name) => shownCall(name, (primaryWindow() || {}).key);
  /** The review mark of a cell while the primary stimulus is on screen, or null. */
  const markOf = (name) => { const marks = marked(); return marks && name && marks[name] && primaryShown() ? marks[name] : null; };
  const number = (text) => { const v = parseFloat(text); return isFinite(v) ? v : NaN; };
  /** A number as somebody types it into a box: digits with at most one decimal mark (a point or a comma) and a sign before
   *  them. null: the box is empty. NaN: anything else, also "1e1" and "0x10", which a program reads as numbers and a person does not. */
  function typedNumber(text) { const t = String(text).trim().replace(',', '.').replace('−', '-'); return t === '' ? null : /^[-+]?(\d+\.?\d*|\.\d+)$/.test(t) ? Number(t) : NaN; }
  /** A whole number as typed (typedNumber), for a box that counts frames; NaN for anything else, also "2.5", never rounded. */
  const typedWhole = (text) => { const v = typedNumber(text); return Number.isInteger(v) ? v : NaN; };

  /** Sets a canvas to the size it is shown at (sharp on high-density screens) and returns a context in CSS pixels. */
  function fit(canvas) {
    const dpr = window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
    const ctx = canvas.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h, dpr };
  }
  /** Round positions for axis ticks. */
  function ticks(lo, hi, about) {
    if (!(hi > lo)) return { values: [lo], step: 1 };
    const raw = (hi - lo) / Math.max(1, about), power = Math.pow(10, Math.floor(Math.log10(raw))), unit = raw / power;
    const step = (unit <= 1.5 ? 1 : unit <= 3 ? 2 : unit <= 7 ? 5 : 10) * power, values = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) values.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return { values, step };
  }
  /** A tick's number, with as many decimals as its spacing needs and a true minus sign (−, not a hyphen). */
  const tickText = (v, step) => v.toFixed(step >= 1 ? 0 : Math.min(6, Math.ceil(-Math.log10(step) - 1e-9))).replace('-', '−');
  /** The value below which this share of the finite numbers lies. */
  function share(rows, part, absolute) {
    const all = [];
    for (const row of rows) for (const v of row) if (finite(v)) all.push(absolute ? Math.abs(v) : v);
    if (!all.length) return 0;
    all.sort((a, b) => a - b);
    return all[Math.max(0, Math.min(all.length - 1, Math.round(part * (all.length - 1))))];
  }
  function tip(text, x, y) {
    const box = $('tip');
    if (!text) { box.hidden = true; return; }
    box.textContent = text; box.hidden = false;
    const w = box.offsetWidth; box.style.left = Math.max(4, Math.min(window.innerWidth - w - 4, x + 14)) + 'px'; box.style.top = (y + 16) + 'px';
  }

  // ---------------------------------------------------------------- editable draft, applied results, and protocol exchange
  const clone = (v) => JSON.parse(JSON.stringify(v));
  const fmt = (v) => finite(Number(v)) && v !== '' && v !== null ? String(+Number(v).toPrecision(7)) : (v == null ? '' : String(v));
  /** The frame rate as shown: four significant figures when it has more, so a rate read from the file
   *  (1.0000537066691941) reads 1.000. The stored value keeps every digit. */
  const rateText = (value) => { const n = Number(value); if (value == null || value === '' || !finite(n)) return value == null ? '' : String(value);
    return String(+n.toPrecision(4)) === String(n) ? String(n) : n >= 1000 ? String(+n.toPrecision(4)) : n.toPrecision(4); };
  /** The samples in a row the engine needs for a length given in seconds: the first sample plus every interval that fits
   *  (ceil(length / median interval) + 1, as the engine does). 1 s is 2 samples at exactly 1 frame per second and 3 at
   *  0.99995 s per frame. null when the frame times are not known. */
  function samplesFor(seconds, times) {
    if (!finite(seconds) || seconds <= 0 || !times || times.length < 2) return null;
    const gaps = []; for (let i = 1; i < times.length; i++) gaps.push(times[i] - times[i - 1]);
    if (!gaps.every((g) => finite(g) && g > 0)) return null;
    gaps.sort((a, b) => a - b); const m = gaps.length >> 1, step = gaps.length % 2 ? gaps[m] : (gaps[m - 1] + gaps[m]) / 2, ratio = seconds / step;
    return Math.ceil(ratio - 16 * Number.EPSILON * Math.max(1, Math.abs(ratio))) + 1;
  }
  function exactInput(box, value, scale) {
    const shown = box.id === 'sampling_hz' ? rateText(value) : scale ? fmt(Number(value) * scale) : fmt(value);
    box.value = shown; box.dataset.shown = shown; box._exact = value; box._scale = scale || 1;
    box.title = value == null ? '' : 'Stored value: ' + value + '. Untouched fields keep their full precision.';
  }
  function readInput(box, optional) {
    if (box.value === box.dataset.shown) return box._exact;
    const text = box.value.trim(); if (optional && text === '') return null;
    const value = text === '' ? NaN : Number(text);
    if (!Number.isFinite(value)) return text;
    // A typed time that reads the same as a frame's time is that frame: 113.091 means the frame stored as 113.09100770950317,
    // so a retyped value never differs from an untouched one by a hidden fraction.
    if (TIME_KEYS.includes(box.dataset.key)) { const shown = fmt(value); for (const t of editedTimes() || []) if (fmt(t) === shown) return t; }
    return value / box._scale;
  }
  const TIME_KEYS = ['Baseline_start_s', 'Baseline_end_s', 'Stimulus_onset_s', 'Start_s', 'End_s'];
  const FIELD_TITLES = { Baseline_start_s: 'baseline start', Baseline_end_s: 'baseline end', Stimulus_onset_s: 'stimulus time', Start_s: 'response window start', End_s: 'response window end',
    Floor_local_fraction: 'minimum rise', SD_multiplier: 'baseline SD multiplier', Min_duration_s: 'minimum duration', Consecutive_samples: 'consecutive samples', Max_gap_s: 'maximum gap',
    Min_baseline_samples: 'minimum baseline samples', Noise_floor_local_fraction: 'noise floor', Expected_interval_s: 'expected interval' };
  const seconds = (v) => fmt(v) + ' s';
  /** Frame times of the settings being edited: from the frame rate, or the timestamps as the plugin read them. Null when unknown. */
  function draftTimes(timestamps, hz) {
    const n = S.info ? S.info.frames : 0; if (!n) return null;
    if (!timestamps) { const rate = hz === null || hz === '' || typeof hz === 'boolean' ? NaN : Number(hz); return Number.isFinite(rate) && rate > 0 ? Array.from({ length: n }, (_, i) => i / rate) : null; }
    if (!S.results) { const t = S.state && S.state.frameTimes; return t && t.length >= n ? t : null; }
    const applied = S.results.appliedDraft && S.results.appliedDraft.timestamps, t = S.results.frameTimes || S.results.time;
    return applied && JSON.stringify(applied) === JSON.stringify(timestamps) && t && t.length >= n ? t : null;
  }
  const editedTimes = () => (S.draft && $('sampling_hz') ? draftTimes(S.draft.timestamps, readInput($('sampling_hz'), false)) : null);
  /** Which frames a window uses, by the engine's own rule: baseline start ≤ t < baseline end, and scoring start ≤ t ≤ scoring end. */
  function windowFrames(row, times) {
    const n = S.info ? S.info.frames : 0; if (!n || !times || times.length < n) return null;
    const num = (v) => (v === null || v === '' || typeof v === 'boolean' ? NaN : Number(v));
    const b0 = num(row.Baseline_start_s), b1 = num(row.Baseline_end_s), on = num(row.Stimulus_onset_s), s0 = num(row.Start_s), s1 = num(row.End_s);
    const tol = 8 * Number.EPSILON * Math.max(1, ...[b0, b1, s0, s1].filter(Number.isFinite).map(Math.abs));
    const span = (test) => { let first = 0, last = 0, count = 0; for (let i = 0; i < n; i++) if (test(times[i])) { if (!first) first = i + 1; last = i + 1; count++; } return { first, last, count }; };
    let before = 0, at = false; if (Number.isFinite(on)) for (let i = 0; i < n; i++) { if (times[i] < on - tol) before = i + 1; else { at = Math.abs(times[i] - on) <= tol; break; } }
    return { n, times, stimulus: Number.isFinite(on) ? { before, at } : null,
      baseline: Number.isFinite(b0) && Number.isFinite(b1) ? span((t) => t >= b0 - tol && t < b1 - tol) : null,
      scoring: Number.isFinite(s0) && Number.isFinite(s1) ? span((t) => t >= s0 - tol && t <= s1 + tol) : null };
  }
  /** The timing is entered in frames: one fact, "the stimulus was added after frame N". F₀ frames, the baseline and the response
   *  window follow from it until a row is set separately. The seconds the engine needs are taken from the recording's own frame
   *  times, never from a typed round number: the stimulus is the time of the first frame after it; a baseline runs up to, not
   *  including, the frame after its last frame; the response window runs from its first frame to its last, both included.
   *  The hidden boxes named window-<i>-<key> hold those seconds exactly as before, so the protocol file is unchanged. */
  // F₀, and a baseline of its own: the last 20 frames before the stimulus, or all of them when there are fewer; frames in which the
  // field moved are left out when at least 3 frames remain after them
  const frameRule = (N) => { let a = Math.max(1, N - 19); const moved = movedWithin(a, N); if (moved.count && N - moved.after + 1 >= 3) a = moved.after; return [a, N]; };
  const wholeIn = (box) => { const text = box ? String(box.value).trim() : '', v = Number(text); return text !== '' && Number.isInteger(v) ? v : NaN; };
  const frameBox = (id, index) => $(id + '-' + index);
  function applyFrames(index) {
    const times = editedTimes(), n = S.info ? S.info.frames : 0, i = String(index), get = (id) => frameBox(id, i), legacy = Number($('protocol-mode').value) === 1;
    if (!times || !get('stim-frame')) return;
    const put = (key, value) => exactInput($('window-' + i + '-' + key), value);
    const N = wholeIn(get('stim-frame')), okN = N >= 1 && N < n;
    if (i === '0' && okN && $('f0-cell').dataset.linked === '1') { const f = frameRule(N); exactInput($('baseline_start'), f[0]); exactInput($('baseline_end'), f[1]); }
    if (okN && get('base-cell').dataset.linked === '1') {
      const own = [wholeIn($('baseline_start')), wholeIn($('baseline_end'))], f = i === '0' && own[0] >= 1 && own[1] >= own[0] && own[0] <= N ? [own[0], Math.min(own[1], N)] : frameRule(N);
      get('base-from').value = f[0]; get('base-to').value = f[1];
    }
    if (okN && get('resp-cell').dataset.linked === '1') { get('resp-from').value = N + 1; get('resp-to').value = n; }
    if (legacy && okN) { get('base-to').value = N; get('resp-from').value = N + 1; }
    const a = wholeIn(get('base-from')), b = wholeIn(get('base-to')), c = wholeIn(get('resp-from')), d = String(get('resp-to').value).trim() === '' ? n : wholeIn(get('resp-to'));
    const at = (frame) => (frame >= 1 && frame <= n ? times[frame - 1] : null);
    put('Stimulus_onset_s', okN ? times[N] : null);
    put('Baseline_start_s', at(a)); put('Baseline_end_s', b >= 1 && b < n ? times[b] : null);
    put('Start_s', at(c)); put('End_s', at(d));
  }
  /** Frames in which the field was displaced (found by the plugin when the recording was opened). */
  const movedRuns = () => (S.info && S.info.motion && S.info.motion.moved > 0 ? S.info.motion.runs : []);
  const movedText = () => 'frames ' + movedRuns().map((r) => (r[0] === r[1] ? r[0] : r[0] + '–' + r[1])).join(' and ');
  const isMoved = (frame) => movedRuns().some((r) => frame >= r[0] && frame <= r[1]);
  /** How many frames of a..b moved, and the first frame after the last of them inside a..b (0 when none moved). */
  function movedWithin(a, b) { let count = 0, last = 0; for (let f = a; f <= b; f++) if (isMoved(f)) { count++; last = f; } return { count, after: count ? last + 1 : 0 }; }
  function applyAllFrames() { for (const section of document.querySelectorAll('.window-editor')) applyFrames(section.dataset.index); }
  /** Fills the frame boxes from the seconds of a loaded window, by the engine's own rule, and decides which rows follow the stimulus. */
  function showFrames(index, row) {
    const i = String(index), get = (id) => frameBox(id, i), used = windowFrames(row, editedTimes()), n = S.info.frames;
    const N = used && used.stimulus ? used.stimulus.before : NaN, okN = N >= 1 && N < n, has = (s) => s && s.count > 0;
    get('stim-frame').closest('.window-editor').dataset.framesShown = used ? '1' : '';   // with recorded timestamps the frame times arrive a moment later
    get('stim-frame').value = Number.isFinite(N) ? N : '';
    get('base-from').value = used && has(used.baseline) ? used.baseline.first : ''; get('base-to').value = used && has(used.baseline) ? used.baseline.last : '';
    get('resp-from').value = used && has(used.scoring) ? used.scoring.first : ''; get('resp-to').value = used && has(used.scoring) ? used.scoring.last : '';
    const rule = okN ? frameRule(N) : null, same = (pair, a, b) => !!pair && pair[0] === a && pair[1] === b, last = used ? used.times[n - 1] : NaN, end = Number(row.End_s);
    const f0 = [wholeIn($('baseline_start')), wholeIn($('baseline_end'))];
    const set = (cellId, key, follows) => { const cell = i === '0' && cellId === 'f0-cell' ? $('f0-cell') : get(cellId); cell.dataset.linked = follows && !S.unlinked.has(key) ? '1' : ''; };
    if (i === '0') set('f0-cell', 'f0', !okN || same(rule, f0[0], f0[1]));
    const baseRule = okN ? (i === '0' && f0[0] >= 1 && f0[1] >= f0[0] && f0[0] <= N ? [f0[0], Math.min(f0[1], N)] : rule) : null;
    set('base-cell', 'base-cell-' + i, !okN || !has(used.baseline) && row.Baseline_start_s == null || has(used.baseline) && same(baseRule, used.baseline.first, used.baseline.last));
    set('resp-cell', 'resp-cell-' + i, !okN || !has(used.scoring) && row.Start_s == null || has(used.scoring) && used.scoring.first === N + 1 && used.scoring.last === n && !(Number.isFinite(end) && end > last + 1e-9));
  }
  /** What each row that follows the stimulus says, and the same windows in seconds on one line. */
  function paintFrames(index) {
    const i = String(index), get = (id) => frameBox(id, i), times = editedTimes(), n = S.info.frames; if (!get('stim-frame')) return;
    const pair = (from, to) => { const a = wholeIn(from), b = wholeIn(to); return Number.isFinite(a) && Number.isFinite(b) ? [a, b] : null; }, text = (p) => (p ? p[0] + '–' + p[1] : 'not set');
    const f0 = pair($('baseline_start'), $('baseline_end')), base = pair(get('base-from'), get('base-to')), resp = pair(get('resp-from'), get('resp-to')), N = wholeIn(get('stim-frame'));
    const chip = (cell, words, hint) => { const button = cell.querySelector('.link-toggle'); button.textContent = words; button.title = hint + ' Click to set it yourself.'; };
    if (i === '0') chip($('f0-cell'), text(f0), 'Follows the stimulus: the frames before it, or the last 20 of a long wait.');
    chip(get('base-cell'), i === '0' && base && f0 && base[0] === f0[0] && base[1] === f0[1] ? 'same as F₀' : text(base), 'Follows the stimulus.');
    chip(get('resp-cell'), resp ? resp[0] + '–' + resp[1] + (resp[1] === n ? ' (end)' : '') : 'not set', 'Follows the stimulus: from the first frame after it to the last frame.');
    const s = (frame) => (times && frame >= 1 && frame <= n ? fmt(+times[frame - 1].toFixed(2)) : '?');
    $('frames-preview-' + i).textContent = !times ? 'Enter the frame rate to set the timing.' : !(N >= 1 && N < n) ? '' :
      'In seconds: baseline ' + (base ? s(base[0]) + '–' + s(base[1]) : '?') + ' · stimulus ' + s(N + 1) + ' · response ' + (resp ? s(resp[0]) + '–' + s(resp[1]) : '?');
  }
  const FILES = ['timestamps', 'sensitivity', 'negative_control'];
  const NUMBERS = ['Baseline_start_s', 'Baseline_end_s', 'Stimulus_onset_s', 'Start_s', 'End_s', 'Floor_local_fraction', 'SD_multiplier', 'Min_duration_s', 'Consecutive_samples', 'Max_gap_s', 'Min_baseline_samples', 'Noise_floor_local_fraction', 'Expected_interval_s'];
  function expandWindow(row, legacy) {
    const w = clone(row);
    if (legacy) Object.assign(w, { Baseline_start_s: Number(w.Start_s) - Number(w.Baseline_s), Baseline_end_s: w.Start_s, Stimulus_onset_s: w.Start_s, Duration_mode: 'samples', Min_duration_s: 1, Min_baseline_samples: 3, Noise_floor_local_fraction: 0, Expected_interval_s: null });
    return w;
  }
  function windowInput(index, key, title, opts) {
    const options = opts || {}, id = 'window-' + index + '-' + key;
    const label = document.createElement('label'); label.className = 'row';
    const name = document.createElement('span'); name.textContent = title; label.appendChild(name);
    const box = document.createElement(options.choices ? 'select' : 'input'); box.id = id; box.dataset.window = index; box.dataset.key = key;
    if (options.choices) for (const [value, text] of options.choices) { const o = document.createElement('option'); o.value = value; o.textContent = text; box.appendChild(o); }
    else if (key === 'Stimulus') box.className = 'wide'; else box.inputMode = 'decimal';
    label.appendChild(box); return label;
  }
  function renderWindows() {
    // The positive-control tick box is moved between the two windows; keep references, because a detached element cannot be found by id.
    S.validationShown = null;
    const container = $('windows-editor'), keep = [S.controlLabel || (S.controlLabel = $('control').closest('label')), S.controlHelp || (S.controlHelp = $('control-help'))]; keep.forEach((e) => e.remove()); container.innerHTML = '';
    const legacy = Number($('protocol-mode').value) === 1;
    S.draft.protocol.windows.forEach((row, index) => {
      const w = expandWindow(row, legacy), section = document.createElement('div'); section.className = 'window-editor'; section.dataset.index = index;
      const heading = document.createElement('h3'); heading.textContent = index === 0 ? 'Primary stimulus' : 'Positive control'; section.appendChild(heading);
      const add = (target, key, title, opts) => { const label = windowInput(index, key, title, opts); target.appendChild(label); const box = label.lastChild;
        if (box.tagName === 'SELECT' || key === 'Stimulus') box.value = w[key] == null ? '' : String(w[key]); else exactInput(box, w[key], opts && opts.scale); return box; };
      const frames = (id, label) => { const box = document.createElement('input'); box.id = id; box.inputMode = 'numeric'; box.dataset.frames = '1'; box.dataset.window = index; box.setAttribute('aria-label', label); return box; };
      const line = (title, hint) => { const r = document.createElement('div'); r.className = 'row'; const name = document.createElement('span'); name.textContent = title; name.title = hint || ''; r.appendChild(name); section.appendChild(r); return r; };
      // A row that follows the stimulus shows its frames as a dashed chip. A click shows two boxes instead; "auto" puts it back.
      const derived = (title, hint, cellId, key, first, last) => {
        const r = line(title, hint), cell = document.createElement('span'); cell.className = 'cell derived'; cell.id = cellId; cell.dataset.linked = '1'; r.appendChild(cell);
        const chip = document.createElement('button'); chip.type = 'button'; chip.className = 'link-toggle'; cell.appendChild(chip);
        const pair = document.createElement('span'); pair.className = 'pair'; const to = document.createElement('i'); to.textContent = 'to'; pair.append(first, to, last);
        const back = document.createElement('button'); back.type = 'button'; back.className = 'relink'; back.textContent = 'auto'; back.title = 'Follow the stimulus again'; pair.appendChild(back); cell.appendChild(pair);
        chip.addEventListener('click', () => { S.unlinked.add(key); cell.dataset.linked = ''; first.focus(); first.select(); });
        back.addEventListener('click', () => { S.unlinked.delete(key); cell.dataset.linked = '1'; applyFrames(index); draftChanged(); });
      };
      const holder = document.createElement('span'); holder.className = 'with-button'; line('Stimulus after frame', 'The last frame before the stimulus was added. The rows below follow it.').appendChild(holder);
      const stimBox = frames('stim-frame-' + index, 'Stimulus added after frame'); holder.appendChild(stimBox);
      const fromFrame = document.createElement('button'); fromFrame.type = 'button'; fromFrame.className = 'from-frame'; fromFrame.textContent = 'Use frame';
      fromFrame.addEventListener('click', () => { stimBox.value = S.frame; applyFrames(index); draftChanged(); });
      holder.appendChild(fromFrame);
      const suggestion = document.createElement('p'); suggestion.className = 'fine suggestion'; suggestion.id = 'suggestion-' + index; suggestion.hidden = true; section.appendChild(suggestion);
      if (index === 0) {
        const first = frames('baseline_start', 'First F₀ frame'), last = frames('baseline_end', 'Last F₀ frame'); exactInput(first, S.draft.values.baseline_start); exactInput(last, S.draft.values.baseline_end);
        derived('F₀ frames', 'These frames set the zero of ΔF/F₀', 'f0-cell', 'f0', first, last);
      }
      derived('Baseline', 'The frames before the stimulus that set the noise level a response must exceed', 'base-cell-' + index, 'base-cell-' + index, frames('base-from-' + index, 'First baseline frame'), frames('base-to-' + index, 'Last baseline frame'));
      derived('Response window', 'The frames in which a response is looked for', 'resp-cell-' + index, 'resp-cell-' + index, frames('resp-from-' + index, 'First frame of the response window'), frames('resp-to-' + index, 'Last frame of the response window'));
      if (index === 0 && S.info && S.info.fijiFrames) { const chosen = document.createElement('p'); chosen.className = 'fine fiji-frames'; chosen.id = 'fiji-frames'; chosen.textContent = 'Cells chosen in Fiji on frames ' + S.info.fijiFrames.first + '–' + S.info.fijiFrames.last + ', before the response.'; chosen.title = 'Blind selection uses the frames before any cell responds. They need not end where the stimulus was added.'; section.appendChild(chosen); }
      // the same timing in seconds, exactly as the engine gets it: filled from the frames, not shown
      const exact = document.createElement('div'); exact.hidden = true; exact.className = 'exact'; for (const key of TIME_KEYS) add(exact, key, key); section.appendChild(exact);
      const preview = document.createElement('p'); preview.className = 'frames-preview'; preview.id = 'frames-preview-' + index; section.appendChild(preview);
      const details = document.createElement('details'); details.className = 'rule'; const summaryTitle = document.createElement('summary'); summaryTitle.title = 'The response rule. Click to adjust it.';
      const summary = document.createElement('span'); summary.className = 'rule-summary'; summary.id = 'rule-summary-' + index; summaryTitle.appendChild(summary); details.appendChild(summaryTitle);
      add(details, 'Floor_local_fraction', 'Minimum rise (%)', { scale: 100 }); add(details, 'SD_multiplier', 'Baseline SD multiplier');
      add(details, 'Duration_mode', 'Duration measured in', { choices: [['seconds', 'Seconds'], ['samples', 'Consecutive samples']] });
      add(details, 'Min_duration_s', 'Minimum duration (s)'); add(details, 'Consecutive_samples', 'Consecutive samples');
      add(details, 'Min_baseline_samples', 'Baseline samples').previousSibling.title = 'The fewest baseline samples a cell needs for a call'; add(details, 'Max_gap_s', 'Maximum gap (s)').placeholder = 'auto';
      add(details, 'Expected_interval_s', 'Expected interval (s)').placeholder = 'auto'; add(details, 'Noise_floor_local_fraction', 'Noise floor (%)', { scale: 100 });
      add(details, 'Stimulus', 'Name in the output');
      section.appendChild(details); container.appendChild(section); showFrames(index, w);
      if (index === 0) keep.forEach((e) => container.appendChild(e));
    });
    if (!S.draft.protocol.windows.length) keep.forEach((e) => container.appendChild(e));
    frameButtons();
    $('control').checked = S.draft.protocol.windows.length > 1;
    $('control-help').hidden = !$('control').checked;
    updateProtocolControls();
  }
  function updateProtocolControls() {
    const legacy = Number($('protocol-mode').value) === 1, busy = !!(S.state && S.state.busy) || S.saving;
    $('protocol-name').disabled = legacy || busy; $('protocol-validation').disabled = legacy || busy;
    $('field-overrides').disabled = legacy || busy;
    for (const key of ['sensitivity', 'negative_control']) $(key + '-file').disabled = legacy || busy;
    for (const section of document.querySelectorAll('.window-editor')) {
      const i = section.dataset.index, box = (key) => $('window-' + i + '-' + key), mode = box('Duration_mode').value;
      for (const field of section.querySelectorAll('input,select')) {
        const key = field.dataset.key;
        field.disabled = busy || (legacy && (['Baseline_end_s', 'Start_s', 'Duration_mode', 'Min_duration_s', 'Min_baseline_samples', 'Noise_floor_local_fraction', 'Expected_interval_s'].includes(key) || (i === '1' && key === 'Stimulus'))) || (!legacy && ((key === 'Min_duration_s' && mode !== 'seconds') || (key === 'Consecutive_samples' && mode !== 'samples')));
      }
      for (const button of section.querySelectorAll('button')) button.disabled = busy || button.dataset.fixed === '1';
      const row = {}; for (const key of TIME_KEYS) row[key] = readInput(box(key), false);
      // the boxes were empty until the frame times arrived: check again at once, so "enter the frame" does not stay up until the next poll
      if (!section.dataset.framesShown && editedTimes()) { showFrames(i, row); if (section.dataset.framesShown) queueMicrotask(() => { if (S.state) showState(); }); }
      const used = windowFrames(row, editedTimes());
      for (const id of ['base-to', 'resp-from']) frameBox(id, i).disabled = busy || legacy;
      paintFrames(i);
      // where the starting stimulus frame came from, until it is changed
      const note = $('suggestion-' + i), first = S.initialDraft && S.initialDraft.protocol.windows[i], from = S.info && S.info.stimulusFrom;
      const untouched = i === '0' && !S.results && !!first && typeof row.Stimulus_onset_s === 'number' && Number(first.Stimulus_onset_s) === row.Stimulus_onset_s;
      const words = !untouched ? '' : from === 'last-session' ? 'From your last session on this recording.' : from === 'last-recording' ? 'Same as your last recording.' : from === 'fiji' ? 'A placeholder from the frames used in Fiji. Set it to when the stimulus was added.' : 'A suggestion, not a measured time.';
      if (note.dataset.shown !== words) {
        note.dataset.shown = words; note.hidden = !words; note.replaceChildren(words);
        if (words && from === 'last-session' && S.info.fijiDraft) { const use = document.createElement('button'); use.type = 'button'; use.className = 'ghost repair'; use.textContent = 'Use Fiji’s values';
          use.addEventListener('click', () => { S.unlinked.clear(); S.initialDraft = clone(S.info.fijiDraft); S.info.stimulusFrom = 'fiji'; fillDraft(S.info.fijiDraft); $('confirmed').checked = false; draftChanged(); }); note.append(' ', use); }
      }
      if (i === '0') shadeTimeline(used);
      // the label says what is applied: samples in a row, also when the length was typed in seconds
      const inRow = (n) => n + (String(n) === '1' ? ' sample' : ' samples in a row'), needed = mode === 'seconds' ? samplesFor(Number(box('Min_duration_s').value), editedTimes()) : null;
      $('rule-summary-' + i).textContent = 'Rule: rise ≥ ' + box('Floor_local_fraction').value + '% and ≥ ' + box('SD_multiplier').value + ' × SD, for ' + (mode === 'seconds' ? box('Min_duration_s').value + ' s' + (needed ? ', here ' + inRow(needed) : '') : inRow(box('Consecutive_samples').value)) + (legacy ? ' (legacy)' : '');
    }
    $('rule-help').textContent = legacy ? 'Legacy rule: the baseline ends at the response start; the control is KCl. Sample counts depend on the frame rate.' : 'Starting values. Check against your controls.';
  }
  /** Under the frame slider: grey for the noise baseline, a mark for the stimulus, orange for the response window of the primary stimulus. */
  function shadeTimeline(used) {
    const strip = $('time-windows'); if (!strip) return; strip.replaceChildren(); if (!used || used.n < 2) return;
    const at = (frame) => (frame - 1) / (used.n - 1) * 100, piece = (className, from, to, title) => { const e = document.createElement('i'); e.className = className; e.style.left = at(from) + '%'; e.style.width = Math.max(0.4, at(to) - at(from)) + '%'; e.title = title; strip.appendChild(e); };
    if (used.baseline && used.baseline.count) piece('baseline', used.baseline.first, used.baseline.last, 'Baseline');
    if (used.scoring && used.scoring.count) piece('scoring', used.scoring.first, used.scoring.last, 'Response window');
    if (used.stimulus && used.stimulus.before < used.n) { const e = document.createElement('i'); e.className = 'stimulus'; e.style.left = at(used.stimulus.at ? used.stimulus.before + 1 : Math.max(1, used.stimulus.before + 0.5)) + '%'; e.title = 'Stimulus added'; strip.appendChild(e); }
  }
  const frameTime = () => (S.info ? (editedTimes() || frameTimes())[S.frame - 1] : NaN);   // of the frame on screen, in the time base being edited
  function frameButtons() {
    const t = frameTime();
    for (const button of document.querySelectorAll('.from-frame')) { button.textContent = S.info ? 'Use frame ' + S.frame : 'Use frame'; button.title = 'Set it to the movie frame on screen' + (finite(t) ? ' (frame ' + S.frame + ', ' + fmt(+t.toFixed(2)) + ' s)' : '') + '.'; }
  }
  function collectDraft() {
    const draft = clone(S.draft), p = draft.protocol, legacy = Number($('protocol-mode').value) === 1;
    for (const key of ['sampling_hz', 'baseline_start', 'baseline_end']) draft.values[key] = readInput($(key), false);
    for (const key of ['condition', 'biological_unit']) draft.values[key] = $(key).value;
    p.schema_version = legacy ? 1 : 2; p.engine_version = legacy ? 'legacy-1' : '2.0.0'; p.signal_reference = legacy ? 'fixed_reference_fractional_dff0' : 'fixed_reference_fractional_change';
    p.name = legacy ? 'Legacy consecutive-sample protocol' : $('protocol-name').value;
    p.validation_status = legacy ? 'legacy_unvalidated' : $('protocol-validation').value;
    p.windows = Array.from(document.querySelectorAll('.window-editor')).map((section, index) => {
      const w = expandWindow(draft.protocol.windows[index], legacy), get = (k) => $('window-' + index + '-' + k);
      w.Role = index === 0 ? 'Stimulus' : 'Control'; w.Stimulus = get('Stimulus').value; w.Duration_mode = get('Duration_mode').value;
      for (const key of NUMBERS) w[key] = readInput(get(key), ['Max_gap_s', 'Expected_interval_s'].includes(key));
      if (legacy) { w.Baseline_end_s = w.Start_s; w.Stimulus_onset_s = w.Start_s; }
      w.Baseline_s = Number(w.Baseline_end_s) - Number(w.Baseline_start_s);
      if (!Number.isFinite(w.Baseline_s)) w.Baseline_s = '';
      if (legacy) return Object.fromEntries(['Role', 'Stimulus', 'Start_s', 'End_s', 'Baseline_s', 'Floor_local_fraction', 'SD_multiplier', 'Consecutive_samples', 'Max_gap_s'].map((k) => [k, w[k]]));
      Object.assign(w, { Protocol_version: 2, Protocol_name: p.name, Engine_version: p.engine_version, Validation_status: p.validation_status, Signal_reference: p.signal_reference, SD_estimator: 'sample_sd', Threshold_operator: '>=', Coverage_policy: 'complete_windows', Response_direction: 'increase', Smoothing: 'none' });
      return w;
    });
    // a setting of the analysis that the engine never sees: kept only while ticked, as the session keeps it
    if ($('automatic-marks').checked) draft.values.automatic_marks = true; else delete draft.values.automatic_marks;
    const overrides = $('field-overrides').value.trim();
    if (legacy) p.field_overrides = null;
    else { try { p.field_overrides = overrides ? JSON.parse(overrides) : null; } catch (_) { p.field_overrides = overrides; } }
    return draft;
  }
  function fillDraft(draft) {
    S.draft = clone(draft);
    exactInput($('sampling_hz'), draft.values.sampling_hz);   // the F₀ frame boxes are made with the primary stimulus's rows
    for (const key of ['condition', 'biological_unit']) $(key).value = draft.values[key] || '';
    $('automatic-marks').checked = draft.values.automatic_marks === true;
    $('protocol-mode').value = draft.protocol.schema_version;
    $('protocol-name').value = draft.protocol.name; $('protocol-validation').value = draft.protocol.validation_status === 'legacy_unvalidated' ? 'user_defined' : draft.protocol.validation_status;
    $('field-overrides').value = draft.protocol.field_overrides ? JSON.stringify(draft.protocol.field_overrides, null, 2) : '';
    for (const key of FILES) fileNote(key);
    renderWindows();
  }
  function fileNote(key) { const file = S.draft[key]; $(key + '-name').textContent = file ? file.name : 'No file selected'; $(key + '-clear').hidden = !file; }
  // Draft validation reports errors without changing settings, saved calls or the applied movie clock.
  function validateDraft(draft) {
    // Every message sits under the field it is about, with its repair. Only a message with no field, or whose field is folded away
    // (Advanced), is listed in the Run box, which otherwise just says how many there are.
    const errors = [], warnings = [], notes = [], fields = new Set();
    const add = (field, message) => { errors.push(message); if (field) fields.add(field); notes.push({ beside: field, message, level: 'error' }); };
    // Two time fields that disagree: both are marked, and one message sits beside the field to change, with the repair.
    const conflict = (mark, beside, message, repair) => { errors.push(message); mark.forEach((f) => fields.add(f)); notes.push({ beside, message, repair, level: 'error' }); };
    const warn = (message, repair, beside) => { warnings.push(message); notes.push({ beside, message, repair, level: 'warning' }); };
    const times = draftTimes(draft.timestamps, draft.values.sampling_hz);
    const numeric = (v) => v !== null && v !== '' && typeof v !== 'boolean' && Number.isFinite(Number(v));
    const value = draft.values, n = S.info.frames;
    if (!draft.timestamps && (!numeric(value.sampling_hz) || Number(value.sampling_hz) <= 0)) add('sampling_hz', 'Frame rate must be a positive number, or supply valid timestamps.');
    if (!numeric(value.baseline_start) || !Number.isInteger(Number(value.baseline_start)) || Number(value.baseline_start) < 1) add('baseline_start', 'F₀: the first frame must be a whole number from 1.');
    if (!numeric(value.baseline_end) || !Number.isInteger(Number(value.baseline_end)) || Number(value.baseline_end) < Number(value.baseline_start) || Number(value.baseline_end) >= n) add('baseline_end', 'F₀: the last frame must follow the first and be before frame ' + n + '.');
    const legacy = draft.protocol.schema_version === 1, names = new Set();
    draft.protocol.windows.forEach((original, index) => {
      const row = expandWindow(original, legacy), prefix = 'window-' + index + '-', label = index ? 'Control: ' : 'Stimulus: ';
      const required = NUMBERS.filter((key) => !['Max_gap_s', 'Expected_interval_s'].includes(key));
      const box = (key) => $(prefix + key), v = (key) => Number(row[key]), has = (key) => numeric(row[key]), fr = (id) => frameBox(id, index);
      for (const key of required) if (!has(key) && !TIME_KEYS.includes(key)) add(prefix + key, label + FIELD_TITLES[key] + ' needs a number.');
      const name = String(row.Stimulus || '').trim().toLowerCase();
      if (!name || names.has(name)) add(prefix + 'Stimulus', label + 'use a nonempty, distinct name.'); names.add(name);
      const used = windowFrames(row, times), last = used ? used.times[used.n - 1] : NaN, N = wholeIn(fr('stim-frame')), okN = N >= 1 && N < n;
      const relink = (cellId) => () => { S.unlinked.delete(cellId === 'f0-cell' ? 'f0' : cellId + '-' + index); (cellId === 'f0-cell' ? $('f0-cell') : fr(cellId)).dataset.linked = '1'; applyFrames(index); draftChanged(); };
      const setFrames = (changes) => () => { for (const [id, frame] of changes) fr(id).value = frame; applyFrames(index); draftChanged(); };
      const ids = (...names) => names.map((id) => id + '-' + index);
      if (!times || !fr('stim-frame')) { /* no frame times: the frame rate message says so */ }
      else if (!okN) conflict(ids('stim-frame'), 'stim-frame-' + index, label + (String(fr('stim-frame').value).trim() === '' ? 'enter the frame after which it was added.' : 'the stimulus must come after a frame from 1 to ' + (n - 1) + '.'));
      else {
        const a = wholeIn(fr('base-from')), b = wholeIn(fr('base-to')), c = wholeIn(fr('resp-from')), d = String(fr('resp-to').value).trim() === '' ? n : wholeIn(fr('resp-to')), rule = frameRule(N);
        if (!(a >= 1 && b >= a)) conflict(ids('base-from', 'base-to'), 'base-from-' + index, label + 'the baseline needs a first and a last frame.', { text: 'Use ' + rule[0] + '–' + rule[1], run: relink('base-cell') });
        else if (b > N) conflict(ids('base-to'), 'base-to-' + index, label + (Math.max(a, N + 1) === b ? 'baseline frame ' + b + ' comes' : 'baseline frames ' + Math.max(a, N + 1) + '–' + b + ' come') + ' after the stimulus.', a > N ? { text: 'Use ' + rule[0] + '–' + rule[1], run: relink('base-cell') } : { text: 'Use ' + a + '–' + N, run: setFrames([['base-to', N]]) });
        else if (has('Min_baseline_samples') && used && used.baseline && used.baseline.count < v('Min_baseline_samples')) {
          // for the primary stimulus this leaves every cell without a call, so the run is stopped; for the control it only loses the control
          const words = label + 'the baseline holds ' + plural(used.baseline.count, 'frame') + ' but at least ' + v('Min_baseline_samples') + ' are needed. No cell would get a call.', fix = N >= v('Min_baseline_samples') ? { text: 'Use ' + rule[0] + '–' + rule[1], run: relink('base-cell') } : null;
          if (index) warn(words, fix, 'base-from-' + index); else conflict(ids('base-from', 'base-to'), 'base-from-' + index, words, fix);
        }
        if (!(c >= 1 && d >= c && d <= n)) conflict(ids('resp-from', 'resp-to'), 'resp-from-' + index, label + 'the response window needs a first and a last frame, up to ' + n + '.', { text: 'Use ' + (N + 1) + '–' + n, run: relink('resp-cell') });
        else if (c <= N) conflict(ids('resp-from'), 'resp-from-' + index, label + 'the response window starts at frame ' + c + ', before the stimulus.', { text: 'Start at ' + (N + 1), run: setFrames([['resp-from', N + 1]]) });
        else if (has('End_s') && v('End_s') > last + 1e-9) {   // a loaded setting that runs past the recording
          if (legacy) warn(label + 'the response window ends after the last frame. Cells may get no call.', { text: 'End at frame ' + n, run: setFrames([['resp-to', n]]) }, 'resp-to-' + index);
          else conflict(ids('resp-to'), 'resp-to-' + index, label + 'the response window ends after the last frame (' + seconds(+last.toFixed(2)) + '). No cell would get a call.', { text: 'End at frame ' + n, run: setFrames([['resp-to', n]]) });
        }
        // A field that sits a few pixels from where the cells were chosen: said once, as the plugin words it (no repair: nothing is wrong with the frames chosen)
        if (!index && S.info.motion && S.info.motionNote && !S.info.motion.moved && S.info.motion.slight > 0) warn(S.info.motionNote + ' The outlines are slightly off the cells there.', null, 'fiji-frames');
        if (!index && S.info.motion && S.info.motion.moved && S.info.motion.smeared) warn('The cells were chosen on frames ' + S.info.motion.smeared[0] + '–' + S.info.motion.smeared[1] + ', which include the movement, so their picture is smeared.', null, 'fiji-frames');
        if (!index && movedRuns().length) {
          // The field moved: in those frames the outlines are not on the same cells. Baselines should leave them out.
          const f0 = [Number(value.baseline_start), Number(value.baseline_end)], size = 'The field moved (about ' + S.info.motion.largest + ' px) during ' + movedText() + '. ';
          const inF0 = Number.isInteger(f0[0]) && Number.isInteger(f0[1]) ? movedWithin(f0[0], f0[1]) : { count: 0 }, inBase = movedWithin(a, b), inResp = movedWithin(c, d);
          const sameAsF0 = a === f0[0] && b === f0[1];
          if (inF0.count) warn(size + 'F₀ frames ' + f0[0] + '–' + f0[1] + ' include ' + plural(inF0.count, 'moved frame') + '.', inF0.after <= f0[1] ? { text: 'Use ' + inF0.after + '–' + f0[1], run: () => { S.unlinked.add('f0'); $('f0-cell').dataset.linked = ''; exactInput($('baseline_start'), inF0.after); applyFrames(index); draftChanged(); } } : null, 'f0-cell');
          if (inBase.count && !(sameAsF0 && inF0.count && fr('base-cell').dataset.linked === '1')) warn(size + 'Baseline frames ' + a + '–' + b + ' include ' + plural(inBase.count, 'moved frame') + '.', inBase.after <= b ? { text: 'Use ' + inBase.after + '–' + b, run: () => { S.unlinked.add('base-cell-' + index); fr('base-cell').dataset.linked = ''; fr('base-from').value = inBase.after; applyFrames(index); draftChanged(); } } : null, 'base-from-' + index);
          if (inResp.count) warn(size + 'The response window includes ' + plural(inResp.count, 'moved frame') + ': values there are not of the same cells.', null, 'resp-from-' + index);
        }
        if (!index) {
          // F₀ frames after the stimulus, and cells that were chosen on frames after it
          const f0 = [Number(value.baseline_start), Number(value.baseline_end)], chosen = S.info.fijiFrames;
          if (Number.isInteger(f0[0]) && Number.isInteger(f0[1]) && f0[1] >= f0[0] && f0[1] > N)
            warn((Math.max(f0[0], N + 1) === f0[1] ? 'F₀ frame ' + f0[1] + ' comes' : 'F₀ frames ' + Math.max(f0[0], N + 1) + '–' + f0[1] + ' come') + ' after the stimulus.', { text: 'Use ' + rule[0] + '–' + rule[1], run: relink('f0-cell') }, 'f0-cell');

        }
      }
      for (const key of ['Floor_local_fraction', 'SD_multiplier', 'Noise_floor_local_fraction']) if (numeric(row[key]) && Number(row[key]) < 0) add(prefix + key, label + 'thresholds cannot be negative.');
      if (numeric(row.Floor_local_fraction) && numeric(row.SD_multiplier) && Number(row.Floor_local_fraction) === 0 && Number(row.SD_multiplier) === 0) add(prefix + 'Floor_local_fraction', label + 'specify an amplitude or noise threshold.');
      if (legacy && numeric(row.Floor_local_fraction) && Number(row.Floor_local_fraction) <= 0) add(prefix + 'Floor_local_fraction', label + 'legacy mode requires a positive minimum rise.');
      if (numeric(row.Min_duration_s) && Number(row.Min_duration_s) <= 0) add(prefix + 'Min_duration_s', label + 'duration must be positive.');
      for (const [key, minimum] of [['Consecutive_samples', 1], ['Min_baseline_samples', 2]]) if (numeric(row[key]) && (!Number.isInteger(Number(row[key])) || Number(row[key]) < minimum)) add(prefix + key, label + FIELD_TITLES[key] + ' must be a whole number ≥ ' + minimum + '.');
      for (const key of ['Max_gap_s', 'Expected_interval_s']) if (row[key] !== null && row[key] !== '' && (!numeric(row[key]) || Number(row[key]) <= 0)) add(prefix + key, label + FIELD_TITLES[key] + ' must be positive or blank.');
      if (!draft.timestamps && numeric(value.sampling_hz) && Number(value.sampling_hz) > 0) {
        const interval = 1 / Number(value.sampling_hz);
        if (!legacy && row.Duration_mode === 'seconds' && numeric(row.Min_duration_s) && Number(row.Min_duration_s) < interval) warn(label + 'minimum duration (' + row.Min_duration_s + ' s) is shorter than the frame interval (' + fmt(interval) + ' s). Timing is coarsely resolved.', null, prefix + 'Min_duration_s');
      }
    });
    if (!String(draft.protocol.name || '').trim()) add('protocol-name', 'Enter a protocol name.');
    if (draft.protocol.field_overrides !== null && !Array.isArray(draft.protocol.field_overrides)) add('field-overrides', 'Field overrides must be a valid JSON array.');
    const signature = JSON.stringify(draft);
    if (signature === S.draftAck && S.state && S.state.validation) {
      // One message per problem: the plugin's general checks are shown only for what this page did not already explain.
      if (!errors.length) for (const message of S.state.validation.errors) { errors.push(message); notes.push({ beside: null, message, level: 'error' }); }
      const mine = [/F₀ baseline includes/, /response end exceeds/, /holds \d+ samples?/, /^Check acquisition timing/];
      for (const message of S.state.validation.warnings) if (!times || !mine.some((p) => p.test(message))) { warnings.push(message); notes.push({ beside: null, message, level: 'warning' }); }
    }
    for (const box of document.querySelectorAll('#side [aria-invalid]')) box.removeAttribute('aria-invalid');
    for (const id of fields) if ($(id)) $(id).setAttribute('aria-invalid', 'true');
    // Rebuilt only when what is shown changes, so a button is never replaced under the mouse by a routine refresh.
    const onPage = (id) => { const field = id && $(id); return !!field && (field.closest('.row') || field).offsetParent !== null; };   // its row, since a row that follows the stimulus hides its own boxes
    for (const note of notes) note.inline = onPage(note.beside);
    const shown = JSON.stringify(notes.map((e) => [e.beside, e.message, e.level, e.inline, e.repair && e.repair.text]));
    if (shown !== S.validationShown) {
      S.validationShown = shown;
      const action = (repair) => { const button = document.createElement('button'); button.type = 'button'; button.className = 'ghost repair'; button.textContent = repair.text; button.addEventListener('click', repair.run); return button; };
      for (const old of document.querySelectorAll('#side .field-error, #side .field-warning')) old.remove();
      const said = new Set();
      for (const entry of notes) {
        if (!entry.inline || said.has(entry.beside + '|' + entry.message)) continue; said.add(entry.beside + '|' + entry.message);
        const field = $(entry.beside), after = field.closest('.row') || field;
        const note = document.createElement('div'); note.className = entry.level === 'error' ? 'field-error' : 'field-warning'; note.setAttribute('role', entry.level === 'error' ? 'alert' : 'status');
        const text = document.createElement('span'); text.textContent = entry.message; note.appendChild(text);
        if (entry.repair) note.appendChild(action(entry.repair));
        let last = after; while (last.nextElementSibling && /^field-(error|warning)$/.test(last.nextElementSibling.className)) last = last.nextElementSibling;
        last.after(note);
      }
      const box = $('validation'), loose = notes.filter((e) => !e.inline), above = { error: notes.filter((e) => e.inline && e.level === 'error').length, warning: notes.filter((e) => e.inline && e.level === 'warning').length };
      box.hidden = !notes.length; box.replaceChildren();
      if (above.error || above.warning) {   // the messages themselves are beside their fields; this takes you to the first
        const jump = document.createElement('button'); jump.type = 'button'; jump.id = 'validation-jump'; jump.className = 'ghost jump' + (above.error ? '' : ' warning');
        jump.textContent = above.error ? (above.error === 1 ? 'Fix 1 problem above' : 'Fix ' + above.error + ' problems above') : plural(above.warning, 'warning') + ' above';
        jump.addEventListener('click', () => { const first = document.querySelector(above.error ? '#side .field-error' : '#side .field-warning'); if (first) first.scrollIntoView({ block: 'center', behavior: 'smooth' }); });
        box.appendChild(jump);
      }
      for (const [level, className] of [['error', ''], ['warning', 'warning']]) {
        const messages = loose.filter((e) => e.level === level); if (!messages.length) continue;
        const list = document.createElement('ul'); list.className = className;
        const once = new Map(); for (const entry of messages) once.set(String(entry.message).split('. ')[0], entry);
        for (const entry of once.values()) { const li = document.createElement('li'); li.textContent = entry.message; if (entry.repair) { li.append(' '); li.appendChild(action(entry.repair)); } list.appendChild(li); }
        box.appendChild(list);
      }
    }
    return { errors, warnings };
  }
  /** Cells measured since 2026.11.20.5 are named by a long stable ID. Show them as "Cell N" (their measured order), with the
   *  ID in tooltips; the tables and exports keep the IDs. Short names, such as Cell001, are shown as they are. */
  function shortName(name) { if (!name) return ''; const i = S.info.names.indexOf(name); return name.length > 14 && i >= 0 ? 'Cell ' + (i + 1) : name; }
  const fullName = (name) => shortName(name) === name ? name : shortName(name) + ' (' + name + ')';
  function qualityFor(name) {
    if (!S.results || !name) return [];
    const q = S.quality[name] || {}, row = callRow(name), flags = new Set();
    for (const text of [q.QC_flags, row.QC_flags]) for (const flag of String(text || '').split(/[;|,]/)) if (flag.trim() && !RECORDING_FLAGS.has(flag.trim())) flags.add(flag.trim());
    if (belowZero(name)) flags.add('TRACE_BELOW_MINUS_ONE');
    for (const flag of checkFlags(name)) flags.add(CHECK + codeOf(flag));   // a check's text has commas, so it is never split like the engine's flags
    return Array.from(flags);
  }
  // ---------------------------------------------------------------- the checks on the calls (results.checks), made in Fiji after the engine
  /** A check flag is keyed by its code (POSSIBLE_PASSING_OBJECT, for every rule that finds one); its own text ("possible passing
   *  object, frames 201 to 210") is the cell's. The labels come with the checks; these are the same words for an analysis without them. */
  const CHECK = 'check:', isCheck = (flag) => flag.startsWith(CHECK);
  const CHECK_LABELS = { POSSIBLE_PASSING_OBJECT: 'Possible passing object', RISE_ON_EDGE: 'Rise on the edge', FIELD_MOVED: 'Field moved', SHORT_SEPARATE_CROSSINGS: 'Short separate crossings' };
  const codeOf = (flag) => flag.code || flag.rule;
  const checkLabel = (code) => ((S.results && S.results.checks && S.results.checks.labels) || CHECK_LABELS)[code];
  /** The checks of one cell, or null (an analysis made before the checks has none). */
  const checksOf = (name) => { const c = S.results && S.results.checks; return (c && c.cells && name && c.cells[name]) || null; };
  /** The checks of one cell as the page shows them: they are about the calls for the primary stimulus, so a control window shows
   *  none, as in the Fiji window. Every flag on screen, in the flagged list, its filters and counts and the cell panel, comes from here. */
  const shownChecks = (name) => primaryShown() ? checksOf(name) : null;
  /** The flags in a cell's checks, one by one: { text, rule }. */
  const flagsIn = (c) => c && Array.isArray(c.flags) ? c.flags.filter((f) => f && f.text) : [];
  const checkFlags = (name) => flagsIn(shownChecks(name));
  /** The reason the mark form starts with on a responder, from its first check flag: none for a cell that did not respond. */
  const REASON_OF_CODE = { POSSIBLE_PASSING_OBJECT: 'Passing object', RISE_ON_EDGE: 'Passing object', FIELD_MOVED: 'Movement', SHORT_SEPARATE_CROSSINGS: 'Noise' };
  function startingReason(name, automatic) {
    if (automatic !== 'Responded') return '';
    const reasons = (S.results.review || {}).reasons || [];
    for (const flag of flagsIn(checksOf(name))) { const reason = REASON_OF_CODE[codeOf(flag)]; if (reason && reasons.includes(reason)) return reason; }
    return '';
  }
  /** A cell whose ΔF/F₀ falls below −1: its corrected fluorescence goes below zero, the mark of an F₀ near zero. Flagged
   *  "Below zero"; such a cell is still drawn, but no automatic axis is taken from it (plotCells). */
  const belowZero = (name) => (S.results.traces[name] || []).some((v) => finite(v) && v < -1);
  // Set on every cell by the engine whenever the window has the matching Window_warning: a property of the protocol and frame
  // rate, shown once as a notice for the recording instead of flagging every cell.
  const RECORDING_FLAGS = new Set(['DURATION_BELOW_SAMPLE_INTERVAL']);
  const FLAG_TEXT = {
    BASELINE_SD_AT_LEAST_ONE_THIRD_LOCAL_F0: 'Baseline noise is at least a third of the local F₀, so the fractional rise threshold is dominated by noise',
    DURATION_BELOW_SAMPLE_INTERVAL: 'Minimum response duration is shorter than one frame interval; timing is coarsely resolved',
    reference_below_3_baseline_sd_review: 'Low F₀ relative to baseline noise; review fractional amplitudes and the response call',
    TRACE_BELOW_MINUS_ONE: 'ΔF/F₀ falls below −1; corrected fluorescence becomes negative relative to its positive reference. Review background and baseline stability'
  };
  const plural = (n, word) => n + ' ' + word + (n === 1 ? '' : 's');
  /** How many cells of the selected window got a call (responded or did not). Zero means the run finished without a single call. */
  const called = () => { const c = S.results ? ((S.results.countsByWindow || {})[S.window] || S.results.counts || {}) : {}; const v = number(c.Evaluable); return isFinite(v) ? v : 0; };
  /** Why cells are unscorable, in the engine's words: the reason most of them share. */
  function unscorableReason() {
    const seen = new Map(); for (const name of S.info.names) { const row = callRow(name); if (row.Status === 'Unscorable' && row.Reason) seen.set(row.Reason, (seen.get(row.Reason) || 0) + 1); }
    let best = ''; for (const [reason, n] of seen) if (!best || n > seen.get(best)) best = reason; return best;
  }
  const flagText = (flag) => isCheck(flag) ? flagShort(flag) : FLAG_TEXT[flag] || flag.replace(/_/g, ' ');
  /** Each flag in two or three words, for the list of flagged cells and the tables; FLAG_TEXT has the full sentence. */
  const FLAG_SHORT = { BASELINE_SD_AT_LEAST_ONE_THIRD_LOCAL_F0: 'Noisy baseline', reference_below_3_baseline_sd_review: 'Low F₀', TRACE_BELOW_MINUS_ONE: 'Below zero', DURATION_BELOW_SAMPLE_INTERVAL: 'Duration under a frame' };
  const flagShort = (flag) => isCheck(flag) ? checkLabel(flag.slice(CHECK.length)) || flag.slice(CHECK.length) : FLAG_SHORT[flag] || checkLabel(flag) || flagText(flag);
  /** What a cell is flagged for, in the list of flagged cells: the engine's flags in their short words, then each check's own text. */
  const flaggedFor = (name) => qualityFor(name).filter((flag) => !isCheck(flag)).map(flagShort).concat(checkFlags(name).map((f) => f.text));
  function showQuality() {
    if (!S.results) return;
    const summary = $('quality-summary'), list = $('quality-list'), counts = new Map(); let affected = 0, responders = 0;
    for (const name of S.info.names) { const flags = qualityFor(name); if (flags.length) { affected++; if (callOf(name) === 'Responded') responders++; }
      for (const flag of flags) { const item = counts.get(flag) || { total: 0, responders: 0 }; item.total++; if (callOf(name) === 'Responded') item.responders++; counts.set(flag, item); }
    }
    const notices = [], window = selectedCall(), baseline = S.results.baselineFrames;
    if (baseline && S.results.time[baseline[1] - 1] > number(window.Stimulus_onset_s || window.Start_s) + 1e-9) notices.push('F₀ frames ' + baseline[0] + '–' + baseline[1] + ' of these results include frames after the stimulus.');
    if (window.Window_warning) notices.push(window.Window_warning);
    if (!window.Window_warning && S.info.names.some((name) => { const q = S.quality[name] || {}, row = callRow(name); return [q.QC_flags, row.QC_flags].some((t) => String(t || '').includes('DURATION_BELOW_SAMPLE_INTERVAL')); }))
      notices.push('Applies to the whole recording: ' + FLAG_TEXT.DURATION_BELOW_SAMPLE_INTERVAL + '.');
    // the note of the checks, once for the analysis, whatever window is shown: why the cells were not checked for a passing object,
    // or that the field's movement was not saved (a cell's own reason is in its panel)
    const checks = S.results.checks, note = (checks && checks.note) || ''; if (note) notices.push(note);
    if (S.qualityFlag && !counts.has(S.qualityFlag)) S.qualityFlag = null;
    $('quality-review').hidden = !affected && !notices.length;
    // when that note is all there is to say, the chip names it and the note leads the list
    const onlyNote = !affected && notices.length === 1 && !!note;
    summary.textContent = affected ? plural(affected, 'cell') + ' flagged' : onlyNote ? (checks.ringMeasured ? 'Checks note' : 'Passing objects not checked') : 'Timing note';
    $('quality-lead').textContent = affected ? plural(affected, 'cell') + ' flagged' + (called() ? ', including ' + plural(responders, 'responder') : '') + '.' : onlyNote ? note : 'Review the applied timing.';
    // Every flagged cell, with its call and what it is flagged for. A click opens the cell; the arrows then step through
    // the flagged cells. With one flag picked below, only its cells are listed.
    const flagged = S.info.names.filter((name) => { const flags = qualityFor(name); return flags.length && (!S.qualityFlag || flags.includes(S.qualityFlag)); });
    const cells = $('table-flagged'); cells.replaceChildren(); $('quality-cells').hidden = !flagged.length;
    if (flagged.length) {
      // the call of the window shown: after review on the primary stimulus while marks apply, named as in the metric tables
      const head = cells.createTHead().insertRow(); for (const title of ['Cell', primaryShown() ? CALL_COLUMN.name : 'Call', 'Flagged for']) { const th = document.createElement('th'); th.textContent = title; head.appendChild(th); }
      const body = cells.createTBody();
      for (const name of flagged) { const tr = body.insertRow(), call = callOf(name) || '', mark = markOf(name); tr.dataset.cell = name; tr.classList.toggle('on', name === S.cell); tr.title = 'Open ' + shortName(name);
        tr.insertCell().textContent = shortName(name); const said = tr.insertCell(); said.innerHTML = '<span class="dot" style="background:' + CALL[kind(call)] + '"></span>' + escape(call); if (mark) said.title = 'Automatic call: ' + mark.automatic;
        tr.insertCell().textContent = flaggedFor(name).join(' · '); }
    }
    list.replaceChildren();
    for (const [flag, count] of counts) {
      // each flag is a filter of its own, so a small group (for example the cells whose ΔF/F₀ falls below −1) can be isolated
      const li = document.createElement('li'), link = document.createElement('button'); link.type = 'button'; link.className = 'flag-filter' + (S.qualityFlag === flag ? ' on' : '');
      link.textContent = flagShort(flag) + ': ' + plural(count.total, 'cell') + (called() ? ', ' + plural(count.responders, 'responder') : ''); link.title = S.qualityFlag === flag ? 'Show all flagged cells again' : 'Show only these cells';
      link.addEventListener('click', () => { S.qualityFlag = S.qualityFlag === flag ? null : flag; S.group = 5; $('group').value = 5; S.keep = null; showQuality(); listCells(); showTab('review'); });
      li.appendChild(link);
      if (flagShort(flag) !== flagText(flag)) { const what = document.createElement('span'); what.className = 'flag-what'; what.textContent = flagText(flag) + '.'; li.appendChild(what); }
      list.appendChild(li);
    }
    // the notices; a note that already leads the list is not said again under it
    for (const message of new Set(notices)) { if (onlyNote && message === note) continue; const li = document.createElement('li'); li.textContent = message; list.appendChild(li); }
    list.hidden = !list.children.length;
  }
  // ---------------------------------------------------------------- the saved columns in words: one dictionary for every table on screen
  const DFF = 'ΔF/F₀', SEC = 's', PCT = '%', GRAY = 'gray levels', AREA = 'ΔF/F₀ × s', RATE = 'ΔF/F₀ per s', FRAMES = 'frames', CELLS = 'cells', PIXELS = 'pixels';
  /** Decimals on screen by unit, the same down a whole column: time 2, ΔF/F₀ 3, percentages 1, gray levels 2. */
  const DECIMALS = { [SEC]: 2, [DFF]: 3, [PCT]: 1, [GRAY]: 2, [AREA]: 2, [RATE]: 3, [FRAMES]: 0, [CELLS]: 0, [PIXELS]: 0 };
  /** Units written after a column's name ("Peak, ΔF/F₀"). A percentage carries its "%" in every cell; a count needs none. */
  const NAMED_UNITS = new Set([SEC, DFF, GRAY, AREA, RATE, PIXELS]);
  /** What each saved column is: [saved name, name on screen, unit, where it is measured, more]. more: alias (the saved column
   *  it equals by definition in this plugin: the two are shown once), response (the timing or shape of a response: shown for a
   *  cell that responded only, and summarised over those cells), fraction (saved as a fraction, shown in %), decimals, about.
   *  Two tables can use one saved name for two things (Status, Baseline_SD), so the words are looked up by table: "cell" is
   *  metrics.csv and cell_quality.csv, "window" the calls, counts and window metrics; "any" holds the rest. The saved files keep
   *  their own names, values and every column; a column missing here is shown under its saved name in words (columnWords). */
  const COLUMN_WORDS = {
    any: [
      ['Cell', 'Cell'], ['Field', 'Recording'], ['Role', 'Role'], ['Stimulus', 'Stimulus'], ['Population', 'Population'], ['Reason', 'Reason'], ['QC_flags', 'Flags', '', '', { flags: true }],
      ['Condition', 'Condition'], ['Biological_unit', 'Biological unit'], ['Group', 'Group'], ['Replicate', 'Replicate'], ['Selection_basis', 'Chosen'], ['Note', 'Note'],
      ['Raw_QC_available', 'Raw checks'], ['Signal_kind', 'Signal'], ['Input_mode', 'Input'], ['Normalized_units', 'Normalized units'], ['Normalization_reference', 'Reference'],
      ['Responded', 'Responded', CELLS], ['Evaluable', 'Evaluable', CELLS], ['Percent_responding', '% responding', PCT],
      ['ROI_ID', 'Outline ID'], ['Accepted', 'Accepted'], ['Included_in_measurement', 'Measured'], ['Exclusion_reason', 'Why left out'], ['ROI_sha256', 'Outline hash'],
      ['Area_pixels', 'Area', PIXELS], ['Background_pixels', 'Background area', PIXELS], ['Near_edge', 'Near the edge'],
      ['Accepted_ROIs', 'Accepted outlines', CELLS], ['Premeasurement_excluded', 'Left out before measuring', CELLS], ['Registered_for_analysis', 'Registered', CELLS],
      ['Normalization_valid', 'Normalization valid', CELLS], ['Normalization_excluded', 'Normalization excluded', CELLS], ['Fields', 'Recordings', '', '', { decimals: 0 }],
      ['Replicates_evaluable', 'Replicates', '', '', { decimals: 0 }], ['Mean_replicate_percent', 'Mean of replicates', PCT], ['SD_replicate_percent', 'SD of replicates', PCT], ['Cell_pooled_percent', 'Cells pooled', PCT],
    ],
    cell: [
      ['Status', 'Normalization', '', '', { about: 'Whether F₀ can normalize this cell; the response call is under Call' }], ['Reason', 'Normalization reason'],
      ['F0', 'F₀', GRAY, 'F₀ frames'], ['Reference_value', 'Reference value', GRAY, 'F₀ frames', { alias: 'F0' }],
      ['Baseline_SD_raw', 'Raw baseline SD', GRAY, 'F₀ frames'], ['Baseline_SD_input', 'Input baseline SD', GRAY, 'F₀ frames', { alias: 'Baseline_SD_raw' }],
      ['Baseline_samples', 'F₀ frames', FRAMES], ['Normalization_minimum_reference', 'Smallest F₀ allowed', GRAY], ['Metric_status', 'Metrics'],
      ['Peak_dFF0', 'Peak', DFF, 'whole recording'], ['Response_Amplitude', 'Peak above baseline', DFF, 'whole recording', { alias: 'Peak_dFF0' }],
      ['Time_to_Peak', 'Time to peak', SEC, 'whole recording, from its start', { response: true }],
      ['Time_to_25_Peak', 'Time to 25% of peak', SEC, 'whole recording, from its start', { response: true }],
      ['Time_to_50_Peak', 'Time to 50% of peak', SEC, 'whole recording, from its start', { response: true }],
      ['Time_to_75_Peak', 'Time to 75% of peak', SEC, 'whole recording, from its start', { response: true }],
      ['Rise_Time', 'Rise time', SEC, 'whole recording, 10% to 90% of peak', { response: true }], ['Rise_Rate_10_90_dFF0_per_s', 'Rise rate', RATE, 'whole recording, 10% to 90% of peak', { response: true }],
      ['AUC', 'AUC', AREA, 'whole recording, signed'],
      ['FWHM', 'FWHM', SEC, 'whole recording, width at half of peak', { response: true }], ['FWHM_Censored', 'FWHM cut off', '', 'still above half at the end'],
      ['FWHM_Lower_Bound', 'FWHM, lower bound', SEC, 'whole recording, still above half at the end', { response: true }], ['Half_Width', 'Half-width', SEC, 'whole recording, FWHM ÷ 2', { response: true }],
      ['Baseline_SD', 'Baseline SD', DFF, 'F₀ frames'], ['SNR', 'Signal-to-noise', '', 'whole recording, peak ÷ baseline SD', { decimals: 2 }],
    ],
    window: [
      ['Status', 'Call'],
      ['Baseline_start_s', 'Baseline start', SEC], ['Baseline_end_s', 'Baseline end', SEC], ['Stimulus_onset_s', 'Stimulus at', SEC], ['Start_s', 'Window start', SEC], ['End_s', 'Window end', SEC],
      ['Floor_local_fraction', 'Minimum rise', PCT, '', { fraction: true }], ['SD_multiplier', 'SD multiple', '', '', { decimals: 1 }], ['Duration_mode', 'Duration in'], ['Min_duration_s', 'Shortest response', SEC],
      ['Consecutive_samples', 'Frames in a row', FRAMES], ['Max_gap_s', 'Largest gap', SEC], ['Min_baseline_samples', 'Fewest baseline frames', FRAMES], ['Noise_floor_local_fraction', 'Noise floor', PCT, '', { fraction: true }],
      ['Expected_interval_s', 'Expected frame interval', SEC, '', { decimals: 6 }], ['Baseline_s', 'Baseline length', SEC], ['Protocol_version', 'Protocol version', '', '', { decimals: 0 }], ['Protocol_name', 'Protocol'],
      ['Engine_version', 'Engine version'], ['Validation_status', 'Protocol status'], ['Signal_reference', 'Signal reference'], ['SD_estimator', 'SD estimator'], ['Threshold_operator', 'Threshold test'],
      ['Coverage_policy', 'Coverage rule'], ['Response_direction', 'Direction'], ['Smoothing', 'Smoothing'], ['Protocol_hash', 'Protocol hash'], ['Applied_protocol_hash', 'Applied protocol hash'],
      ['Resolved_interval_s', 'Frame interval', SEC, '', { decimals: 6 }], ['Gap_mode', 'Gap rule'], ['Nominal_min_samples', 'Frames required', FRAMES], ['Response_reference', 'Response reference'],
      ['Window_clipped', 'Window cut short'], ['Window_warning', 'Window warning'], ['Effective_end_s', 'Window end used', SEC], ['Field_reason', 'Recording note'],
      ['Requested_response_span_s', 'Window length set', SEC], ['Response_observed_span_s', 'Window length recorded', SEC], ['Baseline_observed_span_s', 'Baseline length recorded', SEC],
      ['Baseline_mean', 'Baseline mean', PCT, 'response baseline, of the local level', { fraction: true }],
      ['Baseline_SD', 'Baseline SD', PCT, 'response baseline, of the local level', { fraction: true, about: 'Not the SD of the F₀ frames in the metrics' }],
      ['Baseline_SD_global', 'Global baseline SD', PCT, '', { fraction: true }], ['Baseline_samples', 'Baseline frames', FRAMES], ['Response_samples', 'Window frames', FRAMES],
      ['Local_F0_relative', 'Local level ÷ F₀', '', '', { decimals: 3 }], ['Local_reference_relative', 'Local reference ÷ F₀', '', '', { decimals: 3, alias: 'Local_F0_relative' }],
      ['Input_reference_value', 'Input reference', '', '', { decimals: 3 }], ['Local_reference_value', 'Local reference', '', '', { decimals: 3 }],
      ['Threshold_increase', 'Required rise', PCT, 'of the local level', { fraction: true }], ['Peak_increase', 'Largest rise', PCT, 'response window, of the local level', { fraction: true }],
      ['Floor_threshold', 'Minimum rise threshold', PCT, '', { fraction: true, alias: 'Floor_local_fraction' }], ['Noise_threshold', 'Noise threshold', PCT, '', { fraction: true }], ['Noise_SD_used', 'Noise SD used', PCT, '', { fraction: true }],
      ['Threshold_source', 'Required rise set by', '', '', { values: { Floor: 'Minimum rise' }, noise: true }], ['Margin', 'Largest ÷ required', '', '', { decimals: 2 }],
      ['Onset_s', 'Onset', SEC, 'from recording start', { response: true, about: 'The first frame of the held rise' }], ['Run_start_s', 'Run start', SEC, 'from recording start', { response: true, alias: 'Onset_s' }],
      ['Confirmation_s', 'Response confirmed', SEC, 'from recording start', { response: true }], ['Run_end_s', 'Response end', SEC, 'from recording start', { response: true }],
      ['Run_duration_s', 'Response duration', SEC, 'above the required rise', { response: true }], ['Run_samples', 'Frames held', FRAMES, 'above the required rise', { response: true }],
      ['Response_latency_s', 'Onset latency', SEC, 'from the stimulus', { response: true }],
      ['Peak_s', 'Largest rise at', SEC, 'from recording start', { response: true, about: 'The highest frame of the response window; it comes before the onset when that frame alone was not held' }],
      ['Peak_latency_s', 'Largest rise latency', SEC, 'from the stimulus', { response: true }], ['Window_peak_local_fraction', 'Window peak', PCT, 'response window', { fraction: true, alias: 'Peak_increase' }],
      ['Window_AUC_fraction_s', 'Window AUC', AREA, 'response window, signed'], ['Window_AUC_start_s', 'Window AUC start', SEC], ['Window_AUC_end_s', 'Window AUC end', SEC],
      ['Window_AUC_covered_s', 'Window AUC covered', SEC], ['Window_AUC_coverage_fraction', 'Window AUC coverage', PCT, '', { fraction: true }], ['AUC_method', 'AUC method'],
      ['Eligible', 'Eligible'], ['Measurements_available', 'Measured'], ['Response_window_complete', 'Window complete'],
      ['Raw_F0', 'Raw F₀', '', '', { decimals: 3 }], ['Raw_baseline_SD', 'Raw baseline SD', '', '', { decimals: 3 }], ['Raw_baseline_samples', 'Raw baseline frames', FRAMES], ['Local_F0_raw', 'Raw local F₀', '', '', { decimals: 3 }],
      ['Registered', 'Registered', CELLS], ['Nonresponders', 'Did not respond', CELLS], ['Excluded', 'Excluded', CELLS], ['Unscorable', 'Unscorable', CELLS],
      ['Control_name', 'Control'], ['Control_effective_end_s', 'Control window end used', SEC], ['KCl_effective_end_s', 'KCl window end used', SEC], ['Control_window_warning', 'Control window warning'],
      ['KCl_window_warning', 'KCl window warning'], ['Control_field_reason', 'Control note'], ['KCl_field_reason', 'KCl note'],
    ],
  };
  const WORDS = Object.fromEntries(Object.entries(COLUMN_WORDS).map(([family, list]) => [family, new Map(list.map(([saved, name, unit, measured, more]) => [saved, Object.assign({ saved, name, unit: unit || '', measured: measured || '' }, more)]))]));
  const FAMILY = { metrics: 'cell', quality: 'cell', calls: 'window', windowMetrics: 'window', counts: 'window', sensitivityCalls: 'window', sensitivityCounts: 'window', negativeControlCalls: 'window', negativeControlCounts: 'window' };
  /** A saved column of a table in words. One the dictionary does not know keeps its saved name, with spaces for the underscores. */
  function columnWords(kind, saved) {
    const known = (WORDS[FAMILY[kind]] && WORDS[FAMILY[kind]].get(saved)) || WORDS.any.get(saved); if (known) return known;
    const timed = /_s$/.test(saved); return { saved, name: saved.replace(/_s$/, '').replace(/_/g, ' '), unit: timed ? SEC : '', measured: '' };
  }
  /** The saved columns of a table that are the same as this one by definition. */
  const twinsOf = (kind, saved) => Array.from((WORDS[FAMILY[kind]] || new Map()).values()).filter((w) => w.alias === saved).map((w) => w.saved);
  const TABLE_NAMES = { metrics: 'Whole-recording metrics', quality: 'Cell quality', registry: 'Registered cells', windowMetrics: 'Response-window metrics', replicates: 'Biological-replicate counts', groups: 'Condition counts',
    unassigned: 'Unassigned biological units', cohort: 'Cohort flow', recordingMetadata: 'Recording details', sensitivityCalls: 'Sensitivity: cell calls', sensitivityCounts: 'Sensitivity: responder counts',
    negativeControlCalls: 'Negative controls: cell calls', negativeControlCounts: 'Negative controls: responder counts' };
  /** The response call of a cell, put beside it in the per-cell metric tables (it is not a column of those files). */
  const CALL_COLUMN = { get name() { return marked() ? 'Call after review' : 'Call'; }, call: true,
    get about() { return marked() ? 'The response call for the primary stimulus after review, from the cell response calls and the review marks' : 'The response call for the primary stimulus, from the cell response calls'; } };
  /** The columns a table shows until "Show all saved columns" is ticked, in this order. */
  const SHORT = {
    calls: ['Cell', 'Role', 'Stimulus', 'Status', 'Reason', 'Threshold_increase', 'Threshold_source', 'Peak_increase', 'Peak_s', 'Onset_s', 'Run_samples', 'Response_latency_s', 'Run_duration_s', 'Window_AUC_fraction_s', 'Baseline_SD', 'QC_flags'],
    windowMetrics: ['Cell', 'Role', 'Stimulus', 'Status', 'Reason', 'Peak_increase', 'Peak_latency_s', 'Onset_s', 'Response_latency_s', 'Run_duration_s', 'Window_AUC_fraction_s', 'QC_flags'],
    metrics: ['Cell', CALL_COLUMN, 'Status', 'Peak_dFF0', 'Time_to_Peak', 'Rise_Time', 'FWHM', 'AUC', 'SNR', 'Baseline_SD', 'F0'],
    quality: ['Cell', CALL_COLUMN, 'Status', 'Reason', 'F0', 'Baseline_SD_raw', 'Baseline_samples', 'QC_flags'],
    counts: ['Stimulus', 'Role', 'Population', 'Registered', 'Responded', 'Nonresponders', 'Evaluable', 'Excluded', 'Unscorable', 'Percent_responding', 'Start_s', 'End_s', 'Effective_end_s', 'Window_warning', 'Floor_local_fraction', 'SD_multiplier', 'Consecutive_samples'],
  };
  /** A number as the tables write it: so many decimals, the true minus, never "−0.000", and "%" after a percentage. */
  function decimal(v, decimals, percent) {
    if (!finite(v)) return '–';
    let text = v.toFixed(decimals); if (/^-[0.]*$/.test(text)) text = text.slice(1);
    return text.replace('-', '−') + (percent ? '%' : '');
  }
  /** A share of cells as the page writes it everywhere: one decimal and "%"; a dash when no cell was evaluable. */
  const shareText = (part, whole) => whole ? decimal(100 * part / whole, 1, true) : '–';
  const primaryWindow = () => { const windows = (S.results && S.results.windows) || []; return windows.find((w) => w.role === 'Stimulus') || windows[0]; };
  /** The rows of the summary, always these and in this order, so that the summaries of different recordings line up: [table,
   *  saved column, then the Metric, Unit and Measured text of the saved summary]. Columns that are the same by definition
   *  (alias) are not rows of their own: "Saved as" names them. The response is timed in the response window; the timing of
   *  the whole recording's peak is named for the whole recording. The screen takes its words from the column dictionary; the
   *  saved summary keeps the words and the saved units (fractions stay fractions) it has always had, with two names made
   *  unique (window area, raw baseline SD) and one unit for both areas, so that summaries saved by any release stack. */
  const SUMMARY_ROWS = [['metrics', 'Peak_dFF0', 'Peak', DFF, 'whole recording'], ['windowMetrics', 'Peak_increase', 'Peak rise', 'fraction of baseline', 'response window'],
    ['metrics', 'SNR', 'Signal-to-noise', 'peak ÷ baseline SD', 'whole recording'], ['metrics', 'AUC', 'Area under the curve', AREA, 'whole recording, signed'],
    ['windowMetrics', 'Window_AUC_fraction_s', 'Window area under the curve', AREA, 'response window, signed'],
    ['metrics', 'Baseline_SD', 'Baseline SD', DFF, 'F₀ frames'], ['metrics', 'F0', 'F₀', GRAY, 'F₀ frames'], ['metrics', 'Baseline_SD_raw', 'Raw baseline SD', GRAY, 'F₀ frames'],
    ['windowMetrics', 'Onset_s', 'Onset time', SEC, 'from recording start'], ['windowMetrics', 'Response_latency_s', 'Latency', SEC, 'from the stimulus'],
    ['windowMetrics', 'Confirmation_s', 'Response confirmed', SEC, 'from recording start'], ['windowMetrics', 'Run_end_s', 'Response end', SEC, 'from recording start'],
    ['windowMetrics', 'Run_duration_s', 'Response duration', SEC, 'above the threshold'], ['windowMetrics', 'Peak_latency_s', 'Peak latency', SEC, 'from the stimulus'],
    ['metrics', 'Time_to_Peak', 'Peak time', SEC, 'from recording start'], ['metrics', 'Time_to_25_Peak', '25 % of peak reached', SEC, 'from recording start'],
    ['metrics', 'Time_to_50_Peak', '50 % of peak reached', SEC, 'from recording start'], ['metrics', 'Time_to_75_Peak', '75 % of peak reached', SEC, 'from recording start'],
    ['metrics', 'Rise_Time', 'Rise time', SEC, '10 % to 90 % of peak'], ['metrics', 'Rise_Rate_10_90_dFF0_per_s', 'Rise rate', RATE, '10 % to 90 % of peak'],
    ['metrics', 'FWHM', 'FWHM', SEC, 'width at half of peak'], ['metrics', 'FWHM_Lower_Bound', 'FWHM, lower bound', SEC, 'still above half at the end'], ['metrics', 'Half_Width', 'Half-width', SEC, 'FWHM ÷ 2']];
  /** Mean, SD, SEM and number of cells of every metric: for all cells that were not excluded, by response call, and by how the
   *  cells were chosen when that was recorded. Worked out here from the saved per-cell tables, which are not changed.
   *  - The rows are SUMMARY_ROWS; the saved file has every one of them (no cells where a column was not saved), the screen
   *    those that were saved. A response-window column comes from the window metrics, or from the calls when there are none.
   *  - The timing and shape of a response are summarised over the cells that responded. On screen such a row gives its value
   *    once, under Responded: All cells would repeat it and Did not respond has no cells.
   *  - On screen the two kinds of choosing are columns only when both have cells; n has a column of its own; a fraction is in %.
   *  - The saved file holds every value as saved (a fraction stays a fraction), at full precision. */
  function metricSummary() {
    const r = S.results, tags = r.chosenBlind || null, primary = primaryWindow(), status = primaryCall;
    const left = (n) => status(n) === 'Excluded' || (S.quality[n] || {}).Status === 'Excluded';
    const groups = [['All cells', (n) => !left(n)], ['Responded', (n) => status(n) === 'Responded'], ['Did not respond', (n) => status(n) === 'Did not respond']];
    if (tags) groups.push(['Chosen blind', (n) => tags[n] === true && !left(n)], ['Added after the response', (n) => tags[n] === false && !left(n)]);
    const both = !!tags && S.info.names.some((n) => tags[n] === true) && S.info.names.some((n) => tags[n] === false), shown = groups.filter((g, i) => i < 3 || both);
    /** The saved table that holds a column, and its rows for the primary stimulus. */
    const source = (kind, saved) => {
      for (const key of kind === 'windowMetrics' ? ['windowMetrics', 'calls'] : [kind]) {
        const data = r.tables[key]; if (!data || !data.header.includes(saved) || !data.header.includes('Cell')) continue;
        const role = data.header.indexOf('Role');
        return { data, rows: data.rows.filter((row) => role < 0 || row[role] === ((primary && primary.role) || 'Stimulus')) };
      }
      return null;
    };
    const long = [['Metric', 'Unit', 'Measured', 'Saved_as', 'Group', 'Mean', 'SD', 'SEM', 'N_cells']], wide = [];
    for (const [kind, saved, ...file] of SUMMARY_ROWS) {
      const words = columnWords(kind, saved), found = source(kind, saved), scale = words.fraction ? 100 : 1, decimals = words.decimals != null ? words.decimals : DECIMALS[words.unit];
      const names = found ? [saved].concat(twinsOf(kind, saved).filter((twin) => found.data.header.includes(twin))) : [saved], values = new Map();
      if (found) { const cell = found.data.header.indexOf('Cell'), c = found.data.header.indexOf(saved);
        for (const row of found.rows) { const v = row[c] === '' || row[c] == null ? NaN : Number(row[c]); if (isFinite(v)) values.set(row[cell], v); } }
      const stats = groups.map(([group, belongs]) => {
        const xs = []; for (const [cell, x] of values) if (belongs(cell) && (!words.response || status(cell) === 'Responded')) xs.push(x);
        const k = xs.length, mean = k ? xs.reduce((a, b) => a + b, 0) / k : NaN, sd = k > 1 ? Math.sqrt(xs.reduce((a, b) => a + (b - mean) * (b - mean), 0) / (k - 1)) : NaN;
        return { group, k, mean, sd, sem: k > 1 ? sd / Math.sqrt(k) : NaN };
      });
      for (const s of stats) long.push([...file, found ? names.join('; ') : '', s.group, s.k ? s.mean : '', s.k > 1 ? s.sd : '', s.k > 1 ? s.sem : '', s.k]);
      if (!found) continue;
      const line = [{ text: words.name, title: 'Saved as ' + names.join(', ') + (words.fraction ? ', as fractions' : '') }, { text: words.unit }, { text: words.measured }];
      for (const s of stats.filter((s) => shown.some((g) => g[0] === s.group))) {
        if (words.response && (s.group === 'All cells' || s.group === 'Did not respond')) { line.push({ text: '', key: '', num: true }, { text: '', key: '', num: true }); continue; }
        const mean = scale * s.mean, sd = scale * s.sd, sem = scale * s.sem;
        line.push(s.k ? { text: decimal(mean, decimals, false) + (s.k > 1 ? ' ± ' + decimal(sem, decimals, false) : '') + (words.unit === PCT ? '%' : ''), title: 'Mean ' + mean + (s.k > 1 ? ', SD ' + sd + ', SEM ' + sem : ''), key: mean, num: true } : { text: '–', key: '', num: true });
        line.push({ text: String(s.k), key: s.k, num: true });
      }
      wide.push(line);
    }
    const columns = [{ name: 'Metric' }, { name: 'Unit' }, { name: 'Measured' }].concat(...shown.map(([group]) => [{ name: group, title: 'Mean ± SEM of the cells in this group' }, { name: 'n', title: 'Cells in ' + group.toLowerCase() }]));
    return { wide: { header: columns.map((c) => c.name), columns: columns.map((c, index) => Object.assign({ given: true, index }, c)), rows: wide }, long };
  }
  const csvText = (rows) => rows.map((row) => row.map((v) => { const s = v == null ? '' : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }).join(',')).join('\n') + '\n';
  /** The table on the Metrics tab as CSV rows: the summary in long form (one row per metric and group), any other table as saved. */
  function metricRows() { const key = $('metrics-table').value; if (key === 'summary') return metricSummary().long; return savedRows(key, true); }
  /** A saved table as CSV rows, with its saved names and every column; coded: with the codes of the call checks in its QC_flags (withCodes). */
  const savedRows = (key, coded) => { let data = (S.results && S.results.tables[key]) || { header: [], rows: [] }; if (coded) data = withCodes(key, data).data; return [data.header].concat(data.rows); };
  /** Whether the table on the Metrics tab, as exported, has codes of the call checks in its QC_flags. */
  const metricCoded = () => { const key = $('metrics-table').value; return key !== 'summary' && !!S.results && withCodes(key, S.results.tables[key]).added; };
  const METRIC_HELP = {
    summary: 'Mean ± SEM. Onset, latency, time to peak, rise and width are for cells that responded. Point at a metric for its saved columns; the saved file also has the SD.',
    metrics: 'Measured over the whole recording. AUC is signed and turns negative when the signal drifts down. Times are from the start of the recording; rise time, from 10% to 90% of the peak, can follow noise. For the response itself, choose Response-window metrics.',
    windowMetrics: 'Measured in the response window against the baseline just before the stimulus. Latencies are from the stimulus, to the nearest frame.',
  };
  function showMetricTable() {
    const key = $('metrics-table').value; $('metric-help').textContent = METRIC_HELP[key] || ''; $('all-metric-columns-label').hidden = !SHORT[key];
    if (key === 'summary') table('table-metrics', metricSummary().wide); else table('table-metrics', withCodes(key, S.results.tables[key]).data, key);
  }

  function draftChanged() {
    if (!S.draft) return;
    S.revision++; S.note = null;
    clearTimeout(S.draftTimer); S.draftTimer = setTimeout(() => syncDraft().catch(draftError), 350);
    updateProtocolControls(); showState();
  }
  function draftError(error) { S.note = error.message; showState(); }
  async function syncDraft() {
    clearTimeout(S.draftTimer);
    const draft = collectDraft(), signature = JSON.stringify(draft), revision = S.revision;
    S.draftQueue = S.draftQueue.catch(() => {}).then(async () => {
      const state = await api('draft', { draft, revision });
      if (Number(state.draftRevision) === revision) S.draftAck = signature;
      acceptState(state); showState(); return state;
    });
    return S.draftQueue;
  }
  function acceptState(state) {
    // A poll begun before a draft save must not erase the pending warning or export provenance.
    if (S.state && (Number(state.draftRevision) < Number(S.state.draftRevision) || Number(state.version) < Number(S.state.version))) return false;
    S.state = state; if (Buf.version !== (state.frameVersion || '')) resetFrames(state.frameVersion || ''); return true;
  }
  function hasPending() { return !!S.results && (!!(S.state && S.state.settingsPending) || JSON.stringify(collectDraft()) !== S.draftAck); }
  async function start() {
    S.info = await api('session');
    $('save-picture-folder').disabled = $('save-picture').disabled = !(S.info.hasPicture || S.info.hasMovie); $('save-frames').disabled = !S.info.hasMovie;
    const prefs = S.info.prefs || {}; for (const key of Object.keys(Layout.apply)) if (prefs[key]) Layout.apply[key](prefs[key]);
    // How pictures were last saved. The page's lists show it from the start, so that the figure editor, which shows those lists, never reads a format that is not the one in use.
    for (const key of Object.keys(Save)) { const kept = prefs['step4.' + key]; if (kept && kept !== String(Save[key])) Fig.earlier.save[key] = Save[key] = key === 'dpi' ? Number(kept) : kept; }
    showSaveOptions();
    for (const kind of ['plot', 'heatmap']) if (prefs['step4.figure.' + kind]) figureRestore(kind, prefs['step4.figure.' + kind]);
    if (prefs['step4.outline']) { S.outlineWidth = Number(prefs['step4.outline']) || 1; $('outline-width').value = prefs['step4.outline']; }
    if (prefs['step4.plain'] === '1') { S.outlineColours = 'plain'; $('outline-colours').value = 'plain'; } if (prefs['step4.numbers'] === '1') { S.numberCells = true; $('number-cells').checked = true; }
    if (prefs['step4.viewhint']) Pic.hinted = Number(prefs['step4.viewhint']) || 0;
    if (prefs['step4.every']) { try { window.localStorage.setItem('step4.every', prefs['step4.every']); } catch (blocked) { } }
    if (prefs['step4.help']) { const on = prefs['step4.help'] === '1'; document.body.classList.toggle('show-help', on); $('toggle-help').textContent = on ? 'Hide explanations' : 'Show explanations'; }
    $('recording').textContent = S.info.recording; $('recording').title = S.info.recording;
    document.title = S.info.names.length + ' cells · Step 4 · ' + S.info.recording + ' · ' + S.info.key.slice(0, 6); $('version').textContent = 'Ca²⁺ Imaging (Wang Lab) ' + S.info.pluginVersion;
    S.initialDraft = clone(S.info.draft);
    S.revision = Math.max(0, Number(S.info.draftRevision) || 0); fillDraft(S.info.draft); S.draftAck = JSON.stringify(collectDraft());
    $('time').max = Math.max(1, S.info.frames);
    for (const option of Array.from($('view').options)) if ((option.value === 'picture' && !S.info.hasPicture) || (option.value === 'background' && !S.info.hasBackground) || (!['picture', 'background'].includes(option.value) && !S.info.hasMovie)) option.remove();
    if (S.info.hasBackground) { const areas = new Image(); areas.onload = () => { Pic.background = areas; if (S.view === 'background') requestFrame(); }; areas.src = 'api/background'; }
    S.view = S.info.hasMovie ? 'raw' : 'picture'; $('view').value = S.view;
    S.paths = S.info.rois.map((points) => { const p = new Path2D(); points.forEach((q, i) => i ? p.lineTo(q[0], q[1]) : p.moveTo(q[0], q[1])); p.closePath(); return p; });
    S.bounds = S.info.rois.map((points) => { let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const q of points) { x0 = Math.min(x0, q[0]); y0 = Math.min(y0, q[1]); x1 = Math.max(x1, q[0]); y1 = Math.max(y1, q[1]); } return { x0, y0, x1, y1 }; });
    // for the cells' numbers: the middle of each outline (the centre of its area; of its bounds when it has none) and how wide a typical cell is
    S.centres = S.info.rois.map((points, i) => { let area = 0, x = 0, y = 0; points.forEach((p, k) => { const q = points[(k + 1) % points.length], cross = p[0] * q[1] - q[0] * p[1]; area += cross; x += (p[0] + q[0]) * cross; y += (p[1] + q[1]) * cross; });
      const b = S.bounds[i]; return Math.abs(area) > 1e-9 ? [x / (3 * area), y / (3 * area)] : [(b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2]; });
    { const sizes = S.bounds.map((b) => Math.sqrt((b.x1 - b.x0) * (b.y1 - b.y0))).filter((size) => finite(size) && size > 0).sort((a, b) => a - b); S.cellSize = sizes.length ? sizes[sizes.length >> 1] : 20; }
    if (S.info.hasPicture) { const image = new Image(); image.onload = () => { Pic.base = image; if (S.view === 'picture') { Pic.shown = image; movieMessage(''); } drawPicture(); }; image.onerror = () => movieMessage('The static overview could not be loaded. Choose Recording to view the movie.', true); image.src = 'api/picture'; }
    $('picture-card').hidden = !(S.info.hasPicture || S.info.hasMovie);
    listCells(); setFrame(1); await poll(); setInterval(poll, 600);
  }
  async function poll() {
    if (S.polling || S.closed) return; S.polling = true;
    try {
      const state = await api('state'); if (S.closed) return; S.lost = 0; $('lost').classList.remove('on');
      if (S.info && state.key !== S.info.key) { location.reload(); return; }
      const wasBusy = S.state && S.state.busy; if (!acceptState(state)) return;
      if (wasBusy && !state.busy) S.note = null;
      // only a newer analysis reloads the results: a mark's own answer brings its review and tally (applyReview), so a poll
      // that crossed it, or one made while it is on its way, does not load the whole results again
      if (state.version > 0 && state.version > S.version && !S.reviewing) { await loadResults(); S.version = state.version; S.finishing = !!state.busy; }
      // results read while the run was still finishing (Fiji's load() moves the version on before the finish message and its
      // tone): once the run is over, the status line is what loadResults would have made of it
      else if (S.finishing && !state.busy) { S.finishing = false; S.said = saidNow(); }
      askAroundAgain();   // Fiji answers again: what did not come while it was not connected is asked again
      showState(); updateTransport();
    } catch (problem) { if (++S.lost >= 3) $('lost').classList.add('on'); }
    finally { S.polling = false; }
  }
  function showState() {
    const state = S.state; if (!state) return;
    const engine = $('engine'); engine.textContent = state.engineReady ? 'Analysis engine ready' : 'First run: set up analysis engine'; engine.className = 'pill ' + (state.engineReady ? 'good' : 'warn'); engine.hidden = !!state.engineReady;
    if (S.reapplyFrames && editedTimes() && JSON.stringify(collectDraft()) === S.draftAck) { S.reapplyFrames = false; const before = JSON.stringify(collectDraft()); applyAllFrames(); if (JSON.stringify(collectDraft()) !== before) { draftChanged(); return; } }
    const validation = validateDraft(collectDraft());
    $('analyze').disabled = state.busy || S.saving || !state.engineReady || validation.errors.length > 0 || JSON.stringify(collectDraft()) !== S.draftAck || S.closed; $('log').disabled = !state.logAvailable;
    $('cancel').hidden = !state.busy; $('cancel').disabled = !state.cancellable; $('cancel').title = state.busy && !state.cancellable ? 'The checks run to the end.' : '';
    $('setup').hidden = state.engineReady; $('setup').disabled = state.busy; $('setup').classList.toggle('primary', !state.engineReady);
    $('analyze').textContent = state.busy ? 'Analysis in progress…' : S.results ? (hasPending() ? 'Apply changes & reanalyze' : 'Run analysis again') : 'Run analysis';
    for (const box of document.querySelectorAll('#side input,#side select,#side textarea,#protocol-export')) box.disabled = state.busy || S.saving;
    updateProtocolControls();
    $('sampling_hz').disabled = state.busy || S.saving || !!S.draft.timestamps;
    $('close-page').disabled = state.busy || S.saving || S.closed;
    $('discard').hidden = JSON.stringify(collectDraft()) === JSON.stringify((S.results && S.results.appliedDraft) || S.initialDraft);
    $('discard').disabled = state.busy || S.saving;
    const progress = $('progress'); progress.hidden = !state.busy; progress.classList.toggle('busy', state.busy && state.installing && !(state.percent > 0));
    $('bar').style.width = Math.max(3, state.percent || 0) + '%';
    const status = $('status');
    if (state.busy) { status.className = 'status'; status.textContent = (state.percent > 0 ? state.percent + '% · ' : '') + state.stage + ' · ' + state.elapsed + ' s elapsed'; }
    else if (S.note) { status.className = 'status bad'; status.textContent = S.note; }
    else if (S.imported && S.imported.revision === S.revision) { status.className = 'status good'; status.textContent = S.imported.text; }
    else if (!S.results && !state.error && state.tone !== 'bad') { status.className = 'status'; status.textContent = state.engineReady ? (state.message && !/^Check acquisition timing/.test(state.message) ? state.message : 'Check the timing, tick the box, then run.') : 'First run: choose “Set up analysis engine” (needs internet).'; }
    else if (S.results && state.tone === 'good' && !called()) { status.className = 'status warn'; status.textContent = 'Finished, but no cell got a call' + (unscorableReason() ? ': ' + unscorableReason().replace(/\.$/, '').toLowerCase() : '') + '.'; }
    else if (S.results && state.tone === 'good' && state.folder) { status.className = 'status good'; status.textContent = finishLine(); }
    else { status.className = 'status ' + (state.tone === 'good' ? 'good' : state.tone === 'bad' ? 'bad' : state.tone === 'warn' ? 'warn' : ''); status.textContent = state.message; }
    status.title = state.message || ''; status.classList.toggle('full', !state.busy && !!(S.carried || S.said) && status.classList.contains('good'));   // carried marks and what a mark changed are read in full
    if (state.busy && S.marking) closeMark();
    showMark();
    // what the tick box confirms, in the user's own numbers
    { const d = collectDraft(), hz = Number(d.values.sampling_hz), first = d.protocol.windows[0], used = first ? windowFrames(expandWindow(first, d.protocol.schema_version === 1), editedTimes()) : null, st = used && used.stimulus;
      $('confirm-text').textContent = 'Timing checked' + (d.timestamps || (isFinite(hz) && hz > 0) ? ': ' + (d.timestamps ? 'recorded timestamps' : fmt(+hz.toPrecision(4)) + (fmt(+hz.toPrecision(4)) === '1' ? ' frame/s' : ' frames/s')) + (st ? ', stimulus ' + (st.before >= 1 && st.before < used.n ? 'after frame ' + st.before : 'outside the recording') : '') : '');
      // the recording's own frame rate, when a different one is typed
      const own = Number(S.info.recordedHz), note = $('rate-note'), known = !d.timestamps && isFinite(own) && own > 0 && isFinite(hz) && hz > 0, differs = known && Math.abs(hz - own) > 1e-9 * own;
      const text = differs ? 'The file says ' + fmt(+own.toPrecision(6)) + ' (' + fmt(+(1 / own).toPrecision(6)) + ' s per frame). ' : known ? 'From the file' + (rateText(own) === String(own) ? '.' : ': ' + fmt(+own.toPrecision(8)) + '.') : '';
      if (note.dataset.shown !== text) { note.dataset.shown = text; note.hidden = !text; note.replaceChildren(text);
        if (differs) { const use = document.createElement('button'); use.type = 'button'; use.className = 'ghost repair'; use.textContent = 'Use it'; use.addEventListener('click', () => { exactInput($('sampling_hz'), S.info.recordedHz); applyAllFrames(); draftChanged(); }); note.appendChild(use); } } }
    $('open-folder').disabled = !state.folder; $('open-folder').hidden = !state.folder; $('open-folder').title = state.folder ? 'Open ' + String(state.folder).split(/[\\/]/).pop() : '';
    const pending = hasPending(), notice = $('applied-note'); notice.hidden = !S.results || !pending; notice.classList.toggle('pending', pending);
    notice.textContent = pending ? 'Settings changed. Everything shown is still the last completed analysis.' : '';
    $('workspace-title').textContent = S.results ? 'Review the responses' : 'Review the recording';
    $('review-stage').textContent = S.results ? (pending ? 'Changes pending' : 'Completed analysis') : 'Preview · no analysis yet';
    $('workflow-help').textContent = S.results ? 'Select a cell or response window. Use the same timeline for the movie and traces.' : 'Drag the frame slider or press Play. You can review the movie before running analysis.';
  }
  /** What Fiji says after a good finish of the analysis shown, while marks apply or when its run made an automatic mark (the
   *  message then names the cells marked, until the next run); null otherwise (finishLine says "Done."). */
  function saidNow() { const now = S.state, r = S.results; return now && r && (marked() || automaticMade(now.message)) && now.tone === 'good' && now.folder === r.folder ? (now.message || null) : null; }
  /** Whether the analysis shown was run with "Mark possible passing objects as not responders" ticked. */
  const automaticRun = () => { const d = S.results && S.results.appliedDraft; return !!(d && d.values && d.values.automatic_marks === true); };
  /** Whether that run made at least one automatic mark: one is among the marks, or Fiji's message names the cells it marked,
   *  as it does until the next run, also once they are removed. */
  const automaticMade = (message) => automaticRun() && (Object.values(marked() || {}).some(automaticMark) || /(^|\s)Marked automatically: /.test(message || ''));
  /** The status line after a good finish: "Done." and the marks carried over from the earlier analysis; after a mark, what Fiji
   *  said then. Both are Fiji's own words; the carried marks stay until the next run. */
  function finishLine() {
    const carried = S.carried ? S.carried.text : '';
    if (S.said) return S.said + (carried && !S.said.includes(carried) ? ' ' + carried : '');
    return 'Done.' + (carried ? ' ' + carried : '');
  }
  async function send(path, body) {
    S.note = null;
    try { acceptState(await api(path, body || {})); showState(); }
    catch (problem) { S.note = problem.message; showState(); }
  }
  // ---------------------------------------------------------------- saved results; selecting a window never recalculates a call
  async function loadResults() {
    S.results = await api('results'); const r = S.results;
    Around.generation++; Around.got.clear(); Around.asking.clear();   // "Around the cell" is asked for again, for these results
    // the marks carried over when this analysis was run: Fiji gives them for this analysis until the next run, marks or none
    const carried = r.review && r.review.carried; S.carried = carried && carried.text ? { folder: r.folder, text: carried.text } : null;
    // while marks apply after a good finish of this analysis, the status line is what Fiji says now, as after a mark (applyReview)
    // and as in the Fiji window; a page opened on another analysis, a failed or empty run, or no marks: "Done." (finishLine)
    S.said = saidNow();
    S.markSaid = null; S.keep = null; closeMark();
    const records = (data) => Object.fromEntries(((data && data.rows) || []).map((row) => [row[data.header.indexOf('Cell')], Object.fromEntries(data.header.map((key, i) => [key, row[i]]))]));
    S.quality = records(r.tables.quality);
    const windows = r.windows && r.windows.length ? r.windows : [{ key: 'Stimulus::' + (r.counts.Stimulus || 'Stimulus'), role: 'Stimulus', stimulus: r.counts.Stimulus || 'Stimulus' }];
    const selector = $('response-window'); selector.innerHTML = '';
    for (const w of windows) { const o = document.createElement('option'); o.value = w.key; o.textContent = w.stimulus + (w.role === 'Control' ? ' · control' : ' · primary'); selector.appendChild(o); }
    if (!windows.some((w) => w.key === S.window)) S.window = windows[0].key; selector.value = S.window;
    document.body.classList.add('has-results');
    for (const id of ['result-tiles', 'plot-card', 'group-label', 'window-select-label', 'every-label']) $(id).hidden = false;   // ("Thresholds": drawPlot)
    $('window-select-label').hidden = !(S.results.windows && S.results.windows.length > 1);
    for (const tab of document.querySelectorAll('.tab')) tab.disabled = false;
    const tables = $('metrics-table'); tables.innerHTML = '';
    { const o = document.createElement('option'); o.value = 'summary'; o.textContent = 'Summary: mean ± SEM of each metric'; tables.appendChild(o); }
    for (const key of Object.keys(r.tables || {})) if (!['calls', 'counts'].includes(key) && r.tables[key].header.length) { const o = document.createElement('option'); o.value = key; o.textContent = TABLE_NAMES[key] || key.replace(/_/g, ' '); tables.appendChild(o); }
    tables.value = 'summary'; showMetricTable();
    S.withoutBelowZero = false;   // every cell is counted until the switch says otherwise
    showCalls(); showCounts(); showTally(); showMethods();
    for (const id of ['save-plot', 'save-heatmap', 'save-plot-folder', 'save-heatmap-folder', 'edit-plot', 'edit-heatmap', 'data-save', 'data-download', 'metrics-save', 'metrics-download', 'calls-save', 'calls-download', 'counts-save', 'counts-download']) $(id).disabled = false;
    for (const id of ['plot-saved', 'heatmap-saved', 'data-saved', 'metrics-saved', 'calls-saved', 'counts-saved']) $(id).textContent = '';
    if (!$('panel-data').hidden) showData();
    for (const option of $('signal').options) option.disabled = option.value !== 'dff' && !hasTable(option.value);
    if ($('signal').selectedOptions[0].disabled) { S.signal = Data.signal = 'dff'; $('signal').value = 'dff'; $('data-signal').value = 'dff'; }
    Plot.x0 = Plot.x1 = null; updateWindow(); setFrame(S.frame); showState();
  }
  /** % responding to the primary stimulus as worked out in Fiji (results.tally), by row key: all, blind, later, quiet.
   *  With the switch on, the cells whose corrected fluorescence went below zero are left out of every row. While review marks
   *  apply the figures are those after review; before: the responders before review, changed: the calls a mark changed. */
  function tallyRows() {
    const t = S.results && S.results.tally; if (!t || !t.rows) return null;
    const rows = {}, reviewed = !!(t.reviewed && marked()), side = (from) => S.withoutBelowZero ? from.without : from.counted;
    for (const row of t.rows) { const after = reviewed && row.reviewed ? row.reviewed : null, c = side(after || row);
      rows[row.key] = { label: row.label, responded: c[0], evaluable: c[1], share: shareText(c[0], c[1]), reviewed: !!after, before: side(row)[0], changed: after ? after.changed : 0 }; }
    return rows;
  }
  /** The same tally as a table on the Responder counts tab: every row with all its cells and without those below zero. */
  function showTally() {
    // the engine's counts under the tally are as saved, before review: while marks apply, a line above them says so
    $('counts-before').hidden = !marked();
    const t = S.results && S.results.tally; $('tally-box').hidden = !t; if (!t) return;
    const share = (c) => ({ text: shareText(c[0], c[1]), key: c[1] ? c[0] / c[1] : '', num: true }), count = (n) => ({ text: String(n), key: n, num: true });
    // six columns, so the table fits its tab; the saved file has every number in a column of its own. While review marks
    // apply, two more: the share after review and how many calls the marks changed.
    // After review follows "Leave out … below zero" as the headline does: the figure written there is the one in this column
    const reviewed = !!(t.reviewed && marked()), of = (c) => c[1] ? { text: shareText(c[0], c[1]) + ' (' + c[0] + ' of ' + c[1] + ')', key: c[0] / c[1], num: true } : { text: '–', key: '', num: true };
    const after = (row) => S.withoutBelowZero ? row.reviewed.without : row.reviewed.counted;
    const header = ['Cells', '% responding', 'Responded', 'Evaluable', 'Cells below zero', 'Without them'].concat(reviewed ? ['After review', 'Calls changed'] : []);
    table('table-tally', { header, columns: header.map((name, index) => ({ name, index, given: true })),
      rows: t.rows.map((row) => [row.label, share(row.counted), count(row.counted[0]), count(row.counted[1]), count(row.counted[2] - row.without[2]), of(row.without)]
        .concat(reviewed ? (row.reviewed ? [of(after(row)), count(row.reviewed.changed)] : [{ text: '–', key: '', num: true }, { text: '–', key: '', num: true }]) : [])) });
  }
  /** The headline numbers for the selected response window (the primary stimulus or its control), from the saved counts. */
  function showTiles() {
    const r = S.results; if (!r) return;
    const windows = r.windows && r.windows.length ? r.windows : [], chosen = windows.find((w) => w.key === S.window) || windows[0] || {};
    const counts = (r.countsByWindow && r.countsByWindow[S.window]) || r.counts || {}, whole = (key) => { const v = number(counts[key]); return isFinite(v) ? String(Math.round(v)) : '–'; }, percent = number(counts.Percent_responding);
    const label = escape(counts.Stimulus || chosen.stimulus || 'Stimulus') + (chosen.role === 'Control' ? ' (control)' : '');
    const all = label + ': ' + whole('Responded') + ' of ' + whole('Evaluable') + ' evaluable cells';
    // The primary stimulus: all cells lead; under them the cells chosen blind and those added later, then the quiet cells.
    // Every figure comes from the tally worked out in Fiji, with or without the cells that went below zero.
    const tally = chosen.role !== 'Control' ? tallyRows() : null, below = tally ? r.tally.belowZero.length : 0, said = (row) => row.share + ' (' + row.responded + ' of ' + row.evaluable + (row.changed > 0 ? '; ' + row.before + ' before review' : '') + ')';
    $('below-zero-label').hidden = !below; $('below-zero').checked = !!S.withoutBelowZero;
    $('below-zero-text').textContent = 'Leave out ' + below + (below === 1 ? ' cell' : ' cells') + ' below zero';
    const reason = unscorableReason();
    if (!called()) { $('tile-responding').innerHTML = '<b>–</b><span>' + label + ': no cell got a call' + (reason ? ' (' + escape(reason.replace(/\.$/, '').toLowerCase()) + ')' : '') + '</span>'; }
    else if (tally && tally.all) {
      const lines = [];
      if (tally.blind && tally.later && tally.blind.evaluable && tally.later.evaluable) lines.push(['Chosen blind: ' + said(tally.blind) + ' · added after the response was shown: ' + said(tally.later), 'Cells chosen before the response was seen, and cells added once it was shown']);
      if (tally.quiet && tally.quiet.evaluable) lines.push([tally.quiet.label + ': ' + said(tally.quiet), 'Cells that need a rise of ' + Math.round(100 * r.tally.quietRise) + ' % or less to respond. Noisier cells are set aside whether they responded or not']);
      // with review marks the figure is the one after review, and the count before review is beside it (as in the Fiji window)
      $('tile-responding').innerHTML = '<b>' + tally.all.share + '</b><span>' + label + ': ' + tally.all.responded + ' of ' + tally.all.evaluable + ' evaluable cells' + (tally.all.reviewed ? ' after review; ' + tally.all.before + ' before review' : '') + (S.withoutBelowZero && below ? ' · ' + below + ' below zero left out' : '') + '</span>'
        + lines.map(([text, title]) => '<span class="tile-line" title="' + escape(title) + '">' + escape(text) + '</span>').join('');
    }
    else $('tile-responding').innerHTML = '<b>' + decimal(percent, 1, true) + '</b><span>' + all + '</span>';
    $('tile-registered').firstChild.textContent = whole('Registered'); $('tile-excluded').firstChild.textContent = whole('Excluded'); $('tile-unscorable').firstChild.textContent = whole('Unscorable');
    for (const [id, key] of [['tile-excluded', 'Excluded'], ['tile-unscorable', 'Unscorable']]) $(id).classList.toggle('zero', whole(key) === '0');
  }
  function updateWindow() {
    // the flagged list first: a filter whose flag this window does not have is let go before Show counts its group
    showTiles(); showQuality();
    const group = $('group'); group.innerHTML = '';
    // while review marks apply to the window on screen, the two calls are those after review, and the marked cells are a group of their own
    const reviewing = !!marked() && primaryShown(); if (!reviewing && S.group === 6) S.group = 0;
    GROUPS.forEach((name, i) => { if (i === 6 && !reviewing) return; const o = document.createElement('option'); o.value = i;
      o.textContent = name + (reviewing && (i === 1 || i === 2) ? ' after review' : '') + ' (' + S.info.names.filter((n) => inGroup(callOf(n), i, n)).length + ')'; group.appendChild(o); });
    group.value = S.group; Plot.dirty = true; listCells(); drawHeatmap();
  }
  /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. A cell that a mark or its
   *  removal took out of the group (S.keep, refreshReview) stays in the list, in its place, while it is selected: the arrows go
   *  on from it to the next cell of the group. It is not counted or drawn with the group. */
  function listCells() {
    const names = S.info.names;
    S.listed = S.results ? names.filter((n) => inGroup(callOf(n), S.group, n) || n === S.keep) : names.slice();
    offerCells();
    pick(S.listed.includes(S.cell) ? S.cell : null);   // the selection never changes by itself: a cell that is not listed is no longer selected
  }
  /** The Cell box's list and the arrows, from the cells listed. */
  function offerCells() {
    const offered = $('cell-names'); offered.replaceChildren();
    for (const name of S.listed) { const o = document.createElement('option'); o.value = shortName(name); offered.appendChild(o); }
    $('previous').disabled = $('next').disabled = S.listed.length < 1;
  }
  /** Show is set to another group, by the person or by a flag: the cell a mark kept listed is let go first. */
  function showGroup(group) { if (S.group !== group) { S.group = group; $('group').value = group; S.keep = null; listCells(); } }
  /** Selects a cell everywhere: the list, the picture, the plot, the heatmap and the tables. */
  function pick(name) {
    // the selection moves away from a cell kept listed after its mark: it leaves the list when it is not in the group
    if (S.keep && name !== S.keep) { const kept = S.keep; S.keep = null; if (S.results && !inGroup(callOf(kept), S.group, kept)) { S.listed = S.listed.filter((n) => n !== kept); offerCells(); } }
    if (name && S.results && !S.listed.includes(name)) { S.group = 0; $('group').value = 0; S.cell = name; listCells(); return; }
    S.cell = name || null; $('cell').value = name ? shortName(name) : '';
    const verdict = $('verdict'), row = name && S.results ? calls()[name] : null;
    if (!name) verdict.textContent = !S.results ? '' : S.listed.length ? 'No cell selected.' : 'No cells in this group.';
    else if (!row) verdict.innerHTML = S.results ? '<b>' + escape(shortName(name)) + '</b> · no response call was saved for this cell' : '';
    else {
      // a responder: what the call rests on, the rise held long enough, in place of the engine's "sustained threshold crossing".
      // A marked cell leads with its call after review and the reason, then its automatic call; the rest describes the automatic call.
      const basis = row.Status === 'Responded' ? heldRise(row) : '', mark = markOf(name), call = callOf(name);
      const parts = ['<span class="dot" style="background:' + CALL[kind(call)] + '"></span><b>' + escape(shortName(name)) + ': ' + escape(call || '') + (mark ? ' after review' : '') + '</b>'];
      if (mark) parts.push('<span title="' + attribute(markedAt(mark)) + '">' + escape(mark.reason.toLowerCase() + (mark.note ? ': ' + mark.note : '')) + '</span>',
        ...(automaticMark(mark) ? ['marked automatically'] : []), ...(mark.carriedFrom ? ['<span title="' + attribute(mark.carriedFrom) + '">from an earlier analysis</span>'] : []), 'automatic call ' + escape(mark.automatic));
      parts.push(escape(basis || (row.Reason || '').toLowerCase()));
      const rise = number(row.Peak_increase), at = number(row.Peak_s), need = number(row.Threshold_increase), setBy = { noise: ', set by its noise', floor: ', the minimum rise' }[requiredBy(row)] || '';
      // a responder's largest rise that lies outside the rise that held says so, so that its time does not read against the onset
      const held = row.Status === 'Responded' ? [heldStart(row), number(row.Run_end_s)] : [NaN, NaN], apart = !isFinite(at) ? '' : at < held[0] - 1e-9 ? ', not held' : at > held[1] + 1e-9 ? ', after the held rise' : '';
      if (isFinite(rise)) parts.push('largest rise ' + decimal(100 * rise, 1, true) + (isFinite(at) ? ' at ' + decimal(at, 1) + ' s' + apart : '') + (isFinite(need) ? ' (needed ' + decimal(100 * need, 1, true) + setBy + ')' : ''));
      verdict.innerHTML = parts.filter(Boolean).join(' · ');
      verdict.title = 'Rise and threshold are shown as percentages of local prestimulus fluorescence. The saved call table stores the corresponding fractions.';
    }
    const flags = qualityFor(name).filter((flag) => !isCheck(flag)); $('cell-quality').hidden = !flags.length; $('cell-quality').textContent = flags.map(flagText).join(' · '); $('cell-quality').title = $('cell-quality').textContent;
    // a refusal belongs to its cell: it goes when another cell is chosen, as in the Fiji window
    if (S.markSaid && S.markSaid.cell !== name) S.markSaid = null;
    showChecks(name); if (S.marking && S.marking.cell !== name) closeMark(); showMark();
    if (name && S.results) askAround(name); showAroundButton();
    $('zoom-cell').disabled = !name || !(S.info.hasPicture || S.info.hasMovie);
    for (const row of document.querySelectorAll('tbody tr')) row.classList.toggle('on', row.dataset.cell === name);
    if (Pic.zoom > 1 && name) centreOn(S.info.names.indexOf(name));
    Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap();
  }
  /** The rise a response call rests on, from the saved call: "rise held for 2 frames from 121.0 s, 102.0 s after the stimulus"
   *  (onset, run length and latency). Empty when the call saved no onset. */
  const heldStart = (row) => number(row.Onset_s !== undefined && row.Onset_s !== '' ? row.Onset_s : row.Run_start_s);
  function heldRise(row) {
    const onset = heldStart(row), frames = number(row.Run_samples), after = number(row.Response_latency_s);
    if (!isFinite(onset)) return '';
    return 'rise held ' + (isFinite(frames) ? 'for ' + plural(Math.round(frames), 'frame') + ' ' : '') + 'from ' + decimal(onset, 1) + ' s' + (isFinite(after) ? ', ' + decimal(after, 1) + ' s after the stimulus' : '');
  }
  const escape = (text) => String(text).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const attribute = (text) => escape(text).replace(/"/g, '&quot;');
  /** When a mark was made, for its tooltip: "Marked 2026-11-24 10:05". */
  const markedAt = (mark) => mark.time ? 'Marked ' + String(mark.time).replace('T', ' ').replace(/:\d\d(\.\d+)?([+-]\d\d:?\d\d|Z)?$/, '') : '';

  // ---------------------------------------------------------------- the selected cell's checks, and its review mark
  /** The first frame a check flag names: the first of its frames as Fiji gives them, or, from an answer without them, the first
   *  frame its own text names ("frames 9 to 12", "frame 7"); null when it names none, as short separate crossings. */
  const flagFirst = (flag) => {
    if (flag.frames !== undefined) return Array.isArray(flag.frames) ? flag.frames[0] : null;
    const named = /\bframes? (\d+)/.exec(flag.text || ''); return named ? Number(named[1]) : null;
  };
  /** The selected cell's check flags, under its call. A flag goes to the first frame it names, or, when it names none, to the
   *  first frame of the cell's event (from a still picture to the recording, when there is one); without either it is text. An
   *  evaluable cell that the checks could not look at, in an analysis whose rings were measured, says why. */
  function showChecks(name) {
    const box = $('cell-checks'); box.replaceChildren();
    const c = shownChecks(name), checks = S.results && S.results.checks, items = [];
    if (c) {
      const frameOf = (frame) => Number.isInteger(frame) && frame >= 1 && frame <= S.info.frames ? frame : null;
      for (const flag of checkFlags(name)) {
        const frame = frameOf(flagFirst(flag)) || frameOf(c.firstFrame), item = document.createElement(frame ? 'button' : 'span'); item.textContent = flag.text;
        if (frame) { item.type = 'button'; item.className = 'jump'; item.title = 'Go to frame ' + frame + (c.detail ? '\n' + c.detail : ''); item.addEventListener('click', () => goToFrame(frame)); }
        else if (c.detail) item.title = c.detail;
        items.push(item);
      }
      if (checks.ringMeasured && !c.checked && c.why && ['Responded', 'Did not respond'].includes(automaticCall(name))) { const item = document.createElement('span'); item.textContent = 'passing object not checked: ' + c.why; items.push(item); }
    }
    items.forEach((item, i) => { if (i) box.append(' · '); box.appendChild(item); });
    box.hidden = !items.length;
  }
  /** A check flag's frame: the movie goes there, and from the colour overview or the background areas to the recording. */
  function goToFrame(frame) {
    if (still() && S.info.hasMovie && Array.from($('view').options).some((o) => o.value === 'raw')) { play(false); S.view = 'raw'; $('view').value = 'raw'; S.viewChosen = true; }
    setFrame(frame, true);
  }
  /** What the mark button does for a cell: null where there is nothing to mark (no call to review, a control window, an analysis
   *  running, a review file that is not applied); else its words and, for a new mark, the call it gives. */
  function markAction(name) {
    const review = S.results && S.results.review;
    if (!review || review.problem || !name || !primaryShown() || !S.state || S.state.busy || S.closed) return null;
    if (markOf(name)) return { remove: true, label: 'Remove mark' };
    const automatic = automaticCall(name), choices = review.calls || [];
    if (!choices.includes(automatic)) return null;
    return { call: choices.find((c) => c !== automatic), automatic, label: automatic === 'Responded' ? 'Mark as not a responder' : 'Mark as a responder' };
  }
  /** The mark button, the open form and what Fiji answered, for the selected cell. Called on every state poll: it never rebuilds the form. */
  function showMark() {
    const name = S.cell, action = S.results ? markAction(name) : null, button = $('mark'), form = $('mark-form'), said = $('mark-said');
    if (S.marking && (!action || action.remove)) closeMark();
    button.hidden = !action || !!S.marking; if (action) button.textContent = action.label; button.disabled = S.reviewing;
    if (S.marking) { const other = $('mark-reason').value === 'Other'; $('mark-save').disabled = S.reviewing || !$('mark-reason').value; $('mark-note').required = other; $('mark-note').placeholder = other ? '' : 'Optional'; }
    // a refusal, under the cell it is for; else why the review file is not used, as Fiji gives it. Neither on a control window.
    const problem = S.results && S.results.review && S.results.review.problem, own = !!S.markSaid && S.markSaid.cell === name && primaryShown();
    const text = own ? S.markSaid.text : problem && name && primaryShown() ? problem : '';
    said.hidden = !text; said.textContent = text; said.className = 'mark-said ' + (own ? 'bad' : 'warn');
    form.hidden = !S.marking;
  }
  /** Opens the mark form in place, on the selected cell: the reason it starts with comes from the cell's check flag on a responder. */
  function openMark() {
    const name = S.cell, action = markAction(name); if (!action || action.remove) return;
    const review = S.results.review, reason = $('mark-reason'), note = $('mark-note');
    reason.replaceChildren(); const none = document.createElement('option'); none.value = ''; none.hidden = true; reason.appendChild(none);
    for (const text of review.reasons || []) { const o = document.createElement('option'); o.value = o.textContent = text; reason.appendChild(o); }
    reason.value = startingReason(name, action.automatic); note.value = ''; note.maxLength = review.noteMax || 120;
    S.marking = { cell: name, call: action.call }; S.markSaid = null; showMark(); Plot.dirty = true; drawPlot();
    reason.focus();
  }
  function closeMark() {
    if (!S.marking && $('mark-form').hidden) return;
    const had = document.activeElement && $('mark-form').contains(document.activeElement);
    S.marking = null; S.markSaid = null; $('mark-form').hidden = true; showMark(); Plot.dirty = true; drawPlot();
    if (had && !$('mark').hidden) $('mark').focus();
  }
  /** The note as Fiji keeps it: without the spaces at either end. */
  const tidy = (text) => String(text).replace(/^[\s\u001C-\u001F]+|[\s\u001C-\u001F]+$/g, '');
  /** What is wrong with a reason and note before they are sent, in Fiji's own words; '' when nothing is (Fiji has the last word). */
  function markProblem(reason, typed, review) {
    const note = tidy(typed), most = review.noteMax || 120;
    if (/[\x00-\x1F\x7F\x85\u2028\u2029]/.test(note)) return 'The note must be one line.';
    if (/\p{Cf}/u.test(note)) return 'The note must be plain text, without hidden characters.';
    if (note.length > most) return 'Keep the note to ' + most + ' characters or fewer.';
    if (note && '=+-@'.includes(note[0])) return 'The note cannot start with =, +, - or @.';
    if (reason === 'Other' && !note) return 'Other needs a short note.';
    return '';
  }
  async function saveMark() {
    const marking = S.marking; if (!marking || S.reviewing) return;
    const review = S.results.review, reason = $('mark-reason').value, note = $('mark-note').value; if (!reason) return;
    const problem = markProblem(reason, note, review); if (problem) { S.markSaid = { cell: marking.cell, text: problem }; showMark(); $('mark-note').focus(); return; }
    if (await sendReview({ folder: review.folder, cell: marking.cell, call: marking.call, reason, note })) $('mark').focus();
  }
  async function removeMark() {
    const name = S.cell, review = S.results && S.results.review; if (!name || !review || S.reviewing) return;
    if (await sendReview({ folder: review.folder, cell: name, remove: true })) $('mark').focus();
  }
  /** Sends a mark or its removal to Fiji. A refusal is shown as Fiji gives it, beside the form. */
  async function sendReview(body) {
    S.reviewing = true; showMark();
    try { applyReview(await api('review', body)); return true; }
    catch (problem) { S.markSaid = { cell: body.cell, text: /fetch/i.test(problem.message) ? 'Fiji is not connected.' : problem.message }; return false; }
    finally { S.reviewing = false; showMark(); }
  }
  /** Fiji's answer to a mark: its review and tally replace those of the results, and everything is drawn again in place, without
   *  loading the results anew (which would reset the plot's zoom, the below-zero switch, the table chosen and the notes). The
   *  version it answers with is the page's from now on, so that a poll from before the mark does not load the results again. */
  function applyReview(answer) {
    const r = S.results; if (!r || answer.folder !== r.folder) return;
    r.review = answer.review; r.tally = answer.tally;
    S.version = answer.version; if (S.state) S.state.version = answer.version;
    // what Fiji said is the status line while marks apply or after a run that marked automatically; else "Done." again (finishLine)
    S.said = marked() || automaticMade(answer.message) ? (answer.message || null) : null; S.marking = null; S.markSaid = null; $('mark-form').hidden = true;
    refreshReview(); showState();
  }
  /** Draws again, from the review and tally on hand, everything that counts, groups, colours or lists cells by call. Show stays
   *  on its group and the selected cell stays selected, also when it leaves that group: it stays listed until another cell is
   *  chosen (listCells). Only Marked, once no mark is left, gives way to All cells (updateWindow). */
  function refreshReview() {
    if (S.cell) S.keep = S.cell;
    showCalls(); showTally(); showMethods(); showMetricTable();
    updateWindow();
  }
  /** The Methods sentence of the review, to copy, with how many cells were marked; nothing without marks that changed a call.
   *  Called whenever the review may have changed (results loaded, a mark), so "Copied" from an earlier sentence goes. */
  function showMethods() {
    const review = S.results && S.results.review, sentence = (review && review.methods) || '', note = $('calls-saved');
    $('review-count').hidden = $('copy-methods').hidden = !sentence;
    if (sentence) { $('review-count').textContent = plural(review.changed || 0, 'cell') + ' marked'; $('copy-methods').title = sentence; }
    if (['Copied', 'Not copied'].includes(note.textContent)) { note.textContent = ''; note.className = 'save-note'; }
  }
  async function copyMethods() {
    const sentence = (S.results && S.results.review && S.results.review.methods) || '', note = $('calls-saved'); if (!sentence) return;
    let copied = false;
    try { await navigator.clipboard.writeText(sentence); copied = true; }
    catch (blocked) {   // a browser that keeps the clipboard from the page: the old way, through a selected box
      const box = document.createElement('textarea'); box.value = sentence; box.setAttribute('readonly', ''); box.style.position = 'fixed'; box.style.opacity = '0'; document.body.appendChild(box); box.select();
      try { copied = document.execCommand('copy'); } catch (failed) { copied = false; } box.remove();
    }
    note.className = 'save-note' + (copied ? '' : ' bad'); note.textContent = copied ? 'Copied' : 'Not copied';
  }
  /** While review marks apply: for each row of the calls table, the mark of its cell when the row is the primary stimulus's. */
  function reviewOfRows(data) {
    const marks = marked(); if (!marks || !data || !data.header) return null;
    const cell = data.header.indexOf('Cell'), primary = primaryRows(data);
    return data.rows.map((row) => cell >= 0 && primary(row) ? marks[row[cell]] || null : null);
  }
  /** Whether a row of a table of calls is the primary stimulus's, by its Role and Stimulus. */
  function primaryRows(data) {
    const role = data.header.indexOf('Role'), stimulus = data.header.indexOf('Stimulus'), w = primaryWindow();
    return (row) => !w || ((role < 0 ? '' : row[role]) || 'Stimulus') + '::' + ((stimulus < 0 ? '' : row[stimulus]) || 'Stimulus') === w.key;
  }
  /** A saved table with the codes of the call checks in its QC_flags, after the engine's own, each after the separator of that
   *  table: "; " in the calls and the response-window metrics, ";" in the per-cell tables, as the engine writes them. The checks
   *  are for the primary stimulus, so a row of a control window gets none. added: whether any code was put in. The saved file
   *  is not changed. */
  const CODED = ['calls', 'windowMetrics', 'quality', 'metrics'];
  function withCodes(kind, data) {
    const checks = S.results && S.results.checks, plain = { data, added: false };
    if (!CODED.includes(kind) || !checks || !checks.cells || !data || !data.header) return plain;
    const flags = data.header.indexOf('QC_flags'), cell = data.header.indexOf('Cell'); if (flags < 0 || cell < 0) return plain;
    const perCell = FAMILY[kind] === 'cell', join = perCell ? ';' : checks.join || '; ', primary = perCell ? () => true : primaryRows(data);
    let added = false;
    const rows = data.rows.map((row) => {
      const c = primary(row) ? checks.cells[row[cell]] : null, codes = c && Array.isArray(c.codes) ? c.codes.filter(Boolean) : []; if (!codes.length) return row;
      const own = row[flags] == null ? '' : String(row[flags]), out = fullRow(row, data.header.length); out[flags] = (own ? [own] : []).concat(codes).join(join); added = true; return out;
    });
    return added ? { data: Object.assign({}, data, { rows }), added } : plain;
  }
  /** What the record beside an Export says when the codes of the call checks were put into its QC_flags. */
  const codesRecord = (added) => added ? { qc_flags: 'QC_flags includes the codes of the call checks, after the engine\'s codes.' } : {};
  /** A saved row with every column of its header, for columns put after it. */
  const fullRow = (row, n) => { const own = row.slice(0, n); while (own.length < n) own.push(''); return own; };
  /** The plot opens on the average alone and stays on what was chosen: picking a cell never changes what the plot shows.
   *  A cell's own trace appears only under "One cell + average". */
  function setPlotMode(mode) {
    // "Around the cell" is in gray levels: the Y axis typed in the figure editor for the other plots is kept for them, and its
    // own for it, each put back when its kind of plot is shown again
    const from = yKind(S.plotMode), to = yKind(mode);
    if (from !== to) { Fig.yKept[from] = Object.fromEntries(Y_TYPED.map((key) => [key, Fig.plot[key]])); const back = Fig.yKept[to] || {}; for (const key of Y_TYPED) Fig.plot[key] = back[key] !== undefined ? back[key] : STYLE.plot[key]; }
    S.plotMode = mode; Plot.hover = null;
    for (const button of document.querySelectorAll('[data-plot]')) { const on = button.dataset.plot === mode; button.classList.toggle('on', on); button.setAttribute('aria-pressed', on ? 'true' : 'false'); }
    showY(); Plot.dirty = true; drawPlot();
  }
  /** The Y axis boxes under the plot show the limits typed for what the plot shows (yKey). */
  const showY = () => { const now = Plot.y[yKey()] || {}; $('y-min').value = finite(now.lo) ? now.lo : ''; $('y-max').value = finite(now.hi) ? now.hi : ''; };
  /** A cell chosen by the person (list, arrows, a click on the picture, a trace, a table row or the heatmap), or none. The plot
   *  follows: with a cell it shows that cell over the average, with none the average alone. "Responders vs non-responders" stays,
   *  and so does "Around the cell" while a cell is selected. */
  function chooseCell(name) {
    pick(name || null);
    if (S.results && S.cell && S.plotMode === 'mean') setPlotMode('cell'); else if (S.results && !S.cell && (S.plotMode === 'cell' || S.plotMode === 'around')) setPlotMode('mean');
  }
  /** The cell meant by what was typed in the Cell box: its shown name ("Cell 12"), its saved name, or its number alone.
   *  null for an empty box (no cell); undefined when nothing in the Show list matches. */
  function cellTyped(text, numberToo) {
    const typed = String(text || '').trim().toLowerCase(); if (!typed) return null;
    let hit = S.listed.find((n) => shortName(n).toLowerCase() === typed || n.toLowerCase() === typed);
    if (!hit && numberToo && /^\d+$/.test(typed)) hit = S.listed.find((n) => cellNumber(n) === String(Number(typed)));
    return hit;
  }
  function stepCell(by) { if (!S.listed.length) return; const i = S.listed.indexOf(S.cell), next = i < 0 ? (by > 0 ? 0 : S.listed.length - 1) : ((i + by) % S.listed.length + S.listed.length) % S.listed.length; chooseCell(S.listed[next]); }

  // ---------------------------------------------------------------- one acquisition-frame cursor for the recording and saved traces
  function frameTimes() {
    if (S.results) return S.results.frameTimes || S.results.time || [];
    // Null is an explicit invalid/unknown draft time base; do not resurrect the initial FPS.
    if (S.state) return S.state.frameTimes || [];
    return (S.info && S.info.frameTimes) || [];
  }
  function nearestFrame(t) { const time = frameTimes(); let best = 0, gap = Infinity; for (let i = 0; i < time.length; i++) { const d = Math.abs(time[i] - t); if (d < gap) { gap = d; best = i; } } return best + 1; }
  function updateTransport() {
    if (!S.info) return;
    const available = S.info.frames > 1 && (S.info.hasMovie || !!S.results), staticView = still();
    $('time').disabled = !available || (staticView && !S.results);
    $('play').disabled = !available || staticView;
    $('frame-back').disabled = $('time').disabled || S.frame <= 1;
    $('frame-forward').disabled = $('time').disabled || S.frame >= S.info.frames;
    $('speed').disabled = $('play').disabled;
    const t = frameTimes()[S.frame - 1], provisional = !S.results && (!S.state || S.state.timingProvisional !== false);
    const label = (finite(t) ? (provisional ? '≈ ' : '') + t.toFixed(2) + ' s · ' : '') + 'Frame ' + S.frame + ' / ' + S.info.frames;
    $('time-label').textContent = label; $('time').setAttribute('aria-valuetext', label); frameButtons();
    $('playback-note').textContent = S.view === 'background' ? 'Background areas: a still picture. Blue marks the pixels whose mean is taken off each cell.' : staticView ? 'Static colour overview — this image does not move. Choose “Recording” in View to play or scrub the movie.' : !S.info.hasMovie ? 'The recording is unavailable in this session. Reopen it in Fiji to view its frames.' : frameTimes().length ? (provisional ? 'Preview timing is provisional. Playback follows the displayed frame times; confirm timing before analysis.' : 'Playback follows the applied frame times. Drag the slider to pause and inspect a frame.') : 'Frame-by-frame preview. No valid time base is available; 1× plays 20 frames per second.';
    $('playback-note').classList.toggle('help', !staticView && S.info.hasMovie && frameTimes().length > 0);
    $('playback-note').classList.toggle('help', /^Playback follows the applied frame times/.test($('playback-note').textContent));
  }
  function setFrame(wanted, manual) {
    if (!S.info) return;
    if (manual) play(false);
    const before = S.frame; S.frame = Math.max(1, Math.min(Math.max(1, S.info.frames), Math.round(wanted))); $('time').value = S.frame;
    if (S.frame !== before) { Plot.movedAt = performance.now(); clearTimeout(Plot.rest); Plot.rest = setTimeout(() => { if (S.results) drawPlot(); }, LINGER + 40); }   // then the plot is redrawn without the line
    updateTransport(); drawPlot(); requestFrame();
  }
  function seek(t) { setFrame(nearestFrame(t), true); }
  /** Frame times for the clock: the displayed time base, or 20 frames per second at 1× when there is none. */
  function clockTimes() {
    const n = S.info.frames, t = frameTimes(), real = t.length >= n && t.slice(0, n).every(finite) && t.slice(1, n).every((v, i) => v > t[i]);
    Clock.real = real; return real ? t.slice(0, n) : Array.from({ length: n }, (_, i) => i / 20);
  }
  function play(on) {
    if (S.timer) cancelAnimationFrame(S.timer); S.timer = null; Plot.playhead = null; Clock.waitingSince = 0;
    if (on && S.info && S.info.hasMovie && !still() && S.info.frames > 1) {
      Buf.failed.clear();   // pressing Play retries frames that failed before
      if (S.frame >= S.info.frames) setFrame(1);
      const times = clockTimes(), gaps = times.slice(1).map((v, i) => v - times[i]).sort((a, b) => a - b), typical = gaps[Math.floor(gaps.length / 2)] || 0.05;
      // above 30 frames per second, every Nth frame is shown instead of slowing down
      Clock.times = times; Clock.step = Math.max(1, Math.ceil(Number($('speed').value) / typical / MAX_FPS));
      Clock.t = times[S.frame - 1]; Clock.last = performance.now(); S.timer = requestAnimationFrame(tick);
    }
    $('play').textContent = S.timer ? '❚❚ Pause' : '▶ Play'; $('play').setAttribute('aria-pressed', String(!!S.timer));
    if (S.info && S.info.hasMovie) { frameNote(); pumpFrames(); } if (S.results) drawPlot();
  }
  function nextFrame(f) { return f >= S.info.frames ? null : Math.min(S.info.frames, f + Clock.step); }
  function tick(now) {
    if (!S.timer) return;
    S.timer = requestAnimationFrame(tick);
    const times = Clock.times, speed = Number($('speed').value), n = S.info.frames;
    Clock.t += Math.min(0.1, Math.max(0, (now - Clock.last) / 1000)) * speed; Clock.last = now;   // a hidden tab does not jump ahead
    let target = S.frame, waiting = false;
    for (let k = nextFrame(S.frame); k && times[k - 1] <= Clock.t; k = nextFrame(k)) {
      if (bitmapFor(k, S.view)) { target = k; continue; }
      if (Buf.failed.has(frameKey(k, S.view))) { frameFailed(k, S.view); return; }
      Clock.t = times[k - 1]; waiting = true; break;   // the playhead never runs ahead of the picture
    }
    if (waiting && !Clock.waitingSince) Clock.waitingSince = now; else if (!waiting && Clock.waitingSince) { Clock.waitingSince = 0; frameNote(); }
    if (target !== S.frame) setFrame(target);
    else if (S.frame === n) { play(false); return; }
    else if (waiting) { if (now - Clock.waitingSince > 250) frameNote(); }
    const ahead = nextFrame(S.frame);   // a long gap between recorded frames waits at most about 1 s
    if (ahead && !waiting && times[ahead - 1] - Clock.t > speed) Clock.t = times[ahead - 1] - speed;
    if (S.results && Clock.real) { Plot.playhead = Math.min(Clock.t, times[n - 1]); drawPlot(); }
    pumpFrames();
  }
  function movieMessage(message, error) { const box = $('movie-message'); box.textContent = message; box.hidden = !message; box.classList.toggle('error', !!error); }
  // ---------------------------------------------------------------- movie loading: buffered, latest requested frame wins, failures stay visible
  const frameKey = (f, mode) => Buf.version + '|' + mode + '|' + f;
  function resetFrames(version) {
    for (const entry of Buf.bitmaps.values()) if (entry.bitmap !== Pic.shown && entry.bitmap.close) entry.bitmap.close();
    Buf.version = version; Buf.blobs.clear(); Buf.bitmaps.clear(); Buf.failed.clear(); Buf.bytes = 0;
    if (S.info && S.info.hasMovie && !still()) requestFrame();
  }
  /** Frames are fetched at the size they are drawn, or at full resolution when zoomed in or when that is nearly as large. */
  function frameWidth() {
    const full = S.info.width, canvas = $('picture');
    if (Pic.zoom > 1.01 || !canvas.clientWidth) return full;
    const need = full * place(canvas.clientWidth, canvas.clientHeight).z * (window.devicePixelRatio || 1);
    return need >= full * 0.75 ? full : Math.min(full, Math.max(160, Math.ceil(need / 160) * 160));
  }
  function bitmapFor(f, mode) { const entry = Buf.bitmaps.get(frameKey(f, mode)); return entry && entry.w >= frameWidth() ? entry.bitmap : null; }
  function fetchFrame(f, mode, w) {
    const key = frameKey(f, mode), busy = Buf.inflight.get(key); if (busy && busy.w >= w) return busy.promise;
    const job = { w };
    job.promise = (async () => {
      const response = await fetch('api/frame?f=' + f + '&mode=' + mode + '&w=' + w + '&v=' + encodeURIComponent(Buf.version));
      if (response.status === 410) Buf.closed = (await response.json().catch(() => ({}))).error || '';   // the recording was closed in Fiji
      if (!response.ok) throw new Error('The recording frame could not be loaded.');
      Buf.closed = '';
      const blob = await response.blob();
      if (key.startsWith(Buf.version + '|')) {
        const old = Buf.blobs.get(key); if (old) Buf.bytes -= old.blob.size;
        Buf.blobs.set(key, { f, w, blob }); Buf.bytes += blob.size; Buf.failed.delete(key); trimFrames();
      }
      return blob;
    })();
    job.promise.catch(() => { Buf.failed.set(key, performance.now()); }).finally(() => { if (Buf.inflight.get(key) === job) Buf.inflight.delete(key); pumpFrames(); });
    Buf.inflight.set(key, job); return job.promise;
  }
  function decodeFrame(f, mode) {
    const key = frameKey(f, mode), w = frameWidth(), have = Buf.bitmaps.get(key); if (have && have.w >= w) return Promise.resolve(have.bitmap);
    const busy = Buf.decoding.get(key); if (busy && busy.w >= w) return busy.promise;
    const stored = Buf.blobs.get(key), job = { w: stored && stored.w >= w ? stored.w : w };
    job.promise = (stored && stored.w >= w ? Promise.resolve(stored.blob) : fetchFrame(f, mode, w)).then((blob) => createImageBitmap(blob)).then((bitmap) => {
      const old = Buf.bitmaps.get(key); if (old && old.bitmap !== Pic.shown && old.bitmap.close) old.bitmap.close();
      if (key.startsWith(Buf.version + '|')) Buf.bitmaps.set(key, { f, w: job.w, bitmap }); return bitmap;
    });
    job.promise.catch(() => {}).finally(() => { if (Buf.decoding.get(key) === job) Buf.decoding.delete(key); });
    Buf.decoding.set(key, job); return job.promise;
  }
  /** Keeps the buffer filled: the frame on screen first, then about 3 s of playback ahead (a few frames when paused). */
  // A frame that just failed is left alone for 2 s, the one on screen included: every finished request calls this again, so without the wait a
  // frame that cannot load would be asked for without pause. Moving the slider or pressing Play asks for it again at once.
  function pumpFrames() {
    if (!S.info || !S.info.hasMovie || still()) return;
    const mode = S.view, w = frameWidth(), speed = Number($('speed').value), wanted = [S.frame];
    const ahead = S.timer ? Math.min(120, Math.ceil(3 * speed / Math.max(1e-6, (Clock.times[Math.min(Clock.times.length - 1, S.frame)] - Clock.times[S.frame - 1]) || 0.05) / Clock.step)) : 3;
    for (let k = nextFrame(S.frame), i = 0; k && i < ahead; k = nextFrame(k), i++) wanted.push(k);
    if (!S.timer && S.frame > 1) wanted.push(S.frame - 1);
    wanted.forEach((f, i) => {
      const key = frameKey(f, mode), stored = Buf.blobs.get(key), failed = Buf.failed.get(key);
      if ((!stored || stored.w < w) && Buf.inflight.size < PARALLEL && !Buf.inflight.has(key) && !(failed && performance.now() - failed < 2000)) fetchFrame(f, mode, w).catch(() => {});
      else if (stored && stored.w >= w && i < 8 && !bitmapFor(f, mode)) decodeFrame(f, mode).catch(() => {});
    });
    trimFrames();
  }
  /** Decoded frames are kept near the playhead only; compressed ones until the memory budget is reached, farthest first. */
  function trimFrames() {
    for (const [key, entry] of Buf.bitmaps) if (Math.abs(entry.f - S.frame) > DECODED_AROUND && entry.bitmap !== Pic.shown) { if (entry.bitmap.close) entry.bitmap.close(); Buf.bitmaps.delete(key); }
    if (Buf.bytes <= MAX_BYTES) return;
    const far = Array.from(Buf.blobs.entries()).sort((a, b) => Math.abs(b[1].f - S.frame) - Math.abs(a[1].f - S.frame));
    for (const [key, entry] of far) { if (Buf.bytes <= MAX_BYTES * 0.8) break; Buf.blobs.delete(key); Buf.bytes -= entry.blob.size; }
  }
  /** What the colours of the view on screen mean: [colour, words]. "As recorded" is plain gray values and has none. */
  function viewColours() {
    const rest = S.results ? S.results.baselineFrames : (S.state && S.state.rest), frames = rest ? ' (frames ' + rest[0] + '–' + rest[1] + ')' : '';
    return S.view === 'background' ? [['#2896ff', 'pixels taken as the cells’ background']] : S.view === 'picture' ? [['#00e000', 'cells'], ['#ffe600', 'very bright (precipitate)'], ['#ff40ff', 'response']]
      : S.view === 'magenta' ? [['#00e000', 'cells, background removed'], ['#ff40ff', 'rise above baseline' + frames]]
      : S.view === 'change' ? [['#ffffff', 'rise above baseline' + frames]] : S.view === 'flat' ? [['#ffffff', 'background removed']] : S.view === 'green' ? [['#00e000', 'cells, background removed']] : [];
  }
  /** What the outlines on the picture mean, as a key: [colour, words, true, dashed]. By call: the three calls. Plain, or before
   *  an analysis: "your cells", and the dashed outline of an excluded or unscorable cell when one is drawn. None when the
   *  outlines are off. The one place for these words: the chips under the picture, the key of a saved picture, and the file
   *  kept beside it. */
  function outlineKey() {
    if (!S.outlines) return [];
    if (S.results && S.outlineColours === 'call') return [[ON_PICTURE.responded, 'responded', true, false], [ON_PICTURE.not, 'did not respond', true, false], [ON_PICTURE.other, 'excluded or unscorable', true, true]];
    const dashes = S.info.names.some((name, i) => { const look = outlineLook(i); return !!look && look.dashed; });
    return [['#ffffff', 'your cells', true, false]].concat(dashes ? [['#ffffff', 'excluded or unscorable', true, true]] : []);
  }
  /** Under the picture: what the colours of the view on screen mean (filled swatches) and, once there are results, what the
   *  outlines mean (lines, dashed where the outlines are). Each row is rebuilt only when it changes. */
  function pictureChips() {
    const fill = (box, entries) => { const shown = JSON.stringify(entries); if (box.dataset.shown === shown) return; box.dataset.shown = shown; box.replaceChildren();
      for (const [colour, words, line, dashed] of entries) { const chip = document.createElement('span'), mark = document.createElement('i'); chip.className = 'chip';
        if (dashed) { mark.className = 'dash'; mark.style.borderTopColor = colour; } else { if (!line) mark.className = 'fill'; mark.style.background = colour; }
        chip.append(mark, words); box.appendChild(chip); } };
    fill($('view-key'), viewColours()); fill($('call-legend'), S.results ? outlineKey() : []); $('call-legend').hidden = !$('call-legend').firstChild;
  }
  function frameNote() {
    if (still() || !S.info.hasMovie || !Pic.shownFrame || Pic.error) return;
    const f = Pic.shownFrame, rest = S.results ? S.results.baselineFrames : (S.state && S.state.rest), ordinal = (n) => n + (n % 10 === 2 && n % 100 !== 12 ? 'nd' : n % 10 === 3 && n % 100 !== 13 ? 'rd' : n % 10 === 1 && n % 100 !== 11 ? 'st' : 'th');
    const frames = rest ? rest[0] + '–' + rest[1] : '(preview)';
    $('picture-note').textContent = 'Frame ' + f + (isMoved(f) ? ' · the field is displaced in this frame (' + movedText() + ' moved)' : '') + (S.timer && Clock.step > 1 ? ' · playing every ' + ordinal(Clock.step) + ' frame' : '') + (S.timer && Clock.waitingSince && performance.now() - Clock.waitingSince > 250 ? ' · buffering' : '');
  }
  function showFrame(f, mode, bitmap) {
    if (S.view !== mode || S.frame !== f) return;
    Pic.shown = bitmap; Pic.shownFrame = f; Pic.error = null; movieMessage(''); frameNote(); drawPicture();
  }
  function frameFailed(f, mode) {
    if (S.view !== mode || S.frame !== f && !S.timer) return;
    Pic.error = 'The recording frame could not be loaded.'; play(false);
    movieMessage(Buf.closed || 'Frame ' + f + ' could not load. Keep the recording open in Fiji, then move the slider to retry.', true);
    $('picture-note').textContent = Pic.shownFrame ? 'Last displayed frame: ' + Pic.shownFrame : 'No movie frame loaded';
  }
  function requestFrame() {
    if (!(S.info.hasPicture || S.info.hasMovie)) return;
    if (still() || !S.info.hasMovie) {
      const image = S.view === 'background' ? Pic.background : Pic.base; Pic.shown = image || Pic.base; Pic.shownFrame = null;
      $('picture-note').textContent = ''; movieMessage(image ? '' : 'Loading the picture…'); drawPicture(); return;
    }
    const f = S.frame, mode = S.view, ready = bitmapFor(f, mode);
    if (ready) showFrame(f, mode, ready);
    else {
      const any = Buf.bitmaps.get(frameKey(f, mode));   // a smaller copy stands in while the full-size one loads
      if (any) showFrame(f, mode, any.bitmap); else if (!S.timer && (Pic.shownFrame !== f || Pic.error)) movieMessage('Loading frame ' + f + '…');
      if (Pic.error) Buf.failed.delete(frameKey(f, mode));
      decodeFrame(f, mode).then((bitmap) => showFrame(f, mode, bitmap), () => frameFailed(f, mode));
    }
    pumpFrames();
  }
  // ---------------------------------------------------------------- the picture: its view, and the cells drawn on it
  /** The view of the picture, in one place. Pic.zoom says how far it is zoomed in: 1 is the whole picture, fitted into the
   *  room beside the controls (viewRoom); ZOOM_MOST is the most. Pic.cx and Pic.cy are the point of the recording, in its
   *  own pixels, wanted in the middle of that room. Every way of zooming and moving goes through setView: the buttons and
   *  the keys (in the steps of ZOOM_STEPS; ZOOM_KEYS says which key steps which way, 0 for the whole picture), the wheel, a
   *  pinch, a drag, a rectangle, the small map, "Zoom to cell" and the cell that is picked. */
  const ZOOM_MOST = 16, ZOOM_STEPS = [1, 1.5, 2, 3, 4, 6, 8, 12, 16], ZOOM_KEYS = { '+': 1, '=': 1, '-': -1, '_': -1, '0': 0 };
  /** The room the picture is shown in: the panel less a strip TOOLS pixels high along the top, where the controls lie in one
   *  row in every view and at every zoom. The picture is drawn in the room only (drawPicture), so nothing of it, no cell
   *  and no outline, ever lies under the controls. */
  const TOOLS = 36;
  /** Where the picture lies in a panel of w × h: its left top corner (x, y) and the pixels of the screen to one of its own (z);
   *  with these, the top of the room it is shown in (roomTop), the middle of that room (midX, midY) and its size (wide, high).
   *  Along each side: while the picture fits into the room, it sits in the middle of it. Larger, the point asked for is in
   *  the middle as far as this limit allows: an edge of the picture comes into the room by no more than the picture is
   *  larger than the room, and never past the middle. So the picture cannot be dragged out of view, yet a cell at its edge
   *  can be brought clear of the controls and of the map. */
  function place(w, h) {
    const iw = S.info.width, ih = S.info.height, wide = w, high = Math.max(1, h - TOOLS), fitted = Math.max(0.01, Math.min((wide - 8) / iw, (high - 8) / ih)), z = fitted * Pic.zoom;
    const along = (from, free, size, centre) => { if (size <= free - 8 + 1e-6) return from + (free - size) / 2;
      const slack = Math.min(free / 2 - 4, size - free + 8); return Math.max(from + free - 4 - size - slack, Math.min(from + 4 + slack, from + free / 2 - centre * z)); };
    return { x: along(0, wide, iw * z, Pic.cx == null ? iw / 2 : Pic.cx), y: along(TOOLS, high, ih * z, Pic.cy == null ? ih / 2 : Pic.cy), z, fitted, roomTop: TOOLS, midX: wide / 2, midY: TOOLS + high / 2, wide, high };
  }
  /** Sets the view and draws it. The zoom is kept from 1 to ZOOM_MOST and the point in the middle where place() allows it,
   *  and that point is remembered as it came out, so a drag against an edge turns back at once. Zoomed in, the frame on screen
   *  is asked for at its full size. Without a point, the middle of the picture. */
  function setView(zoom, cx, cy) {
    if (!S.info) return;
    const canvas = $('picture'), w = canvas.clientWidth, h = canvas.clientHeight, movie = S.info.hasMovie && !still(), had = movie ? frameWidth() : 0;
    Pic.zoom = Math.max(1, Math.min(ZOOM_MOST, zoom)); Pic.cx = cx == null ? S.info.width / 2 : cx; Pic.cy = cy == null ? S.info.height / 2 : cy;
    if (w && h) { const at = place(w, h); Pic.cx = (at.midX - at.x) / at.z; Pic.cy = (at.midY - at.y) / at.z; }
    drawPicture(); if (movie && frameWidth() > had) requestFrame();
  }
  /** Zooms while the point of the picture at (px, py) of the panel stays where it is, as far as place() allows: the point under the pointer. */
  function zoomAt(zoom, px, py) {
    if (!S.info) return;
    const canvas = $('picture'), before = place(canvas.clientWidth, canvas.clientHeight), z = before.fitted * Math.max(1, Math.min(ZOOM_MOST, zoom));
    setView(zoom, (px - before.x) / before.z + (before.midX - px) / z, (py - before.y) / before.z + (before.midY - py) / z);
  }
  /** One step in (by > 0) or out on ZOOM_STEPS, from whatever the zoom is now; about (px, py), or the middle of the view. */
  function zoomStep(by, px, py) {
    if (!S.info) return;
    const canvas = $('picture'), at = place(canvas.clientWidth, canvas.clientHeight), next = by > 0 ? ZOOM_STEPS.find((z) => z > Pic.zoom * 1.02) : ZOOM_STEPS.slice().reverse().find((z) => z < Pic.zoom / 1.02);
    zoomAt(next || (by > 0 ? ZOOM_MOST : 1), px == null ? at.midX : px, py == null ? at.midY : py);
  }
  /** The part of the recording that is in view, under the strip of the controls, in the recording's pixels and never beyond the picture: { x, y, width, height }. */
  function viewArea() {
    const canvas = $('picture'), w = canvas.clientWidth, h = canvas.clientHeight, at = place(w, h), x = Math.max(0, -at.x / at.z), y = Math.max(0, (at.roomTop - at.y) / at.z);
    return { x, y, width: Math.max(0, Math.min(S.info.width, (w - at.x) / at.z) - x), height: Math.max(0, Math.min(S.info.height, (h - at.y) / at.z) - y) };
  }
  /** The part of the recording a saved picture or frame shows. null: the whole picture ("Area" says so, or all of it is on
   *  screen). With "Area on screen": the part that is in the panel (viewArea) in whole pixels of the recording, every pixel
   *  counted of which any part shows, and the zoom it is seen at. */
  function savedArea() {
    if (S.saveArea !== 'screen') return null;
    const part = viewArea(), edge = 1e-6, x = Math.floor(part.x + edge), y = Math.floor(part.y + edge);
    const width = Math.min(S.info.width, Math.ceil(part.x + part.width - edge)) - x, height = Math.min(S.info.height, Math.ceil(part.y + part.height - edge)) - y;
    return width < 1 || height < 1 || (width === S.info.width && height === S.info.height) ? null : { x, y, width, height, zoom: Pic.zoom };
  }
  /** "Area", in the picture's Export menu and in Save frames: the whole picture ('whole') or the area on screen ('screen'). */
  function setSaveArea(value) { S.saveArea = value === 'screen' ? 'screen' : 'whole'; for (const id of ['save-area', 'frames-area']) $(id).value = S.saveArea; }
  /** Zooms to a rectangle drawn on the picture ([x0, y0, x1, y1] in the panel's pixels): the part of the recording under it
   *  then fills the room beside the controls as the whole picture does at 1×. A rectangle too small to be meant, or off the picture, does
   *  nothing; says whether it zoomed. */
  function zoomToBand(band) {
    const canvas = $('picture'), at = place(canvas.clientWidth, canvas.clientHeight), within = (v, most) => Math.max(0, Math.min(most, v));
    const x0 = within((Math.min(band[0], band[2]) - at.x) / at.z, S.info.width), x1 = within((Math.max(band[0], band[2]) - at.x) / at.z, S.info.width);
    const y0 = within((Math.min(band[1], band[3]) - at.y) / at.z, S.info.height), y1 = within((Math.max(band[1], band[3]) - at.y) / at.z, S.info.height);
    if ((x1 - x0) * at.z < 6 || (y1 - y0) * at.z < 6) return false;
    setView(Math.min((at.wide - 8) / (x1 - x0), (at.high - 8) / (y1 - y0)) / at.fitted, (x0 + x1) / 2, (y0 + y1) / 2); return true;
  }
  /** "Zoom to region" on or off. While it is on, a drag on the picture draws a rectangle instead of moving the picture. */
  function setRegion(on) {
    Pic.region = on; if (!on) Pic.band = null;
    const button = $('zoom-region'); button.classList.toggle('on', on); button.setAttribute('aria-pressed', on ? 'true' : 'false'); drawPicture();
  }
  /** What the controls on the picture say about its view: the zoom level, which buttons can do anything, the cursor (a hand
   *  where a drag moves the picture, a cross where it draws a rectangle) and the small map; and, where a picture is saved,
   *  that "Area" now makes a difference. */
  function showView(image, w, h) {
    const zoomed = Pic.zoom > 1, text = (Pic.zoom >= 10 ? Math.round(Pic.zoom) : +Pic.zoom.toFixed(1)) + '×', stage = $('picture-stage');
    if ($('zoom-level').textContent !== text) $('zoom-level').textContent = text;
    $('zoom-out').disabled = $('zoom-all').disabled = !zoomed; $('zoom-in').disabled = Pic.zoom >= ZOOM_MOST;
    stage.classList.toggle('zoomed', zoomed); stage.classList.toggle('region', Pic.region || !!Pic.band); stage.classList.toggle('moving', zoomed && !!Pic.drag && Pic.moved && !Pic.drag.band);
    for (const id of ['save-area', 'frames-area']) $(id).closest('label').classList.toggle('matters', zoomed);   // zoomed in, "Area" decides what is saved: its row stands out
    drawMap(image, w, h);
  }
  /** The small map, shown while the picture is zoomed in: the whole field as the panel shows it now, dimmed outside a frame
   *  around the part that is in view. A click or a drag in it moves the view there (wire). Its longer side is about a
   *  quarter of the panel's shorter one. It lies on the bottom right corner, faint until the pointer is on it, and on the
   *  bottom left while the selected cell would lie under it. Pic.map keeps the frame as drawn, in the map's own pixels. */
  function drawMap(image, w, h) {
    const map = $('view-map'), shown = Pic.zoom > 1 && !!image; map.hidden = !shown; Pic.map = null; if (!shown) return;
    const iw = S.info.width, ih = S.info.height, long = Math.round(Math.max(64, Math.min(128, 0.26 * Math.min(w, h))));
    const mw = iw >= ih ? long : Math.max(1, Math.round(long * iw / ih)), mh = iw >= ih ? Math.max(1, Math.round(long * ih / iw)) : long;
    if (map.style.width !== mw + 'px') map.style.width = mw + 'px'; if (map.style.height !== mh + 'px') map.style.height = mh + 'px';
    const selected = S.cell ? S.info.names.indexOf(S.cell) : -1, at = place(w, h), b = selected >= 0 ? S.bounds[selected] : null, gap = 8 + 4;
    map.classList.toggle('left', !!b && at.x + b.x1 * at.z > w - gap - mw && at.y + b.y1 * at.z > h - gap - mh && at.x + b.x0 * at.z < w && at.y + b.y0 * at.z < h);
    const { ctx } = fit(map), area = viewArea(), x = area.x * mw / iw, y = area.y * mh / ih, fw = area.width * mw / iw, fh = area.height * mh / ih;
    ctx.imageSmoothingEnabled = true; ctx.drawImage(image, 0, 0, mw, mh);
    ctx.fillStyle = 'rgba(16,22,28,0.55)'; ctx.fillRect(0, 0, mw, y); ctx.fillRect(0, y + fh, mw, mh - y - fh); ctx.fillRect(0, y, x, fh); ctx.fillRect(x + fw, y, mw - x - fw, fh);
    ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 1.5; ctx.strokeRect(x + 0.75, y + 0.75, Math.max(1, fw - 1.5), Math.max(1, fh - 1.5));
    Pic.map = { x, y, width: fw, height: fh };
  }
  /** Zoomed in this far, the recording's pixels are drawn as squares on screen, not smoothed. (A saved area has its own rule: paintView.) */
  const pixelsShown = (zoom) => zoom > 3;
  /** The pens the cells are drawn with, in pixels: on screen, where a pixel is one of the screen at any zoom; and in a saved
   *  picture, where it is one of a recording 1388 wide and grows with a wider one. For each kind of outline, its dark halo
   *  and its line (plain: a cell without a call, or any cell under "Plain"); the dashes of a cell that is excluded or
   *  unscorable; and how high a cell's number is at the least and at the most. */
  const PENS = { screen: { halo: { plain: 2.2, responded: 2.8, not: 2.2, other: 2.2 }, line: { plain: 1.3, responded: 1.5, not: 1.1, other: 1.1 }, dash: [4, 3], number: [9, 13] },
    saved: { halo: { plain: 3.2, responded: 3.2, not: 3.2, other: 3.2 }, line: { plain: 1.8, responded: 1.8, not: 1.8, other: 1.8 }, dash: [5, 4], number: [12, 22] } };
  /** How a cell's outline is drawn, in one place, for the picture on screen and for a saved one. null: it is not drawn ("All
   *  outlines" is off, or the cell is not in the group under Show). In the colour of its call; or white, as before an
   *  analysis, when "Outlines" says Plain. A cell that is excluded or unscorable keeps its dashes either way. */
  function outlineLook(i) {
    const name = S.info.names[i], known = !!S.results, call = callOf(name); if (!S.outlines || (known && !inGroup(call, S.group, name))) return null;
    const which = kind(call), plain = !known || S.outlineColours === 'plain';
    return { pen: plain ? 'plain' : which, colour: plain ? '#ffffff' : ON_PICTURE[which], dashed: known && which === 'other' };
  }
  /** A cell's number, as the page counts cells: the number at the end of its shown name ("Cell 12": 12; "Cell007": 7), which
   *  also opens the cell when typed into the Cell box. A name that ends in no number is written whole. */
  const cellNumber = (name) => { const m = shortName(name).match(/(\d+)\s*$/); return m ? String(Number(m[1])) : shortName(name); };
  /** Draws the cells on a picture: every outline that is shown and, with "Number cells", its number, white on a dark halo, a
   *  little over half as high as a typical cell is wide (within the pen's limits), in the cell or beside it (placeLabels).
   *  One drawing for the screen and for a saved picture. ctx is in the recording's pixels, with round line joins; unit is
   *  the length there of one pixel of the pen (PENS); "except" is a cell left out (the selected one, which the screen draws
   *  in its own way); "part" is the part of the recording the picture shows ({ x, y, width, height }; all of it when not
   *  given): a cell whose middle lies in it is numbered, and the numbers stay inside it. Leaves no dashes set. */
  function paintCells(ctx, pen, unit, except, part) {
    const drawn = [], box = part || { x: 0, y: 0, width: S.info.width, height: S.info.height }, within = (c) => c[0] >= box.x && c[0] <= box.x + box.width && c[1] >= box.y && c[1] <= box.y + box.height;
    S.paths.forEach((path, i) => { const look = i === except ? null : outlineLook(i); if (!look) return; drawn.push(i);
      ctx.setLineDash([]); ctx.strokeStyle = 'rgba(0,0,0,0.55)'; ctx.lineWidth = pen.halo[look.pen] * unit * S.outlineWidth; ctx.stroke(path);
      ctx.setLineDash(look.dashed ? pen.dash.map((length) => length * unit) : []); ctx.strokeStyle = look.colour; ctx.lineWidth = pen.line[look.pen] * unit * S.outlineWidth; ctx.stroke(path); });
    ctx.setLineDash([]);
    if (!S.numberCells || !drawn.length) return;
    const high = Math.max(pen.number[0], Math.min(pen.number[1], 0.55 * S.cellSize / unit)) * unit;
    ctx.save(); ctx.font = '600 ' + high + 'px ' + getComputedStyle(document.body).fontFamily; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.strokeStyle = 'rgba(0,0,0,0.8)'; ctx.lineWidth = 0.3 * high; ctx.fillStyle = '#ffffff';
    const reach = (i) => { const b = S.bounds[i]; return b ? [(b.x1 - b.x0) / 2, (b.y1 - b.y0) / 2] : [0, 0]; };
    const labels = drawn.filter((i) => S.centres[i].every(finite) && within(S.centres[i])).map((i) => { const text = cellNumber(S.info.names[i]); return { text, x: S.centres[i][0], y: S.centres[i][1], wide: ctx.measureText(text).width + 0.3 * high, home: S.centres[i], reach: reach(i) }; });
    placeLabels(labels, high, box);
    for (const label of labels) { ctx.strokeText(label.text, label.x, label.y); ctx.fillText(label.text, label.x, label.y); }
    ctx.restore();
  }
  /** A cell is large enough to show its number in its middle, with its outline still to be seen around it. */
  const LABEL_ROOM = 1.25;
  const holdsLabel = (label, high) => 2 * label.reach[0] >= LABEL_ROOM * label.wide && 2 * label.reach[1] >= LABEL_ROOM * high;
  /** Places the numbers of the cells so that none lies on another and each stays inside the picture ("box": the part of the
   *  recording it shows). Each label: its cell's middle (home), the half width and half height of its outline (reach), its
   *  own width; "high" is the height of them all. A number goes in the middle of a cell that holds it (holdsLabel); for a
   *  smaller cell, or where that place is taken, beside the cell: right, left, above, below, then a corner. Of these it
   *  takes the first that lies on no number placed before it and on no other cell; failing that, the first on no number;
   *  failing that, the one that overlaps numbers least. Cells that hold their number go first, and keep its middle while no
   *  number is there. One rule on screen and saved. */
  function placeLabels(labels, high, box) {
    const gap = 0.15 * high, placed = [];
    const inside = (x, y, label) => [Math.max(box.x + label.wide / 2, Math.min(box.x + box.width - label.wide / 2, x)), Math.max(box.y + high / 2, Math.min(box.y + box.height - high / 2, y))];
    const over = (x, y, label, other, wide, tall) => Math.max(0, (label.wide + wide) / 2 - Math.abs(x - other[0])) * Math.max(0, (high + tall) / 2 - Math.abs(y - other[1]));
    const onLabels = (x, y, label) => placed.reduce((sum, o) => sum + over(x, y, label, [o.x, o.y], o.wide, high), 0);
    const onCells = (x, y, label) => labels.reduce((sum, o) => sum + (o === label ? 0 : over(x, y, label, o.home, 2 * o.reach[0], 2 * o.reach[1])), 0);
    for (const label of labels.slice().sort((a, b) => holdsLabel(b, high) - holdsLabel(a, high))) {
      const [cx, cy] = label.home, dx = label.reach[0] + gap + label.wide / 2, dy = label.reach[1] + gap + high / 2;
      const spots = (holdsLabel(label, high) ? [[cx, cy]] : []).concat([[cx + dx, cy], [cx - dx, cy], [cx, cy - dy], [cx, cy + dy], [cx + dx, cy - dy], [cx - dx, cy - dy], [cx + dx, cy + dy], [cx - dx, cy + dy]])
        .map(([x, y]) => { const at = inside(x, y, label); return { x: at[0], y: at[1], labels: onLabels(at[0], at[1], label), cells: onCells(at[0], at[1], label) }; });
      const best = (holdsLabel(label, high) && spots[0].labels <= 1e-9 ? spots[0] : null) || spots.find((s) => s.labels <= 1e-9 && s.cells <= 1e-9) || spots.find((s) => s.labels <= 1e-9) || spots.reduce((a, b) => (b.labels < a.labels - 1e-9 ? b : a));
      label.x = best.x; label.y = best.y; placed.push(label);
    }
  }
  function drawPicture() {
    const canvas = $('picture'); if (!S.info || !(S.info.hasPicture || S.info.hasMovie) || !canvas.clientWidth) return;
    const { ctx, w, h } = fit(canvas); ctx.fillStyle = '#10161C'; ctx.fillRect(0, 0, w, h);
    const image = Pic.shown || Pic.base; pictureChips(); showView(image, w, h); if (!image) return;
    const at = place(w, h); ctx.save(); ctx.beginPath(); ctx.rect(0, at.roomTop, w, h - at.roomTop); ctx.clip();   // in the room only: the strip of the controls stays clear
    ctx.imageSmoothingEnabled = !pixelsShown(Pic.zoom); ctx.drawImage(image, at.x, at.y, S.info.width * at.z, S.info.height * at.z);
    const selected = S.cell ? S.info.names.indexOf(S.cell) : -1;
    ctx.save(); ctx.translate(at.x, at.y); ctx.scale(at.z, at.z); ctx.lineJoin = 'round';
    paintCells(ctx, PENS.screen, 1 / at.z, selected);
    if (selected >= 0) for (const [colour, width] of [['#000000', 5.5], ['#5BC0FF', 2.6]]) { ctx.strokeStyle = colour; ctx.lineWidth = width * S.outlineWidth / at.z; ctx.stroke(S.paths[selected]); }
    ctx.restore();
    if (selected >= 0) {
      const b = S.bounds[selected], name = shortName(S.cell); ctx.font = '600 12px ' + getComputedStyle(document.body).fontFamily;
      const tw = ctx.measureText(name).width, x = Math.max(4, Math.min(w - tw - 12, at.x + b.x1 * at.z + 6)), y = Math.max(at.roomTop + 16, Math.min(h - 8, at.y + b.y0 * at.z + 12));
      ctx.fillStyle = 'rgba(0,0,0,0.75)'; ctx.beginPath(); ctx.roundRect(x - 4, y - 13, tw + 8, 18, 6); ctx.fill(); ctx.fillStyle = '#fff'; ctx.fillText(name, x, y);
    }
    if (Pic.band) {   // the rectangle being drawn for "Zoom to region"
      const [x0, y0, x1, y1] = Pic.band, x = Math.min(x0, x1) + 0.5, y = Math.min(y0, y1) + 0.5, wide = Math.abs(x1 - x0), high = Math.abs(y1 - y0);
      ctx.fillStyle = 'rgba(91,192,255,0.16)'; ctx.fillRect(x, y, wide, high); ctx.strokeStyle = 'rgba(0,0,0,0.6)'; ctx.lineWidth = 3; ctx.strokeRect(x, y, wide, high);
      ctx.setLineDash([5, 4]); ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 1.2; ctx.strokeRect(x, y, wide, high); ctx.setLineDash([]);
    }
    ctx.restore();
  }
  /** A point of the pointer in the panel's own pixels. */
  const panelPoint = (event) => { const box = $('picture').getBoundingClientRect(); return { x: event.clientX - box.left, y: event.clientY - box.top }; };
  function cellAt(event) {
    const canvas = $('picture'), at = place(canvas.clientWidth, canvas.clientHeight), p = panelPoint(event), x = (p.x - at.x) / at.z, y = (p.y - at.y) / at.z, ctx = canvas.getContext('2d');
    if (p.y < at.roomTop) return -1;   // the strip of the controls shows no cell
    ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0);
    let found = -1; for (let i = 0; i < S.paths.length && found < 0; i++) if (ctx.isPointInPath(S.paths[i], x, y)) found = i;
    ctx.restore(); return found;
  }
  /** Brings a cell to the middle of the view, at the zoom there is. */
  function centreOn(index) { if (index < 0) return; const b = S.bounds[index]; setView(Pic.zoom, (b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2); }
  function zoomToCell() {
    const index = S.cell ? S.info.names.indexOf(S.cell) : -1; if (index < 0) return;
    const canvas = $('picture'), b = S.bounds[index], wanted = Math.min(canvas.clientWidth / (Math.max(4, b.x1 - b.x0) * 6), canvas.clientHeight / (Math.max(4, b.y1 - b.y0) * 6));
    setView(Math.max(1.01, wanted / place(canvas.clientWidth, canvas.clientHeight).fitted), (b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2);
  }

  function selectedCall() { const own = callRow(S.cell); return Object.keys(own).length ? own : (Object.values(calls())[0] || {}); }
  function samplingSupport() {
    const time = S.results.time, intervals = time.slice(1).map((t, i) => t - time[i]).filter((d) => d > 0).sort((a, b) => a - b);
    const typical = intervals.length ? intervals[Math.floor(intervals.length / 2)] : 1;
    const saved = number(selectedCall().Max_gap_s), gap = saved > 0 ? saved : 1.5 * typical;
    const bounds = time.map((t, i) => [i === 0 ? t : t - (time[i] - time[i - 1] <= gap ? (time[i] - time[i - 1]) / 2 : typical / 2), i === time.length - 1 ? t : t + (time[i + 1] - time[i] <= gap ? (time[i + 1] - time[i]) / 2 : typical / 2)]);
    return { typical, gap, bounds };
  }
  function responseSpans() {
    const row = selectedCall(), start = number(row.Start_s), end = number(row.End_s), before = number(row.Baseline_start_s), baselineEnd = number(row.Baseline_end_s);
    return { baseline: finite(before) && finite(baselineEnd) ? [before, baselineEnd] : [start - number(row.Baseline_s), start], response: [start, end], run: [number(row.Run_start_s), number(row.Run_end_s)] };
  }
  /** What set a cell's required rise: its baseline noise ("noise") or the minimum rise ("floor"). The required rise is always
   *  the larger of the noise threshold and the minimum rise, and the call table says which (Threshold_source); a table saved
   *  without that column is read by the engine's own rule. null: the call has no thresholds. */
  function requiredBy(row) {
    if (row.Threshold_source === 'Baseline noise') return 'noise'; if (row.Threshold_source === 'Floor') return 'floor';
    const floor = number(row.Floor_threshold), noise = number(row.Noise_threshold); return finite(floor) && finite(noise) ? (noise > floor ? 'noise' : 'floor') : null;
  }
  /** The levels the selected cell's call was judged against, on the ΔF/F₀ axis, each a line of the one-cell plot. First the
   *  required rise. With "all", also the threshold that did not set it, under its own name (the minimum rise where the noise
   *  set the required rise, the noise threshold where the minimum rise did), and the mean of the baseline. Levels with one
   *  value are one line, named as the first of them. "colour" and "width" name the entries of the style table a level is drawn by. */
  function referenceLevels(all) {
    const row = selectedCall(), baseline = number(row.Baseline_mean); if (!finite(baseline) || baseline <= -1) return [];
    const level = (key, name, colour, rise) => ({ key, name, colour, width: key === 'required' ? 'requiredWidth' : 'levelWidth', value: baseline + (1 + baseline) * rise });
    const levels = [level('required', 'required rise', 'colourLevelRequired', number(row.Threshold_increase))];
    if (all) levels.push(requiredBy(row) === 'noise' ? level('floor', 'minimum rise', 'colourLevelFloor', number(row.Floor_threshold)) : level('noise', 'noise threshold', 'colourLevelNoise', number(row.Noise_threshold)),
      level('baseline', 'baseline mean', 'colourLevelBaseline', 0));
    return levels.filter((one, i) => finite(one.value) && !levels.slice(0, i).some((before) => before.value === one.value));
  }
  // ---------------------------------------------------------------- around the cell: the traces the passing object check used
  /** What Fiji gave for "Around the cell" (api/around), by cell, for the analysis shown: asked for one cell at a time, as it is
   *  selected, and kept for a few cells. generation: the results it belongs to; an answer for earlier results is dropped. */
  const Around = { got: new Map(), asking: new Set(), generation: 0 };
  const AROUND_CONTROL = 'Checks are made for the primary stimulus.', AROUND_KEPT = 12;
  function askAround(name) {
    const r = S.results; if (!r || !name || Around.asking.has(name)) return;
    const had = Around.got.get(name); if (had && !had.failed) return;
    const generation = Around.generation; Around.asking.add(name);
    api('around?cell=' + encodeURIComponent(name) + (r.folder ? '&folder=' + encodeURIComponent(r.folder) : ''))
      .catch((problem) => { const lost = /fetch/i.test(problem.message); return { why: lost ? 'Fiji is not connected.' : problem.message, failed: true, lost }; })
      .then((answer) => {
        if (generation !== Around.generation) return; Around.asking.delete(name);
        Around.got.delete(name); Around.got.set(name, answer); while (Around.got.size > AROUND_KEPT) Around.got.delete(Around.got.keys().next().value);
        if (name === S.cell) { showAroundButton(); Plot.dirty = true; drawPlot(); }
      });
  }
  /** Once Fiji answers a poll again: the selected cell's answer that did not come because Fiji was not connected is asked again. */
  function askAroundAgain() {
    const got = S.cell && Around.got.get(S.cell); if (!got || !got.lost || Around.asking.has(S.cell)) return;
    Around.got.delete(S.cell); askAround(S.cell); showAroundButton(); Plot.dirty = true; drawPlot();
  }
  /** "Around the cell" for the selected cell on the window shown: { data } to draw, { why } when there is nothing to show, and
   *  {} while Fiji's answer is on its way. */
  function aroundNow() {
    if (!S.results) return { why: '' }; if (!S.cell) return { why: 'No cell selected.' }; if (!primaryShown()) return { why: AROUND_CONTROL };
    const got = Around.got.get(S.cell); if (!got) return {};
    return got.why || !Array.isArray(got.time) ? { why: got.why || '' } : { data: got };
  }
  /** The colour of each part of a flag on "Around the cell", the same in the Fiji window: a style entry. A colour met a second
   *  time for one cell, or a part without one, takes the next of colourAgain. */
  const PART_COLOUR = { before: 'colourBefore', rim: 'colourBefore', after: 'colourAfter', ring: 'colourAfter', core: 'colourCore' };
  /** A zone line's colour in a look. */
  const zoneInk = (look, line) => line.colour === 'colourAgain' ? look.colourAgain[Math.min(line.again, look.colourAgain.length - 1)] : look[line.colour];
  /** The lines of "Around the cell" from one cell's data: its time; the cell's rise; the zones the flag rests on, each by its
   *  place, in the colour of its part (PART_COLOUR); every other usable near-ring sector; and the flagged frames (band, frames
   *  from 1) as a span of time, from half the median frame step before the first to half of it after the last, as in the Fiji
   *  window. Any mix of parts is drawn. */
  function aroundParts(a) {
    const zones = new Map((a.zones || []).map((z) => [z.zone, z])), time = a.time;
    const trace = (zone) => { const z = zones.get(zone); return z && z.usable !== false && Array.isArray(z.trace) ? z.trace : null; };
    const placeOf = (zone, given) => given || (zones.get(zone) || {}).place || 'zone ' + zone;
    const used = new Set(); let again = 0;
    const colourOf = (part) => { const own = PART_COLOUR[part]; if (own && !used.has(own)) { used.add(own); return { colour: own }; } return { colour: 'colourAgain', again: again++ }; };
    const lines = (a.lines || []).filter((line) => line && trace(line.zone)).map((line) => Object.assign({ zone: line.zone, part: line.part, frame: line.frame == null ? null : line.frame, place: placeOf(line.zone, line.place), values: trace(line.zone) }, colourOf(line.part)));
    const others = (a.others || []).filter((zone) => trace(zone)).map((zone) => ({ zone, place: placeOf(zone), values: trace(zone) }));
    const frames = Array.isArray(a.band) && a.band.length === 2 && Number.isInteger(a.band[0]) && Number.isInteger(a.band[1]) && a.band[0] >= 1 && a.band[1] >= a.band[0] && a.band[1] <= time.length ? a.band : null;
    // the median frame step, as the Fiji window takes it: the middle one of the steps in order, the upper middle of an even count
    const steps = time.slice(1).map((t, i) => t - time[i]).filter((d) => d > 0).sort((x, y) => x - y), half = (steps.length ? steps[steps.length >> 1] : 0) / 2;
    const band = frames ? [time[frames[0] - 1] - half, time[frames[1] - 1] + half] : null;
    return { time, cell: a.rise, lines, others, frames, band, rule: a.rule || '', text: a.text || '' };
  }
  /** The "Around the cell" button: off, with why in its tooltip, when there is nothing to show; off without a tooltip while an
   *  answer is on its way. */
  function showAroundButton() {
    const button = document.querySelector('[data-plot="around"]'), now = aroundNow();
    if (!now.data && now.why === undefined) { button.disabled = true; button.title = ''; return; }
    button.disabled = !now.data; button.title = now.data ? button.dataset.title : now.why || button.dataset.title;
  }
  /** Whether "Around the cell" is on screen with something to draw. */
  const aroundDrawn = () => S.plotMode === 'around' && !!aroundNow().data;
  /** The plot's Y axis title, as the painters write it when none is typed: the signal, or the gray levels of "Around the cell". */
  const plotYTitle = () => S.plotMode === 'around' ? 'ΔF (gray levels)' : signalName();
  /** What the Y axis limits typed under the plot belong to: the signal on screen, or "Around the cell", which is in gray levels. */
  const yKey = () => S.plotMode === 'around' ? 'around' : S.signal;
  /** The kind of Y axis a plot has, for what the figure editor keeps (Y_TYPED): the signal's, or the gray levels of "Around the cell". */
  const yKind = (mode) => mode === 'around' ? 'around' : 'signal';
  const Y_TYPED = ['yTitle', 'yMin', 'yMax', 'yStep'];
  /** The style entries of "Around the cell" alone: the record beside a figure of another plot leaves them out. */
  const AROUND_STYLE = ['zoneWidth', 'colourBefore', 'colourAfter', 'colourCore', 'colourAgain', 'otherWidth', 'colourOtherSectors', 'flagged', 'colourFlagged'];
  // ---------------------------------------------------------------- style: how the plot and the heatmap look
  const PT = 96 / 72;   // pixels in a point
  /** The fonts a figure can be written in, in the order the editor offers them: [key, name, family]. The one list for the editor
   *  and for the painters. The page's own font, last, has no family here: it is whatever this computer shows the page in. */
  const FONTS = [['arial', 'Arial', 'Arial, Helvetica, sans-serif'], ['helvetica', 'Helvetica', 'Helvetica, Arial, sans-serif'], ['times', 'Times New Roman', '"Times New Roman", Times, serif'], ['page', 'As on the page', null]];
  /** The style table: the one place that says how the plot and the heatmap look before anybody edits them. A saved figure
   *  and the page are drawn from these entries by the same painters; PAGE_STYLE lists the few the page sets differently.
   *  The figure editor shows some entries as boxes (FIGURE_FIELDS); any other entry can be offered by adding a row there.
   *  Units as in a graphing program: sizes in inches; text, line widths, dashes and tick lengths in points (STYLE_POINTS);
   *  opacity from 0 to 1. null: the figure's own text (for the title: none), or the automatic value. A saved figure is
   *  written in a font with a name (Arial), not in whatever font this computer shows the page in. */
  const STYLE_SHARED = { title: null, xTitle: null, yTitle: null, font: 'arial', text: 9, titleText: 10, ink: '#000000', background: 'white' };
  const STYLE = {
    plot: Object.assign({}, STYLE_SHARED, { width: 7, height: 4,
      // axes: left and bottom only, ticks pointing out, no grid, a dotted line at zero
      xMin: null, xMax: null, xStep: null, yMin: null, yMax: null, yStep: null, frame: 'open', grid: false, axisWidth: 0.75, tickLength: 3, colourGrid: '#dee4ea', colourZero: '#9aa5b1', zeroDash: [0.75, 2.25],
      // a group's mean in the group's colour on its SEM band; a white casing keeps it clear of the single cells under it
      meanWidth: 1.8, colourMean: '#0072b2', colourResponded: CALL.responded.toLowerCase(), colourNot: CALL.not.toLowerCase(), sem: 'band', semOpacity: 0.25, casing: 0.7, colourCasing: '#ffffff',
      // the selected cell
      cellWidth: 1.3, colourCell: '#1f2933',
      // single cells behind the mean: thin, neutral and faint (opacity null: by their number, traceFade); responders after them, a little wider and stronger
      traceWidth: 0.6, colourTrace: CALL.other.toLowerCase(), traceOpacity: null, responderWidth: 0.9, responderOpacity: null,
      // marks: the legend; the stimulus line, with its label above the plot; the baseline band; the response window, a thin bar above the plot
      legend: 'top', stimulus: true, stimulusText: 'stimulus', stimulusWidth: 0.75, stimulusDash: [3, 2.25], colourStimulus: '#b45f06',
      baseline: true, colourBaseline: '#edf0f4', window: true, windowBar: 2.25, colourWindow: '#e9a84c',
      // the selected cell's call. Its required rise, a dashed line (levels); the "Thresholds" box of the page adds the other
      // threshold and the baseline mean, thinner (referenceLevels)
      levels: true, levelWidth: 0.75, requiredWidth: 1.35, levelDash: [3.75, 3], colourLevelBaseline: '#7c8794', colourLevelFloor: '#6a3d9a', colourLevelNoise: '#a46ea0', colourLevelRequired: CALL.responded.toLowerCase(),
      // and the run that counted, a band: never narrower than runWidth, with a tick at the top of the plot where it starts, so
      // that a run of two or three frames can be seen
      run: true, colourRun: '#d5eedd', runWidth: 2.25, runTick: 4.5,
      // "Around the cell", in the Fiji window's colours: the zones a flag rests on by their part (PART_COLOUR), before and rim,
      // after and ring, the core, and a colour met a second time in the next of colourAgain; none is a colour of a call. The
      // other near-ring sectors thin and gray, and the flagged frames as a band
      zoneWidth: 1.3, colourBefore: '#0072b2', colourAfter: '#aa3377', colourCore: '#009e73', colourAgain: ['#e69f00', '#56b4e9', '#cc79a7', '#882255'], otherWidth: 0.6, colourOtherSectors: '#aeb8c2',
      flagged: true, colourFlagged: '#e2d6f3' }),
    heatmap: Object.assign({}, STYLE_SHARED, { width: 7, height: 5, barTitle: null,
      // time: the whole recording, or the limits typed in the editor. The colour range: automatic (heatmapRange), or the two numbers typed there
      xMin: null, xMax: null, xStep: null, barMin: null, barMax: null,
      // the colour bar: narrow, at the right, as tall as the picture; a row without a trace, and a pause in the recording
      bar: true, barWidth: 8, colourMissing: '#e5e5e5',
      // beside the rows: a strip with each row's call and a mark for a cell below zero, each a solid bar over its run of rows and
      // never shorter than stripLeast, so a group of one cell can be found; between the groups of rows, a white line
      strip: true, stripWidth: 3, flagWidth: 1.5, stripLeast: 3, colourOther: CALL.other.toLowerCase(), colourFlag: '#1f2933', divider: 1.5, colourDivider: '#ffffff',
      // above the picture: the stimulus label, the baseline and the response window as thin bars. Over it: the stimulus line, dark on a
      // palette with white in it and light on the others
      marks: true, colourBaseline: '#b8c2cc', colourOnLight: 'rgba(31,41,51,0.55)', colourOnDark: 'rgba(255,255,255,0.75)',
      // the cell selected on the page
      colourSelected: '#1f2933' })
  };
  /** What the heatmap draws as the plot does (its axis, the marks of time, the colours of the calls) is the plot's entry, written once. */
  for (const key of ['axisWidth', 'tickLength', 'stimulusText', 'stimulusWidth', 'stimulusDash', 'colourStimulus', 'colourWindow', 'windowBar', 'colourResponded', 'colourNot']) STYLE.heatmap[key] = STYLE.plot[key];
  /** The entries that are lengths in points; a painter gets them in pixels. */
  const STYLE_POINTS = ['text', 'titleText', 'axisWidth', 'tickLength', 'zeroDash', 'meanWidth', 'casing', 'cellWidth', 'traceWidth', 'responderWidth', 'stimulusWidth', 'stimulusDash', 'windowBar', 'levelWidth', 'requiredWidth', 'levelDash',
    'runWidth', 'runTick', 'zoneWidth', 'otherWidth', 'barWidth', 'stripWidth', 'flagWidth', 'stripLeast', 'divider'];
  /** What the page sets differently from a saved figure, the plot and the heatmap alike: its ink and its own font. "page" also
   *  tells a painter to keep what it drew for the pointer (Plot.geom and the rest) and to mark the selected cell. */
  const PAGE_STYLE = { plot: { page: true, ink: '#3e4c59', font: 'page' } }; PAGE_STYLE.heatmap = PAGE_STYLE.plot;
  /** A style as a painter reads it, a "look": the same entries with points turned into pixels, the font as a family and the
   *  background as a fill (null: none); px and titlePx are the two text sizes. One function for the page and for a saved
   *  figure, for the plot and the heatmap. */
  function lookOf(style) {
    const look = Object.assign({}, style);
    for (const key of STYLE_POINTS) if (look[key] != null) look[key] = Array.isArray(look[key]) ? look[key].map((v) => v * PT) : look[key] * PT;
    return Object.assign(look, { px: look.text, titlePx: look.titleText, family: (FONTS.find((font) => font[0] === style.font) || [])[2] || null, background: style.background === 'transparent' ? null : '#fff' });
  }
  const PAGE_LOOK = lookOf(Object.assign({}, STYLE.plot, PAGE_STYLE.plot)), PAGE_HEAT = lookOf(Object.assign({}, STYLE.heatmap, PAGE_STYLE.heatmap));
  /** How strongly single cells are drawn when the style leaves it automatic: fainter the more there are, so a hundred cells
   *  stay a light cloud under the mean and a handful stay readable. n is every single cell drawn, both kinds; a responder gets
   *  half as much again, so it shows in the cloud without outshining the mean. */
  const traceFade = (n) => Math.max(0.06, Math.min(0.8, 1.5 / Math.sqrt(Math.max(1, n))));
  const responderFade = (n) => Math.min(0.9, 1.5 * traceFade(n));
  // ---------------------------------------------------------------- the traces
  /** ΔF/F₀ of the raw values (no background taken off), worked out here for display: each cell's raw trace against its own mean
   *  over the F₀ frames of the analysis. The response calls do not use it. */
  function rawDff() {
    const r = S.results; if (!r || !r.raw || !Object.keys(r.raw).length) return null;
    if (!r.rawdff) {
      const frames = r.baselineFrames || [1, 1], out = {};
      for (const name of Object.keys(r.raw)) {
        const y = r.raw[name]; let sum = 0, k = 0; for (let i = frames[0] - 1; i < frames[1] && i < y.length; i++) if (finite(y[i])) { sum += y[i]; k++; }
        const f0 = k ? sum / k : NaN; out[name] = y.map((v) => finite(v) && f0 > 0 ? (v - f0) / f0 : NaN);
      }
      r.rawdff = out;
    }
    return r.rawdff;
  }
  const hasTable = (key) => { if (key === 'rawdff') rawDff(); return !!(S.results && S.results[key] && Object.keys(S.results[key]).length); };
  /** The traces the plot draws: ΔF/F₀ as the analysis saved it, or the measurement's gray values (background corrected, or raw). */
  const shownTraces = () => S.signal !== 'dff' && hasTable(S.signal) ? S.results[S.signal] : S.results.traces;
  const isDff = () => !(S.signal !== 'dff' && hasTable(S.signal));
  const signalName = () => isDff() ? 'ΔF/F₀' : S.signal === 'raw' ? 'Raw gray value (not corrected)' : S.signal === 'rawdff' ? 'ΔF/F₀ of raw values (no background removed)' : 'Gray levels above background';
  /** Mean and SEM of some cells' saved traces, frame by frame (a display summary; the tables are not changed). traces: cell → values (when not given, the signal on screen). */
  function groupMean(names, traces) {
    const r = S.results, T = traces || shownTraces(), n = r.time.length, mean = new Array(n), sem = new Array(n), used = names.filter((name) => T[name]);
    for (let i = 0; i < n; i++) {
      let sum = 0, squares = 0, k = 0;
      for (const name of used) { const v = T[name][i]; if (finite(v)) { sum += v; squares += v * v; k++; } }
      mean[i] = k ? sum / k : NaN; sem[i] = k > 1 ? Math.sqrt(Math.max(0, (squares - sum * sum / k) / (k - 1)) / k) : NaN;
    }
    return { mean, sem, n: used.length };
  }
  /** The cells the "Leave out … below zero" switch leaves out of the headline, and so of the plot: those the tally lists, while
   *  the switch is on and the primary stimulus is selected (the switch belongs to it). */
  function leftOut() {
    const r = S.results, control = ((r.windows || []).find((w) => w.key === S.window) || {}).role === 'Control';
    return new Set(S.withoutBelowZero && r.tally && !control ? r.tally.belowZero || [] : []);
  }
  /** The cells a plot uses. The one place that decides it, for the single traces, the means and the automatic range alike:
   *    shown       the cells of a group that have a trace of the signal on screen, without the cells left out by the
   *                "Leave out … below zero" switch. These are drawn and averaged.
   *    steady      those of them an automatic range is taken from: neither excluded nor below zero. A cell that is off for
   *                the whole recording is still drawn but cannot set a limit. All the shown cells when none is steady.
   *  group: an index of GROUPS (when not given, the group under Show); traces: cell → values (when not given, the signal on screen). */
  function plotCells(group, traces) {
    const T = traces || shownTraces(), which = group == null ? S.group : group, out = leftOut();
    const shown = S.info.names.filter((name) => T[name] && !out.has(name) && inGroup(callOf(name), which, name));
    const steady = shown.filter((name) => callOf(name) !== 'Excluded' && !belowZero(name));
    return { shown, steady: steady.length ? steady : shown };
  }
  /** The mean lines of the plot, one per group: its name, the mean and SEM frame by frame, the number of cells, and the entry
   *  of the style table it is coloured by. "Average" and "One cell + average" show the cells under Show; "Responders vs
   *  non-responders" shows every cell by its call. All cells in ΔF/F₀ with nobody left out is the mean the analysis saved;
   *  any other mean is averaged here from the saved traces. */
  function plotSeries() {
    const r = S.results; if (!r) return [];
    if (S.plotMode === 'calls') {
      const cells = plotCells(0).shown;
      return [['Responded', 'responded', 'colourResponded'], ['Did not respond', 'did not respond', 'colourNot']]
        .map(([status, name, colour]) => Object.assign({ name, short: name, colour }, groupMean(cells.filter((n) => callOf(n) === status)))).filter((s) => s.n > 0);
    }
    const cells = plotCells().shown, name = GROUPS[S.group].toLowerCase(), saved = S.group === 0 && isDff() && cells.length === Object.keys(r.traces).length;
    return [Object.assign({ name, short: 'mean of ' + name, colour: 'colourMean' }, saved ? { mean: r.mean, sem: r.sem, n: cells.length } : groupMean(cells))].filter((s) => s.n > 0);
  }
  /** Whether a mean gets its SEM band: when the look draws bands and there is more than one cell (one cell has no SEM). */
  const bandDrawn = (series, look) => look.sem !== 'none' && series.n > 1;
  /** A name with the number of cells it stands for, as every key writes it: "excluded (n = 7)". */
  const counted = (text, n) => text + ' (n = ' + n + ')';
  /** What a mean line is called, in one place, for the legend and for the file kept beside a figure: "all cells ± SEM (n = 113)";
   *  without "± SEM" where no band is drawn, as for one cell: "responded (n = 1)". */
  const meanLabel = (series, band) => counted(series.name + (band ? ' ± SEM' : ''), series.n);
  const shade = (hex, alpha) => { const c = rgb(hex); return 'rgba(' + c[0] + ',' + c[1] + ',' + c[2] + ',' + alpha + ')'; };
  /** Whether "Leave out … below zero" takes cells out of the means the plot shows: a cell left out that would be in them. */
  const meansLeaveOut = () => { const T = shownTraces(), group = S.plotMode === 'calls' ? 0 : S.group; return Array.from(leftOut()).some((name) => T[name] && inGroup(callOf(name), group, name)); };
  /** What a saved plot is called, in its file name; a plot whose means leave out the cells below zero says so ("_without_below_zero"),
   *  so that a file named mean_trace always holds the mean of every cell of its group, as the release saved it. */
  const plotKind = () => S.plotMode === 'around' ? 'around_the_cell' : (S.plotMode === 'cell' ? 'traces' : S.plotMode === 'calls' ? 'responders_vs_nonresponders' : 'mean_trace' + (S.group ? '_' + GROUPS[S.group].toLowerCase().replace(/[^a-z0-9]+/g, '_') : '')) + (isDff() ? '' : S.signal === 'raw' ? '_raw_gray_values' : S.signal === 'rawdff' ? '_dFF0_of_raw' : '_gray_levels')
    + (meansLeaveOut() ? '_without_below_zero' : '');
  /** The numbers behind the plot as a table: time, each mean line with its SEM and number of cells (named for the cells below zero
   *  when they are left out of it: mean_of_all_cells_without_below_zero), and the selected cell. */
  function plotNumbers() {
    const r = S.results, series = plotSeries(), own = S.plotMode === 'cell' && S.cell && shownTraces()[S.cell] ? shownTraces()[S.cell] : null;
    const clean = (text) => String(text).replace(/[^A-Za-z0-9]+/g, '_').replace(/^_+|_+$/g, ''), without = meansLeaveOut() ? '_without_below_zero' : '';
    // "Around the cell": time, the cell's rise, each zone of its flag and each other ring sector, in gray levels, under their places
    const a = S.plotMode === 'around' ? aroundNow().data : null;
    if (a) { const p = aroundParts(a), columns = [p.cell].concat(p.lines.map((line) => line.values), p.others.map((other) => other.values)), value = (v) => finite(v) ? v : '';
      return [['Time_s', S.cell].concat(p.lines.map((line) => clean(line.place)), p.others.map((other) => clean(other.place))).join(',')].concat(p.time.map((t, i) => [t].concat(columns.map((c) => value(c && c[i]))).join(','))).join('\n') + '\n'; }
    const head = ['Time_s']; for (const s of series) head.push(clean(s.short) + without, clean(s.short) + without + '_SEM', clean(s.short) + without + '_n_cells'); if (own) head.push(S.cell);
    const rows = [head.join(',')];
    for (let i = 0; i < r.time.length; i++) {
      const row = [r.time[i]]; for (const s of series) row.push(finite(s.mean[i]) ? s.mean[i] : '', finite(s.sem[i]) ? s.sem[i] : '', s.n);
      if (own) row.push(finite(own[i]) ? own[i] : ''); rows.push(row.join(','));
    }
    return rows.join('\n') + '\n';
  }
  const realSpan = (span) => !!span && finite(span[0]) && finite(span[1]);
  /** An axis whose ends can be typed: its automatic ends with each typed one in its place. null when that leaves no range (a
   *  minimum at or above the maximum). The one rule for the painters, which then keep the automatic ends, and for the boxes
   *  the limits are typed into, which say so (figureRefusal; the page's two Y boxes). */
  const typedRange = (auto, min, max) => { const lo = finite(min) ? min : auto[0], hi = finite(max) ? max : auto[1]; return hi > lo ? [lo, hi] : null; };
  /** The recording's first and last frame time. */
  const wholeTime = () => { const time = S.results.time; return [time[0], time[time.length - 1]]; };
  /** The time the plot shows. auto: the whole recording, so the axis ends where the data end, or the page's zoom. x0, x1: that,
   *  or the limits typed in the figure editor. full0 and full1 are the recording's first and last frame. */
  function plotTimes(look) {
    const [full0, full1] = wholeTime(), auto = [Plot.x0 === null ? full0 : Plot.x0, Plot.x1 === null ? full1 : Plot.x1], [x0, x1] = typedRange(auto, look.xMin, look.xMax) || auto;
    return { x0, x1, full0, full1, auto };
  }
  /** Whether the plot is asked for a cell it cannot draw: "One cell + average" with a selected cell that has no trace of the signal on screen. */
  const traceMissing = () => S.plotMode === 'cell' && !!S.cell && !shownTraces()[S.cell];
  /** What the plot draws in the time it shows (x0 to x1), gathered once for the painter and its legend:
   *    series    the mean lines            traces  the single cells drawn behind them      own  the selected cell's trace
   *    spans     the baseline, the response window and the run that counted: each only when it is drawn, null otherwise
   *    stimulus  when it was added; null when it is not marked
   *    levels    the selected cell's threshold lines (referenceLevels): the required rise unless the look leaves it out, and the
   *              others when "Thresholds" is ticked. They and the run belong to one cell's call. */
  function plotParts(look, x0, x1) {
    const T = shownTraces(), oneCell = S.plotMode === 'cell', all = responseSpans(), onset = number(selectedCall().Stimulus_onset_s), shown = (span) => realSpan(span) && span[1] > x0 && span[0] < x1;
    const stimulus = look.stimulus && finite(onset) && onset >= x0 && onset <= x1 ? onset : null;
    // "Around the cell": the cell, the zones of its flag, the other ring sectors and the flagged frames, and the stimulus; nothing else
    if (S.plotMode === 'around') { const a = aroundNow().data, around = a ? aroundParts(a) : null; if (around && !(look.flagged && around.band && around.band[1] > x0 && around.band[0] < x1)) around.band = null;
      return { series: [], traces: [], own: null, spans: { baseline: null, response: null, run: null }, stimulus, levels: [], around }; }
    const called = oneCell && isDff() && !!S.cell && Object.keys(callRow(S.cell)).length > 0;
    return { series: plotSeries(), traces: S.every ? plotCells().shown : [], own: oneCell && S.cell && T[S.cell] ? T[S.cell] : null,
      spans: { baseline: look.baseline && shown(all.baseline) ? all.baseline : null, response: look.window && shown(all.response) ? all.response : null,
        run: called && look.run && realSpan(all.run) && all.run[1] >= x0 && all.run[0] <= x1 ? all.run : null },   // (a run of one sample has no length and is drawn all the same)
      stimulus, levels: called ? referenceLevels(S.thresholds).filter((level) => look.levels || level.key !== 'required') : [], around: null };
  }
  /** The legend, in one place: one entry for each thing the plot draws and none for anything else, each with what its swatch
   *  needs to show the thing as drawn. In order: the selected cell; each group once (its mean line, on its band when a band
   *  is drawn); "single cells", which are drawn in two ways and so have a swatch of two lines; the selected cell's threshold
   *  lines; then the bands, "baseline", "response window" and the "response" of a cell that responded. The stimulus is named
   *  above its own line and has no entry. */
  function plotLegend(parts, look) {
    const items = [], around = parts.around;
    if (around) {   // the cell; each zone of its flag by its place; the other ring sectors once; the flagged frames
      items.push({ kind: 'line', colour: look.colourCell, width: look.cellWidth, text: shortName(S.cell) });
      for (const line of around.lines) items.push({ kind: 'line', colour: zoneInk(look, line), width: look.zoneWidth, text: line.place });
      if (around.others.length) items.push({ kind: 'line', colour: look.colourOtherSectors, width: look.otherWidth, text: 'other ring sectors' });
      if (around.band) items.push({ kind: 'box', colour: look.colourFlagged, outline: true, text: 'flagged frames' });
      return items;
    }
    if (parts.own) items.push({ kind: 'line', colour: look.colourCell, width: look.cellWidth, text: shortName(S.cell) });
    for (const s of parts.series) { const band = bandDrawn(s, look); items.push({ kind: 'line', colour: look[s.colour], width: look.meanWidth, band: band ? shade(look[s.colour], look.semOpacity) : null, text: meanLabel(s, band) }); }
    if (parts.traces.length) {   // one entry; its swatch is a line of each kind that is drawn, the responders' under the others'
      const yes = parts.traces.filter((name) => callOf(name) === 'Responded').length, lines = [];
      if (yes < parts.traces.length) lines.push({ colour: look.colourTrace, width: look.traceWidth }); if (yes) lines.push({ colour: look.colourResponded, width: look.responderWidth });
      items.push({ kind: 'line', colour: lines[0].colour, width: lines[0].width, second: lines[1] || null, text: 'single cells' });
    }
    for (const level of parts.levels) items.push({ kind: 'line', colour: look[level.colour], width: look[level.width], dash: look.levelDash, text: level.name });
    if (parts.spans.baseline) items.push({ kind: 'box', colour: look.colourBaseline, outline: true, text: 'baseline' });
    if (parts.spans.response) items.push({ kind: 'box', colour: look.colourWindow, height: look.windowBar, text: 'response window' });
    if (parts.spans.run) items.push({ kind: 'box', colour: look.colourRun, outline: true, text: 'response' });
    return items;
  }
  /** A trace's own extent: its 1st and 99th percentile (what an automatic range is taken from, so that a spike does not count)
   *  and its smallest and largest value. Kept per trace, since the traces of a result do not change. null: no number in it. */
  const Extents = new WeakMap();
  function extent(trace) {
    if (!Extents.has(trace)) { const v = trace.filter(finite).sort((a, b) => a - b), at = (part) => v[Math.round(part * (v.length - 1))];
      Extents.set(trace, v.length ? { low: at(0.01), high: at(0.99), min: v[0], max: v[v.length - 1] } : null); }
    return Extents.get(trace);
  }
  /** The automatic Y range, from what the plot draws:
   *    the means    everything drawn as a mean, with its SEM where a band is drawn, and zero where the axis has a zero;
   *    the bulk     with single cells drawn, from the 2nd percentile of the cells' lows to the 95th percentile of their highs.
   *                 Each cell counts by its own extent, over the cells that are neither excluded nor below zero (plotCells),
   *                 never all samples pooled: one cell that is off for the whole recording holds more samples than any
   *                 pooled trim removes.
   *  The range holds both, but reaches no further beyond the means than one and a half times their own height on either side:
   *  the means then fill at least a quarter of the axis, whatever a few cells do, and a trace that leaves the range ends at the
   *  plot's edge. Always inside the range: the selected cell's whole trace and its thresholds.
   *  Then 5 % of room, and outward to the next half step of the ticks: the numbers stay round, and no more than half a step
   *  lies empty beyond the data. */
  function plotRange(parts, look) {
    const T = shownTraces(), least = (values) => values.reduce((a, b) => Math.min(a, b), Infinity), most = (values) => values.reduce((a, b) => Math.max(a, b), -Infinity), drawn = [];
    for (const s of parts.series) { const band = bandDrawn(s, look); s.mean.forEach((v, i) => { if (finite(v)) { const e = band && finite(s.sem[i]) ? s.sem[i] : 0; drawn.push(v - e, v + e); } }); }
    if (drawn.length && (isDff() || S.signal === 'rawdff')) drawn.push(0);
    const means = drawn.length ? [least(drawn), most(drawn)] : null; let range = means;
    if (parts.traces.length) {
      const ends = (end) => plotCells().steady.map((name) => extent(T[name])).filter(Boolean).map((e) => e[end]), lows = ends('low'), highs = ends('high');
      if (lows.length) { const bulk = [share([lows], 0.02), share([highs], 0.95)], reach = means && means[1] > means[0] ? 1.5 * (means[1] - means[0]) : Infinity;
        range = means ? [Math.max(Math.min(means[0], bulk[0]), means[0] - reach), Math.min(Math.max(means[1], bulk[1]), means[1] + reach)] : bulk; }
    }
    const values = range ? range.slice() : [], own = parts.own && extent(parts.own); if (own) values.push(own.min, own.max);
    for (const level of parts.levels) values.push(level.value);
    // "Around the cell": every line it draws, whole, and zero, the cell's level over its F₀ frames
    if (parts.around) for (const trace of [parts.around.cell].concat(parts.around.lines.map((line) => line.values), parts.around.others.map((other) => other.values))) { const e = trace && extent(trace); if (e) values.push(e.min, e.max, 0); }
    let lo = least(values), hi = most(values);
    if (!values.length) { lo = 0; hi = 1; } else if (hi - lo < 0.001) hi = lo + 1;
    const room = (hi - lo) * 0.05; lo -= room; hi += room;
    for (let step = 0, i = 0; i < 4; i++) {   // rounding widens the range, which can call for a wider tick step: until the step stays
      const next = ticks(lo, hi, 5).step; if (next === step) break; step = next;
      lo = +(Math.floor(2 * lo / step + 1e-9) * step / 2).toPrecision(12); hi = +(Math.ceil(2 * hi / step - 1e-9) * step / 2).toPrecision(12);
    }
    return [lo, hi];
  }
  /** The plot's Y range. auto: the automatic range (plotRange), with the limits typed under "Y axis" on the page in their place
   *  (either end alone, too; they belong to the signal on screen). lo, hi: that, or the limits typed in the figure editor. */
  function plotLimits(parts, look) {
    const range = plotRange(parts, look), typed = Plot.y[yKey()] || {}, auto = typedRange(range, typed.lo, typed.hi) || range, [lo, hi] = typedRange(auto, look.yMin, look.yMax) || auto;
    return { lo, hi, auto };
  }
  /** Places legend entries: in rows from the left, wrapped by whole entries, with one row height and even gaps; or stacked in
   *  one column. As few rows as the width allows, and those rows as even as they can be: two entries and two, or three and
   *  two, never three with the fourth left alone. No entry is wider than the width: the words of one that would be go on two
   *  lines or more (legendLines), and its row is that much taller. Returns each entry's place and lines and the size of the
   *  block; drawLegend draws it. ctx carries the font. Lengths are written for 12-pixel text, as in the painters. */
  function legendLayout(ctx, items, width, look, stacked) {
    const u = (n) => n * look.px / 12, rowHeight = u(18), lineHeight = u(15), gap = u(14), lines = items.map((item) => legendLines(ctx, item.text, width - u(21)));
    const wide = lines.map((words) => u(21) + words.reduce((most, line) => Math.max(most, ctx.measureText(line).width), 0)), tall = lines.map((words) => rowHeight + (words.length - 1) * lineHeight);
    /** The rows, filled from the left: an entry starts a new row when it would pass the width, or when its row holds "most" entries already. */
    const fill = (most, sizes) => { const rows = []; for (const w of sizes || wide) { const row = rows[rows.length - 1]; if (row && !stacked && row.n < most && row.width + gap + w <= width) { row.n++; row.width += gap + w; } else rows.push({ n: 1, width: w }); } return rows; };
    let rows = fill(Infinity);
    // the same number of rows with the entries shared out evenly, when they fit that way; else with one more to a row, and so on
    if (!stacked) for (let most = Math.ceil(wide.length / rows.length); most < wide.length; most++) { const even = fill(most); if (even.length === rows.length) { rows = even; break; } }
    // entries of different widths can still leave the last one alone ("response window" under two long entries): filled from
    // the end instead, the same number of rows may end on two
    if (!stacked && rows.length > 1 && rows[rows.length - 1].n === 1) { const back = fill(Infinity, wide.slice().reverse()).reverse(); if (back.length === rows.length && back[back.length - 1].n > 1) rows = back; }
    const places = []; let i = 0, y = 0;
    for (const row of rows) { let x = 0, high = 0; for (let k = 0; k < row.n; k++, i++) { places.push({ item: items[i], lines: lines[i], text: lines[i].join(' '), x, y }); x += wide[i] + gap; high = Math.max(high, tall[i]); } y += high; }
    return { places, rowHeight, lineHeight, width: rows.reduce((most, row) => Math.max(most, row.width), 0), height: y };
  }
  /** The words of a legend entry as lines no wider than "room": one line when they fit; else broken between words, never inside
   *  "± SEM" or the number of cells, "(n = 116)", so that nothing of the entry is lost. Only a word wider than the room on its
   *  own is shortened (fitted). */
  function legendLines(ctx, text, room) {
    if (ctx.measureText(text).width <= room) return [text];
    const lines = [];
    for (const word of text.match(/± \S+|\(n = [^)]*\)|\S+/g) || [text]) {
      const longer = lines.length ? lines[lines.length - 1] + ' ' + word : null;
      if (longer !== null && ctx.measureText(longer).width <= room) lines[lines.length - 1] = longer; else lines.push(fitted(ctx, word, room));
    }
    return lines;
  }
  /** Draws a legend placed by legendLayout with its left top corner at (left, top). A swatch shows the thing as it is drawn:
   *  a line at its own width and dashes, on its band when it has one (two lines, one above the other, for an entry that
   *  stands for lines of two kinds); a box in a band's colour, as thin as a bar when it stands for one. The swatch stands
   *  beside the first line of an entry's words. */
  function drawLegend(ctx, layout, left, top, look) {
    const u = (n) => n * look.px / 12; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.lineCap = 'butt';
    for (const { item, lines, x, y } of layout.places) {
      const sx = left + x, sy = top + y + layout.rowHeight / 2;
      if (item.band) { ctx.fillStyle = item.band; ctx.fillRect(sx, sy - u(5), u(16), u(10)); }
      if (item.kind === 'line') { const rule = (line, y) => { ctx.strokeStyle = line.colour; ctx.lineWidth = line.width; ctx.setLineDash(line.dash || []); ctx.beginPath(); ctx.moveTo(sx, y); ctx.lineTo(sx + u(16), y); ctx.stroke(); ctx.setLineDash([]); };
        if (item.second) { rule(item, sy - u(2.5)); rule(item.second, sy + u(2.5)); } else rule(item, sy); }
      else { const high = item.height || u(10); ctx.fillStyle = item.colour; ctx.fillRect(sx + u(2), sy - high / 2, u(14), high);
        if (item.outline) { ctx.strokeStyle = shade(look.ink, 0.3); ctx.lineWidth = 1; ctx.strokeRect(sx + u(2) + 0.5, sy - high / 2 + 0.5, u(14) - 1, high - 1); } }
      ctx.fillStyle = look.ink; lines.forEach((line, k) => writeText(ctx, line, sx + u(21), sy + k * layout.lineHeight, look));
    }
  }
  /** A text as it fits into a width: the text itself, or its beginning with an ellipsis in place of the rest. The one rule for
   *  text on one line that does not fit, wherever a painter writes words that somebody typed: titles, the stimulus label, a
   *  word of a legend entry wider than the legend. ctx carries the font. */
  function fitted(ctx, text, room) {
    if (ctx.measureText(text).width <= room) return text;
    let n = text.length; while (n > 0 && ctx.measureText(text.slice(0, n).trimEnd() + '…').width > room) n--;
    return text.slice(0, n).trimEnd() + '…';
  }
  /** A title centred at "at" along a side of the drawing, as it fits there without running past either end of the side. */
  const centred = (ctx, text, at, side, look) => fitted(ctx, text, 2 * Math.min(at, side - at) - look.px / 2);
  /** Writes a title in a figure. A font chosen in the editor may have no small lowered or raised digits (the ₀ of ΔF/F₀, the ²⁺ of
   *  Ca²⁺) and would borrow them from another font at the wrong size and place: they are drawn from its ordinary digits instead.
   *  The page's own font has them, so the page writes its text as it always did. */
  function writeText(ctx, text, x, y, look) {
    if (!look.family || !/[₀-₉₊₋⁰¹²³⁴-⁹⁺⁻]/.test(text)) { ctx.fillText(text, x, y); return; }
    const font = ctx.font, size = parseFloat(/([\d.]+)px/.exec(font)[1]), small = font.replace(/[\d.]+px/, (size * 0.68) + 'px'), runs = [];
    for (const letter of text) {
      const low = '₀₁₂₃₄₅₆₇₈₉₊₋'.indexOf(letter), high = '⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻'.indexOf(letter), level = low >= 0 ? -1 : high >= 0 ? 1 : 0, shown = level ? '0123456789+−'[Math.max(low, high)] : letter;
      if (runs.length && runs[runs.length - 1].level === level) runs[runs.length - 1].text += shown; else runs.push({ level, text: shown });
    }
    let width = 0; for (const run of runs) { ctx.font = run.level ? small : font; run.width = ctx.measureText(run.text).width; width += run.width; }
    const align = ctx.textAlign; let at = align === 'center' ? x - width / 2 : align === 'right' ? x - width : x; ctx.textAlign = 'left';
    for (const run of runs) { ctx.font = run.level ? small : font; ctx.fillText(run.text, at, y - (run.level > 0 ? 0.33 : run.level < 0 ? -0.2 : 0) * size); at += run.width; }
    ctx.font = font; ctx.textAlign = align;
  }
  /** Ticks at a spacing typed in the figure editor; the automatic ones when none is typed or it would give more than TICKS_MOST
   *  (the editor refuses such a spacing: figureRefusal). */
  const TICKS_MOST = 60;
  const steppedTicks = (lo, hi, step, auto) => { if (!finite(step) || step <= 0 || (hi - lo) / step > TICKS_MOST) return auto;
    const values = []; for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) values.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12)); return { values, step }; };
  /** The numbers along a time axis W pixels long by a look, one rule for the time course and the heatmap, on the page and in a
   *  figure: one about every 90 pixels of 12-pixel text, and no more than TIME_TICKS, so that the two pictures of one recording
   *  are numbered alike (a spacing typed in the editor is kept). */
  const TIME_TICKS = 6;
  const timeTicks = (x0, x1, W, look) => steppedTicks(x0, x1, look.xStep, ticks(x0, x1, Math.max(2, Math.min(TIME_TICKS, Math.floor(W / (90 * look.px / 12))))));
  /** Strokes one trace: its samples joined in time, with a break where a sample is missing or the recording paused. */
  function strokeTrace(ctx, time, values, X, Y, gapLimit) {
    ctx.beginPath(); let pen = false;
    for (let i = 0; i < time.length && i < values.length; i++) {
      if (finite(values[i])) { if (pen && time[i] - time[i - 1] <= gapLimit) ctx.lineTo(X(time[i]), Y(values[i])); else ctx.moveTo(X(time[i]), Y(values[i])); pen = true; } else pen = false;
    }
    ctx.stroke();
  }
  /** Writes the stimulus label above its line at x, kept inside the width of the plot area (L, W) and no wider than it. Besides its axis titles and
   *  tick labels, this one word is all a painter writes on a plot; every other name is an entry of the legend. */
  function stimulusLabel(ctx, look, x, L, W, y) {
    const text = fitted(ctx, look.stimulusText, W), half = ctx.measureText(text).width / 2; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillStyle = look.colourStimulus;
    ctx.fillText(text, Math.max(L + half, Math.min(L + W - half, x)), y);
  }
  /** Where the parts of the plot go in a w × h area, in one step, from the top: the legend (rows starting on the Y axis, or a
   *  column at the right); a strip for what marks time without lying on the traces (the stimulus label, and under it the
   *  response-window bar, clear of the plot's top edge); the plot area, with a left margin that fits its tick labels.
   *  ctx carries the font; items are the legend's entries and tickLabels the Y axis numbers.
   *  Gaps and margins are written for 12-pixel text and grow with the text; line widths and tick lengths come from the look.
   *  A column of legend entries takes at most a third of the width. The plot area gets what the rest leaves, and never less
   *  than 10 pixels each way, so a figure too small for its text still draws; the editor says when it is (FIT). */
  function plotLayout(ctx, w, h, look, parts, items, tickLabels) {
    const u = (n) => n * look.px / 12, widest = (texts) => texts.reduce((most, text) => Math.max(most, ctx.measureText(text).width), 0);
    const column = look.legend === 'right' && items.length ? legendLayout(ctx, items, w / 3, look, true) : null;
    const L = Math.ceil(u(6) + look.px + u(6) + widest(tickLabels) + u(4) + look.tickLength);
    const R = Math.ceil(u(14) + (column ? u(12) + column.width + u(8) : 0)), W = Math.max(10, w - L - R);
    const rows = look.legend === 'top' && items.length ? legendLayout(ctx, items, w - L - u(6), look, false) : null, top = u(8), under = top + (rows ? rows.height + u(4) : 0);
    const labelled = parts.stimulus !== null && !!look.stimulusText, barred = !!parts.spans.response;
    const labels = labelled ? u(15) : 0, T = Math.round(under + labels + (barred ? look.windowBar + u(4) : labelled ? u(2) : 0));
    const below = look.tickLength + u(3) + look.px + u(5), H = Math.max(10, Math.floor(h - T - below - look.px - u(5)));
    return { L, T, W, H, labelY: under + u(7), lineTop: labelled ? under + labels : T, barY: barred ? T - u(4) - look.windowBar : null, titleY: T + H + below,
      legend: column ? { layout: column, x: w - u(8) - column.width, y: T } : rows ? { layout: rows, x: L, y: top } : null };
  }
  /** Everything about the time course that is settled before a line is drawn, in a w × h area by a look: the time and the range
   *  shown, what is drawn in them (plotParts), the legend's entries, the Y axis numbers and where the parts go (g: plotLayout).
   *  ctx is only measured with, in the look's font. The painter draws from this; the figure editor asks it how much room
   *  the plot would get (figureFit). */
  function plotPlan(ctx, w, h, look) {
    const { x0, x1, full0, full1 } = plotTimes(look), parts = plotParts(look, x0, x1), { lo, hi } = plotLimits(parts, look);
    parts.levels = parts.levels.filter((level) => level.value >= lo && level.value <= hi);   // a threshold that typed limits leave outside the plot is not drawn, and so not listed
    const items = plotLegend(parts, look);
    ctx.font = look.px + 'px ' + (look.family || getComputedStyle(document.body).fontFamily);
    // about five numbers on the Y axis, or as many as its height holds at a line and a half each (a spacing typed in the editor is kept)
    const laid = (about) => { const ys = steppedTicks(lo, hi, look.yStep, ticks(lo, hi, about)); return { ys, g: plotLayout(ctx, w, h, look, parts, items, ys.values.map((v) => tickText(v, ys.step))) }; };
    let plan = laid(5); const held = Math.floor(plan.g.H / (look.px * 1.5)); if (held < 5) plan = laid(Math.max(1, held));
    return Object.assign(plan, { x0, x1, full0, full1, parts, lo, hi, items });
  }
  /** Draws the time course into a w × h area by a look: the page's (PAGE_LOOK) or a saved figure's, the same drawing for both.
   *  On the page this is everything but the time cursor and the hovered trace, drawn once and reused while the cursor moves.
   *  The only words it writes are the axis titles, the tick labels, the stimulus label and the legend.
   *  Returns where things went (plotLayout), the axis ranges, the time ticks, the legend's entries and the means drawn. */
  function paintTraces(ctx, w, h, look) {
    const u = (n) => n * look.px / 12, r = S.results, shown = shownTraces(), time = r.time, gapLimit = samplingSupport().gap;
    const { x0, x1, full0, full1, parts, lo, hi, items, ys, g } = plotPlan(ctx, w, h, look), spans = parts.spans, { L, T, W, H } = g;
    const xs = timeTicks(x0, x1, W, look);
    const X = (t) => L + (t - x0) / (x1 - x0) * W, Y = (v) => T + (hi - v) / (hi - lo) * H, crisp = (at) => Math.round(at) + 0.5;
    if (look.background) { ctx.fillStyle = look.background; ctx.fillRect(0, 0, w, h); }
    const inside = (draw) => { ctx.save(); ctx.beginPath(); ctx.rect(L, T, W + 1, H + 1); ctx.clip(); draw(); ctx.restore(); };
    const across = (y, colour, width, dash) => { ctx.strokeStyle = colour; ctx.lineWidth = width; ctx.lineCap = 'butt'; ctx.setLineDash(dash || []); ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + W, y); ctx.stroke(); ctx.setLineDash([]); };
    /** A span of time as a band over the whole height of the plot, at least "least" wide (from where it starts, kept inside the plot). */
    const band = (span, colour, least) => { const b = Math.min(L + W, Math.max(X(span[1]), X(span[0]) + least)), a = Math.max(L, Math.min(X(span[0]), b - least)); ctx.fillStyle = colour; ctx.fillRect(a, T, b - a, H); };
    const line = (values) => strokeTrace(ctx, time, values, X, Y, gapLimit);
    // behind the data: the baseline band and the run that counted; grid lines when wanted; a dotted line at zero
    const around = parts.around;
    inside(() => {
      if (spans.baseline) band(spans.baseline, look.colourBaseline, 0);
      if (spans.run) band(spans.run, look.colourRun, look.runWidth);
      if (around && around.band) band(around.band, look.colourFlagged, look.runWidth);
      if (look.grid) for (const v of ys.values) across(crisp(Y(v)), look.colourGrid, look.axisWidth);
      // (not under a threshold line that lies on zero, the baseline mean of ΔF/F₀: two dash patterns on one line read as neither)
      if (lo < 0 && hi > 0 && (isDff() || S.signal === 'rawdff' || around) && !parts.levels.some((level) => Math.abs(Y(level.value) - Y(0)) < 1)) across(crisp(Y(0)), look.colourZero, look.axisWidth, look.zeroDash);
    });
    // marks of time, in the strip above the plot: the response window, a thin bar; the stimulus, a thin dashed line with its
    // label above it. Nothing of this lies on the traces.
    if (spans.response) { const a = Math.max(L, X(spans.response[0])), b = Math.min(L + W, X(spans.response[1])); ctx.fillStyle = look.colourWindow; ctx.fillRect(a, g.barY, b - a, look.windowBar); }
    if (parts.stimulus !== null) {
      const sx = crisp(X(parts.stimulus)); ctx.strokeStyle = look.colourStimulus; ctx.lineWidth = look.stimulusWidth; ctx.lineCap = 'butt'; ctx.setLineDash(look.stimulusDash);
      ctx.beginPath(); ctx.moveTo(sx, g.lineTop); ctx.lineTo(sx, T + H); ctx.stroke(); ctx.setLineDash([]);
      if (look.stimulusText) stimulusLabel(ctx, look, sx, L, W, g.labelY);
    }
    // the data: single cells (the others first, faint; then the responders, a little wider and stronger, so they are not
    // hidden), SEM bands, the selected cell's thresholds, the means, and the selected cell on top; last, a tick at the top
    // of the plot where the run that counted starts
    inside(() => {
      ctx.lineJoin = 'round'; ctx.lineCap = 'round';
      if (parts.traces.length) {
        const responded = (name) => callOf(name) === 'Responded', yes = parts.traces.filter(responded), no = parts.traces.filter((name) => !responded(name));
        const all = parts.traces.length;   // both kinds faded by every trace drawn: a lone responder among a hundred cells is as faint as they are, only a little stronger
        ctx.strokeStyle = look.colourTrace; ctx.lineWidth = look.traceWidth; ctx.globalAlpha = look.traceOpacity == null ? traceFade(all) : look.traceOpacity; for (const name of no) line(shown[name]);
        ctx.strokeStyle = look.colourResponded; ctx.lineWidth = look.responderWidth; ctx.globalAlpha = look.responderOpacity == null ? responderFade(all) : look.responderOpacity; for (const name of yes) line(shown[name]);
        ctx.globalAlpha = 1;
      }
      for (const s of parts.series) if (bandDrawn(s, look)) {   // one shape per unbroken stretch, so a band has no seams
        ctx.fillStyle = shade(look[s.colour], look.semOpacity); let run = [];
        const fill = () => { if (run.length > 1) { ctx.beginPath(); run.forEach((i, k) => { const x = X(time[i]), y = Y(s.mean[i] + s.sem[i]); if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
          for (let k = run.length - 1; k >= 0; k--) ctx.lineTo(X(time[run[k]]), Y(s.mean[run[k]] - s.sem[run[k]])); ctx.closePath(); ctx.fill(); } run = []; };
        for (let i = 0; i < time.length; i++) { const has = finite(s.mean[i]) && finite(s.sem[i]); if (!has || (run.length && time[i] - time[i - 1] > gapLimit)) fill(); if (has) run.push(i); }
        fill();
      }
      for (const level of parts.levels) across(Y(level.value), look[level.colour], look[level.width], look.levelDash);
      ctx.lineCap = 'round';
      const cased = (values, colour, width, casing) => { if (casing) { ctx.strokeStyle = look.colourCasing; ctx.lineWidth = width + 2 * look.casing; line(values); } ctx.strokeStyle = colour; ctx.lineWidth = width; line(values); };
      for (const s of parts.series) cased(s.mean, look[s.colour], look.meanWidth, parts.traces.length > 0);
      if (parts.own) cased(parts.own, look.colourCell, look.cellWidth, true);
      if (around) {   // the other ring sectors, thin, under the zones of the flag; the cell on top, as on the other plots
        const zone = (values) => strokeTrace(ctx, around.time, values, X, Y, gapLimit);
        ctx.strokeStyle = look.colourOtherSectors; ctx.lineWidth = look.otherWidth; for (const other of around.others) zone(other.values);
        for (const one of around.lines) { ctx.strokeStyle = zoneInk(look, one); ctx.lineWidth = look.zoneWidth; zone(one.values); }
        if (around.cell) { ctx.strokeStyle = look.colourCasing; ctx.lineWidth = look.cellWidth + 2 * look.casing; zone(around.cell); ctx.strokeStyle = look.colourCell; ctx.lineWidth = look.cellWidth; zone(around.cell); }
      }
      if (spans.run && spans.run[0] >= x0 && spans.run[0] <= x1) { const sx = X(spans.run[0]); ctx.strokeStyle = look.colourLevelRequired; ctx.lineWidth = look.requiredWidth; ctx.lineCap = 'butt'; ctx.beginPath(); ctx.moveTo(sx, T); ctx.lineTo(sx, T + look.runTick); ctx.stroke(); }
    });
    // axes: left and bottom (or a box), ticks pointing out, numbers with a true minus sign, all in the one ink
    const left = L + 0.5, bottom = T + H + 0.5; ctx.strokeStyle = ctx.fillStyle = look.ink; ctx.lineWidth = look.axisWidth; ctx.lineCap = 'butt'; ctx.lineJoin = 'miter'; ctx.beginPath();
    if (look.frame === 'box') ctx.rect(left, T + 0.5, W, H); else { ctx.moveTo(left, T); ctx.lineTo(left, bottom); ctx.lineTo(L + W, bottom); }
    for (const v of ys.values) { ctx.moveTo(left, crisp(Y(v))); ctx.lineTo(left - look.tickLength, crisp(Y(v))); }
    for (const v of xs.values) { ctx.moveTo(crisp(X(v)), bottom); ctx.lineTo(crisp(X(v)), bottom + look.tickLength); }
    ctx.stroke();
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; for (const v of ys.values) ctx.fillText(tickText(v, ys.step), L - look.tickLength - u(4), crisp(Y(v)));
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; for (const v of xs.values) ctx.fillText(tickText(v, xs.step), crisp(X(v)), T + H + look.tickLength + u(3));
    writeText(ctx, centred(ctx, look.xTitle == null ? 'Time (s)' : look.xTitle, L + W / 2, w, look), L + W / 2, g.titleY, look);
    ctx.save(); ctx.translate(u(6) + look.px / 2, T + H / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'middle'; writeText(ctx, centred(ctx, look.yTitle == null ? plotYTitle() : look.yTitle, T + H / 2, h, look), 0, 0, look); ctx.restore();
    if (g.legend) drawLegend(ctx, g.legend.layout, g.legend.x, g.legend.y, look);
    const drew = Object.assign(g, { x0, x1, lo, hi, full0, full1, items, times: xs.values, plotted: around ? items.map((item) => item.text) : parts.series.map((s) => meanLabel(s, bandDrawn(s, look))), around });
    if (look.page) { Plot.geom = drew; Plot.drawn = parts.traces; Plot.series = parts.series; Plot.around = around; }   // the page's cursor and hover read these; a figure leaves them alone
    return drew;
  }
  function drawPlot() {
    const canvas = $('plot'), r = S.results, around = r && S.plotMode === 'around' ? aroundNow() : null, blocked = !!around && !around.data;
    // "Around the cell" with nothing to show for the cell selected: why, in Fiji's words, where the plot would be, and nothing to save
    $('plot-empty').hidden = !!r && !blocked; $('plot-empty').textContent = blocked ? around.why || '' : EMPTY_PLOT;
    if (r) for (const id of ['save-plot', 'save-plot-folder', 'edit-plot']) $(id).disabled = blocked;
    // What the page says beside the plot, not in it, set before the canvas is measured (the note takes a line from it): that the
    // selected cell has no trace to draw; and "Thresholds" is there only where there are thresholds, with one cell in ΔF/F₀
    // (not a greyed box in the views where it would do nothing). "Around the cell" is in gray levels and draws no other cells.
    const missing = !!r && traceMissing(); $('cell-note').hidden = !missing; if (missing) $('cell-note').textContent = 'No ΔF/F₀ trace' + (hasTable('corrected') ? ': choose gray levels under Signal' : '');
    $('thresholds-label').hidden = !(r && S.plotMode === 'cell' && isDff());
    if (r) { $('every-label').hidden = $('signal-label').hidden = S.plotMode === 'around'; showAroundButton(); }
    segRows();
    if (!canvas.clientWidth) return;
    const { ctx, w, h, dpr } = fit(canvas);
    if (!r || blocked) { ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); if (blocked) { Plot.geom = null; Plot.around = null; Plot.dirty = true; canvas.setAttribute('aria-label', $('plot-empty').textContent); } return; }
    if (Plot.dirty || Plot.layer.width !== canvas.width || Plot.layer.height !== canvas.height) {
      Plot.layer.width = canvas.width; Plot.layer.height = canvas.height; const lctx = Plot.layer.getContext('2d'); lctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      paintTraces(lctx, w, h, PAGE_LOOK); Plot.dirty = false;
      // A limit typed under "Y axis" that the plot does not use marks its box: text that is no number, or a minimum at or above the maximum (typedRange).
      for (const [id, end, why] of [['y-min', 'lo', 'Minimum must be below maximum'], ['y-max', 'hi', 'Maximum must be above minimum']]) { const box = $(id), value = typedNumber(box.value), unused = value !== null && Plot.geom[end] !== value;
        if (unused) { box.setAttribute('aria-invalid', 'true'); box.title = Number.isNaN(value) ? 'Not a number' : why; } else { box.removeAttribute('aria-invalid'); box.removeAttribute('title'); } }
      canvas.setAttribute('aria-label', 'Plot of ' + plotYTitle() + ' against time: ' + Plot.geom.items.map((item) => item.text).join('; '));   // the legend is drawn in the plot: said here for a screen reader
    }
    ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.drawImage(Plot.layer, 0, 0); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const g = Plot.geom, time = r.time, X = (t) => g.L + (t - g.x0) / (g.x1 - g.x0) * g.W, Y = (v) => g.T + (g.hi - v) / (g.hi - g.lo) * g.H;
    ctx.save(); ctx.beginPath(); ctx.rect(g.L, g.T, g.W + 1, g.H + 1); ctx.clip();
    if (Plot.hover && S.plotMode !== 'around' && shownTraces()[Plot.hover] && Plot.hover !== S.cell) {   // the trace under the pointer, in the colour of its call
      ctx.strokeStyle = CALL[kind(callOf(Plot.hover))]; ctx.lineWidth = 2.2; ctx.lineJoin = 'round'; strokeTrace(ctx, time, shownTraces()[Plot.hover], X, Y, samplingSupport().gap);
    }
    // Nothing sits on the traces while nobody is reading it. While the frame moves (a drag in the plot, the slider or the step
    // buttons under the picture, the arrow keys, Play) and for a moment after, a solid line shows where the picture is, with
    // that frame's numbers. Otherwise a dashed line and the numbers follow the pointer over the plot and go when it leaves.
    // "Around the cell" says the cell's rise and that of each zone of its flag, in gray levels.
    const near = S.plotMode === 'around' ? Plot.around : null, lines = near ? near.lines : [];
    const own = near ? near.cell : S.plotMode === 'cell' && S.cell ? shownTraces()[S.cell] : null, moving = Plot.scrubbing || !!S.timer, awake = moving || performance.now() - Plot.movedAt < LINGER;
    Plot.readout = null; Plot.line = false;
    const plotted = (i) => (S.plotMode !== 'cell' ? Plot.series.map((series) => series.mean[i]) : []).concat(own ? [own[i]] : [], lines.map((line) => line.values[i]));
    const numbers = (at, i) => { let text = at.toFixed(1) + ' s'; if (S.plotMode !== 'cell') for (const series of Plot.series) if (finite(series.mean[i])) text += ' · ' + series.short + ' ' + series.mean[i].toFixed(3);
      const digits = near ? 1 : 3; if (own && finite(own[i])) text += ' · ' + shortName(S.cell) + ' ' + own[i].toFixed(digits);
      for (const line of lines) if (finite(line.values[i])) text += ' · ' + line.place + ' ' + line.values[i].toFixed(digits); return text; };
    /** The numbers in a box at the line x: along the bottom of the plot, or along the top when a plotted value lies under the box there. */
    const label = (text, x, i) => {
      ctx.font = '12px ' + getComputedStyle(document.body).fontFamily; const tw = ctx.measureText(text).width + 10, bx = Math.max(g.L, Math.min(g.L + g.W - tw, x - tw / 2));
      const covers = (top) => plotted(i).some((v) => finite(v) && Y(v) > top - 4 && Y(v) < top + 22), bottom = g.T + g.H - 22, by = covers(bottom) && !covers(g.T + 8) ? g.T + 8 : bottom;
      ctx.fillStyle = '#1F2933'; ctx.beginPath(); ctx.roundRect(bx, by, tw, 18, 6); ctx.fill(); ctx.fillStyle = '#fff'; ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; ctx.fillText(text, bx + 5, by + 9); Plot.readout = text;
    };
    const t = Plot.playhead != null ? Plot.playhead : S.frame <= time.length ? time[S.frame - 1] : NaN;   // glides with the playhead while playing
    if (awake && isFinite(t) && t >= g.x0 && t <= g.x1) {
      const x = Math.round(X(t)) + 0.5; ctx.strokeStyle = 'rgba(31,41,51,0.8)'; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(x, g.T); ctx.lineTo(x, g.T + g.H); ctx.stroke(); Plot.line = true;
      if (own && finite(own[S.frame - 1])) { const y = Math.max(g.T, Math.min(g.T + g.H, Y(own[S.frame - 1]))); ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(x, y, 5, 0, 7); ctx.fill(); ctx.fillStyle = '#1F2933'; ctx.beginPath(); ctx.arc(x, y, 3.2, 0, 7); ctx.fill(); }
      if (moving || Plot.at == null) label(numbers(t, S.frame - 1), x, S.frame - 1);   // with the pointer on the plot, its own numbers are enough
    }
    const at = moving ? null : Plot.at;   // the sample under the pointer
    if (at != null && finite(time[at]) && time[at] >= g.x0 && time[at] <= g.x1) {
      const x = Math.round(X(time[at])) + 0.5; ctx.strokeStyle = 'rgba(31,41,51,0.4)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(x, g.T); ctx.lineTo(x, g.T + g.H); ctx.stroke(); ctx.setLineDash([]);
      label(numbers(time[at], at), x, at);
    }
    ctx.restore();
  }
  /** The plot buttons that do not fit on one line go on to the next, whole: each line of them gets its own rounded ends. */
  function segRows() {
    const buttons = Array.from(document.querySelectorAll('.plot-head .seg button'));
    buttons.forEach((button, i) => { const top = button.offsetTop;
      button.classList.toggle('row-start', i > 0 && buttons[i - 1].offsetTop !== top); button.classList.toggle('row-end', i < buttons.length - 1 && buttons[i + 1].offsetTop !== top); });
  }
  /** The drawn trace nearest to the mouse, when it is within a few pixels. */
  function traceAt(event) {
    const g = Plot.geom, r = S.results; if (!g || !r) return null;
    const box = $('plot').getBoundingClientRect(), mx = event.clientX - box.left, my = event.clientY - box.top;
    if (mx < g.L || mx > g.L + g.W || my < g.T || my > g.T + g.H) return null;
    const t = g.x0 + (mx - g.L) / g.W * (g.x1 - g.x0); let i = 0, gap = Infinity; for (let k = 0; k < r.time.length; k++) { const d = Math.abs(r.time[k] - t); if (d < gap) { gap = d; i = k; } }
    // ("Around the cell" draws one cell and its zones, in gray levels: no trace to pick, only the time under the pointer)
    const candidates = S.plotMode === 'around' ? [] : Plot.drawn.slice(); if (S.plotMode !== 'around' && S.cell && shownTraces()[S.cell] && !candidates.includes(S.cell)) candidates.push(S.cell);
    let best = null, near = 7;
    for (const name of candidates) { const v = shownTraces()[name][i]; if (!finite(v)) continue; const d = Math.abs(g.T + (g.hi - v) / (g.hi - g.lo) * g.H - my); if (d < near) { near = d; best = name; } }
    return { name: best, index: i, time: t, mx, my };
  }

  // ---------------------------------------------------------------- heatmap
  // Every cell in one picture, one row each: ΔF/F₀ as colour against time. Blue–white–red has white at zero and scales its two
  // sides separately; the other palettes, those of the SimpleCa²⁺ app, start at zero. Beside the rows a strip gives each row's
  // response call, and the key says what its colours mean and how many cells each holds; at the right, a colour bar. Nothing
  // else is written beside the rows but the names of the cells, when they are asked for. One painter for the page and for a
  // saved figure.
  const PALETTES = {
    plasma: ['#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786', '#d8576b', '#ed7953', '#fb9f3a', '#fdca26', '#f0f921'],
    viridis: ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'],
    magma: ['#000004', '#180f3d', '#440f76', '#721f81', '#9e2f7f', '#cd4071', '#f1605d', '#fd9668', '#feca8d', '#fcfdbf'],
    inferno: ['#000004', '#1b0c41', '#4a0c6b', '#781c6d', '#a52c60', '#cf4446', '#ed6925', '#fb9b06', '#f7d13d', '#fcffa4'],
    cividis: ['#00224e', '#123570', '#3b496c', '#575d6d', '#707173', '#8a8678', '#a59c74', '#c3b369', '#e1cc55', '#fee838'],
    signed: ['#3B4CC0', '#FFFFFF', '#B40426']
  };
  const rgb = (hex) => [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)];
  /** The colours along a palette, in RAMP steps: an even number, so that the middle colour (the white of blue–white–red) is one of them. */
  const RAMP = 256;
  function ramp(name) {
    const stops = PALETTES[name].map(rgb), out = new Uint8Array((RAMP + 1) * 3);
    for (let i = 0; i <= RAMP; i++) { const at = i / RAMP * (stops.length - 1), k = Math.min(stops.length - 2, Math.floor(at)), f = at - k;
      for (let c = 0; c < 3; c++) out[i * 3 + c] = Math.round(stops[k][c] + (stops[k + 1][c] - stops[k][c]) * f); }
    return out;
  }
  /** The largest value after the baseline frames and when it happens: the app's rule for the "Time to peak" and "Peak amplitude" orders. */
  function peakAfterBaseline(trace, time, lastBaselineFrame) {
    let best = -Infinity, at = Infinity; if (!trace) return { value: -Infinity, time: Infinity };
    for (let i = Math.max(0, lastBaselineFrame); i < trace.length; i++) if (finite(trace[i]) && trace[i] > best) { best = trace[i]; at = time[i]; }
    return { value: best, time: at };
  }
  /** The entry of the style table that colours a row's response call, by the kind of call ("kind"). */
  const CALL_STYLE = { responded: 'colourResponded', not: 'colourNot', other: 'colourOther' };
  /** The cells the heatmap uses. The one place that decides it, for the rows, their groups and the automatic colour range:
   *    rows    every measured cell in the order chosen under Sort, without the excluded cells when they are hidden. In every
   *            order but the original one, a cell without a trace comes after those with one, and a cell below zero last
   *            among those (its ΔF/F₀ is divided by an F₀ near zero, so its size and timing rank nothing).
   *            "Response call" and "Onset time" keep the calls together (responded, did not respond, the rest), and the rule
   *            holds inside each group; the largest rise in the response window (the call table's Peak_increase) is on top.
   *            "Onset time" puts the responders in the order of their onsets and leaves the rest as "Response call" has them.
   *    groups  those groups as runs of rows, { kind, first, n }; none in the other orders, where the calls are mixed.
   *    steady  the cells an automatic colour range is taken from: neither excluded nor below zero (plotCells). */
  function heatmapCells() {
    const r = S.results, all = S.info.names, place = new Map(all.map((name, i) => [name, i])), original = (a, b) => place.get(a) - place.get(b);
    const rows = all.filter((name) => !(S.hideExcluded && callOf(name) === 'Excluded')), which = (name) => kind(callOf(name)), grouped = S.order === 'call' || S.order === 'onset';
    if (S.order !== 'orig') {
      const rank = (name) => (grouped ? Object.keys(CALL_STYLE).indexOf(which(name)) : 0), bare = (name) => (r.traces[name] ? 0 : 1), below = (name) => (belowZero(name) ? 1 : 0);
      const onset = (name) => { const v = number(callRow(name).Onset_s); return S.order === 'onset' && which(name) === 'responded' && finite(v) ? v : Infinity; };
      const rise = (name) => { const v = number(callRow(name).Peak_increase); return finite(v) ? v : -Infinity; };
      const last = r.baselineFrames && finite(r.baselineFrames[1]) ? r.baselineFrames[1] : 0, peaks = grouped ? null : new Map(rows.map((name) => [name, peakAfterBaseline(r.traces[name], r.time, last)]));
      const key = grouped ? (a, b) => rise(b) - rise(a) : S.order === 'tpeak' ? (a, b) => peaks.get(a).time - peaks.get(b).time : (a, b) => peaks.get(b).value - peaks.get(a).value;
      // (two cells without the number differ by "not a number", which counts as a tie and passes to the next key)
      rows.sort((a, b) => rank(a) - rank(b) || onset(a) - onset(b) || bare(a) - bare(b) || below(a) - below(b) || key(a, b) || original(a, b));
    }
    const groups = [];
    if (grouped) rows.forEach((name, c) => { const now = groups[groups.length - 1]; if (now && now.kind === which(name)) now.n++; else groups.push({ kind: which(name), first: c, n: 1 }); });
    return { rows, groups, steady: plotCells(0, r.traces).steady };
  }
  /** The heatmap's colour range, from the values of the steady cells alone: the cells left out of it (excluded, or below zero) are
   *  the few that would otherwise set it, since one of them can hold more extreme samples than a trim of all cells removes.
   *  "1st to 99th percentile": of those values. "Smallest to largest": their smallest and largest. Each end is rounded outward
   *  to a tidy number of its own, so the two sides of blue–white–red scale separately; a palette that starts at zero starts at
   *  zero. A limit typed in the figure editor replaces its end (and moves the other end when that would not lie beyond it);
   *  two typed limits that are not a range are not used. */
  function heatmapRange(steady, look) {
    const T = S.results.traces, values = Float64Array.from(steady.flatMap((name) => T[name].filter(finite))).sort(), whole = S.range === 'full';
    const at = (part) => (values.length ? values[Math.round(part * (values.length - 1))] : 0), bottom = at(whole ? 0 : 0.01), top = at(whole ? 1 : 0.99);
    const tidy = (v) => { const step = ticks(0, v, 5).step; return +(Math.ceil(v / step - 1e-9) * step).toPrecision(12); };
    let low = S.palette === 'signed' && bottom < 0 ? -tidy(-bottom) : 0, high = top > 0 ? tidy(top) : 1;
    const span = high - low, a = look.barMin, b = look.barMax;
    if (finite(a) && finite(b)) { if (b > a) { low = a; high = b; } }
    else if (finite(a)) { low = a; if (!(high > low)) high = low + span; }
    else if (finite(b)) { high = b; if (!(high > low)) low = high - span; }
    return { low, high };
  }
  /** What the heatmap draws, gathered once for the painter and its key:
   *    names, groups   the rows and their groups (heatmapCells)       flagged  the rows of cells below zero
   *    x0, x1          the time shown: the whole recording, or the limits typed in the figure editor
   *    low, high       the colour range (heatmapRange)                 colourAt the colour of a value
   *    below, above    whether a value that is drawn lies beyond that end of the range: the colour bar then says ≤ or ≥
   *    calls           the kinds of call among the rows, with what the key calls each and its number of rows
   *    spans, stimulus the baseline, the response window and when the stimulus was added (in results saved without that time,
   *                    the start of the response window)
   *    small           the picture itself at one pixel per sample and cell; a missing value is left clear. */
  function heatmapData(look) {
    const r = S.results, time = r.time, cells = heatmapCells(), names = cells.rows, signed = S.palette === 'signed', spans = responseSpans(), onset = number(selectedCall().Stimulus_onset_s);
    const [x0, x1] = typedRange(wholeTime(), look.xMin, look.xMax) || wholeTime();
    const { low, high } = heatmapRange(cells.steady, look), colours = ramp(S.palette);
    // Where a value lies along the palette, from 0 to 1. Blue–white–red: white at zero, the full blue at the lower limit and
    // the full red at the upper one, each side on its own scale. The others: evenly from the lower limit to the upper.
    const along = signed ? (v) => { const c = Math.max(low, Math.min(high, v)); return c >= 0 ? 0.5 + (high > 0 ? 0.5 * c / high : 0) : 0.5 - (low < 0 ? 0.5 * c / low : 0); }
      : (v) => Math.max(0, Math.min(1, (v - low) / (high - low)));
    const colourAt = (v) => { const k = Math.round(along(v) * RAMP) * 3; return [colours[k], colours[k + 1], colours[k + 2]]; };
    const small = document.createElement('canvas'); small.width = Math.max(1, time.length); small.height = Math.max(1, names.length);
    const sctx = small.getContext('2d'), data = sctx.createImageData(small.width, small.height), slack = (high - low) * 1e-9; let below = false, above = false;
    names.forEach((name, c) => { const trace = r.traces[name]; if (!trace) return;
      for (let f = 0; f < time.length; f++) {
        const v = trace[f]; if (!finite(v)) continue; const o = (c * time.length + f) * 4, k = Math.round(along(v) * RAMP) * 3;
        data.data[o] = colours[k]; data.data[o + 1] = colours[k + 1]; data.data[o + 2] = colours[k + 2]; data.data[o + 3] = 255;
        if (time[f] >= x0 && time[f] <= x1) { if (v < low - slack) below = true; else if (v > high + slack) above = true; }
      } });
    sctx.putImageData(data, 0, 0);
    const said = { Excluded: 'excluded', Unscorable: 'unscorable' }, others = Array.from(new Set(names.filter((name) => kind(callOf(name)) === 'other').map((name) => said[callOf(name)] || 'no call'))).sort().join(' or ');
    const calls = [['responded', 'responded'], ['not', 'did not respond'], ['other', others]].map(([k, text]) => ({ kind: k, text, n: names.filter((name) => kind(callOf(name)) === k).length })).filter((call) => call.n > 0);
    return { names, groups: cells.groups, flagged: new Set(names.filter(belowZero)), signed, x0, x1, low, high, below, above, colourAt, calls, spans, stimulus: finite(onset) ? onset : spans.response[0], small };
  }
  /** Whether a span of time lies, at least in part, in the time the heatmap shows. */
  const spanShown = (d, span) => realSpan(span) && span[1] > d.x0 && span[0] < d.x1;
  /** The heatmap's key, in one place: what the colours of the strip mean (the calls that are drawn, each with its number of
   *  cells), the mark of a cell below zero, and the two bars above the picture, each only when it is drawn. Entries as
   *  drawLegend draws them. */
  function heatmapKey(d, look) {
    const items = [];
    if (look.strip) { for (const call of d.calls) items.push({ kind: 'box', colour: look[CALL_STYLE[call.kind]], outline: true, text: counted(call.text, call.n) });
      if (d.flagged.size) items.push({ kind: 'box', colour: look.colourFlag, text: 'below zero' }); }
    if (look.marks && spanShown(d, d.spans.baseline)) items.push({ kind: 'box', colour: look.colourBaseline, height: look.windowBar, text: 'baseline' });
    if (look.marks && spanShown(d, d.spans.response)) items.push({ kind: 'box', colour: look.colourWindow, height: look.windowBar, text: 'response window' });
    return items;
  }
  /** Where the parts of the heatmap go in a w × h area, in one step, sized from the text they hold. From the top: the key (rows
   *  starting at the left edge, since most of it is about what stands left of the picture); a strip for the stimulus label and the colour bar's title, with the bars of the baseline
   *  and the response window under it; the picture. Left of the picture, outward (each as its distance from the picture's edge):
   *  the strip of calls, the mark of a cell below zero, the cell names with their ticks, the axis title. Right of it: the colour bar, as tall as the picture, with its numbers: both ends, each with ≤ or ≥ when values
   *  lie beyond it, and the round values between that are clear of them. ctx gets the look's font, for measuring here and for
   *  the painter's own text. Gaps are written for 12-pixel text and grow with the text; widths come from the look. */
  function heatmapLayout(ctx, w, h, look, d) {
    ctx.font = look.px + 'px ' + (look.family || getComputedStyle(document.body).fontFamily);
    const u = (n) => n * look.px / 12, wide = (text) => ctx.measureText(text).width, widest = (texts) => texts.reduce((most, text) => Math.max(most, wide(text)), 0);
    const yTitle = look.yTitle == null ? 'Cells' : look.yTitle, barTitle = fitted(ctx, look.barTitle == null ? 'ΔF/F₀' : look.barTitle, w / 4), g = {}; let at = 0;   // (the colour bar's title makes room for itself, up to a quarter of the width)
    if (look.strip) { at += u(2) + look.stripWidth; g.strip = at; if (d.flagged.size) { at += u(1) + look.flagWidth; g.flag = at; } }
    if (S.rowNames) { at += u(1); g.tick = at; at += look.tickLength + u(2); g.names = at; at += widest(d.names.map(shortName)); }
    g.beside = at;   // how far the strip and the names reach: a click there picks the row too
    if (yTitle) { at += u(6) + look.px / 2; g.title = at; at += look.px / 2; }
    const L = Math.ceil(at + u(6)), items = heatmapKey(d, look), key = items.length ? legendLayout(ctx, items, w - u(14), look, false) : null, top = u(8), under = top + (key ? key.height + u(4) : 0);
    const stimulus = look.marks && finite(d.stimulus) && d.stimulus >= d.x0 && d.stimulus <= d.x1;
    const labelled = stimulus && !!look.stimulusText, barred = look.marks && (spanShown(d, d.spans.baseline) || spanShown(d, d.spans.response)), titled = look.bar && !!barTitle;
    const labels = labelled ? u(15) : 0, T = Math.round(under + Math.max(labels + (barred ? look.windowBar + u(4) : labelled ? u(2) : 0), titled ? u(9) + look.px : 0));
    const below = look.tickLength + u(3) + look.px + u(5), H = Math.max(10, Math.floor(h - T - below - look.px - u(5))), Y = (v) => T + (d.high - v) / (d.high - d.low) * H;
    // the colour bar's numbers: an end is written in full, whatever the spacing of the round values between
    const between = ticks(d.low, d.high, Math.max(1, Math.floor(H / u(60)))), edge = (d.high - d.low) * 1e-9, clear = look.px * 2;   // a round number less than a line of text from an end is not written
    const exact = (v) => { let digits = between.step >= 1 ? 0 : Math.min(6, Math.ceil(-Math.log10(between.step) - 1e-9)); while (digits < 6 && Math.abs(+v.toFixed(digits) - v) > 1e-9 * Math.max(1, Math.abs(v))) digits++; return v.toFixed(digits).replace('-', '−'); };
    const numbers = [{ value: d.high, text: (d.above ? '≥ ' : '') + exact(d.high) }, { value: d.low, text: (d.below ? '≤ ' : '') + exact(d.low) }]
      .concat(between.values.filter((v) => v > d.low + edge && v < d.high - edge).map((value) => ({ value, text: Y(value) - T < clear || T + H - Y(value) < clear ? '' : tickText(value, between.step) })));
    const R = Math.ceil(look.bar ? u(12) + Math.max(look.barWidth + look.tickLength + u(3) + widest(numbers.map((n) => n.text)), titled ? wide(barTitle) : 0) + u(6) : u(14)), W = Math.max(10, w - L - R);
    return Object.assign(g, { L, T, W, H, yTitle, barTitle: titled ? barTitle : '', items, key: key ? { layout: key, x: u(8), y: top } : null, stimulus, labelY: under + u(7), lineTop: labelled ? under + labels : T,
      barsY: barred ? T - u(4) - look.windowBar : null, titleY: T + H + below, bar: look.bar ? { x: L + W + u(12), numbers } : null });
  }
  /** Draws the heatmap into a w × h area by a look: the page's (PAGE_HEAT) or a saved figure's, the same drawing for both.
   *  Returns where things went (heatmapLayout), the rows, the time shown and its ticks, and the colour range. */
  function paintHeatmap(ctx, w, h, look, d) {
    const u = (n) => n * look.px / 12, names = d.names, support = samplingSupport();
    const g = heatmapLayout(ctx, w, h, look, d), { L, T, W, H } = g, rowH = H / Math.max(1, names.length), rowY = (c) => T + c * rowH;
    const X = (t) => L + (t - d.x0) / (d.x1 - d.x0) * W, crisp = (at) => Math.round(at) + 0.5;
    if (look.background) { ctx.fillStyle = look.background; ctx.fillRect(0, 0, w, h); }
    // the picture, stretched without smoothing: each sample as wide as the frame times say, so a pause in the recording stays
    // blank; white lines between the groups
    ctx.fillStyle = look.colourMissing; ctx.fillRect(L, T, W, H);
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W, H); ctx.clip(); ctx.imageSmoothingEnabled = false;
    if (names.length) support.bounds.forEach((span, f) => { const x = X(span[0]), width = X(span[1]) - x; if (width > 0 && x < L + W && x + width > L) ctx.drawImage(d.small, f, 0, 1, names.length, x, T, width, H); });
    ctx.fillStyle = look.colourDivider; for (const group of d.groups.slice(1)) ctx.fillRect(L, rowY(group.first) - look.divider / 2, W, look.divider);
    ctx.restore();
    // marks of time: the baseline and the response window as thin bars above the picture; the stimulus as a thin dashed line,
    // its label above
    if (g.barsY !== null) for (const [span, colour] of [[d.spans.baseline, look.colourBaseline], [d.spans.response, look.colourWindow]]) {
      if (!realSpan(span)) continue; const a = Math.max(L, X(span[0])), b = Math.min(L + W, X(span[1])); if (b > a) { ctx.fillStyle = colour; ctx.fillRect(a, g.barsY, b - a, look.windowBar); }
    }
    if (g.stimulus) {
      const sx = crisp(X(d.stimulus)), dashed = (from, to, colour) => { ctx.strokeStyle = colour; ctx.beginPath(); ctx.moveTo(sx, from); ctx.lineTo(sx, to); ctx.stroke(); };
      ctx.lineWidth = look.stimulusWidth; ctx.lineCap = 'butt'; ctx.setLineDash(look.stimulusDash);
      dashed(g.lineTop, T, look.colourStimulus); dashed(T, T + H, d.signed ? look.colourOnLight : look.colourOnDark); ctx.setLineDash([]);
      if (look.stimulusText) stimulusLabel(ctx, look, sx, L, W, g.labelY);
    }
    // beside the rows: the calls as a strip, and a mark for the cells below zero. Each is one solid bar for each run of rows
    // alike, with no gaps between its rows, at least stripLeast high (centred on its rows, inside the picture's height); the
    // longer runs first, so that a short one lies on top of its neighbours.
    ctx.lineCap = 'butt'; ctx.lineJoin = 'miter'; ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    if (look.strip) {
      const runs = (of) => { const out = []; names.forEach((name, c) => { const k = of(name), last = out[out.length - 1]; if (k == null) return; if (last && last.k === k && last.end === c) last.end = c + 1; else out.push({ k, start: c, end: c + 1 }); });
        return out.sort((a, b) => (b.end - b.start) - (a.end - a.start)); };
      const bar = (run, x, wide, colour) => { const y0 = rowY(run.start), y1 = rowY(run.end), high = Math.min(H, Math.max(look.stripLeast, y1 - y0)), y = Math.max(T, Math.min(T + H - high, (y0 + y1 - high) / 2));
        ctx.fillStyle = colour; ctx.fillRect(x, y, wide, high); };
      for (const run of runs((name) => kind(callOf(name)))) bar(run, L - g.strip, look.stripWidth, look[CALL_STYLE[run.k]]);
      for (const run of runs((name) => (d.flagged.has(name) ? 'below' : null))) bar(run, L - g.flag, look.flagWidth, look.colourFlag);
    }
    // cell names: every row when the rows are tall enough for the text, else every Nth, each with a tick to its row
    const every = Math.max(1, Math.ceil(look.px * 1.1 / rowH)), named = [];
    if (S.rowNames) {
      ctx.strokeStyle = ctx.fillStyle = look.ink; ctx.lineWidth = look.axisWidth; ctx.beginPath();
      for (let c = 0; c < names.length; c += every) { const y = rowY(c + 0.5); ctx.moveTo(L - g.tick, crisp(y)); ctx.lineTo(L - g.tick - look.tickLength, crisp(y)); ctx.fillText(shortName(names[c]), L - g.names, y); named.push(names[c]); }
      ctx.stroke();
    }
    // on the page, the selected cell: an arrow beside its row, and a frame around it when the row is tall enough to keep its colours
    const chosen = look.page ? names.indexOf(S.cell) : -1;
    if (chosen >= 0) {
      const y = rowY(chosen), mid = y + rowH / 2; ctx.strokeStyle = ctx.fillStyle = look.colourSelected;
      if (rowH >= u(8)) { ctx.lineWidth = look.divider; ctx.strokeRect(L + look.divider / 2, y + look.divider / 2, W - look.divider, rowH - look.divider); }
      ctx.beginPath(); ctx.moveTo(L + W + u(2), mid); ctx.lineTo(L + W + u(9), mid - u(5)); ctx.lineTo(L + W + u(9), mid + u(5)); ctx.closePath(); ctx.fill();
    }
    // the time axis: a line under the picture, ticks pointing out, numbers with a true minus sign; then the two titles
    const bottom = T + H + 0.5, xs = timeTicks(d.x0, d.x1, W, look);
    ctx.strokeStyle = ctx.fillStyle = look.ink; ctx.lineWidth = look.axisWidth; ctx.beginPath(); ctx.moveTo(L, bottom); ctx.lineTo(L + W, bottom);
    for (const v of xs.values) { ctx.moveTo(crisp(X(v)), bottom); ctx.lineTo(crisp(X(v)), bottom + look.tickLength); }
    ctx.stroke();
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; for (const v of xs.values) ctx.fillText(tickText(v, xs.step), crisp(X(v)), T + H + look.tickLength + u(3));
    writeText(ctx, centred(ctx, look.xTitle == null ? 'Time (s)' : look.xTitle, L + W / 2, w, look), L + W / 2, g.titleY, look);
    if (g.yTitle) { ctx.save(); ctx.translate(L - g.title, T + H / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'middle'; writeText(ctx, centred(ctx, g.yTitle, T + H / 2, h, look), 0, 0, look); ctx.restore(); }
    // the colour bar: every colour of the range from the upper limit down, a tick at each number, the title above
    if (g.bar) {
      const x = g.bar.x, Y = (v) => T + (d.high - v) / (d.high - d.low) * H, right = x + look.barWidth;
      for (let i = 0; i < H; i++) { const c = d.colourAt(d.high - (i + 0.5) / H * (d.high - d.low)); ctx.fillStyle = 'rgb(' + c[0] + ',' + c[1] + ',' + c[2] + ')'; ctx.fillRect(x, T + i, look.barWidth, Math.min(1.5, H - i)); }
      ctx.strokeStyle = ctx.fillStyle = look.ink; ctx.lineWidth = look.axisWidth; ctx.strokeRect(x, T, look.barWidth, H);
      ctx.beginPath(); for (const n of g.bar.numbers) { const y = Y(n.value); ctx.moveTo(right, y); ctx.lineTo(right + look.tickLength, y); } ctx.stroke();
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; for (const n of g.bar.numbers) if (n.text) ctx.fillText(n.text, right + look.tickLength + u(3), Y(n.value));
      if (g.barTitle) writeText(ctx, g.barTitle, x, g.labelY, look);
    }
    if (g.key) drawLegend(ctx, g.key.layout, g.key.x, g.key.y, look);
    return Object.assign(g, { rows: names.length, names, named, groups: d.groups, bounds: support.bounds, x0: d.x0, x1: d.x1, times: xs.values, low: d.low, high: d.high, below: d.below, above: d.above });
  }
  function drawHeatmap() {
    const canvas = $('heatmap'), box = $('heatmap-box'), r = S.results; $('heatmap-empty').hidden = !!r; if (!r || $('panel-heatmap').hidden) return;
    canvas.style.width = Math.max(320, box.clientWidth - 2) + 'px'; canvas.style.height = Math.max(220, box.clientHeight - 2) + 'px';
    const { ctx, w, h } = fit(canvas), g = Heat.geom = paintHeatmap(ctx, w, h, PAGE_HEAT, heatmapData(PAGE_HEAT));
    // the key and the colour bar are drawn in the picture: said here for a screen reader
    canvas.setAttribute('aria-label', 'Heatmap of ΔF/F₀, one row for each of ' + plural(g.rows, 'cell') + ' against time. Colours from ' + g.bar.numbers[1].text + ' to ' + g.bar.numbers[0].text + '. ' + g.items.map((item) => item.text).join('; '));
  }
  function heatmapAt(event) {
    const g = Heat.geom; if (!S.results || !g) return null;
    const box = $('heatmap').getBoundingClientRect(), y = event.clientY - box.top, x = event.clientX - box.left, row = Math.floor((y - g.T) / g.H * g.rows);
    if (row < 0 || row >= g.rows || x < g.L - g.beside || x > g.L + g.W) return null;
    const wanted = g.x0 + Math.max(0, Math.min(1, (x - g.L) / g.W)) * (g.x1 - g.x0), f = nearestFrame(wanted) - 1;
    return { name: g.names[row], frame: f, inside: x >= g.L, gap: wanted < g.bounds[f][0] || wanted > g.bounds[f][1] };
  }

  // ---------------------------------------------------------------- tables: the saved tables in words, rounded for reading only
  /** The columns of a saved table as they are shown: its short set (every column with "Show all saved columns"), the Cell
   *  first, a column that is the same as another by definition once, the response call beside the cell in the per-cell metric
   *  tables, and the recording, role and stimulus only when they differ between rows. Each says how its values are written. */
  function tableColumns(kind, data, everything) {
    const has = (name) => data.header.includes(name);
    const describe = (saved, index) => {
      const words = columnWords(kind, saved), values = data.rows.map((row) => row[index]).filter((v) => v !== '' && v != null && v !== 'NA');
      const numeric = values.length > 0 && values.every((v) => isFinite(Number(v)));
      const decimals = words.decimals != null ? words.decimals : DECIMALS[words.unit] != null ? DECIMALS[words.unit] : values.every((v) => Number.isInteger(Number(v))) ? 0 : 3;
      // "Baseline noise" names the noise the row's Noise_source gives, when the table has it (as the Fiji window says it)
      return Object.assign({}, words, { index, numeric, decimals }, words.noise ? { noiseIndex: data.header.indexOf('Noise_source') } : {});
    };
    const all = data.header.map(describe), joined = FAMILY[kind] === 'cell' && has('Cell');
    let columns = everything || !SHORT[kind] ? all.slice() : SHORT[kind].map((name) => name === CALL_COLUMN ? (joined ? CALL_COLUMN : null) : all.find((c) => c.saved === name)).filter(Boolean);
    if (!everything) columns = columns.filter((c) => !c.alias || !has(c.alias));
    if (joined && !columns.includes(CALL_COLUMN)) columns.splice(columns.findIndex((c) => c.saved === 'Cell') + 1, 0, CALL_COLUMN);
    columns = [...columns.filter((c) => c.saved === 'Cell'), ...columns.filter((c) => c.saved !== 'Cell')];
    // the recording's name, the role and the stimulus repeat on every row of a one-window run: left out so the numbers fit
    if (data.rows.length > 1) columns = columns.filter((c) => !['Field', 'Role', 'Stimulus'].includes(c.saved) || new Set(data.rows.map((row) => row[c.index])).size > 1);
    return columns;
  }
  /** The response call of a row: its own in a table of calls, the primary stimulus's call of its cell in a per-cell table. */
  function callReader(kind, data) {
    const status = data.header.indexOf('Status'), cell = data.header.indexOf('Cell');
    if (FAMILY[kind] === 'window' && status >= 0) return (row) => row[status];
    if (FAMILY[kind] === 'cell' && cell >= 0 && S.results) return (row) => primaryCall(row[cell]);
    return () => undefined;
  }
  /** A saved value as a table shows it: { text, title, key (for sorting; '' sorts last), num (right-aligned) }. Numbers have the
   *  decimals of their column, the true minus, no "−0.000" and "%" after a percentage; TRUE and FALSE read yes and no; a missing
   *  value is one dash. The timing of a cell that did not respond is a dash too, with the saved value in its tooltip. */
  function shownCell(column, row, callOfRow) {
    if (column.call) { const call = callOfRow(row) || ''; return { text: call || '–', key: call }; }
    const raw = row[column.index] == null ? '' : String(row[column.index]);
    if (column.saved === 'Cell') return { text: shortName(raw) || '–', title: raw, key: shortName(raw) };
    if (raw === '' || raw === 'NA') return { text: '–', key: '', num: column.numeric };
    if (column.response && callOfRow(row) !== 'Responded') return { text: '–', title: 'Saved: ' + raw + ' · ' + (callOfRow(row) || 'no call'), key: '', num: column.numeric };
    if (column.numeric) { const v = Number(raw); return { text: decimal(column.fraction ? 100 * v : v, column.decimals, column.unit === PCT), title: raw, key: v, num: true }; }
    if (column.flags) return { text: raw.split(/[;|,]/).map((flag) => flag.trim()).filter(Boolean).map(flagShort).join(' · ') || '–', title: raw, key: raw };
    if (/^(true|false)$/i.test(raw)) return { text: /^true$/i.test(raw) ? 'yes' : 'no', title: raw, key: raw.toLowerCase() };
    if (column.noise && raw === 'Baseline noise' && column.noiseIndex >= 0) { const from = row[column.noiseIndex] == null ? '' : String(row[column.noiseIndex]);
      if (from !== '' && from !== 'NA') return { text: from, title: 'Saved: ' + raw + ' · Noise_source: ' + from, key: from }; }
    if (column.values && column.values[raw]) return { text: column.values[raw], title: 'Saved: ' + raw, key: column.values[raw] };
    return { text: raw, key: raw };
  }
  /** A column's name on screen, with its unit ("Peak, ΔF/F₀"), and what pointing at it says: what it is, and its saved name. */
  const headerText = (c) => c.given ? c.name : c.name + (NAMED_UNITS.has(c.unit) ? ', ' + c.unit : '');
  const headerTitle = (c) => c.given ? c.title || '' : [c.name + (c.measured ? ', ' + c.measured : ''), c.about,
    c.call ? '' : 'Saved as ' + c.saved + (c.fraction ? ', a fraction' : '') + (c.alias ? ', the same as ' + c.alias + ' by definition' : '')].filter(Boolean).join('. ') + '.';
  /** A table with one row for each window (the responder counts), turned: a row for each column, a column for each window. */
  function turned(columns, rows, callOfRow) {
    const head = columns.find((c) => c.saved === 'Stimulus'), role = columns.find((c) => c.saved === 'Role');
    const out = [{ name: '', given: true, index: 0 }].concat(rows.map((row, k) => ({ name: head ? String(row[head.index] || '–') + (role && row[role.index] === 'Control' ? ' (control)' : '') : 'Value', given: true, index: k + 1 })));
    return { columns: out, rows: columns.filter((c) => c !== head).map((c) => [{ text: headerText(c), title: headerTitle(c) }].concat(rows.map((row) => shownCell(c, row, callOfRow)))) };
  }
  /** Draws a table. A saved table (kind: which one) keeps its saved names and values; on screen its columns are named in words
   *  with their units, its numbers are rounded by their unit, and pointing at a header or a number shows what was saved. A table
   *  worked out here brings its columns (data.columns, given: its cells are { text, title, key } or plain text). The first column
   *  (data.frozen: that many) stays in place when the table scrolls sideways. A click on a header sorts by it, blanks last. */
  function table(id, data, kind) {
    const element = $(id); element.innerHTML = ''; if (!data || !data.header) return;
    const box = { 'table-calls': 'all-call-columns', 'table-metrics': 'all-metric-columns', 'table-counts': 'all-count-columns' }[id], everything = !!(box && $(box).checked);
    const callOfRow = kind ? callReader(kind, data) : () => undefined, flip = kind === 'counts' && data.rows.length >= 1 && data.rows.length <= 4;
    let columns = (data.columns || tableColumns(kind, data, everything)).concat(data.more || []), rows = data.rows;   // more: columns put after the saved ones
    if (flip) ({ columns, rows } = turned(columns, rows, callOfRow));
    const given = (v) => v && typeof v === 'object' ? v : { text: v == null || v === '' ? '–' : String(v), key: v == null ? '' : v };
    const shown = rows.map((row) => ({ row, cells: columns.map((c) => c.given ? given(row[c.index]) : shownCell(c, row, callOfRow)) }));
    const frozen = data.frozen || 1, cellColumn = kind && !flip ? data.header.indexOf('Cell') : -1;
    const head = element.createTHead().insertRow(), body = element.createTBody(); let sorted = -1, direction = 1;
    /** Frozen columns after the first stand where the ones before them end. */
    const pin = () => { if (frozen < 2) return; let left = 0;
      for (let k = 0; k < frozen && head.cells[k]; k++) { for (const tr of element.rows) if (tr.cells[k]) tr.cells[k].style.left = left + 'px'; left += head.cells[k].getBoundingClientRect().width; } };
    const fill = (list) => { body.innerHTML = '';
      if (!list.length) { const td = body.insertRow().insertCell(); td.colSpan = Math.max(1, columns.length); td.className = 'empty-table'; td.textContent = 'No rows were saved for this table.'; return; }
      for (const { row, cells } of list) { const tr = body.insertRow(); if (cellColumn >= 0) { tr.dataset.cell = row[cellColumn]; tr.classList.toggle('on', row[cellColumn] === S.cell); }
        cells.forEach((cell, k) => { const td = tr.insertCell(); td.textContent = cell.text; if (cell.title) td.title = cell.title; if (cell.num) td.className = 'num'; if (k < frozen) td.classList.add('cell-column'); }); }
      pin();
    };
    const blank = (key) => key === '' || key == null || (typeof key === 'number' && !isFinite(key));
    const compare = (a, b) => blank(a) || blank(b) ? (blank(a) ? 1 : 0) - (blank(b) ? 1 : 0) : direction * (typeof a === 'number' && typeof b === 'number' ? a - b : String(a).localeCompare(String(b), undefined, { numeric: true }));
    columns.forEach((c, k) => { const th = document.createElement('th'), title = headerTitle(c); th.textContent = headerText(c); if (title) th.title = title; if (k < frozen) th.classList.add('cell-column'); head.appendChild(th);
      if (flip) { th.classList.add('fixed'); return; }
      th.addEventListener('click', () => { direction = sorted === k ? -direction : 1; sorted = k; for (const other of head.children) other.classList.remove('up', 'down'); th.classList.add(direction > 0 ? 'up' : 'down');
        fill(shown.slice().sort((a, b) => compare(a.cells[k].key, b.cells[k].key))); }); });
    body.addEventListener('click', (event) => { const tr = event.target.closest('tr'); if (tr && tr.dataset.cell) chooseCell(tr.dataset.cell); });
    fill(shown);
  }
  /** The calls table as saved. While review marks apply, two columns on its right: the reviewed call and the review reason with
   *  its note, on the primary stimulus's row of each marked cell; empty on every other row. The saved Status and Reason stay. */
  function showCalls() {
    const data = withCodes('calls', S.results.tables.calls).data, review = reviewOfRows(data);
    if (!review) { table('table-calls', data, 'calls'); return; }
    const n = data.header.length, blank = { text: '', key: '' };
    table('table-calls', Object.assign({}, data, {
      rows: data.rows.map((row, i) => { const m = review[i], why = m ? reviewReason(m) : ''; return fullRow(row, n).concat(m ? [{ text: m.reviewed, key: m.reviewed }, { text: why, key: why }] : [blank, blank]); }),
      more: [{ name: 'Reviewed call', given: true, index: n, title: 'The call after review. The automatic call is in Call.' }, { name: 'Review reason', given: true, index: n + 1, title: 'Why the call was changed, and the note' }] }), 'calls');
  }
  /** A mark's reason as the Review reason column says it: "Passing object, automatic" for a mark made by the option, then the note. */
  const automaticMark = (mark) => mark.markedBy === 'automatic';
  const reviewReason = (mark) => mark.reason + (automaticMark(mark) ? ', automatic' : '') + (mark.note ? ': ' + mark.note : '');
  /** The calls table as a CSV: as saved, with the codes of the call checks in QC_flags (withCodes), and while review marks apply
   *  with Reviewed_call, Review_reason and Review_note on the right, then Marked_by when a mark was made automatically. */
  function callsRows() {
    const data = withCodes('calls', S.results.tables.calls).data, rows = [data.header].concat(data.rows), review = reviewOfRows(data); if (!review) return rows;
    const n = data.header.length, by = review.some((m) => m && automaticMark(m));
    return [fullRow(rows[0], n).concat(['Reviewed_call', 'Review_reason', 'Review_note'], by ? ['Marked_by'] : [])]
      .concat(data.rows.map((row, i) => { const m = review[i]; return fullRow(row, n).concat(m ? [m.reviewed, m.reason, m.note || ''] : ['', '', ''], by ? [m ? m.markedBy || 'manual' : ''] : []); }));
  }
  /** Whether the calls table, as exported, has codes of the call checks in its QC_flags. */
  const callsCoded = () => withCodes('calls', S.results.tables.calls).added;
  const showCounts = () => table('table-counts', S.results.tables.counts, 'counts');

  // ---------------------------------------------------------------- the processed data, cell by cell
  const Data = { signal: 'dff', first: 0, page: 50 };
  const dataTraces = () => Data.signal !== 'dff' && hasTable(Data.signal) ? S.results[Data.signal] : S.results.traces;
  /** The whole table as CSV text: time, then one column per cell under its full ID. */
  function dataCsv() {
    const r = S.results, traces = dataTraces(), names = S.info.names.filter((n) => traces[n]);
    const rows = [['Frame', 'Time_s'].concat(names).join(',')];
    for (let i = 0; i < r.time.length; i++) rows.push([i + 1, r.time[i]].concat(names.map((n) => finite(traces[n][i]) ? traces[n][i] : '')).join(','));
    return rows.join('\n') + '\n';
  }
  /** Fifty frames at a time on screen (the saved file has them all); cells are headed "Cell 1", with the ID in the tooltip.
   *  Frame and time stay in place when the table scrolls sideways. */
  function showData() {
    const r = S.results; if (!r) return;
    for (const option of $('data-signal').options) option.disabled = option.value !== 'dff' && !hasTable(option.value);
    if ($('data-signal').selectedOptions[0].disabled) { Data.signal = 'dff'; $('data-signal').value = 'dff'; }
    const traces = dataTraces(), names = S.info.names.filter((n) => traces[n]), total = r.time.length;
    Data.first = Math.max(0, Math.min(Data.first, Math.floor((total - 1) / Data.page) * Data.page)); const last = Math.min(total, Data.first + Data.page);
    const rows = []; for (let i = Data.first; i < last; i++) rows.push([i + 1, r.time[i]].concat(names.map((n) => finite(traces[n][i]) ? traces[n][i] : '')));
    // the saved names are those of the file (dataCsv): Frame, Time_s and each cell's full ID
    const unit = Data.signal === 'raw' || Data.signal === 'corrected' ? GRAY : DFF;
    const columns = [{ saved: 'Frame', name: 'Frame', unit: '', decimals: 0 }, { saved: 'Time_s', name: 'Time', unit: SEC, decimals: DECIMALS[SEC] }]
      .concat(names.map((n) => ({ saved: n, name: shortName(n), unit: '', measured: unit, decimals: DECIMALS[unit] }))).map((c, index) => Object.assign(c, { index, numeric: true, measured: c.measured || '' }));
    table('table-data', { header: columns.map((c) => c.saved), columns, rows, frozen: 2 });
    $('data-range').textContent = 'Frames ' + (Data.first + 1) + '–' + last + ' of ' + total + ' · ' + names.length + ' cells';
    $('data-back').disabled = Data.first === 0; $('data-forward').disabled = last >= total;
  }
  const dataKind = () => 'processed_data_' + (Data.signal === 'raw' && hasTable('raw') ? 'raw_gray_values' : Data.signal === 'rawdff' && hasTable('rawdff') ? 'dFF0_of_raw_values' : Data.signal === 'corrected' && hasTable('corrected') ? 'gray_levels_above_background' : 'dFF0');
  function showTab(name) {
    for (const tab of document.querySelectorAll('.tab')) { tab.classList.toggle('on', tab.dataset.tab === name); tab.setAttribute('aria-selected', tab.dataset.tab === name ? 'true' : 'false'); }
    for (const panel of document.querySelectorAll('.panel')) panel.hidden = panel.id !== 'panel-' + name;
    if (name === 'heatmap') drawHeatmap(); if (name === 'review') { Plot.dirty = true; drawPlot(); drawPicture(); } if (name === 'data') showData();
    const on = document.querySelector('tbody tr.on'); if (on && name !== 'review' && name !== 'heatmap') on.scrollIntoView({ block: 'nearest' });
  }
  function download(blob, name) { const link = document.createElement('a'); link.href = URL.createObjectURL(blob); link.download = name; link.click(); setTimeout(() => URL.revokeObjectURL(link.href), 2000); }
  /** A file name that says what it is: the recording, what was saved, the cell when it matters, and when. */
  function exportName(kind, extension) {
    const clean = (text) => String(text || '').replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 80);
    const now = new Date(), two = (n) => String(n).padStart(2, '0');
    const stamp = now.getFullYear() + two(now.getMonth() + 1) + two(now.getDate()) + '-' + two(now.getHours()) + two(now.getMinutes()) + two(now.getSeconds());
    return [clean(S.info.recording), kind, /^(traces|around_the_cell)/.test(kind) && S.cell ? (shortName(S.cell) === S.cell ? clean(S.cell) : clean(shortName(S.cell)) + '_' + clean(S.cell).slice(0, 13)) : '', stamp].filter(Boolean).join('_') + '.' + extension;
  }
  /** What the page showed, kept beside a saved file. @param followsReview whether what is saved follows the review marks: a figure
   *  or picture of the primary stimulus, the calls table and the metric summary do; the counts, the processed data, the per-cell
   *  metric tables and anything of a control window do not. */
  function exportMetadata(followsReview) {
    return { schema: 'simpleca.view-export', schema_version: 1, recording: S.info.recording, plugin_version: S.info.pluginVersion,
      appliedFingerprint: S.results.appliedFingerprint || S.state.appliedFingerprint, appliedProtocol: S.results.appliedProtocol || S.results.protocol,
      manifest: S.results.manifest, settingsPending: hasPending(), displayed_window: S.window, cell: S.cell, group: GROUPS[S.group],
      frame: S.frame, time_s: frameTimes()[S.frame - 1], view: S.view, results_folder: S.results.folder || (S.state && S.state.folder) || null,
      time_range: Plot.geom ? [Plot.geom.x0, Plot.geom.x1] : [Plot.x0, Plot.x1], time_zoomed: Plot.x0 !== null, palette: S.palette, range: S.range, order: S.order,
      plot: S.plotMode, signal: S.plotMode === 'around' ? 'dF (gray levels)' : S.results && isDff() ? 'dF/F0' : S.signal, y_axis: Plot.geom ? [Plot.geom.lo, Plot.geom.hi] : null, y_axis_typed: Plot.y[yKey()] || null, plotted: Plot.geom ? Plot.geom.plotted : [], all_traces_drawn: S.every, thresholds_drawn: S.thresholds, cells_left_out: Array.from(leftOut()),
      // while review marks apply, the calls of the primary stimulus are those after review, as the Fiji window's record says too
      ...(marked() && followsReview ? { calls: 'after review', review_file: S.results.review.file, review_sha256: S.results.review.sha256 } : {}) };
  }
  /** What a saved figure itself shows, in place of what the page shows: its time and Y range and the means it drew, as the
   *  figure's own drawing (Fig.drawn, made by picture()) has them; the limits typed in the editor are under picture.figure. */
  function figureShows(kind) {
    const drew = Fig.drawn; if (!drew || drew.kind !== kind) return {};
    const typedTime = Fig[kind].xMin != null || Fig[kind].xMax != null;
    return kind === 'plot' ? Object.assign({ time_range: [drew.x0, drew.x1], time_zoomed: Plot.x0 !== null || typedTime, y_axis: [drew.lo, drew.hi], plotted: drew.plotted }, drew.around ? { around_the_cell: aroundRecord(drew.around) } : {})
      : { time_range: [drew.x0, drew.x1], time_zoomed: typedTime };
  }
  /** What a saved "Around the cell" shows, for the file kept beside it: the cell, its flag, each zone drawn (by number, place and
   *  part, with the frame where a sector was seen), the other ring sectors, and the flagged frames, as frames from 1 and as the
   *  times of the first and the last of them. */
  const aroundRecord = (around) => ({ cell: S.cell, rule: around.rule, flag: around.text, units: 'gray levels',
    lines: around.lines.map((line) => ({ zone: line.zone, place: line.place, part: line.part, frame: line.frame })), other_ring_sectors: around.others.map((other) => ({ zone: other.zone, place: other.place })),
    flagged_frames: around.band ? around.frames : null, flagged_frames_s: around.band ? [around.time[around.frames[0] - 1], around.time[around.frames[1] - 1]] : null });
  /** How saved pictures are written. Format, resolution and TIFF compression are one setting for every saved picture, the
   *  figures and the pictures of the cells alike; size is the width of a picture of the cells (a figure has its own). Each has
   *  a list on the page (saveBoxes), and showSaveOptions makes those lists show what is set here. */
  const Save = { format: 'png', size: '8x6', dpi: 300, compression: 'lzw' };
  const inches = () => Save.size === 'screen' ? null : Save.size.split('x').map(Number);
  const saveBoxes = { format: 'save-format', size: 'save-size', dpi: 'save-dpi', compression: 'save-compression' };
  function showSaveOptions() {
    for (const key of Object.keys(saveBoxes)) $(saveBoxes[key]).value = String(Save[key]);
    $('save-compression-label').hidden = Save.format !== 'tiff';
    const size = inches(), wide = S.info ? Math.round(S.info.width * pictureScale()) : 0;
    $('save-pixels').textContent = wide ? 'Pictures and frames ' + wide + ' px wide' + (size && wide < Math.round(size[0] * Save.dpi) ? ' (4× the recording, the most it allows).' : '.') : '';
  }
  // ---------------------------------------------------------------- figures: how a saved plot or heatmap looks
  /** What belongs to one recording or one signal and is not carried to the next: the title, the Y axis title and the axis numbers. */
  const FIGURE_OWN = ['title', 'yTitle', 'xMin', 'xMax', 'xStep', 'yMin', 'yMax', 'yStep', 'barMin', 'barMax'];
  const FIGURE_MOST = 80e6;   // pixels in one saved figure: beyond this a browser may not draw it
  /** The room a figure's plot area (the heatmap's picture) should have. said: the share of the figure; under it, the editor
   *  says beside the text size that the text is too large for the figure. least: its width and its height in heights of the
   *  figure's text; under it the text and the legend have left no plot to speak of. A number that would do that is not taken
   *  (FIGURE_ROOM), and a figure that got there another way is not saved. */
  const FIT = { said: 0.4, least: 4 };
  /** The numbers that decide whether the plot has room, each with what the editor says when it would have none. */
  const FIGURE_ROOM = { text: 'Too large for this figure size', titleText: 'Too large for this figure size', width: 'Too small for this text size', height: 'Too small for this text size' };
  /** The two figures as edited (each starts as its style), which one the editor shows, and what the last drawing placed where.
   *  earlier: what start() took from the preferences, for each figure and for the file settings (Save), so that the editor
   *  can say while any of it is in use. rows: the editor's number boxes, each with the line it speaks on (figureSay). */
  const Fig = { plot: Object.assign({}, STYLE.plot), heatmap: Object.assign({}, STYLE.heatmap), open: null, queued: false, drawn: null, earlier: { plot: {}, heatmap: {}, save: {} }, rows: {}, yKept: {} };
  /** The boxes of the figure editor, in order: [section, key, label, type, more]; key is an entry of the style table. Types:
   *  number (unit, with the least and the most it takes, or "above": the number it must be greater than; an empty box is
   *  "auto" where the default is null; scale: what the box shows for 1 of the entry, 100 for an opacity typed in %), text
   *  (the figure's own text until one is typed), choice, tick, colour, and page (a list
   *  or a tick box that lives on the page, such as the heatmap's colours or its cell names, shown here too so one place has
   *  everything: its key is the id of the page's control). "when" leaves out a box for something the plot on screen does not show. */
  const FIGURE_SIZE = [['Size', 'width', 'Width', 'number', { unit: 'in', least: 1, most: 20 }], ['Size', 'height', 'Height', 'number', { unit: 'in', least: 1, most: 20 }],
    ['File', 'save-dpi', 'Resolution', 'page'], ['File', 'save-format', 'Format', 'page'], ['File', 'save-compression', 'TIFF compression', 'page']];
  /** What the editor says beside the name of a section: the file settings are not this figure's alone (Save). */
  const FIGURE_SECTIONS = { File: 'for every saved picture' };
  const FIGURE_TEXT = [['Text', 'font', 'Font', 'choice', { options: FONTS.map((font) => font.slice(0, 2)) }],
    ['Text', 'text', 'Text size', 'number', { unit: 'pt', least: 4, most: 48 }], ['Text', 'titleText', 'Title size', 'number', { unit: 'pt', least: 4, most: 72 }],
    ['Text', 'ink', 'Text and axes', 'colour'], ['Text', 'background', 'Background', 'choice', { options: [['white', 'White'], ['transparent', 'Transparent']] }]];
  const cellDrawn = () => (S.plotMode === 'cell' && !!S.cell && !!shownTraces()[S.cell]) || aroundDrawn(), byCall = () => S.plotMode === 'calls', cellCalled = () => S.plotMode === 'cell' && isDff();
  const meansDrawn = () => S.plotMode !== 'around', everyDrawn = () => S.every && meansDrawn();
  const FIGURE_FIELDS = {
    plot: FIGURE_SIZE.concat([['Text', 'title', 'Title', 'text'], ['Text', 'xTitle', 'X axis title', 'text'], ['Text', 'yTitle', 'Y axis title', 'text']], FIGURE_TEXT, [
      ['Axes', 'xMin', 'X min', 'number', {}], ['Axes', 'xMax', 'X max', 'number', {}], ['Axes', 'xStep', 'X ticks every', 'number', { above: 0 }],
      ['Axes', 'yMin', 'Y min', 'number', {}], ['Axes', 'yMax', 'Y max', 'number', {}], ['Axes', 'yStep', 'Y ticks every', 'number', { above: 0 }],
      ['Axes', 'frame', 'Frame', 'choice', { options: [['open', 'Left and bottom'], ['box', 'Box']] }], ['Axes', 'grid', 'Grid lines', 'tick'], ['Axes', 'axisWidth', 'Axis thickness', 'number', { unit: 'pt', least: 0.25, most: 6 }],
      ['Lines', 'meanWidth', 'Mean line', 'number', { unit: 'pt', least: 0.25, most: 10, when: meansDrawn }], ['Lines', 'colourMean', 'Mean colour', 'colour', { when: () => S.plotMode !== 'calls' && meansDrawn() }],
      ['Lines', 'sem', 'SEM', 'choice', { options: [['band', 'Shaded band'], ['none', 'Not shown']], when: meansDrawn }],
      ['Lines', 'cellWidth', 'Cell line', 'number', { unit: 'pt', least: 0.25, most: 10, when: cellDrawn }], ['Lines', 'colourCell', 'Cell colour', 'colour', { when: cellDrawn }],
      ['Lines', 'zoneWidth', 'Zone lines', 'number', { unit: 'pt', least: 0.25, most: 10, when: aroundDrawn }],
      ['Lines', 'traceWidth', 'Trace line', 'number', { unit: 'pt', least: 0.1, most: 6, when: everyDrawn }], ['Lines', 'colourTrace', 'Trace colour', 'colour', { when: everyDrawn }],
      ['Lines', 'traceOpacity', 'Trace opacity', 'number', { unit: '%', least: 1, most: 100, scale: 100, when: everyDrawn }],
      ['Lines', 'responderWidth', 'Responder line', 'number', { unit: 'pt', least: 0.1, most: 6, when: everyDrawn }], ['Lines', 'responderOpacity', 'Responder opacity', 'number', { unit: '%', least: 1, most: 100, scale: 100, when: everyDrawn }],
      ['Lines', 'colourResponded', 'Responded', 'colour', { when: () => everyDrawn() || byCall() }], ['Lines', 'colourNot', 'Did not respond', 'colour', { when: byCall }],
      ['Marks', 'legend', 'Legend', 'choice', { options: [['top', 'Above the plot'], ['right', 'Right of the plot'], ['none', 'Not shown']] }], ['Marks', 'stimulus', 'Stimulus line', 'tick'], ['Marks', 'stimulusText', 'Stimulus label', 'text'],
      ['Marks', 'baseline', 'Baseline', 'tick', { when: meansDrawn }], ['Marks', 'window', 'Response window', 'tick', { when: meansDrawn }], ['Marks', 'flagged', 'Flagged frames', 'tick', { when: aroundDrawn }],
      ['Marks', 'run', 'Response', 'tick', { when: cellCalled }], ['Marks', 'levels', 'Required rise', 'tick', { when: cellCalled }], ['Marks', 'thresholds', 'Thresholds', 'page', { when: cellCalled }]]),
    heatmap: FIGURE_SIZE.concat([['Text', 'title', 'Title', 'text'], ['Text', 'xTitle', 'X axis title', 'text'], ['Text', 'yTitle', 'Y axis title', 'text'], ['Text', 'barTitle', 'Colour bar title', 'text']], FIGURE_TEXT, [
      ['Axes', 'xMin', 'X min', 'number', {}], ['Axes', 'xMax', 'X max', 'number', {}], ['Axes', 'xStep', 'X ticks every', 'number', { above: 0 }],
      ['Colours', 'palette', 'Colours', 'page'], ['Colours', 'range', 'Range', 'page'], ['Colours', 'barMin', 'Scale from', 'number', {}], ['Colours', 'barMax', 'Scale to', 'number', {}], ['Colours', 'bar', 'Colour bar', 'tick'],
      ['Rows', 'order', 'Sort', 'page'], ['Rows', 'row-names', 'Cell names', 'page'], ['Rows', 'hide-excluded', 'Hide excluded cells', 'page'],
      ['Marks', 'marks', 'Windows and stimulus', 'tick'], ['Marks', 'strip', 'Response calls', 'tick']])
  };
  /** The boxes of an axis: its two limits and the spacing of its ticks (the colour scale has no spacing). */
  const FIGURE_AXES = { x: ['xMin', 'xMax', 'xStep'], y: ['yMin', 'yMax', 'yStep'], colour: ['barMin', 'barMax'] };
  /** Whether a value may be kept for a box by the box's own limits: the one rule for what is typed and for what comes back from the preferences. */
  function figureValid(kind, key, value) {
    const field = FIGURE_FIELDS[kind].find((f) => f[1] === key), type = field && field[3], more = (field && field[4]) || {};
    if (type === 'number') { if (value === null) return STYLE[kind][key] === null; if (!finite(value)) return false; const v = boxUnits(more, value);
      return (more.least == null || v >= more.least) && (more.most == null || v <= more.most) && (more.above == null || v > more.above); }
    if (type === 'text') return value === null || (typeof value === 'string' && value.length <= 80);
    if (type === 'tick') return typeof value === 'boolean'; if (type === 'colour') return typeof value === 'string' && /^#[0-9a-f]{6}$/i.test(value);
    return type === 'choice' && more.options.some((option) => option[0] === value);
  }
  /** A number of the style table as its box shows it, and a number typed in a box as the style table keeps it (more.scale). */
  const boxUnits = (more, value) => more.scale ? +(value * more.scale).toPrecision(12) : value;
  const styleUnits = (more, value) => more.scale && finite(value) ? value / more.scale : value;
  /** The ends of a figure's axes before a limit is typed for them in the editor, as the painters find them: the time shown and,
   *  for the plot, the Y range. The colour scale has none: a limit typed for one of its ends moves the other (heatmapRange). */
  function figureAuto(kind) {
    if (kind === 'heatmap') return { x: wholeTime() };
    const look = lookOf(Fig.plot), times = plotTimes(look); return { x: times.auto, y: plotLimits(plotParts(look, times.x0, times.x1), look).auto };
  }
  /** Why a figure cannot use a number, in the words the editor puts under its box; '' when it can. value: the number, null for
   *  an empty box, NaN for text that is no number. First the box's own limits (figureValid). Then the room for the plot: a
   *  size of the text or of the figure that would leave it under the least (FIT), and with less than it has now. Then the rules
   *  of an axis, which span its boxes: its two limits must leave a range, each with the other one or, where that is not
   *  typed, with the automatic end; and its ticks may not be closer than the painters draw them. The painters go by the same
   *  two rules (typedRange, steppedTicks): a number that passed here and no longer does, because the plot under the figure
   *  has changed, is not drawn either, and the editor marks it. */
  function figureRefusal(kind, key, value) {
    const more = FIGURE_FIELDS[kind].find((f) => f[1] === key)[4] || {}, fig = Fig[kind], said = (v) => String(+v.toPrecision(6)).replace('-', '−');
    if (!figureValid(kind, key, value)) return more.least != null ? more.least + ' to ' + more.most + ' ' + more.unit : more.above != null ? 'Above ' + more.above : 'Not a number';
    if (FIGURE_ROOM[key] && value !== fig[key]) { const room = figureFit(kind, { [key]: value }).room; return room < FIT.least && room < figureFit(kind).room ? FIGURE_ROOM[key] : ''; }
    const axis = value === null ? null : Object.keys(FIGURE_AXES).find((name) => FIGURE_AXES[name].includes(key)); if (!axis) return '';
    const [min, max, step] = FIGURE_AXES[axis], auto = figureAuto(kind)[axis] || [-Infinity, Infinity], now = (k) => (k === key ? value : fig[k]), range = typedRange(auto, now(min), now(max));
    if (key === min) return range ? '' : finite(fig[max]) ? 'Minimum must be below maximum' : 'Must be below ' + said(auto[1]);
    if (key === max) return range ? '' : finite(fig[min]) ? 'Maximum must be above minimum' : 'Must be above ' + said(auto[0]);
    const [lo, hi] = range || auto, least = (hi - lo) / TICKS_MOST, digit = Math.pow(10, Math.floor(Math.log10(least)) - 1);   // (the least spacing is said to two figures, rounded up)
    return value < least ? 'At least ' + said(Math.ceil(least / digit - 1e-9) * digit) : '';
  }
  /** A figure's settings as kept in the preferences: only what differs from the defaults, without what belongs to one recording. */
  const figureKept = (kind) => { const kept = {}; for (const key of Object.keys(STYLE[kind])) if (!FIGURE_OWN.includes(key) && Fig[kind][key] !== STYLE[kind][key]) kept[key] = Fig[kind][key]; return kept; };
  /** Takes kept settings back: known boxes and valid values only, so an old or damaged preference cannot break a figure. */
  function figureRestore(kind, text) {
    let kept; try { kept = JSON.parse(text); } catch (broken) { return; } if (!kept || typeof kept !== 'object') return;
    for (const key of Object.keys(STYLE[kind])) if (!FIGURE_OWN.includes(key) && kept[key] !== undefined && kept[key] !== STYLE[kind][key] && figureValid(kind, key, kept[key])) Fig.earlier[kind][key] = Fig[kind][key] = kept[key];
  }
  /** Whether a figure still uses a setting kept from an earlier session: one of its own, or one of the file settings the editor shows. */
  const earlierInUse = (kind) => Object.keys(Fig.earlier[kind]).some((key) => Fig[kind][key] === Fig.earlier[kind][key]) ||
    ['format', 'dpi', 'compression'].some((key) => Save[key] === Fig.earlier.save[key] && (key !== 'compression' || Save.format === 'tiff'));
  /** A figure measured out from its settings: its size in pixels at 96 to the inch, its look, its title and the height of the
   *  row the title takes (none without a title). */
  function figureFrame(fig) { const look = lookOf(fig), title = fig.title || ''; return { w: Math.round(fig.width * 96), h: Math.round(fig.height * 96), look, title, top: title ? 6 * look.px / 12 + look.titlePx * 1.5 : 0 }; }
  /** What FIT is held against, from a figure's frame and where its parts went (g): how much of the figure the plot area takes
   *  (share), how much of what the title leaves (body), and its shorter side in heights of the text (room). */
  const fitOf = (frame, g) => ({ share: g.W * g.H / (frame.w * frame.h), body: g.W * g.H / (frame.w * Math.max(1, frame.h - frame.top)), room: Math.min(g.W, g.H) / frame.look.px });
  /** The fit a figure would have with some settings changed (none: as it is): laid out as the painters lay it out, nothing drawn. */
  function figureFit(kind, changes) {
    const frame = figureFrame(Object.assign({}, Fig[kind], changes)), ctx = document.createElement('canvas').getContext('2d'), { w, h, look, top } = frame;
    return fitOf(frame, kind === 'heatmap' ? heatmapLayout(ctx, w, h - top, look, heatmapData(look)) : plotPlan(ctx, w, h - top, look).g);
  }
  /** The plot or the heatmap as a figure of the size set in the editor: a title when one was typed, then the drawing by the
   *  page's own painter. The preview and the saved file are this one drawing at different numbers of pixels per inch. */
  function figureCanvas(kind, perInch) {
    const fig = Fig[kind], frame = figureFrame(fig), { w, h, look, title, top: y } = frame;
    const canvas = document.createElement('canvas'); canvas.width = Math.max(1, Math.round(fig.width * perInch)); canvas.height = Math.max(1, Math.round(fig.height * perInch));
    if (canvas.width * canvas.height > FIGURE_MOST) throw new Error('This figure would be ' + canvas.width + ' × ' + canvas.height + ' pixels. Choose a smaller size or resolution.');
    const ctx = canvas.getContext('2d'), u = (n) => n * look.px / 12, family = look.family || getComputedStyle(document.body).fontFamily; ctx.setTransform(canvas.width / w, 0, 0, canvas.height / h, 0, 0);
    if (look.background) { ctx.fillStyle = look.background; ctx.fillRect(0, 0, w, h); }
    ctx.save(); ctx.translate(0, y); look.background = null;   // under the title's row, written last; the whole figure has its background already
    const drew = kind === 'heatmap' ? paintHeatmap(ctx, w, h - y, look, heatmapData(look)) : paintTraces(ctx, w, h - y, look); ctx.restore();
    Fig.drawn = Object.assign({ kind, title, top: y }, fitOf(frame, drew), drew);
    // the title starts where the legend or the key under it does: on the plot's Y axis, or at the heatmap's left edge. One too long for the figure is shortened.
    if (title) { const x = drew.key ? drew.key.x : drew.L; ctx.fillStyle = look.ink; ctx.font = '600 ' + look.titlePx + 'px ' + family; ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; writeText(ctx, fitted(ctx, title, w - x - u(6)), x, u(6) + look.titlePx * 0.75, look); }
    return canvas;
  }
  /** The pictures the Export menus and the editor save. "figure" says which figure's settings made them, for the file kept
   *  beside. A figure whose plot has no room (FIT.least) is not saved; the editor shows it and says so. */
  const figurePicture = (kind) => Object.assign(() => { const canvas = figureCanvas(kind, Save.dpi); if (Fig.drawn.room < FIT.least) throw new Error('The text is too large for this figure size. Change it in Edit figure.'); return canvas.toDataURL('image/png'); }, { figure: kind });
  const plotPicture = figurePicture('plot'), heatmapPicture = figurePicture('heatmap');
  /** What a text box of the editor shows before anything is typed: the figure's own text. A figure has no title of its own. */
  const figureOwnText = (kind, key) => key === 'xTitle' ? 'Time (s)' : key === 'yTitle' ? (kind === 'plot' ? plotYTitle() : 'Cells') : key === 'barTitle' ? 'ΔF/F₀' : '';
  /** Says what a number box of the editor has to say, on the line under it: why the number in it is not used ("why"; the box is
   *  then marked), or else a remark on a number that is used ("remark"). */
  function figureSay(key, words) {
    const row = Fig.rows[key]; if (!row) return; Object.assign(row, words); const text = row.why || row.remark;
    row.note.textContent = text; row.note.hidden = !text; if (row.why) row.box.setAttribute('aria-invalid', 'true'); else row.box.removeAttribute('aria-invalid');
  }
  /** Draws the preview: the saved picture itself, fitted into the stage. Asked for on every change, drawn once per screen refresh.
   *  Then what the drawing showed about its fit: text that crowds the plot, beside the size of that text; and no saving of a
   *  figure with too many pixels or with no room for its plot. */
  function previewFigure() {
    if (!Fig.open || Fig.queued) return; Fig.queued = true;
    requestAnimationFrame(() => {
      Fig.queued = false; const kind = Fig.open; if (!kind || !S.results) return;
      const fig = Fig[kind], stage = $('figure-stage'), shown = $('figure-preview'), wide = Math.round(fig.width * Save.dpi), high = Math.round(fig.height * Save.dpi);
      const scale = Math.max(0.05, Math.min((stage.clientWidth - 32) / (fig.width * 96), (stage.clientHeight - 32) / (fig.height * 96), 2)), made = figureCanvas(kind, 96 * scale * (window.devicePixelRatio || 1));
      shown.width = made.width; shown.height = made.height; shown.style.width = Math.round(fig.width * 96 * scale) + 'px'; shown.style.height = Math.round(fig.height * 96 * scale) + 'px';
      const ctx = shown.getContext('2d'); ctx.clearRect(0, 0, shown.width, shown.height); ctx.drawImage(made, 0, 0);
      const d = Fig.drawn, note = $('figure-pixels'), most = wide * high > FIGURE_MOST, none = d.room < FIT.least; note.classList.toggle('bad', most || none);
      note.textContent = wide + ' × ' + high + ' pixels' + (most ? ': too large to save. Choose a smaller size or resolution.' : none ? ': no room for the ' + kind + '. Choose smaller text or a larger size.' : ' · shown at ' + Math.round(100 * scale) + ' %');
      $('figure-save').disabled = $('figure-download').disabled = most || none;
      // the title's size is the one to change when the plot would have its room without the title; else the text size
      const crowding = d.share >= FIT.said ? null : d.body < FIT.said ? 'text' : 'titleText';
      for (const key of ['text', 'titleText']) figureSay(key, { remark: key === crowding ? 'Too large for this figure size' : '' });
    });
  }
  /** Rows that depend on another: TIFF compression belongs to TIFF. */
  function figureRows() { const box = $('figure-save-compression'); if (box) box.closest('.row').hidden = Save.format !== 'tiff'; }
  /** Builds the editor's boxes for the plot or the heatmap from FIGURE_FIELDS and opens it. */
  function openFigure(kind) {
    if (!S.results) return; Fig.open = kind; Fig.rows = {}; const fig = Fig[kind], holder = $('figure-controls'), numbers = []; holder.replaceChildren();
    $('figure-name').textContent = kind === 'plot' ? 'Plot figure' : 'Heatmap figure'; $('figure-saved').textContent = '';
    const earlier = document.createElement('p'); earlier.className = 'fine'; earlier.id = 'figure-earlier'; earlier.textContent = 'Settings kept from an earlier session are in use.'; holder.appendChild(earlier);
    // After anything is set, the rules that span boxes are looked at again. A box still holding a number the figure could not
    // use hands it over now if it can, and so on until none does. Then every number box says where it stands: also one whose
    // number the figure holds and can no longer use, because its axis has changed under it.
    const review = () => { while (numbers.some((box) => box.take())) { } for (const box of numbers) box.mark(); earlier.hidden = !earlierInUse(kind); };
    const changed = () => { review(); store('step4.figure.' + kind, JSON.stringify(figureKept(kind))); previewFigure(); };
    let section = null;
    for (const [group, key, label, type, more] of FIGURE_FIELDS[kind]) {
      const options = more || {}; if (options.when && !options.when()) continue;
      if (group !== section) { section = group; const head = document.createElement('h3'); head.textContent = group; holder.appendChild(head);
        if (FIGURE_SECTIONS[group]) { const about = document.createElement('span'); about.textContent = FIGURE_SECTIONS[group]; head.appendChild(about); } }
      const row = document.createElement('label'), name = document.createElement('span'); row.className = 'row'; name.textContent = label; row.appendChild(name); let box;
      if (type === 'page') {   // the same list or tick box as on the page: choosing here chooses there
        const source = $(key), ticked = source.type === 'checkbox', value = ticked ? 'checked' : 'value'; box = document.createElement(ticked ? 'input' : 'select'); if (ticked) box.type = 'checkbox';
        for (const o of source.options || []) { const copy = document.createElement('option'); copy.value = o.value; copy.textContent = o.textContent; box.appendChild(copy); }
        box[value] = source[value]; box.addEventListener('change', () => { source[value] = box[value]; source.dispatchEvent(new Event('change', { bubbles: true })); figureRows(); changed(); });
      } else if (type === 'choice') {
        box = document.createElement('select'); for (const [value, text] of options.options) { const o = document.createElement('option'); o.value = value; o.textContent = text; box.appendChild(o); }
        box.value = fig[key]; box.addEventListener('change', () => { fig[key] = box.value; changed(); });
      } else if (type === 'tick') { box = document.createElement('input'); box.type = 'checkbox'; box.checked = !!fig[key]; box.addEventListener('change', () => { fig[key] = box.checked; changed(); });
      } else if (type === 'colour') {   // shown while it is picked, kept when it is chosen
        box = document.createElement('input'); box.type = 'color'; box.value = fig[key]; box.addEventListener('input', () => { fig[key] = box.value; previewFigure(); }); box.addEventListener('change', changed);
      } else if (type === 'text') {   // shown while it is typed, kept when the box is left
        box = document.createElement('input'); box.type = 'text'; box.maxLength = 80; box.spellcheck = false; box.value = fig[key] == null ? figureOwnText(kind, key) : fig[key]; if (!figureOwnText(kind, key)) box.placeholder = key === 'title' ? 'No title' : 'No label';
        box.addEventListener('input', () => { fig[key] = box.value; previewFigure(); }); box.addEventListener('change', changed);
      } else {   // a number; an empty box is "auto" where the figure can choose. While it is typed the box says whether the
        // figure can use it (figureRefusal). It is taken when the box is left or Enter is pressed, and only whole: "40" is
        // never the 4 it passes through. A number the figure cannot use stays in the box, marked, with the reason under it,
        // and the figure keeps the number it had.
        box = document.createElement('input'); box.type = 'text'; box.inputMode = 'decimal'; if (STYLE[kind][key] === null) box.placeholder = 'auto';
        const note = document.createElement('i'); note.className = 'range-note'; note.hidden = true; Fig.rows[key] = { box, note, why: '', remark: '' };
        const show = () => { box.value = fig[key] == null ? '' : String(boxUnits(options, fig[key])); }, judged = () => { const value = styleUnits(options, typedNumber(box.value)); return { value, why: figureRefusal(kind, key, value) }; };
        const take = () => { const { value, why } = judged(); if (why || value === fig[key]) return false; fig[key] = value; return true; };
        numbers.push({ take, mark: () => { const why = judged().why; figureSay(key, { why }); if (!why) show(); } }); show();
        box.addEventListener('input', () => figureSay(key, { why: judged().why }));
        box.addEventListener('change', () => { if (take()) changed(); }); box.addEventListener('keydown', (event) => { if (event.key === 'Enter' && take()) changed(); });
      }
      box.id = 'figure-' + key;
      if (options.unit) { const cell = document.createElement('span'), unit = document.createElement('i'); cell.className = 'unit'; unit.textContent = options.unit; cell.append(box, unit); row.appendChild(cell); } else row.appendChild(box);
      if (Fig.rows[key]) { const note = Fig.rows[key].note; note.id = box.id + '-note'; box.setAttribute('aria-describedby', note.id); row.appendChild(note); }
      holder.appendChild(row);
    }
    figureRows(); review(); figureBehind(true); $('figure-editor').hidden = false; holder.scrollTop = 0; previewFigure(); $('figure-close').focus();
  }
  /** While the editor is open the page behind it takes no clicks, keys or focus. */
  function figureBehind(off) { for (const part of document.body.children) if (part.id !== 'figure-editor' && part.tagName !== 'SCRIPT') part.inert = off; }
  function closeFigure() { const kind = Fig.open; if (!kind) return; $('figure-editor').hidden = true; Fig.open = null; figureBehind(false); const menu = $(kind === 'plot' ? 'export-plot' : 'export-heatmap'); if (menu) menu.querySelector('summary').focus(); }
  /** "Reset to defaults": the built-in look, for every box the editor shows. The figure's own settings go back to the style
   *  table and are no longer kept. The lists and tick boxes that live on the page go back to how the page opens, the file
   *  settings that all saved pictures share among them. */
  function resetFigure() {
    const kind = Fig.open; if (!kind) return; Fig[kind] = Object.assign({}, STYLE[kind]); Fig.earlier[kind] = {}; Fig.earlier.save = {}; store('step4.figure.' + kind, null); if (kind === 'plot') Fig.yKept = {};
    for (const [, key, , type, more] of FIGURE_FIELDS[kind]) {
      if (type !== 'page' || (more && more.when && !more.when())) continue;
      const box = $(key), ticked = box.type === 'checkbox', value = ticked ? 'checked' : 'value', first = ticked ? box.defaultChecked : (Array.from(box.options).find((option) => option.defaultSelected) || box.options[0]).value;
      if (box[value] !== first) { box[value] = first; box.dispatchEvent(new Event('change', { bubbles: true })); }
    }
    openFigure(kind);
  }
  /** A saved picture or frame is as wide as Options says (inches × dpi), up to 4 times its recorded pixels; its shape is the recording's. */
  const pictureScale = () => { const size = inches(), iw = S.info.width || 1; return size ? Math.max(0.25, Math.min(4, size[0] * Save.dpi / iw)) : (Save.dpi >= 600 ? 3 : Save.dpi >= 300 ? 2 : 1); };
  const pictureOptions = () => ({ format: Save.format, dpi: Save.dpi, compression: Save.compression });
  /** What a saved heatmap shows that its settings alone do not say: the cell of each row from the top, and the colour range used. */
  function heatmapRecord() { const cells = heatmapCells(), range = heatmapRange(cells.steady, Fig.heatmap); return { rows: cells.rows, colour_range: [range.low, range.high], cell_names: S.rowNames, excluded_hidden: S.hideExcluded }; }
  /** A figure's settings, for the record beside it: every one, but for a plot other than "Around the cell" without that plot's
   *  own (AROUND_STYLE), so that its record is as it was before that plot. */
  const figureSettings = (kind) => kind === 'plot' && S.plotMode !== 'around' ? Object.fromEntries(Object.entries(Fig.plot).filter(([key]) => !AROUND_STYLE.includes(key))) : Fig[kind];
  /** How a saved picture was made, for the file kept beside it: a figure's size and every setting of the editor; or, for a
   *  picture of the cells, the size from Options, the area it shows (picture.area; picture.barred: with a scale bar) and
   *  how the cells were drawn. */
  const pictureRecord = (picture) => Object.assign(picture.figure ? { size: Fig[picture.figure].width + 'x' + Fig[picture.figure].height, figure: Object.assign({ kind: picture.figure }, figureSettings(picture.figure)) } : Object.assign({ size: Save.size }, areaRecord(picture.area, picture.barred), outlineRecord()),
    picture.figure === 'heatmap' ? { heatmap: heatmapRecord() } : {}, pictureOptions());
  /** The file kept beside a saved picture, made once the picture is drawn: what the page showed (exportMetadata; before an
   *  analysis, the view and the frame), then what the picture itself shows in place of that (a figure's own ranges and means,
   *  figureShows; picture.shows: the two frames of a pair, not the frame on screen), then how the picture was made. */
  const savedRecord = (picture) => Object.assign(S.results ? exportMetadata(primaryShown()) : viewMetadata(), picture.figure ? figureShows(picture.figure) : {}, picture.shows, { picture: pictureRecord(picture) });
  const viewMetadata = () => ({ schema: 'simpleca.view-export', schema_version: 1, recording: S.info.recording, plugin_version: S.info.pluginVersion, view: S.view, frame: S.frame, outlines: S.outlines });
  /** What the colours of the picture mean in the view on screen, then what the outlines mean (outlineKey). Each: [colour, text, is an outline, dashed]. */
  const pictureKey = () => (viewColours().length ? viewColours() : [['#ffffff', 'as recorded']]).concat(outlineKey());
  /** What a saved picture or frame shows of the cells, for the file kept beside it: whether the outlines were drawn, their
   *  colours by what they mean (as its key says it), and whether the cells were numbered. */
  const outlineRecord = () => ({ outlines_shown: S.outlines, outline_colours: Object.fromEntries(outlineKey().map((entry) => [entry[1].replace(/ /g, '_'), entry[0]])), cells_numbered: S.outlines && S.numberCells });
  /** The length of the scale bar of a picture "microns" wide: the longest of 1, 2, 5, 10, 20, 50 … µm that is at most a fifth of it. */
  function scaleBarLength(microns) { const most = microns / 5, power = Math.pow(10, Math.floor(Math.log10(most) + 1e-9)); return (most >= 5 * power ? 5 : most >= 2 * power ? 2 : 1) * power; }
  /** The scale bar of a saved picture of an area (of the whole picture when none is given), in µm; null when the recording does not say how large a pixel is. */
  const scaleBar = (area) => (S.info.micronsPerPixel > 0 ? scaleBarLength((area ? area.width : S.info.width) * S.info.micronsPerPixel) : null);
  /** The area of a saved picture, for the file kept beside it: the whole picture; or the area's left top corner and its size
   *  in pixels of the recording and, where the recording says how large a pixel is, in micrometres, with the zoom it was
   *  seen at. With it the length of the picture's scale bar: a picture of an area has one, and so has one that is "barred"
   *  (a pair of frames). */
  function areaRecord(area, barred) {
    if (!area) return Object.assign({ area: 'whole picture' }, barred ? { scale_bar_um: scaleBar(null) } : null);
    const pixel = S.info.micronsPerPixel, box = (by) => ({ left: +(area.x * by).toFixed(3), top: +(area.y * by).toFixed(3), width: +(area.width * by).toFixed(3), height: +(area.height * by).toFixed(3) });
    return { area: 'area on screen', area_pixels: box(1), area_um: pixel > 0 ? box(pixel) : null, zoom: +area.zoom.toFixed(2), scale_bar_um: scaleBar(area) };
  }
  const VIEW_NAMES = { picture: 'colour overview', raw: 'as recorded', flat: 'background removed', green: 'GFP', magenta: 'response in magenta', change: 'change from baseline' };
  /** The view on screen at frame f, at its full size: a frame of the recording, fetched whole, or the view's still picture. */
  async function viewImage(f) {
    if (still()) { const image = S.view === 'background' ? Pic.background : Pic.base; if (!image) throw new Error('The picture is not loaded yet.'); return image; }
    const reply = await fetch('api/frame?f=' + f + '&mode=' + S.view + '&v=' + encodeURIComponent(Buf.version || ''), { cache: 'no-store' }); if (!reply.ok) throw new Error(reply.status === 410 && (await reply.json().catch(() => ({}))).error || 'Frame ' + f + ' could not be read.');
    return createImageBitmap(await reply.blob());
  }
  /** How a saved picture of the cells is measured out, in one place, for a picture, a frame of a set and a pair of frames.
   *  Lengths are in pixels of the recording as a picture of the whole field has them; "scale" turns them into pixels of the
   *  file (pictureScale). A picture of an area is as wide as one of the whole field: the area is enlarged, while its labels,
   *  its key and the pens of its outlines keep the size they have in every saved picture.
   *    part        the part of the recording shown (all of it when no area is given); cut: whether that is an area
   *    wide, high  its size as drawn; by: how far it is enlarged to get there
   *    unit, size, gap, font   what labels are made of: a pen's pixel, the height of text, a row of the key, the page's font */
  function sheet(area) {
    const iw = S.info.width, scale = pictureScale(), unit = Math.max(1, iw / 1388), part = area || { x: 0, y: 0, width: iw, height: S.info.height };
    return { scale, part, cut: !!area, by: iw / part.width, wide: iw, high: area ? Math.round(iw * scale * part.height / part.width) / scale : part.height,
      unit, size: Math.round(20 * unit), gap: Math.round(30 * unit), font: getComputedStyle(document.body).fontFamily };
  }
  /** A black canvas for a sheet of wide × high, to draw on in the sheet's own lengths. */
  function sheetCanvas(at, wide, high) {
    const canvas = document.createElement('canvas'); canvas.width = Math.round(wide * at.scale); canvas.height = Math.round(high * at.scale);
    if (canvas.width * canvas.height > FIGURE_MOST) throw new Error('This picture would be ' + canvas.width + ' × ' + canvas.height + ' pixels. Choose a smaller size or resolution in Options.');
    const ctx = canvas.getContext('2d'); ctx.scale(at.scale, at.scale); ctx.fillStyle = '#000'; ctx.fillRect(0, 0, wide, high); ctx.lineJoin = 'round';
    return { canvas, ctx };
  }
  /** Draws one picture's part of the recording with its top left corner at (x, y) of a sheet: the recording's own pixels,
   *  then the cells as the page draws them (paintCells: the outlines shown, in their colours or plain, and their numbers)
   *  with the pens of a saved picture. The whole picture is drawn as it always was. An area ends cleanly at its edges, and
   *  its pixels are squares, not smoothed, wherever the saved file enlarges them (more of its pixels than the recording's
   *  across the area), whatever the zoom on screen. "framed": the whole picture is cut at its edges too, as a panel beside
   *  another must be. */
  function paintView(ctx, image, at, x, y, framed) {
    const part = at.part, fx = image.width / S.info.width, fy = image.height / S.info.height;
    ctx.save(); ctx.translate(x, y);
    if (at.cut || framed) { ctx.beginPath(); ctx.rect(0, 0, at.wide, at.high); ctx.clip(); }
    ctx.imageSmoothingEnabled = !(at.cut && at.wide * at.scale > part.width * fx + 1e-6);
    ctx.drawImage(image, part.x * fx, part.y * fy, part.width * fx, part.height * fy, 0, 0, at.wide, at.high);
    ctx.scale(at.by, at.by); ctx.translate(-part.x, -part.y); paintCells(ctx, PENS.saved, at.unit / at.by, -1, at.cut ? part : undefined);
    ctx.restore();
  }
  /** The scale bar of a sheet, measured: its length in µm (scaleBar), how long that is on the sheet, the words written over
   *  it, their font, and the width of the two together. null where the recording does not say how large a pixel is. */
  function sheetBar(at) {
    const length = scaleBar(at.cut ? at.part : null); if (!length) return null;
    const measure = document.createElement('canvas').getContext('2d'), font = '600 ' + Math.round(at.size * 1.1) + 'px ' + at.font, text = +length.toPrecision(2) + ' µm', bar = length / S.info.micronsPerPixel * at.by;
    measure.font = font; return { length, bar, text, font, wide: Math.max(bar, measure.measureText(text).width) };
  }
  /** Draws a scale bar (sheetBar) in white: its length in words, their middle at "middle", over a bar of that length; the
   *  two end at "right". */
  function paintScaleBar(ctx, at, bar, right, middle) {
    ctx.font = bar.font; ctx.textBaseline = 'middle'; ctx.textAlign = 'center'; ctx.fillStyle = '#fff';
    ctx.fillText(bar.text, right - bar.wide / 2, middle); ctx.fillRect(right - bar.wide / 2 - bar.bar / 2, middle + at.size * 0.67, bar.bar, 6 * at.unit); ctx.textAlign = 'left';
  }
  /** A key laid out for a strip "wide" across: its entries in rows from the left, a new row where one would not fit. A
   *  scale bar ("bar", sheetBar), when the picture has one, takes the right end of the first row, and the entries stay
   *  clear of it: in a strip too narrow for the first entry beside the bar, the entries start on the row under it.
   *  With the places of the entries, the strip's height: the rows, and a line under them for the recording's name. */
  function keyLayout(at, key, wide, bar) {
    const measure = document.createElement('canvas').getContext('2d'), left = 14 * at.unit, beside = bar ? bar.wide + left + at.size : 0, fits = (x, width, row) => x + width <= wide - 8 - (row ? 0 : beside); measure.font = at.size + 'px ' + at.font;
    const places = []; let x = left, row = 0;
    for (const item of key) { const width = at.size * 1.4 + measure.measureText(item[1]).width + at.size * 1.4; if (!fits(x, width, row) && (x > left || (bar && !row))) { x = left; row++; } places.push([x, row]); x += width; }
    return { key, places, bar, wide, height: Math.round((row + 1) * at.gap + at.gap * 1.5) };
  }
  /** Draws a key (keyLayout) under a picture that ends at y: a filled square for a colour of the view, an open one for an
   *  outline, dashed where dashes alone tell it from another outline of its colour; under them the recording's name; and
   *  on the right, within the height of the first row, the picture's scale bar when it has one. */
  function paintKey(ctx, at, laid, y) {
    const { unit, size, gap } = at, key = laid.key; ctx.font = size + 'px ' + at.font; ctx.textBaseline = 'middle';
    key.forEach((item, i) => { const kx = laid.places[i][0], ky = y + gap * 0.9 + laid.places[i][1] * gap;
      if (item[2]) { ctx.setLineDash(item[3] && key.some((other) => other !== item && other[2] && other[0] === item[0]) ? [4 * unit, 3 * unit] : []); ctx.strokeStyle = item[0]; ctx.lineWidth = 2 * unit; ctx.strokeRect(kx, ky - size * 0.45, size * 0.9, size * 0.9); ctx.setLineDash([]); }
      else { ctx.fillStyle = item[0]; ctx.fillRect(kx, ky - size * 0.45, size * 0.9, size * 0.9); }
      ctx.fillStyle = '#fff'; ctx.fillText(item[1], kx + size * 1.4, ky); });
    ctx.fillStyle = '#AEB8C2'; ctx.font = Math.round(size * 0.8) + 'px ' + at.font; ctx.fillText(S.info.recording, 14 * unit, y + laid.height - gap * 0.55);
    if (laid.bar) paintScaleBar(ctx, at, laid.bar, laid.wide - 14 * unit, y + gap * 0.9 - size * 0.3);
  }
  /** A frame's number and time, as its label says them: "Frame 22 · 10.5 s". */
  const frameText = (f) => { const time = frameTimes()[f - 1]; return 'Frame ' + f + (finite(time) ? ' · ' + time.toFixed(1) + ' s' : ''); };
  /** The labels of a frame of a set, white on dark boxes: on its top left corner its number, its time and the view; on its
   *  bottom right corner its scale bar ("bar", sheetBar), when it has one. */
  function paintFrameLabels(ctx, at, f, bar) {
    const { unit, size } = at, text = frameText(f) + ' · ' + (VIEW_NAMES[S.view] || S.view);
    ctx.font = '600 ' + Math.round(size * 1.1) + 'px ' + at.font; ctx.textBaseline = 'middle'; const tw = ctx.measureText(text).width;
    ctx.fillStyle = 'rgba(0,0,0,0.7)'; ctx.fillRect(10 * unit, 10 * unit, tw + 20 * unit, size * 1.9); ctx.fillStyle = '#fff'; ctx.fillText(text, 20 * unit, 10 * unit + size * 0.98);
    if (!bar) return;
    const high = size * 1.65 + 14 * unit, top = at.high - 10 * unit - high;
    ctx.fillStyle = 'rgba(0,0,0,0.7)'; ctx.fillRect(at.wide - 30 * unit - bar.wide, top, bar.wide + 20 * unit, high); paintScaleBar(ctx, at, bar, at.wide - 20 * unit, top + size * 0.98);
  }
  /** A saved picture of the cells in the view on screen (at the frame on screen when none is given): the whole picture, or
   *  the area given (savedArea), with the cells on it (paintView). Under a picture, a key to the colours. A frame of a set
   *  ("labelled") has its number, time and view on its corner instead. A picture or frame of an area has a scale bar: at
   *  the right of the key, or on the frame's bottom right corner. One ordinary colour picture. */
  async function viewCanvas(frame, labelled, area) {
    const f = frame || S.frame, image = await viewImage(f), at = sheet(area), bar = at.cut ? sheetBar(at) : null, key = labelled ? null : keyLayout(at, pictureKey(), at.wide, bar);
    const { canvas, ctx } = sheetCanvas(at, at.wide, at.high + (key ? key.height : 0));
    paintView(ctx, image, at, 0, 0);
    if (key) paintKey(ctx, at, key, at.high); else paintFrameLabels(ctx, at, f, bar);
    return canvas;
  }
  /** How "Baseline and peak" is measured out: two panels of one size side by side, a strip over them that says what each is,
   *  and the key under both, with the scale bar. Every edge lies on a whole pixel of the file, so each panel is as sharp as a frame saved alone. */
  function pairLayout(area) {
    const at = sheet(area), whole = (v) => Math.round(v * at.scale) / at.scale, head = whole(at.size * 2.2), second = whole(at.wide + 12 * at.unit), wide = second + at.wide, key = keyLayout(at, pictureKey(), wide, sheetBar(at));
    return { at, head, second, wide, key, high: head + at.high + key.height };
  }
  /** Baseline and peak: two frames of the view on screen in one picture, for showing a response. Each panel is drawn as a
   *  saved frame is (paintView: the same area, outlines and numbers), and the plugin shows every frame of a view between the
   *  same gray levels, so the two can be compared. Over each panel what it is, with its frame and time; under both the key,
   *  once, with the scale bar at its right. */
  async function pairCanvas(baseline, peak, area) {
    const images = [await viewImage(baseline), await viewImage(peak)], laid = pairLayout(area), at = laid.at, { canvas, ctx } = sheetCanvas(at, laid.wide, laid.high), left = 14 * at.unit;
    [['baseline', baseline, 0], ['peak', peak, laid.second]].forEach(([name, f, x], i) => {
      paintView(ctx, images[i], at, x, laid.head, true);
      ctx.textBaseline = 'middle'; ctx.font = '600 ' + Math.round(at.size * 1.1) + 'px ' + at.font; ctx.fillStyle = '#fff'; ctx.fillText(name, x + left, laid.head / 2);
      const named = ctx.measureText(name).width; ctx.font = at.size + 'px ' + at.font; ctx.fillStyle = '#AEB8C2'; ctx.fillText(frameText(f), x + left + named + at.size * 0.8, laid.head / 2); });
    paintKey(ctx, at, laid.key, laid.head + at.high);
    return canvas;
  }
  /** The two frames "Baseline and peak" opens on. Baseline: the last of the F₀ frames. Peak: the selected cell's peak in the
   *  response window, as its call has it; with no such cell, the frame of the response window where the mean ΔF/F₀ of the
   *  cells that responded is highest (of the cells that did not, when none responded). Before an analysis, the frame on
   *  screen for both. */
  function pairFrames() {
    const r = S.results; if (!r) return { baseline: S.frame, peak: S.frame };
    const baseline = r.baselineFrames ? r.baselineFrames[1] : S.frame, own = number(callRow(S.cell).Peak_s); if (finite(own)) return { baseline, peak: nearestFrame(own) };
    const names = S.info.names.filter((name) => r.traces[name]), yes = names.filter((name) => callOf(name) === 'Responded'), mean = groupMean(yes.length ? yes : names.filter((name) => callOf(name) === 'Did not respond'), r.traces).mean;
    const [from, to] = responseSpans().response; let peak = S.frame, top = -Infinity;
    r.time.forEach((t, i) => { if (t >= from - 1e-9 && t <= to + 1e-9 && finite(mean[i]) && mean[i] > top) { top = mean[i]; peak = i + 1; } });
    return { baseline, peak };
  }
  /** Where n frames of w x h go in one picture: a grid about as wide as it is tall, each tile at most its own size and the
   *  whole at most COMBINED_MAX pixels wide. Returns { cols, rows, tw, th, gap, width, height }. */
  const COMBINED_MAX = 6000;
  function combinedLayout(n, w, h) {
    const cols = Math.max(1, Math.ceil(Math.sqrt(n * h / w) * w / h)), rows = Math.ceil(n / cols), gap = Math.max(2, Math.round(w / 200));
    const tw = Math.max(1, Math.min(w, Math.floor((COMBINED_MAX - (cols + 1) * gap) / cols))), th = Math.max(1, Math.round(tw * h / w));
    return { cols, rows, tw, th, gap, width: cols * tw + (cols + 1) * gap, height: rows * th + (rows + 1) * gap };
  }
  /** What a saved picture of the cells is called in its file's name: the view; for a view of the recording the frame on screen,
   *  or the frames given ("baseline_frame5_peak_frame12"); and "area" when it shows the area on screen. */
  const viewWord = () => ({ picture: 'colour_overview', background: 'background_areas', raw: 'as_recorded', flat: 'background_removed', green: 'gfp', magenta: 'response_in_magenta', change: 'change_from_baseline' }[S.view] || S.view);
  const viewKind = (area, frames) => 'picture_' + viewWord() + (still() ? '' : '_' + (frames || 'frame' + S.frame)) + (area ? '_area' : '');
  const badSave = (error) => /fetch/i.test(error.message) ? 'Fiji is not connected. Keep Fiji open, or reopen the analysis from Fiji.' : error.message;
  /** Download through the browser: the same file Save writes (format, size, resolution), plus how it was made. */
  async function save(picture, kind, noteId) {
    const note = $(noteId); note.className = 'save-note'; note.textContent = 'Preparing…';
    try {
      const png = await picture(), metadata = JSON.stringify(savedRecord(picture), null, 2);   // the record of what was drawn: after the drawing
      const reply = await api('convert', Object.assign({ png }, pictureOptions())), name = exportName(kind, reply.extension);
      const bytes = Uint8Array.from(atob(reply.base64), (c) => c.charCodeAt(0));
      download(new Blob([bytes], { type: reply.extension === 'tif' ? 'image/tiff' : 'image/png' }), name);
      download(new Blob([metadata], { type: 'application/json' }), name.replace(/\.(png|tif)$/, '.metadata.json')); note.textContent = '';
    } catch (error) { note.className = 'save-note bad'; note.textContent = 'Not downloaded: ' + badSave(error); }
  }
  /** Saves the picture, the numbers behind it and how it was made into the "plots" folder beside the measurement tables,
   *  and says where. Nothing goes through the browser's Downloads folder. */
  async function saveToFolder(picture, kind, noteId, csv) {
    const note = $(noteId); note.className = 'save-note'; note.textContent = 'Saving…';
    try {
      const png = await picture(), metadata = savedRecord(picture), name = exportName(kind, 'png').replace(/\.png$/, '');   // the record of what was drawn: after the drawing
      const reply = await api('export', Object.assign({ name, png, csv: csv || null, metadata }, pictureOptions()));
      note.textContent = 'Saved ' + reply.saved[0] + (csv ? ' + CSV' : '') + ' ';
      const open = document.createElement('button'); open.type = 'button'; open.textContent = 'Open folder';
      open.addEventListener('click', () => api('open', { what: 'plots' }).catch((error) => { note.className = 'save-note bad'; note.textContent = error.message; }));
      note.appendChild(open);
    } catch (error) { note.className = 'save-note bad'; note.textContent = 'Not saved: ' + badSave(error); }
  }
  const stored = (key) => { try { return window.localStorage.getItem(key); } catch (blocked) { return null; } };
  // The page's address is new every time it opens, so the browser's own storage would forget; Fiji keeps these too.
  const store = (key, value) => { try { if (value === null) window.localStorage.removeItem(key); else window.localStorage.setItem(key, value); } catch (blocked) { } api('prefs', { key, value }).catch(() => { }); };
  const Layout = { apply: {} };
  /** A draggable edge: "read" turns the pointer position into a size, "apply" sets it (null resets). The size is remembered. */
  function dragEdge(handle, read, apply, key) {
    Layout.apply[key] = apply; const saved = stored(key); if (saved) apply(saved);
    handle.addEventListener('pointerdown', (event) => {
      event.preventDefault(); try { handle.setPointerCapture(event.pointerId); } catch (unsupported) { } handle.classList.add('on'); document.body.classList.add('dragging');
      const move = (e) => { const value = read(e); if (value) { apply(value); store(key, value); redraw(); } };
      const done = () => { handle.classList.remove('on'); document.body.classList.remove('dragging'); handle.removeEventListener('pointermove', move); handle.removeEventListener('pointerup', done); handle.removeEventListener('pointercancel', done); redraw(); };
      handle.addEventListener('pointermove', move); handle.addEventListener('pointerup', done); handle.addEventListener('pointercancel', done);
    });
    handle.addEventListener('dblclick', () => { apply(null); store(key, null); redraw(); });
  }
  function wire() {
    const helpOn = stored('step4.help') === '1'; document.body.classList.toggle('show-help', helpOn); $('toggle-help').textContent = helpOn ? 'Hide explanations' : 'Show explanations';
    $('toggle-help').addEventListener('click', () => { const on = document.body.classList.toggle('show-help'); store('step4.help', on ? '1' : '0'); $('toggle-help').textContent = on ? 'Hide explanations' : 'Show explanations'; redraw(); });
    dragEdge($('resize-side'), (e) => Math.round(Math.max(260, Math.min(window.innerWidth * 0.6, e.clientX - $('side').getBoundingClientRect().left - 10))) + 'px',
      (value) => { if (value) document.documentElement.style.setProperty('--side-w', value); else document.documentElement.style.removeProperty('--side-w'); }, 'step4.side');
    dragEdge($('resize-height'), (e) => Math.round(Math.max(240, Math.min(1600, e.clientY - document.querySelector('.split').getBoundingClientRect().top - 6))) + 'px',
      (value) => { const split = document.querySelector('.split'); split.classList.toggle('sized', !!value); if (value) split.style.setProperty('--split-h', value); else split.style.removeProperty('--split-h'); }, 'step4.height');
    dragEdge($('resize-split'), (e) => { const box = document.querySelector('.split').getBoundingClientRect(); return Math.max(20, Math.min(80, (e.clientX - box.left - 7) / box.width * 100)).toFixed(1) + '%'; },
      (value) => { const split = document.querySelector('.split'); if (value) split.style.setProperty('--pic-w', value); else split.style.removeProperty('--pic-w'); }, 'step4.split');
    for (const button of document.querySelectorAll('[data-plot]')) button.addEventListener('click', () => setPlotMode(button.dataset.plot));
    { const button = document.querySelector('[data-plot="around"]'); button.dataset.title = button.title; }   // its tooltip says why while it is off
    for (const [id, end] of [['y-min', 'lo'], ['y-max', 'hi']]) $(id).addEventListener('input', () => {
      const value = typedNumber($(id).value), now = Plot.y[yKey()] || (Plot.y[yKey()] = {});
      now[end] = finite(value) ? value : null; Plot.dirty = true; drawPlot();
    });
    $('y-auto').addEventListener('click', () => { delete Plot.y[yKey()]; showY(); Plot.dirty = true; drawPlot(); });
    const setSignal = (value) => { if (value !== S.signal) { if (yKind(S.plotMode) === 'signal') Fig.plot.yTitle = Fig.plot.yMin = Fig.plot.yMax = Fig.plot.yStep = null; else delete Fig.yKept.signal; }   // typed for the other signal
      S.signal = Data.signal = value; $('signal').value = value; $('data-signal').value = value; showY(); Plot.dirty = true; drawPlot(); if (S.results && !$('panel-data').hidden) showData(); };
    $('signal').addEventListener('change', () => setSignal($('signal').value));
    $('save-plot-folder').addEventListener('click', () => saveToFolder(plotPicture, plotKind(), 'plot-saved', plotNumbers()));
    $('save-heatmap-folder').addEventListener('click', () => saveToFolder(heatmapPicture, 'heatmap', 'heatmap-saved'));
    // A picture of the cells: the frame on screen, and the area "Area" says, both as they are when its entry is clicked.
    const cellsPicture = (how) => () => { const area = savedArea(), picture = async () => (await viewCanvas(null, false, area)).toDataURL('image/png'); picture.area = area; how(picture, viewKind(area), 'picture-saved'); };
    $('save-picture-folder').addEventListener('click', cellsPicture(saveToFolder)); $('save-picture').addEventListener('click', cellsPicture(save));
    $('save-area').addEventListener('change', () => setSaveArea($('save-area').value));
    $('outline-width').addEventListener('change', () => { S.outlineWidth = Number($('outline-width').value) || 1; store('step4.outline', $('outline-width').value); drawPicture(); });
    $('outline-colours').addEventListener('change', () => { S.outlineColours = $('outline-colours').value === 'plain' ? 'plain' : 'call'; store('step4.plain', S.outlineColours === 'plain' ? '1' : '0'); drawPicture(); });
    $('number-cells').addEventListener('change', () => { S.numberCells = $('number-cells').checked; store('step4.numbers', S.numberCells ? '1' : '0'); drawPicture(); });
    // One Export menu per panel holds its save buttons. Choosing an entry closes the menu; so does a click elsewhere or Escape.
    // A panel opened from an entry (Options, Save frames) is placed against the menu's own button, which stays on screen:
    // its right edge on the button's, inside the window; above the button in the lower half of the window, under it in the upper.
    const anchorOf = (button) => button.closest('details.export') || button;
    const openBeside = (panel, entry) => { panel.hidden = false; const box = anchorOf(entry).getBoundingClientRect(), wide = panel.offsetWidth;
      panel.style.left = Math.max(8, Math.min(window.innerWidth - wide - 16, box.right - wide)) + 'px'; panel.style.top = ''; panel.style.bottom = '';
      if (box.top > window.innerHeight / 2) panel.style.bottom = (window.innerHeight - box.top + 6) + 'px'; else panel.style.top = (box.bottom + 6) + 'px'; };
    const menus = Array.from(document.querySelectorAll('details.export, details.quality-chip'));
    for (const menu of menus) {
      menu.addEventListener('click', (event) => { if (event.target.closest('.export-menu button')) menu.open = false; }, true);
      menu.addEventListener('toggle', () => {
        if (!menu.open) return;
        for (const other of menus) if (other !== menu) other.open = false;
        const pop = menu.querySelector('.quality-pop'); if (!pop) return;   // the chip can sit anywhere in its row: its list is kept inside the work area
        const area = document.querySelector('main').getBoundingClientRect(), chip = menu.getBoundingClientRect(), width = Math.min(520, area.width - 16);
        pop.style.width = width + 'px'; pop.style.left = Math.round(Math.max(area.left + 8, Math.min(chip.left, area.right - 8 - width)) - chip.left) + 'px';
      });
    }
    // the path as the click travelled: an entry that redraws itself when clicked is no longer inside the menu by now
    document.addEventListener('click', (event) => { const path = event.composedPath(); for (const menu of menus) if (menu.open && !path.includes(menu)) menu.open = false; });
    document.addEventListener('keydown', (event) => { if (event.key === 'Escape') for (const menu of menus) if (menu.open) { menu.open = false; menu.querySelector('summary').focus(); event.preventDefault(); } });   // (an Escape that closed a menu does nothing else)
    // Save frames…, for a view of the recording, in two ways (the two buttons at its top). Every Nth frame: each its own
    // labelled picture, in a folder of their own, with one picture of them all. Baseline and peak: those two frames side by
    // side in one picture (pairCanvas). Its two boxes open on the frames suggested (pairFrames); a frame typed over a
    // suggestion is kept, for as long as the suggestion stays the same. A box that does not hold a frame of the recording
    // keeps what was typed, is marked, and says beside it which frames there are; nothing is saved meanwhile, and the
    // panel opens on the suggestion again. "Every" goes by the same rule: a number of frames it cannot use, also "2.5", is
    // not rounded to one it can.
    const framesPanel = $('frames-options'), framesEvery = () => typedWhole($('frames-every').value), framesWay = () => framesPanel.querySelector('[data-frames].on').dataset.frames;
    const fileWords = () => (Save.format === 'tiff' ? 'TIFF ' + Save.compression.toUpperCase() : 'PNG') + ', ' + Save.dpi + ' dpi.';
    const Pair = { suggested: {}, typed: {} }, pairBoxes = { baseline: 'frames-baseline', peak: 'frames-peak' }, pairFrame = (key) => { const frame = typedWhole($(pairBoxes[key]).value); return frame >= 1 && frame <= S.info.frames ? frame : NaN; };
    /** Marks a box of the panel that holds no number it can use, and says beside it which numbers there are. */
    const framesMark = (id, bad) => { const box = $(id), range = $(id + '-range'); if (bad) box.setAttribute('aria-invalid', 'true'); else box.removeAttribute('aria-invalid'); range.hidden = !bad; range.textContent = '1 to ' + S.info.frames; };
    const framesCount = () => {
      const pair = framesWay() === 'pair', go = $('frames-go');
      $('frames-every-row').hidden = pair; $('frames-baseline-row').hidden = $('frames-peak-row').hidden = $('frames-download').hidden = !pair;
      go.textContent = pair ? 'Save picture' : 'Save frames'; go.title = pair ? 'Saves the picture and how it was made into the plots folder' : 'Saves the frames into a folder of their own inside the plots folder';
      if (!pair) { const every = framesEvery(), ok = every >= 1 && every <= S.info.frames; go.disabled = !ok; framesMark('frames-every', !ok);
        $('frames-count').textContent = ok ? plural(Math.ceil(S.info.frames / every), 'picture') + (Math.ceil(S.info.frames / every) > 1 ? ' and one with all of them' : '') + ', ' + fileWords() : ''; return; }
      let ok = true;
      for (const key of Object.keys(pairBoxes)) { const bad = !finite(pairFrame(key)); ok = ok && !bad; framesMark(pairBoxes[key], bad); }
      go.disabled = $('frames-download').disabled = !ok;
      const laid = ok ? pairLayout(savedArea()) : null; $('frames-count').textContent = ok ? 'One picture, ' + Math.round(laid.wide * laid.at.scale) + ' × ' + Math.round(laid.high * laid.at.scale) + ' px, ' + fileWords() : '';
    };
    $('frames-every').addEventListener('input', framesCount); $('frames-cancel').addEventListener('click', () => { framesPanel.hidden = true; });
    for (const key of Object.keys(pairBoxes)) $(pairBoxes[key]).addEventListener('input', () => { const typed = $(pairBoxes[key]).value; if (finite(pairFrame(key)) && typed.trim() !== String(Pair.suggested[key])) Pair.typed[key] = typed; else delete Pair.typed[key]; framesCount(); });
    for (const button of framesPanel.querySelectorAll('[data-frames]')) button.addEventListener('click', () => {
      for (const other of framesPanel.querySelectorAll('[data-frames]')) { other.classList.toggle('on', other === button); other.setAttribute('aria-pressed', other === button ? 'true' : 'false'); } framesCount(); });
    $('frames-area').addEventListener('change', () => { setSaveArea($('frames-area').value); framesCount(); });
    framesPanel.addEventListener('keydown', (event) => { if (event.key === 'Enter' && event.target.tagName === 'INPUT' && !$('frames-go').disabled) $('frames-go').click(); if (event.key === 'Escape') { framesPanel.hidden = true; event.preventDefault(); } });   // (an Escape that closed the panel does nothing else)
    $('save-frames').addEventListener('click', (event) => {
      const note = $('picture-saved'); note.className = 'save-note';
      if (still()) { note.className = 'save-note bad'; note.textContent = 'Choose a recording view first: this view is one still picture.'; return; }
      if (!framesPanel.hidden) { framesPanel.hidden = true; return; }
      $('frames-every').value = stored('step4.every') || '50';
      const suggested = pairFrames(); for (const key of Object.keys(pairBoxes)) { if (Pair.suggested[key] !== suggested[key]) delete Pair.typed[key]; $(pairBoxes[key]).value = Pair.typed[key] === undefined ? suggested[key] : Pair.typed[key]; } Pair.suggested = suggested;
      framesCount(); openBeside(framesPanel, $('save-frames'));
      const first = $(framesWay() === 'pair' ? 'frames-baseline' : 'frames-every'); first.focus(); first.select(); event.stopPropagation();
    });
    document.addEventListener('click', (event) => { if (!framesPanel.hidden && !framesPanel.contains(event.target) && event.target.id !== 'save-frames') framesPanel.hidden = true; });
    // The same area for every frame of a set, also when the view is moved while they are saved: the area is taken once, as saving begins.
    const saveEvery = async () => {
      const note = $('picture-saved'), every = framesEvery(), area = savedArea(); note.className = 'save-note'; framesPanel.hidden = true;
      store('step4.every', String(every));
      const frames = []; for (let f = 1; f <= S.info.frames; f += every) frames.push(f);
      const stamp = exportName('x', 'png').match(/(\d{8}-\d{6})\.png$/)[1], folder = 'frames_' + viewWord() + (area ? '_area' : '') + '_every_' + every + '_' + stamp, width = String(S.info.frames).length;
      $('save-frames').disabled = true;
      try {
        let reply = null, combined = null, layout = null;
        for (let i = 0; i < frames.length; i++) {
          note.textContent = 'Saving frame ' + frames[i] + ' (' + (i + 1) + ' of ' + frames.length + ')…';
          const canvas = await viewCanvas(frames[i], true, area);
          reply = await api('export', Object.assign({ folder, name: 'frame_' + String(frames[i]).padStart(width, '0'), png: canvas.toDataURL('image/png') }, pictureOptions()));
          // the same labelled frame, as one tile of the combined picture
          if (!combined) { layout = combinedLayout(frames.length, canvas.width, canvas.height); combined = document.createElement('canvas'); combined.width = layout.width; combined.height = layout.height; const c = combined.getContext('2d'); c.fillStyle = '#000'; c.fillRect(0, 0, layout.width, layout.height); c.imageSmoothingQuality = 'high'; }
          combined.getContext('2d').drawImage(canvas, layout.gap + (i % layout.cols) * (layout.tw + layout.gap), layout.gap + Math.floor(i / layout.cols) * (layout.th + layout.gap), layout.tw, layout.th);
        }
        if (combined && frames.length > 1) { note.textContent = 'Saving the combined picture…'; await api('export', Object.assign({ folder, name: 'all_frames_combined', png: combined.toDataURL('image/png') }, pictureOptions())); }
        await api('export', { folder, name: 'about_these_frames', metadata: Object.assign(viewMetadata(), { frames, every, combined_picture: combined && frames.length > 1 ? { file: 'all_frames_combined', columns: layout.cols, rows: layout.rows, tile_pixels: [layout.tw, layout.th] } : null, picture: pictureOptions() }, areaRecord(area), outlineRecord()) });
        note.textContent = 'Saved ' + plural(frames.length, 'frame') + (frames.length > 1 ? ' and one combined picture' : '') + ' in plots/' + reply.folderName + ' ';
        const open = document.createElement('button'); open.type = 'button'; open.textContent = 'Open plots folder'; open.addEventListener('click', () => api('open', { what: 'plots' }).catch(() => { })); note.appendChild(open);
      } catch (error) { note.className = 'save-note bad'; note.textContent = 'Stopped: ' + badSave(error); }
      finally { $('save-frames').disabled = false; }
    };
    // The pair is one picture: saved or downloaded as a picture of the cells is ("how": saveToFolder or save), in the format and
    // resolution of Options. Its record names both frames, in place of the frame on screen, and the gray levels the view's
    // frames are shown between.
    const savePair = async (how) => {
      const note = $('picture-saved'), baseline = pairFrame('baseline'), peak = pairFrame('peak'), area = savedArea(), times = frameTimes(), named = (f) => ({ frame: f, time_s: finite(times[f - 1]) ? times[f - 1] : null });
      framesPanel.hidden = true; note.className = 'save-note'; note.textContent = 'Preparing…';
      try {
        const shown = await api('range?mode=' + S.view), picture = async () => (await pairCanvas(baseline, peak, area)).toDataURL('image/png');
        picture.area = area; picture.barred = true; picture.shows = { frame: undefined, time_s: undefined, baseline: named(baseline), peak: named(peak), brightness_range: shown.ranges };
        await how(picture, viewKind(area, 'baseline_frame' + baseline + '_peak_frame' + peak), 'picture-saved');
      } catch (error) { note.className = 'save-note bad'; note.textContent = 'Not saved: ' + badSave(error); }
    };
    $('frames-go').addEventListener('click', () => (framesWay() === 'pair' ? savePair(saveToFolder) : saveEvery())); $('frames-download').addEventListener('click', () => savePair(save));
    // how pictures are saved: format, size, resolution, TIFF compression (kept with Fiji's preferences)
    for (const key of Object.keys(saveBoxes)) $(saveBoxes[key]).addEventListener('change', () => { Save[key] = key === 'dpi' ? Number($(saveBoxes[key]).value) : $(saveBoxes[key]).value; store('step4.' + key, String(Save[key])); showSaveOptions(); });
    for (const button of document.querySelectorAll('.save-options-toggle')) button.addEventListener('click', (event) => {
      const panel = $('save-options'), open = panel.hidden || panel.dataset.by !== String(Array.from(document.querySelectorAll('.save-options-toggle')).indexOf(button));
      panel.hidden = !open; if (!open) return;
      panel.dataset.by = String(Array.from(document.querySelectorAll('.save-options-toggle')).indexOf(button)); showSaveOptions(); openBeside(panel, button);
      event.stopPropagation();
    });
    document.addEventListener('click', (event) => { const panel = $('save-options'); if (!panel.hidden && !panel.contains(event.target) && !event.target.classList.contains('save-options-toggle')) panel.hidden = true; });
    // the figure editor, opened from the plot's and the heatmap's Export menu: its preview is the picture that gets saved
    $('edit-plot').addEventListener('click', () => openFigure('plot')); $('edit-heatmap').addEventListener('click', () => openFigure('heatmap'));
    $('figure-close').addEventListener('click', closeFigure); $('figure-reset').addEventListener('click', resetFigure);
    $('figure-save').addEventListener('click', () => Fig.open === 'plot' ? saveToFolder(plotPicture, plotKind(), 'figure-saved', plotNumbers()) : saveToFolder(heatmapPicture, 'heatmap', 'figure-saved'));
    $('figure-download').addEventListener('click', () => Fig.open === 'plot' ? save(plotPicture, plotKind(), 'figure-saved') : save(heatmapPicture, 'heatmap', 'figure-saved'));
    $('figure-editor').addEventListener('mousedown', (event) => { if (event.target === $('figure-editor')) closeFigure(); });
    document.addEventListener('keydown', (event) => { if (Fig.open && event.key === 'Escape') { closeFigure(); event.preventDefault(); event.stopPropagation(); } }, true);   // before the page's own keys
    window.addEventListener('resize', previewFigure);
    // Every table tab has one Export menu: a CSV into the plots folder, with the record of how the page was set beside it
    // (exportMetadata), or a CSV through the browser. [tab, file name, CSV text, whether the table follows the review marks:
    // the calls table with its reviewed columns and the metric summary do, the saved tables as they are do not; whether the
    // codes of the call checks were put into its QC_flags]
    const tableExports = [['metrics', () => $('metrics-table').value === 'summary' ? 'metrics_summary' : 'table_' + $('metrics-table').value, () => csvText(metricRows()), () => $('metrics-table').value === 'summary', metricCoded],
      ['calls', () => 'responder_calls', () => csvText(callsRows()), () => true, callsCoded], ['counts', () => 'responder_counts', () => csvText(savedRows('counts')), () => false, () => false], ['data', dataKind, dataCsv, () => false, () => false]];
    for (const [tab, kind, csv, followsReview, coded] of tableExports) {
      $(tab + '-download').addEventListener('click', () => download(new Blob([csv()], { type: 'text/csv' }), exportName(kind(), 'csv')));
      $(tab + '-save').addEventListener('click', async () => {
        const note = $(tab + '-saved'); note.className = 'save-note'; note.textContent = 'Saving…';
        try {
          const reply = await api('export', { name: exportName(kind(), 'csv').replace(/\.csv$/, ''), csv: csv(), metadata: Object.assign(exportMetadata(followsReview()), codesRecord(coded())) });
          note.textContent = 'Saved in ' + reply.folderName + ': ' + reply.saved[0] + ' ';
          const open = document.createElement('button'); open.type = 'button'; open.textContent = 'Open folder'; open.addEventListener('click', () => api('open', { what: 'plots' }).catch(() => { })); note.appendChild(open);
        } catch (error) { note.className = 'save-note bad'; note.textContent = 'Not saved: ' + (/fetch/i.test(error.message) ? 'Fiji is not connected.' : error.message); }
      });
    }
    $('data-signal').addEventListener('change', () => setSignal($('data-signal').value));
    $('data-back').addEventListener('click', () => { Data.first -= Data.page; showData(); }); $('data-forward').addEventListener('click', () => { Data.first += Data.page; showData(); });
    for (const tab of document.querySelectorAll('.tab')) tab.addEventListener('click', () => showTab(tab.dataset.tab));
    $('toggle-settings').addEventListener('click', () => { const hidden = document.body.classList.toggle('no-side'), button = $('toggle-settings'); button.textContent = hidden ? '›' : '‹'; button.setAttribute('aria-label', hidden ? 'Show settings' : 'Hide settings'); button.title = hidden ? 'Show the settings' : 'Hide the settings'; redraw(); });
    $('side').addEventListener('input', (event) => {
      const box = event.target;
      if (!S.draft || box.type === 'file' || ['confirmed', 'control', 'protocol-mode'].includes(box.id)) return;
      if (box.dataset.shown !== undefined) box.classList.toggle('changed', box.value !== box.dataset.shown);
      if (box.dataset.frames) applyFrames(box.dataset.window); else if (box.id === 'sampling_hz') applyAllFrames();   // seconds always come from the frames
      draftChanged();
    });
    $('control').addEventListener('change', () => {
      const draft = collectDraft();
      if ($('control').checked) {
        const legacy = draft.protocol.schema_version === 1, first = expandWindow(draft.protocol.windows[0], legacy);
        draft.protocol.windows.push(Object.assign({}, first, { Role: 'Control', Stimulus: legacy ? 'KCl' : 'Positive control', Baseline_start_s: null, Baseline_end_s: null, Baseline_s: null, Stimulus_onset_s: null, Start_s: null, End_s: null }));
      } else draft.protocol.windows = draft.protocol.windows.slice(0, 1);
      fillDraft(draft); draftChanged();
    });
    $('protocol-mode').addEventListener('change', () => {
      const wanted = Number($('protocol-mode').value), previous = S.draft.protocol.schema_version;
      $('protocol-mode').value = previous; const draft = collectDraft();
      if (wanted === 1 && (draft.protocol.field_overrides || draft.sensitivity || draft.negative_control)) { S.note = 'Remove field timing overrides and diagnostic files before selecting the legacy rule. They require a version 2 protocol.'; showState(); return; }
      S.protocolDrafts[previous] = clone(draft.protocol);
      if (S.protocolDrafts[wanted]) draft.protocol = clone(S.protocolDrafts[wanted]);
      else {
        draft.protocol.windows = draft.protocol.windows.map((w, i) => {
          const row = expandWindow(w, previous === 1);
          if (wanted === 1) { row.Duration_mode = 'samples'; row.Baseline_end_s = row.Start_s; row.Stimulus_onset_s = row.Start_s; row.Baseline_s = row.Start_s - row.Baseline_start_s; if (i === 1) row.Stimulus = 'KCl'; }
          else row.Duration_mode = 'seconds';
          return row;
        });
        draft.protocol.schema_version = wanted; draft.protocol.validation_status = wanted === 1 ? 'legacy_unvalidated' : 'user_defined';
        draft.protocol.name = wanted === 1 ? 'Legacy consecutive-sample protocol' : 'Custom response protocol';
      }
      fillDraft(draft); draftChanged();
    });
    for (const key of FILES) {
      $(key + '-file').addEventListener('change', async () => {
        const file = $(key + '-file').files[0]; if (!file) return;
        try { if (file.size > 8 * 1024 * 1024) throw new Error('Choose a CSV smaller than 8 MB.'); const content = await file.text(); S.draft[key] = { name: file.name, content }; fileNote(key); if (key === 'timestamps') S.reapplyFrames = true; draftChanged(); }
        catch (error) { draftError(error); }
      });
      $(key + '-clear').addEventListener('click', () => { S.draft[key] = null; $(key + '-file').value = ''; fileNote(key); if (key === 'timestamps') S.reapplyFrames = true; draftChanged(); });
    }
    $('protocol-import').addEventListener('change', async () => {
      const file = $('protocol-import').files[0]; if (!file) return;
      try { let parsed; try { parsed = JSON.parse(await file.text()); } catch (bad) { throw new Error(file.name + ' is not a valid JSON file (' + bad.message + '). Choose a protocol file.'); }
        const timingOf = (protocol) => JSON.stringify([(protocol.windows || []).map((w) => [w.Role].concat(TIME_KEYS.map((k) => (w[k] == null || w[k] === '' ? null : Number(w[k]))))), protocol.field_overrides || null]);
        const result = await api('protocol', { protocol: parsed }); const draft = collectDraft(), before = timingOf(draft.protocol); draft.protocol = result.protocol; S.protocolDrafts = {}; S.note = null; S.unlinked.clear(); fillDraft(draft);
        const changed = timingOf(collectDraft().protocol) !== before; if (changed) $('confirmed').checked = false; draftChanged();
        S.imported = { revision: S.revision, text: 'Imported “' + (result.protocol.name || file.name) + '”. ' + (changed ? 'Its timing replaced yours: check it, then tick the box.' : 'The timing is unchanged.') }; showState(); }
      catch (error) { draftError(error); }
      finally { $('protocol-import').value = ''; }
    });
    $('protocol-export').addEventListener('click', async () => {
      try { const result = await api('protocol', { protocol: collectDraft().protocol }); download(new Blob([JSON.stringify(result.protocol, null, 2)], { type: 'application/json' }), exportName('response-protocol', 'json')); }
      catch (error) { draftError(error); }
    });
    $('analyze').addEventListener('click', async () => {
      if (validateDraft(collectDraft()).errors.length) { $('validation').scrollIntoView({ block: 'nearest' }); return; }
      if (!$('confirmed').checked) { S.note = 'Tick “Timing checked” first.'; showState(); $('confirmed').focus(); return; }
      if (!S.state.engineReady) { S.note = 'Choose “Set up analysis engine” once, then run the analysis.'; showState(); $('setup').focus(); return; }
      S.saving = true; S.note = null; showState();
      try { await syncDraft(); await send('analyze', { draft: collectDraft(), confirmed: true }); }
      catch (error) { draftError(error); }
      finally { S.saving = false; showState(); }
    });
    $('setup').addEventListener('click', () => send('setup')); $('cancel').addEventListener('click', () => send('cancel'));
    const openFile = async (what) => { $('open-feedback').textContent = 'Opening ' + what + ' on the Fiji computer…';
      try { const result = await api('open', { what }); $('open-feedback').textContent = 'Opened ' + result.opened + ' on the Fiji computer.'; }
      catch (error) { $('open-feedback').textContent = error.message; }
    };
    $('log').addEventListener('click', () => openFile('log'));
    $('open-folder').addEventListener('click', () => openFile('folder'));
    $('side').addEventListener('toggle', () => { if (S.state) showState(); }, true);   // a folded section that opens or closes changes where its messages can sit
    $('close-page').addEventListener('click', async () => {
      try { await api('close', {}); S.closed = true; play(false); $('lost').textContent = 'This analysis page is closed. Saved results remain in the results folder. You can close this tab.'; $('lost').classList.add('on'); showState(); }
      catch (error) { draftError(error); }
    });
    $('discard').addEventListener('click', () => { S.note = null; S.unlinked.clear(); S.protocolDrafts = {}; fillDraft((S.results && S.results.appliedDraft) || S.initialDraft); $('confirmed').checked = !!(S.results && S.results.appliedDraft); for (const key of FILES) $(key + '-file').value = ''; draftChanged(); });
    // a flagged cell is opened from its row: the Show list goes to the flagged cells, so the arrows step through them
    $('table-flagged').addEventListener('click', (event) => { const tr = event.target.closest('tr'); if (!tr || !tr.dataset.cell) return;
      const name = tr.dataset.cell; if (S.group !== 5) showGroup(5); chooseCell(name); showTab('review'); $('quality-review').open = false; });
    $('group').addEventListener('change', () => { if (Number($('group').value) !== 5 && S.qualityFlag) { S.qualityFlag = null; showQuality(); } });
    // The review mark of the selected cell: one button, and in its place a form. Enter saves; Esc cancels and keeps the cell.
    // A key held down repeats: only its first press counts, so a held Enter never saves and then removes, or opens and then saves.
    $('mark').addEventListener('click', () => { const action = markAction(S.cell); if (action && action.remove) removeMark(); else openMark(); });
    $('mark').addEventListener('keydown', (event) => { if ((event.key === 'Enter' || event.key === ' ') && event.repeat) event.preventDefault(); });
    $('mark-form').addEventListener('submit', (event) => { event.preventDefault(); saveMark(); });
    $('mark-form').addEventListener('keydown', (event) => {
      if (event.key === 'Enter' && event.repeat) event.preventDefault();
      else if (event.key === 'Escape') { event.preventDefault(); event.stopPropagation(); closeMark(); }
      else if (event.key === 'Enter' && event.target.tagName === 'SELECT') { event.preventDefault(); saveMark(); }
      else if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') event.stopPropagation();   // they move in the box, not to another cell
    });
    $('mark-reason').addEventListener('change', () => { S.markSaid = null; showMark(); });
    $('mark-note').addEventListener('input', () => { if (S.markSaid) { S.markSaid = null; showMark(); } });
    $('mark-cancel').addEventListener('click', closeMark);
    $('copy-methods').addEventListener('click', copyMethods);
    $('all-call-columns').addEventListener('change', showCalls); $('all-count-columns').addEventListener('change', showCounts); $('all-metric-columns').addEventListener('change', showMetricTable);
    $('response-window').addEventListener('change', () => { S.window = $('response-window').value; S.keep = null; updateWindow(); });
    $('metrics-table').addEventListener('change', showMetricTable);
    $('below-zero').addEventListener('change', () => { S.withoutBelowZero = $('below-zero').checked; showTiles(); showTally(); Plot.dirty = true; drawPlot(); });   // the plot and After review leave the same cells out
    $('group').addEventListener('change', () => { S.group = Number($('group').value); S.keep = null; listCells(); });
    // The Cell box. A full name, typed or picked from its list, opens that cell at once; a number alone opens it on Enter.
    // An empty box is no cell, and its ✕ empties it. Text that names no cell in the Show list is put back to the cell that is open.
    // Esc in the box undoes the typing and leaves the box (the browser's own Esc, which would empty a search box, is stopped).
    { const box = $('cell'), settle = (numberToo, leave) => { const hit = cellTyped(box.value, numberToo);
        if (hit !== undefined) { if (hit !== S.cell) chooseCell(hit); if (leave) box.blur(); } else if (leave) box.value = S.cell ? shortName(S.cell) : ''; };
      box.addEventListener('input', () => settle(false, false));
      box.addEventListener('change', () => settle(true, true));
      box.addEventListener('keydown', (event) => { if (event.key === 'Enter') { settle(true, true); event.preventDefault(); } else if (event.key === 'Escape') { event.preventDefault(); box.value = S.cell ? shortName(S.cell) : ''; box.blur(); } });
      box.addEventListener('focus', () => box.select()); }
    $('previous').addEventListener('click', () => stepCell(-1)); $('next').addEventListener('click', () => stepCell(1));
    $('every').addEventListener('change', () => { S.every = $('every').checked; Plot.dirty = true; drawPlot(); });
    $('thresholds').addEventListener('change', () => { S.thresholds = $('thresholds').checked; Plot.dirty = true; drawPlot(); });
    $('view').addEventListener('change', () => { play(false); S.view = $('view').value; S.viewChosen = true; updateTransport(); requestFrame(); });
    $('outlines').addEventListener('change', () => { S.outlines = $('outlines').checked; $('outline-colours').disabled = $('number-cells').disabled = !S.outlines; drawPicture(); });   // colours and numbers belong to outlines that are drawn
    $('zoom-cell').addEventListener('click', zoomToCell); $('zoom-all').addEventListener('click', () => setView(1));
    $('zoom-in').addEventListener('click', () => zoomStep(1)); $('zoom-out').addEventListener('click', () => zoomStep(-1)); $('zoom-region').addEventListener('click', () => setRegion(!Pic.region));
    $('play').addEventListener('click', () => play(!S.timer)); $('speed').addEventListener('change', () => { if (S.timer) play(true); });
    $('time').addEventListener('input', () => setFrame(Number($('time').value), true));
    $('frame-back').addEventListener('click', () => setFrame(S.frame - 1, true)); $('frame-forward').addEventListener('click', () => setFrame(S.frame + 1, true));
    // the heatmap's lists and tick boxes: how they stand on the page is how the heatmap opens, so the page has one place for these defaults
    for (const [id, key] of [['order', 'order'], ['palette', 'palette'], ['range', 'range'], ['row-names', 'rowNames'], ['hide-excluded', 'hideExcluded']]) {
      const box = $(id), read = () => { S[key] = box.type === 'checkbox' ? box.checked : box.value; }; read(); box.addEventListener('change', () => { read(); drawHeatmap(); });
    }
    $('save-plot').addEventListener('click', () => save(plotPicture, plotKind(), 'plot-saved')); $('save-heatmap').addEventListener('click', () => save(heatmapPicture, 'heatmap', 'heatmap-saved'));

    // The picture. Zooming is at the pointer: scrolling over the picture scrolls the page, and with Ctrl or Cmd held it zooms. A
    // pinch on a trackpad arrives as wheel events marked Ctrl with no Ctrl key held, in small steps that give the change of scale
    // as exp(−deltaY / 100): taken at that, the picture follows the fingers. Safari sends a pinch as gesture events with the scale itself.
    const picture = $('picture'), stage = $('picture-stage');
    window.addEventListener('keydown', (event) => { if (event.key === 'Control') Pic.ctrl = true; }, true); window.addEventListener('keyup', (event) => { if (event.key === 'Control') Pic.ctrl = false; }, true); window.addEventListener('blur', () => { Pic.ctrl = false; });
    picture.addEventListener('wheel', (event) => { if (!(event.ctrlKey || event.metaKey)) return; event.preventDefault(); const p = panelPoint(event);
      const pinch = event.ctrlKey && !event.metaKey && !Pic.ctrl && event.deltaMode === 0 && Math.abs(event.deltaY) < 50;
      zoomAt(Pic.zoom * (pinch ? Math.exp(-event.deltaY / 100) : Math.pow(1.0015, -event.deltaY)), p.x, p.y); }, { passive: false });
    picture.addEventListener('gesturestart', (event) => { event.preventDefault(); Pic.pinch = Pic.zoom; });
    picture.addEventListener('gesturechange', (event) => { event.preventDefault(); if (Pic.pinch) { const p = panelPoint(event); zoomAt(Pic.pinch * event.scale, p.x, p.y); } });
    picture.addEventListener('gestureend', (event) => { event.preventDefault(); Pic.pinch = null; });
    // A press that moves more than 3 pixels is a drag, never a click. Zoomed in, it moves the picture; with "Zoom to region" on,
    // or Shift held, it draws a rectangle and the view zooms to it. What was open when a click (or the first click of a double-click) began is kept.
    picture.addEventListener('mousedown', (event) => { if (event.button) return; const p = panelPoint(event);
      Pic.drag = { x: event.clientX, y: event.clientY, band: Pic.region || event.shiftKey ? [p.x, p.y] : null }; Pic.moved = false; if (event.detail < 2) Pic.held = { cell: S.cell, keep: S.keep, mode: S.plotMode }; });
    window.addEventListener('mousemove', (event) => { const drag = Pic.drag; if (!drag) return;
      if (!Pic.moved && Math.abs(event.clientX - drag.x) + Math.abs(event.clientY - drag.y) <= 3) return;
      Pic.moved = true; tip('');
      if (drag.band) { const p = panelPoint(event); Pic.band = [drag.band[0], drag.band[1], Math.max(0, Math.min(picture.clientWidth, p.x)), Math.max(0, Math.min(picture.clientHeight, p.y))]; drawPicture(); return; }
      const z = place(picture.clientWidth, picture.clientHeight).z, dx = event.clientX - drag.x, dy = event.clientY - drag.y; drag.x = event.clientX; drag.y = event.clientY;
      if (Pic.zoom > 1) setView(Pic.zoom, Pic.cx - dx / z, Pic.cy - dy / z); });
    window.addEventListener('mouseup', () => { const drag = Pic.drag; Pic.drag = null; Plot.scrubbing = false; if (!drag || !Pic.moved) return;
      const band = Pic.band; Pic.band = null; if (band && zoomToBand(band)) setRegion(false); else drawPicture(); });   // the button turns itself off once the view has zoomed
    // A click on a cell opens it; a click on empty space shows no cell. A drag is not a click.
    picture.addEventListener('click', (event) => { if (Pic.moved) return; const i = cellAt(event), name = i >= 0 ? S.info.names[i] : null; if (name || S.cell) chooseCell(name); });   // (the open cell clicked again: the zoomed view centres on it)
    // A double-click shows the whole picture. On empty space it keeps the cell that was open: its first click had cleared it.
    picture.addEventListener('dblclick', (event) => { const held = Pic.held;
      if (held && held.cell && !S.cell && cellAt(event) < 0) {
        // a cell a mark kept listed (S.keep) comes back to the list in its place, so Show stays on its group
        if (held.keep === held.cell && !S.listed.includes(held.cell)) { S.keep = held.cell; S.listed = S.info.names.filter((n) => n === held.cell || S.listed.includes(n)); offerCells(); }
        pick(held.cell); if (S.plotMode !== held.mode) setPlotMode(held.mode); }
      setView(1); });
    picture.addEventListener('mousemove', (event) => { Pic.pointer = panelPoint(event); if (Pic.drag) return; const i = cellAt(event); tip(i >= 0 ? fullName(S.info.names[i]) + (callOf(S.info.names[i]) ? ' · ' + callOf(S.info.names[i]) : '') : '', event.clientX, event.clientY); });
    picture.addEventListener('mouseleave', () => { Pic.pointer = null; tip(''); });
    // The first times the pointer comes onto the picture, a line beside the controls says how to zoom and move with the pointer alone.
    // It goes when the pointer leaves, or after 8 s. A pointer that only crosses the picture does not use up one of the times.
    { const hint = $('view-hint'); let shown = null, read = null; hint.textContent = 'Pinch or ' + (/Mac|iPhone|iPad/.test(navigator.platform) ? 'Cmd' : 'Ctrl') + ' + scroll to zoom. Drag to move.';
      stage.addEventListener('mouseenter', () => { Pic.over = true; if (!S.info || Pic.hinted >= HINT_TIMES) return; hint.hidden = false; clearTimeout(shown); clearTimeout(read);
        read = setTimeout(() => { Pic.hinted++; store('step4.viewhint', String(Pic.hinted)); }, HINT_READ); shown = setTimeout(() => { hint.hidden = true; }, 8000); });
      stage.addEventListener('mouseleave', () => { Pic.over = false; hint.hidden = true; clearTimeout(shown); clearTimeout(read); }); }
    // The small map: a press or a drag in it puts that point of the recording in the middle of the view.
    { const map = $('view-map'), move = (event) => { const box = map.getBoundingClientRect(); setView(Pic.zoom, (event.clientX - box.left) / box.width * S.info.width, (event.clientY - box.top) / box.height * S.info.height); };
      map.addEventListener('pointerdown', (event) => { if (event.button) return; event.preventDefault(); try { map.setPointerCapture(event.pointerId); } catch (unsupported) { } move(event); });
      map.addEventListener('pointermove', (event) => { if (event.buttons & 1) move(event); }); }

    const plot = $('plot');
    plot.addEventListener('pointermove', (event) => {
      if (Plot.scrubbing) {
        const g = Plot.geom, box = plot.getBoundingClientRect(); if (g) seek(g.x0 + Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W)) * (g.x1 - g.x0)); tip(''); return;
      }
      const at = traceAt(event), sample = at ? at.index : null, moved = sample !== Plot.at; Plot.at = sample;   // the numbers follow the pointer
      if (!at) { tip(''); if (Plot.hover || moved) { Plot.hover = null; drawPlot(); } return; }
      if (at.name !== Plot.hover || moved) { Plot.hover = at.name; drawPlot(); }
      tip(at.name ? fullName(at.name) + ' · ' + shownTraces()[at.name][at.index].toFixed(3) + ' at ' + S.results.time[at.index].toFixed(1) + ' s' + (callOf(at.name) ? ' · ' + callOf(at.name) : '') : '', event.clientX, event.clientY);
    });
    plot.addEventListener('pointerleave', () => { tip(''); if (Plot.hover || Plot.at != null) { Plot.hover = null; Plot.at = null; drawPlot(); } });
    plot.addEventListener('pointerdown', (event) => { const at = traceAt(event); if (!at) return; event.preventDefault(); if (at.name && at.name !== S.cell) chooseCell(at.name); Plot.scrubbing = true; plot.setPointerCapture(event.pointerId); seek(at.time); });
    for (const type of ['pointerup', 'pointercancel', 'lostpointercapture']) plot.addEventListener(type, (event) => { if (!Plot.scrubbing) return; Plot.scrubbing = false; const at = traceAt(event); Plot.at = at ? at.index : null; drawPlot(); });   // the drag is over: the numbers go back to the pointer
    plot.addEventListener('wheel', (event) => { const g = Plot.geom; if (!g || !(event.ctrlKey || event.metaKey)) return; event.preventDefault(); const box = plot.getBoundingClientRect(), share01 = Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W));
      const centre = g.x0 + share01 * (g.x1 - g.x0), factor = Math.pow(1.0015, event.deltaY), span = Math.min(g.full1 - g.full0, Math.max((g.full1 - g.full0) / 60, (g.x1 - g.x0) * factor));
      let a = centre - share01 * span, b = a + span; if (a < g.full0) { b += g.full0 - a; a = g.full0; } if (b > g.full1) { a -= b - g.full1; b = g.full1; }
      Plot.x0 = Math.max(g.full0, a); Plot.x1 = Math.min(g.full1, b); if (Plot.x1 - Plot.x0 >= g.full1 - g.full0 - 1e-9) Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); }, { passive: false });
    plot.addEventListener('dblclick', () => { Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); });

    const heatmap = $('heatmap');
    heatmap.addEventListener('click', (event) => { const at = heatmapAt(event); if (at) { chooseCell(at.name); if (at.inside) setFrame(at.frame + 1, true); showTab('review'); } });
    heatmap.addEventListener('mousemove', (event) => { const at = heatmapAt(event); if (!at) { tip(''); return; } const trace = S.results.traces[at.name], v = trace ? trace[at.frame] : null;
      tip(fullName(at.name) + (callOf(at.name) ? ' · ' + callOf(at.name) : '') + (belowZero(at.name) ? ' · below zero' : '') + (at.inside ? ' · ' + (at.gap ? 'acquisition gap; nearest frame' : finite(v) ? v.toFixed(3) : 'no value') + ' at ' + S.results.time[at.frame].toFixed(1) + ' s' : ''), event.clientX, event.clientY); });
    heatmap.addEventListener('mouseleave', () => tip(''));

    document.addEventListener('keydown', (event) => { const tag = (event.target.tagName || '').toLowerCase(), field = tag === 'input' || tag === 'select' || tag === 'textarea'; if ($('panel-review').hidden || Fig.open) return;
      // With the pointer over the picture, + and − zoom about it and 0 shows the whole picture: not while typing or choosing from a list, and not with Ctrl or Cmd, which is the browser's own zoom.
      if (Pic.over && Object.prototype.hasOwnProperty.call(ZOOM_KEYS, event.key) && !(field && !['checkbox', 'radio', 'range'].includes(event.target.type)) && !event.ctrlKey && !event.metaKey && !event.altKey) {
        const by = ZOOM_KEYS[event.key], p = Pic.pointer || {}; if (by) zoomStep(by, p.x, p.y); else setView(1); event.preventDefault(); return; }
      if (field) return;
      if (event.key === 'Escape' && (Pic.region || Pic.band) && !event.defaultPrevented) { Pic.drag = null; setRegion(false); event.preventDefault(); }   // Esc leaves "Zoom to region" and keeps the cell
      else if (event.key === 'Escape' && S.cell && !event.defaultPrevented) { chooseCell(null); event.preventDefault(); }
      else if (event.key === 'ArrowLeft') { stepCell(-1); event.preventDefault(); } else if (event.key === 'ArrowRight') { stepCell(1); event.preventDefault(); } else if (event.key === ' ' && S.info && S.info.hasMovie && S.view !== 'picture' && tag !== 'button') { play(!S.timer); event.preventDefault(); } });
    window.addEventListener('resize', redraw);
    if (window.ResizeObserver) { new ResizeObserver(redraw).observe($('panel-review')); new ResizeObserver(() => drawHeatmap()).observe($('heatmap-box'));
      // the plot also changes size when the cell card above it does (a refusal, the mark form, the mark button going while a run
      // starts). Drawing sets only the canvas's pixels, not its size on the page, so this never calls itself again.
      new ResizeObserver(() => { Plot.dirty = true; drawPlot(); }).observe($('plot')); }
  }
  function redraw() { Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); }

  wire();
  start().catch((problem) => { $('lost').classList.add('on'); $('status').textContent = problem.message; });
  window.calciumStep4 = { S, Plot, Pic, Heat, Buf, Clock, Fig, STYLE, shownCall, startingReason, markProblem, refreshReview, loadResults, callsRows, outlineLook, place, tableColumns, shownCell, setView, viewArea, savedArea, viewCanvas, pairCanvas, scaleBarLength, pictureKey, cellNumber, plotCells, figureCanvas, figureKept, figureRestore, openFigure, closeFigure, rateText, samplesFor, cellTyped, cellAt, resetFrames: () => resetFrames(Buf.version), pick, setFrame, showGroup, showTab, collectDraft, frameTimes, nearestFrame, samplingSupport, referenceLevels, legendLayout, exportMetadata, heatmapAt, heatmapData, pictureRecord, decimal, columnWords, metricSummary, showCalls, showCounts, showMetricTable, plotNumbers, plotKind, COLUMN_WORDS, SUMMARY_ROWS, withCodes, metricRows, Around, aroundNow, aroundParts, zoneInk, savedRecord, plotPicture };   // for the browser test
})();

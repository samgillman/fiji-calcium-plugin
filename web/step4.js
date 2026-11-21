// Step 4 page. Everything shown here is read from the tables the analysis saved, or from the recording's own frames.
// The page draws; it does not calculate results. Settings go back to the plugin, which runs the same engine as the Fiji window.
'use strict';
(function () {
  const $ = (id) => document.getElementById(id);
  const CALL = { responded: '#177465', not: '#C8660E', other: '#8A94A0' };
  const ON_PICTURE = { responded: '#3ddbb5', not: '#ffa24a', other: '#aeb8c2' };
  const GROUPS = ['All cells', 'Responded', 'Did not respond', 'Excluded', 'Unscorable', 'Needs quality review'];
  /** The colour overview and the background areas are still pictures; every other view is a frame of the recording. */
  const still = () => S.view === 'picture' || S.view === 'background';
  const kind = (call) => call === 'Responded' ? 'responded' : call === 'Did not respond' ? 'not' : 'other';
  const inGroup = (call, group, name) => group === 0 || (group === 1 && call === 'Responded') || (group === 2 && call === 'Did not respond') ||
    (group === 3 && call === 'Excluded') || (group === 4 && call === 'Unscorable') || (group === 5 && qualityFor(name).length > 0 && (!S.qualityFlag || qualityFor(name).includes(S.qualityFlag)));
  const S = { unlinked: new Set(), validationShown: null, exportScale: null, imported: null, info: null, state: null, results: null, version: -1, group: 0, cell: null, listed: [], frame: 1, view: 'raw', viewChosen: false,
    every: false, plotMode: 'mean', signal: 'dff', outlines: true, outlineWidth: 1, order: 'call', palette: 'plasma', range: 'robust', tall: false, timer: null, lost: 0, note: null, polling: false, paths: [], boxes: [], draft: null, draftAck: '', revision: 0, draftTimer: null, draftQueue: Promise.resolve(), window: null, saving: false, initialDraft: null, protocolDrafts: {}, quality: {}, metrics: {}, closed: false };
  const Plot = { dirty: true, layer: document.createElement('canvas'), x0: null, x1: null, y: {}, hover: null, geom: null, scrubbing: false, drawn: [], series: [] };
  const Pic = { base: null, shown: null, zoom: 1, panX: 0, panY: 0, drag: null, moved: false, shownFrame: null, error: null };
  // Movie frames: compressed images are kept for the session (about 40 MB for 300 full-size frames, much less at display size),
  // decoded ones only near the playhead. Keys carry the server's frame version, which changes with the baseline behind
  // "Change from baseline" and the display range, so a stale image is never reused.
  const Buf = { version: null, blobs: new Map(), bitmaps: new Map(), inflight: new Map(), decoding: new Map(), failed: new Map(), bytes: 0 };
  const PARALLEL = 4, DECODED_AROUND = 30, MAX_BYTES = 90e6, MAX_FPS = 30;
  Object.defineProperty(Pic, 'loading', { get: () => Buf.inflight.size > 0 || Buf.decoding.size > 0 });
  // The playback clock runs in recording seconds; frames are shown when the playhead reaches their timestamps.
  const Clock = { t: 0, last: 0, step: 1, times: [], real: false, waitingSince: 0 };
  const Heat = { geom: null, names: [] };

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
  const callOf = (name) => callRow(name).Status;
  const number = (text) => { const v = parseFloat(text); return isFinite(v) ? v : NaN; };

  /** Sets a canvas to the size it is shown at (sharp on high-density screens) and returns a context in CSS pixels. */
  function fit(canvas) {
    const dpr = (S && S.exportScale) || window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight;   // exports draw at a fixed minimum resolution
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
  const tickText = (v, step) => v.toFixed(step >= 1 ? 0 : Math.min(6, Math.ceil(-Math.log10(step) - 1e-9)));
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
      add(details, 'Min_baseline_samples', 'Baseline samples').previousSibling.title = 'The fewest baseline samples a cell needs to be scored'; add(details, 'Max_gap_s', 'Maximum gap (s)').placeholder = 'auto';
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
  /** Under the frame slider: grey for the noise baseline, a mark for the stimulus, orange for the scored frames of the primary stimulus. */
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
    const overrides = $('field-overrides').value.trim();
    if (legacy) p.field_overrides = null;
    else { try { p.field_overrides = overrides ? JSON.parse(overrides) : null; } catch (_) { p.field_overrides = overrides; } }
    return draft;
  }
  function fillDraft(draft) {
    S.draft = clone(draft);
    exactInput($('sampling_hz'), draft.values.sampling_hz);   // the F₀ frame boxes are made with the primary stimulus's rows
    for (const key of ['condition', 'biological_unit']) $(key).value = draft.values[key] || '';
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
    if ((S.results.traces[name] || []).some((v) => finite(v) && v < -1)) flags.add('TRACE_BELOW_MINUS_ONE');
    return Array.from(flags);
  }
  // Set on every cell by the engine whenever the window has the matching Window_warning: a property of the protocol and frame
  // rate, shown once as a notice for the recording instead of flagging every cell.
  const RECORDING_FLAGS = new Set(['DURATION_BELOW_SAMPLE_INTERVAL']);
  const FLAG_TEXT = {
    BASELINE_SD_AT_LEAST_ONE_THIRD_LOCAL_F0: 'Baseline noise is at least a third of the local F₀, so the fractional rise threshold is dominated by noise',
    DURATION_BELOW_SAMPLE_INTERVAL: 'Minimum response duration is shorter than one frame interval; timing is coarsely resolved',
    reference_below_3_baseline_sd_review: 'Low F₀ relative to baseline noise; review fractional amplitudes and the response call',
    TRACE_BELOW_MINUS_ONE: 'ΔF/F₀ falls below −1; corrected fluorescence becomes negative relative to its positive reference. Review background and baseline stability',
    RISE_START_BELOW_BASELINE_SD: 'The legacy 10% rise level is below baseline SD; the rise-time estimate is sensitive to baseline noise'
  };
  const plural = (n, word) => n + ' ' + word + (n === 1 ? '' : 's');
  /** How many cells of the selected window got a call (responded or did not). Zero means the run finished without a single call. */
  const called = () => { const c = S.results ? ((S.results.countsByWindow || {})[S.window] || S.results.counts || {}) : {}; const v = number(c.Evaluable); return isFinite(v) ? v : 0; };
  /** Why cells are unscorable, in the engine's words: the reason most of them share. */
  function unscorableReason() {
    const seen = new Map(); for (const name of S.info.names) { const row = callRow(name); if (row.Status === 'Unscorable' && row.Reason) seen.set(row.Reason, (seen.get(row.Reason) || 0) + 1); }
    let best = ''; for (const [reason, n] of seen) if (!best || n > seen.get(best)) best = reason; return best;
  }
  const flagText = (flag) => FLAG_TEXT[flag] || flag.replace(/_/g, ' ');
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
    if (S.qualityFlag && !counts.has(S.qualityFlag)) S.qualityFlag = null;
    $('quality-review').hidden = !affected && !notices.length;
    summary.textContent = affected ? plural(affected, 'cell') + ' flagged' : 'Timing note';
    $('quality-lead').textContent = affected ? plural(affected, 'cell') + ' flagged' + (called() ? ', including ' + plural(responders, 'responder') : '') + '.' : 'Review the applied timing.';
    list.replaceChildren();
    for (const [flag, count] of counts) {
      // each flag is a filter of its own, so a small group (for example the cells whose ΔF/F₀ falls below −1) can be isolated
      const li = document.createElement('li'), link = document.createElement('button'); link.type = 'button'; link.className = 'flag-filter' + (S.qualityFlag === flag ? ' on' : '');
      link.textContent = plural(count.total, 'cell') + ' (' + plural(count.responders, 'responder') + '): ' + flagText(flag) + '.'; link.title = S.qualityFlag === flag ? 'Show all flagged cells again' : 'Show only these cells';
      link.addEventListener('click', () => { S.qualityFlag = S.qualityFlag === flag ? null : flag; S.group = 5; $('group').value = 5; showQuality(); listCells(); showTab('review'); });
      li.appendChild(link); list.appendChild(li);
    }
    for (const message of new Set(notices)) { const li = document.createElement('li'); li.textContent = message; list.appendChild(li); }
    $('quality-filter').disabled = !affected;
  }
  /** What each saved metric column is: [column, name, unit, where it is measured or counted from, response].
   *  response: the column describes the timing or shape of a response, so it is summarised only over cells that responded
   *  (for a cell that did not respond it is the timing of noise). The order here is the order in the summary. A numeric
   *  column that is not listed is still shown, under its saved name. */
  const METRICS = [
    ['Peak_dFF0', 'Peak', 'ΔF/F₀', 'whole recording'],
    ['Response_Amplitude', 'Peak above baseline', 'ΔF/F₀', 'whole recording'],
    ['Peak_increase', 'Peak rise', 'fraction of baseline', 'response window'],
    ['Window_peak_local_fraction', 'Peak rise', 'fraction of baseline', 'response window'],
    ['SNR', 'Signal-to-noise', 'peak ÷ baseline SD', 'whole recording'],
    ['AUC', 'Area under the curve', 'ΔF/F₀ × s', 'whole recording, signed'],
    ['Window_AUC_fraction_s', 'Area under the curve', 'fraction × s', 'response window, signed'],
    ['Baseline_SD', 'Baseline SD', 'ΔF/F₀', 'F₀ frames'],
    ['F0', 'F₀', 'gray levels', 'F₀ frames'],
    ['Baseline_SD_raw', 'Baseline SD', 'gray levels', 'F₀ frames'],
    ['Onset_s', 'Onset time', 's', 'from recording start', true],
    ['Run_start_s', 'Onset time', 's', 'from recording start', true],
    ['Response_latency_s', 'Latency', 's', 'from the stimulus', true],
    ['Confirmation_s', 'Response confirmed', 's', 'from recording start', true],
    ['Run_end_s', 'Response end', 's', 'from recording start', true],
    ['Run_duration_s', 'Response duration', 's', 'above the threshold', true],
    ['Time_to_Peak', 'Peak time', 's', 'from recording start', true],
    ['Peak_latency_s', 'Peak latency', 's', 'from the stimulus', true],
    ['Time_to_25_Peak', '25 % of peak reached', 's', 'from recording start', true],
    ['Time_to_50_Peak', '50 % of peak reached', 's', 'from recording start', true],
    ['Time_to_75_Peak', '75 % of peak reached', 's', 'from recording start', true],
    ['Rise_Time', 'Rise time', 's', '10 % to 90 % of peak', true],
    ['Rise_Rate_10_90_dFF0_per_s', 'Rise rate', 'ΔF/F₀ per s', '10 % to 90 % of peak', true],
    ['FWHM', 'FWHM', 's', 'width at half of peak', true],
    ['FWHM_Lower_Bound', 'FWHM, lower bound', 's', 'still above half at the end', true],
    ['Half_Width', 'Half-width', 's', 'FWHM ÷ 2', true],
  ];
  const METRIC = new Map(METRICS.map((m, order) => [m[0], { name: m[1], unit: m[2], measured: m[3], response: !!m[4], order }]));
  /** Columns that say how the analysis was set up or how much of a window had data, not what a cell did: left out of the summary. */
  const NOT_METRICS = new Set(['Reference_value', 'Baseline_SD_input', 'Baseline_samples', 'Normalization_minimum_reference', 'Baseline_start_s', 'Baseline_end_s', 'Stimulus_onset_s', 'Start_s', 'End_s', 'Window_AUC_start_s', 'Window_AUC_end_s', 'Window_AUC_covered_s', 'Window_AUC_coverage_fraction', 'Run_samples', 'Eligible']);
  /** Mean, SD, SEM and number of cells of every metric: for all cells, by response call, and by how the cells were chosen
   *  when that was recorded. Worked out here from the saved per-cell tables, which are not changed.
   *  - Each metric is one row. Saved columns that hold the same number for every cell in this analysis (which ones do
   *    depends on the settings) share a row, and "Saved as" names them all.
   *  - The timing and shape of a response are summarised over the cells that responded. */
  function metricSummary() {
    const r = S.results, tags = r.chosenBlind || null, primary = (r.windows || []).find((w) => w.role === 'Stimulus') || (r.windows || [])[0];
    const called = (primary && r.callsByWindow && r.callsByWindow[primary.key]) || r.calls || {}, status = (n) => (called[n] || {}).Status;
    const groups = [['All cells', () => true], ['Responded', (n) => status(n) === 'Responded'], ['Did not respond', (n) => status(n) === 'Did not respond']];
    if (tags) groups.push(['Chosen blind', (n) => tags[n] === true], ['Added after the response', (n) => tags[n] === false]);
    // every numeric column of the two per-cell tables, as cell -> number
    const columns = [];
    for (const key of ['metrics', 'windowMetrics']) {
      const data = r.tables[key]; if (!data || !data.rows.length) continue;
      const cell = data.header.indexOf('Cell'), role = data.header.indexOf('Role'); if (cell < 0) continue;
      const rows = data.rows.filter((row) => role < 0 || row[role] === 'Stimulus');
      data.header.forEach((name, c) => {
        if (c === cell || NOT_METRICS.has(name)) return;
        const pairs = rows.map((row) => [row[cell], row[c]]).filter((pair) => pair[1] !== '' && pair[1] != null && pair[1] !== 'NA');
        if (!pairs.length || !pairs.every((pair) => isFinite(Number(pair[1])))) return;
        const known = METRIC.get(name) || { name, unit: '', measured: key === 'windowMetrics' ? 'response window' : '', response: false, order: METRICS.length + columns.length };
        columns.push({ saved: [name], values: new Map(pairs.map((pair) => [pair[0], Number(pair[1])])), ...known });
      });
    }
    columns.sort((a, b) => a.order - b.order);
    // columns with the same number for every cell become one metric
    const same = (a, b) => { if (a.values.size !== b.values.size) return false; for (const [cell, x] of a.values) { const y = b.values.get(cell); if (y === undefined || Math.abs(x - y) > 1e-9 * Math.max(1, Math.abs(x), Math.abs(y))) return false; } return true; };
    const metrics = [];
    for (const column of columns) { const twin = metrics.find((m) => same(m, column)); if (twin) twin.saved.push(column.saved[0]); else metrics.push(column); }
    const long = [['Metric', 'Unit', 'Measured', 'Saved_as', 'Group', 'Mean', 'SD', 'SEM', 'N_cells']], wide = [];
    for (const m of metrics) {
      const line = [m.name, m.unit, m.measured];
      for (const [group, belongs] of groups) {
        const xs = []; for (const [cell, x] of m.values) if (belongs(cell) && (!m.response || status(cell) === 'Responded')) xs.push(x);
        const k = xs.length, mean = k ? xs.reduce((a, b) => a + b, 0) / k : NaN, sd = k > 1 ? Math.sqrt(xs.reduce((a, b) => a + (b - mean) * (b - mean), 0) / (k - 1)) : NaN, sem = k > 1 ? sd / Math.sqrt(k) : NaN;
        long.push([m.name, m.unit, m.measured, m.saved.join('; '), group, k ? mean : '', k > 1 ? sd : '', k > 1 ? sem : '', k]);
        line.push(k ? (+mean.toPrecision(4)) + (k > 1 ? ' ± ' + (+sem.toPrecision(3)) : '') + '  (n = ' + k + ')' : '–');
      }
      line.push(m.saved.join(', ')); wide.push(line);
    }
    return { wide: { header: ['Metric', 'Unit', 'Measured'].concat(groups.map((g) => g[0]), ['Saved as']), rows: wide, computed: true }, long };
  }
  const csvText = (rows) => rows.map((row) => row.map((v) => { const s = v == null ? '' : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }).join(',')).join('\n') + '\n';
  /** The table on the Metrics tab as CSV rows: the summary in long form (one row per metric and group), any other table as saved. */
  function metricRows() { const key = $('metrics-table').value; if (key === 'summary') return metricSummary().long; const data = S.results.tables[key] || { header: [], rows: [] }; return [data.header].concat(data.rows); }
  function showMetricTable() {
    const key = $('metrics-table').value;
    if (key === 'summary') { $('metric-help').textContent = 'Mean ± SEM (n). Onset, latency, peak time, rise and width are given for cells that responded. Saved columns with the same numbers share a row. The saved file also has the SD.'; table('table-metrics', metricSummary().wide); return; }
    $('metric-help').textContent = key === 'metrics' ? 'Legacy metrics: AUC is the signed net integral over the whole recording; negative drift can make it negative. Peak time is a recording timestamp. Rise time is the 10–90% crossing interval after the display baseline and can be noise-sensitive. For stimulus-specific values, choose Response-window metrics (including Peak_latency_s and Window_AUC_fraction_s).' : key === 'windowMetrics' ? 'Response-window metrics use the baseline before the stimulus. Peak_latency_s is measured from the stimulus (the first frame after it), to the nearest frame; Window_AUC_fraction_s covers the response window.' : '';
    table('table-metrics', S.results.tables[key]);
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
    for (const key of ['format', 'size', 'compression']) if (prefs['step4.' + key]) Save[key] = prefs['step4.' + key]; if (prefs['step4.dpi']) Save.dpi = Number(prefs['step4.dpi']);
    if (prefs['step4.outline']) { S.outlineWidth = Number(prefs['step4.outline']) || 1; $('outline-width').value = prefs['step4.outline']; }
    if (prefs['step4.every']) { try { window.localStorage.setItem('step4.every', prefs['step4.every']); } catch (blocked) { } }
    if (prefs['step4.help']) { const on = prefs['step4.help'] === '1'; document.body.classList.toggle('show-help', on); $('toggle-help').textContent = on ? 'Hide explanations' : 'Show explanations'; }
    $('recording').textContent = S.info.recording; $('recording').title = S.info.recording;
    document.title = S.info.names.length + ' cells · SimpleCa²⁺ · ' + S.info.recording + ' · ' + S.info.key.slice(0, 6); $('version').textContent = 'Calcium imaging plugin ' + S.info.pluginVersion;
    S.initialDraft = clone(S.info.draft);
    S.revision = Math.max(0, Number(S.info.draftRevision) || 0); fillDraft(S.info.draft); S.draftAck = JSON.stringify(collectDraft());
    $('time').max = Math.max(1, S.info.frames);
    for (const option of Array.from($('view').options)) if ((option.value === 'picture' && !S.info.hasPicture) || (option.value === 'background' && !S.info.hasBackground) || (!['picture', 'background'].includes(option.value) && !S.info.hasMovie)) option.remove();
    if (S.info.hasBackground) { const areas = new Image(); areas.onload = () => { Pic.background = areas; if (S.view === 'background') requestFrame(); }; areas.src = 'api/background'; }
    S.view = S.info.hasMovie ? 'raw' : 'picture'; $('view').value = S.view;
    S.paths = S.info.rois.map((points) => { const p = new Path2D(); points.forEach((q, i) => i ? p.lineTo(q[0], q[1]) : p.moveTo(q[0], q[1])); p.closePath(); return p; });
    S.bounds = S.info.rois.map((points) => { let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const q of points) { x0 = Math.min(x0, q[0]); y0 = Math.min(y0, q[1]); x1 = Math.max(x1, q[0]); y1 = Math.max(y1, q[1]); } return { x0, y0, x1, y1 }; });
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
      if (state.version > 0 && state.version !== S.version) { await loadResults(); S.version = state.version; }
      showState(); updateTransport();
    } catch (problem) { if (++S.lost >= 3) $('lost').classList.add('on'); }
    finally { S.polling = false; }
  }
  function showState() {
    const state = S.state; if (!state) return;
    const engine = $('engine'); engine.textContent = state.engineReady ? 'Analysis engine ready' : 'First run: set up analysis engine'; engine.className = 'pill ' + (state.engineReady ? 'good' : 'warn'); engine.hidden = !!state.engineReady;
    if (S.reapplyFrames && editedTimes() && JSON.stringify(collectDraft()) === S.draftAck) { S.reapplyFrames = false; const before = JSON.stringify(collectDraft()); applyAllFrames(); if (JSON.stringify(collectDraft()) !== before) { draftChanged(); return; } }
    const validation = validateDraft(collectDraft());
    $('analyze').disabled = state.busy || S.saving || !state.engineReady || validation.errors.length > 0 || JSON.stringify(collectDraft()) !== S.draftAck || S.closed; $('cancel').disabled = !state.busy; $('log').disabled = !state.logAvailable;
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
    else if (S.results && state.tone === 'good' && state.folder) { status.className = 'status good'; status.textContent = 'Done.'; }
    else { status.className = 'status ' + (state.tone === 'good' ? 'good' : state.tone === 'bad' ? 'bad' : state.tone === 'warn' ? 'warn' : ''); status.textContent = state.message; }
    status.title = state.message || '';
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
  async function send(path, body) {
    S.note = null;
    try { acceptState(await api(path, body || {})); showState(); }
    catch (problem) { S.note = problem.message; showState(); }
  }
  // ---------------------------------------------------------------- saved results; selecting a window never recalculates a call
  async function loadResults() {
    S.results = await api('results'); const r = S.results;
    const records = (data) => Object.fromEntries(((data && data.rows) || []).map((row) => [row[data.header.indexOf('Cell')], Object.fromEntries(data.header.map((key, i) => [key, row[i]]))]));
    S.quality = records(r.tables.quality); S.metrics = records(r.tables.metrics);
    const windows = r.windows && r.windows.length ? r.windows : [{ key: 'Stimulus::' + (r.counts.Stimulus || 'Stimulus'), role: 'Stimulus', stimulus: r.counts.Stimulus || 'Stimulus' }];
    const selector = $('response-window'); selector.innerHTML = '';
    for (const w of windows) { const o = document.createElement('option'); o.value = w.key; o.textContent = w.stimulus + (w.role === 'Control' ? ' · control' : ' · primary'); selector.appendChild(o); }
    if (!windows.some((w) => w.key === S.window)) S.window = windows[0].key; selector.value = S.window;
    document.body.classList.add('has-results');
    for (const id of ['result-tiles', 'plot-card', 'group-label', 'window-select-label', 'every-label', 'call-legend']) $(id).hidden = false;
    $('window-select-label').hidden = !(S.results.windows && S.results.windows.length > 1);
    for (const tab of document.querySelectorAll('.tab')) tab.disabled = false;
    const tables = $('metrics-table'); tables.innerHTML = '';
    const titles = { metrics: 'Legacy whole-recording metrics', quality: 'Cell quality', registry: 'Registered cells', windowMetrics: 'Response-window metrics', replicates: 'Biological-replicate counts', groups: 'Condition counts', unassigned: 'Unassigned biological units', cohort: 'Cohort flow', recordingMetadata: 'Recording details', sensitivityCalls: 'Sensitivity: cell calls', sensitivityCounts: 'Sensitivity: responder counts', negativeControlCalls: 'Negative controls: cell calls', negativeControlCounts: 'Negative controls: responder counts' };
    { const o = document.createElement('option'); o.value = 'summary'; o.textContent = 'Summary: mean ± SEM of each metric'; tables.appendChild(o); }
    for (const key of Object.keys(r.tables || {})) if (!['calls', 'counts'].includes(key) && r.tables[key].header.length) { const o = document.createElement('option'); o.value = key; o.textContent = titles[key] || key.replace(/_/g, ' '); tables.appendChild(o); }
    tables.value = 'summary'; showMetricTable();
    S.withoutBelowZero = false;   // every cell is counted until the switch says otherwise
    table('table-calls', r.tables.calls); table('table-counts', r.tables.counts); showTally();
    for (const id of ['save-plot', 'save-heatmap', 'save-plot-folder', 'save-heatmap-folder', 'data-save', 'data-download', 'metrics-save', 'metrics-download']) $(id).disabled = false;
    for (const id of ['plot-saved', 'heatmap-saved', 'data-saved', 'metrics-saved']) $(id).textContent = '';
    if (!$('panel-data').hidden) showData();
    for (const option of $('signal').options) option.disabled = option.value !== 'dff' && !hasTable(option.value);
    if ($('signal').selectedOptions[0].disabled) { S.signal = Data.signal = 'dff'; $('signal').value = 'dff'; $('data-signal').value = 'dff'; }
    Plot.x0 = Plot.x1 = null; updateWindow(); setFrame(S.frame); showState();
  }
  /** % responding to the primary stimulus as worked out in Fiji (results.tally), by row key: all, blind, later, quiet.
   *  With the switch on, the cells whose corrected fluorescence went below zero are left out of every row. */
  function tallyRows() {
    const t = S.results && S.results.tally; if (!t || !t.rows) return null;
    const rows = {};
    for (const row of t.rows) { const c = S.withoutBelowZero ? row.without : row.counted; rows[row.key] = { label: row.label, responded: c[0], evaluable: c[1], share: c[1] ? Math.round(100 * c[0] / c[1]) + '%' : '–' }; }
    return rows;
  }
  /** The same tally as a table on the Responder counts tab: every row with all its cells and without those below zero. */
  function showTally() {
    const t = S.results && S.results.tally; $('tally-box').hidden = !t; if (!t) return;
    const percent = (c) => c[1] ? +(100 * c[0] / c[1]).toFixed(1) : '';
    // six columns, so the table fits its tab; the saved file has every number in a column of its own
    table('table-tally', { computed: true, header: ['Cells', '% responding', 'Responded', 'Evaluable', 'Cells below zero', 'Without them'],
      rows: t.rows.map((row) => [row.label, percent(row.counted), row.counted[0], row.counted[1], row.counted[2] - row.without[2], row.without[1] ? percent(row.without) + ' % (' + row.without[0] + ' of ' + row.without[1] + ')' : '–']) });
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
    const tally = chosen.role !== 'Control' ? tallyRows() : null, below = tally ? r.tally.belowZero.length : 0, said = (row) => row.share + ' (' + row.responded + ' of ' + row.evaluable + ')';
    $('below-zero-label').hidden = !below; $('below-zero').checked = !!S.withoutBelowZero;
    $('below-zero-text').textContent = 'Leave out ' + below + (below === 1 ? ' cell' : ' cells') + ' below zero';
    const reason = unscorableReason();
    if (!called()) { $('tile-responding').innerHTML = '<b>–</b><span>' + label + ': no cell got a call' + (reason ? ' (' + escape(reason.replace(/\.$/, '').toLowerCase()) + ')' : '') + '</span>'; }
    else if (tally && tally.all) {
      const lines = [];
      if (tally.blind && tally.later && tally.blind.evaluable && tally.later.evaluable) lines.push(['Chosen blind: ' + said(tally.blind) + ' · added after the response was shown: ' + said(tally.later), 'Cells chosen before the response was seen, and cells added once it was shown']);
      if (tally.quiet && tally.quiet.evaluable) lines.push([tally.quiet.label + ': ' + said(tally.quiet), 'Cells that need a rise of ' + Math.round(100 * r.tally.quietRise) + ' % or less to respond. Noisier cells are set aside whether they responded or not']);
      $('tile-responding').innerHTML = '<b>' + tally.all.share + '</b><span>' + label + ': ' + tally.all.responded + ' of ' + tally.all.evaluable + ' evaluable cells' + (S.withoutBelowZero && below ? ' · ' + below + ' below zero left out' : '') + '</span>'
        + lines.map(([text, title]) => '<span class="tile-line" title="' + escape(title) + '">' + escape(text) + '</span>').join('');
    }
    else $('tile-responding').innerHTML = '<b>' + (isFinite(percent) ? Math.round(percent) + '%' : '–') + '</b><span>' + all + '</span>';
    $('tile-registered').firstChild.textContent = whole('Registered'); $('tile-excluded').firstChild.textContent = whole('Excluded'); $('tile-unscorable').firstChild.textContent = whole('Unscorable');
    for (const [id, key] of [['tile-excluded', 'Excluded'], ['tile-unscorable', 'Unscorable']]) $(id).classList.toggle('zero', whole(key) === '0');
  }
  function updateWindow() {
    showTiles();
    const group = $('group'); group.innerHTML = '';
    GROUPS.forEach((name, i) => { const o = document.createElement('option'); o.value = i; o.textContent = name + ' (' + S.info.names.filter((n) => inGroup(callOf(n), i, n)).length + ')'; group.appendChild(o); });
    group.value = S.group; showQuality(); Plot.dirty = true; listCells(); drawHeatmap();
  }
  /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. */
  function listCells() {
    const names = S.info.names;
    S.listed = S.results ? names.filter((n) => inGroup(callOf(n), S.group, n)) : names.slice();
    const select = $('cell'); select.innerHTML = '';
    const none = document.createElement('option'); none.value = ''; none.textContent = 'None'; select.appendChild(none);
    for (const name of S.listed) { const o = document.createElement('option'); o.value = name; o.textContent = shortName(name); o.title = name; select.appendChild(o); }
    $('previous').disabled = $('next').disabled = S.listed.length < 1;
    pick(S.listed.includes(S.cell) ? S.cell : null);   // the selection never changes by itself: a cell that is not listed is no longer selected
  }
  function showGroup(group) { if (S.group !== group) { S.group = group; $('group').value = group; listCells(); } }
  /** Selects a cell everywhere: the list, the picture, the plot, the heatmap and the tables. */
  function pick(name) {
    if (name && S.results && !S.listed.includes(name)) { S.group = 0; $('group').value = 0; S.cell = name; listCells(); return; }
    S.cell = name || null; $('cell').value = name || '';
    const verdict = $('verdict'), row = name && S.results ? calls()[name] : null;
    if (!name) verdict.textContent = !S.results ? '' : S.listed.length ? 'No cell selected.' : 'No cells in this group.';
    else if (!row) verdict.innerHTML = S.results ? '<b>' + escape(shortName(name)) + '</b> · no response call was saved for this cell' : '';
    else {
      const parts = ['<span class="dot" style="background:' + CALL[kind(row.Status)] + '"></span><b>' + escape(shortName(name)) + ': ' + escape(row.Status || '') + '</b>', escape((row.Reason || '').toLowerCase())];
      const rise = number(row.Peak_increase), at = number(row.Peak_s), need = number(row.Threshold_increase);
      if (isFinite(rise)) parts.push('largest rise ' + (100 * rise).toFixed(1) + '%' + (isFinite(at) ? ' at ' + at.toFixed(1) + ' s' : '') + (isFinite(need) ? ' (needed ' + (100 * need).toFixed(1) + '%)' : ''));
      verdict.innerHTML = parts.filter(Boolean).join(' · ');
      verdict.title = 'Rise and threshold are shown as percentages of local prestimulus fluorescence. The saved call table stores the corresponding fractions.';
    }
    const flags = qualityFor(name); $('cell-quality').hidden = !flags.length; $('cell-quality').textContent = flags.map(flagText).join(' · '); $('cell-quality').title = $('cell-quality').textContent;
    $('zoom-cell').disabled = !name || !(S.info.hasPicture || S.info.hasMovie);
    for (const row of document.querySelectorAll('tbody tr')) row.classList.toggle('on', row.dataset.cell === name);
    if (Pic.zoom > 1 && name) centreOn(S.info.names.indexOf(name));
    Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); legend();
  }
  const escape = (text) => String(text).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  /** The plot opens on the average alone and stays on what was chosen: picking a cell never changes what the plot shows.
   *  A cell's own trace appears only under "One cell + average". */
  function setPlotMode(mode) {
    S.plotMode = mode;
    for (const button of document.querySelectorAll('[data-plot]')) { const on = button.dataset.plot === mode; button.classList.toggle('on', on); button.setAttribute('aria-pressed', on ? 'true' : 'false'); }
    Plot.dirty = true; drawPlot(); legend();
  }
  /** A cell chosen by the person (list, arrows, a click on the picture, a trace, a table row or the heatmap), or none. The plot
   *  follows: with a cell it shows that cell over the average, with none the average alone. "Responders vs non-responders" stays. */
  function chooseCell(name) {
    pick(name || null);
    if (S.results && S.cell && S.plotMode === 'mean') setPlotMode('cell'); else if (S.results && !S.cell && S.plotMode === 'cell') setPlotMode('mean');
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
    S.frame = Math.max(1, Math.min(Math.max(1, S.info.frames), Math.round(wanted))); $('time').value = S.frame;
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
      if (!response.ok) throw new Error('The recording frame could not be loaded.');
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
  function viewKey() {
    const box = $('view-key'), shown = viewColours().map((c) => c.join('|')).join(';'); if (!box || box.dataset.shown === shown) return;
    box.dataset.shown = shown; box.replaceChildren();
    for (const [colour, words] of viewColours()) { const chip = document.createElement('span'); chip.className = 'chip'; const dot = document.createElement('i'); dot.className = 'fill'; dot.style.background = colour; chip.append(dot, words); box.appendChild(chip); }
  }
  function frameNote() {
    if (still() || !S.info.hasMovie || !Pic.shownFrame || Pic.error) return;
    const f = Pic.shownFrame, rest = S.results ? S.results.baselineFrames : (S.state && S.state.rest), ordinal = (n) => n + (n % 10 === 2 && n % 100 !== 12 ? 'nd' : n % 10 === 3 && n % 100 !== 13 ? 'rd' : n % 10 === 1 && n % 100 !== 11 ? 'st' : 'th');
    const frames = rest ? rest[0] + '–' + rest[1] : '(preview)';
    $('picture-note').textContent = 'Frame ' + f + (isMoved(f) ? ' · the field is displaced in this frame (' + movedText() + ' moved)' : '') + (S.timer && Clock.step > 1 ? ' · playing every ' + ordinal(Clock.step) + ' frame' : '') + (S.timer && Clock.waitingSince && performance.now() - Clock.waitingSince > 250 ? ' · buffering' : '');
    viewKey();
  }
  function showFrame(f, mode, bitmap) {
    if (S.view !== mode || S.frame !== f) return;
    Pic.shown = bitmap; Pic.shownFrame = f; Pic.error = null; movieMessage(''); frameNote(); drawPicture();
  }
  function frameFailed(f, mode) {
    if (S.view !== mode || S.frame !== f && !S.timer) return;
    Pic.error = 'The recording frame could not be loaded.'; play(false);
    movieMessage('Frame ' + f + ' could not load. Keep the recording open in Fiji, then move the slider to retry.', true);
    $('picture-note').textContent = Pic.shownFrame ? 'Last displayed frame: ' + Pic.shownFrame : 'No movie frame loaded';
  }
  function requestFrame() {
    if (!(S.info.hasPicture || S.info.hasMovie)) return;
    if (still() || !S.info.hasMovie) {
      const image = S.view === 'background' ? Pic.background : Pic.base; Pic.shown = image || Pic.base; Pic.shownFrame = null;
      $('picture-note').textContent = ''; viewKey(); movieMessage(image ? '' : 'Loading the picture…'); drawPicture(); return;
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
  function place(w, h) {
    const iw = S.info.width, ih = S.info.height, fitted = Math.max(0.01, Math.min((w - 8) / iw, (h - 8) / ih)), z = fitted * Pic.zoom;
    return { x: (w - iw * z) / 2 + Pic.panX, y: (h - ih * z) / 2 + Pic.panY, z, fitted };
  }
  function drawPicture() {
    const canvas = $('picture'); if (!S.info || !(S.info.hasPicture || S.info.hasMovie) || !canvas.clientWidth) return;
    const { ctx, w, h } = fit(canvas); ctx.fillStyle = '#10161C'; ctx.fillRect(0, 0, w, h);
    const image = Pic.shown || Pic.base; if (!image) return;
    const at = place(w, h); ctx.imageSmoothingEnabled = Pic.zoom <= 3; ctx.drawImage(image, at.x, at.y, S.info.width * at.z, S.info.height * at.z);
    const selected = S.cell ? S.info.names.indexOf(S.cell) : -1, known = !!S.results;
    ctx.save(); ctx.translate(at.x, at.y); ctx.scale(at.z, at.z); ctx.lineJoin = 'round';
    const stroke = (path, colour, width, dash) => { ctx.setLineDash(dash ? [4 / at.z, 3 / at.z] : []); ctx.strokeStyle = colour; ctx.lineWidth = width * S.outlineWidth / at.z; ctx.stroke(path); };
    if (S.outlines) S.paths.forEach((path, i) => {
      if (i === selected) return;
      const call = callOf(S.info.names[i]); if (known && !inGroup(call, S.group, S.info.names[i])) return;
      const which = kind(call);
      stroke(path, 'rgba(0,0,0,0.55)', which === 'responded' ? 2.8 : 2.2, false);
      if (!known) stroke(path, '#ffffff', 1.3, false); else stroke(path, ON_PICTURE[which], which === 'responded' ? 1.5 : 1.1, which === 'other');
    });
    if (selected >= 0) { stroke(S.paths[selected], '#000000', 5.5, false); stroke(S.paths[selected], '#5BC0FF', 2.6, false); }
    ctx.restore();
    if (selected >= 0) {
      const b = S.bounds[selected], name = shortName(S.cell); ctx.font = '600 12px ' + getComputedStyle(document.body).fontFamily;
      const tw = ctx.measureText(name).width, x = Math.max(4, Math.min(w - tw - 12, at.x + b.x1 * at.z + 6)), y = Math.max(16, Math.min(h - 8, at.y + b.y0 * at.z + 12));
      ctx.fillStyle = 'rgba(0,0,0,0.75)'; ctx.beginPath(); ctx.roundRect(x - 4, y - 13, tw + 8, 18, 6); ctx.fill(); ctx.fillStyle = '#fff'; ctx.fillText(name, x, y);
    }
  }
  function cellAt(event) {
    const canvas = $('picture'), box = canvas.getBoundingClientRect(), at = place(canvas.clientWidth, canvas.clientHeight);
    const x = (event.clientX - box.left - at.x) / at.z, y = (event.clientY - box.top - at.y) / at.z, ctx = canvas.getContext('2d');
    ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0);
    let found = -1; for (let i = 0; i < S.paths.length && found < 0; i++) if (ctx.isPointInPath(S.paths[i], x, y)) found = i;
    ctx.restore(); return found;
  }
  function centreOn(index) {
    if (index < 0) return; const canvas = $('picture'), b = S.bounds[index]; Pic.panX = Pic.panY = 0;
    const at = place(canvas.clientWidth, canvas.clientHeight);
    Pic.panX = canvas.clientWidth / 2 - (at.x + (b.x0 + b.x1) / 2 * at.z); Pic.panY = canvas.clientHeight / 2 - (at.y + (b.y0 + b.y1) / 2 * at.z);
  }
  function zoomToCell() {
    const index = S.cell ? S.info.names.indexOf(S.cell) : -1; if (index < 0) return;
    const canvas = $('picture'), b = S.bounds[index], at = place(canvas.clientWidth, canvas.clientHeight);
    const wanted = Math.min(canvas.clientWidth / (Math.max(4, b.x1 - b.x0) * 6), canvas.clientHeight / (Math.max(4, b.y1 - b.y0) * 6));
    Pic.zoom = Math.max(1.01, Math.min(16, wanted / at.fitted)); centreOn(index); drawPicture();
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
  function referenceLevels() {
    const row = selectedCall(), baseline = number(row.Baseline_mean); if (!finite(baseline) || baseline <= -1) return [];
    const levels = [{ name: 'response baseline', value: baseline, colour: '#7c8794' }];
    for (const [key, name, colour] of [['Floor_threshold', 'minimum rise', '#ad7423'], ['Noise_threshold', 'noise threshold', '#a46ea0'], ['Threshold_increase', 'required rise', '#177465']]) {
      const value = number(row[key]); if (finite(value)) levels.push({ name, value: baseline + (1 + baseline) * value, colour });
    }
    return levels;
  }
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
  /** Mean and SEM of some cells' saved traces, frame by frame (a display summary; the tables are not changed). */
  function groupMean(names) {
    const r = S.results, T = shownTraces(), n = r.time.length, mean = new Array(n), sem = new Array(n), used = names.filter((name) => T[name]);
    for (let i = 0; i < n; i++) {
      let sum = 0, squares = 0, k = 0;
      for (const name of used) { const v = T[name][i]; if (finite(v)) { sum += v; squares += v * v; k++; } }
      mean[i] = k ? sum / k : NaN; sem[i] = k > 1 ? Math.sqrt(Math.max(0, (squares - sum * sum / k) / (k - 1)) / k) : NaN;
    }
    return { mean, sem, n: used.length };
  }
  /** The mean lines of the plot. "Selected cell and mean" and "Mean of the cells shown" with every cell shown use the mean the
   *  analysis saved; a group of cells, or responders against non-responders, is averaged here from the saved traces. */
  function plotSeries() {
    const r = S.results; if (!r) return [];
    if (S.plotMode === 'calls') {
      const names = S.info.names.filter((n) => shownTraces()[n]);
      return [['Responded', 'responded', CALL.responded], ['Did not respond', 'did not respond', CALL.not]].map(([status, text, colour]) => {
        const g = groupMean(names.filter((n) => callOf(n) === status)); return { label: text + ' ± SEM (n = ' + g.n + ')', short: text, colour, mean: g.mean, sem: g.sem, n: g.n };
      }).filter((s) => s.n > 0);
    }
    const all = isDff() ? { mean: r.mean, sem: r.sem, n: Object.keys(r.traces).length } : groupMean(S.info.names);
    if (S.plotMode === 'mean') {
      const g = S.group === 0 ? all : groupMean(S.listed), who = S.group === 0 ? 'all cells' : GROUPS[S.group].toLowerCase();
      return [{ label: who + ' ± SEM (n = ' + g.n + ')', short: 'mean of ' + who, colour: '#0072B2', mean: g.mean, sem: g.sem, n: g.n }];
    }
    return [{ label: 'mean ± SEM', short: 'mean of all cells', colour: '#0072B2', mean: all.mean, sem: all.sem, n: all.n }];
  }
  const shade = (hex, alpha) => { const c = rgb(hex); return 'rgba(' + c[0] + ',' + c[1] + ',' + c[2] + ',' + alpha + ')'; };
  const plotKind = () => (S.plotMode === 'cell' ? 'traces' : S.plotMode === 'calls' ? 'responders_vs_nonresponders' : 'mean_trace' + (S.group ? '_' + GROUPS[S.group].toLowerCase().replace(/[^a-z0-9]+/g, '_') : '')) + (isDff() ? '' : S.signal === 'raw' ? '_raw_gray_values' : S.signal === 'rawdff' ? '_dFF0_of_raw' : '_gray_levels');
  /** The numbers behind the plot as a table: time, each mean line with its SEM and number of cells, and the selected cell. */
  function plotNumbers() {
    const r = S.results, series = plotSeries(), own = S.plotMode === 'cell' && S.cell && shownTraces()[S.cell] ? shownTraces()[S.cell] : null;
    const clean = (text) => String(text).replace(/[^A-Za-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
    const head = ['Time_s']; for (const s of series) head.push(clean(s.short), clean(s.short) + '_SEM', clean(s.short) + '_n_cells'); if (own) head.push(S.cell);
    const rows = [head.join(',')];
    for (let i = 0; i < r.time.length; i++) {
      const row = [r.time[i]]; for (const s of series) row.push(finite(s.mean[i]) ? s.mean[i] : '', finite(s.sem[i]) ? s.sem[i] : '', s.n);
      if (own) row.push(finite(own[i]) ? own[i] : ''); rows.push(row.join(','));
    }
    return rows.join('\n') + '\n';
  }
  /** What the plot's legend says, as a list: lines, bands and notes. The page shows it above the plot; a saved picture has it drawn in. */
  function legendItems() {
    const r = S.results; if (!r) return [];
    const oneCell = S.plotMode === 'cell', items = plotSeries().map((s) => ({ kind: 'line', colour: s.colour, text: s.label }));
    if (oneCell && S.cell) items.push(shownTraces()[S.cell] ? { kind: 'line', colour: '#1F2933', text: shortName(S.cell) } : { kind: 'note', colour: '#9A5B00', text: shortName(S.cell) + ': no ΔF/F₀ trace' + (hasTable('corrected') ? ': choose gray levels under Signal' : '') });
    if (S.every) {
      const drawn = S.listed.filter((n) => shownTraces()[n]), yes = drawn.filter((n) => callOf(n) === 'Responded').length, no = drawn.filter((n) => callOf(n) === 'Did not respond').length;
      if (yes) items.push({ kind: 'line', colour: CALL.responded, text: 'responded (' + yes + ')' }); if (no) items.push({ kind: 'line', colour: CALL.not, text: 'did not respond (' + no + ')' });
      if (drawn.length - yes - no) items.push({ kind: 'line', colour: CALL.other, text: 'no call (' + (drawn.length - yes - no) + ')' });
    }
    const spans = responseSpans(), real = (span) => span && finite(span[0]) && finite(span[1]);
    if (real(spans.baseline)) items.push({ kind: 'band', colour: '#EDF0F4', text: 'baseline' }); if (real(spans.response)) items.push({ kind: 'band', colour: '#F5C97A', text: 'response window' });
    if (oneCell && isDff()) {
      if (real(spans.run)) items.push({ kind: 'band', colour: '#D5EEDD', text: 'response' });
      for (const level of referenceLevels().slice(1)) items.push({ kind: 'line', colour: level.colour, text: level.name });
    }
    return items;
  }
  function legend() {
    $('legend').innerHTML = legendItems().map((item) => item.kind === 'note' ? '<span class="chip" style="color:' + item.colour + '">' + escape(item.text) + '</span>'
      : '<span class="chip"><i' + (item.kind === 'band' ? ' class="band"' : '') + ' style="background:' + item.colour + '"></i>' + escape(item.text) + '</span>').join('');
  }
  /** Everything except the time cursor and the hovered trace, drawn once and reused while the cursor moves. */
  function paintTraces(ctx, w, h) {
    const r = S.results, shown = shownTraces(), time = r.time, L = 54, T = 10, W = Math.max(10, w - L - 14), H = Math.max(10, h - T - 38);
    const full0 = time[0], full1 = time[time.length - 1], x0 = Plot.x0 === null ? full0 : Plot.x0, x1 = Plot.x1 === null ? full1 : Plot.x1;
    const drawn = S.every ? S.listed.filter((n) => shown[n]) : []; Plot.drawn = drawn;
    const oneCell = S.plotMode === 'cell', series = plotSeries(); Plot.series = series;
    const spans = responseSpans(), gapLimit = samplingSupport().gap, levels = oneCell && isDff() ? referenceLevels() : [];
    if (!oneCell || !isDff()) spans.run = null;                                  // the qualifying run and the thresholds belong to one cell
    const values = (isDff() || S.signal === 'rawdff' ? [0] : []).concat(levels.map((l) => l.value)); for (const s of series) s.mean.forEach((v, i) => { if (finite(v)) { values.push(v); if (finite(s.sem[i])) values.push(v - s.sem[i], v + s.sem[i]); } });
    const own = oneCell && S.cell ? shown[S.cell] : null; if (own) for (const v of own) if (finite(v)) values.push(v);
    // with every cell drawn, the axis covers 99 % of all values (plus the mean and the selected cell); rarer extremes run off the edge
    if (drawn.length) { const rows = drawn.map((n) => shown[n]); values.push(share(rows, 0.005, false), share(rows, 0.995, false)); }
    let lo = Math.min.apply(null, values), hi = Math.max.apply(null, values); if (hi - lo < 0.001) hi = lo + 1; const pad = (hi - lo) * 0.08; lo -= pad; hi += pad;
    // limits typed under "Y axis" replace the automatic ones (either end alone, too); they belong to the signal on screen
    const typed = Plot.y[S.signal] || {}; if (finite(typed.lo)) lo = typed.lo; if (finite(typed.hi)) hi = typed.hi; if (!(hi > lo)) hi = lo + Math.max(1e-6, Math.abs(lo) * 0.1);
    const X = (t) => L + (t - x0) / (x1 - x0) * W, Y = (v) => T + (hi - v) / (hi - lo) * H;
    Plot.geom = { L, T, W, H, x0, x1, lo, hi, full0, full1 };
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); ctx.font = '12px ' + getComputedStyle(document.body).fontFamily;
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W, H); ctx.clip();
    for (const [span, colour, name] of [[spans.baseline, '#EDF0F4', 'response baseline'], [spans.response, '#FFF4E3', 'response window'], [spans.run, '#D5EEDD', 'qualifying run']]) {
      if (!span || !finite(span[0]) || !finite(span[1])) continue;
      const a = Math.max(L, X(span[0])), b = Math.min(L + W, X(span[1])); if (b <= a) continue;
      const strip = name === 'response window';
      ctx.fillStyle = strip ? '#F5C97A' : colour; ctx.fillRect(a, T, b - a, strip ? 5 : H);
      if (!strip && ctx.measureText(name).width + 8 <= b - a) { ctx.fillStyle = '#65717D'; ctx.fillText(name, a + 4, T + 14); }
    }
    const added = number(selectedCall().Stimulus_onset_s);
    if (finite(added) && added >= x0 && added <= x1) {
      const sx = Math.round(X(added)) + 0.5; ctx.strokeStyle = '#B45F06'; ctx.lineWidth = 1.4; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(sx, T); ctx.lineTo(sx, T + H); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = '#8A4B08'; ctx.textAlign = sx > L + W - 70 ? 'right' : 'left'; ctx.fillText('stimulus', sx + (sx > L + W - 70 ? -4 : 4), T + 20); ctx.textAlign = 'left';
    }
    ctx.restore();
    const ys = ticks(lo, hi, 5), xs = ticks(x0, x1, Math.max(2, Math.floor(W / 90)));
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const v of ys.values) { const y = Math.round(Y(v)) + 0.5; ctx.strokeStyle = v === 0 ? '#AEB8C2' : '#DEE4EA'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + W, y); ctx.stroke(); ctx.fillStyle = '#65717D'; ctx.fillText(tickText(v, ys.step), L - 8, y); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const v of xs.values) { const x = Math.round(X(v)) + 0.5; ctx.strokeStyle = '#DEE4EA'; ctx.beginPath(); ctx.moveTo(x, T + H); ctx.lineTo(x, T + H + 4); ctx.stroke(); ctx.fillStyle = '#65717D'; ctx.fillText(tickText(v, xs.step), x, T + H + 7); }
    ctx.fillText('Time (s)', L + W / 2, h - 15);
    ctx.save(); ctx.translate(13, T + H / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'middle'; ctx.fillText(signalName(), 0, 0); ctx.restore();
    ctx.strokeStyle = '#C5CED8'; ctx.strokeRect(L + 0.5, T + 0.5, W, H);
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W + 1, H + 1); ctx.clip(); ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    const line = (y) => { ctx.beginPath(); let pen = false; for (let i = 0; i < time.length && i < y.length; i++) { if (finite(y[i])) { const px = X(time[i]), py = Y(y[i]); if (pen && time[i] - time[i - 1] <= gapLimit) ctx.lineTo(px, py); else ctx.moveTo(px, py); pen = true; } else pen = false; } ctx.stroke(); };
    if (drawn.length) {
      ctx.globalAlpha = drawn.length > 120 ? 0.22 : drawn.length > 40 ? 0.32 : 0.48; ctx.lineWidth = 0.9;
      // non-responders first, so the responders are not hidden under them
      for (const responders of [false, true]) for (const name of drawn) { const call = callOf(name); if ((call === 'Responded') === responders) { ctx.strokeStyle = CALL[kind(call)]; line(shown[name]); } }
      ctx.globalAlpha = 1;
    }
    for (const s of series) {
      ctx.fillStyle = shade(s.colour, drawn.length ? 0.28 : 0.20);
      for (let i = 1; i < time.length; i++) if (time[i] - time[i - 1] <= gapLimit && [s.mean[i - 1], s.mean[i], s.sem[i - 1], s.sem[i]].every(finite)) {
        ctx.beginPath(); ctx.moveTo(X(time[i - 1]), Y(s.mean[i - 1] - s.sem[i - 1])); ctx.lineTo(X(time[i]), Y(s.mean[i] - s.sem[i])); ctx.lineTo(X(time[i]), Y(s.mean[i] + s.sem[i])); ctx.lineTo(X(time[i - 1]), Y(s.mean[i - 1] + s.sem[i - 1])); ctx.closePath(); ctx.fill();
      }
    }
    for (const level of levels) { ctx.strokeStyle = level.colour; ctx.lineWidth = level.name === 'required rise' ? 1.8 : 1; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(L, Y(level.value)); ctx.lineTo(L + W, Y(level.value)); ctx.stroke(); } ctx.setLineDash([]);
    for (const s of series) { if (drawn.length) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 4.2; line(s.mean); } ctx.strokeStyle = s.colour; ctx.lineWidth = 2.4; line(s.mean); }
    if (own) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 3.6; line(own); ctx.strokeStyle = '#1F2933'; ctx.lineWidth = 1.7; line(own); }
    ctx.restore();
  }
  function drawPlot() {
    const canvas = $('plot'), r = S.results; $('plot-empty').hidden = !!r;
    if (!canvas.clientWidth) return;
    const { ctx, w, h, dpr } = fit(canvas);
    if (!r) { ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); return; }
    if (Plot.dirty || Plot.layer.width !== canvas.width || Plot.layer.height !== canvas.height) {
      Plot.layer.width = canvas.width; Plot.layer.height = canvas.height; const lctx = Plot.layer.getContext('2d'); lctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      paintTraces(lctx, w, h); Plot.dirty = false;
    }
    ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.drawImage(Plot.layer, 0, 0); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const g = Plot.geom, time = r.time, X = (t) => g.L + (t - g.x0) / (g.x1 - g.x0) * g.W, Y = (v) => g.T + (g.hi - v) / (g.hi - g.lo) * g.H;
    ctx.save(); ctx.beginPath(); ctx.rect(g.L, g.T, g.W + 1, g.H + 1); ctx.clip();
    if (Plot.hover && shownTraces()[Plot.hover] && Plot.hover !== S.cell) {
      const y = shownTraces()[Plot.hover], gapLimit = samplingSupport().gap; ctx.strokeStyle = CALL[kind(callOf(Plot.hover))]; ctx.lineWidth = 2.2; ctx.lineJoin = 'round'; ctx.beginPath(); let pen = false;
      for (let i = 0; i < time.length; i++) { if (finite(y[i])) { if (pen && time[i] - time[i - 1] <= gapLimit) ctx.lineTo(X(time[i]), Y(y[i])); else ctx.moveTo(X(time[i]), Y(y[i])); pen = true; } else pen = false; } ctx.stroke();
    }
    const t = Plot.playhead != null ? Plot.playhead : S.frame <= time.length ? time[S.frame - 1] : NaN;   // glides with the playhead while playing
    if (isFinite(t) && t >= g.x0 && t <= g.x1) {
      const x = Math.round(X(t)) + 0.5; ctx.strokeStyle = 'rgba(31,41,51,0.8)'; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(x, g.T); ctx.lineTo(x, g.T + g.H); ctx.stroke();
      const own = S.plotMode === 'cell' && S.cell ? shownTraces()[S.cell] : null; let text = t.toFixed(1) + ' s';
      if (S.plotMode !== 'cell') for (const s of Plot.series) if (finite(s.mean[S.frame - 1])) text += ' · ' + s.short + ' ' + s.mean[S.frame - 1].toFixed(3);
      if (own && finite(own[S.frame - 1])) { const y = Math.max(g.T, Math.min(g.T + g.H, Y(own[S.frame - 1]))); ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(x, y, 5, 0, 7); ctx.fill(); ctx.fillStyle = '#1F2933'; ctx.beginPath(); ctx.arc(x, y, 3.2, 0, 7); ctx.fill(); text += ' · ' + shortName(S.cell) + ' ' + own[S.frame - 1].toFixed(3); }
      ctx.font = '12px ' + getComputedStyle(document.body).fontFamily; const tw = ctx.measureText(text).width + 10, bx = Math.max(g.L, Math.min(g.L + g.W - tw, x - tw / 2));
      ctx.fillStyle = '#1F2933'; ctx.beginPath(); ctx.roundRect(bx, g.T + g.H - 22, tw, 18, 6); ctx.fill(); ctx.fillStyle = '#fff'; ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; ctx.fillText(text, bx + 5, g.T + g.H - 13);
    }
    ctx.restore();
  }
  /** The drawn trace nearest to the mouse, when it is within a few pixels. */
  function traceAt(event) {
    const g = Plot.geom, r = S.results; if (!g || !r) return null;
    const box = $('plot').getBoundingClientRect(), mx = event.clientX - box.left, my = event.clientY - box.top;
    if (mx < g.L || mx > g.L + g.W || my < g.T || my > g.T + g.H) return null;
    const t = g.x0 + (mx - g.L) / g.W * (g.x1 - g.x0); let i = 0, gap = Infinity; for (let k = 0; k < r.time.length; k++) { const d = Math.abs(r.time[k] - t); if (d < gap) { gap = d; i = k; } }
    const candidates = Plot.drawn.slice(); if (S.cell && shownTraces()[S.cell] && !candidates.includes(S.cell)) candidates.push(S.cell);
    let best = null, near = 7;
    for (const name of candidates) { const v = shownTraces()[name][i]; if (!finite(v)) continue; const d = Math.abs(g.T + (g.hi - v) / (g.hi - g.lo) * g.H - my); if (d < near) { near = d; best = name; } }
    return { name: best, index: i, time: t, mx, my };
  }

  // ---------------------------------------------------------------- heatmap
  // Drawn the way the SimpleCa²⁺ app draws it: every cell in one picture, classic axes, a colour bar with ticks, the app's
  // blue–white–red for signed values (or one of its palettes from zero up), and its "peak after the baseline" ordering.
  const PALETTES = {
    plasma: ['#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786', '#d8576b', '#ed7953', '#fb9f3a', '#fdca26', '#f0f921'],
    viridis: ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'],
    magma: ['#000004', '#180f3d', '#440f76', '#721f81', '#9e2f7f', '#cd4071', '#f1605d', '#fd9668', '#feca8d', '#fcfdbf'],
    inferno: ['#000004', '#1b0c41', '#4a0c6b', '#781c6d', '#a52c60', '#cf4446', '#ed6925', '#fb9b06', '#f7d13d', '#fcffa4'],
    cividis: ['#00224e', '#123570', '#3b496c', '#575d6d', '#707173', '#8a8678', '#a59c74', '#c3b369', '#e1cc55', '#fee838'],
    signed: ['#3B4CC0', '#FFFFFF', '#B40426']
  };
  const rgb = (hex) => [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)];
  /** 256 colours along a palette. */
  function ramp(name) {
    const stops = PALETTES[name].map(rgb), out = new Uint8Array(256 * 3);
    for (let i = 0; i < 256; i++) { const at = i / 255 * (stops.length - 1), k = Math.min(stops.length - 2, Math.floor(at)), f = at - k;
      for (let c = 0; c < 3; c++) out[i * 3 + c] = Math.round(stops[k][c] + (stops[k + 1][c] - stops[k][c]) * f); }
    return out;
  }
  /** The largest value after the baseline frames and when it happens: the app's rule for ordering heatmap rows. */
  function peakAfterBaseline(trace, time, lastBaselineFrame) {
    let best = -Infinity, at = Infinity; if (!trace) return { value: -Infinity, time: Infinity };
    for (let i = Math.max(0, lastBaselineFrame); i < trace.length; i++) if (finite(trace[i]) && trace[i] > best) { best = trace[i]; at = time[i]; }
    return { value: best, time: at };
  }
  const TOP = 30, BOTTOM = 46, LEFT = 66, RIGHT = 92, ROW = 15;
  function drawHeatmap() {
    const canvas = $('heatmap'), box = $('heatmap-box'), r = S.results; $('heatmap-empty').hidden = !!r; if (!r || $('panel-heatmap').hidden) return;
    const all = S.info.names, time = r.time, rank = (n) => callOf(n) === 'Responded' ? 0 : callOf(n) === 'Did not respond' ? 1 : 2;
    const lastBaseline = r.baselineFrames && finite(r.baselineFrames[1]) ? r.baselineFrames[1] : 0;
    const peaks = {}; for (const n of all) peaks[n] = peakAfterBaseline(r.traces[n], time, lastBaseline);
    const names = all.slice(), original = (a, b) => all.indexOf(a) - all.indexOf(b);
    if (S.order === 'call') names.sort((a, b) => rank(a) - rank(b) || peaks[b].value - peaks[a].value || original(a, b));
    else if (S.order === 'tpeak') names.sort((a, b) => peaks[a].time - peaks[b].time || original(a, b));
    else if (S.order === 'amp') names.sort((a, b) => peaks[b].value - peaks[a].value || original(a, b));
    Heat.names = names;
    const rows = names.map((n) => r.traces[n] || null), have = rows.filter(Boolean), signed = S.palette === 'signed';
    // the colour range: the largest value, or a limit that covers 99 % of the values so one extreme cell cannot wash out the rest
    const largest = share(have, 1, true), most = S.range === 'full' ? largest : share(have, 0.99, true);
    const step0 = ticks(0, most > 0 ? most : 1, 5).step, limit = Math.max(step0, Math.ceil(most / step0 - 1e-9) * step0);
    const colours = ramp(S.palette), missing = [229, 229, 229];
    const callKey = [[CALL.responded, 'responded'], [CALL.not, 'did not respond'], [CALL.other, 'excluded or unscorable']];
    $('heatmap-legend').innerHTML = '<span class="chip">ΔF/F₀ <i class="scale" style="background:linear-gradient(to right,' + PALETTES[S.palette].join(',') + ')"></i>' + (signed ? '−' + fmt(limit) + ' to ' + fmt(limit) : '0 to ' + fmt(limit) + ' (below 0: darkest colour)') + '</span>' +
      '<span class="chip">strip beside the Cell axis:</span>' + callKey.map(([colour, text]) => '<span class="chip"><i class="band" style="background:' + colour + '"></i>' + text + '</span>').join('') +
      '<span class="chip"><i class="band" style="background:#e5e5e5"></i>no ΔF/F₀ trace</span><span class="chip">◀ selected cell · dashed line: stimulus · top strip: baseline and response window</span>';
    const tall = S.tall, width = Math.max(320, box.clientWidth - 2), height = tall ? TOP + names.length * ROW + BOTTOM : Math.max(220, box.clientHeight - 2);
    box.classList.toggle('fit', !tall); canvas.style.width = width + 'px'; canvas.style.height = height + 'px';
    const { ctx, w, h } = fit(canvas), W = Math.max(10, w - LEFT - RIGHT), H = Math.max(10, h - TOP - BOTTOM), rowH = H / names.length;
    const support = samplingSupport(), spans = responseSpans();
    Heat.geom = { L: LEFT, T: TOP, W, H, rows: names.length, bounds: support.bounds };
    // one pixel per sample and cell, then stretched without smoothing: sharp edges at any size
    const small = document.createElement('canvas'); small.width = time.length; small.height = names.length;
    const sctx = small.getContext('2d'), data = sctx.createImageData(time.length, names.length);
    rows.forEach((trace, c) => { for (let f = 0; f < time.length; f++) {
      const v = trace ? trace[f] : null, o = (c * time.length + f) * 4; let colour = missing;
      if (finite(v)) { const t = signed ? (Math.max(-1, Math.min(1, v / limit)) + 1) / 2 : Math.max(0, Math.min(1, v / limit)), k = Math.round(t * 255) * 3; colour = [colours[k], colours[k + 1], colours[k + 2]]; }
      data.data[o] = colour[0]; data.data[o + 1] = colour[1]; data.data[o + 2] = colour[2]; data.data[o + 3] = 255; } });
    sctx.putImageData(data, 0, 0);
    const family = getComputedStyle(document.body).fontFamily;
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); ctx.imageSmoothingEnabled = false;
    const X = (t) => LEFT + (t - time[0]) / (time[time.length - 1] - time[0]) * W;
    ctx.fillStyle = '#e5e5e5'; ctx.fillRect(LEFT, TOP, W, H);
    support.bounds.forEach((span, f) => { const x = X(span[0]), width = Math.max(0, X(span[1]) - x); if (width) ctx.drawImage(small, f, 0, 1, names.length, x, TOP, width, H); });
    // Actual timestamps determine sample widths; unsampled gaps remain blank.
    // where the baseline and the stimulus are
    ctx.font = '11px ' + family; ctx.textBaseline = 'alphabetic'; ctx.textAlign = 'left';
    for (const [span, colour, name] of [[spans.baseline, '#B8C2CC', 'response baseline'], [spans.response, '#E9A84C', 'response window']]) {
      if (!span || !finite(span[0]) || !finite(span[1])) continue; const a = Math.max(LEFT, X(span[0])), b = Math.min(LEFT + W, X(span[1])); if (b <= a) continue;
      ctx.fillStyle = colour; ctx.fillRect(a, TOP - 7, b - a, 4); if (ctx.measureText(name).width + 4 < b - a) { ctx.fillStyle = '#65717D'; ctx.fillText(name, a + 1, TOP - 11); }
    }
    if (finite(spans.response[0]) && spans.response[0] > time[0]) { const x = Math.round(X(spans.response[0])) + 0.5; ctx.strokeStyle = signed ? 'rgba(31,41,51,0.55)' : 'rgba(255,255,255,0.75)'; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(x, TOP); ctx.lineTo(x, TOP + H); ctx.stroke(); ctx.setLineDash([]); }
    // the response call of each row, and lines between the groups when the rows are sorted by call
    names.forEach((name, c) => { ctx.fillStyle = CALL[kind(callOf(name))]; ctx.fillRect(LEFT - 7, TOP + c * rowH, 4, Math.max(1, rowH - (rowH > 3 ? 0.5 : 0))); });
    if (W > 760) { ctx.font = '11px ' + family; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; let kx = LEFT;   // the same key inside the picture, for saved files
      for (const [colour, text] of callKey) { ctx.fillStyle = colour; ctx.fillRect(kx, h - 14, 10, 8); ctx.fillStyle = '#65717D'; ctx.fillText(text, kx + 14, h - 10); kx += 24 + ctx.measureText(text).width; } ctx.textBaseline = 'alphabetic'; }
    if (S.order === 'call') names.forEach((name, c) => { if (c > 0 && rank(name) !== rank(names[c - 1])) { ctx.fillStyle = '#fff'; ctx.fillRect(LEFT, TOP + c * rowH - 1, W, 2); } });
    const chosen = names.indexOf(S.cell);
    if (chosen >= 0) { const y = TOP + chosen * rowH;
      // thin rows are marked beside the picture, never drawn over; tall rows get a frame
      if (rowH >= 8) { ctx.strokeStyle = '#1F2933'; ctx.lineWidth = 1.5; ctx.strokeRect(LEFT + 0.75, y + 0.75, W - 1.5, rowH - 1.5); }
      ctx.fillStyle = '#1F2933'; ctx.beginPath(); ctx.moveTo(LEFT + W + 2, y + rowH / 2); ctx.lineTo(LEFT + W + 9, y + rowH / 2 - 5); ctx.lineTo(LEFT + W + 9, y + rowH / 2 + 5); ctx.closePath(); ctx.fill(); }
    // classic axes: a line on the left and below, ticks outside
    ctx.strokeStyle = '#333'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(LEFT - 9.5, TOP); ctx.lineTo(LEFT - 9.5, TOP + H + 0.5); ctx.lineTo(LEFT + W, TOP + H + 0.5); ctx.stroke();
    ctx.fillStyle = '#333'; ctx.font = '12px ' + family; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    const xs = ticks(time[0], time[time.length - 1], Math.max(2, Math.floor(W / 110)));
    for (const v of xs.values) { const x = Math.round(X(v)) + 0.5; ctx.beginPath(); ctx.moveTo(x, TOP + H); ctx.lineTo(x, TOP + H + 5); ctx.stroke(); ctx.fillText(tickText(v, xs.step), x, TOP + H + 8); }
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    if (tall) { ctx.font = '11px ' + family; names.forEach((name, c) => { ctx.fillStyle = name === S.cell ? '#1F2933' : '#65717D'; ctx.fillText(shortName(name), LEFT - 13, TOP + (c + 0.5) * rowH); }); }
    else { const every = ticks(0, names.length, Math.max(2, Math.floor(H / 60))).step; for (let c = every; c <= names.length; c += every) { const y = Math.round(TOP + (c - 0.5) * rowH) + 0.5; ctx.beginPath(); ctx.moveTo(LEFT - 9.5, y); ctx.lineTo(LEFT - 14, y); ctx.stroke(); ctx.fillText(String(c), LEFT - 17, y); } }
    ctx.font = '600 13px ' + family; ctx.fillStyle = '#1F2933'; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic'; ctx.fillText('Time (s)', LEFT + W / 2, h - 10);
    if (!tall) { ctx.save(); ctx.translate(16, TOP + H / 2); ctx.rotate(-Math.PI / 2); ctx.fillText('Cell', 0, 0); ctx.restore(); }
    // colour bar with ticks, as in the app
    const barX = LEFT + W + 26, barH = Math.min(190, H), barY = TOP + (tall ? 0 : Math.max(0, (H - barH) / 2)), low = signed ? -limit : 0;
    for (let i = 0; i < barH; i++) { const k = Math.round((1 - i / (barH - 1)) * 255) * 3; ctx.fillStyle = 'rgb(' + colours[k] + ',' + colours[k + 1] + ',' + colours[k + 2] + ')'; ctx.fillRect(barX, barY + i, 13, 1.5); }
    ctx.strokeStyle = '#4d4d4d'; ctx.lineWidth = 0.6; ctx.strokeRect(barX + 0.3, barY + 0.3, 12.4, barH - 0.6);
    ctx.font = '12px ' + family; ctx.fillStyle = '#333'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    const bs = ticks(low, limit, signed ? 6 : 5);
    for (const v of bs.values) { const y = barY + (1 - (v - low) / (limit - low)) * (barH - 1); ctx.fillRect(barX + 13, Math.round(y), 4, 1);
      ctx.fillText((v === limit && largest > limit ? '≥ ' : v === low && signed && largest > limit ? '≤ ' : '') + tickText(v, bs.step).replace('-', '−'), barX + 20, y); }
    ctx.font = '600 13px ' + family; ctx.fillStyle = '#1F2933'; ctx.textBaseline = 'alphabetic'; ctx.fillText('ΔF/F₀', barX - 2, barY - 9);
  }
  function heatmapAt(event) {
    const g = Heat.geom; if (!S.results || !g) return null;
    const box = $('heatmap').getBoundingClientRect(), y = event.clientY - box.top, x = event.clientX - box.left, row = Math.floor((y - TOP) / g.H * g.rows);
    if (row < 0 || row >= Heat.names.length || x < LEFT - 12 || x > LEFT + g.W) return null;
    const time = S.results.time, wanted = time[0] + Math.max(0, Math.min(1, (x - LEFT) / g.W)) * (time[time.length - 1] - time[0]), f = nearestFrame(wanted) - 1;
    return { name: Heat.names[row], frame: f, inside: x >= LEFT, gap: wanted < g.bounds[f][0] || wanted > g.bounds[f][1] };
  }

  // ---------------------------------------------------------------- tables (the saved tables themselves; long decimals are shortened for reading only)
  function table(id, data) {
    const element = $(id); element.innerHTML = ''; if (!data || !data.header) return;
    const core = ['Cell', 'Role', 'Stimulus', 'Status', 'Reason', 'QC_flags', 'Baseline_mean', 'Baseline_SD', 'Threshold_increase', 'Peak_increase', 'Response_latency_s', 'Peak_latency_s', 'Run_duration_s', 'Run_samples', 'Window_peak_local_fraction', 'Window_AUC_fraction_s'];
    let columns = data.header.map((name, index) => ({ name, index }));
    if (id === 'table-calls' && !$('all-call-columns').checked) columns = core.filter((name) => data.header.includes(name)).map((name) => ({ name, index: data.header.indexOf(name) }));
    columns = [...columns.filter((c) => c.name === 'Cell'), ...columns.filter((c) => c.name !== 'Cell')];
    // the recording's name, the role and the stimulus repeat on every row of a one-window run: leave them out so the numbers fit
    if (data.rows.length > 1) columns = columns.filter((c) => !['Field', 'Role', 'Stimulus'].includes(c.name) || new Set(data.rows.map((row) => row[c.index])).size > 1);
    // a saved table keeps its saved column names; pointing at one says what the column is
    const frozen = (name) => name === 'Cell' || (!!data.computed && name === data.header[0]);
    const about = (name) => { const m = METRIC.get(name); return data.computed ? '' : (m ? m.name + (m.unit ? ' (' + m.unit + ')' : '') + (m.measured ? ', ' + m.measured : '') + '. ' : '') + 'Saved column: ' + name; };
    const cellColumn = data.header.indexOf('Cell'), numeric = data.header.map((_, c) => data.rows.some((row) => row[c] !== '' && row[c] !== undefined) && data.rows.every((row) => row[c] === '' || row[c] === undefined || isFinite(Number(row[c]))));
    const head = element.createTHead().insertRow(), body = element.createTBody(); let sorted = -1, direction = 1;
    const fill = (rows) => { body.innerHTML = '';
      if (!rows.length) { const td = body.insertRow().insertCell(); td.colSpan = Math.max(1, columns.length); td.className = 'empty-table'; td.textContent = 'No rows were saved for this table.'; return; }
      for (const row of rows) { const tr = body.insertRow(); if (cellColumn >= 0) { tr.dataset.cell = row[cellColumn]; tr.classList.toggle('on', row[cellColumn] === S.cell); }
        for (const { name, index: c } of columns) { const value = row[c] == null ? '' : row[c], td = tr.insertCell();
          if (numeric[c] && value !== '') { const v = Number(value); td.className = 'num'; td.textContent = Number.isInteger(v) ? String(v) : String(+v.toPrecision(4)); td.title = value; }
          else if (name === 'Cell') { td.textContent = shortName(value); td.title = value; } else td.textContent = value;
          if (frozen(name)) td.classList.add('cell-column');
        }
      }
    };
    columns.forEach(({ name, index: c }) => { const th = document.createElement('th'); th.textContent = name; th.title = about(name); if (frozen(name)) th.classList.add('cell-column'); head.appendChild(th);
      th.addEventListener('click', () => { direction = sorted === c ? -direction : 1; sorted = c; for (const other of head.children) other.classList.remove('up', 'down'); th.classList.add(direction > 0 ? 'up' : 'down');
        fill(data.rows.slice().sort((a, b) => direction * (numeric[c] ? (Number(a[c] === '' ? NaN : a[c]) || 0) - (Number(b[c] === '' ? NaN : b[c]) || 0) : String(a[c]).localeCompare(String(b[c]))))); }); });
    body.addEventListener('click', (event) => { const tr = event.target.closest('tr'); if (tr && tr.dataset.cell) chooseCell(tr.dataset.cell); });
    fill(data.rows);
  }

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
  /** Fifty frames at a time on screen (the saved file has them all); cells are headed "Cell 1", with the ID in the tooltip. */
  function showData() {
    const r = S.results; if (!r) return;
    for (const option of $('data-signal').options) option.disabled = option.value !== 'dff' && !hasTable(option.value);
    if ($('data-signal').selectedOptions[0].disabled) { Data.signal = 'dff'; $('data-signal').value = 'dff'; }
    const traces = dataTraces(), names = S.info.names.filter((n) => traces[n]), total = r.time.length;
    Data.first = Math.max(0, Math.min(Data.first, Math.floor((total - 1) / Data.page) * Data.page)); const last = Math.min(total, Data.first + Data.page);
    const rows = []; for (let i = Data.first; i < last; i++) rows.push([i + 1, r.time[i]].concat(names.map((n) => finite(traces[n][i]) ? traces[n][i] : '')));
    table('table-data', { header: ['Frame', 'Time (s)'].concat(names.map(shortName)), rows });
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
    return [clean(S.info.recording), kind, /^traces/.test(kind) && S.cell ? (shortName(S.cell) === S.cell ? clean(S.cell) : clean(shortName(S.cell)) + '_' + clean(S.cell).slice(0, 13)) : '', stamp].filter(Boolean).join('_') + '.' + extension;
  }
  function exportMetadata() {
    return { schema: 'simpleca.view-export', schema_version: 1, recording: S.info.recording, plugin_version: S.info.pluginVersion,
      appliedFingerprint: S.results.appliedFingerprint || S.state.appliedFingerprint, appliedProtocol: S.results.appliedProtocol || S.results.protocol,
      manifest: S.results.manifest, settingsPending: hasPending(), displayed_window: S.window, cell: S.cell, group: GROUPS[S.group],
      frame: S.frame, time_s: frameTimes()[S.frame - 1], view: S.view, results_folder: S.results.folder || (S.state && S.state.folder) || null,
      time_range: Plot.geom ? [Plot.geom.x0, Plot.geom.x1] : [Plot.x0, Plot.x1], time_zoomed: Plot.x0 !== null, palette: S.palette, range: S.range, order: S.order,
      plot: S.plotMode, signal: S.results && isDff() ? 'dF/F0' : S.signal, y_axis: Plot.geom ? [Plot.geom.lo, Plot.geom.hi] : null, y_axis_typed: Plot.y[S.signal] || null, plotted: (Plot.series || []).map((s) => s.label), all_traces_drawn: S.every };
  }
  const Save = { format: 'png', size: '8x6', dpi: 300, compression: 'lzw' };
  const inches = () => Save.size === 'screen' ? null : Save.size.split('x').map(Number);
  /** The plot as a picture of the chosen size and resolution: its title and legend are drawn in, the time cursor is left out. */
  function plotPicture() {
    const screen = $('plot'), size = inches() || [screen.clientWidth / 96, screen.clientHeight / 96], w = Math.round(size[0] * 96), h = Math.round(size[1] * 96), k = Save.dpi / 96;
    const canvas = document.createElement('canvas'); canvas.width = Math.round(size[0] * Save.dpi); canvas.height = Math.round(size[1] * Save.dpi);
    const ctx = canvas.getContext('2d'), font = getComputedStyle(document.body).fontFamily; ctx.setTransform(k, 0, 0, k, 0, 0); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h);
    ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; ctx.fillStyle = '#1F2933'; ctx.font = '600 13px ' + font; ctx.fillText(S.info.recording, 12, 16);
    ctx.font = '12px ' + font; let x = 12, y = 38;
    for (const item of legendItems()) {
      const width = (item.kind === 'note' ? 0 : 24) + ctx.measureText(item.text).width + 16; if (x + width > w - 8 && x > 12) { x = 12; y += 18; }
      if (item.kind === 'line') { ctx.fillStyle = item.colour; ctx.fillRect(x, y - 1.5, 18, 3); }
      else if (item.kind === 'band') { ctx.fillStyle = item.colour; ctx.fillRect(x, y - 5, 14, 10); ctx.strokeStyle = '#C5CED8'; ctx.lineWidth = 1; ctx.strokeRect(x + 0.5, y - 4.5, 14, 10); }
      ctx.fillStyle = item.kind === 'note' ? item.colour : '#3E4C59'; ctx.fillText(item.text, x + (item.kind === 'note' ? 0 : 24), y); x += width;
    }
    const head = y + 14; ctx.translate(0, head); paintTraces(ctx, w, h - head);
    Plot.dirty = true; drawPlot();   // the page's own plot is measured again for the screen
    return canvas.toDataURL('image/png');
  }
  /** A saved picture or frame is as wide as Options says (inches × dpi), up to 4 times its recorded pixels; its shape is the recording's. */
  const pictureScale = () => { const size = inches(), iw = S.info.width || 1; return size ? Math.max(0.25, Math.min(4, size[0] * Save.dpi / iw)) : (Save.dpi >= 600 ? 3 : Save.dpi >= 300 ? 2 : 1); };
  /** The heatmap at the chosen width and resolution; its shape stays as on screen. */
  function heatmapPicture() {
    const canvas = $('heatmap'), size = inches(); S.exportScale = size ? size[0] * Save.dpi / Math.max(1, canvas.clientWidth) : Save.dpi / 96; drawHeatmap();
    try { return canvas.toDataURL('image/png'); } finally { S.exportScale = null; drawHeatmap(); }
  }
  const pictureOptions = () => ({ format: Save.format, dpi: Save.dpi, compression: Save.compression });
  const viewMetadata = () => ({ schema: 'simpleca.view-export', schema_version: 1, recording: S.info.recording, plugin_version: S.info.pluginVersion, view: S.view, frame: S.frame, outlines: S.outlines });
  /** What the colours of the picture mean in the view on screen, then what the outlines mean. Each: [colour, text, is an outline]. */
  function pictureKey() {
    const colours = viewColours().length ? viewColours() : [['#ffffff', 'as recorded']];
    const lines = !S.outlines ? [] : S.results ? [[ON_PICTURE.responded, 'responded', true], [ON_PICTURE.not, 'did not respond', true], [ON_PICTURE.other, 'excluded or unscorable', true]] : [['#ffffff', 'your cells', true]];
    return colours.concat(lines);
  }
  /** The whole picture in the view on screen (a movie frame is fetched at full size), with the outlines shown on the page and a
   *  key to the colours underneath. One ordinary colour picture. */
  const VIEW_NAMES = { picture: 'colour overview', raw: 'as recorded', flat: 'background removed', green: 'GFP', magenta: 'response in magenta', change: 'change from baseline' };
  async function viewCanvas(frame, labelled) {
    const f = frame || S.frame; let image = S.view === 'background' ? Pic.background : Pic.base;
    if (!still()) { const reply = await fetch('api/frame?f=' + f + '&mode=' + S.view + '&v=' + encodeURIComponent(Buf.version || ''), { cache: 'no-store' }); if (!reply.ok) throw new Error('Frame ' + f + ' could not be read.'); image = await createImageBitmap(await reply.blob()); }
    if (!image) throw new Error('The picture is not loaded yet.');
    const iw = S.info.width, ih = S.info.height, scale = pictureScale(), key = labelled ? [] : pictureKey(), font = getComputedStyle(document.body).fontFamily;
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
    if (labelled) {   // a frame of a set: its number, time and view in the corner
      const time = frameTimes()[f - 1], text = 'Frame ' + f + (finite(time) ? ' · ' + time.toFixed(1) + ' s' : '') + ' · ' + (VIEW_NAMES[S.view] || S.view);
      ctx.font = '600 ' + Math.round(size * 1.1) + 'px ' + font; const tw = ctx.measureText(text).width;
      ctx.fillStyle = 'rgba(0,0,0,0.7)'; ctx.fillRect(10 * unit, 10 * unit, tw + 20 * unit, size * 1.9); ctx.fillStyle = '#fff'; ctx.fillText(text, 20 * unit, 10 * unit + size * 0.98);
    } else { ctx.fillStyle = '#AEB8C2'; ctx.font = Math.round(size * 0.8) + 'px ' + font; ctx.fillText(S.info.recording, 14 * unit, ih + keyHeight - gap * 0.55); }
    return canvas;
  }
  async function viewPicture(frame, labelled) { return (await viewCanvas(frame, labelled)).toDataURL('image/png'); }
  /** Where n frames of w x h go in one picture: a grid about as wide as it is tall, each tile at most its own size and the
   *  whole at most COMBINED_MAX pixels wide. Returns { cols, rows, tw, th, gap, width, height }. */
  const COMBINED_MAX = 6000;
  function combinedLayout(n, w, h) {
    const cols = Math.max(1, Math.ceil(Math.sqrt(n * h / w) * w / h)), rows = Math.ceil(n / cols), gap = Math.max(2, Math.round(w / 200));
    const tw = Math.max(1, Math.min(w, Math.floor((COMBINED_MAX - (cols + 1) * gap) / cols))), th = Math.max(1, Math.round(tw * h / w));
    return { cols, rows, tw, th, gap, width: cols * tw + (cols + 1) * gap, height: rows * th + (rows + 1) * gap };
  }
  const viewKind = () => 'picture_' + ({ picture: 'colour_overview', background: 'background_areas', raw: 'as_recorded', flat: 'background_removed', green: 'gfp', magenta: 'response_in_magenta', change: 'change_from_baseline' }[S.view] || S.view) + (still() ? '' : '_frame' + S.frame);
  const badSave = (error) => /fetch/i.test(error.message) ? 'Fiji is not connected. Keep Fiji open, or reopen the analysis from Fiji.' : error.message;
  /** Download through the browser: the same file Save writes (format, size, resolution), plus how it was made. */
  async function save(picture, kind, noteId) {
    const note = $(noteId); note.className = 'save-note'; note.textContent = 'Preparing…';
    try {
      const metadata = JSON.stringify(Object.assign(S.results ? exportMetadata() : viewMetadata(), { picture: Object.assign({ size: Save.size }, pictureOptions()) }), null, 2), png = await picture();
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
      const metadata = Object.assign(S.results ? exportMetadata() : viewMetadata(), { picture: Object.assign({ size: Save.size }, pictureOptions()) }), name = exportName(kind, 'png').replace(/\.png$/, '');
      const reply = await api('export', Object.assign({ name, png: await picture(), csv: csv || null, metadata }, pictureOptions()));
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
    const showY = () => { const now = Plot.y[S.signal] || {}; $('y-min').value = finite(now.lo) ? now.lo : ''; $('y-max').value = finite(now.hi) ? now.hi : ''; };
    for (const [id, end] of [['y-min', 'lo'], ['y-max', 'hi']]) $(id).addEventListener('input', () => {
      const text = $(id).value.trim(), value = Number(text), now = Plot.y[S.signal] || (Plot.y[S.signal] = {});
      now[end] = text !== '' && isFinite(value) ? value : null; Plot.dirty = true; drawPlot();
    });
    $('y-auto').addEventListener('click', () => { delete Plot.y[S.signal]; showY(); Plot.dirty = true; drawPlot(); });
    const setSignal = (value) => { S.signal = Data.signal = value; $('signal').value = value; $('data-signal').value = value; showY(); Plot.dirty = true; drawPlot(); legend(); if (S.results && !$('panel-data').hidden) showData(); };
    $('signal').addEventListener('change', () => setSignal($('signal').value));
    $('save-plot-folder').addEventListener('click', () => saveToFolder(plotPicture, plotKind(), 'plot-saved', plotNumbers()));
    $('save-heatmap-folder').addEventListener('click', () => saveToFolder(heatmapPicture, 'heatmap', 'heatmap-saved'));
    $('save-picture-folder').addEventListener('click', () => saveToFolder(viewPicture, viewKind(), 'picture-saved'));
    $('save-picture').addEventListener('click', () => save(viewPicture, viewKind(), 'picture-saved'));
    $('outline-width').addEventListener('change', () => { S.outlineWidth = Number($('outline-width').value) || 1; store('step4.outline', $('outline-width').value); drawPicture(); });
    // One Export menu per panel holds its save buttons. Choosing an entry closes the menu; so does a click elsewhere or Escape.
    // A panel opened from an entry (Options, Save frames) is placed against the menu's own button, which stays on screen.
    const anchorOf = (button) => button.closest('details.export') || button;
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
    // every Nth frame of the view on screen, each its own labelled picture, in a folder of their own
    const framesPanel = $('frames-options'), framesEvery = () => Math.round(Number($('frames-every').value));
    const framesCount = () => { const every = framesEvery(), ok = every >= 1 && every <= S.info.frames; $('frames-go').disabled = !ok;
      $('frames-count').textContent = ok ? plural(Math.ceil(S.info.frames / every), 'picture') + (Math.ceil(S.info.frames / every) > 1 ? ' and one with all of them' : '') + ', ' + (Save.format === 'tiff' ? 'TIFF ' + Save.compression.toUpperCase() : 'PNG') + ', ' + Save.dpi + ' dpi.' : 'A whole number from 1 to ' + S.info.frames + '.'; };
    $('frames-every').addEventListener('input', framesCount); $('frames-cancel').addEventListener('click', () => { framesPanel.hidden = true; });
    $('frames-every').addEventListener('keydown', (event) => { if (event.key === 'Enter' && !$('frames-go').disabled) $('frames-go').click(); if (event.key === 'Escape') framesPanel.hidden = true; });
    $('save-frames').addEventListener('click', (event) => {
      const note = $('picture-saved'); note.className = 'save-note';
      if (still()) { note.className = 'save-note bad'; note.textContent = 'Choose a recording view first: this view is one still picture.'; return; }
      if (!framesPanel.hidden) { framesPanel.hidden = true; return; }
      $('frames-every').value = stored('step4.every') || '50'; $('frames-every').max = S.info.frames; framesCount(); framesPanel.hidden = false;
      const box = anchorOf($('save-frames')).getBoundingClientRect(); framesPanel.style.left = Math.max(8, Math.min(window.innerWidth - 316, box.right - 300)) + 'px'; framesPanel.style.top = ''; framesPanel.style.bottom = '';
      if (box.top > window.innerHeight / 2) framesPanel.style.bottom = (window.innerHeight - box.top + 6) + 'px'; else framesPanel.style.top = (box.bottom + 6) + 'px';
      $('frames-every').focus(); $('frames-every').select(); event.stopPropagation();
    });
    document.addEventListener('click', (event) => { if (!framesPanel.hidden && !framesPanel.contains(event.target) && event.target.id !== 'save-frames') framesPanel.hidden = true; });
    // every Nth frame of the view on screen, each its own labelled picture, in a folder of their own
    $('frames-go').addEventListener('click', async () => {
      const note = $('picture-saved'), every = framesEvery(); note.className = 'save-note'; framesPanel.hidden = true;
      store('step4.every', String(every));
      const frames = []; for (let f = 1; f <= S.info.frames; f += every) frames.push(f);
      const stamp = exportName('x', 'png').match(/(\d{8}-\d{6})\.png$/)[1], folder = 'frames_' + viewKind().replace(/^picture_/, '').replace(/_frame\d+$/, '') + '_every_' + every + '_' + stamp, width = String(S.info.frames).length;
      $('save-frames').disabled = true;
      try {
        let reply = null, combined = null, layout = null;
        for (let i = 0; i < frames.length; i++) {
          note.textContent = 'Saving frame ' + frames[i] + ' (' + (i + 1) + ' of ' + frames.length + ')…';
          const canvas = await viewCanvas(frames[i], true);
          reply = await api('export', Object.assign({ folder, name: 'frame_' + String(frames[i]).padStart(width, '0'), png: canvas.toDataURL('image/png') }, pictureOptions()));
          // the same labelled frame, as one tile of the combined picture
          if (!combined) { layout = combinedLayout(frames.length, canvas.width, canvas.height); combined = document.createElement('canvas'); combined.width = layout.width; combined.height = layout.height; const c = combined.getContext('2d'); c.fillStyle = '#000'; c.fillRect(0, 0, layout.width, layout.height); c.imageSmoothingQuality = 'high'; }
          combined.getContext('2d').drawImage(canvas, layout.gap + (i % layout.cols) * (layout.tw + layout.gap), layout.gap + Math.floor(i / layout.cols) * (layout.th + layout.gap), layout.tw, layout.th);
        }
        if (combined && frames.length > 1) { note.textContent = 'Saving the combined picture…'; await api('export', Object.assign({ folder, name: 'all_frames_combined', png: combined.toDataURL('image/png') }, pictureOptions())); }
        await api('export', { folder, name: 'about_these_frames', metadata: Object.assign(viewMetadata(), { frames, every, combined_picture: combined && frames.length > 1 ? { file: 'all_frames_combined', columns: layout.cols, rows: layout.rows, tile_pixels: [layout.tw, layout.th] } : null, outlines_shown: S.outlines, outline_colours: S.results ? { responded: ON_PICTURE.responded, did_not_respond: ON_PICTURE.not, excluded_or_unscorable: ON_PICTURE.other } : { your_cells: '#ffffff' }, picture: pictureOptions() }) });
        note.textContent = 'Saved ' + plural(frames.length, 'frame') + (frames.length > 1 ? ' and one combined picture' : '') + ' in plots/' + reply.folderName + ' ';
        const open = document.createElement('button'); open.type = 'button'; open.textContent = 'Open plots folder'; open.addEventListener('click', () => api('open', { what: 'plots' }).catch(() => { })); note.appendChild(open);
      } catch (error) { note.className = 'save-note bad'; note.textContent = 'Stopped: ' + badSave(error); }
      finally { $('save-frames').disabled = false; }
    });
    // how pictures are saved: format, size, resolution, TIFF compression (kept with Fiji's preferences)
    const saveBoxes = { format: 'save-format', size: 'save-size', dpi: 'save-dpi', compression: 'save-compression' };
    const showSaveOptions = () => {
      for (const key of Object.keys(saveBoxes)) $(saveBoxes[key]).value = String(Save[key]);
      $('save-compression-label').hidden = Save.format !== 'tiff';
      const size = inches(), wide = S.info ? Math.round(S.info.width * pictureScale()) : 0;
      $('save-pixels').textContent = (size ? 'Plot ' + Math.round(size[0] * Save.dpi) + ' × ' + Math.round(size[1] * Save.dpi) + ' px.' : 'Plot as on screen, at ' + Save.dpi + ' dpi.') + (wide ? ' Pictures and frames ' + wide + ' px wide' + (size && wide < Math.round(size[0] * Save.dpi) ? ' (4× the recording, the most it allows).' : '.') : '');
    };
    for (const key of Object.keys(saveBoxes)) $(saveBoxes[key]).addEventListener('change', () => { Save[key] = key === 'dpi' ? Number($(saveBoxes[key]).value) : $(saveBoxes[key]).value; store('step4.' + key, String(Save[key])); showSaveOptions(); });
    for (const button of document.querySelectorAll('.save-options-toggle')) button.addEventListener('click', (event) => {
      const panel = $('save-options'), open = panel.hidden || panel.dataset.by !== String(Array.from(document.querySelectorAll('.save-options-toggle')).indexOf(button));
      panel.hidden = !open; if (!open) return;
      panel.dataset.by = String(Array.from(document.querySelectorAll('.save-options-toggle')).indexOf(button)); showSaveOptions();
      const box = anchorOf(button).getBoundingClientRect(); panel.style.left = Math.max(8, Math.min(window.innerWidth - 316, box.right - 300)) + 'px'; panel.style.top = ''; panel.style.bottom = '';
      if (box.top > window.innerHeight / 2) panel.style.bottom = (window.innerHeight - box.top + 6) + 'px'; else panel.style.top = (box.bottom + 6) + 'px';
      event.stopPropagation();
    });
    document.addEventListener('click', (event) => { const panel = $('save-options'); if (!panel.hidden && !panel.contains(event.target) && !event.target.classList.contains('save-options-toggle')) panel.hidden = true; });
    const metricsName = () => exportName($('metrics-table').value === 'summary' ? 'metrics_summary' : 'table_' + $('metrics-table').value, 'csv');
    $('metrics-download').addEventListener('click', () => download(new Blob([csvText(metricRows())], { type: 'text/csv' }), metricsName()));
    $('metrics-save').addEventListener('click', async () => {
      const note = $('metrics-saved'); note.className = 'save-note'; note.textContent = 'Saving…';
      try {
        const reply = await api('export', { name: metricsName().replace(/\.csv$/, ''), csv: csvText(metricRows()) });
        note.textContent = 'Saved in ' + reply.folderName + ': ' + reply.saved[0] + ' ';
        const open = document.createElement('button'); open.type = 'button'; open.textContent = 'Open folder'; open.addEventListener('click', () => api('open', { what: 'plots' }).catch(() => { })); note.appendChild(open);
      } catch (error) { note.className = 'save-note bad'; note.textContent = 'Not saved: ' + (/fetch/i.test(error.message) ? 'Fiji is not connected.' : error.message); }
    });
    $('data-signal').addEventListener('change', () => setSignal($('data-signal').value));
    $('data-back').addEventListener('click', () => { Data.first -= Data.page; showData(); }); $('data-forward').addEventListener('click', () => { Data.first += Data.page; showData(); });
    $('data-download').addEventListener('click', () => download(new Blob([dataCsv()], { type: 'text/csv' }), exportName(dataKind(), 'csv')));
    $('data-save').addEventListener('click', async () => {
      const note = $('data-saved'); note.className = 'save-note'; note.textContent = 'Saving…';
      try {
        const reply = await api('export', { name: exportName(dataKind(), 'csv').replace(/\.csv$/, ''), csv: dataCsv(), metadata: exportMetadata() });
        note.textContent = 'Saved in ' + reply.folderName + ': ' + reply.saved[0] + ' ';
        const open = document.createElement('button'); open.type = 'button'; open.textContent = 'Open folder'; open.addEventListener('click', () => api('open', { what: 'plots' }).catch(() => { })); note.appendChild(open);
      } catch (error) { note.className = 'save-note bad'; note.textContent = 'Not saved: ' + error.message; }
    });
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
      try { let parsed; try { parsed = JSON.parse(await file.text()); } catch (bad) { throw new Error(file.name + ' is not a valid JSON file (' + bad.message + '). Choose a protocol exported from SimpleCa²⁺ or this page.'); }
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
    $('quality-filter').addEventListener('click', () => { S.qualityFlag = null; S.group = -1; showGroup(5); showQuality(); showTab('review'); });
    $('group').addEventListener('change', () => { if (Number($('group').value) !== 5 && S.qualityFlag) { S.qualityFlag = null; showQuality(); } });
    $('all-call-columns').addEventListener('change', () => table('table-calls', S.results.tables.calls));
    $('response-window').addEventListener('change', () => { S.window = $('response-window').value; updateWindow(); });
    $('metrics-table').addEventListener('change', showMetricTable);
    $('below-zero').addEventListener('change', () => { S.withoutBelowZero = $('below-zero').checked; showTiles(); });
    $('group').addEventListener('change', () => { S.group = Number($('group').value); listCells(); });
    $('cell').addEventListener('change', () => chooseCell($('cell').value || null));
    $('previous').addEventListener('click', () => stepCell(-1)); $('next').addEventListener('click', () => stepCell(1));
    $('every').addEventListener('change', () => { S.every = $('every').checked; Plot.dirty = true; drawPlot(); legend(); });
    $('view').addEventListener('change', () => { play(false); S.view = $('view').value; S.viewChosen = true; updateTransport(); requestFrame(); });
    $('outlines').addEventListener('change', () => { S.outlines = $('outlines').checked; drawPicture(); });
    $('zoom-cell').addEventListener('click', zoomToCell); $('zoom-all').addEventListener('click', () => { Pic.zoom = 1; Pic.panX = Pic.panY = 0; drawPicture(); });
    $('play').addEventListener('click', () => play(!S.timer)); $('speed').addEventListener('change', () => { if (S.timer) play(true); });
    $('time').addEventListener('input', () => setFrame(Number($('time').value), true));
    $('frame-back').addEventListener('click', () => setFrame(S.frame - 1, true)); $('frame-forward').addEventListener('click', () => setFrame(S.frame + 1, true));
    $('order').addEventListener('change', () => { S.order = $('order').value; drawHeatmap(); });
    $('palette').addEventListener('change', () => { S.palette = $('palette').value; drawHeatmap(); });
    $('range').addEventListener('change', () => { S.range = $('range').value; drawHeatmap(); });
    $('tall').addEventListener('change', () => { S.tall = $('tall').checked; drawHeatmap(); });
    $('save-plot').addEventListener('click', () => save(plotPicture, plotKind(), 'plot-saved')); $('save-heatmap').addEventListener('click', () => save(heatmapPicture, 'heatmap', 'heatmap-saved'));

    const picture = $('picture');
    picture.addEventListener('wheel', (event) => { if (!(event.ctrlKey || event.metaKey)) return; event.preventDefault(); const box = picture.getBoundingClientRect(), mx = event.clientX - box.left, my = event.clientY - box.top, before = place(picture.clientWidth, picture.clientHeight);
      const ix = (mx - before.x) / before.z, iy = (my - before.y) / before.z; Pic.zoom = Math.max(1, Math.min(16, Pic.zoom * Math.pow(1.0015, -event.deltaY)));
      if (Pic.zoom === 1) Pic.panX = Pic.panY = 0; else { const after = place(picture.clientWidth, picture.clientHeight); Pic.panX += mx - (after.x + ix * after.z); Pic.panY += my - (after.y + iy * after.z); } drawPicture(); }, { passive: false });
    picture.addEventListener('mousedown', (event) => { Pic.drag = { x: event.clientX, y: event.clientY }; Pic.moved = false; });
    window.addEventListener('mousemove', (event) => { if (!Pic.drag) return; if (Math.abs(event.clientX - Pic.drag.x) + Math.abs(event.clientY - Pic.drag.y) > 3) Pic.moved = true;
      if (Pic.moved) { Pic.panX += event.clientX - Pic.drag.x; Pic.panY += event.clientY - Pic.drag.y; Pic.drag = { x: event.clientX, y: event.clientY }; drawPicture(); } });
    window.addEventListener('mouseup', () => { Pic.drag = null; Plot.scrubbing = false; });
    picture.addEventListener('click', (event) => { if (Pic.moved) return; const i = cellAt(event); if (i >= 0) chooseCell(S.info.names[i]); });
    picture.addEventListener('dblclick', () => { Pic.zoom = 1; Pic.panX = Pic.panY = 0; drawPicture(); });
    picture.addEventListener('mousemove', (event) => { if (Pic.drag) return; const i = cellAt(event); tip(i >= 0 ? fullName(S.info.names[i]) + (callOf(S.info.names[i]) ? ' · ' + callOf(S.info.names[i]) : '') : '', event.clientX, event.clientY); });
    picture.addEventListener('mouseleave', () => tip(''));

    const plot = $('plot');
    plot.addEventListener('pointermove', (event) => {
      if (Plot.scrubbing) {
        const g = Plot.geom, box = plot.getBoundingClientRect(); if (g) seek(g.x0 + Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W)) * (g.x1 - g.x0)); tip(''); return;
      }
      const at = traceAt(event); if (!at) { tip(''); if (Plot.hover) { Plot.hover = null; drawPlot(); } return; }
      if (at.name !== Plot.hover) { Plot.hover = at.name; drawPlot(); }
      tip(at.name ? fullName(at.name) + ' · ' + shownTraces()[at.name][at.index].toFixed(3) + ' at ' + S.results.time[at.index].toFixed(1) + ' s' + (callOf(at.name) ? ' · ' + callOf(at.name) : '') : '', event.clientX, event.clientY);
    });
    plot.addEventListener('pointerleave', () => { tip(''); if (Plot.hover) { Plot.hover = null; drawPlot(); } });
    plot.addEventListener('pointerdown', (event) => { const at = traceAt(event); if (!at) return; event.preventDefault(); if (at.name && at.name !== S.cell) chooseCell(at.name); Plot.scrubbing = true; plot.setPointerCapture(event.pointerId); seek(at.time); });
    for (const type of ['pointerup', 'pointercancel', 'lostpointercapture']) plot.addEventListener(type, () => { Plot.scrubbing = false; });
    plot.addEventListener('wheel', (event) => { const g = Plot.geom; if (!g || !(event.ctrlKey || event.metaKey)) return; event.preventDefault(); const box = plot.getBoundingClientRect(), share01 = Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W));
      const centre = g.x0 + share01 * (g.x1 - g.x0), factor = Math.pow(1.0015, event.deltaY), span = Math.min(g.full1 - g.full0, Math.max((g.full1 - g.full0) / 60, (g.x1 - g.x0) * factor));
      let a = centre - share01 * span, b = a + span; if (a < g.full0) { b += g.full0 - a; a = g.full0; } if (b > g.full1) { a -= b - g.full1; b = g.full1; }
      Plot.x0 = Math.max(g.full0, a); Plot.x1 = Math.min(g.full1, b); if (Plot.x1 - Plot.x0 >= g.full1 - g.full0 - 1e-9) Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); }, { passive: false });
    plot.addEventListener('dblclick', () => { Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); });

    const heatmap = $('heatmap');
    heatmap.addEventListener('click', (event) => { const at = heatmapAt(event); if (at) { chooseCell(at.name); if (at.inside) setFrame(at.frame + 1, true); showTab('review'); } });
    heatmap.addEventListener('mousemove', (event) => { const at = heatmapAt(event); if (!at) { tip(''); return; } const trace = S.results.traces[at.name], v = trace ? trace[at.frame] : null;
      tip(fullName(at.name) + (callOf(at.name) ? ' · ' + callOf(at.name) : '') + (at.inside ? ' · ' + (at.gap ? 'acquisition gap; nearest frame' : finite(v) ? v.toFixed(3) : 'no value') + ' at ' + S.results.time[at.frame].toFixed(1) + ' s' : ''), event.clientX, event.clientY); });
    heatmap.addEventListener('mouseleave', () => tip(''));

    document.addEventListener('keydown', (event) => { const tag = (event.target.tagName || '').toLowerCase(); if (tag === 'input' || tag === 'select' || tag === 'textarea' || $('panel-review').hidden) return;
      if (event.key === 'Escape' && S.cell && !event.defaultPrevented) { chooseCell(null); event.preventDefault(); }
      else if (event.key === 'ArrowLeft') { stepCell(-1); event.preventDefault(); } else if (event.key === 'ArrowRight') { stepCell(1); event.preventDefault(); } else if (event.key === ' ' && S.info && S.info.hasMovie && S.view !== 'picture' && tag !== 'button') { play(!S.timer); event.preventDefault(); } });
    window.addEventListener('resize', redraw);
    if (window.ResizeObserver) { new ResizeObserver(redraw).observe($('panel-review')); new ResizeObserver(() => drawHeatmap()).observe($('heatmap-box')); }
  }
  function redraw() { Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); }

  wire();
  start().catch((problem) => { $('lost').classList.add('on'); $('status').textContent = problem.message; });
  window.calciumStep4 = { S, Plot, Pic, Heat, Buf, Clock, rateText, samplesFor, resetFrames: () => resetFrames(Buf.version), pick, setFrame, showGroup, showTab, collectDraft, frameTimes, nearestFrame, samplingSupport, referenceLevels, exportMetadata, heatmapAt };   // for the browser test
})();

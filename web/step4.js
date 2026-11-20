// Step 4 page. Everything shown here is read from the tables the analysis saved, or from the recording's own frames.
// The page draws; it does not calculate results. Settings go back to the plugin, which runs the same engine as the Fiji window.
'use strict';
(function () {
  const $ = (id) => document.getElementById(id);
  const CALL = { responded: '#177465', not: '#925223', other: '#8A94A0' };
  const ON_PICTURE = { responded: '#3ddbb5', not: '#ffa24a', other: '#aeb8c2' };
  const GROUPS = ['All cells', 'Responded', 'Did not respond', 'Excluded', 'Unscorable', 'Needs quality review'];
  const kind = (call) => call === 'Responded' ? 'responded' : call === 'Did not respond' ? 'not' : 'other';
  const inGroup = (call, group, name) => group === 0 || (group === 1 && call === 'Responded') || (group === 2 && call === 'Did not respond') ||
    (group === 3 && call === 'Excluded') || (group === 4 && call === 'Unscorable') || (group === 5 && qualityFor(name).length > 0);
  const S = { info: null, state: null, results: null, version: -1, group: 0, cell: null, listed: [], frame: 1, view: 'raw', viewChosen: false,
    every: true, outlines: true, order: 'call', palette: 'signed', range: 'robust', tall: false, timer: null, lost: 0, note: null, polling: false, paths: [], boxes: [], draft: null, draftAck: '', revision: 0, draftTimer: null, draftQueue: Promise.resolve(), window: null, saving: false, initialDraft: null, protocolDrafts: {}, quality: {}, metrics: {}, closed: false };
  const Plot = { dirty: true, layer: document.createElement('canvas'), x0: null, x1: null, hover: null, geom: null, scrubbing: false, drawn: [] };
  const Pic = { base: null, shown: null, zoom: 1, panX: 0, panY: 0, pending: null, loading: false, drag: null, moved: false, shownFrame: null, error: null };
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
  function exactInput(box, value, scale) {
    const shown = box.id === 'sampling_hz' ? String(value == null ? '' : value) : scale ? fmt(Number(value) * scale) : fmt(value);
    box.value = shown; box.dataset.shown = shown; box._exact = value; box._scale = scale || 1;
    box.title = value == null ? '' : 'Stored value: ' + value + '. Untouched fields keep their full precision.';
  }
  function readInput(box, optional) {
    if (box.value === box.dataset.shown) return box._exact;
    const text = box.value.trim(); if (optional && text === '') return null;
    const value = text === '' ? NaN : Number(text);
    return Number.isFinite(value) ? value / box._scale : text;
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
    const container = $('windows-editor'), keep = [S.controlLabel || (S.controlLabel = $('control').closest('label')), S.controlHelp || (S.controlHelp = $('control-help'))]; keep.forEach((e) => e.remove()); container.innerHTML = '';
    const legacy = Number($('protocol-mode').value) === 1;
    S.draft.protocol.windows.forEach((row, index) => {
      const w = expandWindow(row, legacy), section = document.createElement('div'); section.className = 'window-editor'; section.dataset.index = index;
      const heading = document.createElement('h3'); heading.textContent = index === 0 ? 'Primary stimulus' : 'Positive control'; section.appendChild(heading);
      const add = (target, key, title, opts) => { const label = windowInput(index, key, title, opts); target.appendChild(label); const box = label.lastChild;
        if (box.tagName === 'SELECT' || key === 'Stimulus') box.value = w[key] == null ? '' : String(w[key]); else exactInput(box, w[key], opts && opts.scale); return box; };
      add(section, 'Stimulus', 'Name'); add(section, 'Stimulus_onset_s', 'Stimulus added at (s)');
      add(section, 'Baseline_start_s', 'Response baseline from (s)'); add(section, 'Baseline_end_s', 'Response baseline to (s)');
      add(section, 'Start_s', 'Score response from (s)'); add(section, 'End_s', 'Score response to (s)');
      const summary = document.createElement('p'); summary.className = 'rule-summary'; summary.id = 'rule-summary-' + index; section.appendChild(summary);
      const details = document.createElement('details'); const summaryTitle = document.createElement('summary'); summaryTitle.textContent = 'Adjust the response rule'; details.appendChild(summaryTitle);
      add(details, 'Floor_local_fraction', 'Minimum rise (%)', { scale: 100 }); add(details, 'SD_multiplier', 'Baseline SD multiplier');
      add(details, 'Duration_mode', 'Persistence measured in', { choices: [['seconds', 'Seconds'], ['samples', 'Consecutive samples']] });
      add(details, 'Min_duration_s', 'Minimum duration (s)'); add(details, 'Consecutive_samples', 'Consecutive samples');
      add(details, 'Min_baseline_samples', 'Minimum baseline samples'); add(details, 'Max_gap_s', 'Maximum gap (s; blank = auto)');
      add(details, 'Expected_interval_s', 'Expected interval (s; blank = auto)'); add(details, 'Noise_floor_local_fraction', 'Noise floor (%)', { scale: 100 });
      section.appendChild(details); container.appendChild(section);
      if (index === 0) keep.forEach((e) => container.appendChild(e));
    });
    if (!S.draft.protocol.windows.length) keep.forEach((e) => container.appendChild(e));
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
        field.disabled = busy || (legacy && (['Baseline_end_s', 'Stimulus_onset_s', 'Duration_mode', 'Min_duration_s', 'Min_baseline_samples', 'Noise_floor_local_fraction', 'Expected_interval_s'].includes(key) || (i === '1' && key === 'Stimulus'))) || (!legacy && ((key === 'Min_duration_s' && mode !== 'seconds') || (key === 'Consecutive_samples' && mode !== 'samples')));
      }
      $('rule-summary-' + i).textContent = 'Rise ≥ ' + box('Floor_local_fraction').value + '% and ≥ ' + box('SD_multiplier').value + ' × baseline SD, for ' + (mode === 'seconds' ? box('Min_duration_s').value + ' s.' : box('Consecutive_samples').value + ' consecutive samples.') + (legacy ? ' Legacy rule.' : '');
    }
    $('rule-help').textContent = legacy ? 'Legacy reproduction: the baseline ends at the response start; the optional control is KCl. Sample counts depend on the acquisition rate.' : 'Starter settings are provisional. Check them against your controls before using responder percentages.';
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
    for (const key of ['sampling_hz', 'baseline_start', 'baseline_end']) exactInput($(key), draft.values[key]);
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
    const errors = [], warnings = [], fields = new Set(), add = (field, message) => { errors.push(message); if (field) fields.add(field); };
    const numeric = (v) => v !== null && v !== '' && typeof v !== 'boolean' && Number.isFinite(Number(v));
    const value = draft.values, n = S.info.frames;
    if (!draft.timestamps && (!numeric(value.sampling_hz) || Number(value.sampling_hz) <= 0)) add('sampling_hz', 'Frame rate must be a positive number, or supply valid timestamps.');
    if (!numeric(value.baseline_start) || !Number.isInteger(Number(value.baseline_start)) || Number(value.baseline_start) < 1) add('baseline_start', 'First baseline frame must be a whole number from 1.');
    if (!numeric(value.baseline_end) || !Number.isInteger(Number(value.baseline_end)) || Number(value.baseline_end) < Number(value.baseline_start) || Number(value.baseline_end) >= n) add('baseline_end', 'Last baseline frame must follow the first and be before frame ' + n + '.');
    const legacy = draft.protocol.schema_version === 1, names = new Set();
    draft.protocol.windows.forEach((original, index) => {
      const row = expandWindow(original, legacy), prefix = 'window-' + index + '-', label = index ? 'Control: ' : 'Stimulus: ';
      const required = NUMBERS.filter((key) => !['Max_gap_s', 'Expected_interval_s'].includes(key));
      for (const key of required) if (!numeric(row[key])) add(prefix + key, label + key.replace(/_/g, ' ') + ' needs a number.');
      const name = String(row.Stimulus || '').trim().toLowerCase();
      if (!name || names.has(name)) add(prefix + 'Stimulus', label + 'use a nonempty, distinct name.'); names.add(name);
      if (numeric(row.Baseline_start_s) && Number(row.Baseline_start_s) < 0) add(prefix + 'Baseline_start_s', label + 'baseline cannot start before 0 s.');
      for (const [left, right, strict, message] of [
        ['Baseline_start_s', 'Baseline_end_s', true, 'baseline end must follow baseline start.'],
        ['Baseline_end_s', 'Stimulus_onset_s', false, 'baseline must finish at or before stimulus addition.'],
        ['Stimulus_onset_s', 'Start_s', false, 'scoring cannot start before stimulus addition.'],
        ['Start_s', 'End_s', true, 'response end must follow response start.']]) {
        if (numeric(row[left]) && numeric(row[right]) && (strict ? Number(row[right]) <= Number(row[left]) : Number(row[right]) < Number(row[left]))) add(prefix + right, label + message);
      }
      for (const key of ['Floor_local_fraction', 'SD_multiplier', 'Noise_floor_local_fraction']) if (numeric(row[key]) && Number(row[key]) < 0) add(prefix + key, label + 'thresholds cannot be negative.');
      if (numeric(row.Floor_local_fraction) && numeric(row.SD_multiplier) && Number(row.Floor_local_fraction) === 0 && Number(row.SD_multiplier) === 0) add(prefix + 'Floor_local_fraction', label + 'specify an amplitude or noise threshold.');
      if (legacy && numeric(row.Floor_local_fraction) && Number(row.Floor_local_fraction) <= 0) add(prefix + 'Floor_local_fraction', label + 'legacy mode requires a positive minimum rise.');
      if (numeric(row.Min_duration_s) && Number(row.Min_duration_s) <= 0) add(prefix + 'Min_duration_s', label + 'duration must be positive.');
      for (const [key, minimum] of [['Consecutive_samples', 1], ['Min_baseline_samples', 2]]) if (numeric(row[key]) && (!Number.isInteger(Number(row[key])) || Number(row[key]) < minimum)) add(prefix + key, label + key.replace(/_/g, ' ') + ' must be a whole number ≥ ' + minimum + '.');
      for (const key of ['Max_gap_s', 'Expected_interval_s']) if (row[key] !== null && row[key] !== '' && (!numeric(row[key]) || Number(row[key]) <= 0)) add(prefix + key, label + key.replace(/_/g, ' ') + ' must be positive or blank.');
      if (!draft.timestamps && numeric(value.sampling_hz) && Number(value.sampling_hz) > 0) {
        const interval = 1 / Number(value.sampling_hz), end = (n - 1) * interval;
        if (numeric(row.Start_s) && Number(row.Start_s) >= end) add(prefix + 'Start_s', label + 'response must start before the recording ends at ' + fmt(end) + ' s.');
        if (numeric(row.End_s) && Number(row.End_s) > end) warnings.push(label + 'response end exceeds this recording; complete-window rules may mark cells unscorable.');
        if (!legacy && row.Duration_mode === 'seconds' && numeric(row.Min_duration_s) && Number(row.Min_duration_s) < interval) warnings.push(label + 'minimum duration (' + row.Min_duration_s + ' s) is shorter than the frame interval (' + fmt(interval) + ' s). Timing is coarsely resolved.');
        if (!index && numeric(row.Stimulus_onset_s) && numeric(value.baseline_end) && (Number(value.baseline_end) - 1) * interval >= Number(row.Stimulus_onset_s)) warnings.push('The trace F₀ baseline includes the stimulus onset or later frames. Review the resting baseline.');
      }
    });
    if (!String(draft.protocol.name || '').trim()) add('protocol-name', 'Enter a protocol name.');
    if (draft.protocol.field_overrides !== null && !Array.isArray(draft.protocol.field_overrides)) add('field-overrides', 'Field overrides must be a valid JSON array.');
    const signature = JSON.stringify(draft);
    if (signature === S.draftAck && S.state && S.state.validation) {
      errors.push(...S.state.validation.errors); warnings.push(...S.state.validation.warnings);
    }
    for (const box of document.querySelectorAll('#side [aria-invalid]')) box.removeAttribute('aria-invalid');
    for (const id of fields) if ($(id)) $(id).setAttribute('aria-invalid', 'true');
    const box = $('validation'); box.hidden = !errors.length && !warnings.length; box.replaceChildren();
    for (const [messages, className] of [[errors, ''], [warnings, 'warning']]) if (messages.length) {
      const list = document.createElement('ul'); list.className = className;
      for (const message of new Set(messages)) { const li = document.createElement('li'); li.textContent = message; list.appendChild(li); }
      box.appendChild(list);
    }
    return { errors, warnings };
  }
  function qualityFor(name) {
    if (!S.results || !name) return [];
    const q = S.quality[name] || {}, row = callRow(name), flags = new Set();
    for (const text of [q.QC_flags, row.QC_flags]) for (const flag of String(text || '').split(/[;|,]/)) if (flag.trim()) flags.add(flag.trim());
    if ((S.results.traces[name] || []).some((v) => finite(v) && v < -1)) flags.add('TRACE_BELOW_MINUS_ONE');
    const metric = S.metrics[name] || {}, peak = number(metric.Peak_dFF0), noise = number(metric.Baseline_SD);
    if (finite(peak) && finite(noise) && peak > 0 && 0.1 * peak < noise && finite(number(metric.Rise_Time))) flags.add('RISE_START_BELOW_BASELINE_SD');
    return Array.from(flags);
  }
  const FLAG_TEXT = {
    DURATION_BELOW_SAMPLE_INTERVAL: 'Minimum response duration is shorter than one frame interval; timing is coarsely resolved',
    reference_below_3_baseline_sd_review: 'Low F₀ relative to baseline noise; review fractional amplitudes and the response call',
    TRACE_BELOW_MINUS_ONE: 'ΔF/F₀ falls below −1; corrected fluorescence becomes negative relative to its positive reference. Review background and baseline stability',
    RISE_START_BELOW_BASELINE_SD: 'The legacy 10% rise level is below baseline SD; the rise-time estimate is sensitive to baseline noise'
  };
  const flagText = (flag) => FLAG_TEXT[flag] || flag.replace(/_/g, ' ');
  function showQuality() {
    if (!S.results) return;
    const summary = $('quality-summary'), list = $('quality-list'), counts = new Map(); let affected = 0, responders = 0;
    for (const name of S.info.names) { const flags = qualityFor(name); if (flags.length) { affected++; if (callOf(name) === 'Responded') responders++; }
      for (const flag of flags) { const item = counts.get(flag) || { total: 0, responders: 0 }; item.total++; if (callOf(name) === 'Responded') item.responders++; counts.set(flag, item); }
    }
    const notices = [], window = selectedCall(), baseline = S.results.baselineFrames;
    if (baseline && S.results.time[baseline[1] - 1] >= number(window.Stimulus_onset_s || window.Start_s)) notices.push('The saved trace F₀ baseline includes the stimulus onset or later frames. These results keep their original settings; review timing before rerunning.');
    if (window.Window_warning) notices.push(window.Window_warning);
    $('quality-review').hidden = !affected && !notices.length;
    summary.textContent = affected ? 'Quality review: ' + affected + ' cells flagged, including ' + responders + ' responders' : 'Review the applied timing';
    list.replaceChildren();
    for (const [flag, count] of counts) { const li = document.createElement('li'); li.textContent = count.total + ' cells (' + count.responders + ' responders): ' + flagText(flag) + '.'; list.appendChild(li); }
    for (const message of new Set(notices)) { const li = document.createElement('li'); li.textContent = message; list.appendChild(li); }
    $('quality-filter').disabled = !affected;
  }
  function showMetricTable() {
    const key = $('metrics-table').value;
    $('metric-help').textContent = key === 'metrics' ? 'Legacy metrics: AUC is the signed net integral over the whole recording; negative drift can make it negative. Peak time is a recording timestamp. Rise time is the 10–90% crossing interval after the display baseline and can be noise-sensitive. For stimulus-specific values, choose Response-window metrics (including Peak_latency_s and Window_AUC_fraction_s).' : key === 'windowMetrics' ? 'Stimulus-window metrics use the local prestimulus reference. Peak_latency_s is relative to stimulus onset; Window_AUC_fraction_s covers the declared response window.' : '';
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
    S.state = state; return true;
  }
  function hasPending() { return !!S.results && (!!(S.state && S.state.settingsPending) || JSON.stringify(collectDraft()) !== S.draftAck); }
  async function start() {
    S.info = await api('session');
    $('recording').textContent = S.info.recording; $('recording').title = S.info.recording;
    document.title = S.info.names.length + ' cells · SimpleCa²⁺ · ' + S.info.recording + ' · ' + S.info.key.slice(0, 6); $('version').textContent = 'Calcium imaging plugin ' + S.info.pluginVersion;
    S.initialDraft = clone(S.info.draft);
    S.revision = Math.max(0, Number(S.info.draftRevision) || 0); fillDraft(S.info.draft); S.draftAck = JSON.stringify(collectDraft());
    $('time').max = Math.max(1, S.info.frames);
    for (const option of Array.from($('view').options)) if ((option.value === 'picture' && !S.info.hasPicture) || (option.value !== 'picture' && !S.info.hasMovie)) option.remove();
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
    const engine = $('engine'); engine.textContent = state.engineReady ? 'Analysis engine ready' : 'First run: set up analysis engine'; engine.className = 'pill ' + (state.engineReady ? 'good' : 'warn');
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
    else if (!S.results && !state.error && state.tone !== 'bad') { status.className = 'status'; status.textContent = state.engineReady ? 'Check the settings above, confirm the timing, then run. Movie playback works now.' : 'Movie playback is ready. Before the first analysis, choose “Set up analysis engine” below (requires internet).'; }
    else { status.className = 'status ' + (state.tone === 'good' ? 'good' : state.tone === 'bad' ? 'bad' : state.tone === 'warn' ? 'warn' : ''); status.textContent = state.message; }
    $('open-folder').disabled = !state.folder; $('open-folder').hidden = !state.folder;
    const pending = hasPending(), notice = $('applied-note'); notice.hidden = !S.results; notice.classList.toggle('pending', pending);
    notice.textContent = pending ? 'Settings changed. The results, movie timing and exports still use the last completed analysis. Choose “Apply changes & reanalyze” to update them.' : 'Showing the completed analysis. Response calls and counts come from the saved results.';
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
    const counts = r.counts || {}, whole = (key) => { const v = number(counts[key]); return isFinite(v) ? String(Math.round(v)) : '–'; }, percent = number(counts.Percent_responding);
    $('tile-responding').innerHTML = '<b>' + (isFinite(percent) ? Math.round(percent) + '%' : '–') + '</b><span>' + escape(counts.Stimulus || windows[0].stimulus) + ': ' + whole('Responded') + ' of ' + whole('Evaluable') + ' evaluable cells</span>';
    $('tile-registered').firstChild.textContent = whole('Registered'); $('tile-excluded').firstChild.textContent = whole('Excluded'); $('tile-unscorable').firstChild.textContent = whole('Unscorable');
    document.body.classList.add('has-results');
    for (const id of ['result-tiles', 'plot-card', 'group-label', 'window-select-label', 'every-label', 'call-legend']) $(id).hidden = false;
    for (const tab of document.querySelectorAll('.tab')) tab.disabled = false;
    const tables = $('metrics-table'); tables.innerHTML = '';
    const titles = { metrics: 'Legacy whole-recording metrics', quality: 'Cell quality', registry: 'Registered cells', windowMetrics: 'Response-window metrics', replicates: 'Biological-replicate counts', groups: 'Condition counts', unassigned: 'Unassigned biological units', cohort: 'Cohort flow', recordingMetadata: 'Recording details', sensitivityCalls: 'Sensitivity: cell calls', sensitivityCounts: 'Sensitivity: responder counts', negativeControlCalls: 'Negative controls: cell calls', negativeControlCounts: 'Negative controls: responder counts' };
    for (const key of Object.keys(r.tables || {})) if (!['calls', 'counts'].includes(key) && r.tables[key].header.length) { const o = document.createElement('option'); o.value = key; o.textContent = titles[key] || key.replace(/_/g, ' '); tables.appendChild(o); }
    if (r.tables.windowMetrics && r.tables.windowMetrics.rows.length) tables.value = 'windowMetrics'; else if (r.tables.metrics) tables.value = 'metrics'; showMetricTable();
    table('table-calls', r.tables.calls); table('table-counts', r.tables.counts);
    for (const id of ['save-plot', 'save-heatmap']) $(id).disabled = false;
    Plot.x0 = Plot.x1 = null; updateWindow(); setFrame(S.frame); showState();
  }
  function updateWindow() {
    const group = $('group'); group.innerHTML = '';
    GROUPS.forEach((name, i) => { const o = document.createElement('option'); o.value = i; o.textContent = name + ' (' + S.info.names.filter((n) => inGroup(callOf(n), i, n)).length + ')'; group.appendChild(o); });
    group.value = S.group; showQuality(); Plot.dirty = true; listCells(); drawHeatmap();
  }
  /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. */
  function listCells() {
    const names = S.info.names;
    S.listed = S.results ? names.filter((n) => inGroup(callOf(n), S.group, n)) : names.slice();
    const select = $('cell'); select.innerHTML = '';
    for (const name of S.listed) { const o = document.createElement('option'); o.value = name; o.textContent = name; select.appendChild(o); }
    $('previous').disabled = $('next').disabled = S.listed.length < 2;
    for (const tile of document.querySelectorAll('.tile')) tile.classList.toggle('on', !!S.results && ((tile.id === 'tile-responding' && S.group === 1) || (tile.id === 'tile-registered' && S.group === 0) || (tile.id === 'tile-excluded' && S.group === 3) || (tile.id === 'tile-unscorable' && S.group === 4)));
    pick(S.listed.includes(S.cell) ? S.cell : (S.listed[0] || null));
  }
  function showGroup(group) { if (S.group !== group) { S.group = group; $('group').value = group; listCells(); } }
  /** Selects a cell everywhere: the list, the picture, the plot, the heatmap and the tables. */
  function pick(name) {
    if (name && S.results && !S.listed.includes(name)) { S.group = 0; $('group').value = 0; S.cell = name; listCells(); return; }
    S.cell = name; if (name) $('cell').value = name;
    const verdict = $('verdict'), row = name && S.results ? calls()[name] : null;
    if (!name) verdict.textContent = S.results ? 'No cells in this group.' : '';
    else if (!row) verdict.innerHTML = S.results ? '<b>' + escape(name) + '</b> · no response call was saved for this cell' : '';
    else {
      const parts = ['<span class="dot" style="background:' + CALL[kind(row.Status)] + '"></span><b>' + escape(name) + ': ' + escape(row.Status || '') + '</b>', escape((row.Reason || '').toLowerCase())];
      const rise = number(row.Peak_increase), at = number(row.Peak_s), need = number(row.Threshold_increase);
      if (isFinite(rise)) parts.push('largest rise ' + (100 * rise).toFixed(1) + '%' + (isFinite(at) ? ' at ' + at.toFixed(1) + ' s' : '') + (isFinite(need) ? ' (needed ' + (100 * need).toFixed(1) + '%)' : ''));
      verdict.innerHTML = parts.filter(Boolean).join(' · ');
      verdict.title = 'Rise and threshold are shown as percentages of local prestimulus fluorescence. The saved call table stores the corresponding fractions.';
    }
    const flags = qualityFor(name); $('cell-quality').hidden = !flags.length; $('cell-quality').textContent = flags.map(flagText).join(' · ');
    $('zoom-cell').disabled = !name || !(S.info.hasPicture || S.info.hasMovie);
    for (const row of document.querySelectorAll('tbody tr')) row.classList.toggle('on', row.dataset.cell === name);
    if (Pic.zoom > 1 && name) centreOn(S.info.names.indexOf(name));
    Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); legend();
  }
  const escape = (text) => String(text).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  function stepCell(by) { if (!S.listed.length) return; const i = S.listed.indexOf(S.cell); pick(S.listed[((i + by) % S.listed.length + S.listed.length) % S.listed.length]); }

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
    const available = S.info.frames > 1 && (S.info.hasMovie || !!S.results), staticView = S.view === 'picture';
    $('time').disabled = !available || (staticView && !S.results);
    $('play').disabled = !available || staticView;
    $('frame-back').disabled = $('time').disabled || S.frame <= 1;
    $('frame-forward').disabled = $('time').disabled || S.frame >= S.info.frames;
    $('speed').disabled = $('play').disabled;
    const t = frameTimes()[S.frame - 1], provisional = !S.results && (!S.state || S.state.timingProvisional !== false);
    const label = (finite(t) ? (provisional ? '≈ ' : '') + t.toFixed(2) + ' s · ' : '') + 'Frame ' + S.frame + ' / ' + S.info.frames;
    $('time-label').textContent = label; $('time').setAttribute('aria-valuetext', label);
    $('playback-note').textContent = staticView ? 'Static colour overview — this image does not move. Choose “Recording” in View to play or scrub the movie.' : !S.info.hasMovie ? 'The recording is unavailable in this session. Reopen it in Fiji to view its frames.' : frameTimes().length ? (provisional ? 'Preview timing is provisional. Playback follows the displayed frame times; confirm timing before analysis.' : 'Playback follows the applied frame times. Drag the slider to pause and inspect a frame.') : 'Frame-by-frame preview. No valid time base is available; 1× plays 20 frames per second.';
  }
  function setFrame(wanted, manual) {
    if (!S.info) return;
    if (manual) play(false);
    S.frame = Math.max(1, Math.min(Math.max(1, S.info.frames), Math.round(wanted))); $('time').value = S.frame;
    updateTransport(); drawPlot(); requestFrame();
  }
  function seek(t) { setFrame(nearestFrame(t), true); }
  function play(on) {
    if (S.timer) clearTimeout(S.timer); S.timer = null;
    if (on && S.info && S.info.hasMovie && S.view !== 'picture' && S.info.frames > 1) {
      if (S.frame >= S.info.frames) setFrame(1);
      const schedule = () => {
        const times = frameTimes(), delta = times[S.frame] - times[S.frame - 1], delay = (finite(delta) && delta > 0 ? delta * 1000 : 50) / Number($('speed').value);
        S.timer = setTimeout(() => {
          const advance = () => {
            if (!S.timer) return;
            if (Pic.loading) { S.timer = setTimeout(advance, 20); return; }
            if (S.frame >= S.info.frames) { play(false); return; }
            setFrame(S.frame + 1); schedule();
          }; advance();
        }, Math.max(10, delay));
      }; schedule();
    }
    $('play').textContent = S.timer ? '❚❚ Pause' : '▶ Play'; $('play').setAttribute('aria-pressed', String(!!S.timer));
  }
  function movieMessage(message, error) { const box = $('movie-message'); box.textContent = message; box.hidden = !message; box.classList.toggle('error', !!error); }
  // ---------------------------------------------------------------- movie loading: latest requested frame wins, and failures stay visible
  function requestFrame() {
    if (!(S.info.hasPicture || S.info.hasMovie)) return;
    if (S.view === 'picture' || !S.info.hasMovie) {
      Pic.pending = null; Pic.shown = Pic.base; Pic.shownFrame = null;
      $('picture-note').textContent = 'Static overview · not a movie frame'; movieMessage(Pic.base ? '' : 'Loading the static overview…'); drawPicture(); return;
    }
    Pic.pending = { f: S.frame, mode: S.view };
    if (Pic.shownFrame !== S.frame || Pic.error) movieMessage('Loading frame ' + S.frame + '…');
    if (!Pic.loading) pump();
  }
  async function pump() {
    Pic.loading = true;
    try {
      while (Pic.pending) {
        const want = Pic.pending; Pic.pending = null;
        try {
          const response = await fetch('api/frame?f=' + want.f + '&mode=' + want.mode, { cache: 'no-store' });
          if (!response.ok) throw new Error('The recording frame could not be loaded.');
          const bitmap = await createImageBitmap(await response.blob());
          if (S.view === want.mode && S.frame === want.f) {
            if (Pic.shown && Pic.shown !== Pic.base && Pic.shown.close) Pic.shown.close();
            Pic.shown = bitmap; Pic.shownFrame = want.f; Pic.error = null; movieMessage('');
            const rest = S.results ? S.results.baselineFrames : (S.state && S.state.rest);
            $('picture-note').textContent = want.mode === 'change' ? 'Frame ' + want.f + ' − mean of baseline frames ' + (rest ? rest[0] + '–' + rest[1] : '(preview)') + ' · display smoothed' : 'Frame ' + want.f + ' · as recorded';
            drawPicture();
          } else if (bitmap.close) bitmap.close();
        } catch (problem) {
          if (S.view === want.mode && S.frame === want.f) {
            Pic.error = problem.message; play(false);
            movieMessage('Frame ' + want.f + ' could not load. Keep the recording open in Fiji, then move the slider to retry.', true);
            $('picture-note').textContent = Pic.shownFrame ? 'Last displayed frame: ' + Pic.shownFrame : 'No movie frame loaded';
          }
        }
      }
    } finally { Pic.loading = false; }
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
    const stroke = (path, colour, width, dash) => { ctx.setLineDash(dash ? [4 / at.z, 3 / at.z] : []); ctx.strokeStyle = colour; ctx.lineWidth = width / at.z; ctx.stroke(path); };
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
      const b = S.bounds[selected], name = S.cell; ctx.font = '600 12px ' + getComputedStyle(document.body).fontFamily;
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
  function legend() {
    const box = $('legend'), r = S.results; if (!r) { box.innerHTML = ''; return; }
    const chip = (colour, text) => '<span class="chip"><i style="background:' + colour + '"></i>' + escape(text) + '</span>';
    const parts = [chip('#0072B2', 'mean ± SEM')];
    if (S.cell) parts.push(r.traces[S.cell] ? chip('#1F2933', S.cell) : '<span class="chip" style="color:#9A5B00">' + escape(S.cell) + ': no trace (see Metrics + quality)</span>');
    if (S.every) {
      const drawn = S.listed.filter((n) => r.traces[n]), yes = drawn.filter((n) => callOf(n) === 'Responded').length, no = drawn.filter((n) => callOf(n) === 'Did not respond').length;
      if (yes) parts.push(chip(CALL.responded, 'responded (' + yes + ')')); if (no) parts.push(chip(CALL.not, 'did not respond (' + no + ')'));
      if (drawn.length - yes - no) parts.push(chip(CALL.other, 'no call (' + (drawn.length - yes - no) + ')'));
    }
    if (referenceLevels().length > 1) parts.push(chip('#177465', 'required rise · dashed'));
    box.innerHTML = parts.join('');
  }
  /** Everything except the time cursor and the hovered trace, drawn once and reused while the cursor moves. */
  function paintTraces(ctx, w, h) {
    const r = S.results, time = r.time, L = 54, T = 10, W = Math.max(10, w - L - 14), H = Math.max(10, h - T - 38);
    const full0 = time[0], full1 = time[time.length - 1], x0 = Plot.x0 === null ? full0 : Plot.x0, x1 = Plot.x1 === null ? full1 : Plot.x1;
    const drawn = S.every ? S.listed.filter((n) => r.traces[n]) : []; Plot.drawn = drawn;
    const spans = responseSpans(), gapLimit = samplingSupport().gap, levels = referenceLevels();
    const values = [0].concat(levels.map((l) => l.value)); r.mean.forEach((v, i) => { if (finite(v)) { values.push(v); if (finite(r.sem[i])) values.push(v - r.sem[i], v + r.sem[i]); } });
    const own = S.cell ? r.traces[S.cell] : null; if (own) for (const v of own) if (finite(v)) values.push(v);
    // with every cell drawn, the axis covers 99 % of all values (plus the mean and the selected cell); rarer extremes run off the edge
    if (drawn.length) { const rows = drawn.map((n) => r.traces[n]); values.push(share(rows, 0.005, false), share(rows, 0.995, false)); }
    let lo = Math.min.apply(null, values), hi = Math.max.apply(null, values); if (hi - lo < 0.001) hi = lo + 1; const pad = (hi - lo) * 0.08; lo -= pad; hi += pad;
    const X = (t) => L + (t - x0) / (x1 - x0) * W, Y = (v) => T + (hi - v) / (hi - lo) * H;
    Plot.geom = { L, T, W, H, x0, x1, lo, hi, full0, full1 };
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); ctx.font = '12px ' + getComputedStyle(document.body).fontFamily;
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W, H); ctx.clip();
    for (const [span, colour, name] of [[spans.baseline, '#EDF0F4', 'response baseline'], [spans.response, '#FFF4E3', 'response window'], [spans.run, '#D5EEDD', 'qualifying run']]) {
      if (!span || !finite(span[0]) || !finite(span[1])) continue;
      const a = Math.max(L, X(span[0])), b = Math.min(L + W, X(span[1])); if (b <= a) continue;
      ctx.fillStyle = colour; ctx.fillRect(a, T, b - a, H);
      if (ctx.measureText(name).width + 8 <= b - a) { ctx.fillStyle = '#65717D'; ctx.fillText(name, a + 4, T + 14); }
    }
    ctx.restore();
    const ys = ticks(lo, hi, 5), xs = ticks(x0, x1, Math.max(2, Math.floor(W / 90)));
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const v of ys.values) { const y = Math.round(Y(v)) + 0.5; ctx.strokeStyle = v === 0 ? '#AEB8C2' : '#DEE4EA'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + W, y); ctx.stroke(); ctx.fillStyle = '#65717D'; ctx.fillText(tickText(v, ys.step), L - 8, y); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const v of xs.values) { const x = Math.round(X(v)) + 0.5; ctx.strokeStyle = '#DEE4EA'; ctx.beginPath(); ctx.moveTo(x, T + H); ctx.lineTo(x, T + H + 4); ctx.stroke(); ctx.fillStyle = '#65717D'; ctx.fillText(tickText(v, xs.step), x, T + H + 7); }
    ctx.fillText('Time (s)', L + W / 2, h - 15);
    ctx.save(); ctx.translate(13, T + H / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'middle'; ctx.fillText('ΔF/F₀', 0, 0); ctx.restore();
    ctx.strokeStyle = '#C5CED8'; ctx.strokeRect(L + 0.5, T + 0.5, W, H);
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W + 1, H + 1); ctx.clip(); ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    const line = (y) => { ctx.beginPath(); let pen = false; for (let i = 0; i < time.length && i < y.length; i++) { if (finite(y[i])) { const px = X(time[i]), py = Y(y[i]); if (pen && time[i] - time[i - 1] <= gapLimit) ctx.lineTo(px, py); else ctx.moveTo(px, py); pen = true; } else pen = false; } ctx.stroke(); };
    if (drawn.length) {
      ctx.globalAlpha = drawn.length > 120 ? 0.22 : drawn.length > 40 ? 0.32 : 0.48; ctx.lineWidth = 0.9;
      // non-responders first, so the responders are not hidden under them
      for (const responders of [false, true]) for (const name of drawn) { const call = callOf(name); if ((call === 'Responded') === responders) { ctx.strokeStyle = CALL[kind(call)]; line(r.traces[name]); } }
      ctx.globalAlpha = 1;
    }
    ctx.fillStyle = drawn.length ? 'rgba(0,114,178,0.28)' : 'rgba(0,114,178,0.20)';
    for (let i = 1; i < time.length; i++) if (time[i] - time[i - 1] <= gapLimit && [r.mean[i - 1], r.mean[i], r.sem[i - 1], r.sem[i]].every(finite)) {
      ctx.beginPath(); ctx.moveTo(X(time[i - 1]), Y(r.mean[i - 1] - r.sem[i - 1])); ctx.lineTo(X(time[i]), Y(r.mean[i] - r.sem[i])); ctx.lineTo(X(time[i]), Y(r.mean[i] + r.sem[i])); ctx.lineTo(X(time[i - 1]), Y(r.mean[i - 1] + r.sem[i - 1])); ctx.closePath(); ctx.fill();
    }
    for (const level of levels) { ctx.strokeStyle = level.colour; ctx.lineWidth = level.name === 'required rise' ? 1.8 : 1; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(L, Y(level.value)); ctx.lineTo(L + W, Y(level.value)); ctx.stroke(); } ctx.setLineDash([]);
    if (drawn.length) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 4.2; line(r.mean); }
    ctx.strokeStyle = '#0072B2'; ctx.lineWidth = 2.4; line(r.mean);
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
    if (Plot.hover && r.traces[Plot.hover] && Plot.hover !== S.cell) {
      const y = r.traces[Plot.hover], gapLimit = samplingSupport().gap; ctx.strokeStyle = CALL[kind(callOf(Plot.hover))]; ctx.lineWidth = 2.2; ctx.lineJoin = 'round'; ctx.beginPath(); let pen = false;
      for (let i = 0; i < time.length; i++) { if (finite(y[i])) { if (pen && time[i] - time[i - 1] <= gapLimit) ctx.lineTo(X(time[i]), Y(y[i])); else ctx.moveTo(X(time[i]), Y(y[i])); pen = true; } else pen = false; } ctx.stroke();
    }
    const t = S.frame <= time.length ? time[S.frame - 1] : NaN;
    if (isFinite(t) && t >= g.x0 && t <= g.x1) {
      const x = Math.round(X(t)) + 0.5; ctx.strokeStyle = 'rgba(31,41,51,0.8)'; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(x, g.T); ctx.lineTo(x, g.T + g.H); ctx.stroke();
      const own = S.cell ? r.traces[S.cell] : null; let text = t.toFixed(1) + ' s';
      if (own && finite(own[S.frame - 1])) { const y = Math.max(g.T, Math.min(g.T + g.H, Y(own[S.frame - 1]))); ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(x, y, 5, 0, 7); ctx.fill(); ctx.fillStyle = '#1F2933'; ctx.beginPath(); ctx.arc(x, y, 3.2, 0, 7); ctx.fill(); text += ' · ' + S.cell + ' ' + own[S.frame - 1].toFixed(3); }
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
    const candidates = Plot.drawn.slice(); if (S.cell && r.traces[S.cell] && !candidates.includes(S.cell)) candidates.push(S.cell);
    let best = null, near = 7;
    for (const name of candidates) { const v = r.traces[name][i]; if (!finite(v)) continue; const d = Math.abs(g.T + (g.hi - v) / (g.hi - g.lo) * g.H - my); if (d < near) { near = d; best = name; } }
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
    if (tall) { ctx.font = '11px ' + family; names.forEach((name, c) => { ctx.fillStyle = name === S.cell ? '#1F2933' : '#65717D'; ctx.fillText(name, LEFT - 13, TOP + (c + 0.5) * rowH); }); }
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
    const labels = { Time_to_Peak: 'Peak timestamp from recording start (s)', AUC: 'Whole-recording signed AUC', Rise_Time: 'Legacy 10–90% rise time (s)' };
    const cellColumn = data.header.indexOf('Cell'), numeric = data.header.map((_, c) => data.rows.some((row) => row[c] !== '' && row[c] !== undefined) && data.rows.every((row) => row[c] === '' || row[c] === undefined || isFinite(Number(row[c]))));
    const head = element.createTHead().insertRow(), body = element.createTBody(); let sorted = -1, direction = 1;
    const fill = (rows) => { body.innerHTML = '';
      if (!rows.length) { const td = body.insertRow().insertCell(); td.colSpan = Math.max(1, columns.length); td.className = 'empty-table'; td.textContent = 'No rows were saved for this table.'; return; }
      for (const row of rows) { const tr = body.insertRow(); if (cellColumn >= 0) { tr.dataset.cell = row[cellColumn]; tr.classList.toggle('on', row[cellColumn] === S.cell); }
        for (const { name, index: c } of columns) { const value = row[c] == null ? '' : row[c], td = tr.insertCell();
          if (numeric[c] && value !== '') { const v = Number(value); td.className = 'num'; td.textContent = Number.isInteger(v) ? String(v) : String(+v.toPrecision(4)); td.title = value; } else td.textContent = value;
          if (name === 'Cell') td.classList.add('cell-column');
        }
      }
    };
    columns.forEach(({ name, index: c }) => { const th = document.createElement('th'); th.textContent = labels[name] || name; th.title = 'Saved column: ' + name; if (name === 'Cell') th.classList.add('cell-column'); head.appendChild(th);
      th.addEventListener('click', () => { direction = sorted === c ? -direction : 1; sorted = c; for (const other of head.children) other.classList.remove('up', 'down'); th.classList.add(direction > 0 ? 'up' : 'down');
        fill(data.rows.slice().sort((a, b) => direction * (numeric[c] ? (Number(a[c] === '' ? NaN : a[c]) || 0) - (Number(b[c] === '' ? NaN : b[c]) || 0) : String(a[c]).localeCompare(String(b[c]))))); }); });
    body.addEventListener('click', (event) => { const tr = event.target.closest('tr'); if (tr && tr.dataset.cell) pick(tr.dataset.cell); });
    fill(data.rows);
  }

  function showTab(name) {
    for (const tab of document.querySelectorAll('.tab')) tab.classList.toggle('on', tab.dataset.tab === name);
    for (const panel of document.querySelectorAll('.panel')) panel.hidden = panel.id !== 'panel-' + name;
    if (name === 'heatmap') drawHeatmap(); if (name === 'review') { Plot.dirty = true; drawPlot(); drawPicture(); }
    const on = document.querySelector('tbody tr.on'); if (on && name !== 'review' && name !== 'heatmap') on.scrollIntoView({ block: 'nearest' });
  }
  function download(blob, name) { const link = document.createElement('a'); link.href = URL.createObjectURL(blob); link.download = name; link.click(); setTimeout(() => URL.revokeObjectURL(link.href), 2000); }
  /** A file name that says what it is: the recording, what was saved, the cell when it matters, and when. */
  function exportName(kind, extension) {
    const clean = (text) => String(text || '').replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 80);
    const now = new Date(), two = (n) => String(n).padStart(2, '0');
    const stamp = now.getFullYear() + two(now.getMonth() + 1) + two(now.getDate()) + '-' + two(now.getHours()) + two(now.getMinutes()) + two(now.getSeconds());
    return [clean(S.info.recording), kind, kind === 'traces' && S.cell ? clean(S.cell) : '', stamp].filter(Boolean).join('_') + '.' + extension;
  }
  function exportMetadata() {
    return { schema: 'simpleca.view-export', schema_version: 1, recording: S.info.recording, plugin_version: S.info.pluginVersion,
      appliedFingerprint: S.results.appliedFingerprint || S.state.appliedFingerprint, appliedProtocol: S.results.appliedProtocol || S.results.protocol,
      manifest: S.results.manifest, settingsPending: hasPending(), displayed_window: S.window, cell: S.cell, group: GROUPS[S.group],
      frame: S.frame, time_s: frameTimes()[S.frame - 1], view: S.view, results_folder: S.results.folder || (S.state && S.state.folder) || null,
      time_range: Plot.geom ? [Plot.geom.x0, Plot.geom.x1] : [Plot.x0, Plot.x1], time_zoomed: Plot.x0 !== null, palette: S.palette, range: S.range, order: S.order };
  }
  function save(canvas, kind, redraw) {
    // Capture provenance now, not after the asynchronous PNG encoder has finished.
    const metadata = JSON.stringify(exportMetadata(), null, 2), name = exportName(kind, 'png');
    // Draw once more at least 1600 pixels wide (and at least twice the screen's density), save that, then draw for the screen again.
    S.exportScale = Math.min(6, Math.max(2, Math.ceil(1600 / Math.max(1, canvas.clientWidth)))); Plot.dirty = true; redraw();
    const restore = () => { S.exportScale = null; Plot.dirty = true; redraw(); };
    try { canvas.toBlob((blob) => { download(blob, name); download(new Blob([metadata], { type: 'application/json' }), name.replace(/\.png$/, '.metadata.json')); }, 'image/png'); }
    finally { restore(); }
  }
  function wire() {
    for (const tab of document.querySelectorAll('.tab')) tab.addEventListener('click', () => showTab(tab.dataset.tab));
    for (const tile of document.querySelectorAll('.tile')) tile.addEventListener('click', () => {
      if (!S.results) return;
      const primary = (S.results.windows || []).find((w) => w.role === 'Stimulus');
      if (primary && S.window !== primary.key) { S.window = primary.key; $('response-window').value = S.window; updateWindow(); }
      showGroup(Number(tile.dataset.group)); showTab('review');
    });
    $('toggle-settings').addEventListener('click', () => { const hidden = document.body.classList.toggle('no-side'); $('toggle-settings').textContent = hidden ? 'Show settings' : 'Hide settings'; redraw(); });
    $('side').addEventListener('input', (event) => {
      const box = event.target;
      if (!S.draft || box.type === 'file' || ['confirmed', 'control', 'protocol-mode'].includes(box.id)) return;
      if (box.dataset.shown !== undefined) box.classList.toggle('changed', box.value !== box.dataset.shown);
      if ($('protocol-mode').value === '1' && box.dataset.key === 'Start_s') for (const key of ['Baseline_end_s', 'Stimulus_onset_s']) exactInput($('window-' + box.dataset.window + '-' + key), readInput(box, false));
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
        try { if (file.size > 8 * 1024 * 1024) throw new Error('Choose a CSV smaller than 8 MB.'); const content = await file.text(); S.draft[key] = { name: file.name, content }; fileNote(key); draftChanged(); }
        catch (error) { draftError(error); }
      });
      $(key + '-clear').addEventListener('click', () => { S.draft[key] = null; $(key + '-file').value = ''; fileNote(key); draftChanged(); });
    }
    $('protocol-import').addEventListener('change', async () => {
      const file = $('protocol-import').files[0]; if (!file) return;
      try { let parsed; try { parsed = JSON.parse(await file.text()); } catch (bad) { throw new Error(file.name + ' is not a valid JSON file (' + bad.message + '). Choose a protocol exported from SimpleCa²⁺ or this page.'); }
        const result = await api('protocol', { protocol: parsed }); const draft = collectDraft(); draft.protocol = result.protocol; S.protocolDrafts = {}; fillDraft(draft); $('confirmed').checked = false; draftChanged(); }
      catch (error) { draftError(error); }
      finally { $('protocol-import').value = ''; }
    });
    $('protocol-export').addEventListener('click', async () => {
      try { const result = await api('protocol', { protocol: collectDraft().protocol }); download(new Blob([JSON.stringify(result.protocol, null, 2)], { type: 'application/json' }), exportName('response-protocol', 'json')); }
      catch (error) { draftError(error); }
    });
    $('analyze').addEventListener('click', async () => {
      if (validateDraft(collectDraft()).errors.length) { $('validation').scrollIntoView({ block: 'nearest' }); return; }
      if (!$('confirmed').checked) { S.note = 'Confirm the frame rate or timestamps, baseline and stimulus timing using the checkbox above.'; showState(); $('confirmed').focus(); return; }
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
    $('close-page').addEventListener('click', async () => {
      try { await api('close', {}); S.closed = true; play(false); $('lost').textContent = 'This analysis page is closed. Saved results remain in the results folder. You can close this tab.'; $('lost').classList.add('on'); showState(); }
      catch (error) { draftError(error); }
    });
    $('discard').addEventListener('click', () => { S.protocolDrafts = {}; fillDraft((S.results && S.results.appliedDraft) || S.initialDraft); $('confirmed').checked = false; for (const key of FILES) $(key + '-file').value = ''; draftChanged(); });
    $('quality-filter').addEventListener('click', () => { showGroup(5); showTab('review'); });
    $('all-call-columns').addEventListener('change', () => table('table-calls', S.results.tables.calls));
    $('response-window').addEventListener('change', () => { S.window = $('response-window').value; updateWindow(); });
    $('metrics-table').addEventListener('change', showMetricTable);
    $('group').addEventListener('change', () => { S.group = Number($('group').value); listCells(); });
    $('cell').addEventListener('change', () => pick($('cell').value));
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
    $('save-plot').addEventListener('click', () => save($('plot'), 'traces', drawPlot)); $('save-heatmap').addEventListener('click', () => save($('heatmap'), 'heatmap', drawHeatmap));

    const picture = $('picture');
    picture.addEventListener('wheel', (event) => { event.preventDefault(); const box = picture.getBoundingClientRect(), mx = event.clientX - box.left, my = event.clientY - box.top, before = place(picture.clientWidth, picture.clientHeight);
      const ix = (mx - before.x) / before.z, iy = (my - before.y) / before.z; Pic.zoom = Math.max(1, Math.min(16, Pic.zoom * Math.pow(1.0015, -event.deltaY)));
      if (Pic.zoom === 1) Pic.panX = Pic.panY = 0; else { const after = place(picture.clientWidth, picture.clientHeight); Pic.panX += mx - (after.x + ix * after.z); Pic.panY += my - (after.y + iy * after.z); } drawPicture(); }, { passive: false });
    picture.addEventListener('mousedown', (event) => { Pic.drag = { x: event.clientX, y: event.clientY }; Pic.moved = false; });
    window.addEventListener('mousemove', (event) => { if (!Pic.drag) return; if (Math.abs(event.clientX - Pic.drag.x) + Math.abs(event.clientY - Pic.drag.y) > 3) Pic.moved = true;
      if (Pic.moved) { Pic.panX += event.clientX - Pic.drag.x; Pic.panY += event.clientY - Pic.drag.y; Pic.drag = { x: event.clientX, y: event.clientY }; drawPicture(); } });
    window.addEventListener('mouseup', () => { Pic.drag = null; Plot.scrubbing = false; });
    picture.addEventListener('click', (event) => { if (Pic.moved) return; const i = cellAt(event); if (i >= 0) pick(S.info.names[i]); });
    picture.addEventListener('dblclick', () => { Pic.zoom = 1; Pic.panX = Pic.panY = 0; drawPicture(); });
    picture.addEventListener('mousemove', (event) => { if (Pic.drag) return; const i = cellAt(event); tip(i >= 0 ? S.info.names[i] + (callOf(S.info.names[i]) ? ' · ' + callOf(S.info.names[i]) : '') : '', event.clientX, event.clientY); });
    picture.addEventListener('mouseleave', () => tip(''));

    const plot = $('plot');
    plot.addEventListener('pointermove', (event) => {
      if (Plot.scrubbing) {
        const g = Plot.geom, box = plot.getBoundingClientRect(); if (g) seek(g.x0 + Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W)) * (g.x1 - g.x0)); tip(''); return;
      }
      const at = traceAt(event); if (!at) { tip(''); if (Plot.hover) { Plot.hover = null; drawPlot(); } return; }
      if (at.name !== Plot.hover) { Plot.hover = at.name; drawPlot(); }
      tip(at.name ? at.name + ' · ' + S.results.traces[at.name][at.index].toFixed(3) + ' at ' + S.results.time[at.index].toFixed(1) + ' s' + (callOf(at.name) ? ' · ' + callOf(at.name) : '') : '', event.clientX, event.clientY);
    });
    plot.addEventListener('pointerleave', () => { tip(''); if (Plot.hover) { Plot.hover = null; drawPlot(); } });
    plot.addEventListener('pointerdown', (event) => { const at = traceAt(event); if (!at) return; event.preventDefault(); if (at.name && at.name !== S.cell) pick(at.name); Plot.scrubbing = true; plot.setPointerCapture(event.pointerId); seek(at.time); });
    for (const type of ['pointerup', 'pointercancel', 'lostpointercapture']) plot.addEventListener(type, () => { Plot.scrubbing = false; });
    plot.addEventListener('wheel', (event) => { const g = Plot.geom; if (!g) return; event.preventDefault(); const box = plot.getBoundingClientRect(), share01 = Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W));
      const centre = g.x0 + share01 * (g.x1 - g.x0), factor = Math.pow(1.0015, event.deltaY), span = Math.min(g.full1 - g.full0, Math.max((g.full1 - g.full0) / 60, (g.x1 - g.x0) * factor));
      let a = centre - share01 * span, b = a + span; if (a < g.full0) { b += g.full0 - a; a = g.full0; } if (b > g.full1) { a -= b - g.full1; b = g.full1; }
      Plot.x0 = Math.max(g.full0, a); Plot.x1 = Math.min(g.full1, b); if (Plot.x1 - Plot.x0 >= g.full1 - g.full0 - 1e-9) Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); }, { passive: false });
    plot.addEventListener('dblclick', () => { Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); });

    const heatmap = $('heatmap');
    heatmap.addEventListener('click', (event) => { const at = heatmapAt(event); if (at) { pick(at.name); if (at.inside) setFrame(at.frame + 1, true); showTab('review'); } });
    heatmap.addEventListener('mousemove', (event) => { const at = heatmapAt(event); if (!at) { tip(''); return; } const trace = S.results.traces[at.name], v = trace ? trace[at.frame] : null;
      tip(at.name + (callOf(at.name) ? ' · ' + callOf(at.name) : '') + (at.inside ? ' · ' + (at.gap ? 'acquisition gap; nearest frame' : finite(v) ? v.toFixed(3) : 'no value') + ' at ' + S.results.time[at.frame].toFixed(1) + ' s' : ''), event.clientX, event.clientY); });
    heatmap.addEventListener('mouseleave', () => tip(''));

    document.addEventListener('keydown', (event) => { const tag = (event.target.tagName || '').toLowerCase(); if (tag === 'input' || tag === 'select' || tag === 'textarea' || $('panel-review').hidden) return;
      if (event.key === 'ArrowLeft') { stepCell(-1); event.preventDefault(); } else if (event.key === 'ArrowRight') { stepCell(1); event.preventDefault(); } else if (event.key === ' ' && S.info && S.info.hasMovie && S.view !== 'picture' && tag !== 'button') { play(!S.timer); event.preventDefault(); } });
    window.addEventListener('resize', redraw);
    if (window.ResizeObserver) { new ResizeObserver(redraw).observe($('panel-review')); new ResizeObserver(() => drawHeatmap()).observe($('heatmap-box')); }
  }
  function redraw() { Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); }

  wire();
  start().catch((problem) => { $('lost').classList.add('on'); $('status').textContent = problem.message; });
  window.calciumStep4 = { S, Plot, Pic, Heat, pick, setFrame, showGroup, showTab, collectDraft, frameTimes, nearestFrame, samplingSupport, referenceLevels, exportMetadata, heatmapAt };   // for the browser test
})();

// Step 4 page. Everything shown here is read from the tables the analysis saved, or from the recording's own frames.
// The page draws; it does not calculate results. Settings go back to the plugin, which runs the same engine as the Fiji window.
'use strict';
(function () {
  const $ = (id) => document.getElementById(id);
  const CALL = { responded: '#177465', not: '#925223', other: '#8A94A0' };
  const ON_PICTURE = { responded: '#3ddbb5', not: '#ffa24a', other: '#aeb8c2' };
  const GROUPS = ['All cells', 'Responded', 'Did not respond', 'Excluded or unscorable'];
  const kind = (call) => call === 'Responded' ? 'responded' : call === 'Did not respond' ? 'not' : 'other';
  const inGroup = (call, group) => group === 0 || (group === 1 && call === 'Responded') || (group === 2 && call === 'Did not respond') ||
    (group === 3 && call !== 'Responded' && call !== 'Did not respond');
  const S = { info: null, state: null, results: null, version: -1, group: 0, cell: null, listed: [], frame: 1, view: 'picture', viewChosen: false,
    every: true, outlines: true, order: 'call', palette: 'signed', range: 'robust', tall: false, timer: null, lost: 0, note: null, polling: false, paths: [], boxes: [] };
  const Plot = { dirty: true, layer: document.createElement('canvas'), x0: null, x1: null, hover: null, geom: null, scrubbing: false, drawn: [] };
  const Pic = { base: null, shown: null, zoom: 1, panX: 0, panY: 0, pending: null, loading: false, drag: null, moved: false };
  const Heat = { geom: null, names: [] };

  async function api(path, body) {
    const options = body === undefined ? { cache: 'no-store' } : { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) };
    const response = await fetch('api/' + path, options);
    const data = await response.json().catch(() => ({}));
    if (!response.ok) throw new Error(data.error || 'The plugin answered ' + response.status);
    return data;
  }
  const finite = (v) => typeof v === 'number' && isFinite(v);
  const callOf = (name) => S.results && S.results.calls[name] ? S.results.calls[name].Status : undefined;
  const number = (text) => { const v = parseFloat(text); return isFinite(v) ? v : NaN; };

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

  // ---------------------------------------------------------------- start, and keeping in step with the plugin
  async function start() {
    S.info = await api('session');
    $('recording').textContent = S.info.recording; $('recording').title = S.info.recording;
    document.title = 'SimpleCa²⁺ · ' + S.info.recording;
    $('version').textContent = 'Calcium imaging plugin ' + S.info.pluginVersion;
    for (const setting of S.info.settings) {
      const box = $(setting.key); if (!box) continue;
      S.boxes.push(box); box.value = setting.shown; box.dataset.shown = setting.shown;
      if (setting.shown !== setting.exact) box.title = 'Exact value: ' + setting.exact + '. It is what the analysis receives unless you edit this box.';
      box.addEventListener('input', () => box.classList.toggle('changed', box.value.trim() !== box.dataset.shown));
    }
    $('time').max = Math.max(2, S.info.frames);
    if (S.info.hasMovie) for (const [value, text] of [['change', 'Movie: change from rest'], ['raw', 'Movie: as recorded']]) { const o = document.createElement('option'); o.value = value; o.textContent = text; $('view').appendChild(o); }
    S.paths = S.info.rois.map((points) => { const p = new Path2D(); points.forEach((q, i) => i ? p.lineTo(q[0], q[1]) : p.moveTo(q[0], q[1])); p.closePath(); return p; });
    S.bounds = S.info.rois.map((points) => { let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9; for (const q of points) { x0 = Math.min(x0, q[0]); y0 = Math.min(y0, q[1]); x1 = Math.max(x1, q[0]); y1 = Math.max(y1, q[1]); } return { x0, y0, x1, y1 }; });
    if (S.info.hasPicture) { const image = new Image(); image.onload = () => { Pic.base = image; if (S.view === 'picture') Pic.shown = image; drawPicture(); }; image.src = 'api/picture'; }
    else $('picture-card').hidden = true;
    listCells(); setFrame(1); await poll(); setInterval(poll, 600);
  }
  async function poll() {
    if (S.polling) return; S.polling = true;
    try {
      const state = await api('state'); S.lost = 0; $('lost').classList.remove('on');
      if (S.info && state.key !== S.info.key) { location.reload(); return; }   // measured again in Fiji: start over with the new cells
      const wasBusy = S.state && S.state.busy; S.state = state;
      if (wasBusy && !state.busy) S.note = null;
      showState();
      if (state.version > 0 && state.version !== S.version) { S.version = state.version; await loadResults(); }
    } catch (problem) { if (++S.lost >= 3) $('lost').classList.add('on'); }
    finally { S.polling = false; }
  }
  function showState() {
    const state = S.state; if (!state) return;
    const engine = $('engine'); engine.textContent = state.engineReady ? 'Analysis engine ready' : 'Analysis engine not set up'; engine.className = 'pill ' + (state.engineReady ? 'good' : 'warn');
    $('analyze').disabled = state.busy; $('setup').disabled = state.busy; $('cancel').disabled = !state.busy; $('log').disabled = !state.logAvailable;
    for (const box of S.boxes) box.disabled = state.busy; $('confirmed').disabled = state.busy; $('kcl').disabled = state.busy;
    const progress = $('progress'); progress.hidden = !state.busy && !(state.percent > 0); progress.classList.toggle('busy', state.busy && state.installing && !(state.percent > 0));
    $('bar').style.width = (state.busy ? Math.max(3, state.percent) : state.percent) + '%';
    const status = $('status');
    if (state.busy) { status.className = 'status'; status.textContent = (state.percent > 0 ? state.percent + ' % · ' : '') + state.stage + ' · ' + state.elapsed + ' s elapsed'; }
    else if (S.note) { status.className = 'status warn'; status.textContent = S.note; }
    else { status.className = 'status ' + (state.tone === 'good' ? 'good' : state.tone === 'bad' ? 'bad' : ''); status.textContent = state.message; }
    $('open-folder').disabled = !state.folder;
  }
  async function send(path, body) {
    S.note = null;
    try { S.state = await api(path, body || {}); showState(); }
    catch (problem) { S.note = problem.message; showState(); }
  }

  // ---------------------------------------------------------------- results
  async function loadResults() {
    S.results = await api('results');
    const r = S.results, names = S.info.names, counts = r.counts || {};
    const whole = (key) => { const v = number(counts[key]); return isFinite(v) ? String(Math.round(v)) : '–'; };
    const percent = number(counts.Percent_responding);
    $('tile-responding').innerHTML = '<b>' + (isFinite(percent) ? Math.round(percent) + ' %' : '–') + '</b><span>responding: ' + whole('Responded') + ' of ' + whole('Evaluable') + ' evaluable cells</span>';
    $('tile-registered').firstChild.textContent = whole('Registered'); $('tile-excluded').firstChild.textContent = whole('Excluded'); $('tile-unscorable').firstChild.textContent = whole('Unscorable');
    const group = $('group'); group.innerHTML = '';
    GROUPS.forEach((name, i) => { const o = document.createElement('option'); o.value = i; o.textContent = name + ' (' + names.filter((n) => inGroup(callOf(n), i)).length + ')'; group.appendChild(o); });
    group.value = S.group;
    table('table-metrics', r.tables.metrics); table('table-calls', r.tables.calls); table('table-counts', r.tables.counts);
    Plot.x0 = Plot.x1 = null;
    if (S.info.hasMovie && !S.viewChosen) { S.view = 'change'; $('view').value = 'change'; }
    for (const id of ['save-plot', 'save-heatmap', 'play']) $(id).disabled = false;
    listCells();
    let peak = 1, top = -Infinity; r.mean.forEach((v, i) => { if (finite(v) && v > top) { top = v; peak = i + 1; } });
    setFrame(peak); drawHeatmap();
  }
  /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. */
  function listCells() {
    const names = S.info.names;
    S.listed = S.results ? names.filter((n) => inGroup(callOf(n), S.group)) : [];
    const select = $('cell'); select.innerHTML = '';
    for (const name of S.listed) { const o = document.createElement('option'); o.value = name; o.textContent = name; select.appendChild(o); }
    $('previous').disabled = $('next').disabled = S.listed.length < 2;
    for (const tile of document.querySelectorAll('.tile')) tile.classList.toggle('on', !!S.results && ((tile.id === 'tile-responding' && S.group === 1) || (tile.id === 'tile-registered' && S.group === 0) || ((tile.id === 'tile-excluded' || tile.id === 'tile-unscorable') && S.group === 3)));
    pick(S.listed.includes(S.cell) ? S.cell : (S.listed[0] || null));
  }
  function showGroup(group) { if (S.group !== group) { S.group = group; $('group').value = group; listCells(); } }
  /** Selects a cell everywhere: the list, the picture, the plot, the heatmap and the tables. */
  function pick(name) {
    if (name && S.results && !S.listed.includes(name)) { S.group = 0; $('group').value = 0; S.cell = name; listCells(); return; }
    S.cell = name; if (name) $('cell').value = name;
    const verdict = $('verdict'), row = name && S.results ? S.results.calls[name] : null;
    if (!name) verdict.textContent = S.results ? 'No cells in this group.' : '';
    else if (!row) verdict.innerHTML = S.results ? '<b>' + name + '</b> · no response call was saved for this cell' : '';
    else {
      const parts = ['<span class="dot" style="background:' + CALL[kind(row.Status)] + '"></span><b>' + name + ': ' + escape(row.Status || '') + '</b>', escape((row.Reason || '').toLowerCase())];
      const rise = number(row.Peak_increase), at = number(row.Peak_s), need = number(row.Threshold_increase);
      if (isFinite(rise)) parts.push('largest rise ' + rise.toFixed(2) + (isFinite(at) ? ' at ' + at.toFixed(0) + ' s' : '') + (isFinite(need) ? ' (needed ' + need.toFixed(2) + ')' : ''));
      verdict.innerHTML = parts.filter(Boolean).join(' · ');
      verdict.title = 'Rise and threshold are fractions of the fluorescence just before the stimulus, as in the Cell response calls table.';
    }
    $('zoom-cell').disabled = !name || !S.info.hasPicture;
    for (const row of document.querySelectorAll('tbody tr')) row.classList.toggle('on', row.dataset.cell === name);
    if (Pic.zoom > 1 && name) centreOn(S.info.names.indexOf(name));
    Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); legend();
  }
  const escape = (text) => String(text).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  function stepCell(by) { if (!S.listed.length) return; const i = S.listed.indexOf(S.cell); pick(S.listed[((i + by) % S.listed.length + S.listed.length) % S.listed.length]); }

  // ---------------------------------------------------------------- time
  function setFrame(wanted) {
    const most = Math.max(2, S.info.frames); S.frame = Math.max(1, Math.min(most, wanted)); $('time').value = S.frame;
    const t = S.results && S.frame <= S.results.time.length ? S.results.time[S.frame - 1] : NaN;
    $('time-label').textContent = (isFinite(t) ? t.toFixed(1) + ' s · ' : '') + 'frame ' + S.frame + ' of ' + S.info.frames;
    drawPlot(); requestFrame();
  }
  function seek(t) { const time = S.results.time; let best = 0, gap = Infinity; for (let i = 0; i < time.length; i++) { const d = Math.abs(time[i] - t); if (d < gap) { gap = d; best = i; } } setFrame(best + 1); }
  function play(on) {
    if (S.timer) { clearInterval(S.timer); S.timer = null; }
    if (on) { if (S.frame >= S.info.frames) setFrame(1); S.timer = setInterval(() => { if (S.frame >= S.info.frames) play(false); else setFrame(S.frame + 1); }, 50); }
    $('play').textContent = S.timer ? '❚❚ Pause' : '▶ Play';
  }

  // ---------------------------------------------------------------- the picture: the colour picture, or the recording at the cursor's moment
  function requestFrame() {
    if (!S.info.hasPicture) return;
    if (S.view === 'picture' || !S.info.hasMovie) { Pic.shown = Pic.base; $('picture-note').textContent = S.info.hasMovie ? 'colour picture (does not change with time)' : ''; drawPicture(); return; }
    Pic.pending = { f: S.frame, mode: S.view }; if (!Pic.loading) pump();
  }
  /** Fetches the newest wanted frame; frames asked for while one is on its way are skipped, so playing never queues up. */
  async function pump() {
    Pic.loading = true;
    while (Pic.pending) {
      const want = Pic.pending; Pic.pending = null;
      try {
        const response = await fetch('api/frame?f=' + want.f + '&mode=' + want.mode, { cache: 'no-store' });
        if (!response.ok) throw new Error('frame');
        const bitmap = await createImageBitmap(await response.blob());
        if (S.view === want.mode) {
          Pic.shown = bitmap; const rest = S.state && S.state.rest ? S.state.rest : null;
          $('picture-note').title = 'Each frame minus the mean of the baseline frames, lightly smoothed for display, or the frame as recorded. Display only.';
          $('picture-note').textContent = want.mode === 'change' ? 'frame ' + want.f + ' − mean of frames ' + (rest ? rest[0] + '–' + rest[1] : 'at rest') + ', smoothed' : 'frame ' + want.f + ' as recorded';
          drawPicture();
        }
      } catch (problem) { /* the next request tries again */ }
    }
    Pic.loading = false;
  }
  function place(w, h) {
    const iw = S.info.width, ih = S.info.height, fitted = Math.max(0.01, Math.min((w - 8) / iw, (h - 8) / ih)), z = fitted * Pic.zoom;
    return { x: (w - iw * z) / 2 + Pic.panX, y: (h - ih * z) / 2 + Pic.panY, z, fitted };
  }
  function drawPicture() {
    const canvas = $('picture'); if (!S.info || !S.info.hasPicture || !canvas.clientWidth) return;
    const { ctx, w, h } = fit(canvas); ctx.fillStyle = '#10161C'; ctx.fillRect(0, 0, w, h);
    const image = Pic.shown || Pic.base; if (!image) return;
    const at = place(w, h); ctx.imageSmoothingEnabled = Pic.zoom <= 3; ctx.drawImage(image, at.x, at.y, S.info.width * at.z, S.info.height * at.z);
    const selected = S.cell ? S.info.names.indexOf(S.cell) : -1, known = !!S.results;
    ctx.save(); ctx.translate(at.x, at.y); ctx.scale(at.z, at.z); ctx.lineJoin = 'round';
    const stroke = (path, colour, width, dash) => { ctx.setLineDash(dash ? [4 / at.z, 3 / at.z] : []); ctx.strokeStyle = colour; ctx.lineWidth = width / at.z; ctx.stroke(path); };
    if (S.outlines) S.paths.forEach((path, i) => {
      if (i === selected) return;
      const call = callOf(S.info.names[i]); if (known && !inGroup(call, S.group)) return;
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
    box.innerHTML = parts.join('');
  }
  /** Everything except the time cursor and the hovered trace, drawn once and reused while the cursor moves. */
  function paintTraces(ctx, w, h) {
    const r = S.results, time = r.time, L = 54, T = 10, W = Math.max(10, w - L - 14), H = Math.max(10, h - T - 38);
    const full0 = time[0], full1 = time[time.length - 1], x0 = Plot.x0 === null ? full0 : Plot.x0, x1 = Plot.x1 === null ? full1 : Plot.x1;
    const drawn = S.every ? S.listed.filter((n) => r.traces[n]) : []; Plot.drawn = drawn;
    const values = [0]; r.mean.forEach((v, i) => { if (finite(v)) { values.push(v); if (finite(r.sem[i])) values.push(v - r.sem[i], v + r.sem[i]); } });
    const own = S.cell ? r.traces[S.cell] : null; if (own) for (const v of own) if (finite(v)) values.push(v);
    // with every cell drawn, the axis covers 99 % of all values (plus the mean and the selected cell); rarer extremes run off the edge
    if (drawn.length) { const rows = drawn.map((n) => r.traces[n]); values.push(share(rows, 0.005, false), share(rows, 0.995, false)); }
    let lo = Math.min.apply(null, values), hi = Math.max.apply(null, values); if (hi - lo < 0.001) hi = lo + 1; const pad = (hi - lo) * 0.08; lo -= pad; hi += pad;
    const X = (t) => L + (t - x0) / (x1 - x0) * W, Y = (v) => T + (hi - v) / (hi - lo) * H;
    Plot.geom = { L, T, W, H, x0, x1, lo, hi, full0, full1 };
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); ctx.font = '12px ' + getComputedStyle(document.body).fontFamily;
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, W, H); ctx.clip();
    for (const [span, colour, name] of [[r.baseline, '#EDF0F4', 'baseline'], [r.stimulus, '#FFF4E3', 'stimulus']]) {
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
    const line = (y) => { ctx.beginPath(); let pen = false; for (let i = 0; i < time.length && i < y.length; i++) { if (finite(y[i])) { const px = X(time[i]), py = Y(y[i]); if (pen) ctx.lineTo(px, py); else ctx.moveTo(px, py); pen = true; } else pen = false; } ctx.stroke(); };
    if (drawn.length) {
      ctx.globalAlpha = drawn.length > 120 ? 0.22 : drawn.length > 40 ? 0.32 : 0.48; ctx.lineWidth = 0.9;
      // non-responders first, so the responders are not hidden under them
      for (const responders of [false, true]) for (const name of drawn) { const call = callOf(name); if ((call === 'Responded') === responders) { ctx.strokeStyle = CALL[kind(call)]; line(r.traces[name]); } }
      ctx.globalAlpha = 1;
    }
    ctx.fillStyle = drawn.length ? 'rgba(0,114,178,0.28)' : 'rgba(0,114,178,0.20)';
    for (let i = 1; i < time.length; i++) if ([r.mean[i - 1], r.mean[i], r.sem[i - 1], r.sem[i]].every(finite)) {
      ctx.beginPath(); ctx.moveTo(X(time[i - 1]), Y(r.mean[i - 1] - r.sem[i - 1])); ctx.lineTo(X(time[i]), Y(r.mean[i] - r.sem[i])); ctx.lineTo(X(time[i]), Y(r.mean[i] + r.sem[i])); ctx.lineTo(X(time[i - 1]), Y(r.mean[i - 1] + r.sem[i - 1])); ctx.closePath(); ctx.fill();
    }
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
      const y = r.traces[Plot.hover]; ctx.strokeStyle = CALL[kind(callOf(Plot.hover))]; ctx.lineWidth = 2.2; ctx.lineJoin = 'round'; ctx.beginPath(); let pen = false;
      for (let i = 0; i < time.length; i++) { if (finite(y[i])) { if (pen) ctx.lineTo(X(time[i]), Y(y[i])); else ctx.moveTo(X(time[i]), Y(y[i])); pen = true; } else pen = false; } ctx.stroke();
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
    Heat.geom = { W, H, rows: names.length };
    // one pixel per sample and cell, then stretched without smoothing: sharp edges at any size
    const small = document.createElement('canvas'); small.width = time.length; small.height = names.length;
    const sctx = small.getContext('2d'), data = sctx.createImageData(time.length, names.length);
    rows.forEach((trace, c) => { for (let f = 0; f < time.length; f++) {
      const v = trace ? trace[f] : null, o = (c * time.length + f) * 4; let colour = missing;
      if (finite(v)) { const t = signed ? (Math.max(-1, Math.min(1, v / limit)) + 1) / 2 : Math.max(0, Math.min(1, v / limit)), k = Math.round(t * 255) * 3; colour = [colours[k], colours[k + 1], colours[k + 2]]; }
      data.data[o] = colour[0]; data.data[o + 1] = colour[1]; data.data[o + 2] = colour[2]; data.data[o + 3] = 255; } });
    sctx.putImageData(data, 0, 0);
    const family = getComputedStyle(document.body).fontFamily;
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h); ctx.imageSmoothingEnabled = false; ctx.drawImage(small, LEFT, TOP, W, H);
    const X = (t) => LEFT + (t - time[0]) / (time[time.length - 1] - time[0]) * W;
    // where the baseline and the stimulus are
    ctx.font = '11px ' + family; ctx.textBaseline = 'alphabetic'; ctx.textAlign = 'left';
    for (const [span, colour, name] of [[r.baseline, '#B8C2CC', 'baseline'], [r.stimulus, '#E9A84C', 'stimulus']]) {
      if (!span || !finite(span[0]) || !finite(span[1])) continue; const a = Math.max(LEFT, X(span[0])), b = Math.min(LEFT + W, X(span[1])); if (b <= a) continue;
      ctx.fillStyle = colour; ctx.fillRect(a, TOP - 7, b - a, 4); if (ctx.measureText(name).width + 4 < b - a) { ctx.fillStyle = '#65717D'; ctx.fillText(name, a + 1, TOP - 11); }
    }
    if (r.stimulus && finite(r.stimulus[0]) && r.stimulus[0] > time[0]) { const x = Math.round(X(r.stimulus[0])) + 0.5; ctx.strokeStyle = signed ? 'rgba(31,41,51,0.55)' : 'rgba(255,255,255,0.75)'; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(x, TOP); ctx.lineTo(x, TOP + H); ctx.stroke(); ctx.setLineDash([]); }
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
    const time = S.results.time, f = Math.max(0, Math.min(time.length - 1, Math.floor((x - LEFT) / g.W * time.length)));
    return { name: Heat.names[row], frame: f, inside: x >= LEFT };
  }

  // ---------------------------------------------------------------- tables (the saved tables themselves; long decimals are shortened for reading only)
  function table(id, data) {
    const element = $(id); element.innerHTML = ''; if (!data || !data.header) return;
    const cellColumn = data.header.indexOf('Cell'), numeric = data.header.map((_, c) => data.rows.some((row) => row[c] !== '' && row[c] !== undefined) && data.rows.every((row) => row[c] === '' || row[c] === undefined || isFinite(Number(row[c]))));
    const head = element.createTHead().insertRow(), body = element.createTBody(); let sorted = -1, direction = 1;
    const fill = (rows) => { body.innerHTML = ''; for (const row of rows) { const tr = body.insertRow(); if (cellColumn >= 0) { tr.dataset.cell = row[cellColumn]; tr.classList.toggle('on', row[cellColumn] === S.cell); }
      row.forEach((value, c) => { const td = tr.insertCell(); if (numeric[c] && value !== '') { const v = Number(value); td.className = 'num'; td.textContent = Number.isInteger(v) ? String(v) : String(+v.toPrecision(4)); td.title = value; } else td.textContent = value; }); } };
    data.header.forEach((name, c) => { const th = document.createElement('th'); th.textContent = name; head.appendChild(th);
      th.addEventListener('click', () => { direction = sorted === c ? -direction : 1; sorted = c; for (const other of head.children) other.className = ''; th.className = direction > 0 ? 'up' : 'down';
        fill(data.rows.slice().sort((a, b) => direction * (numeric[c] ? (Number(a[c] === '' ? NaN : a[c]) || 0) - (Number(b[c] === '' ? NaN : b[c]) || 0) : String(a[c]).localeCompare(String(b[c]))))); }); });
    body.addEventListener('click', (event) => { const tr = event.target.closest('tr'); if (tr && tr.dataset.cell) pick(tr.dataset.cell); });
    fill(data.rows);
  }

  // ---------------------------------------------------------------- controls
  function showTab(name) {
    for (const tab of document.querySelectorAll('.tab')) tab.classList.toggle('on', tab.dataset.tab === name);
    for (const panel of document.querySelectorAll('.panel')) panel.hidden = panel.id !== 'panel-' + name;
    if (name === 'heatmap') drawHeatmap(); if (name === 'review') { Plot.dirty = true; drawPlot(); drawPicture(); }
    const on = document.querySelector('tbody tr.on'); if (on && name !== 'review' && name !== 'heatmap') on.scrollIntoView({ block: 'nearest' });
  }
  function save(canvas, name) { canvas.toBlob((blob) => { const link = document.createElement('a'); link.href = URL.createObjectURL(blob); link.download = name; link.click(); setTimeout(() => URL.revokeObjectURL(link.href), 2000); }, 'image/png'); }
  function wire() {
    for (const tab of document.querySelectorAll('.tab')) tab.addEventListener('click', () => showTab(tab.dataset.tab));
    for (const tile of document.querySelectorAll('.tile')) tile.addEventListener('click', () => { if (S.results) { showGroup(Number(tile.dataset.group)); showTab('review'); } });
    $('toggle-settings').addEventListener('click', () => { const hidden = document.body.classList.toggle('no-side'); $('toggle-settings').textContent = hidden ? 'Show settings' : 'Hide settings'; redraw(); });
    $('kcl').addEventListener('change', () => { $('kcl-row').hidden = !$('kcl').checked; });
    $('analyze').addEventListener('click', () => { const values = {}; for (const box of S.boxes) values[box.id] = box.value; send('analyze', { values, confirmed: $('confirmed').checked, kcl: $('kcl').checked }); });
    $('setup').addEventListener('click', () => send('setup')); $('cancel').addEventListener('click', () => send('cancel'));
    $('log').addEventListener('click', () => api('open', { what: 'log' }).catch((e) => { S.note = e.message; showState(); }));
    $('open-folder').addEventListener('click', () => api('open', { what: 'folder' }).catch((e) => { S.note = e.message; showState(); }));
    $('group').addEventListener('change', () => { S.group = Number($('group').value); listCells(); });
    $('cell').addEventListener('change', () => pick($('cell').value));
    $('previous').addEventListener('click', () => stepCell(-1)); $('next').addEventListener('click', () => stepCell(1));
    $('every').addEventListener('change', () => { S.every = $('every').checked; Plot.dirty = true; drawPlot(); legend(); });
    $('view').addEventListener('change', () => { S.view = $('view').value; S.viewChosen = true; requestFrame(); });
    $('outlines').addEventListener('change', () => { S.outlines = $('outlines').checked; drawPicture(); });
    $('zoom-cell').addEventListener('click', zoomToCell); $('zoom-all').addEventListener('click', () => { Pic.zoom = 1; Pic.panX = Pic.panY = 0; drawPicture(); });
    $('play').addEventListener('click', () => play(!S.timer)); $('time').addEventListener('input', () => setFrame(Number($('time').value)));
    $('order').addEventListener('change', () => { S.order = $('order').value; drawHeatmap(); });
    $('palette').addEventListener('change', () => { S.palette = $('palette').value; drawHeatmap(); });
    $('range').addEventListener('change', () => { S.range = $('range').value; drawHeatmap(); });
    $('tall').addEventListener('change', () => { S.tall = $('tall').checked; drawHeatmap(); });
    $('save-plot').addEventListener('click', () => save($('plot'), 'traces.png')); $('save-heatmap').addEventListener('click', () => save($('heatmap'), 'heatmap.png'));

    const picture = $('picture');
    picture.addEventListener('wheel', (event) => { event.preventDefault(); const box = picture.getBoundingClientRect(), mx = event.clientX - box.left, my = event.clientY - box.top, before = place(picture.clientWidth, picture.clientHeight);
      const ix = (mx - before.x) / before.z, iy = (my - before.y) / before.z; Pic.zoom = Math.max(1, Math.min(16, Pic.zoom * Math.pow(1.0015, -event.deltaY)));
      if (Pic.zoom === 1) Pic.panX = Pic.panY = 0; else { const after = place(picture.clientWidth, picture.clientHeight); Pic.panX += mx - (after.x + ix * after.z); Pic.panY += my - (after.y + iy * after.z); } drawPicture(); }, { passive: false });
    picture.addEventListener('mousedown', (event) => { Pic.drag = { x: event.clientX, y: event.clientY }; Pic.moved = false; });
    window.addEventListener('mousemove', (event) => { if (!Pic.drag) return; if (Math.abs(event.clientX - Pic.drag.x) + Math.abs(event.clientY - Pic.drag.y) > 3) Pic.moved = true;
      if (Pic.moved) { Pic.panX += event.clientX - Pic.drag.x; Pic.panY += event.clientY - Pic.drag.y; Pic.drag = { x: event.clientX, y: event.clientY }; drawPicture(); } });
    window.addEventListener('mouseup', () => { Pic.drag = null; Plot.scrubbing = false; });
    picture.addEventListener('click', (event) => { if (Pic.moved) return; const i = cellAt(event); if (i >= 0 && S.results) pick(S.info.names[i]); });
    picture.addEventListener('dblclick', () => { Pic.zoom = 1; Pic.panX = Pic.panY = 0; drawPicture(); });
    picture.addEventListener('mousemove', (event) => { if (Pic.drag) return; const i = cellAt(event); tip(i >= 0 ? S.info.names[i] + (callOf(S.info.names[i]) ? ' · ' + callOf(S.info.names[i]) : '') : '', event.clientX, event.clientY); });
    picture.addEventListener('mouseleave', () => tip(''));

    const plot = $('plot');
    plot.addEventListener('mousemove', (event) => { const at = traceAt(event); if (!at) { tip(''); if (Plot.hover) { Plot.hover = null; drawPlot(); } return; }
      if (Plot.scrubbing) { seek(at.time); tip(''); return; }
      if (at.name !== Plot.hover) { Plot.hover = at.name; drawPlot(); }
      tip(at.name ? at.name + ' · ' + S.results.traces[at.name][at.index].toFixed(3) + ' at ' + S.results.time[at.index].toFixed(1) + ' s' + (callOf(at.name) ? ' · ' + callOf(at.name) : '') : '', event.clientX, event.clientY); });
    plot.addEventListener('mouseleave', () => { tip(''); if (Plot.hover) { Plot.hover = null; drawPlot(); } });
    plot.addEventListener('mousedown', (event) => { const at = traceAt(event); if (!at) return; if (at.name && at.name !== S.cell) { pick(at.name); return; } Plot.scrubbing = true; seek(at.time); });
    plot.addEventListener('wheel', (event) => { const g = Plot.geom; if (!g) return; event.preventDefault(); const box = plot.getBoundingClientRect(), share01 = Math.max(0, Math.min(1, (event.clientX - box.left - g.L) / g.W));
      const centre = g.x0 + share01 * (g.x1 - g.x0), factor = Math.pow(1.0015, event.deltaY), span = Math.min(g.full1 - g.full0, Math.max((g.full1 - g.full0) / 60, (g.x1 - g.x0) * factor));
      let a = centre - share01 * span, b = a + span; if (a < g.full0) { b += g.full0 - a; a = g.full0; } if (b > g.full1) { a -= b - g.full1; b = g.full1; }
      Plot.x0 = Math.max(g.full0, a); Plot.x1 = Math.min(g.full1, b); if (Plot.x1 - Plot.x0 >= g.full1 - g.full0 - 1e-9) Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); }, { passive: false });
    plot.addEventListener('dblclick', () => { Plot.x0 = Plot.x1 = null; Plot.dirty = true; drawPlot(); });

    const heatmap = $('heatmap');
    heatmap.addEventListener('click', (event) => { const at = heatmapAt(event); if (at) { pick(at.name); showTab('review'); } });
    heatmap.addEventListener('mousemove', (event) => { const at = heatmapAt(event); if (!at) { tip(''); return; } const trace = S.results.traces[at.name], v = trace ? trace[at.frame] : null;
      tip(at.name + (callOf(at.name) ? ' · ' + callOf(at.name) : '') + (at.inside ? ' · ' + (finite(v) ? v.toFixed(3) : 'no value') + ' at ' + S.results.time[at.frame].toFixed(1) + ' s' : ''), event.clientX, event.clientY); });
    heatmap.addEventListener('mouseleave', () => tip(''));

    document.addEventListener('keydown', (event) => { const tag = (event.target.tagName || '').toLowerCase(); if (tag === 'input' || tag === 'select' || tag === 'textarea' || $('panel-review').hidden) return;
      if (event.key === 'ArrowLeft') { stepCell(-1); event.preventDefault(); } else if (event.key === 'ArrowRight') { stepCell(1); event.preventDefault(); } else if (event.key === ' ' && S.results) { play(!S.timer); event.preventDefault(); } });
    window.addEventListener('resize', redraw);
    if (window.ResizeObserver) { new ResizeObserver(redraw).observe($('panel-review')); new ResizeObserver(() => drawHeatmap()).observe($('heatmap-box')); }
  }
  function redraw() { Plot.dirty = true; drawPlot(); drawPicture(); drawHeatmap(); }

  wire();
  start().catch((problem) => { $('lost').classList.add('on'); $('status').textContent = problem.message; });
  window.calciumStep4 = { S, Plot, Pic, Heat, pick, setFrame, showGroup, showTab };   // for the browser test
})();

# Changelog

## 2026.11.20.7

- Fix a false Step 4 warning introduced in 2026.11.20.6: "the response baseline holds 0 samples … every cell will be unscorable" appeared for any baseline. Baseline samples are now counted as the analysis engine counts them (baseline start ≤ frame time < baseline end), and the warning names the frames.
- Step 4 stimulus timing: **Noise baseline to** and **Score response from** now follow **Stimulus added at** by default and are shown as "= stimulus time". Changing the stimulus time, or pressing **Use frame N**, moves them with it, using the identical stored number. **Set separately** gives a field its own value and **Link to stimulus time** reverses it. A protocol whose values differ opens with those fields separate. The protocol file format is unchanged. The legacy rule works the same way, with the stimulus time as the one editable field.
- Step 4: each window shows the frames it will use, for example "F₀: frames 1–20. Noise baseline: frames 1–20 (20 samples). Stimulus: between frames 20 and 21. Scoring: frames 21–300.", and the noise baseline, stimulus and scored frames of the primary stimulus are drawn under the frame slider.
- Step 4: when two time fields disagree, both are marked and one message sits beside the field to change, stating both values and offering the repair, for example "Set to 20 s". A problem is described once: the plugin's general message is not added to the page's specific one.
- Step 4: a typed time that reads the same as a frame's time is stored as that frame's exact time, so a retyped value no longer differs from an untouched one by a hidden fraction of a second.
- Step 4: the warning that the F₀ frames reach the stimulus offers **Use frames 1–N**. The starting stimulus time is labelled as a suggestion until it is changed. "Resting baseline frames" is now **F₀ frames (set the ΔF/F₀ zero)** and "Response baseline" is **Noise baseline**; the fields are in the order they occur in time.
- The plugin window now stays above the picture being worked on only. Any other Fiji window that is clicked, such as Brightness/Contrast or the ROI Manager, comes in front of it.

## 2026.11.20.6

- Step 2 has a **Picture** switch: **All frames**, **Baseline…** and **Response (magenta)**. Magenta can now be turned on after the recording is open; the response is looked for at that moment, exactly as the first screen's option does, and the picture changes in place so outlines and clicks keep working. The first screen's option still works and still starts unticked.
- Fix **Review baseline frames**, which could not show its picture (ImageJ refused to put a flat picture into the colour composite). The baseline view now offers the detected baseline frames when there is a response.
- The selection basis recorded with each measurement follows the picture in use: all frames, baseline frames, or response-assisted. When "look for a response" was ticked but no response was found, it is now recorded as all frames rather than response-assisted.
- The saved `picture_colour.tif` is always the all-frames colour picture, whatever view is on screen, with a magenta channel only when a response was found, as before.
- **‹ Change recording** in the Recording card goes back to the first screen. If the ROI Manager holds outlines, it asks whether to clear or keep them, because the next recording starts with whatever outlines it holds. **‹ Back to cells** and **‹ Back to measuring** move back between steps.
- The plugin window stays above Fiji's own windows only while Fiji is in use, and drops back when another program, such as a browser, is in front.
- Step 4: a stimulus window and a positive-control window that overlap, or a baseline that falls inside the other window's response interval, are now flagged while editing, saying what to change. The analysis engine has always refused these, but the run used to fail with only "exit code 1".
- Step 4: when the analysis engine stops, the page and the window show its own error message instead of only the exit code.
- Step 4: a run in which no cell could be scored for the primary stimulus is reported as a warning naming the commonest reason, not as a green "Complete". The saved results are unchanged.
- Step 4 exports are drawn at least 1600 pixels wide, whatever the window size, and their names include the recording, the cell for trace plots, and the time, for example `Recording_traces_Cell001_20261120-153012.png`. The metadata records the results folder and the time range shown.
- Step 4: a protocol file that is not valid JSON, or whose name was edited in only one place, gets a plain explanation instead of a parser message. After a cancel, the progress stage is cleared. The positive-control tick box now sits above the control's own fields.
- Step 4 quality review: the "minimum response duration is shorter than one frame interval" flag, which the engine puts on every cell, is shown once as a notice for the recording and no longer puts every cell into **Needs quality review**. Each flag in the summary is a filter that lists only its cells. Added plain wording for the baseline-noise flag. The saved tables are unchanged.
- Step 4: the headline tiles follow the response window selected, so the positive control shows its own counts.
- Step 4: cells with long stable IDs are shown as "Cell N" (their measured order) on the picture, plot, heatmap and lists, with the full ID in tooltips. Tables and exports keep the IDs.
- Step 4: **Use frame on screen** under "Stimulus added at" sets the onset from the movie frame shown.
- Step 4: a warning while editing when a response baseline cannot hold the minimum number of baseline samples, which would leave every cell unscorable. Discard and protocol import clear old messages. Field names in number errors use their on-screen labels.
- Step 4 movie playback is faster and steadier. Speeds now run from 0.5× to 50×, with **10×** as the default. Frames are shown on a clock in recording time, so they arrive at even intervals. Above 30 frames per second on screen, every Nth frame is shown and the frame note says so. A long gap between recorded frames waits at most about one second.
- Step 4 frames are loaded up to 4 at a time, about 3 s of playback ahead, and kept for the session (within 90 MB), so scrubbing back over frames already seen makes no new requests. If loading falls behind, the picture holds and the note says "buffering". Frames are fetched at the size they are drawn (full size when zoomed in) and may be cached by the browser. Each frame is keyed to the baseline and display range, so a changed baseline never shows an old image. While playing, the trace cursor glides with the clock. The displayed images and measurements are unchanged.

## 2026.11.20.5

- Validate Step 4 settings while editing, mark invalid fields and disable analysis until errors are corrected. Preserve separate elapsed-time and legacy drafts when switching rules; add Discard pending changes to restore the completed analysis or initial settings.
- Require actual positive-control timing instead of inventing a window beyond the recording. Place the suggested primary onset after the final normalization-baseline frame and show warnings for baseline/onset overlap and incomplete recording coverage.
- Surface saved quality flags, low-reference warnings, coarse temporal resolution and negative corrected-fluorescence traces in an expandable summary and per-cell feedback. Keep response rules and inclusion/exclusion policy unchanged; separate Excluded, Unscorable and Needs quality review filters.
- Give each measurement a separate browser session and preserved input files. Existing pages remain available when another measurement is opened or the cell-selection window is closed, until Close analysis or Fiji exits. Distinguish tabs by cell count and session identity, and give feedback when opening a log or results folder.
- Keep short-window layouts scrollable with readable plots and visible captions. Show core response-call columns first, pin cell names while scrolling, explain empty tables, and label legacy whole-recording metrics separately from stimulus-window metrics.
- Reject invalid frame numbers and display modes at the local frame endpoint. Preserve full frame-rate precision in the settings box and keep completed results and exports tied to their applied settings.

## 2026.11.20.4

- Move Check for updates, Install update, Later, Release notes and the startup-check preference onto Step 1, before choosing or opening a recording. Show update offers and installation progress inline on that screen.
- Keep recording selection responsive during update checks, including offline startup. Pause Next/Cancel only while installing; show retryable errors or the installed version and restart instructions without losing the chosen recording or folder.
- Remove the duplicate update area and delayed automatic check from the cell-selection window. Make Step 1 scrollable on smaller screens, with Next and Cancel always outside the scrolling area.

## 2026.11.20.3

- Bring the portable response-protocol workflow into the local Step 4 browser page while preserving the 2026.11.20.2 movie, interactive traces, heatmap, cell navigation and native Fiji fallback. Browser displays read the shared engine's saved results; they do not recalculate response calls.
- Open the browser in **Recording** view with prominent Play, frame-step, slider and speed controls available before analysis. Playback uses the available frame times; when no valid time base is available, the preview explicitly labels its 20-frames-per-second fallback at 1×. Label the colour overview as static, disable movie playback in that view, and hide result tiles until an analysis has completed.
- Keep stimulus timing and a response-rule summary visible while placing detailed rules, protocol files and diagnostics in expandable sections. Pending edits leave the completed results, their movie timing and exports unchanged until **Apply changes & reanalyze** completes.
- Add versioned protocols shared with SimpleCa²⁺ commit b92d4ce98ddab3fceac352179703c1f48f044c90: independent baseline/onset/scoring windows, one generically named positive control, editable per-window floor/noise/persistence rules, JSON import/export and retained field timing overrides.
- Offer a provisional 0.10/3-SD/1-second preset and preserve the actual historical schema-1 sample rule. Version 2 uses timestamp spans and explicit complete-window eligibility; missing data or inadequate coverage remain unscorable. Review flags do not automatically exclude low-reference or peak-outlier cells.
- Analyze full-precision unlifted corrected fluorescence. Preserve stable ROI identities, accepted/edge-excluded registries and the selection basis; response-assisted selection remains an optional review aid. Older rounded and intensity-lifted exports remain labeled compatibility outputs.
- Add exact timestamp CSVs, input/ROI/protocol snapshots, applied-versus-draft state, transformed local thresholds, qualifying-run evidence, window peak/AUC/latency measurements and plot provenance. Changes to timestamp or diagnostic file contents mark the draft as changed even at the same path.
- Add separately labeled sensitivity and negative-control diagnostics, condition/biological-unit metadata and summaries that distinguish biological units from cell-pooled percentages. Zero eligible cells produce an unavailable percentage.
- Keep legacy mode's inactive fields explicit, preserve imported stimulus/control criteria during timing edits, and reject inconsistent protocol metadata. The managed runtime verifies jsonlite support. Extend independent scientific fixtures, protocol/precision checks, local browser/server checks and native GUI regressions; bundle both the pinned R sources and browser assets in the single plugin file.

## 2026.11.20.2

- Step 4 now opens as a page in your browser, built with the same web technology and styling as the SimpleCa²⁺ app: results first, interactive traces (hover a trace to name it, click to select it, scroll to zoom the time axis), the recording played beside the traces, an ordered heatmap and sortable tables. The page is served by the plugin from your own computer; nothing is sent anywhere and no internet is needed.
- The page drives the same analysis engine with the same settings and checks as before, and reads the tables that engine saves. Nothing is recalculated in the browser. **Use a Fiji window instead** opens the previous Step 4 window, which is also used automatically if a browser page cannot be opened.
- The page is reachable only from this computer and only with a long random key in its address.
- The page's heatmap is drawn the way the SimpleCa²⁺ app draws its own: every cell in one picture, axes with ticks, a colour bar with values, the app's blue–white–red for signed values or one of its palettes (plasma, viridis, magma, inferno, cividis) from zero up, and sorting by response call, time to peak, peak amplitude or original order. The colour range covers 99 % of the values by default, or the largest value as in the app.
- Step 1: choosing **Choose a file from disk...** now opens the file chooser straight away, and a **Choose file…** button with the chosen file's name sits under the list. Before, the chooser only appeared after **Next**, so picking that entry looked as if nothing happened.
- **Next** still asks for the file if none was chosen, and stays on Step 1 if that is cancelled instead of closing the plugin.
- Step 4, movie beside the traces: a time slider, a **Play** button and a line on the plot move together, and the picture can show the recording at that moment. **Movie: change from rest** shows each frame minus the average of the baseline frames (lightly smoothed), so a responding cell lights up as its trace rises; **Movie: as recorded** shows the frame itself. **Zoom to cell** and **All outlines** help follow one cell. The recording is only read. Display only: nothing measured or saved changes.
- Step 4, every trace: all cells are drawn behind the mean, coloured by response call (the SimpleCa²⁺ colours), with the selected cell on top. The y-axis covers 99 % of all values so one extreme cell cannot flatten the rest.
- Step 4, cell groups: list and draw only responders, non-responders, or excluded and unscorable cells, from the **Show** list or by clicking a headline tile. A line under the list says what was decided for the selected cell, with the largest rise and the threshold from the saved table. Left and right arrow keys step through the cells.
- Step 4 heatmap: responders first with the strongest on top (or as measured), a colour scale that covers 99 % of the values instead of the single largest, taller rows with readable names, a strip for each cell's call, and markers for the baseline and the stimulus.
- Step 4 layout: narrower settings column that can be hidden, picture and traces side by side, larger default window. **Also score a separate KCl window** is now **Include a KCl control**, with an explanation.
- Main window: a long recording name or folder name is shortened in the middle (full text in the tooltip), so it no longer pushes the window contents off the right edge.

## 2026.11.20.1

- Redesign every plugin window to match SimpleCa²⁺: the computer's own interface font instead of a fixed one, the SimpleCa²⁺ colours, white cards, and one clearly marked main button per step. Measurements, thresholds and saved tables are unchanged.
- Main window: one card per step. Only one of Steps 2 to 4 is open at a time and finished steps collapse to a one-line summary, so the window fits the screen. The cell count wraps instead of running off the edge, long instructions moved to tooltips, the detection progress area is smaller, and version/update controls sit behind an **Updates** button until an update is found.
- Step 1 is a single styled dialog with a folder **Browse** button.
- Step 4: results first, with headline tiles for percent responding, registered, excluded and unscorable cells. Settings are grouped into **Timing** and **Response rule**, show rounded numbers, and hide the KCl boxes until that option is ticked. A box that still shows its rounded suggestion sends the exact value to the analysis, so results do not change.
- Step 4 plots: axis ticks, units and names, labelled baseline and stimulus bands, a colour key on the heatmap, and previous/next cell buttons.
- The cell review picture is now inside the Step 4 window, in full colour, with outlines styled by response call. Scroll to zoom, drag to move, double-click to see the whole picture.
- Windows no longer cover each other: the colour picture is fitted into the space left of the main window, and the Step 4 window opens in that space too.

## 2026.11.19.8

- Add Step 4, SimpleCa²⁺ analysis inside Fiji: confirmed timing/baseline controls, native traces and heatmap, cell-linked image review, metrics, stimulus/KCl responder calls and CSV/PNG export.
- Reuse the unmodified metric and responder functions pinned from SimpleCa²⁺ commit a472f2460075343d5b97c8bac32ddf3da0c77386. Preserve invalid-baseline cells as excluded records and keep unscorable cells explicit; use the unlifted corrected fluorescence only.
- Bundle the complete analysis source and license inside the plugin, so existing Install update buttons install the integration. Provide one-time, per-user R setup with a verified micromamba launcher, progress/logs, cancellation and retry.
- Save the exact input snapshot, settings, source hashes, R/plugin versions and results together. Preserve completed analyses; partial runs are never displayed as completed results.
- Add real R parity/runtime checks on Windows, macOS and Linux, plus process and native-window regressions.

## 2026.11.19.7

- Make accepted cells clearly distinct: white outlines with a thick black outer border, and a thicker highlight for selected cells. Edge exclusions remain dashed.
- Fix clicking accepted cells so Remove selected cell targets the clicked cell, even with suggestions hidden or no automatic suggestions. Support multiple ROI Manager selections.
- Add Click cells to remove for reviewing many accepted cells directly on the picture. This mode never accepts suggestions; drawing tools turn it off.
- Replace Undo last cell with Undo last change: restore removed cells or undo additions, including Accept all, for up to 20 plugin edits. Preserve ROI names, geometry and order; outside ROI Manager edits invalidate undo instead of being overwritten.
- Add native GUI regressions for actual image clicks, removal, undo, non-cell ROIs, overlap targeting and contrasted outlines.

## 2026.11.19.6

- Fix the visible-window pack exception caused by Groovy treating Rectangle coordinates as doubles. Use explicit integer coordinates and statically compile the window class.
- Start detection independently of cosmetic window updates. Startup failures stop the progress timer and restore the Suggest button; result-display errors cannot leave detection stuck.

- Add Install update to the update prompt and main window. It downloads the plugin inside Fiji, verifies the published SHA-256 checksum and version, backs up the previous file, and replaces it atomically.
- Show installation progress, retryable errors and a restart instruction. Keep the running analysis session open; the new version is used after restarting Fiji.
- Prevent concurrent installers and downgrades. Failed downloads, verification, backups or replacements leave the current plugin in place.
- Provide Install_Calcium_Plugin.groovy for first installation or upgrading versions without the installer: open it in Fiji's Script Editor and press Run. No manual file copying is needed.

## 2026.11.19.5

- Keep controls at their full layout size when the plugin window is made smaller. Scroll vertically with the mouse wheel or scrollbar, and horizontally when the window is too narrow.
- Preserve the window size chosen by the user when detection progress, results, or update messages change.
- Fit the initial window within the screen's width and height and keep a usable minimum window size.
- Add native Windows and Linux GUI checks for resizing, wheel scrolling, unclipped controls, and access to the bottom and right-hand buttons.

## 2026.11.19.4

- Keep a plain-language explanation visible during detection, with total and current-stage time, image dimensions, CPU/GPU, attempt number, and detection settings.
- Show whether saved outlines are being reused or a new detection is running. Add a scrollable Progress details window with recent events and Cellpose messages.
- Keep the main window within the screen height and make its contents scrollable so all controls stay reachable.
- Preserve stage explanations during quiet periods and report whether the Cellpose process is still running. Progress remains tied to reported stages, with no invented time estimate.
- Check for new releases in the background when the plugin window opens, with a non-blocking Download / Later prompt and an option to disable startup checks. Manual checks remain available; installation is still manual.
- Publish a tested release automatically when a higher version is merged to main. Existing releases are never overwritten.
- Install this version once to receive startup prompts for future releases; earlier versions can use Check for updates to obtain it.

## 2026.11.19.3

- Use the latest user-supplied plugin as the base, preserving edge-cell exclusion, matched-image outline reuse, faster CPU detection, and automatic Torch backend selection during setup.
- Integrate percentage stage progress, detection logs, cancellation, and GitHub update checking into that version.
- Prevent failed or interrupted detection from being reused as a completed cached result.
- Keep detection settings and edge-cell metadata in the saved records.

## 2026.11.19.2

- First public GitHub release, with a permanent latest-download link.
- Add an explicit **Check for updates** button and download-page link inside the plugin.
- Add release packaging, checksums, automated regression checks, and a maintainer publishing workflow.
- Include the Step 2 progress changes from 2026-11-19.1 below.

## 2026-11-19.1

- Show percentage stage progress, current stage, elapsed time, and CPU/GPU information when supplied by Cellpose.
- Enable verbose, unbuffered Cellpose logs and add **View detection log**.
- Add cancellation and feedback when a run has no new progress messages or exceeds ten minutes.
- Remove the fixed “about 30 seconds” estimate.
- Show 100% only after the suggestions are read, graded, and saved.
- Prevent a failed rerun from reusing an old outlines ZIP; report failures from the current attempt.

## 2026-11-19

- Original packaged workflow for choosing cells, local background correction, and fluorescence export.

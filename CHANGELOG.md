# Changelog

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

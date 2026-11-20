# Changelog

## 2026.11.20.2

- Add portable, versioned response protocols shared with SimpleCa²⁺ commit b92d4ce98ddab3fceac352179703c1f48f044c90: independent baseline/onset/scoring windows, generic positive controls, editable per-window floor/noise/persistence rules, JSON import/export and retained field timing overrides.
- Offer a provisional 0.10/3-SD/1-second preset and preserve the actual historical schema-1 sample rule for legacy analyses. Use timestamp spans, explicit missing/gap coverage and separately recorded eligibility/QC. Low-reference and zero-SD flags remain warnings.
- Analyze full-precision unlifted corrected fluorescence, preserve stable ROI identities and accepted/edge-excluded registries, record selection basis, and make response-assisted selection an optional review aid. Existing three-decimal and lifted legacy exports remain labeled compatibility outputs.
- Add exact timestamp CSVs, input/ROI/protocol snapshots, applied-versus-draft state, correct transformed local thresholds in native plots, qualifying-run evidence, window/latency/AUC metrics and plot provenance.
- Add separately labeled sensitivity and negative-control diagnostics, condition/biological-unit metadata and unit-aware summaries without treating cells as biological replicates.
- Preserve the 2026.11.20.1 interface theme, embedded cell picture, zoom/navigation, result tiles, collapsible workflow and Cellpose progress/selection behavior. Update the managed runtime to verify jsonlite support and extend scientific, precision, protocol-round-trip and native GUI regressions.

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

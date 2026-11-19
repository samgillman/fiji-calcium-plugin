# Changelog

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

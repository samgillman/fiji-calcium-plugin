# Changelog

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

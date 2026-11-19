# Fiji Calcium Plugin

This repository is the working source for the plugin. Use the current GitHub `main` and inspect existing local changes before starting new edits. Do not restore an older ZIP over the working source.

When the user supplies a newer source file or ZIP, compare it with the current source and preserve its features while integrating requested changes. Uploaded files can share the same date/version header despite differing content. The 2026.11.19.3 release incorporates the user's latest 19 November source, including edge-cell controls, saved-outline reuse, CPU speed settings, and Torch backend selection.

Keep image measurements and scientific thresholds unchanged unless the user requests changes. Test progress, cancellation, successful and failed runs, cache reuse, and update metadata when touching those paths.

Run `python tools/run_tests.py` before publishing. Follow `docs/RELEASING.md`; published versions and tags are immutable. Changes to plugin behaviour require a new version and changelog entry when released. Documentation-only commits do not require a release.

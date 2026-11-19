# Publish an update

The public download is a GitHub Release, not whatever happens to be on `main`. Everyone can use the same latest-download link. Users who subscribe to Releases receive GitHub notifications; installations from version 2026.11.19.4 can show a startup update prompt, and the **Check for updates** button remains available.

## Normal publishing workflow

1. Make the plugin changes in a branch and open a pull request.
2. Choose a higher calendar version in the form `YYYY.MM.DD.N`, for example `2026.11.20.1`. Update `VERSION`, `CalciumRelease.VERSION` in the plugin, the plugin's displayed/logged version strings, and `CHANGELOG.md`. Version validation catches inconsistent plugin version strings.
3. Add a changelog section headed `## VERSION` for the new version. Review the user guide if behaviour changed.
4. Wait for **Checks** to pass, review the changes, and merge to `main`.
5. The merge automatically starts **Publish release**. Wait for the workflow to succeed. It runs regression checks on Linux and Windows, builds the assets, creates a draft, uploads the complete package, and then publishes it as the latest release.

The workflow uses the repository's `GITHUB_TOKEN` with `contents: write`; no personal token or additional secret is needed. A merge with an already-published VERSION passes checks and skips publishing, so documentation-only merges do not make releases. Published releases are never overwritten. Existing drafts or prereleases block publication for manual review; API errors also stop publication. A new version must be higher than the latest release.

Once publishing succeeds, updated installations offer Install update / Later the next time the plugin window opens (unless startup checks are disabled). This is not an instant push to already-open windows. Install update verifies, backs up and replaces the plugin file; restart Fiji after saving work to activate it. Versions before 2026.11.19.6 need the standalone installer run once. The release includes Install_Calcium_Plugin.groovy, generated from the same release/installation code as the plugin, so the bootstrap and in-plugin installer cannot drift.

The ZIP keeps the fixed name `Calcium_Imaging_Fiji_Plugin.zip`, so this URL remains valid:

https://github.com/samgillman/fiji-calcium-plugin/releases/latest/download/Calcium_Imaging_Fiji_Plugin.zip

## Recovery and explicit triggers

The normal merge workflow needs no tag or manual dispatch. If a run fails before publication, fix the cause and rerun it, or use **Actions → Publish release → Run workflow** on `main`. Review an existing draft before retrying; never overwrite a published release.

A matching tag can also trigger the workflow:

```sh
git tag v2026.11.20.1
git push origin v2026.11.20.1
```

The same publishing workflow runs for `v*` tags and checks that the tag matches `VERSION`.

## Local checks and packaging

Requires Java 17 or newer and Python 3.10 or newer. Tests use Groovy 4.0.28, groovy-json 4.0.28, and ImageJ 1.54p. The test runner downloads these version-pinned jars from Maven Central into ignored `build/deps`, or uses matching jars from a `FIJI_JARS` directory you supply.

```sh
python3 tools/run_tests.py
python3 tools/build_release.py
```

With a desktop display available, `python3 tools/run_tests.py --gui` also exercises the actual Install update/Later prompt, offline startup, detection details window, resizing/scrolling, image-canvas cell selection and removal, Undo, and rendered outline contrast with the real plugin controls. Release metadata is simulated. CI runs these GUI checks on Windows and on Linux with a virtual display.

The package, standalone plugin file, standalone installer, version file, checksum file, and release notes are written to ignored `dist/`. Tests use a mock Cellpose process and real ROI ZIPs; do a real Fiji/Cellpose smoke test for changes to detection or image measurements.

## Correct a bad release

Keep existing releases available for reproducibility. Make the fix under a higher version and publish again. Users can download an older release from the Releases page and replace their installed script manually if needed. Do not overwrite an already-published tag or asset.

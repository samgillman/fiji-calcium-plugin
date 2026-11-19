# Publish an update

The public download is a GitHub Release, not whatever happens to be on `main`. Everyone can use the same latest-download link. Users who subscribe to Releases receive GitHub notifications; existing installations can also use **Check for updates**.

## Normal publishing workflow

1. Make the plugin changes in a branch and open a pull request.
2. Choose a higher calendar version in the form `YYYY.MM.DD.N`, for example `2026.11.20.1`. Update `VERSION`, `CalciumRelease.VERSION` in the plugin, the plugin's displayed/logged version strings, and `CHANGELOG.md`. Version validation catches inconsistent plugin version strings.
3. Add a changelog section headed `## VERSION` for the new version. Review the user guide if behaviour changed.
4. Wait for **Checks** to pass, review the changes, and merge to `main`.
5. Go to **Actions → Publish release → Run workflow**, choose `main`, and run it. This is the publishing action.
6. Wait for the workflow to succeed. It runs regression checks on Linux and Windows, builds the assets, creates a draft, uploads the complete package, and then publishes it as the latest release.

The workflow uses the repository's `GITHUB_TOKEN` with `contents: write`; no personal token or additional secret is needed. Only repository users with write access can run it. Running a release twice with the same version fails instead of overwriting an existing release.

The ZIP keeps the fixed name `Calcium_Imaging_Fiji_Plugin.zip`, so this URL remains valid:

https://github.com/samgillman/fiji-calcium-plugin/releases/latest/download/Calcium_Imaging_Fiji_Plugin.zip

## Command-line alternative

After committing and pushing the new version to `main`, push a matching tag:

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

The package, standalone plugin file, version file, checksum file, and release notes are written to ignored `dist/`. Tests use a mock Cellpose process and real ROI ZIPs; do a real Fiji/Cellpose smoke test for changes to detection or image measurements.

## Correct a bad release

Keep existing releases available for reproducibility. Make the fix under a higher version and publish again. Users can download an older release from the Releases page and replace their installed script manually if needed. Do not overwrite an already-published tag or asset.

# Fiji Calcium Plugin

Choose cells in a calcium-imaging recording, review automatic Cellpose outlines, measure fluorescence and local background, and analyze responses with SimpleCa²⁺ inside Fiji.

**[Download the latest plugin ZIP](https://github.com/samgillman/fiji-calcium-plugin/releases/latest/download/Calcium_Imaging_Fiji_Plugin.zip)** · [Release notes and older versions](https://github.com/samgillman/fiji-calcium-plugin/releases) · [Full user guide](docs/USER_GUIDE.txt)

## Install in Fiji

1. **[Download the Fiji installer](https://github.com/samgillman/fiji-calcium-plugin/releases/latest/download/Install_Calcium_Plugin.groovy)**.
2. Open that `.groovy` file in Fiji's Script Editor (or drag it onto Fiji), then press **Run**. It finds this Fiji installation and downloads, verifies and installs the latest plugin. No manual copying is needed.
3. Save your work and restart Fiji. Open **Plugins → Calcium Imaging → Calcium Cells Step By Step**.

The same installer upgrades versions that do not yet have an Install update button. After this one-time setup, use **Install update** inside the plugin for future releases. Installation requires a Fiji folder your account can write to.

Manual installation remains available: extract the plugin ZIP, merge its `scripts` folder into the Fiji folder shown by **File → Show Folder → ImageJ**, and restart Fiji. The plugin belongs at `scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy`. Keep only one copy in Fiji's scripts folders.

Fiji includes Bio-Formats for microscope recordings. Cellpose is optional: you can draw cells by hand, or use the plugin's **Set up Cellpose** button for automatic outlines.

## Use it

1. Choose the original recording and an output folder.
2. Choose cells on the colour image. Draw outlines or click **Suggest cells automatically**, then review the suggestions. Automatic detection shows percentage stage progress, a persistent explanation of the current work, total and stage time, saved-outline reuse, image size, CPU/GPU and detection settings. **Details** shows recent events; **Log** and **Cancel** remain available.
3. Click **Measure and save** to export raw fluorescence, each cell's local background, corrected traces, ROI outlines, and an analysis log.
4. Open **Analyze responses** in Step 4. On first use, click **Set up analysis engine** to install a private R runtime. Confirm the frame rate or load timestamps, choose a versioned response protocol, and review each stimulus/control window before analysis. Review native traces, heatmap, metrics and responder counts; selecting a cell also highlights it on the cell picture inside the same window. See [the integrated analysis guide](docs/SIMPLECA.md).

Accepted cells use white outlines with a black outer border. Click one to select it, then **Remove selected**, or enable **Click cells to remove** to remove cells directly from the picture. **Undo** restores removals and reverses additions. Turn removal mode off to accept suggestions again.

You can resize the plugin window. If its contents do not fit, use the mouse wheel or vertical scrollbar; a horizontal scrollbar appears when needed. Controls keep their full size, and status changes preserve your chosen window size.

This release preserves the latest base version's edge-cell exclusion control, reuse of outlines for matching images, and faster CPU detection settings. **Suggest cells again** forces a fresh detection.

The progress percentage represents detection stages, not remaining time. Cellpose may provide no intermediate updates during inference. Detection can take several minutes; one user-reported Windows run took about 700 seconds. A first run may also download the model.

## Get updates

- **Inside the plugin:** startup checks and **Updates → Check for updates** offer **Install update** or **Later** when a newer release is available. Install update downloads the plugin, verifies its SHA-256 checksum and version, backs up the previous file, and replaces it safely. Save your work and restart Fiji to load the new version; your running session remains open.
- **Earlier versions:** versions before 2026.11.19.6 need the installer above run once to gain one-click installation. The older update button still opens a release page; choose `Install_Calcium_Plugin.groovy` there instead of copying files manually.
- **Backup and recovery:** previous scripts are saved with `.groovy.bak` filenames in `calcium-plugin-backups` inside Fiji's folder. To roll back, close Fiji and copy the chosen backup over the installed `.groovy` script. Existing recordings, Cellpose environments, settings and results are preserved.
- **Preferences:** under **Updates**, turn off **Check for updates when this plugin opens** if preferred. Manual checks remain available.
- **GitHub notifications:** sign in and choose **Watch → Custom → Releases → Apply** on this repository.

Installation happens only after you click Install update or run the standalone installer. Updates are checked and downloaded from GitHub; no recordings, images, traces or results are sent. Failed downloads or verification leave the installed plugin unchanged. If Fiji cannot write to its folder, installation reports the problem and can be retried after correcting folder permissions.

## Scientific scope and validation

The defaults were developed for 1388 × 1040 recordings at 0.645 µm/pixel. Review the user guide and settings for other imaging conditions. Automatic outlines and grades require review. Select identifiable cells independently of their response. The optional magenta display is a review aid and starts disabled. Intensity-lifted legacy tables must not be used for fractional response percentages; native analysis uses full-precision unlifted values.

The original workflow was tested on macOS with Fiji/ImageJ 1.54p. A Windows user reported successful detection, but the newest UI/update features have not been validated in a complete Windows imaging session. Automated tests exercise progress, cancellation, retries, ROI import, and update metadata on simulated detector runs; they do not establish segmentation accuracy.

The integrated analysis supports 64-bit Windows, Intel/Apple Silicon Mac and x64 Linux. CI is configured to install R/jsonlite and check independent synthetic outcomes and parity with the pinned SimpleCa source on Windows, macOS and Linux. Step 4 exchanges protocol JSON with SimpleCa²⁺, retains excluded/unscorable records and reports the evaluable denominator. Its provisional 1-second preset is a starting rule for control-based validation; true legacy sample rules remain available. Completed results keep their applied settings even while a new protocol is being edited. Native analysis handles one recording at a time; it does not include all of the web application's multi-recording plot controls.

## Report a problem

Open an [issue](https://github.com/samgillman/fiji-calcium-plugin/issues) with the plugin version, operating system, Fiji version, steps to reproduce, and the relevant error message. Review logs before sharing them: they can include local filenames and paths. Do not upload confidential recordings.

## Publish an update

Maintainers: see [the release guide](docs/RELEASING.md). Merging a higher version into `main` starts the **Publish release** GitHub Action, which tests and packages the plugin, then publishes versioned downloads and SHA-256 checksums. No manual ZIP upload is needed for subsequent releases.

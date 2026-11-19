# Fiji Calcium Plugin

Choose cells in a calcium-imaging recording, review automatic Cellpose outlines, and export each cell's fluorescence and local background from the original movie.

**[Download the latest plugin ZIP](https://github.com/samgillman/fiji-calcium-plugin/releases/latest/download/Calcium_Imaging_Fiji_Plugin.zip)** · [Release notes and older versions](https://github.com/samgillman/fiji-calcium-plugin/releases) · [Full user guide](docs/USER_GUIDE.txt)

## Install in Fiji

1. Download and extract **Calcium_Imaging_Fiji_Plugin.zip** from the link above. Use this release asset; GitHub's automatically generated “Source code” archive is for development.
2. In Fiji, choose **File → Show Folder → ImageJ**.
3. Copy the ZIP's `scripts` folder into that Fiji folder, merging folders when asked. The plugin belongs at `scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy`.
4. Restart Fiji. Open **Plugins → Calcium Imaging → Calcium Cells Step By Step**.

If you installed an earlier ZIP, replace the existing `.groovy` file in that same location. Keep only one copy of the plugin in Fiji's scripts folders.

Fiji includes Bio-Formats for microscope recordings. Cellpose is optional: you can draw cells by hand, or use the plugin's **Set up Cellpose** button for automatic outlines.

## Use it

1. Choose the original recording and an output folder.
2. Choose cells on the colour image. Draw outlines or click **Suggest cells automatically**, then review the suggestions. Automatic detection shows percentage stage progress, a persistent explanation of the current work, total and stage time, saved-outline reuse, image size, CPU/GPU and detection settings. **Progress details** shows recent events; **View detection log** and **Cancel detection** remain available.
3. Click **Measure and save** to export raw fluorescence, each cell's local background, corrected traces, ROI outlines, and an analysis log.

You can resize the plugin window. If its contents do not fit, use the mouse wheel or vertical scrollbar; a horizontal scrollbar appears when needed. Controls keep their full size, and status changes preserve your chosen window size.

This release preserves the latest base version's edge-cell exclusion control, reuse of outlines for matching images, and faster CPU detection settings. **Suggest cells again** forces a fresh detection.

The progress percentage represents detection stages, not remaining time. Cellpose may provide no intermediate updates during inference. Detection can take several minutes; one user-reported Windows run took about 700 seconds. A first run may also download the model.

## Get updates

- **Inside the plugin:** version 2026.11.19.4 and newer check for a release in the background when the plugin window opens. If a newer release is ready, choose **Download** to open its release page or **Later** to keep working. Turn off **Check for updates when this plugin opens** if preferred; **Check for updates** still works manually. Download the ZIP, replace the plugin file, and restart Fiji after saving your work.
- **One-time upgrade:** earlier versions need this version installed once before they can show startup prompts for future updates.
- **Release notifications:** sign in to GitHub, click **Watch → Custom → Releases → Apply** on this repository. GitHub will notify you about new releases according to your notification settings.
- **Permanent download link:** the “Download the latest plugin ZIP” link above always follows the latest published release.

Updates are not silently installed. The update check makes a request to GitHub and sends no recordings, images, traces, or analysis results. GitHub downloaders are not automatically subscribed to notifications.

## Scientific scope and validation

The defaults were developed for 1388 × 1040 recordings at 0.645 µm/pixel. Review the user guide and settings for other imaging conditions. Automatic outlines and grades require review. The sizes of ΔF/F₀ responses from intensity-lifted, dim cells are not interpretable; see the guide's discussion of the “every cell kept” table.

The original workflow was tested on macOS with Fiji/ImageJ 1.54p. A Windows user reported successful detection, but the newest UI/update features have not been validated in a complete Windows imaging session. Automated tests exercise progress, cancellation, retries, ROI import, and update metadata on simulated detector runs; they do not establish segmentation accuracy.

## Report a problem

Open an [issue](https://github.com/samgillman/fiji-calcium-plugin/issues) with the plugin version, operating system, Fiji version, steps to reproduce, and the relevant error message. Review logs before sharing them: they can include local filenames and paths. Do not upload confidential recordings.

## Publish an update

Maintainers: see [the release guide](docs/RELEASING.md). Merging a higher version into `main` starts the **Publish release** GitHub Action, which tests and packages the plugin, then publishes versioned downloads and SHA-256 checksums. No manual ZIP upload is needed for subsequent releases.

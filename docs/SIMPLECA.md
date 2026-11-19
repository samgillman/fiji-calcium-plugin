# SimpleCa²⁺ inside Fiji

Version 2026.11.19.8 adds native analysis for one recording at a time. Install it using **Install update**, save your work and restart Fiji. Existing installers can install this update: the R analysis source and its license are embedded in the plugin script.

## First use

1. Choose cells and **Measure and save**. If you edit the selected outlines afterward, measure again before opening Step 4.
2. Open **Analyze responses**, then **Set up analysis engine**. This downloads a checksum-verified micromamba launcher and installs R 4.5.2 from conda-forge under your user account. It supports 64-bit Windows, Intel/Apple Silicon Mac and x64 Linux. No terminal, system R installation or administrator access is required.
3. Confirm **frames per second** against the acquisition settings. Recognized timing metadata is suggested; missing metadata leaves this field blank. Frame 1 is time zero. Baseline frame numbers start at 1 and include both endpoints.
4. Enter the stimulus window, its preceding baseline duration, and response thresholds. Optionally enter a separate KCl window. The prestimulus baseline must contain at least three samples; each response window must fit the recording and have enough consecutive samples. Invalid coverage appears as unscorable with a reason.
5. Check the confirmation box and **Analyze responses**. Keep the window open until completion. **Cancel** stops the local job; **View log** shows diagnostics. A failed or cancelled run restores the controls and keeps the previous completed result available.

The first setup needs internet access and disk space for R and downloaded packages; it may take several minutes. Later analyses run locally without internet. Files live under `.calcium-simpleca` in your home folder (`C:\Users\your-name` on typical Windows installations). Setup is separate from Cellpose and does not alter it. Incomplete setup never marks R ready; retry uses the package cache. Partial setup directories may remain after cancellation. Network restrictions must allow GitHub release downloads and conda-forge packages.

## Review and export

- **Traces:** mean ± SEM of cells with valid initial F0, plus the selected cell. Gray marks the normalization baseline and pale orange the stimulus window. SEM is a cell-level summary, not uncertainty across animals.
- **Heatmap:** signed ΔF/F₀, with blue negative, white zero, red positive and gray unavailable. Click a row to inspect that cell's trace.
- **Metrics + quality:** every registered cell, including excluded and unscorable records with reasons. A valid initial F0 does not guarantee an evaluable stimulus window; consult the cell response calls too.
- **Cell response calls / Responder counts:** SimpleCa's stimulus/KCl classifications. Percent responding is Responded / Evaluable, with Registered, Excluded and Unscorable shown separately. The population is the cells you selected; selecting only cells that brighten can bias that population.
- **Cell review image:** click an outline, choose a cell in the trace list, or select its table row. This image is a snapshot, separate from the Step 2 editing controls.

Each successful run saves a new `SimpleCa-date-time-id` folder next to the measurement outputs. It contains the exact unlifted input CSV, normalized traces, mean/SEM, metrics, quality records, responder calls/counts, settings, manifest and analysis log. **Save plots** exports the current trace and full heatmap as PNG. Previous runs are not overwritten. Partial runs have no `COMPLETE` marker and are not loaded as results.

## Calculations and provenance

The bundled `utils.R` and `responders.R` are unchanged from [SimpleCa²⁺ commit a472f246](https://github.com/samgillman/simplecaapp/tree/a472f2460075343d5b97c8bac32ddf3da0c77386). The adapter uses the same raw-fluorescence normalization, metric functions and responder functions. Source hashes, plugin version, R version, exact input checksum and selected settings are saved in the manifest. Runtime dependency records remain in the managed environment's `conda-meta` directory.

The input is `gray_values_background_corrected_for_SimpleCa.csv`: background was subtracted once by Fiji. No second background subtraction or intensity lift is applied. Initial F0 is the mean of finite samples in the selected baseline; nonfinite F0 or F0 ≤ 0.000001 is excluded, matching the pinned SimpleCa version. Missing samples remain missing. Normalized traces are `(F−F0)/F0`.

Metric AUC is the signed net integral over the whole recording. Time-to-peak is the recording time of the peak, not latency after stimulus. Censored widths remain censored. Responder calls use the local prestimulus reference and the existing amplitude/noise and consecutive-sample rules. The legacy intensity-lifted export is still written by Step 3 for compatibility but is never used in Step 4.

This first integration does not reproduce the web application's full multi-recording comparisons, grouping controls or publication figure customization. It adds local analysis and review without changing detection or the original image measurements.

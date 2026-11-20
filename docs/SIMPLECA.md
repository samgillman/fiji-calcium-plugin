# SimpleCa²⁺ inside Fiji

Version 2026.11.19.8 adds native analysis for one recording at a time. Install it using **Install update**, save your work and restart Fiji. Existing installers can install this update: the R analysis source and its license are embedded in the plugin script.

## Where Step 4 opens

**Analyze responses** opens Step 4 as a page in your default browser. The plugin serves that page from your own computer: it is reachable only from this computer, its address carries a long random key, nothing is uploaded and no internet is needed (except for the one-time engine setup). Keep the plugin window in Fiji open while you use the page; closing it stops the page. If you measure again in Fiji, the page reloads with the new cells.

Everything below works the same on the page and in the Fiji window. On the page you can also hover a trace to see which cell it is, click a trace to select that cell, scroll over the plot to zoom the time axis (double-click to reset), sort any table by clicking a column name, and hide the settings column. **Use a Fiji window instead** opens the same analysis in a Fiji window; the plugin also falls back to it if a browser page cannot be opened.

## First use

1. Choose cells and **Measure and save**. If you edit the selected outlines afterward, measure again before opening Step 4.
2. Open **Analyze responses**, then **Set up analysis engine**. This downloads a checksum-verified micromamba launcher and installs R 4.5.2 from conda-forge under your user account. It supports 64-bit Windows, Intel/Apple Silicon Mac and x64 Linux. No terminal, system R installation or administrator access is required.
3. Confirm **frames per second** against the acquisition settings. Recognized timing metadata is suggested; missing metadata leaves this field blank. Frame 1 is time zero. Baseline frame numbers start at 1 and include both endpoints.
4. Enter the stimulus window, its preceding baseline duration, and response thresholds. If KCl was added after the stimulus as a positive control, tick **Include a KCl control** and enter when it was on; each cell is then also scored for that window. The prestimulus baseline must contain at least three samples; each response window must fit the recording and have enough consecutive samples. Invalid coverage appears as unscorable with a reason.
5. Check the confirmation box and **Analyze responses**. Keep the window open until completion. **Cancel** stops the local job; **View log** shows diagnostics. A failed or cancelled run restores the controls and keeps the previous completed result available.

The first setup needs internet access and disk space for R and downloaded packages; it may take several minutes. Later analyses run locally without internet. Files live under `.calcium-simpleca` in your home folder (`C:\Users\your-name` on typical Windows installations). Setup is separate from Cellpose and does not alter it. Incomplete setup never marks R ready; retry uses the package cache. Partial setup directories may remain after cancellation. Network restrictions must allow GitHub release downloads and conda-forge packages.

## Review and export

- **Time and the movie:** the slider under the Traces tab, the **Play** button and a click or drag on the plot all move one time cursor. Choose **Movie: change from rest** above the picture to see each frame minus the average of the baseline frames, lightly smoothed, so a responding cell lights up as its trace rises; **Movie: as recorded** shows the frame itself and **Colour picture** the Step 2 picture. **Zoom to cell** follows the selected cell; untick **All outlines** to see the cells themselves. These views are for looking only: the recording is read, never changed, and nothing shown here is measured or saved.
- **Every trace:** all listed cells are drawn behind the mean, coloured by response call, with the selected cell in black on top. The y-axis covers 99 % of all values (plus the mean and the selected cell), so rarer extremes run off the edge. Untick **Every trace** to see only the mean and the selected cell.
- **Show:** list and draw only the cells with one response call. Clicking a headline tile does the same. The line under the list gives the selected cell's call, the reason, its largest rise and the threshold it needed, from the Cell response calls table (fractions of the fluorescence just before the stimulus). The left and right arrow keys step through the listed cells.
- **Traces:** mean ± SEM of cells with valid initial F0, plus the selected cell. Gray marks the normalization baseline and pale orange the stimulus window. SEM is a cell-level summary, not uncertainty across animals.
- **Heatmap:** signed ΔF/F₀, with blue negative, white zero, red positive and gray unavailable. Rows are ordered responders first with the largest rise on top, then non-responders, then cells without a call (or as measured, from the **Order** list); the strip on the left shows each cell's call. Colours stop at a limit that covers 99 % of the values, so a single extreme cell cannot wash out the others; the key shows that limit. Click a row to inspect that cell's trace.
- **Metrics + quality:** every registered cell, including excluded and unscorable records with reasons. A valid initial F0 does not guarantee an evaluable stimulus window; consult the cell response calls too.
- **Cell response calls / Responder counts:** SimpleCa's stimulus/KCl classifications. Percent responding is Responded / Evaluable, with Registered, Excluded and Unscorable shown separately. The population is the cells you selected; selecting only cells that brighten can bias that population.
- **Cell picture (under the trace):** click an outline, choose a cell in the list or with the previous/next buttons, or select its table row. After an analysis, thick white outlines responded, gray ones did not, and dashed ones were excluded or unscorable. Scroll to zoom, drag to move, double-click to see the whole picture. It is a snapshot, separate from the Step 2 editing controls.
- **Headline tiles:** percent responding with its evaluable denominator, and the registered, excluded and unscorable counts, read from the saved Responder counts table.
- **Settings boxes:** suggested values are shown rounded; hover over a box to see the exact value. A box you have not edited sends that exact value to the analysis.

Each successful run saves a new `SimpleCa-date-time-id` folder next to the measurement outputs. It contains the exact unlifted input CSV, normalized traces, mean/SEM, metrics, quality records, responder calls/counts, settings, manifest and analysis log. **Save plots** exports the current trace and full heatmap as PNG. Previous runs are not overwritten. Partial runs have no `COMPLETE` marker and are not loaded as results.

## Calculations and provenance

The bundled `utils.R` and `responders.R` are unchanged from [SimpleCa²⁺ commit a472f246](https://github.com/samgillman/simplecaapp/tree/a472f2460075343d5b97c8bac32ddf3da0c77386). The adapter uses the same raw-fluorescence normalization, metric functions and responder functions. Source hashes, plugin version, R version, exact input checksum and selected settings are saved in the manifest. Runtime dependency records remain in the managed environment's `conda-meta` directory.

The input is `gray_values_background_corrected_for_SimpleCa.csv`: background was subtracted once by Fiji. No second background subtraction or intensity lift is applied. Initial F0 is the mean of finite samples in the selected baseline; nonfinite F0 or F0 ≤ 0.000001 is excluded, matching the pinned SimpleCa version. Missing samples remain missing. Normalized traces are `(F−F0)/F0`.

Metric AUC is the signed net integral over the whole recording. Time-to-peak is the recording time of the peak, not latency after stimulus. Censored widths remain censored. Responder calls use the local prestimulus reference and the existing amplitude/noise and consecutive-sample rules. The legacy intensity-lifted export is still written by Step 3 for compatibility but is never used in Step 4.

This first integration does not reproduce the web application's full multi-recording comparisons, grouping controls or publication figure customization. It adds local analysis and review without changing detection or the original image measurements.

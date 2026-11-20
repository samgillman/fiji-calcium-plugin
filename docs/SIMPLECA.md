# SimpleCa²⁺ inside Fiji

Step 4 analyzes one measured recording locally. Version 2026.11.20.3 combines the browser/movie workflow from version 2026.11.20.2 with portable response protocols, exact timestamps, persistent ROI identities and explicit analysis provenance. The browser page and native Fiji fallback use the same pinned R engine.

## Where Step 4 opens

**Open analysis page** opens Step 4 in your default browser. Fiji serves the page on your own computer, through a loopback address with a random session key. No recording data is uploaded. Keep Fiji open while using the page. Each measurement gets its own page and private input snapshot; remeasuring or closing the cell-selection window does not shut down earlier review pages. Use **Close analysis** in a page to release its session and temporary inputs, or quit Fiji. Browser tabs begin with the cell count and include a session identifier. Engine setup and update downloads require internet; an already installed engine and the local page can run offline.

The browser adds trace hover/selection and time-axis zoom, sortable tables, and an interactive heatmap. **Use a Fiji window instead** opens the native analysis window, which is also the fallback when a browser cannot be opened. Both views include recording playback, frame scrubbing, linked traces and a heatmap, and preserve the protocol, QC and measurement definitions below.

If you change the accepted cells or measure again, review the new measurement before analyzing it. A completed result remains associated with its own input, ROI and protocol snapshots; draft settings or a newer measurement must not be read as the settings of an older result.

## Preview the recording before analysis

When the recording is available, the browser opens in **Recording** view. The playback controls are at the top of the review area: **Play**, previous/next frame, the frame slider and **Speed**. They work before an analysis has been run and do not require the R engine to be installed. The page says **Preview · no analysis yet**; response-percentage tiles, trace results and result tables appear only after a completed analysis.

Playback follows the available frame times: at **1×**, it schedules frames using their acquisition intervals; other speeds scale those intervals. Before timing has been confirmed, estimated seconds are marked as provisional. With no valid time base, the page explicitly labels a frame-by-frame preview whose **1×** rate is **20 frames per second**. That fallback is a display rate, not an inferred acquisition rate and not timing suitable for analysis. Frame loading can slow playback; the frame/time label identifies what is being reviewed.

**View → Static colour overview** shows the Step 2 overview, which does not move. Play and speed are disabled in that view. Before results exist, the slider and frame-step buttons are also disabled; after analysis, they can inspect the trace cursor while the overview remains static. Choose **Recording** to play or scrub actual movie frames. **Change from baseline** displays baseline-subtracted frames. These viewing choices do not alter measured values or response calls.

## Choose cells and measure

Choose identifiable cells independently of whether they brighten. **Review baseline frames** displays a user-confirmed baseline frame range for reviewing outlines; **Restore all-frames view** restores the original colour picture. The response-assisted magenta display is an optional review aid and starts disabled. The selection basis is recorded. Cellpose suggestions and grades still require review; accepting only brightening cells would change the population whose response percentage is reported.

**Measure and save** writes stable IDs into the accepted ROIs and their measurement columns. An accepted-ROI registry records each ID, geometry hash, measured status, edge exclusion, background pixel count and selection basis. A separate accepted-ROI ZIP includes accepted outlines excluded before measurement. Those exclusions remain visible in the cohort flow; they are not silently relabeled as nonresponders. Editing any accepted outline, edge policy or selection basis requires measuring again before opening Step 4.

The analysis input is the full-precision `gray_values_background_corrected_native.csv`, containing the original cell mean minus its local background once. The older three-decimal CSVs remain available for compatibility. The intensity-lifted `for_SimpleCa_every_cell_kept.csv` is a legacy diagnostic export: an arbitrary offset changes fractional amplitudes and thresholds, so do not use it for ΔF/F₀ responder percentages.

## First use and timing

1. Click **Open analysis page**, then **Set up analysis engine**. A checksum-verified micromamba launcher installs private R 4.5.2 and jsonlite 2.0.0 from conda-forge. Supported runtime targets are 64-bit Windows, Intel/Apple Silicon Mac and x64 Linux. No system R installation or administrator access is required.
2. Confirm **frames per second** against acquisition metadata, or use **Load timestamps CSV**. Constant FPS puts frame 1 at zero; it cannot identify dropped frames. A timestamp file needs a `Time_s` column and exactly one finite, nonnegative, strictly increasing value per measured frame. Optional `Frame` must be 1, 2, … in order. The timestamp file is snapshotted and its checksum is checked before analysis.
3. Set the inclusive first and last **normalization baseline frames**. Suggested onset follows the final baseline frame; always replace it with your actual application timing. The page warns when the fixed F₀ baseline includes onset or later samples. This fixed baseline creates the global trace `(F−F0)/F0`. The local prestimulus baseline in each protocol row is a separate reference for that response call.
4. Choose or load a protocol, inspect each row and any imported field timing overrides, and enter condition and independent biological-unit IDs when known. Check the timing/selection confirmation and click **Run analysis** in the browser, or **Analyze responses** in the native fallback.

Setup needs internet and disk space for R and downloaded packages. It is separate from Cellpose and stored under `.calcium-simpleca` in your home folder. A previous base-only R environment is not reused as protocol-ready: setup verifies the required R and jsonlite versions. Retry uses the package cache, and later analyses can run offline. Cancelled setup may leave partial installation directories, but it never marks them ready.

## Response protocols

Each stimulus has its own response rule: the primary stimulus and, optionally, one named positive control. ATP, capsaicin, KCl, ionomycin or another control name is metadata; the software does not infer stimulation concentration or an optimal threshold from that name. This version supports one primary window and at most one control window. The browser keeps timing and a rule summary visible; expand **Adjust the response rule** for detailed threshold, persistence and noise settings. **Protocol name, files and overrides** and **Experiment details and diagnostics** also expand when needed. In the native fallback, scroll the protocol table horizontally to inspect every setting. Unedited suggested values retain their underlying precision.

| Setting | Meaning |
|---|---|
| Baseline start/end | Local reference interval `[start,end)`, ending at or before application onset |
| Applied onset | Recorded stimulus application time; the reference for latency |
| Response start/end | Closed scoring interval `[start,end]`, beginning at or after onset |
| Minimum rise | Increase above the local baseline. The browser's percentage field displays `10` for 10%; the protocol stores this as the fraction `0.10` |
| Baseline SD multiplier | Sample SD of the local baseline after expressing it as a fraction of the local reference |
| Persistence | Actual timestamp span in `seconds`, or consecutive sample count in `samples` |
| Maximum gap | Largest consecutive timestamp gap allowed; blank uses the engine's declared automatic gap rule |
| Minimum baseline samples | Required finite samples in the local reference interval |
| Noise floor | Optional minimum baseline noise. Browser percentage fields are converted to the protocol's fractional units; this setting does not exclude a cell |
| Expected interval | Optional declared sampling interval used by coverage/gap checks |

The threshold is `max(minimum rise, SD multiplier × local baseline SD)`, with the optional noise floor applied before multiplying. Samples equal to the threshold qualify (`>=`). Only the setting for the selected persistence mode is editable; unused values are retained for replay. Seconds mode uses elapsed span and its inherent minimum of two observations, without adding the sample-mode count as a second requirement. There is no smoothing. The local trace is `(y−b)/(1+b)`, where `y` is the fixed-reference fractional trace and `b` its local baseline mean. Therefore a local threshold `T` is drawn on the global trace at `b+(1+b)×T`.

**Provisional slow evoked: 1 s** starts at a 0.10 fractional floor, 3 baseline SDs and at least a 1-second qualifying span. It is a starting protocol to validate on the experiment's controls, not an established optimum for 50B11 or every indicator. A one-second span at 2 Hz needs at least three qualifying samples: last timestamp minus first timestamp is the duration. Missing values and excessive gaps break the run. Version 2 requires complete window coverage and enough reference samples; unavailable coverage is unscorable with a reason, rather than treated as no response.

**Custom version 2** keeps all settings independently editable. Timing edits retain different stimulus/control threshold, duration and noise rules. **Legacy sample rule (starts at 2)** exports actual schema 1 and calls the historical responder engine, including its historical endpoint handling. It is not version 2 with a two-sample setting. Legacy mode uses an immediately preceding reference and its optional control remains KCl. Version-2-only fields are locked or labeled unused; editing the legacy scoring start also moves its reference end and onset. Historical settings TSV files still follow the same legacy path. New explanatory metrics/QC annotations do not change historical classifications or counts.

The browser's **Import protocol JSON / Export protocol**, or the native **Load protocol / Save protocol**, exchange declarative JSON with SimpleCa²⁺. Files record schema and engine version, name, validation status, each window and optional field overrides. Imported overrides are retained in the browser's expandable protocol settings; the native window also shows a recording match count and offers **View field timing overrides**. Overrides can be explicitly cleared. A matching override changes that recording's timing while preserving its rule. `validated_by_user` is a user declaration, not certification by the program. Unsupported schema, reference, engine or rule semantics are rejected.

## Eligibility, quality and percentages

A cell with no finite positive normalization reference, or F0 ≤ 0.000001, is excluded using the shared normalization rule. A missing local background is a technical exclusion. Missing/insufficient local baseline data, incomplete response coverage or an invalid local reference can make an individual window unscorable. The cell remains in the registry and quality/call exports.

Low F0 relative to raw baseline SD, zero baseline SD and other review flags are warnings. They do not automatically remove cells. There is no automatic peak-outlier filter. Zero SD can leave the amplitude floor as the active requirement; the warning remains available for review. A noise-only rule with zero estimated SD is unscorable unless a justified nonzero measurement-noise floor is supplied.

For each window, percent responding is `100 × Responded / (Responded + Did not respond)`. That denominator is **Evaluable**. Registered, Excluded and Unscorable are exported separately. Zero evaluable cells produces an unavailable percentage, not zero. The optional control-positive population includes only control responders that are also evaluable for the primary stimulus; its numerator is the primary responders within that population. A control nonresponse is not automatically classified as a dead cell or removed from the primary percentage.

`cohort_flow.csv` distinguishes all accepted ROIs, premeasurement exclusions, cells registered for analysis and normalization exclusions. `responder_counts.csv` adds window-specific evaluability. Selection based on response cannot be repaired by a threshold after the fact.

## Review the completed analysis

- **Traces:** fixed-reference ΔF/F₀ for the selected cell and the mean ± SEM of cells. The browser can show every cell behind the mean, coloured by its saved response call. Hover or click to identify a cell; scroll to zoom the time axis and double-click to reset. A view range covering 99% of values is a display choice, not a cell exclusion or a change to peak/AUC measurements.
- **References and qualifying runs:** normalization, local baseline, response windows and qualifying evidence refer to the completed analysis. The local threshold is transformed to the global trace axis as described above. The selected call supplies its status, reason, threshold and run timing; editing a draft must not redraw an old result using new thresholds.
- **Time and movie:** in **Recording** view, the slider, **Play** and the plot cursor move together. **Change from baseline** shows that frame minus the average of the normalization baseline frames, lightly smoothed for viewing. **Static colour overview** shows the Step 2 picture and stays still; movie playback is disabled in that view. The native fallback labels these views **Movie: as recorded**, **Movie: change from rest** and **Static colour overview**. **Zoom to cell** follows a selected cell; **All outlines** controls the overlay. These operations affect display only, not measured fluorescence, normalization or response calls.
- **Cells and response groups:** click a trace, picture outline or heatmap row, choose a cell in the list, or use previous/next navigation. The browser's **Show** filter and headline tiles restrict the displayed group. They do not change the registered/evaluable denominator. Outlines and status colours use the saved calls for the selected stimulus window.
- **Heatmap:** choose signed blue–white–red values or an available nonnegative palette, and order by response call, time to peak, peak amplitude or original cell order. The 99% colour range makes the display less sensitive to extreme values; choose the full range to see the largest magnitude. Missing values and acquisition gaps remain visible. Changing palette, range or order does not drop cells from the result tables.
- **Headline tiles:** the saved primary-window percentage and its evaluable denominator, plus registered, excluded and unscorable counts. The stimulus name identifies that primary summary even when another window is being inspected. The saved tables remain the authoritative record.
- **Quality, ROI registry and call tables:** reasons and warnings remain beside each cell. Whole-recording metrics and stimulus-window metrics are separate. The browser supports table sorting; a selected row can return to that cell's trace. The native fallback retains trace, movie/picture, heatmap and table review without changing Step 2 ROIs.

Editing a draft protocol or settings does not change an already completed result. The browser shows **Changes pending** and offers **Apply changes & reanalyze**; results, movie timing and exports continue to use the last completed analysis until the new run succeeds. Timestamp and diagnostic file content signatures include changes or deletion at the same path. Displays and plot exports remain bound to the completed analysis, with provenance identifying its input and applied protocol. The native **Save plots as PNG** also writes `plot_provenance.json`, containing the manifest checksum, protocol identity, selected cell/window and whether edits were pending. Failed or cancelled reruns keep the preceding completed tables. Partial runs have no `COMPLETE` marker and are not displayed as completed results.

## Metrics, diagnostics and biological units

Legacy **AUC** is the signed whole-recording integral and legacy **Time-to-Peak** is the recording timestamp of that peak. Those columns are not repurposed. New window metrics include local fractional peak, signed local AUC in fraction·seconds, actual covered duration/fraction, response latency from application, peak latency, run duration, run samples and confirmation time. AUC uses observed valid sample intervals and does not interpolate through excessive gaps or invent signal outside the recording. Peak outliers remain visible. Legacy annotations use the actual historically clipped endpoint and explicitly identify the annotation engine.

**Save diagnostic templates** creates a sensitivity grid and a negative-control-window template in the measurement folder. Enter genuine unstimulated/vehicle times in the latter. Loading a grid reruns the shared engine for each declared scenario; loading negative-control windows reports their apparent response fraction. These are separately labeled exploratory diagnostics, not replacements for the applied rule or a claim of statistical specificity. They require version 2; selecting/importing legacy detaches diagnostic file selections and explains why. Files are snapshotted for reproducibility.

A biological unit is the independent culture, animal or preparation. Multiple cells and fields from it are nested observations. The recording's unit and condition are exported with counts; the shared summary function pools fields within an assigned unit and reports unit percentages separately from cell-pooled percentages. Step 4 handles one recording per run, so cross-recording aggregation belongs in SimpleCa²⁺ or a subsequent analysis using the exported metadata. Blank unit IDs create an explicit unassigned record and no invented biological N. Trace SEM remains descriptive cell SEM.

## Files and source provenance

Each successful run gets a new `SimpleCa-date-time-id` folder. It contains the input and optional timestamp/ROI/diagnostic snapshots, settings, `protocol.json`, `applied_protocol.json`, normalized traces, mean trace, metrics, quality, ROI registry, cohort flow, calls/counts, window metrics, replicate/group counts, metadata, manifest and log. Optional diagnostic outputs are separate. Previous runs are not overwritten.

The bundled `utils.R`, `response_protocols.R` and `responders.R` are unchanged from [SimpleCa²⁺ commit b92d4ce9](https://github.com/samgillman/simplecaapp/tree/b92d4ce98ddab3fceac352179703c1f48f044c90). The adapter uses their normalization, metrics, protocol import/validation, response calls and summary functions. Source SHA-256 hashes, plugin/R/jsonlite versions, actual input checksum, applied protocol and units are saved with results. The shared canonical protocol MD5 is an identity for comparison, not a security signature. Managed package records remain in the environment's `conda-meta` directory.

Synthetic tests cover independent hand-calculated response outcomes as well as parity with the shared engine. CI is configured for managed runtime checks on Windows, macOS and Linux and native GUI checks on Windows/Linux. Browser acceptance is a separate Playwright gate described in the release guide; server and JavaScript syntax checks alone are not browser validation. These checks do not establish biological sensitivity/specificity or replace acquisition-specific control validation.

## Step 4 editing and quality review

Invalid frame rates, frame ranges, reversed windows and negative thresholds show inline errors before analysis; submission runs the authoritative validation again. A timestamp CSV takes priority over FPS. Elapsed-time and legacy rule types retain separate drafts when switching. **Discard pending changes** restores the complete applied draft, including uploaded file contents, or the starting draft before the first analysis. Completed results, movie timing and exported plots keep their applied settings until a new analysis succeeds.

A newly enabled positive control has blank timing fields: enter its actual baseline, addition and scoring times. The software cannot infer control timing from the end of the stimulus window. Recording-coverage and baseline-overlap warnings do not silently alter a protocol.

The quality summary and per-cell feedback show saved flags, including `DURATION_BELOW_SAMPLE_INTERVAL` and `reference_below_3_baseline_sd_review`. A separate review indicator identifies saved traces below −1 and legacy rise times whose 10%-peak level is below baseline SD. These display diagnostics do not recalculate response calls or add exclusions. **Needs quality review**, **Excluded** and **Unscorable** are separate filters. A low reference can make fractional values unstable; a value below −1 corresponds to negative background-corrected fluorescence when F₀ is positive. Inspect background and baseline stability before interpreting these traces.

**Response-window metrics** provide the applied per-stimulus peak, latency and signed AUC. **Legacy whole-recording metrics** retain their established definitions: `AUC` integrates the whole recording and can be negative under downward drift; `Time_to_Peak` is the recording timestamp; `Rise_Time` uses 10–90% crossings after the display baseline, which can be sensitive to baseline noise. The interface labels these definitions rather than changing old values. Use `Peak_latency_s` for time from stimulus onset and `Window_AUC_fraction_s` for the declared response window.

The calls table starts with core response/quality columns and a pinned cell-name column. **Show all saved columns** exposes the complete table; saved CSVs always retain every column. Empty tables explain that no rows were saved. Short windows can scroll vertically while plots and captions keep a readable height.

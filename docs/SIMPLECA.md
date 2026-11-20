# SimpleCa²⁺ inside Fiji

Step 4 analyzes one measured recording locally. Version 2026.11.20.2 adds portable response protocols, timestamp support, persistent ROI identities and explicit analysis provenance while retaining the cell picture, result tiles and navigation from version 2026.11.20.1.

## Choose cells and measure

Choose identifiable cells independently of whether they brighten. **Review baseline frames** displays a user-confirmed baseline frame range for reviewing outlines; **Restore all-frames view** restores the original colour picture. The response-assisted magenta display is an optional review aid and starts disabled. The selection basis is recorded. Cellpose suggestions and grades still require review; accepting only brightening cells would change the population whose response percentage is reported.

**Measure and save** writes stable IDs into the accepted ROIs and their measurement columns. An accepted-ROI registry records each ID, geometry hash, measured status, edge exclusion, background pixel count and selection basis. A separate accepted-ROI ZIP includes accepted outlines excluded before measurement. Those exclusions remain visible in the cohort flow; they are not silently relabeled as nonresponders. Editing any accepted outline, edge policy or selection basis requires measuring again before opening Step 4.

The analysis input is the full-precision `gray_values_background_corrected_native.csv`, containing the original cell mean minus its local background once. The older three-decimal CSVs remain available for compatibility. The intensity-lifted `for_SimpleCa_every_cell_kept.csv` is a legacy diagnostic export: an arbitrary offset changes fractional amplitudes and thresholds, so do not use it for ΔF/F₀ responder percentages.

## First use and timing

1. Open **Analyze responses**, then **Set up analysis engine**. A checksum-verified micromamba launcher installs private R 4.5.2 and jsonlite 2.0.0 from conda-forge. Supported runtime targets are 64-bit Windows, Intel/Apple Silicon Mac and x64 Linux. No system R installation or administrator access is required.
2. Confirm **frames per second** against acquisition metadata, or use **Load timestamps CSV**. Constant FPS puts frame 1 at zero; it cannot identify dropped frames. A timestamp file needs a `Time_s` column and exactly one finite, nonnegative, strictly increasing value per measured frame. Optional `Frame` must be 1, 2, … in order. The timestamp file is snapshotted and its checksum is checked before analysis.
3. Set the inclusive first and last **normalization baseline frames**. This fixed baseline creates the global trace `(F−F0)/F0`. The local prestimulus baseline in each protocol row is a separate reference for that response call.
4. Choose or load a protocol, inspect each row and any imported field timing overrides, and enter condition and independent biological-unit IDs when known. Check the timing/selection confirmation and click **Analyze responses**.

Setup needs internet and disk space for R and downloaded packages. It is separate from Cellpose and stored under `.calcium-simpleca` in your home folder. A previous base-only R environment is not reused as protocol-ready: setup verifies the required R and jsonlite versions. Retry uses the package cache, and later analyses can run offline. Cancelled setup may leave partial installation directories, but it never marks them ready.

## Response protocols

The native table exposes a separate rule for the primary stimulus and, optionally, one named positive control. ATP, capsaicin, KCl, ionomycin or another control name is metadata; the software does not infer stimulation concentration or an optimal threshold from that name. This version supports one primary row and at most one control row. Scroll horizontally to inspect every setting. Numeric display is compact, but unedited suggested values retain full precision; tooltips show the underlying value.

| Setting | Meaning |
|---|---|
| Baseline start/end | Local reference interval `[start,end)`, ending at or before application onset |
| Applied onset | Recorded stimulus application time; the reference for latency |
| Response start/end | Closed scoring interval `[start,end]`, beginning at or after onset |
| Minimum rise | Fractional increase above the local baseline; `0.10` means 10% |
| Baseline SD multiplier | Sample SD of the local baseline after expressing it as a fraction of the local reference |
| Persistence | Actual timestamp span in `seconds`, or consecutive sample count in `samples` |
| Maximum gap | Largest consecutive timestamp gap allowed; blank uses the engine's declared automatic gap rule |
| Minimum baseline samples | Required finite samples in the local reference interval |
| Noise floor | Optional minimum baseline noise in fractional units; does not exclude a cell |
| Expected interval | Optional declared sampling interval used by coverage/gap checks |

The threshold is `max(minimum rise, SD multiplier × local baseline SD)`, with the optional noise floor applied before multiplying. Samples equal to the threshold qualify (`>=`). Only the setting for the selected persistence mode is editable; unused values are retained for replay. Seconds mode uses elapsed span and its inherent minimum of two observations, without adding the sample-mode count as a second requirement. There is no smoothing. The local trace is `(y−b)/(1+b)`, where `y` is the fixed-reference fractional trace and `b` its local baseline mean. Therefore a local threshold `T` is drawn on the global trace at `b+(1+b)×T`.

**Provisional slow evoked: 1 s** starts at a 0.10 fractional floor, 3 baseline SDs and at least a 1-second qualifying span. It is a starting protocol to validate on the experiment's controls, not an established optimum for 50B11 or every indicator. A one-second span at 2 Hz needs at least three qualifying samples: last timestamp minus first timestamp is the duration. Missing values and excessive gaps break the run. Version 2 requires complete window coverage and enough reference samples; unavailable coverage is unscorable with a reason, rather than treated as no response.

**Custom version 2** keeps all settings independently editable. Timing edits retain different stimulus/control threshold, duration and noise rules. **Legacy sample rule (starts at 2)** exports actual schema 1 and calls the historical responder engine, including its historical endpoint handling. It is not version 2 with a two-sample setting. Legacy mode uses an immediately preceding reference and its optional control remains KCl. Version-2-only fields are locked or labeled unused; editing the legacy scoring start also moves its reference end and onset. Historical settings TSV files still follow the same legacy path. New explanatory metrics/QC annotations do not change historical classifications or counts.

**Load protocol / Save protocol** exchange declarative JSON with SimpleCa²⁺. Files record schema and engine version, name, validation status, each window and optional field overrides. Imported field overrides are retained, shown with a recording match count and available through **View field timing overrides**; they can be explicitly cleared. A matching override changes that recording's timing while preserving its rule. `validated_by_user` is a user declaration, not certification by the program. Unsupported schema, reference, engine or rule semantics are rejected.

## Eligibility, quality and percentages

A cell with no finite positive normalization reference, or F0 ≤ 0.000001, is excluded using the shared normalization rule. A missing local background is a technical exclusion. Missing/insufficient local baseline data, incomplete response coverage or an invalid local reference can make an individual window unscorable. The cell remains in the registry and quality/call exports.

Low F0 relative to raw baseline SD, zero baseline SD and other review flags are warnings. They do not automatically remove cells. There is no automatic peak-outlier filter. Zero SD can leave the amplitude floor as the active requirement; the warning remains available for review.

For each window, percent responding is `100 × Responded / (Responded + Did not respond)`. That denominator is **Evaluable**. Registered, Excluded and Unscorable are exported separately. Zero evaluable cells produces an unavailable percentage, not zero. The optional control-positive population includes only control responders that are also evaluable for the primary stimulus; its numerator is the primary responders within that population. A control nonresponse is not automatically classified as a dead cell or removed from the primary percentage.

`cohort_flow.csv` distinguishes all accepted ROIs, premeasurement exclusions, cells registered for analysis and normalization exclusions. `responder_counts.csv` adds window-specific evaluability. Selection based on response cannot be repaired by a threshold after the fact.

## Review the completed analysis

- **Traces:** fixed-reference ΔF/F₀ for the selected cell and the mean ± SEM of cells. Shading distinguishes normalization, local baseline, response window and a qualifying run. The selected call supplies the local reference, floor, noise threshold and final threshold, transformed to the plotted global axis. Long gaps are not joined across the applied gap limit. The call reason and qualifying duration/sample count appear below the trace.
- **Heatmap:** signed global ΔF/F₀, blue below zero and red above, with gray unavailable values. Position follows actual timestamps; long acquisition gaps stay blank. Click a row to inspect its trace.
- **Cell picture:** an embedded snapshot with stable measured cell IDs. Thick white outlines responded, gray ones did not, and dashed ones were excluded or unscorable for the selected window. Click a cell, use previous/next, scroll to zoom, drag to move, or double-click to reset. This does not edit Step 2 ROIs.
- **Headline tiles:** the saved primary-window percentage and its evaluable denominator, plus registered, excluded and unscorable counts. The saved tables remain the authoritative record.
- **Quality, ROI registry and call tables:** reasons and warnings remain beside each cell. Whole-recording metrics and stimulus-window metrics have separate tabs and exports.

Editing the draft protocol or settings does not change an already completed result. A visible pending-settings message tells you to rerun. Timestamp and diagnostic file content signatures are checked, including changes or deletion at the same path. Display and **Save plots as PNG** remain bound to the completed analysis; `plot_provenance.json` records its manifest checksum, protocol identity, selected cell/window and whether edits were pending. Failed/cancelled reruns keep the preceding completed tables. Partial runs have no `COMPLETE` marker and are not displayed as completed results.

## Metrics, diagnostics and biological units

Legacy **AUC** is the signed whole-recording integral and legacy **Time-to-Peak** is the recording timestamp of that peak. Those columns are not repurposed. New window metrics include local fractional peak, signed local AUC in fraction·seconds, actual covered duration/fraction, response latency from application, peak latency, run duration, run samples and confirmation time. AUC uses observed valid sample intervals and does not interpolate through excessive gaps or invent signal outside the recording. Peak outliers remain visible. Legacy annotations use the actual historically clipped endpoint and explicitly identify the annotation engine.

**Save diagnostic templates** creates a sensitivity grid and a negative-control-window template in the measurement folder. Enter genuine unstimulated/vehicle times in the latter. Loading a grid reruns the shared engine for each declared scenario; loading negative-control windows reports their apparent response fraction. These are separately labeled exploratory diagnostics, not replacements for the applied rule or a claim of statistical specificity. They require version 2; selecting/importing legacy detaches diagnostic file selections and explains why. Files are snapshotted for reproducibility.

A biological unit is the independent culture, animal or preparation. Multiple cells and fields from it are nested observations. The recording's unit and condition are exported with counts; the shared summary function pools fields within an assigned unit and reports unit percentages separately from cell-pooled percentages. Native Step 4 handles one recording per run, so cross-recording aggregation belongs in SimpleCa²⁺ or a subsequent analysis using the exported metadata. Blank unit IDs create an explicit unassigned record and no invented biological N. Trace SEM remains descriptive cell SEM.

## Files and source provenance

Each successful run gets a new `SimpleCa-date-time-id` folder. It contains the input and optional timestamp/ROI/diagnostic snapshots, settings, `protocol.json`, `applied_protocol.json`, normalized traces, mean trace, metrics, quality, ROI registry, cohort flow, calls/counts, window metrics, replicate/group counts, metadata, manifest and log. Optional diagnostic outputs are separate. Previous runs are not overwritten.

The bundled `utils.R`, `response_protocols.R` and `responders.R` are unchanged from [SimpleCa²⁺ commit b92d4ce9](https://github.com/samgillman/simplecaapp/tree/b92d4ce98ddab3fceac352179703c1f48f044c90). The adapter uses their normalization, metrics, protocol import/validation, response calls and summary functions. Source SHA-256 hashes, plugin/R/jsonlite versions, actual input checksum, applied protocol and units are saved with results. The shared canonical protocol MD5 is an identity for comparison, not a security signature. Managed package records remain in the environment's `conda-meta` directory.

Synthetic tests cover independent hand-calculated response outcomes as well as parity with the shared engine. CI is configured for managed runtime checks on Windows, macOS and Linux and native GUI checks on Windows/Linux. These checks do not establish biological sensitivity/specificity or replace acquisition-specific control validation.

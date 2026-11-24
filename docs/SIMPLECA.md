# SimpleCa²⁺ inside Fiji

Step 4 analyzes one measured recording locally. Version 2026.11.20.3 combines the browser/movie workflow from version 2026.11.20.2 with portable response protocols, exact timestamps, persistent ROI identities and explicit analysis provenance. The browser page and native Fiji fallback use the same pinned R engine.

## Where Step 4 opens

**Open analysis page** opens Step 4 in your default browser. Fiji serves the page on your own computer, through a loopback address with a random session key. No recording data is uploaded. Keep Fiji open while using the page. Each measurement gets its own page and private input snapshot; remeasuring or closing the cell-selection window does not shut down earlier review pages. Use **Close analysis** in a page to release its session and temporary inputs, or quit Fiji. Browser tabs begin with the cell count and include a session identifier. Engine setup and update downloads require internet; an already installed engine and the local page can run offline.

The browser adds trace hover/selection and time-axis zoom, sortable tables, and an interactive heatmap. **Use a Fiji window instead** opens the native analysis window, which is also the fallback when a browser cannot be opened. Both views include recording playback, frame scrubbing, linked traces and a heatmap, and preserve the protocol, QC and measurement definitions below.

If you change the accepted cells or measure again, review the new measurement before analyzing it. A completed result remains associated with its own input, ROI and protocol snapshots; draft settings or a newer measurement must not be read as the settings of an older result.

## Preview the recording before analysis

When the recording is available, the browser opens in **Recording** view. The playback controls are at the top of the review area: **Play**, previous/next frame, the frame slider and **Speed**. They work before an analysis has been run and do not require the R engine to be installed. The page says **Preview · no analysis yet**; response-percentage tiles, trace results and result tables appear only after a completed analysis.

Playback follows the available frame times on a clock in recording seconds: at **1×** a frame is shown when its timestamp is reached, and the other speeds (0.5× to 50×, **10×** by default) scale that clock. Above 30 frames per second on screen, every second, third or later frame is shown, and the frame note says so. A long gap between recorded frames waits at most about one second. Frames are loaded a few seconds ahead and kept, so scrubbing back over frames already seen is immediate; if loading falls behind, the picture holds and the note says **buffering** instead of the clock running ahead. While playing, the trace cursor moves continuously with the clock. Before timing has been confirmed, estimated seconds are marked as provisional. With no valid time base, the page explicitly labels a frame-by-frame preview whose **1×** rate is **20 frames per second**. That fallback is a display rate, not an inferred acquisition rate and not timing suitable for analysis. The frame/time label identifies what is being reviewed.

**View → Static colour overview** shows the Step 2 overview, which does not move. Play and speed are disabled in that view. Before results exist, the slider and frame-step buttons are also disabled; after analysis, they can inspect the trace cursor while the overview remains static. Choose **Recording** to play or scrub actual movie frames. **Change from baseline** displays baseline-subtracted frames. These viewing choices do not alter measured values or response calls.

## Choose cells and measure

Choose identifiable cells independently of whether they brighten. In Step 2, **Picture** switches between **All frames** (the original colour picture), **Before response** (a frame range before the stimulus that you confirm) and **Response** (the response in magenta). Magenta is an optional review aid: it starts off, and choosing it looks for the response at that moment, as the first screen's option does. The picture changes in place, so outlines and clicks keep working. The selection basis is recorded. Cellpose suggestions and grades still require review; accepting only brightening cells would change the population whose response percentage is reported.

**Measure and save** writes stable IDs into the accepted ROIs and their measurement columns. An accepted-ROI registry records each ID, geometry hash, measured status, edge exclusion, background pixel count and selection basis. A separate accepted-ROI ZIP includes accepted outlines excluded before measurement. Those exclusions remain visible in the cohort flow; they are not silently relabeled as nonresponders. Editing any accepted outline, edge policy or selection basis requires measuring again before opening Step 4.

The analysis input is the full-precision `gray_values_background_corrected_native.csv`, containing the original cell mean minus its local background once. The older three-decimal CSVs remain available for compatibility. The intensity-lifted `for_SimpleCa_every_cell_kept.csv` is a legacy diagnostic export: an arbitrary offset changes fractional amplitudes and thresholds, so do not use it for ΔF/F₀ responder percentages.

## First use and timing

1. Click **Open analysis page**, then **Set up analysis engine**. A checksum-verified micromamba launcher installs private R 4.5.2 and jsonlite 2.0.0 from conda-forge. Supported runtime targets are 64-bit Windows, Intel/Apple Silicon Mac and x64 Linux. No system R installation or administrator access is required.
2. Confirm **frames per second** against acquisition metadata, or use **Load timestamps CSV**. Constant FPS puts frame 1 at zero; it cannot identify dropped frames. A timestamp file needs a `Time_s` column and exactly one finite, nonnegative, strictly increasing value per measured frame. Optional `Frame` must be 1, 2, … in order. The timestamp file is snapshotted and its checksum is checked before analysis.
3. Set the inclusive first and last **normalization baseline frames**. Suggested onset follows the final baseline frame; always replace it with your actual application timing. The page warns when the fixed F₀ baseline includes onset or later samples. This fixed baseline creates the global trace `(F−F0)/F0`. The local prestimulus baseline in each protocol row is a separate reference for that response call.
4. Choose or load a protocol, inspect each row and any imported field timing overrides, and enter condition and independent biological-unit IDs when known. Check the timing/selection confirmation and click **Run analysis** in the browser, or **Analyze responses** in the native fallback.

Setup needs internet and disk space for R and downloaded packages. It is separate from Cellpose and stored under `.calcium-simpleca` in your home folder. A previous base-only R environment is not reused as protocol-ready: setup verifies the required R and jsonlite versions. Retry uses the package cache, and later analyses can run offline. Cancelled setup may leave partial installation directories, but it never marks them ready.

## Stimulus timing

Each stimulus has one **Stimulus added at** time. **Noise baseline to** and **Score response from** follow it by default and show "= stimulus time"; **Set separately** gives either its own value, and **Link to stimulus time** reverses that. **Use frame N** sets the stimulus time to the movie frame on screen, and linked fields move with it. A typed time that reads the same as a frame's time is stored as that frame's exact time.

Under the fields, one line states the frames each window will use, by the engine's own rule: the noise baseline takes frames from its start up to, but not including, its end; scoring takes frames from its start to its end inclusive. The same windows are drawn under the frame slider. **F₀ frames** set the zero for ΔF/F₀ and are separate from the noise baseline, which sets how large a rise must be. When two time fields disagree, the message sits beside the field to change and offers the repair.

## Response protocols

Each stimulus has its own response rule: the primary stimulus and, optionally, one named positive control. ATP, capsaicin, KCl, ionomycin or another control name is metadata; the software does not infer stimulation concentration or an optimal threshold from that name. This version supports one primary window and at most one control window. The browser keeps timing and a rule summary visible; expand **Adjust the response rule** for detailed threshold, persistence and noise settings. **Protocol name, files and overrides** and **Experiment details and diagnostics** also expand when needed. In the native fallback, scroll the protocol table horizontally to inspect every setting. Unedited suggested values retain their underlying precision.

| Setting | Meaning |
|---|---|
| Baseline start/end | Local reference interval `[start,end)`, ending at or before application onset |
| Applied onset | Recorded stimulus application time; the reference for latency |
| Response start/end | Closed scoring interval `[start,end]`, beginning at or after onset |
| Minimum rise | Increase above the local baseline. The browser's percentage field displays `10` for 10%; the protocol stores this as the fraction `0.10` |
| Baseline SD multiplier | Multiplies the cell's noise, expressed as a fraction of the local reference (see **Noise** below) |
| Persistence | Actual timestamp span in `seconds`, or consecutive sample count in `samples` |
| Maximum gap | Largest consecutive timestamp gap allowed; blank uses the engine's declared automatic gap rule |
| Minimum baseline samples | Required finite samples in the local reference interval |
| Noise floor | Optional minimum baseline noise. Browser percentage fields are converted to the protocol's fractional units; this setting does not exclude a cell |
| Expected interval | Optional declared sampling interval used by coverage/gap checks |

The threshold is `max(minimum rise, SD multiplier × noise)`, with the optional noise floor applied before multiplying. Samples equal to the threshold qualify (`>=`). Only the setting for the selected persistence mode is editable; unused values are retained for replay. Seconds mode uses elapsed span and its inherent minimum of two observations, without adding the sample-mode count as a second requirement. There is no smoothing. The local trace is `(y−b)/(1+b)`, where `y` is the fixed-reference fractional trace and `b` its local baseline mean. Therefore a local threshold `T` is drawn on the global trace at `b+(1+b)×T`.

**Provisional slow evoked: 2 samples** starts at a 0.10 fractional floor, 3 times the cell's noise and 2 qualifying samples in a row. It is a starting protocol to validate on the experiment's controls, not an established optimum for 50B11 or every indicator. The length can instead be given in seconds: last timestamp minus first timestamp is the duration, so a one-second span at 2 Hz needs at least three qualifying samples. Missing values and excessive gaps break the run. Version 2 requires complete window coverage and enough reference samples; unavailable coverage is unscorable with a reason, rather than treated as no response.

**Custom version 2** keeps all settings independently editable. Timing edits retain different stimulus/control threshold, duration and noise rules. **Legacy sample rule (starts at 2)** exports actual schema 1 and calls the historical responder engine, including its historical endpoint handling. It is not version 2 with a two-sample setting. Legacy mode uses an immediately preceding reference and its optional control remains KCl. Version-2-only fields are locked or labeled unused; editing the legacy scoring start also moves its reference end and onset. Historical settings TSV files still follow the same legacy path. New explanatory metrics/QC annotations do not change historical classifications or counts.

The browser's **Import protocol JSON / Export protocol**, or the native **Load protocol / Save protocol**, exchange declarative JSON with SimpleCa²⁺. Files record schema and engine version, name, validation status, the noise estimator, each window and optional field overrides. Files written by the plugin name the estimator at the top as well (`sd_estimator`); SimpleCa²⁺ reads the copies in the windows. Imported overrides are retained in the browser's expandable protocol settings; the native window also shows a recording match count and offers **View field timing overrides**. Overrides can be explicitly cleared. A matching override changes that recording's timing while preserving its rule. `validated_by_user` is a user declaration, not certification by the program. Unsupported schema, reference, engine or rule semantics are rejected, and so is an estimator recorded with another engine version than the earliest that defines it. The browser's import makes one allowance: when a file names an estimator at the top (`sd_estimator`) and its windows say `sample_sd` with engine version 2.0.0, the name at the top stands, because the page writes its own settings that way. **Load protocol** refuses such a file. SimpleCa²⁺ reads a protocol that uses the whole-trace estimator from its engine 2.1.0 on (pull request 14 of the app, not merged when this was written); the app with engine 2.0.0 refuses it as an unsupported engine version. Protocols on the baseline SD go both ways.

### Noise

A protocol names what a cell's noise is taken from (`SD_estimator`). Both estimates are fractions of the local reference, like the rise they are compared with.

| Estimator | Noise | Recorded with engine |
|---|---|---|
| `larger_of_baseline_sd_and_trace_noise` | The larger of the sample SD of the local baseline and the whole-trace noise. Every protocol the plugin makes uses it, including a version 2 rule made from a legacy one | 2.1.0 |
| `sample_sd` | The sample SD of the local baseline alone. The only estimator before engine 2.1.0, and the legacy rule's | 2.0.0 |

Whole-trace noise is 1.4826 × the median absolute deviation, about their median, of the differences between consecutive samples over the whole recording, divided by √2. It measures frame-to-frame noise, so a response moves it little. A difference counts when both samples are finite and no further apart than the maximum gap. A baseline of a few frames can come out quieter than the cell is; taking the larger of the two keeps such a cell from being called on two chance samples. A cell with fewer usable differences than the minimum baseline samples keeps its baseline SD and is flagged `WHOLE_TRACE_NOISE_UNAVAILABLE`.

An imported protocol and every saved analysis keep the estimator they were saved with. Settings remembered for a recording keep an estimator they name. Settings remembered before the plugin carried engine 2.1.0 name none, because there was only one: their timing and thresholds are filled in, noise is taken as a new protocol takes it, and the status line says "Noise is now the larger of baseline SD and whole-trace noise." The native window shows the noise in use under the protocol name. The browser page does not yet show or change it; when remembered settings that name the baseline SD are filled in, its status line ends "Noise: baseline SD, as saved." for as long as the settings use it, and the finish message of an analysis on the baseline SD ends the same way. A protocol records the earliest engine release that defines its rule, not the engine that wrote or scores it, so a `sample_sd` protocol keeps its hash and gives the same tables under engine 2.1.0 as under 2.0.0. Any other pairing of estimator and engine version is refused.

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

## Checks on the calls

After the engine has finished, the plugin checks every scored cell and writes `<recording>_call_checks_<analysis folder>.csv` beside `<recording>_percent_responding_<analysis folder>.csv`. A check is a flag for review. It never changes a call, a count or a percentage, and the analysis folder itself is not written to. The same table is made from the browser page and from the Fiji window. The page receives the checks with the results; it does not show them yet. Once the engine has finished, the analysis can no longer be cancelled: the checks run to the end.

### What Measure and save records around each cell

In the same pass over the frames as the cells and their backgrounds, the mean raw gray value of 25 fixed sets of pixels per cell is measured:

| Zone | Pixels |
|---|---|
| Near ring, 8 sectors | Outside the outline, more than 0.645 µm and up to 6.45 µm from it |
| Far ring, 8 sectors | Outside the outline, more than 6.45 µm and up to 14.19 µm from it |
| Core | Inside the outline, deeper than the median depth |
| Rim, 8 sectors | The other pixels inside |

Sectors are 45° wide and centred on right, lower right, bottom, lower left, left, upper left, top and upper right, seen from the mean x and y of the outline's pixels. Distances are Euclidean between pixel centres: outside, to the nearest pixel of the outline; inside, to the nearest pixel that is not in the outline, with the picture's border counted as outside. Both rings leave out every pixel within 1.29 µm of another measured outline and keep bright pixels. A zone with fewer than 12 pixels is not used, whatever the reason: a sector mostly covered by other outlines, one cut by the picture's border, or a small rim. The outlines' pixels are the same ImageJ masks the cell values are measured with.

The widths were set as 1, 10, 22 and 2 pixels on recordings at 0.645 µm per pixel. They are kept in micrometres and turned into pixels with the recording's pixel size, to a hundredth of a pixel. A recording without a pixel size in micrometres takes the pixel values as they are, and so does one whose pixels are not square, within 1 %; `what_was_done.txt` says which. The unit is read as micrometres by the same test as the 100 µm scale bar.

| File | Content |
|---|---|
| `<recording>_ring_and_core_gray_values.csv.gz` | One row per frame. `Frame`, then for every cell 25 columns `<cell>\|near_right` … `\|near_upper_right`, `\|far_right` … `\|far_upper_right`, `\|core`, `\|rim_right` … `\|rim_upper_right`. Two decimals; empty where the zone has fewer than 12 pixels |
| `<recording>_ring_and_core_pixels.csv` | One row per cell: `Cell`; `ROI_sha256`, the outline as the ROI registry names it; `Area_pixels`; the pixel count of each zone; `Near_ring_sectors_usable`; `Raw_table_sha256`, the hash of `<recording>_gray_values_raw.csv`; `Ring_and_core_values_sha256`, the hash of the values table above; `Pixel_size_um`; the three ring widths in pixels; and `Other_outline_margin_px`, the margin left around other outlines |
| `<recording>_field_movement.csv` | `Frame`, `Displacement_px`, `Shift_x_px`, `Shift_y_px`: the shift of the field at every frame as the movement check of Step 1 measures it, and how far that is from where the field sits for most of the recording. The displacement is the middle value of five frames, with a jolt of one or two frames kept. Written only when Step 1 could measure the movement, which needs a picture at least 196 pixels wide and high and 3 frames or more |

`what_was_done.txt` gains one line, "Around each cell: …". No other file of Steps 1 to 3 changes. When one of the three files cannot be written, none of them is left, every other file is saved as it would be without them, and the line says "Around each cell: the rings and the core were not saved: <file> could not be written. The cell measurements are unchanged." The exception goes to Fiji's Log.

Step 4 uses these files only when they name the cells, the outlines and the raw table of the measurement being analysed, and when the values table is the one the pixels table names. When an analysis runs and they are missing or belong to another measurement, they are measured then with the same code and saved beside the tables, provided the recording is open in Fiji and gives the saved raw gray values on its first, middle and last frame. In practice this is when the rings could not be saved at Measure, or the files were taken away or replaced since: Step 4 runs on the measurement just made. A recording measured with an earlier version has no ring files and must be measured again to be checked. Otherwise the passing-object check is not made and every row says `not checked` with the reason. The field movement is not measured again at this point for files that belong to the measurement.

### The rules

Only the event of the largest rise is looked at. The F0 frames are the analysis's own normalization baseline frames.

- **Traces**, in gray levels: the cell is raw minus its background minus its mean over the F0 frames. A zone is its mean minus its own F0 mean minus the change of the cell's background. Both are averaged over 3 frames. A cell or zone with a value that is not a number is not used.
- **Event**: the frames around the peak of the call where the averaged cell trace is at or above half its value at that peak. A dip of up to 3 frames does not end it. The cell's rise is its value at the peak.
- **Noise of a zone**: the larger of the SD of its averaged trace over the F0 frames, and 1.4826 × the median absolute deviation of the frame-to-frame differences of its unaveraged trace over the whole recording, divided by √2 and by √3.
- **Passes through**: a near- or far-ring sector rises in the 10 s before the event and a sector at least 90° away rises in the 10 s after it. Each rise is at least the cell's own rise and at least 5 × that sector's noise, is gone while the cell is up, with its mean over the event at most half its peak, and is at a frame where the cell is still below half its rise. Text: `possible passing object, frames 201 to 210`, the frames of the event.
- **On the edge**: over the event the core rises less than 0.25 × the whole outline, and the near-ring sector outside the rim sector that rises most rises at least as much as that rim sector and at least 5 × its own noise. No neighbouring sector stands in for a sector that is not usable, and without a usable core the rule is not applied. Text: `rise on the lower right edge, not in the cell, frames 92 to 149`: the frames of the event, widened on each side for as long as that ring sector stays above half its mean over the event. The text never says "object": in the recordings the rule was set on, most such rises were a neighbour without an outline.
- **Field moved**: when one of the two rules applies, the field's shifts from 10 s before the event to 10 s after it are looked at. The field counts as moved when at one of those frames it is 8 pixels or more from where it was over the F0 frames, with the displacement as `field_movement.csv` gives it, or, when no frame is, when it moved 8 pixels or more between two of those frames, each taken as the middle value of five frames. When the field moved within the F0 frames themselves by 8 pixels or more, it is counted from where it sits for most of the recording instead. The row then reads `field moved, frames 28 to 35` instead of the rule's text: the longest stretch of displaced frames, or the two closest frames it moved that far between. Its detail says how far, every displaced stretch, and what would otherwise have been written.
- **Short separate crossings**: for a responder, the stretches of samples at or above its required rise inside the response window are counted from the engine's own columns, `Threshold_increase`, `Baseline_mean`, `Local_reference_relative`, `Start_s`, `Effective_end_s` and `Max_gap_s`, and the saved normalized trace. Three or more stretches, none longer than 2 samples, give `short separate crossings`.

A stretch of one frame reads `frame 28`. The 10 s span was set as 10 frames at 1 frame per second and is turned into frames with the median frame interval of the analysis. The 3-frame average and the 3-frame dip are counted in frames.

### The checks table

One row per measured cell.

| Column | Content |
|---|---|
| `Cell` | The cell's ID |
| `Call` | The engine's status, unchanged |
| `Flag` | The flag's text; empty when there is none. Two flags are joined with "; " |
| `Rule` | `passes through`, `on the edge`, `field moved` or `short separate crossings`, in the order of `Flag` |
| `Event_first_frame`, `Event_last_frame` | The event of the largest rise, counted from 1 |
| `Detail` | Which sides and frames, or the crossings, for each flag; empty without a flag |
| `Near_ring_sectors_usable` | Of 8 |
| `Area_pixels` | Pixels in the outline |
| `Ring_measured` | `yes` or `no` |
| `Field_movement_measured` | `yes` or `no`. Without the movement a moving field is not told from an object, and the analysis says so |
| `Passing_object_check` | `checked` or `not checked` |
| `Not_checked_because` | Empty when checked. For one cell: `no call`, `no peak time`, `missing gray values`, `no rise at the peak` or `no ring around the cell`, the last when no ring sector has 12 pixels. For the whole analysis: why the ring and core could not be used, for example `ring and core were not measured and the recording is not open` |
| `Crossings`, `Longest_crossing_samples` | The separate crossings of the required rise in the response window and the longest in samples, for every scored cell |

When the calls cannot be checked at all, no table is written and the message says why, for example "The calls were not checked: the calls table has no Cell or Status column." The details go to the analysis log.

### What the check cannot see

- It looks only at the event of the largest rise. In 53 of 330 responders of the nine test recordings the samples that made the call lay outside that event.
- By eye, 14 capsaicin and KCl responders in those recordings were objects. The check flagged 9: 8 of 9 in outlines of 900 pixels or more, 1 of 5 in outlines under 700 pixels.
- It misses most objects that cross in 2 frames, and objects that come to rest inside the outline and stay.
- Each near-ring sector lost to a neighbouring outline lowers the chance of a flag. `Near_ring_sectors_usable` and `Area_pixels` are in the table for that reason.
- On the edge cannot tell an object resting on the edge from a neighbour without an outline that responds beside it.
- A flag is not proof of an object. Of 296 ATP responders in three recordings 2 were flagged and neither was an object. The thresholds were set and tested on the same nine recordings.

A sentence for a Methods section: "Responder calls were checked for bright objects passing through an outline, using the fluorescence in sectors of two rings around each outline and in its core and rim; flagged cells were reviewed in the movie. The check examines only the largest rise of each cell, does not detect most objects that cross within two frames, and does not detect objects that remain inside the outline."

## Review marks

A cell the engine scored for the primary stimulus can be marked the other way after review: `Responded` to `Did not respond` or the reverse, with a reason: `Noise`, `Movement`, `Passing object`, `Bad outline` or `Other`. `Other` needs a note. A note is one line of plain text, at most 120 characters, without hidden characters, and cannot start with =, +, - or @, which a spreadsheet would read as a formula; spaces at either end are taken off. A mark can be removed. Excluded and unscorable cells cannot be marked. The control window is never reviewed.

| File, beside the percent table | What it holds |
|---|---|
| `<recording>_reviewed_calls_<analysis folder>.csv` | One row per mark: `Cell`, `Automatic_call`, `Reviewed_call`, `Reason`, `Note`; `Time`, ISO 8601 to the second with the offset; `Analysis`, the analysis folder's name; `Responder_calls_sha256`, the hash of that folder's `responder_calls.csv`; `ROI_sha256`, the outline as its ROI registry names it; `Carried_from`, empty for a mark made in this analysis, else the analysis the mark was made in. Present only while there is at least one mark. |
| `<recording>_percent_responding_<analysis folder>.csv` | With marks: the same rows and columns, then `responded_after_review`, `percent_responding_after_review`, `responded_without_them_after_review`, `percent_without_them_after_review`, `cells_reclassified`. Without marks: as before, byte for byte. |

- **One place applies the marks**: `CalciumAnalysisSettings.tally`. The percent table, the finish message, the Fiji window's headline and groups, and the `tally` and `review` the Step 4 page receives all come from it, `review.after` with the call of every cell after review. A mark moves a cell between the two calls only, so the registered and evaluable counts of every group stay as scored; `cells_reclassified` counts the cells of the row whose call a mark changed.
- **Saved together**: the tally is worked out from the marks first, then the percent table and the review file are each written beside their place and moved in. When the percent table cannot be replaced, for example while it is open in another program, nothing changes and the mark is refused: "Close <file> if it is open in another program, then try again." When the review file cannot follow, the percent table is put back.
- **Kept in line**: when the results are read, a percent table whose after-review columns no longer match the marks, because the review file was taken away or no longer matches, is written again. A table that is right, or has no marks and no such columns, is not touched.
- **The status line follows the marks**: after every mark the finish message is made again from the tally, and the part about carried marks goes once they have all been removed.
- **Never rewritten**: `responder_calls.csv`, the manifest and every other file in the analysis folder. The review file is the record of the review.
- **Checked**: a review file whose `Analysis` or `Responder_calls_sha256` does not match its analysis, or whose rows do not fit its calls, is not applied, and the answer says so. A new mark is refused until it is moved away. A byte order mark at the start and blank lines at the end are ignored, and an empty file holds no marks.
- **Carried over**: when a new analysis of the same recording finishes, the most recent earlier analysis made by this version is found: one with a checks file or a review file beside the tables. Its marks are written for the new analysis for every cell whose automatic call and `ROI_sha256` are unchanged. The others are dropped and listed with the reason: `now scored as did not respond`, `now scored as responded`, `now excluded`, `now unscorable`, `no longer measured`, `outline changed`. A cell no longer measured is named by its number in the earlier analysis, "Cell 12 of the earlier analysis". A carried mark keeps its reason, note and time. When that analysis has no marks, nothing is carried, so marks removed or dropped there do not come back from an older one. When its review file does not match it, nothing is carried and the message says so.
- **Methods sentence**, built from the marks that changed a call of all cells, for example: "3 of 143 cells were reclassified after manual review: 2 from responder to non-responder (1 noise, 1 passing object) and 1 from non-responder to responder (other)."
- **Plots saved from the Fiji window** name, in `plot_provenance.json`, `calls: after review`, `review_file` and `review_sha256` when marks apply. The window's groups then read "Responded after review" and "Did not respond after review". The window reads the marks when it loads results.

The Step 4 page sets and removes marks with `POST api/review`: `cell`, `call`, `reason` and `note`, or `cell` and `remove: true`, and `folder` to name the analysis shown. A body that is not a JSON object is refused with "The request is not valid." The page's control is not built yet.

## Files and source provenance

Each successful run gets a new `SimpleCa-date-time-id` folder. It contains the input and optional timestamp/ROI/diagnostic snapshots, settings, `protocol.json`, `applied_protocol.json`, normalized traces, mean trace, metrics, quality, ROI registry, cohort flow, calls/counts, window metrics, replicate/group counts, metadata, manifest and log. Optional diagnostic outputs are separate. Previous runs are not overwritten.

The bundled `utils.R`, `response_protocols.R` and `responders.R` are unchanged from [SimpleCa²⁺ commit 67efdd24](https://github.com/samgillman/simplecaapp/tree/67efdd241c3f9b25df3755ee359c153e45aae847), response engine 2.1.0. The adapter uses their normalization, metrics, protocol import/validation, response calls and summary functions. Source SHA-256 hashes, plugin/R/jsonlite versions, actual input checksum, applied protocol and units are saved with results. `settings.tsv` and `analysis_manifest.csv` name the noise estimator (`sd_estimator`; `sample_sd` for the legacy rule); the manifest also gives the engine's own sentence for it (`noise_rule`), the engine that scored the run (`engine_version`) and the release the protocol is recorded with (`protocol_engine_version`). With the whole-trace estimator, `responder_calls.csv`, `sensitivity_calls.csv` and `negative_control_calls.csv` have three more columns, placed last so that every other column keeps its position: `Whole_trace_noise`, `Whole_trace_noise_samples` and `Noise_source` (`Baseline SD`, `Whole-trace noise` or `Noise floor`; blank for a cell that was not scored). Read `Noise_source` together with `Threshold_source`, which names the term that set the threshold: `Floor` is the minimum rise, and `Baseline noise` is the SD multiplier times the noise, whichever estimate the noise is. A threshold was set by the whole-trace noise when `Threshold_source` is `Baseline noise` and `Noise_source` is `Whole-trace noise`; the native window's line under a cell's plot names the estimate. A `sample_sd` run has the same columns as before. The shared canonical protocol MD5 is an identity for comparison, not a security signature. Managed package records remain in the environment's `conda-meta` directory.

Synthetic tests cover independent hand-calculated response outcomes as well as parity with the shared engine. CI is configured for managed runtime checks on Windows, macOS and Linux and native GUI checks on Windows/Linux. Browser acceptance is a separate Playwright gate described in the release guide; server and JavaScript syntax checks alone are not browser validation. These checks do not establish biological sensitivity/specificity or replace acquisition-specific control validation.

## Step 4 editing and quality review

Invalid frame rates, frame ranges, reversed windows and negative thresholds show inline errors before analysis; submission runs the authoritative validation again. A timestamp CSV takes priority over FPS. Elapsed-time and legacy rule types retain separate drafts when switching. **Discard pending changes** restores the complete applied draft, including uploaded file contents, or the starting draft before the first analysis. Completed results, movie timing and exported plots keep their applied settings until a new analysis succeeds.

A newly enabled positive control has blank timing fields: enter its actual baseline, addition and scoring times. The software cannot infer control timing from the end of the stimulus window. Recording-coverage and baseline-overlap warnings do not silently alter a protocol.

The quality summary and per-cell feedback show saved flags, including `DURATION_BELOW_SAMPLE_INTERVAL` and `reference_below_3_baseline_sd_review`. A separate review indicator identifies saved traces below −1 and legacy rise times whose 10%-peak level is below baseline SD. These display diagnostics do not recalculate response calls or add exclusions. **Needs quality review**, **Excluded** and **Unscorable** are separate filters. A low reference can make fractional values unstable; a value below −1 corresponds to negative background-corrected fluorescence when F₀ is positive. Inspect background and baseline stability before interpreting these traces.

**Response-window metrics** provide the applied per-stimulus peak, latency and signed AUC. **Legacy whole-recording metrics** retain their established definitions: `AUC` integrates the whole recording and can be negative under downward drift; `Time_to_Peak` is the recording timestamp; `Rise_Time` uses 10–90% crossings after the display baseline, which can be sensitive to baseline noise. The interface labels these definitions rather than changing old values. Use `Peak_latency_s` for time from stimulus onset and `Window_AUC_fraction_s` for the declared response window.

The calls table starts with core response/quality columns and a pinned cell-name column. **Show all saved columns** exposes the complete table; saved CSVs always retain every column. Empty tables explain that no rows were saved. Short windows can scroll vertically while plots and captions keep a readable height.

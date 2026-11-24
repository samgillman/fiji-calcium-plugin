# Changelog

## Unreleased

Step 4 page: how the time course and its legend look, on the page and in a saved figure. Display only: no measurement, response call or saved table changes.

- **A cleaner plot.** Left and bottom axes only, ticks pointing out, no box and no grid lines; one thin dotted line at ΔF/F₀ = 0. Tick numbers use a true minus sign. The numbers on the Y axis are round, and the axis ends on a tick or half a step beyond the last one, so little of the plot is empty; the time axis ends where the data end. The editor still offers a box and grid lines.
- **No extra words on a plot.** Besides its axis titles and tick numbers, a plot carries one word, "stimulus", above the stimulus line. Everything else that is drawn is named once, in the legend. This holds on the page and in a saved figure, for the time course and the heatmap.
- **One look for the page and the saved figure.** Both are drawn by the same code from one table of defaults. A saved figure has black text and axes; the page keeps its dark gray.
- **A saved plot is 7 × 4 in and has no title** until you type one in Edit figure. The heatmap figure has no title of its own either. The recording's name is still in the file name and in the file kept beside the figure.
- **The Y axis is no longer stretched by a few cells.** The range holds every mean ± SEM that is drawn. With All traces on it also takes in the bulk of the single cells: from the 2nd percentile of their lows to the 95th percentile of their highs, each cell by its own 1st to 99th percentile, over cells that are neither excluded nor below zero. It reaches no further beyond the means than one and a half times their own height, so the means fill at least a quarter of the axis. A trace that runs past the range ends at the plot's edge. The selected cell is always shown whole.
- **All traces:** cells that did not respond are thin, neutral gray and faint, fainter the more cells there are. Responders are drawn after them in the responder colour, a little wider and stronger. The mean is drawn last, on its white casing.
- **The mean** has its SEM band at 25 %. A group of one cell has no band and no "± SEM".
- **The stimulus label and the response window sit above the plot,** not on the traces: the stimulus is a thin dashed line with its label above it, the response window a thin bar clear of the plot's top edge.
- **Legend:** one entry for each thing drawn, and no entry for anything that is not drawn. Each group once, with the mean line as its swatch: "all cells ± SEM (n = 113)", or "responded (n = 1)" for one cell. With All traces on, one entry "single cells". "baseline", "response window" and, for a cell that responded, "response" are small entries. A swatch shows what is drawn: the line at its width, dashed when the line is dashed, on its band when a band is drawn.
- **One cell: one threshold line** (#31). "One cell + average" draws the required rise, with one legend entry. The required rise is always either the minimum rise or the noise threshold, so the old plot drew one of those lines hidden under it and still listed it. The new box **Thresholds**, off until you tick it, adds the threshold that did not set the required rise ("minimum rise" or "noise threshold", as the call table's `Threshold_source` says) and the baseline mean, each with its own entry. Two levels with the same value are one line. The saved figure follows the box.
- **The verdict above the plot says what set the required rise:** "needed 30.3%, set by its noise" or "needed 10.0%, the minimum rise".
- **The response band is never too thin to see:** at least 3 pixels wide (2.25 pt in a figure), with a small tick at the top of the plot where the run starts, so a run of two or three frames shows.
- **A selected cell without a ΔF/F₀ trace** is said in the cell panel above the plot, no longer in the legend, and a saved figure does not carry the note.
- **Legend placement:** above the plot, starting on the Y axis, with even gaps. When it is too wide for one row it wraps by whole entries into rows of even length (two and two, three and two), so no entry is left alone on the last row. Edit figure › Legend has a new choice, **Right of the plot**, with the entries stacked.
- **"One cell + average" shows the average of the cells under Show,** as "Average" does, and the legend names the group and its number of cells. It used to show the average of all cells whatever Show said.
- **"Leave out … below zero" now reaches the plot:** those cells are not drawn and not averaged while the box is ticked, and the legend's n says so. The file kept beside a saved figure lists them (`cells_left_out`).
- Edit figure: a colour box for the single cells; "Did not respond" is offered with Responders vs non-responders, where it colours that mean. With one cell: **Response** (the band), **Required rise** (the line; this box was called Thresholds) and **Thresholds** (the page's box), each taking its entry out of the legend with its mark.
- The file kept beside a saved plot says whether Thresholds was ticked (`thresholds_drawn`).
- The thresholds of "One cell + average" are drawn only for a selected cell that has a call.

Step 4 page: the heatmap, on the page and in a saved figure. Display only.

- **Blue–white–red is the default,** with white at ΔF/F₀ = 0. Each side has its own scale, so a small negative side does not flatten the red side. The other palettes are still in the list and still start at zero.
- **The colour range is no longer set by a few cells.** It comes from the cells that are neither excluded nor below zero: the 1st to the 99th percentile of their values, rounded outward to tidy numbers. The other choice under Range, **Smallest to largest**, uses the same cells. On the recording that prompted this the range was 0 to 3; it is now −0.6 to 0.2.
- **The colour bar says when values lie beyond it:** "≤" at the lower end and "≥" at the upper end, each only when a value that is drawn is beyond that end, with limits you typed too. Both ends are always written in full. A palette that starts at zero shows "≤ 0" in the saved figure as well.
- **"Scale to" has a partner, "Scale from":** Edit figure takes a minimum and a maximum for the colours.
- **Sort:** "Response call, strongest on top" ranks by the largest rise in the response window. Cells below zero come last within their group (in "Time to peak" and "Peak amplitude", after the other cells), and a dark mark beside their rows shows them. New order **Onset time**: responders in the order they began to respond; the rest as by response call.
- **Rows:** the axis no longer prints row numbers that read like cell numbers; it is titled "Cells". The strip of colours beside the rows shows the groups, and nothing else is drawn or written there. **Cell names** prints the real names beside the rows, every few rows when the rows are thin, and the heatmap stays in view with its time axis. New box **Hide excluded cells.** Both boxes are in Edit figure too, so a saved figure carries them.
- **Key:** drawn in the heatmap, on the page and in the figure: the calls that are present, each with its number of cells ("did not respond (n = 115)"), "below zero", "baseline" and "response window". The stimulus is labelled above the picture, at the time saved as the stimulus.
- **Edit figure:** X min, X max and X ticks every, for the time axis.
- **Layout:** the picture fills the width; the colour bar is narrow, at the right, as tall as the picture. Margins follow the text. A saved heatmap is 7 × 5 in.
- The file kept beside a saved heatmap lists the cell of each row from the top and the colour range used.
- The Fiji window's heatmap takes its colour limit from cells that are not below zero and puts cells below zero last in their group, as the page does.

Step 4 page: zooming and moving the picture of the cells, plain outlines and numbered cells. Display only.

- **Zoom and move controls on the picture.** On its top left corner: **+** and **−**, the zoom level between them ("2.5×"), **Whole picture**, and **Zoom to region**. They are small buttons with tooltips. The whole picture is fitted beside them, never under them: they lie along the top, or down the left side where that leaves the picture larger.
- **Zoom to region:** press the button and drag a rectangle on the picture; the view zooms to it and the button turns itself off. Esc cancels. Shift + drag does the same without the button. It works in every view and while the movie plays.
- **Drag moves the picture** whenever it is zoomed in, under a hand that closes while it drags, with a cell selected too. A drag never selects or unselects a cell. The picture cannot be dragged out of view: no edge of it passes the middle, so a cell at the edge can still be brought to the middle.
- **Keys:** with the pointer over the picture, **+** and **−** zoom about the pointer and **0** shows the whole picture.
- **A pinch on a trackpad follows the fingers.** Chrome and Edge report a pinch as small scroll steps; these were taken at a scroll wheel's pace, about a tenth of the pinch. Ctrl or Cmd + scroll zooms as before, and plain scrolling still scrolls the page. Safari's pinch is taken too.
- **A small map** on the picture's bottom right corner while it is zoomed in: the whole field, with a frame around the part on screen. A click or a drag in it moves the view.
- **A hint** beside the controls the first three times the pointer rests on the picture: "Pinch or Ctrl + scroll to zoom. Drag to move." (Cmd on a Mac.) It goes when the pointer leaves, and after the third time it is not shown again.
- Zooming in on a paused movie fetches the frame at full size at once. Where the frame was shown reduced, it used to stay coarse until the frame changed.
- At 1× a drag no longer moves the picture: there is nothing to bring into view. Zoom steps are 1, 1.5, 2, 3, 4, 6, 8, 12 and 16×.
- **Outlines: By call or Plain** (#30), beside the outline thickness. By call is as before. Plain draws every outline white on its dark halo, all of one width, and keeps the dashes of excluded and unscorable cells. The key then reads "your cells", and names the dashed outline while one is drawn.
- **Number cells** (#30), a box beside it, off until ticked: each outline that is drawn gets its number, the number of its name on the page ("Cell 12": 12), white on a dark halo in the middle of the cell. Numbers that would lie on one another are moved apart, and none runs off the picture. The selected cell keeps its name tag.
- **Saved pictures and frames follow the screen:** Save to plots folder, Download a copy and Save frames… draw the outlines plain and the cells numbered when the page does. The file kept beside them gives the outline colours and says whether the cells were numbered (`outlines_shown`, `outline_colours`, `cells_numbered`). With By call and no numbers, saved pictures and frames are pixel for pixel what they were.
- Both choices are kept for the next recording, as the outline thickness is. Show still chooses which cells are drawn; without All outlines there are no numbers.
- The Cell box and the numbers on the picture count cells the same way.

Step 4 page: saving the area on screen, and a baseline and peak pair. Display only.

- **Area: Whole picture or Area on screen,** in the picture's Export menu and in Save frames…. Whole picture is as before, however the picture is zoomed, and is what the page opens with: saved pictures and frames are then pixel for pixel what they were. Area on screen saves the part of the recording that is in the picture panel. While the picture is zoomed in, the Area row stands out in both places, since it then decides what is saved. At 1× the two are the same picture.
- **A picture of an area** is as wide as a picture of the whole field (Options sets that width) and has the area's shape. It is drawn from the recording's own pixels of that area, not from the screen. Zoomed in beyond 3×, where the page shows the pixels as squares, the saved picture does too; up to 3× it is smoothed as the page is.
- **Outlines, numbers, labels and the key are not stretched with the area.** They are drawn with the pens and at the text sizes of every saved picture. Outlines end cleanly at the edge of the area. A cell is numbered when its middle lies in the area, and every number stays inside the picture. The selected cell is not marked in a saved picture, as before.
- **A scale bar** on a picture or frame of an area: a round length, at most a fifth of the area's width (20 µm, 50 µm, 100 µm), written over the bar. A picture has it at the right of its key, under the picture, where it covers nothing; a frame of a set, which has no key, on its bottom right corner. It is drawn when the recording says how large its pixels are in micrometres. A picture of the whole field has none, as before.
- **The file name says "area",** and so does the folder of a set of frames. The file kept beside the picture gives the area in pixels of the recording and in µm (`area_pixels`, `area_um`: left, top, width, height), the zoom (`zoom`) and the length of the scale bar (`scale_bar_um`). For the whole field it says `"area": "whole picture"`.
- **Save frames… uses one area for every frame:** the area on screen when saving begins, also when the picture is moved or zoomed while the frames are saved.
- **Save frames… › Baseline and peak** (new, beside Every Nth frame, which is still what the page opens with): one picture with the baseline frame on the left and the peak frame on the right, for showing a response. Over each panel "baseline" or "peak" with its frame and time; under the pair the key, once, with one scale bar at its right. Each panel is drawn as a saved frame is: the view on screen, the same outlines and numbers, the same area.
- **Its two frames are filled in and can be changed.** Baseline: the last of the F₀ frames. Peak: the selected cell's peak in the response window; with no cell selected, the frame of the response window where the cells that responded are brightest on average (the cells that got a call, where none responded). Before an analysis both are the frame on screen. A frame you type is kept while that suggestion stays the same. A box that does not hold a frame of the recording keeps what was typed, is marked, and says beside it which frames there are ("1 to 300"); nothing is saved until it does.
- **Both panels have the same brightness.** The plugin shows every frame of a view between the same gray levels, so the two panels can be compared. The file kept beside the pair gives those levels (`brightness_range`) and names both frames with their times (`baseline`, `peak`) in place of the frame on screen.
- **Save picture** puts the pair into the plots folder; **Download a copy** hands it to the browser. Format and resolution are those of Options. Each panel is as large as a frame saved alone, so the pair is twice as wide; the panel says the pair's size in pixels before you save. The file name says baseline and peak with their frame numbers.
- The page's server tells the page the size of a pixel (`micronsPerPixel` in the session) and the gray levels a view is shown between (`api/range`).

Step 4 page: Edit figure, the editor of a saved plot or heatmap. Display only.

- **A number is taken whole.** A number box takes what you typed when you press Enter or leave the box, not letter by letter: typing 40 into Width used to make the figure 4 in wide, and to keep that width for the next recording.
- **A number the figure cannot use stays in its box,** marked, and the line under the box says what the box takes: "1 to 20 in", "4 to 48 pt", "Above 0", "Not a number". The figure keeps the number it had, and nothing is kept for the next recording. It used to be put back to the old number without a word. Only plain numbers count: "1e1" and "0x10" were read as 10 and 16; a decimal comma is still a decimal point.
- **Axis limits must leave a range.** A minimum at or above the maximum is not taken: "Minimum must be below maximum", "Maximum must be above minimum". One limit alone is held against the automatic other end, and the line names it: "Must be below 298.969" on a recording that ends there. This holds for the plot's X and Y limits, the heatmap's X limits and its Scale from and Scale to. Before, such X limits were ignored without a word and such Y limits drew an empty plot. A number that waits in its box is taken as soon as the other end allows it. A limit the figure holds and can no longer use, because the other end has moved, is marked and not drawn.
- **Tick spacing:** "every 0" is not taken ("Above 0"), nor a spacing that would give more than 60 ticks ("At least 5"). That spacing used to stay in its box while the automatic ticks were drawn.
- **Text that does not fit.** A text size, title size, width or height that would leave the plot no room (less than four heights of the text, either way) is not taken: "Too large for this figure size", "Too small for this text size". Text size 30 on a 7 × 3 in figure used to fill the figure with the legend. Text that leaves the plot under 40 % of the figure is taken, and the line under Text size says "Too large for this figure size" (under Title size, when the title is what takes the room). A figure whose plot has no room for another reason, such as a longer legend than its kept settings were made for, is shown with the reason beside its size in pixels and cannot be saved, from the editor or from the Export menu.
- **Small figures.** A title too long for the figure is shortened, with an ellipsis for what is left out. So are an axis title, the stimulus label and the colour bar's title, which takes at most a quarter of the width. A legend entry too long for its row goes on two lines instead, broken between words and never inside "± SEM" or "(n = 115)", so nothing of it is lost: it used to read "did not respond ± S…". A legend right of the plot takes at most a third of the width and stays clear of the figure's edge, which it used to touch. Legend rows of different lengths no longer leave the last entry alone on its row ("baseline" and "response window" share the last row under two long entries). A low plot has fewer numbers on its Y axis, so they never run into one another.
- **Single cells:** with All traces on, Lines has Trace line, Trace colour and Trace opacity (the cells that did not respond) and Responder line and Responder opacity (the responders drawn over them). Opacity is typed in percent, 1 to 100; an empty box is auto, fainter the more cells there are.
- **Reset to defaults** (the button was called Reset) puts every box the editor shows back to its built-in value: the figure's own, and also Resolution, Format and TIFF compression, the heatmap's Colours, Range, Sort, Cell names and Hide excluded cells, and Thresholds. It used to leave those as they were. Nothing of the figure stays kept.
- **"Settings kept from an earlier session are in use."** The editor says so in one line, above its boxes, for as long as one of them is.
- **File, for every saved picture:** Resolution, Format and TIFF compression stand under a heading of their own that says so. They are shared with pictures and frames, as before.
- **The editor opens on the format in use.** With TIFF kept from an earlier session the Format box read PNG while a TIFF was saved. The lists of Options show what was kept from the moment the page opens.
- **Fonts:** Arial, Helvetica, Times New Roman, and As on the page. A saved figure is written in Arial until you choose another; it used to be written in whatever font this computer shows the page in. Courier is no longer offered, and a kept Courier becomes Arial.
- A colour or a title is shown in the preview while it is picked or typed, and kept when it is chosen or the box is left.
- **Save frames…** goes by the same rule: Every, Baseline frame and Peak frame take a whole number of frames the recording has, as typed. "2.5" in Every used to be saved as every 3rd frame, and "1e1" was read as 10 in all three. A number they cannot use stays in its box, marked, with "1 to 300" beside it, and nothing is saved. Every is a plain box now, as the other two.
- The plot's own **Y axis** boxes on the page go by the same rule as the editor's limits: a limit that leaves no range, or text that is no number, is not used, and its box is marked. Y min 5 with Y max 1 used to draw an empty plot from 5 to 5.5.

Step 4 page: the tables and the cell panel. Display only: no measurement, response call, saved column or saved value changes.

- **Columns in words, with their units.** "Required rise", "Largest rise", "Onset, s", "Peak, ΔF/F₀", "F₀, gray levels" in place of `Threshold_increase`, `Peak_increase`, `Onset_s`, `Peak_dFF0`, `F0`. Pointing at a header says what the column measures and its saved name; pointing at a number gives the saved value. A column the page does not know keeps its saved name, with spaces.
- **Numbers rounded by their unit:** time 2 decimals, ΔF/F₀ 3, percentages 1 with "%", gray levels 2. Rises saved as fractions read as percentages ("30.3%", as in the verdict). The true minus sign, and no "−0.000" or "0.9999". One dash for a missing value. TRUE and FALSE read yes and no; flags read as on the warnings button ("Duration under a frame"). "Required rise set by" reads "Baseline noise" or "Minimum rise" (saved: `Floor`).
- **Percent responding has one decimal everywhere:** the headline, the quiet line, the tally and the counts table. The Fiji window's headline too ("66.7 %"; it said "67 %").
- **Columns that are the same by definition are shown once:** `Run_start_s` with `Onset_s`, `Response_Amplitude` with `Peak_dFF0`, `Reference_value` with `F0`, `Baseline_SD_input` with `Baseline_SD_raw`, `Window_peak_local_fraction` with `Peak_increase`. Show all saved columns shows both, and the tooltip says which is the same as which.
- **Each table starts with its core columns** (calls, window metrics, whole-recording metrics, cell quality, counts); **Show all saved columns** shows every one. The whole-recording metrics and the cell quality have the cell's **Call** beside the cell, and the saved `Status` there is headed **Normalization**.
- **The response timing of a cell that did not respond is a dash,** its saved value in the tooltip ("Saved: 3 · Did not respond"). This holds for the onset, latencies, run, the time of the largest rise, and the whole-recording time to peak, rise time and widths.
- **The time of the largest rise and the onset are told apart:** "Largest rise at" (the highest frame of the window) and "Onset" (the first frame of the held rise); the tooltip says the largest rise can come before the onset when that frame alone was not held. The window's Baseline SD says it is not the SD of the F₀ frames in the metrics.
- **The verdict says what a response rests on:** "Cell 69: Responded · rise held for 2 frames from 121.0 s, 102.0 s after the stimulus · largest rise 40.6% at 92.0 s (needed 30.3%, set by its noise)". A cell that did not respond keeps its reason.
- **The first column stays in place** when a table scrolls sideways, with the stripe and highlight of its row; on Processed data, Frame and Time both do. A click on a header sorts by it, with blanks last both ways.
- **Responder counts** (one row per window) are shown turned, a row for each column and a column for each window ("Capsaicin", "Ionomycin (control)"), so they need no sideways scroll.
- **Export on the calls and counts tabs,** as on the other tabs. The menus say what they do: "Save CSV to plots folder" and "Download CSV". A save from the metrics tab now has the file of how the page was set beside it, as the others do.
- **Summary of the metrics:** the same 23 rows for every recording, each once under its own name, with its unit and where it was measured; mean ± SEM with the decimals of its unit and n in a column of its own. All cells leaves out excluded cells. Onset, latency, time to peak, rise and width are for cells that responded. The window rows come from the window metrics (or the calls when there are none). Chosen blind and Added after the response are columns only when both have cells. The saved summary keeps its columns (Metric, Unit, Measured, Saved_as, Group, Mean, SD, SEM, N_cells) and has every row; a metric not saved has N_cells 0, and rises are in % as shown.
- **"Legacy" is no longer shown:** the table is "Whole-recording metrics", on the page and in the Fiji window.
- **The page and its windows are named for the plugin:** "Ca²⁺ Imaging (Wang Lab)" and "Step 4" in the page header and the browser tab (with a recording open, the tab reads "116 cells · Step 4 · …"), "Step 4 · Response protocols" for the Fiji protocol window, "after measuring" under Step 4 in the main window. A protocol file that is not valid JSON says "Choose a protocol file." Saved file and folder names are unchanged.

## 2026.11.22.2

The plugin is renamed **Ca²⁺ Imaging (Wang Lab)**. No measurement, response call or saved table changes.

- **Fiji's menu:** Plugins › Ca²⁺ Imaging (Wang Lab), one entry at the top of the Plugins menu, instead of Plugins › Calcium Imaging › Calcium Cells Step By Step. The windows, the installer and the Step 4 page carry the new name.
- **The file** is `scripts/Plugins/Ca²⁺_Imaging_(Wang_Lab).groovy` in Fiji's folder. Where a system cannot store ² or ⁺ in a file name it is `Ca2+_Imaging_(Wang_Lab).groovy`, and the menu reads Ca2+ Imaging (Wang Lab).
- **Updating from an earlier version:** its Install update button works as before and puts this version under Calcium Imaging. When this version first starts it moves itself to the new place and says so on the first screen; restart Fiji once more to start it from there. The old file goes to `calcium-plugin-backups`. An older copy left next to a newer one is put there too, so the menu shows the plugin once. A folder that is a link to elsewhere is never touched.
- **The ZIP** keeps plain-letter names (some Windows unzip tools garble others): it holds `Ca2+_Imaging_(Wang_Lab).groovy`, which the plugin renames when it first starts. Download names are unchanged, so earlier versions keep finding updates.
- Unchanged: the results folder (Documents › Calcium imaging results), settings, and `cellpose_python.txt`, which is still read from the Calcium Imaging folder and now also from next to the plugin.

## 2026.11.22.1

Step 4 page: a cell can be unselected with the mouse. No measurement, response call or saved table changes.

- **A click on empty space in the picture shows no cell**, as Esc already did. A click on a cell still opens it, and a click on the open cell still brings a zoomed view back to it. Dragging the picture to move it is not a click and keeps the cell.
- **A double-click on empty space shows the whole picture and keeps the cell.** Its first click would otherwise have cleared it.
- **The Cell box has a ✕** (shown on hover, beside the list's arrow) that empties it: no cell. The box looks as before.
- Esc while typing in the box still puts back the open cell and leaves the box; the browser's own Esc, which empties a search box, is stopped so the box never shows nothing while a cell is open.

## 2026.11.21.7

Step 2 of the plugin's window: the four actions no longer move down when the "Next:" sentence gets longer. No measurement, response call or saved table changes.

- **The Next line and the line under Picture each keep room for two lines.** In 2026.11.21.6 the block of four actions sat 17 or 34 pixels lower whenever the Next sentence took two or three lines instead of one. It now stays where it is in every state.
- **The two Next sentences that took three lines are shorter:** "click the dashed outlines that are real cells (43 left), Add cells that respond… (no outline yet), then Continue to measure." and "these 177 cells are measured. Add cells that respond… (no outline yet) and measure again, or go on to Analyze responses."
- **"Hidden on the picture now: …" is under the Show switches,** where things are hidden and shown, not at the end of the Next sentence. It takes no room while everything is shown.
- In a window made narrower than usual a sentence can still need a third line; it gets it, and what is under it moves down by that line.

## 2026.11.21.6

The plugin's window in Fiji, polished: Step 2 "Choose the cells" most of all. No measurement, response call or saved table changes; what is recorded about how cells were chosen is unchanged.

- **Text uses the card's width.** Every sentence used to wrap at one fixed width, short of the card, and left a strip empty on the right. It now fills the card and follows the window when you resize it.
- **The four actions keep their places:** Suggest cells and Accept likely on one line, Add cells that respond… and Continue to measure on the next. The next one to press is blue; one that cannot be pressed now is gray instead of gone. Each column is as wide as its longest label, so nothing moves when a label or the blue button changes. Before, the blue button sat in a different place in each of five states.
- **Picture is one joined switch:** All frames | Before response | Response, like "Average | One cell + average" on the Step 4 page. The line under it still says what the choice means (blind, magenta). The names lost their brackets and the ellipsis; Before response still asks which frames.
- **The counts come straight after the actions,** then the tools for your cells on two lines: Remove selected, Undo and Draw by hand; then "Click cells to remove" and "Leave out cells within 40 px of the edge". The report of the last detection, with Details and Log, comes after them.
- **Show: one switch per kind on the picture** (your cells, likely a cell, unsure, background, precipitate?, responds), each with a tick and a sample of how it is drawn. They were legend entries that looked like plain labels. "Show suggestions" and "Show background", which did the same job a second time, are gone; "Draw by hand, what is shown" is now "Draw by hand".
- **No empty box.** The trace frame takes its room when there is a trace to show. After Add cells that respond it used to stand empty, saying "Traces are hidden while choosing cells blind" although the picture was the response.
- One left edge for text, buttons and ticks (buttons sat 6 pixels in; Details and Log another 2 pixels off). One tick box everywhere in the window. "Next:" is the window's own blue.
- The recording's name is cut in its middle to the width there is (it was cut to 40 characters); hover for the whole name.
- The window opens with room for its scrollbar, so one that has to scroll down has no bar along the bottom.

## 2026.11.21.5

The line on the Step 4 plot that marks the picture's frame is shown only while the frame moves. No measurement, response call or saved table changes.

- **Nothing is drawn on the traces while nothing moves.** The solid vertical line used to stay at the frame the picture shows.
- **While the frame moves, the line shows where the picture is, with that frame's numbers:** a drag in the plot, the slider or the step buttons under the picture, the arrow keys, Play. About a second and a half after the frame stops, the line goes.
- With the pointer over the plot, the dashed line and the numbers under the pointer show as before.

## 2026.11.21.4

The numbers on the Step 4 plot follow the pointer. No measurement, response call or saved table changes.

- **No box sits on the traces any more.** The dark box with the time and the values used to stay at the frame shown in the picture (at 0.0 s when nothing had been moved), over the traces, whatever was clicked.
- **Over the plot, the numbers are those under the pointer:** the time, each plotted mean and the selected cell, beside a dashed line. They move with the pointer and go when it leaves. Looking does not move the recording.
- **While you drag in the plot or the movie plays,** the numbers are the frame's, at the solid line, as before. When the drag ends they go back to the pointer.
- The solid line still marks the frame the picture shows.
- The box moves to the top of the plot when a plotted value would be under it at the bottom.

## 2026.11.21.3

A figure editor for the plot and the heatmap on the Step 4 page. No measurement, response call or saved table changes.

- **Edit figure… in the plot's and the heatmap's Export menu** opens the figure with everything that can be set beside it. The preview is the picture that gets saved, drawn again at every change: what you see is the file.
- **Size** in inches (1 to 20 each way), with resolution and format; the pixel size is shown. A figure too large for a browser to draw (over 80 million pixels) cannot be saved, and the editor says so.
- **Text:** title, X and Y axis titles (the heatmap's colour bar title too), font (the page's, Arial, Times New Roman, Courier New), text and title size in points, one colour for text and axes, white or transparent background.
- **Plot axes:** minimum, maximum and tick spacing for X and Y (empty is automatic), a box or left-and-bottom frame, grid lines on or off, axis thickness.
- **Plot lines:** thickness and colour of the mean, the selected cell and all traces; colours for responded and did not respond; SEM band on or off. Only what the plot on screen shows is offered.
- **Plot marks:** legend, stimulus line and its label, baseline, response window, thresholds.
- **Heatmap:** colours, range and sort (the page's own lists, shown here too), a colour scale typed in so several recordings can share one, colour bar, window and stimulus marks, the response call strip. A heatmap figure now carries its title and a key to the call strip.
- **Kept for the next recording:** size, font, sizes, colours, thicknesses and what is shown. The title, the Y axis title and the axis numbers belong to one recording and are not carried over. Reset returns one figure to how it was.
- "ΔF/F₀" and "Ca²⁺" are drawn correctly in every font: a font without small lowered or raised digits no longer borrows them from another.
- A plot or heatmap saved without opening the editor is 8 × 6 in, as before. The heatmap used to keep its shape on screen; it now has the size set in the editor.
- The file kept beside each saved figure (`.metadata.json`) records every setting of the editor.
- Options… stays in the picture's Export menu for pictures and frames. Format and resolution are shared with the figures.
- Not in this build: SVG or PDF files, error bars other than the SEM band, a legend inside the plot.

## 2026.11.21.2

Reviewing flagged cells on the Step 4 page. No measurement, response call or saved table changes.

- **The flag list names every flagged cell.** One row per cell, with its call and what it is flagged for in two or three words: "Below zero", "Low F₀", "Noisy baseline". Before, the list gave only a total per flag, and those overlap (on P2 Run1: 8 + 4 + 4 for 9 cells).
- **A click on a row opens that cell.** Show goes to the flagged cells, the list closes, and the ← and → keys step through them.
- Under the rows each flag is explained once; a click on one still shows only its cells. The bullets that sat on the last line of each sentence are gone. So is "Show cells needing review": the Show list does that.
- **The Cell box takes typing.** Type a cell's name or pick it from the list; a number alone and Enter opens that cell ("12" opens Cell 12). An empty box is no cell. It replaces a list with an entry for every cell (177 on P2 Run1).

## 2026.11.21.1

**This build changes results.** Four choices on how responses are counted and shown, each made after it was measured on the P2 Run1 and P3 Run2 analyses.

- **A response must hold for 2 samples in a row.** The length is counted in samples, not seconds. "For 1 s" was applied as 2 samples at exactly 1 frame per second and as 3 at the file's own 0.99995 s per frame, so two recordings of one series were held to different rules:
  - P2 Run1 (3 samples): 83 % responding; with 2 samples about 85 % (146 of 172).
  - P3 Run2 (2 samples): 65 %; with 3 samples about 61 % (142 of 231).
  New analyses use 2 samples in a row. Settings remembered for a recording that still say "for 1 s" move to it once and the page says so; a length you set by hand is kept. The rule line says what is applied: "for 2 samples in a row", or for a length in seconds "for 1 s, here 3 samples in a row". **Running an earlier analysis again can change its calls.**
- **All cells lead the headline**, on every tab and in the line Fiji shows when the analysis finishes. The cells chosen blind and those added after the response was shown are the second line. Before, the blind figure led and the Responder counts tab showed only the all-cells one.
- **A line for quiet cells.** The % among cells that need a rise of 30 % or less to respond, with noisier cells set aside whether they responded or not: 84 % (113 of 134) on P2 Run1 and 68 % (121 of 179) on P3 Run2. Setting aside only the noisy cells that did not respond would lift the figure to 87 % and 72 %, because noisy responders stay in; that is not done.
- **A switch for cells below zero.** A cell whose background-corrected fluorescence went below zero (ΔF/F₀ under −1 at some frame) stays counted. "Leave out N cells below zero" recalculates every percentage without them: 8 cells on P2 Run1 and 18 on P3 Run2, at most 1 point of % responding.
- One tally, worked out in Fiji, feeds the headline, a table at the top of the Responder counts tab and a saved file, `<recording>_percent_responding_<analysis>.csv`: each row with every cell and without the cells below zero. It replaces `_responders_by_how_chosen_`, whose rows are all in it.

## 2026.11.20.23

The Step 4 page only. No measurement, response call or saved per-cell table changes.

- **The metrics summary lists each metric once**, with a plain name, its unit and where it is measured: over the whole recording or the response window, from recording start or from the stimulus. Saved columns that hold the same number for every cell share one row, and "Saved as" names them. Which columns do is read from the numbers, because it depends on the analysis. On the P2 Run1 and P3 Run2 analyses of 20 Nov:
  - Peak_dFF0 = Response_Amplitude, Peak_increase = Window_peak_local_fraction and Onset_s = Run_start_s in both;
  - all four peak columns are the same only in P2 Run1;
  - AUC and Window_AUC_fraction_s are close but not the same (-18.996 against -18.991 on P2 Run1), so they stay two rows.
- **The timing and shape of a response are summarised for cells that responded.** Onset, latency, peak time, rise and width of a cell that did not respond are the timing of noise: on P2 Run1 the 28 such cells averaged an FWHM of 2.1 s against 18.7 s for the cells that responded. "Did not respond" shows "–" in those rows, and "All cells" counts the responders there (its n says so). Peak, signal-to-noise, area and baseline are still given for every cell.
- The summary's exported file (`metrics_summary`) gains the columns Unit, Measured and Saved_as. Two columns that say how much of the response window had data (Window_AUC_covered_s, Window_AUC_coverage_fraction) are no longer listed as metrics. The per-cell tables on disk are unchanged.
- In a saved table, pointing at a column name says what the column is. The metric names stay in view when the summary scrolls sideways.
- **Settings labels never wrap.** A label sits beside its control at the panel's usual width and above it when the panel is narrower than 344 pixels. Before, "Stimulus added after frame" took four lines in a narrow panel. Three labels are shorter: "Stimulus after frame", "Duration measured in", "Baseline samples".
- **Frames per second is shown to four significant figures** (1.000, not 1.0000537066691941). The stored value keeps every digit, and the line under the box says when the value is the file's own.
- **The plot zooms with Ctrl or ⌘ + scroll**, like the picture. A plain scroll over it moves the page.

## 2026.11.20.22

The Step 4 page only. No measurement, response call or saved table changes.

- **One shape and one size for the page's controls.** Every button, list and field is the same height (28 pixels) with the same text size, and every box (cards, tiles, menus, buttons, fields) has the same corner. Before, the controls came in eleven heights from 24 to 34 pixels, and the corners in four sizes plus pills and a round dot. What that changes on screen:
  - the "cells flagged" list opens from a button like Export (it was a pill);
  - the outline-thickness list looks like the other lists (it had no border);
  - the headline box has the same edge as the other tiles (it had a thick blue edge that curved);
  - the dot beside a cell's call is a square, like the swatches in the legend;
  - Play, the frame arrows, the cell arrows, "Open results folder" and the export buttons no longer differ in height.
- **The results row is level.** The sentence with the percentage has a row of its own. Under it the counts, "cells flagged" and "Open results folder" sit on one line at the same height. "Cells flagged" used to sit 10 pixels lower than its neighbours, on a leftover second row beside a taller tile, because a spacing rule meant for the settings also reached it; the Export buttons were 5 pixels low for the same reason.
- **The settings arrow has a strip of its own.** With the settings hidden, the arrow sat on the left border of the results box, against the percentage. It is now in the strip between the settings and the work area, above the drag grip. With the settings hidden that strip stays, in the page margin, so the arrow is never on the results. It also stays in view when the page is scrolled.

## 2026.11.20.21

The movie views, the picture's controls and the export of frames. No measurement, response call or saved table changes.

- **"Response in magenta" shows the response again, and "Change from baseline" uses the same measure.** Since 2026.11.20.15 magenta only started about 2.8 gray levels above "no change" on dim recordings, where cells rise by 1 to 3, so almost nothing was painted. Both views now show how far a pixel rose above what the background around it did, averaged over 9 × 9 pixels and 5 frames, on a scale set by the recording's own noise:
  - P2 Run1 ATP: 97 of 177 measured cells are magenta at the peak and 1 before the stimulus;
  - P3 Run2 ATP: 9 to 12 % of the picture during the response, about 1 % before and after;
  - two recordings without a response (C3 Run2, P2 Run2 capsaicin): never more than 1.8 % of the picture.
  The fading of the field and bright objects that fade do not count as a response. Display only.
- **GFP view.** The recording with the background removed, in green: the green of "Response in magenta" without the magenta.
- **No cell selected is a choice.** The Cell list starts with **None** and the page opens that way; Esc deselects. The plot follows the selection: with a cell it shows that cell over the average, with none the average alone ("Responders vs non-responders" stays when you chose it).
- **Outline thickness sits beside "All outlines"** in the picture's toolbar, where the outlines are switched on and off. It is no longer inside Export → Options.
- **Save frames also writes one picture with all the saved frames** (`all_frames_combined.png`), laid out in a grid from the same labelled frames as the single files.
- **The settings panel folds away from its own edge** (the small arrow between the panel and the work area). The "Hide settings" button in the header is gone.
- Scrolling over the picture scrolls the page; hold Ctrl or ⌘ (or pinch) to zoom the picture.
- "Did not respond" is orange on the plot as on the picture (it was brown on the plot).
- Export menus say where the file goes: "Save to plots folder" and "Download a copy". "Results folder" reads "Open results folder".
- Fiji: the question before "Accept likely" says about 7 in 10 of those outlines were real cells (it said 8 in 10; measured 71 % with the grading of 2026.11.20.20), and the line under the buttons is short when some pictures reuse saved outlines.

## 2026.11.20.20

Fewer cells are left for you to pick up by hand. Measured by running this build on the four audited ATP recordings (837 audited cells; cells chosen blind on the frames the plugin offers). No measurement or response threshold changes.

- **"Accept likely" takes more of the cells.** An outline is now also graded "likely a cell" when it is cell-sized (500 pixels or more), stands out from its surroundings by more than 0.3 gray levels, is not bright on average (8 gray levels or less) and has less than 60 % very bright pixels. Until now a bright dot inside (2 % of its pixels) made a cell "unsure", and 25 % bright pixels made it "very bright".
  - Accept likely took the outlines of 232 audited cells; now 360.
  - Of the outlines it accepts, 73 % were on audited cells; now 71 %.
  - Unchanged: an outline brighter than 8 gray levels on average stays "very bright", and a small one or one touching the picture's border is not accepted in bulk.
- **A third blind picture when Cellpose uses a graphics chip.** Besides the pictures shown 0–2 and 0–8 gray levels, Cellpose also gets the one shown 0–4, and its outlines are added by the same rule (less than 10 % shared with every outline kept so far). 567 audited cells get an outline instead of 547, for 6 % more outlines and about 18 seconds. Without a graphics chip the third picture is not run, because each picture takes minutes there.
- Together, audited cells that end in a "likely" outline went from 243 (29 %) to 380 (45 %). Still open: 270 audited cells (32 %) get no outline at all, and 89 (11 %) sit in outlines hidden by the edge rule.
- The grading rule was chosen on these same four recordings, so expect a little less on others. P3 Run2 is the weakest: 56 % of its accepted outlines are on audited cells (58 % before).

## 2026.11.20.19

- **Responding cells are no longer walked past.** Blind selection alone leaves cells that only show when they respond without an outline: on P2 Run1 ATP it gave 65 cells, and 92 more were found once the response was shown. In a recording with a response:
  - after the blind cells are accepted, the highlighted next step in Step 2 is **Add cells that respond…** (it was Continue to measure);
  - **Continue to measure** with blind cells only asks once: "Add cells that respond" or "Measure these only";
  - after a blind-only measurement, Step 3 has **Add cells that respond…** beside Back to cells, and the sensitivity check's line points to it.
  The blind cells stay recorded as chosen blind, and Step 4 still gives % responding for them alone, beside the added cells and all cells (P2 Run1 ATP: 78 % of 65 blind, 96 % of 89 added, 88 % of 154). A recording with no response found is unchanged.

Step 4 page: six things moved to where you look for them. No measurement, threshold or saved file changes.

- **Playback is directly under the picture.** Play, the frame slider, previous/next frame and Speed sit in the picture's own card, under the picture they control.
- **One cell filter.** The result figures at the top (% responding, registered, excluded, unscorable) are figures only. **Show**, beside **Cell**, is the one place that chooses which cells are listed, drawn and averaged.
- **The selected cell's call is on the plot card.** "Cell 12: Responded · largest rise …" and its quality flags sit at the top of the plot they describe, instead of above both panels.
- **One Export menu per panel**, in the panel's corner: Save, Download, Save frames… and Options… for the picture; Save, Download and Options… for the plot and the heatmap; Save and Download for the tables. Choosing an entry closes the menu.
- **Quality review is a chip** in the row of result figures ("19 cells flagged"). Its list opens over the page instead of pushing the data down.
- **A warning or error sits under the setting it is about**, with its one-click repair. The Run box only says how many there are ("1 warning above", "Fix 2 problems above") and takes you there. A message about a setting folded away in Advanced is still listed in the Run box.
- The picture fills its card in a short window (there was an empty strip under it), the two cards are a little taller now that playback sits inside the picture's card, and "1 frames/s" reads "1 frame/s".
- "Enter the frame after which it was added" no longer flashes for a moment when the page opens on a recording with recorded frame times.

## 2026.11.20.18

- **Blind selection uses the frames before the response, and the plugin offers them.** Until now it asked for "the last frame before the stimulus" and offered frame 20 when it found no response. Twenty frames are too few to see faint cells: on one recording the two blind pictures outlined 70 % of the audited cells on frames 1–20 and 79 % on frames 1–50, and the cells gained responded as often as the ones found before. The window offered now:
  - with a response found: up to 25 frames before the rise (it was 9). The rise the plugin finds is the field's; on P3 Run2, C2 Run1 and P2 Run1 the first cells start their main rise 9 to 15 frames before it;
  - with no response found: up to frame 60;
  - without the frames at the start in which the field moved or was still settling. A "rise" that is only the field moving counts as no response.
  You confirm or change the frames as before. They need not end where the stimulus was added; Step 4 asks for that separately.
- Wording follows: **Before the response only (blind)** on the first screen, **Before response (blind)…** in Step 2, "Frames before the response" in the dialog. The saved record reads `Before the response, blind (frames a-b)` (it was `Pre-stimulus, blind (frames a-b)`).
- Step 4 page: the frames the cells were chosen on are shown as "Cells chosen in Fiji on frames a–b, before the response." The page no longer warns that they run past the stimulus, and the starting stimulus frame taken from them is called a placeholder.
- Offered on the 18 Nov recordings: P3 Run2 1–97, C2 Run1 1–58, P2 Run1 1–54, C3 Run2 1–60, P2 Run2 capsaicin 1–60, C1 Run3 KCl 13–60 (frames 1–12 moved; the rise near frame 6 is that movement).
- Checked against the audited cells with margins of 15, 20 and 25 frames. Cells covered by the two blind pictures with the frames now offered: P3 Run2 66 %, C2 Run1 67 %, P2 Run1 56 %, C3 Run2 80 %, P2 Run2 capsaicin 95 %, C1 Run3 KCl 69 %. On frames 1–20 they were 42 % (P2 Run1), 70 % (C3 Run2), 77 % (P2 Run2 capsaicin) and 66 % (C1 Run3 KCl); with frames set by hand up to 9 before the rise, 68 %, 68 % and 59 % on the first three. The share of found cells with a clear rise is the same at every margin. No cell starts its main rise inside the frames: the first does so 10 to 16 frames after them (with a 15-frame margin it was 1 frame on P2 Run1, with 20 it was 5). Not changed by any of this: cells that are missed still respond more often than cells that are found (P3 Run2 71 % against 56 %), and longer windows give more suggestions that are not cells.
- The frames end before the main response, not before every response: a few cells respond briefly earlier (1 to 4 audited cells per recording at frames 22–60). The dialog says so; end the frames earlier to leave those out.
- Step 4 page: a movie frame that cannot load is no longer asked for again without pause (it was requested hundreds of times a second until the slider moved). It is left alone for 2 seconds; moving the slider or pressing Play asks again at once.

## 2026.11.20.17

- **Movement is measured against the frames the cells are chosen on.** Build .16 compared every frame with where the field sits for most of the recording, so on a recording that drifts slowly it called the first hundred frames "moved", although those are the frames the outlines were drawn on. While cells are chosen blind, the message now reads, for example, "The field sits about 4 px from where the cells are chosen (frames 1–20) from frame 139 on."
- The shift is estimated between the 4-pixel steps of the check, with the smooth lamp light and haze taken off first, so a slow drift is no longer reported as a sudden step.
- Worded by size: from 2.5 pixels the field "sits about N px from…" or "drifted by about N px"; from 8 pixels it "moved" and the outlines are no longer on the same cells. Under 2.5 pixels nothing is said. Only movement of 8 pixels or more is left out of the resting picture and offered a repair on the Step 4 page.
- When the frames the cells are chosen on moved themselves, the message says that their picture is smeared (C1 Run3 KCl: 18 px during frames 1–6, settling by frame 15, chosen on frames 1–20).
- A run of noisy frames is reported as one piece ("34 frames between 2 and 119") instead of a list of fragments.
- Measured with it: P2 Run2 capsaicin drifts about 4 px over the recording, C2 Run1 about 6 px with a jolt at frame 8, P3 Run2 and C3 Run2 under 2 px.

## 2026.11.20.16

- **Movement check.** When a recording is opened the plugin measures, frame by frame, how far the field is displaced from where it sits for most of the recording (on a picture reduced four times; shifts of 4 pixels or more count; about 3 s for 300 frames). Nothing is corrected. When frames moved, the Recording card, the log and `what_was_done.txt` say which ("The field moved (about 16 px) during frames 1–10").
- Step 4 page: when F₀, the baseline or the response window include moved frames, the page says so and offers the frames after the movement in one click ("Use 11–20"). The automatic F₀ rule leaves moved frames out when at least 3 frames remain. The note under the movie says when the frame on screen is a displaced one.
- The resting picture behind "Response in magenta" and "Change from baseline" is made without the moved frames, so a movement at the start no longer leaves edges on every later frame.
- **Blind detection: outlines of the first picture that are never shown no longer hide outlines of the second.** The merge compared every outline of picture 2 with every outline of picture 1, including those the size rule drops afterwards. On one recording a single 6445-pixel outline hid six cells and an 11-pixel speck hid a 1731-pixel one. Now only picture-1 outlines that are shown (150 to 6000 pixels, not a strip on another outline) count. Re-run on saved outlines: 7 cell-sized outlines recovered on P2 Run2 capsaicin, 4 on P3 Run2, 4 on C2 Run1; none lost.
- Very bright suggestions ("precipitate?") are drawn as a red 2-pixel line. The orange 1-pixel line ran along the bright edge it traced and could not be seen.
- Not in this build: measurements of moved frames are still taken at the outlines' fixed places (no registration), and Fiji's own pictures (all frames, before the stimulus, response) still average moved frames. The message says so.

## 2026.11.20.15

Step 4 page, stimulus timing:

- **One input, in frames: "Stimulus added after frame".** F₀ frames, the baseline and the response window follow it and are shown beside it as dashed chips ("1–20", "same as F₀", "21–300 (end)"). A click on a chip gives that row two frame boxes of its own; **auto** puts it back. **Use frame N** sets the stimulus to the movie frame on screen.
- Windows are defined by frame. The seconds sent to the analysis are taken from the recording's own frame times (the stimulus is the time of the first frame after it; a baseline runs up to, not including, the frame after its last frame; the response window runs to its last frame, both ends included), never from a typed round number. The same entry gives the same frames whatever the frame interval, and "to the end" always means the last frame. The protocol file format is unchanged.
- New starting rule: F₀ is the last 20 frames before the stimulus (all of them when there are fewer), and the baseline is the same frames as F₀. Before, the starting baseline was every frame from the start of the recording.
- The frames used to choose cells in Fiji are shown as their own line ("Cells chosen in Fiji on frames 1–113") and no longer decide the stimulus once an analysis has been run: the next recording starts from the stimulus frame of your last analysis ("Same as your last recording"). When the last session's settings for this recording are filled in, the page says so and offers **Use Fiji's values**.
- When the frames used to choose cells run past the stimulus, the page says so: "Cells were chosen on frames 1–74. Frames 21–74 are after the stimulus, so these cells were not chosen blind to it."
- A row set separately that is wrong gets one message beside it, naming the frames, with a one-click repair ("baseline frames 10–12 come after the stimulus. Use 1–9"). A setting that guarantees no call for any cell (a baseline of fewer frames than the rule needs, a response window past the last frame) now stops the run instead of warning that cells "may" be unscorable. A control that was just ticked asks for its one number in one message.
- The tick box reads back what it confirms ("Timing checked: 1 frame/s, stimulus after frame 20"). When the typed frame rate differs from the recording's own, the page shows the file's rate with **Use it**.
- A protocol import says what it did, and unticks the box only when the timing changed. A run in which no cell got a call says so in the headline and the status line, without "0 of 0 evaluable cells" or "including 0 responders".

Step 4 page, layout and wording:

- The settings column is one **Timing** card that fits without the 300-pixel scrolling box: every control the same size on one right edge; the rule is one line that opens its settings; stimulus name, rule type, protocol file, recorded timestamps, diagnostic files, **Show explanations** and **View log** are under **Advanced**; condition and biological unit are under **Labels for the output**. The Run box is the tick, the button and one status line.
- Shorter words throughout: "response window" (no "strip along the top"), "baseline", "Background removed", "Raw gray levels", "Done." The legend, menus, notes and messages lost their parentheses and second sentences. One name per thing: baseline and response window, in the settings, the legend and the messages.
- The buttons under the movie (Save picture, Download, Save frames…) are dark with light text; they looked switched off. **Save frames…** opens a small panel on the page instead of the browser's pop-up boxes. Pictures and frames follow **Options** (width in inches × dpi, up to four times the recorded pixels), and the picture has its own Options button; outline thickness moved there.
- **Y axis** boxes beside the plot buttons set the plot's vertical limits (empty = automatic), on screen and in saved plots.
- One **Signal** setting for the plot and the Processed data table. "All traces" sits with the plot buttons. The Results folder button is beside the results; "Analysis engine ready" shows only when the engine is not ready.
- **Background areas** in the View menu: the pixels taken as each cell's background, in blue, on the all-frames picture.
- Each view has colour swatches under the picture (green and magenta for "Response in magenta", and so on) instead of a sentence.
- **"Response in magenta" and "Change from baseline" no longer show noise as a response.** A change is drawn at full strength only from twice what a resting frame shows by itself, so a recording in which nothing responds stays dark. A brightening beside an equal darkening (a bright object that shifted) is taken out. Both views count from each frame's own middle value. Display only.

Choosing cells:

- **Outlines too thin to be a cell are no longer accepted by themselves.** Cellpose sometimes draws a strip along the rim of a cell it has already outlined (seen as cells 78 and 79 on one cell). An outline under 4 pixels wide on average, or with roundness under 0.2, is left out when it lies along another outline, and otherwise kept as "unsure" (never accepted without a click). The log gives both counts.
- **Cell IDs are the same in every session.** A new outline's ID is worked out from the outline itself (its place and pixels), so tables from two sessions on the same outlines can be matched by cell. Before, every Fiji session gave the same cells new random IDs.

Also in this build: the linked-timing work of the cloud session is merged (its seconds boxes are replaced by the frame editor above; its counting of baseline samples and its windows under the frame slider are kept). The plugin window stays an ordinary window.

## 2026.11.20.14

- **Save options for pictures** on the Step 4 page (Options, beside Save plot and Save heatmap): format **PNG or TIFF**, size (6 × 4, 7 × 5, 8 × 6, 10 × 7.5 or 12 × 8 inches, or as on screen), resolution **150, 300 or 600 dpi**, and TIFF compression **LZW, ZIP or none**. The resolution is written into the file, so it opens at the chosen size. A saved plot now carries its title and legend and leaves out the time cursor. Download gives the same file as Save. The choices are remembered. PDF and SVG are not offered yet.
- **The colour picture is one colour photo.** `..._picture_colour.tif` used to be a three-channel stack that only ImageJ shows in colour; it is now an ordinary colour image with the same display ranges. A second file, `..._picture_colour_cells_outlined.png`, has the measured cells outlined in white, a key to the colours (green cells, yellow very bright objects, magenta response) and a 100 µm bar.
- **Save picture** on the Step 4 page saves the whole picture in the view on screen (the colour overview or a movie frame), with the outlines shown and a key to the colours and the outline colours, in the chosen format and resolution.
- **Save frames…** on the Step 4 page writes every Nth frame (you type N; 50 to start with) of the movie view on screen as its own picture, in a folder of its own inside `plots`. Each carries its frame number, time and view in the corner and the outlines shown on the page, in the format, resolution and compression chosen under Options (for example LZW TIFF), with a small file saying what the outline colours mean.
- A small menu beside "All outlines" sets the thickness of the cell outlines (thin, normal, thick, very thick), on screen and in saved pictures and frames.
- **The heatmap has a legend**: the colour scale with its range, the strip beside the Cell axis (responded, did not respond, excluded or unscorable), the gray rows, the selected-cell marker and the stimulus line. The call key is also drawn into the saved heatmap.
- **ΔF/F₀ of raw values** in the Signal menus (plot and Processed data): each cell's raw trace against its own mean over the F₀ frames, with no background removed. For display and export; the response calls still use the background-corrected ΔF/F₀.

## 2026.11.20.13

- **Processed data tab** on the Step 4 page, as in SimpleCa: one row per frame, one column per cell, for ΔF/F₀, gray levels above background or raw gray values, 50 frames at a time on screen, with **Save table** (the whole table as CSV into the `plots` folder) and **Download CSV**.

- Step 4 page, smaller fixes: the page says when its settings were filled in from Fiji or from the last session on this recording; "Saved in plots…" notes clear when a new analysis finishes; a save that cannot reach Fiji says so; a one-cell plot's file name carries the cell ("Cell_12" and the start of its ID) whatever signal is shown.

- **Metrics summary** on the Step 4 page, as in SimpleCa: the Metrics tab opens on mean ± SEM and number of cells for every metric, for all cells, by response call, and by how the cells were chosen. **Save table** and **Download CSV** work for it (with SD) and for every other table on that tab.
- The top of the Step 4 page is one compact row (the heading goes once there are results, the tiles are a single line), and a handle under the movie and the plot sets their height. The routine playback hint is with the explanations.

Step 2, from a look at the window after suggestions were accepted:

- **Add cells that respond…** and **Continue to measure** sit under Suggest and Accept, so the step the "Next:" line names is on screen without scrolling.
- When the outlines are ready, the progress block shrinks to its result line; stages, timings and settings stay under Details. The result is no longer repeated further down.
- The empty trace frame is hidden while cells are chosen blind.
- The count names the orange outlines "very bright (precipitate?)", as the legend does, and the legend has no gap when "responds" is hidden.

## 2026.11.20.12

- The heatmap opens in **plasma**, the colours SimpleCa uses, instead of blue–white–red. Values below zero then take the darkest colour; choose Blue–white–red to see them.
- The Step 4 plot opens on the **Average** alone, chosen with three buttons (Average, One cell + average, Responders vs non-responders) instead of a menu. In 2026.11.20.11 the average was there but behind a menu that still opened on one cell, and it was missed.

## 2026.11.20.11

Step 4 page:

- **The plot opens on the average alone.** Three buttons choose what it shows: **Average** (mean ± SEM of all cells, or of the group chosen under Show, and nothing else; the default), **One cell + average** (the selected cell with its thresholds, as before) and **Responders vs non-responders**. Picking a cell never changes what the plot shows. "All traces" is off until ticked.
- **Corrected or uncorrected signal.** Signal switches the plot between ΔF/F₀, gray levels above background and raw gray values (not background corrected), read from the measurement's own tables. Excluded cells, which have no ΔF/F₀ trace, can be seen this way.
- **Save plot** writes the picture, the numbers behind it (CSV) and how it was made into a `plots` folder beside the tables and says where, with an Open folder link. Download still goes through the browser.
- **Two more movie views.** "Recording, background removed" takes the smooth background off every frame, and "Response in magenta" shows that in green with the rise above the baseline frames in magenta, so responding cells turn pink. Each frame's own median counts as "no change", so the slow fade of the field does not hide the response. Display only.
- **Resizable panes.** Drag the edge between the settings and the work area, and between the movie and the plot; double-click to reset. The sizes are remembered.
- **More room for the data.** The movie and plot get at least 54 % of the window's height (the page scrolls when needed), the tiles are smaller, a zero count takes little room, the green "Showing the completed analysis" banner appears only when settings were changed afterwards, and the Run box is shorter so more of the settings is visible.
- **The headline shows how the cells were chosen.** When cells were chosen blind and others added afterwards, the first tile leads with % responding of the blind cells, with the added cells and all cells beside it.
- **Fewer false alarms.** The note that the 10 % rise level lies below the baseline SD concerns the rise-time metric, not the call; it no longer counts as a quality flag (on P3 Run2 it flagged 149 of 238 cells). A warning worded two ways is shown once.
- The plot's legend names everything drawn (baseline band, response window, the run that counted, each dashed level), the response window is a strip along the top instead of a fill, and a dashed line marks when the stimulus was added.
- Tables show "Cell 1" instead of the long ID (the ID is in the tooltip) and leave out columns that say the same on every row.
- Explanations under the settings and plots are hidden by default; "Show explanations" in the header brings them back.
- "Discard pending changes" keeps the confirmation tick when it returns to the settings of the completed analysis. Cancel shows only while an analysis runs. The response-window chooser shows only when there are two windows.

Plugin:

- The ROI registry records each cell's own selection basis (chosen blind, or added afterwards), and the measurement's basis says "Blind on frames a-b, then response-assisted: n chosen blind, m added afterwards". Before, every cell of such a measurement was listed as response-assisted.
- Cells count as chosen blind when they were held, unchanged, at the moment blind choosing ended: when the response or the all-frames picture was first shown.

## 2026.11.20.10

- **Step 2 says what to do next.** A "Next:" line at the top of Step 2 always names the one thing to do, from what has been done so far: Suggest cells, Accept likely, click the dashed outlines that are cells, then **Add cells that respond…** or **Continue to measure**.
- **Add cells that respond…** (blind only, beside Continue to measure) shows the response and looks for the cells that only show when they respond, without a second press. Switching Picture to Response (magenta) after blind suggestions does the same. Before, the outlines on screen stayed the blind ones until Suggest cells again was pressed, which looked finished and was not.
- Switching to another picture and suggesting again reuses saved outlines of that picture; a new Cellpose run is made only when the outlines on screen came from the same kind of picture.
- Legend kinds that were switched off come back on when new suggestions arrive, and the Next line says which kinds are hidden.
- With the response shown, the button reads **Accept responding and likely (n)**.
- The result line for blind suggestions is short again ("382 suggestions from two pictures; Cellpose took 20 s and 19 s (GPU)"), so it is no longer cut off. The record keeps both commands and timings.
- The question shown before the response appears now says what happens to the cells already chosen and to the ones added afterwards.

## 2026.11.20.9

- **Responding cells are no longer graded away.** An outline that rises on its own at the stimulus is now graded **responds** even when part of it is very bright (precipitate does not respond). Before, such cells were graded unsure or "precipitate?" and left out by **Accept likely**. On P3 Run2 ATP, 17 of 58 clearly responding areas had no chosen cell after Accept likely; with the new rule 6 are left (4 at the picture's edge, 1 whose rise is no larger than its surroundings, 1 real precipitate). The rise must also be at least 5 % of the outline's brightness. Checked on this one recording only; blind grades are unchanged.
- **Blind and response-assisted choosing can be used together.** Choose cells blind first, then switch to **Response (magenta)** and add what is only visible when it responds. Each cell is recorded as chosen blind, or as added (or reshaped) after the response was shown (`..._cells_chosen_blind_or_after_response.csv`). When Step 4 finishes it gives % responding for each group and for all cells (`..._responders_by_how_chosen_<run>.csv`). Only the blind group is free of selection on the response.
- **Step 4 starts from what you confirmed in Fiji.** The frames confirmed for blind selection become the baseline, and the stimulus is taken to follow their last frame (a range longer than 30 frames gives the last 20 as the resting baseline). Without them, the frames before the rise the plugin found are used, as before. The page says where its starting values come from.
- Step 4 remembers the settings last used for a recording (`step4_last_settings.json` beside its tables) and fills them in the next time that recording is opened with the same frames from Fiji.
- Fix: the Step 4 warning "the response baseline holds 0 samples … every cell will be unscorable" appeared on every recording, whatever the baseline. The samples were miscounted; the analysis itself was not affected.
- The frames dialog for blind selection says what the frames are used for, and offers the longer range when the range typed ends 20 or more frames before the one the plugin found (fewer frames give a noisier picture and fewer cells).
- **The legend is a switch.** Click a legend entry (your cells, background, responds, likely a cell, unsure, precipitate?) to hide or show that kind on the picture; double-click to see only that kind, and again for everything. Hidden suggestions cannot be accepted by a click.
- The plugin window is an ordinary window again: clicking a picture or another program puts it behind. It no longer stays above Fiji's other windows.

## 2026.11.20.8

- **Blind detection uses two pictures.** Cellpose now runs twice on the pre-stimulus picture: once shown from 0 to 2 gray levels and once from 0 to 8, brighter pixels capped. Every outline from the 0-2 picture is kept. An outline from the 0-8 picture is added when it shares less than 10 % of the smaller outline's area with every 0-2 outline. Neither picture holds anything of the response, and the grades are unchanged.
- Why: 0-2 alone was best only on P3 Run2, where it was chosen. On four ATP recordings of 18 Nov 2026 (real Cellpose), of the response-assisted outlines that could be cells, the two pictures together found 68 / 95 / 55 / 62 % (C2 Run1 / C3 Run2 / P2 Run1 / P3 Run2). That compares with 46 / 94 / 34 / 50 % for 0-2 alone and 57 / 86 / 45 / 36 % for 0-8 alone. Of the audited cells it found 167 / 128 / 114 / 156 of 240 / 164 / 210 / 223. A cap taken from each picture's own 90th percentile was never better than the best fixed cap. The cost is 25-55 % more suggestions to review on C2, P2 and P3.
- The cells blind detection misses still rise clearly more often than the ones it finds: by about 10-12 points on P2 Run1 and P3 Run2 under every cap tried, and by 4 points on C2 Run1 with the two pictures.
- Each picture has its own folder in `automatic_outlines_blind/`, and saved outlines are reused per picture. `detection_picture.txt` records how many outlines came from each picture, the evidence table has a `found_on_picture_gray_levels` column, and `..._suggestions_rois.zip` holds every suggestion, named as in the table. The record (`what_was_done.txt`) gives both Cellpose commands and the merge counts.
- The progress bar gives each Cellpose run half of the detection part (0-42 % and 42-84 %), so it keeps moving during the second run.
- Response-assisted and all-frames detection, the sensitivity check, measurements and Step 4 are unchanged.

## 2026.11.20.7

- **Blind cell selection for % responding.** Step 1 now asks how cells are chosen: **Before the stimulus only (blind)**, the default, **All frames**, or **All frames, with the response in magenta**. The magenta choice stays off by default.
- Blind: you confirm the last frame before the stimulus. The plugin offers the frames that end well before the rise it finds, and asks you to type it when no response is found. The picture shows those frames only, with nothing in the magenta channel. Cellpose gets that picture only, never the response (saved in `automatic_outlines_blind/` with `detection_picture.txt`). Suggestions are graded on shape, size and brightness in that picture, so no suggestion is graded "responds". Hovering shows what a cell is but no trace. Changing the frames after the stimulus changes nothing about blind detection.
- The blind picture for Cellpose is shown from 0 to 2 gray levels with brighter pixels capped, so very bright objects cannot flatten dim cells. On P3 Run2 ATP (real Cellpose), of 260 response-assisted outlines that could be cells, the pre-stimulus picture found 94 at the old 0–8 range and 131 at 0–2. Of the 182 cells accepted there on the magenta picture, it found 59 and 91. A square root on the levels, or dividing by the lamp light, found fewer. Not yet tried on other recordings.
- Step 2's **Baseline…** view is replaced by **Before stimulus (blind)…**. Before, it changed only the picture on screen while Cellpose still got the all-frames picture. It now switches the whole selection to blind. Suggestions made with the response are dropped from view when switching to blind.
- The recorded selection basis is `Pre-stimulus, blind (frames a-b)` only if nothing of the response (magenta, a trace, the analysis page) had been shown before the outlines were final. Outlines changed after that are recorded as `Pre-stimulus picture (frames a-b), not blind: …`. Showing magenta from the blind view asks first.
- **Sensitivity check** (Step 3 tick box, blind only, on by default): after measuring, response-assisted detection runs once for comparison, reusing saved outlines when the picture matches. The report gives how many of its outlines that could be cells the blind detection also found, by grade, and how many coincide with measured cells. It is added to `what_was_done.txt`, with `..._sensitivity_check.csv` and the missed outlines as `..._sensitivity_check_missed_RoiSet.zip`. The comparison outlines are never added to the cells or measured.
- The plugin window no longer floats over other programs when the plugin is started while another program is in front. It stays above Fiji's own windows only while Fiji is in use, and drops back within about two seconds when it is not.
- A copy newer than the latest release (a test build) now says so on the first screen instead of "You are up to date".
- Measurements, background rules, response thresholds and Step 4 are unchanged. The response-assisted and all-frames detection pictures are unchanged (checked byte for byte against the saved P3 Run2 picture).
- Fix a false Step 4 warning introduced in 2026.11.20.6: "the response baseline holds 0 samples … every cell will be unscorable" appeared for any baseline. Baseline samples are now counted as the analysis engine counts them (baseline start ≤ frame time < baseline end), and the warning names the frames.
- Step 4 stimulus timing: **Noise baseline to** and **Score response from** now follow **Stimulus added at** by default and are shown as "= stimulus time". Changing the stimulus time, or pressing **Use frame N**, moves them with it, using the identical stored number. **Set separately** gives a field its own value and **Link to stimulus time** reverses it. A protocol whose values differ opens with those fields separate. The protocol file format is unchanged. The legacy rule works the same way, with the stimulus time as the one editable field.
- Step 4: each window shows the frames it will use, for example "F₀: frames 1–20. Noise baseline: frames 1–20 (20 samples). Stimulus: between frames 20 and 21. Scoring: frames 21–300.", and the noise baseline, stimulus and scored frames of the primary stimulus are drawn under the frame slider.
- Step 4: when two time fields disagree, both are marked and one message sits beside the field to change, stating both values and offering the repair, for example "Set to 20 s". A problem is described once: the plugin's general message is not added to the page's specific one.
- Step 4: a typed time that reads the same as a frame's time is stored as that frame's exact time, so a retyped value no longer differs from an untouched one by a hidden fraction of a second.
- Step 4: the warning that the F₀ frames reach the stimulus offers **Use frames 1–N**. The starting stimulus time is labelled as a suggestion until it is changed. "Resting baseline frames" is now **F₀ frames (set the ΔF/F₀ zero)** and "Response baseline" is **Noise baseline**; the fields are in the order they occur in time.


## 2026.11.20.6

- Step 2 has a **Picture** switch: **All frames**, **Baseline…** and **Response (magenta)**. Magenta can now be turned on after the recording is open; the response is looked for at that moment, exactly as the first screen's option does, and the picture changes in place so outlines and clicks keep working. The first screen's option still works and still starts unticked.
- Fix **Review baseline frames**, which could not show its picture (ImageJ refused to put a flat picture into the colour composite). The baseline view now offers the detected baseline frames when there is a response.
- The selection basis recorded with each measurement follows the picture in use: all frames, baseline frames, or response-assisted. When "look for a response" was ticked but no response was found, it is now recorded as all frames rather than response-assisted.
- The saved `picture_colour.tif` is always the all-frames colour picture, whatever view is on screen, with a magenta channel only when a response was found, as before.
- **‹ Change recording** in the Recording card goes back to the first screen. If the ROI Manager holds outlines, it asks whether to clear or keep them, because the next recording starts with whatever outlines it holds. **‹ Back to cells** and **‹ Back to measuring** move back between steps.
- The plugin window stays above Fiji's own windows only while Fiji is in use, and drops back when another program, such as a browser, is in front.
- Step 4: a stimulus window and a positive-control window that overlap, or a baseline that falls inside the other window's response interval, are now flagged while editing, saying what to change. The analysis engine has always refused these, but the run used to fail with only "exit code 1".
- Step 4: when the analysis engine stops, the page and the window show its own error message instead of only the exit code.
- Step 4: a run in which no cell could be scored for the primary stimulus is reported as a warning naming the commonest reason, not as a green "Complete". The saved results are unchanged.
- Step 4 exports are drawn at least 1600 pixels wide, whatever the window size, and their names include the recording, the cell for trace plots, and the time, for example `Recording_traces_Cell001_20261120-153012.png`. The metadata records the results folder and the time range shown.
- Step 4: a protocol file that is not valid JSON, or whose name was edited in only one place, gets a plain explanation instead of a parser message. After a cancel, the progress stage is cleared. The positive-control tick box now sits above the control's own fields.
- Step 4 quality review: the "minimum response duration is shorter than one frame interval" flag, which the engine puts on every cell, is shown once as a notice for the recording and no longer puts every cell into **Needs quality review**. Each flag in the summary is a filter that lists only its cells. Added plain wording for the baseline-noise flag. The saved tables are unchanged.
- Step 4: the headline tiles follow the response window selected, so the positive control shows its own counts.
- Step 4: cells with long stable IDs are shown as "Cell N" (their measured order) on the picture, plot, heatmap and lists, with the full ID in tooltips. Tables and exports keep the IDs.
- Step 4: **Use frame on screen** under "Stimulus added at" sets the onset from the movie frame shown.
- Step 4: a warning while editing when a response baseline cannot hold the minimum number of baseline samples, which would leave every cell unscorable. Discard and protocol import clear old messages. Field names in number errors use their on-screen labels.
- Step 4 movie playback is faster and steadier. Speeds now run from 0.5× to 50×, with **10×** as the default. Frames are shown on a clock in recording time, so they arrive at even intervals. Above 30 frames per second on screen, every Nth frame is shown and the frame note says so. A long gap between recorded frames waits at most about one second.
- Step 4 frames are loaded up to 4 at a time, about 3 s of playback ahead, and kept for the session (within 90 MB), so scrubbing back over frames already seen makes no new requests. If loading falls behind, the picture holds and the note says "buffering". Frames are fetched at the size they are drawn (full size when zoomed in) and may be cached by the browser. Each frame is keyed to the baseline and display range, so a changed baseline never shows an old image. While playing, the trace cursor glides with the clock. The displayed images and measurements are unchanged.

## 2026.11.20.5

- Validate Step 4 settings while editing, mark invalid fields and disable analysis until errors are corrected. Preserve separate elapsed-time and legacy drafts when switching rules; add Discard pending changes to restore the completed analysis or initial settings.
- Require actual positive-control timing instead of inventing a window beyond the recording. Place the suggested primary onset after the final normalization-baseline frame and show warnings for baseline/onset overlap and incomplete recording coverage.
- Surface saved quality flags, low-reference warnings, coarse temporal resolution and negative corrected-fluorescence traces in an expandable summary and per-cell feedback. Keep response rules and inclusion/exclusion policy unchanged; separate Excluded, Unscorable and Needs quality review filters.
- Give each measurement a separate browser session and preserved input files. Existing pages remain available when another measurement is opened or the cell-selection window is closed, until Close analysis or Fiji exits. Distinguish tabs by cell count and session identity, and give feedback when opening a log or results folder.
- Keep short-window layouts scrollable with readable plots and visible captions. Show core response-call columns first, pin cell names while scrolling, explain empty tables, and label legacy whole-recording metrics separately from stimulus-window metrics.
- Reject invalid frame numbers and display modes at the local frame endpoint. Preserve full frame-rate precision in the settings box and keep completed results and exports tied to their applied settings.

## 2026.11.20.4

- Move Check for updates, Install update, Later, Release notes and the startup-check preference onto Step 1, before choosing or opening a recording. Show update offers and installation progress inline on that screen.
- Keep recording selection responsive during update checks, including offline startup. Pause Next/Cancel only while installing; show retryable errors or the installed version and restart instructions without losing the chosen recording or folder.
- Remove the duplicate update area and delayed automatic check from the cell-selection window. Make Step 1 scrollable on smaller screens, with Next and Cancel always outside the scrolling area.

## 2026.11.20.3

- Bring the portable response-protocol workflow into the local Step 4 browser page while preserving the 2026.11.20.2 movie, interactive traces, heatmap, cell navigation and native Fiji fallback. Browser displays read the shared engine's saved results; they do not recalculate response calls.
- Open the browser in **Recording** view with prominent Play, frame-step, slider and speed controls available before analysis. Playback uses the available frame times; when no valid time base is available, the preview explicitly labels its 20-frames-per-second fallback at 1×. Label the colour overview as static, disable movie playback in that view, and hide result tiles until an analysis has completed.
- Keep stimulus timing and a response-rule summary visible while placing detailed rules, protocol files and diagnostics in expandable sections. Pending edits leave the completed results, their movie timing and exports unchanged until **Apply changes & reanalyze** completes.
- Add versioned protocols shared with SimpleCa²⁺ commit b92d4ce98ddab3fceac352179703c1f48f044c90: independent baseline/onset/scoring windows, one generically named positive control, editable per-window floor/noise/persistence rules, JSON import/export and retained field timing overrides.
- Offer a provisional 0.10/3-SD/1-second preset and preserve the actual historical schema-1 sample rule. Version 2 uses timestamp spans and explicit complete-window eligibility; missing data or inadequate coverage remain unscorable. Review flags do not automatically exclude low-reference or peak-outlier cells.
- Analyze full-precision unlifted corrected fluorescence. Preserve stable ROI identities, accepted/edge-excluded registries and the selection basis; response-assisted selection remains an optional review aid. Older rounded and intensity-lifted exports remain labeled compatibility outputs.
- Add exact timestamp CSVs, input/ROI/protocol snapshots, applied-versus-draft state, transformed local thresholds, qualifying-run evidence, window peak/AUC/latency measurements and plot provenance. Changes to timestamp or diagnostic file contents mark the draft as changed even at the same path.
- Add separately labeled sensitivity and negative-control diagnostics, condition/biological-unit metadata and summaries that distinguish biological units from cell-pooled percentages. Zero eligible cells produce an unavailable percentage.
- Keep legacy mode's inactive fields explicit, preserve imported stimulus/control criteria during timing edits, and reject inconsistent protocol metadata. The managed runtime verifies jsonlite support. Extend independent scientific fixtures, protocol/precision checks, local browser/server checks and native GUI regressions; bundle both the pinned R sources and browser assets in the single plugin file.

## 2026.11.20.2

- Step 4 now opens as a page in your browser, built with the same web technology and styling as the SimpleCa²⁺ app: results first, interactive traces (hover a trace to name it, click to select it, scroll to zoom the time axis), the recording played beside the traces, an ordered heatmap and sortable tables. The page is served by the plugin from your own computer; nothing is sent anywhere and no internet is needed.
- The page drives the same analysis engine with the same settings and checks as before, and reads the tables that engine saves. Nothing is recalculated in the browser. **Use a Fiji window instead** opens the previous Step 4 window, which is also used automatically if a browser page cannot be opened.
- The page is reachable only from this computer and only with a long random key in its address.
- The page's heatmap is drawn the way the SimpleCa²⁺ app draws its own: every cell in one picture, axes with ticks, a colour bar with values, the app's blue–white–red for signed values or one of its palettes (plasma, viridis, magma, inferno, cividis) from zero up, and sorting by response call, time to peak, peak amplitude or original order. The colour range covers 99 % of the values by default, or the largest value as in the app.
- Step 1: choosing **Choose a file from disk...** now opens the file chooser straight away, and a **Choose file…** button with the chosen file's name sits under the list. Before, the chooser only appeared after **Next**, so picking that entry looked as if nothing happened.
- **Next** still asks for the file if none was chosen, and stays on Step 1 if that is cancelled instead of closing the plugin.
- Step 4, movie beside the traces: a time slider, a **Play** button and a line on the plot move together, and the picture can show the recording at that moment. **Movie: change from rest** shows each frame minus the average of the baseline frames (lightly smoothed), so a responding cell lights up as its trace rises; **Movie: as recorded** shows the frame itself. **Zoom to cell** and **All outlines** help follow one cell. The recording is only read. Display only: nothing measured or saved changes.
- Step 4, every trace: all cells are drawn behind the mean, coloured by response call (the SimpleCa²⁺ colours), with the selected cell on top. The y-axis covers 99 % of all values so one extreme cell cannot flatten the rest.
- Step 4, cell groups: list and draw only responders, non-responders, or excluded and unscorable cells, from the **Show** list or by clicking a headline tile. A line under the list says what was decided for the selected cell, with the largest rise and the threshold from the saved table. Left and right arrow keys step through the cells.
- Step 4 heatmap: responders first with the strongest on top (or as measured), a colour scale that covers 99 % of the values instead of the single largest, taller rows with readable names, a strip for each cell's call, and markers for the baseline and the stimulus.
- Step 4 layout: narrower settings column that can be hidden, picture and traces side by side, larger default window. **Also score a separate KCl window** is now **Include a KCl control**, with an explanation.
- Main window: a long recording name or folder name is shortened in the middle (full text in the tooltip), so it no longer pushes the window contents off the right edge.

## 2026.11.20.1

- Redesign every plugin window to match SimpleCa²⁺: the computer's own interface font instead of a fixed one, the SimpleCa²⁺ colours, white cards, and one clearly marked main button per step. Measurements, thresholds and saved tables are unchanged.
- Main window: one card per step. Only one of Steps 2 to 4 is open at a time and finished steps collapse to a one-line summary, so the window fits the screen. The cell count wraps instead of running off the edge, long instructions moved to tooltips, the detection progress area is smaller, and version/update controls sit behind an **Updates** button until an update is found.
- Step 1 is a single styled dialog with a folder **Browse** button.
- Step 4: results first, with headline tiles for percent responding, registered, excluded and unscorable cells. Settings are grouped into **Timing** and **Response rule**, show rounded numbers, and hide the KCl boxes until that option is ticked. A box that still shows its rounded suggestion sends the exact value to the analysis, so results do not change.
- Step 4 plots: axis ticks, units and names, labelled baseline and stimulus bands, a colour key on the heatmap, and previous/next cell buttons.
- The cell review picture is now inside the Step 4 window, in full colour, with outlines styled by response call. Scroll to zoom, drag to move, double-click to see the whole picture.
- Windows no longer cover each other: the colour picture is fitted into the space left of the main window, and the Step 4 window opens in that space too.

## 2026.11.19.8

- Add Step 4, SimpleCa²⁺ analysis inside Fiji: confirmed timing/baseline controls, native traces and heatmap, cell-linked image review, metrics, stimulus/KCl responder calls and CSV/PNG export.
- Reuse the unmodified metric and responder functions pinned from SimpleCa²⁺ commit a472f2460075343d5b97c8bac32ddf3da0c77386. Preserve invalid-baseline cells as excluded records and keep unscorable cells explicit; use the unlifted corrected fluorescence only.
- Bundle the complete analysis source and license inside the plugin, so existing Install update buttons install the integration. Provide one-time, per-user R setup with a verified micromamba launcher, progress/logs, cancellation and retry.
- Save the exact input snapshot, settings, source hashes, R/plugin versions and results together. Preserve completed analyses; partial runs are never displayed as completed results.
- Add real R parity/runtime checks on Windows, macOS and Linux, plus process and native-window regressions.

## 2026.11.19.7

- Make accepted cells clearly distinct: white outlines with a thick black outer border, and a thicker highlight for selected cells. Edge exclusions remain dashed.
- Fix clicking accepted cells so Remove selected cell targets the clicked cell, even with suggestions hidden or no automatic suggestions. Support multiple ROI Manager selections.
- Add Click cells to remove for reviewing many accepted cells directly on the picture. This mode never accepts suggestions; drawing tools turn it off.
- Replace Undo last cell with Undo last change: restore removed cells or undo additions, including Accept all, for up to 20 plugin edits. Preserve ROI names, geometry and order; outside ROI Manager edits invalidate undo instead of being overwritten.
- Add native GUI regressions for actual image clicks, removal, undo, non-cell ROIs, overlap targeting and contrasted outlines.

## 2026.11.19.6

- Fix the visible-window pack exception caused by Groovy treating Rectangle coordinates as doubles. Use explicit integer coordinates and statically compile the window class.
- Start detection independently of cosmetic window updates. Startup failures stop the progress timer and restore the Suggest button; result-display errors cannot leave detection stuck.

- Add Install update to the update prompt and main window. It downloads the plugin inside Fiji, verifies the published SHA-256 checksum and version, backs up the previous file, and replaces it atomically.
- Show installation progress, retryable errors and a restart instruction. Keep the running analysis session open; the new version is used after restarting Fiji.
- Prevent concurrent installers and downgrades. Failed downloads, verification, backups or replacements leave the current plugin in place.
- Provide Install_Calcium_Plugin.groovy for first installation or upgrading versions without the installer: open it in Fiji's Script Editor and press Run. No manual file copying is needed.

## 2026.11.19.5

- Keep controls at their full layout size when the plugin window is made smaller. Scroll vertically with the mouse wheel or scrollbar, and horizontally when the window is too narrow.
- Preserve the window size chosen by the user when detection progress, results, or update messages change.
- Fit the initial window within the screen's width and height and keep a usable minimum window size.
- Add native Windows and Linux GUI checks for resizing, wheel scrolling, unclipped controls, and access to the bottom and right-hand buttons.

## 2026.11.19.4

- Keep a plain-language explanation visible during detection, with total and current-stage time, image dimensions, CPU/GPU, attempt number, and detection settings.
- Show whether saved outlines are being reused or a new detection is running. Add a scrollable Progress details window with recent events and Cellpose messages.
- Keep the main window within the screen height and make its contents scrollable so all controls stay reachable.
- Preserve stage explanations during quiet periods and report whether the Cellpose process is still running. Progress remains tied to reported stages, with no invented time estimate.
- Check for new releases in the background when the plugin window opens, with a non-blocking Download / Later prompt and an option to disable startup checks. Manual checks remain available; installation is still manual.
- Publish a tested release automatically when a higher version is merged to main. Existing releases are never overwritten.
- Install this version once to receive startup prompts for future releases; earlier versions can use Check for updates to obtain it.

## 2026.11.19.3

- Use the latest user-supplied plugin as the base, preserving edge-cell exclusion, matched-image outline reuse, faster CPU detection, and automatic Torch backend selection during setup.
- Integrate percentage stage progress, detection logs, cancellation, and GitHub update checking into that version.
- Prevent failed or interrupted detection from being reused as a completed cached result.
- Keep detection settings and edge-cell metadata in the saved records.

## 2026.11.19.2

- First public GitHub release, with a permanent latest-download link.
- Add an explicit **Check for updates** button and download-page link inside the plugin.
- Add release packaging, checksums, automated regression checks, and a maintainer publishing workflow.
- Include the Step 2 progress changes from 2026-11-19.1 below.

## 2026-11-19.1

- Show percentage stage progress, current stage, elapsed time, and CPU/GPU information when supplied by Cellpose.
- Enable verbose, unbuffered Cellpose logs and add **View detection log**.
- Add cancellation and feedback when a run has no new progress messages or exceeds ten minutes.
- Remove the fixed “about 30 seconds” estimate.
- Show 100% only after the suggestions are read, graded, and saved.
- Prevent a failed rerun from reusing an old outlines ZIP; report failures from the current attempt.

## 2026-11-19

- Original packaged workflow for choosing cells, local background correction, and fluorescence export.

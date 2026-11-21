// Calcium imaging: cells step by step                                                    version 2026.11.20.19
// One plugin for the whole job on one recording:
//   1. choose the recording (already open in Fiji, or a file) and a folder for the results
//   2. the plugin makes the cells visible: all-frames picture, response (found automatically, if there is one), and a
//      colour picture: green = cells, yellow = very bright objects (precipitate), magenta = brightened after the stimulus
//   3. you choose the cells on the colour picture: either draw around a cell and press t, or press "Suggest cells
//      automatically" (Cellpose, with stage progress) and click inside each suggested outline that is a real cell.
//      While you work, the plugin finds the background next to every cell and shows it in blue.
//   4. press "Measure and save": gray values of every cell and of its background are measured on the ORIGINAL movie,
//      and three tables are saved: raw, background, corrected (raw - background). The corrected table is for SimpleCa2+.
//      A background-corrected copy of the movie can be opened too (and saved, if wanted).
// The recording itself is never changed. Everything the plugin did is written to a log file next to the tables.
//
// Install: in Fiji choose File > Show Folder > ImageJ. In the folder that opens go into scripts > Plugins, make a
//          folder called "Calcium Imaging", copy this file into it, and restart Fiji.
// Run:     Plugins > Calcium Imaging > Calcium Cells Step By Step
// Try without installing: drag this file onto the Fiji bar and press Run in the editor that opens.
// Needs:   Fiji (it includes Bio-Formats, which opens .zvi/.czi files). Cellpose is optional: only for the
//          "Suggest cells automatically" button. See README.txt.
// The settings below are in pixels and were chosen on 1388 x 1040 recordings at 0.645 micrometres per pixel.

import groovy.transform.CompileStatic
import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.CompositeImage
import ij.Prefs
import ij.WindowManager
import ij.gui.GenericDialog
import ij.gui.Overlay
import ij.gui.Roi
import ij.gui.ShapeRoi
import ij.gui.ImageRoi
import ij.io.OpenDialog
import ij.io.RoiEncoder
import ij.plugin.ZProjector
import ij.plugin.RoiEnlarger
import ij.plugin.filter.RankFilters
import ij.plugin.frame.RoiManager
import ij.process.ByteProcessor
import ij.process.ColorProcessor
import ij.process.FloatProcessor
import ij.process.ImageProcessor
import ij.process.LUT
import java.awt.Color
import java.awt.Rectangle
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

@CompileStatic
class CalciumSteps {
    // ---- settings (in pixels; chosen on 1388 x 1040 recordings at 0.645 micrometres per pixel) ----
    static double BALL_RADIUS = 150      // Process > Subtract Background, rolling ball radius (picture)
    static double BLUR_PICTURE = 1.5     // Process > Filters > Gaussian Blur (picture)
    static double SHOW_MAX = 8           // all-frames picture is shown 0 to this many gray levels
    static double BRIGHT_MIN = 10, BRIGHT_MAX = 30   // red channel: very bright objects (precipitate) turn yellow
    static double BLUR_RESPONSE = 2, RESP_MIN = 1.00, RESP_MAX = 1.10
    static int GAP = 6                   // background band starts this many pixels outside the ROI
    static int REACH = 40                // ... and ends this many pixels outside the ROI (about 26 micrometres)
    static int EDGE_MARGIN = 40          // cells closer than this to the picture's border can be left out (tick box): their background band is cut off there
    static double KEEP_DARKEST = 0.25    // the darkest quarter of the allowed band pixels is the background
    static double BRIGHT_ABOVE = 8       // pixels brighter than this in the all-frames picture are never background
    static int BRIGHT_GROW = 4           // ... widened by this many pixels
    static int ROI_MARGIN = 4            // other ROIs, widened by this many pixels, are never background
    static int MIN_PIXELS = 150          // fewer usable pixels: the band is doubled once; still fewer: no background
    static int SUGGEST_MIN_AREA = 150, SUGGEST_MAX_AREA = 6000   // automatic outlines outside this size (pixels) are not shown
    static double SUGGEST_BRIGHT = 12    // automatic outlines brighter than this on average are shown in red: probably precipitate
    // Blind selection: cells are found and accepted on frames before the response only, so whether a cell is chosen cannot depend on its response.
    // Detection pictures for Cellpose: the picture of the frames before the response shown twice, from 0 to BLIND_CAP and from 0 to BLIND_CAP_HIGH gray
    // levels, brighter pixels capped. The low cap keeps precipitate and rounded cells from flattening the dim cells; the high cap
    // keeps the brighter cells whole. All outlines from the first picture are kept; an outline from the second is added when it
    // shares less than BLIND_MERGE_OVERLAP of the smaller outline's area with every outline from the first.
    // Tested on four ATP recordings of 18 Nov 2026 (real Cellpose, 20 Nov): of the response-assisted outlines that could be cells,
    // the two pictures together found 68 / 95 / 55 / 62 % (C2 Run1 / C3 Run2 / P2 Run1 / P3 Run2), against 46 / 94 / 34 / 50 % with
    // cap 2 alone and 57 / 86 / 45 / 36 % with cap 8 alone; no single cap was best on all four, and a cap taken from each picture's
    // own 90th percentile was never better than the best fixed cap. Cost: 25-55 % more suggestions to review.
    static double BLIND_CAP = 2
    static double BLIND_CAP_HIGH = 8
    static double BLIND_MERGE_OVERLAP = 0.10
    // -------------------------------------------------------------------

    ImagePlus movie
    int n, w, h
    float[] rest                 // average of all frames (raw)
    float[] flat                 // the all-frames picture (lamp light removed): used to rank pixels by how empty they are
    float[] respDiff             // mean of the response frames minus mean of the frames before (raw gray levels); null without a response
    int lastClear = 0, lastDim = 0, lastLifted = 0; double lastMeanRise = 0
    int[] lastBackgroundPixels = new int[0]
    static double KEEP_LEVEL = 1.0       // "every cell kept" table: a cell whose lowest 20-frame average (before the response) is below this is lifted to it
    ImagePlus allPic, respPic, colourPic
    /** The all-frames picture with the pixels used as the cells' background in blue (no outlines): made at measuring, shown by the Step 4 page. */
    ImagePlus backgroundPic
    boolean[] bright             // widened bright objects
    boolean found = false; int onset, p0, p1, b0, b1; double share; String why = ""
    /** Shift of every frame (CalciumMotion.shifts), or null when the movement could not be checked. */
    double[][] motionShifts
    /** What the shifts mean now: counted from the frames the cells are chosen on while choosing blind, else from where the field
     *  sits for most of the recording. Null when the movement could not be checked. */
    Map getMotion() { motionShifts == null ? null : CalciumMotion.summary(motionShifts, blind ? blindFirst : 0, blind ? blindLast : 0) }
    /** Says the movement in the log once per reference, so the record matches what the Recording card shows. */
    String motionSaid = null
    void sayMotion() {
        String now = CalciumMotion.sentence(getMotion())
        if (now && now != motionSaid) say("Movement: " + now + " Outlines stay where they were drawn; nothing was corrected.")
        motionSaid = now
    }
    boolean looked = false             // whether the response has been looked for, at the start or later from Step 2
    List<String> log = []
    Roi[] alsoAvoid = new Roi[0]         // outlines that are not measured (cells at the edge that are left out) but must never serve as background
    String edgeNote = ""                 // what was done about cells at the picture's edge, for the record of steps
    boolean blind = false                // cells are being chosen blind: picture, suggestions and grades come from frames blindFirst..blindLast only
    int blindFirst = 0, blindLast = 0
    float[] blindPic                     // frames blindFirst..blindLast, processed like the all-frames picture (shown, and used to grade suggestions)

    void say(String s) { log.add(s) }

    /** Steps 1-2: pictures. The movie is only read. */
    void prepare(ImagePlus imp, boolean lookForResponse) {
        movie = imp; n = imp.getStackSize(); w = imp.getWidth(); h = imp.getHeight()
        say("Calcium imaging plugin, version 2026.11.20.19")
        say("Recording: " + imp.getTitle() + " (" + n + " frames, " + w + " x " + h + ", " + imp.getBitDepth() + "-bit)")
        // Did the field move? Nothing is corrected; moved frames are named so that they can be kept out of baselines.
        try { motionShifts = CalciumMotion.shifts(imp.getStack(), n) } catch (Throwable ignored) { motionShifts = null }
        ImagePlus sum = ZProjector.run(imp, "sum")                                   // Image > Stacks > Z Project, Sum Slices
        FloatProcessor avg = (FloatProcessor) sum.getProcessor().convertToFloat()
        avg.multiply(1.0d / n)                                                       // Process > Math > Divide by the number of frames
        rest = (float[]) ((float[]) avg.getPixels()).clone()
        allPic = new ImagePlus("all_frames_" + imp.getTitle(), avg)
        allPic.setCalibration(imp.getCalibration().copy())
        IJ.run(allPic, "Subtract Background...", "rolling=" + (int) BALL_RADIUS)     // Process > Subtract Background
        IJ.run(allPic, "Gaussian Blur...", "sigma=" + BLUR_PICTURE)                  // Process > Filters > Gaussian Blur
        allPic.setDisplayRange(0, SHOW_MAX)                                          // Image > Adjust > Brightness/Contrast > Set
        say("All-frames picture: Z Project Sum Slices 1-" + n + ", Divide " + n + ", Subtract Background rolling=" + (int) BALL_RADIUS + ", Gaussian Blur sigma=" + BLUR_PICTURE + ", shown 0-" + SHOW_MAX)
        float[] a = (float[]) allPic.getProcessor().getPixels()
        flat = (float[]) a.clone()
        ByteProcessor bp = new ByteProcessor(w, h)
        byte[] bpx = (byte[]) bp.getPixels()
        for (int i = 0; i < a.length; i++) if (a[i] > BRIGHT_ABOVE) bpx[i] = (byte) 255
        new RankFilters().rank(bp, (double) BRIGHT_GROW, RankFilters.MAX)            // Process > Filters > Maximum
        bright = new boolean[w * h]
        for (int i = 0; i < bpx.length; i++) bright[i] = bpx[i] != 0
        if (lookForResponse) findResponse()
        colourPic = colourComposite(true)
        say("Colour picture: red = all-frames picture shown " + BRIGHT_MIN + "-" + BRIGHT_MAX + " (precipitate turns yellow), green = all-frames picture 0-" + SHOW_MAX + (found ? ", magenta = response" : ""))
    }

    /** Looks for the response and, when there is one, makes the response picture. The movie is only read. */
    void findResponse() {
        looked = true
        findResponseFrames()
        if (found) {
            ImagePlus pk = ZProjector.run(movie, "sum", p0, p1), pre = ZProjector.run(movie, "sum", b0, b1)
            float[] x = (float[]) pk.getProcessor().convertToFloat().getPixels(), y = (float[]) pre.getProcessor().convertToFloat().getPixels()
            float f = (float) ((p1 - p0 + 1) / (double) (b1 - b0 + 1))
            float[] r = new float[x.length]
            for (int i = 0; i < r.length; i++) r[i] = y[i] > 0 ? (float) (x[i] / y[i] / f) : 1f      // Process > Image Calculator, Divide
            respDiff = new float[x.length]
            for (int i = 0; i < r.length; i++) respDiff[i] = (float) (x[i] / (p1 - p0 + 1) - y[i] / (b1 - b0 + 1))
            respPic = new ImagePlus("response_" + movie.getTitle(), new FloatProcessor(w, h, r))
            IJ.run(respPic, "Gaussian Blur...", "sigma=" + BLUR_RESPONSE)
            respPic.setDisplayRange(RESP_MIN, RESP_MAX)
            say("Response: brightness rises near frame " + onset + " in " + Math.round(100 * share) + " % of the field. Response picture = frames " + p0 + "-" + p1 + " divided by frames " + b0 + "-" + b1 + ", Gaussian Blur sigma=" + BLUR_RESPONSE + ", shown " + RESP_MIN + "-" + RESP_MAX)
        } else say("Response: none found (" + why + "). The colour picture has no magenta.")
    }

    /** The colour picture as a Fiji composite: red and green from the all-frames picture, magenta from the response.
     *  On screen there is always a third channel (blank, and switched off, until a response is shown), so the picture can
     *  change view later without replacing its window. The saved copy has a magenta channel only when there is a response. */
    CompositeImage colourComposite(boolean forDisplay) {
        ImageStack st = new ImageStack(w, h)                                         // Image > Color > Merge Channels
        st.addSlice("very bright", allPic.getProcessor().duplicate())
        st.addSlice("all frames", allPic.getProcessor().duplicate())
        boolean third = found || forDisplay
        if (third) st.addSlice("response", found ? respPic.getProcessor().duplicate() : new FloatProcessor(w, h))
        ImagePlus tmp = new ImagePlus("colour_" + movie.getTitle(), st)
        tmp.setDimensions(st.getSize(), 1, 1)
        CompositeImage ci = new CompositeImage(tmp, IJ.COMPOSITE)
        ci.setChannelLut(LUT.createLutFromColor(Color.red), 1); ci.setC(1); ci.setDisplayRange(BRIGHT_MIN, BRIGHT_MAX)
        ci.setChannelLut(LUT.createLutFromColor(Color.green), 2); ci.setC(2); ci.setDisplayRange(0, SHOW_MAX)
        if (third) { ci.setChannelLut(LUT.createLutFromColor(Color.magenta), 3); ci.setC(3); ci.setDisplayRange(RESP_MIN, RESP_MAX) }
        if (third && !found) ci.setActiveChannels("110")
        ci.setC(1); ci.setCalibration(movie.getCalibration().copy())
        return ci
    }
    /** The all-frames colour picture as one ordinary colour photo, so every program shows it in colour (the window's own
     *  picture is a three-channel stack, which only ImageJ shows coloured). Green = all-frames picture 0..SHOW_MAX; red is
     *  added where it is BRIGHT_MIN..BRIGHT_MAX, so very bright objects turn yellow; magenta = the response picture
     *  RESP_MIN..RESP_MAX when a response was found. The same display ranges as the window. */
    ImagePlus colourPhoto() {
        float[] g = (float[]) allPic.getProcessor().getPixels()
        float[] m = found ? (float[]) respPic.getProcessor().getPixels() : null
        int[] rgb = new int[w * h]
        for (int i = 0; i < rgb.length; i++) {
            int gv = (int) Math.max(0, Math.min(255, Math.round(g[i] / SHOW_MAX * 255)))
            int rv = (int) Math.max(0, Math.min(255, Math.round((g[i] - BRIGHT_MIN) / (BRIGHT_MAX - BRIGHT_MIN) * 255)))
            int mv = m == null ? 0 : (int) Math.max(0, Math.min(255, Math.round((m[i] - RESP_MIN) / (RESP_MAX - RESP_MIN) * 255)))
            rgb[i] = (Math.min(255, rv + mv) << 16) | (gv << 8) | mv
        }
        ImagePlus photo = new ImagePlus("colour_" + movie.getTitle(), new ColorProcessor(w, h, rgb))
        photo.setCalibration(movie.getCalibration().copy())
        return photo
    }
    /** The colour photo with the measured cells outlined in white and, below it, a key to the colours and a scale bar. */
    ImagePlus colourPhotoWithCells(Roi[] cells) {
        int keyHeight = 96
        ColorProcessor sheet = new ColorProcessor(w, h + keyHeight); sheet.setColor(Color.black); sheet.fill()
        sheet.insert(colourPhoto().getProcessor(), 0, 0)
        ImagePlus outlined = new ImagePlus("outlined", sheet); Overlay ov = new Overlay()
        for (Roi c : cells) { Roi r = (Roi) c.clone(); r.setPosition(0); r.setStrokeColor(Color.white); r.setStrokeWidth(2); r.setFillColor(null); ov.add(r) }
        outlined.setOverlay(ov)
        ImageProcessor ip = outlined.flatten().getProcessor()
        ip.setAntialiasedText(true); ip.setFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 20))
        List<Color> colours = [Color.green, Color.yellow]; List<String> texts = ["cells: average of all " + n + " frames", "very bright objects (precipitate)"]
        if (found) { colours << Color.magenta; texts << ("brightened after the stimulus (near frame " + onset + ")") }
        int x = 14, y = h + 12
        for (int i = 0; i < texts.size(); i++) {
            int width = 26 + ip.getStringWidth(texts[i]) + 26
            if (x + width > w - 8 && x > 14) { x = 14; y += 28 }
            ip.setColor(colours[i]); ip.fillRect(x, y + 2, 18, 18); ip.setColor(Color.white); ip.drawString(texts[i], x + 26, y + 22); x += width
        }
        x = 14; y += 28
        ip.setColor(Color.white); ip.setLineWidth(2); ip.drawRect(x, y + 2, 18, 18); ip.drawString("outlines: the " + cells.length + " cells measured", x + 26, y + 22)
        double pixel = movie.getCalibration().pixelWidth; String unit = movie.getCalibration().getUnit()
        if (pixel > 0 && unit != null && (unit.startsWith("micro") || unit == "um" || unit == "µm" || unit == "\u00B5m")) {   // a 100 micrometre bar at the right
            int bar = (int) Math.round(100d / pixel)
            if (bar < w - 40) { ip.fillRect(w - 14 - bar, y + 8, bar, 6); String label = "100 \u00B5m"; ip.drawString(label, w - 14 - bar + (int) ((bar - ip.getStringWidth(label)) / 2), y + 40) }
        }
        return new ImagePlus("colour picture with the cells", ip)
    }
    /** Shows the response in magenta on the colour picture, looking for it first if needed. False when there is none. */
    boolean showResponse() {
        if (!looked) findResponse()
        if (!found || colourPic.getNChannels() < 3) return false
        float[] target = (float[]) colourPic.getStack().getPixels(3)
        System.arraycopy((float[]) respPic.getProcessor().getPixels(), 0, target, 0, target.length)
        ((CompositeImage) colourPic).setActiveChannels("111"); colourPic.updateAndDraw()
        return true
    }
    void hideResponse() { if (colourPic.getNChannels() >= 3) { ((CompositeImage) colourPic).setActiveChannels("110"); colourPic.updateAndDraw() } }
    /** A picture from frames first to last only, made like the all-frames picture: Z Project Sum Slices, Divide by the number
     *  of frames, Subtract Background, Gaussian Blur. */
    float[] framesPicture(int first, int last) {
        ImagePlus sum = ZProjector.run(movie, "sum", first, last)
        FloatProcessor avg = (FloatProcessor) sum.getProcessor().convertToFloat()
        avg.multiply(1.0d / (last - first + 1))
        ImagePlus picture = new ImagePlus("before_stimulus_" + movie.getTitle(), avg)
        IJ.run(picture, "Subtract Background...", "rolling=" + (int) BALL_RADIUS)
        IJ.run(picture, "Gaussian Blur...", "sigma=" + BLUR_PICTURE)
        return (float[]) picture.getProcessor().getPixels()
    }
    /** Blind selection uses the frames before the RESPONSE. That is not the same as before the stimulus: a drug can reach the
     *  chamber long before any cell responds, and more frames give a clearer picture. Measured 20 Nov 2026: the two blind
     *  pictures outline 70 % of one recording's audited cells on frames 1-20 and 79 % on frames 1-50, and the cells gained
     *  respond as often as the ones found before. The window the plugin offers, for the user to confirm or change:
     *   first frame: 1, or the frame after a movement of the field at the start (CalciumMotion.MOVED pixels or more);
     *   last frame, response found: BLIND_MARGIN frames before the rise. The rise the plugin finds is the field's (about the
     *     median cell); the first cells start their main rise 9-15 frames before it (P3 Run2 12, C2 Run1 9, P2 Run1 14-15;
     *     measured 20 Nov 2026 on the audited traces). A margin of 15 left 1-6 frames, 20 leaves 5-11, 25 leaves 10-16.
     *     A rise too early for that keeps 20 frames when at least 9 frames then remain before the rise, and says so;
     *   last frame, no response found (or a "rise" that is the field moving): BLIND_CEILING, or the last frame of a shorter recording.
     *  Only the time of the rise is used, never which cells rose. A few cells can still respond briefly inside the window
     *  (1 to 4 audited cells per recording at frames 22-60 on the 18 Nov recordings): it ends before the main response. */
    static final int BLIND_MARGIN = 25, BLIND_CEILING = 60
    /** Why the offered window is what it is, in a sentence for the dialog. */
    String blindWindowWhy = ""
    /** True when the "rise" the plugin found lies in or right after frames in which the field moved: then it is the movement. */
    boolean riseIsMovement() {
        if (!found || motionShifts == null) return false
        for (Object r : (List) CalciumMotion.summary(motionShifts).get('runs')) { List run = (List) r; if (onset >= ((Number) run.get(0)).intValue() && onset <= ((Number) run.get(1)).intValue() + 5) return true }
        return false
    }
    int[] proposedBlindWindow() {
        int first = 1; String moved = ""
        if (motionShifts != null) {
            List runs = (List) CalciumMotion.summary(motionShifts).get('runs')
            if (!runs.isEmpty() && ((Number) ((List) runs.get(0)).get(0)).intValue() == 1) {
                int end = ((Number) ((List) runs.get(0)).get(1)).intValue()
                if (end < n - 10) first = end + 1
            }
        }
        boolean rise = found && !riseIsMovement()
        int last
        int full = onset - BLIND_MARGIN, least = Math.max(first, Math.min(b1, onset - 1))
        if (rise && full >= Math.min(least, first + 19)) {
            last = full
            blindWindowWhy = "The cells start to rise near frame " + onset + ". The frames offered end " + BLIND_MARGIN + " frames before that: the first cells rise 10 to 15 frames earlier."
        } else if (rise) {
            last = Math.min(least, first + 19)
            blindWindowWhy = "The cells start to rise near frame " + onset + ", which leaves few frames before the response. Check these frames."
        } else {
            last = Math.max(Math.min(n, BLIND_CEILING), Math.min(n, first + 19))
            blindWindowWhy = (found ? "The only rise found is near frame " + onset + ", where the field moved, so no response was found." : "No response was found.") + " The frames offered go up to frame " + last + ": more frames give a clearer picture."
        }
        // The field may still be settling after the movement: frames at the start of the window that sit away from the rest of it are left out too.
        if (motionShifts != null) for (int round = 0; round < 3; round++) {
            Map inside = CalciumMotion.summary(motionShifts, first, last); int until = 0
            for (String kind : ['runs', 'slightRuns']) for (Object r : (List) inside.get(kind)) { List run = (List) r; int a = ((Number) run.get(0)).intValue(), b = ((Number) run.get(1)).intValue(); if (a <= first && b >= first && b < last - 9) until = Math.max(until, b) }
            if (until == 0) break
            first = until + 1
        }
        if (first > 1) moved = " Frames 1-" + (first - 1) + " are left out: the field moved in them."
        blindWindowWhy += moved
        return [first, last] as int[]
    }
    /** The last frame the plugin offers for blind selection (see proposedBlindWindow). */
    int estimatedLastBeforeStimulus() { return proposedBlindWindow()[1] }
    /** Blind selection from here on: the picture shows frames first..last only, the response channel is emptied, and new
     *  suggestions are found and graded on that picture without the response. Measurements are not affected. */
    void startBlind(int first, int last) { startBlind(first, last, framesPicture(first, last)) }
    /** The same with the picture already made (framesPicture can run off the window's thread; this part changes the window). */
    void startBlind(int first, int last, float[] picture) {
        blindPic = picture; blindFirst = first; blindLast = last; blind = true
        if (colourPic != null) {
            if (colourPic.getNChannels() >= 3) Arrays.fill((float[]) colourPic.getStack().getPixels(3), 0f)   // nothing of the response stays in the window
            showPicture(blindPic); hideResponse()
        }
        log.removeAll { String s -> s.startsWith("Blind selection:") }
        say("Blind selection: cells are chosen on frames " + first + "-" + last + " only (before the response" + (found && !riseIsMovement() ? "; the rise starts near frame " + onset : "; no response was found") + "). Picture: Z Project Sum Slices " + first + "-" + last + ", Divide " + (last - first + 1) + ", Subtract Background rolling=" + (int) BALL_RADIUS + ", Gaussian Blur sigma=" + BLUR_PICTURE + "; the colour picture shows it in red and green, with nothing in magenta. Cellpose gets this picture only (" + blindRecipe() + "), never the response; suggestions are graded on shape, size and brightness in it; traces are not shown while choosing.")
    }
    void stopBlind() { blind = false }
    String blindBasis() { return "Before the response, blind (frames " + blindFirst + "-" + blindLast + ")" }
    String blindRecipe() {
        return "shown to it twice, from 0 to " + levels(BLIND_CAP) + " and from 0 to " + levels(BLIND_CAP_HIGH) + " gray levels with brighter pixels capped; all outlines from the first picture are kept, and an outline from the second is added when it shares less than " + (int) Math.round(BLIND_MERGE_OVERLAP * 100) + " % of the smaller outline's area with every outline from the first"
    }
    static String levels(double v) { return v == Math.rint(v) ? String.valueOf((long) v) : String.valueOf(v) }
    /** Shows a picture (the all-frames picture, or a baseline picture) in the red and green channels of the colour picture. */
    void showPicture(float[] picture) {
        for (int c = 1; c <= 2; c++) { float[] target = (float[]) colourPic.getStack().getPixels(c); System.arraycopy(picture, 0, target, 0, target.length) }
        colourPic.updateAndDraw()
    }
    void showAllFrames() { showPicture((float[]) allPic.getProcessor().getPixels()) }

    /** Finds the frames of the response without outlines: 24 x 24 pixel squares, largest step up per square, at least 5 % of the squares stepping within the same 15 frames. */
    void findResponseFrames() {
        found = false; why = ""
        if (n < 20) { why = "fewer than 20 frames"; return }
        int T = 24, tw = (int) (w / T), th = (int) (h / T), m = tw * th
        float[][] X = new float[n][m]
        ImageStack st = movie.getStack()
        for (int f = 0; f < n; f++) {
            ImageProcessor ip = st.getProcessor(f + 1); float[] row = X[f]
            for (int ty = 0; ty < th; ty++) for (int tx = 0; tx < tw; tx++) {
                double s = 0
                for (int y = ty * T; y < ty * T + T; y++) { int o = y * w + tx * T; for (int x = 0; x < T; x++) s += ip.getf(o + x) }
                row[ty * tw + tx] = (float) (s / (T * T))
            }
        }
        int[] at = new int[m]; boolean[] resp = new boolean[m]; int nResp = 0
        double[] cs = new double[n + 1]; double[] d = new double[n - 1]
        for (int i = 0; i < m; i++) {
            cs[0] = 0
            for (int f = 0; f < n; f++) { cs[f + 1] = cs[f] + X[f][i]; if (f > 0) d[f - 1] = Math.abs(X[f][i] - X[f - 1][i]) }
            double[] ds = d.clone(); Arrays.sort(ds)
            double med = (ds.length % 2 == 1) ? ds[(int) (ds.length / 2)] : 0.5d * (ds[(int) (ds.length / 2) - 1] + ds[(int) (ds.length / 2)])
            double noise = 1.4826d * med / Math.sqrt(2.0d), best = -1e9; int bestAt = 0
            for (int f = 3; f < n - 12; f++) {
                int lo = Math.max(0, f - 15), hi = Math.max(1, f - 2)
                double s = (cs[f + 12] - cs[f + 2]) / 10.0d - (cs[hi] - cs[lo]) / (hi - lo)
                if (s > best) { best = s; bestAt = f }
            }
            at[i] = bestAt + 1
            if (best > 5 * noise / Math.sqrt(3.0d) && best > 0.3d) { resp[i] = true; nResp++ }
        }
        int need = Math.max(5, (int) (0.02d * m))
        if (nResp < need) { why = "only " + nResp + " of " + m + " squares step up"; return }
        int[] hist = new int[n + 16]
        for (int i = 0; i < m; i++) if (resp[i]) hist[at[i]]++
        int bestCount = -1, f0 = 1
        for (int f = 1; f < n; f++) { int c = 0; for (int k = f; k < f + 15; k++) c += hist[k]; if (c > bestCount) { bestCount = c; f0 = f } }
        if (bestCount < Math.max(need, 0.3d * nResp)) { why = "the steps are spread over the recording"; return }
        int[] on = new int[bestCount]; double[] mean = new double[n]; int k = 0
        for (int i = 0; i < m; i++) if (resp[i] && at[i] >= f0 && at[i] < f0 + 15) { on[k++] = at[i]; for (int f = 0; f < n; f++) mean[f] += X[f][i] }
        Arrays.sort(on)
        int on0 = (bestCount % 2 == 1) ? on[(int) ((bestCount - 1) / 2)] : (int) Math.floor((on[(int) (bestCount / 2) - 1] + on[(int) (bestCount / 2)]) / 2.0d)
        int peak = on0; double top = -1e30
        for (int f = on0 - 1; f < Math.min(n, on0 + 60); f++) {
            double s = 0; for (int j = f - 2; j <= f + 2; j++) if (j >= 0 && j < n) s += mean[j]
            if (s > top) { top = s; peak = f + 1 }
        }
        int gapBefore = on0 > 40 ? 9 : 3                     // the frames before end well ahead of the rise, so early responders are not in them
        onset = on0; b1 = Math.max(2, on0 - gapBefore); b0 = Math.max(1, b1 - 19); p0 = Math.max(on0 + 1, peak - 9); p1 = Math.min(n, peak + 10); share = bestCount / (double) m
        if (share < 0.05d) { why = "only " + Math.round(100 * share) + " % of the field changes, near frame " + on0; return }
        found = true
    }

    /** Pixel indices inside a ROI, clipped to the picture. */
    int[] pixelsOf(Roi roi) {
        Rectangle r = roi.getBounds(); ImageProcessor mask = roi.getMask()
        int rx = (int) r.getX(), ry = (int) r.getY(), rw = (int) r.getWidth(), rh = (int) r.getHeight()
        int x0 = Math.max(0, rx), y0 = Math.max(0, ry), x1 = Math.min(w, rx + rw), y1 = Math.min(h, ry + rh)
        int[] out = new int[Math.max(0, (x1 - x0)) * Math.max(0, (y1 - y0))]; int k = 0
        for (int y = y0; y < y1; y++) for (int x = x0; x < x1; x++)
            if (mask == null || mask.get(x - rx, y - ry) != 0) out[k++] = y * w + x
        return Arrays.copyOf(out, k)
    }

    /** True when an outline comes within EDGE_MARGIN pixels of the picture's border. There its background band is cut off on one side
     *  and the lamp is weaker, so the background comes from better-lit glass than the cell sits on. */
    boolean nearEdge(Roi roi) {
        Rectangle b = roi.getBounds()
        return b.getX() < EDGE_MARGIN || b.getY() < EDGE_MARGIN || b.getX() + b.getWidth() > w - EDGE_MARGIN || b.getY() + b.getHeight() > h - EDGE_MARGIN
    }

    /** Pixels that belong to a ROI (widened): never background. */
    boolean[] takenBy(Roi[] cells) {
        boolean[] taken = new boolean[w * h]
        for (Roi c : cells) { Roi big = RoiEnlarger.enlarge(c, ROI_MARGIN); for (int i : pixelsOf(big)) taken[i] = true }
        for (Roi c : alsoAvoid) { Roi big = RoiEnlarger.enlarge(c, ROI_MARGIN); for (int i : pixelsOf(big)) taken[i] = true }
        return taken
    }

    /** Background pixels for one target: band GAP..REACH outside it, minus ROIs and bright objects, darkest quarter in the all-frames picture. */
    int[] backgroundOf(Roi target, boolean[] taken, boolean[] scratch) {
        int reach = REACH; int[] pick = new int[0]
        for (int attempt = 0; attempt < 2 && pick.length == 0; attempt++) {
            int[] innerPx = pixelsOf(RoiEnlarger.enlarge(target, GAP))
            for (int i : innerPx) scratch[i] = true
            int[] outer = pixelsOf(RoiEnlarger.enlarge(target, reach))
            int[] cand = new int[outer.length]; int k = 0
            for (int i : outer) if (!scratch[i] && !taken[i] && !bright[i]) cand[k++] = i
            for (int i : innerPx) scratch[i] = false
            if (k >= MIN_PIXELS) {
                float[] v = new float[k]
                for (int j = 0; j < k; j++) v[j] = flat[cand[j]]      // darkest in the all-frames picture = emptiest glass, on every side of the cell
                float[] s = v.clone(); Arrays.sort(s)
                float thr = s[Math.max(0, (int) (KEEP_DARKEST * k) - 1)]
                int[] p = new int[k]; int q = 0
                for (int j = 0; j < k; j++) if (v[j] <= thr) p[q++] = cand[j]
                pick = Arrays.copyOf(p, q)
            }
            reach = reach * 2
        }
        return pick
    }

    /** Step 3: the background of every cell. Returns one array of pixel indices per cell (empty = none found). */
    int[][] backgrounds(Roi[] cells) {
        boolean[] taken = takenBy(cells), scratch = new boolean[w * h]
        int[][] out = new int[cells.length][]
        for (int c = 0; c < cells.length; c++) out[c] = backgroundOf(cells[c], taken, scratch)
        return out
    }

    /** The background-corrected movie (32-bit copy; the original is untouched).
     *  "Empty glass" = the background pixels of every cell, plus the same kind of pixels found around a grid of points over
     *  the rest of the field. For every frame a smooth background picture is made from the empty glass only (near glass
     *  counts more: Gaussian weights, sigma MAP_SIGMA pixels) and subtracted from the frame. */
    static double MAP_SIGMA = 25, MAP_WIDE = 100, MAP_EPS = 0.02
    ImagePlus correctedMovie(Roi[] cells) {
        boolean[] taken = takenBy(cells), scratch = new boolean[w * h], empty = new boolean[w * h]
        for (Roi c : cells) for (int i : backgroundOf(c, taken, scratch)) empty[i] = true
        for (int gy = 24; gy < h; gy += 48) for (int gx = 24; gx < w; gx += 48)
            for (int i : backgroundOf(new ij.gui.OvalRoi(gx - 2, gy - 2, 4, 4), taken, scratch)) empty[i] = true
        int S = 4, sw = (int) (w / S), sh = (int) (h / S)
        float[] e = new float[sw * sh]
        for (int y = 0; y < sh * S; y++) { int row = ((int) (y / S)) * sw; for (int x = 0; x < sw * S; x++) if (empty[y * w + x]) e[row + (int) (x / S)] += 1f }
        FloatProcessor e1 = new FloatProcessor(sw, sh, (float[]) e.clone()), e2 = new FloatProcessor(sw, sh, (float[]) e.clone())
        e1.blurGaussian(MAP_SIGMA / S); e2.blurGaussian(MAP_WIDE / S)
        float[] d1 = (float[]) e1.getPixels(), d2 = (float[]) e2.getPixels()
        ImageStack src = movie.getStack(), out = new ImageStack(w, h)
        for (int f = 1; f <= n; f++) {
            ImageProcessor ip = src.getProcessor(f); float[] s = new float[sw * sh]
            for (int y = 0; y < sh * S; y++) { int row = ((int) (y / S)) * sw, o = y * w; for (int x = 0; x < sw * S; x++) if (empty[o + x]) s[row + (int) (x / S)] += ip.getf(o + x) }
            FloatProcessor s1 = new FloatProcessor(sw, sh, (float[]) s.clone()), s2 = new FloatProcessor(sw, sh, s)
            s1.blurGaussian(MAP_SIGMA / S); s2.blurGaussian(MAP_WIDE / S)
            float[] a1 = (float[]) s1.getPixels(), a2 = (float[]) s2.getPixels(); float[] b = new float[sw * sh]
            for (int i = 0; i < b.length; i++) { double den = d1[i] + MAP_EPS * d2[i]; b[i] = den > 1e-9 ? (float) ((a1[i] + MAP_EPS * a2[i]) / den) : 0f }
            FloatProcessor bs = new FloatProcessor(sw, sh, b); bs.setInterpolationMethod(ImageProcessor.BILINEAR)
            float[] big = (float[]) bs.resize(w, h).getPixels(); float[] o32 = new float[w * h]
            for (int i = 0; i < o32.length; i++) o32[i] = ip.getf(i) - big[i]
            out.addSlice(src.getSliceLabel(f), new FloatProcessor(w, h, o32))
        }
        String title = movie.getTitle(); int dot = title.lastIndexOf("."); if (dot > 0) title = title.substring(0, dot)
        ImagePlus res = new ImagePlus(title + "_background_corrected", out)
        res.setCalibration(movie.getCalibration().copy()); res.setDisplayRange(0, SHOW_MAX)
        say("Background-corrected movie (32-bit copy): for every frame, a smooth background picture made only from empty glass (the background pixels of the cells plus the same kind of pixels around a 48-pixel grid; Gaussian weights, sigma " + MAP_SIGMA + " px) was subtracted. Shown 0-" + SHOW_MAX + ". The tables use each cell's own background pixels, so Multi Measure on this movie gives close but not identical numbers.")
        return res
    }

    /** Evidence about a suggested outline, all measured, nothing guessed.
     *  Returns [grade, area, mean, brightShare, contrast, riseInside, riseAround, touchesEdge].
     *  grade 0 = responds on its own, 1 = likely a cell, 2 = unsure, 3 = probably precipitate, 4 = wrong size (not shown).
     *  On four reviewed recordings (1396 outlines) the share accepted as cells was 83 %, 80 %, 41 %, 12 % and 2 %.
     *  (Those shares were measured before 2026.11.20.9, when grade 0 did not yet include responding outlines with very bright parts.)
     *  Blind (cells chosen before the response): brightness comes from the picture of those frames and the response is not used, so grade 0 never occurs. */
    double[] grade(Roi roi) { return grade(roi, blind) }
    double[] grade(Roi roi, boolean blindGrade) {
        float[] pic = blindGrade && blindPic != null ? blindPic : flat
        float[] rise = blindGrade ? null : respDiff
        int[] px = pixelsOf(roi); int a = px.length
        if (a == 0) return [4d, 0d, 0d, 0d, 0d, 0d, 0d, 0d] as double[]
        boolean[] near = new boolean[w * h]; for (int i : pixelsOf(RoiEnlarger.enlarge(roi, 3))) near[i] = true
        double s = 0, sr = 0; int nb = 0; boolean edge = false
        for (int i : px) { s += pic[i]; if (pic[i] > BRIGHT_ABOVE) nb++; if (rise != null) sr += rise[i]; int x = i % w, y = (int) (i / w); if (x <= 1 || y <= 1 || x >= w - 2 || y >= h - 2) edge = true }
        double mean = s / a, brightShare = nb / (double) a, riseIn = sr / a
        double ringSum = 0, ringRise = 0; int nr = 0
        for (int i : pixelsOf(RoiEnlarger.enlarge(roi, 12))) if (!near[i]) { ringSum += pic[i]; if (rise != null) ringRise += rise[i]; nr++ }
        double contrast = nr > 0 ? mean - ringSum / nr : 0, riseRing = nr > 0 ? ringRise / nr : 0
        boolean own = rise != null && riseIn > 0.5d && riseIn > 1.25d * Math.max(riseRing, 0.05d)
        int g = 2
        boolean precip = mean > SUGGEST_BRIGHT || brightShare >= 0.25d
        if (precip) g = 3
        else {
            if (brightShare < 0.02d && contrast > 0.2d && !edge && mean <= BRIGHT_ABOVE) g = 1
            if (own && brightShare < 0.10d && mean <= BRIGHT_ABOVE) g = 0
        }
        // An outline that rises on its own at the stimulus is a cell even when part of it is very bright: precipitate does not
        // respond. Before 2026.11.20.9 such cells were graded unsure or precipitate and left out by "Accept likely" (on P3 Run2,
        // 17 of 58 clearly responding areas had no chosen cell; with this rule 6, of which 4 lie at the picture's edge). The rise must also be at least 5 % of the outline's brightness, so
        // flicker of a very bright lump does not count.
        if (g != 0 && own && riseIn >= 0.05d * mean) g = 0
        if (a < SUGGEST_MIN_AREA || a > SUGGEST_MAX_AREA) g = 4
        return [(double) g, (double) a, mean, brightShare, contrast, riseIn, riseRing, edge ? 1d : 0d] as double[]
    }

    /** Gray value over time of any outline, minus the background next to it (the same rule as for the tables). */
    double[] traceOf(Roi roi, Roi[] cells) {
        Roi[] all = new Roi[cells.length + 1]; System.arraycopy(cells, 0, all, 0, cells.length); all[cells.length] = roi
        int[] px = pixelsOf(roi), bg = backgroundOf(roi, takenBy(all), new boolean[w * h])
        double[] out = new double[n]; ImageStack st = movie.getStack()
        for (int f = 0; f < n; f++) {
            ImageProcessor ip = st.getProcessor(f + 1); double s = 0, b = 0
            for (int i : px) s += ip.getf(i)
            for (int i : bg) b += ip.getf(i)
            out[f] = (px.length > 0 ? s / px.length : 0) - (bg.length > 0 ? b / bg.length : 0)
        }
        return out
    }

    /** Area in pixels and mean brightness (all-frames picture) of a ROI. */
    double[] describe(Roi roi) {
        int[] px = pixelsOf(roi); double s = 0
        for (int i : px) s += flat[i]
        return [(double) px.length, px.length > 0 ? s / px.length : 0d] as double[]
    }

    /** The picture given to the automatic outliner: green = all-frames picture 0-8, magenta = response (if any). */
    /** True when two colour pictures are the same, allowing a difference of 2 in a channel (rounding can differ between computers). */
    static boolean samePicture(ImagePlus a, ImagePlus b) {
        if (a == null || b == null || a.getWidth() != b.getWidth() || a.getHeight() != b.getHeight() || a.getBitDepth() != 24 || b.getBitDepth() != 24) return false
        int[] p = (int[]) a.getProcessor().getPixels(), q = (int[]) b.getProcessor().getPixels()
        for (int i = 0; i < p.length; i++) {
            int x = p[i], y = q[i]
            if (x != y && (Math.abs(((x >> 16) & 255) - ((y >> 16) & 255)) > 2 || Math.abs(((x >> 8) & 255) - ((y >> 8) & 255)) > 2 || Math.abs((x & 255) - (y & 255)) > 2)) return false
        }
        return true
    }

    /** Blind: the picture of the frames before the response alone (see blindDetectionPicture). Otherwise the all-frames picture, with the response when one was found. */
    ImagePlus outlinePicture() { return blind ? blindOutlinePicture() : outlinePicture(found) }
    /** The response-assisted picture, used once by the sensitivity check after blind selection. */
    ImagePlus responseOutlinePicture() { return outlinePicture(true) }
    /** A picture Cellpose gets in blind mode, in 0..1: the frames before the response only, never the response, from 0 to cap gray levels
     *  with brighter pixels capped (BLIND_CAP and BLIND_CAP_HIGH; see blindOutlinesToAdd). For detection only; measurements are
     *  made on the raw movie. */
    float[] blindDetectionPicture() { return blindDetectionPicture(BLIND_CAP) }
    float[] blindDetectionPicture(double cap) {
        float[] v = new float[blindPic.length]
        for (int i = 0; i < v.length; i++) v[i] = (float) Math.max(0d, Math.min(1d, blindPic[i] / cap))
        return v
    }
    ImagePlus blindOutlinePicture() { return blindOutlinePicture(BLIND_CAP) }
    ImagePlus blindOutlinePicture(double cap) {
        float[] g = blindDetectionPicture(cap); int[] rgb = new int[w * h]
        for (int i = 0; i < rgb.length; i++) rgb[i] = ((int) Math.round(g[i] * 255)) << 8
        return new ImagePlus("for_automatic_outlines_blind", new ColorProcessor(w, h, rgb))
    }
    /** An outline too thin to be a cell body. Cellpose sometimes draws a strip along the rim of a cell it has already outlined
     *  (seen 20 Nov 2026: a 3-pixel-wide crescent under a cell, accepted as a second cell). "Thin" is a mean width under
     *  THIN_WIDTH pixels (2 x area / perimeter, which is the radius of a round outline) or a roundness under THIN_ROUNDNESS
     *  (4 pi area / perimeter^2). Measured on the reviewed and audited outlines of four recordings: none is this thin. */
    static final double THIN_WIDTH = 4d, THIN_ROUNDNESS = 0.2d
    static boolean tooThin(Roi roi) {
        ij.process.FloatPolygon p = roi.getFloatPolygon(); int n = p.npoints
        if (n < 3) return true
        double area = 0d, length = 0d
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n
            area += (double) p.xpoints[i] * (double) p.ypoints[j] - (double) p.xpoints[j] * (double) p.ypoints[i]
            length += Math.hypot((double) p.xpoints[j] - (double) p.xpoints[i], (double) p.ypoints[j] - (double) p.ypoints[i])
        }
        area = Math.abs(area) / 2d
        if (length <= 0d || area <= 0d) return true
        return 2d * area / length < THIN_WIDTH || 4d * Math.PI * area / (length * length) < THIN_ROUNDNESS
    }
    /** A thin outline lying along another outline's edge: at least THIN_ON_RIM of its pixels are within THIN_REACH pixels of
     *  the other outline. Such a strip is the other cell's rim, not a cell of its own. */
    static final int THIN_REACH = 5
    static final double THIN_ON_RIM = 0.6d
    static boolean alongAnother(Roi thin, List<Roi> others) {
        Rectangle tb = thin.getBounds(); Rectangle grown = new Rectangle((int) tb.getX() - THIN_REACH, (int) tb.getY() - THIN_REACH, (int) tb.getWidth() + 2 * THIN_REACH, (int) tb.getHeight() + 2 * THIN_REACH)
        java.awt.Point[] points = thin.getContainedPoints()
        if (points.length == 0) return false
        int[][] reach = [[0, 0], [THIN_REACH, 0], [-THIN_REACH, 0], [0, THIN_REACH], [0, -THIN_REACH], [3, 3], [3, -3], [-3, 3], [-3, -3]] as int[][]
        for (Roi other : others) {
            if (other.is(thin) || !other.getBounds().intersects(grown)) continue
            int near = 0
            for (java.awt.Point q : points) {
                int x = (int) q.getX(), y = (int) q.getY()
                for (int[] step : reach) if (other.contains(x + step[0], y + step[1])) { near++; break }
            }
            if (near >= THIN_ON_RIM * points.length) return true
        }
        return false
    }
    /** Blind detection runs Cellpose on two pictures (0 to BLIND_CAP and 0 to BLIND_CAP_HIGH gray levels). Every outline from the
     *  first is kept. An outline from the second is added when, for every outline from the first, the pixels they share are fewer
     *  than BLIND_MERGE_OVERLAP of the smaller one's pixels. Returns the positions (in second) of the outlines to add. */
    static List<Integer> blindOutlinesToAdd(List<Roi> first, List<Roi> second, int width) {
        List<Set<Integer>> firstPixels = new ArrayList<Set<Integer>>()
        for (Roi r : first) firstPixels.add(pixelSet(r, width))
        List<Integer> add = new ArrayList<Integer>()
        for (int j = 0; j < second.size(); j++) {
            Roi s = second.get(j); Set<Integer> sp = pixelSet(s, width); Rectangle sb = s.getBounds()
            boolean same = false
            for (int i = 0; i < first.size() && !same; i++) {
                if (!first.get(i).getBounds().intersects(sb)) continue
                Set<Integer> fp = firstPixels.get(i)
                int shared = 0
                for (Integer p : sp) if (fp.contains(p)) shared++
                same = shared >= BLIND_MERGE_OVERLAP * Math.min(sp.size(), fp.size())
            }
            if (!same) add.add(j)
        }
        return add
    }
    static Set<Integer> pixelSet(Roi r, int width) {
        Set<Integer> s = new HashSet<Integer>()
        for (java.awt.Point q : r.getContainedPoints()) s.add(((int) q.getY()) * width + (int) q.getX())
        return s
    }
    ImagePlus outlinePicture(boolean withResponse) {
        float[] g = (float[]) allPic.getProcessor().getPixels()
        float[] m = withResponse && found ? (float[]) respPic.getProcessor().getPixels() : null
        int[] rgb = new int[w * h]
        for (int i = 0; i < rgb.length; i++) {
            int gv = (int) Math.max(0, Math.min(255, Math.round(g[i] / SHOW_MAX * 255)))
            int mv = m == null ? 0 : (int) Math.max(0, Math.min(255, Math.round((m[i] - RESP_MIN) / (RESP_MAX - RESP_MIN) * 255)))
            rgb[i] = (mv << 16) | (gv << 8) | mv
        }
        return new ImagePlus("for_automatic_outlines", new ColorProcessor(w, h, rgb))
    }

    /** A picture-sized overlay that paints the given pixels translucent blue. */
    ImageRoi blueLayer(int[][] bg, int only) {
        ColorProcessor cp = new ColorProcessor(w, h); int[] px = (int[]) cp.getPixels()
        ByteProcessor alpha = new ByteProcessor(w, h); byte[] al = (byte[]) alpha.getPixels()
        for (int c = 0; c < bg.length; c++) {
            if (only >= 0 && c != only) continue
            for (int i : bg[c]) { px[i] = 0x2896ff; al[i] = (byte) 255 }
        }
        java.awt.image.BufferedImage bi = new java.awt.image.BufferedImage(w, h, java.awt.image.BufferedImage.TYPE_INT_ARGB)
        int[] argb = new int[w * h]
        for (int i = 0; i < argb.length; i++) argb[i] = al[i] != 0 ? ((0x99 << 24) | px[i]) : 0
        bi.setRGB(0, 0, w, h, argb, 0, w)
        ImageRoi ir = new ImageRoi(0, 0, bi); ir.setZeroTransparent(false)
        return ir
    }

    /** Step 4: measure on the original movie and save. Returns a one-line summary. */
    String measureAndSave(Roi[] cells, String[] names, String prefix) {
        int nc = cells.length
        int[][] cellPx = new int[nc][]; for (int c = 0; c < nc; c++) cellPx[c] = pixelsOf(cells[c])
        int[][] bg = backgrounds(cells)
        lastBackgroundPixels = new int[nc]
        for (int c = 0; c < nc; c++) lastBackgroundPixels[c] = bg[c].length
        double[][] raw = new double[n][nc], back = new double[n][nc]
        ImageStack st = movie.getStack()
        for (int f = 0; f < n; f++) {
            ImageProcessor ip = st.getProcessor(f + 1)
            for (int c = 0; c < nc; c++) {
                double s = 0; for (int i : cellPx[c]) s += ip.getf(i); raw[f][c] = cellPx[c].length > 0 ? s / cellPx[c].length : Double.NaN
                s = 0; for (int i : bg[c]) s += ip.getf(i); back[f][c] = bg[c].length > 0 ? s / bg[c].length : Double.NaN
            }
        }
        StringBuilder a = new StringBuilder("Frame"), b = new StringBuilder("Frame"), d = new StringBuilder("Frame")
        for (String nm : names) { a.append(",").append(nm); b.append(",").append(nm); d.append(",").append(nm) }
        a.append("\n"); b.append("\n"); d.append("\n")
        for (int f = 0; f < n; f++) {
            a.append(f + 1); b.append(f + 1); d.append(f + 1)
            for (int c = 0; c < nc; c++) {
                a.append(",").append(String.format(Locale.US, "%.3f", raw[f][c]))
                if (bg[c].length > 0) { b.append(",").append(String.format(Locale.US, "%.3f", back[f][c])); d.append(",").append(String.format(Locale.US, "%.3f", raw[f][c] - back[f][c])) }
                else { b.append(","); d.append(",") }
            }
            a.append("\n"); b.append("\n"); d.append("\n")
        }
        new File(prefix + "gray_values_raw.csv").text = a.toString()
        new File(prefix + "background_next_to_each_cell.csv").text = b.toString()
        new File(prefix + "gray_values_background_corrected_for_SimpleCa.csv").text = d.toString()
        // Native analysis preserves the measured doubles. Legacy three-decimal exports remain compatible.
        new File(prefix + "gray_values_background_corrected_native.csv").setText(CalciumMeasurement.fullPrecision(raw, back, names), "UTF-8")
        // Compatibility only: an arbitrary additive offset changes fractional responses and must not
        // be used to rescue a denominator or classify responses with a fractional amplitude threshold.
        double[] lift = new double[nc]; lastLifted = 0
        for (int c = 0; c < nc; c++) {
            if (bg[c].length == 0) continue
            double lowest = Double.MAX_VALUE; int win = Math.min(20, n), upTo = found ? Math.max(win, b1) : n     // up to the last frame before the response (what any baseline uses)
            for (int f = 0; f + win <= upTo; f++) { double v = 0; for (int j = 0; j < win; j++) v += raw[f + j][c] - back[f + j][c]; v /= win; if (v < lowest) lowest = v }
            if (lowest < KEEP_LEVEL) { lift[c] = KEEP_LEVEL - lowest; lastLifted++ }
        }
        StringBuilder k = new StringBuilder("Frame"); for (String nm : names) k.append(",").append(nm); k.append("\n")
        for (int f = 0; f < n; f++) {
            k.append(f + 1)
            for (int c = 0; c < nc; c++) { if (bg[c].length > 0) k.append(",").append(String.format(Locale.US, "%.3f", raw[f][c] - back[f][c] + lift[c])); else k.append(",") }
            k.append("\n")
        }
        new File(prefix + "for_SimpleCa_every_cell_kept.csv").text = k.toString()
        ZipOutputStream z = new ZipOutputStream(new FileOutputStream(prefix + "cells_RoiSet.zip"))
        for (int c = 0; c < nc; c++) { Roi r = (Roi) cells[c].clone(); r.setName(names[c]); r.setPosition(0); z.putNextEntry(new ZipEntry(names[c] + ".roi")); z.write(RoiEncoder.saveAsByteArray(r)); z.closeEntry() }
        z.close()
        int none = 0; for (int c = 0; c < nc; c++) if (bg[c].length == 0) none++
        // per-cell summary: level before the response (or frames 1-20), peak, rise
        int s0 = found ? b0 : 1, s1 = found ? b1 : Math.min(20, n)
        StringBuilder sm = new StringBuilder("cell,x,y,area_pixels,background_pixels,raw_level,background_level,cell_above_background,peak_above_background,rise_gray_levels,peak_frame,dF_over_F0,rise_over_noise,clear_rise,not_brighter_than_background_frames_1_20,lifted_by_gray_levels_in_every_cell_kept_table,within_" + EDGE_MARGIN + "_px_of_picture_edge\n")
        lastClear = 0; lastMeanRise = 0; lastDim = 0
        for (int c = 0; c < nc; c++) {
            Rectangle r = cells[c].getBounds(); double rl = 0, bl = 0
            for (int f = s0 - 1; f < s1; f++) { rl += raw[f][c]; bl += back[f][c] }
            rl /= (s1 - s0 + 1); bl /= (s1 - s0 + 1)
            double F0 = rl - bl, pk = -1e30; int pf = 0
            double ss = 0; for (int f = s0 - 1; f < s1; f++) { double dv = (raw[f][c] - back[f][c]) - F0; ss += dv * dv }
            double sd = Math.sqrt(ss / Math.max(1, s1 - s0))
            for (int f = s1; f + 4 < n; f++) { double v = 0; for (int j = 0; j < 5; j++) v += raw[f + j][c] - back[f + j][c]; v /= 5; if (v > pk) { pk = v; pf = f + 3 } }
            sm.append(names[c]).append(",").append((int) Math.round(r.getCenterX())).append(",").append((int) Math.round(r.getCenterY())).append(",").append(cellPx[c].length).append(",").append(bg[c].length).append(",")
              .append(String.format(Locale.US, "%.2f,%.2f,%.2f,%.2f,%.2f,%d,%s", rl, bl, F0, pk, pk - F0, pf, F0 > 0.5d ? String.format(Locale.US, "%.3f", (pk - F0) / F0) : ""))
            double first = 0; int nf = Math.min(20, n); for (int f = 0; f < nf; f++) first += raw[f][c] - back[f][c]
            boolean dim = first / nf <= 0d            // not brighter than its background in frames 1-20: no dF/F0; SimpleCa2+ leaves such cells out
            if (dim) lastDim++
            boolean clear = sd > 0 && (pk - F0) / sd >= 5.0d
            if (clear) lastClear++
            lastMeanRise += (pk - F0) / nc
            sm.append(String.format(Locale.US, ",%.1f,%s,%s,%.2f,%s", sd > 0 ? (pk - F0) / sd : 0d, clear ? "yes" : "no", dim ? "yes" : "no", lift[c], nearEdge(cells[c]) ? "yes" : "no")).append("\n")
        }
        new File(prefix + "summary_per_cell.csv").text = sm.toString()
        // check picture: all-frames picture, cells yellow, background blue
        ImageProcessor base = allPic.getProcessor().duplicate(); base.setMinAndMax(0, SHOW_MAX)
        ImagePlus chk = new ImagePlus("check", base.convertToByte(true).convertToRGB())
        Overlay ov = new Overlay(); ov.add(blueLayer(bg, -1))
        for (Roi c : cells) { Roi r = (Roi) c.clone(); r.setPosition(0); r.setStrokeColor(Color.yellow); r.setFillColor(null); ov.add(r) }
        chk.setOverlay(ov); ImagePlus flat = chk.flatten()
        IJ.saveAs(flat, "PNG", prefix + "cells_and_background_areas.png")
        ImagePlus areas = new ImagePlus("background areas", base.convertToByte(true).convertToRGB()); Overlay only = new Overlay(); only.add(blueLayer(bg, -1)); areas.setOverlay(only); backgroundPic = areas.flatten()
        IJ.saveAsTiff(allPic, prefix + "picture_all_frames.tif")
        // the all-frames colour picture, whatever view is on screen: one colour photo, and the same with the cells outlined and a key
        IJ.saveAsTiff(colourPhoto(), prefix + "picture_colour.tif")
        IJ.saveAs(colourPhotoWithCells(cells), "PNG", prefix + "picture_colour_cells_outlined.png")
        say("Legacy screening rise (not a stimulus responder call): highest 5-frame average after frames " + s0 + "-" + s1 + " at least 5 baseline standard deviations above their mean: " + lastClear + " of " + nc + " cells. Use Step 4 for protocol-defined responder counts. Mean rise " + String.format(Locale.US, "%.2f", lastMeanRise) + " gray levels.")
        if (lastDim > 0) say(lastDim + " cells are not brighter than their own background in frames 1-20. Step 4 retains their registered records and reports technical exclusions separately from nonresponders.")
        say("Compatibility export for_SimpleCa_every_cell_kept.csv adds constants to " + lastLifted + " dim traces. These offsets alter fractional response thresholds. Do not use this export for fractional amplitudes or responder percentages. Step 4 uses full-precision, unlifted fluorescence.")
        String summary = nc + " cells measured on " + n + " frames; background found for " + (nc - none) + (none > 0 ? " (" + none + " cells had too little empty space nearby: their corrected column is empty)" : "")
        say("Background for each cell: band from " + GAP + " to " + REACH + " pixels outside the ROI, without other ROIs (+" + ROI_MARGIN + " px) and without pixels brighter than " + BRIGHT_ABOVE + " gray levels in the all-frames picture (+" + BRIGHT_GROW + " px); the darkest " + Math.round(100 * KEEP_DARKEST) + " % of those pixels are used, ranked in the all-frames picture, where the uneven lamp light is already removed (rule changed 19 Nov: ranking on raw brightness picked one side of each cell).")
        log.removeAll { String s -> s.startsWith("Cells at the picture's edge") }      // only the latest measurement's choice stays in the record
        if (edgeNote != "") say(edgeNote)
        say("Measured on the original movie: mean gray value of every cell and of its background, every frame. corrected = raw - background.")
        say(summary)
        new File(prefix + "what_was_done.txt").text = log.join("\n") + "\n"
        return summary
    }
}

/** Sensitivity check after blind selection. Response-assisted detection is run once, for comparison only, and each of its
 *  outlines is looked up among the blind suggestions. Found = a blind outline contains its centre, or it contains the centre of
 *  a blind outline (the rule the window uses to match accepted cells to suggestions). The comparison outlines are never added
 *  to the cells or measured. */
class CalciumSensitivity {
    static final List<String> GRADES = ["responds on its own", "likely a cell", "unsure", "very bright: probably precipitate"]
    static Set<Long> pointsOf(Roi r) { Set<Long> s = new HashSet<Long>(); r.getContainedPoints().each { java.awt.Point p -> s.add((((long) (int) p.getY()) << 32) | (((int) p.getX()) & 0xffffffffL)) }; s }
    /** Both lists hold suggestions as the window keeps them: maps with roi, grade, cx, cy and edge. */
    static Map compare(List assisted, List blindSuggestions, Roi[] measured, boolean leaveOutEdge) {
        def counted = assisted.findAll { !(leaveOutEdge && it.edge) }
        def rows = []
        counted.eachWithIndex { a, int i ->
            Roi r = (Roi) a.roi; Rectangle rb = r.getBounds()
            boolean found = blindSuggestions.any { b -> ((Roi) b.roi).contains((int) a.cx, (int) a.cy) || r.contains((int) b.cx, (int) b.cy) }
            double best = 0; Set<Long> mine = null
            blindSuggestions.each { b ->
                Roi o = (Roi) b.roi
                if (!o.getBounds().intersects(rb)) return
                if (mine == null) mine = pointsOf(r)
                Set<Long> theirs = pointsOf(o); int both = 0; theirs.each { if (mine.contains(it)) both++ }
                double iou = both / (double) (mine.size() + theirs.size() - both); if (iou > best) best = iou
            }
            boolean inCells = measured.any { Roi c -> c.contains((int) a.cx, (int) a.cy) || r.contains((int) c.getBounds().getCenterX(), (int) c.getBounds().getCenterY()) }
            rows << [n: i + 1, grade: (int) a.grade, x: (int) a.cx, y: (int) a.cy, area: (int) (a.ev != null ? a.ev[1] : r.getContainedPoints().length), found: found, overlap: best, inCells: inCells, roi: r]
        }
        def cellLike = rows.findAll { it.grade <= 2 }
        return [rows: rows, cellLike: cellLike.size(), found: cellLike.count { it.found }, bright: rows.count { it.grade == 3 }, brightFound: rows.count { it.grade == 3 && it.found },
                byGrade: (0..3).collect { int g -> [rows.count { it.grade == g }, rows.count { it.grade == g && it.found }] },
                inCells: cellLike.count { it.inCells }, blindShown: blindSuggestions.count { !(leaveOutEdge && it.edge) }, assistedShown: assisted.size(), edgeLeftOut: assisted.size() - counted.size()]
    }
    static String pct(int a, int b) { b > 0 ? Math.round(100d * a / b) + " %" : "-" }
    /** Plain-language report, for the record of steps and the result line. */
    static String report(Map m, String how) {
        int t = (int) m.cellLike, f = (int) m.found
        def by = (0..2).collect { int g -> def c = m.byGrade[g]; GRADES[g] + " " + c[1] + " of " + c[0] }.join("; ")
        return "Sensitivity check (response-assisted detection run once for comparison; " + how + "; its outlines were not added to the cells and were not measured): " +
            t + " response-assisted outlines could be cells (graded responds, likely or unsure" + (((int) m.edgeLeftOut) > 0 ? ", not at the edge" : "") + "). " +
            "The blind detection also found " + f + " of them (" + pct(f, t) + ") and missed " + (t - f) + ". By grade: " + by + ". " +
            "Very bright (probably precipitate): " + m.brightFound + " of " + m.bright + " found. " +
            m.inCells + " of the " + t + " coincide with a measured cell. Blind detection made " + m.blindShown + " suggestions. " +
            "Found = a blind outline contains its centre or it contains the centre of a blind outline. " +
            "How to read it: if most are found, the blind picture picks up dim cells well. If many are missed, blind detection misses cells that are dim at rest; " +
            "% responding is still unbiased if those cells respond at the same rate as the others, but that assumption then carries weight and should be stated with the result."
    }
    static void save(Map m, String prefix) {
        def csv = new StringBuilder("outline,grade_with_response,x,y,area_pixels,found_by_blind_detection,best_overlap_with_a_blind_outline,coincides_with_a_measured_cell\n")
        m.rows.each { r -> csv.append(String.format(Locale.US, "%d,%s,%d,%d,%d,%s,%.2f,%s\n", r.n, GRADES[(int) r.grade], r.x, r.y, r.area, r.found ? "yes" : "no", r.overlap, r.inCells ? "yes" : "no")) }
        new File(prefix + "sensitivity_check.csv").setText(csv.toString(), "UTF-8")
        new ZipOutputStream(new FileOutputStream(prefix + "sensitivity_check_missed_RoiSet.zip")).withCloseable { ZipOutputStream z ->
            m.rows.findAll { !it.found && (int) it.grade <= 2 }.each { r ->   // the missed ones that could be cells, as counted in the report
                Roi o = (Roi) ((Roi) r.roi).clone(); String name = String.format(Locale.US, "missed_%03d_%s", r.n, GRADES[(int) r.grade].replaceAll(/[^A-Za-z]+/, "_"))
                o.setName(name); o.setPosition(0); z.putNextEntry(new ZipEntry(name + ".roi")); z.write(RoiEncoder.saveAsByteArray(o)); z.closeEntry()
            }
        }
    }
}

/** Release lookup sends no recording data. Installation requires an explicit button click. */
class CalciumRelease {
    static final String VERSION = "2026.11.20.19"
    static final String REPOSITORY = "https://github.com/samgillman/fiji-calcium-plugin"
    static final String API = "https://api.github.com/repos/samgillman/fiji-calcium-plugin/releases/latest"
    static final String ASSET = "Calcium_Imaging_Fiji_Plugin.zip"
    static final String SCRIPT = "Calcium_Cells_Step_By_Step.groovy"
    static final String CHECKSUMS = "SHA256SUMS.txt"

    static List<Integer> versionParts(String version) {
        String clean = version.replaceFirst(/^v/, "")
        if (!(clean ==~ /\d{4}\.\d{1,2}\.\d{1,2}\.\d{1,4}/)) throw new IOException("Unrecognised release version")
        clean.split(/\./).collect { Integer.parseInt(it) }
    }
    static boolean newer(String candidate, String installed) {
        def a = versionParts(candidate), b = versionParts(installed)
        for (int i = 0; i < a.size(); i++) { if (a[i] != b[i]) return a[i] > b[i] }
        false
    }
    static Map describe(Map release) {
        if (release.draft || release.prerelease) throw new IOException("No stable release found")
        String tag = release.tag_name ?: ""
        versionParts(tag)
        String page = REPOSITORY + "/releases/tag/" + tag
        if (release.html_url != page) throw new IOException("Unexpected release address")
        def asset = release.assets?.find { it.name == ASSET && it.state == "uploaded" && it.size > 0 }
        if (asset == null || asset.browser_download_url != REPOSITORY + "/releases/download/" + tag + "/" + ASSET)
            throw new IOException("This release does not have a ready plugin download yet")
        def ready = { String name -> release.assets?.find { it.name == name && it.state == "uploaded" && it.size > 0 &&
            it.browser_download_url == REPOSITORY + "/releases/download/" + tag + "/" + name } }
        def script = ready(SCRIPT), sums = ready(CHECKSUMS)
        [version: tag.replaceFirst(/^v/, ""), newer: newer(tag, VERSION), page: page,
         installable: script != null && sums != null, script: script?.browser_download_url,
         checksums: sums?.browser_download_url, size: script?.size]
    }
    static Map latest() {
        def connection = (java.net.HttpURLConnection) new URL(API).openConnection()
        connection.setConnectTimeout(8000); connection.setReadTimeout(8000)
        connection.setRequestProperty("Accept", "application/vnd.github+json")
        connection.setRequestProperty("User-Agent", "Fiji-Calcium-Plugin/" + VERSION)
        try {
            int code = connection.getResponseCode()
            if (code == 404) throw new IOException("No public release is available yet")
            if (code == 403 || code == 429) throw new IOException("GitHub temporarily limited update checks; try again later")
            if (code != 200) throw new IOException("GitHub returned HTTP " + code)
            byte[] bytes = connection.getInputStream().withCloseable { input ->
                def buffer = new ByteArrayOutputStream(); byte[] block = new byte[8192]; int n
                while ((n = input.read(block)) != -1) {
                    if (buffer.size() + n > 1048576) throw new IOException("The update response was too large")
                    buffer.write(block, 0, n)
                }
                buffer.toByteArray()
            }
            def release = new groovy.json.JsonSlurper().parseText(new String(bytes, "UTF-8"))
            if (!(release instanceof Map)) throw new IOException("Unexpected update response")
            return describe(release)
        } finally { connection.disconnect() }
    }
}

/** Verify the complete download before touching the installed file; commit with one atomic rename. */
class CalciumInstaller {
    File applicationDirectory
    Closure download = { String url, int limit, Closure progress -> fetch(url, limit, progress) }
    Closure replaceFile = { java.nio.file.Path staged, java.nio.file.Path destination ->
        java.nio.file.Files.move(staged, destination, java.nio.file.StandardCopyOption.ATOMIC_MOVE, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    }
    static final int MAX_SCRIPT = 2 * 1024 * 1024
    static String escape(String value) { value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;") }
    static String digest(byte[] bytes) {
        java.security.MessageDigest.getInstance("SHA-256").digest(bytes).collect { String.format("%02x", it & 255) }.join("")
    }
    static String scriptVersion(byte[] bytes) {
        def match = new String(bytes, "UTF-8") =~ /static final String VERSION = "(\d{4}\.\d{1,2}\.\d{1,2}\.\d{1,4})"/
        match.find() ? match.group(1) : null
    }
    static byte[] fetch(String address, int limit, Closure progress) {
        URL url = new URL(address)
        // GitHub release downloads redirect to GitHub's asset host. Every hop must stay on HTTPS.
        for (int redirects = 0; redirects < 6; redirects++) {
            if (url.protocol != "https") throw new IOException("The download must use HTTPS")
            def connection = (java.net.HttpURLConnection) url.openConnection()
            connection.setInstanceFollowRedirects(false)
            connection.setConnectTimeout(10000); connection.setReadTimeout(20000)
            connection.setRequestProperty("User-Agent", "Fiji-Calcium-Plugin/" + CalciumRelease.VERSION)
            try {
                int code = connection.responseCode
                if (code in [301,302,303,307,308]) {
                    String location = connection.getHeaderField("Location")
                    if (!location) throw new IOException("The download redirect is missing its address")
                    url = new URL(url, location); continue
                }
                if (code != 200) throw new IOException("Download returned HTTP " + code)
                long total = connection.getContentLengthLong()
                if (total > limit) throw new IOException("The download exceeds the expected size")
                return connection.inputStream.withCloseable { input ->
                    def output = new ByteArrayOutputStream(); byte[] buffer = new byte[8192]; int count
                    while ((count = input.read(buffer)) != -1) {
                        if (output.size() + count > limit) throw new IOException("The download exceeds the expected size")
                        output.write(buffer, 0, count); progress.call(output.size(), total)
                    }
                    if (total >= 0 && output.size() != total) throw new IOException("The download was incomplete")
                    output.toByteArray()
                }
            } finally { connection.disconnect() }
        }
        throw new IOException("Too many download redirects")
    }
    Map install(Map release, Closure progress = { int p, String message -> }) {
        String version = release.version ?: ""
        CalciumRelease.versionParts(version)
        String base = CalciumRelease.REPOSITORY + "/releases/download/v" + version + "/"
        if (!release.installable || release.script != base + CalciumRelease.SCRIPT || release.checksums != base + CalciumRelease.CHECKSUMS)
            throw new IOException("This release does not have verified installation files")
        if (applicationDirectory == null || !applicationDirectory.isDirectory()) throw new IOException("Could not locate this Fiji installation")
        File app = applicationDirectory.canonicalFile
        File target = new File(app, "scripts/Plugins/Calcium Imaging/" + CalciumRelease.SCRIPT)
        if (!target.canonicalPath.startsWith(app.path + File.separator) || java.nio.file.Files.isSymbolicLink(target.toPath()))
            throw new IOException("The plugin installation path points outside this Fiji folder or uses a symbolic link")
        progress.call(5, "Downloading the release checksum...")
        byte[] manifest = download.call(release.checksums, 65536, { a,b -> })
        def checks = new String(manifest, "UTF-8").readLines().findAll { it ==~ /[a-fA-F0-9]{64}  Calcium_Cells_Step_By_Step\.groovy/ }
        if (checks.size() != 1) throw new IOException("The release checksum is missing or ambiguous")
        String expected = checks[0].substring(0,64).toLowerCase(Locale.ROOT)
        progress.call(15, "Downloading plugin " + version + "...")
        byte[] bytes = download.call(release.script, MAX_SCRIPT, { count,total ->
            long size = total > 0 ? total : (release.size ?: 0) as long
            if (size > 0) progress.call(15 + (int) Math.min(65L, (long) (65L * count / size)), "Downloading plugin " + version + "...")
        })
        progress.call(85, "Verifying the download...")
        if (digest(bytes) != expected) throw new IOException("Checksum verification failed; the installed plugin was not changed")
        if (scriptVersion(bytes) != version || !new String(bytes,"UTF-8").contains('static final String REPOSITORY = "' + CalciumRelease.REPOSITORY + '"'))
            throw new IOException("The downloaded plugin does not match the selected release")
        def lockPath = new File(app, ".calcium-plugin-update.lock").toPath()
        java.nio.channels.FileChannel channel = null
        java.nio.channels.FileLock lock = null
        java.nio.file.Path staged = null
        try {
            channel = java.nio.channels.FileChannel.open(lockPath, java.nio.file.StandardOpenOption.CREATE, java.nio.file.StandardOpenOption.WRITE)
            try { lock = channel.tryLock() } catch (java.nio.channels.OverlappingFileLockException ignored) { }
            if (lock == null) throw new IOException("Another plugin installation is in progress; wait for it to finish")
            byte[] previous = target.isFile() ? java.nio.file.Files.readAllBytes(target.toPath()) : null
            String diskVersion = previous == null ? null : scriptVersion(previous)
            if (diskVersion && CalciumRelease.newer(diskVersion, version)) throw new IOException("A newer plugin is already installed; check for updates again")
            if (previous != null && digest(previous) == expected) return [version:version, target:target, backup:null, already:true]
            java.nio.file.Files.createDirectories(target.parentFile.toPath())
            staged = java.nio.file.Files.createTempFile(target.parentFile.toPath(), ".calcium-update-", ".tmp")
            java.nio.channels.FileChannel.open(staged, java.nio.file.StandardOpenOption.WRITE).withCloseable { output ->
                def buffer = java.nio.ByteBuffer.wrap(bytes)
                while (buffer.hasRemaining()) output.write(buffer)
                output.force(true)
            }
            File backup = null
            if (previous != null) {
                progress.call(92, "Backing up the installed plugin...")
                File backups = new File(app, "calcium-plugin-backups")
                java.nio.file.Files.createDirectories(backups.toPath())
                backup = new File(backups, "Calcium-" + (diskVersion ?: "previous") + "-" + UUID.randomUUID() + ".groovy.bak")
                java.nio.file.Files.write(backup.toPath(), previous, java.nio.file.StandardOpenOption.CREATE_NEW, java.nio.file.StandardOpenOption.WRITE)
            }
            // Do not overwrite an external edit made while this update was being staged.
            if (previous == null ? target.exists() : (!target.isFile() || digest(java.nio.file.Files.readAllBytes(target.toPath())) != digest(previous)))
                throw new IOException("The installed file changed during the update; please try again")
            progress.call(98, "Installing the verified plugin...")
            replaceFile.call(staged, target.toPath())
            progress.call(100, "Installed " + version + ". Save your work and restart Fiji.")
            return [version:version, target:target, backup:backup, already:false]
        } catch (java.nio.file.AccessDeniedException denied) {
            throw new IOException("Fiji cannot write to its installation folder. Use a Fiji folder your account can modify, then try again.", denied)
        } catch (java.nio.file.AtomicMoveNotSupportedException unsupported) {
            throw new IOException("This folder does not support safe replacement. The existing plugin was kept.", unsupported)
        } finally {
            if (staged != null) java.nio.file.Files.deleteIfExists(staged)
            if (lock != null) lock.release()
            if (channel != null) channel.close()
        }
    }
}

/** Background requests with results delivered on Swing's event thread. */
class CalciumUpdateCheck {
    Closure fetch = { -> CalciumRelease.latest() }
    String installed = CalciumRelease.VERSION
    volatile boolean enabled = true, closed = false
    java.util.concurrent.atomic.AtomicBoolean busy = new java.util.concurrent.atomic.AtomicBoolean(false)
    Set<String> offered = new HashSet<String>() // accessed only on the event thread

    boolean check(boolean automatic, Closure receive) {
        if (closed || (automatic && !enabled) || !busy.compareAndSet(false, true)) return false
        Thread.startDaemon("Calcium plugin update check") {
            Map result
            try {
                result = new LinkedHashMap(fetch.call())
                result.newer = CalciumRelease.newer(result.version, installed)
            } catch (Exception problem) { result = [error: problem.message ?: problem.toString()] }
            final Map answer = result
            javax.swing.SwingUtilities.invokeLater {
                busy.set(false)
                if (!closed && (!automatic || enabled)) {
                    answer.prompt = answer.newer && (!automatic || offered.add(answer.version))
                    receive.call(answer)
                }
            }
        }
        true
    }
    void close() { closed = true }
}

/** Stage progress comes from actual work/log events, never from elapsed time.
 *  Cellpose's single-image inference does not expose a reliable overall percentage. */
class CalciumOutlineTask {
    Closure progress
    volatile boolean cancelled = false
    volatile boolean cancellable = true
    boolean succeeded = false
    volatile Process process
    volatile File logFile
    int percent = 0
    String stage = "Starting", detail = "", device = ""
    boolean downloading = false
    String source = "Checking for saved outlines", settings = "", imageSize = ""
    long began = System.nanoTime(), stageStarted = began, lastEvent = began
    List<String> history = []

    synchronized void remember(String message) {
        lastEvent = System.nanoTime()
        long seconds = (long) ((lastEvent - began) / 1000000000L)
        history.add(String.format(Locale.US, "%d:%02d  %s", (long) (seconds / 60), seconds % 60, message))
        if (history.size() > 120) history.remove(0)
    }
    synchronized Map snapshot() {
        [percent: percent, stage: stage, detail: detail, at: lastEvent, stageAt: stageStarted,
         context: [source, imageSize, device, settings].findAll { it }.join(" | "),
         history: history.join("\n")]
    }

    void checkCancelled() {
        if (cancelled || Thread.currentThread().isInterrupted()) throw new java.util.concurrent.CancellationException("Detection cancelled")
    }
    // Blind detection runs Cellpose twice; each run's 0-84 % is drawn in its own part of the bar (partFrom + value x partShare).
    int partFrom = 0; double partShare = 1; int lastAsked = 0
    synchronized void report(int value, String what, String extra = "") {
        checkCancelled()
        lastAsked = value
        percent = Math.max(percent, Math.min(99, (int) Math.round(partFrom + value * partShare)))
        if (stage != what) stageStarted = System.nanoTime()
        if (stage != what || detail != extra) remember(what + (extra ? ": " + extra : ""))
        stage = what; detail = extra
        progress.call(percent, stage, detail)
    }
    void cancel() { cancelled = true } // the worker stops the process, keeping the UI responsive

    void logLine(String raw) {
        String line = raw.replaceAll(/\u001B\[[0-?]*[ -\/]*[@-~]/, "").trim()
        if (line.isEmpty()) return
        remember("Cellpose: " + line)
        String low = line.toLowerCase(Locale.ROOT)
        if (low.contains("using cpu")) {
            device = "CPU"
        } else if (low.contains("using gpu")) {
            device = "GPU"
        } else if (low.contains("downloading:")) {
            downloading = true; report(15, "Downloading the Cellpose model", "First use needs a large download. See the detection log for details.")
        } else if (low.contains("loading model")) {
            downloading = false; report(30, "Loading the Cellpose model", "Preparing the model in memory" + (device ? " (" + device + ")." : "."))
        } else if (low.contains("network run in")) {
            downloading = false; report(70, "Building cell outlines", "The image analysis has finished; Cellpose is creating the outlines.")
        } else if (low.contains("masks created") || low.contains("saving rois")) {
            downloading = false; report(80, "Saving Cellpose outlines", "Writing the detected shapes to an outlines file for Fiji to read.")
        } else if (low.contains(">>>> completed in")) {
            downloading = false; report(84, "Finishing Cellpose", "Waiting for Cellpose to close its output files.")
        } else {
            def match = (line =~ /(\d{1,3})%\s*\|/)
            if (match.find()) {
                int p = Math.max(0, Math.min(100, Integer.parseInt(match.group(1))))
                if (downloading) report(15 + (int) (10 * p / 100), "Downloading the Cellpose model", "Model download: " + p + "% (this is one part of detection).")
                // The CLI's outer progress bar counts images, not tiles. Its 100% is NOT our 100%.
                else if (stage == "Loading the Cellpose model" && line.contains("0/1")) report(40, "Finding cells in the image", "Cellpose is locating cells. Outlines appear after this stage; intermediate progress is not reported.")
            }
        }
        // Raw log messages belong in the history; retain the useful stage explanation.
        progress.call(percent, stage, detail)
    }

    /** A regular file captures stdout AND stderr. Poll it so silent processes remain cancellable. */
    int run(List command, File directory, File output) {
        checkCancelled()
        logFile = output
        def builder = new ProcessBuilder(command as String[])
        builder.directory(directory); builder.redirectErrorStream(true); builder.redirectOutput(output)
        builder.environment().put("PYTHONUNBUFFERED", "1")
        Process child = builder.start(); process = child
        RandomAccessFile reader = null
        def pending = new StringBuilder()
        try {
            reader = new RandomAccessFile(output, "r")
            boolean finished = false
            while (!finished) {
                checkCancelled()
                finished = child.waitFor(250, java.util.concurrent.TimeUnit.MILLISECONDS)
                // Read a bounded snapshot, including carriage-return-only download updates.
                long end = reader.length()
                while (reader.getFilePointer() < end) {
                    checkCancelled()
                    int b = reader.read()
                    if (b < 0) break
                    if (b == 10 || b == 13) {
                        if (pending.length() > 0) { logLine(pending.toString()); pending.setLength(0) }
                    } else if (pending.length() < 4096) pending.append((char) b)
                }
            }
            if (pending.length() > 0) logLine(pending.toString())
            checkCancelled()
            return child.exitValue()
        } finally {
            if (reader != null) reader.close()
            if (child.isAlive()) {
                child.destroy()
                if (!child.waitFor(2, java.util.concurrent.TimeUnit.SECONDS)) {
                    child.destroyForcibly(); child.waitFor(2, java.util.concurrent.TimeUnit.SECONDS)
                }
            }
            process = null
        }
    }
}

/** One look for every plugin window: the computer's own interface font and the SimpleCa²⁺ colours. */
class CalciumTheme {
    static final Color BLUE = new Color(0x0072B2), BLUE_DARK = new Color(0x004D7A), BLUE_LIGHT = new Color(0xE8F4F8)
    static final Color PAGE = new Color(0xF7F9FB), CARD = Color.white, LINE = new Color(0xDEE4EA), MUTED = new Color(0x65717D), INK = new Color(0x1F2933)
    static final Color GOOD = new Color(0x1E7E34), WARN = new Color(0x9A5B00), BAD = new Color(0xB42318), ACCENT = new Color(0xD55E00), FAINT = new Color(0x9AA5B1)

    /** The look and feel's own label font (San Francisco on a Mac, Segoe UI on Windows), never a fixed family. */
    static java.awt.Font font(int style = java.awt.Font.PLAIN, float delta = 0f) {
        java.awt.Font f = javax.swing.UIManager.getFont("Label.font") ?: new java.awt.Font(java.awt.Font.DIALOG, java.awt.Font.PLAIN, 13)
        float size = (float) Math.max(12f, Math.min(14f, f.getSize2D()))
        return f.deriveFont(style, (float) (size + delta))
    }
    static String hex(Color c) { String.format("#%06x", c.getRGB() & 0xFFFFFF) }
    /** A long name without spaces cannot wrap and would push a window wider than the screen: keep its start and end. */
    static String middle(String text, int most) {
        if (text == null || text.length() <= most) return text ?: ""
        int head = (int) Math.ceil((most - 1) * 0.6d), tail = most - 1 - head
        return text.substring(0, head) + "…" + text.substring(text.length() - tail)
    }
    static void smooth(java.awt.Graphics2D g) {
        g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
        def hints = null
        try { hints = java.awt.Toolkit.getDefaultToolkit().getDesktopProperty("awt.font.desktophints") } catch (Throwable ignored) { }
        if (hints instanceof Map) g.addRenderingHints((Map) hints)
        else g.setRenderingHint(java.awt.RenderingHints.KEY_TEXT_ANTIALIASING, java.awt.RenderingHints.VALUE_TEXT_ANTIALIAS_ON)
    }
    static javax.swing.JButton button(String text, boolean primary = false) { def b = new javax.swing.JButton(text); style(b, primary); b }
    static void style(javax.swing.AbstractButton b, boolean primary) { b.setUI(new CalciumButtonUI(primary: primary)) }
    /** Makes a button the filled blue one (the next thing to press), or an ordinary one again. */
    static void primary(javax.swing.AbstractButton b, boolean on) { if (b.getUI() instanceof CalciumButtonUI && ((CalciumButtonUI) b.getUI()).primary != on) { ((CalciumButtonUI) b.getUI()).primary = on; b.setFont(font(on ? java.awt.Font.BOLD : java.awt.Font.PLAIN)); b.revalidate(); b.repaint() } }
    /** Marks a button as the current option of a group, or not. */
    static void choose(javax.swing.AbstractButton b, boolean on) { if (b.getUI() instanceof CalciumButtonUI) { ((CalciumButtonUI) b.getUI()).chosen = on; b.repaint() } }
    static javax.swing.JLabel label(String text, Color colour = INK, int style = java.awt.Font.PLAIN, float delta = 0f) {
        def l = new javax.swing.JLabel(text); l.setFont(font(style, delta)); l.setForeground(colour); l.setAlignmentX(0f); l
    }
    static void plain(javax.swing.JComponent c) { c.setFont(font()); c.setOpaque(false); c.setAlignmentX(0f) }
    /** Text that fits the given width, ending in an ellipsis when it had to be cut. */
    static String fit(String text, java.awt.FontMetrics fm, int width) {
        if (text == null || width <= 0) return ""
        if (fm.stringWidth(text) <= width) return text
        String s = text
        while (s.length() > 1 && fm.stringWidth(s + "…") > width) s = s.substring(0, s.length() - 1)
        return s.length() > 1 ? s + "…" : ""
    }
    /** Round positions for the ticks of a plot axis. */
    static List<Double> ticks(double lo, double hi, int about) {
        if (!(hi > lo) || !Double.isFinite(lo) || !Double.isFinite(hi)) return [lo]
        double raw = (hi - lo) / Math.max(1, about), power = Math.pow(10d, Math.floor(Math.log10(raw))), unit = raw / power
        double step = (unit <= 1.5d ? 1d : unit <= 3d ? 2d : unit <= 7d ? 5d : 10d) * power
        List<Double> out = []
        for (double v = Math.ceil(lo / step - 1e-9d) * step; v <= hi + step * 1e-9d; v += step) out << (Math.abs(v) < step * 1e-9d ? 0d : v)
        return out
    }
    static String tick(double v, double step) {
        int decimals = step >= 1d ? 0 : (int) Math.min(6, Math.ceil(-Math.log10(step) - 1e-9d))
        return String.format(Locale.US, "%." + decimals + "f", v)
    }
    /** A number for a settings box: at most three decimals, no trailing zeros. The exact value is kept separately. */
    static String shortNumber(double v) {
        if (!Double.isFinite(v)) return ""
        def format = new java.text.DecimalFormat("0.###", java.text.DecimalFormatSymbols.getInstance(Locale.US))
        return format.format(v)
    }
    /** Puts a picture window in the free space to the left of a control window, so neither covers the other. */
    static void beside(ImagePlus picture, java.awt.Window panel) {
        try {
            def window = picture?.getWindow(); if (window == null || panel == null) return
            def configuration = panel.getGraphicsConfiguration(); def bounds = configuration.getBounds()
            def insets = java.awt.Toolkit.getDefaultToolkit().getScreenInsets(configuration)
            int left = (int) bounds.getX() + insets.left + 8, top = (int) bounds.getY() + insets.top + 8
            int freeWidth = panel.getX() - left - 8, freeHeight = (int) bounds.getHeight() - insets.top - insets.bottom - 16
            int chromeWidth = Math.max(0, window.getWidth() - picture.getCanvas().getWidth()), chromeHeight = Math.max(0, window.getHeight() - picture.getCanvas().getHeight())
            if (freeWidth - chromeWidth < 240 || freeHeight - chromeHeight < 180) return
            double zoom = Math.min(1d, Math.min((freeWidth - chromeWidth) / (double) picture.getWidth(), (freeHeight - chromeHeight) / (double) picture.getHeight()))
            window.setLocationAndSize(left, top, (int) (picture.getWidth() * zoom) + chromeWidth, (int) (picture.getHeight() * zoom) + chromeHeight)
        } catch (Throwable ignored) { } // Tidy window placement is a convenience; never stop the plugin for it.
    }
}

/** Flat, rounded buttons that look the same on Mac, Windows and Linux. A primary button is filled blue. */
class CalciumButtonUI extends javax.swing.plaf.basic.BasicButtonUI {
    boolean primary = false, chosen = false   // chosen: the current option of a group, drawn light blue
    void installUI(javax.swing.JComponent c) {
        super.installUI(c)
        def b = (javax.swing.AbstractButton) c
        b.setOpaque(false); b.setContentAreaFilled(false); b.setBorderPainted(false); b.setFocusPainted(false); b.setRolloverEnabled(true)
        b.setMargin(new java.awt.Insets(0, 0, 0, 0)); b.setBorder(javax.swing.BorderFactory.createEmptyBorder(5, 12, 5, 12))
        b.setFont(CalciumTheme.font(primary ? java.awt.Font.BOLD : java.awt.Font.PLAIN))
    }
    void paint(java.awt.Graphics g0, javax.swing.JComponent c) {
        def b = (javax.swing.AbstractButton) c; def m = b.getModel()
        def g = (java.awt.Graphics2D) g0.create()
        try {
            CalciumTheme.smooth(g)
            int w = c.getWidth(), h = c.getHeight()
            Color fill, line, ink
            if (!m.isEnabled()) { fill = primary ? new Color(0xB7D3E3) : CalciumTheme.PAGE; line = primary ? fill : CalciumTheme.LINE; ink = primary ? Color.white : CalciumTheme.FAINT }
            else if (primary) { fill = m.isArmed() && m.isPressed() ? CalciumTheme.BLUE_DARK : m.isRollover() ? new Color(0x0064A0) : CalciumTheme.BLUE; line = fill; ink = Color.white }
            else if (chosen) { fill = new Color(0xD3E6F0); line = CalciumTheme.BLUE; ink = CalciumTheme.BLUE_DARK }
            else { fill = m.isArmed() && m.isPressed() ? new Color(0xD3E6F0) : m.isRollover() ? CalciumTheme.BLUE_LIGHT : Color.white; line = m.isRollover() ? CalciumTheme.BLUE : new Color(0xC5CED8); ink = CalciumTheme.INK }
            g.setColor(fill); g.fillRoundRect(0, 0, w - 1, h - 1, 10, 10)
            g.setColor(line); g.drawRoundRect(0, 0, w - 1, h - 1, 10, 10)
            if (b.isFocusOwner() && m.isEnabled()) { g.setColor(new Color(0, 114, 178, 110)); g.setStroke(new java.awt.BasicStroke(2f)); g.drawRoundRect(1, 1, w - 3, h - 3, 9, 9) }
            g.setFont(b.getFont()); def fm = g.getFontMetrics(); String text = b.getText() ?: ""
            g.setColor(ink); g.drawString(text, (int) ((w - fm.stringWidth(text)) / 2), (int) ((h - fm.getHeight()) / 2) + fm.getAscent())
        } finally { g.dispose() }
    }
}

/** The clickable title line of a card: number, title, a short summary and an open/closed arrow. */
class CalciumCardHeader extends javax.swing.JComponent {
    CalciumCard card; boolean hover = false
    CalciumCardHeader(CalciumCard owner) {
        card = owner; setAlignmentX(0f); setOpaque(false)
        addMouseListener(new java.awt.event.MouseAdapter() {
            void mouseEntered(java.awt.event.MouseEvent e) { hover = true; repaint() }
            void mouseExited(java.awt.event.MouseEvent e) { hover = false; repaint() }
            void mouseReleased(java.awt.event.MouseEvent e) { if (card.collapsible && javax.swing.SwingUtilities.isLeftMouseButton(e) && contains(e.getPoint())) card.onToggle.call(card) }
        })
    }
    java.awt.Dimension getPreferredSize() { new java.awt.Dimension(220, 38) }
    java.awt.Dimension getMinimumSize() { new java.awt.Dimension(120, 38) }
    java.awt.Dimension getMaximumSize() { new java.awt.Dimension(Integer.MAX_VALUE, 38) }
    protected void paintComponent(java.awt.Graphics g0) {
        def g = (java.awt.Graphics2D) g0.create()
        try {
            CalciumTheme.smooth(g)
            int w = getWidth(), h = getHeight(), x = 12
            if (hover && card.collapsible) { g.setColor(CalciumTheme.BLUE_LIGHT); g.fillRoundRect(2, 2, w - 5, h - 4, 10, 10) }
            if (card.badge) {
                g.setColor(card.done ? CalciumTheme.GOOD : card.expanded ? CalciumTheme.BLUE : CalciumTheme.FAINT); g.fillOval(x, (int) ((h - 22) / 2), 22, 22)
                g.setColor(Color.white)
                if (card.done) { g.setStroke(new java.awt.BasicStroke(2f, java.awt.BasicStroke.CAP_ROUND, java.awt.BasicStroke.JOIN_ROUND)); int cy = (int) (h / 2); g.drawPolyline([x + 6, x + 10, x + 16] as int[], [cy, cy + 4, cy - 4] as int[], 3) }
                else { g.setFont(CalciumTheme.font(java.awt.Font.BOLD, -1f)); def fm = g.getFontMetrics(); g.drawString(card.badge, x + (int) ((22 - fm.stringWidth(card.badge)) / 2), (int) ((h - fm.getHeight()) / 2) + fm.getAscent()) }
                x += 30
            }
            g.setFont(CalciumTheme.font(java.awt.Font.BOLD, 1f)); def fm = g.getFontMetrics(); int baseline = (int) ((h - fm.getHeight()) / 2) + fm.getAscent()
            g.setColor(CalciumTheme.INK); g.drawString(card.title, x, baseline); x += fm.stringWidth(card.title) + 12
            int right = w - 12
            if (card.collapsible) {
                g.setColor(CalciumTheme.MUTED); g.setStroke(new java.awt.BasicStroke(1.6f, java.awt.BasicStroke.CAP_ROUND, java.awt.BasicStroke.JOIN_ROUND)); int cx = w - 20, cy = (int) (h / 2)
                if (card.expanded) g.drawPolyline([cx - 4, cx, cx + 4] as int[], [cy - 2, cy + 2, cy - 2] as int[], 3) else g.drawPolyline([cx - 2, cx + 2, cx - 2] as int[], [cy - 4, cy, cy + 4] as int[], 3)
                right = w - 32
            }
            if (card.summary && !(card.collapsible && card.expanded && !card.summaryWhenOpen)) {
                g.setFont(CalciumTheme.font()); fm = g.getFontMetrics(); String s = CalciumTheme.fit(card.summary, fm, right - x)
                g.setColor(CalciumTheme.MUTED); g.drawString(s, right - fm.stringWidth(s), baseline)
            }
        } finally { g.dispose() }
    }
}

/** A white card with a numbered title. Clicking the title opens or closes what is inside. */
class CalciumCard extends javax.swing.JPanel {
    String badge = "", title = "", summary = ""
    boolean collapsible = true, expanded = true, done = false, summaryWhenOpen = false
    javax.swing.JPanel body = new javax.swing.JPanel()
    CalciumCardHeader header
    Closure onToggle = { CalciumCard card -> card.setExpanded(!card.expanded) }
    CalciumCard(String badge, String title, boolean collapsible = true) {
        this.badge = badge ?: ""; this.title = title; this.collapsible = collapsible
        setLayout(new java.awt.BorderLayout()); setOpaque(false); setAlignmentX(0f)
        setBorder(javax.swing.BorderFactory.createEmptyBorder(1, 1, 1, 1))
        body.setLayout(new javax.swing.BoxLayout(body, javax.swing.BoxLayout.Y_AXIS)); body.setOpaque(false)
        body.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 12, 10, 12))
        header = new CalciumCardHeader(this); header.setToolTipText(collapsible ? "Click to open or close this step" : null)
        add(header, java.awt.BorderLayout.NORTH); add(body, java.awt.BorderLayout.CENTER)
    }
    void setExpanded(boolean value) { expanded = value; body.setVisible(value); revalidate(); repaint() }
    void setSummary(String value) { summary = value ?: ""; header?.repaint() }
    void setDone(boolean value) { done = value; header?.repaint() }
    java.awt.Dimension getMaximumSize() { new java.awt.Dimension(Integer.MAX_VALUE, (int) getPreferredSize().height) }
    protected void paintComponent(java.awt.Graphics g0) {
        def g = (java.awt.Graphics2D) g0.create()
        try {
            CalciumTheme.smooth(g)
            g.setColor(CalciumTheme.CARD); g.fillRoundRect(0, 0, getWidth() - 1, getHeight() - 1, 12, 12)
            g.setColor(CalciumTheme.LINE); g.drawRoundRect(0, 0, getWidth() - 1, getHeight() - 1, 12, 12)
        } finally { g.dispose() }
    }
}

/** One compact Step 2 panel, updated only on Swing's event thread: a slim bar, the stage in plain words, and buttons for the details. */
class CalciumOutlineProgress extends javax.swing.JPanel {
    javax.swing.JProgressBar bar = new javax.swing.JProgressBar(0, 100)
    javax.swing.JLabel text = new javax.swing.JLabel()
    javax.swing.JButton cancelButton = new javax.swing.JButton("Cancel detection")
    javax.swing.JButton logButton = new javax.swing.JButton("View detection log")
    javax.swing.JButton detailsButton = new javax.swing.JButton("Progress details")
    javax.swing.JTextArea detailsText = new javax.swing.JTextArea(18, 64)
    javax.swing.JDialog detailsWindow
    javax.swing.Timer clock
    java.util.concurrent.atomic.AtomicReference latest = new java.util.concurrent.atomic.AtomicReference()
    CalciumOutlineTask task
    long started, stopped
    boolean running = false

    static String escape(String value) { value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;") }
    CalciumOutlineProgress() {
        setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.Y_AXIS)); setAlignmentX(0f)
        bar.setString("0% - Stage progress"); bar.setAlignmentX(0f)   // the percentage is written next to the stage name below, not over the bar
        bar.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, 14))
        bar.setToolTipText("Progress through detection stages, not a time estimate. It can pause while Cellpose runs.")
        setOpaque(false)
        text.setAlignmentX(0f); text.setFont(CalciumTheme.font()); text.setForeground(CalciumTheme.INK); text.setVerticalAlignment(javax.swing.SwingConstants.TOP)
        // Reserve space for a quiet/long-running message without growing beyond the screen mid-run.
        text.setPreferredSize(new java.awt.Dimension(330, 118))
        text.setMinimumSize(new java.awt.Dimension(330, 118))
        text.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, 118))
        text.setBorder(javax.swing.BorderFactory.createEmptyBorder(6, 0, 0, 0))
        def buttons = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 4, 2)); buttons.setAlignmentX(0f); buttons.setOpaque(false)
        cancelButton.setText("Cancel"); detailsButton.setText("Details"); logButton.setText("Log")
        cancelButton.setToolTipText("Stop the detection"); detailsButton.setToolTipText("Everything the detection has reported so far"); logButton.setToolTipText("Open Cellpose's own messages")
        [cancelButton, detailsButton, logButton].each { CalciumTheme.style(it, false); buttons.add(it) }
        add(bar); add(text); add(buttons)
        detailsText.setEditable(false); detailsText.setLineWrap(true); detailsText.setWrapStyleWord(true)
        detailsText.setFont(text.getFont()); detailsText.setMargin(new java.awt.Insets(8, 8, 8, 8))
        ((javax.swing.text.DefaultCaret) detailsText.getCaret()).setUpdatePolicy(javax.swing.text.DefaultCaret.NEVER_UPDATE)
        detailsButton.addActionListener({ e ->
            if (detailsWindow == null) {
                detailsWindow = new javax.swing.JDialog(javax.swing.SwingUtilities.getWindowAncestor(this), "Detection progress details", java.awt.Dialog.ModalityType.MODELESS)
                detailsWindow.add(new javax.swing.JScrollPane(detailsText)); detailsWindow.pack()
                detailsWindow.setLocationRelativeTo(this)
            }
            render(); detailsWindow.setVisible(true)
        } as java.awt.event.ActionListener)
        cancelButton.addActionListener({ e ->
            if (task != null && task.cancellable) task.cancel()
            cancelButton.setEnabled(false); cancelButton.setText("Cancelling...")
        } as java.awt.event.ActionListener)
        logButton.addActionListener({ e ->
            File log = task?.logFile
            if (log != null && log.isFile()) {
                try { java.awt.Desktop.getDesktop().open(log) }
                catch (Throwable t) { IJ.showMessage("Detection log", "Open this file to see Cellpose's messages:\n" + log.getPath()) }
            }
        } as java.awt.event.ActionListener)
        clock = new javax.swing.Timer(500, { e -> render() } as java.awt.event.ActionListener)
        setVisible(false)
    }
    CalciumOutlineTask begin() {
        started = System.nanoTime(); stopped = 0; running = true
        task = new CalciumOutlineTask(progress: { int pct, String stage, String detail ->
            latest.set(task.snapshot())
        })
        task.began = started
        task.report(0, "Preparing cell detection", "Checking Cellpose and looking for saved outlines.")
        bar.setValue(0); cancelButton.setText("Cancel"); cancelButton.setEnabled(true)
        bar.setVisible(true); cancelButton.setVisible(true); roomFor(118)
        setVisible(true); render(); clock.start()
        return task
    }
    void render() {
        def state = latest.get(); if (state == null) return
        long now = running ? System.nanoTime() : stopped
        long seconds = (long) ((now - started) / 1000000000L)
        long quiet = (long) ((now - (long) state.at) / 1000000000L)
        long stageSeconds = (long) ((now - (long) (state.stageAt ?: now)) / 1000000000L)
        bar.setValue((int) state.percent)
        bar.setString(state.percent + "% - Stage progress")
        String elapsed = String.format(Locale.US, "%d:%02d total | %d:%02d in this stage", (long) (seconds / 60), seconds % 60, (long) (stageSeconds / 60), stageSeconds % 60)
        String message = state.detail
        String quietMessage = ""
        if (running && quiet >= 60) quietMessage = "No new progress message for " + quiet + " s. " +
            (task.process?.isAlive() ? "Cellpose is still running; no measurable progress reported." : "Waiting for the current stage to report progress.")
        if (running && seconds >= 600) quietMessage += " Over 10 minutes so far; check the log."
        String context = state.context ?: ""
        String muted = CalciumTheme.hex(CalciumTheme.MUTED)
        text.setText("<html><body style='width:300px'><b>" + escape(state.stage) + "</b> &nbsp;<font color='" + muted + "'>" + state.percent + " %</font><br><font color='" + muted + "'>" + escape(elapsed) +
            (context ? " &nbsp;·&nbsp; " + escape(context) : "") + "</font><br>" + escape(message) +
            (quietMessage ? "<br><font color='#865000'>" + escape(quietMessage) + "</font>" : "") + "</body></html>")
        String details = state.stage + " (" + state.percent + "%)\n" + elapsed + "\n" + context + "\n" + message + "\n" + quietMessage +
            "\n\nPercentages mark completed stages, not remaining time. Cellpose does not report a reliable percentage while finding cells in one image.\n\nRecent events (up to 120; full output is in View detection log):\n" + (state.history ?: "")
        if (detailsText.getText() != details) {
            int position = detailsText.getCaretPosition()
            detailsText.setText(details); detailsText.setCaretPosition(Math.min(position, details.length()))
        }
        logButton.setEnabled(task?.logFile?.isFile() ?: false)
        // Before Cellpose starts there may be an installation or a modal file chooser.
        cancelButton.setEnabled(running && task.cancellable && !task.cancelled)
    }
    void finish(boolean ok, String message) {
        stopped = System.nanoTime(); running = false; clock.stop()
        String finalStage = ok ? "Cell outlines ready" : task.cancelled ? "Detection cancelled" : "Detection did not finish"
        task.remember(finalStage + ": " + message)
        latest.set(task.snapshot() + [percent: ok ? 100 : task.percent, stage: finalStage, detail: message, stageAt: stopped])
        render(); cancelButton.setEnabled(false)
        if (!ok) bar.setString(task.percent + "% - " + (task.cancelled ? "Cancelled" : "Stopped"))
        else {
            // Done: only the result stays (the stages, timings and settings are under Details), so the next step is not pushed away.
            bar.setVisible(false); cancelButton.setVisible(false); roomFor(40)
            text.setText("<html><body style='width:300px'>" + escape(message) + "</body></html>")
        }
    }
    void roomFor(int height) {
        text.setPreferredSize(new java.awt.Dimension(330, height)); text.setMinimumSize(new java.awt.Dimension(330, height)); text.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, height))
        revalidate()
    }
    void close() { clock.stop(); if (running && task != null) task.cancel(); detailsWindow?.dispose() }
}

/** A display failure must not prevent detection or leave a running clock without a worker. */
class CalciumDetectionRunner {
    Closure launch = { Closure job -> Thread.start(job) }
    void safelyRefresh(Closure refreshWindow) {
        try { refreshWindow.call() } catch (Throwable problem) { IJ.log("Calcium window refresh failed: " + problem) }
    }
    void complete(javax.swing.JButton button, CalciumOutlineProgress panel, String original, boolean ok, String message, Closure after, Closure refreshWindow) {
        try { panel.finish(ok, message) }
        catch (Throwable problem) { IJ.log("Calcium progress display failed: " + problem) }
        finally {
            panel.clock.stop(); panel.running = false
            button.setEnabled(true); button.setText(original)
        }
        try { after.call(message) } catch (Throwable problem) { IJ.log("Calcium result display failed: " + problem) }
        safelyRefresh(refreshWindow)
    }
    void start(javax.swing.JButton button, CalciumOutlineProgress panel, Closure detect, Closure after, Closure refreshWindow) {
        String original = button.getText()
        button.setEnabled(false); button.setText("Looking for cell outlines...")
        try {
            def task = panel.begin()
            // Start the worker before doing any cosmetic window refresh.
            launch.call({ ->
                String message = "Detection failed"
                try { message = detect.call(task) }
                catch (java.util.concurrent.CancellationException ignored) { message = "Cancelled. You can draw cells by hand or try detection again." }
                catch (Throwable problem) { message = "Failed: " + problem; IJ.log("Calcium imaging plugin: " + problem) }
                finally {
                    final String result = message
                    javax.swing.SwingUtilities.invokeLater { complete(button, panel, original, task.succeeded, result, after, refreshWindow) }
                }
            })
        } catch (Throwable problem) {
            complete(button, panel, original, false, "Could not start detection: " + problem, after, refreshWindow)
            return
        }
        safelyRefresh(refreshWindow)
    }
}

/** A cancellable local job. Processes and their descendants are stopped together. */
class CalciumLocalJob {
    volatile boolean cancelled = false
    volatile Process process
    File logFile
    Closure progress = { int p, String message -> }
    void check() { if (cancelled) throw new java.util.concurrent.CancellationException('Cancelled') }
    void cancel() {
        cancelled = true
        def p = process
        if (p != null) {
            try { p.toHandle().descendants().forEach { it.destroyForcibly() } } catch (Throwable ignored) { }
            p.destroyForcibly()
        }
    }
    /** R's last error message in a log ("Error in f(x) :" followed by the message, or "Error: message"), or null. */
    static String engineError(File log) {
        try {
            if (log == null || !log.isFile()) return null
            // Windows line endings leave a blank line after every line of the log, so blank lines are skipped
            List<String> lines = log.readLines('UTF-8').findAll { it.trim() }
            int at = lines.findLastIndexOf { it ==~ /^Error( in .*)?\s*:.*/ }
            if (at < 0) return null
            String first = lines[at].replaceFirst(/^Error( in .*?)?\s*:\s*/, '').trim()
            List<String> rest = lines.drop(at + 1).takeWhile { it.startsWith(' ') && !it.trim().startsWith('Calls:') }.collect { it.trim() }
            String message = ([first] + rest).findAll { it }.join(' ').trim()
            return message ? message.take(400) : null
        } catch (Throwable ignored) { return null }
    }
    void run(List<String> command, File directory, int timeoutSeconds = 900) {
        check(); logFile.parentFile.mkdirs()
        def pb = new ProcessBuilder(command as String[]).directory(directory).redirectErrorStream(true)
        // An unrelated R installation or project .Rprofile must not alter this analysis.
        ['R_HOME','R_LIBS','R_LIBS_USER','R_LIBS_SITE','R_PROFILE','R_ENVIRON'].each { pb.environment().remove(it) }
        Process p = pb.start(); process = p
        def readerError = new java.util.concurrent.atomic.AtomicReference<Throwable>()
        def reader = Thread.startDaemon('SimpleCa output') {
            try {
                logFile.withWriterAppend('UTF-8') { log ->
                    p.inputStream.withReader('UTF-8') { input ->
                        String line
                        while ((line = input.readLine()) != null) {
                            log.write(line.take(8192) + '\n'); log.flush()
                            if (line.startsWith('SIMPLECA\t')) {
                                def fields = line.split('\t',3)
                                if (fields.length == 3) progress.call(Integer.parseInt(fields[1]),fields[2])
                            }
                        }
                    }
                }
            } catch (Throwable problem) { readerError.set(problem) }
        }
        long deadline = System.nanoTime() + timeoutSeconds * 1000000000L
        try {
            while (!p.waitFor(200, java.util.concurrent.TimeUnit.MILLISECONDS)) {
                check()
                if (System.nanoTime() > deadline) throw new IOException('This operation timed out. See the analysis log and retry.')
            }
            reader.join(3000); check()
            if (readerError.get() != null) throw new IOException('Could not read analysis progress',readerError.get())
            if (p.exitValue() != 0) { String said = engineError(logFile); throw new IOException(said ? 'The analysis engine stopped: ' + said : 'Analysis engine returned exit code '+p.exitValue()+'. See '+logFile.name+'.') }
        } finally {
            if (p.isAlive()) { try { p.toHandle().descendants().forEach { it.destroyForcibly() } } catch (Throwable ignored) { }; p.destroyForcibly() }
            process = null
        }
    }
}

/** Per-user R runtime. No system installation, shell initialization or changes to Cellpose. */
class CalciumRRuntime {
    File home = new File(System.getProperty('user.home'), '.calcium-simpleca')
    static final String R_VERSION = '4.5.2'
    static final String JSONLITE_VERSION = '2.0.0'
    static final String MAMBA_VERSION = '2.10.0-0'
    static final Map HASHES = [
        'linux-64':'701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd',
        'osx-64':'8421462a1ef1df9f7acd114f9e2a24c1332ce4f4e29f44465b4acb44e2755172',
        'osx-arm64':'c547639f5caf742236641c362746abb4a3505774753ecd08fe341643670c4c01',
        'win-64':'f8b4284bd4eb1134122eaf93b3ea6c7be567e596de5829e1a6f478e34e41974a']
    static String platform(String os = System.getProperty('os.name'), String arch = System.getProperty('os.arch')) {
        boolean arm = arch.toLowerCase(Locale.ROOT) in ['aarch64','arm64']
        boolean x64 = arch.toLowerCase(Locale.ROOT) in ['amd64','x86_64']
        if (os.toLowerCase(Locale.ROOT).contains('mac') && (arm || x64)) return arm ? 'osx-arm64' : 'osx-64'
        if (os.toLowerCase(Locale.ROOT).contains('win') && x64) return 'win-64'
        if (os.toLowerCase(Locale.ROOT).contains('linux') && x64) return 'linux-64'
        throw new IOException('Automatic R setup supports 64-bit Windows, Intel/Apple Silicon Mac, and x64 Linux.')
    }
    File executable() { new File(home,'bin/'+platform()+'/micromamba'+(platform().startsWith('win')?'.exe':'')) }
    File environment() {
        File marker = new File(home,'ready-protocol-v2-'+platform()+'.txt')
        if (!marker.isFile()) return null
        String name = marker.getText('UTF-8').trim()
        if (!(name ==~ /r-4\.5\.2-jsonlite-2\.0\.0-[a-f0-9-]+/)) return null
        File env = new File(home,name)
        executable().isFile() && new File(env,'conda-meta/history').isFile() ? env : null
    }
    List<String> command(File env) {
        def base=[executable().path,'--no-rc','-r',new File(home,'packages').path,'run','-p',env.path]
        if(platform().startsWith('win'))return base+['Rscript','--vanilla']
        File rhome=new File(env,'lib/R')
        // Conda's Unix R shell wrapper embeds unquoted prefixes. The native entry point supports spaces.
        base+['-e','R_HOME='+rhome.path,'-e','R_SHARE_DIR='+new File(rhome,'share').path,
            new File(rhome,'bin/exec/R').path,'--vanilla','--slave']
    }
    List<String> scriptCommand(File env,File script,List<String> args=[]) {
        return command(env) + (platform().startsWith('win') ? [script.path] : ['--file='+script.path,'--args']) + args
    }
    void verifyExecutable() {
        if (!executable().isFile() || CalciumInstaller.digest(executable().bytes) != HASHES[platform()]) throw new IOException('The R launcher failed verification. Run setup again.')
    }
    void probe(File env,CalciumLocalJob job) {
        // Windows launchers can strip quotes inside inline -e expressions. Use the same file path as analysis.
        File script=File.createTempFile('verify-r-', '.R', home)
        try {
            script.setText('stopifnot(as.character(getRversion()) == "'+R_VERSION+'", requireNamespace("jsonlite", quietly=TRUE), as.character(packageVersion("jsonlite")) == "'+JSONLITE_VERSION+'")\ncat("R ready\\n")\n','UTF-8')
            job.run(scriptCommand(env,script),home,60)
        } finally {script.delete()}
    }
    File setup(CalciumLocalJob job) {
        home.mkdirs(); job.check()
        def channel = java.nio.channels.FileChannel.open(new File(home,'setup.lock').toPath(), java.nio.file.StandardOpenOption.CREATE,java.nio.file.StandardOpenOption.WRITE)
        def lock = null
        try {
            try { lock = channel.tryLock() } catch (java.nio.channels.OverlappingFileLockException ignored) { }
            if (lock == null) throw new IOException('R setup is already running in another window.')
            File ready = environment()
            if (ready != null) {
                try { verifyExecutable(); probe(ready,job); job.progress(100,'R is ready'); return ready }
                catch (java.util.concurrent.CancellationException cancelled) { throw cancelled }
                catch (Exception broken) { new File(home,'ready-protocol-v2-'+platform()+'.txt').delete() }
            }
            String platform = platform()
            job.progress(5,'Downloading the verified R installer')
            int lastDownload = -1
            byte[] binary = CalciumInstaller.fetch('https://github.com/mamba-org/micromamba-releases/releases/download/'+MAMBA_VERSION+'/micromamba-'+platform,32*1024*1024,{ n,total ->
                job.check(); int percent=total>0?(int)(100L*n/total):0; if(percent!=lastDownload){lastDownload=percent; job.progress(5+(int)(percent/5),'Downloading R setup: '+percent+'%')}
            })
            if (CalciumInstaller.digest(binary) != HASHES[platform]) throw new IOException('R installer checksum verification failed')
            executable().parentFile.mkdirs(); File temp = File.createTempFile('mamba-', '.tmp', executable().parentFile); temp.bytes=binary
            try { java.nio.file.Files.move(temp.toPath(),executable().toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING) } finally { temp.delete() }
            if (!platform.startsWith('win') && !executable().setExecutable(true)) throw new IOException('Could not enable the R installer')
            job.check(); verifyExecutable()
            File env = new File(home,'r-'+R_VERSION+'-jsonlite-'+JSONLITE_VERSION+'-'+UUID.randomUUID())
            job.progress(30,'Installing R (first run only). Download details are in the setup log.')
            job.run([executable().path,'--no-rc','-r',new File(home,'packages').path,'create','--yes','--override-channels','--strict-channel-priority','-c','https://conda.anaconda.org/conda-forge','-p',env.path,'r-base='+R_VERSION,'r-jsonlite='+JSONLITE_VERSION],home,1800)
            job.progress(90,'Checking the installed analysis engine')
            probe(env,job)
            job.check()
            File markerTemp=File.createTempFile('ready-', '.tmp', home); markerTemp.setText(env.name,'UTF-8')
            java.nio.file.Files.move(markerTemp.toPath(),new File(home,'ready-protocol-v2-'+platform+'.txt').toPath(),java.nio.file.StandardCopyOption.ATOMIC_MOVE,java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            job.progress(100,'R is ready. You can now analyze responses.'); env
        } finally { if(lock!=null) lock.release(); channel.close() }
    }
}

/** Engine files travel inside the verified plugin script, including when updated by older installers. */
class CalciumSimpleCaEngine {
    static File unpack(File home) {
        byte[] bytes = Base64.getDecoder().decode(CalciumSimpleCaBundle.DATA)
        if(CalciumInstaller.digest(bytes)!=CalciumSimpleCaBundle.HASH) throw new IOException('Bundled analysis code failed verification')
        File dir = new File(home,'engine/'+CalciumSimpleCaBundle.HASH); dir.mkdirs()
        def zip = new java.util.zip.ZipInputStream(new ByteArrayInputStream(bytes))
        try {
            java.util.zip.ZipEntry entry
            while((entry=zip.nextEntry)!=null) {
                File target = new File(dir,entry.name)
                if(!target.canonicalPath.startsWith(dir.canonicalPath+File.separator) || entry.isDirectory()) throw new IOException('Invalid engine bundle path')
                def out = new ByteArrayOutputStream(); byte[] buffer=new byte[8192]; int count
                while((count=zip.read(buffer))!=-1) { if(out.size()+count>1024*1024) throw new IOException('Engine file too large'); out.write(buffer,0,count) }
                byte[] content=out.toByteArray(); target.parentFile.mkdirs()
                if(!target.isFile() || CalciumInstaller.digest(target.bytes)!=CalciumInstaller.digest(content)) {
                    File temp=File.createTempFile('engine-', '.tmp',target.parentFile); temp.bytes=content
                    java.nio.file.Files.move(temp.toPath(),target.toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING)
                }
            }
        } finally { zip.close() }
        dir
    }
    static File analyze(CalciumRRuntime runtime, File input, Map settings, File outputRoot, CalciumLocalJob job) {
        File env=runtime.environment()
        if(env==null) throw new IOException('Set up the analysis engine first.')
        runtime.verifyExecutable(); job.check()
        File engine=unpack(runtime.home)
        File output=new File(outputRoot,'SimpleCa-'+new java.text.SimpleDateFormat('yyyyMMdd-HHmmss').format(new Date())+'-'+UUID.randomUUID().toString().take(8))
        output.mkdirs(); job.logFile=new File(output,'analysis.log')
        // Snapshot the exact measurement file; later measurements cannot change this run.
        File snapshot=new File(output,'input_corrected_fluorescence.csv'); java.nio.file.Files.copy(input.toPath(),snapshot.toPath())
        settings=new LinkedHashMap(settings)
        settings.plugin_version=CalciumRelease.VERSION; settings.engine_sha256=CalciumSimpleCaBundle.HASH
        settings.input_sha256=CalciumInstaller.digest(snapshot.bytes)
        if(settings.expected_input_sha256 && settings.input_sha256!=settings.expected_input_sha256) throw new IOException('The measurement file changed. Close this review window, measure again and reopen Step 4.')
        settings.measurement_input_sha256=settings.input_sha256
        if(settings.timestamps_file) {
            File times=new File(settings.timestamps_file.toString())
            File savedTimes=new File(output,'timestamps.csv');java.nio.file.Files.copy(times.toPath(),savedTimes.toPath())
            if(settings.timestamps_sha256 && CalciumInstaller.digest(savedTimes.bytes)!=settings.timestamps_sha256)throw new IOException('The timestamps changed. Load and confirm the timing again.')
            double[] clock=CalciumMeasurement.readTimes(savedTimes,CalciumCsv.read(snapshot).size()-1)
            File original=new File(output,'measured_fluorescence.csv');java.nio.file.Files.copy(snapshot.toPath(),original.toPath())
            CalciumMeasurement.attachTimes(original,snapshot,clock)
            settings.timestamps_file='timestamps.csv';settings.timestamps_sha256=CalciumInstaller.digest(savedTimes.bytes)
            settings.time_source='explicit_timestamps';settings.input_sha256=CalciumInstaller.digest(snapshot.bytes)
        } else settings.time_source='constant_fps'
        if(settings.roi_registry_file) {
            File source=new File(settings.roi_registry_file.toString()),target=new File(output,'roi_registry.csv')
            java.nio.file.Files.copy(source.toPath(),target.toPath())
            if(settings.roi_registry_sha256 && CalciumInstaller.digest(target.bytes)!=settings.roi_registry_sha256)throw new IOException('The ROI registry changed. Measure and save again.')
            settings.roi_registry_file='roi_registry.csv';settings.roi_registry_sha256=CalciumInstaller.digest(target.bytes)
        }
        if(settings.roi_set_file) {
            File source=new File(settings.roi_set_file.toString()),target=new File(output,'accepted_rois.zip')
            java.nio.file.Files.copy(source.toPath(),target.toPath())
            if(settings.roi_set_sha256 && CalciumInstaller.digest(target.bytes)!=settings.roi_set_sha256)throw new IOException('The accepted ROI snapshot changed. Measure and save again.')
            settings.roi_set_file='accepted_rois.zip';settings.roi_set_sha256=CalciumInstaller.digest(target.bytes)
        }
        if(settings.protocol_json) {
            String json=settings.remove('protocol_json').toString()
            new File(output,'protocol_input.json').setText(json,'UTF-8')
            settings.protocol_file='protocol_input.json'
            settings.protocol_input_sha256=CalciumInstaller.digest(json.getBytes('UTF-8'))
        }
        ['sensitivity_file':'sensitivity_grid.csv','negative_control_file':'negative_control_windows.csv'].each{key,name->
            if(settings[key]) {
                File source=new File(settings[key].toString()),target=new File(output,name)
                java.nio.file.Files.copy(source.toPath(),target.toPath())
                String hash=CalciumInstaller.digest(target.bytes)
                if(settings[key+'_sha256']&&settings[key+'_sha256']!=hash)throw new IOException('The selected diagnostic file changed. Load it again before analysis.')
                settings[key]=name;settings[key+'_sha256']=hash
            }
        }
        File config=new File(output,'settings.tsv')
        config.withWriter('UTF-8') { writer -> settings.each { k,v ->
            if((''+k+v).contains('\n') || (''+k+v).contains('\r') || (''+k+v).contains('\t')) throw new IOException('Setting values must fit on one line')
            writer.write(k+'\t'+v+'\n')
        } }
        job.run(runtime.scriptCommand(env,new File(engine,'runner.R'),[engine.path,snapshot.path,config.path,output.path]),output)
        job.check()
        if(!new File(output,'COMPLETE').isFile()) throw new IOException('Analysis did not finish writing its outputs; see the log')
        output
    }
}

class CalciumCsv {
    static String encode(Object value) {
        String text=value==null?'':value.toString()
        text.contains(',')||text.contains('"')||text.contains('\n')||text.contains('\r') ? '"'+text.replace('"','""')+'"' : text
    }
    static void write(File file,List rows) {
        file.withWriter('UTF-8'){writer->rows.each{row->writer.write(row.collect{encode(it)}.join(',')+'\n')}}
    }
    static List<List<String>> read(File file) {parse(file.getText('UTF-8'))}
    static List<List<String>> parse(String data) {
        def rows=[]; def row=[]; def cell=new StringBuilder(); boolean quoted=false
        for(int i=0;i<data.length();i++) {
            char ch=data.charAt(i)
            if(ch=='"') {
                if(quoted && i+1<data.length() && data.charAt(i+1)=='"') {cell.append('"'); i++} else quoted=!quoted
            } else if(!quoted && (ch==',' || ch=='\n')) {
                row << cell.toString(); cell.setLength(0)
                if(ch=='\n') {rows<<row; row=[]}
            } else if(ch!='\r' || quoted) cell.append(ch)
        }
        if(quoted) throw new IOException('Incomplete CSV')
        if(cell.length()>0 || !row.isEmpty()) {row<<cell.toString(); rows<<row}
        if(rows.size()>1 && rows.any{it.size()!=rows[0].size()}) throw new IOException('Inconsistent CSV columns')
        rows
    }
    static double number(String s) { try { s ? Double.parseDouble(s) : Double.NaN } catch(Exception ignored) {Double.NaN} }
}

/** Measurement transport preserves doubles; display/export rounding never enters native analysis. */
class CalciumMeasurement {
    static String fullPrecision(double[][] raw,double[][] background,String[] names) {
        def out=new StringBuilder('Frame');names.each{out.append(',').append(CalciumCsv.encode(it))};out.append('\n')
        for(int f=0;f<raw.length;f++) {
            out.append(f+1)
            for(int c=0;c<names.length;c++) {
                double value=raw[f][c]-background[f][c]
                out.append(',');if(Double.isFinite(value))out.append(Double.toString(value))
            }
            out.append('\n')
        }
        out.toString()
    }
    static String idFor(Roi roi) {
        String id=roi.getProperty('calcium.roi.id')
        if(!id || !(id==~/Cell_[a-f0-9]{32}/)) {id=assignId(roi)}
        id
    }
    /** A new outline's ID is worked out from the outline itself (where it is and which pixels it covers), so the same outline
     *  gets the same ID in every session and tables from two sessions can be matched by cell. The ID then stays with the outline,
     *  also when it is reshaped. A second, identical outline in one set gets a random ID. */
    static String assignId(Roi roi,boolean random=false) {
        String id=null
        if(!random){
            try{
                java.security.MessageDigest md=java.security.MessageDigest.getInstance('SHA-256')
                java.awt.Rectangle b=roi.getBounds();md.update((b.x+','+b.y+','+b.width+','+b.height+';').getBytes('UTF-8'))
                def mask=roi.getMask();if(mask!=null)md.update((byte[])mask.getPixels())
                id='Cell_'+md.digest().encodeHex().toString().substring(0,32)
            }catch(Exception ignored){id=null}
        }
        if(id==null)id='Cell_'+UUID.randomUUID().toString().replace('-','')
        roi.setProperty('calcium.roi.id',id);id
    }
    static void ensureIds(Roi[] rois) {
        Set<String> seen=new HashSet<>()
        rois.eachWithIndex {roi,i->
            String id=idFor(roi)
            if(seen.contains(id)) {
                // Roi.clone shares mutable properties. A serialization round trip detaches them.
                Roi copy=new ij.io.RoiDecoder(RoiEncoder.saveAsByteArray(roi),'cell.roi').getRoi()
                id=assignId(copy,true);rois[i]=copy
            }
            seen.add(id)
        }
    }
    static List registry(Roi[] accepted,Roi[] measured,int[] backgroundPixels,CalciumSteps steps,String selection,Map selectionByCell=null) {
        def measuredIds=measured.collect{idFor(it)}
        def rows=[['Cell','ROI_ID','Accepted','Included_in_measurement','Exclusion_reason','ROI_sha256','Area_pixels','Background_pixels','Near_edge','Selection_basis']]
        accepted.each {roi->
            String id=idFor(roi);int index=measuredIds.indexOf(id);boolean included=index>=0
            Roi canonical=(Roi)roi.clone();canonical.setName(id);canonical.setPosition(0)
            rows<<[id,id,true,included,included?'':'Edge exclusion policy',CalciumInstaller.digest(RoiEncoder.saveAsByteArray(canonical)),roi.getContainedPoints().length,
                included?backgroundPixels[index]:'',steps.nearEdge(roi),selectionByCell?.get(id)?:selection]
        }
        rows
    }
    static void writeRois(File file,Roi[] rois) {
        def zip=new ZipOutputStream(new FileOutputStream(file))
        try{rois.each{roi->Roi copy=(Roi)roi.clone();String id=idFor(copy);copy.setName(id);copy.setPosition(0);zip.putNextEntry(new ZipEntry(id+'.roi'));zip.write(RoiEncoder.saveAsByteArray(copy));zip.closeEntry()}}finally{zip.close()}
    }
    static double[] readTimes(File file,int frames) {
        readTimesRows(CalciumCsv.read(file),frames)
    }
    static double[] readTimesRows(List rows,int frames) {
        int column=rows?rows[0].indexOf('Time_s'):-1
        if(column<0 || rows.size()!=frames+1)throw new IOException('Timestamps need a Time_s column and exactly one row per measured frame.')
        int frameColumn=rows[0].indexOf('Frame')
        double[] time=new double[frames]
        rows.drop(1).eachWithIndex{row,i->
            if(frameColumn>=0 && CalciumCsv.number(row[frameColumn])!=i+1)throw new IOException('Timestamp Frame values must match the measured frames, starting at 1.')
            time[i]=CalciumCsv.number(row[column])
        }
        for(int i=0;i<time.length;i++)if(!Double.isFinite(time[i])||time[i]<0||(i>0&&time[i]<=time[i-1]))throw new IOException('Timestamps must be finite, nonnegative, strictly increasing seconds.')
        time
    }
    static void attachTimes(File input,File output,double[] time) {
        def rows=CalciumCsv.read(input)
        if(rows.size()!=time.length+1 || rows[0].contains('Time_s'))throw new IOException('Timestamp rows do not match the fluorescence input.')
        def timed=[rows[0].take(1)+['Time_s']+rows[0].drop(1)]
        rows.drop(1).eachWithIndex{row,i->timed<<row.take(1)+[Double.toString(time[i])]+row.drop(1)}
        CalciumCsv.write(output,timed)
    }
}

/** Native Fiji trace view; missing samples remain gaps. */
/** Did the field move during the recording? Each frame is compared with a reference picture (20 frames spread over the second
 *  half), on a picture reduced four times from which the smooth lamp light and haze are taken off, so only structures the size
 *  of cells are matched. The shift at which a frame differs least from the reference is found in whole reduced pixels within
 *  +-12 (48 pixels) and refined between them from the neighbouring values. Shifts are then counted from where the field sits
 *  for most of the recording, or from the frames the cells were chosen on. Nothing is corrected: the check says which frames
 *  are displaced and by how much. Under SLIGHT pixels nothing is reported; from MOVED pixels on the outlines are no longer on
 *  the same cells (a cell body here is 20 to 40 pixels across). */
@CompileStatic
class CalciumMotion {
    static final int REDUCE = 4, REACH = 12, LOCAL = 7
    static final double SLIGHT = 2.5d, MOVED = 8d
    /** A frame reduced REDUCE times (block means) minus its own local mean over 15 x 15 reduced pixels. */
    static float[] small(ImageProcessor ip, int sw, int sh) {
        float[] blocks = new float[sw * sh]; int w = ip.getWidth()
        Object pixels = ip.getPixels(); short[] s16 = pixels instanceof short[] ? (short[]) pixels : null; float[] f32 = pixels instanceof float[] ? (float[]) pixels : null
        for (int y = 0; y < sh; y++) for (int x = 0; x < sw; x++) {
            double sum = 0d
            for (int dy = 0; dy < REDUCE; dy++) { int row = (y * REDUCE + dy) * w + x * REDUCE
                for (int dx = 0; dx < REDUCE; dx++) sum += s16 != null ? (double) (s16[row + dx] & 0xffff) : f32 != null ? (double) f32[row + dx] : (double) ip.getf(row + dx) }
            blocks[y * sw + x] = (float) (sum / (REDUCE * REDUCE))
        }
        // local mean from a summed-area table; the box is cut at the picture's edges
        double[] table = new double[(sw + 1) * (sh + 1)]
        for (int y = 0; y < sh; y++) { double run = 0d; for (int x = 0; x < sw; x++) { run += blocks[y * sw + x]; table[(y + 1) * (sw + 1) + x + 1] = table[y * (sw + 1) + x + 1] + run } }
        float[] out = new float[sw * sh]
        for (int y = 0; y < sh; y++) { int y0 = Math.max(0, y - LOCAL), y1 = Math.min(sh, y + LOCAL + 1)
            for (int x = 0; x < sw; x++) { int x0 = Math.max(0, x - LOCAL), x1 = Math.min(sw, x + LOCAL + 1)
                double sum = table[y1 * (sw + 1) + x1] - table[y0 * (sw + 1) + x1] - table[y1 * (sw + 1) + x0] + table[y0 * (sw + 1) + x0]
                out[y * sw + x] = (float) (blocks[y * sw + x] - sum / ((y1 - y0) * (x1 - x0))) } }
        return out
    }
    /** The shift (x, y, in reduced pixels, with a fraction) at which frame "a" differs least from "ref": mean absolute difference
     *  on every other pixel for every whole shift, then the position of the minimum between whole shifts from its two neighbours. */
    static double[] shift(float[] a, float[] ref, int sw, int sh) {
        int side = 2 * REACH + 1; double[] error = new double[side * side]
        double best = Double.MAX_VALUE; int bx = 0, by = 0
        for (int dy = -REACH; dy <= REACH; dy++) for (int dx = -REACH; dx <= REACH; dx++) {
            double sum = 0d; int count = 0
            for (int y = REACH; y < sh - REACH; y += 2) { int row = y * sw, moved = (y + dy) * sw + dx
                for (int x = REACH; x < sw - REACH; x += 2) { sum += Math.abs((double) a[moved + x] - (double) ref[row + x]); count++ } }
            double v = count > 0 ? sum / count : Double.MAX_VALUE; error[(dy + REACH) * side + dx + REACH] = v
            // ties and near-ties go to the smaller shift, so noise on a featureless field reads as "not moved"
            if (v < best * 0.999d || (v <= best && Math.abs(dx) + Math.abs(dy) < Math.abs(bx) + Math.abs(by))) { best = v; bx = dx; by = dy }
        }
        double fx = bx, fy = by; int at = (by + REACH) * side + bx + REACH
        if (bx > -REACH && bx < REACH) { double l = error[at - 1], h = error[at + 1], bend = l - 2d * best + h; if (bend > 0d) fx += Math.max(-0.5d, Math.min(0.5d, 0.5d * (l - h) / bend)) }
        if (by > -REACH && by < REACH) { double l = error[at - side], h = error[at + side], bend = l - 2d * best + h; if (bend > 0d) fy += Math.max(-0.5d, Math.min(0.5d, 0.5d * (l - h) / bend)) }
        return [fx, fy] as double[]
    }
    static double median(double[] values) { if (values.length == 0) return 0d; double[] s = (double[]) values.clone(); Arrays.sort(s); int m = (int) (s.length / 2); return s.length % 2 == 1 ? s[m] : (s[m - 1] + s[m]) / 2d }
    /** Shift of every frame in pixels, [n][2], counted from where the field sits for most of the recording (the median shift). */
    static double[][] shifts(ImageStack stack, int frames) {
        int n = Math.min(frames, stack.getSize()), sw = (int) (stack.getWidth() / REDUCE), sh = (int) (stack.getHeight() / REDUCE)
        double[][] out = new double[Math.max(0, n)][2]
        if (n < 3 || sw <= 4 * REACH || sh <= 4 * REACH) return out
        float[][] reduced = new float[n][]
        for (int f = 0; f < n; f++) reduced[f] = small(stack.getProcessor(f + 1), sw, sh)
        float[] ref = new float[sw * sh]; int used = 0
        for (int k = 0; k < 20; k++) { int f = (int) Math.round((n - 1) * (0.5d + 0.5d * k / 19d)); float[] r = reduced[f]; for (int i = 0; i < ref.length; i++) ref[i] = (float) (ref[i] + r[i]); used++ }
        for (int i = 0; i < ref.length; i++) ref[i] = (float) (ref[i] / used)
        double[] xs = new double[n], ys = new double[n]
        for (int f = 0; f < n; f++) { double[] s = shift(reduced[f], ref, sw, sh); xs[f] = s[0] * REDUCE; ys[f] = s[1] * REDUCE }
        double mx = median(xs), my = median(ys)
        for (int f = 0; f < n; f++) { out[f][0] = xs[f] - mx; out[f][1] = ys[f] - my }
        return out
    }
    static Map summary(double[][] shifts) { summary(shifts, 0, 0) }
    /** What the shifts mean, counted from the frames refFirst..refLast (the frames the cells were chosen on) when given, else
     *  from where the field sits for most of the recording. "runs": frames displaced by MOVED pixels or more; "slightRuns":
     *  by SLIGHT up to MOVED. Frames are counted from 1, both ends included. A frame's displacement is the middle value of
     *  five frames around it, so single noisy frames do not break a run. */
    static Map summary(double[][] shifts, int refFirst, int refLast) {
        int n = shifts.length; boolean chosen = refFirst >= 1 && refLast >= refFirst && refFirst <= n
        double rx = 0d, ry = 0d, within = 0d
        if (chosen) {
            int a = refFirst - 1, b = Math.min(n, refLast); double[] xs = new double[b - a], ys = new double[b - a]
            for (int f = a; f < b; f++) { xs[f - a] = shifts[f][0]; ys[f - a] = shifts[f][1] }
            double cx = median(xs), cy = median(ys)
            // how far the field moves inside those frames themselves (the second-largest distance, so one noisy frame does not count)
            double[] far = new double[b - a]; for (int f = a; f < b; f++) far[f - a] = Math.hypot(shifts[f][0] - cx, shifts[f][1] - cy)
            Arrays.sort(far); within = far.length >= 2 ? far[far.length - 2] : far[0]
            // a picture made from frames that moved is smeared and has no one place: then count from most of the recording instead
            if (within < MOVED) { rx = cx; ry = cy }
        }
        boolean steady = !chosen || within < MOVED
        double[] raw = new double[n], size = new double[n]
        for (int f = 0; f < n; f++) raw[f] = Math.hypot(shifts[f][0] - rx, shifts[f][1] - ry)
        for (int f = 0; f < n; f++) { int a = Math.max(0, f - 2), b = Math.min(n, f + 3); double[] part = new double[b - a]; for (int k = a; k < b; k++) part[k - a] = raw[k]; size[f] = Math.max(median(part), n < 5 ? raw[f] : 0d) }
        // a jolt of one or two frames is kept: the middle value would hide it
        for (int f = 0; f < n; f++) if (raw[f] >= MOVED && size[f] < MOVED) size[f] = raw[f]
        List<List<Integer>> runs = new ArrayList<List<Integer>>(), slightRuns = new ArrayList<List<Integer>>()
        double largest = 0d, slightLargest = 0d; int moved = 0, slight = 0, start = -1, slightStart = -1
        for (int f = 0; f <= n; f++) {
            boolean off = f < n && size[f] >= MOVED, little = f < n && size[f] >= SLIGHT && size[f] < MOVED
            if (off) { moved++; largest = Math.max(largest, size[f]); if (start < 0) start = f } else if (start >= 0) { runs.add([start + 1, f] as List<Integer>); start = -1 }
            if (little) { slight++; slightLargest = Math.max(slightLargest, size[f]); if (slightStart < 0) slightStart = f } else if (slightStart >= 0) { slightRuns.add([slightStart + 1, f] as List<Integer>); slightStart = -1 }
        }
        return [moved: moved, largest: (int) Math.round(largest), runs: runs, slight: slight, slightLargest: (int) Math.round(slightLargest), slightRuns: slightRuns, frames: n,
            reference: chosen && steady ? [refFirst, Math.min(n, refLast)] : null, smeared: chosen && !steady ? [refFirst, Math.min(n, refLast)] : null]
    }
    /** "frames 1-10", "frames 1-10 and 40-42", or for many pieces "34 frames between 2 and 119". */
    static String framesText(List runs) {
        if (runs == null || runs.isEmpty()) return ""
        List<String> parts = new ArrayList<String>(); int count = 0
        for (Object r : runs) { List run = (List) r; int a = ((Number) run.get(0)).intValue(), b = ((Number) run.get(1)).intValue(); count += b - a + 1; parts.add(a == b ? "" + a : a + "–" + b) }
        if (parts.size() > 3) { List first = (List) runs.get(0), last = (List) runs.get(runs.size() - 1); return count + " frames between " + first.get(0) + " and " + last.get(1) }
        return (count == 1 ? "frame " : "frames ") + parts.join(" and ")
    }
    /** One or two short sentences for the Recording card and the log; empty when the field stayed put. */
    static String sentence(Map m) {
        if (m == null) return ""
        List<String> parts = new ArrayList<String>()
        int moved = ((Number) m.get('moved')).intValue(), slight = ((Number) m.get('slight')).intValue(), n = ((Number) m.get('frames')).intValue()
        List reference = (List) m.get('reference')
        String where = reference != null ? " from where the cells are chosen (frames " + reference.get(0) + "–" + reference.get(1) + ")" : ""
        if (moved > 0) parts.add("The field moved (about " + m.get('largest') + " px" + where + ") during " + framesText((List) m.get('runs')) + ".")
        if (slight > 0) {
            List runs = (List) m.get('slightRuns'); List only = runs.size() == 1 ? (List) runs.get(0) : null
            String when = only != null && ((Number) only.get(1)).intValue() == n && ((Number) only.get(0)).intValue() > 1 ? "from frame " + only.get(0) + " on" : "during " + framesText(runs)
            parts.add(reference != null ? "The field sits about " + m.get('slightLargest') + " px" + where + " " + when + "." : "The field drifted by about " + m.get('slightLargest') + " px " + when + ".")
        }
        List smeared = (List) m.get('smeared')
        if (smeared != null) parts.add("The cells are chosen on frames " + smeared.get(0) + "\u2013" + smeared.get(1) + ", which include that movement, so their picture is smeared.")
        return parts.join(" ")
    }
}

/** Fast pixel work for showing movie frames next to the traces. Display only: nothing here is measured or saved. */
@CompileStatic
class CalciumFrameView {
    /** Mean of frames first..last (counted from 1, both included), pixel by pixel: the resting picture. */
    static float[] mean(ImageStack stack, int first, int last) { mean(stack, first, last, null) }
    /** The same without the frames in "skip" (frames in which the field moved), unless that would leave none. */
    static float[] mean(ImageStack stack, int first, int last, Set<Integer> skip) {
        int n = stack.getWidth() * stack.getHeight(), count = 0; double[] sum = new double[n]
        int from = Math.max(1, first), to = Math.min(stack.getSize(), last), kept = 0
        if (skip != null) for (int f = from; f <= to; f++) if (!skip.contains(f)) kept++
        for (int f = from; f <= to; f++) {
            if (skip != null && kept > 0 && skip.contains(f)) continue
            ImageProcessor ip = stack.getProcessor(f)
            for (int i = 0; i < n; i++) sum[i] += ip.getf(i)
            count++
        }
        float[] out = new float[n]
        if (count > 0) for (int i = 0; i < n; i++) out[i] = (float) (sum[i] / count)
        return out
    }
    /** The gray value (minus the resting picture, when given) below which this share of a spread sample of the pixels lies. */
    static double level(ImageProcessor ip, float[] rest, double share) {
        int n = ip.getWidth() * ip.getHeight(), step = Math.max(1, (int) (n / 200000)), k = 0
        float[] sample = new float[(int) ((n + step - 1) / step)]
        for (int i = 0; i < n && k < sample.length; i += step) sample[k++] = (float) (ip.getf(i) - (rest != null ? rest[i] : 0f))
        if (k == 0) return 0d
        Arrays.sort(sample, 0, k)
        return sample[(int) Math.max(0L, Math.min((long) (k - 1), Math.round(share * (k - 1))))]
    }
    /** Writes one frame (minus the resting picture, when given) into rgb, black at lo and white at hi. "a" and "b" are work
     *  buffers of the frame's size. With smooth, each pixel is first averaged with its eight neighbours (edges repeat), which
     *  calms the camera noise of a difference picture. */
    static void paint(ImageProcessor ip, float[] rest, double lo, double hi, int[] rgb, float[] a, float[] b, boolean smooth) { paint(ip, rest, lo, hi, rgb, a, b, smooth, Float.NaN) }
    /** With "none" given (the value that means no change in this frame), rims are taken off the difference first: see withoutRims. */
    static void paint(ImageProcessor ip, float[] rest, double lo, double hi, int[] rgb, float[] a, float[] b, boolean smooth, float none) {
        int n = rgb.length, width = ip.getWidth(), height = ip.getHeight()
        Object pixels = ip.getPixels()
        if (pixels instanceof short[]) { short[] p = (short[]) pixels; for (int i = 0; i < n; i++) a[i] = (float) (p[i] & 0xffff) }
        else if (pixels instanceof byte[]) { byte[] p = (byte[]) pixels; for (int i = 0; i < n; i++) a[i] = (float) (p[i] & 0xff) }
        else if (pixels instanceof float[]) System.arraycopy((float[]) pixels, 0, a, 0, n)
        else for (int i = 0; i < n; i++) a[i] = ip.getf(i)
        if (rest != null) for (int i = 0; i < n; i++) a[i] = (float) (a[i] - rest[i])
        float[] shown = a
        if (smooth && width > 2 && height > 2) {
            for (int y = 0; y < height; y++) {
                int row = y * width
                b[row] = (float) (a[row] + a[row] + a[row + 1]); b[row + width - 1] = (float) (a[row + width - 2] + a[row + width - 1] + a[row + width - 1])
                for (int i = row + 1; i < row + width - 1; i++) b[i] = (float) (a[i - 1] + a[i] + a[i + 1])
            }
            float ninth = (float) (1d / 9d)
            for (int i = 0; i < width; i++) { a[i] = (float) ((b[i] + b[i] + b[i + width]) * ninth); int j = n - width + i; a[j] = (float) ((b[j - width] + b[j] + b[j]) * ninth) }
            // rows 1 .. height-2: b holds the three-across sums, a receives the nine-pixel mean (a's own rows are no longer needed)
            for (int i = width; i < n - width; i++) a[i] = (float) ((b[i - width] + b[i] + b[i + width]) * ninth)
        }
        if (!Float.isNaN(none)) withoutRims(a, b, width, height, none)
        float low = (float) lo, scale = (float) (255d / Math.max(1e-9d, hi - lo))
        for (int i = 0; i < n; i++) {
            int g = (int) ((shown[i] - low) * scale)
            if (g < 0) g = 0 else if (g > 255) g = 255
            rgb[i] = (g << 16) | (g << 8) | g
        }
    }
    /** Takes the rims of shifted objects off a difference picture. A bright object that moves by a pixel gets brighter on one
     *  side and darker on the other; a cell that responds only gets brighter. Each pixel of "a" (a smoothed difference; "none" is
     *  the value that means no change) loses the size of the strongest darkening within about 10 pixels of it, beyond what the
     *  frame's own noise explains: block means of 4 x 4 pixels, the lowest over 5 x 5 blocks, less 3 noise SDs of the block means.
     *  So a rise beside an equal fall disappears and a rise on its own stays. "work" is a buffer of the frame's size. Display only. */
    static void withoutRims(float[] a, float[] work, int width, int height, float none) {
        int bw = (int) ((width + 3) / 4), bh = (int) ((height + 3) / 4), blocks = bw * bh
        if (blocks < 25 || work.length < 3 * blocks) return
        // work[0..blocks): block means, work[blocks..2*blocks): counts then the lowest nearby mean, work[2*blocks..3*blocks): sort
        for (int k = 0; k < 2 * blocks; k++) work[k] = 0f
        for (int y = 0; y < height; y++) { int by = y >> 2, row = y * width; for (int x = 0; x < width; x++) { int k = by * bw + (x >> 2); work[k] = (float) (work[k] + a[row + x] - none); work[blocks + k] = (float) (work[blocks + k] + 1f) } }
        for (int k = 0; k < blocks; k++) { work[k] = work[blocks + k] > 0f ? (float) (work[k] / work[blocks + k]) : 0f; work[2 * blocks + k] = Math.abs(work[k]) }
        Arrays.sort(work, 2 * blocks, 3 * blocks)
        float allowed = (float) (3d * 1.4826d * work[2 * blocks + (int) (blocks / 2)])   // 3 SDs, from the median size of a block mean
        for (int by = 0; by < bh; by++) for (int bx = 0; bx < bw; bx++) {
            float lowest = 0f
            for (int dy = -2; dy <= 2; dy++) { int yy = by + dy; if (yy < 0 || yy >= bh) continue
                for (int dx = -2; dx <= 2; dx++) { int xx = bx + dx; if (xx < 0 || xx >= bw) continue; float v = work[yy * bw + xx]; if (v < lowest) lowest = v } }
            float fall = (float) (lowest + allowed)
            work[blocks + by * bw + bx] = fall < 0f ? fall : 0f
        }
        for (int y = 0; y < height; y++) { int by = y >> 2, row = y * width; for (int x = 0; x < width; x++) { float fall = work[blocks + by * bw + (x >> 2)], v = a[row + x]; if (fall < 0f && v > none) { v = (float) (v + fall); a[row + x] = v < none ? none : v } } }
    }
    /** A colour frame: "gray" (a frame as paint() left it) becomes the green channel and "rise" (the smoothed rise above the
     *  resting picture, also from paint()) the magenta one, so a cell that responds turns from green to pink or white. */
    static void tint(int[] gray, int[] rise) {
        for (int i = 0; i < gray.length; i++) { int g = gray[i] & 255, m = rise[i] & 255; gray[i] = (m << 16) | (g << 8) | m }
    }
    /** The value below which this share of the finite numbers lies (0 when there are none). */
    static double share(double[][] rows, double share, boolean absolute) {
        int n = 0
        for (double[] row : rows) for (double v : row) if (!Double.isNaN(v) && !Double.isInfinite(v)) n++
        if (n == 0) return 0d
        double[] all = new double[n]; int k = 0
        for (double[] row : rows) for (double v : row) if (!Double.isNaN(v) && !Double.isInfinite(v)) all[k++] = absolute ? Math.abs(v) : v
        Arrays.sort(all)
        return all[(int) Math.max(0L, Math.min((long) (n - 1), Math.round(share * (n - 1))))]
    }
}

/** How each response call is drawn, everywhere in Step 4. The two call colours are the ones the SimpleCa²⁺ app uses. */
class CalciumCalls {
    static final Color RESPONDED = new Color(0x177465), NOT = new Color(0x925223), OTHER = new Color(0x8A94A0)
    static Color of(String call) { call == 'Responded' ? RESPONDED : call == 'Did not respond' ? NOT : OTHER }
    static final List<String> GROUPS = ['All cells', 'Responded', 'Did not respond', 'Excluded or unscorable']
    static boolean inGroup(String call, int group) { group <= 0 || (group == 1 && call == 'Responded') || (group == 2 && call == 'Did not respond') || (group == 3 && call != 'Responded' && call != 'Did not respond') }
}

class CalciumAnalysisPlot extends javax.swing.JPanel {
    double[] time=[], mean=[], sem=[], selected=[]
    String caption='Run the analysis to see traces', cell=''
    double baselineStart=0, baselineEnd=0, stimulusStart=0, stimulusEnd=0
    // every cell's trace behind the mean, coloured by its response call; "moment" is the time shown on the picture
    double[][] all=[]; List<String> allCalls=[]; boolean showAll=true; String allLabel='all cells'
    double moment=Double.NaN; Closure seek={double t->}
    Map responseCall=[:]
    private java.awt.image.BufferedImage layer; private String layerKey=''
    private int l=62,t=40,w=1,h=1; private double lo=0,hi=1
    CalciumAnalysisPlot() {
        setPreferredSize(new java.awt.Dimension(520,320)); setMinimumSize(new java.awt.Dimension(240,180)); setBackground(Color.white)
        def mouse=new java.awt.event.MouseAdapter(){
            void mousePressed(java.awt.event.MouseEvent e){move(e)}
            void mouseDragged(java.awt.event.MouseEvent e){move(e)}
            void move(java.awt.event.MouseEvent e){if(time.length>1&&w>1&&javax.swing.SwingUtilities.isLeftMouseButton(e))seek.call(time[0]+Math.max(0d,Math.min(1d,(e.x-l)/(double)w))*(time[-1]-time[0]))}
        }
        addMouseListener(mouse);addMouseMotionListener(mouse)
        setToolTipText('Click or drag to move through time; the picture follows')
    }
    /** Call after changing what is drawn (traces, selection, bands); moving the time cursor needs only repaint(). */
    void changed(){layerKey='';repaint()}
    Map referenceValues() {
        double base=CalciumCsv.number(responseCall.Baseline_mean?.toString())
        def values=[baseline:base]
        ['floor':'Floor_threshold','noise':'Noise_threshold','threshold':'Threshold_increase'].each{k,v->
            double fraction=CalciumCsv.number(responseCall[v]?.toString());values[k]=Double.isFinite(base)&&Double.isFinite(fraction)?base+(1d+base)*fraction:Double.NaN
        }
        values
    }
    double maxGap(){double gap=CalciumCsv.number(responseCall.Max_gap_s?.toString());if(Double.isFinite(gap)&&gap>0)return gap;if(time.length<2)return Double.NaN;def intervals=(1..<time.length).collect{time[it]-time[it-1]}.sort();1.5d*(intervals[(int)(intervals.size()/2)] as double)}
    int X(double x){l+(int)Math.round((x-time[0])/(time[-1]-time[0])*w)}
    int Y(double y){t+(int)Math.round((hi-y)/(hi-lo)*h)}
    protected void paintComponent(java.awt.Graphics graphics) {
        super.paintComponent(graphics)
        def g=(java.awt.Graphics2D)graphics.create()
        try {
            CalciumTheme.smooth(g); g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f)); def fm=g.getFontMetrics()
            if(time.length<2) {g.setColor(CalciumTheme.MUTED);String s=CalciumTheme.fit(caption,fm,getWidth()-24);g.drawString(s,(int)((getWidth()-fm.stringWidth(s))/2),(int)(getHeight()/2));return}
            // Everything except the time cursor is drawn once into a picture and reused, so dragging through time stays smooth with hundreds of cells.
            double scale=Math.max(1d,g.getTransform().getScaleX())
            String key=getWidth()+'x'+getHeight()+'@'+scale
            if(layer==null||layerKey!=key) {
                layer=new java.awt.image.BufferedImage((int)Math.ceil(getWidth()*scale),(int)Math.ceil(getHeight()*scale),java.awt.image.BufferedImage.TYPE_INT_RGB)
                def lg=layer.createGraphics();try{lg.scale(scale,scale);lg.setColor(Color.white);lg.fillRect(0,0,getWidth(),getHeight());still(lg)}finally{lg.dispose()}
                layerKey=key
            }
            g.drawImage(layer,0,0,getWidth(),getHeight(),null)
            if(Double.isFinite(moment)&&moment>=time[0]&&moment<=time[-1]) {
                int x=X(moment);g.setColor(new Color(31,41,51,200));g.setStroke(new java.awt.BasicStroke(1.2f));g.drawLine(x,t,x,t+h)
                int i=0;double best=Double.MAX_VALUE;for(int k=0;k<time.length;k++){double d=Math.abs(time[k]-moment);if(d<best){best=d;i=k}}
                String s=String.format(Locale.US,'%.1f s',time[i])
                if(i<selected.length&&Double.isFinite(selected[i])) {int y=Math.max(t,Math.min(t+h,Y(selected[i])));g.setColor(Color.white);g.fillOval(x-5,y-5,10,10);g.setColor(CalciumTheme.INK);g.fillOval(x-3,y-3,7,7)
                    s+=String.format(Locale.US,' · %s %.3f',cell,selected[i])}
                int sw=fm.stringWidth(s)+10,sx=Math.max(l,Math.min(l+w-sw,x-(int)(sw/2)))
                g.setColor(CalciumTheme.INK);g.fillRoundRect(sx,t+h-fm.getHeight()-6,sw,fm.getHeight()+4,8,8);g.setColor(Color.white);g.drawString(s,sx+5,t+h-8)
            }
        } finally {g.dispose()}
    }
    /** Bands, grid, every cell, the mean and the selected cell. */
    void still(java.awt.Graphics2D g) {
        CalciumTheme.smooth(g); g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f)); def fm=g.getFontMetrics()
        l=62;w=Math.max(1,getWidth()-l-18)
        boolean many=showAll&&all.length>0
        // the key above the plot: it takes as many lines as it needs, and the plot starts below it
        def entries=[[CalciumTheme.BLUE,2.4f,'mean ± SEM']]
        if(cell&&selected.length>0)entries<<[CalciumTheme.INK,1.7f,cell]
        if(cell&&selected.length==0)entries<<[null,0f,cell+': no trace (see Metrics + quality)']
        if(many) {
            int yes=allCalls.count{it=='Responded'} as int,no=allCalls.count{it=='Did not respond'} as int
            if(yes>0)entries<<[CalciumCalls.RESPONDED,1.4f,'responded ('+yes+')']
            if(no>0)entries<<[CalciumCalls.NOT,1.4f,'did not respond ('+no+')']
            if(all.length-yes-no>0)entries<<[CalciumCalls.OTHER,1.4f,(yes+no>0?'no call':allLabel)+' ('+(all.length-yes-no)+')']
        }
        int line=0,run=0;def places=[]
        entries.each{entry->int need=(entry[0]!=null?21:0)+fm.stringWidth(entry[2] as String)+14;if(run>0&&run+need>w){line++;run=0};places<<[run,line];run+=need}
        t=30+16*(line+1);h=Math.max(1,getHeight()-t-98)
        def values=[]; mean.eachWithIndex {v,i->if(Double.isFinite(v)) {values<<v;if(i<sem.length&&Double.isFinite(sem[i])) {values<<v-sem[i];values<<v+sem[i]}}}
        selected.each {if(Double.isFinite(it)) values<<it}; values<<0d
        Map references=referenceValues();references.values().each{if(Double.isFinite(it as double))values<<it}
        // With every cell drawn, the axis covers 99 % of all values (plus the mean and the selected cell); rarer extremes run off the edge.
        if(many) {values<<CalciumFrameView.share(all,0.005d,false);values<<CalciumFrameView.share(all,0.995d,false)}
        lo=values.min() as double;hi=values.max() as double
        if(hi-lo<0.001) hi=lo+1
        double pad=(hi-lo)*0.08;lo-=pad;hi+=pad
        double t0=time[0],t1=time[-1]
        double localEnd=CalciumCsv.number(responseCall.Baseline_end_s?.toString()),localStart=CalciumCsv.number(responseCall.Baseline_start_s?.toString())
        if(!Double.isFinite(localEnd))localEnd=stimulusStart
        if(!Double.isFinite(localStart))localStart=localEnd-CalciumCsv.number(responseCall.Baseline_s?.toString())
        double runStart=CalciumCsv.number(responseCall.Run_start_s?.toString()),runEnd=CalciumCsv.number(responseCall.Run_end_s?.toString())
        [[baselineStart,baselineEnd,new Color(0xEEF1F4),'normalization'],[localStart,localEnd,new Color(0xDCEBFA),'local baseline'],
         [stimulusStart,stimulusEnd,new Color(0xFFF0DC),'response'],[runStart,runEnd,new Color(0xC3EBC8),'qualifying run']].each { band ->
            if(!Double.isFinite(band[0] as double)||!Double.isFinite(band[1] as double))return
            int x0=Math.max(l,X(band[0] as double)),x1=Math.min(l+w,X(band[1] as double))
            if(x1>x0) {g.setColor(band[2] as Color);g.fillRect(x0,t,x1-x0,h);String name=band[3] as String
                if(fm.stringWidth(name)+8<=x1-x0){g.setColor(CalciumTheme.MUTED);g.drawString(name,x0+4,t+fm.getAscent()+2)}}
        }
        def yTicks=CalciumTheme.ticks(lo,hi,5); double yStep=yTicks.size()>1?yTicks[1]-yTicks[0]:1d
        yTicks.each { double v -> int y=Y(v);g.setColor(v==0d?new Color(0xAEB8C2):CalciumTheme.LINE);g.drawLine(l,y,l+w,y)
            String s=CalciumTheme.tick(v,yStep);g.setColor(CalciumTheme.MUTED);g.drawString(s,l-8-fm.stringWidth(s),y+(int)(fm.getAscent()/2)-1) }
        def xTicks=CalciumTheme.ticks(t0,t1,Math.max(2,(int)(w/90))); double xStep=xTicks.size()>1?xTicks[1]-xTicks[0]:1d
        xTicks.each { double v -> int x=X(v);g.setColor(CalciumTheme.LINE);g.drawLine(x,t+h,x,t+h+4)
            String s=CalciumTheme.tick(v,xStep);g.setColor(CalciumTheme.MUTED);g.drawString(s,x-(int)(fm.stringWidth(s)/2),t+h+6+fm.getAscent()) }
        g.setColor(new Color(0xC5CED8));g.drawRect(l,t,w,h)
        def clip=g.getClip();g.clipRect(l,t,w+1,h+1)
        double gap=maxGap()
        def path={double[] y->def p=new java.awt.geom.Path2D.Double();boolean pen=false
            for(int i=0;i<Math.min(time.length,y.length);i++){if(i>0&&Double.isFinite(gap)&&time[i]-time[i-1]>gap)pen=false;if(Double.isFinite(y[i])){double px=l+(time[i]-t0)/(t1-t0)*w,py=t+(hi-y[i])/(hi-lo)*h;if(pen)p.lineTo(px,py) else p.moveTo(px,py);pen=true} else pen=false};p}
        if(many) {
            int alpha=all.length>120?55:all.length>40?80:120
            g.setStroke(new java.awt.BasicStroke(0.9f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND))
            // non-responders first, so the responders are not hidden under them
            [false,true].each{boolean responders->all.eachWithIndex{double[] y,int i->String call=i<allCalls.size()?allCalls[i]:null
                if((call=='Responded')==responders){Color c=CalciumCalls.of(call);g.setColor(new Color(c.getRed(),c.getGreen(),c.getBlue(),alpha));g.draw(path(y))}}}
        }
        g.setColor(new Color(0,114,178,many?70:50))
        for(int i=1;i<time.length;i++) if((!Double.isFinite(gap)||time[i]-time[i-1]<=gap)&&i<mean.length&&i<sem.length&&[mean[i-1],mean[i],sem[i-1],sem[i]].every{Double.isFinite((double)it)}) {
            g.fillPolygon([X(time[i-1]),X(time[i]),X(time[i]),X(time[i-1])] as int[],[Y(mean[i-1]-sem[i-1]),Y(mean[i]-sem[i]),Y(mean[i]+sem[i]),Y(mean[i-1]+sem[i-1])] as int[],4)
        }
        g.setColor(Color.white);g.setStroke(new java.awt.BasicStroke(4.2f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));if(many)g.draw(path(mean))
        g.setColor(CalciumTheme.BLUE);g.setStroke(new java.awt.BasicStroke(2.4f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));g.draw(path(mean))
        if(selected.length>0){g.setColor(Color.white);g.setStroke(new java.awt.BasicStroke(3.6f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));g.draw(path(selected))
            g.setColor(CalciumTheme.INK);g.setStroke(new java.awt.BasicStroke(1.7f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));g.draw(path(selected))}
        if(Double.isFinite(stimulusStart)&&Double.isFinite(stimulusEnd))references.each{key,value->if(Double.isFinite(value as double)){
            g.setColor(key=='threshold'?CalciumTheme.GOOD:key=='baseline'?CalciumTheme.MUTED:CalciumTheme.WARN)
            g.setStroke(new java.awt.BasicStroke(key=='threshold'?2f:1f,java.awt.BasicStroke.CAP_BUTT,java.awt.BasicStroke.JOIN_MITER,10f,[5f,4f] as float[],0f))
            g.drawLine(X(stimulusStart),Y(value as double),X(stimulusEnd),Y(value as double))}}
        g.setClip(clip);g.setStroke(new java.awt.BasicStroke(1f))
        // the key, then the axis names
        entries.eachWithIndex{entry,int i->int kx=l+(places[i][0] as int),ky=20+16*(places[i][1] as int);String text=entry[2] as String
            if(entry[0]!=null){g.setColor(entry[0] as Color);g.setStroke(new java.awt.BasicStroke(entry[1] as float));g.drawLine(kx,ky,kx+16,ky);g.setStroke(new java.awt.BasicStroke(1f));g.setColor(CalciumTheme.INK);g.drawString(text,kx+21,ky+4)}
            else {g.setColor(CalciumTheme.WARN);g.drawString(CalciumTheme.fit(text,fm,l+w-kx),kx,ky+4)}}
        g.setColor(CalciumTheme.MUTED);String xl='Time (s)';g.drawString(xl,l+(int)((w-fm.stringWidth(xl))/2),getHeight()-60)
        def turned=(java.awt.Graphics2D)g.create();try{turned.rotate(-Math.PI/2);String yl='Global ΔF/F₀';turned.drawString(yl,-(t+(int)((h+fm.stringWidth(yl))/2)),14)}finally{turned.dispose()}
        g.setColor(CalciumTheme.MUTED);g.drawString(CalciumTheme.fit(caption,fm,w),l,getHeight()-37)
        String rule='Local threshold '+(responseCall.Threshold_increase?:'unavailable')+' ('+(responseCall.Threshold_source?:'')+') · span '+(responseCall.Run_duration_s?:'unavailable')+' s · '+(responseCall.Run_samples?:'')+' samples'
        g.drawString(CalciumTheme.fit(rule,fm,w),l,getHeight()-21)
        g.setColor(CalciumTheme.GOOD);g.drawString(CalciumTheme.fit('Green dashed: local threshold transformed to this global axis',fm,w),l,getHeight()-5)
    }
}
class CalciumAnalysisHeatmap extends javax.swing.JPanel {
    static final int TOP=46, BOTTOM=26, LEFT=78, ROW=15
    double[][] traces=[]; List<String> names=[]; Closure choose={int i->}; double[] time=[]
    List<String> calls=[]            // the response call of each row, for the coloured strip on the left
    double stimulusStart=Double.NaN, stimulusEnd=Double.NaN, baselineStart=Double.NaN, baselineEnd=Double.NaN
    String selected=''
    double maxGap=Double.NaN
    CalciumAnalysisHeatmap() {
        setBackground(Color.white);setPreferredSize(new java.awt.Dimension(650,260))
        addMouseListener(new java.awt.event.MouseAdapter(){void mouseClicked(java.awt.event.MouseEvent e){
            if(e.y>=TOP&&e.y<getHeight()-BOTTOM&&!names.isEmpty()) choose.call(Math.min(names.size()-1,(int)((e.y-TOP)*names.size()/(getHeight()-TOP-BOTTOM))))
        }})
        setToolTipText('Click a row to see that cell')
    }
    /** Tall enough for every row to be readable; the tab scrolls. */
    void fitRows(){setPreferredSize(new java.awt.Dimension(650,Math.max(260,names.size()*ROW+TOP+BOTTOM)));revalidate();repaint()}
    /** Colours run from −limit to +limit. The limit covers 99 % of the values, so one extreme cell cannot wash out all the others. */
    double limit() {
        double most=CalciumFrameView.share(traces,0.99d,true)
        if(!(most>0d))return 1d
        double power=Math.pow(10d,Math.floor(Math.log10(most))),unit=most/power
        return (unit<=1d?1d:unit<=1.5d?1.5d:unit<=2d?2d:unit<=3d?3d:unit<=5d?5d:10d)*power
    }
    protected void paintComponent(java.awt.Graphics graphics) {
        super.paintComponent(graphics);def g=(java.awt.Graphics2D)graphics.create()
        try {
            CalciumTheme.smooth(g);g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f));def fm=g.getFontMetrics()
            if(names.isEmpty()||time.length<2){g.setColor(CalciumTheme.MUTED);String s='Run the analysis to see the heatmap';g.drawString(s,(int)((getWidth()-fm.stringWidth(s))/2),(int)(getHeight()/2));return}
            double scale=limit(),largest=CalciumFrameView.share(traces,1d,true)
            g.setColor(CalciumTheme.INK);g.drawString('ΔF/F₀ of every cell',10,19)
            // colour key: blue below zero, white at zero, red above; gray where there is no value
            int kw=120,kx=getWidth()-kw-58
            if(kx>250) {
                for(int i=0;i<kw;i++){double v=(i/(double)(kw-1))*2-1;int fade=(int)(255*(1-Math.abs(v)));g.setColor(v>=0?new Color(255,fade,fade):new Color(fade,fade,255));g.drawLine(kx+i,9,kx+i,20)}
                g.setColor(CalciumTheme.LINE);g.drawRect(kx,9,kw-1,11);g.setColor(CalciumTheme.MUTED)
                String low='−'+CalciumTheme.shortNumber(scale),high=(largest>scale?'≥ ':'+')+CalciumTheme.shortNumber(scale)
                g.drawString(low,kx-6-fm.stringWidth(low),19);g.drawString(high,kx+kw+6,19)
            }
            int w=Math.max(1,getWidth()-LEFT-12),h=Math.max(1,getHeight()-TOP-BOTTOM)
            def X={double v->LEFT+(int)Math.round((v-time[0])/(time[-1]-time[0])*w)}
            // a strip above the rows shows where the baseline and the stimulus are
            [[baselineStart,baselineEnd,new Color(0xC9D1D9),'baseline'],[stimulusStart,stimulusEnd,new Color(0xF2C078),'stimulus']].each{band->
                if(Double.isFinite(band[0] as double)&&Double.isFinite(band[1] as double)){int x0=Math.max(LEFT,X(band[0] as double)),x1=Math.min(LEFT+w,X(band[1] as double))
                    if(x1>x0){g.setColor(band[2] as Color);g.fillRect(x0,TOP-12,x1-x0,6);String s=band[3] as String;if(fm.stringWidth(s)+4<x1-x0){g.setColor(CalciumTheme.MUTED);g.drawString(s,x0+2,TOP-15)}}}}
            g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-2f));def small=g.getFontMetrics()
            def intervals=(1..<time.length).collect{time[it]-time[it-1]}.sort();double typical=intervals[(int)(intervals.size()/2)],gap=Double.isFinite(maxGap)?maxGap:1.5d*typical
            double first=time[0],last=time[-1]
            g.setColor(new Color(0xEEEEEE));g.fillRect(LEFT,TOP,w,h)
            for(int c=0;c<names.size();c++) {
                int y=TOP+(int)(h*c/names.size()),y2=TOP+(int)(h*(c+1)/names.size())
                String call=c<calls.size()?calls[c]:null
                if(call!=null){g.setColor(CalciumCalls.of(call));g.fillRect(4,y,4,Math.max(1,y2-y-1))}
                if(y2-y>=11){boolean chosen=names[c]==selected;g.setFont(CalciumTheme.font(chosen?java.awt.Font.BOLD:java.awt.Font.PLAIN,-2f));g.setColor(chosen?CalciumTheme.INK:CalciumTheme.MUTED)
                    g.drawString(names[c],14,y+(int)((y2-y+small.getAscent())/2)-1)}
                for(int f=0;f<time.length;f++) {
                    double v=traces[c][f];int fade=Double.isFinite(v)?(int)(255*(1-Math.min(1,Math.abs(v)/scale))):0
                    g.setColor(!Double.isFinite(v)?new Color(0xD9DEE3):v>=0?new Color(255,fade,fade):new Color(fade,fade,255))
                    double before=f==0?0d:(time[f]-time[f-1]<=gap?(time[f]-time[f-1])/2d:typical/2d)
                    double after=f==time.length-1?0d:(time[f+1]-time[f]<=gap?(time[f+1]-time[f])/2d:typical/2d)
                    int x=LEFT+(int)(w*(time[f]-before-first)/(last-first)),x2=LEFT+(int)(w*(time[f]+after-first)/(last-first))
                    g.fillRect(x,y,Math.max(1,x2-x),Math.max(1,y2-y))
                }
                if(names[c]==selected){g.setColor(CalciumTheme.INK);g.drawRect(LEFT,y,w-1,Math.max(1,y2-y-1))}
                // a thin line where one group of calls ends and the next begins
                if(c>0&&c<calls.size()&&calls[c]!=calls[c-1]){g.setColor(Color.white);g.fillRect(LEFT,y-1,w,2);g.setColor(new Color(0x65717D));g.drawLine(4,y,LEFT+w,y)}
            }
            g.setColor(CalciumTheme.LINE);g.drawRect(LEFT,TOP,w-1,h-1)
            if(Double.isFinite(stimulusStart)&&stimulusStart>time[0]&&stimulusStart<time[-1]){g.setColor(new Color(31,41,51,150));int x=X(stimulusStart);g.drawLine(x,TOP-12,x,TOP+h)}
            g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f))
            def ticks=CalciumTheme.ticks(time[0],time[-1],Math.max(2,(int)(w/110)));double step=ticks.size()>1?ticks[1]-ticks[0]:1d
            g.setColor(CalciumTheme.MUTED)
            ticks.each{double v->int x=X(v);g.drawLine(x,TOP+h,x,TOP+h+3);String s=CalciumTheme.tick(v,step)+' s';g.drawString(s,Math.max(LEFT,Math.min(getWidth()-fm.stringWidth(s)-4,x-(int)(fm.stringWidth(s)/2))),getHeight()-8)}
        } finally {g.dispose()}
    }
}

/** The measured cells on a picture inside the analysis window: the colour picture, or the movie at the moment the time cursor is on. */
class CalciumCellPicture extends javax.swing.JPanel {
    java.awt.image.BufferedImage image; Roi[] rois=new Roi[0]; List<String> names=[]; Map<String,String> calls=[:]
    int selected=-1; Closure choose={int i->}
    boolean outlines=true, onGray=false   // draw every outline or only the selected one; colours that read on a gray movie frame
    int group=0                           // which calls are drawn (CalciumCalls.GROUPS)
    String note=''                        // what the picture shows, bottom right
    double magnify=1d, panX=0d, panY=0d   // wheel zoom and drag offset, on top of "fit the panel"
    private java.awt.Point dragFrom; private boolean dragged=false
    CalciumCellPicture() {
        setBackground(new Color(0x10161C));setPreferredSize(new java.awt.Dimension(420,320));setMinimumSize(new java.awt.Dimension(160,120))
        setToolTipText('Click a cell to see its trace. Scroll to zoom, drag to move, double-click to see the whole picture.')
        def mouse=new java.awt.event.MouseAdapter(){
            void mousePressed(java.awt.event.MouseEvent e){dragFrom=e.getPoint();dragged=false}
            void mouseDragged(java.awt.event.MouseEvent e){
                if(dragFrom==null||image==null)return
                if(Math.abs(e.x-dragFrom.x)+Math.abs(e.y-dragFrom.y)>3)dragged=true
                if(dragged){panX+=e.x-dragFrom.x;panY+=e.y-dragFrom.y;dragFrom=e.getPoint();repaint()}
            }
            void mouseReleased(java.awt.event.MouseEvent e){
                if(!javax.swing.SwingUtilities.isLeftMouseButton(e)||image==null||dragged)return
                if(e.getClickCount()>=2){whole();return}
                double[] f=frame();int x=(int)Math.floor((e.x-f[0])/f[2]),y=(int)Math.floor((e.y-f[1])/f[2])
                int i=(rois as List).findIndexOf{it.contains(x,y)}
                if(i>=0)choose.call(i)
            }
            void mouseWheelMoved(java.awt.event.MouseWheelEvent e){
                if(image==null)return
                double[] before=frame();double ix=(e.x-before[0])/before[2],iy=(e.y-before[1])/before[2]
                magnify=Math.max(1d,Math.min(16d,magnify*Math.pow(1.2d,-e.getPreciseWheelRotation())))
                if(magnify==1d){panX=0d;panY=0d} else {double[] after=frame();panX+=e.x-(after[0]+ix*after[2]);panY+=e.y-(after[1]+iy*after[2])}
                repaint()
            }
        }
        addMouseListener(mouse);addMouseMotionListener(mouse);addMouseWheelListener(mouse)
    }
    /** Where the picture is drawn: left, top and zoom, keeping its shape and leaving room for the key below. */
    double[] frame() {
        double fit=Math.min((getWidth()-8)/(double)image.getWidth(),(getHeight()-30)/(double)image.getHeight());if(!(fit>0))fit=0.01d
        double zoom=fit*magnify
        return [(getWidth()-image.getWidth()*zoom)/2d+panX,4d+(getHeight()-30-image.getHeight()*zoom)/2d+panY,zoom] as double[]
    }
    void whole(){magnify=1d;panX=0d;panY=0d;repaint()}
    /** Brings a cell into view when the picture is zoomed in. */
    void reveal(int index) {
        if(image==null||index<0||index>=rois.length||magnify==1d)return
        def b=rois[index].getBounds();double[] f=frame()
        panX+=getWidth()/2d-(f[0]+b.getCenterX()*f[2]);panY+=(getHeight()-30)/2d-(f[1]+b.getCenterY()*f[2])
    }
    /** Zooms so the cell and its neighbourhood fill the panel. */
    void zoomTo(int index) {
        if(image==null||index<0||index>=rois.length||getWidth()<=0)return
        def b=rois[index].getBounds();double fit=Math.min((getWidth()-8)/(double)image.getWidth(),(getHeight()-30)/(double)image.getHeight())
        double wanted=Math.min((getWidth()-8)/(double)(b.getWidth()*6d),(getHeight()-30)/(double)(b.getHeight()*6d))
        magnify=Math.max(1.01d,Math.min(16d,wanted/Math.max(1e-6d,fit)));panX=0d;panY=0d;reveal(index);repaint()
    }
    protected void paintComponent(java.awt.Graphics graphics) {
        super.paintComponent(graphics);if(image==null)return
        def g=(java.awt.Graphics2D)graphics.create()
        try {
            CalciumTheme.smooth(g);g.setRenderingHint(java.awt.RenderingHints.KEY_INTERPOLATION,magnify>3d?java.awt.RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR:java.awt.RenderingHints.VALUE_INTERPOLATION_BILINEAR)
            double[] f=frame();g.drawImage(image,(int)f[0],(int)f[1],(int)(image.getWidth()*f[2]),(int)(image.getHeight()*f[2]),null)
            def outline={Roi roi->def p=roi.getPolygon();def path=new java.awt.geom.Path2D.Double()
                for(int i=0;i<p.npoints;i++){double x=f[0]+p.xpoints[i]*f[2],y=f[1]+p.ypoints[i]*f[2];if(i==0)path.moveTo(x,y) else path.lineTo(x,y)};path.closePath();path}
            def dashed=new java.awt.BasicStroke(1.3f,java.awt.BasicStroke.CAP_BUTT,java.awt.BasicStroke.JOIN_ROUND,10f,[4f,3f] as float[],0f)
            Color yes=onGray?new Color(0x3DDBB5):Color.white, no=onGray?new Color(0xFFA24A):new Color(0xAEB8C2), other=new Color(0xAEB8C2)
            if(outlines) rois.eachWithIndex{Roi roi,int i->if(i==selected)return
                String call=calls[i<names.size()?names[i]:'']
                if(!calls.isEmpty()&&!CalciumCalls.inGroup(call,group))return
                def shape=outline(roi)
                g.setColor(new Color(0,0,0,150));g.setStroke(new java.awt.BasicStroke(call=='Responded'?3.2f:2.4f));g.draw(shape)
                if(call==null){g.setColor(Color.white);g.setStroke(new java.awt.BasicStroke(1.3f))}
                else if(call=='Responded'){g.setColor(yes);g.setStroke(new java.awt.BasicStroke(1.8f))}
                else if(call=='Did not respond'){g.setColor(no);g.setStroke(new java.awt.BasicStroke(1.2f))}
                else {g.setColor(other);g.setStroke(dashed)}
                g.draw(shape)}
            if(selected>=0&&selected<rois.length) {
                def shape=outline(rois[selected]);g.setColor(Color.black);g.setStroke(new java.awt.BasicStroke(5.5f));g.draw(shape)
                g.setColor(new Color(0x5BC0FF));g.setStroke(new java.awt.BasicStroke(2.6f));g.draw(shape)
                if(selected<names.size()) {g.setFont(CalciumTheme.font(java.awt.Font.BOLD,-1f));def fm=g.getFontMetrics();def b=shape.getBounds();String s=names[selected]
                    int x=(int)Math.max(2,Math.min(getWidth()-fm.stringWidth(s)-8,b.x+b.width+6)),y=(int)Math.max(fm.getAscent()+4,Math.min(getHeight()-34,b.y+fm.getAscent()))
                    g.setColor(new Color(0,0,0,190));g.fillRoundRect(x-4,y-fm.getAscent()-2,fm.stringWidth(s)+8,fm.getHeight()+3,8,8);g.setColor(Color.white);g.drawString(s,x,y)}
            }
            g.setColor(getBackground());g.fillRect(0,getHeight()-24,getWidth(),24)
            g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-2f));def fm=g.getFontMetrics()
            int right=getWidth()-8
            if(note){g.setColor(new Color(0x8A94A0));g.drawString(note,right-fm.stringWidth(note),getHeight()-9);right-=fm.stringWidth(note)+14}
            int x=8,y=getHeight()-9
            def chip={Color c,String text->if(x+18+fm.stringWidth(text)>right)return
                g.setColor(c);g.setStroke(new java.awt.BasicStroke(2f));g.drawLine(x,y-4,x+12,y-4);g.setColor(new Color(0xC5CED8));g.drawString(text,x+16,y);x+=16+fm.stringWidth(text)+12}
            if(calls.isEmpty()){g.setColor(new Color(0xC5CED8));g.drawString(CalciumTheme.fit('Click a cell to see its trace',fm,right-x),x,y)}
            else {chip(yes,'responded');chip(no,'did not respond');chip(other,'excluded or unscorable (dashed)')}
        } finally {g.dispose()}
    }
}

/** One headline number with a short caption under it. Clicking it can show just those cells. */
class CalciumStat extends javax.swing.JPanel {
    javax.swing.JLabel value=CalciumTheme.label('–',CalciumTheme.INK,java.awt.Font.BOLD,7f), caption=CalciumTheme.label(' ',CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f)
    Closure pressed=null; boolean active=false; private boolean hover=false
    CalciumStat(String text,Color colour=CalciumTheme.INK) {
        setLayout(new javax.swing.BoxLayout(this,javax.swing.BoxLayout.Y_AXIS));setOpaque(false)
        setBorder(javax.swing.BorderFactory.createEmptyBorder(8,12,8,12));value.setForeground(colour);caption.setText(text);add(value);add(caption)
        addMouseListener(new java.awt.event.MouseAdapter(){
            void mouseEntered(java.awt.event.MouseEvent e){hover=true;repaint()}
            void mouseExited(java.awt.event.MouseEvent e){hover=false;repaint()}
            void mouseReleased(java.awt.event.MouseEvent e){if(pressed!=null&&javax.swing.SwingUtilities.isLeftMouseButton(e))pressed.call()}
        })
    }
    void show(String number,String text){value.setText(number);caption.setText(text)}
    protected void paintComponent(java.awt.Graphics g0){def g=(java.awt.Graphics2D)g0.create();try{CalciumTheme.smooth(g)
        g.setColor(pressed!=null&&hover?CalciumTheme.BLUE_LIGHT:CalciumTheme.CARD);g.fillRoundRect(0,0,getWidth()-1,getHeight()-1,12,12)
        g.setColor(active?CalciumTheme.BLUE:CalciumTheme.LINE);g.setStroke(new java.awt.BasicStroke(active?2f:1f));g.drawRoundRect(active?1:0,active?1:0,getWidth()-(active?3:1),getHeight()-(active?3:1),12,12)}finally{g.dispose()}}
}

class CalciumResponseProtocol {
    static final List<String> COLUMNS=['Role','Stimulus','Baseline_start_s','Baseline_end_s','Stimulus_onset_s','Start_s','End_s','Floor_local_fraction','SD_multiplier','Duration_mode','Min_duration_s','Consecutive_samples','Max_gap_s','Min_baseline_samples','Noise_floor_local_fraction','Expected_interval_s']
    static final List<String> LABELS=['Role','Stimulus / control','Baseline start (s)','Baseline end (s)','Applied onset (s)','Response start (s)','Response end (s)','Minimum rise (fraction)','Baseline SD multiplier','Persistence: seconds / samples','Minimum duration (s)','Consecutive samples','Maximum gap (s; blank = auto)','Minimum baseline samples','Noise floor (fraction)','Expected interval (s; blank = auto)']
    static final List<String> LEGACY_EDITABLE=['Stimulus','Baseline_start_s','Start_s','End_s','Floor_local_fraction','SD_multiplier','Consecutive_samples','Max_gap_s']
    static final List<String> LEGACY_UNUSED=['Min_duration_s','Noise_floor_local_fraction','Expected_interval_s']
    static final List<String> LEGACY=['Role','Stimulus','Start_s','End_s','Baseline_s','Floor_local_fraction','SD_multiplier','Consecutive_samples','Max_gap_s']
    static Map window(String role,String stimulus,double baselineStart,double onset,double end) {
        [Role:role,Stimulus:stimulus,Baseline_start_s:baselineStart,Baseline_end_s:onset,Stimulus_onset_s:onset,Start_s:onset,End_s:end,
         Floor_local_fraction:0.10d,SD_multiplier:3d,Duration_mode:'seconds',Min_duration_s:1d,Consecutive_samples:2,Max_gap_s:null,Min_baseline_samples:3,Noise_floor_local_fraction:0d,Expected_interval_s:null]
    }
    static double numeric(Object value,String name,boolean optional=false) {
        if(optional&&(value==null||value.toString().trim()==''))return Double.NaN
        double n=CalciumCsv.number(value?.toString())
        if(!Double.isFinite(n))throw new IllegalArgumentException('Enter a number for '+(COLUMNS.indexOf(name)>=0?LABELS[COLUMNS.indexOf(name)].toLowerCase(Locale.ROOT):name)+'.')
        n
    }
    static Map payload(List<Map> windows,int version,String name,String validation,List overrides=null) {
        if(!(version in [1,2])||windows.size()<1||windows.size()>2)throw new IllegalArgumentException('Use one primary stimulus and at most one positive control.')
        def seen=[];def out=[]
        windows.eachWithIndex {original,index->
            Map row=new LinkedHashMap(original)
            row.Role=index==0?'Stimulus':'Control'
            row.Stimulus=row.Stimulus?.toString()?.trim()
            if(!row.Stimulus||seen.contains(row.Stimulus.toLowerCase(Locale.ROOT)))throw new IllegalArgumentException('Stimulus and control need distinct, nonempty names.')
            seen<<row.Stimulus.toLowerCase(Locale.ROOT)
            COLUMNS.findAll{!(it in ['Role','Stimulus','Duration_mode'])}.each{key->double value=numeric(row[key],key,key in ['Max_gap_s','Expected_interval_s']);row[key]=Double.isFinite(value)?value:null}
            if(row.End_s<=row.Start_s||row.Baseline_end_s<=row.Baseline_start_s||row.Baseline_start_s<0||row.Baseline_end_s>row.Stimulus_onset_s||row.Start_s<row.Stimulus_onset_s)throw new IllegalArgumentException('Baseline must end before onset; the response window must start at or after onset and end later.')
            if(row.Floor_local_fraction<0||row.SD_multiplier<0||(row.Floor_local_fraction==0&&row.SD_multiplier==0)||row.Noise_floor_local_fraction<0)throw new IllegalArgumentException('Choose nonnegative thresholds with an amplitude or noise requirement.')
            if(!(row.Duration_mode in ['seconds','samples'])||row.Min_duration_s<=0||row.Consecutive_samples<1||row.Consecutive_samples!=Math.floor(row.Consecutive_samples as double)||row.Min_baseline_samples<2||row.Min_baseline_samples!=Math.floor(row.Min_baseline_samples as double))throw new IllegalArgumentException('Use seconds or samples, a positive duration and integer sample counts.')
            ['Max_gap_s','Expected_interval_s'].each{key->if(row[key]!=null&&row[key]<=0)throw new IllegalArgumentException(key+' must be positive or blank.')}
            row.Baseline_s=row.Baseline_end_s-row.Baseline_start_s
            if(version==1) {
                if(row.Duration_mode!='samples'||row.Min_duration_s!=1d||row.Min_baseline_samples!=3d||row.Noise_floor_local_fraction!=0d||row.Expected_interval_s!=null)throw new IllegalArgumentException('Legacy mode has fixed sample persistence and three baseline samples; version-2-only duration/noise/interval settings are unsupported.')
                if(row.Baseline_end_s!=row.Start_s||row.Stimulus_onset_s!=row.Start_s||row.Floor_local_fraction<=0||(index==1&&row.Stimulus!='KCl'))throw new IllegalArgumentException('Legacy mode uses a baseline immediately before the response start and KCl as its optional control.')
                out<<row.subMap(LEGACY)
            } else {
                row.putAll([Protocol_version:2,Protocol_name:name,Engine_version:'2.0.0',Validation_status:validation,Signal_reference:'fixed_reference_fractional_change',SD_estimator:'sample_sd',Threshold_operator:'>=',Coverage_policy:'complete_windows',Response_direction:'increase',Smoothing:'none'])
                out<<row
            }
        }
        if(windows.size()==2)crossCheck(windows.collect{Map w->[name:w.Stimulus?.toString()?.trim(),baselineStart:numeric(w.Baseline_start_s,'Baseline_start_s'),baselineEnd:numeric(w.Baseline_end_s,'Baseline_end_s'),
            onset:numeric(w.Stimulus_onset_s,'Stimulus_onset_s'),end:numeric(w.End_s,'End_s')]})
        if(!name?.trim()||!(validation in ['provisional','user_defined','validated_by_user']))throw new IllegalArgumentException('Enter a protocol name and validation status.')
        if(version==1&&overrides)throw new IllegalArgumentException('Legacy protocols cannot contain field timing overrides.')
        [schema:'simpleca.response-protocol',schema_version:version,engine_version:version==2?'2.0.0':'legacy-1',name:version==2?name:'Legacy consecutive-sample protocol',validation_status:version==2?validation:'legacy_unvalidated',signal_reference:version==2?'fixed_reference_fractional_change':'fixed_reference_fractional_dff0',windows:out,field_overrides:overrides?:null]
    }
    /** The engine's two rules for a stimulus and a positive control (responder_validate_protocol in the pinned R source),
     *  checked here too so the page and window can say what to change before anything runs. */
    static void crossCheck(List<Map> w) {
        def at={double v->new java.text.DecimalFormat('0.##',java.text.DecimalFormatSymbols.getInstance(Locale.US)).format(v)}
        if(Math.max(w[0].onset as double,w[1].onset as double)<Math.min(w[0].end as double,w[1].end as double))
            throw new IllegalArgumentException('The '+w[0].name+' and '+w[1].name+' response windows overlap ('+at(w[0].onset as double)+'–'+at(w[0].end as double)+' s and '+at(w[1].onset as double)+'–'+at(w[1].end as double)+' s). End the '+w[0].name+' window at or before '+at(w[1].onset as double)+' s.')
        for(int i=0;i<2;i++){Map own=w[i],other=w[1-i]
            if(Math.max(own.baselineStart as double,other.onset as double)<Math.min(own.baselineEnd as double,other.end as double))
                throw new IllegalArgumentException('The '+own.name+' baseline ('+at(own.baselineStart as double)+'–'+at(own.baselineEnd as double)+' s) falls inside the '+other.name+' response window ('+at(other.onset as double)+'–'+at(other.end as double)+' s). '+
                    (i==1?'End the '+other.name+' window at or before '+at(own.baselineStart as double)+' s, or start the '+own.name+' baseline after '+at(other.end as double)+' s.':'Move the '+own.name+' baseline outside that window.'))}
    }
    static Map read(File file) {parse(file.getText('UTF-8'))}
    static Map parse(Object input) {
        def value
        try { value=new groovy.json.JsonSlurper().parseText(input instanceof String?input:groovy.json.JsonOutput.toJson(input)) }
        catch(groovy.json.JsonException bad) { throw new IOException('This protocol file is not valid JSON. Choose a protocol exported from SimpleCa or this plugin.') }
        if(!(value instanceof Map)||value.schema!='simpleca.response-protocol'||!(value.schema_version in [1,2])||!(value.windows instanceof List)||value.windows.size()<1||value.windows.size()>2)throw new IOException('Unsupported protocol schema. Use a SimpleCa version 1 or 2 protocol with one stimulus and optional control.')
        if(value.signal_reference!=(value.schema_version==2?'fixed_reference_fractional_change':'fixed_reference_fractional_dff0'))throw new IOException('Unsupported signal reference.')
        if(value.engine_version!=(value.schema_version==2?'2.0.0':'legacy-1'))throw new IOException('Unsupported protocol engine version.')
        value.windows.eachWithIndex {row,i->
            if(!(row instanceof Map)||row.Role!=(i==0?'Stimulus':'Control'))throw new IOException('Protocol roles must be Stimulus then optional Control.')
            if(value.schema_version==2) {
                [Signal_reference:'fixed_reference_fractional_change',SD_estimator:'sample_sd',Threshold_operator:'>=',Coverage_policy:'complete_windows',Response_direction:'increase',Smoothing:'none'].each{k,v->if(row[k]!=v)throw new IOException('Unsupported or missing '+k+'.')}
                if(row.Protocol_version!=2)throw new IOException('Protocol versions disagree.')
                [Engine_version:'engine_version',Protocol_name:'name',Validation_status:'validation_status'].each{k,v->if(row[k]!=value[v])throw new IOException(k=='Protocol_name'?'The protocol name at the top of this file ("'+value[v]+'") differs from the name stored with its '+row.Role+' window ("'+row[k]+'"). Change both to match, or import the file unchanged and rename the protocol on the page.':k=='Validation_status'?'The validation status at the top of this file differs from the one stored with its '+row.Role+' window. Change both to match, or import it unchanged and set the status on the page.':'Protocol metadata disagree: '+k)}
                double baseline=numeric(row.Baseline_s,'Baseline_s'),duration=numeric(row.Baseline_end_s,'Baseline_end_s')-numeric(row.Baseline_start_s,'Baseline_start_s')
                if(Math.abs(baseline-duration)>16d*Math.ulp(1d)*Math.max(1d,Math.abs(baseline)))throw new IOException('Baseline_s must match the explicit baseline start and end.')
            } else {
                row.Baseline_start_s=(row.Start_s as double)-(row.Baseline_s as double);row.Baseline_end_s=row.Start_s;row.Stimulus_onset_s=row.Start_s
                row.Duration_mode='samples';row.Min_duration_s=1d;row.Min_baseline_samples=3;row.Noise_floor_local_fraction=0d;row.Expected_interval_s=null
            }
        }
        if(value.field_overrides!=null&&!(value.field_overrides instanceof List))throw new IOException('Field timing overrides must be an array.')
        if(value.field_overrides) {
            if(value.schema_version!=2)throw new IOException('Legacy protocols cannot contain field timing overrides.')
            List timings=['Baseline_start_s','Baseline_end_s','Stimulus_onset_s','Start_s','End_s'];Set seen=new HashSet()
            value.field_overrides.each{override->
                if(!(override instanceof Map)||override.keySet().any{!(it in ['Field','Role','Stimulus']+timings)}||!(override.Field instanceof String)||!override.Field.trim())throw new IOException('Timing overrides require a field identifier and known timing columns.')
                Map original=value.windows.find{it.Role==override.Role}
                if(original==null||!seen.add(override.Field+'\r'+override.Role))throw new IOException('Each field and existing protocol role can have only one timing override.')
                if(override.containsKey('Stimulus')&&override.Stimulus!=original.Stimulus)throw new IOException('Timing override stimulus must match its protocol role.')
                Map combined=new LinkedHashMap(original)
                timings.each{k->if(override[k]!=null){if(!(override[k] instanceof Number))throw new IOException('Timing overrides must be finite numbers or null.');combined[k]=numeric(override[k],k)}}
                payload(value.windows.collect{it.Role==override.Role?combined:it} as List<Map>,2,value.name.toString(),value.validation_status.toString())
            }
        }
        // Validate supported fields without executing content. The R import validates the full contract again.
        payload(value.windows as List<Map>,value.schema_version as int,(value.name?:'Imported protocol').toString(),value.schema_version==1?'user_defined':value.validation_status.toString(),value.field_overrides as List)
        value as Map
    }
    static Map canonical(Object input) {Map value=parse(input);payload(value.windows as List<Map>,value.schema_version as int,value.name.toString(),value.schema_version==1?'user_defined':value.validation_status.toString(),value.field_overrides as List)}
    static String json(Map payload) {groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(payload))}
}

/** Step 4 binds displayed results to their completed protocol snapshot. */
class CalciumAnalysisWindow extends javax.swing.JFrame {
    Map measurement
    CalciumRRuntime runtime=new CalciumRRuntime()
    Closure performAnalysis={rt,input,settings,out,job->CalciumSimpleCaEngine.analyze(rt,input,settings,out,job)}
    CalciumLocalJob job
    boolean busy=false,closed=false,loading=false
    File results,timestampsFile,sensitivityFile,diagnosticFile
    String timestampsHash='',sensitivityHash='',diagnosticHash=''
    double[] explicitTime=[]
    List fieldOverrides=[]
    Map appliedManifest=[:]
    String appliedFingerprint='',submittedFingerprint=''
    List<Map> responseRows=[]
    Map<String,javax.swing.JTextField> fields=[:]
    Map<String,List<String>> exact=[:]   // Untouched rounded suggestions retain their exact source values.
    def setup=CalciumTheme.button('Set up analysis engine (one time)')
    def runButton=CalciumTheme.button('Analyze responses',true)
    def cancelButton=CalciumTheme.button('Cancel')
    def logButton=CalciumTheme.button('View log')
    def folderButton=CalciumTheme.button('Open analysis folder')
    def savePlots=CalciumTheme.button('Save plots as PNG')
    def loadProtocolButton=CalciumTheme.button('Load protocol')
    def saveProtocolButton=CalciumTheme.button('Save protocol')
    def loadTimesButton=CalciumTheme.button('Load timestamps CSV')
    def clearTimesButton=CalciumTheme.button('Use confirmed frame rate')
    def sensitivityButton=CalciumTheme.button('Load sensitivity grid CSV')
    def diagnosticButton=CalciumTheme.button('Load negative-control windows CSV')
    def clearDiagnosticsButton=CalciumTheme.button('Clear diagnostics')
    def templateButton=CalciumTheme.button('Save diagnostic templates')
    def overrideButton=CalciumTheme.button('View field timing overrides')
    def clearOverridesButton=CalciumTheme.button('Clear field timing overrides')
    def overrideNote=new javax.swing.JLabel('No field timing overrides. The recording clock is used directly.')
    def preset=new javax.swing.JComboBox<String>(['Provisional slow evoked: 1 s','Custom version 2','Legacy sample rule (starts at 2)'] as String[])
    def validation=new javax.swing.JComboBox<String>(['provisional','user_defined','validated_by_user'] as String[])
    def settingsButton=CalciumTheme.button('Hide settings')
    def previousCell=CalciumTheme.button('‹'), nextCell=CalciumTheme.button('›')
    def status=new javax.swing.JTextArea('Check acquisition timing and a protocol before analysis.',4,24)
    def appliedNote=new javax.swing.JTextArea('No completed analysis in this window.',3,24)
    def timingNote=new javax.swing.JLabel('Time is reconstructed from confirmed frames per second unless timestamps are loaded.')
    def diagnosticNote=new javax.swing.JLabel('Optional diagnostics are saved separately; they never replace the applied responder counts.')
    def progress=new javax.swing.JProgressBar(0,100)
    def confirmed=new javax.swing.JCheckBox('I checked the timing, baseline windows and cell selection',false)
    def kcl=new javax.swing.JCheckBox('Also score a separate positive control',false)
    def cells=new javax.swing.JComboBox<String>()
    def stimuli=new javax.swing.JComboBox<String>()
    def group=new javax.swing.JComboBox<String>(['All cells'] as String[])
    def allTraces=new javax.swing.JCheckBox('Every trace',true)
    def cellStatus=CalciumTheme.label(' ',CalciumTheme.MUTED)
    def heatmapOrder=new javax.swing.JComboBox<String>(['Responders first, strongest on top','As measured (original ROI order)'] as String[])
    def plot=new CalciumAnalysisPlot()
    def heatmap=new CalciumAnalysisHeatmap()
    def cellPicture=new CalciumCellPicture()
    def responding=new CalciumStat('responding',CalciumTheme.BLUE),registered=new CalciumStat('cells registered'),excluded=new CalciumStat('excluded'),unscorable=new CalciumStat('unscorable')
    def metrics=new javax.swing.JTable(),calls=new javax.swing.JTable(),counts=new javax.swing.JTable()
    def qc=new javax.swing.JTable(),registry=new javax.swing.JTable(),windowMetrics=new javax.swing.JTable(),replicates=new javax.swing.JTable(),sensitivity=new javax.swing.JTable(),diagnostics=new javax.swing.JTable()
    def windowTable=new javax.swing.JTable()
    def tabs=new javax.swing.JTabbedPane()
    javax.swing.table.DefaultTableModel windowModel
    Map<String,double[]> traces=[:]
    Map<String,Map<String,String>> callRows=[:]   // cell -> its row of the saved Cell response calls table
    // the movie beside the traces
    ImagePlus movie; java.awt.image.BufferedImage snapshot, frameImage; float[] rest, work, work2; double[] rawRange, changeRange
    int frame=1, restFirst=1, restLast=1
    def view=new javax.swing.JComboBox<String>(['Static colour overview'] as String[])
    def outlinesBox=new javax.swing.JCheckBox('All outlines',true)
    def zoomButton=CalciumTheme.button('Zoom to cell'), wholeButton=CalciumTheme.button('Whole picture')
    def playButton=CalciumTheme.button('▶ Play')
    def slider=new javax.swing.JSlider(1,2,1)
    def timeLabel=CalciumTheme.label(' ',CalciumTheme.MUTED)
    javax.swing.Timer clock, player, inputWatch
    private boolean listing=false

    long started
    String stage=''
    CalciumAnalysisWindow(Map data) {
        super('Step 4 · SimpleCa²⁺ response protocols');measurement=new LinkedHashMap(data);loading=true
        if(data.picture!=null) {
            snapshot=(data.picture as ImagePlus).getBufferedImage()
            measurement.picture=new ImagePlus('Analysis cells · '+data.recording,new ij.process.ColorProcessor(snapshot))
            cellPicture.image=snapshot;cellPicture.rois=(data.rois?:new Roi[0]) as Roi[];cellPicture.names=(data.names?:[]) as List<String>
        }
        // A matching original recording remains read-only while movie frames are displayed.
        if(data.movie instanceof ImagePlus&&snapshot!=null&&(data.movie as ImagePlus).width==snapshot.width&&(data.movie as ImagePlus).height==snapshot.height&&(data.movie as ImagePlus).stackSize>=(data.frames as int)) {
            movie=data.movie as ImagePlus;view.addItem('Movie: change from rest');view.addItem('Movie: as recorded')
        }
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
        def root=new javax.swing.JPanel(new java.awt.BorderLayout(12,12));root.setBorder(javax.swing.BorderFactory.createEmptyBorder(12,12,12,12));root.setBackground(CalciumTheme.PAGE)
        def row={java.awt.Component... parts->def p=new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT,6,3));p.setOpaque(false);p.setAlignmentX(0f);parts.each{p.add(it)};p.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE,(int)p.preferredSize.height));p}
        def config=new javax.swing.JPanel();config.setLayout(new javax.swing.BoxLayout(config,javax.swing.BoxLayout.Y_AXIS));config.setBackground(CalciumTheme.PAGE);config.setBorder(javax.swing.BorderFactory.createEmptyBorder(10,10,10,10))
        int n=data.frames as int,end=Math.min(n-2,Math.max(3,(data.baselineEnd?:20) as int));double rate=(data.hz?:Double.NaN) as double
        restFirst=(data.baselineStart?:1) as int;restLast=end
        def box={String key,Object value,String tip->
            String full=value==null?'':''+value,shown=(value instanceof Double||value instanceof BigDecimal||value instanceof Float)?CalciumTheme.shortNumber(value as double):full
            def f=new javax.swing.JTextField(shown,6);f.setFont(CalciumTheme.font());f.setToolTipText(shown!=full?tip+' (exact value: '+full+')':tip)
            if(shown!=full)exact[key]=[shown,full];fields[key]=f;f}
        def timing=new CalciumCard('','Acquisition and normalization',false)
        timing.body.add(row(CalciumTheme.label('Frames per second'),box('sampling_hz',Double.isFinite(rate)&&rate>0?rate:'','Confirm acquisition rate; load timestamps for irregular sampling')))
        timing.body.add(row(CalciumTheme.label('Baseline frames'),box('baseline_start',data.baselineStart?:1,'First normalization frame'),CalciumTheme.label('to'),box('baseline_end',end,'Last normalization frame')))
        timing.body.add(row(loadTimesButton));timing.body.add(row(clearTimesButton))
        def identity=new CalciumCard('','Protocol identity',false)
        def protocolName=box('protocol_name','Slow evoked response — provisional','Name recorded with each completed analysis');protocolName.columns=26
        identity.body.add(row(CalciumTheme.label('Name'),protocolName));identity.body.add(row(preset,validation));identity.body.add(row(loadProtocolButton,saveProtocolButton))
        identity.body.add(CalciumTheme.label('Imported settings remain editable for each stimulus and control.',CalciumTheme.MUTED))
        config.add(identity);config.add(javax.swing.Box.createVerticalStrut(8))
        def windows=new CalciumCard('','Stimulus and control windows',false)
        def protocolOwner=this
        windowModel=new javax.swing.table.DefaultTableModel(CalciumResponseProtocol.LABELS as Object[],0){boolean isCellEditable(int r,int c){protocolOwner.editableWindowCell(r,c)}}
        windowTable.setModel(windowModel);windowTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);windowTable.rowHeight=24;windowTable.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f));windowTable.tableHeader.font=CalciumTheme.font(java.awt.Font.BOLD,-1f)
        windowTable.setDefaultRenderer(Object,new javax.swing.table.DefaultTableCellRenderer(){
            java.awt.Component getTableCellRendererComponent(javax.swing.JTable table,Object value,boolean selected,boolean focus,int r,int c){
                def label=(javax.swing.JLabel)super.getTableCellRendererComponent(table,value,selected,focus,r,c)
                boolean legacy=protocolOwner.preset.selectedIndex==2;String key=CalciumResponseProtocol.COLUMNS[c],mode=table.model.getValueAt(r,CalciumResponseProtocol.COLUMNS.indexOf('Duration_mode'))?.toString()
                boolean unused=legacy&&key in CalciumResponseProtocol.LEGACY_UNUSED||key=='Consecutive_samples'&&mode=='seconds'||key=='Min_duration_s'&&mode=='samples'
                label.text=unused?'Not used ('+(legacy?'legacy':mode)+')':value instanceof Number?CalciumTheme.shortNumber(value as double):(value==null?'':value.toString())
                label.toolTipText=(value==null?'':value.toString())+(unused?' — retained for replay; the selected persistence mode does not use this setting':legacy&&!table.model.isCellEditable(r,c)?' — fixed by the legacy engine':'')
                if(!selected)label.foreground=!table.model.isCellEditable(r,c)?CalciumTheme.MUTED:CalciumTheme.INK;label
            }})
        for(int i=0;i<windowTable.columnCount;i++)windowTable.columnModel.getColumn(i).preferredWidth=i==0?80:170
        windowTable.columnModel.getColumn(CalciumResponseProtocol.COLUMNS.indexOf('Duration_mode')).cellEditor=new javax.swing.DefaultCellEditor(new javax.swing.JComboBox<String>(['seconds','samples'] as String[]))
        double onset=Double.isFinite(rate)&&rate>0?end/rate:10d,finish=Double.isFinite(rate)&&rate>0?(n-1)/rate:20d
        addWindow(CalciumResponseProtocol.window('Stimulus','Stimulus',0d,onset,finish))
        def windowScroll=new javax.swing.JScrollPane(windowTable);windowScroll.preferredSize=new java.awt.Dimension(720,110);windows.body.add(windowScroll);windows.body.add(row(kcl))
        windows.body.add(CalciumTheme.label('Baseline [start,end), applied onset, response [start,end]. Scroll horizontally to edit every rule.',CalciumTheme.MUTED))
        windows.body.add(CalciumTheme.label('A 1 s span is last minus first qualifying timestamp; at 2 Hz it needs at least 3 samples.',CalciumTheme.MUTED))
        windows.body.add(row(overrideButton,clearOverridesButton));windows.body.add(overrideNote);config.add(windows);config.add(javax.swing.Box.createVerticalStrut(8))
        def metadata=new CalciumCard('','Recording metadata',false)
        def condition=box('condition','','Condition or treatment group');condition.columns=20
        def unit=box('biological_unit','','Independent culture, animal or preparation; leave blank when unknown');unit.columns=20
        metadata.body.add(row(CalciumTheme.label('Condition / group'),condition,CalciumTheme.label('Biological unit'),unit))
        metadata.body.add(CalciumTheme.label('A biological unit is an independent culture, animal or preparation; cells and fields are not replicates.',CalciumTheme.MUTED))
        metadata.body.add(timingNote);config.add(metadata);config.add(javax.swing.Box.createVerticalStrut(8))
        def diagnosticCard=new CalciumCard('','Optional diagnostics',true)
        diagnosticCard.body.add(row(sensitivityButton,diagnosticButton));diagnosticCard.body.add(row(clearDiagnosticsButton,templateButton));diagnosticCard.body.add(diagnosticNote);config.add(diagnosticCard)
        [confirmed,kcl].each{CalciumTheme.plain(it)};confirmed.text='<html>I checked the timing and baselines<br>and reviewed cell selection</html>'
        [status,appliedNote].each{it.editable=false;it.lineWrap=true;it.wrapStyleWord=true;it.opaque=false;it.focusable=false;it.border=null;it.font=CalciumTheme.font();it.foreground=CalciumTheme.INK;it.alignmentX=0f;it.maximumSize=new java.awt.Dimension(Integer.MAX_VALUE,(int)it.preferredSize.height)}
        progress.setStringPainted(true);progress.alignmentX=0f;progress.maximumSize=new java.awt.Dimension(Integer.MAX_VALUE,24)
        def side=new javax.swing.JPanel();side.layout=new javax.swing.BoxLayout(side,javax.swing.BoxLayout.Y_AXIS);side.opaque=false;side.preferredSize=new java.awt.Dimension(305,640)
        def recording=CalciumTheme.label(''+(data.recording?:''),CalciumTheme.MUTED);recording.toolTipText=''+(data.recording?:'')
        def source=CalciumTheme.label("<html><body style='width:245px'>Uses full-precision background-corrected fluorescence. Background is subtracted once; no intensity lift. QC warnings remain reviewable.</body></html>",CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f)
        [CalciumTheme.label('Analyze responses',CalciumTheme.INK,java.awt.Font.BOLD,5f),recording,javax.swing.Box.createVerticalStrut(8),timing,javax.swing.Box.createVerticalStrut(8),confirmed,row(runButton,cancelButton),progress,status,appliedNote,row(setup),row(logButton),source,javax.swing.Box.createVerticalGlue()].each{side.add(it)}
        root.add(side,java.awt.BorderLayout.WEST)
        def stats=new javax.swing.JPanel(new java.awt.GridBagLayout());stats.opaque=false
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->def c=new java.awt.GridBagConstraints();c.gridx=i;c.fill=java.awt.GridBagConstraints.BOTH;c.weightx=i==0?2.2d:1d;c.insets=new java.awt.Insets(0,0,0,i<3?8:0);tile.preferredSize=new java.awt.Dimension(i==0?220:100,(int)tile.preferredSize.height);stats.add(tile,c)}
        // the Traces tab: which cells, the picture at the cursor's moment, the traces, and the time slider under both
        [cells,group,view,heatmapOrder,stimuli].each{it.setFont(CalciumTheme.font())}
        previousCell.setToolTipText('Previous cell (left arrow)');nextCell.setToolTipText('Next cell (right arrow)')
        group.setToolTipText('List and draw only the cells with this response call');allTraces.setToolTipText('Draw every listed cell behind the mean, coloured by its response call')
        def white={javax.swing.JPanel p->p.setOpaque(true);p.setBackground(Color.white);p}
        def pickers=white(row(CalciumTheme.label('Show'),group,CalciumTheme.label('Cell'),cells,previousCell,nextCell,CalciumTheme.label('Window'),stimuli,allTraces))
        cellStatus.setBorder(javax.swing.BorderFactory.createEmptyBorder(0,8,6,8))
        def head=new javax.swing.JPanel(new java.awt.BorderLayout());head.setBackground(Color.white);head.add(pickers,java.awt.BorderLayout.NORTH);head.add(cellStatus,java.awt.BorderLayout.SOUTH)
        def tracePanel=new javax.swing.JPanel(new java.awt.BorderLayout());tracePanel.setBackground(Color.white);tracePanel.add(head,java.awt.BorderLayout.NORTH)
        if(cellPicture.image!=null) {
            view.setToolTipText('What the picture shows. "Change from rest" is each frame minus the average of the baseline frames, so a responding cell lights up.')
            outlinesBox.setToolTipText('Untick to see the cells themselves: only the selected cell keeps its outline')
            def pictureTools=row(view,outlinesBox,zoomButton,wholeButton);pictureTools.setOpaque(true);pictureTools.setBackground(cellPicture.getBackground());outlinesBox.setForeground(new Color(0xC5CED8))
            def picturePane=new javax.swing.JPanel(new java.awt.BorderLayout());picturePane.add(cellPicture);picturePane.add(pictureTools,java.awt.BorderLayout.NORTH)
            def split=new javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT,picturePane,plot);split.setResizeWeight(0.5d);split.setBorder(null);split.setContinuousLayout(true);tracePanel.add(split)
            // Give the picture as much width as its shape can use at this height (at most 62 %), until the divider is dragged by hand.
            boolean dragged=false
            if(split.getUI() instanceof javax.swing.plaf.basic.BasicSplitPaneUI)(split.getUI() as javax.swing.plaf.basic.BasicSplitPaneUI).getDivider().addMouseListener(new java.awt.event.MouseAdapter(){void mousePressed(java.awt.event.MouseEvent e){dragged=true}})
            split.addComponentListener(new java.awt.event.ComponentAdapter(){void componentResized(java.awt.event.ComponentEvent e){
                if(dragged||split.getWidth()<=0||split.getHeight()<=0)return
                int ideal=(int)((split.getHeight()-pictureTools.getHeight()-30)*snapshot.getWidth()/(double)snapshot.getHeight())+8
                split.setDividerLocation((int)Math.max(240,Math.min(split.getWidth()*0.62d,ideal)))
            }})
        } else tracePanel.add(plot)
        slider.setMaximum(Math.max(2,n));slider.setOpaque(false);slider.setToolTipText('Move through the recording; the picture and the line on the plot follow')
        def timeRow=new javax.swing.JPanel(new java.awt.BorderLayout(8,0));timeRow.setBackground(Color.white);timeRow.setBorder(javax.swing.BorderFactory.createEmptyBorder(4,8,6,8))
        timeLabel.setPreferredSize(new java.awt.Dimension(170,(int)timeLabel.getPreferredSize().height))
        timeRow.add(playButton,java.awt.BorderLayout.WEST);timeRow.add(slider);timeRow.add(timeLabel,java.awt.BorderLayout.EAST)
        tracePanel.add(timeRow,java.awt.BorderLayout.SOUTH)
        tabs.setFont(CalciumTheme.font())
        def heatmapPanel=new javax.swing.JPanel(new java.awt.BorderLayout());heatmapPanel.setBackground(Color.white)
        heatmapPanel.add(white(row(CalciumTheme.label('Order'),heatmapOrder,CalciumTheme.label('strip on the left: response call · click a row to see that cell',CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f))),java.awt.BorderLayout.NORTH)
        def heatmapScroll=new javax.swing.JScrollPane(heatmap);heatmapScroll.setBorder(null);heatmapScroll.getVerticalScrollBar().setUnitIncrement(CalciumAnalysisHeatmap.ROW*3);heatmapPanel.add(heatmapScroll)
        tabs.addTab('Protocol',new javax.swing.JScrollPane(config));tabs.addTab('Traces',tracePanel);tabs.addTab('Heatmap',heatmapPanel)
        def tables=['Legacy whole-recording metrics':metrics,'Cell quality':qc,'ROI registry':registry,'Cell response calls':calls,'Responder counts':counts,'Stimulus-window metrics':windowMetrics,'Biological-unit summaries':replicates,'Sensitivity (diagnostic)':sensitivity,'Negative controls (diagnostic)':diagnostics]
        tables.each{title,table->table.autoResizeMode=javax.swing.JTable.AUTO_RESIZE_OFF;table.autoCreateRowSorter=true;table.selectionMode=javax.swing.ListSelectionModel.SINGLE_SELECTION;table.font=CalciumTheme.font(java.awt.Font.PLAIN,-1f);table.rowHeight=22;table.gridColor=CalciumTheme.LINE;table.tableHeader.font=CalciumTheme.font(java.awt.Font.BOLD,-1f);tabs.addTab(title,new javax.swing.JScrollPane(table))}
        def main=new javax.swing.JPanel(new java.awt.BorderLayout(0,10));main.opaque=false;main.add(stats,java.awt.BorderLayout.NORTH);main.add(tabs);main.add(row(folderButton,savePlots,settingsButton,CalciumTheme.label('Cell N is not biological replicate N.',CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f)),java.awt.BorderLayout.SOUTH);root.add(main)
        setContentPane(new javax.swing.JScrollPane(root));root.minimumSize=new java.awt.Dimension(860,660);root.preferredSize=new java.awt.Dimension(1120,760)
        setSize(1180,830);minimumSize=new java.awt.Dimension(360,300);def bounds=graphicsConfiguration.bounds;setSize(Math.min(width,(int)bounds.width-50),Math.min(height,(int)bounds.height-80))
        cancelButton.enabled=false;folderButton.enabled=false;savePlots.enabled=false;logButton.enabled=false;previousCell.enabled=false;nextCell.enabled=false
        setup.addActionListener({beginSetup()} as java.awt.event.ActionListener);runButton.addActionListener({beginAnalysis()} as java.awt.event.ActionListener)
        cancelButton.addActionListener({job?.cancel();stage='Cancelling...';status.text=stage} as java.awt.event.ActionListener)
        logButton.addActionListener({if(job?.logFile?.isFile())IJ.open(job.logFile.path)} as java.awt.event.ActionListener)
        folderButton.addActionListener({if(results!=null)java.awt.Desktop.desktop.open(results)} as java.awt.event.ActionListener)
        settingsButton.addActionListener({side.visible=!side.visible;settingsButton.text=side.visible?'Hide settings':'Show settings';root.revalidate();root.repaint()} as java.awt.event.ActionListener)
        savePlots.addActionListener({exportPlots()} as java.awt.event.ActionListener)
        loadProtocolButton.addActionListener({chooseFile('Load response protocol',false){loadProtocol(it)}} as java.awt.event.ActionListener)
        saveProtocolButton.addActionListener({chooseFile('Save response protocol',true){saveProtocol(it)}} as java.awt.event.ActionListener)
        loadTimesButton.addActionListener({chooseFile('Load per-frame timestamps',false){loadTimes(it)}} as java.awt.event.ActionListener)
        clearTimesButton.addActionListener({timestampsFile=null;timestampsHash='';explicitTime=new double[0];fields.sampling_hz.enabled=true;timingNote.text='Time source: confirmed constant FPS. Dropped frames cannot be inferred.';changed()} as java.awt.event.ActionListener)
        sensitivityButton.addActionListener({chooseFile('Sensitivity grid CSV',false){loadDiagnostic(it,false)}} as java.awt.event.ActionListener)
        diagnosticButton.addActionListener({chooseFile('Negative-control windows CSV',false){loadDiagnostic(it,true)}} as java.awt.event.ActionListener)
        clearDiagnosticsButton.addActionListener({sensitivityFile=null;diagnosticFile=null;sensitivityHash='';diagnosticHash='';diagnosticNote.text='No optional diagnostic files selected.';changed()} as java.awt.event.ActionListener)
        overrideButton.addActionListener({showOverrides()} as java.awt.event.ActionListener)
        clearOverridesButton.addActionListener({fieldOverrides=[];updateOverrideNote();changed()} as java.awt.event.ActionListener)
        templateButton.addActionListener({saveDiagnosticTemplates()} as java.awt.event.ActionListener)
        preset.addActionListener({applyPreset()} as java.awt.event.ActionListener);validation.addActionListener({changed()} as java.awt.event.ActionListener)
        kcl.addActionListener({if(!loading){if(kcl.selected&&windowModel.rowCount==1){Map primary=windowRows()[0];double baseEnd=(primary.End_s as double)+(primary.Baseline_end_s as double)-(primary.Baseline_start_s as double)+1d;Map control=CalciumResponseProtocol.window('Control','KCl',primary.End_s as double,baseEnd,baseEnd+10d);if(preset.selectedIndex==2)control.Duration_mode='samples';addWindow(control)};if(!kcl.selected&&windowModel.rowCount>1)windowModel.removeRow(1);changed()}} as java.awt.event.ActionListener)
        fields.values().each{field->field.document.addDocumentListener([insertUpdate:{e->changed()},removeUpdate:{e->changed()},changedUpdate:{e->changed()}] as javax.swing.event.DocumentListener)}
        windowModel.addTableModelListener({e->if(!loading&&!busy){if(preset.selectedIndex==2&&e.column==CalciumResponseProtocol.COLUMNS.indexOf('Start_s')&&e.firstRow>=0){loading=true;def start=windowModel.getValueAt(e.firstRow,e.column);windowModel.setValueAt(start,e.firstRow,3);windowModel.setValueAt(start,e.firstRow,4);loading=false};if(preset.selectedIndex==0){loading=true;preset.selectedIndex=1;validation.selectedItem='user_defined';loading=false};windowTable.repaint();changed()}} as javax.swing.event.TableModelListener)
        stimuli.addActionListener({if(!listing)selectWindow()} as java.awt.event.ActionListener)
        cells.addActionListener({if(!listing)selectCell(cells.selectedItem as String,true)} as java.awt.event.ActionListener)
        def step={int by->int count=cells.getItemCount();if(count>0)cells.setSelectedIndex((int)(((cells.getSelectedIndex()+by)%count+count)%count))}
        previousCell.addActionListener({step(-1)} as java.awt.event.ActionListener);nextCell.addActionListener({step(1)} as java.awt.event.ActionListener)
        group.addActionListener({if(!listing)listCells()} as java.awt.event.ActionListener)
        allTraces.addActionListener({plot.showAll=allTraces.selected;plot.changed()} as java.awt.event.ActionListener)
        heatmapOrder.addActionListener({orderHeatmap()} as java.awt.event.ActionListener)
        def primaryGroup={int which->String primary=responseRows.find{it.Role=='Stimulus'}?.Stimulus;if(primary!=null&&stimuli.selectedItem!=primary)stimuli.selectedItem=primary;showGroup(which)}
        responding.pressed={->primaryGroup(1)};registered.pressed={->primaryGroup(0)};excluded.pressed={->primaryGroup(3)};unscorable.pressed={->primaryGroup(3)}
        registered.setToolTipText('Click to list all cells');[excluded,unscorable].each{it.setToolTipText('Click to list excluded or unscorable cells')}
        heatmap.choose={int i->showGroup(0);cells.setSelectedItem(heatmap.names[i]);tabs.selectedIndex=1}
        cellPicture.choose={int i->if(!busy&&!closed&&cells.getItemCount()>0){String name=measurement.names[i];if((cells.getModel() as javax.swing.DefaultComboBoxModel).getIndexOf(name)<0)showGroup(0);cells.setSelectedItem(name)}}
        [metrics,qc,calls,windowMetrics].each{table->table.selectionModel.addListSelectionListener({e->
            if(!e.valueIsAdjusting&&table.selectedRow>=0) {
                int col=(0..<table.columnCount).find{table.getColumnName(it)=='Cell'} ?: 0
                String name=table.getValueAt(table.selectedRow,col).toString()
                if((cells.getModel() as javax.swing.DefaultComboBoxModel).getIndexOf(name)<0)showGroup(0)
                cells.setSelectedItem(name)
            }
        } as javax.swing.event.ListSelectionListener)}
        // time: the slider, the plot and the Play button all move the same time cursor
        plot.seek={double t->stopPlaying();if(plot.time.length>1){int best=0;double gap=Double.MAX_VALUE;plot.time.eachWithIndex{double v,int i->double d=Math.abs(v-t);if(d<gap){gap=d;best=i}};showFrame(best+1)}}
        slider.addChangeListener({e->if(slider.getValue()!=frame){stopPlaying();showFrame(slider.getValue())}} as javax.swing.event.ChangeListener)
        player=new javax.swing.Timer(50,{e->if(frame>=slider.getMaximum())stopPlaying() else showFrame(frame+1)} as java.awt.event.ActionListener)
        playButton.toolTipText='Display playback at 20 frames/s; cursor positions use the recording timestamps.'
        playButton.enabled=movie!=null&&view.selectedIndex>0
        playButton.addActionListener({if(player.isRunning())stopPlaying() else {if(frame>=slider.getMaximum())showFrame(1);player.start();playButton.setText('❚❚ Pause')}} as java.awt.event.ActionListener)
        view.addActionListener({if(view.selectedIndex==0)stopPlaying();playButton.enabled=movie!=null&&view.selectedIndex>0;showFrame(frame)} as java.awt.event.ActionListener)
        outlinesBox.addActionListener({cellPicture.outlines=outlinesBox.selected;cellPicture.repaint()} as java.awt.event.ActionListener)
        zoomButton.addActionListener({cellPicture.zoomTo(cellPicture.selected)} as java.awt.event.ActionListener)
        wholeButton.addActionListener({cellPicture.whole()} as java.awt.event.ActionListener)
        // left and right arrows step through the listed cells, unless a text box or list is using them
        def keys=getRootPane()
        keys.registerKeyboardAction({e->if(previousCell.enabled)step(-1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_LEFT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        keys.registerKeyboardAction({e->if(nextCell.enabled)step(1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_RIGHT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        clock=new javax.swing.Timer(500,{e->if(busy)status.text=stage+' · '+(long)((System.currentTimeMillis()-started)/1000)+' s elapsed'} as java.awt.event.ActionListener)
        addWindowListener(new java.awt.event.WindowAdapter(){void windowClosed(java.awt.event.WindowEvent e){closed=true;job?.cancel();clock.stop();player?.stop();inputWatch?.stop();(measurement.picture as ImagePlus)?.close()}})
        inputWatch=new javax.swing.Timer(1500,{e->if(!busy&&!closed&&results!=null&&(timestampsFile||sensitivityFile||diagnosticFile))changed()} as java.awt.event.ActionListener);inputWatch.start()
        loading=false;[allTraces,outlinesBox].each{CalciumTheme.plain(it)};updatePreviewTiming();if(movie!=null)view.selectedIndex=2;showFrame(1)
    }
    void stopPlaying(){player.stop();playButton.setText('▶ Play')}
    /** The resting picture and the display ranges for the movie views. Display only. */
    void prepareMovie() {
        def stack=movie.getStack();int n=measurement.frames as int
        rest=CalciumFrameView.mean(stack,restFirst,restLast)
        // the frame where the cells are brightest on average sets the display range, so it does not flicker from frame to frame
        int peak=Math.max(1,Math.min(n,(restLast+n).intdiv(2)))
        if(plot.mean.length>0){double top=Double.NEGATIVE_INFINITY;plot.mean.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top&&i<n){top=v;peak=i+1}}}
        def brightest=stack.getProcessor(peak);def resting=new ij.process.FloatProcessor(stack.getWidth(),stack.getHeight(),rest)
        double rise=CalciumFrameView.level(brightest,rest,0.9995d);if(!(rise>0d))rise=1d
        changeRange=[-0.1d*rise,rise] as double[]
        double low=CalciumFrameView.level(resting,null,0.005d),high=Math.max(CalciumFrameView.level(resting,null,0.9995d),CalciumFrameView.level(brightest,null,0.9995d));if(!(high>low))high=low+1d
        rawRange=[low,high] as double[]
        frameImage=new java.awt.image.BufferedImage(stack.getWidth(),stack.getHeight(),java.awt.image.BufferedImage.TYPE_INT_RGB);work=new float[stack.getWidth()*stack.getHeight()];work2=new float[work.length]
    }
    /** Moves the time cursor: the line on the plot, the slider and, in a movie view, the picture. */
    void showFrame(int wanted) {
        frame=Math.max(1,Math.min(slider.getMaximum(),wanted))
        if(slider.getValue()!=frame)slider.setValue(frame)
        double t=frame<=plot.time.length?plot.time[frame-1]:Double.NaN
        plot.moment=t;plot.repaint()
        timeLabel.setText((Double.isFinite(t)?String.format(Locale.US,'%.1f s · ',t):'')+'frame '+frame+' of '+slider.getMaximum())
        int mode=view.getSelectedIndex()
        if(movie!=null&&mode>0) {
            if(rest==null)prepareMovie()
            double[] range=mode==1?changeRange:rawRange
            CalciumFrameView.paint(movie.getStack().getProcessor(frame),mode==1?rest:null,range[0],range[1],((java.awt.image.DataBufferInt)frameImage.getRaster().getDataBuffer()).getData(),work,work2,mode==1)
            cellPicture.image=frameImage;cellPicture.onGray=true
            cellPicture.note=mode==1?'frame '+frame+' minus the mean of frames '+restFirst+'–'+restLast+', lightly smoothed':'frame '+frame+' as recorded'
        } else if(snapshot!=null) {cellPicture.image=snapshot;cellPicture.onGray=false;cellPicture.note=movie!=null?'colour picture (does not change with time)':''}
        cellPicture.repaint()
    }
    /** Preview time is provisional until results are saved; later draft edits cannot change the applied clock. */
    void updatePreviewTiming(){if(results!=null)return;String typed=fields.sampling_hz?.text?.trim();double rate=CalciumCsv.number(exact.sampling_hz!=null&&exact.sampling_hz[0]==typed?exact.sampling_hz[1]:typed);plot.time=explicitTime.length>0?explicitTime.clone():Double.isFinite(rate)&&rate>0?(0..<(measurement.frames as int)).collect{it/rate} as double[]:new double[0];plot.changed();showFrame(frame)}
    /** Uses the free space to the left of the Step 1-3 window, so the two never cover each other. */
    void besideControls(java.awt.Window panel) {
        try {
            def configuration=panel.getGraphicsConfiguration();def bounds=configuration.getBounds();def insets=java.awt.Toolkit.getDefaultToolkit().getScreenInsets(configuration)
            int left=(int)bounds.getX()+insets.left+8,top=(int)bounds.getY()+insets.top+8,free=panel.getX()-left-8,tall=(int)bounds.getHeight()-insets.top-insets.bottom-16
            if(free>=760)setBounds(left,top,Math.min(free,1240),Math.min(tall,820))
        } catch(Throwable ignored) { }
    }
    boolean editableWindowCell(int r,int c){String key=CalciumResponseProtocol.COLUMNS[c];if(c==0)return false;if(preset.selectedIndex==2)return key in CalciumResponseProtocol.LEGACY_EDITABLE&&!(r==1&&key=='Stimulus');String mode=windowModel.getValueAt(r,CalciumResponseProtocol.COLUMNS.indexOf('Duration_mode'))?.toString();if(key=='Min_duration_s')return mode=='seconds';if(key=='Consecutive_samples')return mode=='samples';true}
    void addWindow(Map row){windowModel.addRow(CalciumResponseProtocol.COLUMNS.collect{row[it]} as Object[])}
    List<Map> windowRows(){(0..<windowModel.rowCount).collect{r->CalciumResponseProtocol.COLUMNS.collectEntries{key->[key,windowModel.getValueAt(r,CalciumResponseProtocol.COLUMNS.indexOf(key))]}}}
    void stopEditing(){if(windowTable.isEditing()&&!windowTable.cellEditor.stopCellEditing())throw new IllegalArgumentException('Finish the window setting edit first.')}
    Map protocol(){stopEditing();CalciumResponseProtocol.payload(windowRows(),preset.selectedIndex==2?1:2,fields.protocol_name.text.trim(),validation.selectedItem as String,fieldOverrides)}
    Map fileSignature(File file,String loaded){if(file==null)return null;String current;try{current=file.isFile()?CalciumInstaller.digest(file.bytes):'missing'}catch(Throwable ignored){current='unreadable'};[path:file.absolutePath,loaded_sha256:loaded,current_sha256:current]}
    String fingerprint(){groovy.json.JsonOutput.toJson([fields:fields.collectEntries{k,v->[k,v.text]},windows:windowRows(),preset:preset.selectedIndex,validation:validation.selectedItem,overrides:fieldOverrides,timestamps:fileSignature(timestampsFile,timestampsHash),sensitivity:fileSignature(sensitivityFile,sensitivityHash),diagnostic:fileSignature(diagnosticFile,diagnosticHash)])}
    void changed(){if(loading||busy)return;if(results==null){updatePreviewTiming();return};if(results!=null){boolean dirty=fingerprint()!=appliedFingerprint;appliedNote.text=dirty?'Settings changed — rerun to apply. Display and plot export still show the completed run.':'Displayed results match the applied settings.';appliedNote.foreground=dirty?new Color(155,75,0):new Color(20,95,45)}}
    void applyPreset(){if(loading)return;loading=true
        boolean legacy=preset.selectedIndex==2
        if(preset.selectedIndex!=1){fields.protocol_name.text=legacy?'Legacy consecutive-sample protocol':'Slow evoked response — provisional';validation.selectedItem=legacy?'user_defined':'provisional';for(int r=0;r<windowModel.rowCount;r++){['Floor_local_fraction':0.10d,'SD_multiplier':3d,'Duration_mode':legacy?'samples':'seconds','Min_duration_s':1d,'Consecutive_samples':2].each{k,v->windowModel.setValueAt(v,r,CalciumResponseProtocol.COLUMNS.indexOf(k))};if(legacy){['Min_baseline_samples':3,'Noise_floor_local_fraction':0d,'Expected_interval_s':null].each{k,v->windowModel.setValueAt(v,r,CalciumResponseProtocol.COLUMNS.indexOf(k))};windowModel.setValueAt(windowModel.getValueAt(r,5),r,3);windowModel.setValueAt(windowModel.getValueAt(r,5),r,4);if(r==1)windowModel.setValueAt('KCl',r,1)}}}
        loading=false;windowTable.repaint();updateDiagnosticAvailability();changed()
    }
    void updateDiagnosticAvailability(){boolean available=preset.selectedIndex!=2&&!busy;sensitivityButton.enabled=available;diagnosticButton.enabled=available;templateButton.enabled=available;if(preset.selectedIndex==2){sensitivityFile=null;diagnosticFile=null;sensitivityHash='';diagnosticHash='';diagnosticNote.text='Optional sensitivity and negative-control diagnostics require version 2. Previously selected diagnostic files are detached.'}}
    void chooseFile(String title,boolean save,Closure action){try{def chooser=new javax.swing.JFileChooser(measurement.folder as File);chooser.dialogTitle=title;int answer=save?chooser.showSaveDialog(this):chooser.showOpenDialog(this);if(answer==javax.swing.JFileChooser.APPROVE_OPTION)action.call(chooser.selectedFile)}catch(Throwable t){status.text=t.message?:t.toString()}}
    void saveProtocol(File file){file.setText(CalciumResponseProtocol.json(protocol())+'\n','UTF-8');status.text='Saved draft protocol: '+file.name}
    void loadProtocol(File file){Map p=CalciumResponseProtocol.read(file);loading=true;try{preset.selectedIndex=p.schema_version==1?2:1;fields.protocol_name.text=p.name;validation.selectedItem=p.schema_version==1?'user_defined':p.validation_status;fieldOverrides=p.field_overrides?:[];windowModel.rowCount=0;p.windows.each{addWindow(it as Map)};kcl.selected=p.windows.size()==2;confirmed.selected=false}finally{loading=false};updateOverrideNote();updateDiagnosticAvailability();changed();status.text='Loaded '+file.name+(fieldOverrides?' · inspect the retained field timing overrides.':'')}
    void updateOverrideNote(){def active=fieldOverrides.findAll{it.Field==measurement.recording};overrideNote.text=fieldOverrides?'Retained overrides: '+fieldOverrides.size()+'. Recording '+measurement.recording+' uses '+active.size()+' override row(s); inspect before confirming.':'No field timing overrides. Table timings are used directly.';overrideNote.foreground=active?new Color(155,75,0):Color.darkGray}
    void showOverrides(){def columns=['Field','Role','Stimulus','Baseline_start_s','Baseline_end_s','Stimulus_onset_s','Start_s','End_s'];def model=new javax.swing.table.DefaultTableModel(fieldOverrides.collect{r->columns.collect{r[it]} as Object[]} as Object[][],columns as Object[]){boolean isCellEditable(int r,int c){false}};def table=new javax.swing.JTable(model);table.autoResizeMode=javax.swing.JTable.AUTO_RESIZE_OFF;for(int c=0;c<table.columnCount;c++)table.columnModel.getColumn(c).preferredWidth=140;def scroll=new javax.swing.JScrollPane(table);scroll.preferredSize=new java.awt.Dimension(850,240);javax.swing.JOptionPane.showMessageDialog(this,scroll,'Imported field timing overrides — recording '+measurement.recording,javax.swing.JOptionPane.INFORMATION_MESSAGE)}
    void loadDiagnostic(File file,boolean negative){String hash=CalciumInstaller.digest(file.bytes);if(negative){diagnosticFile=file;diagnosticHash=hash}else{sensitivityFile=file;sensitivityHash=hash};diagnosticNote.text=(negative?'Negative-control windows: ':'Sensitivity: ')+file.name+'; saved separately from applied results.';changed()}
    void loadTimes(File file){double[] times=CalciumMeasurement.readTimes(file,measurement.frames as int);timestampsFile=file;timestampsHash=CalciumInstaller.digest(file.bytes);explicitTime=times;fields.sampling_hz.enabled=false;timingNote.text='Time source: '+file.name+' · '+times.length+' timestamps · '+times[0]+' to '+times[-1]+' s';confirmed.selected=false;changed()}
    void saveDiagnosticTemplates(){File dir=measurement.folder as File;dir.mkdirs();CalciumCsv.write(new File(dir,'sensitivity_grid_template.csv'),[['Floor_local_fraction','SD_multiplier','Min_duration_s'],[0.10,3,1],[0.20,3,1],[0.10,5,1],[0.10,3,0.5],[0.10,3,2]]);CalciumCsv.write(new File(dir,'negative_control_windows_template.csv'),[['Diagnostic_window','Baseline_start_s','Baseline_end_s','Stimulus_onset_s','Start_s','End_s'],['Edit to a genuine unstimulated or vehicle window','','','','','']]);status.text='Saved diagnostic templates beside measurement outputs. Fill in control windows before loading.'}
    Map settings(){
        Map p=protocol();Map values=fields.collectEntries{k,v->String typed=v.text.trim();[k,exact[k]!=null&&exact[k][0]==typed?exact[k][1]:typed]};Map files=[:]
        if(timestampsFile){files.timestamps_file=timestampsFile.path;files.timestamps_sha256=timestampsHash}
        if(sensitivityFile){files.sensitivity_file=sensitivityFile.path;files.sensitivity_file_sha256=sensitivityHash}
        if(diagnosticFile){files.negative_control_file=diagnosticFile.path;files.negative_control_file_sha256=diagnosticHash}
        CalciumAnalysisSettings.check(values,p,confirmed.selected,measurement,files)
    }
    void launch(boolean installing,Closure action){
        if(busy)return
        busy=true;started=System.currentTimeMillis();stage=installing?'Preparing R setup':'Preparing analysis'
        setup.enabled=false;runButton.enabled=false;cancelButton.enabled=true;progress.value=0;progress.indeterminate=installing
        fields.values().each{it.enabled=false};windowTable.enabled=false;[confirmed,kcl,preset,validation,loadProtocolButton,saveProtocolButton,loadTimesButton,clearTimesButton,sensitivityButton,diagnosticButton,clearDiagnosticsButton,templateButton,overrideButton,clearOverridesButton].each{it.enabled=false};savePlots.enabled=false
        job=new CalciumLocalJob(logFile:new File(runtime.home,installing?'setup.log':'analysis.log'))
        job.progress={int p,String message->javax.swing.SwingUtilities.invokeLater{if(!closed){progress.value=p;progress.setString(p+'% · '+message);stage=message;logButton.enabled=job.logFile.isFile()}}}
        clock.start()
        Thread.startDaemon('SimpleCa analysis'){
            File output=null;String error=null
            try{output=action.call(job)}catch(java.util.concurrent.CancellationException ignored){error='Cancelled. Completed analyses are kept; partial output is not used.'}catch(Throwable t){error=t.message?:t.toString()}
            final File answer=output;final String failure=error
            javax.swing.SwingUtilities.invokeLater{
                busy=false;clock.stop();if(closed)return
                progress.indeterminate=false;cancelButton.enabled=false;setup.enabled=true;runButton.enabled=true;fields.values().each{it.enabled=true};fields.sampling_hz.enabled=timestampsFile==null;windowTable.enabled=true
                [confirmed,kcl,preset,validation,loadProtocolButton,saveProtocolButton,loadTimesButton,clearTimesButton,sensitivityButton,diagnosticButton,clearDiagnosticsButton,templateButton,overrideButton,clearOverridesButton].each{it.enabled=true};logButton.enabled=job.logFile.isFile();savePlots.enabled=results!=null;updateDiagnosticAvailability()
                if(failure){status.text=failure;status.foreground=CalciumTheme.INK;progress.string='Not completed';changed();return}
                if(installing){progress.value=100;status.text='R and protocol support are ready. Check settings, then analyze.';return}
                try{loadResults(answer);appliedFingerprint=submittedFingerprint;changed();String empty=CalciumAnalysisSettings.noEvaluableCells(answer);status.text=empty?:'Complete. Protocol, inputs and tables saved in '+answer.name+'.'+CalciumAnalysisSettings.byHowChosen(answer,measurement);status.foreground=empty?CalciumTheme.WARN:CalciumTheme.INK;progress.value=100;progress.string='Complete';tabs.selectedIndex=1}catch(Throwable t){status.text='Could not display results: '+t.message;progress.string='Display error; see saved outputs'}
            }
        }
    }
    void beginSetup(){launch(true){task->runtime.setup(task)}}
    void beginAnalysis(){try{Map s=settings();submittedFingerprint=fingerprint();launch(false){task->performAnalysis.call(runtime,measurement.csv,s,measurement.folder,task)}}catch(Throwable t){status.text=t.message?:t.toString()}}
    List<Map> table(javax.swing.JTable target,File file){
        if(!file.isFile()){target.model=new javax.swing.table.DefaultTableModel();return []}
        def rows=CalciumCsv.read(file);if(!rows)return []
        target.model=new javax.swing.table.DefaultTableModel(rows.drop(1).collect{it as Object[]} as Object[][],rows[0] as Object[]){boolean isCellEditable(int r,int c){false}}
        for(int i=0;i<target.columnCount;i++)target.columnModel.getColumn(i).preferredWidth=160
        rows.drop(1).collect{row->rows[0].withIndex().collectEntries{name,i->[name,row[i]]}}
    }
    /** The headline numbers and the per-cell calls, read from the tables the analysis just saved. Display only. */
    void summarize(File dir) {
        def countRows=CalciumCsv.read(new File(dir,'responder_counts.csv')),callRows=CalciumCsv.read(new File(dir,'responder_calls.csv'))
        def column={List rows,String name->(rows[0] as List).findIndexOf{it==name}}
        int role=column(countRows,'Role')
        def row=countRows.drop(1).find{role<0||it[role]=='Stimulus'}
        def cell={String name->int i=column(countRows,name);row!=null&&i>=0?''+row[i]:''}
        def whole={String name->double v=CalciumCsv.number(cell(name));Double.isFinite(v)?''+(long)Math.round(v):'–'}
        double percent=CalciumCsv.number(cell('Percent_responding'))
        responding.show(Double.isFinite(percent)?String.format(Locale.US,'%.0f %%',percent):'–',cell('Stimulus')+': '+whole('Responded')+' of '+whole('Evaluable')+' evaluable cells')
        responding.setToolTipText('Primary stimulus summary: '+cell('Stimulus')+'. Selecting another trace window does not change this saved primary summary.')
        registered.show(whole('Registered'),'cells registered');excluded.show(whole('Excluded'),'excluded');unscorable.show(whole('Unscorable'),'unscorable')
        int callRole=column(callRows,'Role'),callCell=column(callRows,'Cell'),callStatus=column(callRows,'Status');Map<String,String> found=[:]
        if(callCell>=0&&callStatus>=0)callRows.drop(1).each{if(callRole<0||it[callRole]=='Stimulus')found[''+it[callCell]]=''+it[callStatus]}
        cellPicture.calls=found;cellPicture.repaint()
    }
    String callOf(String name){callRows[name]?.Status}
    /** All linked views use one saved stimulus/control window, never draft settings. */
    void selectWindow(){
        callRows=responseRows.findAll{it.Stimulus==stimuli.selectedItem}.collectEntries{[(it.Cell):it]}
        cellPicture.calls=callRows.collectEntries{k,v->[k,v.Status]}
        Map first=callRows.values().find{true}?:[:];plot.stimulusStart=CalciumCsv.number(first.Start_s);plot.stimulusEnd=CalciumCsv.number(first.Effective_end_s?:first.End_s)
        listing=true;try{int keep=Math.max(0,group.selectedIndex);group.removeAllItems();CalciumCalls.GROUPS.eachWithIndex{String name,int i->group.addItem(name+' ('+(measurement.names as List<String>).count{callRows.isEmpty()||CalciumCalls.inGroup(callOf(it),i)}+')')};group.selectedIndex=callRows.isEmpty()?0:Math.min(keep,group.itemCount-1)}finally{listing=false}
        listCells();orderHeatmap();plot.changed()
    }
    /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. */
    void listCells() {
        int which=Math.max(0,group.getSelectedIndex());String keep=cells.getSelectedItem() as String
        List<String> listed=(measurement.names as List<String>).findAll{callRows.isEmpty()||CalciumCalls.inGroup(callOf(it),which)}
        listing=true
        try {cells.removeAllItems();listed.each{cells.addItem(it)};if(listed.contains(keep))cells.setSelectedItem(keep)} finally {listing=false}
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->tile.active=!callRows.isEmpty()&&((i==0&&which==1)||(i==1&&which==0)||(i>=2&&which==3));tile.repaint()}
        cellPicture.group=which
        List<String> drawn=listed.findAll{traces.containsKey(it)}
        plot.all=drawn.collect{traces[it]} as double[][];plot.allCalls=drawn.collect{callOf(it)};plot.allLabel=which==0?'all cells':CalciumCalls.GROUPS[which].toLowerCase()
        previousCell.enabled=listed.size()>1;nextCell.enabled=listed.size()>1
        if(cells.getItemCount()>0)selectCell(cells.getSelectedItem() as String,true)
        else {plot.cell='';plot.selected=new double[0];plot.responseCall=[:];cellStatus.setText('No cells in this group.');cellPicture.selected=-1;plot.changed();cellPicture.repaint()}
    }
    void showGroup(int which){if(which<group.getItemCount()&&group.getSelectedIndex()!=which)group.setSelectedIndex(which)}
    /** Rows of the heatmap: as measured, or responders first with the strongest on top (then non-responders, then cells without a call). */
    void orderHeatmap() {
        List<String> names=(measurement.names as List<String>).collect{it}
        if(heatmapOrder.getSelectedIndex()==0&&!callRows.isEmpty()) {
            def rank={String name->String call=callOf(name);call=='Responded'?0:call=='Did not respond'?1:2}
            def peak={String name->double v=CalciumCsv.number(callRows[name]?.Peak_increase);Double.isFinite(v)?v:Double.NEGATIVE_INFINITY}
            names=names.sort(false){a,b->rank(a)<=>rank(b)?:peak(b)<=>peak(a)?:(measurement.names as List).indexOf(a)<=>(measurement.names as List).indexOf(b)}
        }
        heatmap.names=names;heatmap.time=plot.time;heatmap.calls=callRows.isEmpty() ? [] : names.collect{callOf(it)?:''}
        heatmap.traces=names.collect{traces[it]?:new double[plot.time.length].collect{Double.NaN} as double[]} as double[][]
        heatmap.baselineStart=plot.baselineStart;heatmap.baselineEnd=plot.baselineEnd;heatmap.stimulusStart=plot.stimulusStart;heatmap.stimulusEnd=plot.stimulusEnd
        heatmap.maxGap=plot.maxGap();heatmap.fitRows()
    }
    void loadResults(File dir){
        if(!new File(dir,'COMPLETE').isFile())throw new IOException('This analysis has no COMPLETE marker.')
        results=dir;table(metrics,new File(dir,'metrics.csv'));table(qc,new File(dir,'cell_quality.csv'));table(registry,new File(dir,'roi_registry.csv'));responseRows=table(calls,new File(dir,'responder_calls.csv'));table(counts,new File(dir,'responder_counts.csv'));table(windowMetrics,new File(dir,'window_metrics.csv'));table(replicates,new File(dir,'replicate_counts.csv'));table(sensitivity,new File(dir,'sensitivity_counts.csv'));table(diagnostics,new File(dir,'negative_control_counts.csv'))
        def normalized=CalciumCsv.read(new File(dir,'normalized_traces.csv')),avg=CalciumCsv.read(new File(dir,'mean_trace.csv'))
        plot.time=normalized.drop(1).collect{CalciumCsv.number(it[0])} as double[]
        traces=[:];normalized[0].drop(1).eachWithIndex{name,i->traces[name]=normalized.drop(1).collect{CalciumCsv.number(it[i+1])} as double[]}
        plot.mean=avg.drop(1).collect{CalciumCsv.number(it[1])} as double[];plot.sem=avg.drop(1).collect{CalciumCsv.number(it[2])} as double[]
        appliedManifest=CalciumCsv.read(new File(dir,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
        int first=Math.max(0,(CalciumCsv.number(appliedManifest.baseline_start) as int)-1),last=Math.max(0,(CalciumCsv.number(appliedManifest.baseline_end) as int)-1)
        plot.baselineStart=plot.time[Math.min(first,plot.time.length-1)];plot.baselineEnd=plot.time[Math.min(last,plot.time.length-1)]
        listing=true;try{stimuli.removeAllItems();responseRows.collect{it.Stimulus}.unique().each{stimuli.addItem(it)}}finally{listing=false}
        try{summarize(dir)}catch(Throwable ignored){}
        restFirst=first+1;restLast=last+1;rest=null
        selectWindow();folderButton.enabled=true;savePlots.enabled=true
        int peak=1;double top=Double.NEGATIVE_INFINITY;plot.mean.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top){top=v;peak=i+1}}
        showFrame(peak)
        if(!appliedFingerprint)appliedFingerprint=fingerprint();changed()
    }
    void selectCell(String name,boolean highlight) {
        if(name==null)return
        plot.cell=name;plot.selected=traces[name]?:new double[0]
        Map call=callRows[name]?:[:];plot.responseCall=call
        plot.stimulusStart=CalciumCsv.number(call.Start_s);plot.stimulusEnd=CalciumCsv.number(call.Effective_end_s?:call.End_s)
        plot.caption=name+(traces.containsKey(name)?' · '+(call.Status?:'')+' · '+(call.Reason?:''):' · unavailable: inspect Cell quality for the exclusion reason')
        plot.changed()
        int index=measurement.names.findIndexOf{it==name}
        if(highlight&&index>=0&&measurement.picture!=null){Roi roi=(Roi)measurement.rois[index].clone();roi.setStrokeColor(Color.white);roi.setStrokeWidth(4);(measurement.picture as ImagePlus).setRoi(roi)}
        cellPicture.selected=index;cellPicture.reveal(cellPicture.selected);cellPicture.repaint()
        zoomButton.enabled=cellPicture.selected>=0
        heatmap.selected=name;heatmap.repaint()
        // one line saying what the analysis decided for this cell, in the words and numbers of the saved table
        Map<String,String> row=callRows[name]
        if(row==null) cellStatus.setText(callRows.isEmpty()?' ':'<html><b>'+name+'</b> · no response call was saved for this cell</html>')
        else {
            String colour=CalciumTheme.hex(CalciumCalls.of(row.Status));def number={String key,String format->double v=CalciumCsv.number(row[key]);Double.isFinite(v)?String.format(Locale.US,format,v):null}
            def parts=["<font color='"+colour+"'>●</font> <b>"+name+': '+CalciumOutlineProgress.escape(row.Status?:'')+'</b>',CalciumOutlineProgress.escape((row.Reason?:'').toLowerCase())]
            String rise=number('Peak_increase','%.2f'),at=number('Peak_s','%.0f'),need=number('Threshold_increase','%.2f')
            if(rise!=null)parts<<'largest rise '+rise+(at!=null?' at '+at+' s':'')+(need!=null?' (needed '+need+')':'')
            cellStatus.setText('<html>'+parts.findAll{it}.join(' · ')+'</html>')
            cellStatus.setToolTipText('Rise and threshold are fractions of the fluorescence just before the stimulus, as in the Cell response calls table.')
        }
    }
    void exportPlots(){
        if(results==null)return
        changed()
        try{
            ['trace':plot,'heatmap':heatmap].each{name,panel->int w=Math.max(650,panel.width),h=Math.max(260,panel.height);def image=new java.awt.image.BufferedImage(w,h,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics();try{panel.paint(g)}finally{g.dispose()};javax.imageio.ImageIO.write(image,'png',new File(results,name+'.png'))}
            new File(results,'plot_provenance.json').setText(CalciumResponseProtocol.json([result_folder:results.name,cell:cells.selectedItem,stimulus:stimuli.selectedItem,protocol_hash:appliedManifest.protocol_hash,manifest_sha256:CalciumInstaller.digest(new File(results,'analysis_manifest.csv').bytes),settings_pending:fingerprint()!=appliedFingerprint,source:'Completed analysis snapshot; display is unsmoothed'])+'\n','UTF-8')
            status.text='Saved completed-run plots and plot_provenance.json in '+results.name
        }catch(Throwable t){status.text='Could not save plots: '+t.message}
    }
}

/** Shared settings and portable protocol validation for both Step 4 views. */
/** Writes a picture as PNG or TIFF with its resolution recorded in the file, so it opens at the intended size in a layout
 *  program or a journal's checker. TIFF can be LZW- or ZIP-compressed; both are lossless. */
class CalciumPictureFile {
    static final List FORMATS=['png','tiff'],COMPRESSIONS=['lzw','zip','none']
    static String extension(String format){format=='tiff'?'tif':'png'}
    static byte[] encode(java.awt.image.BufferedImage source,String format,int dpi,String compression) {
        if(!(format in FORMATS))throw new IllegalArgumentException('Choose PNG or TIFF.')
        if(!(compression in COMPRESSIONS))throw new IllegalArgumentException('TIFF compression is LZW, ZIP or none.')
        if(dpi<72||dpi>1200)throw new IllegalArgumentException('Resolution must be 72 to 1200 dpi.')
        def bytes=new ByteArrayOutputStream(1<<20);def stream=javax.imageio.ImageIO.createImageOutputStream(bytes)
        def writer=javax.imageio.ImageIO.getImageWritersByFormatName(format).next();def param=writer.getDefaultWriteParam()
        try {
            writer.setOutput(stream)
            if(format=='png') {
                def meta=writer.getDefaultImageMetadata(javax.imageio.ImageTypeSpecifier.createFromRenderedImage(source),param)
                def root=new javax.imageio.metadata.IIOMetadataNode('javax_imageio_png_1.0'),size=new javax.imageio.metadata.IIOMetadataNode('pHYs')
                String perMetre=''+(Math.round(dpi/0.0254d) as long)
                size.setAttribute('pixelsPerUnitXAxis',perMetre);size.setAttribute('pixelsPerUnitYAxis',perMetre);size.setAttribute('unitSpecifier','meter');root.appendChild(size)
                meta.mergeTree('javax_imageio_png_1.0',root)
                writer.write(null,new javax.imageio.IIOImage(source,null,meta),param)
            } else {
                // TIFF without transparency: on white, three colour channels
                def flat=new java.awt.image.BufferedImage(source.width,source.height,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=flat.createGraphics()
                g.setColor(java.awt.Color.WHITE);g.fillRect(0,0,flat.width,flat.height);g.drawImage(source,0,0,null);g.dispose()
                if(compression!='none'){param.setCompressionMode(javax.imageio.ImageWriteParam.MODE_EXPLICIT);param.setCompressionType(compression=='zip'?'Deflate':'LZW')}
                def meta=writer.getDefaultImageMetadata(javax.imageio.ImageTypeSpecifier.createFromRenderedImage(flat),param)
                def directory=javax.imageio.plugins.tiff.TIFFDirectory.createFromMetadata(meta),tags=javax.imageio.plugins.tiff.BaselineTIFFTagSet.getInstance()
                long[][] resolution=[[dpi as long,1L] as long[]] as long[][]
                directory.addTIFFField(new javax.imageio.plugins.tiff.TIFFField(tags.getTag(javax.imageio.plugins.tiff.BaselineTIFFTagSet.TAG_X_RESOLUTION),javax.imageio.plugins.tiff.TIFFTag.TIFF_RATIONAL,1,resolution))
                directory.addTIFFField(new javax.imageio.plugins.tiff.TIFFField(tags.getTag(javax.imageio.plugins.tiff.BaselineTIFFTagSet.TAG_Y_RESOLUTION),javax.imageio.plugins.tiff.TIFFTag.TIFF_RATIONAL,1,resolution))
                directory.addTIFFField(new javax.imageio.plugins.tiff.TIFFField(tags.getTag(javax.imageio.plugins.tiff.BaselineTIFFTagSet.TAG_RESOLUTION_UNIT),javax.imageio.plugins.tiff.TIFFTag.TIFF_SHORT,1,[javax.imageio.plugins.tiff.BaselineTIFFTagSet.RESOLUTION_UNIT_INCH as char] as char[]))
                writer.write(null,new javax.imageio.IIOImage(flat,null,directory.getAsMetadata()),param)
            }
        } finally {stream.close();writer.dispose()}
        bytes.toByteArray()
    }
    /** The picture the page drew (PNG data), re-written in the chosen format. Returns [bytes, extension]. */
    static List fromPage(String png,Object format,Object dpi,Object compression) {
        byte[] sent=Base64.decoder.decode(png.replaceFirst(/^data:image\/png;base64,/,''))
        if(sent.length<8||sent[1]!=(byte)0x50||sent[2]!=(byte)0x4E||sent[3]!=(byte)0x47)throw new IllegalArgumentException('That is not a PNG picture.')
        String kind=format==null?'png':(''+format).toLowerCase()
        if(kind=='png'&&dpi==null)return [sent,'png']   // as drawn, nothing to add
        def image=javax.imageio.ImageIO.read(new ByteArrayInputStream(sent));if(image==null)throw new IllegalArgumentException('The picture could not be read.')
        int dots=dpi instanceof Number?(dpi as int):(dpi==null?300:Integer.parseInt((''+dpi).trim()))
        [encode(image,kind,dots,compression==null?'lzw':(''+compression).toLowerCase()),extension(kind)]
    }
}
class CalciumAnalysisSettings {
    static final List VALUE_KEYS=['sampling_hz','baseline_start','baseline_end','condition','biological_unit']
    static final Map FILE_KEYS=[timestamps:'timestamps_file',sensitivity:'sensitivity_file',negative_control:'negative_control_file']
    static Map<String,Object> suggested(Map data) {
        int n=data.frames as int,end=Math.min(n-2,Math.max(3,(data.baselineEnd?:20) as int))
        double rate=(data.hz?:Double.NaN) as double
        [sampling_hz:Double.isFinite(rate)&&rate>0?rate:'',baseline_start:data.baselineStart?:1,baseline_end:end,condition:'',biological_unit:'']
    }
    /** Step 4 starts from what was confirmed in Fiji, so nothing is typed twice. With frames confirmed for blind selection
     *  ("frames before the stimulus", a to b): the stimulus is taken to follow frame b, and the resting baseline is the last 20
     *  of those frames, or all of them when there are fewer (close to the stimulus, so fading over a long wait does not count).
     *  Without them: the 20 frames before the rise the plugin found, or frames 1-20 when it found none. */
    static Map framesFromFiji(def steps) {
        int n=steps.n as int,a=steps.blindFirst as int,b=steps.blindLast as int
        if(b>0&&b>=a){int first=Math.max(a,b-19)
            return [first:first,last:b,chosen:[first:a,last:b],note:'Baseline frames '+first+'–'+b+' and the stimulus (after frame '+b+') are placeholders taken from the frames you confirmed in Fiji ('+a+'–'+b+'), which are the frames before the response. Set when the stimulus was added, then run.']}
        if(steps.found)return [first:steps.b0 as int,last:steps.b1 as int,note:'Baseline frames '+steps.b0+'–'+steps.b1+' and the stimulus time come from the rise the plugin found near frame '+steps.onset+'. Check them, then run.']
        [first:1,last:Math.min(20,n-1),note:'No response was found, so the baseline is a guess (frames 1–'+Math.min(20,n-1)+'). Set the stimulus time, then run.']
    }
    /** The settings last used for this recording are kept beside its tables and filled in again the next time Step 4 opens
     *  for the same recording with the same frames from Fiji. Uploaded files are not kept. */
    /** After which frame the stimulus was added in the last analysis that was run, on any recording: the starting value for the
     *  next recording, so the one fact is asked once per series. 0 when no analysis was run yet. */
    static int lastStimulusFrame(){try{return (int)Math.rint(Prefs.get('calcium.step4.stimulus_after_frame',0d))}catch(Exception ignored){return 0}}
    static void rememberStimulusFrame(Map data,Map draft){
        try{
            if(draft?.timestamps!=null)return
            double hz=CalciumCsv.number(''+draft.values.sampling_hz);def w=(draft.protocol.windows as List)[0] as Map;def at=w.Stimulus_onset_s!=null?w.Stimulus_onset_s:w.Start_s
            if(!(hz>0)||at==null)return
            double onset=CalciumCsv.number(''+at);int n=data.frames as int,before=0
            for(int k=0;k<n;k++)if(k/hz<onset-1e-9)before=k+1
            if(before>=1&&before<n){Prefs.set('calcium.step4.stimulus_after_frame',before as double)}
        }catch(Exception ignored){}
    }
    static File rememberedFile(Map data){data.folder instanceof File?new File(data.folder as File,'step4_last_settings.json'):null}
    /** The baseline frames that came from Fiji itself (not from the last recording's stimulus frame): what remembered settings are filed under. */
    static int fijiFrame(Map data,String end){def own=data['fijiBaseline'+end];(own!=null?own:data['baseline'+end]) as int}
    static void remember(Map data,Map draft) {
        File file=rememberedFile(data);if(file==null||!(draft?.values instanceof Map)||!(draft?.protocol instanceof Map))return
        try{file.setText(groovy.json.JsonOutput.toJson([recording:''+data.recording,frames:data.frames as int,baseline_start:fijiFrame(data,'Start'),baseline_end:fijiFrame(data,'End'),
            values:draft.values,protocol:draft.protocol]),'UTF-8')}catch(Exception ignored){}
        rememberStimulusFrame(data,draft)
    }
    static Map remembered(Map data) {
        File file=rememberedFile(data);if(file==null||!file.isFile())return null
        try{def saved=new groovy.json.JsonSlurper().parseText(file.getText('UTF-8'))
            if(!(saved instanceof Map)||saved.recording!=''+data.recording||(saved.frames as int)!=(data.frames as int)||(saved.baseline_start as int)!=fijiFrame(data,'Start')||(saved.baseline_end as int)!=fijiFrame(data,'End'))return null
            if(!(saved.values instanceof Map)||!(saved.protocol instanceof Map))return null
            Map candidate=[values:saved.values,protocol:saved.protocol,timestamps:null,sensitivity:null,negative_control:null]
            return diagnose(candidate,data).errors?null:copy(candidate)
        }catch(Exception ignored){return null}
    }
    /** When blind frames were used: % responding for the cells chosen blind and for the cells added after the response was shown.
     *  Only the blind group is free of selection on the response. The table goes beside the measurement tables; the sentence
     *  returned is shown on the page when the analysis finishes. */
    static String byHowChosen(File dir,Map measurement) {
        try {
            if(!(measurement.chosenBlind instanceof Map))return ''
            Map tags=measurement.chosenBlind as Map
            def calls=CalciumCsv.read(new File(dir,'responder_calls.csv'));List h=calls[0];int cell=h.indexOf('Cell'),role=h.indexOf('Role'),status=h.indexOf('Status')
            if(cell<0||status<0)return ''
            Map tally=[blind:[0,0,0],later:[0,0,0]]   // responded, evaluable, registered
            calls.drop(1).findAll{role<0||it[role]=='Stimulus'}.each{row->
                if(!tags.containsKey(row[cell]))return
                List t=tally[tags[row[cell]]?'blind':'later'];t[2]++
                if(row[status]=='Responded'){t[0]++;t[1]++}else if(row[status]=='Did not respond')t[1]++
            }
            List all=(0..2).collect{int i->(tally.blind[i] as int)+(tally.later[i] as int)}
            def pct={List t->(t[1] as int)>0?String.format(Locale.US,'%.1f',100d*(t[0] as int)/(t[1] as int)):''}
            def line={String label,List t->label+': '+t[0]+' of '+t[1]+' responded'+((t[1] as int)>0?' ('+pct(t)+' %)':'')}
            File out=new File(measurement.folder as File,measurement.recording+'_responders_by_how_chosen_'+dir.name+'.csv')
            CalciumCsv.write(out,[['chosen','registered','evaluable','responded','percent_responding'],['blind']+[tally.blind[2],tally.blind[1],tally.blind[0],pct(tally.blind)],
                ['after the response was shown']+[tally.later[2],tally.later[1],tally.later[0],pct(tally.later)],['all cells']+[all[2],all[1],all[0],pct(all)]])
            if((tally.later[2] as int)==0)return ' All '+tally.blind[2]+' cells were chosen blind: '+tally.blind[0]+' of '+tally.blind[1]+' responded ('+pct(tally.blind)+' %).'
            return ' '+line('Chosen blind',tally.blind)+'. '+line('Added after the response was shown',tally.later)+'. '+line('All cells',all)+'.'
        }catch(Throwable ignored){return ''}
    }
    /** How the page was last arranged (pane sizes, explanations on or off). The page's address changes every time, so the
     *  browser cannot remember these itself; they are kept with Fiji's preferences. */
    static final Map PAGE_PREFS=['step4.side':/\d{2,4}px/,'step4.split':/\d{1,2}(\.\d)?%/,'step4.height':/\d{3,4}px/,'step4.help':/[01]/,
        'step4.every':/\d{1,4}/,'step4.outline':/\d(\.\d)?/,'step4.format':/(png|tiff)/,'step4.size':/(screen|\d{1,2}(\.\d)?x\d{1,2}(\.\d)?)/,'step4.dpi':/(150|300|600)/,'step4.compression':/(lzw|zip|none)/]
    static Map pagePrefs(){PAGE_PREFS.keySet().collectEntries{String k->[k,Prefs.get('calcium.'+k,'')]}.findAll{k,v->v&&(v==~PAGE_PREFS[k])}}
    static Map setPagePref(Object key,Object value) {
        if(!(key instanceof String)||!PAGE_PREFS.containsKey(key))throw new IllegalArgumentException('Unknown page preference.')
        if(value!=null&&!(value instanceof String&&(value==~PAGE_PREFS[key])))throw new IllegalArgumentException('That value cannot be kept.')
        Prefs.set('calcium.'+key,value==null?'':value as String);[kept:key]
    }
    static String full(Object value){value==null?'':''+value}
    static String shown(Object value){(value instanceof Double||value instanceof BigDecimal||value instanceof Float)?CalciumTheme.shortNumber(value as double):full(value)}
    static String identity(Map data) {
        String source=data.csv instanceof File?(data.csv as File).canonicalPath:''
        CalciumInstaller.digest(groovy.json.JsonOutput.toJson([generation:data.measurementId,key:data.key,recording:data.recording,input:source,input_sha256:data.inputHash,registry_sha256:data.registryHash,rois_sha256:data.acceptedRoisHash]).getBytes('UTF-8'))
    }
    static Map draft(Map data) {
        Map values=suggested(data).collectEntries{k,v->[k,full(v)]};double hz=CalciumCsv.number(values.sampling_hz)
        double onset=Double.isFinite(hz)&&hz>0?(values.baseline_end as int)/hz:10d,end=Double.isFinite(hz)&&hz>0?((data.frames as int)-1)/hz:20d
        double from=Double.isFinite(hz)&&hz>0?Math.max(0,(values.baseline_start as int)-1)/hz:0d
        [values:values,protocol:CalciumResponseProtocol.payload([CalciumResponseProtocol.window('Stimulus','Stimulus',from,onset,end)],2,'Slow evoked response — provisional','provisional'),timestamps:null,sensitivity:null,negative_control:null]
    }
    /** Read-only feedback for incomplete drafts; the same authoritative checks also run on submission. */
    static Map diagnose(Map draft,Map measurement) {
        List errors=[],warnings=[]
        try {
            if(!(draft?.values instanceof Map)||!(draft?.protocol instanceof Map))throw new IllegalArgumentException('Complete the recording settings and response protocol.')
            Map flags=[:],uploads=[:]
            FILE_KEYS.each{k,v->uploads[k]=upload(draft[k],k);if(uploads[k])flags[v]=true}
            check(draft.values as Map,draft.protocol as Map,true,measurement,flags)
            double[] times=uploads.timestamps?CalciumMeasurement.readTimesRows(CalciumCsv.parse(uploads.timestamps.content as String),measurement.frames as int):
                ((0..<(measurement.frames as int)).collect{it/CalciumCsv.number(draft.values.sampling_hz.toString())} as double[])
            Map protocol=CalciumResponseProtocol.canonical(draft.protocol)
            protocol.windows.each { original ->
                Map row=new LinkedHashMap(original)
                (protocol.field_overrides?:[]).findAll{it.Field==measurement.recording&&it.Role==row.Role}.each{override->override.each{k,v->if(v!=null)row[k]=v}}
                double onset=(row.Stimulus_onset_s?:row.Start_s) as double,start=row.Start_s as double,end=row.End_s as double
                if(start>=times[-1]||end<=times[0])errors<<row.Stimulus+': the response window must overlap recorded time.'
                else if(end>times[-1])warnings<<row.Stimulus+': response end exceeds the recording; complete-window rules may mark cells unscorable.'
                int last=CalciumCsv.number(draft.values.baseline_end.toString()) as int
                // same sample selection as the engine (responders.R): from ≤ t < to
                double from=(row.Baseline_start_s==null?times[0]:row.Baseline_start_s) as double,to=(row.Baseline_end_s==null?onset:row.Baseline_end_s) as double
                int needed=(row.Min_baseline_samples==null?3:row.Min_baseline_samples) as int,available=0
                double edge=8d*Math.ulp(1d)*Math.max(1d,[from,to,start,end].collect{Math.abs(it as double)}.max() as double)
                for(double t:times)if(t>=from-edge&&t<to-edge)available++
                if(available<needed)warnings<<row.Stimulus+': the response baseline ('+CalciumTheme.shortNumber(from)+'–'+CalciumTheme.shortNumber(to)+' s) holds '+available+' samples but at least '+needed+' are needed, so every cell will be unscorable.'
                if(row.Role=='Stimulus'&&times[last-1]>onset+edge)warnings<<'The trace F₀ baseline includes frames after the stimulus. Choose frames up to the stimulus if this is not intentional.'
            }
        }catch(Exception problem){errors<< (problem.message?:problem.toString())}
        [errors:errors,warnings:warnings]
    }
    /** A finished run where no cell could be scored for the primary stimulus is not a success: says why, from the saved
     *  tables, or returns null when at least one cell was evaluable. Display only; the saved results are unchanged. */
    static String noEvaluableCells(File dir) {
        try {
            def counts=CalciumCsv.read(new File(dir,'responder_counts.csv'));List head=counts[0]
            def row=counts.drop(1).find{head.indexOf('Role')<0||it[head.indexOf('Role')]=='Stimulus'};if(row==null)return null
            double evaluable=CalciumCsv.number(row[head.indexOf('Evaluable')]);if(!(evaluable==0d))return null
            String stimulus=head.indexOf('Stimulus')>=0?row[head.indexOf('Stimulus')]:'the stimulus'
            def calls=CalciumCsv.read(new File(dir,'responder_calls.csv'));List h=calls[0];int role=h.indexOf('Role'),reason=h.indexOf('Reason')
            def reasons=calls.drop(1).findAll{role<0||it[role]=='Stimulus'}.collect{reason>=0?it[reason]:''}.findAll{it}.countBy{it}.sort{-it.value}
            String why=reasons?reasons.keySet().first()+' ('+reasons.values().first()+' cells)':'see the Cell response calls table'
            return 'Finished, but no cell could be scored for '+stimulus+': '+why+'. Check the settings; the results were saved for the record.'
        } catch(Throwable ignored) { return null }
    }
    /** A copy is necessary: neither an incoming HTTP map nor later edits can mutate a submitted run. */
    static Map copy(Map draft){new groovy.json.JsonSlurper().parseText(groovy.json.JsonOutput.toJson(draft)) as Map}
    static String fingerprint(Map draft) {
        def sorted;sorted={v->v instanceof Map?new TreeMap(v.collectEntries{k,x->[''+k,sorted(x)]}):v instanceof Collection?v.collect{sorted(it)}:v}
        CalciumInstaller.digest(groovy.json.JsonOutput.toJson(sorted(draft)).getBytes('UTF-8'))
    }
    /** Acquisition/source checks are shared with the native window. Protocol response math stays in the pinned R engine. */
    static Map check(Map typed,Map protocol,boolean confirmed,Map measurement,Map files=[:]) {
        if(!confirmed)throw new IllegalArgumentException('Confirm the frame rate/timestamps, baseline windows and cell selection first.')
        Map values=VALUE_KEYS.collectEntries{k->[k,full(typed[k]).trim()]}
        double first=CalciumResponseProtocol.numeric(values.baseline_start,'normalization first frame'),last=CalciumResponseProtocol.numeric(values.baseline_end,'normalization last frame')
        if(first<1||first!=Math.floor(first)||last!=Math.floor(last)||last<first||last>=(measurement.frames as int))throw new IllegalArgumentException('Check the normalization baseline frame range.')
        if(files.timestamps_file){values.sampling_hz=''}else{double hz=CalciumResponseProtocol.numeric(values.sampling_hz,'frame rate');if(hz<=0)throw new IllegalArgumentException('Use a positive frame rate.')}
        if(protocol==null)throw new IllegalArgumentException('A portable response protocol is required; choose a version 2 or explicit legacy protocol.')
        Map canonical=CalciumResponseProtocol.canonical(protocol)
        if(canonical.schema_version==1&&(files.sensitivity_file||files.negative_control_file))throw new IllegalArgumentException('Sensitivity and negative-control diagnostics require a version 2 protocol.')
        values.protocol_json=groovy.json.JsonOutput.toJson(canonical);values.recording=measurement.recording;values.input_mode='background_corrected_fluorescence';values.selection_basis=measurement.selectionBasis?:'Not recorded'
        if(measurement.inputHash)values.expected_input_sha256=measurement.inputHash
        if(measurement.registry){values.roi_registry_file=(measurement.registry as File).path;values.roi_registry_sha256=measurement.registryHash}
        if(measurement.acceptedRois){values.roi_set_file=(measurement.acceptedRois as File).path;values.roi_set_sha256=measurement.acceptedRoisHash}
        values.putAll(files)
        values
    }
    /** Optional browser file inputs contain text, never server filesystem paths. */
    static Map upload(Object item,String kind) {
        if(item==null)return null
        if(!(item instanceof Map)||!(item.content instanceof String)||!(item.name instanceof String)||!item.name.trim())throw new IllegalArgumentException('Choose a '+kind+' text file with a name and content.')
        if((item.content as String).getBytes('UTF-8').length>(1<<24))throw new IllegalArgumentException('The '+kind+' file is too large (maximum 16 MiB).')
        [name:item.name.toString(),content:item.content.toString()]
    }
}

/** One measurement and its immutable completed-run display, shared with the local browser. */
class CalciumAnalysisSession {
    Map measurement
    final String measurementIdentity
    CalciumRRuntime runtime=new CalciumRRuntime()
    Closure performAnalysis={rt,input,settings,out,job->CalciumSimpleCaEngine.analyze(rt,input,settings,out,job)}
    Closure open={File target->java.awt.Desktop.getDesktop().open(target)}
    CalciumLocalJob job
    volatile boolean busy=false,installing=false
    volatile int percent=0,version=0
    volatile String stage='',message='Check acquisition timing and a protocol before analysis.',tone='plain'
    volatile long started=0L
    volatile File results
    volatile String resultsJson='{}',appliedFingerprint='',submittedFingerprint=''
    Map draft;long draftRevision=-1L
    String validationFingerprint='';Map draftValidation=[:]
    double[] appliedTime=new double[0]
    String appliedTimeSource=''
    ImagePlus movie;byte[] picture,backgroundPicture;int width=0,height=0
    float[] rest,work,work2,background;int[] tintBuffer;double[] rawRange,changeRange,flatRange,meanTrace=new double[0];double riseTop=1d
    int restFirst=1,restLast=1,frameEpoch=0
    java.awt.image.BufferedImage frameImage
    /** Whether the field moved (CalciumMotion.summary): from Fiji, or checked here when the session was started without it. */
    Map motion
    /** Where the starting stimulus frame came from ('fiji', 'last-recording' or 'last-session'), and Fiji's own values when the last session's were used. */
    String stimulusFrom='fiji';Map fijiDraft
    /** Names what a frame looks like: the baseline behind "Change from baseline" and the display range. The page keys its buffer by it. */
    synchronized String frameVersion(){restFirst+'-'+restLast+'.'+frameEpoch}
    CalciumAnalysisSession(Map data) {
        measurement=new LinkedHashMap(data);measurementIdentity=CalciumAnalysisSettings.identity(data);draft=CalciumAnalysisSettings.draft(data)
        stimulusFrom=data.stimulusFrom?(''+data.stimulusFrom):'fiji'
        if(data.framesNote)message=''+data.framesNote
        Map before=CalciumAnalysisSettings.remembered(data)
        restFirst=(data.baselineStart?:1) as int;restLast=draft.values.baseline_end as int
        if(before!=null){
            double a=CalciumCsv.number(''+before.values.baseline_start),b=CalciumCsv.number(''+before.values.baseline_end)
            if(Double.isFinite(a)&&Double.isFinite(b)&&a==Math.rint(a)&&b==Math.rint(b)&&a>=1&&b>=a&&b<=(data.frames as int)){
                fijiDraft=CalciumAnalysisSettings.copy(draft);stimulusFrom='last-session';draft=before;restFirst=a as int;restLast=b as int;message='Your last settings for this recording are filled in.'}
        }
        if(data.backgroundPicture instanceof ImagePlus){def out=new ByteArrayOutputStream();javax.imageio.ImageIO.write((data.backgroundPicture as ImagePlus).getBufferedImage(),'png',out);backgroundPicture=out.toByteArray()}
        if(data.picture!=null){def image=(data.picture as ImagePlus).getBufferedImage();width=image.width;height=image.height;def out=new ByteArrayOutputStream();javax.imageio.ImageIO.write(image,'png',out);picture=out.toByteArray()}
        // The original recording must match the measured picture, frame count and ROI coordinates. It is read only.
        if(data.movie instanceof ImagePlus&&picture!=null&&(data.movie as ImagePlus).width==width&&(data.movie as ImagePlus).height==height&&(data.movie as ImagePlus).stackSize>=(data.frames as int))movie=data.movie as ImagePlus
        if(data.motion instanceof Map)motion=data.motion as Map
        else if(movie!=null)try{Map chosen=data.fijiFrames instanceof Map?data.fijiFrames as Map:null;motion=CalciumMotion.summary(CalciumMotion.shifts(movie.getStack(),data.frames as int),chosen?chosen.first as int:0,chosen?chosen.last as int:0)}catch(Throwable ignored){motion=null}
    }
    /** The frames (counted from 1) in which the field was displaced. */
    Set<Integer> movedSet(){Set<Integer> out=new HashSet<Integer>();if(motion!=null)(motion.runs as List).each{r->int a=(r as List)[0] as int,b=(r as List)[1] as int;for(int f=a;f<=b;f++)out.add(f)};out}
    static String json(Object value) {
        def clean;clean={v->v instanceof Map?v.collectEntries{k,x->[''+k,clean(x)]}:v instanceof Collection?v.collect{clean(it)}:(v instanceof double[]||v instanceof float[]||v instanceof int[])?(v as List).collect{clean(it)}:(v instanceof Double||v instanceof Float)?(Double.isFinite(v as double)?v:null):v instanceof File?v.name:v}
        groovy.json.JsonOutput.toJson(clean(value))
    }
    /** Before analysis the recording is scrubbable using preview timing; saved results always keep their applied clock. */
    synchronized Map timing() {
        if(results!=null)return [frameTimes:appliedTime,timeSource:appliedTimeSource,timingProvisional:false]
        try {
            if(draft.timestamps!=null){Map item=CalciumAnalysisSettings.upload(draft.timestamps,'timestamps');double[] times=CalciumMeasurement.readTimesRows(CalciumCsv.parse(item.content as String),measurement.frames as int);return [frameTimes:times,timeSource:'explicit_timestamps',timingProvisional:true]}
            double hz=CalciumCsv.number(draft.values?.sampling_hz?.toString())
            if(Double.isFinite(hz)&&hz>0)return [frameTimes:(0..<(measurement.frames as int)).collect{it/hz},timeSource:'constant_fps',timingProvisional:true]
        }catch(Throwable ignored){}
        [frameTimes:null,timeSource:'unknown',timingProvisional:true]
    }
    synchronized Map described() {
        [recording:''+(measurement.recording?:''),key:measurementIdentity,frames:measurement.frames as int,names:(measurement.names?:[]) as List,width:width,height:height,
         hasPicture:picture!=null,hasBackground:backgroundPicture!=null&&picture!=null,hasMovie:movie!=null,pluginVersion:CalciumRelease.VERSION,prefs:CalciumAnalysisSettings.pagePrefs(),
         recordedHz:(measurement.hz instanceof Number&&Double.isFinite(measurement.hz as double)&&(measurement.hz as double)>0)?measurement.hz:null,fijiFrames:measurement.fijiFrames instanceof Map?measurement.fijiFrames:null,stimulusFrom:stimulusFrom,fijiDraft:fijiDraft,motion:motion,motionNote:CalciumMotion.sentence(motion),selectionBasis:measurement.selectionBasis?:'Not recorded',
         rois:((measurement.rois?:new Roi[0]) as List).collect{Roi roi->def p=roi.getPolygon();(0..<p.npoints).collect{int i->[p.xpoints[i],p.ypoints[i]]}},
         draft:CalciumAnalysisSettings.copy(draft),draftRevision:draftRevision,
         settings:CalciumAnalysisSettings.suggested(measurement).collect{k,v->[key:k,shown:CalciumAnalysisSettings.shown(v),exact:CalciumAnalysisSettings.full(v)]}]+timing()
    }
    synchronized Map state() {
        boolean ready=false;try{ready=runtime.environment()!=null}catch(Throwable ignored){}
        String fingerprint=CalciumAnalysisSettings.fingerprint(draft)
        [key:measurementIdentity,busy:busy,installing:installing,percent:percent,stage:stage,message:message,tone:tone,elapsed:busy?(long)((System.currentTimeMillis()-started)/1000L):0L,
         version:version,folder:results?.name,engineReady:ready,logAvailable:job?.logFile?.isFile()?:false,rest:[restFirst,restLast],frameVersion:frameVersion(),
         settingsPending:results!=null&&fingerprint!=appliedFingerprint,draftFingerprint:fingerprint,appliedFingerprint:appliedFingerprint,draftRevision:draftRevision,validation:validation(fingerprint)]+timing()
    }
    synchronized Map validation(String fingerprint) {
        if(validationFingerprint!=fingerprint){draftValidation=CalciumAnalysisSettings.diagnose(draft,measurement);validationFingerprint=fingerprint}
        draftValidation
    }
    /** Drafts may contain incomplete text. They receive feedback without running analysis or writing uploads. */
    synchronized Map updateDraft(Map value,Object revision=null) {
        if(value==null)throw new IllegalArgumentException('Provide the pending draft settings.')
        if(revision!=null){if(!(revision instanceof Number)||!Double.isFinite(revision as double)||(revision as double)!=(revision as long)||revision<0)throw new IllegalArgumentException('Draft revision must be a nonnegative integer.');if((revision as long)<draftRevision)return state();draftRevision=revision as long}
        draft=CalciumAnalysisSettings.copy(value)
        Map now=state()
        if(!(draftValidation?.errors))CalciumAnalysisSettings.remember(measurement,draft)   // kept for the next time this recording is opened
        now
    }
    void launch(boolean setup,Closure action,String fingerprint=null,Map submitted=null) {
        synchronized(this){if(busy)throw new IllegalStateException('A job is already running.');busy=true;if(fingerprint!=null)submittedFingerprint=fingerprint}
        installing=setup;started=System.currentTimeMillis();percent=0;stage=setup?'Preparing R setup':'Preparing analysis';message=stage;tone='plain'
        final CalciumLocalJob mine=new CalciumLocalJob(logFile:new File(runtime.home,setup?'setup.log':'analysis-'+measurementIdentity.substring(0,12)+'.log'));job=mine
        mine.progress={int p,String text->percent=p;stage=text}
        Thread.startDaemon('SimpleCa analysis for the Step 4 page') {
            try {File output=action.call(mine) as File
                if(setup){percent=100;message='R is ready. Check the settings, then run analysis.';tone='good'}
                else {load(output,fingerprint,submitted);percent=100;String empty=CalciumAnalysisSettings.noEvaluableCells(output)
                    if(empty){message=empty;tone='warn'} else {message='Complete. Protocol, input snapshots and tables saved in '+output.name+'.'+CalciumAnalysisSettings.byHowChosen(output,measurement);tone='good'}}
            }catch(java.util.concurrent.CancellationException ignored){message='Cancelled. Completed analyses are kept; partial output is not used.';tone='plain';stage='';percent=0}
            catch(Throwable t){message=t.message?:t.toString();tone='bad'}
            finally{busy=false}
        }
    }
    /** Saves what the page shows (the picture, the numbers behind it, and how it was made) in a "plots" folder beside the
     *  measurement tables, so nothing has to be fished out of the browser's Downloads. Only that folder is written. */
    Map export(Map body) {
        if(!(measurement.folder instanceof File))throw new IllegalStateException('This measurement has no folder to save into.')
        // a set of pictures (the frames of a movie view) goes into a folder of its own inside "plots"
        String sub=body.folder instanceof String?(body.folder as String).replaceAll(/[^A-Za-z0-9._-]+/,'_').replaceAll(/^[._]+/,''):''
        File dir=sub.isEmpty()?plotsFolder():new File(plotsFolder(),sub.length()>120?sub.substring(0,120):sub);if(!dir.isDirectory()&&!dir.mkdirs())throw new IOException('Could not create the folder '+dir.path)
        String stem=(''+(body.name?:'plot')).replaceAll(/[^A-Za-z0-9._-]+/,'_').replaceAll(/^[._]+/,'');if(stem.isEmpty())stem='plot';if(stem.length()>160)stem=stem.substring(0,160)
        List saved=[]
        if(body.png instanceof String){
            List picture=CalciumPictureFile.fromPage(body.png as String,body.format,body.dpi,body.compression)
            new File(dir,stem+'.'+picture[1]).bytes=picture[0] as byte[];saved<<stem+'.'+picture[1]
        }
        if(body.csv instanceof String&&!(body.csv as String).isEmpty()){new File(dir,stem+'.csv').setText(body.csv as String,'UTF-8');saved<<stem+'.csv'}
        if(body.metadata instanceof Map){new File(dir,stem+'.metadata.json').setText(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(body.metadata)),'UTF-8');saved<<stem+'.metadata.json'}
        if(saved.isEmpty())throw new IllegalArgumentException('There was nothing to save.')
        [folder:dir.path,folderName:dir.name,saved:saved]
    }
    File plotsFolder(){new File(measurement.folder as File,'plots')}
    void setup(){launch(true){task->runtime.setup(task)}}
    /** Validate an exact copy, write private uploads, then let the engine snapshot them into the completed run. */
    synchronized void analyze(Map incoming,boolean confirmed) {
        if(busy)throw new IllegalStateException('A job is already running.')
        if(!confirmed)throw new IllegalArgumentException('Confirm the frame rate/timestamps, baseline windows and cell selection first.')
        if(!(incoming instanceof Map)||!(incoming.values instanceof Map)||!(incoming.protocol instanceof Map))throw new IllegalArgumentException('A complete draft with a portable response protocol is required.')
        Map submitted=CalciumAnalysisSettings.copy(incoming),uploads=[:]
        Map validation=CalciumAnalysisSettings.diagnose(submitted,measurement)
        if(validation.errors)throw new IllegalArgumentException(validation.errors.join(' '))
        CalciumAnalysisSettings.remember(measurement,submitted)
        CalciumAnalysisSettings.FILE_KEYS.each{k,v->uploads[k]=CalciumAnalysisSettings.upload(submitted[k],k)}
        if(uploads.timestamps)CalciumMeasurement.readTimesRows(CalciumCsv.parse(uploads.timestamps.content as String),measurement.frames as int)
        // Check everything independent of filesystem writes before creating the request snapshot.
        Map flags=[:];CalciumAnalysisSettings.FILE_KEYS.each{k,v->if(uploads[k])flags[v]=true}
        Map checked=CalciumAnalysisSettings.check(submitted.values as Map,submitted.protocol as Map,confirmed,measurement,flags)
        File request=new File(runtime.home,'requests/'+UUID.randomUUID().toString());if(!request.mkdirs())throw new IOException('Could not create the analysis input snapshot folder.')
        try {
            CalciumAnalysisSettings.FILE_KEYS.each{k,v->if(uploads[k]){
                File file=new File(request,k+'.csv');file.setText(uploads[k].content as String,'UTF-8');checked[v]=file.path
                checked[k=='timestamps'?'timestamps_sha256':v+'_sha256']=CalciumInstaller.digest(file.bytes)
            }}
            String fingerprint=CalciumAnalysisSettings.fingerprint(submitted)
            checked.draft_sha256=fingerprint
            draft=submitted
            launch(false,{task->try{performAnalysis.call(runtime,measurement.csv,checked,measurement.folder,task)}finally{request.deleteDir()}},fingerprint,submitted)
        }catch(Throwable t){request.deleteDir();throw t}
    }
    void cancel(){job?.cancel();if(busy)stage='Cancelling...'}
    static Map table(File dir,String name){File file=new File(dir,name);if(!file.isFile())return [header:[],rows:[]];def rows=CalciumCsv.read(file);[header:rows?rows[0]:[],rows:rows?rows.drop(1):[]]}
    static List<Map> records(Map table){(table.rows as List).collect{line->table.header.withIndex().collectEntries{h,i->[h,i<line.size()?line[i]:'']}}}
    /** Draw only completed files. No response classifications or scientific metrics are recomputed in Java. */
    void load(File dir,String fingerprint=null,Map submitted=null) {
        if(!new File(dir,'COMPLETE').isFile())throw new IOException('This analysis has no COMPLETE marker.')
        def normalized=CalciumCsv.read(new File(dir,'normalized_traces.csv')),average=CalciumCsv.read(new File(dir,'mean_trace.csv'))
        def numbers={List rows,int column->rows.drop(1).collect{CalciumCsv.number(it[column])}}
        List time=numbers(normalized,0);if(time.size()!=(measurement.frames as int))throw new IOException('Saved timestamps do not match the measured recording frames.')
        for(int i=0;i<time.size();i++)if(!Double.isFinite(time[i] as double)||time[i]<0||(i>0&&time[i]<=time[i-1]))throw new IOException('Saved timestamps are invalid.')
        Map traces=[:];normalized[0].drop(1).eachWithIndex{name,i->traces[''+name]=numbers(normalized,i+1)}
        Map manifest=CalciumCsv.read(new File(dir,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
        int first=CalciumResponseProtocol.numeric(manifest.baseline_start,'saved baseline first frame') as int,last=CalciumResponseProtocol.numeric(manifest.baseline_end,'saved baseline last frame') as int
        if(first<1||last<first||last>time.size())throw new IOException('Saved baseline frames do not match the recording.')
        Map tables=[:]
        [metrics:'metrics.csv',calls:'responder_calls.csv',counts:'responder_counts.csv',quality:'cell_quality.csv',registry:'roi_registry.csv',windowMetrics:'window_metrics.csv',replicates:'replicate_counts.csv',groups:'group_counts.csv',unassigned:'unassigned_biological_units.csv',cohort:'cohort_flow.csv',recordingMetadata:'recording_metadata.csv',sensitivityCalls:'sensitivity_calls.csv',sensitivityCounts:'sensitivity_counts.csv',negativeControlCalls:'negative_control_calls.csv',negativeControlCounts:'negative_control_counts.csv'].each{k,file->tables[k]=table(dir,file)}
        Map callsByWindow=[:],countsByWindow=[:];List windows=[]
        records(tables.calls).each{row->String role=row.Role?:'Stimulus',stimulus=row.Stimulus?:'Stimulus',key=role+'::'+stimulus
            if(!callsByWindow.containsKey(key)){callsByWindow[key]=[:];windows<<[key:key,role:role,stimulus:stimulus]};callsByWindow[key][row.Cell]=row}
        records(tables.counts).each{row->String key=(row.Role?:'Stimulus')+'::'+(row.Stimulus?:'Stimulus');if(!countsByWindow.containsKey(key))countsByWindow[key]=row}
        Map primary=windows.find{it.role=='Stimulus'}?:windows.find{true}?:[:];Map calls=callsByWindow[primary.key]?:[:],counts=countsByWindow[primary.key]?:[:]
        Map firstCall=calls.values().find{true}?:[:]
        def protocol={String name->File file=new File(dir,name);file.isFile()?new groovy.json.JsonSlurper().parse(file,'UTF-8'):null}
        String applied=fingerprint?:manifest.draft_sha256?:''
        List mean=numbers(average,1)
        // The plot can also show the gray values themselves: corrected (cell minus its background) and raw (not corrected).
        // They come from the measurement's own tables and are offered only when those hold exactly the analysed cells and frames.
        def grayTable={File file->Map out=[:]
            try{if(file!=null&&file.isFile()){def rows=CalciumCsv.read(file)
                if(rows.size()==time.size()+1&&(rows[0].drop(1).collect{''+it} as Set).containsAll(traces.keySet())&&(rows[0].drop(1).collect{''+it} as Set)==((measurement.names?:[]).collect{''+it} as Set))rows[0].drop(1).eachWithIndex{name,i->out[''+name]=numbers(rows,i+1)}}}catch(Throwable ignored){out=[:]}
            out}
        File tableFolder=measurement.folder instanceof File?measurement.folder as File:null
        Map corrected=grayTable(measurement.csv instanceof File?measurement.csv as File:null)
        Map raw=grayTable(tableFolder==null?null:new File(tableFolder,measurement.recording+'_gray_values_raw.csv'))
        Map loaded=[folder:dir.name,time:time,frameTimes:time,timeSource:manifest.time_source?:'constant_fps',timingProvisional:false,mean:mean,sem:numbers(average,2),traces:traces,calls:calls,counts:counts,windows:windows,callsByWindow:callsByWindow,countsByWindow:countsByWindow,
            baseline:[time[first-1],time[last-1]],baselineFrames:[first,last],stimulus:[CalciumCsv.number(firstCall.Start_s?:manifest.stimulus_start),CalciumCsv.number(firstCall.Effective_end_s?:firstCall.End_s?:manifest.stimulus_end)],
            corrected:corrected,raw:raw,chosenBlind:(measurement.chosenBlind instanceof Map?measurement.chosenBlind:null),appliedDraft:submitted==null?null:CalciumAnalysisSettings.copy(submitted),protocol:protocol('protocol.json'),appliedProtocol:protocol('applied_protocol.json'),manifest:manifest,appliedFingerprint:applied,tables:tables]
        String text=json(loaded)
        synchronized(this){meanTrace=mean as double[];restFirst=first;restLast=last;rest=null;frameEpoch++;appliedTime=time as double[];appliedTimeSource=loaded.timeSource;appliedFingerprint=applied;results=dir;resultsJson=text;version++}
    }
    /** One frame of the recording as a JPEG: as recorded, or minus the mean of the baseline frames and lightly smoothed. Display only.
     *  A width below the recording's gives a smaller copy for playback; it is painted at full size first, so it looks the same. */
    byte[] frame(int wanted,String mode,int size=0) {
        java.awt.image.BufferedImage copy
        synchronized(this){copy=scaled(paintFrame(wanted,mode),size)}
        // encoding is the slow part and needs no shared state, so several frames can be encoded at once
        def out=new ByteArrayOutputStream(1<<16);javax.imageio.ImageIO.write(copy,'jpg',out);return out.toByteArray()
    }
    static java.awt.image.BufferedImage scaled(java.awt.image.BufferedImage image,int size) {
        int w=image.width,h=image.height
        if(size<=0||size>=w){def copy=new java.awt.image.BufferedImage(w,h,java.awt.image.BufferedImage.TYPE_INT_RGB);image.copyData(copy.getRaster());return copy}
        def current=image
        while(true){   // halve in steps, then the remainder, so fine detail is averaged rather than skipped
            int next=Math.max(size,(current.width/2) as int),nh=Math.max(1,Math.round(h*next/(double)w) as int)
            def smaller=new java.awt.image.BufferedImage(next,nh,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=smaller.createGraphics()
            g.setRenderingHint(java.awt.RenderingHints.KEY_INTERPOLATION,java.awt.RenderingHints.VALUE_INTERPOLATION_BILINEAR);g.drawImage(current,0,0,next,nh,null);g.dispose()
            current=smaller;if(next==size)return current
        }
    }
    private java.awt.image.BufferedImage paintFrame(int wanted,String mode) {
        if(movie==null)throw new FileNotFoundException('No recording is attached to this measurement.')
        int n=measurement.frames as int
        if(wanted<1||wanted>n)throw new IllegalArgumentException('Frame must be between 1 and '+n+'.')
        if(!(mode in ['raw','flat','magenta','change']))throw new IllegalArgumentException('Frame mode must be raw, flat, magenta or change.')
        def stack=movie.getStack();int f=wanted
        if(rest==null) {
            rest=CalciumFrameView.mean(stack,restFirst,restLast,movedSet())   // moved frames would smear the resting picture
            // the frame where the cells are brightest on average sets the display range, so it does not flicker from frame to frame
            int peak=Math.max(1,Math.min(n,(restLast+n).intdiv(2)));double top=Double.NEGATIVE_INFINITY
            meanTrace.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top&&i<n){top=v;peak=i+1}}
            def brightest=stack.getProcessor(peak);def resting=new ij.process.FloatProcessor(width,height,rest)
            double rise=CalciumFrameView.level(brightest,rest,0.9995d);if(!(rise>0d))rise=1d
            // What a resting frame shows by itself (the last baseline frame against the resting picture): camera noise. A change is
            // drawn at full strength only from twice that, so a recording in which nothing responds stays dark instead of showing
            // its noise as if it were a response.
            def quiet=stack.getProcessor(Math.max(1,Math.min(n,restLast)));double quietMiddle=CalciumFrameView.level(quiet,rest,0.5d)
            double quietSpread=CalciumFrameView.level(quiet,rest,0.99d)-quietMiddle,quietTop=CalciumFrameView.level(quiet,rest,0.9995d)-quietMiddle
            rise=Math.max(rise,2d*quietTop)
            changeRange=[-0.1d*rise,rise] as double[]
            double low=CalciumFrameView.level(resting,null,0.005d),high=Math.max(CalciumFrameView.level(resting,null,0.9995d),CalciumFrameView.level(brightest,null,0.9995d));if(!(high>low))high=low+1d
            rawRange=[low,high] as double[]
            // "Background removed" and "response in magenta": the smooth background of the resting picture (Process > Subtract
            // Background, rolling ball, as for the plugin's all-frames picture) is taken off every frame, so the lamp's uneven
            // light goes and the cells stand out. Shown from 0 to the level 97 % of the brightest frame's pixels stay under, so
            // very bright objects do not set the range. Magenta is the rise above the resting picture. The whole field fades
            // during a recording, so each frame's own middle value (the median over the picture) is taken as "no change"; magenta
            // starts a little above it and is full at 0.6 of the level 99 % of the brightest frame's pixels stay under. Display only.
            def ground=new ij.process.FloatProcessor(width,height,(float[])rest.clone())
            new ij.plugin.filter.BackgroundSubtracter().rollingBallBackground(ground,CalciumSteps.BALL_RADIUS,true,false,false,true,true)
            background=(float[])ground.getPixels()
            double flatTop=CalciumFrameView.level(brightest,background,0.97d);if(!(flatTop>0d))flatTop=1d
            flatRange=[0d,flatTop] as double[]
            riseTop=CalciumFrameView.level(brightest,rest,0.99d)-CalciumFrameView.level(brightest,rest,0.5d);if(!(riseTop>0d))riseTop=rise
            riseTop=Math.max(riseTop,2d*quietSpread)
            frameImage=new java.awt.image.BufferedImage(width,height,java.awt.image.BufferedImage.TYPE_INT_RGB);work=new float[width*height];work2=new float[width*height]
        }
        int[] pixels=((java.awt.image.DataBufferInt)frameImage.getRaster().getDataBuffer()).getData()
        if(mode=='flat'){CalciumFrameView.paint(stack.getProcessor(f),background,flatRange[0],flatRange[1],pixels,work,work2,false);return frameImage}
        if(mode=='magenta'){
            if(tintBuffer==null||tintBuffer.length!=pixels.length)tintBuffer=new int[pixels.length]
            double middle=CalciumFrameView.level(stack.getProcessor(f),rest,0.5d)
            CalciumFrameView.paint(stack.getProcessor(f),rest,middle+0.15d*riseTop,middle+0.6d*riseTop,tintBuffer,work,work2,true,(float)middle)
            CalciumFrameView.paint(stack.getProcessor(f),background,flatRange[0],flatRange[1],pixels,work,work2,false)
            CalciumFrameView.tint(pixels,tintBuffer);return frameImage
        }
        boolean change=mode!='raw';double[] range=change?changeRange:rawRange
        // "Change from baseline" counts from each frame's own middle value too (the whole field fades), and loses the rims of shifted objects
        if(change){double middle=CalciumFrameView.level(stack.getProcessor(f),rest,0.5d);CalciumFrameView.paint(stack.getProcessor(f),rest,middle+range[0],middle+range[1],pixels,work,work2,true,(float)middle)}
        else CalciumFrameView.paint(stack.getProcessor(f),null,range[0],range[1],pixels,work,work2,false)
        frameImage
    }
}

/** Serves the Step 4 page to a browser on this computer only. Nothing is reachable from other computers, and every
 *  address carries a random key, so no other program or web page can read or drive the analysis. */
class CalciumWebServer {
    volatile CalciumAnalysisSession session
    def server; def executor; int port=0
    boolean closed=false
    Closure onStop={ -> }
    String token=UUID.randomUUID().toString().replace('-','')+UUID.randomUUID().toString().replace('-','')
    Map<String,byte[]> files=[:]
    /** The page's own files, from the copy embedded in this plugin (or from a folder while developing the page). */
    static Map<String,byte[]> page() {
        Map<String,byte[]> out=[:]
        String folder=System.getProperty('calcium.web.dir')
        if(folder){new File(folder).listFiles().findAll{it.isFile()}.each{out[it.name]=it.bytes};return out}
        byte[] data=Base64.getDecoder().decode(CalciumWebBundle.CONTENT)
        if(CalciumInstaller.digest(data)!=CalciumWebBundle.SHA256)throw new IOException('The embedded Step 4 page is damaged; reinstall the plugin.')
        def zip=new java.util.zip.ZipInputStream(new ByteArrayInputStream(data))
        try {def entry;while((entry=zip.getNextEntry())!=null)if(!entry.isDirectory())out[entry.getName()]=read(zip,1<<24)} finally {zip.close()}
        return out
    }
    /** Reads to the end of a stream (or of one ZIP entry) without closing it, up to a limit. */
    static byte[] read(InputStream input,int limit) {
        def out=new ByteArrayOutputStream();byte[] buffer=new byte[8192];int count
        while((count=input.read(buffer))>0){out.write(buffer,0,count);if(out.size()>limit)throw new IOException('Too much data')}
        return out.toByteArray()
    }
    /** Release the controls only after the browser launch succeeds, so they cannot cover the review page. */
    static void openBrowser(java.awt.Window controls,java.net.URI address,Closure browse={java.net.URI target->java.awt.Desktop.desktop.browse(target)}){browse.call(address);controls.setAlwaysOnTop(false);controls.toBack()}
    String url(){'http://127.0.0.1:'+port+'/'+token+'/'}
    /** Throws when this Java cannot run the small built-in web server; the caller then opens the Fiji window instead. */
    void start() {
        files=page()
        // Looked up by name, so the plugin still loads on a Java that leaves this optional part out.
        def serverClass=Class.forName('com.sun.net.httpserver.HttpServer'),handlerClass=Class.forName('com.sun.net.httpserver.HttpHandler')
        server=serverClass.create(new InetSocketAddress(InetAddress.getLoopbackAddress(),0),0)
        server.createContext('/',{exchange->handle(exchange)}.asType(handlerClass))
        executor=java.util.concurrent.Executors.newCachedThreadPool({Runnable r->def t=new Thread(r,'Calcium Step 4 page');t.setDaemon(true);t} as java.util.concurrent.ThreadFactory)
        server.setExecutor(executor)
        server.start();port=server.getAddress().getPort()
    }
    synchronized void stop(){if(closed)return;closed=true;try{server?.stop(0);executor?.shutdownNow()}finally{onStop.call()}}
    static final Map TYPES=[html:'text/html; charset=utf-8',css:'text/css; charset=utf-8',js:'text/javascript; charset=utf-8',json:'application/json; charset=utf-8',png:'image/png',jpg:'image/jpeg',txt:'text/plain; charset=utf-8']
    void reply(exchange,int status,String type,byte[] body,String cache='no-store') {
        exchange.getResponseHeaders().set('Content-Type',TYPES[type]?:type);exchange.getResponseHeaders().set('Cache-Control',cache);exchange.getResponseHeaders().set('X-Content-Type-Options','nosniff')
        exchange.sendResponseHeaders(status,body.length==0?-1L:(long)body.length)
        if(body.length>0){def out=exchange.getResponseBody();out.write(body);out.close()}
    }
    void replyJson(exchange,int status,Object value){reply(exchange,status,'json',(value instanceof String?value:CalciumAnalysisSession.json(value)).getBytes('UTF-8'))}
    void handle(exchange) {
        try {
            String host=exchange.getRequestHeaders().getFirst('Host')?:'',path=exchange.getRequestURI().getPath(),prefix='/'+token+'/'
            if(!(host.startsWith('127.0.0.1:')||host.startsWith('localhost:'))||!path.startsWith(prefix)){reply(exchange,404,'txt','Not found'.bytes);return}
            String what=path.substring(prefix.length()),method=exchange.getRequestMethod();def now=session
            Map query=[:];(exchange.getRequestURI().getRawQuery()?:'').split('&').each{String pair->int cut=pair.indexOf('=');if(cut>0)query[pair.substring(0,cut)]=URLDecoder.decode(pair.substring(cut+1),'UTF-8')}
            if(method=='GET') {
                if(what==''||what=='index.html'){reply(exchange,200,'html',files['step4.html']);return}
                if(what=='app.css'){reply(exchange,200,'css',files['step4.css']);return}
                if(what=='app.js'){reply(exchange,200,'js',files['step4.js']);return}
                if(now==null){replyJson(exchange,409,[error:'Measure and save in Fiji first.']);return}
                switch(what) {
                    case 'api/session':replyJson(exchange,200,now.described());return
                    case 'api/state':replyJson(exchange,200,now.state());return
                    case 'api/results':replyJson(exchange,200,now.resultsJson);return
                    case 'api/picture':if(now.picture==null)reply(exchange,404,'txt','No picture'.bytes) else reply(exchange,200,'png',now.picture);return
                    case 'api/background':if(now.backgroundPicture==null)reply(exchange,404,'txt','No background picture'.bytes) else reply(exchange,200,'png',now.backgroundPicture);return
                    case 'api/frame':
                        if(!query.containsKey('f')||!query.containsKey('mode'))throw new IllegalArgumentException('Ask for a frame number and a display mode.')
                        String frame=query.f,size=query.w?:'0'
                        if(!(frame==~ /[0-9]+/))throw new IllegalArgumentException('Frame must be a whole number.')
                        if(!(size==~ /[0-9]{1,5}/))throw new IllegalArgumentException('Width must be a whole number of pixels.')
                        String version=now.frameVersion();byte[] image=now.frame(frame as int,query.mode,size as int)
                        // a frame of a given version never changes, so the browser may keep it; an outdated request is not kept
                        reply(exchange,200,'jpg',image,query.v==version?'private, max-age=3600':'no-store');return
                }
            } else if(method=='POST'&&now!=null) {
                byte[] sent=read(exchange.getRequestBody(),1<<26)
                Map body=sent.length>0?(new groovy.json.JsonSlurper().parseText(new String(sent,'UTF-8')) as Map):[:]
                switch(what) {
                    case 'api/draft':replyJson(exchange,200,now.updateDraft(body.draft as Map,body.revision));return
                    case 'api/protocol':replyJson(exchange,200,[protocol:CalciumResponseProtocol.canonical(body.protocol)]);return
                    case 'api/analyze':now.analyze(body.draft as Map,body.confirmed==true);replyJson(exchange,200,now.state());return
                    case 'api/setup':now.setup();replyJson(exchange,200,now.state());return
                    case 'api/cancel':now.cancel();replyJson(exchange,200,now.state());return
                    case 'api/export':replyJson(exchange,200,now.export(body));return
                    case 'api/prefs':replyJson(exchange,200,CalciumAnalysisSettings.setPagePref(body.key,body.value));return
                    case 'api/convert':   // the same file as Save writes, handed back for the browser's own download
                        if(!(body.png instanceof String))throw new IllegalArgumentException('There is no picture to convert.')
                        List converted=CalciumPictureFile.fromPage(body.png as String,body.format,body.dpi,body.compression)
                        replyJson(exchange,200,[extension:converted[1],base64:Base64.encoder.encodeToString(converted[0] as byte[])]);return
                    case 'api/close':
                        if(now.busy)throw new IllegalStateException('Cancel or finish the analysis before closing this page.')
                        replyJson(exchange,200,[closed:true]);Thread.startDaemon('Close calcium analysis page'){stop()};return
                    case 'api/open':
                        File target=body.what=='log'?now.job?.logFile:body.what=='plots'?(now.measurement.folder instanceof File?now.plotsFolder():null):now.results
                        if(target==null||!target.exists()){replyJson(exchange,404,[error:body.what=='log'?'There is no log yet.':body.what=='plots'?'No plot has been saved yet.':'There are no saved results yet.']);return}
                        now.open.call(target);replyJson(exchange,200,[opened:target.name]);return
                }
            }
            reply(exchange,404,'txt','Not found'.bytes)
        } catch(IllegalArgumentException|IllegalStateException|IOException refused){try{replyJson(exchange,400,[error:refused.message?:'That is not possible right now.'])}catch(Throwable ignored){}}
        catch(Throwable problem){try{replyJson(exchange,500,[error:problem.message?:problem.toString()])}catch(Throwable ignored){}}
        finally {exchange.close()}
    }
}

/** Separate browser sessions for separate measurements. Closing the selection window does not end review. */
class CalciumWebPages {
    Map<String,CalciumWebServer> pages=[:]
    File snapshotRoot=new File(new CalciumRRuntime().home,'browser-sessions')
    synchronized CalciumWebServer page(Map measurement) {
        String identity=CalciumAnalysisSettings.identity(measurement)
        if(pages[identity]!=null&&!pages[identity].closed)return pages[identity]
        File folder=new File(snapshotRoot,UUID.randomUUID().toString())
        if(!folder.mkdirs())throw new IOException('Could not preserve the measurement for browser review.')
        folder.deleteOnExit()
        try {
            Map captured=new LinkedHashMap(measurement)
            ['csv','registry','acceptedRois'].each{key->if(measurement[key] instanceof File){
                File source=measurement[key] as File,target=new File(folder,source.name)
                java.nio.file.Files.copy(source.toPath(),target.toPath());target.deleteOnExit();captured[key]=target
            }}
            def opened=new CalciumWebServer(session:new CalciumAnalysisSession(captured))
            opened.onStop={ -> synchronized(this){pages.remove(identity)};folder.deleteDir() }
            opened.start();pages[identity]=opened;opened
        }catch(Throwable problem){folder.deleteDir();throw problem}
    }
    void closeAll(){List active;synchronized(this){active=new ArrayList(pages.values())};active.each{it.session.cancel();it.stop()}}
}

// BEGIN GENERATED SIMPLECA BUNDLE
class CalciumSimpleCaBundle {
    static final String HASH = "5ef618cd736e4736b8b3500dd971b44bdefda7fa564609b48afadb6067f0bc27"
    static final String DATA = "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"
}
// END GENERATED SIMPLECA BUNDLE

// BEGIN GENERATED WEB BUNDLE
class CalciumWebBundle {
    static final String SHA256 = "b58bce179731c5cdfe12ceaf3ee2a8ddf6e041ba831cdacad47bf35aac9680f9"
    static final List<String> PARTS = [
        "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"]
    static final String CONTENT = PARTS.join("")
}
// END GENERATED WEB BUNDLE

/** Cell edits keep ROI order/names, and can be undone without overwriting outside edits. */
class CalciumCellEdits {
    RoiManager manager
    ImagePlus picture
    List history = []

    List snapshot() {
        // Assign identities inside the edit transaction. A later display refresh must
        // not change its fingerprint and make an ordinary addition impossible to undo.
        Roi[] original=manager.getRoisAsArray(),stable=original.collect{it} as Roi[]
        CalciumMeasurement.ensureIds(stable)
        for(int i=0;i<original.length;i++)if(!original[i].is(stable[i]))manager.setRoi(stable[i],i)
        (0..<manager.getCount()).collect { int i ->
            Roi copy=new ij.io.RoiDecoder(RoiEncoder.saveAsByteArray(manager.getRoi(i)),'snapshot.roi').getRoi()
            [roi:copy,name:manager.getName(i)]
        }
    }
    List fingerprint(List state) {
        state.collect { [it.name, Base64.getEncoder().encodeToString(RoiEncoder.saveAsByteArray(it.roi))] }
    }
    void change(String label, Closure action) {
        def before = snapshot()
        if (!history.isEmpty() && fingerprint(before) != history.last().after) history.clear()
        action.call()
        def after = fingerprint(snapshot())
        if (fingerprint(before) != after) {
            history << [before: before, after: after, label: label]
            if (history.size() > 20) history.remove(0)
        }
    }
    int hit(int x, int y) {
        def hits = (0..<manager.getCount()).findAll { int i ->
            def r = manager.getRoi(i); r.isArea() && r.contains(x, y)
        }
        hits.isEmpty() ? -1 : hits.min { int i ->
            def b = manager.getRoi(i).getBounds(); b.getWidth() * b.getHeight()
        } as int
    }
    boolean selectAt(int x, int y) {
        int i = hit(x, y)
        if (i < 0) { manager.deselect(); picture.deleteRoi(); return false }
        manager.select(picture, i)
        // The plugin's overlay supplies the highlight; do not leave ImageJ's yellow editing handles on top.
        picture.deleteRoi()
        true
    }
    int remove(List indexes) {
        def chosen = indexes.unique().findAll { int i -> i >= 0 && i < manager.getCount() && manager.getRoi(i).isArea() }.sort()
        if (chosen.isEmpty()) return 0
        change('remove ' + chosen.size() + (chosen.size() == 1 ? ' cell' : ' cells')) {
            chosen.reverseEach { int i -> manager.delete(i) }
            manager.deselect(); picture.deleteRoi()
        }
        chosen.size()
    }
    String undo() {
        if (history.isEmpty()) return 'No cell changes to undo.'
        def entry = history.last()
        if (fingerprint(snapshot()) != entry.after) {
            history.clear()
            return 'The ROI Manager changed outside this window. Undo was cleared to keep those changes.'
        }
        history.remove(history.size() - 1)
        manager.reset()
        entry.before.each { saved ->
            manager.addRoi((Roi) saved.roi.clone())
            manager.rename(manager.getCount() - 1, saved.name)
        }
        manager.deselect(); picture.deleteRoi()
        'Undid: ' + entry.label + '.'
    }
}

/** White outlines with a black border remain visible over bright, dark and coloured image regions. */
class CalciumCellAppearance {
    static void add(Overlay overlay, Roi cell, boolean edge, boolean selected) {
        [Color.black, Color.white].eachWithIndex { Color colour, int layer ->
            def outline = (Roi) cell.clone()
            float width = (float) (layer == 0 ? (selected ? 8 : 6) : (selected ? 4 : 2.5))
            outline.setFillColor(null)
            outline.setStrokeColor(colour)
            outline.setStroke(edge ? new java.awt.BasicStroke(width, java.awt.BasicStroke.CAP_BUTT,
                java.awt.BasicStroke.JOIN_ROUND, 10f, [6f, 5f] as float[], 0f) : new java.awt.BasicStroke(width))
            overlay.add(outline)
        }
    }
}

/** Keep the full layout inside the viewport instead of compressing or wrapping its rows. */
class CalciumScrollPanel extends javax.swing.JPanel implements javax.swing.Scrollable {
    CalciumScrollPanel() { setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.Y_AXIS)) }
    java.awt.Dimension getPreferredScrollableViewportSize() { getPreferredSize() }
    boolean getScrollableTracksViewportWidth() {
        getParent() instanceof javax.swing.JViewport && getParent().getWidth() >= getPreferredSize().width
    }
    boolean getScrollableTracksViewportHeight() { false }
    int getScrollableUnitIncrement(java.awt.Rectangle visible, int orientation, int direction) { 24 }
    int getScrollableBlockIncrement(java.awt.Rectangle visible, int orientation, int direction) {
        Math.max(24, (int) (orientation == javax.swing.SwingConstants.VERTICAL ? visible.height : visible.width) - 24)
    }
}

/** Fit the initial window on screen; later updates preserve the user's chosen size. */
@CompileStatic
class CalciumWindow extends javax.swing.JFrame {
    CalciumWindow(String title) { super(title) }
    @Override void pack() {
        super.pack()
        def configuration = getGraphicsConfiguration()
        def bounds = configuration.getBounds()
        def insets = java.awt.Toolkit.getDefaultToolkit().getScreenInsets(configuration)
        int screenTop = (int) bounds.getY(), screenHeight = (int) bounds.getHeight()
        int screenWidth = (int) bounds.getWidth()
        int available = screenHeight - insets.top - insets.bottom - 60
        int availableWidth = screenWidth - insets.left - insets.right - 24
        setMinimumSize(new java.awt.Dimension(Math.min(320, availableWidth), Math.min(240, available)))
        setSize(Math.min(getWidth(), availableWidth), Math.min(getHeight(), available))
        if (isShowing()) setLocation(getX(), (int) Math.max(screenTop + insets.top,
            Math.min(getY(), screenTop + screenHeight - insets.bottom - getHeight())))
    }
    void refreshContents() {
        getContentPane().revalidate()
        validate()
        repaint()
    }
}

/** Small plot of one outline's gray value over time (background removed). */
class CalciumTracePanel extends javax.swing.JPanel {
    double[] vals = null; String caption = "Move the mouse over a cell or a suggestion to see its trace"; int band0 = 0, band1 = 0
    CalciumTracePanel() { setPreferredSize(new java.awt.Dimension(330, 118)); setBackground(Color.white); setBorder(javax.swing.BorderFactory.createLineBorder(CalciumTheme.LINE)) }
    void show(double[] v, String cap, int b0, int b1) { vals = v; caption = cap; band0 = b0; band1 = b1; repaint() }
    protected void paintComponent(java.awt.Graphics g0) {
        super.paintComponent(g0)
        def g = (java.awt.Graphics2D) g0.create()
        try {
            CalciumTheme.smooth(g)
            int W = getWidth(), H = getHeight(), L = 38, T = 22, B = 20, R = 10
            g.setFont(CalciumTheme.font(java.awt.Font.PLAIN, -2f)); def fm = g.getFontMetrics()
            g.setColor(CalciumTheme.MUTED); g.drawString(CalciumTheme.fit(caption, fm, W - 12), 8, 14)
            if (vals == null || vals.length < 2) return
            double lo = vals.min(), hi = vals.max(); if (hi - lo < 1) { hi = lo + 1 }
            double pad = (hi - lo) * 0.06d; lo -= pad; hi += pad
            def X = { int f -> L + (int) ((W - L - R) * f / (double) (vals.length - 1)) }
            def Y = { double v -> T + (int) ((H - T - B) * (1 - (v - lo) / (hi - lo))) }
            if (band1 > band0) { g.setColor(new Color(0xFBE6FB)); g.fillRect(X(band0 - 1), T, Math.max(2, X(band1 - 1) - X(band0 - 1)), H - T - B) }
            def ticks = CalciumTheme.ticks(lo, hi, 3); double step = ticks.size() > 1 ? ticks[1] - ticks[0] : 1d
            ticks.each { double v -> int y = Y(v); g.setColor(v == 0d ? new Color(0xAEB8C2) : CalciumTheme.LINE); g.drawLine(L, y, W - R, y)
                String s = CalciumTheme.tick(v, step); g.setColor(CalciumTheme.MUTED); g.drawString(s, L - 5 - fm.stringWidth(s), y + (int) (fm.getAscent() / 2) - 1) }
            g.setColor(new Color(0xC5CED8)); g.drawRect(L, T, W - L - R, H - T - B)
            g.setColor(CalciumTheme.BLUE); g.setStroke(new java.awt.BasicStroke(1.4f, java.awt.BasicStroke.CAP_ROUND, java.awt.BasicStroke.JOIN_ROUND))
            for (int f = 1; f < vals.length; f++) g.drawLine(X(f - 1), Y(vals[f - 1]), X(f), Y(vals[f]))
            g.setColor(CalciumTheme.MUTED); g.drawString("frame 1", L, H - 6)
            String last = "" + vals.length; g.drawString(last, W - R - fm.stringWidth(last), H - 6)
            if (band1 > band0) {
                String s = "response frames"; int x = Math.max(L + fm.stringWidth("frame 1") + 8, Math.min(W - R - fm.stringWidth(last) - fm.stringWidth(s) - 8, X(band0 - 1)))
                g.setColor(new Color(0x9C4F9C)); g.drawString(s, x, H - 6)
            }
        } finally { g.dispose() }
    }
}

/** Update controls are available before a recording is opened. All UI callbacks stay on the event thread. */
class CalciumUpdatesPanel extends CalciumCard {
    def updateBtn = CalciumTheme.button('Check for updates')
    def installBtn = CalciumTheme.button('Install update', true)
    def laterBtn = CalciumTheme.button('Later')
    def releaseBtn = CalciumTheme.button('Release notes')
    def autoUpdateBox = new javax.swing.JCheckBox('Check for updates when this plugin opens', ij.Prefs.get('calcium.steps.check_updates', true))
    def updateNote = CalciumTheme.label('Check for an update before opening a recording.', CalciumTheme.MUTED)
    def installBar = new javax.swing.JProgressBar(0, 100)
    def updateCheck = new CalciumUpdateCheck()
    def installBusy = new java.util.concurrent.atomic.AtomicBoolean(false)
    Map selectedRelease
    String installedVersion
    String releasePage = CalciumRelease.REPOSITORY + '/releases/latest'
    boolean closed = false
    Closure onBusy = { boolean busy -> }
    Closure installRelease = { Map found, Closure progress ->
        String directory = ij.IJ.getDirectory('imagej')
        new CalciumInstaller(applicationDirectory: directory ? new File(directory) : null).install(found, progress)
    }
    Closure browse = { String page -> java.awt.Desktop.getDesktop().browse(new java.net.URI(page)) }
    CalciumUpdatesPanel() {
        super('', 'Plugin updates', false)
        def buttons = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 6, 0))
        buttons.setOpaque(false); buttons.setAlignmentX(0f)
        [updateBtn, installBtn, laterBtn, releaseBtn].each { buttons.add(it) }
        [installBtn, laterBtn].each { it.setVisible(false) }
        CalciumTheme.plain(autoUpdateBox)
        updateCheck.enabled = autoUpdateBox.isSelected()
        installBar.setStringPainted(true); installBar.setVisible(false); installBar.setAlignmentX(0f)
        installBar.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, 22))
        [buttons, javax.swing.Box.createVerticalStrut(6), autoUpdateBox, updateNote, installBar].each { body.add(it) }
        updateBtn.addActionListener({ e -> checkUpdates(false) } as java.awt.event.ActionListener)
        installBtn.addActionListener({ e -> startInstall(selectedRelease) } as java.awt.event.ActionListener)
        laterBtn.addActionListener({ e ->
            installBtn.setVisible(false); laterBtn.setVisible(false)
            message('Update postponed. You can continue or check again when ready.')
        } as java.awt.event.ActionListener)
        releaseBtn.addActionListener({ e ->
            try { browse.call(releasePage) }
            catch (Exception problem) { message('Open this page in your browser: '+CalciumInstaller.escape(releasePage)) }
        } as java.awt.event.ActionListener)
        autoUpdateBox.addActionListener({ e ->
            updateCheck.enabled = autoUpdateBox.isSelected()
            ij.Prefs.set('calcium.steps.check_updates', updateCheck.enabled)
            // A disabled automatic request suppresses its callback. Restore the manual control now.
            if (!updateCheck.enabled && !installBusy.get() && installedVersion == null) {
                updateBtn.setEnabled(true); updateBtn.setText('Check for updates')
            }
        } as java.awt.event.ActionListener)
    }
    void message(String text) {
        updateNote.setText("<html><body style='width:300px'>"+text+'</body></html>')
        revalidate(); repaint()
    }
    void checkUpdates(boolean automatic) {
        if (closed || installBusy.get() || installedVersion != null) return
        boolean accepted = updateCheck.check(automatic, { Map found ->
            if (closed || installBusy.get() || installedVersion != null) return
            updateBtn.setEnabled(true); updateBtn.setText('Check for updates')
            if (automatic && found.error) return // Offline startup never interrupts choosing a recording.
            if (found.error) {
                message('Could not check for updates: '+CalciumInstaller.escape(found.error)+'. You can retry or open Release notes.')
                return
            }
            releasePage = found.page; selectedRelease = found
            installBtn.setVisible(found.newer && found.installable); laterBtn.setVisible(found.newer)
            boolean aheadOfRelease = false
            try { aheadOfRelease = !found.newer && CalciumRelease.newer(CalciumRelease.VERSION, String.valueOf(found.version)) } catch (Exception ignored) { }
            message(found.newer ? '<b>Update available: '+CalciumInstaller.escape(found.version)+'.</b> '+
                (found.installable ? 'Install it here, then restart Fiji to use it.' : 'Open Release notes for installation instructions.') :
                aheadOfRelease ? 'This copy ('+CalciumRelease.VERSION+') is newer than the latest release ('+CalciumInstaller.escape(found.version)+'): a test build that has not been released.' :
                'You are up to date (installed '+CalciumRelease.VERSION+'; latest release '+CalciumInstaller.escape(found.version)+').')
        })
        if (accepted) {
            updateBtn.setEnabled(false); updateBtn.setText('Checking...')
            if (!automatic) message('Checking the latest GitHub release...')
        }
    }
    void startInstall(Map found) {
        if (closed || installedVersion != null || found == null || !found.installable || !installBusy.compareAndSet(false, true)) return
        installBtn.setEnabled(false); updateBtn.setEnabled(false); laterBtn.setEnabled(false)
        autoUpdateBox.setEnabled(false); onBusy.call(true)
        installBar.setValue(0); installBar.setString('Preparing installation...'); installBar.setVisible(true)
        message('Downloading and verifying the update. Please keep this window open until installation finishes.')
        Thread.start('Calcium plugin installation') {
            Map installed = null; String problem = null
            try {
                installed = installRelease.call(found, { int percent, String text ->
                    javax.swing.SwingUtilities.invokeLater { if (!closed) { installBar.setValue(percent); installBar.setString(text) } }
                })
                if (installed == null) throw new IOException('The installer did not confirm completion')
            } catch (Exception failure) { problem = failure.message ?: failure.toString() }
            final Map completed = installed; final String error = problem
            javax.swing.SwingUtilities.invokeLater {
                installBusy.set(false)
                if (completed != null) { installedVersion = completed.version; updateCheck.close() }
                if (closed) return
                onBusy.call(false); autoUpdateBox.setEnabled(true); laterBtn.setEnabled(true)
                if (error != null) {
                    installBar.setString('Update was not installed')
                    installBtn.setEnabled(true); updateBtn.setEnabled(true)
                    message('Could not install: '+CalciumInstaller.escape(error)+'<br>You can retry Install update or open Release notes.')
                } else {
                    installBar.setValue(100); installBar.setString('Installed '+completed.version)
                    installBtn.setText('Installed — restart Fiji'); installBtn.setEnabled(false); laterBtn.setVisible(false)
                    message('<b>Installed '+CalciumInstaller.escape(completed.version)+'.</b> Save your work, then restart Fiji to use it. This session continues with the version already loaded.'+
                        (completed.backup ? '<br>The previous file was backed up in calcium-plugin-backups inside Fiji’s folder.' : ''))
                    updateNote.setToolTipText('Installed file: '+completed.target+(completed.backup ? '; backup: '+completed.backup : ''))
                }
            }
        }
    }
    void close() { closed = true; updateCheck.close() }
}

/** Step 1: which recording, whether to look for a response, and where the results go. */
class CalciumStartDialog {
    static Closure updatesFactory = { -> new CalciumUpdatesPanel() }
    static void fit(javax.swing.JDialog dialog) {
        dialog.pack()
        def config = dialog.getGraphicsConfiguration(), bounds = config.getBounds()
        def insets = java.awt.Toolkit.getDefaultToolkit().getScreenInsets(config)
        dialog.setSize(Math.min(dialog.getWidth(), (int) bounds.getWidth() - insets.left - insets.right - 40),
            Math.min(dialog.getHeight(), (int) bounds.getHeight() - insets.top - insets.bottom - 60))
    }
    /** Blind selection: which frames come before the response. Offers the plugin's window (CalciumSteps.proposedBlindWindow)
     *  with its reason, for the user to confirm or change. Returns [first, last], or null when cancelled. */
    static int[] askBlindFrames(def steps) {
        int n = (int) steps.n; int[] offer = (int[]) steps.proposedBlindWindow(); String why = "" + steps.blindWindowWhy
        boolean rise = steps.found && !steps.riseIsMovement()
        int first = steps.blind ? (int) steps.blindFirst : offer[0], last = steps.blind ? (int) steps.blindLast : offer[1]
        while (true) {
            def gd = new GenericDialog("Frames before the response")
            gd.addMessage(why.replace(". ", ".\n"))
            gd.addNumericField("First frame", first, 0)
            gd.addNumericField("Last frame before the response", last, 0)
            gd.addMessage("Cells are found and chosen on these frames only; measurements use every frame.\nA few cells may respond briefly inside them: end the frames earlier to leave those out.\nStep 4 asks separately when the stimulus was added.")
            gd.showDialog(); if (gd.wasCanceled()) return null
            double a = gd.getNextNumber(), b = gd.getNextNumber()
            if (!Double.isFinite(a) || !Double.isFinite(b) || a != Math.floor(a) || b != Math.floor(b) || a < 1 || b < a || b > n) { IJ.error("Calcium imaging", "Choose whole frame numbers from 1 to " + n + ", the first no later than the last."); continue }
            first = (int) a; last = (int) b
            if (rise && last >= (int) steps.onset && !IJ.showMessageWithCancel("Calcium imaging", "Frame " + last + " is at or after the rise near frame " + steps.onset + ".\nFrames in which cells respond are not blind: cells that respond are brighter in them.\n\nUse these frames anyway?")) continue
            if (offer[1] - last >= 20 && IJ.showMessageWithCancel("Calcium imaging", "The plugin offers frames up to " + offer[1] + ": " + (rise ? "the cells start to rise near frame " + steps.onset : "no response was found") + ".\nYou chose frames " + first + "-" + last + ". With more frames the picture is clearer and more cells are found.\n\nOK: use frames " + first + "-" + offer[1] + ".    Cancel: keep " + first + "-" + last + ".")) last = offer[1]
            return [first, last] as int[]
        }
    }
    /** The file and folder choosers; tests replace them so no native dialog opens. Each returns a path, or null when cancelled. */
    static Closure chooseFile = { -> def od = new ij.io.OpenDialog("Choose the recording (.zvi, .czi, .tif ...)"); od.getFileName() == null ? null : od.getPath() }
    static Closure chooseFolder = { String from -> ij.io.DirectoryChooser.setDefaultDirectory(from); new ij.io.DirectoryChooser("Folder for the results").getDirectory() }
    /** Returns [recording, look, folder, file], or null when the dialog was cancelled. "file" is the path chosen for the
     *  from-disk entry (the last label unless fromDisk names another); it is null for a recording that is already open. */
    static Map ask(List<String> labels, String preselected, String defaultFolder, String fromDisk = null) {
        Map answer = null
        String diskLabel = fromDisk ?: (labels ? labels[-1] : null)
        def show = {
            def dialog = new javax.swing.JDialog((java.awt.Frame) null, "Calcium imaging: cells step by step", true)
            def content = new javax.swing.JPanel(); content.setLayout(new javax.swing.BoxLayout(content, javax.swing.BoxLayout.Y_AXIS))
            content.setBorder(javax.swing.BorderFactory.createEmptyBorder(18, 20, 14, 20)); content.setBackground(CalciumTheme.PAGE)
            def hint = { String text -> CalciumTheme.label(text, CalciumTheme.MUTED, java.awt.Font.PLAIN, -1f) }
            content.add(CalciumTheme.label("Step 1 of 4", CalciumTheme.MUTED, java.awt.Font.PLAIN, -1f))
            content.add(CalciumTheme.label("Choose the recording", CalciumTheme.INK, java.awt.Font.BOLD, 5f)); content.add(javax.swing.Box.createVerticalStrut(12))
            def recording = new CalciumCard("", "Recording", false)
            def choice = new javax.swing.JComboBox<String>(labels as String[]); choice.setSelectedItem(preselected); choice.setFont(CalciumTheme.font()); choice.setAlignmentX(0f)
            choice.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) choice.getPreferredSize().height))
            // How cells are chosen. Blind is the default: for % responding, whether a cell is chosen must not depend on its response.
            def blindChoice = new javax.swing.JRadioButton("Before the response only (blind): for % responding", true)
            def allChoice = new javax.swing.JRadioButton("All frames", false)
            def look = new javax.swing.JRadioButton("All frames, with the response in magenta (not for % responding)", false)
            blindChoice.setToolTipText("<html>Cellpose and you see only the frames before any cell responds, and traces are hidden while choosing,<br>so a cell is chosen the same way whether or not it responds. The plugin offers the frames; you confirm them.</html>")
            allChoice.setToolTipText("<html>The average of all frames. It includes the frames after the stimulus,<br>so cells that stay bright after it are a little easier to see.</html>")
            look.setToolTipText("<html>Cells that brighten after the stimulus turn magenta and are easier to find.<br>Dim cells are then chosen only if they respond, which raises % responding.<br>Use it for the size or timing of responses, not for % responding.</html>")
            def choiceGroup = new javax.swing.ButtonGroup(); [blindChoice, allChoice, look].each { choiceGroup.add(it); CalciumTheme.plain(it) }
            // Choosing "from disk" asks for the file here, straight away, and shows its name; Next is not needed to get the file chooser.
            String file = null
            def fileButton = CalciumTheme.button("Choose file…"), fileName = CalciumTheme.label("No file chosen yet", CalciumTheme.MUTED)
            def fileRow = new javax.swing.JPanel(new java.awt.BorderLayout(8, 0)); fileRow.setOpaque(false); fileRow.setAlignmentX(0f)
            fileRow.add(fileButton, java.awt.BorderLayout.WEST); fileRow.add(fileName, java.awt.BorderLayout.CENTER)
            fileRow.setBorder(javax.swing.BorderFactory.createEmptyBorder(6, 0, 0, 0))
            fileName.setPreferredSize(new java.awt.Dimension(200, (int) fileName.getPreferredSize().height))
            fileRow.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) fileRow.getPreferredSize().height))
            boolean choosing = false
            def pickFile = { ->
                if (choosing) return
                choosing = true
                try {
                    String chosen = chooseFile.call()
                    if (chosen) { file = chosen; fileName.setText(new File(chosen).getName()); fileName.setForeground(CalciumTheme.INK); fileName.setToolTipText(chosen); fileButton.setText("Change file…") }
                } finally { choosing = false }
            }
            def showFileRow = { -> fileRow.setVisible(choice.getSelectedItem() == diskLabel); fit(dialog) }
            fileButton.addActionListener({ e -> pickFile() } as java.awt.event.ActionListener)
            choice.addActionListener({ e -> showFileRow(); if (dialog.isVisible() && choice.getSelectedItem() == diskLabel) javax.swing.SwingUtilities.invokeLater { pickFile() } } as java.awt.event.ActionListener)
            fileRow.setVisible(choice.getSelectedItem() == diskLabel)
            [choice, fileRow, javax.swing.Box.createVerticalStrut(4), hint("Pick the original movie, not a processed copy."), javax.swing.Box.createVerticalStrut(10),
             CalciumTheme.label("Choose cells on:", CalciumTheme.INK), blindChoice, allChoice, look].each { recording.body.add(it) }
            def results = new CalciumCard("", "Folder for the results", false)
            def folder = new javax.swing.JTextField(defaultFolder, 30); folder.setFont(CalciumTheme.font()); def browse = CalciumTheme.button("Browse…")
            def line = new javax.swing.JPanel(new java.awt.BorderLayout(6, 0)); line.setOpaque(false); line.setAlignmentX(0f); line.add(folder, java.awt.BorderLayout.CENTER); line.add(browse, java.awt.BorderLayout.EAST)
            line.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) line.getPreferredSize().height))
            [line, javax.swing.Box.createVerticalStrut(4), hint("A subfolder with the recording's name is made inside this folder.")].each { results.body.add(it) }
            def next = CalciumTheme.button("Next", true), cancel = CalciumTheme.button("Cancel")
            def buttons = new javax.swing.JPanel(new java.awt.BorderLayout()); buttons.setOpaque(false); buttons.setAlignmentX(0f)
            def right = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT, 6, 0)); right.setOpaque(false); right.add(cancel); right.add(next)
            buttons.add(CalciumTheme.label("Version " + CalciumRelease.VERSION, CalciumTheme.FAINT, java.awt.Font.PLAIN, -1f), java.awt.BorderLayout.WEST); buttons.add(right, java.awt.BorderLayout.EAST)
            def updates = updatesFactory.call()
            updates.onBusy = { boolean busy -> next.setEnabled(!busy); cancel.setEnabled(!busy) }
            def dismiss = { -> if (!updates.installBusy.get()) { updates.close(); dialog.dispose() } }
            [recording, javax.swing.Box.createVerticalStrut(10), results, javax.swing.Box.createVerticalStrut(10), updates].each { content.add(it) }
            browse.addActionListener({ e ->
                String chosen = chooseFolder.call(folder.getText())
                if (chosen) folder.setText(chosen)
            } as java.awt.event.ActionListener)
            next.addActionListener({ e ->
                if (updates.installBusy.get()) return
                boolean disk = choice.getSelectedItem() == diskLabel
                if (disk && file == null) pickFile()          // Next without a file still asks for one,
                if (disk && file == null) return              // and stays here if that was cancelled
                answer = [recording: choice.getSelectedItem() as String, look: look.isSelected(), blind: blindChoice.isSelected(), folder: folder.getText(), file: disk ? file : null]; dismiss()
            } as java.awt.event.ActionListener)
            cancel.addActionListener({ e -> dismiss() } as java.awt.event.ActionListener)
            dialog.getRootPane().setDefaultButton(next)
            dialog.getRootPane().registerKeyboardAction({ e -> dismiss() } as java.awt.event.ActionListener, javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ESCAPE, 0), javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
            def scroll = new javax.swing.JScrollPane(content)
            scroll.setBorder(javax.swing.BorderFactory.createEmptyBorder()); scroll.getViewport().setBackground(CalciumTheme.PAGE)
            scroll.getVerticalScrollBar().setUnitIncrement(24); scroll.setWheelScrollingEnabled(true)
            def page = new javax.swing.JPanel(new java.awt.BorderLayout()); page.setBackground(CalciumTheme.PAGE)
            buttons.setBorder(javax.swing.BorderFactory.createEmptyBorder(4, 20, 14, 20))
            page.add(scroll, java.awt.BorderLayout.CENTER); page.add(buttons, java.awt.BorderLayout.SOUTH)
            dialog.setContentPane(page); dialog.setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE)
            dialog.addWindowListener(new java.awt.event.WindowAdapter() {
                void windowClosing(java.awt.event.WindowEvent event) { dismiss() }
                void windowClosed(java.awt.event.WindowEvent event) { updates.close() }
                void windowOpened(java.awt.event.WindowEvent event) { updates.checkUpdates(true) }
            })
            fit(dialog); dialog.setResizable(true); dialog.setLocationRelativeTo(null); dialog.setVisible(true)
        }
        if (javax.swing.SwingUtilities.isEventDispatchThread()) show.call() else javax.swing.SwingUtilities.invokeAndWait(show as Runnable)
        return answer
    }
}

// ---------------------------------------------------------------------------------------------------------------
// The part with windows. (When the variable CALCIUM_TEST is set, only the class above is returned, for testing.)
if (binding.hasVariable("CALCIUM_TEST")) return CalciumSteps
boolean testGui = binding.hasVariable("CALCIUM_TEST_GUI")

// ---- Step 1: which recording, where to save ----
def movies = []
def idList = WindowManager.getIDList()
if (idList != null) idList.each { id -> def im = WindowManager.getImage(id); if (im != null && im.getStackSize() > 1) movies << im }
def fromDisk = "Choose a file from disk..."
def labels = []
movies.eachWithIndex { im, i -> labels << ("" + (i + 1) + ":  " + im.getTitle() + "   (" + im.getStackSize() + " frames, " + im.getBitDepth() + "-bit)") }
labels << fromDisk
def front = WindowManager.getCurrentImage()
def pre = fromDisk
movies.eachWithIndex { im, i -> if (im == front) pre = labels[i] }
if (pre == fromDisk && movies.size() > 0) pre = labels[0]
def defaultDir = Prefs.get("calcium.steps.folder", System.getProperty("user.home") + File.separator + "Documents" + File.separator + "Calcium imaging results")

def start = CalciumStartDialog.ask(labels, pre, defaultDir, fromDisk)
if (start == null) return
def picked = start.recording
boolean look = start.look
boolean blindStart = start.blind == true
def outRoot = start.folder
if (outRoot == null || outRoot.trim().isEmpty()) outRoot = defaultDir
Prefs.set("calcium.steps.folder", outRoot)

ImagePlus imp = null
if (picked == fromDisk) {
    String recordingPath = start.file      // chosen in the Step 1 dialog
    if (recordingPath == null) return
    IJ.run("Bio-Formats Importer", "open=[" + recordingPath + "] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT")
    imp = WindowManager.getCurrentImage()
} else imp = movies[labels.indexOf(picked)]
if (imp == null || imp.getStackSize() < 2) { IJ.error("Calcium imaging", "No movie was opened."); return }

IJ.showStatus("Making the cells visible...")
if (imp.getBitDepth() == 32) IJ.log("Calcium imaging plugin: this movie is 32-bit, so it looks like a processed copy. Background is meant to be taken on the original movie.")
def steps = new CalciumSteps()
steps.prepare(imp, look || blindStart)   // blind: the response is looked for only to offer when the stimulus came; it is never shown while choosing
if (blindStart) {
    int[] frames = CalciumStartDialog.askBlindFrames(steps)
    if (frames == null) { IJ.showStatus("Calcium imaging: cancelled"); return }
    steps.startBlind(frames[0], frames[1])
}
steps.colourPic.show()
IJ.run(steps.colourPic, "Set... ", "zoom=75")
def rm = testGui ? new RoiManager(false) : RoiManager.getRoiManager()
if (!testGui) { rm.runCommand(steps.colourPic, "Show None"); IJ.setTool("oval") } // The plugin draws its own contrasted overlay.
def base = imp.getTitle(); int dot = base.lastIndexOf("."); if (dot > 0) base = base.substring(0, dot)
def outDir = new File(outRoot, base)

// ---- the window that stays open for steps 2 and 3 ----
def GRADE_TEXT = ["responds on its own", "likely a cell", "unsure", "very bright: probably precipitate"]
def GRADE_SHARE = ["about 8 in 10 such outlines were real cells in tests", "about 8 in 10 such outlines were real cells in tests", "about 4 in 10 such outlines were real cells in tests", "about 1 in 10 such outlines were real cells in tests"]
def CYAN = new Color(0, 220, 255), GREY = new Color(205, 205, 205), RED = new Color(255, 70, 70), BLUE = new Color(40, 150, 255)
def DASH = new java.awt.BasicStroke(1f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_MITER, 10f, [4f, 4f] as float[], 0f)
def small = CalciumTheme.font(), bold = CalciumTheme.font(java.awt.Font.BOLD, 1f)
def MUTED = CalciumTheme.hex(CalciumTheme.MUTED)
def lineIcon = { Color c, int width, boolean dashed, boolean filled ->
    def bi = new java.awt.image.BufferedImage(30, 14, java.awt.image.BufferedImage.TYPE_INT_ARGB); def g = bi.createGraphics()
    g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
    g.setColor(new Color(24, 30, 36)); g.fillRoundRect(0, 0, 30, 14, 6, 6); g.setColor(c)
    if (filled) g.fillRoundRect(5, 3, 20, 8, 3, 3)
    else { g.setStroke(dashed ? DASH : new java.awt.BasicStroke((float) width)); g.drawOval(5, 2, 20, 10) }
    g.dispose(); new javax.swing.ImageIcon(bi)
}
// The legend is also the switch for what is drawn on the picture: click an entry to hide or show that kind (your cells,
// the background, or one grade of suggestion); double-click an entry to see only that kind, and again to see everything.
def hiddenKinds = [] as Set          // 'cells', 'background', 'g0' responds, 'g1' likely, 'g2' unsure, 'g3' precipitate?
def legendEntries = [:]
Closure legendClicked = null         // set below, once the picture can be redrawn
def legendRow = { javax.swing.Icon ic, String s, String kind = null ->
    def l = new javax.swing.JLabel(s, ic, javax.swing.SwingConstants.LEFT); l.setFont(CalciumTheme.font(java.awt.Font.PLAIN, -1f)); l.setForeground(CalciumTheme.MUTED); l.setIconTextGap(6); l.setAlignmentX(0f)
    if (kind != null) {
        legendEntries[kind] = [label: l, text: s]
        l.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR))
        l.setToolTipText("Click to hide or show these on the picture. Double-click to see only these (double-click again for everything).")
        l.addMouseListener(new java.awt.event.MouseAdapter() { void mouseClicked(java.awt.event.MouseEvent e) { if (legendClicked != null) legendClicked.call(kind, e.getClickCount() >= 2) } })
    }
    l
}
def section = { String number, String title, boolean collapsible = true -> new CalciumCard(number, title, collapsible) }
def row = { java.awt.Component... cs -> def p = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 6, 3)); p.setOpaque(false); cs.each { p.add(it) }; p.setAlignmentX(0f); p }
def note = { String s -> def l = new javax.swing.JLabel("<html><body style='width:300px'>" + s + "</body></html>"); l.setFont(small); l.setForeground(CalciumTheme.INK); l.setAlignmentX(0f); l }
def hint = { String s -> def l = note(s); l.setForeground(CalciumTheme.MUTED); l }
def btn = { String s -> CalciumTheme.button(s) }
def box = { String s, boolean on -> def b = new javax.swing.JCheckBox(s, on); b.setFont(small); b.setOpaque(false); b.setAlignmentX(0f); b }

def frame = new CalciumWindow("Calcium imaging: cells step by step")
frame.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
def root = new CalciumScrollPanel(); root.setBorder(javax.swing.BorderFactory.createEmptyBorder(10, 10, 10, 10)); root.setBackground(CalciumTheme.PAGE); root.setOpaque(true)

// One card per step. Step 1 is already done; of steps 2 to 4 only one is open at a time, so the window fits the screen.
def s1 = section("1", "Recording", false); s1.setDone(true)
def recordingNote = note("&nbsp;")
// "The field moved (about 16 px) during frames 1–10.", or nothing (also for a window built without the movement check)
// For Step 4 the movement is counted from the frames the cells were chosen on whenever such frames were confirmed, also after cells were added with the response.
def motionForPage = { chosen -> try { return chosen != null && steps.motionShifts != null ? CalciumMotion.summary(steps.motionShifts, chosen.first as int, chosen.last as int) : steps.motion } catch (Throwable ignored) { return null } }
def movedNote = { -> try { steps.sayMotion(); return CalciumMotion.sentence(steps.motion) } catch (Throwable ignored) { return "" } }
def describeRecording = { ->
    recordingNote.setText("<html><body style='width:300px'><b>" + CalciumTheme.middle(imp.getTitle(), 40) + "</b><br><font color='" + MUTED + "'>" + steps.n + " frames · " + (steps.blind ? "cells chosen blind on frames " + steps.blindFirst + "–" + steps.blindLast : steps.found ? "response near frame " + steps.onset + " (magenta)" : (steps.looked ? "no response found (no magenta)" : "response not looked for")) + "</font>" + (movedNote() ? "<br><font color='#9A5B00'>" + movedNote() + "</font>" : "") + "</body></html>")
    recordingNote.setToolTipText(imp.getTitle() + ". " + (steps.blind ? "Cells are chosen on frames " + steps.blindFirst + "–" + steps.blindLast + ", before the response" + (steps.found ? "; the cells start to rise near frame " + steps.onset : "") + ". The response is not shown while choosing." : steps.found ? "Brightness rises near frame " + steps.onset + "; the response can be shown in magenta." : (steps.looked ? "No response was found, so there is no magenta. Every cell is still shown in green." : "The response has not been looked for. Choose Response (magenta) in Step 2 to look for it.")))
}
describeRecording()
def changeBtn = btn("‹ Change recording")
changeBtn.setToolTipText("Back to the first screen, to choose another recording or results folder")
s1.body.add(recordingNote); s1.body.add(row(changeBtn))

def s2 = section("2", "Choose the cells")
String selectionBasis = steps.blind ? steps.blindBasis() : steps.found ? 'Response-assisted review' : 'All-frames structural review'
// What the colour picture shows while cells are chosen. The choice is recorded with the measurement as the selection basis.
def viewAllBtn = btn('All frames'), viewBaselineBtn = btn('Before response (blind)…'), viewResponseBtn = btn('Response (magenta)')
viewAllBtn.setToolTipText('The average of all frames with the background removed: cells by their structure')
viewBaselineBtn.setToolTipText('<html>Only frames before the response. Suggestions are made on this picture too, traces are hidden<br>and grades use shape, size and brightness only, so whether a cell is chosen cannot depend on its response.</html>')
def blindNote = { -> 'Blind: frames ' + steps.blindFirst + '–' + steps.blindLast + ' only, before the response. Traces are hidden; grades use shape, size and brightness.' }
viewResponseBtn.setToolTipText('<html>Cells that brighten after the stimulus turn magenta.<br>Choosing cells with this on is recorded as response-assisted selection,<br>which can raise percent responding.</html>')
def selectionNote = hint('&nbsp;')
def showView = { String which, String text ->
    [all: viewAllBtn, baseline: viewBaselineBtn, response: viewResponseBtn].each { key, button -> CalciumTheme.choose(button, key == which) }
    selectionNote.setText("<html><body style='width:300px'>" + text + "</body></html>")
}
def suggestBtn = CalciumTheme.button("Suggest cells automatically", true)
suggestBtn.setToolTipText("Cellpose outlines the cells; then click inside each outline that is a real cell")
def acceptAllBtn = btn("Accept likely"); acceptAllBtn.setEnabled(false)
acceptAllBtn.setToolTipText("Accept every suggestion graded 'responds' or 'likely a cell'")
def ovalBtn = btn("Oval"), freeBtn = btn("Freehand"), addBtn = btn("Add outline  (t)")
addBtn.setToolTipText("Draw around a cell on the picture, then add it here or press t")
def removeBtn = btn("Remove selected"), undoBtn = btn("Undo"); undoBtn.setEnabled(false)
removeBtn.setToolTipText("Click a white outline to select it, then remove it"); undoBtn.setToolTipText("Undo the last change to your cells")
def removeMode = box("Click cells to remove", false)
def editHint = hint("")   // messages about removing and undoing appear here; the how-to is on the buttons' tooltips
removeMode.setToolTipText("While this is on, clicking a cell on the picture removes it")
def showSug = box("Show suggestions", true), showBg = box("Show background", true)
showBg.setToolTipText("The pixels used as each cell's background, in blue")
def edgeBox = box("Leave out cells within " + CalciumSteps.EDGE_MARGIN + " px of the edge", true)
edgeBox.setToolTipText("<html>At the edge a cell's background can only be taken from the inner, better-lit side,<br>so the cell comes out too dim. Untick to suggest, accept and measure edge cells like the others.</html>")
def status = new javax.swing.JLabel("<html><body style='width:300px'><b>Cells chosen: 0</b></body></html>"); status.setFont(small); status.setForeground(CalciumTheme.INK); status.setAlignmentX(0f)
def info = note("&nbsp;")
def trace = new CalciumTracePanel(); trace.setAlignmentX(0f); trace.setFont(small)
if (steps.blind) { trace.show(null, "Traces are hidden while choosing cells blind", 0, 0); trace.setVisible(false) }   // no empty frame while choosing blind
// Step 2 leads through its own order: suggest, accept the likely ones, check the dashed ones, add the cells that only show when
// they respond, continue. This line always names the one thing to do next (set in refresh, from what has been done so far).
def nextLine = note("&nbsp;")
def sayNext = { String text -> nextLine.setText("<html><body style='width:300px'><b><font color='#1F6FB2'>Next:</font></b> " + text + "</body></html>") }
def addRespondersBtn = btn("Add cells that respond…")
addRespondersBtn.setToolTipText("<html>Shows the response (magenta) and looks for the cells that only show when they respond.<br>Your cells so far stay recorded as chosen blind; the ones you add are counted apart in Step 4.</html>")
addRespondersBtn.setVisible(steps.blind == true)
boolean detecting = false, autoSuggestAsked = false
s2.body.add(nextLine); s2.body.add(javax.swing.Box.createVerticalStrut(4))
s2.body.add(row(CalciumTheme.label('Picture:', CalciumTheme.MUTED), viewAllBtn, viewBaselineBtn, viewResponseBtn)); s2.body.add(selectionNote)
def cpLabel = hint("&nbsp;")
def setupBtn = btn("Set up Cellpose (one time)"); def setupRow = row(setupBtn)
def detectionProgress = new CalciumOutlineProgress()
def continueBtn = btn("Continue to measure  ›")
showView(steps.blind ? 'baseline' : steps.found ? 'response' : 'all', steps.blind ? blindNote() : steps.found ? 'Magenta: cells that brightened near frame ' + steps.onset + '. Picking cells this way is recorded as response-assisted.' : 'Choosing cells on the average of all frames.')
s2.body.add(row(suggestBtn, acceptAllBtn)); s2.body.add(row(addRespondersBtn, continueBtn)); s2.body.add(detectionProgress); s2.body.add(cpLabel); s2.body.add(setupRow)
// Drawing by hand and what is shown are used now and then: one line opens them, so the sequence above stays short.
def moreRows = new javax.swing.JPanel(); moreRows.setLayout(new javax.swing.BoxLayout(moreRows, javax.swing.BoxLayout.Y_AXIS)); moreRows.setOpaque(false); moreRows.setAlignmentX(0f)
moreRows.add(row(CalciumTheme.label("Draw:", CalciumTheme.MUTED), ovalBtn, freeBtn, addBtn)); moreRows.add(row(showSug, showBg)); moreRows.add(row(edgeBox))
def moreBtn = btn("Draw by hand, what is shown  ▸"); moreBtn.setToolTipText("Oval, Freehand and Add outline; show or hide suggestions and background; cells at the edge")
moreRows.setVisible(false)
moreBtn.addActionListener { moreRows.setVisible(!moreRows.isVisible()); moreBtn.setText("Draw by hand, what is shown  " + (moreRows.isVisible() ? "▾" : "▸")); frame.revalidate(); frame.repaint() }
s2.body.add(javax.swing.Box.createVerticalStrut(4)); s2.body.add(status)
s2.body.add(row(removeBtn, undoBtn, removeMode)); s2.body.add(editHint)
s2.body.add(row(moreBtn)); s2.body.add(moreRows)
s2.body.add(info); s2.body.add(trace); s2.body.add(javax.swing.Box.createVerticalStrut(6))
def legend = new javax.swing.JPanel(new java.awt.GridLayout(0, 3, 6, 2)); legend.setAlignmentX(0f); legend.setOpaque(false)
def acceptedIcon = new java.awt.image.BufferedImage(30, 14, java.awt.image.BufferedImage.TYPE_INT_ARGB)
def acceptedGraphics = acceptedIcon.createGraphics()
acceptedGraphics.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
acceptedGraphics.setColor(new Color(24, 30, 36)); acceptedGraphics.fillRoundRect(0, 0, 30, 14, 6, 6)
acceptedGraphics.setColor(Color.black); acceptedGraphics.setStroke(new java.awt.BasicStroke(5f)); acceptedGraphics.drawOval(6, 3, 18, 8)
acceptedGraphics.setColor(Color.white); acceptedGraphics.setStroke(new java.awt.BasicStroke(2f)); acceptedGraphics.drawOval(6, 3, 18, 8); acceptedGraphics.dispose()
legend.add(legendRow(new javax.swing.ImageIcon(acceptedIcon), "your cells", "cells"))
def respondsLegend = legendRow(lineIcon(CYAN, 2, false, false), "responds", "g0"); respondsLegend.setVisible(!steps.blind)   // no such grade while choosing blind
legend.add(legendRow(lineIcon(CYAN, 1, false, false), "likely a cell", "g1"))
legend.add(legendRow(lineIcon(GREY, 1, true, false), "unsure", "g2"))
legend.add(legendRow(lineIcon(BLUE, 1, false, true), "background", "background"))
legend.add(legendRow(lineIcon(RED, 2, false, false), "precipitate?", "g3"))
legend.add(respondsLegend)   // last, so the grid has no hole when it is hidden
legend.setToolTipText("White: your cells. Blue: the background used. The others are suggestions, by grade. Click an entry to hide or show it on the picture; double-click to see only that one.")
def legendHint = CalciumTheme.label("Click an entry to hide or show it on the picture; double-click to see only that one.", CalciumTheme.MUTED)
legendHint.setFont(CalciumTheme.font(java.awt.Font.PLAIN, -2f)); legendHint.setAlignmentX(0f)
legend.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) legend.getPreferredSize().height))
s2.body.add(legend)   // how to use the legend is in its tooltip

def s3 = section("3", "Measure and save")
def makeMovie = box("Also open a background-corrected movie", true)
def saveMovie = box("Save that movie as a TIFF too (large)", false)
saveMovie.setToolTipText("About 1.7 GB for 300 frames of 1388 x 1040")
def sensBox = box("Sensitivity check after measuring (blind only)", true); sensBox.setVisible(steps.blind == true)
sensBox.setToolTipText("<html>Runs response-assisted detection once, for comparison, and reports how many of its outlines<br>the blind detection also found. A second Cellpose run: about 20 s with a graphics chip, a few minutes<br>without; saved outlines are reused. Its outlines are never added to your cells or measured.</html>")
def go = CalciumTheme.button("Measure and save", true)
def savedIn = "..." + File.separator + CalciumTheme.middle(outDir.getParentFile().getName(), 24) + File.separator + "<br>" + CalciumTheme.middle(outDir.getName(), 44)
def result = note("Saves the tables in: " + savedIn); result.setToolTipText(outDir.getPath())
def openBtn = btn("Open results folder"); openBtn.setVisible(false)
def closeBtn = btn("Close")
closeBtn.setToolTipText("Close this window; the pictures and saved results stay")
def backToCellsBtn = btn('‹ Back to cells')
s3.body.add(makeMovie); s3.body.add(saveMovie); s3.body.add(sensBox); s3.body.add(row(go, openBtn)); s3.body.add(result); s3.body.add(row(backToCellsBtn))
def s4 = section("4", "Analyze responses")
def analyzeBtn = CalciumTheme.button("Open analysis page", true); analyzeBtn.setEnabled(false)
analyzeBtn.setToolTipText("Opens the analysis page in your browser. It runs on this computer only; nothing is sent anywhere.")
def classicBtn = btn("Use a Fiji window instead"); classicBtn.setEnabled(false)
classicBtn.setToolTipText("The same analysis in a Fiji window, for when a browser is not wanted")
analyzeBtn.addPropertyChangeListener("enabled", { classicBtn.setEnabled(analyzeBtn.isEnabled()) } as java.beans.PropertyChangeListener)
def analysisNote = hint("Measure and save first, then open the analysis page to review the recording and run analysis. Keep Fiji open. One-time engine setup needs internet.")
Map analysisMeasurement = null
def backToMeasureBtn = btn('‹ Back to measuring')
s4.body.add(row(analyzeBtn)); s4.body.add(analysisNote); s4.body.add(row(classicBtn)); s4.body.add(row(backToMeasureBtn))
def stepCards = [s2, s3, s4]
def openStep = { CalciumCard card -> stepCards.each { it.setExpanded(it.is(card)) }; frame.refreshContents() }
stepCards.each { it.onToggle = { CalciumCard card -> if (card.expanded) { card.setExpanded(false); frame.refreshContents() } else openStep(card) } }
s3.setExpanded(false); s4.setExpanded(false); s3.setSummary("after choosing cells"); s4.setSummary("SimpleCa²⁺, after measuring")
continueBtn.addActionListener({ e -> openStep(s3) } as java.awt.event.ActionListener)
backToCellsBtn.addActionListener({ e -> openStep(s2) } as java.awt.event.ActionListener)
backToMeasureBtn.addActionListener({ e -> openStep(s3) } as java.awt.event.ActionListener)
[s1, s2, s3, s4].each { root.add(it); root.add(javax.swing.Box.createVerticalStrut(8)) }
root.add(row(closeBtn, CalciumTheme.label("Version " + CalciumRelease.VERSION, CalciumTheme.FAINT, java.awt.Font.PLAIN, -1f)))
def scroll = new javax.swing.JScrollPane(root, javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED, javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED)
scroll.setBorder(javax.swing.BorderFactory.createEmptyBorder()); scroll.getViewport().setBackground(CalciumTheme.PAGE)
scroll.setWheelScrollingEnabled(true)
scroll.getVerticalScrollBar().setUnitIncrement(24)
scroll.getHorizontalScrollBar().setUnitIncrement(24)
frame.getContentPane().add(scroll)
frame.pack()
// An ordinary window: clicking a picture or another program puts this window behind it, like any other window.
// (It used to stay above all of Fiji's windows, and could not be sent back.) It opens beside the picture, so both stay visible.
if (!testGui) { def scr = java.awt.Toolkit.getDefaultToolkit().getScreenSize(); frame.setLocation((int) Math.max(0, scr.width - frame.getWidth() - 12), 40) }
if (!testGui) CalciumTheme.beside(steps.colourPic, frame)   // the colour picture takes the free space on the left; this window stays on the right

// ---- what the window does ----
def lastKey = ""; int[][] lastBg = new int[0][]
def proposals = []        // automatic outlines: [roi, grade (0 responds, 1 likely, 2 unsure, 3 bright), ev (the measured evidence), cx, cy]
def cellEdits = new CalciumCellEdits(manager: rm, picture: steps.colourPic)
def measuredKey = { Roi[] rois -> rois.collect { CalciumInstaller.digest(RoiEncoder.saveAsByteArray(it)) }.join(':') }
// Blind selection stays blind only while nothing of the response has been shown: a trace, the magenta picture or Step 4.
// The outlines held at that moment were chosen blind; outlines changed afterwards are recorded as not blind.
boolean responseSeen = steps.found && !steps.blind && look      // magenta shown from the start
String blindKey = null
// Per cell: the outlines held when the response was first shown. A cell with one of these outlines, unchanged, was chosen blind;
// a cell added or reshaped afterwards was not. Lets blind and response-assisted choosing be used together and counted apart.
Set blindCellKeys = null             // set when blind choosing ends; null while it has never ended
Set notBlindKeys = [] as Set         // cells added outside blind choosing, kept apart if blind choosing is started again
def cellKey = { Roi r -> r.getBounds().toString() + '#' + r.getContainedPoints().length }
Closure endBlindChoosing = null      // set below, once the outlines can be read
Closure reveal = null                                             // set below, once the outlines can be read
def basisFor = { Roi[] all ->
    if (!steps.blind) {
        // cells chosen blind first and others added afterwards: say so, with the counts (the registry has it cell by cell)
        int blindly = steps.blindLast > 0 && blindCellKeys != null ? (all.count { blindCellKeys.contains(cellKey(it)) } as int) : 0
        return blindly == 0 ? selectionBasis : 'Blind on frames ' + steps.blindFirst + '-' + steps.blindLast + ', then ' + (selectionBasis.startsWith('Response') ? 'response-assisted' : 'all frames') + ': ' + blindly + ' chosen blind, ' + (all.length - blindly) + ' added afterwards'
    }
    (!responseSeen || measuredKey(all) == blindKey) ? steps.blindBasis() :
        'Picture of frames ' + steps.blindFirst + '-' + steps.blindLast + ' (before the response), not blind: the response had been shown before these outlines were final'
}
// Step 4 opens as a page in the browser, served from this computer by the plugin. The Fiji window does the same analysis
// and is used when asked for, or when the page cannot be served or opened here.
CalciumWebPages stepFourPages = new CalciumWebPages()
def openFijiWindow = { -> def window = new CalciumAnalysisWindow(analysisMeasurement); window.besideControls(frame); window.setVisible(true) }
analyzeBtn.addActionListener({ e ->
    if (analysisMeasurement == null) return
    reveal()
    if (testGui) { openFijiWindow(); return }
    try {
        def stepFourPage = stepFourPages.page(analysisMeasurement)
        CalciumWebServer.openBrowser(frame,new java.net.URI(stepFourPage.url()))
        s4.setSummary("open in your browser")
    } catch (Throwable problem) {
        IJ.log("Calcium imaging plugin: the analysis page could not be opened in a browser (" + problem + "). Using the Fiji window instead.")
        openFijiWindow()
    }
} as java.awt.event.ActionListener)
classicBtn.addActionListener({ e -> if (analysisMeasurement != null) { reveal(); openFijiWindow() } } as java.awt.event.ActionListener)
// Browser pages keep their own measurement snapshots until Close analysis or Fiji exits.
def cellRois = { ->
    Roi[] accepted=rm.getRoisAsArray()
    Roi[] unique=accepted.collect{it} as Roi[]
    CalciumMeasurement.ensureIds(unique)
    for(int i=0;i<accepted.length;i++)if(!accepted[i].is(unique[i]))rm.setRoi(unique[i],i)
    def all = unique as List
    all.findAll { it != null && it.isArea() }.collect { def r = (Roi) it.clone(); r.setPosition(0); r } as Roi[]
}
endBlindChoosing = { -> if (steps.blindLast > 0 && blindCellKeys == null) blindCellKeys = (cellRois().collect { cellKey(it) } as Set) - notBlindKeys }
reveal = { -> if (!responseSeen) { responseSeen = true; blindKey = measuredKey(cellRois()); endBlindChoosing() } }
// the outlines in the ROI Manager, split into the cells that count and (when the box is ticked) the ones at the picture's edge
def splitCells = { ->
    Roi[] all = cellRois(); boolean cut = edgeBox.isSelected()
    def keptIdx = (0..<all.length).findAll { !(cut && steps.nearEdge(all[it])) }
    [all: all, keptIdx: keptIdx, cells: keptIdx.collect { all[it] } as Roi[], edge: (0..<all.length).findAll { !keptIdx.contains(it) }.collect { all[it] } as Roi[]]
}
def unusedProposals = { Roi[] cells ->
    proposals.findAll { p -> !(edgeBox.isSelected() && p.edge) && !cells.any { c -> c.contains(p.cx, p.cy) || p.roi.contains((int) c.getBounds().getCenterX(), (int) c.getBounds().getCenterY()) } }
}
def evidenceText = { p ->
    double[] e = p.ev; def parts = []
    parts << String.format(Locale.US, "%.0f um across", 2 * Math.sqrt(e[1] / Math.PI) * steps.movie.getCalibration().pixelWidth)
    if (steps.respDiff != null && !p.blind && !steps.blind) parts << String.format(Locale.US, "rise after the stimulus %.1f gray levels inside, %.1f around it", e[5], e[6])
    parts << String.format(Locale.US, "%.1f gray levels brighter than its surroundings", e[4])
    if (e[3] >= 0.02) parts << String.format(Locale.US, "%.0f %% of it is very bright", 100 * e[3])
    if (e[7] > 0) parts << "touches the edge of the picture"
    if (e[1] > 2300) parts << "large: check it is one cell, not two"
    "<b>Suggestion: " + GRADE_TEXT[p.grade] + "</b> (" + GRADE_SHARE[p.grade] + ")<br>" + parts.join("; ") + "."
}
def refresh = { boolean force ->
    try {
        def sp = splitCells(); Roi[] all = (Roi[]) sp.all; Roi[] cells = (Roi[]) sp.cells; Roi[] edgeCells = (Roi[]) sp.edge
        def key = all.collect { it.getBounds().toString() + it.getContainedPoints().length }.join("|") + "#" + Arrays.toString(rm.getSelectedIndexes()) + "#" + showBg.isSelected() + "#" + proposals.size() + "#" + showSug.isSelected() + "#" + edgeBox.isSelected() + '#'+selectionBasis + '#' + responseSeen + (analysisMeasurement!=null ? '#'+measuredKey(all) : '')
        if (!force && key == lastKey) return
        lastKey = key
        boolean measuredNow = false   // the cells on the picture are the ones in the saved tables
        if (analysisMeasurement != null) {
            boolean current = measuredKey(all)+'#edge='+edgeBox.isSelected()+'#selection='+basisFor(all) == analysisMeasurement.key
            measuredNow = current
            analyzeBtn.setEnabled(current && go.isEnabled())
            s3.setDone(current); s4.setSummary(current ? "ready" : "measure again first")
            analysisNote.setText("<html><body style='width:300px'>" + (current ? "Review the recording, check timing, and run analysis in your browser. Keep Fiji open." : "Cell outlines changed. Measure and save again before analysis.") + "</body></html>")
        }
        steps.alsoAvoid = edgeCells
        lastBg = steps.backgrounds(cells)
        int none = lastBg.count { it.length == 0 } as int
        def ov = new Overlay()
        if (showBg.isSelected() && cells.length > 0 && !hiddenKinds.contains('background')) {
            int sel = rm.getSelectedIndex(); int only = -1
            if (sel >= 0 && rm.getSelectedIndexes().length == 1) { def allR = rm.getRoisAsArray(); int k = -1; for (int i = 0; i <= sel && i < allR.length; i++) if (allR[i].isArea()) k++; only = sp.keptIdx.indexOf(k) }
            ov.add(steps.blueLayer(lastBg, only))
        }
        def left = unusedProposals(all)
        if (showSug.isSelected()) left.each { p ->
            if (hiddenKinds.contains('g' + Math.min(3, (int) p.grade))) return
            def r = (Roi) p.roi.clone()
            if (p.grade == 0) { r.setStrokeColor(CYAN); r.setStrokeWidth(2) }
            else if (p.grade == 1) { r.setStrokeColor(CYAN); r.setStrokeWidth(1) }
            else if (p.grade == 2) { r.setStrokeColor(GREY); r.setStroke(DASH) }
            else { r.setStrokeColor(RED); r.setStrokeWidth(2) }
            ov.add(r)
        }
        def selectedAreas = []; int areaIndex = 0
        rm.getRoisAsArray().eachWithIndex { r, i -> if (r.isArea()) { if (rm.isSelected(i)) selectedAreas << areaIndex; areaIndex++ } }
        if (!hiddenKinds.contains('cells')) all.eachWithIndex { c, i -> CalciumCellAppearance.add(ov, c, !sp.keptIdx.contains(i), selectedAreas.contains(i)) }
        steps.colourPic.setOverlay(ov)
        int nLikely = left.count { it.grade <= 1 } as int
        int edgeSug = edgeBox.isSelected() ? (proposals.count { it.edge } as int) : 0
        status.setText("<html><body style='width:300px'><b>Cells chosen: " + cells.length + "</b>" + (edgeCells.length > 0 ? " &nbsp;(+" + edgeCells.length + " at the edge, left out)" : "") + (none > 0 ? " &nbsp;(" + none + " without background)" : "") + (proposals.size() > 0 ? "<br><font color='" + MUTED + "'>Suggestions left: " + nLikely + " likely, " + left.count { it.grade == 2 } + " unsure, " + left.count { it.grade == 3 } + " very bright (precipitate?)" + (edgeSug > 0 ? "; " + edgeSug + " at the edge not shown" : "") + "</font>" : "") + "</body></html>")
        s2.setSummary(cells.length + (cells.length == 1 ? " cell chosen" : " cells chosen")); s2.setDone(cells.length > 0 && !s2.expanded)
        boolean withResponse = !steps.blind && steps.found && selectionBasis.startsWith('Response')
        acceptAllBtn.setText((withResponse ? "Accept responding and likely (" : "Accept likely (") + nLikely + ")"); acceptAllBtn.setEnabled(nLikely > 0); acceptAllBtn.setVisible(nLikely > 0)
        int nUnsure = left.count { it.grade == 2 && !(edgeBox.isSelected() && it.edge) } as int
        def hiddenNames = legendEntries.findAll { k, v -> hiddenKinds.contains(k) }.collect { k, v -> v.text }
        String hiddenNote = hiddenNames ? " <font color='#9A5B00'>Hidden on the picture now: " + hiddenNames.join(", ") + " (click them in the legend to show them).</font>" : ""
        addRespondersBtn.setVisible(steps.blind == true); addRespondersBtn.setEnabled(steps.blind == true && all.length > 0 && !detecting)
        // one blue button: the next thing to press (cosmetic, so a window without these buttons carries on)
        try {
            String blue = detecting ? '' : measuredNow ? '' : (proposals.isEmpty() && all.length == 0) || (!steps.blind && proposals.any { it.blind }) ? 'suggest' : nLikely > 0 ? 'accept' : all.length > 0 ? 'continue' : 'suggest'
            CalciumTheme.primary(suggestBtn, blue == 'suggest'); CalciumTheme.primary(acceptAllBtn, blue == 'accept'); CalciumTheme.primary(continueBtn, blue == 'continue')
        } catch (Throwable ignored) { }
        if (detecting) sayNext("wait while the cell outlines are made (see the bar below).")
        else if (measuredNow) sayNext("these " + cells.length + " cells are measured. Go on to <b>Analyze responses</b>, or change the cells here and measure again." + hiddenNote)
        else if (proposals.isEmpty() && all.length == 0) sayNext("press <b>Suggest cells automatically</b>, or draw the cells yourself." + hiddenNote)
        else if (!steps.blind && proposals.any { it.blind }) sayNext("press <b>Suggest cells again</b>: the outlines on the picture were made before the response, not with it." + hiddenNote)
        else if (nLikely > 0) sayNext("press <b>" + (withResponse ? "Accept responding and likely" : "Accept likely") + " (" + nLikely + ")</b>." + hiddenNote)
        else if (steps.blind) sayNext((nUnsure > 0 ? "click the dashed outlines that are real cells (" + nUnsure + " left). Then " : "") + "<b>Add cells that respond…</b> or <b>Continue to measure</b>." + hiddenNote)
        else sayNext((nUnsure > 0 ? "click the dashed outlines that are real cells (" + nUnsure + " left), then " : "") + "<b>Continue to measure</b>." + hiddenNote)
        undoBtn.setEnabled(!cellEdits.history.isEmpty())
    } catch (Throwable t) { status.setText("" + t) }
}
def paintLegend = { -> legendEntries.each { k, v -> boolean off = hiddenKinds.contains(k); v.label.setText(off ? "<html><strike>" + v.text + "</strike></html>" : v.text); v.label.setForeground(off ? new Color(0xB0B7C0) : CalciumTheme.MUTED) } }
legendClicked = { String kind, boolean only ->
    def kinds = legendEntries.keySet()
    if (only) {
        if (!hiddenKinds.remove(kind)) hiddenKinds.add(kind)                      // take back what the first click of the double-click did
        boolean alreadyOnly = !hiddenKinds.contains(kind) && hiddenKinds.size() == kinds.size() - 1
        hiddenKinds.clear(); if (!alreadyOnly) hiddenKinds.addAll(kinds.findAll { it != kind })
    } else if (!hiddenKinds.remove(kind)) hiddenKinds.add(kind)
    paintLegend(); refresh(true)
}
def viewButtons = [viewAllBtn, viewBaselineBtn, viewResponseBtn]
def busyView = { boolean busy -> viewButtons.each { it.setEnabled(!busy) } }
def currentView = { -> steps.blind ? 'baseline' : selectionBasis.startsWith('Response') ? 'response' : 'all' }
def blindControls = { boolean on -> respondsLegend.setVisible(!on); sensBox.setVisible(on); trace.setVisible(!on); frame.refreshContents() }
viewAllBtn.addActionListener({
    if (steps.blind) endBlindChoosing()          // the all-frames picture holds the response frames: cells added from here on are not blind
    steps.stopBlind(); steps.showAllFrames(); steps.hideResponse()
    selectionBasis = 'All-frames structural review'
    showView('all', 'Choosing cells on the average of all frames.'); blindControls(false); describeRecording(); refresh(true)
} as java.awt.event.ActionListener)
viewBaselineBtn.addActionListener({
    int[] range = CalciumStartDialog.askBlindFrames(steps); if (range == null) return
    busyView(true); selectionNote.setText("<html><body style='width:300px'>Making the picture of frames " + range[0] + "–" + range[1] + "…</body></html>")
    Thread.start('Calcium pre-stimulus picture') {
        float[] picture = null; String problem = null
        try { picture = steps.framesPicture(range[0], range[1]) } catch (Throwable t) { problem = t.message ?: t.toString() }
        javax.swing.SwingUtilities.invokeLater {
            busyView(false)
            if (picture == null) { showView(currentView(), 'Could not make the picture of the frames before the response: ' + problem); return }
            // Blind choosing starts again. If the response has not been shown, the cells added from now on are blind again; the
            // ones added in between are remembered as not blind. Once the response has been seen, nothing added is blind.
            if (!responseSeen && blindCellKeys != null) { notBlindKeys = (cellRois().collect { cellKey(it) } as Set) - blindCellKeys; blindCellKeys = null }
            steps.startBlind(range[0], range[1], picture)
            // suggestions made with the response (or on all frames) would show it: only blind suggestions stay
            proposals.removeAll { !it.blind }; if (proposals.isEmpty()) suggestBtn.setText("Suggest cells automatically")
            trace.show(null, "Traces are hidden while choosing cells blind", 0, 0); info.setText("<html><body style='width:300px'>&nbsp;</body></html>")
            selectionBasis = steps.blindBasis()
            showView('baseline', blindNote()); blindControls(true); describeRecording(); refresh(true)
        }
    }
} as java.awt.event.ActionListener)
viewResponseBtn.addActionListener({
    if (steps.blind && !testGui && javax.swing.JOptionPane.showConfirmDialog(frame, "<html><body style='width:320px'>The cells you have now stay recorded as <b>chosen blind</b>, as long as you do not change them.<br><br>The response is shown next and the plugin looks for the cells that only show when they respond. Cells you add from now on are recorded as <b>added after the response was shown</b>. Step 4 gives % responding for both groups.</body></html>", "Add cells that respond", javax.swing.JOptionPane.OK_CANCEL_OPTION) != javax.swing.JOptionPane.OK_OPTION) return
    boolean first = !steps.looked
    busyView(true); if (first) selectionNote.setText("<html><body style='width:300px'>Looking for the response in all " + steps.n + " frames…</body></html>")
    Thread.start('Calcium response picture') {
        String problem = null
        try { if (!steps.looked) steps.findResponse() } catch (Throwable t) { problem = t.message ?: t.toString() }
        javax.swing.SwingUtilities.invokeLater {
            busyView(false); describeRecording()
            String previous = currentView()
            if (problem != null) { showView(previous, 'Could not look for the response: ' + problem); return }
            if (!steps.found) { showView(previous, 'No response found (' + steps.why + '), so there is no magenta. The picture is unchanged.'); return }
            reveal(); steps.stopBlind(); steps.showAllFrames(); steps.showResponse()
            selectionBasis = 'Response-assisted review'
            showView('response', 'Magenta: cells that brightened near frame ' + steps.onset + '. Picking cells this way is recorded as response-assisted.'); blindControls(false); describeRecording(); refresh(true)
            // The outlines on the picture were made before the response. Look for the cells that respond straight away, so
            // nobody has to know that "Suggest cells again" is the next press. (Saved outlines of the same picture are reused.)
            if (proposals.any { it.blind }) { autoSuggestAsked = true; if (!testGui) suggestBtn.doClick() }
        }
    }
} as java.awt.event.ActionListener)
addRespondersBtn.addActionListener({ e -> viewResponseBtn.doClick() } as java.awt.event.ActionListener)
def addCell = { Roi r0 -> def r = (Roi) r0.clone(); r.setPosition(0); r.setStrokeColor(null); rm.addRoi(r) }
def acceptAt = { int x, int y ->
    def hit = unusedProposals(cellRois()).findAll { !hiddenKinds.contains('g' + Math.min(3, (int) it.grade)) && it.roi.contains(x, y) }
    if (hit.isEmpty()) return false
    cellEdits.change('add cell') { addCell(hit.min { it.roi.getBounds().width * it.roi.getBounds().height }.roi) }
    rm.deselect(); steps.colourPic.deleteRoi()
    refresh(true); return true
}
def acceptAll = { -> def todo = unusedProposals(cellRois()).findAll { it.grade <= 1 }; cellEdits.change('accept all likely cells') { todo.each { addCell(it.roi) } }; rm.deselect(); steps.colourPic.deleteRoi(); refresh(true); todo.size() }
def lastHover = null
def hoverAt = { int x, int y ->
    Roi[] cells = cellRois()
    def p = showSug.isSelected() ? unusedProposals(cells).findAll { !hiddenKinds.contains('g' + Math.min(3, (int) it.grade)) && it.roi.contains(x, y) }.min { it.roi.getBounds().width * it.roi.getBounds().height } : null
    def target = null; def text = "&nbsp;"; def cap = ""
    if (p != null) { target = p.roi; text = evidenceText(p); cap = "this suggestion: gray levels above its background, over time" }
    else {
        int ci = -1; cells.eachWithIndex { c, i -> if (ci < 0 && c.contains(x, y)) ci = i }
        if (ci >= 0) { target = cells[ci]; text = "<b>Your cell " + (ci + 1) + "</b>"; cap = "cell " + (ci + 1) + ": gray levels above its background, over time" }
    }
    if (target == null || target.is(lastHover)) { if (target == null && lastHover != null) { lastHover = null; info.setText("<html><body style='width:300px'>&nbsp;</body></html>") }; return text }
    lastHover = target
    info.setText("<html><body style='width:300px'>" + text + "</body></html>")
    if (steps.blind) trace.show(null, "Traces are hidden while choosing cells blind", 0, 0)   // a trace would show the response
    else { reveal(); trace.show(steps.traceOf(target, cells), cap, steps.found ? steps.p0 : 0, steps.found ? steps.p1 : 0) }
    return text
}
def findCellposePython = { ->
    def home = System.getProperty("user.home"); def cands = [Prefs.get("calcium.steps.cellpose_python", "")]
    try {                                    // a text file next to the plugin may hold the path (one line)
        def nf = new File(IJ.getDirectory("imagej"), "scripts/Plugins/Calcium Imaging/cellpose_python.txt")
        if (nf.exists()) cands << nf.text.trim()
    } catch (Throwable ignore) { }
    ["cellpose-env", ".venvs/cellpose", "miniconda3/envs/cellpose", "anaconda3/envs/cellpose", "miniforge3/envs/cellpose", "mambaforge/envs/cellpose"].each { d ->
        cands << home + "/" + d + "/bin/python"; cands << home + "/" + d + "/Scripts/python.exe"; cands << home + "/" + d + "/python.exe"
    }
    cands.find { it && new File(it).isFile() && new File(it).canExecute() }
}
// One-click setup of Cellpose in the user's home folder. Everything goes into two folders of its own
// (cellpose-env and cellpose-env-tools); deleting them removes it again. Returns "" when it worked, else what went wrong.
def setupCellpose = { File envDir, File toolsDir, boolean offline, Closure progress ->
    def logFile = new File(toolsDir, "setup_log.txt"); toolsDir.mkdirs()
    def run = { List cmd, String what ->
        progress(what); logFile << ("\n== " + what + "\n" + cmd.join(" ") + "\n")
        def pb = new ProcessBuilder(cmd as String[]); pb.redirectErrorStream(true); pb.redirectOutput(ProcessBuilder.Redirect.appendTo(logFile))
        if (offline) pb.environment().put("UV_OFFLINE", "1")
        else { pb.environment().put("UV_CACHE_DIR", new File(toolsDir, "cache").getPath()); pb.environment().put("UV_PYTHON_INSTALL_DIR", new File(toolsDir, "python").getPath()) }   // keep everything inside the two folders
        pb.start().waitFor()
    }
    boolean win = System.getProperty("os.name").toLowerCase().contains("win"), mac = System.getProperty("os.name").toLowerCase().contains("mac")
    def arch = System.getProperty("os.arch").toLowerCase(); boolean arm = arch.contains("aarch64") || arch.contains("arm")
    def home = System.getProperty("user.home")
    // 1. uv: a small free program that installs Python packages (https://github.com/astral-sh/uv)
    def uvName = win ? "uv.exe" : "uv"
    def uv = [new File(toolsDir, uvName), new File(home, ".local/bin/" + uvName), new File(home, ".cargo/bin/" + uvName), new File("/opt/homebrew/bin/uv"), new File("/usr/local/bin/uv")].find { it.isFile() && it.canExecute() && (it.getParentFile() == toolsDir || System.getProperty("calcium.steps.own_uv_only") == null) }
    if (uv == null) {
        if (offline) return "uv is not on this computer and downloading is switched off"
        def asset = win ? "uv-x86_64-pc-windows-msvc.zip" : mac ? (arm ? "uv-aarch64-apple-darwin.tar.gz" : "uv-x86_64-apple-darwin.tar.gz") : (arm ? "uv-aarch64-unknown-linux-gnu.tar.gz" : "uv-x86_64-unknown-linux-gnu.tar.gz")
        def url = "https://github.com/astral-sh/uv/releases/latest/download/" + asset
        progress("Downloading the installer uv (about 17 MB)..."); logFile << ("\n== download " + url + "\n")
        def pack = new File(toolsDir, asset)
        try { new URL(url).withInputStream { i -> pack.withOutputStream { o -> o << i } } } catch (Throwable t) { return "could not download " + url + " (" + t + ")" }
        // compare with the checksum that is published next to the file
        try {
            def want = new URL(url + ".sha256").getText().trim().split(/\s+/)[0].toLowerCase()
            def md = java.security.MessageDigest.getInstance("SHA-256"); pack.withInputStream { i -> byte[] buf = new byte[1 << 16]; int k; while ((k = i.read(buf)) > 0) md.update(buf, 0, k) }
            def have = md.digest().collect { String.format("%02x", it) }.join("")
            logFile << ("checksum published " + want + "\nchecksum of file   " + have + "\n")
            if (want.length() == 64 && want != have) { pack.delete(); return "the downloaded installer does not match its published checksum; nothing was installed" }
        } catch (Throwable t) { logFile << ("checksum could not be compared: " + t + "\n") }
        if (win) {
            def zf = new java.util.zip.ZipFile(pack)
            zf.entries().each { e -> if (!e.isDirectory()) { def out = new File(toolsDir, new File(e.getName()).getName()); out.bytes = zf.getInputStream(e).bytes } }
            zf.close()
        } else if (run(["tar", "-xzf", pack.getPath(), "-C", toolsDir.getPath(), "--strip-components", "1"], "Unpacking uv") != 0) return "could not unpack " + pack.getPath()
        uv = new File(toolsDir, uvName); uv.setExecutable(true)
        if (!uv.isFile()) return "uv was not found after unpacking"
    }
    // 2. a Python of its own with Cellpose in it
    if (run([uv.getPath(), "venv", "--python", "3.11", envDir.getPath()], "Making a Python environment for Cellpose (downloads Python 3.11 if needed)...") != 0) return "could not make the Python environment; see " + logFile.getPath()
    def py = new File(envDir, win ? "Scripts/python.exe" : "bin/python")
    // "packaging" is used by a part of Cellpose (fastremap) but is not installed with it automatically
    // On Windows and Linux "--torch-backend auto" picks the version that can use an NVIDIA graphics card when there is one
    // (a much bigger download in that case); if that fails, the plain install is tried.
    int inst = mac ? 1 : run([uv.getPath(), "pip", "install", "--python", py.getPath(), "--torch-backend", "auto", "cellpose", "packaging"], "Downloading and installing Cellpose (the long part)...")
    if (inst != 0 && run([uv.getPath(), "pip", "install", "--python", py.getPath(), "cellpose", "packaging"], "Downloading and installing Cellpose (the long part)...") != 0) return "could not install Cellpose; see " + logFile.getPath()
    if (run([py.getPath(), "-m", "cellpose", "--version"], "Checking that Cellpose starts...") != 0) return "Cellpose was installed but does not start; see " + logFile.getPath()
    if (!offline) run([uv.getPath(), "cache", "clean"], "Tidying up...")   // the installed copy does not need the download cache
    return ""
}
def cellposeEnv = new File(System.getProperty("user.home"), "cellpose-env"), cellposeTools = new File(System.getProperty("user.home"), "cellpose-env-tools")
def SETUP_TEXT = "Cellpose is the free program that suggests the cell outlines. It is separate from Fiji.\n\n" +
    "Setting it up downloads from the internet:\n" +
    "   - uv, a small installer (about 17 MB, from github.com/astral-sh/uv), unless it is already here\n" +
    "   - Python 3.11 and Cellpose with the packages it needs (about 250 MB on a Mac, more on\n" +
    "     Linux; from github.com/astral-sh/python-build-standalone and pypi.org)\n" +
    "and the first time you press \"Suggest cells\" Cellpose downloads its model (about 1.2 GB).\n\n" +
    "Everything goes into two folders in your home folder:  cellpose-env  and  cellpose-env-tools\n" +
    "(the model goes into a third one, called  .cellpose ). Nothing else on the computer is changed.\n" +
    "It needs about 1 GB of disk space, plus the model. To remove it later, delete those folders.\n\n" +
    "It takes a few minutes, depending on the internet connection. You can keep choosing cells\n" +
    "by hand meanwhile."

def cellposeBusy = false
def showCellposeState = { String busyText ->
    def found = findCellposePython() != null
    cpLabel.setText("<html><body style='width:300px'><font color='" + MUTED + "'>" + (busyText != null ? busyText : found ? "Cellpose: found on this computer." : "Cellpose (the program that suggests outlines) is not installed on this computer yet.") + "</font></body></html>")
    cpLabel.setVisible(!found || busyText != null)   // nothing to say once Cellpose is installed and idle
    setupRow.setVisible(!found || busyText != null); setupBtn.setEnabled(busyText == null)
}
// asks first (unless already asked), then installs; returns "" when Cellpose is ready
def runSetup = { boolean ask, Closure detectionUpdate = null ->
    if (ask && javax.swing.JOptionPane.showConfirmDialog(frame, SETUP_TEXT + "\n\nSet up Cellpose now?", "Set up Cellpose", javax.swing.JOptionPane.OK_CANCEL_OPTION, javax.swing.JOptionPane.QUESTION_MESSAGE) != javax.swing.JOptionPane.OK_OPTION) return "cancelled"
    cellposeBusy = true
    def problem = ""
    try {
        problem = setupCellpose(cellposeEnv, cellposeTools, false, { String s ->
            if (detectionUpdate != null) detectionUpdate.call(s)
            javax.swing.SwingUtilities.invokeLater { showCellposeState("Setting up Cellpose: " + s) }
        })
    } catch (Throwable t) { problem = t.toString() }
    cellposeBusy = false
    steps.say(problem == "" ? "Cellpose was set up by the plugin in " + cellposeEnv.getPath() : "Cellpose setup failed: " + problem)
    javax.swing.SwingUtilities.invokeLater {
        showCellposeState(problem == "" ? null : "Cellpose could not be set up: " + problem + "  (details in " + new File(cellposeTools, "setup_log.txt").getPath() + ")")
        if (problem != "") setupBtn.setEnabled(true)
    }
    return problem
}

def forceFresh = false     // set when "Suggest cells again" is pressed: do not reuse saved outlines
// comparison (the sensitivity check): [picture, dir, into] runs Cellpose on that picture into that folder, grades the outlines
// with the response, and puts them in "into" instead of the suggestions. Nothing is shown and nothing is accepted.
def suggest = { CalciumOutlineTask task = new CalciumOutlineTask(progress: { p, s, d -> }), Map comparison = null ->
    task.report(0, "Checking Cellpose", "Finding the installed Python environment.")
    def py = findCellposePython()
    if (py == null && cellposeBusy) return "Cellpose is still being set up; try again when it says it is found"
    if (py == null && comparison != null) return "Cellpose was not found"
    if (py == null && !testGui) {
        def choices = ["Set up Cellpose now", "I already have it: show where", "Cancel"] as Object[]
        int ans = javax.swing.JOptionPane.showOptionDialog(frame, "Cellpose was not found on this computer.\n\n" + SETUP_TEXT, "Cellpose is needed for suggestions", javax.swing.JOptionPane.DEFAULT_OPTION, javax.swing.JOptionPane.QUESTION_MESSAGE, null, choices, choices[0])
        if (ans == 0) {
            javax.swing.SwingUtilities.invokeLater { suggestBtn.setText("Setting up Cellpose first (a few minutes)...") }
            task.cancellable = false
            def problem
            try { problem = runSetup(false, { String s -> task.report(0, "Setting up Cellpose", s) }) }
            finally { task.cancellable = true }
            task.checkCancelled()
            if (problem != "") return "Cellpose could not be set up: " + problem
            py = findCellposePython()
            javax.swing.SwingUtilities.invokeLater { suggestBtn.setText("Looking for cell outlines...") }
        } else if (ans == 1) {
            def g2 = new GenericDialog("Where is Cellpose?")
            g2.addMessage("Give the 'python' file of the environment in which Cellpose is installed\n(for example  ~/cellpose-env/bin/python ).")
            g2.addFileField("Python with Cellpose:", "", 40)
            g2.showDialog(); if (g2.wasCanceled()) return "cancelled"
            py = g2.getNextString(); if (!new File(py).canExecute()) { IJ.error("Calcium imaging", "That file cannot be run:\n" + py); return "that python file cannot be run" }
            Prefs.set("calcium.steps.cellpose_python", py)
        } else return "cancelled"
    }
    if (py == null) return "Cellpose was not found"
    task.report(5, "Preparing the detection image", "Checking whether saved outlines can be reused.")
    boolean blindRun = comparison == null && steps.blind == true
    def dir = comparison != null ? (File) comparison.dir : new File(outDir, blindRun ? "automatic_outlines_blind" : "automatic_outlines"); dir.mkdirs()
    // Blind: Cellpose runs on two pictures, 0 to BLIND_CAP and 0 to BLIND_CAP_HIGH gray levels, each in its own folder, and the
    // outlines are combined (CalciumSteps.blindOutlinesToAdd). Otherwise one picture, in dir.
    def runs = blindRun ? (0..1).collect { int i ->
            double cap = i == 0 ? CalciumSteps.BLIND_CAP : CalciumSteps.BLIND_CAP_HIGH; String lv = CalciumSteps.levels(cap)
            [label: "picture " + (i + 1) + " (0-" + lv + " gray levels)", part: "Blind picture " + (i + 1) + " of 2 (0-" + lv + " gray levels)", shown: "0 to " + lv,
             dir: new File(dir, "picture_" + (i + 1) + "_shown_0_to_" + lv + "_gray_levels"), pic: steps.blindOutlinePicture(cap)] }
        : [[label: "", part: "", shown: "", dir: dir, pic: comparison != null ? (ImagePlus) comparison.picture : steps.outlinePicture()]]
    String pictureNote = blindRun ? steps.blindBasis() + ". Pictures for Cellpose: " + steps.blindRecipe() + "; no response channel. Folders: " + runs.collect { ((File) it.dir).getName() }.join(", ") + ".\n" : ""
    if (blindRun) new File(dir, "detection_picture.txt").setText(pictureNote, "UTF-8")
    // One Cellpose run on one picture, in that picture's folder. Saved outlines are reused when the saved picture matches.
    // Gives [zip, how, cmd, pending], or [error] when Cellpose made no outlines.
    def detect = { Map run ->
        File d = (File) run.dir; d.mkdirs()
        ImagePlus pic = (ImagePlus) run.pic
        String part = run.part ? run.part + "; " : ""
        task.settings = run.part; task.device = ""; task.downloading = false
        def tif = new File(d, base + "_for_automatic_outlines.tif")
        def zip = new File(d, base + "_for_automatic_outlines_rois.zip")
        def optFile = new File(d, "cellpose_options_used.txt")
        def pending = new File(d, "cellpose_run_incomplete.txt")
        task.imageSize = pic.getWidth() + " x " + pic.getHeight() + " pixels"
        def opts = []; def cmd = []
        boolean reuse = false
        if (!(forceFresh && comparison == null) && !pending.exists() && tif.isFile() && zip.isFile() && zip.length() > 0) {
            ImagePlus saved = null
            try { saved = IJ.openImage(tif.getPath()); reuse = CalciumSteps.samePicture(saved, pic) }
            catch (Throwable ignore) { }
            finally { if (saved != null) saved.close() }
        }
        String how
        if (reuse) {
            task.source = "Saved outlines reused (no new Cellpose run)"
            how = "the outlines made earlier for this recording were used again" + (optFile.exists() ? "" : " (made by an earlier version of the plugin)")
            cmd = ["(no new run; earlier settings: " + (optFile.exists() ? optFile.text.trim() : "full size") + ")"]
            def previousLog = new File(d, "cellpose_output.txt")
            if (previousLog.isFile()) task.logFile = previousLog
            task.report(84, "Reusing saved cell outlines", "The detection image matches this recording; a new Cellpose run is not needed.")
        } else {
            task.source = "New detection"
            task.checkCancelled()
            // A failed/cancelled attempt may leave a partial ZIP. Never treat it as a reusable result.
            pending.text = "Detection has not completed successfully. Do not reuse its outlines.\n"
            def removeOldZip = { -> if (zip.exists() && !zip.delete()) throw new IOException("Cannot replace " + zip.getPath()) }
            removeOldZip()
            if (optFile.exists() && !optFile.delete()) throw new IOException("Cannot replace " + optFile.getPath())
            if (!IJ.saveAsTiff(pic, tif.getPath())) throw new IOException("Could not save the detection image")
            long started = System.nanoTime()
            def key = "calcium.steps.cellpose_gpu." + Integer.toHexString(py.hashCode())
            String known = Prefs.get(key, "")
            if (System.getProperty("calcium.steps.pretend_no_gpu") != null) known = "CPU"
            if (known == "") {
                task.report(8, "Checking CPU / graphics support", "Checking whether this Cellpose installation can use a graphics processor.")
                // A script file avoids Windows command-line quoting differences for Python -c.
                def probe = new File(d, "cellpose_device_check.py")
                probe.text = "import torch\nprint('CHIP=GPU' if (torch.cuda.is_available() or (hasattr(torch.backends, 'mps') and torch.backends.mps.is_available())) else 'CHIP=CPU')\n"
                try {
                    def probeLog = new File(d, "cellpose_device_check.txt")
                    int probeCode = task.run([py, "-u", probe.getAbsolutePath()], d, probeLog)
                    String answer = probeLog.getText("UTF-8")
                    if (probeCode == 0) known = answer.contains("CHIP=GPU") ? "GPU" : answer.contains("CHIP=CPU") ? "CPU" : ""
                    if (known != "") Prefs.set(key, known)
                } catch (java.util.concurrent.CancellationException cancelled) { throw cancelled }
                catch (Exception ignored) { } // Keep the uploaded version's GPU-first behaviour if the check is inconclusive.
            }
            task.checkCancelled()
            task.device = known
            boolean gpu = known != "CPU"
            task.settings = part + "Attempt 1 of 2; " + (gpu ? "full-size detection" : "faster CPU setting (diameter 60)")
            opts = ["--pretrained_model", "cpsam_v2", "--use_gpu"] + (gpu ? [] : ["--diameter", "60"]) + ["--save_rois", "--no_npy", "--cellprob_threshold", "-2", "--flow_threshold", "0.8", "--verbose"]
            cmd = [py, "-u", "-m", "cellpose", "--image_path", tif.getPath()] + opts
            task.report(10, "Starting Cellpose (attempt 1 of 2)", (run.label ? "Blind " + run.label + ". " : "") + (gpu ? "Using full-size detection. First use may download the model." : "CPU detection uses the faster diameter-60 setting from this version; it may still take several minutes."))
            int code = task.run(cmd, d, new File(d, "cellpose_output.txt"))
            if (code != 0 || !zip.isFile() || zip.length() == 0) {
                task.checkCancelled(); removeOldZip()
                task.downloading = false; task.device = "CPU"
                task.settings = part + "Attempt 2 of 2; default model, full-size detection"
                task.report(task.lastAsked, "Retrying Cellpose on the CPU (attempt 2 of 2)", "The first attempt did not produce outlines. Trying this version's default model and full-size image.")
                opts = ["--save_rois", "--no_npy", "--cellprob_threshold", "-2", "--flow_threshold", "0.8", "--verbose"]
                cmd = [py, "-u", "-m", "cellpose", "--image_path", tif.getPath()] + opts
                code = task.run(cmd, d, new File(d, "cellpose_output_second_try.txt"))
            }
            task.checkCancelled()
            if (code != 0 || !zip.isFile() || zip.length() == 0) return [error: "Cellpose made no outlines" + (run.label ? " on blind " + run.label : "") + " (exit code " + code + "). Open the detection log: " + task.logFile.getPath()]
            optFile.text = opts.join(" ")
            how = "Cellpose took " + (long) ((System.nanoTime() - started) / 1000000000L) + " s" + (task.device ? " (" + task.device + ")" : "") +
                (opts.contains("--diameter") ? "; used the faster CPU diameter-60 setting" : "; used full-size detection")
        }
        return [zip: zip, how: how, cmd: cmd, pending: pending]
    }
    def results = []
    try {
        for (int i = 0; i < runs.size(); i++) {
            if (runs.size() > 1) { task.partFrom = 42 * i; task.partShare = 0.5 }   // two runs: 0-42 % and 42-84 % of the bar
            def done = detect(runs[i])
            if (done.error) return done.error
            results << done
        }
    } finally { task.partFrom = 0; task.partShare = 1 }
    String how = runs.size() == 1 ? results[0].how : (0..<runs.size()).collect { runs[it].label + ": " + results[it].how }.join("; ")
    String cmdText = runs.size() == 1 ? results[0].cmd.join(" ") : (0..<runs.size()).collect { runs[it].label + ": " + results[it].cmd.join(" ") }.join("; ")
    task.report(85, "Reading the cell outlines", "Opening the outlines and checking their measurements.")
    def readOutlines = { File z ->
        def list = []
        def zf = new java.util.zip.ZipFile(z)
        try {
            zf.entries().toList().findAll { !it.isDirectory() && it.getName().endsWith(".roi") }.sort { it.getName() }.each { e ->
                task.checkCancelled()
                def r = new ij.io.RoiDecoder(zf.getInputStream(e).withCloseable { it.bytes }, e.getName()).getRoi()
                if (r != null && r.isArea()) list << r
            }
        } finally { zf.close() }
        return list
    }
    def outlines = readOutlines((File) results[0].zip).collect { [roi: it, from: 0] }
    def gradeOne = { r -> comparison != null ? steps.grade(r, false) : steps.grade(r) }
    String mergeNote = ""
    if (blindRun) {
        // Only an outline of picture 1 that will be shown may stand in for an outline of picture 2. One that is too small or too
        // large for a cell, or a thin strip along another outline, is never shown and must not hide a good outline (seen 20 Nov
        // 2026 on P2 Run2: one 6445-pixel outline hid six cells of picture 2, and an 11-pixel speck hid one of 1731 pixels).
        outlines.each { it.ev = gradeOne(it.roi) }
        def first = outlines.collect { it.roi }, second = readOutlines((File) results[1].zip)
        def shownFirst = outlines.findAll { it.ev[0] < 4 && !(CalciumSteps.tooThin(it.roi) && CalciumSteps.alongAnother(it.roi, first)) }.collect { it.roi }
        def add = CalciumSteps.blindOutlinesToAdd(shownFirst, second, ((ImagePlus) runs[0].pic).getWidth())
        add.each { int j -> outlines << [roi: second[j], from: 1] }
        int unshown = first.size() - shownFirst.size()
        mergeNote = "all " + first.size() + " outlines from " + runs[0].label + " kept; " + add.size() + " of " + second.size() + " from " + runs[1].label +
            " added (the other " + (second.size() - add.size()) + " share at least " + (int) Math.round(CalciumSteps.BLIND_MERGE_OVERLAP * 100) + " % of the smaller outline's area with an outline from picture 1" +
            (unshown > 0 ? "; " + unshown + " outlines of picture 1 that are not shown, for their size or as a strip on another outline, were not counted" : "") + ")"
        new File(dir, "detection_picture.txt").setText(pictureNote + "Combined: " + mergeNote + ".\n", "UTF-8")
    }
    def got = []; int hidden = 0, strips = 0, thinKept = 0
    List everyOutline = outlines.collect { it.roi }
    outlines.eachWithIndex { o, i ->
        task.checkCancelled()
        def r = o.roi
        double[] ev = o.ev != null ? (double[]) o.ev : gradeOne(r)
        boolean thin = ev[0] < 4 && CalciumSteps.tooThin(r)
        if (ev[0] >= 4) hidden++
        else if (thin && CalciumSteps.alongAnother(r, everyOutline)) strips++   // the rim of a cell that has its own outline
        else {
            // a thin outline standing alone stays, as "unsure": it is never accepted without a click
            int grade = (int) ev[0]; if (thin) { thinKept++; grade = Math.max(grade, 2); ev[0] = grade }
            def cen = r.getContourCentroid(); got << [roi:r, grade:grade, ev:ev, cx:(int) cen[0], cy:(int) cen[1], edge:steps.nearEdge(r), blind:blindRun, picture:runs[o.from].shown, thin:thin]
        }
        task.report(85 + (int) (12 * (i + 1) / outlines.size()), "Checking the suggested cells", "Checked " + (i + 1) + " of " + outlines.size() + " outlines.")
    }
    task.report(98, "Saving the suggestion measurements", "Saving the evidence behind each outline's grade.")
    def csv = new StringBuilder("suggestion,grade,area_pixels,mean_brightness,share_very_bright,brighter_than_surroundings,rise_inside,rise_around,touches_edge,x,y" + (blindRun ? ",found_on_picture_gray_levels" : "") + "\n")
    got.eachWithIndex { p, i -> csv.append(String.format(Locale.US, "%d,%s,%.0f,%.2f,%.3f,%.2f,%.2f,%.2f,%s,%d,%d", i + 1, GRADE_TEXT[p.grade], p.ev[1], p.ev[2], p.ev[3], p.ev[4], p.ev[5], p.ev[6], p.ev[7] > 0 ? "yes" : "no", p.cx, p.cy) + (blindRun ? "," + p.picture : "") + "\n") }
    new File(dir, base + "_suggestions_with_evidence.csv").text = csv.toString()
    if (blindRun) {   // the combined suggestions, named as in the evidence table, so they can be opened in the ROI Manager
        def zout = new java.util.zip.ZipOutputStream(new FileOutputStream(new File(dir, base + "_suggestions_rois.zip")))
        try {
            got.eachWithIndex { p, i ->
                def copy = p.roi.clone(); copy.setName(String.format(Locale.US, "suggestion_%03d", i + 1))
                zout.putNextEntry(new java.util.zip.ZipEntry(copy.getName() + ".roi")); zout.write(ij.io.RoiEncoder.saveAsByteArray(copy)); zout.closeEntry()
            }
        } finally { zout.close() }
    }
    def finishRecord = { -> results.each { File pending = (File) it.pending; if (pending.exists() && !pending.delete()) throw new IOException("Could not finish the detection record: " + pending.getPath()) } }
    if (comparison != null) {
        task.checkCancelled(); comparison.into.addAll(got); task.succeeded = true
        finishRecord()
        return got.size() + " response-assisted outlines; " + how
    }
    if (blindRun) steps.say("Automatic outlines are blind: Cellpose got two pictures of frames " + steps.blindFirst + "-" + steps.blindLast + " only (" + steps.blindRecipe() + "), saved as " + ((File) results[0].zip).getName().replace("_rois.zip", ".tif") + " in " + dir.getName() + "/" + runs.collect { ((File) it.dir).getName() }.join(" and " + dir.getName() + "/") + "; " + mergeNote + ". The outlines were graded on the pre-stimulus picture without the response.")
    steps.say("Automatic outlines (suggestions only): Cellpose, command: " + cmdText + " (" + how + "). " + got.size() + " outlines shown: " + got.count { it.grade == 0 } + " respond on their own, " + got.count { it.grade == 1 } + " likely cells, " + got.count { it.grade == 2 } + " unsure, " + got.count { it.grade == 3 } + " very bright; " + hidden + " left out for their size; " + strips + " thin strips along another outline's edge left out, " + thinKept + " thin outlines kept as unsure (narrower than " + (int) CalciumSteps.THIN_WIDTH + " pixels on average, or roundness under " + CalciumSteps.THIN_ROUNDNESS + "); " + got.count { it.edge } + " lie within " + CalciumSteps.EDGE_MARGIN + " pixels of the picture's edge" + (edgeBox.isSelected() ? " and are not shown (box ticked)" : "") + ". Grades are measured (brightness, share of very bright pixels, contrast" + (blindRun ? "; the response is not used" : ", rise after the stimulus inside against around") + "). Only the outlines the user accepted became cells.")
    task.checkCancelled()
    javax.swing.SwingUtilities.invokeAndWait { proposals.clear(); proposals.addAll(got) }
    finishRecord()
    task.succeeded = true
    if (!blindRun) return got.size() + " suggestions; " + how
    task.settings = "Two pictures of the frames before the response (0-" + CalciumSteps.levels(CalciumSteps.BLIND_CAP) + " and 0-" + CalciumSteps.levels(CalciumSteps.BLIND_CAP_HIGH) + " gray levels)"
    def took = results.collect { (it.how =~ /^Cellpose took (\d+) s/) }.collect { it.find() ? it.group(1) : null }
    String shortHow = took.every { it != null } ? "Cellpose took " + took.join(" s and ") + " s" + (task.device ? " (" + task.device + ")" : "") :
        results.every { it.how.startsWith("the outlines made earlier") } ? "the outlines made earlier for this recording were used again" : how
    return got.size() + " suggestions from two pictures; " + shortHow
}


def detectionRunner = new CalciumDetectionRunner()
suggestBtn.addActionListener({ e ->
    // "again" asks for a new Cellpose run only when the outlines on screen came from the same kind of picture as now
    forceFresh = proposals.any { (it.blind == true) == (steps.blind == true) }
    detecting = true
    detectionRunner.start(suggestBtn, detectionProgress,
        { task -> suggest(task) },
        { String message ->
            detecting = false
            if (hiddenKinds.removeAll(['g0', 'g1', 'g2', 'g3'])) paintLegend()      // new suggestions are never made invisible by an old legend click
            suggestBtn.setText(proposals.isEmpty() ? "Suggest cells automatically" : "Suggest cells again")
            // a finished run is already reported under the button; only a problem is repeated here, where the eye is
            info.setText("<html><body style='width:300px'>" + (message ==~ /\d+ suggestions.*/ ? "&nbsp;" : "Suggestions: " + CalciumOutlineProgress.escape(message)) + "</body></html>")
            if (!cellposeBusy) showCellposeState(null)
            refresh(true)
        },
        { -> frame.refreshContents() })
    try { refresh(true) } catch (Throwable ignore) { }       // only so the Next line says to wait; the worker is already started
} as java.awt.event.ActionListener)
acceptAllBtn.addActionListener({ e ->
    int nL = unusedProposals(cellRois()).count { it.grade <= 1 } as int
    if (nL == 0) return
    int ans = javax.swing.JOptionPane.showConfirmDialog(frame, "Accept all " + nL + " outlines graded " + (steps.blind ? "" : "\"responds\" or ") + "\"likely a cell\"?\n\nIn tests about 8 in 10 of these were real cells, so about 2 in 10 are not.\nLook through them afterwards and remove the wrong ones.", "Accept all likely ones", javax.swing.JOptionPane.OK_CANCEL_OPTION)
    if (ans == javax.swing.JOptionPane.OK_OPTION) acceptAll()
} as java.awt.event.ActionListener)
def updateEditHint = { String message -> editHint.setText("<html><body style='width:300px'>" + message + "</body></html>") }
def leaveRemoveMode = { -> removeMode.setSelected(false); updateEditHint('Click a white outline to select it, then Remove. Or turn on <b>Click cells to remove</b>.') }
ovalBtn.addActionListener({ e -> leaveRemoveMode(); IJ.setTool("oval") } as java.awt.event.ActionListener)
freeBtn.addActionListener({ e -> leaveRemoveMode(); IJ.setTool("freehand") } as java.awt.event.ActionListener)
addBtn.addActionListener({ e -> def r = steps.colourPic.getRoi()
    if (r != null && r.isArea() && edgeBox.isSelected() && steps.nearEdge(r)) { info.setText("<html><body style='width:300px'><b>Not added:</b> this outline is within " + CalciumSteps.EDGE_MARGIN + " pixels of the picture's edge, where a cell's background can only come from one side. Untick \"Leave out cells at the picture's edge\" to add it anyway.</body></html>"); return }
    if (r != null && r.isArea()) { cellEdits.change('add drawn cell') { addCell(r) }; rm.deselect(); steps.colourPic.deleteRoi(); refresh(true) } else info.setText("<html><body style='width:300px'>Draw around a cell on the colour picture first.</body></html>") } as java.awt.event.ActionListener)
def removedMessage = { int count -> updateEditHint('Removed ' + count + (count == 1 ? ' cell.' : ' cells.') + ' <b>Undo last change</b> restores them.' + (removeMode.isSelected() ? '<br><b>Removal mode is ON:</b> click another white outline to remove it. Turn it off to accept suggestions.' : '')) }
def removeSelected = { -> int count = cellEdits.remove(rm.getSelectedIndexes() as List); refresh(true); if (count > 0) removedMessage(count); count > 0 }
def undoLast = { -> updateEditHint(cellEdits.undo() + (removeMode.isSelected() ? '<br><b>Removal mode is ON.</b> Turn it off to accept suggestions.' : '')); refresh(true) }
removeBtn.addActionListener({ e -> if (!removeSelected()) updateEditHint('Click inside a white cell outline first, or turn on <b>Click cells to remove</b>. You can also select cells in the ROI Manager list.') } as java.awt.event.ActionListener)
undoBtn.addActionListener({ e -> undoLast() } as java.awt.event.ActionListener)
removeMode.addActionListener({ e ->
    if (removeMode.isSelected()) { IJ.setTool('oval'); steps.colourPic.deleteRoi(); updateEditHint('<b>Removal mode is ON:</b> click inside each white cell outline to remove it. Undo restores mistakes. Turn this off to accept suggestions.') }
    else leaveRemoveMode()
} as java.awt.event.ActionListener)
def clickCellAt = { int x, int y ->
    if (removeMode.isSelected()) {
        int count = cellEdits.remove([cellEdits.hit(x, y)])
        steps.colourPic.deleteRoi(); refresh(true)
        if (count > 0) removedMessage(count)
        return count > 0
    }
    if (cellEdits.selectAt(x, y)) { refresh(true); updateEditHint('<b>Cell selected</b> (thicker white outline). Press <b>Remove selected cell</b> to remove it.'); return true }
    boolean added = showSug.isSelected() && acceptAt(x, y)
    refresh(true)
    if (added) updateEditHint('<b>Cell added:</b> its outline is now white with a black border. Undo last change reverses this.')
    added
}
showSug.addActionListener({ e -> refresh(true) } as java.awt.event.ActionListener)
edgeBox.addActionListener({ e -> refresh(true) } as java.awt.event.ActionListener)
showBg.addActionListener({ e -> refresh(true) } as java.awt.event.ActionListener)
if (steps.colourPic.getCanvas() != null) {
    def canvas = steps.colourPic.getCanvas(); int px0 = 0, py0 = 0; long lastMove = 0
    def cellMouse = new java.awt.event.MouseAdapter() {
        void mousePressed(java.awt.event.MouseEvent e) { px0 = e.getX(); py0 = e.getY() }
        void mouseReleased(java.awt.event.MouseEvent e) {
            if (!javax.swing.SwingUtilities.isLeftMouseButton(e) || e.isPopupTrigger() || e.getClickCount() > 1) return
            if (Math.abs(e.getX() - px0) > 3 || Math.abs(e.getY() - py0) > 3) return                 // a drag is drawing, not a click
            int tool = ij.gui.Toolbar.getToolId()
            if (tool == ij.gui.Toolbar.POLYGON || tool == ij.gui.Toolbar.POLYLINE || tool == ij.gui.Toolbar.ANGLE || tool == ij.gui.Toolbar.MAGNIFIER || tool == ij.gui.Toolbar.HAND) return
            int x = canvas.offScreenX(e.getX()), y = canvas.offScreenY(e.getY())
            javax.swing.SwingUtilities.invokeLater { try { clickCellAt(x, y) } catch (Throwable t) { status.setText("" + t) } }
        }
    }
    def cellMotion = new java.awt.event.MouseMotionAdapter() {
        void mouseMoved(java.awt.event.MouseEvent e) {
            long now = System.currentTimeMillis(); if (now - lastMove < 120) return; lastMove = now
            int x = canvas.offScreenX(e.getX()), y = canvas.offScreenY(e.getY())
            javax.swing.SwingUtilities.invokeLater { try { hoverAt(x, y) } catch (Throwable t) { } }
        }
    }
    canvas.addMouseListener(cellMouse); canvas.addMouseMotionListener(cellMotion)
    frame.addWindowListener(new java.awt.event.WindowAdapter() {
        void windowClosed(java.awt.event.WindowEvent e) { canvas.removeMouseListener(cellMouse); canvas.removeMouseMotionListener(cellMotion) }
    })
}
showCellposeState(null)
setupBtn.addActionListener({ e -> setupBtn.setEnabled(false); Thread.start { if (runSetup(true) == "cancelled") javax.swing.SwingUtilities.invokeLater { showCellposeState(null) } } } as java.awt.event.ActionListener)
def timer = new javax.swing.Timer(600, { e -> refresh(false) } as java.awt.event.ActionListener)
closeBtn.addActionListener({ e -> timer.stop(); detectionProgress.close(); frame.dispose() } as java.awt.event.ActionListener)
// Back to the first screen. The plugin starts again from its menu command; the outlines can be kept or cleared first,
// because the next recording uses whatever outlines the ROI Manager holds.
Closure relaunch = { ->
    String command = 'Calcium Cells Step By Step'
    if (ij.Menus.getCommands()?.containsKey(command)) Thread.start('Calcium plugin restart') { IJ.run(command) }
    else IJ.showMessage('Calcium imaging', 'To choose another recording, run the plugin again:\nPlugins > Calcium Imaging > Calcium Cells Step By Step')
}
def changeRecording = { int choice ->   // 0 clear outlines, 1 keep them
    if (choice == 0) rm.reset()
    timer.stop(); detectionProgress.close(); frame.dispose()
    steps.colourPic.changes = false; steps.colourPic.close()
    if (!testGui) relaunch()
}
changeBtn.addActionListener({ e ->
    int count = rm.getCount()
    if (count == 0) { changeRecording(1); return }
    String unmeasured = analysisMeasurement == null ? ' They have not been measured yet.' : ''
    def options = ['Clear outlines and go back', 'Keep outlines and go back', 'Cancel'] as Object[]
    int answer = javax.swing.JOptionPane.showOptionDialog(frame, '<html><body style="width:320px">The ROI Manager holds ' + count + ' outline' + (count == 1 ? '' : 's') + '.' + unmeasured +
        ' The next recording starts with whatever outlines it holds, so clear them unless they belong to the next recording too.</body></html>',
        'Change recording', javax.swing.JOptionPane.DEFAULT_OPTION, javax.swing.JOptionPane.QUESTION_MESSAGE, null, options, options[0])
    if (answer == 0 || answer == 1) changeRecording(answer)
} as java.awt.event.ActionListener)
openBtn.addActionListener({ e -> try { java.awt.Desktop.getDesktop().open(outDir) } catch (Throwable t) { IJ.log("Results folder: " + outDir.getPath()) } } as java.awt.event.ActionListener)
frame.addWindowListener(new java.awt.event.WindowAdapter() { void windowClosing(java.awt.event.WindowEvent e) { timer.stop(); detectionProgress.close() } })
// Sensitivity check after blind selection: response-assisted detection once, compared with the blind suggestions. Returns the
// comparison and its report, or a message saying why it did not run. Its outlines are never added to the cells or measured.
def runSensitivity = { Roi[] cells, String prefix, Closure tell ->
    if (!steps.found) return [text: 'Sensitivity check not run: no response was found, so there is no response-assisted picture to compare with.']
    def assisted = [], blindSuggestions = []
    javax.swing.SwingUtilities.invokeAndWait { blindSuggestions.addAll(proposals.findAll { it.blind }) }
    def task = new CalciumOutlineTask(progress: { p, s, d -> tell('Sensitivity check: ' + s.toLowerCase() + '…') })
    String how
    try { how = suggest(task, [picture: steps.responseOutlinePicture(), dir: new File(outDir, 'automatic_outlines'), into: assisted]) }
    catch (Throwable t) { return [text: 'Sensitivity check failed: ' + t] }
    if (!task.succeeded) return [text: 'Sensitivity check not completed: ' + how]
    def m = CalciumSensitivity.compare(assisted, blindSuggestions, cells, edgeBox.isSelected())
    CalciumSensitivity.save(m, prefix)
    String text = CalciumSensitivity.report(m, how)
    steps.say(text + ' Files: ' + new File(prefix + 'sensitivity_check.csv').getName() + ' (one row per response-assisted outline) and ' + new File(prefix + 'sensitivity_check_missed_RoiSet.zip').getName() + ' (the missed outlines that could be cells, to open in the ROI Manager).')
    new File(prefix + 'what_was_done.txt').setText(steps.log.join('\n') + '\n', 'UTF-8')
    [m: m, text: text]
}
go.addActionListener({ e ->
    def sp = splitCells(); Roi[] cells = (Roi[]) sp.cells; Roi[] edgeCells = (Roi[]) sp.edge
    final String measuredSelectionBasis=basisFor(sp.all as Roi[])
    final boolean measuredExcludeEdges=edgeBox.isSelected()
    final boolean checkSensitivity = steps.blind && sensBox.isSelected()
    steps.alsoAvoid = edgeCells
    steps.edgeNote = edgeBox.isSelected() ? "Cells at the picture's edge: outlines within " + CalciumSteps.EDGE_MARGIN + " pixels of the border are left out, because their background can only be taken from the inner, better-lit side (" + edgeCells.length + " outlines in the ROI Manager left out; they are not in the tables or the saved outlines)." : "Cells at the picture's edge were NOT left out (box unticked); the summary table marks them."
    if (cells.length == 0) { result.setText("<html><body style='width:300px'>There are no cells yet. Choose cells in step 2 first.</body></html>"); return }
    go.setEnabled(false); go.setText("Working..."); analyzeBtn.setEnabled(false); analysisMeasurement = null
    Thread.start {
        try {
            outDir.mkdirs()
            def names = cells.collect { CalciumMeasurement.idFor(it) } as String[]
            ImagePlus corrected = null
            if (makeMovie.isSelected()) {
                result.setText("<html><body style='width:300px'>Making the background-corrected movie...</body></html>")
                corrected = steps.correctedMovie(cells)
                if (saveMovie.isSelected()) IJ.saveAsTiff(corrected, outDir.getPath() + File.separator + base + "_background_corrected_movie.tif")
            }
            result.setText("<html><body style='width:300px'>Measuring " + cells.length + " cells on " + steps.n + " frames...</body></html>")
            steps.log.removeAll{it.startsWith('Cell selection basis:') || it.startsWith('Sensitivity check') || it.startsWith('Chosen blind or later:')};steps.say('Cell selection basis: '+measuredSelectionBasis+'. Accepted ROI registry includes premeasurement exclusions.')
            Map chosenBlind = null
            if (steps.blindLast > 0 && (steps.blind || blindCellKeys != null)) {
                chosenBlind = [:]
                cells.eachWithIndex { Roi c, int i -> chosenBlind[names[i]] = blindCellKeys != null ? blindCellKeys.contains(cellKey(c)) : !notBlindKeys.contains(cellKey(c)) }
                int nBlind = chosenBlind.values().count { it } as int
                CalciumCsv.write(new File(outDir, base + '_cells_chosen_blind_or_after_response.csv'), [['cell', 'chosen']] + names.collect { [it, chosenBlind[it] ? 'blind (frames ' + steps.blindFirst + '-' + steps.blindLast + ')' : 'after the response was shown'] })
                steps.say('Chosen blind or later: ' + nBlind + ' of the ' + cells.length + ' cells were chosen blind on frames ' + steps.blindFirst + '-' + steps.blindLast + ' and are unchanged since; ' + (cells.length - nBlind) + ' were added or reshaped after the response was shown. Step 4 gives % responding for each group; only the blind group is free of selection on the response. List: ' + base + '_cells_chosen_blind_or_after_response.csv')
            }
            def summary = steps.measureAndSave(cells, names, outDir.getPath() + File.separator + base + "_")
            if (corrected != null && !testGui) { corrected.show(); corrected.setSlice(Math.min(steps.n, steps.found ? steps.p0 + 5 : 1)); corrected.setDisplayRange(0, CalciumSteps.SHOW_MAX); corrected.updateAndDraw() }
            IJ.log(summary); IJ.log("Saved in: " + outDir.getPath()); IJ.log("Step 4 uses full-precision unlifted fluorescence: " + base + "_gray_values_background_corrected_native.csv")
            File registryFile=new File(outDir,base+'_roi_registry.csv')
            CalciumCsv.write(registryFile,CalciumMeasurement.registry(sp.all as Roi[],cells,steps.lastBackgroundPixels,steps,measuredSelectionBasis,
                chosenBlind==null?null:chosenBlind.collectEntries{id,blindly->[id,blindly?steps.blindBasis():'Added after the response or all frames were shown (blind frames '+steps.blindFirst+'-'+steps.blindLast+')']}))
            File acceptedRoisFile=new File(outDir,base+'_accepted_rois.zip');CalciumMeasurement.writeRois(acceptedRoisFile,sp.all as Roi[])
            double interval = imp.getCalibration().frameInterval
            String unit = imp.getCalibration().getTimeUnit().toLowerCase(Locale.ROOT)
            if (unit in ['ms','msec','millisecond']) interval /= 1000d
            else if (unit in ['min','minute']) interval *= 60d
            else if (!(unit in ['sec','s','second'])) interval = Double.NaN
            Map handed = CalciumAnalysisSettings.framesFromFiji(steps)
            final Map fromFiji = handed
            String stimulusFrom = 'fiji'
            int usual = CalciumAnalysisSettings.lastStimulusFrame()
            if (usual >= 1 && usual < steps.n - 1) {
                handed = [first: Math.max(1, usual - 19), last: usual, chosen: handed.chosen, note: 'The stimulus is set after frame ' + usual + ', as in your last analysis. Check it, then run.']
                stimulusFrom = 'last-recording'
            }
            final Map measured = [csv:new File(outDir,base + "_gray_values_background_corrected_native.csv"),folder:outDir,
                measurementId:UUID.randomUUID().toString(),recording:base,rois:cells.collect{(Roi)it.clone()} as Roi[],names:names,frames:steps.n,picture:steps.colourPic,backgroundPicture:steps.backgroundPic,motion:motionForPage(handed.chosen),movie:imp,
                hz:interval > 0 ? 1d/interval : Double.NaN,baselineStart:handed.first,baselineEnd:handed.last,framesNote:handed.note,fijiFrames:handed.chosen,stimulusFrom:stimulusFrom,fijiBaselineStart:fromFiji.first,fijiBaselineEnd:fromFiji.last,chosenBlind:chosenBlind,key:measuredKey(sp.all as Roi[])+'#edge='+measuredExcludeEdges+'#selection='+measuredSelectionBasis,
                registry:registryFile,registryHash:CalciumInstaller.digest(registryFile.bytes),selectionBasis:measuredSelectionBasis,
                acceptedRois:acceptedRoisFile,acceptedRoisHash:CalciumInstaller.digest(acceptedRoisFile.bytes)]
            measured.inputHash=CalciumInstaller.digest(measured.csv.bytes)
            String sensitivityHtml = ""
            if (checkSensitivity) {
                def check = runSensitivity(cells, outDir.getPath() + File.separator + base + "_", { String s -> javax.swing.SwingUtilities.invokeLater { result.setText("<html><body style='width:300px'>" + CalciumOutlineProgress.escape(s) + "</body></html>") } })
                def m = check.m
                sensitivityHtml = m != null ? "Sensitivity check: the blind detection also found <b>" + m.found + " of " + m.cellLike + "</b> response-assisted outlines that could be cells (" + CalciumSensitivity.pct((int) m.found, (int) m.cellLike) + "); " + (m.cellLike - m.found) + " missed. Details in what_was_done.txt.<br>" :
                    "<font color='#b05000'>" + CalciumOutlineProgress.escape((String) check.text) + "</font><br>"
                IJ.log((String) check.text)
            }
            def html = "<b>Done.</b> " + cells.length + " cells measured on " + steps.n + " frames.<br>" + sensitivityHtml +
                "Legacy 5-SD screening rise: <b>" + steps.lastClear + " of " + cells.length + "</b>. Protocol responder counts are in Step 4.<br>" +
                (edgeCells.length > 0 ? edgeCells.length + " outlines at the picture's edge were left out (dashed).<br>" : "") +
                (steps.lastDim > 0 ? "<font color='#b05000'><b>" + steps.lastDim + " cells are not brighter than their background at rest.</b> Review their quality status in Step 4.</font><br>" : "") +
                "<b>Continue to Step 4 to analyze responses inside Fiji.</b> It uses unlifted fluorescence and reports excluded and unscorable cells explicitly.<br>Saved in: " + savedIn
            javax.swing.SwingUtilities.invokeLater { analysisMeasurement = measured; analyzeBtn.setEnabled(true); result.setText("<html><body style='width:300px'>" + html + "</body></html>"); go.setText("Measure again"); go.setEnabled(true); s3.setSummary(cells.length + " cells measured"); s4.setExpanded(true); refresh(true); openBtn.setVisible(true); frame.refreshContents() }
        } catch (Throwable t) { go.setEnabled(true); go.setText("Measure and save"); result.setText("<html><body style='width:300px'>Could not measure: " + t + "</body></html>"); IJ.log("Calcium imaging plugin: " + t) }
    }
} as java.awt.event.ActionListener)
if (testGui) { binding.setVariable("hooks", [nextLine: nextLine, addRespondersBtn: addRespondersBtn, acceptAllBtn: acceptAllBtn, autoSuggestAsked: { -> autoSuggestAsked }, legendClicked: { String k, boolean o -> legendClicked.call(k, o) }, hiddenKinds: hiddenKinds, cellRois: cellRois, sensBox: sensBox, basisFor: basisFor, responseSeen: { -> responseSeen }, reveal: { -> reveal() }, runSensitivity: runSensitivity, trace: trace, respondsLegend: respondsLegend, viewAllBtn: viewAllBtn, viewBaselineBtn: viewBaselineBtn, viewResponseBtn: viewResponseBtn, selectionNote: selectionNote, selectionBasis: { -> selectionBasis }, changeBtn: changeBtn, changeRecording: changeRecording, continueBtn: continueBtn, backToCellsBtn: backToCellsBtn, backToMeasureBtn: backToMeasureBtn, cards: [s1, s2, s3, s4], recordingNote: recordingNote, detectionProgress: detectionProgress, suggestBtn: suggestBtn, edgeBox: edgeBox, addBtn: addBtn, setForceFresh: { boolean v -> forceFresh = v }, setupCellpose: setupCellpose, findCellposePython: findCellposePython, setupBtn: setupBtn, setupRow: setupRow, cpLabel: cpLabel, showCellposeState: showCellposeState, setupText: SETUP_TEXT, saveMovie: saveMovie, makeMovie: makeMovie, rm: rm, refresh: refresh, go: go, status: status, steps: steps, outDir: outDir, showBg: showBg, frame: frame, suggest: suggest, acceptAt: acceptAt, acceptAll: acceptAll, hoverAt: hoverAt, proposals: proposals, result: result, info: info, removeSelected: removeSelected, undoLast: undoLast, root: root]); refresh(true); return }
frame.setVisible(true)
timer.start()
refresh(true)

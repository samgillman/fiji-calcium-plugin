// Ca²⁺ Imaging (Wang Lab)                                                                  version 2026.11.22.2
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
// Install: run Install_Calcium_Plugin.groovy from the release page in Fiji, or see README.md for copying by hand.
// Run:     Plugins > Ca²⁺ Imaging (Wang Lab)
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
    // the wider "likely" rule of 2026.11.20.20 (see grade): at least this many pixels, this much brighter than the surroundings, and less than this share of very bright pixels
    static int LIKELY_MIN_AREA = 500
    static double LIKELY_CONTRAST = 0.3, LIKELY_BRIGHT_SHARE = 0.6
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
    // With a graphics chip a third picture, 0 to BLIND_CAP_MID, is run as well and its outlines are added by the same rule against
    // the outlines kept so far. Measured 20 Nov 2026 on four audited recordings: 68.3 % of the audited cells get an outline against
    // 65.5 % with two pictures, for 6 % more outlines. Without a graphics chip each picture takes minutes, so it is left out there.
    static double BLIND_CAP_MID = 4
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
        say(CalciumRelease.NAME + ", version 2026.11.22.2")
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
        // Wider "likely" (2026.11.20.20). A cell-sized outline that stands out from its surroundings and is not bright on average is
        // likely a cell even with bright parts inside. Measured 20 Nov 2026 on the blind suggestions of four audited recordings
        // (871 outlines, 427 on audited cells): the rule above alone marks 312 outlines likely, 73 % of them on audited cells, which is
        // 53 % of the outlines that are cells; with this rule 472, 72 % and 80 %. It re-grades outlines that were unsure (a bright dot
        // inside, or contrast between 0.2 and 0.3 is NOT enough: contrast must exceed LIKELY_CONTRAST) and outlines that were "very
        // bright" only for their share of bright pixels. An outline brighter than BRIGHT_ABOVE on average stays as it was.
        if (g >= 2 && mean <= BRIGHT_ABOVE && a >= LIKELY_MIN_AREA && contrast > LIKELY_CONTRAST && brightShare < LIKELY_BRIGHT_SHARE && !edge) g = 1
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
        new File(prefix + "what_was_done.txt").setText(log.join("\n") + "\n", "UTF-8")
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
    static final String VERSION = "2026.11.22.2"
    static final String REPOSITORY = "https://github.com/samgillman/fiji-calcium-plugin"
    static final String API = "https://api.github.com/repos/samgillman/fiji-calcium-plugin/releases/latest"
    // Download names. They stay as they are: installed copies of earlier versions look for exactly these.
    static final String ASSET = "Calcium_Imaging_Fiji_Plugin.zip"
    static final String SCRIPT = "Calcium_Cells_Step_By_Step.groovy"
    static final String CHECKSUMS = "SHA256SUMS.txt"
    /** The name in Fiji's Plugins menu and on the windows. Fiji takes a menu entry from the file name, with _ shown as a space. */
    static final String NAME = "Ca\u00b2\u207a Imaging (Wang Lab)"
    static final String PLACE = "scripts/Plugins/Ca\u00b2\u207a_Imaging_(Wang_Lab).groovy"
    /** The ZIP's place, since ZIP names stay plain letters; elsewhere moved to PLACE when the plugin starts. Kept where the
     *  system cannot store ² or ⁺ in a file name: the menu entry then reads Ca2+ Imaging (Wang Lab). */
    static final String PLAIN_PLACE = "scripts/Plugins/Ca2+_Imaging_(Wang_Lab).groovy"
    /** Versions up to 2026.11.22.1: Plugins > Calcium Imaging > Calcium Cells Step By Step. */
    static final String OLD_PLACE = "scripts/Plugins/Calcium Imaging/" + SCRIPT
    static final List<String> PLACES = [PLACE, PLAIN_PLACE, OLD_PLACE].asImmutable()
    /** The first version that knows these places. An older copy is never moved, only put aside. */
    static final String FIRST_PLACED = "2026.11.22.2"
    static String command(String place) { new File(place).name.replaceFirst(/\.groovy$/, "").replace("_", " ") }

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
        File app = home()
        File target = inside(app, place(app))
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
        java.nio.file.Path staged = null
        try {
            return (Map) locked(app) {
                List<File> copies = elsewhere(app)
                byte[] previous = target.isFile() ? java.nio.file.Files.readAllBytes(target.toPath()) : null
                String diskVersion = previous == null ? null : scriptVersion(previous)
                if (([diskVersion] + copies.collect { versionOf(it) }).any { it && CalciumRelease.newer(it, version) })
                    throw new IOException("A newer plugin is already installed; check for updates again")
                boolean already = previous != null && digest(previous) == expected
                File backup = null
                if (!already) {
                    java.nio.file.Files.createDirectories(target.parentFile.toPath())
                    staged = java.nio.file.Files.createTempFile(target.parentFile.toPath(), ".calcium-update-", ".tmp")
                    java.nio.channels.FileChannel.open(staged, java.nio.file.StandardOpenOption.WRITE).withCloseable { output ->
                        def buffer = java.nio.ByteBuffer.wrap(bytes)
                        while (buffer.hasRemaining()) output.write(buffer)
                        output.force(true)
                    }
                    if (previous != null) {
                        progress.call(92, "Backing up the installed plugin...")
                        backup = backupFor(app, previous)
                        java.nio.file.Files.write(backup.toPath(), previous, java.nio.file.StandardOpenOption.CREATE_NEW, java.nio.file.StandardOpenOption.WRITE)
                    }
                    // Do not overwrite an external edit made while this update was being staged.
                    if (previous == null ? target.exists() : (!target.isFile() || digest(java.nio.file.Files.readAllBytes(target.toPath())) != digest(previous)))
                        throw new IOException("The installed file changed during the update; please try again")
                    progress.call(98, "Installing the verified plugin...")
                    replaceFile.call(staged, target.toPath())
                }
                // A copy left at an earlier place would be a second menu entry.
                Map tidy = putAside(app, copies)
                if (backup == null && copies) backup = tidy.aside[copies[0]]
                progress.call(100, "Installed " + version + ". Save your work and restart Fiji.")
                [version:version, target:target, backup:backup, already:already && !tidy.aside, aside:tidy.aside, left:tidy.left]
            }
        } catch (java.nio.file.AccessDeniedException denied) {
            throw new IOException("Fiji cannot write to its installation folder. Use a Fiji folder your account can modify, then try again.", denied)
        } catch (java.nio.file.AtomicMoveNotSupportedException unsupported) {
            throw new IOException("This folder does not support safe replacement. The existing plugin was kept.", unsupported)
        } finally {
            if (staged != null) java.nio.file.Files.deleteIfExists(staged)
        }
    }

    /** Versions up to 2026.11.22.1 were installed under Plugins > Calcium Imaging, and their updater still puts a new
     *  version there; the ZIP puts it at PLAIN_PLACE. A copy at another place is moved here when it is version FIRST_PLACED or later and newer than the
     *  copy here, or there is none. Copies that are not newer than the one here go to calcium-plugin-backups. An older
     *  version with nothing here is left where it is. Returns [moved: the copy moved here or null, backup, aside, left, target]. */
    Map settle() {
        File app = home()
        File target = inside(app, place(app))
        Map none = [moved:null, backup:null, aside:[:], left:[], target:target]
        if (!elsewhere(app)) return none
        try {
            return (Map) locked(app) {
                List<File> copies = elsewhere(app)
                String here = target.isFile() ? versionOf(target) : null
                if (target.isFile() && here == null) return none   // a copy here without a version: change nothing
                File best = null; String now = here
                copies.each { File copy ->
                    String v = versionOf(copy)
                    if (v && !CalciumRelease.newer(CalciumRelease.FIRST_PLACED, v) && (now == null || CalciumRelease.newer(v, now))) { best = copy; now = v }
                }
                File backup = null
                if (best != null) {
                    if (target.isFile()) {
                        byte[] there = java.nio.file.Files.readAllBytes(target.toPath())
                        backup = backupFor(app, there)
                        java.nio.file.Files.write(backup.toPath(), there, java.nio.file.StandardOpenOption.CREATE_NEW, java.nio.file.StandardOpenOption.WRITE)
                    }
                    java.nio.file.Files.createDirectories(target.parentFile.toPath())
                    replaceFile.call(best.toPath(), target.toPath())
                    removeIfEmpty(app, best.parentFile)
                }
                Map tidy = now == null ? [aside:[:], left:[]] :
                    putAside(app, copies.findAll { File copy -> String v = versionOf(copy); copy != best && v && !CalciumRelease.newer(v, now) })
                [moved:best, backup:backup, aside:tidy.aside, left:tidy.left, target:target]
            }
        } catch (java.nio.file.AtomicMoveNotSupportedException unsupported) {
            throw new IOException("This folder does not support safe replacement. The plugin was left where it was.", unsupported)
        }
    }

    /** Fiji's folder, checked. */
    File home() {
        if (applicationDirectory == null || !applicationDirectory.isDirectory()) throw new IOException("Could not locate this Fiji installation")
        applicationDirectory.canonicalFile
    }
    /** Where the plugin belongs in Fiji's folder: PLACE, or PLAIN_PLACE where the system cannot store ² or ⁺ in a file name. */
    static File place(File app) { [CalciumRelease.PLACE, CalciumRelease.PLAIN_PLACE].collect { new File(app, it) }.find { nameable(it) } }
    /** Copies at the plugin's other places that may be moved: inside Fiji's folder and not symbolic links. */
    static List<File> elsewhere(File app) {
        File here = place(app)
        CalciumRelease.PLACES.collect { new File(app, it) }.findAll { File f -> f != here && nameable(f) && f.isFile() && within(app, f) }
    }
    static boolean nameable(File file) {
        try { file.toPath(); return true } catch (java.nio.file.InvalidPathException ignored) { return false }
    }
    /** Inside Fiji's folder with no link on the way: what exists of the path is really where its name says. On Windows
     *  a canonical path does not resolve links, and Files.isSymbolicLink misses junctions; a real path resolves both. */
    static boolean within(File app, File file) {
        def base = app.toPath().toAbsolutePath().normalize(), path = file.toPath().toAbsolutePath().normalize()
        if (path == base || !path.startsWith(base)) return false
        def part = path
        while (!java.nio.file.Files.exists(part, java.nio.file.LinkOption.NOFOLLOW_LINKS)) part = part.parent
        try { return part.toRealPath() == base.toRealPath().resolve(base.relativize(part)) } catch (IOException ignored) { return false }
    }
    static File inside(File app, File file) {
        if (!within(app, file)) throw new IOException("The plugin installation path points outside this Fiji folder or uses a symbolic link")
        file
    }
    static String versionOf(File file) {
        try { return scriptVersion(java.nio.file.Files.readAllBytes(file.toPath())) } catch (IOException ignored) { return null }
    }
    /** One installation at a time, also across Fiji windows. */
    static Object locked(File app, Closure work) {
        def channel = java.nio.channels.FileChannel.open(new File(app, ".calcium-plugin-update.lock").toPath(), java.nio.file.StandardOpenOption.CREATE, java.nio.file.StandardOpenOption.WRITE)
        try {
            java.nio.channels.FileLock lock = null
            try { lock = channel.tryLock() } catch (java.nio.channels.OverlappingFileLockException ignored) { }
            if (lock == null) throw new IOException("Another plugin installation is in progress; wait for it to finish")
            try { return work.call() } finally { lock.release() }
        } finally { channel.close() }
    }
    /** A new file name in calcium-plugin-backups inside Fiji's folder. */
    static File backupFor(File app, byte[] bytes) {
        File backups = new File(app, "calcium-plugin-backups")
        java.nio.file.Files.createDirectories(backups.toPath())
        new File(backups, "Calcium-" + (scriptVersion(bytes) ?: "previous") + "-" + UUID.randomUUID() + ".groovy.bak")
    }
    /** Moves copies into calcium-plugin-backups, so that the menu shows the plugin once. A folder left empty is removed. */
    static Map putAside(File app, List<File> copies) {
        Map aside = [:]; List<File> left = []
        copies.each { File copy ->
            try {
                File backup = backupFor(app, java.nio.file.Files.readAllBytes(copy.toPath()))
                java.nio.file.Files.move(copy.toPath(), backup.toPath())
                aside[copy] = backup
                removeIfEmpty(app, copy.parentFile)
            } catch (IOException ignored) { left << copy }
        }
        [aside:aside, left:left]
    }
    /** An earlier version's folder, such as Calcium Imaging, once nothing but Finder's .DS_Store is left in it. */
    static void removeIfEmpty(File app, File folder) {
        if (folder == new File(app, "scripts/Plugins") || !(folder.list() as List)?.every { it == ".DS_Store" }) return
        new File(folder, ".DS_Store").delete(); folder.delete()
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
    static void choose(javax.swing.AbstractButton b, boolean on) {
        if (!(b.getUI() instanceof CalciumButtonUI)) return
        ((CalciumButtonUI) b.getUI()).chosen = on; if (on && b.getParent() instanceof CalciumSegments) ((CalciumSegments) b.getParent()).raise(b); b.repaint()
    }
    /** Controls side by side, the first at the left edge so that it lines up with the text above and below it. */
    static javax.swing.JPanel row(java.awt.Component... parts) {
        def p = new javax.swing.JPanel(new CalciumRowLayout()); p.setOpaque(false); p.setAlignmentX(0f)
        parts.each { p.add(it) }
        p
    }
    /** A long name without spaces, cut in the middle so that it fits a width: its start and its end stay readable. */
    static String fitMiddle(String text, java.awt.FontMetrics fm, int width) {
        if (text == null || width <= 0) return ""
        if (fm.stringWidth(text) <= width) return text
        for (int most = text.length() - 1; most >= 3; most--) { String cut = middle(text, most); if (fm.stringWidth(cut) <= width) return cut }
        return "…"
    }
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
    String segment = null                     // 'first', 'middle' or 'last' part of a CalciumSegments: square where it meets its neighbour
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
            if (segment) {
                // one part of a joined switch: round only at the ends of the switch, and the edge between two parts stays gray
                // unless one of them is the chosen one (which is drawn last, see CalciumSegments)
                def shape = new java.awt.geom.Area(new java.awt.geom.RoundRectangle2D.Float(0f, 0f, (float) (w - 1), (float) (h - 1), 10f, 10f))
                int half = w.intdiv(2)
                if (segment != 'first') shape.add(new java.awt.geom.Area(new java.awt.geom.Rectangle2D.Float(0f, 0f, (float) half, (float) (h - 1))))
                if (segment != 'last') shape.add(new java.awt.geom.Area(new java.awt.geom.Rectangle2D.Float((float) half, 0f, (float) (w - 1 - half), (float) (h - 1))))
                if (!chosen && m.isEnabled()) line = new Color(0xC5CED8)
                g.setColor(fill); g.fill(shape); g.setColor(line); g.draw(shape)
            } else {
                g.setColor(fill); g.fillRoundRect(0, 0, w - 1, h - 1, 10, 10)
                g.setColor(line); g.drawRoundRect(0, 0, w - 1, h - 1, 10, 10)
            }
            if (b.isFocusOwner() && m.isEnabled() && !segment) { g.setColor(new Color(0, 114, 178, 110)); g.setStroke(new java.awt.BasicStroke(2f)); g.drawRoundRect(1, 1, w - 3, h - 3, 9, 9) }
            g.setFont(b.getFont()); def fm = g.getFontMetrics(); String text = b.getText() ?: ""
            g.setColor(ink); g.drawString(text, (int) ((w - fm.stringWidth(text)) / 2), (int) ((h - fm.getHeight()) / 2) + fm.getAscent())
        } finally { g.dispose() }
    }
}

/** Text that wraps at the width it is given, so it fills its card and follows the window when that is resized. (A label with a
 *  width written into its text wraps there whatever the card's width is.) say() takes the words without html tags. */
class CalciumText extends javax.swing.JLabel {
    static final int NARROW = 200, EXPECTED = 440   // the width asked for (text never makes a card wider than its controls do), and the width a card usually has
    int laidOut = -1
    int lines = 0   // room kept for this many lines even when the words take fewer: what is under the text then stays where it is
    CalciumText(String words = "&nbsp;", Color colour = CalciumTheme.INK) {
        setFont(CalciumTheme.font()); setForeground(colour); setAlignmentX(0f); setVerticalAlignment(javax.swing.SwingConstants.TOP); say(words)
        // the height depends on the width, which is known only once the text has been laid out: lay out once more when it changes
        addComponentListener(new java.awt.event.ComponentAdapter() { void componentResized(java.awt.event.ComponentEvent e) { if (getWidth() != laidOut) { laidOut = getWidth(); revalidate() } } })
    }
    void say(String words) { setText("<html>" + words + "</html>") }
    /** The height the text takes when wrapped at a width. */
    private int tall(int width) {
        def view = getClientProperty(javax.swing.plaf.basic.BasicHTML.propertyKey)
        if (!(view instanceof javax.swing.text.View)) return (int) super.getPreferredSize().height
        def insets = getInsets()
        ((javax.swing.text.View) view).setSize((float) Math.max(1, width - insets.left - insets.right), 0f)
        return (int) Math.ceil(((javax.swing.text.View) view).getPreferredSpan(javax.swing.text.View.Y_AXIS)) + insets.top + insets.bottom
    }
    private int room() {
        if (getWidth() > 0) return getWidth()
        def p = getParent(); if (p != null && p.getWidth() > 0) { def i = p.getInsets(); return p.getWidth() - i.left - i.right }
        return EXPECTED
    }
    private boolean wordless() { String t = getText(); t == null || t.replaceAll("<[^>]*>", "").trim().isEmpty() }
    java.awt.Dimension getPreferredSize() {
        if (isPreferredSizeSet()) return super.getPreferredSize()
        def insets = getInsets(); int kept = lines > 0 ? lines * getFontMetrics(getFont()).getHeight() + insets.top + insets.bottom : 0
        new java.awt.Dimension(NARROW, Math.max(kept, wordless() ? 0 : tall(room())))
    }
    java.awt.Dimension getMinimumSize() { getPreferredSize() }
    java.awt.Dimension getMaximumSize() { new java.awt.Dimension(Integer.MAX_VALUE, (int) getPreferredSize().height) }
    protected void paintComponent(java.awt.Graphics g) { if (getWidth() > 0) tall(getWidth()); super.paintComponent(g) }
}

/** One line that may be long and has no spaces, such as a file name: shown whole when it fits, else with its middle left out. */
class CalciumName extends javax.swing.JLabel {
    String whole = ""
    CalciumName(String text) { super(text ?: ""); whole = text ?: ""; setFont(CalciumTheme.font(java.awt.Font.BOLD)); setForeground(CalciumTheme.INK); setAlignmentX(0f); setToolTipText(whole) }
    void setWhole(String text) { whole = text ?: ""; setText(whole); setToolTipText(whole); repaint() }
    java.awt.Dimension getPreferredSize() { def fm = getFontMetrics(getFont()); new java.awt.Dimension(Math.min(CalciumText.NARROW, fm.stringWidth(whole)), fm.getHeight()) }
    java.awt.Dimension getMinimumSize() { getPreferredSize() }
    java.awt.Dimension getMaximumSize() { new java.awt.Dimension(Integer.MAX_VALUE, (int) getPreferredSize().height) }
    protected void paintComponent(java.awt.Graphics g0) {
        def g = (java.awt.Graphics2D) g0.create()
        try { CalciumTheme.smooth(g); g.setFont(getFont()); def fm = g.getFontMetrics(); g.setColor(getForeground()); g.drawString(CalciumTheme.fitMiddle(whole, fm, getWidth()), 0, fm.getAscent()) } finally { g.dispose() }
    }
}

/** Lays controls side by side, 6 pixels apart, starting at the left edge; hidden ones take no room. */
class CalciumRowLayout implements java.awt.LayoutManager {
    static final int GAP = 6, PAD = 3
    void addLayoutComponent(String name, java.awt.Component c) { }
    void removeLayoutComponent(java.awt.Component c) { }
    java.awt.Dimension preferredLayoutSize(java.awt.Container p) {
        int w = 0, h = 0, n = 0
        for (java.awt.Component c : p.getComponents()) if (c.isVisible()) { def d = c.getPreferredSize(); w += (int) d.width + (n++ > 0 ? GAP : 0); h = Math.max(h, (int) d.height) }
        def i = p.getInsets(); new java.awt.Dimension(w + i.left + i.right, h + 2 * PAD + i.top + i.bottom)
    }
    java.awt.Dimension minimumLayoutSize(java.awt.Container p) { preferredLayoutSize(p) }
    void layoutContainer(java.awt.Container p) {
        def i = p.getInsets(); int x = i.left, h = p.getHeight() - i.top - i.bottom
        for (java.awt.Component c : p.getComponents()) if (c.isVisible()) { def d = c.getPreferredSize(); c.setBounds(x, i.top + (int) ((h - (int) d.height) / 2), (int) d.width, (int) d.height); x += (int) d.width + GAP }
    }
}

/** Choices of one thing joined into one control, like the plot's "Average | One cell + average" on the Step 4 page. Each part
 *  keeps its place whichever is chosen; the chosen one is tinted (CalciumTheme.choose). */
class CalciumSegments extends javax.swing.JPanel {
    List<javax.swing.AbstractButton> parts = []
    CalciumSegments(javax.swing.AbstractButton... buttons) {
        setLayout(null); setOpaque(false); setAlignmentX(0f)
        buttons.eachWithIndex { b, i -> parts << b; if (b.getUI() instanceof CalciumButtonUI) ((CalciumButtonUI) b.getUI()).segment = i == 0 ? 'first' : i == buttons.length - 1 ? 'last' : 'middle'; add(b) }
    }
    boolean isOptimizedDrawingEnabled() { false }   // neighbours share one edge
    java.awt.Dimension getPreferredSize() { int w = 0, h = 0; parts.each { def d = it.getPreferredSize(); w += (int) d.width; h = Math.max(h, (int) d.height) }; new java.awt.Dimension(w - (parts.size() - 1), h) }
    java.awt.Dimension getMinimumSize() { getPreferredSize() }
    java.awt.Dimension getMaximumSize() { getPreferredSize() }
    void doLayout() { int x = 0, h = getHeight(); parts.each { def d = it.getPreferredSize(); it.setBounds(x, 0, (int) d.width, h); x += (int) d.width - 1 } }
    /** The chosen part is drawn last, so where two parts meet the edge seen is the chosen one's. */
    void raise(javax.swing.AbstractButton part) { if (part.getParent() == this) setComponentZOrder(part, 0); repaint() }
}

/** The tick box of this window: the same for a tick with words and for a Show switch, and it starts at the very left of its
 *  control, so it lines up with the text and the buttons above and below. */
class CalciumTick implements javax.swing.Icon {
    boolean on; CalciumTick(boolean on) { this.on = on }
    int getIconWidth() { 15 }
    int getIconHeight() { 15 }
    void paintIcon(java.awt.Component c, java.awt.Graphics g0, int x, int y) { draw((java.awt.Graphics2D) g0, x, y, on, c != null && !c.isEnabled(), false) }
    static void draw(java.awt.Graphics2D g0, int x, int y, boolean on, boolean gray, boolean hover) {
        def g = (java.awt.Graphics2D) g0.create()
        try {
            CalciumTheme.smooth(g); Color blue = gray ? new Color(0xB7D3E3) : CalciumTheme.BLUE
            g.setColor(on ? blue : gray ? CalciumTheme.PAGE : Color.white); g.fillRoundRect(x, y, 14, 14, 5, 5)
            g.setColor(on ? blue : hover ? CalciumTheme.BLUE : new Color(0xC5CED8)); g.drawRoundRect(x, y, 14, 14, 5, 5)
            if (on) { g.setColor(Color.white); g.setStroke(new java.awt.BasicStroke(1.8f, java.awt.BasicStroke.CAP_ROUND, java.awt.BasicStroke.JOIN_ROUND)); g.drawPolyline([x + 3, x + 6, x + 11] as int[], [y + 7, y + 10, y + 4] as int[], 3) }
        } finally { g.dispose() }
    }
    /** Gives a tick with words this box. */
    static void dress(javax.swing.JCheckBox box) {
        box.setIcon(new CalciumTick(false)); box.setSelectedIcon(new CalciumTick(true)); box.setIconTextGap(6); box.setFocusPainted(false)
        box.setBorder(javax.swing.BorderFactory.createEmptyBorder(3, 0, 3, 0)); box.setMargin(new java.awt.Insets(0, 0, 0, 0))
    }
}

/** A switch for one kind of thing on the picture: a tick box, a sample of how it is drawn, its name. A click turns it off or
 *  on; a double-click is passed on too ("only this one"). The owner decides what the clicks do and sets on or off. */
class CalciumSwitch extends javax.swing.JComponent {
    String words; javax.swing.Icon sample; boolean on = true, hover = false
    Closure clicked = { boolean twice -> }
    CalciumSwitch(String words, javax.swing.Icon sample) {
        this.words = words; this.sample = sample
        setFont(CalciumTheme.font(java.awt.Font.PLAIN, -1f)); setAlignmentX(0f); setOpaque(false); setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR))
        addMouseListener(new java.awt.event.MouseAdapter() {
            void mouseClicked(java.awt.event.MouseEvent e) { if (javax.swing.SwingUtilities.isLeftMouseButton(e)) clicked.call(e.getClickCount() >= 2) }
            void mouseEntered(java.awt.event.MouseEvent e) { hover = true; repaint() }
            void mouseExited(java.awt.event.MouseEvent e) { hover = false; repaint() }
        })
    }
    void setOn(boolean value) { if (on != value) { on = value; repaint() } }
    java.awt.Dimension getPreferredSize() { def fm = getFontMetrics(getFont()); new java.awt.Dimension(15 + 6 + (sample != null ? sample.getIconWidth() + 6 : 0) + fm.stringWidth(words) + 2, Math.max(22, fm.getHeight() + 4)) }
    java.awt.Dimension getMinimumSize() { getPreferredSize() }
    protected void paintComponent(java.awt.Graphics g0) {
        def g = (java.awt.Graphics2D) g0.create()
        try {
            CalciumTheme.smooth(g); int h = getHeight(), x = 0
            CalciumTick.draw(g, x, (int) ((h - 15) / 2), on, false, hover)
            x += 15 + 6
            if (sample != null) { sample.paintIcon(this, g, x, (int) ((h - sample.getIconHeight()) / 2)); x += sample.getIconWidth() + 6 }
            g.setFont(getFont()); def fm = g.getFontMetrics(); g.setColor(on ? CalciumTheme.MUTED : CalciumTheme.FAINT)
            g.drawString(words, x, (int) ((h - fm.getHeight()) / 2) + fm.getAscent())
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
    javax.swing.JLabel text = new CalciumText()   // wraps at the card's width
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
        def buttons = CalciumTheme.row()   // its buttons line up with every other row of the card
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
        text.setText("<html><b>" + escape(state.stage) + "</b> &nbsp;<font color='" + muted + "'>" + state.percent + " %</font><br><font color='" + muted + "'>" + escape(elapsed) +
            (context ? " &nbsp;·&nbsp; " + escape(context) : "") + "</font><br>" + escape(message) +
            (quietMessage ? "<br><font color='#865000'>" + escape(quietMessage) + "</font>" : "") + "</html>")
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
            bar.setVisible(false); cancelButton.setVisible(false); roomAsNeeded()
            text.setText("<html>" + escape(message) + "</html>")
        }
    }
    void roomFor(int height) {
        text.setPreferredSize(new java.awt.Dimension(330, height)); text.setMinimumSize(new java.awt.Dimension(330, height)); text.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, height))
        revalidate()
    }
    /** A finished report takes the lines its words need at the card's width. */
    void roomAsNeeded() { text.setPreferredSize(null); text.setMinimumSize(null); text.setMaximumSize(null); revalidate() }
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
                catch (Throwable problem) { message = "Failed: " + problem; IJ.log(CalciumRelease.NAME + ": " + problem) }
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

/** The response as the movie views show it ("Response in magenta", "Change from baseline"): how far each pixel rose above
 *  "no change". One definition for every view that shows a response, made so that a recording in which nothing responds stays
 *  dark and one whose cells rise by a gray level or two shows them. Measured 20 Nov 2026: on P2 Run1 ATP 97 of 177 measured
 *  cells are magenta at the peak and 1 before the stimulus; on two recordings without a response at most 1.8 % of the picture
 *  is ever magenta. (The scale before this, taken from single pixels, started magenta at 2.8 gray levels on that recording,
 *  where the cells rise by 1 to 3.)
 *   - Averaging. The frame is averaged with FRAMES frames before and after it and over a box of (2 x BOX + 1) squared pixels.
 *     Camera noise of about 4 gray levels in one pixel becomes about 0.2.
 *   - "No change" is what the background did. In every block of BLOCK pixels, the middle value of the change over the block's
 *     dimmest BACK_SHARE of pixels (in the resting picture), interpolated between blocks. The whole field fades during a
 *     recording and the lamp is uneven; cells that respond, or fade more than the glass, do not enter it.
 *   - The rims of shifted objects are taken off (CalciumFrameView.withoutRims).
 *   - Scale. Dark up to START_SD noise SDs; full at the level FULL_SHARE of the brightest frame's pixels stay under, and never
 *     closer to the start than SPAN_SD more SDs.
 *  Display only: no measurement uses any of this. */
@CompileStatic
class CalciumRise {
    static final int BOX = 4, FRAMES = 2, BLOCK = 128
    static final double BACK_SHARE = 0.4d, START_SD = 3d, SPAN_SD = 3d, FULL_SHARE = 0.95d
    final int width, height, frames, gw, gh
    final float[] rest
    final boolean[] back
    final double noise, start, full
    private final float[] a, b, grid, sample

    /** rest: the resting picture (mean of frames restFirst..restLast); brightest: the frame where the cells are brightest on average. */
    CalciumRise(ImageStack stack, float[] rest, int restFirst, int restLast, int brightest) {
        width = stack.getWidth(); height = stack.getHeight(); frames = stack.getSize(); this.rest = rest
        int n = width * height
        gw = (int) ((width + BLOCK - 1) / BLOCK); gh = (int) ((height + BLOCK - 1) / BLOCK)
        a = new float[n]; b = new float[n]; grid = new float[gw * gh]; sample = new float[BLOCK * BLOCK]; back = new boolean[n]
        // background pixels: the dimmest share of each block in the smoothed resting picture
        System.arraycopy(rest, 0, a, 0, n); box(a, b, width, height, BOX)
        for (int by = 0; by < gh; by++) for (int bx = 0; bx < gw; bx++) {
            int y1 = Math.min(height, (by + 1) * BLOCK), x1 = Math.min(width, (bx + 1) * BLOCK), k = 0
            for (int y = by * BLOCK; y < y1; y++) for (int x = bx * BLOCK; x < x1; x++) sample[k++] = a[y * width + x]
            if (k == 0) continue
            Arrays.sort(sample, 0, k); float cut = sample[(int) (BACK_SHARE * (k - 1))]
            for (int y = by * BLOCK; y < y1; y++) for (int x = bx * BLOCK; x < x1; x++) if (a[y * width + x] <= cut) back[y * width + x] = true
        }
        // camera noise of one averaged frame, from two neighbouring resting frames (their difference holds no cells, only noise)
        int second = Math.max(1, Math.min(frames, restLast)), first = second > 1 ? second - 1 : Math.min(frames, 2)
        double oneFrame = 0d
        if (first != second) {
            ImageProcessor p1 = stack.getProcessor(first), p2 = stack.getProcessor(second)
            for (int i = 0; i < n; i++) a[i] = (float) (p1.getf(i) - p2.getf(i))
            box(a, b, width, height, BOX)
            oneFrame = spread(a, b) / Math.sqrt(2d)
        }
        int resting = Math.max(1, restLast - restFirst + 1), averaged = Math.min(frames, 2 * FRAMES + 1)
        double sd = oneFrame * Math.sqrt(1d / averaged + 1d / resting)
        noise = sd > 0d ? sd : 1d
        start = START_SD * noise
        rise(stack, Math.max(1, Math.min(frames, brightest)), a)
        full = Math.max(start + SPAN_SD * noise, level(a, b, FULL_SHARE))
    }

    /** The rise of frame f above "no change", in gray levels, written into out (one value per pixel; not the buffer b). */
    void rise(ImageStack stack, int f, float[] out) {
        int n = width * height, from = Math.max(1, f - FRAMES), to = Math.min(frames, f + FRAMES)
        Arrays.fill(out, 0f)
        for (int g = from; g <= to; g++) {
            ImageProcessor ip = stack.getProcessor(g); Object pixels = ip.getPixels()
            if (pixels instanceof short[]) { short[] p = (short[]) pixels; for (int i = 0; i < n; i++) out[i] = (float) (out[i] + (p[i] & 0xffff)) }
            else if (pixels instanceof byte[]) { byte[] p = (byte[]) pixels; for (int i = 0; i < n; i++) out[i] = (float) (out[i] + (p[i] & 0xff)) }
            else if (pixels instanceof float[]) { float[] p = (float[]) pixels; for (int i = 0; i < n; i++) out[i] = (float) (out[i] + p[i]) }
            else for (int i = 0; i < n; i++) out[i] = (float) (out[i] + ip.getf(i))
        }
        float share = (float) (1d / (to - from + 1))
        for (int i = 0; i < n; i++) out[i] = (float) (out[i] * share - rest[i])
        box(out, b, width, height, BOX)
        // "no change": the middle value of the background pixels of each block
        double sum = 0d; int filled = 0
        for (int by = 0; by < gh; by++) for (int bx = 0; bx < gw; bx++) {
            int y1 = Math.min(height, (by + 1) * BLOCK), x1 = Math.min(width, (bx + 1) * BLOCK), k = 0
            for (int y = by * BLOCK; y < y1; y += 2) for (int x = bx * BLOCK; x < x1; x += 2) { int i = y * width + x; if (back[i]) sample[k++] = out[i] }
            if (k > 0) { Arrays.sort(sample, 0, k); grid[by * gw + bx] = sample[(int) (k / 2)]; sum += grid[by * gw + bx]; filled++ } else grid[by * gw + bx] = Float.NaN
        }
        float fallback = filled > 0 ? (float) (sum / filled) : 0f
        for (int k = 0; k < grid.length; k++) if (Float.isNaN(grid[k])) grid[k] = fallback
        for (int y = 0; y < height; y++) {
            double fy = (y + 0.5d) / BLOCK - 0.5d; int y0 = (int) Math.floor(fy); float ty = (float) (fy - y0)
            int ya = Math.max(0, Math.min(gh - 1, y0)) * gw, yb = Math.max(0, Math.min(gh - 1, y0 + 1)) * gw, row = y * width
            for (int x = 0; x < width; x++) {
                double fx = (x + 0.5d) / BLOCK - 0.5d; int x0 = (int) Math.floor(fx); float tx = (float) (fx - x0)
                int xa = Math.max(0, Math.min(gw - 1, x0)), xb = Math.max(0, Math.min(gw - 1, x0 + 1))
                float top = (float) (grid[ya + xa] + tx * (grid[ya + xb] - grid[ya + xa])), bottom = (float) (grid[yb + xa] + tx * (grid[yb + xb] - grid[yb + xa]))
                out[row + x] = (float) (out[row + x] - (top + ty * (bottom - top)))
            }
        }
        CalciumFrameView.withoutRims(out, b, width, height, 0f)
    }

    /** Frame f as the views show it: 0 (no rise) to 255 (full) in each entry of shade. */
    void paint(ImageStack stack, int f, int[] shade) {
        rise(stack, f, a)
        float low = (float) start, scale = (float) (255d / Math.max(1e-9d, full - start))
        for (int i = 0; i < shade.length; i++) { int v = (int) ((a[i] - low) * scale); shade[i] = v < 0 ? 0 : v > 255 ? 255 : v }
    }

    /** Mean over a box of (2r+1) x (2r+1) pixels, in place; the picture's edge repeats. work: a buffer of the same size. */
    static void box(float[] a, float[] work, int width, int height, int r) {
        if (r < 1 || width < 2 || height < 2) return
        float size = (float) (2 * r + 1)
        for (int y = 0; y < height; y++) {
            int row = y * width; double s = 0d
            for (int k = -r; k <= r; k++) s += a[row + Math.min(width - 1, Math.max(0, k))]
            for (int x = 0; x < width; x++) { work[row + x] = (float) (s / size); s += a[row + Math.min(width - 1, x + r + 1)] - a[row + Math.max(0, x - r)] }
        }
        for (int x = 0; x < width; x++) {
            double s = 0d
            for (int k = -r; k <= r; k++) s += work[Math.min(height - 1, Math.max(0, k)) * width + x]
            for (int y = 0; y < height; y++) { a[y * width + x] = (float) (s / size); s += work[Math.min(height - 1, y + r + 1) * width + x] - work[Math.max(0, y - r) * width + x] }
        }
    }
    /** The value below which this share of a spread sample of a lies. work: a buffer at least as long as the sample. */
    static double level(float[] a, float[] work, double share) {
        int step = Math.max(1, (int) (a.length / 200000)), k = 0
        for (int i = 0; i < a.length; i += step) work[k++] = a[i]
        if (k == 0) return 0d
        Arrays.sort(work, 0, k)
        return work[(int) Math.max(0L, Math.min((long) (k - 1), Math.round(share * (k - 1))))]
    }
    /** A robust SD of a (1.4826 x the median distance from the median), from a spread sample. */
    static double spread(float[] a, float[] work) {
        double middle = level(a, work, 0.5d)
        int step = Math.max(1, (int) (a.length / 200000)), k = 0
        for (int i = 0; i < a.length; i += step) work[k++] = (float) Math.abs(a[i] - middle)
        if (k == 0) return 0d
        Arrays.sort(work, 0, k)
        return 1.4826d * work[(int) (k / 2)]
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
    /** A green frame ("GFP"): "gray" becomes the green channel alone. The green of "Response in magenta" without the magenta. */
    static void green(int[] gray) { for (int i = 0; i < gray.length; i++) gray[i] = (gray[i] & 255) << 8 }
    /** A gray frame from shades of 0 to 255 (CalciumRise.paint). */
    static void gray(int[] shade, int[] rgb) { for (int i = 0; i < rgb.length; i++) { int v = shade[i] & 255; rgb[i] = (v << 16) | (v << 8) | v } }
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
    /** A response must hold for 2 samples in a row. The length is counted in samples, not seconds: "for 1 s" was 2 samples at
     *  exactly 1 frame per second and 3 at 0.99995 s per frame, so two recordings of one series were held to different rules. */
    static Map window(String role,String stimulus,double baselineStart,double onset,double end) {
        [Role:role,Stimulus:stimulus,Baseline_start_s:baselineStart,Baseline_end_s:onset,Stimulus_onset_s:onset,Start_s:onset,End_s:end,
         Floor_local_fraction:0.10d,SD_multiplier:3d,Duration_mode:'samples',Min_duration_s:1d,Consecutive_samples:2,Max_gap_s:null,Min_baseline_samples:3,Noise_floor_local_fraction:0d,Expected_interval_s:null]
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
    ImagePlus movie; java.awt.image.BufferedImage snapshot, frameImage; float[] rest, work, work2; double[] rawRange; CalciumRise riseView; int[] riseShade
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
        riseView=new CalciumRise(stack,rest,restFirst,restLast,peak)
        double low=CalciumFrameView.level(resting,null,0.005d),high=Math.max(CalciumFrameView.level(resting,null,0.9995d),CalciumFrameView.level(brightest,null,0.9995d));if(!(high>low))high=low+1d
        rawRange=[low,high] as double[]
        frameImage=new java.awt.image.BufferedImage(stack.getWidth(),stack.getHeight(),java.awt.image.BufferedImage.TYPE_INT_RGB);work=new float[stack.getWidth()*stack.getHeight()];work2=new float[work.length];riseShade=new int[work.length]
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
            int[] shown=((java.awt.image.DataBufferInt)frameImage.getRaster().getDataBuffer()).getData()
            if(mode==1){riseView.paint(movie.getStack(),frame,riseShade);CalciumFrameView.gray(riseShade,shown)}
            else CalciumFrameView.paint(movie.getStack().getProcessor(frame),null,rawRange[0],rawRange[1],shown,work,work2,false)
            cellPicture.image=frameImage;cellPicture.onGray=true
            cellPicture.note=mode==1?'frame '+frame+': rise above the mean of frames '+restFirst+'–'+restLast+', averaged over neighbouring pixels and frames':'frame '+frame+' as recorded'
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
        if(preset.selectedIndex!=1){fields.protocol_name.text=legacy?'Legacy consecutive-sample protocol':'Slow evoked response — provisional';validation.selectedItem=legacy?'user_defined':'provisional';for(int r=0;r<windowModel.rowCount;r++){['Floor_local_fraction':0.10d,'SD_multiplier':3d,'Duration_mode':'samples','Min_duration_s':1d,'Consecutive_samples':2].each{k,v->windowModel.setValueAt(v,r,CalciumResponseProtocol.COLUMNS.indexOf(k))};if(legacy){['Min_baseline_samples':3,'Noise_floor_local_fraction':0d,'Expected_interval_s':null].each{k,v->windowModel.setValueAt(v,r,CalciumResponseProtocol.COLUMNS.indexOf(k))};windowModel.setValueAt(windowModel.getValueAt(r,5),r,3);windowModel.setValueAt(windowModel.getValueAt(r,5),r,4);if(r==1)windowModel.setValueAt('KCl',r,1)}}}
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
                try{loadResults(answer);appliedFingerprint=submittedFingerprint;changed();String empty=CalciumAnalysisSettings.noEvaluableCells(answer);status.text=empty?:'Complete. Protocol, inputs and tables saved in '+answer.name+'.'+CalciumAnalysisSettings.percentResponding(answer,measurement);status.foreground=empty?CalciumTheme.WARN:CalciumTheme.INK;progress.value=100;progress.string='Complete';tabs.selectedIndex=1}catch(Throwable t){status.text='Could not display results: '+t.message;progress.string='Display error; see saved outputs'}
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
        try{file.setText(groovy.json.JsonOutput.toJson([format:REMEMBERED_FORMAT,recording:''+data.recording,frames:data.frames as int,baseline_start:fijiFrame(data,'Start'),baseline_end:fijiFrame(data,'End'),
            values:draft.values,protocol:draft.protocol]),'UTF-8')}catch(Exception ignored){}
        rememberStimulusFrame(data,draft)
    }
    /** 2: the rule's length is counted in samples by default. Settings remembered by an earlier build carry no format. */
    static final int REMEMBERED_FORMAT=2
    /** A rule still on the former default, "for 1 s", becomes the present default, 2 samples in a row. A length set by hand to
     *  anything else is kept. Returns whether the rule changed. */
    static boolean toSampleRule(Map protocol){
        boolean changed=false
        (protocol?.windows instanceof List?protocol.windows as List:[]).each{w->
            if(w instanceof Map&&w.Duration_mode=='seconds'&&CalciumCsv.number(''+w.Min_duration_s)==1d){w.Duration_mode='samples';w.Consecutive_samples=2;changed=true}}
        changed
    }
    /** @param note when given, note.ruleMoved says the remembered rule was on the former default and now counts samples */
    static Map remembered(Map data,Map note=null) {
        File file=rememberedFile(data);if(file==null||!file.isFile())return null
        try{def saved=new groovy.json.JsonSlurper().parseText(file.getText('UTF-8'))
            if(!(saved instanceof Map)||saved.recording!=''+data.recording||(saved.frames as int)!=(data.frames as int)||(saved.baseline_start as int)!=fijiFrame(data,'Start')||(saved.baseline_end as int)!=fijiFrame(data,'End'))return null
            if(!(saved.values instanceof Map)||!(saved.protocol instanceof Map))return null
            boolean moved=!(saved.format instanceof Number&&(saved.format as int)>=REMEMBERED_FORMAT)&&toSampleRule(saved.protocol as Map)
            Map candidate=[values:saved.values,protocol:saved.protocol,timestamps:null,sensitivity:null,negative_control:null]
            if(diagnose(candidate,data).errors)return null
            if(note!=null)note.ruleMoved=moved
            return copy(candidate)
        }catch(Exception ignored){return null}
    }
    /** A cell counts as quiet when the rise it needs to respond is at most this fraction of its baseline. */
    static final double QUIET_RISE=0.30d
    /** % responding to the primary stimulus, worked out in this one place, saved beside the tables and shown the same way on
     *  every tab. Rows: all cells; the cells chosen blind and those added after the response was shown, when that was
     *  recorded; and the quiet cells (required rise at most QUIET_RISE), where noisy cells are set aside whether they responded
     *  or not, so the row is not tilted toward responders. Each row is given twice: with every cell, and without the cells whose
     *  background-corrected fluorescence went below zero (dF/F0 under -1 at some frame).
     *  @param traces dF/F0 by cell when already read; otherwise normalized_traces.csv is read
     *  @return null when the analysis has no calls */
    static Map tally(File dir,Map measurement,Map traces=null) {
        try {
            def calls=CalciumCsv.read(new File(dir,'responder_calls.csv'));List h=calls[0];int cell=h.indexOf('Cell'),role=h.indexOf('Role'),status=h.indexOf('Status'),need=h.indexOf('Threshold_increase')
            if(cell<0||status<0)return null
            Set low=new LinkedHashSet()
            if(traces!=null)traces.each{name,values->if((values as List).any{it instanceof Number&&Double.isFinite(it as double)&&(it as double)<-1d})low<<(''+name)}
            else{File file=new File(dir,'normalized_traces.csv')
                if(file.isFile()){def saved=CalciumCsv.read(file);List names=saved[0]
                    for(int c=1;c<names.size();c++){boolean negative=false;for(int r=1;r<saved.size()&&!negative;r++){double v=CalciumCsv.number(''+saved[r][c]);if(Double.isFinite(v)&&v<-1d)negative=true};if(negative)low<<(''+names[c])}}}
            Map tags=measurement?.chosenBlind instanceof Map?measurement.chosenBlind as Map:null
            def blank={[counted:[0,0,0],without:[0,0,0]]}   // responded, evaluable, registered
            Map rows=[all:blank(),blind:blank(),later:blank(),quiet:blank()];List under=[]
            calls.drop(1).findAll{role<0||it[role]=='Stimulus'}.each{row->
                String name=''+row[cell],call=''+row[status];boolean responded=call=='Responded',evaluable=responded||call=='Did not respond'
                if(tags!=null&&!tags.containsKey(name))return   // a cell that is not in this measurement
                boolean below=evaluable&&low.contains(name)     // an excluded cell is out of the percentages already
                if(below)under<<name
                double rise=need<0?Double.NaN:CalciumCsv.number(''+row[need])
                List into=['all'];if(tags!=null)into<<(tags[name]?'blind':'later');if(evaluable&&Double.isFinite(rise)&&rise<=QUIET_RISE)into<<'quiet'
                into.each{String key->List both=[rows[key].counted];if(!below)both<<rows[key].without   // (a "?[" would read as a safe index)
                    both.each{List t->t[2]++;if(evaluable)t[1]++;if(responded)t[0]++}}
            }
            List out=[[key:'all',label:'All cells']+rows.all]
            if(tags!=null)out+=[[key:'blind',label:'Chosen blind']+rows.blind,[key:'later',label:'Added after the response was shown']+rows.later]
            if(need>=0)out<<([key:'quiet',label:'Quiet enough to detect a '+Math.round(100*QUIET_RISE)+' % rise']+rows.quiet)
            return [quietRise:QUIET_RISE,belowZero:under,tagged:tags!=null,rows:out]
        }catch(Throwable ignored){return null}
    }
    static String percentOf(List t){(t[1] as int)>0?String.format(Locale.US,'%.1f',100d*(t[0] as int)/(t[1] as int)):''}
    /** Saves the tally beside the measurement tables and returns the sentence shown when the analysis finishes: all cells first,
     *  then the two groups when cells were added after the response was shown. '' when it was not recorded how cells were chosen. */
    static String percentResponding(File dir,Map measurement) {
        try {
            Map found=tally(dir,measurement);if(found==null)return ''
            List rows=found.rows as List
            if(measurement.folder instanceof File)CalciumCsv.write(new File(measurement.folder as File,measurement.recording+'_percent_responding_'+dir.name+'.csv'),
                [['cells','registered','evaluable','responded','percent_responding','cells_below_zero','evaluable_without_them','responded_without_them','percent_without_them']]+rows.collect{Map row->
                    List t=row.counted as List,w=row.without as List;[row.label,t[2],t[1],t[0],percentOf(t),(t[2] as int)-(w[2] as int),w[1],w[0],percentOf(w)]})
            if(!found.tagged)return ''
            def of={String key->(rows.find{it.key==key} as Map).counted as List}
            def line={String label,List t->label+': '+t[0]+' of '+t[1]+' responded'+((t[1] as int)>0?' ('+percentOf(t)+' %)':'')}
            List blind=of('blind'),later=of('later'),all=of('all')
            if((later[2] as int)==0)return ' All '+blind[2]+' cells were chosen blind: '+blind[0]+' of '+blind[1]+' responded ('+percentOf(blind)+' %).'
            return ' '+line('All cells',all)+'. '+line('Chosen blind',blind)+'. '+line('Added after the response was shown',later)+'.'
        }catch(Throwable ignored){return ''}
    }
    /** How the page was last arranged (pane sizes, explanations on or off). The page's address changes every time, so the
     *  browser cannot remember these itself; they are kept with Fiji's preferences. */
    static final Map PAGE_PREFS=['step4.side':/\d{2,4}px/,'step4.split':/\d{1,2}(\.\d)?%/,'step4.height':/\d{3,4}px/,'step4.help':/[01]/,
        'step4.every':/\d{1,4}/,'step4.outline':/\d(\.\d)?/,'step4.format':/(png|tiff)/,'step4.size':/(screen|\d{1,2}(\.\d)?x\d{1,2}(\.\d)?)/,'step4.dpi':/(150|300|600)/,'step4.compression':/(lzw|zip|none)/,
        'step4.figure.plot':/\{[^\r\n]{0,1500}\}/,'step4.figure.heatmap':/\{[^\r\n]{0,1500}\}/]   // the figure editor's settings, one line of JSON each; the page checks every value
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
    float[] rest,work,work2,background;int[] tintBuffer;double[] rawRange,flatRange,meanTrace=new double[0];CalciumRise riseView
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
        Map carried=[:];Map before=CalciumAnalysisSettings.remembered(data,carried)
        restFirst=(data.baselineStart?:1) as int;restLast=draft.values.baseline_end as int
        if(before!=null){
            double a=CalciumCsv.number(''+before.values.baseline_start),b=CalciumCsv.number(''+before.values.baseline_end)
            if(Double.isFinite(a)&&Double.isFinite(b)&&a==Math.rint(a)&&b==Math.rint(b)&&a>=1&&b>=a&&b<=(data.frames as int)){
                fijiDraft=CalciumAnalysisSettings.copy(draft);stimulusFrom='last-session';draft=before;restFirst=a as int;restLast=b as int;message='Your last settings for this recording are filled in.'+(carried.ruleMoved?' The rule now counts samples: 2 in a row. It was "for 1 s".':'')}
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
                    if(empty){message=empty;tone='warn'} else {message='Complete. Protocol, input snapshots and tables saved in '+output.name+'.'+CalciumAnalysisSettings.percentResponding(output,measurement);tone='good'}}
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
            corrected:corrected,raw:raw,chosenBlind:(measurement.chosenBlind instanceof Map?measurement.chosenBlind:null),tally:CalciumAnalysisSettings.tally(dir,measurement,traces),appliedDraft:submitted==null?null:CalciumAnalysisSettings.copy(submitted),protocol:protocol('protocol.json'),appliedProtocol:protocol('applied_protocol.json'),manifest:manifest,appliedFingerprint:applied,tables:tables]
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
        if(!(mode in ['raw','flat','green','magenta','change']))throw new IllegalArgumentException('Frame mode must be raw, flat, green, magenta or change.')
        def stack=movie.getStack();int f=wanted
        if(rest==null) {
            rest=CalciumFrameView.mean(stack,restFirst,restLast,movedSet())   // moved frames would smear the resting picture
            // the frame where the cells are brightest on average sets the display range, so it does not flicker from frame to frame
            int peak=Math.max(1,Math.min(n,(restLast+n).intdiv(2)));double top=Double.NEGATIVE_INFINITY
            meanTrace.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top&&i<n){top=v;peak=i+1}}
            def brightest=stack.getProcessor(peak);def resting=new ij.process.FloatProcessor(width,height,rest)
            double low=CalciumFrameView.level(resting,null,0.005d),high=Math.max(CalciumFrameView.level(resting,null,0.9995d),CalciumFrameView.level(brightest,null,0.9995d));if(!(high>low))high=low+1d
            rawRange=[low,high] as double[]
            // "Background removed", "GFP" and the green of "Response in magenta": the smooth background of the resting picture
            // (Process > Subtract Background, rolling ball, as for the plugin's all-frames picture) is taken off every frame, so the
            // lamp's uneven light goes and the cells stand out. Shown from 0 to the level 97 % of the brightest frame's pixels stay
            // under, so very bright objects do not set the range. The response itself ("Response in magenta", "Change from
            // baseline") is CalciumRise: one definition for both. Display only.
            def ground=new ij.process.FloatProcessor(width,height,(float[])rest.clone())
            new ij.plugin.filter.BackgroundSubtracter().rollingBallBackground(ground,CalciumSteps.BALL_RADIUS,true,false,false,true,true)
            background=(float[])ground.getPixels()
            double flatTop=CalciumFrameView.level(brightest,background,0.97d);if(!(flatTop>0d))flatTop=1d
            flatRange=[0d,flatTop] as double[]
            riseView=new CalciumRise(stack,rest,restFirst,restLast,peak)
            frameImage=new java.awt.image.BufferedImage(width,height,java.awt.image.BufferedImage.TYPE_INT_RGB);work=new float[width*height];work2=new float[width*height]
        }
        int[] pixels=((java.awt.image.DataBufferInt)frameImage.getRaster().getDataBuffer()).getData()
        if(mode=='flat'||mode=='green'){
            CalciumFrameView.paint(stack.getProcessor(f),background,flatRange[0],flatRange[1],pixels,work,work2,false)
            if(mode=='green')CalciumFrameView.green(pixels)
            return frameImage
        }
        if(mode=='magenta'||mode=='change'){
            if(tintBuffer==null||tintBuffer.length!=pixels.length)tintBuffer=new int[pixels.length]
            riseView.paint(stack,f,tintBuffer)
            if(mode=='change'){CalciumFrameView.gray(tintBuffer,pixels);return frameImage}
            CalciumFrameView.paint(stack.getProcessor(f),background,flatRange[0],flatRange[1],pixels,work,work2,false)
            CalciumFrameView.tint(pixels,tintBuffer);return frameImage
        }
        CalciumFrameView.paint(stack.getProcessor(f),null,rawRange[0],rawRange[1],pixels,work,work2,false)
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
    static final String SHA256 = "a6f5057bc06f04f8ce80bcf41e816cb87ec93f286cd3fa0198818850975ec56c"
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
    def placeNote = CalciumTheme.label('', CalciumTheme.INK)
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
        placeNote.setVisible(false); placeNote.setAlignmentX(0f); placeNote.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 8, 0))
        [placeNote, buttons, javax.swing.Box.createVerticalStrut(6), autoUpdateBox, updateNote, installBar].each { body.add(it) }
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
    /** What happened to the plugin's file when it started (CalciumInstaller.settle); stays above the update messages. */
    void notice(String text) {
        placeNote.setText("<html><body style='width:300px'>"+text+'</body></html>'); placeNote.setVisible(true)
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
                        (completed.aside ? '<br>After the restart it is at Plugins › '+CalciumInstaller.escape(CalciumRelease.command(completed.target.path))+'.' : '')+
                        (completed.backup ? '<br>The previous file was backed up in calcium-plugin-backups inside Fiji’s folder.' : '')+
                        (completed.left ? '<br>Fiji could not move an older copy: '+CalciumInstaller.escape(completed.left.join(', '))+'. Delete it so that the menu shows the plugin once.' : ''))
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
    /** Shown in the updates card when the plugin's file was moved as it started. */
    static String notice = null
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
            if (!Double.isFinite(a) || !Double.isFinite(b) || a != Math.floor(a) || b != Math.floor(b) || a < 1 || b < a || b > n) { IJ.error(CalciumRelease.NAME, "Choose whole frame numbers from 1 to " + n + ", the first no later than the last."); continue }
            first = (int) a; last = (int) b
            if (rise && last >= (int) steps.onset && !IJ.showMessageWithCancel(CalciumRelease.NAME, "Frame " + last + " is at or after the rise near frame " + steps.onset + ".\nFrames in which cells respond are not blind: cells that respond are brighter in them.\n\nUse these frames anyway?")) continue
            if (offer[1] - last >= 20 && IJ.showMessageWithCancel(CalciumRelease.NAME, "The plugin offers frames up to " + offer[1] + ": " + (rise ? "the cells start to rise near frame " + steps.onset : "no response was found") + ".\nYou chose frames " + first + "-" + last + ". With more frames the picture is clearer and more cells are found.\n\nOK: use frames " + first + "-" + offer[1] + ".    Cancel: keep " + first + "-" + last + ".")) last = offer[1]
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
            def dialog = new javax.swing.JDialog((java.awt.Frame) null, CalciumRelease.NAME, true)
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
            if (notice) updates.notice(notice)
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

// The updater of 2026.11.22.1 and earlier puts the new version under Plugins > Calcium Imaging, and the ZIP at
// Plugins > Ca2+ Imaging (Wang Lab). Move it to its place.
Map placed = null
if (!testGui && IJ.getDirectory("imagej") != null) {
    try { placed = new CalciumInstaller(applicationDirectory: new File(IJ.getDirectory("imagej"))).settle() }
    catch (Exception problem) { IJ.log(CalciumRelease.NAME + ": the plugin could not be moved to Plugins > " + CalciumRelease.NAME + " (" + (problem.message ?: problem) + ")") }
}
if (placed?.moved || placed?.aside) {
    String menu = CalciumRelease.command(placed.target.path)
    CalciumStartDialog.notice = placed.moved ? "<b>The plugin is now at Plugins › " + CalciumInstaller.escape(menu) + ".</b> Restart Fiji before starting it from the menu again." :
        "An older copy of the plugin was moved to calcium-plugin-backups inside Fiji’s folder, so that the menu shows it once."
    IJ.log(CalciumRelease.NAME + ": " + (placed.moved ? "moved " + placed.moved + " to " + placed.target : "") +
        placed.aside.collect { k, v -> "; moved " + k + " to " + v }.join("") + (placed.left ? "; could not move " + placed.left.join(", ") : ""))
}

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
if (imp == null || imp.getStackSize() < 2) { IJ.error(CalciumRelease.NAME, "No movie was opened."); return }

IJ.showStatus("Making the cells visible...")
if (imp.getBitDepth() == 32) IJ.log(CalciumRelease.NAME + ": this movie is 32-bit, so it looks like a processed copy. Background is meant to be taken on the original movie.")
def steps = new CalciumSteps()
steps.prepare(imp, look || blindStart)   // blind: the response is looked for only to offer when the stimulus came; it is never shown while choosing
if (blindStart) {
    int[] frames = CalciumStartDialog.askBlindFrames(steps)
    if (frames == null) { IJ.showStatus(CalciumRelease.NAME + ": cancelled"); return }
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
// likely, unsure and very bright: share on audited cells with the grading of 2026.11.20.20 (71 %, 31 %, 12 % of 505, 225 and 209 blind outlines on four recordings)
def GRADE_SHARE = ["about 8 in 10 such outlines were real cells in tests", "about 7 in 10 such outlines were real cells in tests", "about 3 in 10 such outlines were real cells in tests", "about 1 in 10 such outlines were real cells in tests"]
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
def legendRow = { javax.swing.Icon ic, String s, String kind, String what ->
    def part = new CalciumSwitch(s, ic); legendEntries[kind] = [part: part, text: s]
    part.setToolTipText("<html>" + what + "<br>Click to hide or show these on the picture. Double-click to see only these, and again for everything.</html>")
    part.clicked = { boolean twice -> if (legendClicked != null) legendClicked.call(kind, twice) }
    part
}
def section = { String number, String title, boolean collapsible = true -> new CalciumCard(number, title, collapsible) }
def row = { java.awt.Component... cs -> CalciumTheme.row(cs) }
def note = { String s -> def l = new CalciumText(s); l.setFont(small); l }
def hint = { String s -> def l = note(s); l.setForeground(CalciumTheme.MUTED); l }
def btn = { String s -> CalciumTheme.button(s) }
def box = { String s, boolean on -> def b = new javax.swing.JCheckBox(s, on); b.setFont(small); b.setOpaque(false); b.setAlignmentX(0f); CalciumTick.dress(b); b }

def frame = new CalciumWindow(CalciumRelease.NAME)
frame.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
def root = new CalciumScrollPanel(); root.setBorder(javax.swing.BorderFactory.createEmptyBorder(10, 10, 10, 10)); root.setBackground(CalciumTheme.PAGE); root.setOpaque(true)

// One card per step. Step 1 is already done; of steps 2 to 4 only one is open at a time, so the window fits the screen.
def s1 = section("1", "Recording", false); s1.setDone(true)
def recordingName = new CalciumName(imp.getTitle())
def recordingNote = note("&nbsp;")
// "The field moved (about 16 px) during frames 1–10.", or nothing (also for a window built without the movement check)
// For Step 4 the movement is counted from the frames the cells were chosen on whenever such frames were confirmed, also after cells were added with the response.
def motionForPage = { chosen -> try { return chosen != null && steps.motionShifts != null ? CalciumMotion.summary(steps.motionShifts, chosen.first as int, chosen.last as int) : steps.motion } catch (Throwable ignored) { return null } }
def movedNote = { -> try { steps.sayMotion(); return CalciumMotion.sentence(steps.motion) } catch (Throwable ignored) { return "" } }
def describeRecording = { ->
    recordingName.setWhole(imp.getTitle())
    recordingNote.setText("<html><font color='" + MUTED + "'>" + steps.n + " frames · " + (steps.blind ? "cells chosen blind on frames " + steps.blindFirst + "–" + steps.blindLast : steps.found ? "response near frame " + steps.onset + " (magenta)" : (steps.looked ? "no response found (no magenta)" : "response not looked for")) + "</font>" + (movedNote() ? "<br><font color='#9A5B00'>" + movedNote() + "</font>" : "") + "</html>")
    recordingNote.setToolTipText(imp.getTitle() + ". " + (steps.blind ? "Cells are chosen on frames " + steps.blindFirst + "–" + steps.blindLast + ", before the response" + (steps.found ? "; the cells start to rise near frame " + steps.onset : "") + ". The response is not shown while choosing." : steps.found ? "Brightness rises near frame " + steps.onset + "; the response can be shown in magenta." : (steps.looked ? "No response was found, so there is no magenta. Every cell is still shown in green." : "The response has not been looked for. Choose Response in Step 2 to look for it.")))
}
describeRecording()
def changeBtn = btn("‹ Change recording")
changeBtn.setToolTipText("Back to the first screen, to choose another recording or results folder")
s1.body.add(recordingName); s1.body.add(recordingNote); s1.body.add(row(changeBtn))

def s2 = section("2", "Choose the cells")
String selectionBasis = steps.blind ? steps.blindBasis() : steps.found ? 'Response-assisted review' : 'All-frames structural review'
// What the colour picture shows while cells are chosen. The choice is recorded with the measurement as the selection basis.
def viewAllBtn = btn('All frames'), viewBaselineBtn = btn('Before response'), viewResponseBtn = btn('Response')
def viewSwitch = new CalciumSegments(viewAllBtn, viewBaselineBtn, viewResponseBtn)   // one joined switch; the line under it says what the choice means
viewAllBtn.setToolTipText('The average of all frames with the background removed: cells by their structure')
viewBaselineBtn.setToolTipText('<html>Only frames before the response; you are asked which. Suggestions are made on this picture too, traces are hidden<br>and grades use shape, size and brightness only, so whether a cell is chosen cannot depend on its response.</html>')
def blindNote = { -> 'Blind: frames ' + steps.blindFirst + '–' + steps.blindLast + ' only, before the response. Traces are hidden; grades use shape, size and brightness.' }
viewResponseBtn.setToolTipText('<html>The response in magenta: cells that brighten after the stimulus turn magenta.<br>Choosing cells with this on is recorded as response-assisted selection,<br>which can raise percent responding.</html>')
def selectionNote = hint('&nbsp;'); selectionNote.lines = 2
def showView = { String which, String text ->
    [all: viewAllBtn, baseline: viewBaselineBtn, response: viewResponseBtn].each { key, button -> CalciumTheme.choose(button, key == which) }
    selectionNote.setText("<html>" + text + "</html>")
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
def edgeBox = box("Leave out cells within " + CalciumSteps.EDGE_MARGIN + " px of the edge", true)
edgeBox.setToolTipText("<html>At the edge a cell's background can only be taken from the inner, better-lit side,<br>so the cell comes out too dim. Untick to suggest, accept and measure edge cells like the others.</html>")
def status = new CalciumText("<b>Cells chosen: 0</b>"); status.setFont(small)
def info = note("&nbsp;")
def trace = new CalciumTracePanel(); trace.setAlignmentX(0f); trace.setFont(small)
if (steps.blind) trace.show(null, "Traces are hidden while choosing cells blind", 0, 0)
trace.setVisible(false)   // it takes its room when there is a trace to show: an empty frame only pushes the rest down
// Step 2 leads through its own order: suggest, accept the likely ones, check the dashed ones, add the cells that only show when
// they respond, continue. This line always names the one thing to do next (set in refresh, from what has been done so far).
def nextLine = note("&nbsp;"); nextLine.lines = 2   // one or two lines in every state at the usual width: the four actions under it stay where they are
def sayNext = { String text -> nextLine.setText("<html><b><font color='" + CalciumTheme.hex(CalciumTheme.BLUE) + "'>Next:</font></b> " + text + "</html>") }
def addRespondersBtn = btn("Add cells that respond…")
String respondersTip = "<html>Shows the response in magenta and looks for the cells that only show when they respond.<br>Your cells so far stay recorded as chosen blind; the ones you add are counted apart in Step 4.</html>"
addRespondersBtn.setToolTipText(respondersTip)
boolean detecting = false, autoSuggestAsked = false
s2.body.add(nextLine); s2.body.add(javax.swing.Box.createVerticalStrut(4))
s2.body.add(row(CalciumTheme.label('Picture:', CalciumTheme.MUTED), viewSwitch)); s2.body.add(selectionNote)
def cpLabel = hint("&nbsp;")
def setupBtn = btn("Set up Cellpose (one time)"); def setupRow = row(setupBtn)
def detectionProgress = new CalciumOutlineProgress()
def continueBtn = btn("Continue to measure  ›")
showView(steps.blind ? 'baseline' : steps.found ? 'response' : 'all', steps.blind ? blindNote() : steps.found ? 'Magenta: cells that brightened near frame ' + steps.onset + '. Picking cells this way is recorded as response-assisted.' : 'Choosing cells on the average of all frames.')
// The four things this step asks for keep their places whatever has been done: suggest and accept on one line, add the
// cells that only show when they respond and continue on the next. The next one to press is blue; one that cannot be pressed
// now is gray, never gone. Each column is as wide as the longest label it ever shows, in bold, so nothing moves when a
// label or the blue button changes.
def actionMetrics = suggestBtn.getFontMetrics(CalciumTheme.font(java.awt.Font.BOLD))
int actionHigh = (int) suggestBtn.getPreferredSize().height
def columnWide = { List<String> longest -> (longest.collect { String one -> actionMetrics.stringWidth(one) }.max() as int) + 26 }
int leftWide = columnWide(["Suggest cells automatically", "Suggest cells again", "Add cells that respond…"]), rightWide = columnWide(["Accept responding and likely (8888)", "Continue to measure  ›"])
[suggestBtn, addRespondersBtn].each { it.setPreferredSize(new java.awt.Dimension(leftWide, actionHigh)) }
[acceptAllBtn, continueBtn].each { it.setPreferredSize(new java.awt.Dimension(rightWide, actionHigh)) }
s2.body.add(row(suggestBtn, acceptAllBtn)); s2.body.add(row(addRespondersBtn, continueBtn))
s2.body.add(cpLabel); s2.body.add(setupRow)   // only when Cellpose is missing: it is what Suggest needs
// Drawing by hand is used now and then: one line opens it, so the sequence above stays short.
def moreRows = new javax.swing.JPanel(); moreRows.setLayout(new javax.swing.BoxLayout(moreRows, javax.swing.BoxLayout.Y_AXIS)); moreRows.setOpaque(false); moreRows.setAlignmentX(0f)
moreRows.add(row(CalciumTheme.label("Draw:", CalciumTheme.MUTED), ovalBtn, freeBtn, addBtn))
def moreBtn = btn("Draw by hand  ▸"); moreBtn.setToolTipText("Oval, Freehand and Add outline")
moreRows.setVisible(false)
moreBtn.addActionListener { moreRows.setVisible(!moreRows.isVisible()); moreBtn.setText("Draw by hand  " + (moreRows.isVisible() ? "▾" : "▸")); frame.revalidate(); frame.repaint() }
s2.body.add(javax.swing.Box.createVerticalStrut(4)); s2.body.add(status)
// Your cells: what can be pressed on one line, how clicks and the edge are treated on the next (the counts above say how
// many are left out at the edge). Drawing by hand opens under them.
s2.body.add(row(removeBtn, undoBtn, moreBtn)); s2.body.add(row(removeMode, javax.swing.Box.createHorizontalStrut(8), edgeBox)); s2.body.add(editHint); s2.body.add(moreRows)
// What is drawn on the picture: one switch per kind, each with a sample of how that kind is drawn. (Two ticks, "Show
// suggestions" and "Show background", used to do the same job a second time.)
def legend = new javax.swing.JPanel(new java.awt.GridLayout(0, 3, 10, 0)); legend.setAlignmentX(0f); legend.setOpaque(false)
def acceptedIcon = new java.awt.image.BufferedImage(30, 14, java.awt.image.BufferedImage.TYPE_INT_ARGB)
def acceptedGraphics = acceptedIcon.createGraphics()
acceptedGraphics.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
acceptedGraphics.setColor(new Color(24, 30, 36)); acceptedGraphics.fillRoundRect(0, 0, 30, 14, 6, 6)
acceptedGraphics.setColor(Color.black); acceptedGraphics.setStroke(new java.awt.BasicStroke(5f)); acceptedGraphics.drawOval(6, 3, 18, 8)
acceptedGraphics.setColor(Color.white); acceptedGraphics.setStroke(new java.awt.BasicStroke(2f)); acceptedGraphics.drawOval(6, 3, 18, 8); acceptedGraphics.dispose()
legend.add(legendRow(new javax.swing.ImageIcon(acceptedIcon), "your cells", "cells", "The cells you have chosen: white outlines."))
def respondsLegend = legendRow(lineIcon(CYAN, 2, false, false), "responds", "g0", "Suggestions that rise after the stimulus on their own."); respondsLegend.setVisible(!steps.blind)   // no such grade while choosing blind
legend.add(legendRow(lineIcon(CYAN, 1, false, false), "likely a cell", "g1", "Suggestions graded likely a cell."))
legend.add(legendRow(lineIcon(GREY, 1, true, false), "unsure", "g2", "Suggestions graded unsure: dashed."))
legend.add(legendRow(lineIcon(BLUE, 1, false, true), "background", "background", "The pixels used as each cell's background, in blue."))
legend.add(legendRow(lineIcon(RED, 2, false, false), "precipitate?", "g3", "Suggestions that are very bright: probably precipitate."))
legend.add(respondsLegend)   // last, so the grid has no hole when it is hidden
def showLabel = CalciumTheme.label('Show:', CalciumTheme.MUTED); showLabel.setVerticalAlignment(javax.swing.SwingConstants.TOP); showLabel.setBorder(javax.swing.BorderFactory.createEmptyBorder(3, 0, 0, 0))
def showRow = new javax.swing.JPanel(new java.awt.BorderLayout(8, 0)); showRow.setOpaque(false); showRow.setAlignmentX(0f)
showRow.setBorder(javax.swing.BorderFactory.createEmptyBorder(3, 0, 3, 0))
showRow.add(showLabel, java.awt.BorderLayout.WEST); showRow.add(legend, java.awt.BorderLayout.CENTER)
showRow.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) showRow.getPreferredSize().height))
s2.body.add(showRow)
def hiddenLine = note(""); hiddenLine.setForeground(CalciumTheme.WARN); s2.body.add(hiddenLine)   // "Hidden on the picture now: ...", only while something is
// Under what is looked at all the time: the bar while outlines are made, then the report of how they were made
s2.body.add(detectionProgress)
s2.body.add(javax.swing.Box.createVerticalStrut(4)); s2.body.add(info); s2.body.add(trace)   // last: when a trace appears nothing under it has to move

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
// Blind selection alone leaves cells that only show when they respond without an outline (20 Nov 2026: it outlined 56-67 % of the
// audited cells on three recordings with a response, and the missed ones respond more often). After a blind-only measurement this
// button goes back and adds them in one press.
def addMissedBtn = btn("Add cells that respond…"); addMissedBtn.setVisible(false)
addMissedBtn.setToolTipText("<html>Goes back to the cells, shows the response and looks for the cells that only show when they respond.<br>The cells measured so far stay recorded as chosen blind.</html>")
s3.body.add(makeMovie); s3.body.add(saveMovie); s3.body.add(sensBox); s3.body.add(row(go, openBtn)); s3.body.add(result); s3.body.add(row(backToCellsBtn, addMissedBtn))
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
boolean blindOnlyChosen = false, respondersExplained = false
// (both answers are false in a window built by a test, which has no response and must not open a dialog)
def respondersWaiting = { -> try { return steps.blind == true && steps.found && !steps.riseIsMovement() } catch (Throwable ignored) { return false } }
def mayAsk = { -> try { return !testGui } catch (Throwable ignored) { return false } }
continueBtn.addActionListener({ e ->
    if (mayAsk.call() && !blindOnlyChosen && respondersWaiting.call() && addRespondersBtn.isEnabled()) {
        Object[] choices = ["Add cells that respond", "Measure these only"]
        int answer = javax.swing.JOptionPane.showOptionDialog(frame, "<html><body style='width:330px'>Cells that only show when they respond have no outline yet. On the recordings checked, blind selection alone outlined about 6 in 10 cells, and the missed ones respond more often.<br><br>Add them first? The cells you have stay recorded as <b>chosen blind</b>; the ones added are counted apart in Step 4.</body></html>", "Cells that respond", javax.swing.JOptionPane.DEFAULT_OPTION, javax.swing.JOptionPane.QUESTION_MESSAGE, null, choices, choices[0])
        if (answer == 0) { respondersExplained = true; addRespondersBtn.doClick(); return }
        if (answer != 1) return
        blindOnlyChosen = true
    }
    openStep(s3)
} as java.awt.event.ActionListener)
addMissedBtn.addActionListener({ e -> openStep(s2); respondersExplained = true; addRespondersBtn.doClick() } as java.awt.event.ActionListener)
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
// Wide enough for the sentences to read well. A narrower window still works: the text follows the width, down to the widest row of controls.
int beside = (int) scroll.getVerticalScrollBar().getPreferredSize().width   // room for the scrollbar: a window that scrolls down gets no bar along the bottom
frame.setSize(Math.max(520, frame.getWidth() + beside), frame.getHeight())
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
        IJ.log(CalciumRelease.NAME + ": the analysis page could not be opened in a browser (" + problem + "). Using the Fiji window instead.")
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
        def key = all.collect { it.getBounds().toString() + it.getContainedPoints().length }.join("|") + "#" + Arrays.toString(rm.getSelectedIndexes()) + "#" + hiddenKinds.sort(false).join(",") + "#" + proposals.size() + "#" + edgeBox.isSelected() + '#'+selectionBasis + '#' + responseSeen + (analysisMeasurement!=null ? '#'+measuredKey(all) : '')
        if (!force && key == lastKey) return
        lastKey = key
        boolean measuredNow = false   // the cells on the picture are the ones in the saved tables
        if (analysisMeasurement != null) {
            boolean current = measuredKey(all)+'#edge='+edgeBox.isSelected()+'#selection='+basisFor(all) == analysisMeasurement.key
            measuredNow = current
            analyzeBtn.setEnabled(current && go.isEnabled())
            s3.setDone(current); s4.setSummary(current ? "ready" : "measure again first")
            analysisNote.setText("<html>" + (current ? "Review the recording, check timing, and run analysis in your browser. Keep Fiji open." : "Cell outlines changed. Measure and save again before analysis.") + "</html>")
        }
        steps.alsoAvoid = edgeCells
        lastBg = steps.backgrounds(cells)
        int none = lastBg.count { it.length == 0 } as int
        def ov = new Overlay()
        if (cells.length > 0 && !hiddenKinds.contains('background')) {
            int sel = rm.getSelectedIndex(); int only = -1
            if (sel >= 0 && rm.getSelectedIndexes().length == 1) { def allR = rm.getRoisAsArray(); int k = -1; for (int i = 0; i <= sel && i < allR.length; i++) if (allR[i].isArea()) k++; only = sp.keptIdx.indexOf(k) }
            ov.add(steps.blueLayer(lastBg, only))
        }
        def left = unusedProposals(all)
        left.each { p ->
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
        status.setText("<html><b>Cells chosen: " + cells.length + "</b>" + (edgeCells.length > 0 ? " &nbsp;(+" + edgeCells.length + " at the edge, left out)" : "") + (none > 0 ? " &nbsp;(" + none + " without background)" : "") + (proposals.size() > 0 ? "<br><font color='" + MUTED + "'>Suggestions left: " + nLikely + " likely, " + left.count { it.grade == 2 } + " unsure, " + left.count { it.grade == 3 } + " very bright (precipitate?)" + (edgeSug > 0 ? "; " + edgeSug + " at the edge not shown" : "") + "</font>" : "") + "</html>")
        s2.setSummary(cells.length + (cells.length == 1 ? " cell chosen" : " cells chosen")); s2.setDone(cells.length > 0 && !s2.expanded)
        boolean withResponse = !steps.blind && steps.found && selectionBasis.startsWith('Response')
        acceptAllBtn.setText((withResponse ? "Accept responding and likely (" : "Accept likely (") + nLikely + ")"); acceptAllBtn.setEnabled(nLikely > 0)
        int nUnsure = left.count { it.grade == 2 && !(edgeBox.isSelected() && it.edge) } as int
        def hiddenNames = legendEntries.findAll { k, v -> hiddenKinds.contains(k) }.collect { k, v -> v.text }
        hiddenLine.setText("<html>" + (hiddenNames ? "Hidden on the picture now: " + hiddenNames.join(", ") + "." : "") + "</html>")
        addRespondersBtn.setEnabled(steps.blind == true && all.length > 0 && !detecting)
        addRespondersBtn.setToolTipText(steps.blind == true ? respondersTip : "The response is already shown or the cells are not chosen blind, so there is nothing to add this way.")
        boolean respondersNext = false
        try { respondersNext = respondersWaiting.call() && all.length > 0 && nLikely == 0 && !detecting; addMissedBtn.setVisible(respondersWaiting.call() && measuredNow) } catch (Throwable ignored) { }
        // one blue button: the next thing to press (cosmetic, so a window without these buttons carries on)
        try {
            String blue = detecting ? '' : respondersNext ? 'responders' : measuredNow ? '' : (proposals.isEmpty() && all.length == 0) || (!steps.blind && proposals.any { it.blind }) ? 'suggest' : nLikely > 0 ? 'accept' : all.length > 0 ? 'continue' : 'suggest'
            CalciumTheme.primary(suggestBtn, blue == 'suggest'); CalciumTheme.primary(acceptAllBtn, blue == 'accept'); CalciumTheme.primary(continueBtn, blue == 'continue'); CalciumTheme.primary(addRespondersBtn, blue == 'responders')
        } catch (Throwable ignored) { }
        if (detecting) sayNext("wait while the cell outlines are made (see the bar below).")
        else if (measuredNow && respondersNext) sayNext("these " + cells.length + " cells are measured. <b>Add cells that respond…</b> (no outline yet) and measure again, or go on to <b>Analyze responses</b>.")
        else if (measuredNow) sayNext("these " + cells.length + " cells are measured. Go on to <b>Analyze responses</b>, or change the cells here and measure again.")
        else if (proposals.isEmpty() && all.length == 0) sayNext("press <b>Suggest cells automatically</b>, or draw the cells yourself.")
        else if (!steps.blind && proposals.any { it.blind }) sayNext("press <b>Suggest cells again</b>: the outlines on the picture were made before the response, not with it.")
        else if (nLikely > 0) sayNext("press <b>" + (withResponse ? "Accept responding and likely" : "Accept likely") + " (" + nLikely + ")</b>.")
        else if (steps.blind && respondersNext) sayNext(nUnsure > 0 ? "click the dashed outlines that are real cells (" + nUnsure + " left), <b>Add cells that respond…</b> (no outline yet), then <b>Continue to measure</b>."
            : "press <b>Add cells that respond…</b>: cells that only show when they respond have no outline yet. Then <b>Continue to measure</b>.")
        else if (steps.blind) sayNext((nUnsure > 0 ? "click the dashed outlines that are real cells (" + nUnsure + " left). Then " : "") + "<b>Add cells that respond…</b> or <b>Continue to measure</b>.")
        else sayNext((nUnsure > 0 ? "click the dashed outlines that are real cells (" + nUnsure + " left), then " : "") + "<b>Continue to measure</b>.")
        undoBtn.setEnabled(!cellEdits.history.isEmpty())
    } catch (Throwable t) { status.setText("" + t) }
}
def paintLegend = { -> legendEntries.each { k, v -> v.part.setOn(!hiddenKinds.contains(k)) } }
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
def blindControls = { boolean on ->
    respondsLegend.setVisible(!on); sensBox.setVisible(on)
    if (!on && trace.vals == null) trace.show(null, "Move the mouse over a cell or a suggestion to see its trace", 0, 0)   // no longer "hidden while choosing blind"
    trace.setVisible(!on && trace.vals != null); frame.refreshContents()
}
viewAllBtn.addActionListener({
    if (steps.blind) endBlindChoosing()          // the all-frames picture holds the response frames: cells added from here on are not blind
    steps.stopBlind(); steps.showAllFrames(); steps.hideResponse()
    selectionBasis = 'All-frames structural review'
    showView('all', 'Choosing cells on the average of all frames.'); blindControls(false); describeRecording(); refresh(true)
} as java.awt.event.ActionListener)
viewBaselineBtn.addActionListener({
    int[] range = CalciumStartDialog.askBlindFrames(steps); if (range == null) return
    busyView(true); selectionNote.setText("<html>Making the picture of frames " + range[0] + "–" + range[1] + "…</html>")
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
            trace.show(null, "Traces are hidden while choosing cells blind", 0, 0); info.setText("<html>&nbsp;</html>")
            selectionBasis = steps.blindBasis()
            showView('baseline', blindNote()); blindControls(true); describeRecording(); refresh(true)
        }
    }
} as java.awt.event.ActionListener)
viewResponseBtn.addActionListener({
    boolean explained = respondersExplained; respondersExplained = false   // the question before measuring has already said it
    if (steps.blind && !testGui && !explained && javax.swing.JOptionPane.showConfirmDialog(frame, "<html><body style='width:320px'>The cells you have now stay recorded as <b>chosen blind</b>, as long as you do not change them.<br><br>The response is shown next and the plugin looks for the cells that only show when they respond. Cells you add from now on are recorded as <b>added after the response was shown</b>. Step 4 gives % responding for both groups.</body></html>", "Add cells that respond", javax.swing.JOptionPane.OK_CANCEL_OPTION) != javax.swing.JOptionPane.OK_OPTION) return
    boolean first = !steps.looked
    busyView(true); if (first) selectionNote.setText("<html>Looking for the response in all " + steps.n + " frames…</html>")
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
    def p = unusedProposals(cells).findAll { !hiddenKinds.contains('g' + Math.min(3, (int) it.grade)) && it.roi.contains(x, y) }.min { it.roi.getBounds().width * it.roi.getBounds().height }
    def target = null; def text = "&nbsp;"; def cap = ""
    if (p != null) { target = p.roi; text = evidenceText(p); cap = "this suggestion: gray levels above its background, over time" }
    else {
        int ci = -1; cells.eachWithIndex { c, i -> if (ci < 0 && c.contains(x, y)) ci = i }
        if (ci >= 0) { target = cells[ci]; text = "<b>Your cell " + (ci + 1) + "</b>"; cap = "cell " + (ci + 1) + ": gray levels above its background, over time" }
    }
    if (target == null || target.is(lastHover)) { if (target == null && lastHover != null) { lastHover = null; info.setText("<html>&nbsp;</html>") }; return text }
    lastHover = target
    info.setText("<html>" + text + "</html>")
    if (steps.blind) trace.show(null, "Traces are hidden while choosing cells blind", 0, 0)   // a trace would show the response
    else { reveal(); trace.show(steps.traceOf(target, cells), cap, steps.found ? steps.p0 : 0, steps.found ? steps.p1 : 0); if (!trace.isVisible()) { trace.setVisible(true); frame.refreshContents() } }
    return text
}
def findCellposePython = { ->
    def home = System.getProperty("user.home"); def cands = [Prefs.get("calcium.steps.cellpose_python", "")]
    try {                                    // a text file next to the plugin may hold the path (one line); up to 2026.11.22.1 it was in Calcium Imaging
        ["scripts/Plugins/cellpose_python.txt", "scripts/Plugins/Calcium Imaging/cellpose_python.txt"].each { String place ->
            def nf = new File(IJ.getDirectory("imagej"), place)
            if (nf.exists()) cands << nf.text.trim()
        }
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
    cpLabel.setText("<html><font color='" + MUTED + "'>" + (busyText != null ? busyText : found ? "Cellpose: found on this computer." : "Cellpose (the program that suggests outlines) is not installed on this computer yet.") + "</font></html>")
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
            py = g2.getNextString(); if (!new File(py).canExecute()) { IJ.error(CalciumRelease.NAME, "That file cannot be run:\n" + py); return "that python file cannot be run" }
            Prefs.set("calcium.steps.cellpose_python", py)
        } else return "cancelled"
    }
    if (py == null) return "Cellpose was not found"
    task.report(5, "Preparing the detection image", "Checking whether saved outlines can be reused.")
    boolean blindRun = comparison == null && steps.blind == true
    def dir = comparison != null ? (File) comparison.dir : new File(outDir, blindRun ? "automatic_outlines_blind" : "automatic_outlines"); dir.mkdirs()
    // Blind: Cellpose runs on two pictures, 0 to BLIND_CAP and 0 to BLIND_CAP_HIGH gray levels, each in its own folder, and the
    // outlines are combined (CalciumSteps.blindOutlinesToAdd). Otherwise one picture, in dir.
    def blindPicture = { int i ->
        double cap = i == 0 ? CalciumSteps.BLIND_CAP : i == 1 ? CalciumSteps.BLIND_CAP_HIGH : CalciumSteps.BLIND_CAP_MID; String lv = CalciumSteps.levels(cap)
        [label: "picture " + (i + 1) + " (0-" + lv + " gray levels)", part: "Blind picture " + (i + 1) + (i < 2 ? " of 2" : ", with a graphics chip") + " (0-" + lv + " gray levels)", shown: "0 to " + lv,
         dir: new File(dir, "picture_" + (i + 1) + "_shown_0_to_" + lv + "_gray_levels"), pic: steps.blindOutlinePicture(cap)] }
    def runs = blindRun ? (0..1).collect(blindPicture)
        : [[label: "", part: "", shown: "", dir: dir, pic: comparison != null ? (ImagePlus) comparison.picture : steps.outlinePicture()]]
    // Whether this Cellpose runs on a graphics chip: known from an earlier run, or found out during the first picture.
    def chip = { -> System.getProperty("calcium.steps.pretend_no_gpu") != null ? "CPU" : System.getProperty("calcium.steps.pretend_gpu") != null ? "GPU" : Prefs.get("calcium.steps.cellpose_gpu." + Integer.toHexString(py.hashCode()), "") }
    boolean thirdPlanned = blindRun && chip.call() == "GPU"
    if (thirdPlanned) runs << blindPicture(2)
    def pictureNoteOf = { -> blindRun ? steps.blindBasis() + ". Pictures for Cellpose: " + steps.blindRecipe() + (runs.size() > 2 ? "; a third picture, from 0 to " + CalciumSteps.levels(CalciumSteps.BLIND_CAP_MID) + " gray levels, adds outlines by the same rule against the outlines kept so far (run because Cellpose uses a graphics chip here)" : "") + "; no response channel. Folders: " + runs.collect { ((File) it.dir).getName() }.join(", ") + ".\n" : "" }
    String pictureNote = pictureNoteOf.call()
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
            else if (System.getProperty("calcium.steps.pretend_gpu") != null) known = "GPU"
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
            // the graphics chip did not work for the first two pictures (Cellpose fell back to the CPU): no third slow run
            if (i == 2 && task.device == "CPU") { runs.remove(2); break }
            if (runs.size() == 2) { task.partFrom = 42 * i; task.partShare = 0.5 }   // two runs: 0-42 % and 42-84 % of the bar
            else if (runs.size() == 3) { task.partFrom = 28 * i; task.partShare = 1d / 3 }   // three runs: 0-28, 28-56 and 56-84 %
            def done = detect(runs[i])
            if (done.error) return done.error
            results << done
        }
        // the first run of all on this computer finds out about the graphics chip only now: the third picture follows
        if (blindRun && runs.size() == 2 && task.device == "GPU" && chip.call() == "GPU") {
            runs << blindPicture(2); task.partFrom = 84; task.partShare = 0.01
            def done = detect(runs[2])
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
        def addedSecond = add.collect { int j -> [roi: second[j], from: 1] }
        outlines.addAll(addedSecond)
        int unshown = first.size() - shownFirst.size()
        mergeNote = "all " + first.size() + " outlines from " + runs[0].label + " kept; " + add.size() + " of " + second.size() + " from " + runs[1].label +
            " added (the other " + (second.size() - add.size()) + " share at least " + (int) Math.round(CalciumSteps.BLIND_MERGE_OVERLAP * 100) + " % of the smaller outline's area with an outline from picture 1" +
            (unshown > 0 ? "; " + unshown + " outlines of picture 1 that are not shown, for their size or as a strip on another outline, were not counted" : "") + ")"
        if (runs.size() > 2) {   // the third picture, against every outline that is kept and will be shown
            addedSecond.each { it.ev = gradeOne(it.roi) }
            def soFar = outlines.collect { it.roi }
            def shownSoFar = shownFirst + addedSecond.findAll { it.ev[0] < 4 && !(CalciumSteps.tooThin(it.roi) && CalciumSteps.alongAnother(it.roi, soFar)) }.collect { it.roi }
            def third = readOutlines((File) results[2].zip)
            def addThird = CalciumSteps.blindOutlinesToAdd(shownSoFar, third, ((ImagePlus) runs[0].pic).getWidth())
            addThird.each { int j -> outlines << [roi: third[j], from: 2] }
            mergeNote += "; " + addThird.size() + " of " + third.size() + " from " + runs[2].label + " added by the same rule against the outlines kept so far"
        }
        pictureNote = pictureNoteOf.call()
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
    if (blindRun) steps.say("Automatic outlines are blind: Cellpose got " + (runs.size() > 2 ? "three" : "two") + " pictures of frames " + steps.blindFirst + "-" + steps.blindLast + " only (" + steps.blindRecipe() + "), saved as " + ((File) results[0].zip).getName().replace("_rois.zip", ".tif") + " in " + dir.getName() + "/" + runs.collect { ((File) it.dir).getName() }.join(" and " + dir.getName() + "/") + "; " + mergeNote + ". The outlines were graded on the pre-stimulus picture without the response.")
    steps.say("Automatic outlines (suggestions only): Cellpose, command: " + cmdText + " (" + how + "). " + got.size() + " outlines shown: " + got.count { it.grade == 0 } + " respond on their own, " + got.count { it.grade == 1 } + " likely cells, " + got.count { it.grade == 2 } + " unsure, " + got.count { it.grade == 3 } + " very bright; " + hidden + " left out for their size; " + strips + " thin strips along another outline's edge left out, " + thinKept + " thin outlines kept as unsure (narrower than " + (int) CalciumSteps.THIN_WIDTH + " pixels on average, or roundness under " + CalciumSteps.THIN_ROUNDNESS + "); " + got.count { it.edge } + " lie within " + CalciumSteps.EDGE_MARGIN + " pixels of the picture's edge" + (edgeBox.isSelected() ? " and are not shown (box ticked)" : "") + ". Grades are measured (brightness, share of very bright pixels, contrast" + (blindRun ? "; the response is not used" : ", rise after the stimulus inside against around") + "). Only the outlines the user accepted became cells.")
    task.checkCancelled()
    javax.swing.SwingUtilities.invokeAndWait { proposals.clear(); proposals.addAll(got) }
    finishRecord()
    task.succeeded = true
    if (!blindRun) return got.size() + " suggestions; " + how
    boolean three = runs.size() > 2
    task.settings = (three ? "Three" : "Two") + " pictures of the frames before the response (0-" + CalciumSteps.levels(CalciumSteps.BLIND_CAP) + (three ? ", 0-" : " and 0-") + CalciumSteps.levels(CalciumSteps.BLIND_CAP_HIGH) + (three ? " and 0-" + CalciumSteps.levels(CalciumSteps.BLIND_CAP_MID) : "") + " gray levels)"
    def took = results.collect { (it.how =~ /^Cellpose took (\d+) s/) }.collect { it.find() ? it.group(1) : null }
    String shortHow = took.every { it != null } ? "Cellpose took " + (three ? took[0] + " s, " + took[1] + " s and " + took[2] : took.join(" s and ")) + " s" + (task.device ? " (" + task.device + ")" : "") :
        results.every { it.how.startsWith("the outlines made earlier") } ? "the outlines made earlier for this recording were used again" :
        // some pictures reused and some new: say so in a few words (the full account is in the Log)
        (results.count { it.how.startsWith("the outlines made earlier") } + " of " + results.size() + " pictures used outlines made earlier; the other" + (took.count { it != null } == 1 ? "" : "s") + " took " + took.findAll { it != null }.join(" s and ") + " s" + (task.device ? " (" + task.device + ")" : ""))
    return got.size() + " suggestions from " + (three ? "three" : "two") + " pictures; " + shortHow
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
            info.setText("<html>" + (message ==~ /\d+ suggestions.*/ ? "&nbsp;" : "Suggestions: " + CalciumOutlineProgress.escape(message)) + "</html>")
            if (!cellposeBusy) showCellposeState(null)
            refresh(true)
        },
        { -> frame.refreshContents() })
    try { refresh(true) } catch (Throwable ignore) { }       // only so the Next line says to wait; the worker is already started
} as java.awt.event.ActionListener)
acceptAllBtn.addActionListener({ e ->
    int nL = unusedProposals(cellRois()).count { it.grade <= 1 } as int
    if (nL == 0) return
    int ans = javax.swing.JOptionPane.showConfirmDialog(frame, "Accept all " + nL + " outlines graded " + (steps.blind ? "" : "\"responds\" or ") + "\"likely a cell\"?\n\nIn tests about 7 in 10 of these were real cells, so about 3 in 10 are not.\nLook through them afterwards and remove the wrong ones.", "Accept all likely ones", javax.swing.JOptionPane.OK_CANCEL_OPTION)
    if (ans == javax.swing.JOptionPane.OK_OPTION) acceptAll()
} as java.awt.event.ActionListener)
def updateEditHint = { String message -> editHint.setText("<html>" + message + "</html>") }
def leaveRemoveMode = { -> removeMode.setSelected(false); updateEditHint('Click a white outline to select it, then Remove. Or turn on <b>Click cells to remove</b>.') }
ovalBtn.addActionListener({ e -> leaveRemoveMode(); IJ.setTool("oval") } as java.awt.event.ActionListener)
freeBtn.addActionListener({ e -> leaveRemoveMode(); IJ.setTool("freehand") } as java.awt.event.ActionListener)
addBtn.addActionListener({ e -> def r = steps.colourPic.getRoi()
    if (r != null && r.isArea() && edgeBox.isSelected() && steps.nearEdge(r)) { info.setText("<html><b>Not added:</b> this outline is within " + CalciumSteps.EDGE_MARGIN + " pixels of the picture's edge, where a cell's background can only come from one side. Untick \"Leave out cells at the picture's edge\" to add it anyway.</html>"); return }
    if (r != null && r.isArea()) { cellEdits.change('add drawn cell') { addCell(r) }; rm.deselect(); steps.colourPic.deleteRoi(); refresh(true) } else info.setText("<html>Draw around a cell on the colour picture first.</html>") } as java.awt.event.ActionListener)
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
    boolean added = acceptAt(x, y)   // (a suggestion hidden with its switch is not there to be clicked)
    refresh(true)
    if (added) updateEditHint('<b>Cell added:</b> its outline is now white with a black border. Undo last change reverses this.')
    added
}
edgeBox.addActionListener({ e -> refresh(true) } as java.awt.event.ActionListener)
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
    // The menu command whose file is still there: one moved as the plugin started stays in the menu until Fiji restarts.
    String fiji = IJ.getDirectory('imagej'); def commands = ij.Menus.getCommands()
    String place = fiji == null ? null : CalciumRelease.PLACES.find { commands?.containsKey(CalciumRelease.command(it)) && new File(fiji, it).isFile() }
    if (place) Thread.start('Calcium plugin restart') { IJ.run(CalciumRelease.command(place)) }
    else IJ.showMessage(CalciumRelease.NAME, 'To choose another recording, ' + (placed?.moved ? 'restart Fiji, then ' : '') + 'start the plugin again:\nPlugins > ' +
        (placed?.target ? CalciumRelease.command(placed.target.path) : CalciumRelease.NAME))
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
        ' The next recording starts with whatever outlines it holds, so clear them unless they belong to the next recording too.</html>',
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
    if (cells.length == 0) { result.setText("<html>There are no cells yet. Choose cells in step 2 first.</html>"); return }
    go.setEnabled(false); go.setText("Working..."); analyzeBtn.setEnabled(false); analysisMeasurement = null
    Thread.start {
        try {
            outDir.mkdirs()
            def names = cells.collect { CalciumMeasurement.idFor(it) } as String[]
            ImagePlus corrected = null
            if (makeMovie.isSelected()) {
                result.setText("<html>Making the background-corrected movie...</html>")
                corrected = steps.correctedMovie(cells)
                if (saveMovie.isSelected()) IJ.saveAsTiff(corrected, outDir.getPath() + File.separator + base + "_background_corrected_movie.tif")
            }
            result.setText("<html>Measuring " + cells.length + " cells on " + steps.n + " frames...</html>")
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
                def check = runSensitivity(cells, outDir.getPath() + File.separator + base + "_", { String s -> javax.swing.SwingUtilities.invokeLater { result.setText("<html>" + CalciumOutlineProgress.escape(s) + "</html>") } })
                def m = check.m
                sensitivityHtml = m != null ? "Sensitivity check: the blind detection also found <b>" + m.found + " of " + m.cellLike + "</b> response-assisted outlines that could be cells (" + CalciumSensitivity.pct((int) m.found, (int) m.cellLike) + "); " + (m.cellLike - m.found) + " missed" + ((m.cellLike - m.found) > 0 ? ": <b>Add cells that respond…</b> below adds them" : "") + ". Details in what_was_done.txt.<br>" :
                    "<font color='#b05000'>" + CalciumOutlineProgress.escape((String) check.text) + "</font><br>"
                IJ.log((String) check.text)
            }
            def html = "<b>Done.</b> " + cells.length + " cells measured on " + steps.n + " frames.<br>" + sensitivityHtml +
                "Legacy 5-SD screening rise: <b>" + steps.lastClear + " of " + cells.length + "</b>. Protocol responder counts are in Step 4.<br>" +
                (edgeCells.length > 0 ? edgeCells.length + " outlines at the picture's edge were left out (dashed).<br>" : "") +
                (steps.lastDim > 0 ? "<font color='#b05000'><b>" + steps.lastDim + " cells are not brighter than their background at rest.</b> Review their quality status in Step 4.</font><br>" : "") +
                "<b>Continue to Step 4 to analyze responses inside Fiji.</b> It uses unlifted fluorescence and reports excluded and unscorable cells explicitly.<br>Saved in: " + savedIn
            javax.swing.SwingUtilities.invokeLater { analysisMeasurement = measured; analyzeBtn.setEnabled(true); result.setText("<html>" + html + "</html>"); go.setText("Measure again"); go.setEnabled(true); s3.setSummary(cells.length + " cells measured"); s4.setExpanded(true); refresh(true); openBtn.setVisible(true); frame.refreshContents() }
        } catch (Throwable t) { go.setEnabled(true); go.setText("Measure and save"); result.setText("<html>Could not measure: " + t + "</html>"); IJ.log(CalciumRelease.NAME + ": " + t) }
    }
} as java.awt.event.ActionListener)
if (testGui) { binding.setVariable("hooks", [nextLine: nextLine, addRespondersBtn: addRespondersBtn, addMissedBtn: addMissedBtn, continueBtn: continueBtn, acceptAllBtn: acceptAllBtn, autoSuggestAsked: { -> autoSuggestAsked }, legendClicked: { String k, boolean o -> legendClicked.call(k, o) }, hiddenKinds: hiddenKinds, cellRois: cellRois, sensBox: sensBox, basisFor: basisFor, responseSeen: { -> responseSeen }, reveal: { -> reveal() }, runSensitivity: runSensitivity, trace: trace, respondsLegend: respondsLegend, viewAllBtn: viewAllBtn, viewBaselineBtn: viewBaselineBtn, viewResponseBtn: viewResponseBtn, selectionNote: selectionNote, selectionBasis: { -> selectionBasis }, changeBtn: changeBtn, changeRecording: changeRecording, continueBtn: continueBtn, backToCellsBtn: backToCellsBtn, backToMeasureBtn: backToMeasureBtn, cards: [s1, s2, s3, s4], recordingNote: recordingNote, detectionProgress: detectionProgress, suggestBtn: suggestBtn, edgeBox: edgeBox, addBtn: addBtn, setForceFresh: { boolean v -> forceFresh = v }, setupCellpose: setupCellpose, findCellposePython: findCellposePython, setupBtn: setupBtn, setupRow: setupRow, cpLabel: cpLabel, showCellposeState: showCellposeState, setupText: SETUP_TEXT, saveMovie: saveMovie, makeMovie: makeMovie, rm: rm, refresh: refresh, go: go, status: status, steps: steps, outDir: outDir, viewSwitch: viewSwitch, hiddenLine: hiddenLine, showSwitches: legendEntries, moreBtn: moreBtn, recordingName: recordingName, frame: frame, suggest: suggest, acceptAt: acceptAt, acceptAll: acceptAll, hoverAt: hoverAt, proposals: proposals, result: result, info: info, removeSelected: removeSelected, undoLast: undoLast, root: root]); refresh(true); return }
frame.setVisible(true)
timer.start()
refresh(true)

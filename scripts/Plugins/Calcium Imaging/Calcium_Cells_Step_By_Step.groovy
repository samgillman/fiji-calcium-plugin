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
        double pixel = CalciumFrameView.micronsPerPixel(movie)
        if (pixel > 0) {   // a 100 micrometre bar at the right
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
    /** The width of one pixel in micrometres, as the picture's own calibration gives it; NaN when it is not in micrometres. For scale bars. */
    static double micronsPerPixel(ImagePlus imp) {
        if (imp == null) return Double.NaN
        double pixel = imp.getCalibration().pixelWidth; String unit = imp.getCalibration().getUnit()
        return pixel > 0 && unit != null && (unit.startsWith("micro") || unit == "um" || unit == "µm" || unit == "\u00B5m") ? pixel : Double.NaN
    }
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
    /** A cell whose ΔF/F₀ falls below −1: its corrected fluorescence goes below zero, the mark of an F₀ near zero. */
    static boolean belowZero(double[] trace){for(double v:trace)if(Double.isFinite(v)&&v<-1d)return true;return false}
    /** Colours run from −limit to +limit. The limit covers 99 % of the values of the cells that are not below zero (as on the
     *  Step 4 page; an excluded cell has no values): a few such cells would otherwise set it and wash out all the others. */
    double limit() {
        double[][] steady=traces.findAll{double[] trace->!belowZero(trace)&&trace.any{Double.isFinite(it)}} as double[][]
        double most=CalciumFrameView.share(steady.length>0?steady:traces,0.99d,true)
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
        catch(groovy.json.JsonException bad) { throw new IOException('This protocol file is not valid JSON. Choose a protocol file.') }
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
        super('Step 4 · Response protocols');measurement=new LinkedHashMap(data);loading=true
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
        def tables=['Whole-recording metrics':metrics,'Cell quality':qc,'ROI registry':registry,'Cell response calls':calls,'Responder counts':counts,'Stimulus-window metrics':windowMetrics,'Biological-unit summaries':replicates,'Sensitivity (diagnostic)':sensitivity,'Negative controls (diagnostic)':diagnostics]
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
        responding.show(Double.isFinite(percent)?String.format(Locale.US,'%.1f %%',percent):'–',cell('Stimulus')+': '+whole('Responded')+' of '+whole('Evaluable')+' evaluable cells')
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
    /** Rows of the heatmap: as measured, or responders first with the strongest on top (then non-responders, then cells without a call).
     *  Inside each group a cell below zero comes last, as on the Step 4 page: its rise is divided by an F₀ near zero and ranks nothing. */
    void orderHeatmap() {
        List<String> names=(measurement.names as List<String>).collect{it}
        if(heatmapOrder.getSelectedIndex()==0&&!callRows.isEmpty()) {
            def rank={String name->String call=callOf(name);call=='Responded'?0:call=='Did not respond'?1:2}
            def below={String name->traces[name]!=null&&CalciumAnalysisHeatmap.belowZero(traces[name])?1:0}
            def peak={String name->double v=CalciumCsv.number(callRows[name]?.Peak_increase);Double.isFinite(v)?v:Double.NEGATIVE_INFINITY}
            names=names.sort(false){a,b->rank(a)<=>rank(b)?:below(a)<=>below(b)?:peak(b)<=>peak(a)?:(measurement.names as List).indexOf(a)<=>(measurement.names as List).indexOf(b)}
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
    /** How the page was last arranged (pane sizes, explanations on or off, how the cells are drawn on the picture, how often the
     *  hint about zooming was shown). The page's address changes every time, so the browser cannot remember these itself; they
     *  are kept with Fiji's preferences. */
    static final Map PAGE_PREFS=['step4.side':/\d{2,4}px/,'step4.split':/\d{1,2}(\.\d)?%/,'step4.height':/\d{3,4}px/,'step4.help':/[01]/,
        'step4.every':/\d{1,4}/,'step4.outline':/\d(\.\d)?/,'step4.plain':/[01]/,'step4.numbers':/[01]/,'step4.viewhint':/\d/,'step4.format':/(png|tiff)/,'step4.size':/(screen|\d{1,2}(\.\d)?x\d{1,2}(\.\d)?)/,'step4.dpi':/(150|300|600)/,'step4.compression':/(lzw|zip|none)/,
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
    /** The width of one pixel in micrometres, from the recording (or, without it, the picture the cells were chosen on); NaN when neither says. */
    double micronsPerPixel(){double pixel=CalciumFrameView.micronsPerPixel(movie);Double.isNaN(pixel)&&measurement.picture instanceof ImagePlus?CalciumFrameView.micronsPerPixel(measurement.picture as ImagePlus):pixel}
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
         micronsPerPixel:micronsPerPixel(),   // for the scale bar of a saved picture; null (json) when the recording does not say
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
    static final List FRAME_MODES=['raw','flat','green','magenta','change']
    /** What every frame of the recording is painted from, made once for the frames the baseline is taken from (rest==null asks
     *  for it again): the resting picture, the display ranges, the smooth background and the rise. */
    private void prepareFrames() {
        if(movie==null)throw new FileNotFoundException('No recording is attached to this measurement.')
        if(rest!=null)return
        def stack=movie.getStack();int n=measurement.frames as int
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
    /** The gray levels each view of the recording is shown between, darkest to brightest. One range serves every frame of a
     *  view, so two frames can be compared by eye; the page writes it into the record of a saved pair of frames. raw: the
     *  recording's own levels. flat, green: levels above the smooth background. change: the rise above the baseline.
     *  magenta: both, the cells in green and the rise in magenta. Display only. */
    synchronized List displayRange(String mode) {
        if(!(mode in FRAME_MODES))throw new IllegalArgumentException('Frame mode must be raw, flat, green, magenta or change.')
        prepareFrames()
        Map recorded=[shown:'gray levels',from:rawRange[0],to:rawRange[1]],cells=[shown:'gray levels above the background',from:flatRange[0],to:flatRange[1]],rise=[shown:'rise above the baseline, gray levels',from:riseView.start,to:riseView.full]
        if(mode=='raw')return [recorded]
        if(mode=='change')return [rise]
        return mode=='magenta' ? [cells,rise] : [cells]
    }
    private java.awt.image.BufferedImage paintFrame(int wanted,String mode) {
        if(movie==null)throw new FileNotFoundException('No recording is attached to this measurement.')
        int n=measurement.frames as int
        if(wanted<1||wanted>n)throw new IllegalArgumentException('Frame must be between 1 and '+n+'.')
        if(!(mode in FRAME_MODES))throw new IllegalArgumentException('Frame mode must be raw, flat, green, magenta or change.')
        prepareFrames()
        def stack=movie.getStack();int f=wanted
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
                    case 'api/range':replyJson(exchange,200,[mode:query.mode,ranges:now.displayRange(query.mode as String)]);return   // what a view's frames are shown between
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
    static final String DATA = "UEsDBBQAAAAIAAAAIVwauiPzMQEAAMkBAAAOAAAAcHJvdmVuYW5jZS50c3ZNkDtu3EAQROPlYbTT8+npcSooMGAosA+w6N9oCSxFgqQC3d4TSICyDqpfPdQxL9vDlW+7b+sxn+v+ebmf53b8ul7f5vP+IU+6LteDl7f58Vj4/Xp8ffC2Td/3bWSW+bxIi5bVG5mxpK7OmkqE2mpICj1TDzlrC9Of388vr/9ebsedY8FLza1aipG6xmxBO0msUgsktY5NJEU2KVm7sCByrdnEQyg+GjLwtPuxre/m+/H09xsaGQCBCM178ky1cRY2KqacrWoxwtJiw9ighKjaqhfrUNURUu7Txzk/fvDcmiTjYSpEQ7GaGwgVjJWCdAMjj90IgAO6qY4G4FJUeuUA+CV5+G3b13PV9SccqzUAbzh0onNr2A2lMpliGXDTQDJAGXpiqpEwYRizOtBYqcL0H1BLAwQUAAAACAAAACFc+MYmeYYQAAAdLgAACAAAAHJ1bm5lci5SrVrrc9y2Ef9+fwU0UYZkQ52ltJl2KrMziix13MqyK8XpB1XDgUicDjGPpAhSjzT+3/vbBcDH3cmPtpmRc1xggcW+d4FvxKn+RQuZy7pVzVxcZlqVrV7oTGSyyLpCtroqjciqlRKLplqJdqlErctS5eJSr+pCHcvvvxOm6ppMmflMNrdGvNwjhJUs8yN8hm0jdaHL27dl8SQS8dPF+5NophciLFR52y5DwonETiL+cBYJ01Z1GJw81iprsQdm6FKJXDf4rpqnWOiy7lpxfPlzLIxqW6xrxE/0VXUtjfRTg2hmyQoXulDzWrqtrg6uYxF0tWkbJVcB/W51YeYXQfTFGI0yNfii0rqp2iqr/iv0XDUezR8EnMOcfJ6rQq8s8u+BvARMNeDd6dHZ5UkM7hbzUq6UASgLgw/qiVa9l0Wngiieie3/Aeu4kMYwWpAtZSPB4waod13VKgIGMUsOOjCncQ9aquxDvyHTYAUoy6dXXV3oTEJWoT/FLuiJvCT78V5aJBc6KL7PaclQmnmhTTvg80GiQb684KzsVoS26MqMlDIkoPi3uOe1uroGU80/ZVMSBq2J+arRWWiurjD1+jqKDsVI6+6jneTgTPz2m9jRZr7QpW4VgI7uWhp8Bq9L0KLzIOYjHWKzj15Y67TEIlcL2RVtEgQRb6SJhqIYCKDNCChHNHksoQoDHjn4DHp1Swea7FLHolWPLZ0aHA1N3eiyXYTB5es3785Ojo/+1X6b48/8q4TQ/GRQvSg6s5xnUNiqUGGEM8BI5hk0DWdkJfvDdQyr6Rqj71VCBhqbZfXguZk4kT804NGEoseYtCISPDLPzD0ga8qPpXlO3FQPVofscoASp2aZLKsSGlKk+CumB36MxONVDIbnerEIGRmwGCr/T13m1UNaK/khLSrCzk9Pg9jDj94fD+AUKhfFeVPVdufrmR2ynmRTpWA/rIFOSUnYPAI4iXWn/JVsI7RaSjK9kUYRcQ60k7ihsSb1Hs2UsgZzW/pNrCJEsagap2PYhHjHZG0w0q7r6eBx9QjTMSFNW9PdN9oYkD5dd/axV63whzi4gGMhXV5AT/dgQZk24INYKWm6hggsugpzM1Vm8CyzvO09FMmayfr+Oh75h162c3g71h0IC3I5P6J/5Dnckxs4MqeS/HSvXXSmEj4qzNvo5fdsmSV0hj8P9u03awAAcKk7SXDa4Duggakj6qdFG/GEUTA9F12p7yDjTBUFecZuhTD3oMF52YoC52/FwT5CHi0VWOJ2dE7BEaozdi95ewVqrqN47HLUXQrJhp7+aKCECTBi1WGDG2xPUSRDALpX8HetbNi1gISDoHex4c69rOviibbaQywZm4ceNn2MfvsNFhRaD/MIlS+qW6b2YLT/eVXuOYyJdEXViFaDtFauajrxQjeGxY3ZYUAqjhiuUqYxgFMp5LZhVeaOcMbfSdoG1NoP0EdIHka/AeKhlwd27CV/2d9/SXr2OdqdNxZl1azw41dOToTfWzSyvFVud2+3ARt5uqpyjEBlbmT24bapujJPs6ppWCnSqZa7zThb6cpCL0hvBry9Hm/KP22EzDKFPMpzIPgJ/ITz+VaX3w4qyf6FOE3My9tdO8nZNEt7CEc0LSKxAky/X+67D/aHPPoy2e/5Y5eyugXFaqUuhV0pBs+Qs91Kr2gQf4vz6ZKiALkJeDqg5CR5S15qcxmiMlCPZFu6TUcaMvtoQxYdp5E2LrAuGAxTupcuf4We6EVIo0TmyAAshjeCujKaCOv3prXArF2eu3cQvaDph+tUke20smzTBZMzI1P2mcUoYIDpFDGss4i9UCJn0y4bYNyRkcDJ3cKxKnKCYxcBEu1I80Q7nb8/O7OOy8aEXu2aSqd+YkqOOoic6MfovScdBaStuFtc7HNuVIhSqRxkUxYeBscgHqd+XWZFByj2SJ13pywPIyePGCGvn5ImVEgdgh8HI6n1I6Qc9DGHHIxdf9BqT6vT1JEj9iO7REUUTa3YW8vexdvXPVd4px7tGarpaG1VVA+qIU/cZ7Lh5xCjiKgGV+CBSBcO8PekRseD/NRdB5c5ofzqc+siBFr9sSccH0jkFdx9WbViJdtsOYTWNa36OPNezYoul62cc/gJSWMTNvYZ0aZb1p2VIiNmfb/ngicMKIsO4qlGz5BahBqGLigoyaKCdnplt+poKPsk4/X7p0bfllBH1G2oZdZCHeNe6WuKeM6vx+zJ4+A5h4pd7qy201a77gzMcnZ2nCQPSsRkgTCLYqFXE3kkiacinqR1hObWdA602d3QZSjpFij5J7ex58kun99b+eHd7mUr244Zbo0GehAAfMFmw+AzMuNRrED8kPcofuVNwTkeqmeECahH4Ha62/3HMbglfennysL0LrPQcJgAszw6/vtfL96+P3+Vvj8/+vno9dnRj2cnzF+BnN7+4fwDS4djgK0IprSJWtXtU+q0J6UIGw2OXAxakI9lPRIfL2clZBf0jQKV0ny/8mT3mLR3TV1IwVNtKEHfT2xPwJ6AdSW7QSIfsrB7Wd99wud9I95hQ9Xcg89ypQuN4rVwyW2f3S3hS8UHpWqCUiOD0l5h4D5gc960oMJVk889IXdX8eBFT/fZO7pc5/JV2sgHQKxiBJxOswedBqC76OuXiNYqFuHpuyJ5MGGAOV73MMu1u/gNwylXw6KJVQlvEdXD7juqmojvURQcef0MrBYE70sDBjCESjYuHHzFcLD/3aKo4FB++OF3+sXU0cR9OXruVIh5T17u21xUC/wbxHrNO5Hfu8vY31Vzqv/Cho4Qu8NGhyMvN53hBqKhnPnjfhz8TJFFtl6+ucoKFq7v1wjfr+EUbSNq+9FpxPbQqYnqVV01bd8Amobw9YXGydI34oRqNptz2QbHvDX3XO4ZbMBpG1EPXwElpwxecE6FrKnDP+pRUu42f5ayniSzi21WHcJ7bDMz9+WT+DUo5+62e8QDfZNrSP59jmgxOdFr4QdT1goPLisNJKzZaoyrZrLoqOBJ7aloOfKwMaWsTVUkQ5lldj9kBegibcxjG7XdYm6yPQqr+KexIt6dBtzhWR6872Q9MOGrVmOm9WuRPuvyHlU0ZoeDSJBAT5TF/9jSLOl1Z/6LcwQzqv60TQs2pZwirkDNVJEPqw5pKDsz6gJQvvsFhLmtttHlhtJN+oZ25pREpiq1gcCEQ2Sx2dJWKuP+DBxMnB9I7jK21y2nv8eu3BtzgChJDs5sEumJ2u0bSwN+gVooe0pliczM9rj/R/pit6I9bTLsHs2+Ece+x+U8tm1kwD2ZrNE31FmXXNFjAhwSqlhOFhTO9sStTYgI4ofbo7y6auazbYdba6SF63OiDazdY3xpLqJHrS4yUQtF+baB8qOuXFsh7RBQJog307FguiMSovZrt7Q4X7cn9yJDFx/iwP2gMguL28GRqIPht014pzMh2GAsaTvo9lhnTjw09q0MJkttnGsynSF+ddvE9HmLLeROyZiQILg84qKiAI28wTrvbVlIcFLoW82BPOjZTtnHGstG+Qj7xtSMQXBwqRntlJIDbx3Izz5xkzizcE0F2vZtPxcEIPtbsa0x4KIrR9vRl9+JfuddszbVxQo6ngtLlHOWMLl+4zFgs0HsjWzaJfbQMRLBB9JGQE/hCJRV99Qk2AqWt2q8K8GhjcuKxOjT+meV6SpGNoUvVFThVCFiX3av2fckY4wDhzTR/xkC9cNaecmKlWx1dkeuRE9Rsvg8cq1am+RzHAe5a9dPoOj6jqrkvlhGzLHF0/YF98/cVQgi7M5nC/r4ou/TcCCU8KFPRptkQhcRcT5uGqbcgkhoj7vMFXY7yVDWISOZTO9JnmAkE4xLVSirSsiXQEHPUTMdoEZ41bo6Q7mUCxCVBBejphM7BIHsddoyMIfCJ2x7plYZ39UquhiQN5oLGFSYSCCfcy2sa6QFVK8vKfbThxuFrkjSjC/XkAvlmj3J8x45/iuK4TrZ4uxjqjbsIfjIfRD4xCVBf//iid3l5aNo+s31+PbVhwjB8+PhQKmHOX407nQuXA089QPwmN1qJZunLd7dr+WIpuxqtN5uvwQMdzt8TQqWi0Oza8z8MdQ76gnwXVW7y/wJmCULy4anuXVL2mPkI9iJ1a5CjWDvVJORLbtjg3eJ7xJFzxyYWiL12mEt7LMH/Zoj9UxBQbOF8DdKliPx1fYcA/ExleSfGKbwm9YVom++OTr48m3yJEXzC4+N8nATyfIlDvj/2yx4wzrH2F0pDfXwyCeygGObEtAkbm4jMqgawYx6qYO57nFy9foVp5yU6ve9xE8uTq2KYARds39Pd1+n/+l7ZA7ynu+5XGHjSgS+lxvTlmt5W1bwd5nZXrAbVfLVAfzetGa/bXT+fId9E+2rGuwDWVOvMF72v6x6iO6Yqh1h79TD9ax+2NtF/Xhyms20c4zg0s4Jwli3Pm7jsb82Sn2JPGG0zTLM87x+Bv3/wnC/5EhPvpbhjv6v5fnGqb6I8ZtYU+6jUqR+2eXJG4ritjysCWEuflwzUxt3tLIpglG15Cs1sqCSCkfY0QL/sgXal1vuwgEbzme20iFmSjNHUt7oxxHj6JY5OrQNQ5ozupBkvOjQhsPq4bJbmdCORTN5z49iSIXI97u5EaaR0zXsRehGwQ7EO+5q8vyIbp6KNKJXA83K9oI5BSR35Ue9K7JXojOjVtv2svfkboOD8U25ycPHq+Egj9F19MLcNW1IGd0Yjv8+t/mGC+5vaDi8JGBEDBki41nF5ymJ1CRlRKUoYg8TN/jEe1XK0l4zjF5+bT4gC4a51C2ECttHYc52+hdhiX8P1j8HG559JeNHX/zmi1/gqEfQZKvMC99PSSb3areqvfCNFnh7avwg1VTbJ9cy+4DK52c3P/Cz6RHOqGeZGniAldy+xpc0egYjTpfSLJMtOAQfEOL1Dtan8dzstW1SWdxWkP9ylQRvXv1Aje3eXQz92L9dvj0XronYPsV8/yfpir1rKEXnC7YWOX0UxNPr9mT0O568c0iC8HTvdD96cbp/KBb6kR6z7FNtUBINtsbghwruNYR7uhKP+hvsDJJgKEmHm6okOCpI/6jfdNOyFHJRkWLePPEb0UOcwb0R6KfwGpRaUTh5SumpREIPTBSfispdLmGhiElwxu0xAv5ZPCwramZ7dxwLmznAcbWCM7VGFkHMF/VtldrbCsUGynxwS1Fp/+dhleH1imW3r5KE6wowTdOaOHEXdvR00M9uFNwm3VseiuqGL5Ry133nskrJbLnlWsGuyzxia6orm5nypij18a2zpyTwL+pcm1/ljCK7tqLOSCZQikHGe221d/mKXuTQGfeqrqU2ulD+pj6IoVjVSpfUDEyCPmcXL0Q4fHwnXvFbmdYHzOhQtCpbWlXtFzPuNZS/9XEcNX1M4TOMMsCEW1A27pRgHCkKc8r2UQ83MsohVDXqrtMUrlyCiWQTGjqK3lRQ1wUIoQe8YpSnuMBmw+eeC5/uyS8iZeXOyQFmPa2Ei0NwTiePbPmizEPiwbfGI68ON5rYVgvLNLLP75Ku5Ld+FhZtrj3+vtrJh9cQ4wF60hpTWD6TprUhb3qXbWPMGCUOfHMjXclSL8D6cUlwBnUyVFyTrpLUtrTwj9/SQ9GfTug2YHSptx8HR25p0eNHs/8AUEsDBBQAAAAIAAAAIVx6uMRcdQIAACwEAAAQAAAAdXBzdHJlYW0vTElDRU5TRV1SS2/bMAy+61cQObWA0Q0FdtlNtZVEmF+QlWY5KrYSa3CswFIW9N+PdNJ2GxDAEMnvRaaQGnLX2jFYxlJ/fpvcsY/w0D7C89fnb9CYE6zcMJzMyFhtp5MLwfkRXIDeTnb/BsfJjNF2CRwma8EfoO3NdLQJRA9mfIOznQIC/D4aN7rxCAZa1GE4GXukCf4Qr2ayONyBCcG3ziAfdL69nOwYTSS9gxtsgIfYW1g0d8TicRbprBmYG4F67y24utj7S4TJhji5ljgScGM7XDry8N4e3MndFQg+hw8MSS8BE5DPBE6+cwf62jnW+bIfXOgT6BxR7y8Ri4GK8xYTyvHFTxDsMDBkcOh7zvrpbp4h62daaLyvKFDl2vvTv0lcYIfLNKKknTGdx5XNir9sG6lC4wc/DP5K0Vo/do4She+MaWyZvf9t5yy3244+otWbBTrA+fOq91bozTDA3t4Xhrq4XvNXnInkQ8TDOzPA2U+z3v8xn1B/LaCplnrLlQDZQK2qV5mJDBa8wfciga3U62qjAScUL/UOqiXwcgc/ZJklIH7WSjQNVIrJos6lwJos03yTyXIFL4grK/wDy0JqJNUVkOCdSoqGyAqh0jU++YvMpd4lbCl1SZzLSgGHmist003OFdQbVVeNQPkMaUtZLhWqiEKU+glVsQbiFR/QrHmekxTjG3SvyB+kVb1TcrXWsK7yTGDxRaAz/pKLmxSGSnMuiwQyXvCVmFEVsihGYzd3sF0LKpEex1+qZVVSjLQqtcJngimV/oBuZSMS4Eo2tJClqoqE0ToRUc0kiCvFjYVWDf9cBEfovWnEByFkgufI1RCYIr4PP7E/UEsDBBQAAAAIAAAAIVy59xCppSUAAEGVAAAVAAAAdXBzdHJlYW0vcmVzcG9uZGVycy5S5T1rc9zIcd/5K0Ac7SzOS0qUK67EOqqKJ1FnlleUTErnSjbMGtwFSeR2gRWAFUnb99/Tr5npAQYgKclxnFzV3S1nel49Pd093T2Nb6J3WbU7z5bLaJ4ul3W0qbOouc6iOl1lUVFWq3SZ/zlbRDf5AiqqdJ7VUVpH73OoflluKgBPi0X0uyxtVul6b6vK6nVZLLJqtq7KppyXy+i73ehyU8ybvCxGdZOvNstNPY7qJq2acZQVi3F0kdbZMi+ycXS5LMtqDAPnNf5VwSwA9ipdj7ei4X/mZdFUMNpB9PpwcnY0NgUzGgiKTz5MJq4UxpWyJPoL9F1sVhdZ5c31Nonyeg8qsiqf41+//GW0zIqr5hr/ODiI9idYBDCXeZE3GZRCR2aF2FVT5aub2i4aq/PLaCSd2OJom/r6619pvFRVQNF28ef5dVq5woQwUTflehQfFQ1MOnWDFoCwvdgMtA1bOvqUrtfLu9Eyr5uRQ3oyliWPo2V5lcPmj/aThEZE1Hx3wBvUHYyIwwzHuEUCwEZ5EdUZIHhRw86VKwKtoKBa5MVVNF+W85/2ovdQiMCrTd1El+VyWd5In9DVwNSHaaR/OUI3pjmu7Cmuknqxf1Fn0XdSRV3CX/vqL9ijpgLSGMmAruvtEZEJzBfolDA44rnw30Ai8P/oRfQ08fbu5XVZ4umJ1mWdN/mnzM2RCWQcpav1Mm82cPJ4tojoVXoLyF/RmYDGRVkU2VVK7c9eAVaXTQ6NEBEInUY31+UyM+RdXuIBgE3aUIMaB8hqxvpNXizKG6LbRdqke7TO0Sm2PojiM9nyeByZn9FB5I7zGW7fjMuIxo6KBf2JtCZ7mtfvTz8cjeQIJnzwuKopgQ4yReV4vuKf5ss4EXR9AFz9/uUSmY9HgumyRKIA7GRyGBb55WVWwV/RGs5fWt1ZYF5oD4F5/MJjFANnRbMToCWvj0Towz+quIa+Y/Og08HT5y52PbrWkwGq8xmgDyq74qZ3lq3TKm2yFnKhJ5ywoY26jLJ0fh1dwy7kTR2VNwUgObPwZggzU0VT1QX8HknBOERiL3m+PoXFMHx8L/NX/zhCbG2oIUi9sTjLn+FfEGxFfQnSzk3we1kKtXHM5zUexBmwsnQ5g+mTpAAA4Ulnr2buBEIx8yhv9i/d+ZvJ+cP2IunepLczONk0KB8Zx1mS6ORwVmUwcJQt4TRg2dbPW0rk8uRnLMm1KFsgQYOQhxEu82wJIpdB+Qg2KMzx2DfTaYySPT4/J0mGdAPlDLsjmIUaonZbToiF0mV2SeCG3qTaIRJgvomOi/kSORoOCpCrNZ7ntJFTy4IjjS7KTbHAo1sWy7vo4g7Y3vwaabfCChAnezhvVi7+Jfo22nvD9TuLcnOxzPYy6PZb5JWjfWCEF/VoPsLpabWDubGVTzNZVJPmyxFObgwyGSHKT1kFjB/rEGTXbwIA82W+XmfU2MC+oLnBGbSNDww63A7/CoEs4935o2ye601+dkBu0qpAiQogSCMCZthpDVyvaC5HMUPTpE2nhN1fXMHA+F9k5qR/eSuiM8EUFsfCue0QgqQ2CoAHEEFc5/NrQl704oDpYZcxITRmaQMKcZzFppKF+E0VmG1LTIuRRgJlcUs4GsHYY+mINzSty8LgBptO989hQ2g21Br5tt7lBObFO4tdt7F4ahWYRQlHtSgbqEsL4pOXcDQdW0VmyWcRUMc7BppAloE65FD+XGtESB26MmlxOiZZ2iJZSP9eKcUSMIJr+vXELCZ+nd2Q6kbzrrIsyLSNOhAHOxUMQ79CjAFGZsd7X5bRZXbj0GE4HfBYwpzSQna5Lqo2y6w1NJ5f5Av29BIWYOfP5QRLCc8N0HMeCftJYMfDB87O8YyHRd0su53TNtGeIn+Mlvkqb9RspE0s8sJwPmBZc0JARsJFiM5xgb4p+Ix8lRej8DmU9q9xTjNH2fwLedOmaelrBIoChdn8S7zZHRjez92RstawdP1Q1DBkCiyzR8ie8qAHMqaSi6ssxXKzlHBzCw26qYIFPQdIAzTTxSwv5th15lW/y9KfgjXBQVgmN6bLyPLakLAeRycolD3w4e7dXGu48M5lPngpgx5B6qLUrq5yDxnRW6BvVkOGO6eVemDR6aZgnSVQnokSM9xpVd7sFXRrsdfe6wwuX6aMbshIX3gcRzlroB9nqEpfjYhWrHJ+ZxQDKp7m56QaRJ5ycmehiSJ3mLwAFtvGR7ck8BexgmCqMhDHxad0mQPzTG9AkdqUgO55VszVhWj0+mli2hfZbUM/f7YTkSsykyhMxoJgZVrcjbbdJf2O+EfSmrI/oTd5XSODhmltAGGAHjsTQFiLx/fMC1tgZ3hMYFBguOesEws5PsXKfWAICGgn6ypBCH6smlFAsRmc+gQ7cNP10JnX0Z+zqoTbjLkxIjcGVl2vQGPsWYfcjImvNvVvf1svzGqiJ3YtDlKdK6U5emrxtwxJbTxoZPeDJ7c1AuOzkvnBtFgOnMtdpzNB4Ml8FVlmI2oFmkars8u8qhunj2CLHaGCX1L7HRaHNbYdEIQgi6g/e7eHLrfZukNDuIucgXCbukx5CvVmNdJDTvGMwm/p4TzRROB4hlCCFozU4643/yl1cg70tz+RI+16Ii7T048HLFxOQLsTUXS0BkZnEbuHW41bwGvAhswIA4NiQxkU4ZjXChx2wQCw13YrLWyL00ujFhVZ6DafF3BiMC3SfjEkYZIoNmKPx4pFmybg2I7nsUmfFEDr5DvdwrR9BfwRdU+57LlevOPf6uVsA1sB01ioczavSmJupuMTOP5DYHYgT/DLeJZ3eQBnrzSmA/g1kt3sdmcrPPnf2WtRwgDQv/5eippUA6+57/5rrNNcPvu4AQnU3Fl5SSpcvrASlE/nN9Hvs2wdXed1U1ZoBxKj+cVdk+2CLN3FHyguQJHMK+CuOZAqQKGOhfJknTbXe1H0I2j9Oelv30TPkCuX62YXWuEFoshgZBDXaXGVsULa3JRaRWZL12JdwvWE73XYc3kZfdrfk9ta/E4WN/vEI8XRL/LiF2Qcrkdm5VaMCBArlm3rvenCNbOcy7Q7OIieTQwzqrJmUxWjnk2ZfXr24L1Qm2DtNAZy5wfQplflImOiR0oeMXe1EFbnhktCfLhpyhWIvDnceOM3aQEjsKUZjXLA39nqiTcNbQNJ9Ii2v2nvQESm+3v/DBJOhOUqW+Qg/u0wUzfgi+ipP8CpqBWAKzJfzmlx8WtjZIItJj7j3d20fMejak0/S7ZsGpFRgEqodn7sjkYuUjRoRBrYrGk+juDeXZVro1Ce0zbNy03ReIevInvaulxvlinLcWflFSUzq63AA+Adcz850IwwESWjD7LNIQUewbH/X0EllfRfloKG7Xd23lCsFtFWvE+zK2ALsHHYIaEbJwiYtiuA8ju08Z2UhcV2TUbCcXSESgZew/DvTt9GhcaJBVZuNWwYzV3oeoDVjS/pDPQuA7lXNDOZH3IWtkOCUNt/+hQIG1H5BLBKwsPcQszxpL2nrXyTrkctAnCbTmBSGuMNlVRs2oRYwXXvNkSTjgR3zIa1HGsElghfEhIzUnFG1EyWI4SaTp+dn+/gHKbq7yD9iUInLgXXwb4oJ7z46fTX58QHeJECPpX/00jMjP0Jdc7SkMkbbeK7VrOk9nEAKzCz7sYsyj2c94jM8WOpSEgUoe/B63mRFeUKBBgIO/QL13bxIvU2TQScrAS4dYZOPzgixgtLpq30KkORxNrUDPabDXbqqkkbNe7fO7bAAg+9Kkqgizlv3ChuGV6Qq/vmU1Ui1lIs0XaUjs/LR4zwT6E47RO+ne6DgqDm1MsFebbk9+C/jQUGVHmLky4X7SdmtChvc1c7jt8AxRkqlzq9Tqi1hmKutbwGaB3FEOAZtUD08zi+FGu7l7aqhoZ4Zu4Ea7ggZCOkT0MJv43HwSu6dGh10MwyQKbnxF4dmH6ZrGQhQVoxqH3WwagwK1DjjMkDjslmVSDXGSSmXtIx+EaNCidF635KC58BtPFm8jDi9TX8Aduc2z3ZNkYU5xWVX9MpNz8/dwZRmTcx/das20ZFa6uEG1oDiyBitMJb0x8AkOrTUg/HyumbtDCK3h7ymRIjOOhwFTnXcgoOWofAzCgh9xUznnm5WqGSc1GDfpQtrGZGSjGqwSB2Qe0hQ32KblGAubiDW2W1yxZbowrvQYeHQF6onep+qgz1pgznY82sNRBNc406hQzv2JbTh8gV1b5OcFCKdzmg3TTG8qYmGiddmEBwI9Yp2x/E60y9qD4NUfW440LaKoKKVQkN1oaofK/hozRS61pkmh2Z6Bnrf+DAGd9LRqUzxOWMseJXg+5F0yNpRL50iqj5NXNKOwTaBAhH05iK4CCdJ2ZAtJFLJRRwFU2xNbIPZ2uoQctbusgu08BNsbYxJcFtRe+R/GnAtpkBfe3tbiGfR0P9xAVM4dEzdCAADmX+JgSRL6YtFstftA3uHzU3BW1KZSPM0FNlmeUSvkUZXmYB/a0DxdZoYRajosoQVtme7KRVV6+dunF7YRLch6Zs5MQBUzzbrGqBn45Rdjj1cIxiTCn2WGAV87HvaYEaI+46rDUxvFUGPpAZ4A2qXNIdgANmsHCqOjo3KrqpUqorWse6DXwl3u6jOHaUoCV0AichWzb/0Wob1QsO/RETMQMBuUUHkd+Jv5EUURjcTiZcuvGS174x18Wrqtys+b6I0VFKnRQoC46qY1CrY0nHHem7sL2QunNBE4C9ZuDptMbj15aFj6BmuZGC4CVaEpuwADMrRyzXAvQGf0s1YDPdq1YwMIZRKdgpaQRPu0PV2YrwhI4EM6J0FrFnok7+81mrX6CVUQFkBP/ZjfZtzFq2gmG+i551R1EXagoLZX6lvIx456V7bfSGSYtGBoCjN1ibrYg3GW3gNZngyQh0eUe2LaMTVJtiL2r5Zkg9MKEkKLfgXlKhqRmN3agLvNpUfHvPa7asO2hgT2hCofFs6Zhd/GxwY2PGkwgjo7RO4GY3k9CQFt1WFJtoDZsUoEDTIKuV+rN+YGyTCRBaWWsTkWq2Wjd0blx4J3sHUaxfkXNSIl11EWiYsOqVsWrous5sXIyS14VZgPJFJi1pScoJhbC6U+HcLRSF4bw/QL9WgtKikq3HBfjQaBqx2EF6UbLvxc2AvA7ozGk7gLaU90ctdSvsJaXhdATjNo01zfEuFuzoOXnWlNdSu4MoajQHFYxjijnOAR1Z/AtP4gQva6tWJEHihsqpY6bMnAOIcBSP8OjmIMGGcGNoDWNdQoh3QI3DJXF7IQ2cDA+KviOEY9q067JT6Lq0cvn75jqH0zXi4sijF1g9Y5IqyT0FZT5Z0mxV/a4rOUdSamPJjm/b+NZq7/DQ2syxwQad85KHzqw7J+JtczjyTouaqYf2xL/x0hnwZeXH+exymV75+ure3p5VeURkbQpaD/0Hq8dovOn6/l0DJbycLIIjInd2I8lax9vILIPEmO7sbAHYFPnHDZx2jn+XCdVNxRLbiLP4eYxm5VtSakj2JAndGpcgpynuI34exQMxjrNPzx4U5hh08gzGPpKyMjOV6hmAJVKPaF8Y6wzqGorbMZOQcndjo2Dwrc+PsOxEVtp2Rs2e+R18r4oJ7nNCJ3Xf3UBKCQEkx5BFXvKAkEBvyhKDxwYPZpwjP1BweID7Age7PQuvUMiyDp2tdjThth7YhLq5OFYK2r7AEDLcfUQRUEyzvIvEk4krQH0FDTwp+vsaQAipJDelaDnEJdoReA8JYpz8Y0UxBu6mfWGMZlfyYmZJpRVvaJep4hx/sRgOchwPduxPaNuPgbR7f1I65dTwf/b+hoKEdIeORwj9kQYQjnfEYBv/9CUB30j4H+lNDqsKqzARlIlTKh4dH2nWek9wlA4mBs2/3rB7CiNifKuwCBgHo/Ur2a1XXXVGNlQx5QfFqnYs218UtXrfJEXnwknSnlTo4B/pDUkSOdkeyfv6bGvCb+3NqHVg4a5DcV3XJXBDOJZpXWc1759VRWy87dZACDppCx2AVgD6Q1ZtlEdL99LoNKvLJSxhZuyQeLMBJAVR4FATasZke5qB5gE6yMKtFLBBmDCsoaTAnBWHZvA2mm6e4xnBpUGbeZlWQNh3iF0aLFaEPBwBrNQHo83vk2aqGDhW+Az8xYHm4PjTM7Js3R8PbLEjOJgZwpght+f58RDEEXQLgTP8TDXoBoEnAzRMU37qeraqx3DPyOy73VKcoeqz3DTrTTMrTLQg6tbGGcRhq1sU/eqk+bTK1gRmmiYd99DWUBw1NhddUvUgYdXaFyVRs+xeNvumg4TR10Ikau0i2LdvF1RDqMDr9cgEX7t66qk3Frs37lqVz66W5UW67AQTf9+Shc6Ch3tk4hnqQL1QxziiINTZ66fQ7ZIDTduDMISNdAkCRscFLtYB0b3Bg2h304Ggwb4o2PzrBZdLHWB+U2fdsPOvFWaOL8zsbTUQLP6lEedUqnhyu3/qLGw0sqSDzrxiftfqmGYWrKKORfisCar9AM/1ImCHH17a+kBn30Tv8c1fjtVwCA7fHcMdLQcCAN12DbKCVJTsFn+6R/Eg2OgGuVtkm6aCZnSp3uuZ3eL1676JWYD+ZSJYaE9UdXdrVCWp5VnotUAbKL3KwphEiFXWXJdIrDEseJ39uQS55jg5b/OsKGfZLdaXHBwlF4CjZX6VcySSvMl/Q6pdtkI37Cz9BOIv9eotgZhbfokDNA6C+v3DSzGEwKyAaZ6mN8BrfBKFooseFujVhQlVsShkYulNpwOYQnf6zN3zKyCo2U852Y/iTVGvszmZw2LD0Vg/imIyfSjuZTZB0Hdisy/MNkVOfvTYgdhqPocuaNB0LL086C0I31WLq/qwfg0DlJVnKfpGCUa8+N3NEWU1oCetIg5JqPmhJirKBahXZJWuKEKF+URtYlO3W/JSJC/8byoz6UpnY7dpebe1TQcfx8LiWPCm9Z6S5H5M6dA6vXgM0zHqaRz94UXRunHV7OTnVP6/I0qEDTDsLu3+BzjW5GclwYg13Y+dEBJvyqiihGcsqnh7sqzTyFuxaX4+HEbTxhYNSmGPHxMyY5t3+9StwRFZSIz5o6TnvDcYwxExptD2gPDeO/3gEM60TAHPr5/6S2XC2mHOIPHjH3fg5/65bnbqKw7BPoJKiO2y1UOrf61tARfpneNFJ3L+406rbV/P9vZ7X9cC2O3fVLRQ02JyvQNoIPMKgfM8fOzW79tHK90wLMU40UPuWCW71n1eqMtaDDC2FB7pmClv+hRIx+bujy7OCqdHC1NFgi4PSOPJyKIOfty5HdHPMeDcACs8UCd8I4i7c8wXICkp4cRH81wErWcq4tcdhIe8vet5MhJzSWcNHw3ovhx+6HC9zOd5Y090hiPQGwNvCLtSebvRcWIwUtxO2Vdn5r2Q85Nt95BaCMfx6eEfEe7DyeGPh8eTw+8nR4qX8Jsgn5JwhttasqIhTgO2yMvCk4wdeor3XpsBQepR+gK4mecVmpSw+a6T2mZ4UCpzpPznyCZBOWvYaJPNQdZSOBvURawOkN9Iv6p8INIxGOXk7MO7d29P3x+9mp0d/3ByOAGMHb8/i82G9D+yDFqCgvQev/pwevj++O3J7Pujyds/zs4O37ybHM2OT94fnf54OInd269HUcvnvEjtUHyP4ELBGTwLvgvSsA53PLgXxpAu7p4h2038HmNfYC8tTT+H7UVtkfc8X+Hlg+3x6c2uNbKak1dlWD/wKnZbG1K9J7HbyqF1pyPa7iiYre2/77yg9ZMuPPAxLauG3MfjXtY+jrbfHJ+dHZ/8MPvj8ckroDqgtA9HvYT9+Fe6D3oK1zGACJD3HNUBdu0gIXghcNkFT78ReayvLPpNplWDvr2vxz5txxsgrAwNKkuBgRVFmTomRfUE+ukjnjrXpD9EjqXm6Ocq7GsGTX+fR1nHJ0BMx69mk7cvgWF+f3h2NDk+OeojLX4zKka27wLvpx2Ue5Uq4E8G0UWgg2fOvj/9lFa5NYKj189ff4jN04QOCPcBnv7vR6dv7dJBN1XSlVu+OIj2n/w62FY1mx2+n02ODs/ez97C3+9/d3xq0AqqvOvygc/hA2P9eHT6b7OzN4eTSbvfRwsc3hmy2YnZWVLPiShks95l6CWwav/gd/E0EjX8H30cH0iVQKFSD8g2UGhNCd8YLDbzbKG4/UqZfXpI72s9D/dtiIo1uQ5UuAU7h2clW1FtJ497N+6ZdT1Yu5N/j+fl9z3H7roCvLlzUWD2X/54m8GCBtPHgBuE3NtGmVIDBGTzH7RadJIgkLeOG5G/LtjKml0V6YU67dS0kiaEutWGWukeI+z3x4OTeBKNlMsvUczVd5a/0M8lOxNQFmsZGh+XuscHxgGIgcPXWboYcQTsLj//INS5EmQ+z5LwQM5mH8agPxHFQSg5RNC+LB2hUcJC9lma27At0Tscv9rPKIlH7JKNdF5BRxVySsyiiAlZIgpTSRvimQD3PGKb8V2URv+1gTqK5wR2SsCKme6yyL10B/+xKg2J86PXr49evj/+EQXx6dHZ795OXvUpNY/JwDGYmUfCpvtCqkd+qEonnDoUcvCwmJhgYMF4IFDExlwzSgIZXWC+OyrjZziry/1ZWmw3D87Hgi2g7vEpWRjUdxT29K9jXzsT842KtAAuacE5TCtQU+gdmraDUK1DrWxYhDtSNR6oEAf4h0m/gstJm6bCCxocWWFY6gFynKh4h/Fnv4+mNfkPXccqnqYdS8KV4bARrOsL/Ag9Jh1O5XJ/mO9nZXNxo1G4DlCczXriZTsRzfjVBk2fyKWth0Z5ONhSC9OpXBwoR0HTG1OM9uTJsm/DpIaZgQy49hmhzbyCVSO3MpcNZctlDelpelkKClV7E8Wiu8HnUCBbHjwLL2HJOw1Ch1QvSkMeyjDrdgs9vm5wRNZSk3vGn1nm1Y06aVysHZDtW15ClySQ5IUPogHH7kxsmB+hZAdJuoHers5V6ne1RCvB97dJJ8ZLJ4Rpx7jxlELpbo5u0XDoB8SNzUpUcogQ3NZwYp3tocw6kkxn3AccnBdzGlNlt4+S8xB+1CMtW5l8QTqeIZwOdIBMAkUUduTq/ZDEJ4Odb6H3/DDiXXeP0cvL6One/r/6//yGGQewdOQnIEbvojq/9QLDqTs8weQ85fBKebCWfcows065xIjFi2yeYtAKjAPKYYqv7Xc5Y5TOnyzUYYMg0ICS5cCzr0a8bEDfb3pfCKzVEwECTyiCdV/v0omMAP/VowitWLD+yNFxe5KPSAVuu++J/v2yzEtsbvwb5FwazpEVSLvkC73/j0mY/ga5l/7WKZc+M9PScADL3yAH02DOJZNK36RB6c3N0s1N9GxyrqKpH5i8iTiUpGmatPM2TYYTN1WS58dkbQLwqfodyMwUitdRaZ4mOs+THC9OA6PRAthtp2pKgrmaJo9J1vS1UyoNJj/Cu0fPvWMo1ZEC/MrpjibD+Y4mAwmPJibjkWQqOinjcdTaL5XDIXlUFqSJlwapQwc2JVL0lXIi0RRets6gHnTrK6c74vDqoYOuN92FzNjD/jk5jNyRM1mWzJIDmZb4TNI8qek3IAIo2Aw96nixdhNk3350c51xLbDTTWrzVqG/7Pcvl3vBCfSkeXJZnVyN2wz9SRtGpH5U/X82udPh1VWF6abl5VhtrsOZRBKhlrvha2B0kZfygR2MbuDbdr0XcXRgHaVVBh1yhtQ6B/rDxrCkVYm3iKbU7YvncP1epfzVC25eZSukhE95jSYhnc7BDjaQWGYMSnyTot5irAcUybOYKf43kkDRsfe9pNhpL2wzMZzMWVDCCWT6VF2Vw6U/yYwOaNXRsJIlB9My0Ka2lqGiXgTU+2jVqcFUxJjKs9pgIoILjtznbSDHLiuiFE3vvu3lG1OYE2uGamjLSkdb8E7rrxmSbPwfVexP8SxbYggNXpDMdFLQ1O9qONEXGfBDKBeaRK0L9j4j0ztRyROj+TxxWqYMzvM3+Grvt0VMH+Id9bRRr3FuoUzgjUKUqWMk9WyLAXKRZogIsUS5KFuMH3ONMAIX48Td3OMfMN9N3DM1OxOCosvS4XKp8qrRU8+0mV8bnzn8bG1xazXQYFPbd6F8JZceMBWBKTGN7NynAnSuEhZI2oEhYBqvSGu8xfEcJWWBnuR0m6d0Lq/Kaif5p1wTUhFJ+3e4peCivolIN+e2nUWoj+IQPO2Tu+zYPAzEux37RAbuqDAxmYw6xQ4bM+GXB6osUWf35HDk5uoCxHzkK4jEWkoPDb7V9Az3sM8GaH6ggl6ApGOKpJLZxZ3j03T9pFwSbpxxC/NqxuY7b92O/NSeotIIKdwmHMnuf/VNH7vDoMQCCVngHRgWABo5yiUJ8sgE97wqIqSZeT7B7OThsuHhQuCnjF9OoO0pdZfBKQzHh9zjXv2iK3GJr9CPgg8vb/fgvNOHzkxGKrWxvdcKt3FjnFzi78FfjGsFO+Bcq4+eaejehs4a0SPQVOPvtMdS2yqz3qdWeLgxg9y6Ih2y0rnii3pIrimzfea+7/xfaHZ4EnlAHSsx+0w4bVZ4j+1OeDv9JZvrMq8NbqwduH97bX6+2+HUfKFqibtQ3KtLJQ8kDGc/MlOeVR1T0i3ZkSxApqxFXmY70xnmOFMa5ZoX8NiEetIZPv4IdNVNtKctXUSRjpYS38TFtY78zFBoAJlxjsHWpNstDMW2xwGibYMG8oSbN9qsVThqISemLZfkfFgmL466ck1TmpVq/ONeecaXk7dribknjs02djr/+hMPS2x0iQPRQzZKV2s/7cAuStZwocMcBZp4c+GacgFtVs85dwScFzhN6tpTrhs0uoF8oMeTOrDjyR9eQmdGNmBsygZVVLL80wc4l6C+zu/E5/tPtbm8ovcZbzd3kf5iH1n+9X2HW8/0or/Ev/pwk7wMLEhrZ899yNcn+FvJRqGYMBr0QuxtxH7XYt9O3ym4Jvzfm1CStCbo+R09l/QAjrR7l0WN6nKHbs+e9mZgnQZnS0hL0PpBO9zB87NLjIdXZk08nQdj2ntvy9gx140MROhOTLsr7MavuxdcrWDwUCOp4ME5WJ2m57+M6xSplejAd+L//W/JeJTHvCdTE/NexEET+/qLPJeh2DMPz61YM1LcJELFjOLH5FBrL37GlIS20QbNkNjzAmFN9/1RlMoO6IfhqQo/VrPVQgXmtUZT4Y2tNm0LpBe8GCrXIZBY756o09EWRsA3uk8ZOkhGMQoOdCT67rit9uPbgMvOaqMcgm3O5dT82OFPFx8oV2JJwUZh10GmrXo29pd4AAYUet+Sh5J2chl5nhisYDmr4icWLtgLR2kNfaA9ECYyLOlOMNzh53wpvOvxyK361eZgUXB6kZqOa0HI9tlduLnVqFostAdafUc6WG0/b65RJon+U9BTF3JLDcm0Z1ZwKZGvhZhNTmbd3fz9dzWXUGQ4TUzeb6hZ60cGBNJKCepSaNnRKkkOMORCp6446brLZhqKzzj3HucHIczAEohFfMd+o3RTc7J2VKpim/YY8bYYkNBhN7rdm5bT3Mpl1zcfbn0I5HATlOU0U3yByBYbbieJwUvm/fQMf2zEeOBm+nObbbXuVrYuCdnU2Hjq8aOx5VPk0/QMpeKiHNlOpYX7u9OE2OpikVs10XgUZEljMQ86NWbLXt/ZM2Bq1DDGnjV2XQfcgjqCU/Nf4+Yl42Q7CXsoAqEbhup159LctApdXs1wcKrXi0nFsyvRpZ3A1O7zFK+95OUZbN5K7fOdl4m7FZQrtcE8xpw0s43PVqjHvR8oeFwS0RBmiYd094F4hhhSTIo7+wpDMn32dNfJ/RnsnVOwOcNQN5ZZLupm+M6rCVvhzSXMJvVLh2DMtelsek9H6suHO4dWXLQiIu+NlqSHnHwp8e+2oArzhwcXKli6NRb84Ogo+bq8ftpl8ySxIwC6sdLRXkOfD1yt0cEX2TnE2kwAs9C2AIkX/ibCfARlhV+owgfR9IFX9g+uNsho2LErXiB7CwSN4y6Sr3+YvI7opsUrt32cASR1J8/0JV01QpovYNmXyiSXcvPCEp9bZgt91a4zSjL1CYTMo79MeVXliy+9cQMKW2rZ4+/bz+x9+8ytxjnVrKfJXbmfBa7cLX8YLo59YaQB059mmHdV+Slf4B3eegUUHgkv3K/5EpA44wI6Ed0GtQZEV7W8faUKnDcSJlpTcivREcHQrLEflCBTnEwq0Q/rfTe8AOKFfp4VQDMYhWKbtbyZ/rpN9hw4bAX6zTGt5cI63RW2ad85nknrL6GdJFB5HywfSTKpXf1tVAQcinXkPfTiHE2cjlaJTZ8hY/u9KLU9qfASW2bvY7gi4qzY3svmwrFgA5H/tjMbo/loS1CfksnbKjtu7n+OBLSzVY5DWu+5tE9YtmOg+b5GCrUEpZgaii5SngwKDOLQH1Jyka7Ziamzx/BqzNfJ7DCk8MnvMOShsAHOeIaSQR3V2EYsmXZkAb4nwkRoSxA8jv40nf4JJ02zV7mPbQTKA9ubRbvPjmGyMJA5VsmAhvB30Vj+TlGRoVAV8wn1XSc7zJWNz9GupLiALQDaKi8x9Sxc8viANSwsmPY8eWEihloBeg+SF3YCvdnuP4/rS7yVFzTYNbO2eX4rUuJhZsqQgfPMNTjiaKiwRJHlK6HilfSFXBggj+0OLdmQg0094eJld2naoRzq9YM4qJ3Lg5no1OPupMIrA0zICkVPyRSLipmd2XdBsVa0OgzaYDDIid0ztLrLdB2XdqTAAbmdcnMZCzWxV7K/Lxun+ATLw19ZMhEbbxykMJ+fS/FOp3GXtxtylPTeX4PDd0ZlmoRVPZTPGxYY4Fldvr8jgTWHwmHda0Lja+S9UaBd1+//Qjny31BLAwQUAAAACAAAACFcig319SgRAACFRAAAHQAAAHVwc3RyZWFtL3Jlc3BvbnNlX3Byb3RvY29scy5S3Vvrb9tGEv+ev2LNpC1V0LrYxR0KpyrgNg5wBzct8ugXJyesxZXMhiJVPmyrbf73m5l9L5eS3PZQ3BkIYpO7s7vz+M1jh4/ZD3XT8etSZOxWNG1RVyJnjWg3ddUKtmnqrl7UZTtl7Ef5mp2womXdjWB9tbjh1QrG3xRtVzfFgpePHrMFzlz0XXErjlu+3pQCH3UNX3TP9BrslK35BwF0inVRrRivcrbBV20nqoVg4n5TFouimz6SW8lFMxfVqqjEXFP46pgtYQcd/J5OWHI6fTp9mjxyxuu9z29PvcEVX8Nh74oqr+/ajC3Lum7YjD2dnjzNWFUXcOwZ+yJ7xPb+5H3DkeR8Xec4KYFzw9ptktlXLTw/OYTYsoF94ejTy4yt+AZ+e/n28vKQqcDC+TVvRYn8kSxHQl9cqvPMzRkPoXbLyyJXm+941yOpBJh5WyDfeZkcQgQEKBadyOdF1YkGSBIj8EAT9ivML5YsPSraKSzEp3T0VElkMmGgTJs0+UHJT4uKrfu2Y9eCcYazJMOmyURRS74HtevmbcI+KapPGEq5NTTZp5+yo+Q1aFtf9u283jFU7o8FT68CejPgiV7wParXkDhQ+XjI3vR534LihUTYEuSGtsY3YBBgaa0agIajDt+In/uigXewi0WavKpLAQpo9oO/f2OUo+MN7sN9JqocnwyEOjwRUTUELmieZv8RL8vUbCV2zAmtMCJbNfOMAeUNBxQwtDJAj7LkGzJLXHcC/+TJ677DQ/OYGl1lzFLImxrt6cX55euL9xnsoQHQac/bF4BJddPqV0iTgKVq6rsUqOMD0N1eeAByn7GSX4tSawodvwR46m7S+4k8+yI9QevTwqXxsG+pwzf8VjAAWkUbRGz/ABRUTJFnRPluaElaYIpnngFhpV3w55MfDNBZaDy9DF+iNOhs8L96dzFA1DG0TRW7n/w4AAeYNUAMNfh1sQK8mDdiKRqCdbSTZXEPsGCezZfoGAhY5tKdJGq2VVqc5z0hlWXBQ6XbavYLhLx5WYNPMkuovfYiJUAEeZwDVhZdD+hNT7RWPXn9fA6q3xVodI2dRmiKVvCc2ddm0nPPHZhJnpeAyXoYw7/N5O8Awx2nMZwNriqBQWCSa+NdzOxvrc81DsAeljwLTHdGMTXKrs/v5+B13KVRrwGiqx4sG16BLr88B8EBQ5kowRzhGe6J39Oe4C9D7CKC/lGyETcxWCYy5hAPlOhNqKMyPd1j+cBtmm3GnKojAv1qwMiX1uOOqp/jlZGk4G3fiLWoumMZfwyUUSDuc8AqMiG55LzNtam8uQHDvanLfF4DfphxX88Sox9gyHwl5psa4qotvV3USKUTc4WYeuwrFfzNcwBPs/GkqBYglNaY5+t1XXc3GLzh2woQDN8AcErYR/TEN+jyyUlpYFFgIUx8poEW/nv0MR7ARUI+/XI0nDADtH89H4a1e0KKoySEVs+BOys0ouubCjCfPIbdbl8VP/d2LwOo1isp76GfsiMIGi/Zb78xNBRun8OTI03e8TOnl553fVu1/WYDcX0slNfbk8cE34l2sRKNWcQXwkBcO6QQdSShFO2cw2KXQcThxTOxuCUG/UTJBfXEx0Nt4YlENAmFRhH0wWYYnU+CXdswI/Gcrdyx60XVZpOBE6XdBc7yQXFbNFbz/FEiscvxMfpJCHKa5DiQ0cmGqKx4bLAK/x4ikyYfQBIOHiIPkdRAQ+JYF22rQAcOmhfLpSOA0C5961JThwkG5LPqnRN/qic7wk+CCAoVrRmQUZ6cnaKl8mr7vMc8FswntxiASo6WrBT3qMgB+DF3TsEYIQRCJoM5+uMzNA/HLHDtYGmlnKjVXVODm/MxoZUh56YB4TRbk2kgYkMSDGes70Cjr7eMdxCZACziaEULB8kzd+K+m8Ny/bpqx402YjOL9E+3joF+/3eVT4scxJreYj62Ncy/ctkCCYZNFSZGyq5s7ydKP16e299Bl37BISlkJ+u7Fl5MQOxlvSLdOJmMZFCk81RAwU2gUXSY3BjXhp55vem2OuuxuuvrZ1ejBpjVjfppsQ7Wl3pkElLcwUKpi81y8gIGADPkNuXaVb8WTbEI1Ggfzo/husK/w7A9rpgPSI+jQHsoru7E1J36FbDMVzFSLjXCVSd4uoRoFbDs/iBNcmpxnbacoT5JkijEa2CklChWKVK5NwizkaXWh0aZHnUfpvRyjyphzn4lCb9/b5LtNDgwFnew/CCDJbQblwf0iH3KfIaw33CZGXuqmaHYIRczqToceFO3BaVMcEhMWui0vO9qwJliobgm+fDREaExINJmXMoxKVIzs0vzwk15aWuOmDQ+MVBIuTcJ2jbCI+3NSH4mO8lNshmcxoMBu4N42gxbGew1SJHlEMVJb2A0N2UnIICdQyAE7hpQ8p2EfE02eO+hXyVWnOSnECpMgZ0kwIhaxcTTWAR4cukWfg7gHklSbfESNgP5F/frDu2IeJhOLIoKq6/XpfCDm49+gDOabSD4TmxWoaQ0Mmvgf/VMl9XnZVCXbcFz6QjD8UzEdFsfYpKir3xH40vrJMevP7OkbyGzyAVYs8jxb52l5PPr7RxfBtHPj+EOLMMt4Yy5ZNHcB2RHrCYopAG3xgtpezI1mMLkFJ1OPcNtUREY4pK6yWGoGgA71Bcl5RZC2BXYuhypdyb3a9cAXC3BI6eDytzskMJcxrxSxMytREgHGSlEzKgOkbGwBDGLVCAyFik+zJzag1zFlh5mqvIQcUIyFzBHn4wZrfUvKDcz3j62BXpHWAnlBcpTKB859AFWs/3aoNZqe1+UmDKSpx9erdDorJqG8lfTRvQyLCrOYkC+q2A1jurx4hk7dWE9OmaI67FhY2wIYDIjbXeKZi7gS1e4HinYwUsN86xeYtZTgop1rLurH+IbwSE8heBixC3OQi9+XjnYv2ggFmmQCkYU1vvI8zhkCnNBkj+Dgb+IpnZCNHhb1R3jMoUb2XxQO3fDkTD8HRN5QOJrN6Dxg+IxCirwAT0Zn+o7dGHlRrLH4AdEBfy6FmDwwuYfND9T0auwUZHEFuWkcH01ny87lDzOUqkt5LENV9cUJ/9gn7Ppd3wBMCOe5HUPvncqNi083az5fXoCK123scjNi+NHxrBjNsrW4x2CmQDTzTY9Tjm06aSA0IsbYoT2EI4BSC5U6EccdYmXEyxo4rF3SPwrtLQg6vVD6905oxGYzgjUSVC10XWUfKODosfseU3P2wJiGPR96CE5PLkBkWrqn7VDmozDooYTU6JGfqNAl3FydqoPzNhPqAZfXII8CvVkwIRQPlcFJGajLLr66T1wySlXevzytcCnRLzE6RMzWzL0AvTTclNBhHtCLKq3RS6Dgv3s0Qz+qP3YMOwcqQ6jhBpYqPXKw+apPYyt1+s7IDOIXFPqV/Gdl5jlyYtZ88xW3u1OnbfxiwF/PUocwRW/KERJoSSVsyZund9dzwUnmS7bkyv+MyJFuo20VJzguOhDauPoJk+9iPufMGMDygEKrzP14dLcabDxY0A5RRdbDql4PKABQLkeRdxnpV8ZVNsIKrO6khvyO9OUJ04Q9qGq76qQBUMuy+qaofWEdjVWawuHTXwP5hHS9VtJx38u1cYv3bpCvIDR28HereFWANdUMVuxBsmhEnFbwluSZsmCMfzejNXyZAU7OFTG/L1mEEdiV0TyLsyWvq8AUjGVCzcKkYYWNaKmQPjxtT1aDyJ4gaC1S5X4s4hlubUf8/xPK/74xTDfLQ0MmfwO+iUOm/frXV4RyDbA6SqZnwJYxR+BE3VmfcPtlr2M+7giV56GsvN1zLIn1FZTvtVcGb5CrNFb2M0a62tk1dkPNYrG5v6ovS43IvA8ctkLjJ1LPY9dNGbSBhxd9trKdl0zu/eOcf8zYZ4b413XOAsntLL1dKr5atwVRt3fI9O2pB8hfuJK8Di2iHOJrp1Z6DbBBELHqDSrwfqMZ1ABKGCYpzgatEnRfKdEb434T/YfGgwa6hILQrJW/DwHJyGDUznChmhAaiWIl9JG8D3ZA4VP1MuAxqFGh0VqebSJ5SonAhZygAqqRN85j+SS1PE3GGxCpz/av/TxjzZMbPi2rPmYBf0Pm87DLYCkrvxmpfsxBu4+1q9ggGgH+Ejd8o7umFoD9PZ5YVD/tm66QegzGZrjR9kP5a2mkOHhnRe6eDvon7SH3dUmSQI2PdQ7UctZyfZmOmc2DaLDE5tuItN3a8lSHbNd3Ig1p0pkgVWdBZ/qpOZYHwSgRg4z9dmZ5pdpHgl6y2fxro+xkq6US1JSbf/4xPakUAU8Sixeo786ea9oqXuCSC+9PZZeJdawvWvJQZXdLiuPMO8rU/m2y7TDmnHsZm9/GVmttWNgvlw+tQtrBZo5Tfu+vsGrsbgCAhtQtd3dS0ivFLrhfgQfCXpUpP4SlXLDFyJNfmrrqoQAEbQMXgksRszYm1dvL4ZtJv96/f1LBmoK+/mb3JaN/DE/16QYUP7AVyqWdrB8HOlj6I6TkaQycdt9oNc5O+tq3FOqiGR0qwnCv67v1SEy6osjYEAja2SRviKLQww+6FsAaQoYNJtZLC9WRUcMPkeEh/DZ8i3EeRQOcPMOq6WXYHgtnUAKje5svtl2JCyark7ta4HkelwLcPScsB7zOreb0GeZHQORtEpc/edB1573Mrhxv8Vo2rstI+2QkItJBmnkhnc3Tq/RX6B9RTvXOrRqxKZMk3+/e9d+/u7dr7Ckd8BHsjmIQkUKEWjixBskLd/pmseHU0p524Bd4fZpel4LWeumKbpapb8X4LnVDk0oY3e8qbRPyQDoF3Wubo7evnlx/CXGn26H17vKJC2O5XXN9lsKMa3tLJt6TdaDpwYXg/6nWG6fg728UPZCghizELBSuh0zeigmptJD6Ourxpm8cqryAgd/J9oWZARzHHMxfWRq40+k28t2OkevE80m1j4J067qqP7YiOCSWZrCDnJH8XcHtrcaHqlIAK/E5G9em6vpIVc3p0pJR5Z+sMMfF4E/Mws38hdcBpsdeJ9l7OfGoY6dDrTXue9Q2iDMyCI7jl7JWlUgCsyMVn7UwX69loon4nXiQQO5D6VekYibIIU3Dd9azN7NUwcJD241H+0lN5/9DChSv+XB5PSHQnrr6oMhTTVo2txHdviNUdw0nBWGDaD7Fol+laTXiX2dpJcatJXuWyn2SdOY5jqVwENL/iM46ccOHuBh4VVdcarRLTZdcvDWSvEfS32lAOyGo//UHZis227EM8SaDq9SEStysSg5trbrIboCBGhyIxoi14G142X5LXULluUxkdYDYQXVaQgborxxyrCtg5a4PfkMS6gQ/AGniJgqyFCujHPvBNBqa+eTYiajFQmKfdkx+ly5xWvMLTP4Md1Vdg4aEuPdh8P28YmtqpoC87BtZOLoh9NMIqNu40+H0/YihG00k9z5/d9I7OyddX5+5ycSAbjqo2LPupUBnWHI44z5lz/mTSZPHS2EPLBcdQgIr0XH6Rsg4rJK3Act64M6QbyH3f+JflAd63KPZNhDhHIqpB/E1nY66QM49dFhHqM4cXUFU0EFh+GceROEcYNQziXht1PpjagRsq8qmBVeZL+pN8eluBWl82mWFoiDcKYZWSMcRcVA0nyTYe+uvTKZ1oCgbrC7ojh6bzzZeXGEnaSVjAOwsYhK8OrWaPp/VvjcW+b1CqI6rdL3HvHaNbiomwcVrv94aoxL0pXr/nz4sILKAws1scIL8Xdf8eWwqoos2spaUke87cR6g9lxajJDkxAey89OK8zKu7SHc1UfUpyHu8xzp8pC9ZhviipFRrypX/G7FADqtO+WXxJHELlpIoztKwSptKvrsj07W+d/b/u1pIrluv8AUEsDBBQAAAAIAAAAIVwz3+m0X2AAAApOAQAQAAAAdXBzdHJlYW0vdXRpbHMuUuW963YbV7Im+J9PkYbpY6AMQCBtl23ZVI8sS7Z6ZLlassvTzWJBSSBJZglIoJAJXnzsWecdpl9kfs3/fpR+kokvLvuWCZJy1enutdprWURm7h37Fjsidtz2+9nRLf9l3xWLdbHJnm2rWVOuqvrW0kd772cvV03xMGsuyjo7KxdFRn/r1XYzK+bZ6U2Wr9fjV1n/6qKcXWR50+Szi6LO8sUiW+ezt/l5UQ+yvOKizUVB4JqibrJ6WzbFkOBkZZMtt/SmWjVa3VWkb3WxOMv++7/912xbF9nDh9nft/miPCtnObo+3iNwT6jyapmtaEx5s9pkZ/T/BTW4KKvz7OVPL15kl/liW9R7bz749dcP3mRfjbIzHXo/H2ang6w8y/rvlfW42i4W/XxA3c2KBTV3yuAvitnb7OqioM5vaCDZ6vRvxazJZquqyUuavW2dn9KkzPMm57ZX22a9berx3kVeT/lt2OL1IPvXvYybtBavqcVN0Ww3Vf/Z4xevnw78d1Qfn23yZRGWqjarKwB6lE1QdlFU582Fvtj7DZ3+brvMq9HZpixoHm6yeVmvF/lNtjrjBXxwvllt11k5L6qG5rLY1Fl/W83p72y1oVnfFPk8y2uCU9NKYBmqbL2g5WnKZoG1pdXMrwkLFvkpTdSXMva3xQ3hRUPNbKtmtSUsmA/Ge2vqc3Mz5ZLpRJzX29N+b/pRb5j1MvrnesATvlpSq4QL22YVNAsELLK/rcqK0E4GUNG81LQ0wEn6eEODa7IVrVAFNKP1JwRbYdmuSlrNnIDRgtHCUf+AeDRFNyNG6C01k8/rbL6qPiQIl8XmbLG64gbRgT10ZYpfU+5KNA7uSj0kaNfT2UVOc3mU/XHSXmUpN8h+/dUWzN7QltTyma1wr4eF/Y3+1wFTi+u8bgrX3Gy1WORrGtYRZq5nOFOhC32pNMi+OvLdGhgoRu16vSmr5qzf+0AnsyYocb8Ggko/EjbQYmSzRZFvqLrMQFkRjhOOYEMQICDSqsoX1P+/b7G9sRGqcPnGBOoVD85ty0XR1Nn5OZfZFEuadp5xaYA6cUqQaciEi4S3hJfF5hJ1GScJGpqg1ovlurkhvKPxnI+zx7RSS6INM8bJuvA4u6G1XK23i7yhQVCXV9yYjOP0BmSJHper+ZYaX52iMd7uc7RMZAUTd71elDMCdwNaRqSJBoQWTrdNA1q0Wssk7MIUftNGDH3t8UJfvHeUHeAtyuVBqfeqX3iRacTLq1o/DBIEwgQbCnEJWcxnqw3NDhHf5TLf3IwaplzLgkDNaCOVzQURr6xcrhc3mOnmgsrOiqrm4dM8rIlbFFRu3lyAOmBOi2vsVZkwpsn1OHvK7579/N33D77LF2ejn1EhI6JFO6yoZ5vytCBoM+wzINQGNOZL7NjVVbEZndL+nCelXSekcFYX65yofbG4Ge+d8ZhovpcEcrupi6kMaGpUL1wD/US7Z5hVtGmLvJoyeyBWVCz1J0/k/f9rVg0tOkGzXk4rWQ3fGPqQ12OsG01NsQn6gVWq8L3erolY1vXPOe+Ruk8VsLXOqXg1QDFt6K7CWoyr+C7dVSvoPCryuupMomowlg/K6oNs1u9hgUF7sMZTXmMiQ0aHmDkpTtIGKmmJeDd4wjP+5Dz7b/9vhr+9HSsx4Oq8WZhXzfJFtILchekLoM30a6BNz+0DtxN8g8TlqiOidzbSATVrXSNsGEhrv7k2oyn4l38J5zIY3C0NMYAvs91NRigTNr8LZFeXf7M5v72/eGMI9KjFcKwlHVKPRD5tCBQC8tfkQTgiPGzK84tmZK30bNt0bIehveQiAyNMPcClcU0GPSFPP2/yNVEBIj0k3xXnqw1WPBQ1iLUT4V5uF025pn6BxUMWzWrCj9HpzQh/CU4oHghPB7vEtytqoUsYGQpVI2568LlMzSXJtYsbfHGliIcKEyZuA0D92lfjv4OYL/+l6g0CYuJ3/wEV/On10/HLx98/fU0lVezjOfi2JD5YQAo+BQ0kJiTSCbHgWui17IXmRuhzCfFsXqxJ0qPXYJXEmQjOKXE/zMA4+76cbVYk2q1vwMJWG2JcZxuSmOfl2RlR4ArcegHpjzq/XK4q4m/1BWg7TyCkwC0ke2Z3NPNExurLoYpdry/K6ubDmmT6JUEmnsJ8d52D3hdLkY+rkmSC8R5hw3ZTTGU005LErXAFIIfVbeaor8H2YnHYf2BpWB9FlFKclncK8L1+z7XcExLGw9NSRDf4x74rhN69fDx1izYFQS/rmuimfhXGnFSjLjkOHVH8tJx1jITQvofr6BctX8kzhDI9m9Sunjs61NGc1VPaIrKflRdc7qgF8AOQ4b9PSSTpBxPMb9fA7VEvpFfJ4I79gE6c6DqxFRgREdPdpT3hwYaVbNekxyb57/3scVYvccKcF9TrZVmVNUS+9WpxU62WJdEMOnsRyl6ugGZAegCmHaRbKoDkEZeRVrD6tEDxs3xDIlldC/GoqRuj9YbOFth3FYDdEEsjUWrswHGrNN5XxXw7K/p+6wcnztkMRGJOshx+Z3/IDj4+yD7SVx98kB0efPLZJ59//MdPPgvlkG1z9vmPq+dV0y+q2SGecCAc0tanIw8hvSs6cL8885h8fj36gAh0FrB7dJXqy3nhemDVfrNGQ2IVU+2OfbLM3xZj2ee3I32CPvxZyN7/idMjkWUQfyLxdiKDkA/BnM78BSQ80DcmhOdQPuTN7GJMVQkmFVrilE4LA4pVnvLxHyIt168JpR1JBM2UA8ecOvS2XK/xg1riVzjtooqKj+M9iKqbOTBzU5wTpm1u2pRrmJFgdVmutjj+sfDtSNl7L1/98LPbrIl8zq8BrkUdPfkqroHeNBqeaRpyOq2+8WCqAZxkvQtPqQwOPtF+XHK7zWq1qB8+BMQpdsC0zqt6Wlw3fZutnRSCMYJP/ZDrAFCfj9/jlmXrW8/ow7H1wAqcWAVfvslLkvECfJr1Q5BDaUgreHya9vB7u+zzB+4ZCUKsegmYBiOyx9yjFJeH6A29RotD0XUddZFVIZAMTWQJqzJbVfMSWIE+0YajMTTSEp4Ii4idHzFhzYHpPDQFRERmhY/biid9XpzlxGynLMpM+WMtxbFrp7NisQCqTV5IbTn/1o/rZ9TPFSshWKjYE4pgjOY9nR4TzKmzLFLQNJFAz2OBRO+GgQcdA37yAKQA9Qc/tCfKDDDlulJDBzxCgxYK+HLDbL5Zra3jJw7gPrcaoMgyv+4HHxTQgIiocawYEUIwehAzKrjJq7cBLJJTy6Ie04uLFS/aWbmpm55CAdthskH8ZUYcYFHkqrJY5KzR3FyWrKFwBCPb1nKOLnR5IbEpLLysiiuT7q5W28U80FYRjIIkQeEt+SXtCluoumggufXvgyg8NH7npskdb4APhg5BOT/PF0U+77um0+31m24xkG8s4dSPOjlw58Sm5thOxWI+FGWsICBxcDq6zX8n9jnJrgw3KEtgTpWCJ+hRXpgEGe5k6g7KKsPTGvraV+IXInZZfz18HYyVVlBxQ1rGaXLk2XECTBAPxZN3zBavSETNTSzO+2W7NgNl+ilaIXd8Ri0ZAsnGOoW2/SEWMbHaoRqotsuCjmDWY4c5GN8ZJI+ibyB4fA4ezcfZYrXaRJ+dYOKbFQWXe34ksjyPqTXCTDZHHW4BHk6wfyEQpTv+mCRNOhv1pTYVIDHorGE67Boe0fINTgRyIPNKR/YSwTmeUcFFHK7fO98U60W/99f3j48fXs/L87J5eHLyr3/8bb9nWN8xqBia7AGGpmcIq+gWmEcGZlp6Enl8zCBOTuglynDhZI+CeqlIM+IzHs8OYfNlWVyNsx+xx2HMoDeiWc0XNzULSjBzgFyVGzX7MLML5SKAiHb95tJpwJgHby73hT5xb2TULBbpWv+rG4oN8bhjbbs4hMkfNi2eKfvX/KjajnXBCmDIBVuoc/HUl/I4urVeDvms/AIE/ij78dVPTwdBC8cGLzzlpJ9g3Oj3hiLohG8HyjXX+aJoVDvW5A1JY4TfL/mI5+my1G8how3fJrjOYOTSAyK9nDc0LH6l4hS9U33v/reo0PqqumH9mkxnaM3pW4eiERAE4UGKfp2GH2VdOlNC9CIx63qg42V+RYevpr8oLgsWeUTS0j4dZbZMJojx36Hr0pH9YvUKzgOeUU2XeUXnClraiGPpS6HBisdBQft9bD/2Zf8S0eu1jwq9LrTdnJZEkXxDkV6DaUE3KJ5soQSJXPq3mvY0EeSHD5vVf3z9w0vluwDMcBkinX1JljHDzuM5yRuz2XazoSMlySLb5SnIK9MGPkjZAaDWbYMlKdYs1ojMclYU89N89paAifmwyas5iylmSmhgPDIRyPTxslK5GZNVmzXeU26nxabSdHrSkteyKv65pV8PipriXI476F7ytS4KZoEmFJoZhisx/ZfthHIDFJRH/e5VhR5qoDPko0LAOaTW8THen5zAQkesxxl2UVjkeWrLF2LGFLz4iCrtxWf8OuuLehhFhlFxJhrD+6kfP8y+KRYlm76y86KCSR0MgS2sYjG7LIlprjYf1k4/CSsgLe4GVvE6h8z6/oeAZLrB+eqqYnyhI/lbscSxeW3Oltbsux9//BOrE2cXdvT+uTh9XNfF8nRxA0Cn23IBi9C2oZVuLujHuZ7q0QK6ywBLwuKr1eat2u3OSHCt2cOAKMeGAbmecDvE3GjhnJGRgK3Z44AVO7V0TnQPl6roBAySTquC1xYzkqsiarYoCY1Z/yzOB4UI+zwhtWgneJ4wlZig7EfTSojeifpdZ87KQVuQujgv6rcNUY5XNqP/B6xeSxpsDW2U9MmexkEBNMKKDnARHoOtVX3B54zaTi68srRzgspXmxLikcNg1h/KYPhT7SvCBwL6YxjqqdR4T5uZuokOt692dOi6N9S2ZH+87xQ5tEs350WTFdcEHVUeasGa1uJtkZ2fc//L6kxWqDSnC4YyLxgTWMVdNqw/aVJNB/Qb1gtsD14DlCmWa7znj6h3JNYskcXoxWDglJlj2GLojbBjsdjTMtABt6GzGXBdpo7kzfk8EB5Wb0VOv3nCfI/3sek2ZZh9r7LNuNpeoJorNhvWE7hpLbxqlPcDiQ8Xq6uXq8a5yXhFpKh8e9/Y8mCPFPOHvUBz8T2tUn5eENjAdtHcrJmFcONUer7d5Krm+HwvVjya2kHYvx1BaNAk1DLyswagljE6bTfvOIgPtBO/LgXpSETa5FcgaaqpoWP72Uoq7tflL45UMlrtE7Oaiy+QDaFXw5RdzPIp8GUKEPBygDRh+mXbKEceKU0n/cdP6K3nq+CGf/xkSryyz71NtKJlBcII5iKqPSWmxPhHtJPp25a9I4izzYql2Fw2tl62W75mRwIihStCKKIoM6GWACQ0hZ2oRBygTZZutv4gdJjIK6asMfVt79CpeC9EG1U20/O5F6ocxgBTFnltWs5hNh6PZRF1LDmDkJGkgOjwNWOc0bJ46vesJ9DgGWz+K8CTmUynMHcD/E7JropLjJ8PhcTa2gIQWIPRLRsxNUXTQJSDWd9Ndr4KGE17mhNyzD4kw5goqyMJf3rQQaHLeaarLR4o9KwKHXZvmUNXL4cxVmHLRO5npwWNuRh00vpn9kN0gdBn/1JsViOipdtlOFqVNrhMVYSw0nlhMjhgqmRCXEj3cQTcgYsRPskEZY7+izbIeIACEyxSz5unl/RG6h0fl3SwHSqhOOMRtuin045EPRcttWuqr6Q63eldFBXIvefJWmSSFJEu8MQ504buIsMtugtfNOLG8HwUGcfNMFtue3d5otxGlo0eq2iZ0KfQeNdCcM+o4xUSObI8rwgHn4vhSZiabNGfwbxoQ/6X538isjW7gNgkClA+U4CAelHmx0Ao6V8Vp68GOgUXhPjVKqtvoOjPfinXzLhl2+ZVxb5i0F5UIowt15hJ1555jrLch47YoETuzEWqpAoir8J/iMA8hbndXFvXKDR32//8FxYIBkOly6+ePckOvvj00BxAIVRtSf4vYSonWXWlFidiXdk3T5+9ePzjU+zIgrYXjQEWO5YxiUARWGicqeyTV08+PhQRC4/Ub4ASxyZ2qaPO1uOMDnVLeEdCT5f1zmlueujwtoIzLQ42D+FMWywBP89+WZSnAKONMx0hpvAYRHJEzc3g2Eqn/aF43maA5wdFhDcicvSRJ+gHpkQq7sVUqM6e2BEjuyxghjCH0xoSInvRZf2zBU6SpzfueEnUjMWeKTUxtXUMyYc2LTRDj34zZVlYGve9d3UaSWGzBS1Yn4omMhiz6+Yw1UpwLyB8BAaCa3bpoK0AeYMAHJLMR4JSzgYeEgwaEikGAu+T3wPvk054BBCrMXVikQOLA70RmvNfprfIrb3x+S+qYEqkUq3XEkwZok6sYr0W9VOb+VGhK0OtIp9kvsM3Js55OF0SnXXICXXs0kQUZ1PK8Flk45U9O4PHJb2CxRv8n8VBbKoAObLM3K5/hzh5/ouM36Zfppw/8cyYk+wvZhT6vqxKOgMSyOtiPpLPsplgQKG9UFyvaTeA3C/z83KWHZxln58Os2cvvlW7Pe2MdXlGuxBuc/BAwrQeH5xAJUy4Q73uT64Pzgb8+jB+/fmpvP4keq30fbaZyWCO+1ByfzZ4yH8/GZzsYCjvS8eVOD0UwgSRGYg5EkyVCS7OFqZIJfB18XegeP/gANM7Iu6TtvA+iyXr7GAyYr8JnZuPss/lWZuUeVYeD+b64+oVDcjZouXUsKdT/2I144lfOHh9kt9YDDscT1iBB0Pe57Tg2uEhy0TsxkyiVpMv1yrIYwPTdE4+mXz8yemnk8GXTCL6h+6X+/F5680k2Rk07bzLA8iKNTZtA6nq3+dqaRDKZK9FrtjRiuPN6b7UNvSbHBplLx0fK2B9hs304OTEbTEVGfREpA0QDh3hX+cZoqQrHdIQzD76gkENbc8e6Y9AqvK7WX99lH08eUH/hsP3j64loxCz+ZSWcBNAMKIAyTtzg/zXaIUPJwfdK/y71rrY71hteskE5Ut73Hr6Eq5vsZ+s8C0/gjF4uDLqVpcS1OiY50/+GMyzdsOm9f3sacVHyBlmkDaYyEqrzc2en8U/Tj4NZzHuaDRGXQUrkL42mNq1kVtVG2X87FowcVa5/+CuwBf4KYQhL+PsOZHj4j9mr3547vwi5ZAIY/12WZGM9X2RVwcwOePHYY9PpOwNyw6aJHsx1JcsbVpbZaPKMvgOvZ99ozY3+p/oEyThh9ocqJBY3SDcMssiGrmdkbxO76BrZyGSD4UukAYNdkfTcDDNX9HR/vFk9MXJRwMYHntP0MO//OUgiq/ZsqDeq8pZ0cv+84g9bEneXotYIt6h6M8mr84Jk0QfB0n+//6UpPYif0sdmgmcKQfH3Ey5etilmynXVu3AOmfuqS/BxUbZkpCVqLS9OzgJTfpqYEZFPmwJhNC/czL+FBVojlzjXOhB9imbbNgR4WDyVzFFL1bnB5O+FWa8Y28QV/0BqkDZD+2hOl5uWEFP7WQHcnZ0Lz/Gy8Pxp8nrz/BaXx5MDNofGDQj6Nd8yMn1aF7NLjisgbhXs8zX6iVSE1crxtlj8zfFoeEyX5TipIaaBAgS5QauKEAU0ck6R7HcRcGcFhf5ZbnafInwJzZ3s+JUq57TExx8tzV0q6rCxxFgQycup//mTp3m9FDO3nL4DYfavCiXXi0g1iPaZ0HPavb6oQ3HUUL5bMa7j5ptoCERMxajmhoauOcg23TuWW3RrwXJzIg1gzZ/UxS/WGt2cOQTo8dEncUpz19i07u24AI3eEIlidry/aW1nphH1e3/EepOzRA3MRuhtnGXO4UrOFCU98DuqmkFB92bxUM2fxf74CoaN4TM2e9918Y7d97MpC5JBjhQQtsNw8TY24vUvyAcNlyzdFadB8ocOrVz81/0Q2Xhl5ucUnnbca70gIHlp3XQeTqx4NlBV4WOewHupRFoIpreNpm27DaXtAmmM4klFZdJnTsHcAAVtgf/SAaglAd99xAGQTnuYgex7M8wVW4CRKqN8fGOkJ2o8C0o4csALZIm2P7Yahb7QGQBKGqaFYnmOazE6lJNB0fbrdckIvhNxYZjIri8m0frFX1gKDkJJxckjsNLeysEQSLHXCwHHca4xzUOU8S681MEBLI5DpE9+eZmLKICNzSVKE5FuYMXiDYp4fLb73uUesCrM4B7TzH6wnncJzAeZckcBRvEFHpxdExGXDXYMG7wJfQcK7UIiucyvbJAxP4Hc08av0Sj5XK7zD6YD0JFX9y5hECZqcWJ4uGGdII1C3GENuLtWJfnVXmWDvoPPDXD7OAQNYSld2xBnQoBQYe8/ojhDjP9cwpvWZ5kAxV5p/iKfDok0pqMb5B0hD1G3EmE95bAV2dbRHqpTS7opXRKGp5ISXlzpP1UNScP8kh/yEvPKY+y9/z+lY+ech3532zU+i0UpMDjEDPDinuTprQ1MVEoJWUxZyyWX1EZzqH94V1i8ZT5guRGqz7OfsZ+InERLJCgQlUztFjnTd2wx3rGnhYrETUk1p3kUcBzg8uvC4QxGEc/2yJ2wgl4zl0jv3p7lW98p9S5xnNYDJXIl8eXVNobWmeNSzg3eJX7nVT43lF2CHLEbsmeXtn3INpd6fVEKKSwI8ZskqFDgu6qBvKmlBZxU38fnGiP2hKlKzEQ2tdFtQcZ5BMhUSy+1V6w4vaJuDUspI2zr6HyRtID0X8QNq/EwYGbvijgCbGyQ5kOzIKtDiajgwOqVNJSwRhOkvd5vq7HAc/ZzVQgNIvjixHur4QIWqxgWQOV+rJYTihw4+elZbHZz4kjpx95sLZ/w9qHXNuRY//Wk+Pu+net2W97tzKADdhEX08ATsYP8AP0y8Eb+mYs8H2MLXCUtgBfm4GZOJ7kixmHkPMq/mfdPtBdVWKc1VMpq5dGEo6c8cmvTlw4xMtuOj/LXstP2Y+qfwfN4YjY+bNnE5Z/EROLh1DpzsHQAPGDRp4jErC8LOfbfDHijy2gKQjYqKYaVP2z5riQuHgDJK3IIECCdIQtIDkRFzpPR1DsHeI0b6u7KU9PaWeEVV8//T7T12lF2ax0PFKyqutrh1WiuLP+y8fD7OVjTu1RrYgN09Dh7mPO7OO9IFxcgU+VfsQ+K26dhn661fZ9j2OCSA/BHIu2/d2q2iz+nro6hark18BW5hDwY9ONMYiJQtBfFnldahSbAH7bAypF8Yjus5ET35L6kNd+EvdR3TlPJI3riOPW/VLgvc5Cz+2SqCtB2ds644vtOzihwdf3SOZRBm5b8a4W4zbdsmnj7rmrE9kofG0N3lHno646JjLy39+MivqOya9jz4DlxcCYZBjokAbYYptNSRRYTLHZ5JcoVrCPrI4IS6Pd/2VPWR1G/BB5ierbijIRfp1TR6GnylWJIEmEEJJT52cwDp8W4jtRVuyQl0cOIEp4rrM/c2VYJg1g3/lIPtC9MQhtnMT04e7qKyol4kPPNQgM6MIDzkjywFtlY6KFLnqnDbDwVcWB2xvXwAP2PrBeTa3w1BdOYtddz9LYQ59iyMefXKe5Qqx2vBsd20R5Pp6LPRrtIcjEyBEdgkpRqs5C1QAfQ9U/4s4WkVtGs9xITqDjvz4e/Zd89Mtk9MV0dMIJgqaiwOwoPHXf9VNnqb9Of52yQjQuF/bTXg/8hMrhXd6rDMBCHS9e3+Rm56g4oNPtAu4FrNcSbw6I6eLarJ60AUptLrPH6zXMlDMIkH+WbYkTUIo0ggtz8aLg1sWgn2n8F1TcHCbe2xPt8hQlp2ZKa0VmRDmvNpf7aj4l+sahL/6Npo7oAVJE8VwRW9z3sx8qQe2HTqG+gDca9xSy7WpT0hGKZApJTPAajrlqwKsfEqNWMNB5LtQRZnVVmXITEFnT7oHigyR7EGFN/d8UDk3JsqxRGgBv0CoEAXYHFplAIoHSdXIcIJkKIn4Hnrg7T6S7I2gNBNtZcLww55xb9rja92RNO+P3gxxOU5fDKe5yQPbDuooo+CI4rSnhiKRx+5wtS5LvrPlgI5PfRUYJg18p9hpX6ePQOyfawcGQDmFDUkr/Nt3uItAsqTDriCSXDmrDwQHebn7f+a04XzVhQ2U1W2znBfKwBaLpKpOgMPo1z9lFHZ/dcRmHuwAGj8CJ2HgacbdwnN6Av/SL8fkYruDrYsZ0kFfjQzEngWwCWJbRDI8EWZ1T0uBL27+1+ddrjgSeP93vqdQb+/mxS3RPtUxT9tW6bqaHk8NP6QA5mhwiXUdvj1XlUyUNDkRMGYa6LkdBrMCAnZzhTwIww3hCjQXwDIW8h18QcP7LAUYpUQIhApk2dFcpCd+1G3yMd3SPvsJ9BYFAEm3haw5v20Zt5gg7dTGWsJGEYYWDc8qvsAt0iLTHofHD1zf1+Buq0B84Z+AgLGP8AXalpo3xdcM0MRy7Dk9y3YySJk1Mi2FuLs3/EnJcVkNYXjLbpWowlPPmk1VFaCp+ZYviPJ/d8CG23C6nT6uGZMVXWEiJOmI5hlPBaLvOPilOyKa0suQvEjvUaKygRA8arZYeMZJ5VzsdzZRmTlExjRXG57b44j6088DYpzC0kl8QCNCDryx1FZSnOgH0rteeBNBD13EUeVUSruLL9GAy/WLC4vSUNvG07jnhmOEFTJGadPzQvzn2P8FEpBqrSaxBRRwqoUjwGhEOuaXH4+WnrV3D25ZmfAN82iAOFnr2OnFInqmPIAjjj8+fPQtfPMwWv1wN4c7HJ+WKPYDZPZ6z1jHPTWOuhBPhu+Y7GpKsjuRPw2y+LgmFP55Mbj2Wqiej+YyH/YPb2y9XGpJ8ttRIDVb0JoknBAYTk52hHJERxLxzl1gSsXeBep5vN4X1iONb2dw17zRw5Ztz74DGI+yKGJA0h0fhDFlGKJspetYpE80yTxv9O7RwlSMMfs8HTPiEZ/SexADa7mc+u1k4gxY2PqbOxlMWluJ5k6kecrT9L3CF68GHFGkcCA0szQHGvJ80EDwqps5XY8xWXxD08OFDCcoZcnWnLkP4H9B4mc9I3ipG8MfjAEGJbgTD5ySbm9WsqCVXAp0dawkiGCvT8l+7AzY1AYtlxrRNAbY0e3vOakGHelD7TOH3F7ErcxpTA4hGqPReP//+Ty+ePnk8/fPTV6+f//CS2R8JEGVjihy3IOdFc6/yLUlsW72tSBroGefgAGsf1BiJ6V5Gb0vjsSh+V2qnkP382mv3yroTJ+va1ZlbE2vd1aMwOdB9u9XSq5BIMLVFPLLl1HgpETzq2UWxzINCvY/V6sZRBpyod7pczQvOzyIYtJ98GjqjVr6Yvi0ZrfRJSkIA7O+qPdhHFYEhJgSIAlMFYJE27wJPgmpSgMiYUXPyko3so15chCMEppuCc78x5emx62o0HcTm4P01hRa4XpnN7V06pxAGWe9gMulpkNxBZ1+QwGtJvDjq00Ex+mPWvyjhNs5pAFUfQb/YFwYJ8M63+Wb+JVYf7vcaHRMDZ4+rgTbLyaObm+l6RUz3RqKskbtgikLTy4MvYWHMfDeQj+oq+5jFWNZhv/6G1cjZlVjj1XUeFjoSHGs/1+sifzslZk7d2ky1DeIVHHmQRG6bltEm8myRn0/D+rUF4kkmBlpfpGkGLB/JXt8F9ba6Bt9SPkzr5Up15kfKGtSCqtMwdQlyHPzkU1Ke/femFojersRf2CoT16Mj2t21DrVWnS/XC46Xh8h28UtYJ/omxZkTCLltzV4rmt6zDeaj4o8+5BMjzenp6todhaqcJ40oZG/gMzttl4TlhKZdS9VqzPU6rfg72mZEqguEVXAeASF9MnOq/WjPti53kNSXNv1WAn56OJZoAWVWJgdhZCH/igoFnOQoYCtSRnKDuJQOuzNd/2oJNT1rp6PeZqMB2UcusZGoreu+Lyfp0r36pheA2KzKKevR6+1pw8syZ0GfZ2MqZGV6NumFWphYSNvJMpAhhvtHIM8WW+RXn4G2hJlqe/kCgtHNtD0qqiUQqxUtI4J0XR9BbW7p0LvyJURHEMyoX/zGd4BwtqqKRUc3QcVpwiBVU49Pb6aEslHnHEgmJ/KlYyXF7fauDaKwVBFAZXlPBEs9jACkmyXwlhlEOKpxbtPl/NM7+9CWhXyyzSASMMnncum6vDNT5zCNvPd5LyVXUCuEWrOp6KGIOo9EXfqNJz9NYyqzHyWEGAxbeR2Z2gRSY7ickX6oPbndpEjEeMKPvJluG1gLAg0KKCyANaDbvZ9+fNJjPQ0/suRs0gk1CldokRyI/RaOYBiKfRMoLkGtRpw6Ro1jJCd8mUn2PKejlzQMG8kp7pl4hnxZxBU5uofz4a6aKNfuaUkH1XMWSRCWzHmD6rFg90DPch2pXSpFBbaODV0+F13GvJU7rA4j8RWturP9qa3ttVjB+BqHWX2Jcz4JKIUq2iX5fS49EwMVbFJTLh6eqzz6uQwKci5Pj98+KEdOT41LEX69qK9xvrxGikBvJ3jGynYEWj6UXkr0zoizrVLnziBb59mPnPBzkz3biKkM1EEhvEIt0dkfTjQPNFUsK5C6xvvRC6tPYoX9jlqLbGahYteLhw/xd8rTpVFjFfsSHlE7zIumkeZwmLEmcQrkKSEm95YSDtbTyBsbsKS0F0VODwPrYWQ9Hlov68+gVabeQ50NxB24yjI1U1TWsDf5j+NcSsS5RDnLdESDQUg4+PhVG4wIw7T8cTnMTgZBBWOkltSMRaBfGZEHotBkeBYOPJ6JClgOxXHb6RjK6Bv7kgVvftvzv4IZfH6WXUFts2XMgM6bMFyDvgD3UXYguVwFGVxNvJq6k33QkUew5QTPI3Pxn4TNvhKdhSngNpLkhmC6MnPeGn4zPXyYq4JScg3txCvu7JF28BY88quiGs55M3BEXH/clhQks9NKv/eUqQAPyCLQkcRkfykJM9otBaMKhuSzBTtSF7H5sNoZWnMJTcTeRLPKV+jky9PyfIucubKn9PxGaMcvRpBM4fm0UPEOFn5OpZUrMRjv8VEQ8r7KyrKgScys5hZwud1qsSZwxjMhxlTkOEj95l4OnDbRYtAcjMST/hv7oF1T1Kp5hoU/+brt61sCAt8baw54TkWh0HLtr2jUdaK8Ymzs7ovxgVFzZ1D4E1w7M5sqiX00k5vQULOeQcW8XUqubLY/sLuOapcfI1Bji2MwE2btWclhHnRihg6VyoTE+r//23+VS1qEt2EtaQwkoeYLDnDnthUOc+NcU9BbbNaT139WBQ6BQriVWDMwhhVhHGYQQCEXs4HQ51u1E1+GE98Y5AOSp/a2ZsNQzvlXk7G4pFeWK0SWwRnwfUk5UiP3QM3BNY0MTKaMZqJy8TgApRE8aXMb3gnF/MswCCbMccJxqqx0qLyxgDMfAMRIz4YkkK3KGa6EKOTamsD3V0hKYiIgqvVYrDuhd6D3GYzcDEs2qC28DTS276gSV5dj1PgRRhkCNvm9wHC5EcSs6w4o0UGeBB1d5DMOt6hmnDbru1+GPo+loWNlmj+zxvrsZrwkImlFHV7hTOCLPWdGF0ws6j0Q+MhDLzQKKGN1pCFDU54XxnKwGw2YQvE30aDeJE6OnphqAopo1xBimPvjG/TnjbiXv3EKizfSibwhkfEU3up2p07t7jTj8V/ldWDfc/uMjURuDtwlDgxeOhJT2qHHFlO5+4W/n99kvMRH4OvpUkTO5W0GjKsu+rs5Mui6KJG7WIeyV7lMa7Vg7tFxQ5eUYhEl7vAdUUdR4R0xO3EZDn6Mmwi6I9zndUjt5Dq908IHHiiynbM5ZiP3F7CHes9nW9C9yEuLQXivT13QQVhKtHJhMbfMNiCcQmf9CO4wAUDSo+gpkwWGBPnI+9rIGJ9crID+krvWmJCc3B92ER/FPPcUIjJXDkZPZWqf2VHW1qvUQ780cUoTl7Te0B2JAECOQK4NN5Nue7jxFH+Hlz4HzxzAWT3STsbe6CzUK8b+lnjvt2YtyDOpZ4zAHjNvIsWEmEOF4gtz/hBVPrQsO8AfKIQxfjH3JxQTIcMsmoF7sQJ8XuiVYcrRQR/G0T0h8+Z4KLUfHiUz1R+ctIRHOv9zZCIyAYO2zHRYQ5cB2g1tqKeqgVkxieBxnZ4jhb1BmuhA0xmTDOb7YmNySgAdTUy8UiIVPTupuItayBE52BThqrkb4HQrtS6Bs/dwhHCPetxmFExWmGONdESR7HHFaZYaPVFp3k0n5fAiRCsnszi/duZm1/pQWt6LFJHxPrdehUBePp5qHOM0PT8IZGKNTJ95Y3J6UHcYbQQJ/A7tPoa22mxNsWtDtMNpy4h6MQ1a2NnWOkajDFrxRDFZyeADL6Vnkvdcy1DMfufFlNai1XQduH05Q4SNoNyynAp653qasu2vYpb4D/39X4Vh7f/KYiA9rvYH8Xprnfaqq1633cH2qqSLn77eufqZpB7F1WaOD4pSeY9VFXx2GdXNDRFNf4IJrtfyN0eE0q5qM04XORI1AhIfnlWAr4W28lcRmOkI9OeDTBLb86EeB98l7l+lPv9f4Hfj8fhgIKCec7LUku/sFNzxiZ5LdoTxXrB25pSYQ47TRbTbnoTEVXTe2ro8AueAiu5/yXfjwqKIYKaQr1AnCYy/2uJtsW7Ge97Dy3YDLRT7E+gqMRqq+7Eut0cArI0pGcGKLzlLh+ga+eeBXB7BPQh+Trwaksc69a4Et0ltxLeoRfPbDeZgKkvn80KFUAfI5nTouw9/GS/lRSUHcSmEuUf9G2l4W1Lpq6z++6bpj78X/5r9+WqL6KZi3QmROWbcwRGkrNugOIrWGrQnT9E+08zUXXvEp7INieh7AW8DChhccxZ3PuxAiGMrJsLC/VfP1bNFlOO4X7ggwsVkrl0r116zNEzmK41tpWo86Q4fsollr+KhcycSOs+6vCeWjtfHDZlTqeUfDrn5hx/UH4rYr6d5lfyHbCrg7VtWM6IRrB+xE6OeebMoz2TviWo5nbpqyGdGJP+CxwOuyZK0JA+qlcbsB9m0uMBqXSDR4CVCPOd63qcvcTsz1N6IBSaf/X1b6mUckeYm63/3C0xMVXT97kyOAh9+q1KbSDZ8to5qfziONGsxSgVfLO2sTxis6pqpSBquYOg+xoKklyDjGl5qjNoc3i44ml7V+e3fQyI2xPYysbMH6PF1MMz85hp0gTQpeojDtsrRJq8fa4IwvUSdd8O7ydbRHO6Yl3+KWP1eyk1uJyWu3ImXFu5FS3zFXcREYN2LmGjRDmJiQFJiYu+FmMS0BPyzWyq9W4r8vfQjEi1N26USlxI47TOxG+fqM0gPviJrpYfC+2wBtx637oEALX73JtDR/Y5dEExGbNCP2eQAhCR8YZ5z8mB2l9ipIgvPr0aH4gr/3tvvzoP9P+VYf48p7zrSK7/qhe5Y0Xn+rnGzId88859LDrkbTs5mLiui90AEKtJXQ8i2TSKB1LHq/afwaGbWHJ69yHcx0gR3BUtZ64EKn7pBpwvayM3VKlb/1eM97ZMEdycnsra57A7DpmrJnP1MfVuwPO1rW+8nqNk9pFAddohbkKb32HVF3S7guHJil1tFtCyfwWU9r/ge+0gDD1NKIG3gnhzJZCtZkGH1Gu8JrbBvzjEpnCZ9JwN8V10h15T4J/59bJqN+x+BdeCEj+rzEBsmApOE3aH6ApnHvwHGmdXop+f3wMvAyvjOhpgIg1/QJjXzA9usZOzDtOd81wHLmhD71CalUXLjPXXxMK29DSU1Q7SsCP+YsUAvJ9HWHsB5Vs2znCwvSC3B023WPi8Xm4nVDtEwr4itEpJsUA4GU9wzqIbMpUtEHq4Dn65W7JMT+bYwgWzZaDonJBY8Y5Jqju4dwnDiZFDsu2pKJ2OqsoviKPVwjne8aSEzECD1exUHNE3Eb1vR5Aoup7Te9AAqeceS+MD0OulmDoz8LZ+2vUzzVAQb31iTHwMYVzwA7nHkG2etDjxMOJpq6gUmmFFLw5TA+LvxcmQjuV0wlSqQSuGBN95YkP5AqKYVPeDeeJ7p/S493+UPtnq+QDB6kQRWV5onx0u6gX+1HpXwR99G03oUP7cdQj0uHSlBm7b9ReMZ1ObBES6Qwkoyrsr9Y0vsMqUoUEmdNUUlV38w8efUW5KR37ciN7Bq3akTHkJVsl9WNltBBofrlOdgYYlArRumTZGPx6wXHSMxavhaDiiJ1OeuqJB7dix2dLY6r5DLIgh+QoZbm/LTG01yy94A34N00Ox+/e0w+/rtOf8z53/hEfy1BwFK/DVUj0N3kcDKzPii7wQ4TpPbY3+QkoTxNVB5mzu7OOf45lQ84GkPs97X32K/np5PDy3Hbv/rbxGsAVhB6we9cfaN42elGwUvLqcFZr+5GzaHswcE2/+bFTtx1PCFxn3twXDCUO/4BhSZpieGeYotuu4utNYxNUnmdooEXiSQvdH6b9RWL35NpUmMebIq472dSL1b1ggKJNF+8j4lQr/HlrlTPukDIQf/oX96/uvpW/5/jn8q+tcN4VfWUQ84CfEf9ruFltfqxA7PmWIzEhVhgrKsz1EvqNQX6cdoj0Lq5YzMAkeV2eiqeZW4E29QScfGOKI5m7n+OER8m1X4/+qd8d5xWHDA8jhz0r0PJTsfuwszJjJ1xF13G29PzSVxdORuQpIJpJ+RqODd6UFECvaNArdIxJ5+dBVQBxax27lNQjdJaMtoLrWGC0LQfryZN294J+fZGw/sTebS3GRZ/0e28/sZZBLGRLWSNXJLoQyONYZvWrGa6rjCBopq1bFucpGguqyM9ywsIthLLaGwM9Yg3CqavKBrM6q/YLTlOkq2dp/z71Bu0a6TZsdhxjwHR437bCGr5up3vwOaHn67RJl2V4bBeW50kDo/JqASD8jWpukTwbqAoYuWb1szAsWao4eRRlicJJM22p6SLmcQEqF3DCoaQNeSpwsSAkqG9HIVUB446y2hGpOwKWErLiN1x+C9U2aYbotlvq6Z3yH43V/mi4U85zgyiDpAXG9TXrsMaCgVim4QdMzBnfFQfLi6kDbALa86FbsuQXu9XdahT5M/xGeTFzvwyqlCgzVIUu92Y9kHtSCaAtLEBqMY20DkP5gLie8j8IKjElVnIozjg7lPxSvo2LVQLd/drHscw2R+JItlcGdgTKlo0iABucgQ5r5Gf7HRoWKbBfgN9I6QV2Wz4KJhDhIwV8EqLu1SCWgbHpVEQ95+PQp6rEPgA6cr21bMdWrs7CiUkq+Adv0oR01Vvbb4CX0Kwr34TP0nuPWy7s10gHyxwyXs22ecfYolxqVKvHwg0Bw7DNF8TV3EuIvBb0hq7S/zt8VYeYAtAg9i0L07HFLs2Dpx+Ea+WBEzDxumXcDgj4/Dt8flSWiNIYTgW+A9xwiWImIc0lUVvh6DBrEPHUJ4Pqy9I65ZAlSpohogF/7sVRC36Y2aPMnxaeCc4+/CSRniqIowxDdDfZA+vGHtj8gAb9S7zM3gm39QE8WCViATtxyh2TfCCSqxKne81yWJtWQNe29usP9khdTt3qtBp+xnEDyMAsh0kJ4gtOi+WxPNGoLTUm/gEye/15ePqnfhz8PM6V80hjjQwPcGCWP9SbJexPo+7hTHS3BY6p7zTJXF3yWEWa9TJOkWxTTfxtRF7EfCg8EI5KCEV7Xrq/dQOLyn1mXPoq4KDSBkH60gKqQNcJfM887KNcXfQL2mOgOMJ1rCeOn8YrUh7EaYUAfWM3SJU/Ek8NT3LW0l8VlNDZivubgFtDhnU1Vyz3Du7kKs0Eyp1xxxcvW2V6u7umGoMR7iohCqDFQZxZbA98NAADM/1F5jScAvzOnaX4LC6cNxHGR9CPtDcRyt8q54rXxYt8aUhztLpyep3a0mv4cm2Ct0u/Sw4vAuSk9qT07ETimbdlMIR9RDYiGnXOH3aKm9jtvDikWLzj57P8bfgfR7ieItmGeVge6pKva0pTV/Ow+NrRkVnTb2za7Z9n4EvsDtUybdccDdROwcV6qxvku/bvIudLISrBXq03TjXfFVuKB8SKnoVRA5SfIcYi3RwHz9qHkPvhcmYIhWzduN33Ew3vBvqXAjLa8Tno9aMRh0WIwEBaeniI+NkB3s1BiLW0E8wvtdKig6r9yYcy+qmODyZdah++gPVF1VMzi2JEkPrJOzVbGZwcps4Zhsy8rnJV/1phFlNd9S4L01Oe7AyH5Hs/5G9+hscxvfvq8SZZcapcUggqkLcmWwwt9p+9sH9jCeIXZpCZCA0w8QCofZJoZOQmnzqs7TqgZY5It6hSgLqx1rQ+wAekeL3bJCR8/bOp9hvCfu2XmLLepEdTGiAz8DlmghOl2ju7tXdw/wvWSEyYYOESUd4T9DvxNIFJExMFQaztoyTstSmCwGgbtabd66g/94HEM/YfWARcOCnhHJnF30fZ4y0HuaMxsxEV+88G4BVVHM5Rzs0jKr+KB95LMzjNbnyN5MJHmxcLA0QgiJWl3gJoMQfKC1GGtZvmvmH+H33tbNbAzwWuoE/McXJrJFPJjI45juh+npCI4a1lymULw7Dv2RALNDreLEkNW2wxlJjeEKV1DE20QJTnsqsNRR6GSLnd1hvmfwliprF/wwEvN3NRBGYO5q5FY4XaGcZqB2CUTYZ8/hNvs//jtx2mAL2ZVOyFoHBwzcMl9c5lXjgDeu1aE/Wfiwb3bbrveYqxZTO5rN0wBZGYlrKYhXlICEvcjf8H4xCXuaIIO1FiFEyRM0tHsKinyzuAknGd35Kkhq3vbre9+7u7l5gKHM0uFAoVBrAmy2Up3mc002g4iUWoyP1RTlgGqatFHylTmPMdpzo4Nkm8UK8QEHmrL6Ak+RB9jAq1EASVTQUaOtLfyuk6CrJfiFddd7dVLMkJtckaJiW6nJw6yMGb/3c+Bo/8SpA/+m6sD+ITFS7ZrPJux9Me6DGX9zvsn6FCoMCIZjn7GvXVBiEAu5af/77jE0Gv3tJOGJYe4MqxBpTEBiva61wwHUV9PhBDzRV3YXRKmfZvYcWcv45C+WQtNxii+I5TokMnYOzaZeu7YpZgXuRNXL3F0myD3EYokcSVsdedDO2XrCFl7NdJnj6i9Jd5FvTkuiTCRFkwROqzNyGlRJyzbea6VnixNhSx5S3FefJI2Tm0Y7lXfQmYnShf3HLGF5G8junHS4buDsbKKSGv+2TKA9lyAOhTabif11BVKjDvrDUY8vfNijvITWx5/EtbepghAILpck8yRlnNZ0nGkUeDqoYdY9JEIfGoT8az0dZjwQYNZmwtbncByeN3ApVqGGWoO7J4GplC6hHM0DVfNRDE2muz3V3cAtpyM/doFKQaRVNWktX1eMcQxsUc3LPMIGrmKJacM6f9lOPv7ik1cPXv1lezj53GV3ldfPHjzT1wNzzV1vW2lZ5bKfrCb8LTS90rPyb+U4+0+So9WSFplDjxxmcamhuOiyB0iQHZbTxRKrWi1JDIV7NOeGjfO0IAcXZz2Fv4BbE8vA60xSbhuakTTJOHqvG7jaWYV37Wbb/e2MjWmiRip8rfveSL23Fd4evmJXxiZDSbXbFlqCyxlf6O2MLx+71/TCmTL55lh6Y2K+JWOl976Mr0AfHplJ7LplZg8xw9FKl8aC40qLTaHGW711Ty93Ej8czWXrNy/NAObi+tguH/U5Zn2XWP5XJ/mvwjt/qTZr0DX7pik4wSulhDijvxi4DI9z1DmWjyc+7yLfiYVGJg4Grex+4BHlSuldxnCCI1BtB07zKUbSjQQIPA5DY/vZRNZqwoIN7QxGkEV+nogetvp8B6cObBAU7PPPIZKRiSE7ztvs8/tCxi9qxwLCefoq+/hFJ8z6gs2anRAdpCApiq4Fhwi7hREJpQM60pvs6u5lvilz3oCmR+yYz/dsvun3jl5gfrOPcYWle6uhxe3+OFFiykmlpx8HszfXFNEOeR0pEldnjRuCOturslt45JMP8CeXfxv3IB9MJlLz2quz+9cIx5ogDutswuwOF++ZF/sOLAsUrviomWNjqrbnEm0706Y4XXAHn00cMwnaG2Ayo00jVPaVmzZLk3mmV018bdP3+pspInE6gTo86QAdAlBkcuWTEorfKOAxW8q8pv2/rbV1t1TOAmis8c94rxz1Pz2ZCoJY4lLFkVDR9aVx7FdFXrvMxg6+1Os9l2fcM2DLLhz+yGXzTRg8ZioV+apu8tv34b3ETo2Q9JS29XQwr2jCaUAecY5iRNJZinPmO9Thjyypu9z7+OaVWC/vzmv/8n557V/emXL+HtQZU6DBETtSsO5JFATLgLKN1fliCAxn9HW748h+saCk/B+3R+ttAF9FtzNl6eV+4fVt4S0C/i6c9MNe1iFlUKFUWBEbbltekRsHgs7Oy3q9yG+mjWaIdd3Fi6GJm70OqScx36QpqHUADV+LlIUZ8ZSuemGzhxacRAqEfvZVvT19xO9fyU890Jhgyl/0N8qH3eNvJuOfpV+0g3ZiZ5VHWdmNeEEv4b+jk8Lu1ig5jMaBOOLVEs5NMlmMrAG/x1vVSvF5U9Q7G+QehigUXLso9dfiAMUar/zcZ3yrHwh1e/5NPW6tm/Ipv25Mif7xdcvPnZ+lXYGisInvl9XNmB4C41yObPcX5WJOw+N4Tg2NCd8Ps47Ou77iT5K0pgyPx9yjsL3jk0DoDvE4mIM2QNnX3bCCTt+zr3tSXd2yHK1z9xZIGMPpTXgxXPddOaE/lFV+ycZV9rVanYV37EjkRp0kGvT5QR0E1A1PSpwa3H0N8ca9VLdRnR69KCLVVey6NkIJj+ipWcqP3I4Ci7fcMpEUjU8VP2C7SH4dOTA4zhbcGFBK1kNWNAT21/sdpnZdSaHWnf7BMDucDHadsSysH7J3mg7g+OAhnU5aDtbaf20litcPwn3lkGRBv06oDMZEw+2v7d70qMV3z2fYecWGzMDB/0L5DWM/GNM13nLPRuB+4xhtUtDdT2166ehjkr0wMkDUbgwWwCBZZHnD+U6nfHnnProvA5fNdA9Nw14YH9jWZXRejNS+HqbzMhjFw6Ms1KbcbjXiIukqHaVv7nVLzjtdjDN0p7n4Rp4OOB3X9uz5sLxvNOmqJrVXy3gseEOgIHq/GZ3lMzYteYJ9VhRzWKCji03PVmGa1yDDq4tKIsbRttJFEUOcwp/9or8RRhWl4o1uLU64haYMVzmYtah8PuV4ndznPm4lHQRLQ9ybTIj0LJ//jbbDMr3ZmsNGfR/vG78T3hwJEHL3teMdQeiw93mKbT8sZ3stsGgF5a4KqR3bNO9xHbZa58O295zajI897sbQh4Fdn/M6sZCpuwa3lUu6aw5nqOHZ3IvmKPDcY0V34W7OiNJxJU02TDrnsuytZN1JXjjLIJezN4T1tWzu0bsuY7JLuZl2KnHNbOVO1ibQIe+5c3dzUYaUtE2V0OdO3aZLkLsc0Zpk/F3mvpWKJA4gviUkLLt1RsKehB5jQXTfQ/o87qUq62SKhu8Q7aL2vV1XUd/ayS/Zsp7kUq/FNzpMp66O9h1rKY3riV44tNHYr3FZpGSQXy+KpvAXzozYWgtNbeEy3S+3Td6Yd8AC2QyJi2hw6Z9AhEu+pTKfP1C5mEl0Xi6InIpBERMjaAv8e1uyQx7fsMLGCU5wr5RMnQOClIviU0QkIbvg3J0WIMc9hRFRvsrVNZzvPV/MtgvfDZjCaDLBjJB2Bi7GEqyBJJNlE4enytBzdID7FRwceBP7fJR8HQhuZKDZ0GrbqlltZxcIYdWBJLnk5VaZ1zjW6drKMCz9cjs3TWXUw249Gt9xfslns2LdCEvbcRyJGRvftkHTnz2v4NsoUbS2pMUGp+UStwHhiMC3lRP/gE6obvg8XTdRbC5h0jkEYB/B4nhUvkR8L1LkIUxowNBccSImKw4xS5hn5UcWTI+EOjusRVflXnW5OGi9oGUxSSC5GxVEtDFGldyL2jqaDf3sQNDSq3/uY9VyA7OgFrm4Da6/Le1Uck9o++Zm/4GjD/Wx7Vdr/v6GPBa3dNstGHapr+qKKg1IMXTzin8kVXOF25eXJkbpUB7yDv50hjSyLtLTxuN15EnMRyCn0NM1SNpwSxMccrSKBxYewXdtB/vBuireod4Jz9+/pwMVFRvsMW6JQ0xKkdybGsIRWVUakf1Ma9rBRmRsf7UyJ/fazuXibclv1e/hS28YYgdLUkppw6rV1I0wKL2313lL0jS+mNaLwOnpyt99dlyeWOKzaJL6kzHytytI5aKTfg+3U7HcDjyBdp5w7wGGYgV7HACEI1FlW0FyUlp2UGJRRMiLKMOP88PhWTJESe6Pc131KBfsOfg62YZjvyfRQXg/KO8G5VHyrLiyQzW0DFS3lgjbyAvK+WG6a5uCILToFL3r3B2ELiT/eQ/fHRUfkTwVdXxHyWMEUFo/0yR6LqGKv4QqypLXuiHAL0ZnbgYHeN544urDAvWX/9Q6A3eqVrS4te/T7NwSbctDIanAXIXVySoeXRxczMHzaSiQi2psBdu5K6GWebXNF0HgYatoNHEtd8lgatlDbzemHrYxtZ1oIHDvbiNoFGXBniKS30MZy7NJVq/wr78GaL66qmqclZYg/2iAjrMOmhefR/4i0cjcNs7ebPKrN+z0Vvu7B53vkgNlrifUgf6zbES9GDxAfy4QNUnU7Htp/ek1nVYqlgolPkTTdZvzdXcqi3cM5cB//5wAjo5NFzsjhHvv9gOSqN7uc/9qWwvWRpwdISN2sbQ41q6y1HoqerBxL+16FLffFZzTnYxkkGwNBbI/bwxj0YHOcGG3adoc6rCLQ2k2JcR4O63yanMXV8qy8CZn7PGy6vvcEaZRvuWCZ4cI/mrnTjBBO7fAO/Tw+JYCB/NRcHFKx2KHDkZUyTa8j2I408gVvtirbrz/Js3pGonQ5SLSS775Ez7BdRsJJPKgzThAWrVInJbJ7Ee+XugDTYRskHU7rYWc5vgYJWFVnMWT6CZnEExoqve1MYi70o7ba8NrN94cH78JE8vBJDPpuXAWSWXgiBH14L6AwxSiJjhc71vGhcD7u7Op0P+EIG8ePjwt+Z7nq7q/2DWKntbyVy7qi/0nhcSwBC3ZREG+fAuKLG7n78kgCes46bwvf/xWMzcEdxjgnQwpQGLYg0Igw0y5Mr04lg8nA2TZhbhyErLKgDUydbyfjru1PaCiCaOiLnJoVWBOgcPzeY4jaZgGPouctsIdcE9a9CTIhrlKPa9Dv2v2ALmbRr3PvqFQYFDL1MmiVmZtN54HJ0Z3FaXI6pptqoIK3EGD2wkdmy7Ly9LuyKtZbVNUHBslekJJfbEC8Pl2Vkj0FCAR0hk7lm74GGioqrYNG0zQvSlhE0SgUHDMer63UrG3myAmsqEToMJ25Yj9jyBIFUbNMXpEE2HOUjQXOj6b/haCyBE8SpUYdnX/T7CnzJ89mwwG/2CHVSyMiXlwwb2LHQp67Wh9u/v6I8z0DWoeXEY96/f+ZLoZ9c+6TdurDl0dckgxf9h7J/1sqMgJ/fFZMPdNuRwVpzcmqvuGpPD9gevhwkP/QU8bJhD2nbaEs+PW27Ozkg0lSUqwge+DwrxHJyKSGbvGJTY2R1Z8My3i2tFgEoWa/bf/L7qylYmQXARD59BhsH0lh5t7IX6K9K9/FWXp8SJfEBDbUTSwfnWgRFBBEL0O3Ca5X5sVVChTuQbdsr6GU+DFc2SMavmrdV+IEWYC9mTot+GttxPzFrICM7v15HuxLvZNIvlNiTvriPgQ7/Q/TLePQkOB7sNQCanqlH1uiJDEnrkon0gUbV6LdgnokdRxKq/3dgF1xhCBcXxsagt5HmQfZRrWl/rEGhJF4/BYYy6gUYOB3mC3K2KwDPHBSxVwro/6HPRRG5PKprcLtryrEV+Oia66Q2Y6AY+icDJWCXLZ+Wo8Y03sKbuTWXktqHlymhKXzl/IddIsOSiEfcyZ6dTsnfSENu7Hyeb9Mut1zkiH8kD1Bm4Kui/oUcx5uYrt5cjKj4wWp4EhChYtnSL0qn3rzJIFJm6gN1SN6JH89UmAzBgCksMpGnxeQV0Uk2jNUUnkckeleAIcFKc65R/HKfRhJiF4YfwzyRRCM2IggyBtDdsEWE2pfJk49Dnt182N2SW87WBwr0bZ91RyGFKd/bfFjTazXykl65iH8BjBAvL1PkRkXFMVTg9n1CqLq7jHeCPZOCM7B3d5ILGGdddxZtG5GL15w7kCeRiscc+vpipN3h8IVYqheCJ/fyCuTgwqIvZdR6fv83U/MGAA6FCsU7atg9ryfT8EGhbZf8ZWLcY7+hF+khB8/g7+i9NX624Q85zuKMr8OsWe34aZmwUdscxel4Ls3tOYVk6Xxn336/ysXEDbxeN9ieNUYVzK0u6/a+NxowOxr7TsFJ9SndfbJYyenHNK48S3GzofjceivkAWwp0L36ym+D4E1kdzKBm+cM7Ad2KEDoxnc5rP0RQF1GM+U0feOzJ8iPtxyW+x+ZKiXH36/JvE/6c9hjOZbum5P2agP4PwhZwy9oJDEw/mdx6W5NzxJ+/BpadCyKNZcGfaSBYgukmJDXy8ULyNeDy/PpK5kTEJNTq96fPMyL07g6SMQKDDr0kYCCyb6uRykBl+J1Fm7h6gYmlFfXxbR/nsgRdCcIMgKwCHyfUtMrNeaq/n94JtxSsnpnaBtVJjnhPmptj5wrNv2wpP+JCHpTHdQB+nzWH2+Kcnw+wV7qnkO6LH44HtDiu4a4Nwm1EOSIm62HXMdyV2H+htT/gdt6f5Y5Ozff3vgqmaYQ87VZtJDNpRD/xVFW6eZpihfkAHdu7lLlKQde9y3fwvNDA53P901liumxs30+je4JiZQEukCAdyixsnprRzmBIAwkysrMxAooX2eZzBmYCNpl9ZfuSoFBZOwISBhlYkUIG07SKACn3fJy9IukG4iJQIsj2nUI5RxbSGfw94HfHTiiTswgU9z/rjQ9oo488+HWiCh9CfEf9Rg5tCFPiAeqzNPsr+jmjfjzgyss+/R/Tq4MRl2/HJLhREqLJ566cpFB1UXAjmKvoqKHTkBoxnl4wC/0WlLQLv+C0fdUxjs6vMsaav3PV9cBJqJHo6qKnX9ApGWZBnK9Iv0or5bFJyCpoG56Q7zlS7j1HRWSfxwUk3lztuRm93HTJtF8U6BU59xZJzp8iodLF1XbaKjBEsO2Cqa27C7e10qHS0CwCoK3Z974HQrrYgOIQVd9zCcUmDK5qIO4fkFdE2yzs8nG2QTlsSOznKx27bYNj1yBIbH8P2RR++Y1TLvCrPcJVN6ERkv/eye7pqOFqfQmsFU3SI0zu9M9Jj4Eko2P6WBCXIWTk6I6tPdSmxuvGsqCduepALXqMOdug+O13WLupAsJ2/rPMF9b1wXJn1exZ44I9zR8FD61IjV6b1Lg6UkCL6WzsTa+PCR3PyTiM1dhrGhfvHcbmtuJCdhsIgfDeMMH1HILJrGUp1d6Bvdb9A33YYyq5ws12BKbuksR2hKp1uOpo8ZcWuPnz6kMoq2B/Zr2EkMjkGpplevOnHNhoVgQuzl+TcFxP0I8fQqECyIXRvDdubuMs5oVvne8/Dani3nDs6yu1oIqnr7CjLsx2csEBR9EUljTm2yuoWTTTdMWPYc6NxVP4oft4z9/KMxQ9NfsZ+A5oxysi8d+aptzPMPHIJht68QHtOi5i47E7PBDJ7qlJ7vZCw9VRvV3M6dUe73FMtyZxAoti/tal76vGiukuMsdciNL1IBwQQASXBc0AXgkdeOAGZbtdWHgG/P6VCsvfDHuhujF7JhuO2w10lsLCd+JNsIfw0E11sLzXE7u0NOETgyaLIN+rKJWIkknEgAxibnvP2cqmn+094xcbhRMfL3t5FEKrluKokCbNchAFWAN5my2l3JTuYJijUqIXQNx/e+BzyTifX7aJhHlTbZeiEJHx3lVx7Tr2pGka4ooZTFIFVedwiByzN0moxR+Il7WfkNH+ZPVY/fvM9/7PAWJ2ylZxGW1SX5WZVQZIZLcq3BTvLk+S5PU1CxN5sLlkDfVnW5SkuWZth8m91WN9ciljBZzhn4Il4ur+zt+V6nnBy/8qxcP9KWLe6w7e5d/ihi3PH3z27Dt8nvDq4bDjiua33xkf9hw7W2P7YZoq+TFsq8N/a3LLjm+OMu69M9l+U5wVFHc/z7zzD8++6uZ2L8klP3N101FmiLo+P+Y1e3A2Ucm/s5K/IWQDzNIroiYTP5EE4fFccEUniyFTO3i+zt4iMkXzBHMOoYUL/9J2l8TG8b14JRZCcm12hIIiHkVigO/bcUCBqwL7dWOSYkqndu6fbcj8KiDQFYASs454M5yGhgyprXDOvlDC4IyMGdMv9GC58qePUIcgp8smuwQSxM2zUCLBoECFpUlMRNHHVfxd8bWEsFwoxVlUDt5nT/6EGbe7abZrdk6d5sGPjPOP4WBaf1itcvDqnX5XcsGJBd6AEHMUPbVg+26zqWih7dJ3QjswVzL+Q8A2azDh5xXdbIhgjRESIl4/F7/nOIMyquCbpS8J4p+7LVL8kaXrrdiATv/z1VwtG4MfY9s0Kit7EN1sHeYa3mrrGJb3V29dSM2Swf6QSt5Jmp/UBl3Pf3AeQf6SSeBxpq1K/Zzmgan+JXnSxjaQqlzbv32j9AFPtGjYfX7Zi+PJItnPwcWc9OP9qq+L7a11QxPpOMkWya2A1VyddhCOzBSvH/dhFfbFasAs/4xRH70LmYo4q14XSNxeNuYVLIBhbvhE/3vVq4WPw6ZS0JUZ/I+HNeQULPiFJjdtrq4KvpOI453wzu2DMlrbopMtB31fUydWVu3dU870G5JvfsIuic33igCdYKrzw1l2ZDydUnKuztxBc1DjWKAHG9+HAgBHWdjMVtu7fSoo1neg4JBODnYo/cTm/lltv7fo88XyyC3jcxGAaOoDAD9lAFNL1OwEEd+369SoEmXz6eQ5Rx/q7hASPgeDu0iXa/GyvKeWCSS0pt5pHuWx42odusod+ioatqRgm45ItU5lXu2lrGaS/FM5AS0RbOrscPKOu9Lemn4nrDWzvJTNtDvb3gmbDMFhp6pa0TQvsjL+7yXD5QVuDfJR2MwbkpryVm0ETsvnMLrqQ6sMCOd3qukVW7n/rys5LSwqEPF+MJD1TKBJqTWnuig339IiX5/BFGxsGljM3KR8PVaWFkNJbP1VIMmW3veY8Onua5h3Y7jctcmrkC/C+G068XuBwmjfBnhA3VwiyzArlBhkGFm21oZpKRcfhKvH2MAYeUEvEcwsHzxmWeYOL+68mwnZ3tzmfW4ubGHs8DdyYLFnrnbuMalsPp/mp3kEVwzuWs49HTV7847jQySCyuhlidJQyPO6EgoATtyRU8iRe3rirJjbEqQVd03FhpPg3CQvsL6Z/Dgv8JjgtmquiqJhlFbMtB5C7xVQ3cb38PVqxOriz2gLhaBFPN+X8HJfT5K4NuUlncRNoFgCO1SOnK1wOP/9bzhmqldfx3e+lT6CBhmTCx8KZNXWdXMsJWO5GEgnEgb8CGiZS7xBJsjmrh7ZgcJ0kJtDmXyqj00dhkJ1c1UpKZ4WRupuoWLknC8oZH9rMOWCv9gvu5Ji26KZvR2ielnw1nFEbzHfvLOdzZS9O0xMToGe8r/lZZ4hmJve3A/oVWOdlNNSEMr3I3xmO8uKXsbhB9SLMdLiCqayV83bQ5XszXqd9vvW/kIYjV7tScT+rg3sCapP8w3vWbNH+WATQWz9dP5nX05FxnG/O++714F1lCHYHsSDZHWxUmYrzrJIA+WSgKnuAuQXCQVJKHK3yVkWRMoJ6cZkugSKFnPQ/ARDkG0/6/ShLi0ZSQzBoFRsE6wOek0AcpgAx5NWV1JHK8CKltzcHPpPkMRcB+b85DN5yebxtuKyj9a50cxi9d+Ud5Y7Sk98cxGzr5jBhYwfGsIJ3XIbaeUSdcHJT4QOCXSytw8Kho0xOCMNYA1ZHEGmcj+I3fRQK3vzKZfyLlkMClX/UAvpVG+ijFOhXCdA9fwuUm7h/sYEa7iFahO912rEjZDlvDo5RyK/kzaF7oYvoS+jqRSU+ymjKMxQbwA/G9xywkXSbVk1+m+bRuKmEcIMZEBe8aEZEgOsV8pM8+/m7783B/f6cjs/Cdt2J3KJxT97XydymF/niLPuO/hkRoUASZf8typkDlxfszB9iXi6HPT202TWqXDqs7VTPPqbXw+F/fdwyvka86UWSWCdmRkOdX55NU/zxWZ/nuZTLrTOX32WcPfXlaQrpVOlSOPECefhXnMvJC8/jPSciTc+uLpahL899uB7P9dBN5b04UHvmdhpVmumiOGMduyVGp1cypPAd00EnleodZnEudZ6bsA5eTJ8Y4uKLvnsBijf9WgOeYihAqunP5ZxvagouT3hnbuhw744ErlbO8pR2YN0dEDpmO8h5eitH5sVN+Z3rked0fjDZAUq750c0KXHtjt4kfJCxwPFAiZS94hkXSyiHbYo9sqpkazpKWhfnbN4IdlbD0jvvXbne9Txfs87sApranN2VYG19sClGcEbmhBZe+GPrxYplYAVe6+2YDK6ooest6wvfzoea4EI7DckTXhAwL+letE083sk93RFOEwINdh8CfQnjHgby2JbhZNCaW5qEwNnTVWEJgjbcFJOEvOBU7FjKukXl6BoMpF3oUVRI50uOyIlvoDVCiw9Zzj0iuEtvkGDBxWAUcklr6IhnXQAI2nT+mYUegVHt2YD+OSoYpXa/Rw/TsXE/yiQVRThNeLdLRWOT21bOxENM/BHdoisxhXAokWqg9EndYXYgPuk8nf9+s2annV3T5kgIT9IulZVDjvaUJAOITHvGQuAWHBcLr4QPtpuU5whvq/tIp9Mbq2S/Q8rSIiMtYsYsWY44ItsxOPkRsjf95VlqJ5cDfjly0KYQg5MdQQuOI3K32+xQPNKjwrs4Y8wVZRoeZIdaN7oNt5PGi3t7sNGzqo28IaEUqhURLUxNcDhsfTMVuj/XaI3wCDNlVxh2hQ4FdKfvwhkgonC0jRZeuddSpTlojH2WOwz+z0rgO4rGijrhwDZHQXPRlT9uhIyhfrwpjnZg4C34d4d45ZoxJLlTwuLsw3s78cn3e+R6s0veUrT6zaRGPaA8rjC7o5ePnfcGkvHG+jF1FnzAF3ZJTnPv+aRuTEMAIzmAryPVaasz+Kxn8BHiHKkAyjcHb1St1zoY1HN3P0/2+hvJvCFaRBHPy8YH4GoADsnjrfCQSBQPgfsZUZVNErPtQiqiVUHM0LRZTfG1Jd/ax8NPW9/dt08nu7991q7HcBGtNP0xkcvl7SscPQ4m0y8m3FfcBzito2KPf3rSBqh61eljOPk023nRQr3OI8A/R+IfhlczBU7QtCrD7PXLV+nhINCOAw/tVKphdqaDPpObAYJrILyFEfVfy02EyIYhudXCjLs4Hl4UebMkoYxRm/XmK1F3i1feti5Me807wemmEfesINRGMhJD62bLfgzISeIvAzPyeMvxsO2W5a6aeHfjY3wstCQWNL96wpkOs/TEl9lVXKEpzo47JGC1uMs/9WrCr/zNhLfcm5F2oZ1uzKcrw9+4+l5XPrEU6J6JZbsG5q7IaJ83XHjJpd79eEv1zuOIN/N6LaN7ddwyyNJxMjqDOvkFtUfPqzNf11rxR8iwkCCISQljTtUmLcSrxeXC61+l3lGwEF2G42Mud9KxgBF6BlgjDTkk9cIaA/K4qoug8I1mmEcIzDRO/YbrAeFDZlyOnUIybFxc51wsFnQIDaKnmTMoEYn9TJr8LaxdgY0kSh5o1rAK+bUCtdswYHfsRcpuh49TD49ucJ26va663KB44iSeJKoBYIbXxXl1O3cBVfPYobrDuOvX5fjFjsuVWSF5Xv2dPJ7aRrYyOK2WtTDax5b3LWO7qxihS9Ekyh1Xlo/ZDOLedaw/f/bg2WSQJOL2TF1EjsBXOlx9Vojg1tlAn8crE4S6Rpn6eMnuR7OH8RgtUquTkjvQt1qB2ufTOB4WJXeqxozquKbadOt3Eywn5MKaGOtkrhMVjJPKt8t+VI1zFKc+dzsG+Q5XOd2f+dz3zqa7utadS9NzpThl5i7+FaV4jLNiuiGZItR8K7qzQoo+ML1jih3KQRGcRd4FQLBk3nIi8LfyMjz1hoPFnwgwi4ZlNeIrNpgciL/BLIf0tODLm5zcFNx4siCQDI76v4qFuIhuBIMmasiLfn0czsNJV5HoxXGHNhXvByd2mA63qy22K6oJWifxW9zwKh98yoMIevu238hAloLn/A13APgHm7fRdmmXCdi90fyWSbtabd7CAc9LUddRvzsnMlowGZifLU4Gd1rHPaWD+kEx+mJ3q3wRblzlQfR8v175uY2biBHwjpW+z5zs7Nru6XGrGSdbjFu7P+1SEwt7jbdOLwlU4oK35KJNrTX4vc9Cmn0QAc4+8VMHboZ2k5Y9hSvdPbj49B88+VE3ehJ3o0ef8BJKXTs05+7QrCV0DH5l1PYCCYzo2QesS2By6PKOZP1Xq1OEGRJ5rNUhpGkOP2UtW/PpRP5+9qloWevCqeBm+bSoGkkdE5yS+W7kdgd5a3xsM/N+9k1xBiRzGjLJKLo+mIR4lX2UTcb06g8dQ5by0s+oPL3aXf7TNvxPb4P/WRv+Z7fB/6IN/4sd8HUiXsvZCsKfV9Oz05v51zXiv2yW7MgvjW/SoTUV3X2Xe28ZHLCghfdGXs10aIu902W5a5utkbGj7TvpQQtkme53g/zp5D6QP/sdff7sjj7rinj3CLdJpNGD3zGag/uM5ovfAfiLOwA7+wenWMj7DWRZvUUSj1/II9qGn06QLjHa5fgO15KJfuTt7QoYwwtpQX8y/rwT48HsXFVnTFDz8MuigeyWR1gveX+B5z6T7qUq166IahQ+hd44wy2sDGq3Pyo7F9ZsVpYbDNlKvOzyTlXDsPmn2l6jHq2LX1blnHVqKlSK5yrO19m2ci2G/sfsfRWIx3JYYQsnUCDybfNlPxLjab6dTaXj0UGmOXYlU7PysQN40nbPSmSGXUCSYm2IYStAoQCSdNosvZAE3Bgcc5w3TiwOIB27gidI5hMADb8wgIsCBg4GsHNMEbiPFIV3jSwsPHCWL5wLra9/sEZb0lTAACWlTrWx4euGC7k8Nl4oA9g26uCbD6KC3KJrC2fOq4tlwoDZ/cd5cvE3eyvzcsoK8rgS7FNTZ/l8V3aeuGy1OOxuhpn6KnV7MfU7V47p4E4PJi/6iZ3JNcekNgC+D9uB70w4c61izuzgyydz2q4SWCUEb6Opjsp7I4Xh0W22IC/wtaxBofS4yyIEvt9lDQLX7rIEgee2rUCOoN9hBjIOYWYg2mu3WIDa2OJsQZiwbjtQtHydJqF0wTpMQ3517mEcol2+5+1CT+Ic9KrN9TeoQIALrjuI/SrnjSkJpSwrCfPsQ0zzh+4+1GoeXZhQVurTpMbQUJ8pWQskxv+xz7Cl6bWcspGLcWj1/0QNbKwlDSZgh5a0yEliVjVpZ8bDJMZ2GM7G3Wat0F4yb/ZZAGRp4Dl7Lp/d+EufNJMrBzKp3KBrtqeJGahIcitLbXG5Q7lvhEGB+boisGb0AKWXhPc4iJ0BwQGx8AfJUJ7tnEMdxJ5LlTDlcOiv3DUxrtHwkpjVIspBFmZi7FRc9/nCGK11cnKLxrqlpvZ6pjibIXzW/LJGBewaFWuw9RGpHr+yfHwRclgVl8xtqmnzdJItOeOe5Yp/Hxv/VA6zhQsl92tRiy7R2XX8FznKbxfNdM5Ms5XrM1wQbeunuvAaUiULldtiQbrK0E1Bo5rC9vaDojjZR++lqVeF+FY4JF8RdSRo4r9/xikr2eMLODp3iD7r93h9goTKPZ1y95PbDVKduspx6mXXKyTQlRbYystXycmcuSKcylkLDUOwok99X9MhI9u8ZAW4wkW2fE2wmuGDmDERmbQZZLJ8vV3WfRHk3OuoFVgSvzJTSdj6MOOoOmYQfPU1YYLcyupJ6bw4y5FjQ/PHYCUlu0xKz+uIlHsAnorXnRfAdtW5KK6lQfp3jnrai6lsB+1LSEilC4HBSIer7/dis1BEmnwGwUlImoDNufTilEjAfHS2KYlrILuQ5L8zk+gS66S3cdt4a9fgV0fZ59ag7/erJ/j59QaXStKu4r9jgsvJgz8eZsjZ33tdNIe9ljj/fvatXi2dLeWqYwaKGG8Oe01aqonmFbgyYFtM0cCgX7lBpslBpFZosVLuVHs5gm8pkhwUkkXC3fnDKW/k3T8oPPwTZAWHZEGvgm6M9zSR923c+A6O+5Svkc2uClGJrcvLVRNfuejvevyfwm71bjz/O2C9GLsKI2HZZbHglMxyLp+PL3POGSU9sKRFOcZtnzwH5q9wG+AJ5stkj4zM6jciYu4DOGhvT/NTa3eIbr15ePRGzhOWJTlkgfzB50Pu4JH4mPBHobLmmSgtGUJ/W/C1V4wfEOUE0WRnK+3t823WC1rc/BqXq3nc1AIIH5CfH9ZuAjpwEZHgbHuFWwRo2y5k5k6M95RUL4wROjyVD7LkNW0ol0x52DqI9f6ynXz8xSfPHjxD8pHPJ71bTjY9e5v5t/2k/qDXPmf1AoV9bQX0SNTjv/Hr5OTT+8HUUlyWP2by0deLDkK9b2ifoAIHYMnbPio/OKSltypymOu9pgM5FYUSD2/SAWX0PJl/5qvJ8amHP0gH0yxwmXMvPqnqsbNnxgt+7reK+QOtK3n46Qe7Svtjriv96WRn6c/asD9rw7713Ns7mPz3f/t/vqA2eAFRrDU9Dxyo+MTp74+kx11IIkgZHEPT3dZgenftNv74j+y33RuMQbsNJr34vRuMZ/x/8C47mIxo2Qb/q2836Mx/qiCnP9lC1r/33mMd2qhZjV6uFDPL1f8+uzGdpf8xuzP7Tl17G84yCQWJJlBxZgN6p75yxMDKzaY4p8P0xhs6stliNXtL4hvLTHRWnBZzzRrutpbEEEiWPdZK8l1WoqUbwW/2o4xPWVJONN3s9qhBBOJZdYyLuQkAMrVrLMGsT890dkfAkYAjaJw/6/8HUEsBAhQDFAAAAAgAAAAhXBq6I/MxAQAAyQEAAA4AAAAAAAAAAAAAAKSBAAAAAHByb3ZlbmFuY2UudHN2UEsBAhQDFAAAAAgAAAAhXPjGJnmGEAAAHS4AAAgAAAAAAAAAAAAAAKSBXQEAAHJ1bm5lci5SUEsBAhQDFAAAAAgAAAAhXHq4xFx1AgAALAQAABAAAAAAAAAAAAAAAKSBCRIAAHVwc3RyZWFtL0xJQ0VOU0VQSwECFAMUAAAACAAAACFcufcQqaUlAABBlQAAFQAAAAAAAAAAAAAApIGsFAAAdXBzdHJlYW0vcmVzcG9uZGVycy5SUEsBAhQDFAAAAAgAAAAhXIoN9fUoEQAAhUQAAB0AAAAAAAAAAAAAAKSBhDoAAHVwc3RyZWFtL3Jlc3BvbnNlX3Byb3RvY29scy5SUEsBAhQDFAAAAAgAAAAhXDPf6bRfYAAACk4BABAAAAAAAAAAAAAAAKSB50sAAHVwc3RyZWFtL3V0aWxzLlJQSwUGAAAAAAYABgB8AQAAdKwAAAAA"
}
// END GENERATED SIMPLECA BUNDLE

// BEGIN GENERATED WEB BUNDLE
class CalciumWebBundle {
    static final String SHA256 = "0f4ecd5e6769ac7c6be6bc9a8aa9126b564705e1c908c5f32dbaec1a1864f746"
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
s3.setExpanded(false); s4.setExpanded(false); s3.setSummary("after choosing cells"); s4.setSummary("after measuring")
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

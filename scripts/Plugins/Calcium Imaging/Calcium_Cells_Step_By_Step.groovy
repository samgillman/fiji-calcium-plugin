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
        if(!(value instanceof Map)||value.schema!='simpleca.response-protocol'||!(value.schema_version in [1,2])||!(value.windows instanceof List)||value.windows.size()<1||value.windows.size()>2)throw new IOException('Unsupported protocol schema. Use a version 1 or 2 protocol file with one stimulus and an optional control.')
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
    static final String DATA = "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"
}
// END GENERATED SIMPLECA BUNDLE

// BEGIN GENERATED WEB BUNDLE
class CalciumWebBundle {
    static final String SHA256 = "8d4fd33d47345b92030c017a01db960b2d4fb5653ba1224995d8099078c8ec49"
    static final List<String> PARTS = [
        "UEsDBBQAAAAIAAAAIVwJ3CxD7CsAAP21AAAJAAAAc3RlcDQuY3NzzX1LjyNHmti9f0W4C4KKI5Liu9gsjDDdrWrswNaMIWm9B0EwkmSQzOlkJp2ZrOpSow+++WzAe/Kv2JNhGD74p+iX+HtERMYrSVZ1D3YlSN3MjIzn937Ft38QP9XyICbikGxlX7wtsuJYVl1RPx6kqFZJJruiOiSrNN+KJF+Lapesi4dKbIosKx5EvZPip3R/yOTb5P/9y+//9X+L5HD4usLneymuf/yW/tL/sdN/IYT4aZcc5ELIe1k+imXxQVyvknINg6U4zF7mx65YHuu6yLsiS6u6KzapzNYdsUuoS7EqylyWotcrk3V6rG7hYZqLZVJWNLkVzV5UD0m92skKPruXOK5u34OWtyIv8KutSOEjcUizTBQl/G2VlqsMtuCn9LdmjtZkaASaD37Z662KvC6LrLcTu3S7Ew9pvRPD8eGDqOWHmpb7AyxISGiVStzRZMmzzNL8fSUeyrSuZS5gZTwKvCwlrbLIs0chP6zkoU7pOe578iiqQjzsJDd6pNZruUlzue6LP3z7YlEWRS0+vsDlLrMjrOFqMLgZvRndqge9DOZZw+O7+bvJu7l5vE7K99R48v3N61v6/mGX1tjBO/oHW27L5LE3HeCzm3ev3r0xz4YDfPj93d3k7rV5OKOHs+nN8OZ78/AVPRy+G70aj3mYbVGs8dHdzd14gu0ekjKHB69eT98MBjTBBBu8mYzGwzl/U8rqUORrSR/e3ExmU2wHZ4rfjaaj0Rh/wxHLEp7MX7+avB7wlwy6veFCDMQQzgnPqtwuk+vxsDsZdqfD7qA/HHRum6YjbDqFZsNJrO2oo6ZEwLUQs8OHWwfWFmLEjwyswJM5PKLPNvBwIXqAMJnsVY9VLfdd8Qah44dk9RP9fgdNuuLlT3JbSPGPf37ZFT8Wy6Iu4Nk/yOxe1ukqEX+RRwlvXpdpkgGyJnnVq2SZbniUfZEXC+jinfgB/gbt4I9khT28TQC/12kCSL/GDl6+BeRJAb3+Ih/gJ36ImC9vX3x68QfxERG2V6W/Aeos4O/lWpY9eHQrPr3Y1XsYeVmsH6HZTjKYAWB8hS/V431SblM4XDhWXjju6bfD/mQq7pPymrcDNh9xGDaOn2m4gefLZPV+WxbHfO28nNKBPcjl+7SmLnrVvmAEXwDi1LAnaVLJ9a1Yp9UhSx4XYpNJmDX+v7dOS7lCJFvguMd9jhP+ZZeu1zL/FSZtvsmLXIp/B4SuKGvoFZv1YTYP0AY7gtXAv8mxLvDNi35dHOyveUSYyTbvAWLtAVZWQBckEKNtcoCPCXIOyXpN0x4OANpGA3wWrhpRttmmq81mc8snQyCrW2lY79yq+eECaNbLEsnJR8F7RbRuiCDJDx7U4c0Q/YgK9AgGsIOHMjlYXcDzXPejP5vgZwWyixqW3Z/fOsNMcBiGg14mN7XCGOgRTgEAComy6nCT7NPsUa8GQbHj9kWI1Yz0CqhcASR7k+Ee8AHeEiXuNY9llqWHKq1aVraHaT2k63pHQAqzIt7wMRzWHBRSkCGfE+ODpgQ8bf7lwS5RkdF02tX/9YfzTvte4xz6SCcRAa1urtavNiM59PEFGkJnnwR/h+TU/w4AZrieed9hww4DNQK6/8larmTifQKUudMCs5o2AFfbL4jUVkWWrqGj6Wq2GphxFkC1r/tF3vExjZEIfhbHOoJhHl5FUWhmpoPHqmkSjZ1U6Vqabgfw73hKLS+kEDweA2AAcp9e7JM0j0zahS5vUk8fGTcIhSfvrNRxIjB19EnYR2BRziGT1dOAe4qwNJs9UiySThZn1d9kx2oHczNtBrHNunDZJ7cOuM/IQ9JpnJwZDgT/zjXhIRJ+CaH+27Gq080jsXKJ/ItwtbeU9YOU+mzmDYUj4qYICY7yHRHMxSYtq7q32qUZHl0bs2ulBklaPoWx0E4gNcDvUrNN9WNG/Zb7JItx3BlN4tOLND8ca8LSX1Af+LXTFfToF/r/HsSGX/UTfP/H/LhfytJ9VsmkXO1+VWMv4BVIZmntUvMxzlMd8YwASZ8wz8rITzbYDcRNQ24cOvN2+vbu+/lTqDLzUbUTZopqB/oPTDDUBIfjAUMPb8+mWB1BvK9kBnCrf7FIz79692mVLjPsAegZyHeS5EI1WeIHg+5wOOkOb+ZdkIdgaqpdr9hsKlnT2prZrHZJvpWE+bw+ly6zfBAXHFgFaM72lwSExl6a3wMIrf9Yl0c8UF5J9BUyVNA7ksjLtukQmwg4kNwQPNP3CWK7jwz8XcB698mHniVg2jS4B/igJTC742PmSp8KfpQEMrxRlMD6gvgmCyMtXBK0y9X7M4hI1KyC7aoVMtJIK9CuscdDkTKeAnZyb3QgZqo9ECFJwIBeV4jT0VPGqRQ50JS9Q2uFYhH9ZMUaZDDThlrRNJHALIQhM6hUxgSfNkph7S4sB78GirO2pzRHDqEZl9aqGdBCyhCiIX9hQCyO7AE6n0Z/m4+cpDIs4QXkihDUwGJ/HDnbGgTlKmVm1kyPTxK2Y7CvhATlpOtgjvWiWfpih2BOpHgBR5kAOVl3orzfQfMTBML0rLtroB30b1DAh6d1rljXllyh92ItN8kxsw6xfyhTgJfH1smfnraj94Q83h/lifuGppAT4+v3wSjWHjpg+Obm+/Hd2O/QPHYWYzrd7oqq9noiSDoA6UXd0+3NeRU/6W//IF6T2cm2OiWKSYn3Uh7YvMbQLIqNSOuKAJwsS2pepUQxQvO2/t+O+0MXsD1Ltr1NmgHAw68quZdojJHCoKwnr1mINfHwqkF4QRq0tx+eymvvQahUTZg6+t0wPHibG1HJpgEU6EMyROvZcsdFRGcWpTkhjUOGqOSS0WygNYNDWWxLWVWWOWZ2Rk1Fe1W7mYWROqbuNGN9B0zmPjAANYJ7K7rbZJKai/HAJoFmiP7yWD2acVTP4+lXjmUmyQEnubMqQ8Ft2B9xXyCw9FCjTIFj5sClse8/vZePmzLZy0q1/igGXzVcmIWEHoyBUIkr8t8pMxdOE5h9faxcacOG/lH/ldwTdHNLrde3qPCq0TIJ2rBw1DRRen5cWEHs/xmQGzYQCHHVpx8Vagn5SrLZGnH/IEsUM5KtJGt7IlBpUZSgeMj74h9zAB342RWA1mglQPpAtukky0RSKwM9QbOGAbSC8+NjXqMNfCfZlg4dI+FAMXYls4xt3HUzS/QvwHA8Wf6aDN6bdHss0ZoONEbbzdG6uSseXqJZngaDDgVTJCJffXQvRGSgQPJp16PmxtDjWdOw6zNy4BIgD7fKlgJdjfbz1fRQkwnZL8916ZvsbgYstdFb26anBb/+9JTop/rtL9OtY/RgtAi22FA4sqoPmSadsDGoiVH3S3dmo2FEChtOW3ig1Y9apC+/TRWnPBKcE6TKhGySDHQLelYBqbBxZQMCGzBfRJKNgb6q2yCITFY7A4U0Pv2tlykZ29sxmgINBJi5QtQEzIKJ3MsGuCuxlOh8+02WBaqKeuhmWtWtQSiacA0o+shMv0ZBBVEYERt+WBhi1JooIpyWlwNWpXa8iqs32sRKE0M3gtkKRi2j1H3QSl2beP/3sm4pwGYbJRxbx5vhE1QHtRVWVx4os8UqHHHZDb6MAe8qyVbX/jxED7VAC0MJWpZRoxOjckyaVmd4VRxk3mOa7DNAy+uBHs64qkmqX6gey3x9wl4cHKdPfJOlJYcNgp7GpqeY+KyBwpNMnKbmKG+M/ho4NC7azqnezmR5iUrSLtPbHjF7ta2KGA/aj1hYQoUn1pVRItD6mMusxbT9XGP2QEN+UWQggT7JazaI2i98IFEdZ8mSJn9x97OTbE+4HRvS6dsuov23fB0lkaS76dYPoMMAD9gV6UpWtF6R7AtgURhmUaFHGC2RBXMNJPB3H1AuFiAnPlZaRitJyQPU03yDbNIk7xEv4LH6KBn29FAf221b1gEYS7WI9/KdrSB6RCyko9bhfYdcRM/kywAIgu8g4vTTE0fM8NX54WA4G75t9Hat77lmLg+kuSdt83KfNsqkM8xqNB3JZpi72d2bd4NAIR2PJ4PJum2oS8jM1Wg2Ho3XQc/TN28HGPbxSUR7jliMpuvZ8mbOhlRZrtNVHZFt3L1uQ63PdeNaU4hzu1cag9cYL9OAU5oTQC6zYvXe6K18mo1Ke5kWrRhkaaneaKqFmff+yxEAtn6MbE9MfQtd3lowYohO8kcKCWICDQtJYw5TRwkdub7ObZkCAOD/e4BC8Ax2lEk0LA0+3CcfrgddMdzguhwHJGqfW8TjQ6G191LC9+m9vD3HI6JoxQrtVvbJfBlxXFujrpL8PqmcwUH6gHmjWp/m5DUZWKeIG+waJswOqAP/9OLqAEADkj5CjdKfVmVRgVCYEtG+OmQEMbGXSKcT0JyBqnKMViLWZbIV+0JL7rrz69+KYg/Ik+adrlhlQNDxE9S+4OzWBcnu3LVItuhLdvqD/z2Qdg50I8m38BXSbN421XFkGQD3S2t/uWEf5obqTUv7JbyzvinlFn0jZ/boRR9giKD79Lk8y+OqX7RRDmV378l7aFhZMsAqOSgf0+W8w430mDUGgEudIopKkB/nQpk2Ss+iRje9ou+0QvAEqumxPGF1Blw2RaNu88SYcT2xSrW4iE2+vnszfzuyd8RwGs3qyKeWHiLk+CIXd5wR8HLI763pgONbHl9GzFsYhJ5yCntAsQ42ebPOHl292AIx0+yEtV1NwEud7uVTZeCRb5T6O+jJLhaMzhttmqVYgQAleq1/jXCnU/5NcYU99bTwboHu8IbQkaCc9mfBgu0ZhewePdfwZ37cg4C+WqAJ5JglJT6oFBcFgpZlpOg+g7M1GPevZbe42smk3icOMimYDb3PprUywsSozQu2GdkO/iw5VJJ0Ofrbbcxk2CKd1buuqB3n8NQPo7vAFGCfO9ohuGfnyIB7rN4/QlPEwqj/ofFmhraPYKuOFcxMRXkoXPutl+ZrggYe3niX2j2HIAr0j4fFItnUqqVibi/F7//8Ly/djVTWI/gETjdv+ej/RD96UVPcbV1ajNqKOdBvFzmpZ2m2vkauGSgI7968e/vue+eLmActFmFivmDbgxu8OJXLDRlc63UfMA/eRxD5MoSF3SHBmoVWhC37dxQi2HDli6FBKC2C3yAEPXPqo3D46L40YGY+HVuWXvUl6eysqKMZoE7eK8Gxqsv0INl5oo7w64rUdW3jXaFZIKdHWfpeakeK8eQq1b710IMtuwwGTn52FiROD9pAiNqlHyRsw6oS3witQTFVIndSTg48ZdFGgyWIzZtjluGqr/b8ZU9TMctlqnkvzIXkVtPG0KbRSRufEvtBIM6BvfPH12r3KUuhVM6rhXKoobuAPAK86g4J/OTMqtBOg4e062/SD6xcB04cod/6lMaKlUHme3DAnr6wqdWgTVZugV3WmMPQi2hQNKZKzC+IiW5jYqP2+EeEnn5e9CiAV4fxhpHDcCb/XsqDwMe4IoSG+1Q+9MVrcSgl/hV9hUDQAUbEUsKhoCMzyR5JPmYNhh2RgmyeFeGldrUgmD0kj+Rb3A09S74TAUqxxfFYIVRkZV3DjgFgprULdCqAyxO5TXsWUMLQCl8yCSPjnh5nbDszMKRaq/7eXFB14Ghkb9IPRfmeI2V3ypv2BcNtQ+k3OuLBtamelA8dF9YLdL5v/Xj4CGVr5fcOVrw6jxUBWrWZuDBrKJVrbUq19yASdOeQVrne3MgmhIaZ4dgXtTgMIbYn9tD9g8zVqfq8dLJ2x+jF7FsKp3dJ1dOef2PFusgmZf/ss5GqMeT0FEw69q/BwDNf2W+HI/VWd9E1Vp+6OK52vUShigE3NKLIHrCfyjeGNXYPI8/SYai4Efp7qYf1T4likoaz7mjUHc27/Vc+/Y0d8mklLkrwY+frrKgvy7KwGA2Mvp6sluEuW0YWD2Y5yjwwwxMGT6y40Jb+tGHAtakOvXMyiHCRzGip2k0kbMSKHwuOPaXa2vopw5n6xeTXsvjOHLsyR5RHtLP20RwrfBjt0ao3DyJEBpHhIAGln+wka1ItLjGJeeMY74clkc3m+lwV89YH62Bpf0yBVC/65TGPRA/E1HZbancyV6IZBMvN+kaOz8owDskbbuablSvYDERP54h56aP9wazzNH5sUlbUkintxY2mMAFoAVq7uzeZ4vadcCWYQaxoXQW8dDhVAlQZdvA/736zjWuasjaJkT7ytGbCRGmG+kD7RdU4NxpmOavFe2lR91iI/sktthIUmQoMrK6aVJCuUI8YgM1PnR1hbQn7vSIKR9NUZcOO3GxYlfL47CBT22Z3i8IsAQOcOKYM09ZaFGWTZkS+2matgrnCZIwmU0FHsGwotsknI8r8ZH/QElWvu1BY/zRv/cWE6BToN8PraLgG8PmlsYnH6MwT7NenI6o+/+hbUxcCRhN1tNlr1VgWk25MFu6g1fNmS98hk2NNbi3rJM0qjnVCEaDTMH0LYj69qI57TZL8hV1gxjsR1A/InaMGKtdpTTLPhS6ccJI6x4geacLt9r4bewtwiM9Ah4AH3wFpxnIBZg+8bEMnYulVewq5UV0CnDkNXJH5MLGPIHKkoZ9KNxlF23JmnZZwusJ7raywbVkXm/FmurlpBOabm/l4vtGhUJVl5W7CC/W7wOOFJ8FN2lQxDkOLvg1PQ8WAtUz90jSPTy8uSS10XCl+SqHaeG0rSbQ5hLQdpkpkcswxvFtVesEaMaAHCaX0A7LXoGKjTWSVpYcD6oFKaq848tqkcNuk7mZkeLTS9Oy3Khe7RT8bmW+bMJOo6ND4MdQWeY7gy3SU071EnK9Wsz/tJRbXuLZTFue4uA7Zba8tmXw+GOJjrN3iZ72zwOAZelpjsT5BD9ou1uIIYQJlT4rCR4GS3e9ET9zMVdAoTCU0OGmUiPm+hnN1NEKoFPjooQth4gijx47jNg44p8kgmPpser/THzmxL+G2OQMoyDJQNYvE9gvhucZ068ko0to9b03ehnze7sHGrGtaD2uOj451Cb/giMfW5jXFN7riSuZbDiMPjKBC+JFsjZgUvu2jTdSIbWFM+qn2TRanlpMi8q7fQTwczkRRRjTUeA+21RQ1LbW+6GnMB+YwogVrIiBtL/9gj0UaHU/NVs6goV3KJHYqDn7HNJEoMTE41QQ728CiRUUXKlt6YlTh5cctb3YnLVa44aZU73DbesDbG9PcCClBV0yBonRaZ/EU5PJogk0ABoYANKYAm/KBVL5qJWeRqQUD88FHdFl4FzNun46TPal5ItwiQ9YOLsyQ0mlPmM1E9c7gB6arafiv0BvRpKwUKH+rKmXIgK3sTE/UUA6fBp7sUh1xIbmL1bum725eR/3vgaTtZXV/cmZjHFjU0RqRRiXt0WI4Ycr9hP3I0SyFTVnse7wvbVrnQvRG5M9Fvd7kLMBu332A1nnCifIr4OdLnXwjUbTZiGs8XJgr7zh6p7g8QNU12WpEkgHWOOcHfmcP6MdFrxFn0AKqsWJT7UDq3Mns0BF9/MOnEl6RKZjf92Wy3ZJnUa63EhPK2AfCvmEF08YVjCCJIJN06ZHdluy6piHKeZ6QRpDYUiJHuejgce+hy8yUtQE2LZROYy7Gw5gA7GGDiFhK2NG4ac4EEoLmy93dtlnb1WiRUIivv76NaqdsnUA6znb3qZWU6kSDTRzXXQ8dv7SY3jmtKJI5azR1a8ocj80T7wr9GIC6WUxrhqzn+rzij3ttHlBqjOGiW1XFIrLBkXAWXUYt8ur5fhryt3TVcg7pCsFncvNVp8NF9XSrCbWyoowHelDGG2vkTrN+NYNgA/wGnosQYZ9I9lMsTKygt2c0akdM44mNjPfFs1C8rlvzUMJZcNMFbOy1q253LG+PrgQRCUUQgmIfkOViVMPDTqWXkAKJnGlVIuVpAhx2KsGkmUfMAWCy2J1UkViK9BNi4628kxP51C9i5Qsu5JmNxeZC3tnC8s5z0VN6RlR78Ej3SJHueFeWkBzVWD4Ji+LH5NuzwuVJgdLOL2D2/AMGZ5RFsScoi/OxLtD2qjIN0EzxNWWti80RWDX6B5nVxQwOSH5I6uyK2YgkVgXb3LsKISEfHMYV0QZb8fxfo9BQrtHko1KvJD3Q5VwBK/ZHTORNSHzFWYRiYxBGcuMYA5Uk8LPN8HnjOCn5x2NOeMjRY4BomPJVeALCEtADUwxwWoSPLBcAaG93gMfv5eOuyNQ+NbJ0Q4Riuo/tKbKqH1jyeX/kuJrsGk3mS6cKlP31jRKLr2AVx4NtNAspPkhuK5mdasOn+pbaUQY1bktdYGTWMc+5wDHHAHEFXtgIaP9PDl1DIlbhpgNt0/UHODIPKEdJCeGUTEfAJh0d0PEHIEWC115dIzsn1S1AMTRxN1tTW8zxn52y2loCXFP9yO5+4E1okSXtXXtd+UN50wyDP59V6ScmNjalLTQabZHXFFjWouacmdIqe6HEcJconQxP0rBFKZiY2gMyeKNqhWQIkY1CNtNSSZUKm5Sot3OF5JGTWTayHTOgThtB7ZPVQSMqKnI9nTkirCvBGgFWeyCI/PWxLz9Sy3VDUTtyQjm2V+3fhf34B+DligAyoUPOucYox2qP5UDYxLeUhlr+BO81uBNRxIxyzIJSuIKMtzjoamXtsYujxrWk7IF/5wB7VhhxC30P/mjUUgTzOV59ewPaChAHDZ8kTl4WUTc4F79hD9/Ib1bci7HfOylCMWnu6mYjB8vprZWDFDQ2FNJIpJv1mB2UCIewpP0xA7YDf9IxYJUrqeq0V1hTnqKYKhsU6V23MaVUHJHLLibUi6puo07TQ9vCQlWCmI9TRz0vQOpM+ULfFBGHliBOylKozx1rU1ghElLvFTrkkl6Nn5XZga7kyrvT0yGzga/xYnf7hUFZ1XELR2XkALVVk2ar4guiSCfKhT9Xdfekv/4E6jiDxOK+HGe3Dr7SsBP3mvNWz60gQBU0zaGCvTBTrS290R0tbRlPZ8G4tu3zZg2//76uZ+R7XbWqqHO2zAcVKDqRENG7wfz18C4ygMFps68OpyQ3Kq2mN7aiUGxByc1G5AsHFMHowT/AREaz/nDYHw36w+nClZOZOvfFX2GFPxNREXxJBd8IoQtcmdI6CP+6IpaqrLDeyq5Ak+S3IF0igqFQgqWwqA3qCYkS0pXgTmLClVbdQk86v9IR1nFB3G4TVEQOYgsmzkdeFeQnpISDptJOmsKSyBMmHI5U0RSfjy760ogxJpat2f6qBiLxyxyPkqgriIV/rtmGrhgE/lWfMWdls0aJ4gwoBhS+IHPSmbS6iYoEJttw51ZxNDzy6+EYaX2xUa+/EXP4Dz+jnTNjfWMUR1AoR+yEzRNgLQ/EdxKWmjQUKs8bTlbrdWbaBEzaDoDPEpgK7Ak9+VarIVq2/JNpwQ1sW8B4MmarggMmTykV4DA9XQ/bP1siqF3hnreKrLMe+qEj8dA6v56H6sAvnP0UEEdOEl/etBX03YH9NVKFOi9wMVYio8FqmUXEPJqDCe9sPmuZEZnrSgrlAPmbL/DRGWZMcRcER0ooAkBOdOIyhmBTWYIsxbDHB8AcTKqpHwokX5LdDy9xo15i8poqVYZx22Q8QWEJ2Azw2/UvIHokPZbS/vhy+PLX8Gxid22or7kWud9FB9NP8HevLrbbLGoSVxup+jGDPhsClHzROO88EcCbz2m55PJi5wGdjPjhTpaLYq1I56M/vVC6n9EY2Cs/t3zLhbKitb2tlW4+ozo6SMZUL9dJyjpti2+AFE3qToV0gFAPHtpLpccN05v0b2lPIWYjG7MW4fJyKrhktKZnE7zLKZw94EU0rcVb61ZXsGyfbXF3cXXB2pO4utCq1bCiEQ9IdCNR7XQRfhKP1QxiPMe+HHZebfJd3ZdFg3r+/NCJEehu4T0+7bIUn7w1aiDPPT9u25H+1BCLhUrStP29v//z//r6smQ+Np4q5qvCjH9B+9OvKGmcHOL/fh2ZjPaTqq568PR91B8i/I3CZOXm+g2bM7VMqwE8rRi62H5xxU5lAxTOd7ErByyx46pEKqHzby5GsDDOxrZ0acnUfqacLUobAmoPUkZzEyCKA8q5QHrnWQeFU4XNeBmaixJiscoRN4Ybc2n7Jca2T6NhdLG8Kq/fBSUXd3VBTSnXiJoLXTQpNPedd5IMw6U21z60ke6s2J50p8CJkMUW5brvi4c8K5J1YPsG/UNLiSV5yBbCCs2r8NpE1IcfCvSGXdOlhVlRoHnOia8BQfKQAFfcltCgyJned+iYHatg1CX7jGJ5kU6/RM19yzIZG+KL1OOL93y6Hp9jse0zl231Sz05EMEU03OGCe8EEkH6lgYx5seseqBCQZFTSUVoT1I66M4Kom61ls5hB6TH7IpsjQEGKixOlbQ0UZiOnbfxCqG9RQkCjdfBqL0TRSu9Fq7LrkkLijcM1LVIO3UabQluZF90Mlu9q5tOW2BZemwRn2D7TcjGim991ZtHctFeJuR+yTLCV3WjqmVGVwkFzTtEWIoLfi8fWwHsbBqwMH2coIpubnCKqUqZJXYOHc/a8ML4rcBrFN4LMZ52ApbNJrfvQNjrisiLQAZ0bCJkwnuyGEh8l2gQiCz7FlFCD6/MZsWmR3f6utxeQ9pjL/mQVs8gC6eLwwvdcdvdY83QHg5MlWgQtdm+WugiKtpQh3j9WByJv5BhDlhT7zP+Uc7//6jDP8jUxj687NFihLoaIgynJlPpOwcob/sZ9fGCAnJWAqBmdpES1pbL0isDMbPFIWdeTRS+/fSLFaUNh/qSl//ER7jMGOBUuI10coKjRjflwk2ws9/9XXcT9n3Q8OfZXhZvOI+ddLuLJZb66VZiCFLxTwLz/M2r+evpbVwTZXo31SbLoHaNIWQnJtDOFkyIFMKvulmEwrPwfhK8vBvj6TFGgXot8iZwg2TTtQSFJV2qECiuTRXWtLEXo6ln0zaofux8ZFkPnpK0FvR+ijM2jb1yw84cGhnsnX0Ny6K5FSrNnftYyuLhVt0KVXHwRiW4CJjOC7SyAQ/Him9Yx4K2eHcMbqiaSU+V82gvrWjW4nzx3Uk1OjBDXmKAeOINcJGEbDlYjearp17PciPn56s9h/nMk9Y7maMB2y3b9wTzRXsXYfQ9m0qCOoThbc/Nnc0302CMiPkjxg6aixj1x4fChanG791UPZsqt7HOd/xKfINT6lj1ABURBcC4nnJI6Hx2/xC9YedfJ8ho8Kwgo9BF6yaNOtvoVcsCzmO/vdI/MhnGUdtg7NktgdKw69K+9omtB+5FU7dKusIgYlJQSNLDr0glqcS1dkMpTyuFZHLEJnlBkY51XIpDd+Z8jCek2KaU0diNIzfRqp992J/82ehCg67bSPitMGQwEnvvN1u7EZttxVTsg1TXolp2pJuW61RDamubli1QidSAjW65N5MsBUTMUk8z0VcX2y05G4ygIkzzcu91PVO28Up5VM32K7Nq9OJiTTCTwwGeYfyw3ojwiYoyRj2kaopxFmh5YFcpc9Hf/+f/UNHarGhgtIm6q2Mv8yNe5cRBB3iXsvLul7m60EwreJI/OM9Lo/f1eL38HfhrQ3s+n9N+5pXOgcXqAlYb4aottpS2rXwKr23v42Jma7HWkz2SitQV7vsLWa918bBtiTwBTF/SWntuSadan12gZWb1xlF/9gg1z4kZcxPQN2jC+UKZwyllPnoCR/83IooogunuFNWp9rdL7wETnhYR7PSeR03/Efu8bTvRa7VuR9NJCvpqtIiJJIynCKsntZCUKHEIPOHOLQG+zeaSXbjIkXB3927ybu5sy0V9X3gZtFeMhiKL+IImNspb2TXI/qBZmcrKSHbi+uXrUiYvu03KFuv7nb74pyboTtvWKmGuTNHlrnGktVyla7wGiy67qjgFgUO3ddFrYLH5ukJTdVew7T+nXizHQ4ShOhvTx/+rG6y+fKS9UwN/xrYRV3Vy/L6Xg/YFUDtrh9rTsGJvSX+f1FiaoCsimQn6Zau573wB3WgMxMmRjKMa64H0XM+wKRFy+UaCCKr7du6rd2WFMxtmTzgIAjpXXizIvxq0z/fUSmzpehxEsDSu+Ne6wIFfXILLSqSEbMkSkOrWyT4kFaxCmRYzffdJ/tg43PpWwoTu/YTy+G8sZ+IAbMqdfFdYz3QGSLxbpWWHjCXmh7GTM/CCdScALRbHYxH4iO/vyZnG7s677Coyvf7pszQmWTKy6jRZOJOqwF/EE/SltspQSzm0ZAoglamkENXkveQMX8oD1gSeOuXL6/GnnXFnJ/hrWGjMtWcS0E2B6HZbqoE/Ld8D/N/aP3oFLA5BoKm4qd+SW22VJfsD1aKNXbhCJaeDssqxixajTb8zpUVOXHV44V1TrcGbrUNvvWyCxqYfbW7nhcXKi9mRgwBNf85N8oRK2LJEBQACClrhrEEFEU2KJlqpOSZJ3YKuJQxAW74yPGNQokk+7KTlIaBHJkeArhqnhG6Ct2eUehNetWu/JFnk4j67SFkTSfypxUk6Gn6WA/QJTlLcScx5ubAwDXEPvtTEMvtxFg4lPJBnpjhw5AZderfFtnTndR9zKn42LaErP0lC3Zei6AjXkcEKFDIlgoIHyIfPM2AjjDnrZtoAFowDXRxSSZ4EQky3TbI/z4RzSohhFo7npC+uX5f6shYac+Ec06jPxkqqzyD5toEmSZirYUWtaBeqi20Jbc14QQmfBrhOaS/2pLRN6Yplm57ZxjNlCI3WPHMNVtHr2kLDURCi1zDUiSP3Df7ONqUw5pTMl+61H/Hrt7wda65XOusG9wwyhghFiI1TPaip6+GaTFVuDmD0XximmcB6+XEY9NeIhFa5KV1wgtAqSylYyNCDvaED3DXtQEXeQQGy0/EQ0NFIPRIbReIV73wkIK+lDahtxTJVdyN+egaMuaacFX9TFlj767r3arCW244lst6QyEp/3LBQQfv7VxU6Ve8AI3IUbIhwYCS9ePka/SLcoHrZd65kZle+ruRKO+TKLsaOrcjyJt0if+NgcXWQqCWrIj1YpoCvilM0CwGhK+SHhMJcEuSGMKsMM/lBsk5BDcpZwSZfTk1OXlL4dd00WRvy10hjPIsmlD7I7m+qRjcXZH3uJaHeDU7B3RqNXWsy1TkW1jRP5uCrlhyAfCEdtlyLwzH7FpGsdhr6hq9eDew3baa+L2DHi94qqtaFonFX6F+bovDDGS42/jeR2DEXpjtiQ/MuuP7vk3A+XXoxdiqA02kDP3I7kH1gT0At8gl3tvo3GxiYWD9GGKnHqZ6VljMez9xdswskN2Jy21DPQaIQSk7lgvkzM3c0O/eTOjA6VCYvH0TnHvSX8gCUCnkVzDFd4aBrVLCur+4md6/u7sTgKzGaftUlNMEfU0AgkqC+pWqEqiSh8s5VoMmhtkAlyZoSSJZV0xrcJmE6xi9MDwjdpg0BUgCvwNK5dukUiIejUiypCR0dNaGjp7N5ZvFsnpmHSfYwqrB1GL4RMX3pSnAb+aBi/SmIX8UWV0wIUVRHY29sM78LcFOFxTa9E3Nilom3JFKZK82cMZ0MmJKQSQliB1aHryocEUV4kNuPlYyO+llZ79OBgannp76HAm3olYzo/LwjeOUhllkDbnHMU+iNwyvYECeZ3eh0cWyAMoWO3U5hk3awPWurYkHrHl3CCd1NdePWaUd1tbIh4GNv2FAclsrI8mZL71zZQa2Upc1EcCYBZwg97Mi4qJ+lXGIQjhoYly5GiOn6ZQESCslVqtBQdMJfor5Aa8/6jpiW91bEAdpafo21pPOL51aey39v6cw6S0X0qygpUEO3RG+fyIX1O/KXKfTFN6ayYQD7l2dFf3ndzqlafGYpTqrvH+vySCVh2jN8/UGHeBmfoaIKoutGejcVNJGSXVsFDjH3k4hqIiiuYM3kgOMFyO9GoqwuhUhR8RVHoXSiyGCtjOYdlLa54EjGPj8af+kTitySdnYhTSQG/e5x0gpR1gPf9G2ZCiJ4cGGH5txjjuhYcZ3LRtoXv33eMKd3CHNBASp+tVwY2tZEirJzh50lPp+6uP0FsfRTxVZjyv5MlTENVUYrt9PeNiVuX1qTxGnm1Kwfcs36yfSrTvBZ5JQ+umJca8ZGRHGxzLVuGVHMe1qo+nm/aUfHHqu7mGELjgXDmvQ4sAoKw3pCoIrXVnbWXhEQrEWJof+qpFmCtkweQWn9TBFQ31KyAHepat+wpdt439lw+oB1Qm3v/AYNA2s7N0xPuMt0Cxk4jowTW3CwfYa1EzFIgNqTuRcOmuxB2hChCtJQTElH5fiQ+JfJ5F55jvQsMoBKmPU17uGPMN0OT1XlDeLYXI3VKQ7LG1FgwfZrfeMNL7/T5UxWE8CH5V+50o41Krs/6agyELQyY/D6/b/9d9rub2B3aGxFynk7dIh937KfW8Vk1dYtdLUobXtmw5qzZzphj8sTXf/817/+h5+weVXLw6T/t6pD0YMsA1JSH8FArZ7gxG/VsdjHuS4T2Byr4h3b09nCvZQ5KGI755T71oKpFiSpV7dN5CPGZjwaleARjrnfZAsqgDsTVjVSdmNSoawQT7t22nMs4rZBsCnm0lZqjWeM5jgMeMAN7OFeO4lz52wi2jbRdEVrdwc0RXKtVmFRj/M5eq2avaK4p2Qrl43PwsJql8XyTV5P57NZfDFYf3XdCBme4tx6f11bRZCg/y+UdB32eypBrDWUMVw+Jd/FHjdugMvS5Mxs/RQ5glGmTh/bD71drvsSsPYECc/PJtMwFZTzuLBMEO0t8cZ4bZTphZfVX3Yzc2y2o9ARZGq+kpWM57in8pU2uVNGcao9Nedpul45/uE6jC9Jc7btbwNlgwtTnWdT4B4DNl6Z62sHXfy3P4sVford2NFEJM+Uw0RRePUC79veY/3VSsZ2xIT1+s9xlHvy2+gBhlyEwG3XJzHio7OV9n4atWuTIIgc8zrNyIhwaztErDqATmohO7l0TKOO0GdnLs6jr+Sm4IpwfYX1dGRNYQVSTakZOkJD4Pjxc0NX05Vczy3V0RQgJMe5qsmBHD1j2SapxMv/BFvzEp9B9/8fUEsDBBQAAAAIAAAAIVwXAir1qh0AAMZ3AAAKAAAAc3RlcDQuaHRtbN09TZMbx3V3/YoOLiFLAJZLSZRMLlFF8UNiQpMsLi2VfWE1gAbQ4mAGnh4sFjwxOdlVOdlJnFQ5tyT35BBXfJTv+RH7S/K+eqZ7PgDscmW7fCExM/3x+vX7fq97T/5qmk2K7cqoRbFMRh+d4H8q0en8fs+kPXxh9BT+W5pCq8lC584U93vrYjb4oudfp3pp7vfOrNmssrzoqUmWFiaFZhs7LRb3p+bMTsyAHvrKprawOhm4iU7M/WMcpLBFYkYP9ff/ffF3v1dPl3pu07m68S1AoZ7p8U31/e/UaWFW6tOTI2770Uli07cqN8n9niu2iXELY2DmRW5m93t6tRpOnOtFzSxA5RtMdaHv9vH7kSxvnE23sliTq0minbvfK7IVtFHqZGrP/LtxrtNpbxewJ26l05GHlx5OjmCE+ki5mWT5FDr3lJ1Gj7TG+71ee795nm3CT+N1UWRp+XWRuYJHnCSZM4OVnhsYErYYgKemAD1+UTrVydZZd3LE72k0hNePtbJJwkOZFFZpoOPCTN7ieotF1V/x14v3/ynLFbSafPTRRxHkGYxB4yUEZJ7hMoEMgGhGrxcwEgKr8KOyhUM6Ss2ksLC6IlNP7Hd2qP7WAF7xp8pWJh2qU31mprDDbp1Aj9wstU2BxghA/3aWJQDMPXjGPjHsszxb8nh2Rl+CSaeZcSrNCoV7c2byoSA9WlSit9m6YDrRzk4NLRB/9JTOrR4keoz098DPCAxUAAoddYm31n8C5gDcJL6FE3ik1UTnU0bjKs+KbJIBM8VjQp/F7dFru4R3sBe3y7cES0mAREe0Y09yYGFAP9A+TJalU0+4Nl2tC16RXq6Am+ZvFu9gcny9zKawfVMzsUudyGKnBkC3YzMdb4F/7NIMFiZZIbnS1CUgKw/EzHqayHVhBoBteFzY6dQg26zaOygalHoFc5wS5A4AMco4+AADOrVOC5so2CPc2ZnNl7TLM1ywwilVBmter1bJVjELAjnhoK6ABbthCAPuFM65sek027iBmdoiywNebOJ4ggzT83hkLqR34+xcuBSEJWx2b6SeppNkDfSj1SpzICTPiBjx4wHYk5aMCo/AxyCHc7VZCNFLG7XRgKQpLhRkmXzOlIPFA/EW1ZKBCJjyPKWCsLeJa9JhtsJWOqTDE7deLnW+HT1DyJEHc4IBeAUwASPL51asAe4Raw9hOy2OHJDhxL/r+cYbYrVVoidmQYx+v/dC4GnSXds0X9osyeYW9BEQiy2Cycbllzf45ZAp7wJOkR/6agKyZ52bvjLFZLiHA5igRw9NAsBlQLIOhA+JHqTDcWLdgrCXrpdj2NFsBhw4NSDMYJcLVYEJJAxMOkHCD7dRNm7fNpbiRFoEu/hgeqbTiZk2N27xyehVk29A6nyyc6kP1MPTb9TGFgugd5BU5g3K+2S9JOnNQsj1FWgMUBMbkkzEs0P1tKBVAu5djZUjZo12emYTw4K4h2KRgSQIgt2uoB9ge680+beeTMwKTBowLM76hTkvjuBHc1cDYU5zkjxjIVufA62m3uh5prAhrDgBXjOl4BWtXuszSYzOezVdH+l2z/qvzBJUVqnaYxGFe7aGScnwQxngd55gifauS2OU/T3AvICYjlBDQAcWD+pMJ2sA9HZv9DjRKyf0oo6AU1yRLU+OuF29PZiIz8xcT7Yov4Aw1iQZSR0ZV3U6YgAO4TNWN7AArzgKnaPY4hndUJGZo/QcbAkwRUB35F52um4aE2ny0qMS9zegrhIptO+RJDlISn2jEwt2KyIG3hRr14rxs7JVA+/Q5sw6kYsvq4cuvK+dyd9MDaINDN6fwJOSp64eMrmZvhlv32D3EmrY6/FW4av9GxawkCbtg3Koyc3eln26RKejouC/OX3xPFK3zL8RluySHZWSqb9zWdoHt4FkJ0x5hC+qfQm5sRzEnPMgsWn9+DwCp4MB2zb4iTUJsQRSIpqbOdCGUzdwQTdPUOSAYaMJhhk2HZRt0JLewECf9JRbgQohAwNa6cSZNrus7FgaZ370g7mnPkSo/mDGXG9RSc3al0Q2t46FTsRXIH8eWT1PQbzZCX123UIplO2nJiXLyRZbNYfJUMSHNmz1+YcT8eEkB8v4sNM1CflOJD0HaYoydFDag2zP1pCVSrM33rT8wTDWmOlgtDV6XhPu2iyWiiLZxXDkeTqz0mh7gPuAilSstslCp+DJktGtE9Aa6huTo7hVt0uyd2DcJNs63ZducETwrTIxNBGy+dyrs9pyJZhwugArCmzwBXqoCBeIrwS8YBoM7F4Me+Brb8GzWZBkBcxEfcP2FdYCGMAErc89tQ6MVzMdfWPNBvz6eRu66+Zp8KXm+ObrdIB+0w5Xyzt5+1wubINqTVVUKG8HSMc9cZ4VdaxIr4WwsWugdT004SuOc7DWDpzbkOB2bbCPxuQWDW9eAsUv3tUDOxXCX63TZnwn3K0JmvNJZ/+H9Lm1J7Rhl6GVrR7xV+EA18FfwXKBHea5ca5U0fLk0eRxPNalp93Bqx7BJE/ld4T6Ud2Z38VSwAjrBjOV+Dk1hVqv6iGwQ6SJ95hNOpgZ0Meo+7vArDvgLQOdsVwpO50cURCqHjTEdwMznUuv6lFGbsqS0pWvoSCMaX2NwiRoyKLmay9iqmjAxfv/DaOM9UCos+/AWJQwKD4MOH4maGERm+UyeZZb8HqJt2j9BXq+5eyPcj3HcMZSv42hUGhr5ygDUzBQso3Jh+pRtobtBK1hQXRAJ5jdFMMwuioRP/iFkcUm8JssfwuyYYLCV1P4lkh2dLI45lhR+Z3g64HewUC5xCcl5Auy/ngEu+s7zJJsI8JcphEVRGsrfV6X+BWtkGfUy0Rvh+qnGOoCmZZXE4Gms0aNzQwjPCBFUxRtnnZJAwljhdHfsS7JhYcaAH0ixbyUkb//Hei6igW2pqiHu2tshva1mUqMj6RY9GafmAzHKtAmLClmnRQDeSPdWrgcG6ixnfu4Iew7dF1hPIl2bQxU+mugUnFwq29IGIhFEcAKQ4vrBBywxmqbE0aTza0rTI5qpz4ZyNoEo9e+xSWHNucUOWwZ2H+55IDrFKOBKOqaQ1bf2geV6BKO9fM1aEGwaJlgSmL2rycLW8Uqfagp6ijvkCF93CmE2jdbYaJG+Me/A0NwylKxrQchvBd1kFceHAy7cSyffg1miZ7PCcVH9MKzzDqJp4UdxDbrZFQT16MnMILrshCRi2WzFEESsmVsIbXHmFlJGpQcIEwzMfZlOV44/t8/Pjl6cvH378Xgu/jFr46VLpTLliJT7hJEqJTmeQaNwEfIc7K+1SxZg/hwEwOWgdpQ3BFnUzjbUL22aCdRdycrICs5MbMCg77oDhITmXyC0htsg+FOI61aSmSlBStkQ63DNAsUGutaitN26XO/6y8wPxSnjWp6q9QH8DvVFf/ocWltcMC26NVhqZoqmlsXegAPXrzWYSsH47gwW+RFjnHPGKa2GYKlVfOAhiqWumHUBBN9zS2uOP7SFDl4RzvG/zG3UB8r4ZcrzkQsswdhiiW4Ewa76kwUjd8x1SvWE5g1pqZXnAd/7ZjlZZ7Bvjvp0pijOz2MJAt06km25lEB5xhJ8dJPL6lbVUSWJWiEs/sACB6QSS4tS6FUJcgG7LjHkkjwBbvCjVQYxvQ7NuBvvc4gYTQbSqpVbRZyWcOhqVEZEx09wESLIbLYE4/0UzFNhZkok5T+kzM6nyx6CvlekIMxDFdLEj2nBLheo/+P4esCumazWWvYTiT2awrRE6x/7cr8T07R5T6Zfmi62oKDaihhEQaQx5lyiAOw0LDvXXXx23/qKzZ0zXJVgBmDVqlPlcMonKwqjcnHboJy23ljHgQzCmwkP5yhogJeKFKcfGrBYBQ/BRLL1g2WkvW+lM8EtbpB6oPM9Zux6xG1CxyzhqdRnz9FrdE+93P4JPPmdr5on7hsFE36+wZLdllaDkzeImSb4Bt5zqBn3UJ4kndlQP521aWVL0Ew5G+jRnU7gWx4iQFSQCbkD9ZA4UJfee8E1wlkvG1kFHK96Y1eVU5Me1YADCdY75elSYFFGtlZdxYB/H9T7cnr0E/ifGFsofjhsLJIUV+2YegnWh0X739byhxosgTTA+yPi/f/RqOhaUIuEr/ujb568rILsrJNy3hdfdjOw8oZsveIS8caEE/hgvY+su2YAE8w94SRc9yfrvYVMiK0rew5Zt4LcIRTpZ0yerLwgqTqwthC4jHA1aCW84IzzFrliPI76P58esuP5kUC1ufAb1xGv0RkBkPlYvyBIAGwt2rMrJSNvwNScxEhYMbhEAHcZfLuMh8FNhAzEsIbKZT5/nX7HAE/SEOuWosZ4wV/AryBME0pWCVYb3zBBLZyEyJFNN8wtU1hY9li1+CpW8M7vRH0T7vToWVAfJRm+bI7gXc8/IKGmrztanEbW9Au1Zr5rWjbAlnrl1uyq+4yWQm6vDqZCOHOiE4iQ2yIYQoLTGpoYt9xs7CFGaoXu9EFGkhw7Fp3i+dtYhUnDvAmsHfyHwIINhf+9wHk6TH1bW6xEClkP9HjWKoW8tF+YukmeB5SvNiRes4ziI3TBnGgE99l2XIQ2jONeOfPMrQu2JZo6teaumuoPIobRRpNQkknE52eaRd+QkuC39bgDQYkTQY6IfVRw/BF3K2lIypMF/STZ7a22U6M+J3WjhuCeia0lhpT1eonCa2uMKtaaJneY9HgjtAqTcptmIr8C3UDvHb11mxvhgbIL35VxZ4bICTAamLlh88+XAtW4kzn4cKUdQpbGiTC3uj4D79pyb3RUDbduwibxmsAsr7xMS9g9HEQMefS0QZCY/dpk/liJZpeN0k2zHfhkjZgfFcLu3ELZyYzV6tpGP3NIkM4RO634QhtxH8IpBjby9qRRfa2mTYsf6zLm2LEV6MFVID9AGDEYA7V6cKChfwxN6SCVYp5Y3i4Mm9gO3lGkJ+ekxmgruU0+JleriI+AiVVBGzEj8xFXWWcnYN7QRCOT6GKcnh64tHtch7T2ev6Pvd9ZRkHyjXbG4Q4sG8qORvZS9EIQ/WQySJn1GJGwbM/wjMMlrdPXFH4fQBekivlYO1VnP95lumpL7IubV8qr25kv+pTYUVV6Qxw9DLXqZNyl1brHk17ykK164e27CN5A12q4uKf/4cyEe2kFUoQ2pwBZ8C6hUjp6VHzpqMor3d7gA3ZTHXDm1yvako1J3O9rLvrqaVNyeha6nM2vkCe04/KHmvHax2uuAhQQhxOejMt4W6tTRDTbMrEJv5mWb5pycc2fNYYfeGrXX6szBoXLK2MmZYlN/gQ2mL0tbaFQmKKvzVt3s96I/gH1UxnDeBx99fbvdHt7q8w9q6RbwV24fGtHZPcgll2fP8Mvn8Wfd9tMMbUINErRpo9k7gQ1+Wr4yYhNIVAV6WWZ2+puumWLasumTLRvFAOMAQFE6B9wZYA93caCEOpBvWNSHaDyu1Foh1fxG1r0cvaKN5slMqhlo5oKQeZxdLORAO6mV/n/rVCaF/QRwde6Pegskl9dsqbFzAu6YZskyaAUlayYRQNDXviL+eLAttzWDLT0qRrDyf9rtU84rvBZJHZScXDkWFyV2n0bznLQoV3fQU2EHpZ7abdUD3AesJSF94N9OMspg940gX2jKOFimLGYFPQQJEMwN1A977B66RkA5arGVftvCXKehQD3M1pUR0HguLJYU8WqEI/eGzLXrjZLkZiyhtOefAbTTNTDjvoCN4+BmFoCU5v3c1+Fc3yUQkK2bLzp9FO9ZlncWsxxOVz0ViZJRkqkMFIjfAlfNtauxMi47JYeCSkHliYfrnFAgysuUTmsBYVgfLtMXqdrbbd8AiXXGlTuK8nWNqOAN8dCPaxh+fwVQxEx444IL/kYSIlBIBGKdmZ03AWo+9ABVqwDb8pLEY7cM/vvw/0hT0pB6rf+hIO7DTe8v8wHp53KIUObSuvkCR0DYGh5GNOcAMu3mlPJDVQ+AQrImYYQyr6CstuCELwRrJkTewHaI0jENTASzWu+Y1Bq+Wra2qq/rivCIhC5R9UBTS3YrNzCYzsC+4SvMFiGBCFB5QBtSnFWqgehuQ4fa0EpdERcybVmZvqMdIKsKAp4Lys8aIHj15fEkd9JZHai8YvX+6oqYmHYdV56TF2+SKEkgWVYXS6Rs7Md0RevkW1U24Y8m4QSYlJXMYrE+vYA0leN7xsMrFDvw8oPMdzp5RbVze+/y91+vjHN71M4Wg2aF4H7jiXTWBsoU+sAFhboH406LiPHvBADUatCbYKPA657QkClLrSBw8ZzlIqgfA1DjOMGG322l7WBEya8hLUx/7dJcDjLPs++Dz+IoSRxSA1VHjQUM7SFPX2eBCZ2k4toTNIqIN6xeO3qS/TgjftBZZtPkpQFEPCvr0eBsTExok2EOyiDmAMjw3s7TRG544ALI3S4zSDL9M4NFJMfkC5j+2wfosbnsl2orXnMgxNbER1UeqlHIONEvBd7XK9VLl1pexLM0sol4Y3OdRcuEpZIdMMFSfefr62Oepg7AP2jk42YLrQtOnu8p1qNYCS1+XDDpSMTu081UlsR9KrulCwk4XiL3Xh0KhzApg3XojE9SB4GolXOcfTKnz0i3PmhDqjHWi8JdY5AXZMMhNzAtAY5g2DRCQn2sBvGzas3ik6cx6ozlShL7Xqjb5CkChMC0zdzHV2jcAJWb3hJXH/HW0jqJAd4d2l3Fch4u1An9sK+c/skihSbDJRzQrbDNUzg0aNxnICqUPwJ4OxLAIsETtRCQ4w9CcNf0o9mwfSt4MlhpbbjqJHRRc4bkxBz7D4l2fxUPIkSNGjIjs5srWJ9PmVJvrJatU5Uee1DYBNGqg1rI3ZU+8V1BAGg+KbHRH1/UkaGGm+KGp5GWCvMPJaBTwjNxb3MjCF+JnOIUT3LeBBc2NikUvScniARSGRiF7s+supF0YQR221iHCajaQJ6RicD6y9OJ4rrovFs8Vfo+h8WOQJsvnFP/wLCUe+jAF7vMukxoYOruJG3ouTCN52VES/fiIYYD3GculJB3Bc34Mtwrk/rs18YLykIyCCm/Kh0RCkhENDIdT2aoGPhsOKsF+T694PztE7r+KBPW88PP2GdSFqAHiBlyQsNXllLT439CSXjc4GUz2wu7QbTnz1IT64d7/boL6yS46XWlwJtlNw2fr8AEpLn+O/WKILnjS53gilBDK8MBTnEiYrEzYrX4H4eIpFbXYOC7yqh9lZgyXOplpE7ubCkPA7xN9cwNO7DB53nDthj7P0+dFIgfYJl++BnYF1Nntdz/oxoMOqOMtC3x2nIaRqDNEucUXXVZXW/CoGG550jsowcmbRskIkMrn6yhV5hqfCKANXZCsOOSY6p3dsZKaxrSZVotZJFyQdqqQgdxELRlDckjwG9eyCzuRByHAEGP7aAuPOdSqYxlZDxdXVVTE7UOkSoXLFXelvc56uukoFOYrq2+k9H2+n6pa4TK+lCqUqHotKULp7YdsCHTm/0AP6FBhi4qsuSCDB0wG99HLVG73E4BTd9VOsp4fMBawARsEL+NeiRU6ovkRVrUiE6DYDnZiiqJzML2Gei/e/phIh+B8dErLr4RFPL+CesTtOhEFyeRFE9OimMTb12UuS8R1epAAUXA7h8DTTlj36gBzQrndRlRyLsZ17jf4J6tqq5qixiANwC/alW+reiP/vvHfB5naKRqX82FG/iIPRfwdMbtOZydOsN5IfnX6LPePp5cfhm/+KaiKjCnBOv/roAnlUZWVzEFlA1kuNpe30Z5zAtc2jMynfsE83NttMpHDOh268iDalSuqXglpK11A2gsnyVubDlqShLn7572Sf/fI/dm5/no3XeATkGGMemfrRj6CrHH2hG1Y6SmXXKCBOl6gmuKPIxoNR2l6BhqSPZdql8ATxdfH+X6XIWkwZW9wTs3xmNizgmuIOmWOn3w8tfUU4CVaat+n274KWXESqiqvmRsD9NhMZ9JH/YSULtLnSTJVuLNnWO0Fc0OHX8sicovOq8ehNeGtmxg5LW5TvFYxtMqEji9sr8kOMbkz6SuvrMbv9Qq4rafbHNK47jjr9udjXVwSvy8Tutq1loiub13tPM5CH6vM2nmLokoYwhlDR8Z7YgR9hT/ggoPZm3ODS5nJ5bi28c6B26knS8aPXdCw0VFnSe1DIcdlGMror4IvlE9x5wDef1QK/YQyq2bZXl20jyktQFQDvvTQs4eg4JLZLlvm1XYcs81g+TJbx8VvJh4PuwctprkGwhQu6ukRrA+6yogy7HCLOPMQeOx8m0loh3y3L8Fag3d53270UQqz+gqn4FHYpNZrHrCsqiY9Z127C2MfSkrraz9DtJ0o8y1GJ05UZLmSs4CIhvrIW0OrUW7zNVmKh1HlYMtXlmZXWfC2sytg7jFEnzaO218mw1aKuzq4tIPYF6+HG9X9wXubFXA8nX35RH8zmh3GwJ54P4l8+Dr03fNUSk2omESmj3Z1JlMxfJR0la/mXnMKL7nNC9LQVYQvZPdZ5Yv31qy7GYO1bVHtdVW/SDOLZRzditJwqobYdJc3eP9RFBzjRl/Zq5p1H5fe5dkxJOwRsMM1hkpbo/DBByxWGTKLXKGDLNV1dvoaQ3ZAQQk7n06s0380fXLjSQq5Htrah+gPEZ6sT1XbF4Iv6bcd8AbLci7yS83wcymQhTxckYrAzQHsZysJwEuuJp4+aN5nvl+RCnJEgb72XbK81JpdofLA5huP8UR0ghvx6jCpGwmHMntduFblWk6pa09V5vg7gn8ig4qVcD9dfcknXaU5hYnA7qIqnqseu66ioScPKOtBC85QYd44l0grrNh/KXfRye3AQhac/qIDVe/dYwjzHlAxebMXX36vncn9VIChOjvgKPREkEbBBDXOviuKVb6QSjrdxajXdNRrdRIhSsyxS9rez8l/lYO5/QnXOzVMKXP/cvC0b74R7+fyrLtursGh5vX765MlOc6u0St+Z5tSYdW5MXJ51cG0nHeK2d84/7Y3uqD/8Rn2qum8M+Pz8s97oc2z12Y5WX5zfCbJVX2DzOzuaH986/xyPbR3fwpbwc1fb2+dfQMvb2PILtedMvU8NVVXoDbxNV7aBtmM8fwX/KPjYBccnt8LDXvC0q/EdaDy6U2/TCXAFHV7tg8Wy0L78kwNAJip431xS8LGxtOTdJgD72c++7QL5nQUb4mdPO+8vSfH6oRFeQrR7QS3hHDk4gleAVFeN1u+3Yw6rnSFAZsQSv1CWVqeXkS8s/3WQljMGPslDYRLOw7uFXlHWl9JydDvHPYxrY3kFXm/MQWHJ2WIGmM41UAoxCHn7mxfktnZUDGwfHSSbGLhLSKfgIEfjD0rtq4DH9cspFwx+o9I+sASee/nC5F1V8I+1n6R5DCY85MKpBdgFaMEZ0+C8iqQYykNoSzBV4mM2h1See6BX2ua7S89fpMHZt8uftulLFHRsio2RNCcdagocb7ytRsak67CpjKJegR7wvxAGF5yThckYiCxMmybs/NqyHAld48DojIdhv7k8dhGVpoYtowrVdL00uZ20kCHzIRrKS51uS6qUH83bFJvL83imFZa7sTKp87cTYO2Mz0NRYMI767KKL+Oduhp6KjAOwJBv3Iqkljskm6vGze9ecf2MBhFdey3TPWYUuZit7Nn3QzFONguTh/ZW/TxFbuRmJamm8icUPBpeVlR/NfTyeg9ALZUaXQ6tf/bnV2Vpf7oTrB1qWOAiA77Sw7uuMm+7vkFGmdeLzMMDh63ulnTscreuJYleXiR/QC7do6PtLvvGFfai2ZsKnk0C/9fceOD4VbeCB4rXidemlaxN+A+tyDD85yxqOyXfqr9n0PKRD8+djEOYeLDyOuYWGuF2lbG2r2o8RCj3pb8c2fpHI1v/gEPbqqbb+ECjfAiugoreRNUCfgX+Xtz6VS9iw1GJHwhTKmHyDl9VYNDwiP3i5C85RVCUL5sHG1uGmGVUCN9EHFXvdsXR2QDm2+jF+CAhP5YzHXQIa22TYmBFblCxqKFCsKmZabyVuP2ap4OPAXik74pnNZd1KYb3u3NdRTOtEqxaRtf9T2Wkitv2Cb2+7v+Hr/av/QGBFrFT2PJPWK7q192fYLxlVSiXT/jP2X7H+0NvcTT5e7VH/Fd7/x9QSwMEFAAAAAgAAAAhXGwZMCYPiwEA8SYFAAgAAABzdGVwNC5qc7S9S3Mcd5YvttenSHLUt7LEQqJASuoRIIgBgmCLd0iQA0CtmYvGIBJVCSCbVZWlyiw8xGZE26u527EjHI6wd7b39mYivJy794fQJ/H5nXP+z8wCSKmvIkRUZv7fj/N+rK8nh00xT75M5vlFkSV7V8XitrksZxdJfVldz5LLYlEkZZ0sinycnC+qadJcFkmTn02Kmn/ms3xyW1OJOr8qxoOkWrhii2JULcbUWK9O0Nj5Ip8WdfbZ+npyRJ/RZTJe5Nf1VlI2ybiiJmdVk4zyyWg5yRs0UC8nTZ0lh0XTUDt1clElZ/noXdJU3MN8srwoZ4Pk+rIcXSaL5UwGVVM/STGjTzRAefWi/HOZXJezcXWdfdZb1lSoWZSjprf1WXq+nI2aspolaT95/1mSjKpZ3SSfJ9tJWo77yfZ3NLbRclrMmuyiaPYmBX4+u305xuctW2F359UrqvMew55Xs3Ex3kx6f7fx+99/+fVXvQGmhufdv//66+EePVc0rAXe/P3ON1/uDHvJB9fUm/3Tty93j3442Gs3+GQ8PjvzGjw/zx9/mfsN5sXZ348eBw3+4eDND28PqbHj3s5kkoyKyaSmKr0D0zIenpdj3gDtDq/2bkaTpX7+YVbTfmLr8bRfFOM6+WmZT8rmlqpclcV17wQ9rn/xBe/vqJpUy0VS0ZnCRzoqY94KbODFolrSY07nqsa/tBsljWtOO7Kk7reSAidRppRwbTpiuZygpDoPj1eWfLFuZyoNbWMvaeMOM668vb2d9LTxXvKXvwTv3Xh6bsHe0VFBK3QYJ9wSfkh5t2bJ06S3cE+bXqF4LVGUnlGox7Py+ipnf6Du56a7QYLRzGmDabLc94V8pnaHGHzqnjeS//SfOofWjwo+DgvGw0NxGk/iV3kSVrFHIW76y7Ccd0zikl+hpB6ZF9Ui5QlmE7qpzWXyHU2OPqcPDjNTZJJfoIVWjXLGQ6nToGi/713GQ744y9mknL3DvZnRbhMUSfuD5IpqjHNc+EPAOPq2xKKX03m1aLisPM/OK/O7bggcmQeFSuaRDmpNbW0maxu6c5vJcMBXzBSZlDU3fHwykCO8mVBZnEA6DQQAe/Kwe1nVBTV0nk/qYsDbwdfAvKFTT31fVpNxbV/NJ1XzuhpTi71pkc+opbq8IJhMz+Pzc0CFZUNrUFCNZrEs7POP5bi55GHoi12+rVSsh60EeFlOz4rFLiCF7Q1AfofuLJW6pnEAEDTlFFBHJ1rVDc+ezpZdr3k1oQ4u3JDz5rKWxTirbgr9SYjgvDFV+GFn9I766Q0YuMgSD/XTkd8pv/nHZbGkHt8S7inrIqN1qiZXBbZbgL7dyfzKH0s5K5synzz3O58vqqYi4MUvaXTvPwzMEZSHaQHUoV9GNGfsLTfog9y3tDN8CMflAjVl+Sf5LUZuEcqIVrMpFKektPazq7zu0bBvhmY8Nxvml/Z/CZBq3uWNO5f5mDbTHbpZYec5pSrjnQaHtPiGj7Atd1FUU7s6o8Xy7MxfIODnmexQXSxK2S1/mq9e7v9h74AmuvHVcLhFR5aw+7RW7Ewr8K4o5jVTE6Aqrpmg4I8CjRmunzcE5fFS4Pskp3Z5xLaX71/uH50evXy9h3v9ZCDPB3s7z9Hx0HZMvRB2aIoZt4YVoB85IYX8tuaPRDf8XBF5AlyEHtA9UzoDxlMJoFhd1DhthIRyO+4RoYeGqJjZqACtcklY66ygbgiwEblU08GnQum07n/G4zhyE8ySl1RR7jp3nW6gU4yPr5Ap17f4cV6Vs6aF5Oio8otpOR5TrXR0Q2fvdotGS8izaP5IHfSz5EfMlorJMLglmhX1N654V2k7aZw5nXGapKzNNL+luYwqWvhc1oRouBwdN/nsYlL09fTQuSsmY5ydC6aVHv4XzIYWVF48RC/VjEqfYSKpHYRtiHrBPPhE9c16U1EzypQOB+N7vPQOiptENfOPDrUxahYT7irZxa93xS2KYZz0cU5o4tIbiGx8o8Rirt8J4jYgEnO6CDWfyLPiIscAL9EvzSQ8VHibXOd6bGjJpzkBfHd4ld6spzhK9E1O83gBEqbJHHQoRwwczvLaAsraR0cYLcPn0Y15Nbr1QN6F+c2bY+8r5h6090LwjbwpFotq4eDFBQNVrYltM58EwBh4LRtgR0Er7b5hDc0XWTCG0rwoeCugglb/dUV3QFkAIpSqKR9BWt9ySnyAEIHvinmTnFcCC/Qa0mk9I6iWfDlMXj/jj0+Gw+Scml6ry59Ni9TjkjZzglNNWzwu6zlB2gQl+gMZwJjOOJEvdIgKnKTJLdEE+UIBVX57SdAzS/6hIEAxyhdEe8oYFrQSxL7I3iqqN/zG6JLONTV2XRL5IpRtXTDQOStoLcbS78NdLiZskSnx0N52M9QFCtGeVYA6TU7XhRcG53kGKoC2i7iWsTtBz5bnfIIs/SGbcDapzmqhd17ncyDAs7Kh7QjfEXEzKS8um+AlrxCDCe/leV5ODAFlGrxtsIfDANvtHBD7swcG6MtB8nxv983zveenO8R27ANGP6Ez8Xrnn06f/fMRg/BvhsXX8uYFsyVPhmjqzdmfCVhk4+KclohQ+bwgzJnSVSF+Y1LlzEoOaMrEhW0qgU+rkJnJZHwgvhMiGR/MhNyHD/0tD0LTujM3STic/mX2kaCsA7k1/ZqNCRHIIRNGhfliAk+z4OgA+44uhSkuF0wV0S5O57XbsF3uBVsmNBIAA/8g4nDOV51rCbKl5twtu85LsL+HdNmKaN2/J/KB23RoHJ/pe17fzkaJ5W3zeZmC8gLRNb4VRtc0Us1RpKZ28I1JdeKGeBfGxLi8pxtBc9sEB7NWNxV4qA/EyLwHJXRZgSt9++bwqAfok48LEJHvk95uRdBg1qwd3c6LHhXJ5/NJOWK6e/3PdTWjNmQsm8l/Pnyzn4EfJ5b9/DaVAfIkzQiFTyECazvJsRrJedGMLlNqtFzvJY8SmZnOo+/XJEo/t7VMMxkGkPYzGg61Iicpff+hrzXLc2JDbNnqXR+UN2EB3IE9gNAUrWYMTXHWekdWFkH3ur4m1DVOMCzbBhiIpRnYoiA8M+OR4cUHu5nnIEcxx/SKh9TQ0hEdcCWMlRDkPTBJZf2Ci1I5dxRAt9ce56t8Cm1gmtqnjEs9u/2RqWIMnqZ9fJgJlXwivHFQVstgvz+EnR1QA9SdZVO5eNo/xosTqRZWeHPeKk9tKFN3yEvkKsh0UaEpbhqu8F4/0YLQji/q4gVBhUa+b5ll9daGpn5F497P97cUD33xRbJjWgYOJ6qHjzxWugZlRygMEJgYEyKWyouyUfAOQqnCVSHoTmCFYDOw++IdE1NCsIFqIrw2zYWUy5kTI1RAOAsUSPJFAuAwzfiSCslA3QDAF9N5c5thnJtUVQVwBV99AgBV8nCj2BCM8XB4szF8OHA0zKK6INjEtH+NGcncah0BAVA66Fa0poISCxMw6/E+19BFNisMmHLIF1I+ZPR7ShdmURDAGxVpbwAZUNbz3vzyr/+Gd2s9txWNHFwIPxgwbSbr/3K89ujkafqn8aM/ZU//NP7iL3/K6Hf/8/WMcAptJfbMjKhv985unZDMbgN5BnRC3ESAsngrsLZM4CrtLjB8Cy1Kia6Ffpx9RaurCBciIcK4vnSJO/qRB7HqYLYW1S6HvMzK+iVBxgv63DqhOs/DogFBKsygkXUyDsMJVdITRzKtL/PFHHTxJSHAtXExqyGNIzaOmJNazqF0zu0BIt/gjCe7h4dEu93Q1OMzcV42qXQc4ojxHHdRJahj4shHxVs0cACIjstO+AvgQOpmo0lJwJ+lDIQW4vffF8DXDtjqx2sUTx7QmXlN8DzjDUiv6eJQ5yxL0nKXXD0ueKkF+RT7DXY2txU11tmUHD0LwpobNxGiQXZlPdPe4zEOPX3OiBM7IkqupgM2TamFAfA7pBbmZ4gC3qMSrRstEZcRKGWO+wHLSOcV7SmjZz61N2Aey9G71sbxy3RSgQ6n4wyi2ewg47P0siQSaFLR+tjer/LJkimOSXVi6ZAI8+YM5KnyGion67JO0/wm3TDdgEe4ZljN3+gh3QDBh4fzSUX4kn9OqouNYUoN9vtUZUlAmmqg/XWp73eLsaBfLvUtMfnZV3RXNuiumDdP6Pmx9/x7ev6KnjeGfdo9bnCgE4Tc+0Raxxqmk0KuKo9qVJQTWjYaBPe5lmwUa9+gCTxuUTlqm2b/SD5/wZ/x+tE2v+lrH9l8WV/KPPOzGlf7W78GDW5Ig7uK918qy9qHu7/DO9ozIH2geKgmtDO7NTiIYT6QVD3PRwBmMxbLC/SHxCmZlrNlLagoJRDNojn6dnk7JxK2HwI36u8IAAIEyEBnR9DtKmuqF3TXxymP8juInWU6chTKWfr1wFvMNW+3pRFd1L6HLtaALIAz+r7OgJeDcOaESAvBcSxeAZizcn+lkgymI4DSugtcPiWSrYa4cdHgoNbVZNkUIVBjoXV4NvTUM8+PPwQGvfdX+rbPd+rcEBp9tCQnwHREK+QfBrf1fBlRXOTe9jYO5TO+1NWiSdOcyGNe/5zW7yw8OFTq2N5Dc9WwEa5hLPrAB2pYCDqNURHalP6JOXceMJkz8hokBJ4ibgGIdTv5PO1RoZ4/KUNCUIHsshwTNgI6pENoMKABpyiAwsoioBQ9bYUVmffxQcI18yc3WXV+TmCWUYtUqZvbSUFTOrdQHKvypbcqirbK2axYcEWa+TX9T2Vu6GZvfEk7+CjpzW96fotNxSDoFiW+tgV4qYSJXPuN/yXFuGRdqsiw6ZSCRyrGRsUwEPmYyqOT4kYkDo6mnoActewCs1FMF6cRR3Xlq0bOp42ro0dY6RWcZGIvrhi19nruN5NvTw1J+MgWpzV6S/wyi+fT3/dxzlOwK6ZCDwovraWDMFddRCoL6HmNNI0YXmgMAazKc+IVZ2CJLqAQFIZbJa9TIqhVUrIQNbGvnD4vJ4bcTjeyIf331ZPfD7/++utvNr75cqOvpDJ/yXggzNGOFfSIuFq0j0z/BxAS3VkIifIB9YdTa1aGP27xvZCGzZoQFWNe8BLT8wPdhJlDzWGdaB2lbR8cmH2ZBfuBU41e9Ous73Zwhp2aMSindfB2tt0ClYveeczUkWjb5zAHAPPLgBO7oNp3QUdCkCvMuSivsJUzI11R4Wm5ANbnpsBKmz1g0SPNWEh5Ik5rs7eCt6XN9WRKVymf2eK4rSBPan8wYIMI420kNYjox3bg1HBxk4+ayS19k2NJjJOOj6/gEyqj/Q6zb+i/r6gNlOHSwtE5mZA0wcIcFhkB476b0fluoSkZABSbuhgqAwoINz0dWoLpYP0N0oSFXQ+kM/rFPwx4/zZ5bE8UhuhD04t8rnSRo4lKaFPo0FJFvx168+hRn2sIjuOPx+UJQVD9CUxy4uMCLsxbmKYXPqS5YAhzAWlcv2NwXK8b/+nApzRILmX0xt9hp5Ve9D/8jrb4Kb85np4AMMlPHiodD/0AivbxADeb2Jhtu7JCDQaXzBE4Upja+ZpOhHJ1e28PX756s08vfPrYkgBcpc9oZqOFbfn0vZzNl01K2EepVmjj8klEsAjnJ3iwHAsM4VNEl/f08mdw2gZGGQC1Ke3QFwL8aQCfQORKH5v80Ycs6EGhkPQqmBFCK8K+mRmI9+mUp0HvuJq+k763dQxgE13zTdnwt05g1zv0wPImS9OkHKHhLPlh1lTL0SV9Py+LCe0X4LbKXqEbgIJL4FXWay03MIC32iI1BMxwl86b/fZ2e+L24Lp5+3erERRhG1HpiaAD0x2ugRRkRLDyoppN8MpCRpLvQ27gixnc5bOCBiMJ4321PTCxxaWJhtlRGQousgBZQZDWjiq3hje9WkqVqjk0pgwbT7LhNxsJ7A9qDwAqXqX6WmI4/P3vh998NXyy8fuB6V6weCFjkKmK9GVcnp+DvGe8TpB9aXccNM/ZLXgYIRWpN97XzC4AFMWn/7D3z4fOUsTfwnfFbd+JusxJ9i6AT/GzJhZ0WjGG1UGdMgA+PlEOYMoyK2xLeDAaQ+sGCH3duxSh6NcOmW20nqma6LRuiGg/ZVMt+66YjeXNIe3GklDlKUTMjXlnK+xxuRNH+L14uffq+enRy6NXe2IeE3ezCWMoVWHxK2om7NYvUbCJWDwIKlHrKz4tXMS0bsX4QpDbTva07fg7d8HL+AJihNNJRSt3ajYcNi90wKfLabIoa+7p+Sn1THwJUdELf6yHzxP3gQq+Lmen4yUD5Rn3bBoyL6nMLkYyWjZEr5wqpoZ1jHtr8Dfay29OCaFIU/kNN0XPOnj0dmbX2jZl+nRrbtvbr2hCp+fdk57hY8IfsXY382JEp/PU0D48iELfWoooMAc0WM6xATj8zOIkdc8Qdy88OkZ579oYYoouX+7Fpkd8W6KeTUHZWNSqwgw5ppoSJttFIr5v6aflrJNSslY/deraGySXP4foEfcYl/68AknLPzLV3W0Sgy3QcdYBaoV/tS076MAMxjZ1xHfcUPD6KOS7amnMu7OqmsASKwbSNNa2JFgBNHphuoi7g6LyabKzWOS3GVY2fa+0M0j2D4MkPR0kJW9bSRBF6orqzwAdno7V40RifVZGFehNf8oS8ep6onugJ0tj0dRpSE3UjQpRlGnd9jRH+o5NqKSnri+Zv5dN0IAbUqiSQo1QEqK9UycR06tfBDpH3/ytxkTvnG0AqD1UYHVt6WE2NjP9PA0IMvBc3uHVgsHEHTkS1x2IEKRvttfyzz+yaMxopQ20XNaQJJ7dejyPWl4vlhNC1CFsT375r/8HzfHbxAfoInGAASfrv71y+Ne8p4LxBZUhMMSoIXiLGJlPuZ53MjRdl9fqCx08S6+C+3rVvq5Xd9/Wq36gRT4bJmhtiqllMdqEScRG13dGmTB1mnlfY4RJ3+th8N00WvuNMooMxtRUEF7+/X0sSJZlx2dDDBEdodWT7LycEFZIIzDUpypzK8IOF6Ce5zNRedVG5QWOUdj2bWPOIL9Y4YafEWs5FNZypvwkdluUfspOMlGmLO8C3ZjWS+aaTA/6xL08ekSwyEnTuYIMxQzjw5bRZ2AYoo2VYbLhhEgZudcYKFeQx9w1ARmsZYa/xT6vYVv6rh8dLHSMVD63EkrwhB4bTX2Bmae6xDcRQHiHaQUE5PtkpncKzK4coc0WIqlYyPNe++c5flDgIaRIYq97u/LZkIFh6/UG2sQBSFPZ+wYgkq4ET1bg57e4Ajr5sD8FGu3u6u7u6u7u6rC7bdyOR35/ofbkSCgPNpysCWA1bJNRzqwFGriIcyKpBslDJmwMzQrTvnwMQzExSxWKZv9hlrz45X/8qzU3C8y9jC2XIV9VVKRg+byaQJ3BNFLZgJEpJyopg2qjADk2z4HGJ7cqkVTyrCVIgzSpyd8x03OXk40xN5DjIgyVMFLK7Yn7hahQRABnF0CNU5nTsxoX3EO1b+RVKRtYx9oFYNMp2MJWrF/S7oX5oqENfPtIrV/LdeaXW8HqmXXjRnXV6mAMTWUbgAkRIQbl89g8Dl0fsaUms4dsY86+FGNtee3b8ru1b4kL/C6BMT11XtVu0Y0okE18z42omelWI4iHiBnrRMiPBfIGGxI/i0MyUP2bXR5aRoyXJdxsmIyZPx4aBK7Awt8Fpp6hnZIdmFrhIpXDITgvoC81DZQzqyxTW1eWiohpKJdn5QjMKLkdggqE9GgMT0wLi2KaQ35rLLGnTluAAgdLMXjY90B/7itbCLfsQyD2jRPTcd/b8vfHEqYWkOjtqxSEX2cCogEEqLK8kiEQ0MSlf9LnbrwvlkA9RmMnWz5Xw8YhLxlN0aYHonknkXGybiueMdYtibgbXDkxfmi/wY04zcj95hzBGj5jtVlajmF9OS5kfJ/TM7iuNci35PWWT1aBhr1Vqkpr+QRVo2RoIJ8Y3ENnDRiX6SJIo/ACaDx/NzNgGS19nhQX+ejWLQxRquY2rE2rcdHrZ0YfAq+kmKliko66SFk2sCbSJCuQ8ukMIoMxDrqcKgbl8XgCUupZbzEvmV08CHVMFb/BfVityMFI2yOge/YOJfZFq80nEQjejr8UMnHIO46yQuKfD9fg4EPTNkIlcTSS0hs9x3LB5M3eoRTHP5zMWUBMMs1/PDxZXQwiES4EebvH7unYeIoofO/4FD2rBSiLwI7NQnWMC540HV95OH34NTGx4y/W9YwmYtcVjxv8qB9ckW+x/k+TY3n01LZSB0AD8vtgFXX0brpAFOYUYsmppS3vc1MFHzdU5d+9fsBE9++vKxr1vS80n/sedD7zd03vlZwt22h7xPtb93YXsOTxofeWCNxK52d0CJ+K+KPXJ4zEHezomqCFpc7yj7au3R535t9TJo1T3mPh4ATVm+MjT98qZ86Es7xibY7HHicAIWmXbBQ7bOrud9XpkLrmTZrTOKPvRgJ7Zod3BqhhWz/rHpHX6Mg2uqdt0cuxrIghYV+0kbviddCo4jcwgvXjORNzyvercM0qI501O2pV82JWjEOzH8auB8tZIM9g5MGCG8Ye04rRUfwiEywvAqvwC1NvsGX3LA+4tNGfc089RUuA4HYcqTCi6UIGswCMwHFaAIA8Tfh5U/4Q9P/lr/8Tm1oDIPazP1cEN3pMgPW8nsv6tdIj3hHzO4QlsOnRnryFwil7+LgTI4T5nmhUtsXSORCplmfZmfGWMkxBm3i2JKAh7spZXY4LqZ0OZe9mYFC4qX4sa4lIqrO+UmSW73bMuMeBA0Dn9Ix5nNFfw4Lr0qTnqhNRrlpbOA/4a/46kGlsan9PpegjNs4bMr/domB2JhMlYtCHp10hoptL0UpYl8iflsXi9rCY0KdqQTXTnlqnr4HMqRYgHHy6SNuwkFrJKHePSlixOw5EGALLQBm6H9uXwMvE8gmrZWmyu7DBG7NtCkcfqK5rw+v5pHzbNs3KyXikA7En88m6T6fP4BpkzXN9KVxAGQakoVKEMa3ELdGZZ2cjyxQ+DZ8zZVmY1F1BQw3YVofGXPOQa3yoleL/zli7tSiyjP1p66a96XZ/ZdyHqrbj8T4FTmbyXR1B2VZSQJ/qNlUF0bbPWECTA89HHL4EcScWK4ZmEW4ksdhno3Xp/rM7qRKzuLQyKS+ooaP6ZoXNi0zAhkxpBVHwCa3xFY0G2EVLtJpU6Y1tUZ/bw2sTOh/dVjQ4VbwIyylI2yf9EuveiF2k8zXPywVdSGMt8uABXqBj/DW4I7cvNuTFmQWSeoJ4SIK/Z4aq4NNdcCQIZXy6RbHnw99EPQeCVrnooD5f0uVmRkigSh1wtCgQk9xSSd4YDoVm5vMrNClsl5QEs04/WjSu6RBtPjjMTBSFDDsJPiu4bx86mKU+puF6HWA49O8D7CbMl2jnUoGj50Mm+ekP0HcgcqchHfiHIOTFuKIjD4fKXMhr+92wFtqN5Sy4vGEsMJK+p+xJZPSOixok7kHYTjeXBx23D4q8Dj2BbwPYXYsXxkx80AUR4pc4w2bd6pB3GSTuYfWo7aXUQR/eMVa/bAc0wO4wM9L6LjcNTfKRaona6R5wo7hq31l6gk3GfYqYvdvhZ6lmhrBB0ZMaiDDh7u8IMIYTgkxqz+4Q0nrQWLyMEY4mQDFrVotdPh0/3yuoMdhYrec+SlQCMMcEbTWlHqoAPPgcIApEXB+VXqmMzru1ALj1ELpZ5mYLumgVrKXiKZDOqdA8JMznxwxKEf0GJ7TXBpyYSLcopA0rBwwXTKUY0Q6SFjuLKxAU99DeIEJd/cE94iJ/6KPLcm4gNUF3ojRqcYEP9uFs2TSs/mNAG5C2ROIArlLPFxcIlrOlhSPTeG7ZfVObOQ5DAHON3Ukp8bCAOMomuUUUl2Jynq0EzRh46uME2UWCifSr96LjQnnxDay0umyshYcRZysFzaxm2WR2p9FfLJcahJiLN5WhtnkyiNsCdH65YV8yK9gzFmoQuwMh8URQctVUOkblSXv0uDzlPxF/iTdsOJqanwLQaBQJAzBBiH3/rK9cT89ava1aMR5/TkGShbemjkQlAl0so9opM4GFz6NYcGIdex6zjUHvqekHx4O3e20uob0EgfSjs6nS3acwNgtjuLAhi57KxuriMqwRgf+QU+irzafK9HovPdtwq0NB/7yzUA6mekL63gbpSyIjdCZ8Pf7j3x1akFKMoexHq3TqmZ3l9nX/w/Z13/2VcnYhL16KVd1xz3Eb2P+afTLLq7K5xeOsuMjZogzOmItq4tvn7f/w+tnewd/cAhCHsMOQjMv69nJ4EdrF4U2HFRwXNAZvplZs38Yx4lZasUlsuZbhmqyGZxc9JyQururCzoq8NMTG7HUKTxgU8cxgTVmNKpHX8CtJrwedho9K3jvzC7of9p0r3m8bRNoqXbaQ3sfnZlmnEq7LMx6MjBE3BissBp/cbRg4XGEMKNrx0GPkumUVLRSSKBxUMMHMB2MctpKuV0WNwDfx9yewDv4hUNGgLV9N40OySX5WgI1ZFQ+LvwM18o9sRCCx3he+r0e7HXCMM3m/qilYFKAlFIu9zzBF00c+n9PO7l6Wk7FEJgiYElbirehCFyQbXVblSMBizbieMUOJlWUsr7b7tFShTb3qvbdlycJvCG207VbPGZLb7gI/xWN1IcA8T4CYW4UNgVLdsWRSCWOuLFevxv3VSgc+f/0qq6JiAxiMmicCxG8uCwFqVAu29rocFz21muHlwtohzB0+quur7n1706D7tcQtf+7wAKAVXghgqVNzrE1YGHa4viJKRUA0+8WaUA0i8D4rmutCBezNdWW4iy1xQ1gU5wX1MGJrvGKUI9ponoyLJmf79UKWF47koBPOisQK8EsxYEB4yjBIYbLdsnfD6IhNYbdrqxStrbhsIGMhXHKY6URe8WVDSMboFTdgEJIngtPbR23ZCt8Xk3nYBL/xW1i7pDcs1uAhZHQs93KEWVFdbrYosIppX2wFZBLim/n90WtEEAolQb9K92wMLU1JvVy1Gw3jE6eJD/WhuIUrsc/Aiq1XX51xeYV7Y2TT0en2BZuuUCDAtmAgGBcC7ECXs7rjyyfoV8tFl1Qb5jiiBJzeLugeLW4tgcRg6q2e/sScBzc+/5JpB/1wePmYVSyE8oiy7sIfHmdkQP/H4J2tRJoMhsANWHsTgcwCDUA4cxkzOufU0+QXsg0MgvZe7e0esVXoKrBkwN71MZU4WeGRKR/7CrBiby75KjNhTS/9zcTvyoIpKmjZtcBqpDY2IzLbgLW8CxfdgW4CYDqjyoty1Auxje26t9G7H0fR+52GVoJY1CLt5YsyXxPIYUZ95zSZusehkf1ucdCL+2/ZooMw+FsRBfLF57oR8KmHTtuUQudVWbj5L7zZs+/VPXKsSxZZ2eA14os1zmugETCwGTUxYgGAFAUuEvVWSZOnK7qVPMyXTfUQumaOb4R4Z1mw/oQIS1WNelswSEIJtDOsDaBl4m0SNjLYxYERUt+3+iwbCbaQK+rIelqAD7GMapXcmo9ruDEo2N+Khquym1XjElkLjwyLDJNJNK6v9W0wXl+Qo4MLBkEVWoNQ+d19i4NiYWd4Yw94cxcBV6KBJqbVek3V03ZlkKluL8w3eYvjoXKUvI9ZLz5erfXit+EVBWfzznyKhocTaz/pzVM5Sng/8guiHcKpKAlIVfsd+4CC/uSwj4Sz9q6o51cI0kw0CNExuFEEvFKFQ54GhAqzBmTlAaTR8GISnTFa1vD5lEdhBPD8wRsAT/GTBjCmdgjK3jkGXIK24eCW+KFINMpxOJIIJMMqtviYsykFY+qmuVyzW89AwSJV36YaIoCjQMrVZQjr7LHFNJq17BIfRnXtLGkMaRMeVUSaoEExwRRY6kt2HZ/qyVTaVuDejP3utOV+jL2rqTjQfcy9saXbl8d9CpYZr1UuHVQPr9IPtQrkemZ0ruj9505nZomgw0zNtT/tdHWsmR1FvEvLC4ThvZu2njNlbUtG68IOTvZjLygpwgn3wm18UCwKXtOF1F1pOwMWszsKu99CksadxRzASN8BOR2XcI4GYOOMIU5cSdxve6/yqEZgOaqA3XBDGqoq7DisIUb1qypAzu3AlyLotOe5Rui1rq2+wEiAfy4WFSQR/9///GKdxfWDpKUg9skMC5s+izozAjgDQFZb0Wvwbx2AuM1OiquC3S+M4He6hFffzaiQLBeRoleBQudbf0tYoxRUkb208mu7n0Gtpgrq8GbGVexK2BU4CD0looWwBoPeJDlEeQX8cM5Ow7FmWPvvfGsGbBVnHdOMEnOEw+tHbXRN+u4GVt0xZ5DomSael38uBTY5KdeIUyvcB1OkVAc8QZNrer5tMQEmnZ8iKMxJFMwYaHs4GU1lfJFY9dAau+rUnf6hXUL15DBEdWffhpz1RmNoqGcTZDURUgQgFg6jHao92DQyAa7ZQWpGvbfsg8TBbIogQHsHql4hNJDxWJkg80BWQ2+8tazyZ+C74niuUBcFszKbcMOZFH7kJ/XNYu0bpGYBajHRSe7hJblYjAHkZXAu+FUvCBUBGQKdYBvKoY/VSLmgslFCtnUtDZeKkKGq3u47tVosJg987V3PldJT26HcUzTYMTgtFg0PAs2SIw6vXFIpgUHqz4gLWAq7pKh/CmnUkR7WlaSnFBMSwNVwh/zIByDoIXM68nz8Z0B7kIwx3cFNfQTJqyXbE1kzI3NlZK39jyHB4Q0/pC34g7dqHV+5mt0RHDQtfIe677UXLSNJf9fnYOosiULYTI6OHrXUUhAa5BtF1FgxkEDn1fNeIFhMvRSXTB6F6iU2k+Pjnt5/VlDqz5NBcuxpzHytpI2acXLi0ZzhONpqzddRvA+YqfZa01+h/Ozq/Y6eu1Sjdh1t/YxvWbWsD8uzCctu/RMNhz+fNrCxywRSi3sonfCzgh1uWWoSjSTQ2LoYJTLzjJ0JLOOnPHjnlLr0tt5rF7LtjoajJu9WFO+7eCcrTm3nOK0kF03gomoGlWrZiGhUa3VBPG0l0FAEl1DqAB60zLmvV3MEXfqQ7g4KQ/N88DTaD1ZqNUws0V/RgdiqMbcJeYUx//DUQZfF6J3G97iz++Q7o31pKYMcQn3Q1bbUWs7HBHHfauO7UkjHFKjvVhUM1NOfqDYiImpZo/iDB2lHiBR87EsoEkkhZSdqG52Jlfi4rBHSE+ulQ6BaqL4VlHYKvjvrWGscOP6sIS/LAj4H99dp0yXH99ui9KkrlIa2Hk7On9TNr/XlCDVvsK3sVIgNVOEhFsjsO3uXJ+hUFBxUJ43QkO67FYi4GYh3FY3fDKA9fFaeDFSt32/z+KKn54Z83b3jmfWTW1U+dtCoOl/A9LjT1sfa9bTQaid++1uZ5LggarzyGGpp/CAx2g7VWV+KPfCm5BT/0WDxlTeLPasN7pcGbJ0uRBxVNHi035IceDus5ph3brGRyvWNzWXHZpkvRl8GIzqzJiFluWBt2fsP93AKVEw1jGF4wlROs8QB8lRHmj0C3kfgwziDhIaW6HZxGW8mDGxFaA53SyQz49gGD4vQeO+hy8VKwO8WgR1cyzMY/CKln4fjHsQX1nfRgVl34IQE7O0hTfF/EsPnO9rpJz8hxd/rckQQNK/fpUZIiXoKpfvcMBKJFCqBtGdA1v0j3aQ6QEM5FhBqTIuNsIQNiQluRhbfXWdGLoNpOzAx93fWZ6zpxkscJuWwZVtH+EesvmVfJCuExqIIJow8jGKcEQs9+wMrj4RkwaVBNJIU96aF749LTmRZTZPttvjFDPYFfQ+vggvhGPmtuFBl/CThfzjQi8Rq6gqb1MqJo2he1C2tIEviQNrRTjhENrYGEeDGqqp+mS76hNhlTbOTcV5XzJTtrkVWaxKXcbY6zp1s8tR2tOPCt6ClQzVndo2tqgwxEFfZSXwCO5KFyIFXiROnUVbzYuse3CW8gRnIjicPNzkDLM8GuOI4aFxAnLMoNiygMi+njyQ7ylkLd/7m6ET/rWYNYMp1UcxSaTe07pCtY7Pn7p0KJYR8sp2AcFl/nKKGyrVVNHgZyAMuLiveO9XUcqW2Ygab8stf/9dak0D03HwSrvJJisHRpMgXgHrRbTa2sa2pb3VdWCPRZDXmRJroqipEPN21KRICeyyCxtPqUAbJTqrCuQep5RIozULoD95x8h3JLvMxBwJkTSLOc4wGxZCIU3xeQ96PjKIS+m+zFZtc8/zY46+My7VNJeTEjyFVOjPJrnjlZzAWd6HUDTEEVwDpke9Qq3dIwsGpa2CcsZohGZoHlucuIrjCMyZgI5LJcC4x3gp82JhtaMueYhP+3gHHBmSR0C//9f9kYTT32SlEMt6dRGv/jh2sgiqhtMgrm/y3/yU5fC6ZmdjQvnPyd0xVI5Oiqi7g04TzzS3Ecp/3Rz/1rR+GvOVmuwhIs4yw8jfUN3tzqJ2fOnT4kULMkipXGy6la+KV/JKwi0GosAKeZ5Kj1UkpGdFLECxjc0qH+B92JwSyJSK+5rAaF7g9Jg+q87TINGK20guqtUt2HcFXi2eQab92YbHhWPcDzGH90M3Em8LD/2JBhKpJyCn6MzOTARACEqDZfJ1e6CpJeGk+sDzKi5IgsbRC28OODC7xtQ8DodPlU7NTYKM1pUx6Skg+4O/Wa06Kt9CIKcxYkig09lmc+WHzPyJCCXKoJOtJqpX58QsIppAdtRixbZtFDoPEeOupgWNg6nYXCmJTnhDN2N/4ECRByRtx+uN72nOfg3xYnDom+5JDj8B7cM1W8+oZIaRaw8k6hnKkwM8s8E61fqDWXZXPcV8Wxumyex/t6upJTv3UM76zacv5NOxT35ouAw/WQdtrlTqMta2tfg0JFUdNMFESvtWD1f+1G20jave6djrsM29WBWuQCB1+yLYVxYbZV11HILKM6a08DYFf1LkJ5JvE0WWeJmkrrrotDXR2rPYmnjd+X+M7mMCEojWeaca7gRHxcvBCdhbz41T78CWUeoZB1Yy1BcZhJKKdMoS7xFyefQ6zvJ2+nXYlnKGOaoPl92ZgwRP7gPYcPY9JT12GZbckjDI19QYnV+wZuOX3w0pkvBAvwcAtkDGvdW5E9oWe22O7nsQdYvZCMEaZAyM6lJ+RMw6XIOAvf2WcuU5RpxdHeZB0mNoEdjTEwPv2Lm1JzOcdcpiubkHLlEaadlZWk+qiBAGFJHWrepG2faHknEDQZTHNTzW3s09b4A4/3qIiog3vKtKTX2sSouFxNsyGUABn8H+jsVj3lKAOi7DcN0vyUYXx+fmQm7qjjAgfemb8anbdoohGjgxbM8l+dH/RQ4c4PVgWJzQ/lYy6HdOmtdZiBKOiNgOhe9CyiQ+w7cdd/xQxtoSsUgHWJziWdGtUjrn+iedyogEG3nWJvkMPO8sfXWcHkqg08vuwclrs6a5183CyEarCntFOoKuLRWUC2bMp+LECdr0n6ujaTzpknmhP7tpxqLL8CNl0oHizrqJJHJEbIgfjl+lP24mW/O8f7Go6N1AHoFrRvrv9RtWp1B9fK0tL4Gfaj/vr9bomp34F6uKKQ7tH21nSehz3sPm9SA/a5RzsuvkkV+F7HYNP5ErIiT1+x44wJydukzr8co1K79QmtX88cC8BETYVuCBniA//NlsQcZD8MYYWm10gZJAcRnBxswNWcm4RiMOmOd136757WsNo8OiSWDlI4E6Rsl4LfLfNOURoNPlFcTqvJuXolvOHoF5TnBqOZZAYAvN0XC4Kk+SDDjZdCslqMq0qTlXMkQVmyL5uF9F353XKYnXtNNpC4ZPaKsQgOdFn8dmaZ1zj1G/GiRc0onlDPNz7zqLu91M/M6F97eWsSDgPe5Ke9u9rqxUWncFnS0HsZFdjFfP5ToGWFpGPOp8g3GqUhiEwhvW+GXIU5JezgVXFDOIuT4mLlWBkXUxvr2Z79l9NTfiUgyWmfOJCXJVammofSreIsIj2WKnnXlF7xvzo3Ui3XasDsbMAtwOjEwIjvmVxOi5gmsko/r7G7tCirxhPfASfxglE7i4/0BBlj/teXLP2DnPMiD6HE9+vjFeHnu7AQbnrgLsqLtZwOSmcmQYfgS1Po68bF7IgXOcp/xGqbRNKYnktam5m9VwrLGgOLDpQVsMtricidXZbACl4tWjqhJamWtR8GRCInAlGVm5pSiEEdIMfGh3wSW1Cy5jsLsLb0K0V7zW7DuZcdFx2Gssep5AkZrImEJzUcKRbOkkXK/tFdcZprDXDMkqJ5D5L3sxgfWpb4O+zSqpy+JtrDh2vTXH4R/YZuc5vzRjSnfFVTghkjJA3yFZci40Uj+FgyT6RAzXYrhDg/RpSNjZTZJH2pYkyaoO/Zx6I10VFHkdqJF/MODmTPELebn5rljwC4MU11EBpgCjUbVenpbNVyZT0IHkfzZctzX5MxfvaNHtGyRuR9GsVpKJAcNPNJGx9IIb4yFSFDnqirN0yq3Z0XQkTrwNnQ37oVC+Iqd2UmP8CXxfvTLBVDpjq77X07G02McrCpwwcTNat9lYDig1C1kxuo/2BNuSNXSrduUKo6SyxzjUFpluqfvc6tTq5a6GUpaBt57HGFaVFHaU5HPE4P20M2krHKEyEMy/ZUTvV0Wo86rWkLsmJTebjpSCOkhPbXD4Pglw+HeHLutP6dGHOlTFSmWiP58S2LQ90yJI/MJgYULCfczP8iPytfTETjyUWL1wYJfZQOUMMCROawuaEJxBUL+GMJZDQi3Sa+XLKaICR/w8HBGxlMwhGHdVoTytu8luIwXleHV5ORCn5WX+FaHKz5PhrOkXR5m581Gw4iuAnzAXl75gJN/ft3dO8uzISicWroG5bbg08/0deAi+EsCwPYNuZdY50wruZCOOCA23lZ3fSdMbicaZe/hFSuC9yRrUoCXPnk9VSDjF6CuQcrpIVacyJuyph0rdCmDGwMSLk5VMrrYByddNxtvQzMroqflqWC1azqqzBpLyy1oMPfoN4wY+/4Jkj6nxsLgr/qwsfK7EibHRkU8Ccaa/EOYd37AwpaTOFdBOWZv4S4fuBjdsK65qO5KiSEBXn1J+CWf5HQeJQYSwgFNaURHpRs160PCoCVC25b3ojbInJFJA11avquljs0g1JQwkJN0GF+Zxy9FkOsxAN1RdwmBH3JOYOuGVYxg1AOxCVOZJhZSbkg+BiL8rTvYZiktXujtjBTg1oAwj7MSXPF5+SgcSdZ1h5JDYusOhRVnqHr4gD3BMxrn6Nss4kq2pF0YNpAjqGdm6MT3A9j9xvCjWF4ylysTqYox/aCheAV48DW9nSNDC6Ki7Hxyd4KoeDKSUxapYxPyTjkMPC8iy9j515e+zJdVlvWhtOs1n/grkIZ6O5GenzHR0LFqCuki/WrRzSRtJ6wClDDLVK46rDngbJChd3c0mMBUs8xFWZPHqRnaj6xgvvUjaxCVlgXsaobVTZSppSWtE7KHMx7hD9uSA2u6Qia4osrfPoVq1OcWK/dmQ4wbeVCU7sx4/NbxJUQDRRiV/thy8PjdYepLmXQ4R+5f14S928jB8yGje/Y39cCwADqxMDqh1BkQfRRSF7hVhA9Y0c75SGG0VAxRtEzV4sZ5sKkoKYrkFoC3tMaVpJ65y25hH6BNsTarXVyPYlgTs5eDtOY+gnLVZIjJhwzES5ERYRT9uOJt0cXQMaP9Rl7HARZGH68h3HFf/brBqNtNVQ7jWxr5UtlEyPfUPk/ZOTFUsPjNntCWBDeLdMNFbYbBBOulrVlm/gKdnnVHoTyzvFHpZ4tCvai4KFM+y9ZmRCOQt/bL49mjRnSGyq+bwYb9lmrY0UbNiJ9eFwdf4XbzShWW/37YDkXs7GfLJc5JO0Y/4D9ZnVY3G2bFwWPdRcuThyihaFGv5lyX6lk66WkzGr9WTaGcc/uDGkwOoG/4bHzrcWNCBJaXpw+qlGtaZhdcMbE/Lst0As14d3hD9EUHJkoeQYv0byC6GNW/DSQXNjo6+9d8YvCE5EnIHybrA5kIyXPi/WAUfTFoyZRZvixZ/uvsYjTmOwep6/bnLMwbItooOfo1Y4AR/q2bmxFhE1vQl2wSh/K7jU3ZBKNJIMmq7cU5gPgAbBIMam6lHxsUgIOenUHMWD1FfRCTcKLj7id66RWGuGSZsEZSYS02Ga3/IVJnLO3GG7Sns4L/7qzlauEcPxGa+OO7QdtKVds46j4J/ybmz6aVNMUo56Lea56SN8iO2C+ndDs9+4Ep3ggKPmqTSXg7tA0JvDPTqZ09AmmtjJOdGMeJ/YSUrCT8AoHN48rNTVyBImaRqjCUIrKe2niD43ISC91DTC1wsObm9Ex0pUaFSLEGKJxOLeBGrQ4agHTDutWivbWj0pLy45e5Ie3nZbwDZH4uQMtCFaSKk3gWjg3FsU1ihgo0RZFYQyienUj5zPypFPi5zlETxs9mhvbUx3RJSwgRDJdRfhIP09o1NRIYexkGMjS0NhlAs6NKNmueAAOVr/UxfEzyJn3LAjeujIaiC48KaofyoXI6mJdwwON2rlzcFReLEyG3ucQ0PiwjEtJUnkqH7WInwkMdFdIszBHQJMIh/q8mcXg8BLKpyyyqxjoyYI1CpxaZI50RHjJRvl2gR/SP+XCoeXBB4vxtHixTDZbufW5TQQYWIS7xOC8j8N8uIFeYWYwGYyjjPUQX72TLKItFLp4duBpAzxv42IL+y3R4u92al5xHmUrcLlqWgTWC+GxiiZrwOv8aPEj6IlxpjD8LhL1gssrDnUHs3qWqWDK3vkEaw41FxCIP23LoVGi2DxikVdK9BeERYR0iXxDFqdHJfjI96X/9YNoP+R8huX495JrvrtZceeu5TTD1K7e/Tklk9cyD4qiW60f886uE2fjztbvXduaKt3j8vY/Tvr3juvkNfxPTvXEd2MYxLevQgS7DKUvFjxmz+SX7GLnbxGa0NxUTtv0lEHoWOk3eHCuya6Fn5TXWm8zOsMmc9jyOxwRgchvpqvCXJLWXThw4G2/EGU3QaDEwnUjUJNMpm/PUqwsdRaEcniHfoVUHxVLje8IR7I6yCRLXciIgX5npjIpgnybLNCQOpERTHwXdFsFzS+V1z0q/l1J3j/EIK21oHqiHVgzcg+2qTzjuAQJ6JKammpPM9r9+5bEKir1UnEh6i9ZO0lQTDRSDxFUtRj1jURHoFfKJhVNL4VLdhYgl7JqBkN4BMqnlYsrJ1nPS9G5fktXdiElfsN4D5tkHjS2VWIZuyCZnzU3O+d3be/YuQ6BrbxVmVi7dsfTL1IYndsWOjMGQ83/to10FZkEztEG7nLGA6YwXnj8dVGfAR13Kw5Ipaz24Z5g6ONrYig8vjkrquQrrQ9sMV8owH/yujQfoUidoXlhHHNNXvlWwt0A4r7FeIydRusxLfLCV6KeY6z1HAT7Zy8M4P5NfO2p5KmdDbJZ++iE9lptHOHyY47pR0fv+uK7GtDnkFIub66dsRDupuO6xB5VgTO0REw+MR7ZccXy7ymrSh4jItabSTiBMZc8mW1EByX+17ILuibcSNzfUrlLDnSCJ8E9SvYgk9uQaRVkyuf6V559aOwPmFcNNHXdZgih7YGaowTmjIPTPa63PoiWSsB28U99sH+RXggHkRlzX/vqWmGFFsmc2RdjkZl+3DXXGy9YJGM4D7UV5QeEL4LLNXY7jRb9qblFd22hsM7o3de6DM/Cpoz63V3gAjWN57l5ZwWkiZ7Nimmm55srVcnF4hhQSeEg0WoYIojgLAS55ytadkzvoQMl5oalxJuNp/APegW9kSTvJxl/tVW40sjePGgmhkQQbb28DOp119pvtlpwion9E4L1sCSYCqpTY453rXTNxlOZH2QrFsuBX79EoWa3wt19KfxIxMS4ine/ot69I8II9fsGaBhc9dPOuD63StgDEL7scUCRp3VRNOmmtV1njVIB2UWh00YPsGcdOW6hcakn7m77TvbVjf3eNr+HRv5HnPamVIcBU44eFd1o5mm4uQ0WspcmlaUJzHTlXX5HPYe8LSAbUlHohvTFt1XhAo2bdKlOCjOluVEVZQccMREJ5FTr4YrLIvMjVMxfZ1BMZpoqIKxZ7A+rWDUdHbLsUSIZpmBNDsHDelbL1ezt9hzZzjmnAMATpAGCDeaJ2Tz0iDkEr4CZ0uAOpsDLOO4JeybzdbdWXV+Tuvwlq7vrLFwD0Y0Mm02n4elVF3CxTS/J8XNJYM2VMKasONMvCkcw4q2hY9WX8LIlDNNGyTTlfhDau8cgEINMRSBQTmlbM4j9tHHpvYgKTJ7UOFlT0cVP6TDAWcvM/nlzW/2qLCObQxWbQSkViY3Bzi7krxxxdAcKTepxtLA6jzwRf+o1CzqPN5KziLv74pf1Om77s3dllgdsUhLE3frZYPiSp79lbfplQRivPfeCz5dYxA8SIJ3BrzQ9QWTZbLORZZnOQdpjixgg6FQ94tb7wC+D4k4/mwOJCe3L8dsqygf1AmB6Jm/gDKSlxaUssVCOeMsE6jHAcPvr7cV0Z/mbn8eVO4PVB5hAlKuutRxcxo17p5I6Xzn/GMjfUuCBSZfFTGyd7ndJfUk93dI2wrB60LcV+9oNJ8UC8mtKW5evda6aHby+2J6o1h0uoP1DmJXSaALlI7IealiLmmrhtzjsJTXwqRojFGpCk2vL+GDxQk5MkRd1GFreGgAoPV/kXVMeU3+osvZ/3xdEPaKmm7P6PxplyuKesNDAR4YA89WjE3fJPrztOf7/A2QeoK1LQJ2jSW2QN4HBrbitJ7RDaVySpRtdlawwJ8hcOtgGCWcdYf6xFYs1DBBjS10AkHhnN6kUSmzpcRGHF3II1R5akJ14trJo9kwY9TAWF6OHKPJKRGMV6oRdM5MhHyEQtkSWrnJ3xUcMtAEAmFzlegm/Hk5/biMZyjYRhL8VsLZu61dw1vzsQN/8GdYiPiT13CKvcQutDtk0ndwD8OqQUvsxMDSXeIXNwzTITVEpEvvWQfjVWK1pJSstSjUg6oPCLZl4B0GEe5y8Wi49wfpMyBagm12Y7R4idrIjWfUjd+MO94Cyfk03diIwOfk5ayp/lgW1yDF4ToJ93O24WUbn8v8qhRXdvYvVzLcz04WZfPFfO8Krnus1JI9CkbMJhOg9ZQ4/mbg3gJDqBYLVexN2BYA0rrB7tbyLxPNy9SzTKHFsVEH8MW8404sJfE0FVoRdStqrpKQKiAlXhNl2d/qJCLM6PpcHjjPmGGHKD6rITBOoSXvswKCP3u709U4NykqKwm0ayZ6xzQnpUzzbgTYRnDvUceFc5Q27kd4yQdZ0iCTa+mxzu1zh/Kx9Ee9pJmONDy85wH7wY9up6ZaJnqpMCaPh4+/zjY2ssfD7CtR7NGWjoW9mlQcdBcBfJOXz7ME1LkYWeR18hDNJfsPk1RgsW23WiDrnHMtxRC/SKgBNveoKigTaonyx03XrMgrbsQ9mrM149vL5zX3uBBfFTCISyRzqnkew+HGwBOdiP3SLbsFtyLoURM4p+LfopGSH8iDrlzPZl8pnT5PfB1YNfnm3CQ21QosF3PR/r9k1Tl0dRxYhleGIX4p2rLNRMIBeOHAiFnUi6TjoqscDxVXmoV3TxN1SY9LAByzpHGm4av6DJft5H9aEpZqbhFHU2fvCQttsGFIO4L1OD7xmcefeE20qeMZAzTNci8ObjB047iTaGOQnE/yiw53Ou+iMj0KmPhT9o+7p1ye27JPJ4EIC69YeiMAgmurPFOgw/rx1l8GJ+t9kVWguHFPgBzyYG/3zcHzl/t/OH3xaucPh8yUeGWoFvcp/sjee8fHcqLF/1IsKuN45Sr0jg52dvdOn+29evPj6euX+z8cnr7Z3zPYXNfTi6XENfv+tdxRy0E2hDJZ4RBha2IyPP7yr/+2scnCgVG1WEh2kfPJsqL9G3EAq4uqMGWRXm4gFBSiUtKq5TPRpM+KfMGfs+QFDUJDZNPNfPjM1ny4pbdMxsRW5yVS0C7y69mAzb1nVYLcJVMifUZJfkNFlPaaickhzAbnk6phaNPX26hksVnG4NSn9iDSwhPhWLsjdnzSVwnclTqOsxPzldjGJt8maxu8lEQzIgYcHXjPmv7s1s+WBckTS5TwTk0hLtXscYqwK+AlxEfz1FLMLAmfE3d16+J1qmQcMItF/tI9vJUGVo47KiRFMRHH5chFALWGuCYtMUvZsBccnM8OPbMrFh1d71Id957/cLBz9PLNvp69w53Xb1/tnb7cP9o7+OPOq94Jr4w0g8qnR3v/dASmgg/ms53DvVcv9/dOD5+f7hydvtrbOTzCwT09+v7lwfPTV292d16dvhhuetY8oq0ta2fwn4PwXozN2rAalc+adV1wkdEksq7V9aKdcTUtZwhmgq3i1nsDHtzdM9v0EjgZqbFV3MSamWoWa2a2bIa1Dv2LDsCFduMje/rE04Ai8BuygNEwXtEZ4qu1KCass+fMQ8F6AWNw8jJvJawOXDBfYI0MUKqD6AQs1OtKELG1CjycFWzfYU0LgvECrFjNoZ14lhzIuJGa92JRLWG6DAdsm3iJcHcJhND7TPPmylETDoIvuPg02EDRPUaJeJWwh5vlW0R0ISATAPF7E+pDjHuMlZEGYTGX96IyNtaQCHKCPiY+jLqEQA+DGqJLZrX1owEE4bThzteG6KCLiax8AK3wQvKBBwzMiLGhQZ3Ei3nwS8ITP7sVOCJIEuE7ZcgnkjcnLKyFYJr5wRAhVxxhivWYo2wPdCzoJEd92GgOVzD3vKK6Q5UeYvl+vDQrB/poOUMYV9S3UUEFIvZqsfPe1AOY17AhAL8qCz5F7OE2OeXaO+AqUZDLumCxQDfxz0QKa2EcgeWo8w5aYsuo2RGGzkZD+sEOoWd0xjKWPnfPvIR7OUA6DXp7Ebxlod+wz+SZobkhejormBfqhRkNjxfa1OxEcqYUM1UXcXm4ZhMZaLvBS2hepC2p66S99Nagf0feHImEjskDPm4WZB/jFR8dpk1UupKun66DM0/ctUFQAiGVsM/XlcQOWhSFdS1SNMT8nsE8hbGhM6BIaPItD2cYNAmalSYJzmhUBFeFyx5+/+aAEcwnIBdYWt1aoMLy8U+CvAg6twJMOqKGCt2LUWxCPtEzqVtuz4dsWC7hR1r7xFMPN4rtu7lUQJCDQvhHoaTTbmq8K+S2S8v4edpTQtzmWYRgUTh47yPe4IsCmeBG4pzn5+cCTYllGZgEpwWHUBq2yPW7760h92NeY8sS404E8d52/OiRfMd9N+yVXO8DA8x7fW9gKN8hazGsgZDVdkwEHJHDQWbPd1I2DMCWkF6TTzbDebNB/BbXy/j7x42Py4eD1D5r7XPAZSLFrlUvINejiVolOGPbIrpd6jmFC7dBuR7escaovsEpcyqe+6pHV5fT4th8grkmh2X4ziAZ6bod/5NWy36SwJ3W/UzHLhrultGo7Sqw7fQHAHJA8EzN1A8jU2MSvtLmNzBF0aGFNHs0su4ynsVMZwHPb0eEAMKDWI7l/a9ijC0GWM31akcN92K4XcPq2jAl91H+LF3qGxvhcKd2OJBcbUTkYihn+ZJN3j0L+bO7Owos6ST9lIEBuJVgvfUygOMOPhJ2DAv7ITU9MKZpe/0oexZyoX0zOdU/cHWTdTYSoJtqT62QW98MIA1m42IcSkWKTJOqoRZUHL3W0IgFGreiCH5yL48EwsD0X1J5yC4bVyFtxEEYtaznh57J8MH7gN9Ks/sRA4XZMVmKbDBAH/d7Af+YnM45rXXeaFBAU5SgbgZhBScRJq5YqWq0INK8fLFQuwbCcw2L8YhxuLg0PTMN4dMcWQK3IuHRcAjmJWfPYcw/EIMRHpQlZiVqoG/mYdqLpHZGRN66squRlV5OH1+J7WZ4UBW91+46Rgd7K4o+K6OXDB0grdbM3iOpKq9CRzoO75xxmTDOpDQQHHqDat37OKHkJeQN29qlyL6Pvqd3KSBLXSwaAKqQXpcA+2wSC2kORzvO5e8LdyhgDWuWl7NqrJKwN5esZL6M7oym1MD4QtXyZcuQ7Kwa38ZzeEbvOu0VDNWiy+INEpYAaMqf+EDOvsBsi/AlcCzVsF4/LFraVnkuvWc1CJReRDZcEB+Z9tiCXKI+seELavS5rE1Y8IZmyfc7kuc6rQaVlsFh4dMYzMTVWsXL2axYfH/0+hV6+xb6fdHXbD8cV83DhFNnbD90vP0mBkNU+avjd4QJGST1GUk//O7bdVT/jo0v6lE+L+Rru89wiPEFYysjS0L3sz9X5SztJf/x76wy8YkjVousUCL7PJlQV4xhJAIS45rAILNwLJHELwBUQFk1tFKjs3qKrcdKzOkWURfFjcTk9xyn75bO9uHAAaPUsq4meRNlH7xX7wQKnuNa3a+XRsG2XprfhonvadZrMmHGMhG6hVcQfiARBpLCeRmeEu0jDI9rdk4JaXYK8zAUL75Qzj6mi5Bbq8LdiM3MVq/N6imwYgp7GLG0yPgkEfJQgBGKUJwCVIP5fkyWuZBiWTketoDe5Cdkl5ODtZ18xXCdnzwPwa+2QqZQVK0sPU81r/ZRfobgBkIK+brpSNOIeQTWFvGuPdAxOubUQUXG+PfbCKFYx0nD655+bZ8c11+iiVOikaOeB+07NaPdhrGeca8RjStR+BsVv1blu1pNq6Ge137jf+pEyUGfq8lyOuNMeSqcA3E0LkV0vBDzcFEUiFrWJrixYornL15gVwyoQtz8vV28gTH9211IaHq/o59/ONj5Z7ZRWeS3YjiAAjsHezt+dWSs49j8O0d7wXvYt+PDi4Od13uHfA4kQMAg2d179YrfyC2jXl/+0568kdgUPSOyel6MSoK9tZsGVAEI6z5wnqVjVvAqqyILtClBkR8PHER+MsCQYNbBthIbg8SbWPI4EFc939t9+XqHR/Q+OablOeEEC8e0dPSLWjqmdaJf1Mgxlkk/Y230J5ZDi8oCnLA84Zinrr9l1iciWfjMJJgD5rlelA2dMhsCTibVq4VqSB++LfJ3bmoP+yC73eQI1SyQ2oKR2MPfPcRhcbqjLW4P/tsStQehH4PJ79Nwn5/+sP/yKFAn0SoMcHbkYMhBkG03+3fSdxJmAhWMW/1TS7hvMzmWN5KaQigwl5VKNlZikghvYYwGkKd8UZxk/GeTQHmZ12xYEPSg+kqqWgDy1jgsHP5eXABKTTQrnhYi14Yc1PewGEFBbbUtqSbJAruFEMKXRN6wytQW2dSauHe59s8koMTZcXoHwi3ikSzsZ4mYYpXoGyuDcDiWqdoxpzIxVhOal0aBSBP5HSLw6e0YSFz2TLvnmOBiOAGiA8bwABJu4UXQi1KXbASSHmqCD5vh5PB536roJNgLc1dVBfaLcNWZgpfN5CEGTkes1s6nRbMoR0R811fW//pU0SBeDpKHIlJ5KMQT6CQrhWTGUsRc2s5W8jCf3T7U0GOqCWsyjhQiE0JQfbEK0RGItQmWSc1C1DWdjUjkGvBh2bIXK5mWNZQ8kibTOh6IpBeXyFs6G/0mlbo/cobd4ALtvnn1w+v90x/fHDw/tPpUmgUdfkVfHruEv2zsxV5MPU6jZzId83uTnYb/8pswT43+5i9vq/lykhu3Ze9JWmI9g3TBv/itkS0Zpo1/6P/vhZWlm7IgUuTDycCN38+44R64xWdRAg6kzrGvEsnJgXJ/YEKHvsoPHeOciCrCvjJM8yDzZrknnNzO8rrkce6yi718Rlgfcc7mCnaoB/n1KU0yv8rLCeuE0HJ+re5V2rRksQFvw6vKj/KJA3BIvib6wk+mvwVdvvLnYszTrKVv8y6Rd0FJ8dCzGgyZonnwR2wFyIMkeBDsgSat0o+98rwHV+StIIRThUFqSPi7JHgGHvM6fvPy9OVzFHsjQXYSeuLGdkajYt7IiOxvXR8Wc8D99VShNcgmTnGkwFsK7t1QwVpWwBxEqCI5QWNlVhVDIDD7+Kuv/WFc5vWlv0I71MSp0goYET1aMkLOoOVcvWLPPFV1q8Z+kS9Oi/GFnCGYwrAqFC/8jnXupzTO2l8MG5Yo3AMshl2T0+JG1goVX+m0TRA5KSe74ho4KC7AbtAinlKp05yO5a2efffJ1bDjDI+cdb0KXifmtesurOYPN6zpfXGVGYLVAQgLgYnBV0zxBOdd73l9Wvjn2r2/qxn4Yxc5TpWWPlVSSI/gTAJ9e03xocetf95d6fD5yiqA16dz5E8Y+zXEclLeB7dK/wAL+uBfEG5rXYNpMk7fxB0pEFVFDEmA0GcGyAg14+SrgZEI8JhgMBbJ6Uo5HBBu6MJgBDvEF8Oeht/vGZLPV+MYcG20sYxkA4iWmFcdlXl+oOFgfT7sBafBI0JOF/m1AddWcXX4/I4B+ZVLhtUGaN/fgD+oeBDdI/TSuYUNKc3fvlLqSB4igUOIl2ALwI3k8P7jY8SshZxvEEOapEmONVNH/n6BJzgdv5Btw0NPyfVepMMxO6e51HaMsZGpp94lns69uxl/uVznwUIhvclpU53qcPgZmqW5PDNLEbc7MFaCdWISYbz3qO0WJWL6ePxVq5vHX/0OF/m/R29fDVu9fTX879bb79tz+/3ffm4HCGqD9vnKsTWePKxofYMmTCP5Jph3Vw+m7QPA2Y3h6TdDPisAoBJCgztbCNUnHOVv6M3h6h92GUPzH2FX283CsZ8pFAf4fvz+Nd9o+bti8pIhnFi8y3xyfv/80djpbjGrK8HX/CIZLWFnc27AvpjSyvXjZvNG7bFaAzzlVBCnz0DMmObgP0YvkzN9uWLkd/WycvzfU8HTHzFn9IantWt9WtENT/C//Xvy+L5t8kCtEGk+mBbQ0wLzh/sHjk5HFFe1DpV1bI0Fe4OxhCjApyQe2yHpH2FDu3E2I9ROfOCyhj7z7CD56kGKFGCoIs4rqjtgCsbWFr2Bn1m8CUrabn/0YhZ7JWwO0x9tMFvz1R2rFdGvXnshjIS2scyhlgpOSitYFlFT+sJuULT4G1o1zI478IyuSuXxWlGOeiy9LurGnrDWvLpjF/U0j0Y5E9f7EGkHAX5eaaROeuEt6YqwR70XxTUKRxH1Y5pgZfwwE10s4a/3rbedZFcQIi82UWTibG9ttBFf60aEmW7t+RS9sbcGcRpaJiGMAbz37g7SPUha69f3L1iYwpYnJkb75g031UpjGwxHX/nsfkiE8TtnX9gzp9kms9XD7J65RDunLaNpa8eOw9kLjmOQ6Zbhib7i4HYReWZT3fJ9sA8yOJPvVohI8xCuK7PO/kIoLw0EKSYgp/O4rH5wzgwx/32g5vDRcXvxSafsD3TBzE2n397s99nwfwJquePKmqRN4Z2yKxbJV5QlWnQJWtSqazQpkTjBg49AzKxDl/GE1l9eOeuWyUCWLXkAZoo2uOVcEd7VYWbZk4RYllmsiYK1/mlZwKzl1IC4U+jagg7UHIXwhNeFXZHqrKYNoQZW1RNUGYzPgoCOyhFAaFVvI0bY2Ac19YUAN8u3GqA5CD1GWEdzHwS8k4z4bf0MLC++XwmtJBKCMJtxbcznVXS9mqu9mFRnfD16f+BfEV1yH4btYj+f3YVu7EHwKugBaBW3A37FeImodOMLwqjQrRPoKZUOdMH3JzpYacZZbbdbs9/ua3Fg2c322IK1FpHtoi2bEEFAACA+Yei2mWjcdzXjWDmLJryE5Xy3Of9eQFx99OE3bL/fpKFWghajGL+/5n4JeWj9s2LKMPG/rDzBbgs7qU2dk5BHQVf71Uf3ETRCt00BrzZBN1dffBxh5fatrpaLUXvXAHSTs9vgEIjGCcHUeZqbERGNEM5KZi4Q8s3bNboCHnrrOlaPg+G9cbwB/7Q4l1l/50+4SgDgQNtRlGZUT8glAuboqAWYLGc+n4PMzB638Uk961Ewc4gVTDSaqUfpW1Q+kk/F+FO79GZgMbRttZj9qvYcml7GrIlt2ry27QvzLSLbEADcMeQ2IXQpGkMB3Z/cbKRuossCw7TZ6NY7Tom+ClfGt7i/r3GGTwEbxZcm/40n9bK8uCzisxpBOXhrSOBqP1EPc9MSZ06irUmGvglU5EjMh0CKvLDtaQTrE0zn1y6Toy3n3EGLF1Q0rcKlTmB+J5gN0UMwJe1454dd22NAGwZysxYCsXKzsKW6JYqgty1xhFchJpVRvEs64VUZgV0q2tX0fXc/I8N2dayurZxffIyEA80RjXdZiTaUKuuTz7FOyovSqGfNb6OrMvrAOlRKh+pSeyuvlUmpcP1FxW24FPMq0naLDgBKk48g0FDhLKSd2+qWj27AA1JBK11kqaPinJbHuYh/FFUV6EA7NaJMCcycGWmP2WiJl6pvg8J7ns5zr0vL6TmXDpLwKVbCahJmK9rQZ9lf87Foc436ycs8FXOP/7DbWZFe31HJ9Hjd4majDh1Xa6eCDtv1vP4CTtj0dB4xuaYfZXG13bgUWvWZYP7j+Vkaw5o3Z3+mBeCoGXvUblnUqb4q9NG3xBED8/T4PJ+WMMWC1eaJxLQ0r4wLJBt0Smm2/BmouZqYp0UmadKEltP+8xrgMX2ftKpv8r/iNuCa2nTBcfgLshig8f5Jv6/xMtWPduf1y1f/zGaJymVuqjH1wBjUuxdsW7Vpkob3jAfha1vRfhDzK/9NXczY458a3I2b8b+1KpoQArrVrcrx96gBZwm5E9oOstGdmLQaE6yM4ykDy3r2r2NEOsHlfjejUwmzsDoy4VJvpnqOUCLW55lV5bjLRR17tHuGXikshAbSWOjXju44+CgftmPZKHaWOOHUHB2v2e1UmmJHEi6S5bNb74P4mHLjNvgOP/keRFDTwXVn/bQ2IQ1NbRv5yTuJm/LgeYpTtQG7ELadx82xlT6esq3wZnR4e70wjtRRy1zZ2z0munKTNZFpL813CgosMNgMrOoaOiH1m3MscbAJdPu86Dmdyw/3e+Pd3Hdhv6wT2DW3cp2JVSns9KVtgQHmqy6pHdDRzrNXe6f7at3sX8gfoX9acyStEQj5d5RDQOljDxanQF2L200fiyXGPjq+uIY2WAuNJbkdY7Sy6Vvdrdn3etd7A3FnQTFrvme/8dlazgSM8Q7/YB+Ss9CWr2b4AVkpN8X+7+cTvut2BWjkOfyjNn0555heljy9DmBz6F5JInaBZiugT1Da4vpwOisA074JeaIpceuou1XwqqNeq2MHzI5a1jl8I8S5ko1TJGoljJ5FiDgvFuzZp3trDHpTMYueSWwTCxo1SR5MYSNb1J1Xr04FDfIhlftvvPQwki7uKhzqqizFal9gfKDCWvHsDSQwYAAmtjBQaspJ8tC65wRQ4yEHjmKXz4G13ebgccEUbVgJgYiyr87C1tjPhlazRj3k7KG6BXRd4p+WyM3xuJ44JmLZu5nsttBgJVcWybet2W6gr7Yg4tPn3z2rtkjg9G8znXj809bIvaMbDDkwdIotjAIbFmPJoYYgajQQLSS1o0OwwPneIbhVk1F0mKx1CenjORuyKzTpNhsWGnPHnI5nIdzib0LjYI+hCXmXbgNhz5jAWg2s1Gz5zMBHZ8/qUsg7Gl9oP03No4FlFPjBLQYwchOuCRwAynlAcDHw6NNyBtAksdwe/vKv/zbMhsPhQ/G9gDeMcapxDjMxwaetple+i4UWD2Ky2HBzLlAkQlkIbYYoKhpX+8omBDbtKWG3/i9rx8Ps5AtDtnG07L6pxuG2a8LbRboRxgvkL4ZY661hVWmm4lGpA4Xb4+88v00XTZDDNrFDrDhBampfKKB5hTmzL3tKsCeOupjJwM0qwnlinNeXIkibmeTGOdLV6+ELgDT3aWIYzemEDcRfS2grttd5atd9YzhMvkhQKlmXjwO4WmF72fHTrLEGFBO09KMJjxJG4xK4CP96F4YrjH2iJfoST9CSzPqaiMTZ2BGBCLYuUV7sle27MCj18fBkK8J7uReirI1Bz27VNlhTurQQv4u7ZGL8sLcSbDapcrDG2jpTNq11WPjRYWjx6TFYuNQxC+k1R83KeOg2Vln71fF1ZlNV6SeJbhItAQd9MEHaZBIQUF7nt7V622KKAYJX5yNdEvWYKriVcXnOerfGkZd1wjKm5ZxAKYOKgboC+eTEQOJOoD3BkgN2tuO+jeTNxh11jmU64JMs2VUKps3ABEyLdp0yM9G3qTC9RRDvpE0ie9TD66F4KnFQtSykFZkCAuelZGFLxO07qJ07D6g4cEuvFiO4stb4vYaoi8qpZ5W4WkrAdHZIZJcnR+4J4enY7i0vjm+wbsqBe35k5jRzKWYgktSgCzZUvbVSV6SWaTi6mJ6Wy3xswgZfV7po03xseqbWfiIMoEF92MeDmJBQmNnnEQCisSAGy3BWwZITdI87du7IyUAR7niGMJAIbMlGfaN3IRX6w+vXOwf/fHrw5kewg8fHPceUfaptdkDLraTNjG6nZ10G6Xh5VtstPYFSHcHIVhtSdttMdg/Ya9AY3dJQvZQ8o+WCVf4fZ4XbsQL3KSvyO/tbrb1or8ivMUf1qn+U20R3d5/o7bBqqWKVcGjF3alt6z513TT+q3tVXauH9tuVuitG+us0uqsH+uuUuc484Y5rHCwm3+P7lYdd5zTmfLipj9zqjmY894nHXyXWxh6OQaPL+3dl9fg8T4mvhp/ecHejnkPE7/8mo73H92EDI2fvA9dZa3Cf7OOwotX26LpdEjo9EFrw6JOcBe72QAhbvtclwHkAnJwYehD+ePT9Ocq8ZlSsnJ5lR+iHeGlLV5viWQ+JXJTJG1SV8UNkQjqQPwm/R281JYHnml8HmcCdJQC4F2PKmCU/Vot34tYpbuHuVjJBYMl04U3BzpSjS0vuSdK6sXHJX3OEMEr4xMKW1ygEeEzvyAqASjGxdlPlsWoNj2Clf8ZmQcTTGhbDRPxA34lKB93KififOMFFJEfWJsVmwk44FDIPEAjBUYK5DonJU5v6fVaEM1cKlcMftGMoCEvTCotQ+FtuIylA62MIVA1Fj1iKF+WV0qos37dzR0SHgVIGVmqymezYA3VdLWHhVNDt4YB2GGOkIOYtMevvJubGYQJJQN/AJ5jOVsURBTAzqwSxCQ6Z4rzMr3SGW8lMqFxPoqsBqLa86A8g3hH5IVxa/+RwkAQ5OxI9KDf0a+o1AxrbjyqBhFKScAJJZUdl7VQvViQim38oRH0U0jliLRsJoEccIV+wZ0TGsFrTpK0tNVJtxH4OXDxtn5HFh1pjO9MlkIQAUvJ4ptE0+xr92deJsVe5rSEtpJrC3MnEUD/1InRGDUal/fZFeyKUvj1MBPq0wwfon357phwig+sekCshjjAtA4XuanG5kxMX7K/heLcyTInteazhEhKkLPMaRUlMHU2yCA2BM/3x74wRwcRFW/VcbTqaOM8nddhGEO6Qz/528uABn5MVgVQ7x9X/6NI8hP7AZrXUVTD6PmI4Sq5TJt8mT7DlGJTJSBooMD2dpdyu3IoTWGghSTzrlToSvUdWDMZKhC4dpglK1043jpKaci0gKJOnXTSm6F5Oks1EtJ/tdFFQxPH9FMzFacI1IxR/QrIX/MgQcBEZ2k0cSx3syu8iNe+vziMlwrMkrCzZc0Tq7cU3tIprlOYouPWm1FyIVE62kh54+bjpb5Mhy6Kq62M8y1lIU7MndHr0J4vxJHijE+P1XQY7E2bMJPOxkWc/BCAG2Y8YAMj6i4S9DENfgLnGup3mvD025slr9YgQLpOIIZbknwoYOeFwz2POiXzSDmzonZ1BkmVYDEkh5tMVcXhPEb9sr7RmQPh3RMbY1kMafazpYGE8YnlhccnTBMJajoRlRNu2kH3xQHN+P40/bNpIW8fyBTIZm6HZRelEmzK6p2pqc5LR11HepGoWEIzWng58lLwwqJx1nlou04dcWVu2AYP8mOg6JI6Zh7Zc6DiN99nVgx5tuRYDTgyxupgzL2gBAtA3iFPKdRcsqza9IxMEzvtIDnuvZy4AvzDLvp/v0/z2JZi3fFUFhJcpoq/zNukRjjEzrAbnIgvjrAIP1Q6wiq0UPww46OWMMzR5YM1UvJFo5mHyBFGD3/ARljH0bVYlainlAKkSa1eOiSUckUZWUKKUieOv96lDQX43rZynCKd5U9ssmPBIojfvaLlu0MN4SXcgJdBzxhPJEZ98gIwQ61RmE0tKx23MNb7jNCWv8+Yyq39aNOmq+ukZ0Wzop0/km3vQVtN3iLrOJ5Ea74CEur7vZKzoHdYQ0007AhrOujeMd9qUB9e6wuAywcob2meQpsSCwpWBvXiCciUs61aiIThpBJmOq86wdnQYsJDmmx1aRoPreEmcjX17ElyzB3rPYmRiAjZKxvj3LLLfVMAi5l0cDnRT4S5IYJP/zIx+kIimrAXL6EteO/Gzas9gjRf2w6Lj1ltrw/ehK8m8W2YLnmqT1O16pqSMpIy4QFh8NkOSteVM8vY8cGj78BrgatQmjCiugKNI+Yr4n2KasS8ZColj02T0MiXsCVEG8oN4dGMB7Sa9qkR/y0vW/CG6dHrLBJt8oadFL5J7qSc7fFdMvSy7drxy6syYjAZRWnXKWyEIsefu9PWS//i/+WiYStRDRx0TQT31th4LiUiZnorVHTsOOYR2ea5hlxB9SFxn5DjqiRSEn/Ui9HVNZQLeqiL/j1l46D9XrX3HEpl6mi+AhmN6oZ/d1c22gQiRRjjzsB9iVbGf8rZ8FdW4SKkhOSv6jumi4I2lkOi6GGROMEdug4dQ1MDW1FNgEyw27SPWUfV2wtYbdaJE0w16Bt1l6u+awtgELps1FYvLdvO6IEb3g81Xb+EwFgW7IQh806yBDHrEwx1lmupQv8WFYM9EaJ/LxjbDkGRYPw4uhtGP+n1DBTNl+GEg9Kc1vNT9qK+Mll0IBaaNr7VfRy1fy4urQEeMbbzyyAZOvKWn5srph9ePHw7+NDsxlqYcuvkhH2JnP/qQ7UcfPpTkAw/RUI3AxAYCE9jRn3+aSRn66yuNbfRap7BFFsEzzu95+EddCl8vTpvIK0Iwd5qkkJ0t2HNmYUzowLAJPB2wMq/isFrSkRGUdMs9Dqiz1BFdYMw4Xr+KxNa4DROyWU12UQbw1mbhMcsXiVIyDNquLQ+Cu6MGNIqwbw1tR2vWwMvN4Eyc28EuQ2MM0wkzotSPfwiUPVxlLeExjZotx1wD5KiRXTk+Qaock8zEI3btRbeMXH/Lt+x/vXd08HL39Pu9V2+tQZ+uX3jXIQQkMnVgdCcDiferoaUG6mPG0UTHotxlVNwtT3xL57CBCCw3JwVl3Yrq1W0Jagk/VDb31eHzrBdZslkQ5wSaLW07vIcQcFRMa1k5TovmZeIzolUukiPRZnne1BzzOEugsBCJsqczgtmMdYjTjrbUR7cEfcRlw3BKYniClUOAufMKIcgkP2GWvKhCqQ9WppicD0TIWSQrjJHNcsS2y3ZRui0aJCZ7LWoHqwX987JufCc+J14RdWTZXgZjYaOJbJDc1DoMytg+tBJvySiPcMQjGed9V95+WbssJvNW2hfvYPPVMXkjqFo+maxpVT1oa5OcOJ8go4fkD2NJjRXwdUAYHpRJI+JUNTHEAQKhi1dATreiSueFZ6pBjckCe7lFft7sisYjTlrG38KUZWj7iqXMSId1mMHXJ8jvM5rQZuFsV8sm1SbwuOijuHvkXFiNKSj2TvXtbPQcJWieI6RtTbn8HhJfE9x/8tVQ0flyTkCoMFFR1LzbRtSH/LdIreFcOFVuLOVc2n0O+68z4DcuFbjfDreSY3DeiXNDNWdt9cS9s8ivRZKDgLjahBg1EC+8wFD+8+Gb/axm1F2e38oacGBrWXiV2MuD2RUu9I/LYim5FNyjLqQalBEShyVVKrNJAzbfSQfQRn6dQ6EyL9MeN9ZjVz788kbyIeD6VEbBLcgIDrSgMPcL+6QD3BmBj7dTN23lHJZVVp8bi/bV4ls8K1kb0nn+ArSOQbt5XYL1dYRir5BduLhYzgzMynXTgD2S6VJyu9F5yetCYS/b3xovuoRTjCDfOULwEDmYI6GivVaHmaww2L67FuxbI/LRCvF3AkJBdY2m1FHRfHGGrswdmZNjNlxWk8f4bHluanFWCe2Cwe8f9b3kLeMkfkUjPu2ri9ktA2Hc2hBCw29lDZlO06IPHgRUTPrggb92ZnLUeQMrQm2A1yW6P+FVkzwZ7gCuuN7AxPZqi94ivBQ1QYqSs7ZwEUIGOB9r83KEo7xGaJhznWTjsgYQHgsG8suEHzE71o5gNaSAZJQ1L19XV2Uh2aO4FbUlChuJivuQZ74ozmuXSkse1eSzrbhQxoae6vRVfosI9Ug6dqtZ2Lm2EpL+Z37lf7WZyZD3VyeuOvxJbshZMVyEEXOvZt9OieIOjW5IG4WmpeflBVapGJdNtTBKfPFNEa05N2XMyRdExbLmFDFdG2lFvXEktAA7Jy5rvat3LwnEU32fqZg3rPvEvHtIj/ZlBqbKqWe4BJ1b/ssnUFgzNCQr1U9elBcZ4RCY2rPDmlAbVNT7bQmH8byEFknvOpqFuAN/jeCGQSZVfTNngVhHliXWTnEGMOR0h2iNSP2GmEtk/XLbrBOS5ZZ5iXZK3hANSetvdA53V3EUUFBO429LsrFD88j2KTRnneOqGnSGN/ha6Du1ZHFJeLormmVqDYa4YIQslGXe6AVD2q0m1ZITqfakWNDxSL56XdtSH9r9iPFMHfckr3c1sRwDTPQhr23COA6Az3deSqyYC7IKXRKHhHUiqJLhN9fqWlO/MBIor9gsZk1lqxriUd+bNKLsPXJIRwFeGQSXXzbFNA3qRKfDNITVYSolSc+okXeSzu3Dqjkxke6lpWOKqKOEWdYtlxyI88K1k7rhpgjtj8SEAmPlWzdHIPLm7yHKIUw/yWe5EzmzO1zw1mzO50jzZExv4yYV",
        "KntMX1zc5smKCnLbdoZRMZGca6g35KYULvo//p3gTzFPvuQUbfQhbpXLht8I8qg/y3CQfC2LZEI4xrPp7eb/8f/88j/8v8nLaX6B5tIficEgTHHW95uUXC5KKmxZVFs2ZT55bghlBHcxuFFI4ZgTAY8Gxck0v8HYLPnjqgRk0xB5EMuJ9BC1HFKmdxMSFvGztSEcjm/8kWwMzDQFU7fBbzU3lt+eEzRWVbKiyudaMW4qjwaybHO+JaEj2F6jRT3wXDtqubyAcUWXYkHqPji2XQyCeidOA+u334/aE4qlr/Ok7ZoSOZx6OKC7a3Ov2apfFbk4R0z688usmk0qTnhp/GQA29wAQamh3Jb2wnlT48n3OdAS4S4mXpnNM83XC6h7e0TnrXsVzD027cW0GWACh5vfdFsjt0T20+AEqW8O8ZzOi0eWLapS5b1zSLbqQLw319V4S3UeP8eQpVBGhwoJ4tP0J2cd8zSZZ8BKR1X60/GQ+P+fjjdOQCPMM+xD+BpNEWAkogltezzW3GWkO8zY3vT+wcJ57mZIxV7O2LHmdpDcho83G/S45n0OngN69Ce2jtXm30u7csfKWXozxPiHQCG3wYfboU5sS/qyt/Jmw9UIPtxu2BofnLoAHaAtVKPyH9xqEElr7IZUWaf4XATb03I8lvSqnEpLKQTJSQUPP3Hig5wSR27LPOgKs+BQ3WdgkCmOL8AsrD/Ik+Z2zk78bMAg+Zs4GSkavmODzOmQPWK/C04RfyN/ICYbtg8V9Nat1NhS6Dh9h6zh/eR3ppbgGjpXo0VVYyhzKES/4LVN1njt4R9IT3LbkkfbUpR2Cr/TueQV502ikvrtVr9tyLcN/5uTQqiVGi+ALOVx6bwCea/zszpFv32o9Yq1b2CVcgMN/hNqkD9gHfxntsdKzzI6eo8S+rMBhf/jAREr2a28upVXJ24oVitT/qz7IcORzRCrAqfsR+MbtDjoQ8wLqEl+vh32vVgXaE1MYsTuQ54J6OIHzWfYz2rkg/QtF6gVRms4KYcoti2DMkTBU3k8Dl5+9x1SyG8mj4c+ZdmBYsw8S8DnGFbzy1WwukZxLtGG0pZB1tLGFFCLTwFtX4ukLmXpwgeIaK51VArKTe8s2LPdh3WPhLFsyhGL+XkEI7YpBlt4hhibOasZdo243NBHTZXYbNzcZtZTJ1fTs4dEHDawlKCRE4zyxTiQFn+kGOAzURP7uURVDAN3YxVUQJaln15qXOcU74h+Gw6tfDQSfUglTxBMOIpeYdI8EMYSYysVTtznbSveSYQzuE+wyK96fXMA2i3Ti0ZgEi0ZHsD5WPJdCYqek8L455TvBcuIwDOL0McQsmw2QeyK0iXYZYfz4rzU10SeLOtbvsUt+dMZfVEDzLZcMRSdq/WHtgZSyTXR7xKlS4VAiIdL7ibmy+gOnQTwva4xpnUgErSUYYApvx3WdxMOpKwiZD9a5LMa0kxDuFlebb6ozibFFP1hnI8e6X59t5086XduWD4e626ZVK+IFD65Zb7XnSIRT9okrIGWR0cX6HfMydI90e1Q0a6/A1K8kPDzLI+T3xiQ/IpTcfMyyaeDIh/fgsbb0TxepiEIlm6Z6nvBIVgXy9kmB5ZdzpM8LNuz/QRJdudw3OFEzp0dXlTVmNuHhNkbqwMaD1r1fCKbBgjZnAYepdMDiVkxZgVkyvjjHonpdigx5e0KG7WbphhYJOf3clAJt7EzmagAWUHBRwhwpYe+0UlYBZZ3iQNblyubXwB2EvJQSIur2Dreu58jMa27sgIN6/xKIWN7D/DW9ZsxInIoVuydP0lW7QNgBokj6BYmkQzYB0x8Cy86C9D7HRPB0smwafHmPj5qHyy/WOfCBCXawhY14yZ02V6xrsXvupG8/v5VLGdQtlwQ8Vv/8tf/i01frPLgaZIGGoankiPhVv286EbgCmt/uGYIi2wvrqQZ78fvey1u/qy64RAHRhJDvOXiVtJeVgs64mnv7yRMFOKYD+R3zZ/1ARNlt/e/M6kb1kSPBDULtX7/MbxHOeoUFUgbNrs4vfz5ow+30QZnEHZQuem8tu3xgVxj8f1HNmcPsWmCasWU0EfApKjIKrOXXDJiPBdttlJTTsjU7xjFRy610azIyROUYp4AqM1vD1T798R+bt8SFOgN/M4t1oeFRS64EkSE4gwTOQZcgJvRWQ4dVN3cTorsWsXpluV9Yto3dVk4BtHf7wJ3KevC9bmQbHC4cPjFJ2Lea9kQvWmdLfMxvNHt8YMnCd5hRCKNNAadqvzjjNmeqFIByiSfwxURHwjzypMlvtlyQhAjKK67B52c5eNVA5f6Hc2WU9xaGoC4POmTJ67c9jX49wyAkf/KEdjG8aE1lFCDasAt80P2NDVQfD0Q2di49+u2sItq0W1Vgwruf/1fdqG1SPLRT8uylqiF4lxqzBH9Kn17DEwbMPnkBlzclAHHmBN7o+pGw8QQAZYJxPYIMjV4+uWv/9thJ2n2y1//9ySVrOqcA2dWNP2s69B4K+qtIPOtQqhhquLkmPbv2Vwm6FYsag/OHfUl+0cvGxum6UKiB3JaUCIYlzYqmMQ1E/y2ybeh/c3ZeP4pM2FDQ5NZc8F6v3ruagTA+u/fcLZ7z+GIHA7izsY8+jkelUdDxwVw6PGd/7Y/8wbhu/yQgFx3XQJrQKQ6WC2pWpnwPIsdmZEm2iBSfKJBTWhYjXpgzO2WdbHo1ZIeXYSNgbBp3GFXdPmzU/aNM+Mk5AgAuI7xJdlOxpmhOzIXEWvAgahB33Opp6rsU3qd6JN8NlYvYC4x8JupR5fFND+1bCet6AZ1GLAfOHLiXcyD4N7oFOFvZgwBDTsMWkPWZA1r31Y7HYmruipH+USMPYqFtRrWa4omj55ohYAz3bUJq0Cor3EFEvcaR+F82qSPLn+m0/LWeF6nX/bZ3WDVJ9VIsgMDa4XWRWmogWHW2ZORj7E4HJgFMIc7U8aKWOwNuWvmzbeyYqyYFGdfyY0gWNGUoq6qZSPUpm9i2jMX30ac00MZlpKjJw0jGsZABNW5F9eLzShqyKeLcSBFEQleqJszC/s9jqGKP6D2pKbXOIg4fB00NvODYFdo7nYfObwyvUA5FY907/FAxwlSRpqlT1YoTGXWEm5LBcNfJDZUs41jLOb5ppmnJvMKYqLCLx+rzXsPC19/87/mfe0lqS2RbiTr3F1nuZpt33ml+xnTPTLgp8gfohYx6JbPBdYLjgOyEtvOuISfnwKWY2M3Vw/v7/sK8+UE+OIn7AP7PBJHp4JYUAsc7BBuGe3P8nFLPjnC13upWGj3spyM6dRI5ETn/YJudYmdhHnJAmPLXSGcVlPsTTgXA8jlphHBEgx4cPoADvStvAxwxsVlxWkn53m50O8RKPmBjaX1Wz4e713RB5DoxYwOMPE8BKXhQ69C7eImHzWcoiqNOSyrB/aOe5fkIxZkwBmLl4uKEvPKq5XSaMSl04mSK/raaW/2wEfFW+2iMVOysmCIvuQ1ztUbKiX6fPVO8lF/xnEO0/XjP/1p/WS9n82reeofL+Wd1HhyO7ACZGhQjhQeKP+mIGFLv/kW1pYwAauqTdpyHSKIwkRM1V96+LRCeBLMCJ8iVLUYG9oIMcn/T9u7dcdxZWli7/0rgqjuzgghkATAi8RMQlykRJY0RYk0ibJKhUJzBYAEkMVEJiYjQSDF4lpVsj3usd9mzXIve9mv82yP7RmP/db9Pj+Cv8Tn+/be5xIRCZLVNZrpIjLixLmfffb124/hwsLMjYpXgHgA4vAwQ1fFKCtJ9HHkOU45jj0vpF3O5q+JaL/Bie7wDQmqjB5cGdTuYN73cockL4xi+zbmfLtBienaBvJovF+pgsTmixx4azhCNJ5LC5DEHKvqmevlaJEOFPDi3V41yTBFcaIw23CpbYQa9LNf1xGKI64GyVWiNlnaZCQmY45cATLnji06EdrJO6yejLGT3SfnlN+fT6plP/txdsFIivkotfCYVtOJE1MBxLHlXGFJqd205bDvOw5+drQM3qQtTb/4UcUGhGAlCRUovknTG1rVnC2lvG9Hn13j3M5IpZvZxj/zPw2w0WUcqnqLc2URIuaL6QSYQyAkj2qVZtrTl5gv/NzZHkkMSfqwV5gCOoG1iemNHA0GT+HWioP44lwkOd+CPwjRToIxW6TReHsItr8GOWC/7Ko+/oDVvR4tzTa/x79R7Xg/pBDB0BXpJtuxUeSGA2KvklCDgJ47N4ZegGANHddbgMP7AYI1JTjVwDUGA2Ey+oJ43fdprFMEjpIwHoMEq9zY10H2ga991LUqnFRrquxgGqNEAUwLuNl2t/E3u989xZXd1skSesFjBHuHwWsYCfFVQhsz77BDqFw8SOnVpZdPGGEcwQtbriAJVP6PlmZAWH3321Tg8UjiS35WxGb3G7ZeEmZuaMbmCfxSF68o/J/onJfgpPO+neCFpCi0qf9gw0OSFjt3M2+EE9ac47G6EEsJd4UpYDpcisW87X4wdFNjo0oDw9mQToXHdAjVn/swIY6PIi5FLUykVPmaB/ev1wZi/XftFf6m6WigYWVvQEh34koHY81We+AIaL1e4Y/TCYnO/RPzZVH4ha9pMlIFayNGbagvO3f6n7uhfdRZe1P3XlrI5lkanU03Julcz3eqe9Ne7zpvNMtuM9njez1L2tHTBByRkyEt9iRdUUCb0U9Zoj//dHee7Sg/jI/5gwNsO8fOB+ZCX3bMfCtm0cg8YoKdlP5o5Dil347ms8bW18hgcfySNAQEY0YyDvJD7nvEQUAeZZT05bge+YAAxWWT5iUhiv3ahSm+I0rADrjErmigAH9YtID95pFXqaFRJnoM9Y+vhz+ienAzbuBD/wPRsmAF8MAOhxWw33EZ7iNfQn4l7znoUEB++hKe7nQf2I554ajJ+UQRFOGBH1FrBI3uNnsX9aVxRnvB+nIDFqdqOppsoB1P4wqu6NfuUceCBr9jiKoMTY7cjc1dNxp/4kRMW8HR8bEout3dINs3cQQO/YtbEGo8OtKAFNxN1ojwqxomvZOh4+EXWxsmnfXnyb+SaV7xXs4ibog+HUjlry3jwdVkGsD0vS/VS3He7g4sRWz/31gwujqFdeHIIeT/MiCSjqfZk/Hvx46l8lG67uAVAoPKvGeO+hIFlVh/Eic/R/6p8Whh+JE/EEg9HPeZYOR7mNEZYZnnc054djy5QHK0Q2bqvsQuAjzTZfYTyAtir4m4iM55BFfATDTwFNhPg1OIr69EeGsgDrhP1CGH9sUbC2WmW/BsHi4LF+Hbd61dq0BaiwaI1qHyMg2a+YAD0MeO78IvpZVDNgLGva9BVW8zcghSin+WAWNgkB1SD+5zY+DB1n4paTEGITtGLuUO6Zz6zjadjhNNxvtm10uudICqQnIlhc0wICxLSIW3g7A2ghwBaF24DWsyEGImyIglAC6scnMtI4L/iWuJgAf8tXEwu0p4qoWtc5fTleQtcVKXYKx4XJvO2VOEG/ztFhIvspvyqw2aowkQswAw2oDMmVp10y64HHeh1uMrg4nQ+EJDMDnG3I5l7q+B+lUU5PG0CwQ2mgNhWFxPJUcRmZ2YfgSocsthGucCkoiwsKI9Aujbao/OeipCJaAAXCmGbkvjEahNLH4K9laAtcnfZtGjMoa1cTNXWGpRxXHpt0Fq9qKTxCXOowNoaxY/g7d31+Ot7sfb8CWPTrh74HuVJS+4hdq7LS6yKeAn0SeRsj4tyOducZtv6JbeKwzZqfHRzUbZPwMcCiSlaBIQLCG1XmoN9Novu2lbAB155/2kiCkqphZlEy9bZfcuEgJp7zpIYaUJ8y6a8GcqKPZLg//eWZHopyUap1l+NOo56omS2J2g5mikzkme7Vm9Pn2OfM5qNTdS1oHO84b8Bg1g8oWGSNsdEeE5IlmnkC/apxiARJhHy6AEDClzj2nU2m/nBkv0QjyW7iPHD1Tno7xDMSPTG9QbqaYH2g4tsULlkdtO6sUGRY3qIuKm9GHdjFKcsjwGf4wPmr4N1FDe+gtZA3A8Qd/t2OODCIZ+4o7NMCQcOYsYp8MIZFnT3eAerQiirHyYT0VEhiwCFkfbNAhYZHoDi12u+jYnqDBQ7hSG23uUQuYnHFs/UR2gzp0sXpAbjQWJ+TYztkttO/r9A4rW7q/+gXFQdvgGsJvW1ZgRHDEGmdzoCZ30PFOTTIalIpUM7mlobQNj6tDR8OUwLZZEPsMpscn4dVXd7SvwdAQ8D8w2einzsS64pRKGIjiDzMdMxZ2GXnNs0TWcapfhZwM3hpbrTSSyxc5BDB9lwlnPBRSpruf+wZfuxN+/efDl/dpJel+it9HxWeEXpD15YMujh10ef9ANqJB1Ct5A929K2y2fINk/jkjK9qkA5JrGfeCCUuDYyKysu43nzH/NAxZ+8m20d9Jy4UUhbXQBpQt5wfaNWyzMZ7BaAY7OxAw06LUqYNuF58eMYAR4q7SO0UjRxpVFZ4PIa4DwPt9KL4VPlZaEvPhBi/i3YtBxkQZtjfrO96HvpC1weUMaBkCM6YmNrhzJ5mfEgd8jF4pM4N+AZE1gOHOCr27gfvb9bFyPR4bXxuwQiDWk08nl6Uhw+05HyyBhoZrpbBFNwkedCnTV77pAjlaelLR8J60KrxOK1bpqcMI6hE63VkIYHgSXVCMsMc3w4nbHCfNM7LousMBpgngp8qWAavY4xoyq+p01zhbKr0mZnbXo0POJYDl+GT+mC0fUtKE62oUdHfWPXhDDY1XeJMr92FoaXMgxafmrdOV9ylTEP87rBTWrTROM8Q2h8NBXYfllPqICnyUi+jyQ8I+oIErO2ta37Y2PyNQTI3tvr9G7MkpSwYwNzdbLJPfrvlcKtv0aVCyUPgn3CVZgMyQR9cx7qvEyuMbA0ce3mgIQY2YEkdXNksCtdhklfvni2a+fvwzRwypIJpzwp+rox20NPaoN3EcCd2HhsqIPGE+J28+8I8+OqQ5wHNe0KCJYDMeS6IhaWnydCXkbjGf8LVrB/0psr9AbUq94NJ7TSCvALEmMJuxU34iCOFUjPiXoUWD7DA5XbhfFzcU1gho0G42wpqou0nhSSZw0Oh1rpkYYcZpyW9SjRG4zsP54Ls0mgW/Ic8VOGh852y8NZHwqzp6S/YNIIsk2m9GHUIGyMA0bU0G4GtqrluNY+6xxWyCbgnb5z7cLuaWdL+Avxh0cdaLb0MPIj9Gb8eyigckF7B5hQKOH1j3bg/ezLRmKW8fXuX/rLWASv02xkH8pH1+oTSjI/Fhg8bGw6KaDpYJ+DsybJoa+0oXlL0Lwui9Ne9CkGdjssrqCxKuBn7q+FDv8C3uoj4YR8fAzrG+SE/LOHwmxpNbW6TinMX2rxnA2jra+/nAnUP5SQ4xPWqqWuMZh4HxzgaOIZ26iWPepwCbNVeF36Wg3O0e6aRH4Srgag07jjpOSlnxpaIciqpXvHzQ3asvFzk0ckr3KVtQf4q9xaZUkY30gicHyYg/v9nWrRfKLNKM1NcjyjcSLSwLx0q3+AEnPZU1to4lz1vezLvjvvo4lBLVQCLXW4xsoaVr5EeV1GnPkuRF1V0uyzwlbLope6tPGdZAC496YnOMOYOW52XmUDpq1Abm29vYZ8tmno8mR4F2PprOLk1M6+JOyGeHXcOE6W6sv6kU1BsCwz1ApuBeOB1tLBK2DCg53TCmSpMAKWhU3L2gaHLyI8rGzr/pjVnMqwvZS/vJotljLGEa2sxYgcegZ9NXDp0+Z1SYPDRf7wm0qb/fhtRDu3ITUMACBj7TVcryNFpGxwp0fhUXG1sINSbaRkYbz7xVm+DjJzyvp1Nrva7ikjki6o1feByX+3kk7j5Y0KhGDeUA//tFC6XBmOXqPJ7PZnC8FLmY6Prs4k7TA2bs9QBON3WXzaMmF2o9jRZobDpiJ1fxkZKMjeZ8w43bs6+9fcffRPaCeBUBFmjmAw300G9WaZKI6ShCdcUIWyWyyqg/suD0Uekk8TQyljCfQ534VkBXmZXH/49jgihnkHT3xStFqURhJcd29z7ahzN0wGJcevYb5XAt9KYWonA+FgrjP7suUN93efbN47fYoz4UI271ktmOhR6Vl9ziRe/LmGNC3+LsKcpI6/lMy4uB98gr/OXagVMA/3XZst44XSeuyHyN9jpeyu2iojFN5ukczJ0hUU5MLKdK2vjbncEmSKretESpK//SLrmpTV1cnkguecIF0vDU1c2yl7kd53EhGxSJHpElllGHd9qZ3n/fF7os4+P94ImkN1UvyyWxuPJ1xmvom0T/yK721ukqmF5+UBmOPv2BjakxbVw06ddfW7eXin2azsw1lAWIvf2MTPgk6pm1Yv8aPbUHv48WceL+wwrUkT6ZUxisLBREWZkf4i0j5CTgfDIRJVRzbIbyEAFc9m+aJHGeutHJNDIMPRyphmfOfCVgRDNBKeWvXKGLVuPvDbd20gbm3g2wt3N+YwW0N3FJY/u2t/mZWu9O3uY0/Ikpje3xN3TfymaRAmF84eUz4IpwczYhQ9LPHZ+eLZUggwL5pjoiZEOIkK4SnsbGWPiK0zySvOHl0GD6OyVHAttV42eup64Q9FGMmCDUReF8JjmHAbk74iVicZCcR19Eg/8eGGNK4BxgzM8JFK9OWFmilM48V6iEUDG0GnGsj6fYzLF4voasK0AjKilXFy/PJxdzR1EgVqoWYmLo6MyNUpMHjHogJsq5xTNmjqwCDLORGSu6CYxqZ7KOOPRQZ1ULuGOGLsPqi1nPLrwZE/g4a/7+l4+Lfuske9qLH9+XxZJE8/VKenuDpMDo6kLIyxAjVpn9wWxO4r1kFuE7JFcxMtbOpsMPgqUWbMaDMJUEBpBKpuOqZZ7ZBNGfzeXoobLhEITKsA1FUo2quChGx6K09m4oWBRep9GqtZboeLUAxvpsdjfIz9z8hvOBcH7uZxIs2ygdDyj5ALxkQQBfH/V6EEk3LkFTg6SRHCTrJ1uztCgJLbFot4j5+uHAr7H6619V87BpEEAs1imL5AXnkZqH/YkBdup6GRlRSJtzbOJYKdz8Hs52LAF7NJZGNW0tExDXUUXjOhSq9exNvm3ksoRelGACmeuljSnTFJZkJ0yaBQqrX68KDfFeawtHnaNEFV4sq6mzvzn625p2r6uwNERY3PDddr8nObW4ZwQ+AyJ6oC7z+wAT1BN408aNiT/lX2GSguNDROQoS70mV9IcrYu9vrK5MPkwrkwaa1x9rwCuKJv6QMnLXAs4xXISiDySjDjk5DjZf46ut7TWoU5P8RaW5k6hblMy5LiiTUzEp/DQb8YpDoLvlRy+G0eXEu88x9Ayu0/4QXBkbLztDuFNbnYNadjGCXAwm0ofd2azh5cZB7sQk0kuP8/FZ3hAizbENX3V4LgJL83S8iNV64o8ytRyBJu820RaYHZl9we5Z9TK+6tAO2CYbF37c/LvfHa3/tWJoyBfFdR3CJGkgdhovrA+1ijSlhquvZUQA2jZPxMHSlIE3GuqegOanTtNxp4zBU+0mRGzGV4+ZDfhB5upVPJjNDkUSEm66x3k+hpltCQzSZpH15pOOQsP4aNu7PfSEWLBIlvOXCI3jtRiATyR9Q3Z44Yjp3PttBfxt3qCSNYzxi/G081NDb3vbJjd+iwbfzfCJiATep3MslMuM826k319IZAFO6NW5I+pjOKoQx+zmxVTC0CUhCj8GkicQ1i0v5QV9jjkaBVnKnjx/2cp/EvUxyhrie2j90UIRpGOMqc3iRSifbE5NGiV+3ItwC7PXO8kkDnmED0YCNllmJxVUuQ04YhQZiy4X+xPV+B00Xl8vEhgOjyqAYntjaCkWSkaO3MeuAZSXdo6G1vRYIrt10Hy4jtwGHZa7CJYxzdmE2elytqsM6y1YeJQTV0ksb8JZK+yXbapSkVL/awG/ZpRzBCqmIOixYBqahKIu+lquMKvYV+GYzuV1VYQa/CccwgY0ki2jS7ND3Pdy6u7vmLnFV+GW+LIF/XV9LV82JjLgq52PRu2aGqNL4r8a+3HPmnAkbr9kGiHiNDAaopkWx0BJBZxM8GsIg/I8+gyineRE7fIEzBVPmNqpPG7OcZHv//XfR1IOONIn4yt3wW6bhJIAc/WeeNgRGwRK3Yxh/xvThSn27l+pUoVPh9E6dHC9tIbQxauU8oD3RwuPyCjHwHdYAuwVRRJoxZ934sRjDgI8vSDEw5QmYf7K6PazR7DJnFXz17UywFcjH4jBUDqgolSvR5AdNJoOl17fUIf0aCD0XQCJJZlIwCV+/8d/K7YIwVr2alpg4HqA4vd//F89RjEQrRzfxGqRXNENHqxZfTi/OIhxi10HmqD9hmwS7iPX7sU0nEYz0WgWpH72YgRxMBuH8BZDSAYrKCvOpuJ9HkxCzU1naAKyk9F6VKDPQH0skQkIHI47WRgjzDZzwwSoh4aiZFVZUq8Ei6GzPoWdiGvrZx5BQLEDMLHVhSo+x44kC2qBZioUEDKhUUu98s9laHCUEnTQcI/y2vWYnNnWP/0DV63OtjdNbXJO++whPK16K/d1W3zU5CY3UiLcXHX3qGt5AI/06W3d/LuPnVXlXK8/oEXbh8VHaV1WCAoo3fGbXlSTD16IeGFF0WLepowewFaJ2NBTszTDh084sDrtRxn700lPwVtHRC14lUhbEa8kbaaAtxTfceqPHgo+yZwpraA0n0r4WpKKkMXBCZl/ylyYjVb2xbdNLjJWrL4rs6fffv/Lxy8cGb+9SUO9eR5Mg7rG7d/5CJ9NE1dmOKZxUG1I6aSRRkKO9nqPXgsn51e+weR5CHYTdJ+EXRbSNSDH0KBBMvwZpPjadd5gtMGJNJ3sfCReE1AoNIRQNJDy5+qc0uS9ytbVj8AzXvMLf/p2MgJVLUL2HSdH0h1C7+0ifr0Vvc7fmOPVG3eIF44P1U3+FXrY15bwj49B4LMHaWOOiEUJat6qwnrguuV2Rf4qpBy56SauvWw8YsS5jDb2gpkxMwEwfjh1dxrKWtDlwicMXWiuUPGB4P5FfYh88U9lNEADcQT+5Xh6KC4XvrXZNBC7FWTvhrKzTQmD3HHwaEZixGNHn53wyDOWmwMOVW+4XkDzMJeA+bA9RD2VfKbHODI0ruImiwTQP7G6ynbeSfYZhRYmDBsH96otBSuJtsGGFMBW6M4VUfrcIjusco/Ui9bqHL9tW7q1FrP0Zn/zTmyYPnCUKbvVcYBKDZf73n0tYx6bTgkmZqzo7Njdq7NLTOSRxc7ZbrVR67UuD6GB6CLJh6PxxNQZEUeuaYpu+jHezL57+JtXTkANVlFtzVpKeHFrlsBOnaQ37FglZo29DXDHxL02SAZNblQqAgf+v/0v7v9nz8FnkF15/z/9B240iXzWz69VCKuCB/KcnK1WFqaOgwFKy7F/PxOooPOLs/OAV/bumvuiLYxf6WY+jlJ8Hrd2PbwKL+je5m/WBsE8hsuJX/uiIRwvqI4NpijlAZTapImMP3qh0lMX7UYnEp+P4nR6ra1Wtml+TIEXyGrjR7rZ9xyFZDLDSNy+C3sOW3drc3OzQJYYttTYkVNEzKiDiFqyAVPihYXfuxXMKlBPr7lcwJthEdiPMlNKmrWBD1AeidHCcoYA+de8hXhkXvOwQMzWYQ4bX70uAk3FEh2MYQiAXf51qXIY7Sa+b6/Jx7tehVyHmU9VqxTZbVsx0v1qtIzr8dv4CctFr9oZP5qH/7Ue+zAlYi85cPfk68j50t9Jim91Ma1lmn2CdpETo4MvyUCkWtxA7Vus6L7ZuMbBHS6qo6uKt6uux/Rohz7pnEvmAVnccBXJy6Lhkud5Y9j70WTMVTfm2H+kHTLFcbzVk75+mW3f2Sy6equKLeMFOvakngS62kHZdjBaXI6AEG2Irkp10GANBu8MsJDuCnOM61ZW+zmRNrBK0WzLDpFX3O8bfvt8KYezaO2npLQc4G5zUeDQPMsf8T2eZGj9pVY/ZdUNrlomK6HdTfKcZEVSmLgy80nSQ4IDgShC6Q0tBisVkhOkt5eHmsOr4FKTPG5LjmwQYuoNaXr4FwOmEwA/ALy4xRtkBxfiyy3YGgzc4x3g5VLHguCucQSD+YphEcygfnByuTe5G62hAq3MxI7s+Kg4e/Z61vsDNF946X8cJx4UcQLtkMKnaXeGgwyhlFC7kEtFFSKGMvYPi+g7Ht+QOgs5X6K3kn2h6HhmyrJ4DDuZzxPExiezgzqwv3F/kodNPpkFlwveoZsfwX3EbHkrV2Nky9QMNHA1Ox7BKEhJjSoapDtjqBleUiyldAdu/EJ8fWaTCy4ChTr4JYntE6G4IuappzzSMhOlg75/TYGPO4HZihsCH5vx9z8TMJSQed5UlizCMp11e0j1N7eohpdP3FyO3ZKxoZDLvQGgol6qbPkzcSrO259bL/TZN6Pxyemi6P+E5HMKsHk0euNkiOfQY74AX8Q8y6lNkG196Rvb7H9+xzFwxw0G7pgBx4FFv7sZM+msw3E1d5n5Dv+0yVNgD/w5M6Ik52In2YbuhgqMgH0RwuzlE38kLjmCeAUfJAfD/N/TLnGrKUcrLTi+KV18ohBkrY4Azr+2Lrt9McHUs8+I2BIbkSX8wL/SwUu/4Hx27iTxcQoi9/vZAb2OLzOFNnAPrBxIlKBrqqonkSe9B5xhanJwOVPXsfsPjndIt0DB/haD2Ilp2t9e8ieDHP/2Df8eTQ/dy1+/+BZosbMpYmwiihJkLfHn1+b7s9cFfEeJEnqZPcYVkDfU0Go17UjR12sIySBTfkC+CTzNk+aZ1hjuaTXAVPIW8Y44Ve/ZRqfnQAvTtXNvi4jYbeygPEv2QZCGvrJQQS0VAKbF7ZLLUjr/LiGaTlAIdcQE9mgEOF7tALwX4ls+sJv+0KIWFUALKxRtlT6BZL1KMGrI97Il9jJOrq/J3BJtYuceJ5/omizSI5AMpCFqRteSL+/7g6r86Y6G0qIjRyNsy+TUfuSRvcxSGlFmp0APaNOdsA9YANjz7t/mGX6uc80b6M2IZbWaVA0dKAX77o7AX4BS0Im6ewuXgYwMrBzzDOAvrZ95WfT3YCW5Wf1xc/D6Hh2BqrGbshZ9aJrzXErBiYURfMw1+ogTJ698OZnMNrGLDm/nsuE9kgngzK5gpcI7Y6SaTz6ZxMT9aVCDAWcWREG68y7sdXlih/kDR/maE9raW+kJ9a+TExqfSV/io88kOLdfjUbnYhQSfpxZ2QG4FhC6ofV3K+1YMc0uQnuDCGi3MsYRmM1IBcG8cjvo0gt2KE4D3VGhXJsTJx6a9Qns3e8vaq+dHdcSmi9erOLhXZca/BL3RkMAjwyR7Vhz9RifKiF0YiNlCA0DbWKbiIBXS3+i3ji2jGGZSJl7yWvOXXFV/VqD0awGDqmfgVGeJvZIwzL3umj3rZg+0Bx6AvYXABT9REcfkdu2+Yw8aId7SMSdN91NzsSHVhQtXfTzw5ozsZUhFk0F+gQj2sT+oCf1DOfWdsJh3jJdGbS9ngsdbdwtszxS5u01hWoq0oOvWenVIZGk39IQF6YSh6ou0lOCst36KE1aaV5GVKjB1YhDbavQSnE7kmmS0KTXkRxhak8RrdTIAGNG/EHUdf1U3/rw/WMfu5+Q0ZU35bWXix6znZit8OQ3opjI3spquMn0+riPq80NJmUCIOPdz54/fPHw6dPHTyk0JgWgFTT84Bu5dsD93WJi3Bzo2/vZtuhYu6+iJlm1jnvN1qqrT5fzC3alLdMUnfxJd2v2b5PfM9L6NWvy+i3Iv69H5wsKsqnCEp7zQyapEGMBCF0dwQmfjc5mjsAdXBxBK0gDb3XIrGTH7kY7Ba0jaW6BhUY9a6szBOSeYtZ+Q6shygzvQCei2LFbHa+N/DL7+vFXz75+/PWrhy+e/fr7r5vqjfTGNgXjRytAkss4vfL84Qp8+f0dmpAe/bj7+GUXIXSz5PZ7ZD8N52Kk2QCKlinOj/5gb2s/HfxGlMa++bIDeqQxz643RXTxt4dAMf6LwhTcobPJRCTijU5iEE7eRRvxBx8LTXcmj2xB755wn8JDapDtSaESGGlH9X4/W3tYe13tGr2AJrjBTuZQx1H5Rb8bQHF2GeLRyFfSbhOjUHwg4ojx4BQLNwC4LKhmaeCdVuPE3/SoiKKZWCNDM/URsc/cQ0WNBGqfPTG0yCjq1HvDrvJB29vr/WL7i3t3j48Jri/uZeJQVkXIIe//+D9n0Yf7+3ScDpWawklr3NwcOVJH4GkCGxGA5heOGbsrT8ncHMxBSeGjNTocn4/pvutL3t6UHpmI7ZpUepg2DEZ9uqhWNFxGnc4kr/tRuwkG7NG+bEtEtQTne0WzElZkrR7zv4+oq1HL8aRatOvo6nLzyxNs708aNCOh94ddZ8hxfoKm1gB9wekpBUOY2KHpORKnmDI7qsChulP1aKnRjITFOHU9FHaVfnVj0ZEGP4XPMriCm7/cIFtbwi2RA1grPbqGVG49xCmn/7hwyagwABpJ5Ak5c80WR9UsHOKmI+9eo037IeMKmx0fS4gQCgpegjr1AG4NYx0obM74vA4AjylCiKYkqNS9PkQq6UgA7quN88o8GEna5tYNp10D+9Pgma3TXhI3d/aOqCAtq2RKN23FIB77eO/Z96+ef/vV7q9fPA6wZf7ICYaULDFNbjg20RdTAKL0jsZH5qGPb679glkDgM3fuXz+U/zvfsKTcxM0IYQ0bUkXDNRkNgNTq1Pw1P3Kx0FqvHGD790s4d++7LAgd/rZiY9k2JuNEe73XZOOlcq1j82z/FFj5clM+KxfN/dYjP7xwRsPjteQeLP6UsKZCuzCUvABxcUM2173S+mrbp4NqYt/2zHPCJWTUg1XFyKL4XoMJSt5uYOL8WQh4ZM8eU5K1EjMDqQcDPArnK2mmWNMM0fOFLzK1SQLbYkBGwm7raiqspDiPM0kuLOTKQ9nRuSOQlJEXq0AiEp5ogZ1xOR46kgAfBuBx5h3Y85Wo0cB4wSAOvD+vqbYmOmpXFVpKkI88dATmvyQfUH7qDItjne9obyQvN4HYI/mu7Nz0BAwmzI+OOho1lySJQyz6KgQi9eoMNxKUWVBCY4hCAhWjq90ItUOHKNjoaRgOiljCnUBQIiQis7RYUxawqPBdV/LgAZuTEaOaTiSFIqBU0sIrxzJrPlFZHRuvQpgfi0NceRlEJF0cw5coQa54WUNcTd1j/BEzecdIkEWiSdPxOfmL8aQ0lBC7z2N9oOSMZ9mf5NtbfJAbTN+kL83KSltbYNvnRI0JSp3q1XuFpPhNstttcoROLcWcMvFaQr/3OSXV3PI7uNcfeUjxGClQt0RHFHwybEE238Hngo+Z4rIKRc8EADAdtD799DCXcfqsCmIguTGCBUvcS78HUHi5pF6JfJF/FJhg5Hly4nYUMOJhhB16trkkV6IdQugwDVVJw4xK5QXq1xnrDuiY7U0jh2QkNiODWWHadFjz13eYzc0Yj2Oh2KMUbrfg9p8x5TVjZ0PJdIwHBfv05t4pDDpduIAlOJ8dJ5i9fZq2Hk+Zgwan9bhNBj38+Mtk8N2lEE6vHTb9qJaUEWfivJGfKhE2xh+OLXi2KGN0BQ4Ii/7PXOL/6hD1Fgft3uewrewEV0jgF2NsgLiJr42Mh8yBb3WCqXuHDH3/HEgNu2AjmtodBORWgKorgv4wriiW3CgDxDpGm9pqUjJfeO17e4QJv4R9KvX2PjSgKJe9Z6K+1LMbr7/47/rtWFv4rl517h4vJdnQ0MPtK+lP6axXjLMMcsUXYRC3kRugemMV9OP9dFQj736DBjpyH1zDgx/xxUTFvDyFJlXSMTdtG5Q+QtZENusjvTHrr3Obrrn3swauVBG9DZvLKAQhORGLxqn15blOD3FsjRRpyKeoG2+75gL/bRTKdwycq4m3pajuosoJp7oXR6Cfwnfu0Q6gp+lbDiTtwXzUYKGZhBA/ipgVIjMdJwK8eNpkP/7mXeXIqId8CGgcR1TVe0dugbuVh6Llkyyc/jKjscL4p9OJXeXSlZzVKiCv4hyTKZQZzl69MK9dfvnt8+efffqu2cvd61qOI9Khw6vODz+Kd5nHMHMNWPD8aScA3LTYkLd5VSDR70FTMEvzhxHO1GEyGrBPmoCaleQOg2M18ADzsQ0eDJjAMocMJNw40Xo30BtrgyPRWFtWrUjbpTaIjgUyq8c6svdx89f6rB/9fjHlzLj7kg6eRJKFSktv12HymzTh1wlk16U+myE7ECVNn7uOJZTwMS4nXBSEpfr0J37k4mqbUgRnJRyXmZrv8XyLAQzZC3ZSN5jD7gsSNAdYWSF9XI86t0yGhRsjFtlttW/U2bbZXarzG6XmSvxhXvmfm/d3S+jYcM7orfeG8Ak2NvRfzfcvxv445X9sen+2BTJxyOOzRRSzBR3UQCMLApT+QmmPiKM507u3H327OlLiyc+HZ+c0jott8Bidl5Gkr7fppPxSI+J39CXeCCsKE9VJRFt8uQn7qXdtGdyJHUrsOtUE+TRZVPQrG0gLUx6VVpGCq87FMgykdkkBkjQKYN+zvqdLJeM28khd20Kf/ADtV6yHubdkolDtt/sn/4hOxVCQ1u+myNgBU5dW/lVmS0Lv2F0TvVAql5mMRNTv0/gleU/CZl0g7nUrHu1bkrUbecZEyRkx5YUSGyzMyeYF+WqA5zl7uFvSrz6UTpGMB3ca/nlGFQcK154BCwBjYUGB7RpEF2INiWSs0zJGdsoiZgkjxv9gFLzKTxa52VEoYK7AaLNpxF10C+rmkS2UgeHyfhszMQ+hGqC3vfoZNQg25oLJ+kYcIem4Be5qtW0cTS0XXrczkMBGRK3FffSeaXwpNK5fvZyljYsHhUH9AA+OdHUO653cg0tRz7LvKKgovfaNtLAuw8PSD4XElfrsdntsKEr+gxmxKbajB64rqHTVGc2vmw5Bo9Pw6NTeuQCx4pc2qVshEaoG/wSeEwKf5FFBTb7m3GEMveT++ALeCaMXYU5a7QHp4Wr5CeIY1LRZ/5mTbwuuAHh8D7HKhy7I1Nyu5aK4Rh5FolpfoeF0IyAsN4NnsvAy1tHVXyP4ujJdurHWU/g2xOFPLC4K+Y+uS1tw3KPh+tuKF577KdBm9FWbmtL+Ac1l3HFLHgbmKnyKv5UWpRBusn5iU6PiSb6bXY1kAlC7JacXrfKrmxpTMGOyAOOo3cvUONA37i5X9rHXFFZb+6JUMEyruA0rgDOFT/ZJigzpTuDTOsCkRnIXrqJGw7kRl+68XEf8HkgOTK2gMWu2Qn81SF5CEA/cGuQAXOUgAYTwcIEIQ03rie5JDApHZL7S85JoWQEF0l0gZBe6rfQWp+N4DQET/+aGmu6il0seBsJB+GBi40cYYnqTBzF1AEq+63xh2WXy5m2Dw8uTxCVRNDFHduon1kKxkr616b0o4AR0gIhBDeW/0SCdujupsPlR4L5dIYQiJdVl6//afP5N0pdRFTfaYrXcegFPF1BVs4UFyTxlB9YJIfnwlfAI/idUHKv4PDYkciScxFTRN3gh1dDv/2z5AwktNJKL4O8SpebKDU4ob0jghx1Iq8WfZwSd8bdX1cgRO7fn6KWtcSPUmIZSoj01BDBiTGRgGrYnH2JCb0uoAW7so5v9lhs8Hk+YI93W+bcdcRekQUS2E05U2RHyui6Tg/ZIKo+Mpbi92je3K1YtocL3azW8j9ns3qIjU+LUZFLSj7u+7vqIzcde5WeOzeLiHuX+rju+vdPkg2Of+vOOOf7n/DRMny0XP0RNss53/+ULDEgSqkuBpMoIHcFLI/kTKaRZKJ3EAx/I8O4NGLLOFBxhbWtUBq0ZkqABG+mvZ4vXRdc83+Z9YwO18evpQL+eZy/MHJFLPyJDMVP7rUnMZ8hLmobxCcqrflk+nNME9y4Gt/fD9/flO9lSLqtp4oAmac94d3lBLwCUxQRHiMWA5251qsf8WrZwtwkrntTAeBl1vFUedJwGEUEbHCcpRfKrIpebcJM4IsPRsuZF3RU++JYFMdo8J7nvSCE812X39TD+ahqmHn/UrdOkwyX2VXCt5bZBmiwEliwRo3XIMXK4jSocdHkyMJwB0kVTXgBnZH8snEBwLIC7kym6vo6jGfPT5u9cr+XRcpSrd4ODc8UHGnhTAQMmRtt0FarZPka1mwtpHYAvzKZiJyH/TUzlkb4m0KYl4yfBbZnzTQR80W8M8MNk9v2KEQlK32IRelItyX6Bb42mUpSOFNWp64IymG2xvEFxeBPkUQNybxqeeBwnpKNqgiT7jkHRURzGVZvZcp5zTQRNn0pHONOJk7ifnMKCAqKu825zkLx5oxeL+113I7wNDurtl5wWPct8A/5bEPq43b0B7ddl23BRmVLq0w/jGpbJidGGrufMexUC+uvXPu/k4qtZLG0OzsN8bUIcB6dlEfu5UEgzckJke3siZ1RnDp7fAXwqOxsNL1Q2LXsJYKtxLbceTR6/N0r7H6ski3v3tsm6WDSX+puyhU9Bume/A7byRS9K9ppNH/rn4NMmx7GLihjbv+9HqrZQF96Bitfy0+f5C9gg1mTwwajuJjFytKgSI+5xXzvytGrpfu/K8cgLbf2kxPtb48iOvotsiSX0tiguWlyhLNGHelzOjTmVYfaXfBIBE+rnz2Mer+YzVTT64Z1IJ6NC1IyoCYmZgDUoX1RBeTQNNKSZNSRDjEAdLE+u7NHbuvkBxWywX7EJfdncTdQF45hKAQCE6wDwbm7eYXgbZm9KRKKcbUJ1Q8ryXNfFH1mznr+sb1fNK6sMjmfuF23WtW41j+xmrhfy9X92rIKb1mFy1aFpzY/y9Ud+4R6gryXu5E6+shIcpR0pOsuSdcSz5fJ86CGChg3JhkEnVlF+uvVZr4B4UFiBZpvQ7tpqpjcLeK6W4JCNCz5Ej+X8jNoq4gqE9M+M3EgD+psuoYDLYfAXdqUDO1KN5NJ88iLgiamDBHI1lmIG1utnHjBpiPkNpBp6Q+cNqeKio73Yv+eeki2aK9owoQdn8FGaugV/8WTHXT6iDRiEpRINeauRtwcBSszVw4CRzJxLM2kVN7FjGhQFB/NwMqQEokWSMG988pJ/N7IJpK5rhlcNuqUrFWS3S1I8K2FDAYMbxMbit+o1R2bl8AilZ5mu4Erlwi8nbPqNehwdjQmHAp1Yo1t4FgyHgl6IjA2vaHHVhvrTrjEv2T2KRHsIjSLHTiMPYjBMO0dRLl1++HxhbckM9w//UOPyNMno4Qmb/BRjJth+4vr03CtAB8oeVGuKyVlvJcMyyFJcxPJmS+qVh4mmYmhLzGepgXiqfDSpRIeDKYroS+r7JU6y25Tryqpp6qMDykBvO1oXvOtUAPfCrV//I571P3C3/R/4yt7o9VyAH8Wb8OAq3qBhGoC/tzRtcpR0XkdzwCdErxlv7T9jPw5R5Y/xja+mAFh/FQfErecXlv3XXWebOumdBZZnMWi1za2pT4Ggc8RF0NlfVQvdypqbRw8x7yMz87EEkWOSQGDTeMtXkddUphoch5aspW5ERK6EETkRMCTSQ4c306/729hENUErwxoqJXGmQX+X164xkbzRJ3YqyV9BR4jvupbyywoBPNg5sgfshbSOkYjK4M8qVC0CEJVKMp1RQsxLWZxlIfWQ3ttmGVLD6oBI3SKg2XbOFHxtUDA32sfSS6SclUb3I+ZGqpzzRsknG6TznXth4TMnVWaAZu+y6iuKL1bezOp2Y0brGdIe2BwRFY/deu0uTryFk3d3D/dTqh2uYi6en7q7u0yjo3+osw2+9t3vdLUPeRwi5T3PEO7Y0aujoHQzQYGbUWrNMa3n4lVS0yJZ6fp9x/48JRGSYZJ41Eg65hCcUsXcZNem8SaOb/qFVnzdfRy2PhepVRWcNpRgReq/dsENMM2446l8vH2hzRbnaVqFjeQtq4L8Bq+Mjc7m6zoALNR79kLy7pruShgPL2trqPVeQehxNEhctoBNmAl+oODvuNXleX9MoNqC7VtYIqk0JKFllGh01DoNK1pMzDPjc+jN6mc8DY7XFxl78SonIszWqVSc6RtgYIFj11Ln6FvYp5e2tMlnnJ/uL1/fGmPZb2jD45P7ZWupP9M+7S46vNcvjybicD4eGp3s6A6ogQoAVFFjBZs8v+fXWJLF6EmiL4vsW/gtzk/Oajyrbvl9nbpztdm/84dZJS3Yi/ckuahGgTQNl9hJo9PrRkIEjBqn3aVLKFJab/CKh1f8v0ZVhpWjeNLKek7XTt29/XId9vin6Qq+Of8oIdoq39nGH3hWwDEl/RW/opP9TGadR8Wjcen+tgyOHry11DPoremZj1uGKC9pVb99Ct1WOnSgBvWG26BGpea3FeikZEcsTV3ADQD+UPVtnJfntKeLBfF/GIyGrhb0F1e4ID7ReycJG29VPqfCxcLQ4PcArdiry+myQuejqF/kHwl9TNrG8S9DKy8eyVX56jhqVSJXy1dNEwPprpj71OnQkVkDrSKqki1s3Xriy/EOwDVnCCjm2VcuxRUEdz9Txy7QdcjJHtGFd6nC1N2hOCn02oyi8z29MFCOGjOuG2fgP7S28wRmCOKajcSvtMkfoxEXQOilQWX1tLrxMsvitrTZkPsnswJHELFac6SB2pasnFtMIGTUUU/AY8bKH6c0Wo/f/y9uP3JzEObiZEOBOkUA9vub5c+I/MRfn/BnaZvGFHJvwGZjimJvt3q30q+pRsiv90CILB+6/7Gt5iKQbbnrvRb+5bezP2+V7pF3HcFRJwQ3rfVy1uNXvI3W7oV9fLWil5+0ejlF76XX0S9/CLqpRvH7biXcKbc3t63gDAckG/c+vilsejhccy+eV/f0ruTesuIV+Bi9eifoaeZG1YMJbLrcerl1OVrD71zooVGEltgdnzs7aeSaFA+UxaS2eh1e+Lgg6XW5KyaqUVcBrGlh4IgOV6MaH33OUf82ZBg6lJiLtee+W5Qgcm937fsix17PQ6mFtZXTiBPyWhM7edltVwRsKxxtmlmAglWiBmavbFbOkn3tRMngiolGe0O/7HsvMNm6DPVbvI1GOLx9JeYvlzO+8s+J7OU7L/FSsuM6Fh2SHD4KeywRIBwPZLKGRvSFUPNcim4wluQ4oHW8MCKZANpp9R1jArY/eiKhPDoPZbet4jRQeZHqd1F6zwOvfQOe5hSoDLIhycjsYrVckWIJGnpE4UmjYTkrkzD6I7e9jDjz83Nz93PzwvVUdmem9QzTdnqN7jE2mvCR3XMtJSP2H9TD0s10gCM6Syin5dz6DjVBJgQzJBhkHGIBC6JwoHPJGjXsiLidZ8JHfOb+e+O1ovf1Z/99c2gHD0Dh5xkKDzb29qnrNCoZRgIy9dUm4VrdzYNijFD7TJ6Y0dMU5VTlUbLqI7BwA5CVstSTrfU6q++4M0+cTMD/wsY4t2bY6y1XEPQt1n2AaMyvHpztRYshGXoqV8tEsjoU0l1Og+IBLxXD0dPoT6pvZMwnEtA6hiTovRsBZn0yl/w6cHO2+aqdEJEMcE5Q27zGvk6SV+DUcjSWYvyAVzCVM3DXrEAgxuuVLdca46qjc4XpL96t1MJQKeWRAcwm45M6RoNiMpRc9Sbet7NET9UDu3JmkVnrPK9iNS+VNHkK5wkhpEx3XKTassn4zejaRHYG0n6JS445qKu19CFOC4GK7c8UYhpN6WGNfF0VL0haL/R9ZDyO/gWgznFga1zt4AA5pyWXJEyk3ktOeqU1MstuJPtwcJCSG/OjCOkdGDdpCfqpufKU9WD8eep59272NR1KPCMCOp10i0oyRWIIx/ctwfr/NcbkA8R8Kull/7BfXsgpaUxoegv++dOxKgDTBl+roCXGJMky4QEg3QTdEJi5GfAoDWsAYElIU7aOKDLQigaLZ66b792C5Mz+rxDtqJIuFny/0XyYCxkueXqg2zsEd/C/dp3m4kH6rNwqbFw3AKHurI7EVSGGytawA9JA2MpUN0U6SGV5nwQfXsYrM4gALr6j99/Vv89jkfXhMZuaHJAuMup0ZZVaWSXjfZ3008evZQ6aO30Sq/oOYyNWCT2G0f4JTzEb8rkFDqoqLfuZOYmf88kcvPu5iZD/9g6tUf8eTJaAAL5wlEwTmnuwSEOZkfLgp8/qc7GxGBbCIr+w8n4BCepB2dzx0WEN48UBwAvhbz0Prz3vujaept9oCKis8O2NsOrBaJZJeobjnd6xA7INKryarwfkEkBsJJT/7RB5ZGZP6lt2qDWiI+IqQLdyL53p49yU8IOL+vtOujmItfmsUpwxEeb7VRYalGOCxWaCSntvRqponzIDQY48CBvM0kofTVIGueWWjaeYTtdMjwHU3vmxEt3uxBKgDYoJXvr0SI4wjpDTHVcTSlzPpB/0PN34dBE934uM2UBA45Stm0ykujT3V1SthDNnG4b9owv+jJG+fvK/oj1VR9T9p0/KIBymDWtrw8928PIHsdaMjbRsZ+4f+PM4Xqhy05XDoQ3vLJtkoFTHEXEQ0z4k3EqvT99+Ojx01cvnj37joqtbWarUloxmxzVTy0LKoegEVfYJttudWRYXAK90aLq7DX5N7fxmh9stT/gkdOZeI5FrBMuwDw9Se3qmXCmTGNvZpdqStmC/IMEgtHhWzmHxOlnze2FtZZLj/IsUWpP9bwzbCCOQCxnKrAoK5Njk2osG5la2ciCwOd+KjdgoXO6SjnnQtLEm9h9qZ5mjkfDhHgeLdTgfp2B4aIYItuBoaytELZYI8TlJF8cFrYQxytjzH1g+UgVTxY9iQATIopZntYycWcaQZUwV1uL41BLAY9DoQkMenSHqtQG1s+ecQSIL/f+UrS100QF+Blp0RY1iFUKSKIag/HCbGUIoeS6oydDIpWa63AZ1RrXZUGL7bIzk3kgn0yQN872H3Vh/UxuWz+n+Gg8D+tgYMzoG/QP0TFV893U1IBBWqRA0Kd0AtVqQ3lDeaQzdm4llUu4WjGTbPa37niSqlO549HX1KAmJ2VHokEtZzDO+15wQVImNTrevL6Cr1OLi91oFb4qYk24cbIh2CqpLGVyXW2h3LIokgFQqEy7r/o3C+BaQFXSdDXLo/6ts2ChQW0eUBQ2A1bkSB35nuR7ZW2k8saXS/8lBPOkt9OndpO3JlxWyBHLo4vDUZ7XF2duIAJncIHAO4w0HebeDHfNrL/cx/9GEbKlz1Frzcoe7mhVdtLqVvMZBQa5Mx9kmxAWWh1xss8MKBog+TN/P6S/MRWhX6uvYx+80ECBDUQsP7BraSN+WunToolwsifRZPuZDkM7e3TlH/gbbZ1Hp73VCQbSuMysrO3NRqTm+WzBCY96mFypwNzTjimUnoHy4bmr90helui/G2n8E3/hUfxz3f/030ZFfAXxM19uvaMcn+0XOqgsE5ZxD+u+n/CNNOoKHUl2V8QtIhR0YdxhtSBPKKs98Eci1xLx66I0HaDu4M5CMSsYUhYTQWzl7AOOGivkausrc32fLu73kNBB36jZmTA6fKZhNALnVzc/RJ19YVb00Ud+GRWzc+i3vZMRrHBW2Z8b1tMD18UqoKMo/8mQMHdzXQ2NEbUny6GRGQryulD82nCKEy8Q7zH4Ya/gYRwpRdn048CLPpDSqWVBF38Ub0eXL4suwW1rc+vu1lfdpmfv49SEQQqIRgmSUYr3OOx2BdQp4KOu/rcjLpsC9MHoZDx9XiE0cqhCg3Q5xPjo+DeiR1rWzd65z8HbRKoYrI6gUqaZ8fwf9A24ERl7g8/iCm8B+Evwf5ep/7S6SDR8l/Vpsix/lqNJt3JigUTbyHCXh27ZEhxWEz7+qdQueB3Bv5hRkdcTNCyVNxuaRqhv+2ar3vKO2db3SIOTuLoU3TicmCAibwpI8qaAJN/p3xEE6DuPvtp88sQ92e7fdfQ5EV5t76/SUdnkp/ooH7sbqaVUpRhccFaLsp0jS2/gRDHiqyzN0BbsFrqKbW3S1van65CSLiwYmNfQQIixRbW5Sdhd7CUGT44F3Tm2yy5novXsbtGMy4ti8taJr+Nro4982eUw5ApuF4kus+1WY3qsz70OtU0yMM3qo0JUCNz5W7dKjGGdKD7u/+5GFDPvpJ7BFcarOgS1mLAxMVSWwXlp+Ejm88AGb/uDEbNYU9mOTd/06E/5tWZ4zo53qY9Xyc3lFQsV1B7diVfAvUMFy/BOoUMCd+9jF2IVKbl3DV24fiHubZVb97bL7Tt33Gps3W15OHmrifHkH9ZN3u3STd5qex61al6ljKfHQXfD1zk7bX9Ek6sV1R9UeKUR/OYOm0ZgtXxTLTjRvX7O7x1vN3rjjn2qgLWkrJ4rAbTeIxwIt/++InfBATWYU9ak3AdC2yF7il5jmb78UV8C6uid6mg9wwTS9XBh3fovFUZ1ThuVzYO2JuciP+9fteKWlvJi2YobgpXTxySTrE6hxcx720dxJAE+vR8zGzZz7p4NR30FY0EbYgSB1bqU3Rbaxb0MAFf4xwl3pn9oL5Ad6FgBmNGoTxmkqYGQRsZMYGrKAcwEv7jPAsgMhLGQsak5cd8KzbR7DhpmErbCNzUeduxkHyeFQvG+fgSLvDfbqvtAG3igNBeGnxT9jH7zTVWP6L2fTXNyN5ajhD8womCU67I7oNT+0IdvGY9GcwPDr3B3RRaITXV7tSCsRAiQ6EDwOg0ZQDrzSVzZikF8+Jx0jC9KiOUpfvsMIRItvtITE8xn7hLsPGTtryIrDb5qQGcE/05iTGnHbn5qulmJkSs+dp249/4qCVITvuqrShbLdD+0b8Nf6MXsMjBYuo+fHfzefdQH0GEOn3s7QA/43SDLtYAmV2ZmjLwoIBzDRv6uuV/q6uwc+tWXF+cIS25sGgCLZhGuN5FG4dHhrgDXAlYeT1QFtKVmqoWZsqEJxHvYorwN7Iivjoga0tQZVVi2RJ4wZ5Od0GgYsn+0F4Wwt8ph9qEa2ErkFHqQ7KiGN0/Xov9ddfXqpDp/BfcV0c1KeQHXkL/p6ui2l/YwiRrkzrfJSedkTwz6qGfh6sAU5TpHfrYkxfr9Hbb9YMX7QiGDbIYwzFLdBdhuSP7m21q3usC7hvrabTXed7Vl/EMwMcpLzlhpc+B9yOKM4Uy189Jdiq3ERrPLKL5AFoPxc8Q40LVyhfov8YjLM5oepa8eT4/4wiP0RO/MBP2q9t8bOv3jZj2+7IgVNgZrnw0yMd0qbA4VVfYk1EztoZRIWqQiUYa30d1PG/5+6SfOPuHY8fwCyaeir19cTKMBNl7IaPZT5z7GlNYeS7BXE+BpDIwypBoSc5r12/EHSD+Ur/HftSKg9kzHZxdn/MK95WlcKwRiLamO7lKTy2qZOF4JUKJevdIEsvvUNOGYv1HcRgTCy5B7OpRGqK75rn3+qnbS9SEc7CpfLPiZh7yWyGA5m1ycTTVT3BGgfcR/8QQploJrfYxjwtYliVfocYtFsBl4tMRSxIgfWJlmV8UP81Ey4wEDpCe/hx/4/AlWIPqKK5L4RfBJuu350Ss/DmygmSTFjQp9j0dRocBmydZnvfFZYB04BfKXI6TS9AMbC2Oftb8D70obBzQy5rVuh9b1xGnZXeN19nvk+Tsy9LzSAv3+8799cvPJ+5//mFVXcFum4bnKLNWT2hQ3yBedT2aLfvZEjZzB3BxtXwOecY0ihdSkFt4xbFZuJEtZhFMFlzYL+5M9JP6vrVMTwHA5K2b1HS1UPI96UXaek1BBUq/rTBGOCxPu6MDtUMNlj/NLP4UAvSsIHkxfRTM+/BjReR8fqpZbsXn3szVRZQmu8Ro1NWvmiTuy3DQ+FITSrZzIShY4QPcetDzA5yON85au5jTlffj6OAh+R10kFrNhmuAG5aau239+H9JMKymW2sTY951M0geL3sUUhHNufZgG3mYdr9VTkSlJcWptlQnEYj/Ih/KQsCX5WcrqDEIX3f2+Ba7T9/8zWf93RauzYE/2+FcemixDi/yQyd7YT075i1AumslAfMZTJA1yoy2ijGFcJWlSjBhNSshRKxV4kGmXhBaU8BQL+7jRnSda5hrqxYhMqVFaKLWpcGYalX6vxa6hdkVp5hup2afgY1Y9XQpsq0bVj0K5zSJlKHR+jEemu7Dwi7olERKiCD6ISNDywndvoqwmKjM2xH2pQH4B4SfU8RdFyOcpHjBWiSohR0E9qTkdVcxpCh9LrxV6Dgzze3cdF/n5tqokNPpNwLAh7keEHxpiCOrH4xNFRIbLlLnvj0OYNDNJCZ05Gi8EcARuIaBMA00mKmfvmArn/ZCObzLGbTiLv1YCaE7nou2icFEr2DggHISco4uwbuLOUU5AmshAjy2KJ8JBdL8ORSs+V6WLj6VA8EykQnvy7PtdQrzvETlkgj31MPmjzL4ZTd6MFo75dsyiY6o36tF8fCwJKE/tFYp/0/mjzLSa5rfMD42iu/gj+350mb2Yncm2Xms8czchn8CK4ivgDKGbAaNEnuB+l+R7Hq4gXAUD7+UjLk28Tn1+hGu3WIhVWsK6wJWs9WayUEnZRtHq+pkPsY2ErRQXKLuylBOscZXZPhhmzx/+8vGrl7s/Pn3MLaRX4ugyVFoDfNggShaTpcU5cGfLjtbdKvsAx9g3iuir2ZX7I3/y7S9//eLxqyffPn76NUIPEHEoLlVMYGungjt+RL61OlI8vvlMwRSs6V9PMSuVHreTeXVOAPzz+cz97U4KoIHF63+KtH5DdR4lWeNtZXn6BLBxATAHETb5Fc+v6zKn5dXzZ9+6Dewx6Gfn1eF4sZQZ3oTebSvmpWVG9FTRwTW3w7cYL7A34E4YkEKri4Xbfa4HchU2VznEdgRyAUoCQ5WGiAonxhNQlBYx508qS37KcZVBv/zm4YvHXzPo0rrNcPir3fjXMvmFptz9boccgx9k90qpYJe/tkDpp69dqWBnDIl+CNbm7goN3Ir6g46o4y1CHlVFVNX1+GSav31XJr0uXaeVJ/k8hEzctjsPCW6ugFPHSBesv8JfwGhdcjPUsgWYVOFCUiqczMdHxHqaEfKWW8md6p9G85nWe/Ud4jR1mr6rrvzfgHn1MxYVWkaFlnEhy/CEcDE3QWh6IDhZJQWAH2RwEgKO/j7l3h0gc4bc1r/kJ71fHI1Gt4kAI49/63qLx/eq6s7BFlBd3IOvJVhUatvub9/Zj2dK4i57muQyjsWE2KKRl+Io/PLxd8QMg5Aq4VmHVY05tBjJFMAfrybmaxuj27n/0JiOkuGt0tB3zFKNrfP59sG2H9WLEBX71cOnT0OO1v5i9nR2OZp/VQGR1op/jx3Egu6sNIvUuGh7GEWPP57Jacdkb98pdUCceqvuK30UspmG2WuJemY3c3/64flF+4pOrr1fbB1v37t1K64nmaqD0ek4koNAdsBCTEcXi3k1keueEDC5kSqhTkL/vR8pkpK5y+lJhRRTPrLYMRvV8UKYEKR5sEg6iUOnp+hiThQb7R0r8Rvyro1lF491mknlmxPN7/zsyr73nfD13YsepoXD5CAHpsZsSorJoYpmjuueXNQq8nnvdXH+k9TTdhUPE0lSN7EIrKKxcp9Pj+BnXHG2XYl5owrPUKMHSLw4OydwlnRioHls7bdQw579jEomZ9se6hG9ZedTJ/mlr733i4Pbd4437/pNE7R70rB88Mg/7f1idHS8eXy7V+rQrKD8elRJuPwd+/IHLeS+u1d9cfvw2l2uCg3BPGrI/JUlyBXndBEDFGdgzQti9Rotqx4KCbzN0ZFcWtxOcetNFVsiwZRcMGafaUjfhV+yN3QLlPEH6dQWIZFgeWDvaDFbFha75ZclEZUwZZ8ffvH5vdu95D3lPLy8W906ulelLymv4WV1++6o2kxfmgz7AWoX0XCdl/nF1FSEhADGYmCrDxSxelrN56ggY64XV1ptR8pnkEtSW2KUeIesrAe4oNYWPPQsXiFAZbB5fHM5E94HWT00HanyfrVmfYBUCZVwvHVfXJD0H90ZjY5Am617tk3d713XQXfP9+9k73T4yln/eQzDHbAm85jDiYY0PhvFgGMR4rLydRI87KO7I7FMZC+9OufI9TyIuMBc+/wCLwKbiFmzIAGpkwzxp7AebjBRMfzyBdOzrF07AA2QPeFtyBqBAc0ZNJZVAp44VFY9wIuQyIvaDimiLupRK8TZ06u5rbf7U1c2XPzj2i7Z0R38v5j8RGEic8mApL4B3LcjTbdtSlbpC1NFSyQ2laaMIiE35JWr9WwyPiKhp8M/gexkBxOfRRElpQtyfgypjceHXXg6okDNrCiCX6EB2WxUd70Y9t1fi0smwDL+qra2Ss9PgaBoq6w/XCvuh06auySOJ9WJp1Z3yqgvEYv4bEHEkNU3tBR74uqKmRInCb4BL6BVS6Gv7VkXExTdkxb+n97P4iWdUG5PstL7V1Ne8f6t+9kzWZhmfZpP3K0w4rVCR86rycgx77IvZEIlMLxiUD0jn6bhirGElsZfrLhHD7443D489CT62fSpkA/xtLq1Vd7eKu9saRByKPW1q9UKwaXL/k8c7ZonEtslCsX3t6L3YiM7oAWi5XKUEHoy01N45FVTOSi4aR0oOaGhc4rXNDWQzcQMkHbTeB5oRAiYg5W+8KH+riJ33ijXl15wnUX4psELFVpj4lV6wQb6liDU4FfMNcW/ffGYTwoqy5eBu+rFXAx+ezYnvPRiRHjkRAXAZfKu6OucQQsHB0F5iKHySawKMuWHXHzzUVu3MCQ9FEe0E03udBagqDoEctFCQI3WW+g0eLkaP66ZPhPyqIw2wYqDpOjCv0wg4deepRd1esyUf9T0J5ObqPsb5gB7LR9yM/fsYmdpudRFLR13wMgdsU2N2vlXJHX4pbSqt986AJ2aLdHpVIkCplytr6sm49dqWR5PX3vELVOn9gGAcTJao4HN7MQCQebXfjGTAD5BTQU4xOgyKGvVUzF/DpPeyQh9C3RxoTjEM7nMWkx4iqAVtHw7xJSCDPyWs2CUTdUyt0a3D+/c63ltDvWdRIyKKrGzAD/Z8JCSlHe6FH1oFY8WBulaEEizNWg7NTqWKkk7NKTOqoGDZcHAcQwIZaEKdalaldQ6LzrLUepoKYQ45VykYGrfhkgnI5o/d4SeX/kUsWCzqLSjBlGiJUNac69EPxk1oVy0ZdszvmjHvmkaBTERz45zzldqDVQAjTbzKmVbcW1KSxO1JS1YqEjo1o0dDTAKj3ayh/N5teyPa/4bCheEHdUf4gT0hnaZN8hOuCtYoVYLniQGobTXKAg2+/xKPtLwdZ1/e2YEzewbgyyn9UBjmbDsYlJyf8AhDDYhAQ+VV3/4A6La9rbpKmZ8btAuStFod9BwxyiN8wokoBfwScQdXe2O79KD9PTZs18pOIdbuBWyBY9D2TwfcNvno28eP9z9cB26ZcqOsycmMMNwE62MOwmJqijCOjwdWU5hHMuJoNuA9zHJYyBaI7U/aVJQgiBiQ5OFPXVzNlfyYmpp4uZUyj4dukvzKEqfRJ2hcNyAYD++kPxoJAQwlPQzOqoIGlTUc4OfO3CMGNHHeGH6m4hXprZv6E5nF+Da4TbBzlpQvcduYs9MMgHus2MjDYAI3UzIpdeMET6rEcrb34yDKzahnoQf3U15Vv/L+SKP8TenhuxrQB46iK7qWd+90vzyrBeujCz2XygbOCuudfd41xLDZqouxRRRw9UlpqUMy4eBBtBms/lrTeXKTQIqdDSuzycV/NchQZkblquO7fmMkHZHYtbNpmK51dgRVQpofwysz9ywVCAgxwnMfxg8KFi24Jlc018fHzcd9GJvUHWamDMicN5HX/FX7KLKp0UD7yZCypMKUOjITUwjCklHspPN+ybaSIZ3UiqkvEYc9QVcA9++s6CKiKDTtOPmod2h0JJPCsJm3Ls9fLU/pBM7Iqp34Nn+WhzXE0d26R3o6Ea2JU7t+mhrHzZ6PFjqyCOXdrXnLwnrAlBORG3vZEviz7xeX7ewrWgakK/kNfxOL5D+7bUjr99X3w8xcumse72MLxht4o04YG2q12r+Bpi6xK053pQ6fCSKNWkr4Wp0tcfxIbpyViAl66dVvUsnHnGAibLoBq8WfgZnNN1W3mXsBuJOdUehw/6H3IzuSbx+6Vu/s4ZpqnI5opFYBkFtEPmB1cnhUJZ1vAgZCCXoDHFq7hSeuIvdn+roSB/O5nOyjPyOOz2mgwxL3ZWuuJlR7B/cUdVEknthSpiCSicwt9cFHfRtqPYU7qKRLzZrDrRxXH/NlZOGbuQf21REXuXJ9xLlJ/Ww1hz96dn09dgNrdxWl+5KLxwJCJMFCrgIk1R0fsd+RZVTbbOSiLolQSYHVtX7JVpSZyZRk4SSPdsO39klCvuaq5v2dV7BPV12mUe1G8LOI3/klZFkHDwnJSxFzU4fBWEOOLxTqBqP3MpLPQO5gt//q3/jh2B4e4K0V6r5jrPgYURaqbuoxELnGcVXq8mnvoYclxlYIt36CFAP2w8KqSkpXOQZXgqD4VhaJzUJ8zoVG17rmbsh6B3NWBFzUvLAmLtChmKOuhHtM1UCaFQ3Ia6Gch1R2tW0HD0JsRJvXGltnfQzIrBvYur6ZuhbwS93ab4xYutpHyYDrvCdtNZNir37kv70HdwKUEaslQ1W8JlVwxRLr3FRMMHS60LrbRPYt2p2qWE/nA44XvPkb2VgPBxZArM1Ai7yMnz/x38XKWTXstqxbcBFEp5Vs8erQHUk2j4aJmexSQKS5axWqc7xCktxbykFIMfEcexjqX2s6E7qTTMf47wEpSKTIKp0nUefoZ+IwlqA5WwdAXgYPLtoxqQ09rzGroFeuctddCe1ijIq+Fxyl1727TJ6qcUKDYnpz2GEIEH6SmpLkXdxFl6OFjm+Iwf8CNMLbwAQVHemOEFwy7POPLCn/QNflF1SrMLOdazkrnIrXgfHtMgFynLE4JryyKjC9RsJM3a89ksRTCO0lajWRRcwU8xY+W/he0JgKtG4C4oS9pbaI4L/QUrAysiP3urxYKhiPPetxv99cOdyLurYO4sikeN5K0d3+2EoSA+2tKFg83ogrmlrGgwmizorWPAEdtqjVE8JVGtdaaFZI53fsXcRTP/rMGtx/xPuTV2bLxg6Z+7hlZi7XCsTyTkjqyJzqJeH4IrLGG3MXKFBRscSBOm58f7yxbNfP3/ZfeW0IMCHf+ErC9tXUBYUILvryrr2jjLIbOlqnCpdngwyrVm4fk8gksguTWuQxDeuurN4al1lgDfRl2DdI7xvwwgvDWdbUL+TFmXn7UjT7bbiaoQZe6z7jBzZDU8j8rhufyGw1lJbGdge8EFw2vpAB966JHjJC6p5ROBLoS7QBsguInJhZQaKs5hrShkjDTwQZDRXI7dpCC6g2cIE07abv3jCin2Afor9bO2hnGWNGXhmBr51O+Rrgq0Y6ErYw8Ns7UXwvXkTUmpON4JPzpoqM0RRIhbLpUe7l2MnFY+nQUagJnU6ox+pJ2RqouH0tCUIc3IMzpni9lV7cqXivrqXJgxo12l6OYJS92Nk8FZcguQSRiXfzY70dqPU32uK2YcKNRaOLwIzme3JxEKteq+XWHnmbZNPeE8H4K81CEZLilGh9SiyFu172Vfhq+pFtbio05gJgbNqKP20BK3Ig6Q4AJ4DN83RhjOaHNBCVaJoEnonX8xgoKaMWk1wE1ZMoc6gh0wRwrynhGy/5Tin4ahG6CREFKRIxS+gSrPneu5RINVrKLE1jUNMQ/aunauexQMBriWeuIEtDSauh2mUbj5QLnXguBzPreKHMa1JT98Z1Q7TXxT710xtCEqUTLiVnCIa+mJHyUFQy1LVL9bYA2bKU2okSe5UG1tNg/NA7v9Sd31XaSq8o56vFWw8r3kOy0wQtoHCB4U7JKUbErwyHZGUSzmOZSvYcrgFSE46KKe6/aDTjqOgc4jQKXCrMP6S4xsQa154kxxd+rxYS3MWaF5y11mxDcjOpjVm3a0sP+IC42fhxeMfxLvIXxKkzgAZPVqVvkR89JK4iOMxU3mcL8ybxII1BoyS04H+4//BpWNPtrZuFWtGLo1xXJMSa+oPRckfVkGfT6XiHPlVdLV7EqTVNmYFw/IAubaIAnuDYy9zltuqYZnWM2bwhUZCO0xVA/PlabFYX3KqaujTEQCqJuenVYJvcgiZ++QAr4PCS7wWsBaHAmHYK+XHVvxjO/xgvbpssTeCHI9rGGhBbVWePgh+IiF45RgvRw//76868rqStvDAHI/OWrNbs/VnZMfyZOy7xhAF9k5pW9Z9LylapZJAP11iXsOVmXt+z6Kduti5buYtSdkybPh0+IQSmI70AIDgcHvjy6HJaJojgfM48dMfGJSwCBK44oTrfO2V7vNXfPuKS1RYcltDS65Ca0ds4JWq/xkvrRjBngVxCxoxNQpdrENVHed8BFTckVdxpvg8bjC/GjPyXRavzTLADiEZirGKosJbsXrB2aF+9aZ+5ahieMBTFMaDbZ3bZee+fIUH196RjuERfI6be39Xbfy0uXFvf/3mSYlPCzmhgI7KI3XlajXlK/fHK+gpX4m8c61i8tXRkyebr2bHr/i1a0m+FL1jT9xgXcPJYZAOdCy4EpNYV21+iZFjuti9a4tWH6jnEA1DgVJ7n2xQVdCs5s2Syy4yKt3cmkZ99RJdqq+NHn9upoFsQTd4t8avWEHHoKLQ4i4fioSjFaT6a1lao7TKTxkPXGYm13XtUBoOeApwDcdkZ0+e70NU6nw+yJqZmxBoMbWrVNT2ksBHEPCjnfhw47fNzRhe/92r9T+8Wv9rvsAFYvfcTvYpmyXJDgFMAjgxIdbuVd3bH8Y6Ugp4MncFi0q8LUeT132yewXRk9lSma16g7vmFa6+DxSZypYwLFHmgA3tGnRNvM4Qw1z/WQg5eHJ3w1mwcJfiOFFYN5TIccz3nhTcG+9fMyWI5WXPVEVc91XniwXwPzjpZebLiOpXiqgaWEoYK2DSlqTAbbThHurX8peef86EhkO7L/xUFB0qYRaVEr+b9phW/HeSGyxY4OajagJAFfI57l8xAd3AnxEKA34CHLv5aGvfU6SHU3ou6g3HxFnqaEu/adETBG0aCwRPYfGtprJKrk6eBAnwM0ID4HifHUi0cbnh21uwN8Dd4xCRM9cn97wIelGCsDfjcqMUS1PxDEtVoOhtmXCujK6M1Nnqd06HJp9ZzKqtYQqawcIA3vbF6PiiriZD79TUq+kA9aNUmkoTrIxO6Vgh9KjEYEuMrJH1B6ZlWRtXABvH/ZOJjzuhlE/HUQl8/QCVWImtkMXEFYSgAugZJPc6HXtSlxpL43xZguZAj97KLPCCJhSPhurN3XFkJfQ5SboAo3yP9nhGMQj+0UAp2t9vWnHHZ6MGm9rnKFfEDigQB/evbEyPhnFULSo88wZeXTDApfUzgaUccHMaQ3ZNIEISttt399tksim7yv215b3hPjSvXdcjY6hzETPjC3KPjZTSAuxQfglwMWJOQATp6Xi1yXtRNafSuUGmr0r9Y6tZZsuXgVPFnswIGgo7V7ctRV+ETPg/sQn/8Af2oqG6vFJQ0Lj32t13XaK+X2yozOrXo6M4ymDsleYQBwfdmkKNtkmTr6sRQxB6rjNitN2YNHYist+v4D6QY9L4jxudjEbbi9a7JtjO4o73zlf5FQOkr7aIQkZ1xhH9wJuEL6Srlfi5YBMRRkqtEpEmOLZXiDI+NicJp6i2Cs+Wnkl5GoxavB6IHyoKTQPeyxtEzH/H8h434hJWhFipLxQii+X2EO2xagO4e6lmvQxYOVkwOvpPgAxUHR2NjobRRRQSq8IB2LugWkCbdbw1SB8vJxPZCosb6AhjvKuL6WRU10FTFfwHyYal7rfxfxJJoXlW4/A+6KHHh6+RdGSXrHs0g2JUxc4xRYmPKGwTnOcIOFP/UkPnbZprGlK8q5VpNlew4mUmCVYbOG/4DsauVZB7Fm7wisUE+MwSVHuOxtgcfcKkAORf3E13tel/I1efG0mKVDkRPHTrfapbjU9vrFZtwDB6uMUvzUkhWGl41AYNkcUMbC8lzxeS2jVVxDQJU7gZxJ0L/dmNJJjdhtRikS48bYMEl47E2YcEmVQEOJ5+wAd6kMW/szSUuLZa9JgmdVgRq8Of6bRjGouo0+9qYIXYpjCd23KyggtZUf1b00/JqtozScJ4tWVtXky1OSin6eGZ+0hJbH6iW46M9OvaxTpF7lVzoreQ1hB2LGpeIyiBa+YeZW9lU/uOyk/ro/yySbVIWZ2LBy08q5f9AB8XbA/8LCidlTy5+5Z/0anhRoIbJQkR24ZAuRxStS63nEKIHFuGdGYQaN5QVDCg8LHkO5dCI6Iq8CuNC6vERbW+REZcM8eNRke1z5HGi45fV7oKzEVN9J5Bm+gOpXoNAMTtx/gqf5mVBqBAve2lZLI60ABD17gtNVKqxvfbmrH1wbEAt7Dj4KlpU58YdX2Q4VgALkn/8Pr7QVsPRad6EVIGWbMj53qx1rgO1zw192/WmrnCQtiyoucERxsqfcwQSQnKXIZ9bKCeBy5716XwlDslP5dg5DY36o7AGeX4yNbIwn3KwHwtwu1b+p8Dtkxwr8ywFGXm/Iq5zDSAWJ5bOJWhoLTyDXhAs6aoL70wgf9tZMtx3fUmndxGNfzYvu7V2tv9tKs+IKxUk5QaD2gYyJtfmtVIsRgKTxpkmN5akXsrxTufxiHMsZqLzada8fv9OR5G508SExvQIc/Ra+HmqpiTCjyZqGpDCIIwIL1E47IU1Vzcl+tdHsgbBMtw4f0PhYUK+0hGifrvZ51j5Qe2Wh17aVcipuMVCuF4nExt4MNV+f6m1aXxfFEao4/ZSGgSCYwa0H/+uSa2qJ3giO/l+VYUjrNQ4IuIjPV8H941T4SAEvpToWARH73l5X7p2vb6SjKfCBbUILqdEI9ofZWSNK1ZP8NOJtcScSNdPTuYXXXTjUdepNDcKBaKp5PkUflWNhxYmE9q+AcFMzGkg4hFelTNffsNqn5NNy5W0MyVPXhxMV01ams1NKccKhtIEzuIJUbuBvf1CBGLoB1bqqq4dw+xhKM50O1hlMp5u3/Im66MbFn1+fj1SIKixRXciXWFd5IRmVVyWnrtyPwEfyua169gU4aTkkIhuI1/OIqCcORaFF2TRbOIQ7aBioXclNRIJgL+Y465Vm/nH0bV6++qc3Em89ehzEvO5mLQ3hv6MR3H5G3qlsxnRhg1vW833Lnmlcgl2TkCBAXOnH7sfOqmK38Taco8vmZmY0Behlwn6U3wDXvrds0lMrgh7mpLc6fo73v3CuogB649qOTOgGfxZi9uZ5/OG4x4HKY6aWv2xJpt+bSGPfKj7JJS3J4u29qPoK8I1mn3HyUmYRSVV6/VU6FMbV+i5qta3gLuF02wQQ1I1eAp68GbYdLwwcXkNeeUlbdVIWXw29p2VUcHI85728OU15ZQ4t6d9AhJyfGcC1EH59Xkv8cWACYHpjZvtXBOyxDzZT7gGiv/AU9Wt8tM/vSm58Z/gsYBuUikJ6C5zSZQxnifmYYbbKfLq1iq6Xpj9dAFx4kO3e1KK+7cjM80zKOOsRKF0kitCGws6UXL5I5iSji+mHPkB6PlLAL4qoPnj4VSTgSKodalwLRa6l43wTJ/cGOxfZn4TEwlHlqsGKBaGVz+FwG6XIAfPIJgRThI3UYzQ3MRLWRsDFn4MVJ1XcXw04YIcQS69lA8AaKUyYqC0yWOyKqoC4ES3AijPMzvNLuT/Y1ApoCEM2/vxeKymh/ZVp7ClYiTVy9GHr+IeH+DyLGp1hhVEK9SpcbIA0sCTlmFl5Dc6Edn54tlvHRQ2HeJJqKJXimZdKirZJkcfRVbv1BY/tlO3ujjSPmkdGIpKfcSlHL2ibUgqZjUshFVc6QuZXv7fyEBRuyWfVcPyEaevzHg4LftaBipbKQ1xVbAbivn5lDlc/IliB/E9ZKtZ8AAeNcQTqRkyOfjXTGYQ7LDucJ1Kao9TYArx20nSyqFIQuLKU3pmuiPYNuCBXmuljZW0+S6mvJTItyMJHFHPlIfsURdJ57XdEz1ko4yCBYD5fU2j0DOqqmmRNGye8jYgD0pdnA0lvfcj55ezeEhfvUS8zK+iYQ+3Ri4s9xucvLmfJTvoYy7yt1lv039KR+yYj69d6dgHomKTvYyx26l+IcoUeXPTagaJYo69y83/EvIrbalrYdZOudYK38M/KmyCkr2GzWVocWNzFKt+/K+BuuEfrgVPqTPnHyILYDXtiljB12Na9jRXj6Qfy11cVDBesEWf7up0dUNyo3I00GPP3cv/P5pz+Uf1VVH1uRuYSyuRGQlQegOmbRoEZZ9r7RHbcDc/fooiqbWGsNOYQXuKPOjLRxcaO1Y+BSpZNz7+9gyjkGUMu7Bepoth1lJd3x54Mq7D+4MUXRjh69Z/7r+nfpz8L5gdF/k2XHbvDkyqi/mmvONUJnT6Do0LR28EAhBJhZBQdQkqF5NlLYLMFgabuCaY4rU5GTz/tqRCysXg3iZ3eGpPpdVlRL0AB+dF9mB21ZuM2nv8dI2O2d0PY+SDiFxtnt6U0qva5bgz+QnsuX0F7Pnjjka1+4iy7eQtotTrZUcjsYT1uEeah0bH1FH21fELP1Nblw5CGXFYyt6U5bL/S2rHkttM7ioh9Z+JKuzFgNqqR3TMVXULusNnxtPCtZnQmz7xWw2FKezYK1aEWMkI1Lj/CeY5TvUmvxoNfNgFKybzyi1tR2xky+jqGyGEXozfGQx1x3MJ30MRP46HfNW1LfeO2OVsX0ZjO3La4zttqfb5vXnEwrL6rmtUEEDzBsdsrxYIykUL+fV+blgiAv7qB+UPjcIna6UZQYjgZyGyPpUD7E+7ujBMumqN1eK6Ujz6jjutSY/LCCAcswlya072gLWRwYQzkdWhrVLWcM5H1AN71HS8cnlrAzAmEH7z+dTRd3HK7+lj2cQGBTGGruyn32vOnGC9QuBOdVAH9V+8U8n4dRmYUq9pE9mPA2X5lupYVdzYYBzWQCk4KzVZ5K4iG4yxyrCEYcFkaMjt4NfcGVr0duyYz31qFLMvekoBHoCckmZcm38YDI7fC0cnCj0zddg0WcuycNqrjhrI0XceqpAaxAkDXcO5HZre4P4UQoEX3l/hZAQoRm5y4FWSydDSMplqp9UiSh7ubSNkh7Ci4DqMiVRdbv+/MrRPyTzdVOlKf92sot864tCdMjxszs+ZZv7cbsISmbRsJEbw59iUQsLErqp4EqyLTdcNdtbRRrwp7lpRYXMGrknWCWLw08ST7xsgLtaupJKCPK8meaY5US9WjBXSCmAoSuaDLOynsnjWFeESfRzpMPwrl88chBo4WoTkwGGlcphEFBaQOXLwQ1+Fdz251UdnWMeQi0gW1tE9jWMdC0c2cmcvLTsG490MuYgdVYIIpZ4l5n7aOrjeUkfT+YscIQRa1PEX4Dpc9/tzQP/TN2Wz54lHjtKtGCVnl32p45DocSnv2UzrHNnubWFhVeWx7XE8sASiAruWEnPbgWHz7cIkdryCvVLkaNib0/1gxP+TwfNxM9ekNRlvOlt15HzNz/wZM5mnMKAQA6BoAKrP3h+H48t3qNaDqW/ntaL7D4TdNqAFTANHKc/x57jU2E58DZYlZCPMVoKN3QWvp9FReSZcIcqluEO0EnAO+XD8TiOiJunMpJOHUoNjZ0zwAc3d1E2Kg9jqAkueM9IyLbEd/BcQJFARRKuiyxvG27VcAYTMZgabaAYxEcsXjkQZaA6ubNQrlzKs2rJguFmSaZdt2gdXFO27FFrz7t9jWnaijULh6+bu0syR5tw5GiYuz0RghH9KdOPj1dMP7qA13txGe3Cl+iArg1eDxNR7XyigDl7isI0FtlhGZBBosOP86sTRVnHlb+S8pomPIZtUhin13Asx6ElBkhpcog0bGcUN8FALg0mGyblNZscHqh9S3+rg3YGSR55iXmor0AIMJdQp5B4DJPs5bgB8FvIO/Qv2J9LfIXHwyY0iPQwugbj688TFK5CevFInq3WveMJFm+ZYNFatp1HPNNTxWzkEryAXLDTWcwyrUEMXBv4pHEprVESc4B85VMP3cwmjFkT9aZyMj4wULn/VvD5mgYV3i3WgvlpOhPjgWoqPWc3cWMHNiWAL9lok9ejtCtOJcbBXjKKWyCqgamZa/LfFYxPYCiEl0CNsfWoeeFLXIncHfd3tLh3mnYv02RzDbN5fBdiOBiw+8Yt9OLwNL/pJu93L9f/8DvO0d7fFfuf/a74g3ty80QgH1l/Uw82YZYKz3F4DZzs9+QhDUWIMe0p25NE1JjyitXdMKdjhDM3mR4Wac3CqgZ3tIu6lyJzviyNTD+6o9PfJS6Ld0601b8WXyfb4zxwlIRijjbKhgHhAWj+h7M5kiS4XZeLEOUeFsxGJH4YIWmQ93yKfFnN0FD5xDhm7lHZCLw7Ey51ODuNxbWZN5J3ThJ4hxBAQT+OUn239CjgkGjDIQA5QopwdRE8kmlTkU9ibJ5VPndN6SHFaVMD3Lv1KaqSshXdlHQy+MqklIA/L07y5q1i/YQVgzxt46gFuUaWesK1KTM//10i/jXSheS3x2F4OHFyPcKlUVUvPPfpWHeQlhCZ4vUlHn5FmaN3cLFY9FrHUq6T0tZGLggMUjrdF8LePIP1lWKcuHPlPqlxA2KzreOKsE8D838z246PHEUZ8SJ6y17ijn9JKA4RhvqSqMVevRgdLvJa2oHUA0HKiVB35Z/NogjIhL56OhNRp09Pi4jvvhA0QHHQMERANAWo1dcj6wdTj8peClMprjX6kidA3tWjBSgr3EtYcR9nQqGepMTB6GQ8fe4uuVwfwIy4O+OwlkVogY/cHOr47JX0zb6Nm2ML74Ky28+A+OqQ2XRDzmWZdQa3+3cwa/6FFub79fDeBIT4+2iy+dLmVbkH1qYKGDd+rs+we5HjyY2WWTpQaF9ZLRPXq9CMB0XHcNXppOhczOD01h9PX8PgcKtor+pWPNdpb9w/m8y4EPfJHrJnoP3K2kG2jfHbmqO3fgz1bvAWMtmTguNAgAVeP0JD+Mr6s1XoUr0GqZDjFvNbJC/RtRJ51kCTLDT+WADOKUBFSiSWcK9Gk2OKzKTp2L1TUYKT6ru1n4zPAS8DbPtJFLgCBPN2MJzpvFD3bBqYL4EPMy+cYwCI0SdC7NMWRKJgE8LpSdLF2dmI+Dca+0eEaYUZayagsBA+40AafGKDxZJ3q/RQKaGP7/N/Pju18Hp8SAQg8+SVTPAkvh3AJL4U/JWuqkPa1WnRh7/CY3cHMS7z/R//Xc8a/lLbnW5sJIraD9fQ8NLCpGdwIIH87lZlrVqsUQRFckdenropcCMyEC7aeYwa8kG78wvZX+eQZCOtvA+HOkoD/azVnSye/orKmaNRhIfSXiNYM4I5Xb9wh7ZaFLQh2cW77UNPf5D9V+mIJVMjdfpgpAhafwjFsM+jp3kzIR/TTduJIHQmg3lVYqUAyjoGduXR+ASzQexBBRk1oCfZzv/479//6f91j3UXf1XxgSQJUPXubK6ZNM8008FUMJ4ETl+cmi6BIy662Mr4x4HIPY3souTw4CoDoETtnqoCYreXRpJZsnoLJpfTaEeaXTxOTB2THsWMOBofNY9Ug+zJiuGOTHkmqhm4VJooAPDON/fcpG28//lfv//5f3j/8//4/k//+z/+p3/89//4f77/0/+18f5P/8nN2vs//X/7Nx3HVOtJjBkQf4ikPY+KkuoAONQd+YipdTmjtM3WoyeTWbXIb+Z7vzvq768X51c3+6Or0aHg24ttWLbBjpCTEK7PD1x5AHWiQlgx737Bc3d+BcndnY8VghVz38xNtGpLS1Az9tzEvP/5T+9//vn9z//N+5//2/c//3fvf/5X73/+79///Pd+vnp9gvU9O86lysJrKnrxVL7/0//9/k//4f2f/uP7P/0/flY7vrXE4+jBlzsEet5ALCfr1Af4vRnCqeSTB1lvc2v71u07dz//4t76+7//N71gsJ/QvRUswD7zZ6OtmLXENEV6Hvzci59Rz6M929EGixXFBKFIvXaGxgpNvR6GH0eBAIC4S0Q5kPFL5S02E8WwJZxytfktKBvLPerbNMheGXCvDPnGamvSfryLrpZh5pW9/qt0G1cqSCSChSi16PKp7yV0DncwAEeu3LxIxWB+BnEh5g9Lywyyq5WCS0tX9mfMR3JubQLkCgCTlocvBWXc8Xy33Mfh8X3Zk5t9jGWT9mx39CBTLrpnLvQq9CAeHOcjUU8x5az6G1eHjF1aaRceNozeM0i7CcamYhjMg4UDGJiRC9vut1/96uWr75693DX7dnQZzYEKAHX7haRBkw4NsgQzIMUGCPW50d3djNClFiMYYmV8EKbUwIvnYuaNHL1umCsXfRcArQm3gfugAO5H8NxQ94Ivo2Y9bxSMyg2vGW9xwdZ9k2j2g9dD6rEwdOXu79ApRF5/xtd4vG4eFrHzC2usDmrA1t+Pv1BYqPU3TfeHgKj1VmuSqUGynQRHwBwTjWdilDVdgX+wxPN8c8BUHDRNtiAm+M0hknaF9Go+vUjsAUGr7tSzw4aDpuoYsb64G//epjXtTmTTvtp0mtz99rvHr7heQbvJTl3OLCucN0gnqQGhU5rbBIB7BPovYiQ7D4puYTcvwHJrblLrAjZp2KOYF79BLdz/h4gzjDexLyDAAS+5jxfJG38HbUd+YPH4I6ebH4Cq7SYy1eMUhaUMwfq/pJgp6rGFZNVl0IG6JENrrxNAiCWFD6CtyPuVaxQpkAQVDUA3RfB2ZnbGFpslMi5dUZXRYiO2WX9TZj/SWE3HFK+oauoxcObOxfyFJNarwXJiPA1NQJE4gzXwc2L3UBazlBR0lRzR/01BddzZlr9E93p/J+p1pFP5Ta7lHWPyY1yp3uuRRuaaojJWRLV43Uhr+OGu8Dqb6EZQSaItpkYBkaZpdatC5MLIn1pu1jhBqztEVZY/dTu7sNMZCbRwqHhENaaw4CQtIivHnijwVWMvKEGPZU9SXEaUHnzKm6I4INU1w7JQDRE4GE2Fb0yTq101obq5CWUGJLJQ1BwS7l9mT3lWl6mOdCGOcrFIl8RBS4qpH8C4wpm7zSfF4jeUkm3+RLmdj1CttjQ6aXLCYeAZUuHC05Gn7v5BPyNygkc/QMPEx1eFaAD/KoEj0YAROoAlW+FkJpLpZfZP/5AxeLjyIdVyMXs/isXsPE6sTWeDWtwpSESE6v6oEQNQzetWEWepJGsslfCSmPWYHhbVQhNL8mYyAX+yjGrWGKm8U10zPfJZ4yOfIbNob2jI/wEC2pK88+pJCD004hE097Q/JZ6EUm/tengyVkwPUYQt6vggmKTbpQbS2FwPpyOTqIkxvdsXUjWyqmgu7cK32n9ZnUtp6c4H/JtY8sQy8Xpl3dBA9+geYG373JSHSCsSPGiwU7X1h7acUPcnHniKz1mpvwSTo8+9MoYt9b33ptAgguD6ECpoATV+RLkGqvQANmJEB/Rpy/MaEqhfLTUjmXLGi5kpTo5VEylJ/Dxivca6K39ARMsIQyV/8u1uJxZ87P8FgdIIj3pVqktYWMJPcAQD9ZV4DHS19li3TcO3B/D7gMdVRLJoC0+gQ2T5dizClUsYy2O4a+lcFqykK/3fXP9v2TzgiiuyNvzg04TBvsjvQkNh44dm+q4A8R1RwxJNH97d9oVDGtNkOC8atW/xi1xH+QB6dmri1WPTPJ8u8i8KSnBufn6IHRi2NjGsDdfrjexFF9pfa9oc7fjESZPaMXCbO/ICUQD9TDz+4BEoZM1MZUJ0H4hXxKm5yXGeBo2FJk0SlBgNwzGiGdutb9zouA0dnRQ96Y0bXRHGrVY4L9bcA/FbFC3Nrs2thH/KYNbtK8ITs6kHjbjjMKik3m1dtLgHEg6409wlrOFWc6/BDPhNc70jJhxqnF3EsbLSoNkVE2IrmYHjOHbJdHxTSkd/HGQ2xov8c/XX3MXNGY8jnYUB6nDz8OMg87OxyxZve+VymBi/R8CRuU92XUXfuP9jjw2nRpGzMn8M3qrFxx45RmHAnciDsJEcD8c/oFoEy+puC1+Li6X79ilLYVO+86A1MdPxOIS7iqDYFDst6rIeLRaYloPR8YzyieF2a5Rqkz0Jcu0g1OkBo+Bsbl4Aglh/GYNDGPCzxAUQLaooO+/jsuMGztRLIeGmHAeVnwyiG6KI2QDyxZOlZWCTrJOliagYBeD0aCbaDR7I6tqsF/BYL6xE8wNAOarFgUFDP2vYZ/w9qbesanycYKUwj086EhNxKiZAs0+vtvQC64a/c6u/04T60ytRX3RAcpXevd9/3hHMIEctji/ylEx+dmEKhZgjqIsnTGkUP6QOx3F5kow+BI+qFoL6A43G8EjYXpIyLD8D7Zt6/1E8QUopgapp4MqswKCJmH1REAY61XP/S0zphrXCLeRXs7PzC9cOpYj8aHZ4gax+fdg1C1b0hKWDS62cvmMo/Wwbd7Dq8MucLsXSUlswhPg6VwvvyRMijsl6faTmJbktxjS/sVOJM/QSM5VoV0w/KJEisXbFh4awllhxtnTDSQ7jh9i1Zd1XpUKUcRKtkI96wwLU7hWFV8SJkv58QnxajMjdDEMPaTw5sutFbhSU65/0v4F2x1b4M4QjFup2yk/uZ+7GjKuMs8aiROPqSVNg4MNyJUCl33Qya7It3yUiojiKNUm02flXE1+15+U+E3KRhTTYSrB6deQMrEZdsuaIezda+SzSOVKbYJmA5Qo5SC4Q17nZPOgtAR5geWcsJJnIWhVTLDLbnxjCUVy/bkXkCywjvQWgKlfrbYy0EPsLNNQfbc2Mdk6VGNqOhcDYHTJGtrhL+Gnn0f1RhiZ5n/kmzzRAvfvG8m1KyKpCkTVIPe4XDSa/ltpfH7eSpfjmaolrRKoreK5AWQfFIApCA+nG/vLi/BwQIkUfvsQfvGu6d3Iph97fJJ03WakAnp4Vxq8yYeBYwUnci6uaAZaqEW6pg+Oiv6HsxumCOiZfOJbqislp86st/fszfPkjCgqVAeMmhow3LBhFozpu8nA+rs9JKiOJT5joaqGuRSFMVtERLV1ol49co0jDhcqNbVNmrEiNfqpI3JGI+NjzzR3xfJWjGlKk5jq3iB8tyaJu6dvDyRjANNyi4Yva3Rms8V08s9XhfFZTHb808CAfdwWHuZW+eKvc8CIPvJajY8tR7qNd8p62PPKolvtUZzzxkSHO+exYWdvasGA8SIrEMRrIR5wFzeN4rPGfNYntyhWqhgmEFhoB1VYUMyV9I4ZJoRqIUBrmn5Un1/eBl6lMIxkpLX/jcdnL7DcBtn1d60HwYCySPY2Um6G4E5FwPuSLtja10zGw4g7Ed+6vb5pbSzW0MQbGp5k4PI8VJbsA0ofcjR6GlHO4ChR4mJ3Mxxb9SDXUZSXPq+xoBp211+1blgtds1yxyiP7RxOBDC03HpbdoGObCThDhOQVVTG/mCZfv/C/3RserwbCg3uDwRWx1f4NNqznugo93znJXf4jkD2SNn7pvtcHuBDTVgCQyozvlKdjbj74CRKPxd18mLwyXRfLOWOuUxKGyyN/XsFJZFonfoeI8cNxVFCiIhkrPQOQkHpsvnYfRgyh/iWRaCQPUJBnvBn5xwROwW3oH/PNAqZlpOBtzOFmYw5/y6Gnc6i/MSsgQuaCEra0aOOVCHmBVVT2wRHfkuu2QLBL0f9OodEYJgySf8MQgKM0+Yu1npq3EMocBb8g6EVXNdgF+n/V2LweDc9IVJPG/KZREkSGVKZNyVpFt/a7aFALVa9Njk76UPR4mpTqdjrQMVs6u06veln53zS+KZKr5zqzU9tXPdYJ/vApN2byZbS1OizBDY/2E705z5ue7aXouD7qIm3Sn1ivWaw2Htbeeuh6QT1e4pSl8Y1C2xMotdyHo2hWBOhyxzD3LAwjN0CRlhSkoe5Seyu1B4s5I3/MBXJpudit5dPx0dEI+ECWKdAkjlVYvXHm5BW5jHh2ZudDBoHyREKeMT+M2bm1HZsKg+KreYcpR7H6brLl/BezMQ2mmtO+uaVCqvv0EKyAOArZYixr3k72cZit5TWor742JJmfrSgUmrnx/7P3bsttJFmW6Ht9RSS7exBIgRBJ3cFkypRKqVNTSilNZFZWDoctCxIAhREJsABQJFsts36ceZ156B+YM39w7LxP/0l9yfG19t7u2yMCpFSX7jlmJ2e6RER4+N237+vasXyW7tYxsHTdbBmIKqEgBUtkUtg5QwnPEvO44JxJSBWl7qnbTEV1Ku4/78KczXVZF5JpQrZgBJkI+6kn4m3cfrbpUl+voxPkh1YQiQRIKwWOT2aH1ckT5vjzRRQk2GVB5uPniLGoQNYGzcKZ+yPdBMJ+nM4kiq6kuKk4VdufN5C4A1YMJofDXTmgWKw5qPiqMbD6R62Di1i+14wu79CmvbTAkesh2CTOuY67lsEuL95VZ5JL+XyqMa1hQkdLZIumnfXQgW4vRtXpIuuevxA/DzdaQc7OPZBcOj+GouCcA8Ud08LDTdjN7hYUiSEyE4uPMVYAd6Z32UEo2q8pgRZ8/AxCTnR07zPHIHX1rjsB6eNPbr2KPFrbdXwdkUPv6dMc/l1fb/M8Qo6A9wfiU6Sd00fhc+mhFYmS9Wwx8vcrI+EFkMHm1wd+3eB3FUGtVLlZCXRDPd1YC/KeJnLGF/DXdANf6Yb1bRKsuoX1mwrSakF3a/GsnLhgssLKrT4DK4DKlGXO2epeK1J0C0R0Dg6dcTgrrzNL3LGQe8pEyppKI7wP7IXb6/agNdzQwc6z2AptR9jPW6bEO9KCJDIGtxa2xxcrUGoVbH8W+eHwdev0GlD0dQz6lkvTZaxuAPVL5fFBJgMJLH8Dj38eselqgq8lX1lYoo6YEyN7ytQY6YAshLD4Ijfx4EyW8UaTXay8luT1jYz4dez1XjtbbfI7dJwNDvtTXSqsLokaQZCpKYKIlsvZaVFS5384u+yqhaY4C4yeeIUhhtnsTupLhVkHVHOQaiAVS4YdFTAl2d57r6VBYztUqzJwUtvciYZvKELbtu/KAMrWOY6S8XXz63nYQKeiz5+f91wnKynnyH4C/7ybtKMS2L1nw/pFFFQWHpvWTsvlqycPZSpa1Y7xVQNGv66AaWnMKWOazZpHQvK0qH/Q3uTliiZNdqXmp21ELQXafIJWONEWTXdN8XO6wVvzmhnL/TNbbIU9+vfU/VHEoaM2We09vNahFP5Gzf5dfkb/LmP/sjm1g9TuQKN9rMXYaSil86e93GOcY2KFmdU1kNRux1hgKdIrxGOVAdoX0UDjHh73xbmlYSj3NodwLUwXJ1VgOmouZaxCxsR4TDlysyWKrlOz89ML92blDrhxxFf1EYuOj5ldItt/pWOOHfJGKbG6KPRxbgWJZORYnc26DYgIe9M3rIj4wL+8ajFXBYERHlO5EfnYW5DNxJbb36KP4ymugbjlelQQMMOdv+7FpC6K9UYymJoA0s3ySchSVsfU3hGf8nhEmo+Ob8sTg572TIK+SgmXXW/A5ERNjtqtnSWZVmRR8YYCCyZHV4ebBGl+KlhZ3hDPHqklPQPzQE9qeaGPqukHstB/W3boeNyt2VK37c06EcQ73b4of+iJN4+3cEq5CHM5HVgd7ohYgqAXZ44IJGBZmncVDeDSi8kiOSVRixW+GJkPr2hlNbK3O4gxOtaDXJmU5aJczjj8bc5qI8sfdUiWAVE6BKcv/pXUVxGtTTJlTqOK3tqn2r8qjuejK2hMBM0Fg/gwGV0skqlL/J8klQYUyOKbZajcMRlmmFzKJS5FJshMWAx0YB0zk62Ffiniif7o1ouDsDydTYG4DLryahZHUaSk8WHq9pAUPXw6m885px0m0B4gOhy4oUjPblkcxdCxSzLjcq8XaDrN3zo1mFmHy7kk4VuR7VMtFR53WXYJTLWh/2J2sdjmzFDvLO3DszmN6ONJIJn83Fc4bzNKd/5mPB672IimLTqPp0ZNcvon8yVdovgrEMDRXH17oSjXzsuDvJBaTV0peZLUio0a8/q2W2qr1WXOQCfwANzx5enANcXOKDtbQ2TOPlEt9h6us3BGT8swjXYzxD+jkOIdOE7c3CcHHCOPnKAsyogn54m4t7VjH8dQPMlmYhAb54uR2oOI6zG7HIhLvblwxqw04kt2U97rQG4iMLAeyIxd3IfJj7n9Lt5dHYC52d/vXK2HWju9onMyw//+qE2cntMDWP11tYEOcv7iiwoZhzrvJvxCG7cvpFva1c7BgQMRnF2SRk+gAha3QcUyfiW5PkOBqDE4n56LRC8FvWt1vLiIk89XLJQD0UgFaB3Vhp3wZLmcT4LQEWhCNZ9U65Ows04mUFp3IDlh06CkwErsFNKn/mTxqnpVqnnwMcgNCKQsC0hFmMoYfyZNSVaSVa1pK41CbLaTQ9XI3m/puhCi0HFMBRYSNC/jkgC0Vh1XEwu+HLBInLl+G8xtRLXtKkzhdoFU8+l2V+/8mtOvGisXcEHk/SDQ64LMzZtf9eFOkPDnctNO5aadSZZBEy9Ow0VcpnOeFfjcwy2b79gIBobfcKJyfkbH/ehpdNxXX6PjPr2N5PfX4d/c4+i4T5+j477zOtIfgV2cySc/tLPb7S4+oRNgNffwP+roA1fHhqtPTryF17L8qzGxdHp3EI+PFMXRsVyJH9M6C8OR0vtR64DZ1wUXDQv2HegWDEJOUZdrC54+eflyH1aYUu1GqXWkp2roDLb6W02VgJm1VniRrBytOZa0+sXVDZFmCl8QciULjRMnx+mMIEfwHgx7GkVhQf8lOkCKSoL+j/DgDfv32J8QNS6e0ECpkcGW34b2qvMggU8XftolZlvdFgkrg2zDSM5eXUmEqRy10xl8lYtqzEUKR292oq4viijofOulzmKi3GAyw1aab76TAuOK10tNmw2nGedWYLZPU0CNZ4Bl93slOVVJVBqC5WYxHEw4f88ziqfjSkZK0xw/rq11LcsxJl8YTxFXjubnh4d4AECar3Z51jFBF4EZL3asdHh5NpqDilTTo1F/Orso4QjCKrCawyfAHX754tXfP3sje4avsAngAb6jCWWESomw65gDQ6sVl6NAMiTtTumG+lUa6mMvaKmcpypeeFHpY7MPIB9Jtx9aOKqWyC2CdCqwbk0ODuSd74KtFV0fU/YfGCk0lLda9pez55PLIIxvyg2y6Agr3trZzP1IRcNAEtwAuj54vNb5bmGQLp3if/8/coFJCcKnRqTQ/LPYwTuZelnTsFja5YuVLTQyxFoz8o2vPkMC836EHitCsC4DX6P8Hdf/cqAYEvSKEsVa5k2IQxsLINRH0xHbLhGWh744iRqgFUp1NU9CcenREEeGa08aPgIaCNHZ3NIAiM8PeDCee3lRXB/CHe6mcOkeXnpvIF5j0elHblZcaOuhuh7haZYXVBklU4/GUn7QrQqAVEuwhDjzSVd9uT4ooJjmrcLyQ/VXfMsJXS/uxiffaIjf1lbXabnlzsaluh7e9OAdsKPtlqotpRuZPpLiD8EC6p9II8Ry3kqVy2Cbz7ce3bmzyobAW42C2eEl2u9xUjYf9or7NSPjddLdtWHxq2BKa1Hwh4iNucc5uFU8MkknEbgEhpf5eS6N0J6F++gdNDxfRZVd/hyZ7eV6/GYnw594LJZKe0t75aAIHLdynccnhEuIIdaxRrmS8VMxTEgHhLJT8H4uO0PSrzP1Org3+fENf2zWXb7i1hV38N+XMJiutIN05seHVXlns3d3s3dvs7fRf9jtNPmZTeNnVpmQqFKsWSHkmWzPhu3I3zNEv7iWDvqZdcndrvJTupedUm0Z1vFmFTdsx8Ywq/lRKZBx98iVP7hxb19/alJ9d/pbrTV6+Nl0xXPaAKdl7n6km6US8jKcAjdOQWf3D3Q/po1obI5n8AyOOYY1IqpwOjs/fuflH/FeWiYG5LEcmoH10UlciinTYMXTjl/GQ5cWnmeqWooHgf3IjoA9++yTYDV+yYG423ogWvwL74TFvCkQ4IvPyTVejLWl17HBy18WvlpmwkHDgOAR70USFolpOqqItaBZqk5n54tRz2MgYPsIYOXoQsEW6pFErOnJkjkrlrmauybC1jTcVAYKqzuPUF11zABTwaiaHHqz7zSl2lNqJHkfIW0eyrEPqqr8PYK3Z2Hqack8vaq9/lVfw/kxbs5QyTeUp0Onwt8mW4MDwJMrvt3Tv791a7mq+8tC9/AtVg7x+iWkbdT4dS6mZ0khJNNPTEbYnuqhnznovPcOOjH0EC7jUnD/Pf1s1JUGgYaQKcMX0tZwm227TCLRVDGcDKslzSjJ3mJZNIzdprTTEORN3CFfEivqT6ZM4rswprbrWhHfHn2RgiwPBRZDpCZs2/DjQQMJ0LznUnV58uq8d/SpC0z0dgY0h5Si4YvlZAqgppbp/DzlCfJKXrnJRp/RGe17mG4bUuhFmvQYwIrHA5aBFmM4Cjz6pKeaMeznS+5pjbT/Denv+p/5nwG/SW2M2492F9hgonhvqdIQ+ztIxozK4PMzRV6/+C4IB3/85/8ROKAl/oWdCYYi/rZwFvGAPqpODO6VGWFgYCI/tRidVfOwnmDvo6N1cVYRw5NuzrCQqMyyO8EdpFC3RXV21hO/ZGvLUK3EhZnZ2DSeAYiIiiYTnncW0oUY0QCNUPKjfj+6MswYzWwgo19IPInYEi8knFp8clkxLTPbGQgSVCA6c4fVPAY5SOuS1HIRgX0Oa323iFhsmAjppPlDUnYS3OzVAml1xkgi+HqaIAYMCjDC/I0Vsgmz78J3of2Q4/DTk5fP9vaeATJPyH1gaxenVZDgO3+zMdx4+PABVL1/c/f+xp1HY/75YGtjs3rIPx8dbT54NOKfh8M7Dx7e55/Dh/ce3D/kn6Phg0f37vDP8eGj8Z1K/hweVVtSdrwxfrS12TkQjIsPk/lkOFmw8bt3Nzbv3ZXGH249fCAt3hndffTwkfy5ef/hQ2l8634oIX9ujh+NrMC9wwcP5M/7o6PRPanh8N5wtHWo/Rg92LpnjZ9WxzZu/CdNbz7cGN8ZSi/ubowf3Ncp2Bw/3JQpGG2NH8jEHA3vbjyQp+PN+xv3htrIo/v3penx6Kh6qE+PxsPDsTU9mY5H8+ms0fjhxtFdqfButXGkc/rg4ebRfamlurd1dH9DGh/fvXv3vk76/UdhWDrphxs60w+Gm3di4+PqrjV+NPkQJ31jY2vrrk7k1p17D6TuO4d3H90/4p/3HtwbauMPwmAfyNo+rB7e1wWq7j06eiDdP7pzeOe+zP9o8+jonnZpNHp456E1DrsBXAxC23e+u/v0qTT4nP/xz+/ubtwNe+XgNyr/aXTAMSJ6ynejS0ri+wRdfjFd4oleZ5uBrwu8BJMntLy+g8Sqq1/fA2vP1wceptOogsF0KtWicvrNkx9/olJ1IYnpkA7MLGgeEVOEZCMSpcR+gnyG837YIK/dQlH4hBqcZoiXbHKn2LpHtMvIxgXB4UyDGTwLF5jIM978euTl0qSiL8xnGK1qFAN7+HMgJw+fzOfVVVmykVuSpc80Xy0uxjvsTc25mJIG8sWyksAlsQ8+9V+PbFAUAGvvtzIcn4phpWNqC8PL95l7KnpzJL05Cnf0nfAPehLGtD9BzwEXdZCLFmws8FL7R/AOt58YKx+tF64ARj+OAp5e7aHyOkt+Us2PR4ZLK1rxPECRYl1EndHcOGdno6mmqA9/dxY5sOsana4CY382qt6viQ/GT0D/pGlhGdgvpJ0atiS5xAdP0AfT0pQKqCDGCwQK2Zvnc7dhHJu2nlKucTUTI0tOi/V1ayi3A/+V8Dn227Q42Q6K2oCNtj6pBzsayl3YnYKXb81v3X4Eth6DoBlWBmNvtmUs6rPeSGKmgxCOTfofyjcSjmXIlgtKokv4fVgcldIK4T/qrMfhlXAeE8ltwLzR5Rp+ruWwtrBivd3d+/XlM9zRKTBqUHRqATAdegjF54H1CE/IXsVnNKp0cuxhYWiIuiFMI3wDFil1h2T1kAQdoyOFEO3F3SnoUUumVAYMwVkCrUgY1krtCH1hqaQKYXrCfxIaFT2XjFEld4h9bUkV1FVmNl/2IpwksWkuKYRYvFTOKonDTD/sQGkntm7/SRPGe80CGzeRBM8jMnPojDVWqbANKMVFPN0zl4uyp4A6/EwcGMgTY1s3WpZwL6mgBNcZ+W+YuXEvS+hYuE74lAKHBD2HG4UgwZbGIexSqHPCBL/P/aKyBtfe+E2opOT1FM5k2OZrgREenalDGTfFLAjpYM9TBF0PyRroQaKPehYQuez2XET8yajRukAPKTgBGWjuFxEDjHDOJwvL0V2PP5bLkieF5ywcKpDBt0EGnYe9M9K7Eir3lqH7YZ6dLxe1OM58w8mhnkDHtgCeJZMZR5dBQbUMck99Pi31RWpeT4Qusp2PhaRx0KyQ8NHEwe9ZtOkUsEQC9j51MpLQ+J6zqMoqYZefwnSUmmWSjqsiCRDYQ+3HkS58FSLC4Bo4sEj4dKimgPVJO5moNgi9WDRAx5R+8GXNTbIGLyMBlLv9wPnO+hR3ekpnhAH5MTAloYxYIfHeLEv7+oNh3fGswq4YSKpaa3BRhH1bVjClpp+HbciPaKQZ/1nuws9u9MymAFl2GqGn9raDrgSW7QgqxlSJdzyQKNKerD9Nsbt92Whiaq5gfyWkgHvKrddJGi17SSsoht6pq055+n0fSmvwsfkFw4hfpiulG9N1cADaUeIxEr5wlFc3V59cYR0V9h5MnOE1urJ8BJgCfaZJPmrWNg4y+9AreYR/5iy+mV1ImT4P8ttFspK2zBpWzI1IIRriRen01R/Qsw/FwDE2n/IukiJ9UQ9zqhT7WW9xfWWTTOMLD6FDzwktXLfrr4BegFDV5mPB7AQniP2eNoTq/e20EXw0Hjd/1ppMZL4NetGTCd3uHtQNqlCu+IbTUcXMhmMZTin/qro01rnFXKLtTvYNB6LHm83ijMdnh/ps0PpFeCt8+XpeTebBB09kYIYoL+HYDM3pKomYcS+vTZ1D3poRAQbQLwTjZznRLEvVAik3VDk/hW9AmBZDG+HcL+Ad5OYmnGWhYPzrkJAj3NzyVP6Uxzio8pR/6UOcP33KP+VxaFYbQX1KQOVJbngQtb9cWCmGl9ELspSSIDlF5MpuOcqOyJRUQWrZl3+aCc6nkuA8/JOyNWbkaGrpy6lA02o0+Q6+GPjieosOiqMes5hL7vI6yKvcu1JVT69LwRiU22uDu5q7vNvX27QhAOh9F3ONyoXq4NU1O0kUEViN3scIOhi4C5pxWNhqdG0rSrt/6bKR7t8urnu95+n5NboQ3lzc4mfRawp0lZECE8DmORUCE4cTUYfpO8JenI9ORzHfBPG4wehO6DmCm1o6KG6jEV9vbVOsTY8ehcU8G80RSzM5GSGj8ljZHUVhLNZ2Ad2t1ill9NYGymEt7B15LHnZL7ClCnX4pLQuieEvqnlMVjDj+QoTmnItqxE0Zh7DORatc6tyRbTUmUY66nRkUgUEI2q240/tg6igc/frtpQ+hSb4EoU4gfYJ9j4zhlIYvJEpBdyCCgopnKGvZpgONQ9iAsLgMnBTfqVAIxA7yd3Z77bkH7p/6bRdyu5sww3c87ybHopeyrzDZGf371Jd1MfW16r645Nq+WN2o+ypzkk5LrnLut2uED+yUTPac0Nz7D1vAQQqdTL8OHIMCAQSViNLIYLLVXJ1OFUPykL/lJeEAor5wza8d061LKUbwt5s9DfA4Qh8t3tHbqb/6FHGVXI/Rt/cGEKMlD87CnQaKMsHKB/7koJIKdKtMgGef2hPV9StZxfKoEuF/9rt2/blzIl6leyODk9zXaGf5br5Gg1cMnqMU/JksQx9oAaGSKHDOSOkqaY4PRFgqAhGOP9xMhXEpfigukyXh/IvlUcYOIw5XQ5D21WXUMIYThWT0R8mu5pYUXxV2QdUS5WSXA5Pw1utA0VusffbK6qSjqQmW+qSdnTsvq54tcTUdPmFEQO8TMtCPGjYh+lZNxR402TBEYsOczUtAbp7lfQmKjCZbOmMXyJ+J01M6UWyrorC4VweH4cG03eWpt7dM7ExiR5MgvTSYEtpgZUbTM5EzDXUMx9jI0wrCWNsJU1aaiWTUktPrLpF/T8p+2RZ8xGvhBbEZgRQXeNEwn+B2pLyarkY/SKeFfSCVMLLNy4RqaiyfHNILb4E9aYl8Y//7X9iFv743/6v2LhI6/lMCkiQKQBNIWRaNWqVYi4NMNJSB5UntjXS9YePkgpA4FAj2JjXQkedTaZeqdL147BvL8DIDnH/lhOEVPAWUGti4o4rUa30GnoXP1gx3ab8wL71buo483zkX0bXcaYjZhASrdjIiAJUG/VVwgiwkbcZsCSRgLKwE2WwmChmxT34fbWsyub1V1Nd1KI35Nzs1HQfPbXh7sj7viynkOOVZDoB2trc7OI3ajMxWcVNC9V8WsHtrG+Yb29nJhq77u/LAT6weCe57XliYVcoDVfwknRb/6wuKR+4UnmEoKNxaexSsww44yZ6URm+I8apOP6EQvGLZjJz/sfJTzna2chZb4D322xzRhjkngiyzcE+kPmLjtLIvas8mya/Mp83ltUkxvh9fnZGR7vQNM+c5DKeRm87MpCqKgchWQxo+Du58ol2Umsq/7FWH4EgI92xDfK4zkAceSMJ5z5azLAIgaNwiOlHllR1o39P4Ih5fdmTr0OB27J2zHY54NN1BGhfKHuQiuEZlTbRs3FgffNGm9ibzV54LUyB4h8bj9DNE3xEet1glqI9UIGRMTUo8jWNiLRCxqHu677af38Q95ha8PLfWwcHNaBYy78bnc+hrVmOnp2M8KvsSAAcguJYMoaNegx357QVy8XoUV+QxCCW9L2QaFL5tB5I2iM4IV4z6or9Y1QaCZXrVi9rG8DYQdigJtDNf8xuKcY94WMYLGL3YRY7Il2+VrgXW8hOkauD2gyDDUvtWCy14zo01dhn/yucnk0sd+PPcPACvF95FEbrfb1vFWNMwN3eZ22v2D7mv4//2Z8dFDtpR4VuuFfcb/41H9SLbNWLbDWK3DmgNf9eHkxKK+XYoRXZgwhUhFIfwlRKehsufjcusWQpjJzuh0AKuCtuWUFbfE1nqG1/ijoT7r6z82XaeuhwS5jjQhJBfCxMMz4oOqbFCJfbz9NF4A9huAnPz+OPTvGpZ2ibO4WTIqGZ3B0tS92JdSV9U7+uCmCaOLs1XSb6tu9LH+B+60xn5Kg6SQbV8FPwbp2cYlVy0e/vO0Vyz2uVGaE8pdW1UzOSyTvpmo32QHwvyv33kgTsQIRZU2m9t6Tngcf+sil437Vt/ykMy746InChgmj2pw50K3NUjHLFQFkX05Gp1CAa48a6RG0/zBxyjfUc8IhyC72MAe/FW6Ans9urM62DGE4AxgYKbmGEBkUD5tcyqtdkLmPvM/j1E+YCirjqkylBRWKIqSRlcJIaAxkzozyq29UUCeVQ+i1q21F1Aq6tlCeKLbZPKLxhP6GNEWgMTza9KT4pFAOvH3M2UnU0SHKAMVNRrQh/RzoqRgOpUwZpUhvJq6eQzElmUK6bS61JKoAKILJTzYidzLtQPQVxZ+FRo82jlA0ZGog5n2viimeW1cLElAR/o+mRnI9TnT3/7eiqHLbppywnj9dRK1BwmBkQSOdXTDkrjG7Y5zx1NcVbpk0mrJftTsH82U82M5IRaqoPDug0BYlqQPJpSeAVQ5dnTtORy6nrOiRJuhP39Vj1mXL8MzujMHPPw6fWYictU0c03vlUCNSDbj7uXMzlsF9Dl//CHiSseWE9BjXw69Q7LdhpQgJd37MI+f2FPftFEctv6FdNCO3UDQVss0ZR2vKuGqW4OfUq1nnoErAK5kaQ7KGR7xfPG5lZIfWXTrKPyVlNmoFYizynpupnpk4xIkwWmirqQpXZUDaIoWHsT20tdWsm/+e5b5Kuo6PJc3opMIkkQWu2FddeR8eQms7BErxu++70i5fNPvZM/1+Uov4QZfpwgmEducymWh5JdJD91dRmiVSqtmVxI7kz44wK8zbOydyn7ElphPrFG8veMVnWVUM9QjdTybPIh1UjyouBwEePCAqeyHZSKAlhFVW2Di6zEDg3HPC4pqk/HC0vRmZcwM3gk+YGogtRhMlbRWqNWezE1Uu8s7DvRFAX53I3t6qx7IhgjH3d17S2f1IqW81i25K0tv12WJUabOiTlv/1krJ9ZkLYQmOlhRFbHdH856ePxfdSYyNbrGHcrUC861CH1amD3YU9bN818m6nd6kW82WLFVkhSTEL1ESoPD5+EomUNp0NDT8s20lNOKtAMkNSRMsDPEJ+7PMzsQb+4TyQxui9JVuou4IlqBidz+yjt4r0SuFCj+UXLT7brTd1qmEz1oASsQL80O8/udt4F7pAwAAs8jq2iZBIqR2fyIs29MYtFjX1YiqryXaH8oqiRYQcSAiex30NOpFPJRlBOH3jau4oZMTYYCtzCQ4KBPIEHaS5AuTr/cJU1cjPHEcoG8aPrg7luG1okM3BWgHvhJAnG66WmmC424ucXwuH2FO/ky/J4MtMw1+evhcNPUZzNyTvjVeqjjVyPirkDFOWDzxNPyH6m/Dgnn6zEwWIlvzAsdiNCYHzzpTXc4igjjdwaprHNlYNwwg7YRSg0eMvzjWclKF/XtLh2FU8e+S36b9vSuJmUrlhX3V5EqUYfw77qtz7IarSG8QzMhaISxGfCruQA2865v4Ghzgy3BzLAGqmmhY2IpsZYSzMws0uYWOIyO8VoW4KAJwajnGYZQwXrBpl6ca4oLTsKa+SiMTXIBKa8FW6pyItMsubDwsDLjAGMBpkItIEQDmjk+CGMrqsjjKdNK6m4eR4wtS3OlBa7nEmN9VFIKq/7/ccoRJc2pPZ8eZG6b/sml0/EGWBkii1hW+K+zgoMTz11ocICiMlurL+3/L7MAeNqcVXwCDWPsNnSiWZRlV99UkpO+tQWP3xv/73Tj2ZWoLt2Y9BEbaqIkCF9RIVQLjwwZkqgGXYsJxJXc4uVHyuBm6PWIGcAVbwP1sq4Db4dKBCsyEO2YSqQ4cpuWTZPlDZIuZ+bqwwqR+obtHNRcFAs8MKup5o3vjDuvarvcNR/UZ3YKB+gvq7nt7buzAEdD+hM0tt+dp/qhlC8lz3kV6mJPfO7iPvflHIm5WESDkBQ9ISRagwpX3jCCPZI59ohFmIot6v4W9eiF3ke6yRMst3/+a6DLbHeTLOXlFnJQepG5GDxPLHSxleD4VcryldOVVjlwO9i2sJy5Pu8M9M3q77DSLtXyKHO2pIfDDHg4TrxMe+JWudMPVT9vUVCX1N3fCFyXx/ePZk709N5ntNmttMAlMlYowKcr4i0WNBpeeaToEeFa3JbX+QBq4T7m4Wv4xjHlq8wUIQ+YqW3LV18MYbZcy2tLM9TAKYANx0K62BLMXL/ihe9rBe4VvfjbZstEPDiBwK+MSwH3PS/vunmjVUblNzaL4h5zWyOJ1JpNJAjevqxLEQUVn1JIKDQ3hyeNZovqKz6tzHCanjEVRcVzHh4OFJhdQQ4hIgyTiTDsT8iJNuoD3p348GC52NN631CiT7a7Pn9mLiO8uai78nsK/t2pw8m8I4NsyNsvSW9rvHNnH/EIvi7LWS23XckiEp5V41k3ZK5YqsxiLuyjvNegkoFSFV/HGrsLcvuzW41CAMnDK8cWyoqr67BLjZi3lxfnCO2deswPcMhJtn6RG4dqLGV3dwjfrudhsLhfMlPut9eoZH2i0RdnMB8NcM0PbQdauG31/Po1lLDLtK0cndLLkzW8wmeTpNjR2oJdTsWfMgoS6ZZtwZzEO5+NWllfQAzD7brwIw14W71jyysFzWZbysoKjaCbbsLBtf5Zawpp2+PXXnjQk7NU1n5ixa3zyr8gbrBN2UqfOTm9EkkeehVlmKTie4d3u2ZDuSohmaA5v2z8prfV2+H1TYTPmznDWT/HhgwL92ElAZcBlzffKMt6QkRXZwKbqnuTOoPTCfJ/fF6+lLKlIG+cPvw7n7t84RWsNoGfi4y4g2Y1G/tCWYBaXuUatBDYrnIDi54AVRXnBqzqfei5PC7cnkvZkIpjPgO6U7DGRA2KxoyWZnXvJPSjqlJjgxlzUp7tODR1uJ+g6REFnrJ8wnKbGqGppqSBaV4JlqIuWRxd1aEMYUtR3CSt3/jd+Gn59+6fOxJTX1T7ueNzu0HEg4l7P8XhTcC2TR+wynJyhoZ2MLNdJAPaBMhP9rhja9T5h/lvSAXcRn4S7Fv/334rQRf480BOqo634HnvAW8Ovox4OmzA6qgEgMdxrBJsqCWdCTE9DwZS3QrDxkC+vFYZ814XosK31W6bNuI0TxkOqYkvnDL8XcUCdzYmzYgIcYrmBk+5PagtS26Z8iE04MQIhk/4deTfZNu5sVrIe6u5roMRb0aJImqIvCorwCdBu/E51Ehv3a5EBaLxFBf5TBspb2LHl2ksN+u9Zdp4tc7HOZxJcEHVNyUjNI9Bo+CJg25Ij87MadAeNdtYgxueI30FGBs6VDYzoZ5OaNjLrDC6FONJPRdKDgDrTwxezJjFWAGAl7qKBURrIp7hLc5/LpKzSY7KCa13g5M4rm1Xf8YKepF0vKFeoON/ubYfkhYnXVVXtY8x/JzDTuLv1L5qNrZJhrh7LxHDSe3NqRceaIqjxMRyLY5YyDrCMmzXKEXdVSwKUSqxPBXXn0Vaq3EqA1ewiQHNs0KmEDNI+TK9QqFktehs3MbnEPGbppdTxqyY1NXbbg0usNrbvBrmIRGyvRDU+W2eYTyIO09cJeIvCFA37zjn8CPmK6Z0CrPdZFseB5Q/QeFOubaQ/ph/DGrt9EacVYBqCX9J0ECDK1BjSCfeaGUy5Lp8ezQqzqWzFbGaPc2JVDk69ScyI6mY4xk5OuWp/+UpOpVPexXpOpsnVfxegmpdiWzEtLMkS+fqTTpkaTG0rd8qVW5dTNt2BUX+HgDsyg0JYk4YsTZLpc8ozjhNLQb7q2jJhI9lcIatGeGlrgYQnlT5RfvTrir0SqsrX6rISVv6ndUP+H5Y78/zMzNjMzivSbTPl1/VbK1GiXR54R8a+TovH4T069mB/qZB01/iSFDMbgvuR2JFFAEsIzZHy5MiEICQRjYjh/QjCwMC2amTaga75AascvN/f2BHA0fHdpu8mMQw1uwmUC/6EG0XdElbg5YqeGyonyEqH5HxA6XetCG3z6/PiwhDfVEbybbxWdnvzY9D+2+KNbzw+nmYsnvXwkPrQIBBBXyqTbykD8ZQidu/8ukx4jdeeH7vaqK8yDF2MjydrqNZAza7+WarnMqakiyF7l9G6uPiUNxuyqluS7zPq2Kv/C52THbes/dc9iP8xp3lQ9nFf0Uwlbzw06U9cc95P1sXb2j53P2GVS0WSEyivrAKnSklEVnjkpnSp+xcd5ItXVVkxR/eS6bGVyyeCmoIlhDJgQlfygpqJHVMSgkJv7cnOgt3cz+2o4Y9FUjoOXLO+U2PDTx1Lgt2iVk80wS11qNrTV2Ust+IC36mX2bJ15rtsSm1qB1tymqgiG0Ryw79V0dLIeW9GitQyUmnYOZ7cZcndna4N98xksgRLaFU9ObHD/fTMUbyuv4Acp4Gtwc+PyYNZTYIrrIqbUktiuslNGu2svizJO5thkVTBn76ardYrjiEm9hPP8nMx3n5HJT/uNQ29emwkZPKonNdXf2cn5vDoJp02UiZaPqZ70r1881TAVXqNM/OMJSrgV+pIjKbwLAlizwEYq0NfXX5A5sH4GdPZXZzeIq6nbNu5zpjU4vjmtQdrZKzMbXJO5oFdcn/VgTucz5z9VXpHx2pMcBD8Qqh5L4hi4uUb0hgHgTybhoO4lPLD0COvJJZRPXaKElUO+qKaSUEzzIawMBc7TJHS79UQJguarqSuIWFZK3fSRaolOG5hyYT+M46AnUv4A1kZRag8wAvTrJRMvDKyrGKvQ3/3xAbgTBLfLq2+zV4FX+csC8ROpUs0G4nQhT3CSL2ZzXAiGbzSezWN6PwRQ/SYDnj4/FcTIylejRnJBOlVADuKCUjmPULky8rehAhFG1+C7SWAnqUkrX1PXDbh/q9OHKv0r+9zgMaK3CEwgUwEvOrwqhiN6suK4AcykhjLByCtnTvEwr4G/Xufv01GgUUcWVsCZcuEMCWBkPlO4h2gzTSFnnBBDaVNjCaHKxABDeA7MACZKnRpdjEIjKQo68VSmqBSYTol/5cQK0mpGTN7xQk2aV8b4viNBTuky+NafqfAwiImHdILksmiuiDyfmhw/bBvRUIdO/YJf2jF+2HXYTGw7QfvNITsgSnF2sc/aD7oN1znCS3bEZc6l99EXr550EpdpToKjsGTFTpGDK6nvgD7kZFkjMUOW5t790K1jBgJo+LSiXzBH248PUpqv2otB8f2zpy9+fPJyd1/enIcmDlz51teD1g7G/LsvAok4Dj10HVWfzxgl3so0fvzUs9P9UdaxZxPVS4MzA03mdKmwBG7TYOlse4TFxS1H6vv8yY8vXv7KDKcHtZyVUPEzjsNWizplJSI7buvyctv94fWbPa3nMSFRNSUM8BjTu1pINTO0oFmaJZ6+fvnzj68AXaH9e5w9V/tCMVDE1cCdq8PTUV+oEGoSq4htye9m4ZBXU59w3J+5NBz9K4U6C4rrUb86mVRyf2NC9LdPia59lWR7Ukk8oq773dhaf3HGqUlNTocvsL5tw9EVKG5pemFfoy25jWG/3++3j6Olxl5xY+mvUumDjMWMRBSO4TRwilb8JPmuRJo6H52NiLXjrDi8ggJzuK6OLfPzcOdE8EIF3ovJWSfLONeJFkUKsfkZq7jfeT4ZnTDU/k3oJP41h4LOQVotHTfjEyxC/Rrqd9RX+tdlPA/7kvEb2ps65GN+m3Eu5pDGDBCG0Z9yNedhhmfzcObDFNrcdhaxCpod9DKs0lXIWuq3EYFmSRbcXVTP81Atzxc1GmKWijB3eIstxEZWlPKkA2vXJDUWMIvAXWlQ7Bw6fdlkS4GD66qLlIu9IoANU+8qD+70VQrDoPMJBKRFmcBvYrs0tJJN3GfiKrJ8Hz8BLAmdqaG2lfwGbNgYBMGv+hPlkjSjwyKur6Q3RqTnRzVaqhKUEbtg5mBpp+9guE/x90IQcnkVFKLXWa9O6P/SLT719doJzEkl3mHKBsX7IqKSy+ZUfWM0LAD+v1j743/97xv9jY0NRXf/uzXFqa8MpLM6Hm0Xe29+fsYCz8ON+IxsZ3GleSmmM4efpX2IKFoQC+HFI6BHCjufwmclTUAOfsEwZH6EODEXL+wndiIeKqHAyXJyVt/2ZG2xKUsbOk19Ykx/gyX/mMxuLNEXvAu/bUBj7APZJ4AA6WwnEUOj91EWr/74z/+jo07oR4ItkcWlBYFczHj72qSQFK+UhzJ+d4mw2bJRzh0vfd68OWo9SybX0HbX9XEp7vThqfa3VvKTEwrRa+HxIBfaLzB29eb8BMArP2zbge69vrF9nxqjiCQyoqbH+eallPJVrGpQR9PZxXwMKPujn7dSMuWyZSEjlssX91jf5AnmlOXD9DV2iJ5I+z5IobJNHxebGxvhuADd+0PPGrbzGx+A9eSs//R0r9uyeB+04zzanzyWgzV4Uh0vGrMXKiB7sixv72//U+/gtoa4oLRkLA5/9AHsWzaZLBZFgV3snAiDg/m+dqNhYdyM3v6HEr3+J4me/Nvbk/4SsSiYxUZ/b/8DimZlcGYCGWIU0HS2or3+cvYSoGpPqwV41OZ6qviD/ecf7IdvDxq9aBZZtQGNFDQ+sC4016M+SelSkVqUBwOHJZo7BwfAXVIy008v6uXWNMkGAR2iwZlJ6yjYDiyFXSC6vRhwIYQFLWUoNnLr70kG9sg/InMetvIRNS1hJfSPcOpePfnx2fdvf371Ym+3Lxw1OslV63GW5IHEcG3X29GI35aGxGZGaozk8anBo37MCuPaiM/YTjhWfaJcaD4zUn4jvrKEuLblU5kJ1uxObacnHiP8nULQTGKQEl7pIZXJS5ClTPmhFRw0jpicqX4h6tJ+ZzvuBuEoLIlMrnf1yU9ivhKBmCFkT9hy2KRV/pXdk1F3U6uvoeZgNSbX9NRVtHHBpgXNOPZ2cS4x6l3V2HzGN2Ty",
        "M9WI+meaBhB0nVvHMHc0p+ZG8enAwgMzpj8cQgvwkxrY/cf+YsYTu76N2HEHsNeBkKAQ/k53PC7p2XQZHiJBQYemM/zq2v4J//s70IgV3X2vfprdbo146AwNxI1ELE8rZD1eqOi7Unp5vG8EKJ3v8CZSNXca8cGKSRPUtMR2NTguqEIymibhaZXJLU8yPWXpcPqxxbv0uFrUqJPwn0JZtxNRNMcsCooENKfHnmlQlUV18CwgQoueB1XhZ6ZoldxbWk6opCellU6RRMppXLFw+ySuQMRln/t2cq0Hs/l7wakXk8zhnFFyvvcik8ZTphvDJEDpZos48Ql9OTupJoKs0hXum+pZrVq7IPWP57N/5IZDd5G1tCuRUYQ6Y5qd6A6nssxR2LknC+KbXoSCGJmgPcymaT5Ehgmzh1ljdJVINK0q03ISNaVY+5yAjARmk4YTOnjpg3C8pqP5D3s/AvMBrDm1PoTBhBbHyao1u6EzxXwsOuzAOiXATjiG4V/+WNd5x0UiRUTdHAvJz2YxIbaxKvyKhYpP+5PhgVcH0/pZojOBcPwt/uj2j96NjsLeaCAK8iyFD5j44vEq8X5QF1ABx0d4kpgwo6N9JEZEQ9OCWPW2F9/sFHfbFIR+l2LmP0cJ3rOcSrGVxJiht10Q4Ox+oVn1posnD8YkyxBdZT4QDvPqbBQE0A8K/0g9bEdT7Bgl/BCFsjCYDyb/PBZCn2S0DyY7fKgLcR9q6LGK/dckmeVHJZQ4zolwO/KcGB/+W2bqqGJwPdHFheE7IifdJl1/hTFuSg9kyWyfQM3JffP4Oq2P82LNbns7ogJEu/dDeFoif9RiNF9CEKPfANMf1Ep+F56WmhQdJITK6/VNZJJDVCQoxo5Z+3CRPJdRRIIfk2gKvSOMmUuCNpPAzvFMfp8iDEIpko3gjPElpQYmcENKG9/A9G+RGUyjAEUmnKXqLlQuW7t+CzU72AauNDFhNW+7c6RZ0qXM5oPmWcF0ncfPuoX7oZ4L2gv+Y64N8mMnb3OVqVl8J4g8lO2WCbUg5clksdTJwJLVqa7zzvkKRWO0aUT+xQrql3H59W/uXMQR0atsV3JFeId8OxIROBklwzWyeEWDQtGhL8m6pPjjW5xhIiPzxui8mgn1uMAO0EB2Bg/QrRvMR9TpJEcxtyj+hOJmHRcyHwnXuGV0chW5M6Vu3Vg9HKbFaNlXjarpHKVg2BdYYAzwJeZyOTs+DrdjJ8gJvXpREiZ1Jk9dV3jOFJeEn8bUZkuynDcXIZ88fNoXhDUbj0hfXZZVES09dsUCI9RtrFV42LFQJzsZqRQHXA2HJfW7emV2JDQpDi8czrLNEEYGA7uVbl2InwGgU9Rcya9I1Eu9BmIh4dt4H34VDY6oqoZDfnqmae5iRBQbRoYSZNE60exbj4syPrdUdnCTdAU2DVA90bWvY7+qRq/0xWH+4nEouR4exjup6vZPZoCZfypdLfX5oQJRKT/w0eyLpjbq5iYmt3saWwdOVqtw0JfvOiY6FDs12QFO7LX9lYkc24pfLdvrXdxebmdlu+bdDbuGpI/JkodP301OhuXyXeZKKFzGx5aKxkC+6XSbZAGYOcPhM7jdoPQo0MHQLDjfcDrtwvD3lF1gDNl7XKynd8iSk+6395lDpSbwHCv1RufnGDKfu85Kjq2yc35GEOnACHRklvPxpEYF1x/lIW/qB5FuMKqB/IRZdOvRf7r/yypdKcRXOcx+d7t2ZrVyUsdrJk59rDKw+OjetESSr6XEYCywxeYd2ypzyd+cUdQuInBCUdK0+rvUpTTW6PyVeLWnajWSFRX5JJMSej7DFd1O+mq860gBp85ihZLkr71G4cTbquSbXmTWQ6V/KTcjmhhnR6PFImw/4c0Fz1QSzsXOfy8pBT4yCqU6waYZjzsxaR5SRVZIbHNvQ6ix+mOEr/YI8h+HjIr6Uono9VGP+h7scTpcCZDHOBv77sUBUz7mucVMLbYXMwtFX6enu79Tjl6yrzCIBhozVXLBT5/DNlzZhaT0ePF9XUbFiJ4uPlyfpnZpQ07j97ldfP7ahM4uG0IzIhxkt40KSPv7HTq74YwjbOMtLNmqhGFlpirsReN9myu/pKDJc6JTIGG03f5E/A6k1P7kIG/Cw0L5JOroMxOpPy78b9FpuX7lWis2K+/+81TUm+Ffrx16Phkvr8RfbyFKFkZWJR1PmSyCY4CiWS5lOIzAIT6qR3jFDIs10INic83ZE198HyHc1YyoWpHnEgwoebKZh+jqC1QhTVMkfXiv2Tnm+ZwrKfx9cCZecmM4jmJzrctx6HT78irICfJHfzhZGBKMPhGDaXbmvoqHzpdxFsZGMxZN+Vqag7ettYQr1J9ubWe72VfpSXzvjaN/iaMDaIslO5BvdTbDHopEuNPugZpKZMCDpVSK5H/FbakGJA/OqelH12BKNMI/VsqPe771W+6ztsN+sJ2f3vSpHGO08G9zenNf80ztStkZOTs1AmosuIRCJLsDOUG9QqgVT5JYGaAxiDTWDV5MnPk2Yg7m6oLce+PF0Ww+54YEC/z3b578Coe+589zXj06S32U3odhGSFV24D9RPtiKjDrKw0E4JPtS6W88VP8jl/uPnvqP42+g+H5QRONsLEmZWxmavW7cOmu758ZswaqjW5pVFwVYSpQPZL3Dc0dD4+aDodmd0gNOQua2NyisJDxMVj9DuoTpn6FMst8TmsKPFNBb0UGzagHg+sC8ahJ9Zrimtb97HQRX1kUdbfkOMLiptEIcpadZ4CP0eFzXhydvAvAMQs9cKTVU5Md3rBWNpxdoM7nxdmNb3ek+ZzdxEe/pVZWGaVO5Mve4t3bNMLGwfDcE5/QtBr+eHs8r67eimEEnH7b5w3+yx6ykuHz5xtvZ+O3qOyaetwxzKpKz1kbu3MS2PmTxVuGHb1N0HAiiYTmhGfOLs5QmzqQfqxfiWHPEdDL5NA/nI/mVxLUPps/Cbx9B1x0h3Hs4Y9VGhW8MgkBVUa30m2+agmDsZvwmo8xZBwiDo2uDZ2owHBDYITTTYNgIQyDf6SAKfk5Uc9NCZayDZ2jv8lKxXCTPLxru1ZsPvowGWEjfSx+OpktwyaeL68Ky/aET/EY3/FviStiwGytIlKDrmN/MhPt1KsQslEHakK1cJDcZlMT9GZUn7Lyr1wv84dphIGhEbbsxXQ5+10oGejr4cns6D3cRCR6pNMa+BakwZNZNSxD4UMhxEnNdzKhemmV4qNCX1Gm/24+QpjKz29eahkhueE367Vi1lghe0afUiTGfIY9hRsGAJJCF1BfGPfsvasvthcI6tbGhiXeEycB3smcHslNDYNe8vYY1EMkcluly2ZheXlOq2VguxbmUNIMfRhdntmdJTafQLFHU6QfrkURwr2h8KkNVEfF+C1eMQlu+Pb+PzxZ/0/V+j9urD/qv10/uHU79LXz1pf4h7e3/untrb/lC7CqVF0Ezu5hjqcuyeXh7Ps9gszBLl4w+9vU92HaDUduuAvgoXIrVLiROxcEWeAUhjTkmD8eLZ8HTubXsKNKXDqhulKf/xguqnd8uFl7I23zjuJ5da9+QPybvPNVTaaB+DSf74almw7xPBOt9jm7pTLMcX27Ylbt0YsKnN5t8aISWyt9Z7nYjwsH1GI4JUnNDJuj1m8QJvK78RFHyMsr+yCtzuadblfhezmpq/xf3pr7C0GedUM1jq7svR9Hy2roha3oJbE4CrJhNUDyVKCGHlV9UJB1+Qxd4Pu3H8IGD9UNyEvb3A2K+mz2ENZ4PJmm8lpCHv9OnhoocHV2djIZDZ9PAI92FnbZ0isxmm9xBHaBfbUctbztWX0/zWfLWWCjWiqzV1KTvTzTp9av02o6GY8WWW/sWQ8ECB4Oi59GU5mEcLvr3zg7gbsJx+8qcCniIIRKLjSLkiDuyIprqPMgMOmvf/4pcMQS8HxgndCQvN3+WESGJdnsgTwHBQw7fF/fEimtV2Dl8AX+7VmK4bfj2Qm5zjQWeUJtv06nbHSZWXkrAQBAYjYNL9rXbM28ARmx+7jYjz8Yi+1+bUJektfxzeXmgQ7lH2ezUzDR+j6ibfYsS+0gpfftWZpozWePxNg6Jv5hfTwLdfGz8O+Ps+GoF7VzKQo1jHSy+H48LpWdux0YLCrPpGivuHorMDUrBnky84N8R6FQ+i0fvmV2YP38an83KuhsPtnLZSyiLaQq9TXmLsxVdXLyVgjTW0Yuo6d0FsBFFMaEdFfuVXqmhrm3sHe+JeC3ywmJh6/DDVr3PfpFMQkdqrXliqbjTi+pewJ7FjPpEcyJBQYKTK1qol8VAsQiYpBdD448gUnCFhXVlCqXDESb8R8Gn10+nxz3NfXeaZCeoIc9M+aqG5VbgvvaliZd8qNTMBddpuXGkubqV7U8RWDnomy626DXCBizPplfDR6H5cW/TGjHzSyfG6HNTHLs4Z4kwA6VaVgYMkfHOLtQnX9TXdqb7FJLbivYNR3iofuDus8uEd2If9x4+uikYr1LZ0HqwdbnH2HXu23MR7ZrU37jL+xJmpNsR/4QmBPZjrpwtSDT57N54L9I7YLgzDXEZtt78fw5jTLh+YIP5yPBCxXqLW6cDMbKau81NuQibt7YAdXtqCoV0KLbzEdnIb0CuEBjpSVczz4pKztb7zT7GqwtGlb7LiZ1rGhP9wBuEtL63exytFD/ZRw6uOaq+lEzO4kKCl8zhPlC2dsFI5nh8JG5L+N7mjLGnMvAC5xNjzuSYi/8eHh5H+qfs8mguLOx0fOzGt6e/COy/SULx2QaGIZk3UDdEhkmGdqpo8Y+5WYbpPfqcd+51JSzEjbgrUU28ryreLwuv1KX+RB/W8flSfizU+8/X7hHOphM2nYTXLZI3DCVh/VVBVL4tXDL1C3+Nv3aDy8PkrpXOWtUL2+idqXeK8WYcLggnDgZt0h5ywl0x54j56zrgpTdmDhNOLKwBA7LXrk02bhf267dhbm8JDe6kXeNiecWDeUTWwi36k/xqCK9XVJI8T2USWeX/JsqHPYz3Mp8+Y3vFd4AZOBrGW1Yvq64CN/913+pS2m8XWbMYoo7Ex5CfaoZ+uJCnIwofwEroVKGAcPf7boEBYR7qSWVACDOwt+rwJmBH+tyJs7pEW0/fETSJPqjKiXWOarm8wl8fzUaVCx2Bo4lg/7VZU6MpIqPLHGJP+vPX/z9z2+evX39yyuofzv8CuYzwQPDX7iF5N/qkv/uLke0ol/piyt9cWUvDqu5vsJf4eXBdr29H1/vIgPVw43RfU03JDso5vR0DMfAMjDSEakqDqEIDRf3aXXFScHNWzAeNgaTzpCSKDEQXAlm8Shdlg+80DSdONXnJ0MGCPaJ/aLAEu+quVPdHxNR3kygPcdPKIWm6sBBMVBO535eGmvERxMG5BUnMNhHgFZjQAg5LXyTnD+LNlHcnclU/7KLp84viRuSdTO1aptBcrVJNCS936YzmaKwrRZnowppO+GorB7a7PoF52c4i1AVmPlluFyig7Su7JvXr3/sRrBQudX4zfFsqSnuDNniAtZCqUndvrONuSfmdC7GRv9uzybmrlBkW2xzRJcIydER5v5C81PznsTIcK9iV3jA28itKkuoeh/R4VzE/RD6Nx21nRmMlF3UGLu92opOFjZ+vXl4uPY+szCXXgtKBu7ki2e7qpPy8F5fLp8xaHMiL7Pg8EeW/ZVA0gtLAkuvSdwUFl6QzZcIADFqilurWiwjr065QD1KzWI8quZB7p4PYuZcpIXBYXjvEs3OR+PwyfQIGCUxvEanRxOkJmOaid8F781uglJvnM6jaopF1uRg1fQqZfMNZ+p8EVY5QdHLlzHJG3zgR1neWBEwpiNGiOHUEKu9VIGhuurme2bCpJwqkTZBLYg+Taa5GxGtVpfTAqHo7AzGIJEj/3A+OgfnTEV6r1AZUN7FebcuoLLYDn7gCOIvGPNkGj5+yvcNpyCnhjpNkotZRPD9hbhQ9SShFBmVHrn5ILHN5iP484548Jk0bjm/SgnXr8wnpV/shQ8WlmNY16AUS16afiLyR1lSb3tQpQU0AsUaTShrAwedQP3sOZPrFcdUOpPABckNvYGDLKMdJsZzr1Xny9ma84oejsbV+YmQQEhfxQJMkcsij88lqIVJeW1wHGmPcaxEFAvzhWSAob0onv5dV1KkGUVN0+yS/wa6vpyc8CxOVLSFkfDdbHKkgLkGF68BOBQTRHTQihl+Q5xN9JVn5WTyQXMNRGjoxbkkg85vTIVxLiRLrAMk76nvvkwUiJFwNaIcSEJMCmyQK073AlqZDGNqasmuZZFf/WINdHkNC85+nkMpgc4LGsDpSOJDIrlP3i/D2UhvmHcxNi+j4rsv/tMzcR3aFYraIeXFH7/YH+pDCpOt2pfJ3uh9tNnj1hsUW4iU6xWpJiHOggRnf31RXeSxAdQhrJgTWhDtrfItfnHC2HZWNspBnefxrxVFvcgD96mavJy+PNj2PGz9/nTMjwTgElhLKMKghWgz6oycQbwBO5ChsXWEoLcu2rOney9ev9olScUwYD5sF+AzYVQ/33v2+z1Zc9zGGBryTss0yb9ynrhK6jc0KJ6/frW3q5HfM/XCxB9JU78F9VlcNat7qf/uuWu7ZRucLdM2uGvb4O5D3VKxMuMhxMcNPP4XVvlgq6WPk+n72EWyboHScx543DtZF5xFOvz6LvvVNm9hmpkqjEeKf7A6Ih3DR1ZmfM/9PIAriBMaQGW+J0xjVCNElW4N5uSrr9T+wr9JksxJal/ehMYPCW+yujL1XGUAQaVpf29oV1XIja32/MWzl99zn/4mKaQd6TFPl/3aEsfljTsoW4PL+O73TtxrL5tEuV9by4bC7lyEL1Vjt995orvARMDfAw8l32qyO1NBEQVDQfljZUGTFX+vGO88ubVdrACsiRL6Ppn0+etNfTLx9Neb+hTl11//1D6Nzfk0eqGuOA/g2vD6JaQvnDeBW5d1AyYsDlb4hwfB9/F4PuE5+/vwr6Tekx1z9L6TlYvoy/jxhIseWPj3ofziZlKx0d+6Z9TifjbQl9Yi9Pax/h9hpUZnvqzqTbsuY62aGTBUF6tVApRIESoGMzBonsiv0uFt7/VidIp/dp/9eM3aHFZCynbf0R2WP9lNyIF48Uq5iakuT7MdUIU4O3Sn/VNmp6fjTNSvdboYXGjtXDNdWTXNPtOmFDtNmvnZvY5cy/1ebW3kQm7tN5tIbX3GQsfKVvT+tTDUqc5ZetAcwt+1cF2b0GcrQ8+/P7sDEUAiTuGbCCnxV5rG2KQbd2r033DssmAJfqjnsYhuXlBEg/FCLrvtOyUcOPz4Pke/aq1Zb3bfzR+RVBKFRe0lZNf+WkEDkJwChDOlk4Q1jX2bGyP/hso4E1ji+1VUwvckYtX1EpRH3CWJmDfK242evjGkab31m4O2FJTCqqW/W1pR5Lm4ejHHpivfbGB+PvWfxLI54VE2qjYP4ugon//hfAK/2fnkC+tIBm6e++wXBZYVlRx0ZTBRB/J/FlOWM9zzWPppAi///zQbpyDmskr06UijS2vHxtxDcXHGfrE/WoslAwQtVtQstg2j9omMVT5Zzm7+IF+PtlPyZnbBwlSNsXJxonLdtiLz2cU6FSnxLo+/moXfBRl7fXRJiE7SsR8gdNsDdc1OXzZP7an9IUlmFxnk80oqNJ+cZeTBBCUtfdD9TbFCY4ikYnQcoD/GxcxcMyIyKm4rATyMicxFAaZo/byYqMafzqx0q2rgye+fiVrgMlDyFcarA6a3319hwTowBRqKNG1ZaYS/qJWhUqhFGKQOkWvgbKlZA6CiElAc/KlaPBn6wDAWivn5iVgBzCgvCkFTdPPp0QwmU0jerWrydq+V31Unk6F6mFIXK8FLmQfLGBiw8Drx4qp6mQiUlOT0HAPBncG5o6uuqHOLHf06iL/8Y//Ogah4IS7XXt096ApQaHK5ljpcjHZXMSTUHF7L8Cmqb+kareOxhHrdaLyOpkKxryTLt8MfDAvxM2CMSvTU5qQGOs13fVEvZwAn3+4U6RXdQaUodc9ZyW+0JN64gpJLJK+zSG+6HpQwTRApfBxSPj8+SD+94WmdHnciVLhDponuvrVGcIhjI80aD8XdVGEJ0nfKhF3zpevL7X/4m/2N9UfV+vjg4/1PEaPQr4KrRhsQJg2fc6aUU+tD8VtqYB43qfwJDwF8l+pMrt5qBGgaHMzilbT3hmAVkaqinl7OtarGfR0GvaVrTZqVEynbfDghfv+h9+mb4nFxS0/B14Wvajn7aT46mkAhW25udQ1l3TnFoDef0wIPpjsvwPVhg7d9KVe/kfTAOGu6hKStvUxAMZWm84pkTG11p7n/nbn6SX4XQbSWq+DUHBpONfECFqFnJu5ZsumatwN5gr5lc6jfFdPRIPYpdQhkV68ajgfgANIfsUyXyhWSzei2U9cn58tZ5hJI2IXof5eCLKJr9eVAwq7h1JbwNTVoAhZQSZ31ekwnR7UFMoMRwljCT/Hy1TxsvloW6sNUxZsNZV/ylinx508w6ZaSmFtKwu1O/9rsaqomfl73/rxKRvyjagrp63wxikeiZ4sqmROcLeDsnECbFisetj3Bms26Hurqy9bytrme0LNxNNh6k1yveFW9kg0AJXVyQ9CP+8VzRpm13bVmleWFKBhvlrM03KV+U4XNYq4c8NpQOsEWZ/Oa6dP5RdAuhYGlHJZyeZTPX+wZsiisliejBT0VOAuyQy98f841WyWsdJFrMvv74qya2nSOmrwUzZrSk6pQt+jcXK1oGdj88160aUZ/Dm2Yx8RZWbEvTqtwO+CsbEdvFOHRzAkn8D3EnRBran606aCDoy2nND4/nhlzBPhPaz0MiNMgV8zwjQwEGR/pDxoa7Q4cUSZea0V0BnUugZbAUqzD/IfMXkcV9m20D6ZlkrXUtrEgR+/QILKELKJ3UdiGE0xdypii25xsdJiNdgrxZjQ+X1QnN3BgyjB9EQMW2CnhpoAvceydhnuSuSoiualL4a0P2d1xv5uigjrr4G7/+F//u8el/+paDjLSNMcmpXwg7qHLQZXYJDyMTxTW1vFGqSLVv8Sy8n4guhXbAS5BvfPLEQbReB9qhseYIjhVxkA1nv4dXaznk6UO9GPBcnr3FZ/CVM2QvSpCMoSvvoF/kg6SKKZ8ltck34VxNPpF98Msqp8qgJ06Zxe9Yr0rqZNzurpDYuoQ924fD10eB0FoksVFa3ENO5mT6P7phJ71l3LkDtLOZKX4EkEIuCp2GhehvMbO3F9/MSV/EfaM/XXQ08iyUmFj39t+jszHQFYJCAUSJ7BTJCJQyv0W6ihDJ7v6V3Xp844YShUKxAVjTYo3qEwPmgmfMqC/EwS8IPieRr8Q5v/DHOApfUZ/zN5gQ+KUsUPhWLa1H7p1Y/uTqbYvLcX2ZZ+fSq+y9qt4IFL7GzkUyj6c8hGFsqMNh9WQiVNxJjA3EyaQBkrD3ounv92lNyYClY4nEZLhLMzt5kaG78A/T2bHmxuliD9MJtZV/80yXXomzMO/HMQI3qrJ5Szl5To/63qGX7bAN1oJjv9S/44Ddjnn+eK2djr0Y7T+CEgT/F1z8NUgT+VWk+vBQkR1S5WVxOmBJb1izgXku1oksVsdgBZyQ8IV5eI6Z96M/ZdO/Bat7hTKPDosp/fyIgjJN7iTO1G4K8CHXyVH3vzMgzzF20HID6hhXZjusnGVq2sfRNKHIttOw8PYAnY6TipUFIPi/RTMl2h/FHl4YvjDnFt66cEHKjBPYaTD6rQ6hs9GXANjNw/n8Ec1JnTVPSvZiIV8S6ys3LAEmGSn4Wb20Sb4P+6+ftU/q+aLkUTWMt1rtQysUnmIhK8Iw02IZjK//FSDYQiHd6Z+9gqMuhKS5s9fvrQyaDBi0zVfNTQklPJaL/L4YegApA31D8zUK/V9s5O+yqUEU4OZB+sScJiB4VrQ90fia0QrhkMErl5aCy817MIJZMwDnWS0pstQ3QU1O2Fa84vpz4tRdsayK7Q+4K7qEgwRsT7wnZ0VkwSllrrSJF8r9dDy3lMHtQZiiEe9bsJ+xMXjfSJmbR+RgmBVF+6xE8M95Cr0BC+ik9BpjvNPx16dTM1qyFilqfnghyof3bcYA4SLCOK3CJDC+vvYAnVLz3zRAZkf4xLEI7MoIY1HolnJuxXCteSsDH87WIFMQsarBJ8YfvlcB044vhj4ABKUs5iWR/fhadt4raPR9yY0S3TF7GygTT4u7ocizMh8dhmuoc2tlKA5vP/pMrzd7N8raPbYzsMp4dk+gTd/oH2a2bWnJyPdUuOIrCViGqb9jK7ZFwitKY+7Emxyei4ZZDO5dOQiHjj35lHKeAYM230ZHVdlZOpiWQJbousyXSBsHScWhoZGBIJIyCX43W7twlu+HlMdLEHLxxHFBz1B1tFfmGAVKeelTP8iPJC/3iny8cpiHv1WPwmsgPwVFqtL6OzZ6SChTIWKkKL6B/A+Us7W8FPXa7nGiBqyCXWxABSKcZQTOVKBUbb3QPSIYX3D+pyA+eFer+m6ThD4sJQ8BIiFAH2knLniMETRRNpa1FT4c8GM9Sen6TjuJESrBX7Dy8troDokxbEk3mVQ12V4uDUEcOlHzcOsB2R2Jqmc5xFUxdTv2AC2+k3tWDhI+udLJhJPKZ7DQjJ7rzTgczxTWwUWj7qtk2ra/lHodab83YvOwXPv0yyZx3SlTS8MeojQyExnOTCaJWosErNqYeqSpelwLPpBVBt68tShGVyR7gLVhQhAS7KDJZw+gQ6zzGBaaWRIQQQs8iX0Wik3YBpBsNv30VMup26lM1ySR+/q1qCaMqFtbwnprW0ARyGv6luhlhL9pt0Ws4038pRv9lYRcxtN+riZpLz5daT16fMoyWV9+LpW67c+jg0Au7jvAJjybD4HMs+ei+cR6nEoYltWKxQh//ov/rlWz4hIrmi/eEqMVjBSCO4hBYZicu4CnPs58oqcaa2xeXLPE5LLtHaHhfWuTicnV3rN9vVXuFNZzenZ+XI03IWJoYyLyHuiDz/o5yy9jfYBB0WHXnAo+VTexr7Z4P/PB30buwkv3CKwF8nfuCu544F/csJuWEdTke1YgIDt0hB2arSssXvhmAE7CH/TD/kEYDOh5BVAh/IaCzMyUsxNmkNxv1gI5rmGnyvm4DICrNbCul2t1QlyA165dVMwgVYiSaJhWFSe3F1dTyHxmfo/t3ylMzBX+UknPqIZhM7U7hGFxqmd9l5O5TGOrv6TASMK8dQIsxRRo+GISpbB8cbQRahvBzEsBT5WC41vpYigjEuKTmE042gIe9DrqYShQCqvhVhSMUF+JkhR/WQJFTxt4zd5MRLAgQoq+7N/WSjCwcvt9s04mb7XNzNB4Lu/scFDnvGHTD/Ax2N3bnBQzU0Ln55OhkOi9dvLJ0hriTcYaniOjTcSWHCs7xjbcCh/6ypdhPVGvvhzaH17RXiBv5rs6kb/wb24LbIknnJOG3epxXRjYp8RHigIG9PzRV1JvmBQ7ZrM/ZqFWU5iSJ/XyRDSA6aCXh7kB+lRr1KJMOknCeeC4AYxzHM6E3VsGTW03SzKdLshQDJEHP4v6Nli1qKx0fj1NoFSj0aZKXLiXddy8cbQdVHGxtPWryuWW++V1QHEeQQpLg6q/ybkYp6F4ZoOJcG1GwVezkA2gJPWmUARcxt4D9gFuOYtCJt06ZOLcCD7FScmm6hS4EdSKOPKctFMmkUVVTJKmNktfM6vl9map5oXyIzNA3fIXcic2yu6QYQSJYVDy5K/vphaAr+oNsmzJai7HzSVhCkpF0xRFt9fxffOLGywLBLbTxwwaiuf0mWr67+Pbn9owfIUdqJqE5O1G6HqKhenWpszYuLMsLt7RkpJXozKwip9lYVITjUqVu29QbZZC0XWtmN8YwyPJNcLKxijEBGBdUJmZT7CQ0ly5S11oVKtUEqsrZD/w8BK6qloWa4hQs8McwdhqqqhFFxo5AoyvV1+PuWWZmXmDLSU5Q1v+pgAJG0Nf0q/VIYJv8McjGqoFhJfH19GCI6vUv4PrZTIw4FFumzDx5xMqRYV9/VzIl9y9uwTyVmw6qtMtpEUccskUAycNGGQM4LbZHdEWGXV8CyQuTickcX76B2hgXQiKGokcXh8NKKEoaGV89EYzrYW2U0jdlQmmMCCAwtnNuSRFNh4bLhkwj8Ko1XXAXHqcAGEyQRfLXXS1bobBmbpkc0zRO8CGq0DZURyNpOK4ESH53YzjDV4FQ3Wt5/O33NWV3rXjq+w4xC4o1BMfQm3ThsuPTMYUL1Ew8PF8sl0AiP6bCpClN4X6naWfZp5rAnlKKxt3eh8isRtLpG3Vz+vlOe41gWyw3Vkxtb5BJTaUm+lV2eGcRqhYlZIX/FKg+Xp+F2zYBS00uWX9VM8djzceH/jnkMcL2WPHjETEyMsAkNzZ6tLdVBdp5cVlsDbvHSu4tsCV0S+o/W6fgRFn/Tv66LUHJ/DMDFH4RIL2+sNFpWZwbpxUJIZxMRX1G15o+RNlE35Sn7ZO0lUlYm+bbpL61Q35rHyn+fCb6tys1FBtiAiQ0qVdRFSmFGAVc6djOXG3MuGqeVBEF6AuSgx6EzSiwKQh1KjZizfqgYyJBErBi/0tey5TCTv0fMLDHyDq9pmvS0ox4cVSDFrBjIfcmla91rugIhcpAI8+5DkdgIZsa5wdw8cr8aL+MPoM+V6QAliIKzEky+QR7RBQgDQT6vPajOXpUraleoFAAmI3upZJ6ZVt0kkzXbaGMXfpbQ3jmaEEeS43emdofXm7/20Wn1eKuwsMsA0Io3M9P5JiSTI4zttLIg4p8SMCrE+vUqjahqV5xsc947klBz2Bebn2x3uEpqso+fFkAoO3UD6StxxMZUpGjvLa5cZ/vZjAHgqfdCt8zkflWVK/F/s4eMb0WLIE0Y8+E8Za4CQgUKxcc4oZ08NgGfQRMVzVmy8W2FvBRhiEpQlO2i+OzI4AZV0kK8ThVOCIjBwn4FhFsfy3fnkJPMz7KgvXObDV1fw0rSSh2ODjcB12ua7hed6/de9O1uu23gzq+ZmOzKmashfcRG/E8hVP2+KdkECZ6pdSXMhhc1j66lmuSoTPJzWAJefBgZbU/boAOTRxCTsnR9spuTRdm0xh82kApnPkFlzV6t2z8jcqoUzy3MHS3YnvaPrmrWtD93rWi92TW1wnXGZinTF9+nEufTZzvSDpKx6wsyYXrwM7HsveWjKeaI7pvo4hCagOQGnuxSjHoROOpPTGI5Wqb9P0lDURh2RoKn/JXWMCwiq8EuGkShw4vSbmozVc1aBHuFnomAwvCOgK1NXUuHenUBo0FZzQStCos/FoAh8pjirUINirTc7yHg6NgoEpeRceQ6UIfOtZDQPvNicE2UUMftefBOLR0obKnBMpaHTiaGcxAI5AwO9gAmz7DJBwSdPYFXU1Q+FtoCGlm6zJeHMewbIyd7OXbCt06ln0lWCzFNB2oFL3F1DirX7tyf+JBLkMBlflcnLR9qBVaImV8Smmbs1JsFLeK5umPvEgV4Jp4QpaHEe7UbhwmdGWhQSlqKxOJIzQMMpeM0CNsM/ANojCNNkeh4vbnzEPpFIa++xOmkgfG3ET/PcriIO7+7EHIT5+dZKWk4sU7OndHzO+dOQYWTSnLOnSKCre4HTjH6wXK0j7XW35E3kx3Qn+pRNvqgsVjWtWJ6SR+EzuggFQUZEcY1u8+tax1mzfOC7qXk7sPfgme+WNQXd0GKJ6XN5jU5HikA7T/BRVQYbNZCcgiB1pDaSYZBX9nwU10uFvtn5/Eg4fA0qmyCVN+QOvun7CKDR0XvBxTCsVC38WF+OJD0K34a5EJ5kxUSmTyfTs3Ph5CRZScyXyBJCUTTULfXBpYF0icdwDLXfdtTgB+t24NHs7Oqa1ZWvKF6FghETdiZ/6dN8eWf+NwedLTO+yNLdhgL7rO0gTrL+1o+biSffaaCtkcTsK4mmi1Xo8AF4D2c6VlVSe/1BDLxW2cfi8PwQiRtdTtWMt9w2Yly6RJSBpxWuvhZ+JtFhid4V16++W2pPZDkI8R4kSTUiGNPG2TLOPmsNZ3EBdflmrUrE+pLN6usVqzH39c9eKftAFinuos+eVw0H/HjDbMoR6m6vOCvC8MsJxS381Z87Dq3ri3aIBSgqMRPhV9gOUXSf8byLV2QMUpqAjwnUa/qZ+6p9JkLbc52Gz11KqeimlWxyFN0b51Tn7IaFl2DTlZOlDi9prpJVgAbYP2e62Pa2cnKXLyVkFRGr8mxxNjo54Q5IetK2mS1SBEVmzfGmnIFbBxfy0iwp1wBFryi1RSPNMlpoXqlZqSMRKorNsZ0d5v8D1lqX1QSRBu5lE/Eyxf/pdd4vfqlvh4R8qZLGUuF34/BDm64ikRgsLk8Dpbr94oVUCVzh1r0F1uMZ3Kd4aufM/ibykHrszwBfs3Z3Y01tVNb2FDIRa7uLhhkxtqCB9fz4XQZ27AdsAY9IbLowbzAGI4rNy0XJhc29oECm4NC+cTOIa8USJuwNbrBcqyfCn39w+F4g6oqOZj3UoO2VIfRte7yDjeA2sCWHWl7Ho07QsaRhjRwq9uI6nnf0dRTMJGFaZs6D6kSmGUURqXgl+Y9MR9bpmPzk2Ml3My9V3kAWGIxjGd8NEUDv+l6KFRMj+H85H+aCoaIJaPUp6DrVwFMhOQvK2JWuD9tVhoNjy45BS7zgp+3meHFIGl3y1dINT7pu6evUzpnCzFJQXA6X4MiQkiCLyudquc6Y2M7sqh/ZKzkdg1rX0HbqD6yj200lKLttYWo0oGJZJW0enuuvLyKqbY0Mso6kPJ2r62zwJbw0RTNR40VW1BBaZx75Whp31CP52+OlQgpHvIP2Boyi/yZ1OFOcUTMxuvKCnbGziLd0IsmIMJjXC5w9Szp77ZlHVfmZx1edbX5cl1tcZ/RLYYJLHnr2sSm1Wl76ZCX3b6Ew8gPOCEoasZKv7C3N+EqWOI86ocg8pyRrhQl/OFochQej4SFAj7SClq7jjePsWzQacF1Q9b3ofzJpKCmi5Pl3o3cI/MRh7HptrWrFnUpdz7O2KMke92ZIzLfRwlY4ZTR09B04dh6dL8rMoiCXv3MwgY/y2chlpTlk9yIWN1SGTNe4IFWjQZ71tpsXdHCz8bgriV1SGtD5MssCCtNMEAhEHS6BVCgTk33ms4LzxNfL6vhVTIa5+/TNi5/2Ol15NQknlTt0PN7O8+ZxQpJbwI3m+WQhuG550uWnNgRxNc0mgksYvVbgYic6k6ZeX5L0rcvvQfydfKzYPdTQZT31hKKL89NAua/cssdFX3szgkP8chbjHxVW/fB8crJcD6yRuKEmFGbBFq+Hau3VHbSi6x/QAIBvFL1SwqVqXkaC0AKvWcJFmwIaIohUKamNJBG9KKZMSw9k8yxlkmvnnUYp8TntQb3ocqhNu6iziolk6smnxHZYnc4Ef/y0vqfnmLeaO8nn751oMWr64loyAFX3bjcj1NQC9bwWXaaPr1dpkylsqqKVLfkyBXQUL7+KikWCFlERbaA+Ufscf1HlX9c95ybHmt4wMcR/gtIQXJ8kvY+lUJ9u96eqwghcjUubhwKmINJY/AYckdWwq+mQiYPlPkQAdz/C7Nh8OVWdgBfMF0tVxiTlmya7Zzf/NKXbJ3fl1A2fefx07ksGCs7QkAnjqOkNEf61rGOUAktN/fWv/1KIb9D5GQ7dXYW0ofGeUdKsF14T2xpqVp1FT4AYR93J4zx9RqoGF9zIdDWBZJDltILDTmRqWf6x9z7a8t5HdyWPWJ55qrgdaqWTfRkffbtTwNf7cXHHMpnp0zt8uhWebprpJxvG62icKS1AzhKaOZu45i6zimuZy/jY2/I/1X1qZQnNLC6utOKVufRh8SdiTVSP6Ssh8+QUw0EnoAzMGTEengFXJuEqApMgD8DTs04ODViJy1rmnCjGr+/pDls6GAgPyEROcNEPlKsaXjGyLSZEcfKV5DCRonNmqpHOxQSEAhx1AsdQ+RM+PAcC+/6W4ItMEzq7eCX5LQxe8a3c3ngJHMZnhrro4jx/8GnA9MggVkxcnxqe7eb0CNdIFwPKbYmZ1bQGGkmdefduE8qHaIMWataSYbCX/Cm5bq+jmPpuJGGik6VuiNLyYeLxdsyOeVjN53Cwpr6jUoe48LCLDmrLdp3KWl5Af+SjGrMdLxsAGz7m4Gr41OeZOZnLUvL54VbMXx6kkKpLXGJtJcQdrRcd2ttiWwa1dKCFRk7W6goSedGsQDqXMhyGbzGJutf9vPZq8wrf6fMlvLLtYCg6ba039bCgjwm+tn6w0EHwCDmFKVOE+F7bDqxtWnWMpOfvMsV8gDhrdqFmylemksqSTG8nh33dKdW0OrlaTHQDZhGQvFa6Gk+Zas98mC1sw6WeJR0rq5y35LEWtpCJZiOMXErqmna4pqslSQNaiWQplKSd1WTe0zQl2sOUYsb6GpnJ2nnPdj9n93P2fspP/LiZsnvAKUsPevV98tgPsXaOuvnG6Gsar4/2YJAf0dRHH/vh2/e31r9L9vBanut6omw9WaTm9nf9ckx3Vwxpt3XEzjHVL3dqXPn6Jo2V85NSf/52dKWpXMOdY5mZoA/paQ4sLRdu92rxbjQ8aCOXv4VqyKYZvVAk4rJruJ6Pi+xxQRzzvxnzPzgfVonbAlKvQVynTjrSkPELDZZPDl9+vVx7pflsfHGO0vXAaZzEJG26CodXcWKvdEI1iBJqMjKYk6VB67kG3OUjqslaXsGMzvqtaz17S4uX3yxxC73VzkXSj8v02XQ5nwQ208+kZAliyi9Wv88/9zcPIura7eJ22N6dt1CvycuNsCpdS9JtXc/2LHOyyxvyRzXMAt0IFsNuN3Sed3jtdHI0DxOxRp5dyB2l6QVVK+E8bvWKez2mn9gK/3dvo/jjP/+v4n//36cxgKdaijsviO1YGmx6UbL576q5mPBKbTXxe+rBrc8DMx3aPAuXx7z4LAAqQdW9JbBPiesTl+vAcN9DUDere1zcI76cPN9yz4UdtwfbWXBjbfZqJ0GQGZlPVGdbon0j6QdZmGqKtuPJB14SgYSEWdwWA0S0baWkrylTWXXFy0QcfiuRjUJN+UWi84sdjJ7ILlbyqLP602jOEIXi2wLiR21J+Fl4jH+UexpkMhKmpr3Crqa4z7Yfq2vO1g20oTF52+a+i/ospBf4EmGepvTKHEbMGiUZCblGVyNNKrBjnSkzvSABaZtlLhTHezpaEq4smviQmZy5QWMyuwXioKqw+X9BGc21mp9DS0O2SHtqkOcBt730jq5Uo+jeqb/juVsTXnFN2y6FL0ENwqs0Qukc8ylMpzKb3qNZ9o4enzpDi5fhLs/WpwO2VioSjhyjehsevD0/HcQ9VopR8ZPfKOlKw0zvrNhbPVXulIdXRpmxBQbFLQ6ij1Dl8Kq/nD0PxYflna4Gouv7q5b3mktVS1gQTb2UpVLVYjFWJi8X7S9R0NRZkv1vvEFHeP+3sjEHGFQJqF0+xEyduZOJl2f+YPW410JX2BP8HXuw1e2tmnNZy4jgK7P9uxfPfnn76smPioofmTz1iqF/MaO8IAVfDHJGIRzeE2ZVd8AFEpaId8cYaHj5989/QuRMEACmS8zE3PID4CTJ0+iPgIYlqoOCaExPkucNyPks3DjCeox7haorxqCiInBVxpHXDn/o/GgZ+M2hEJieERaOdqGO2cYFy9GpFlfTo3SAUFJilsb+0PDLspsuMwZNc0eLa7NE78a8eGGFf5oceUiJgT6AaYbHkDWsCPh2QheuBwTXhDFdjZY+kpsV5JiiYDXANVYXFYK/MRdlpzqb3OZ0PR7vEHwAIvN/OJ0NR/ypI8CzD3wwmh6Fdz+/eQH4j0CNpsvyu/NxXxlygdsiGsxRYHCxutPZOjXLnWTHZUf6s/ct42MwZBE7knzhAZgCgIyIbWIR6zQ9clG+mxANQ4YnbRyezA6Vj71OFeN5V0yrx0jrWQpy8mmqWkmXmdtsxLZTiPicEltIrNy0C4tAWHlJLfILwcOHjBjrRF1o4B3WeNbXCkyFPJMw3qzmKOTCq97rSbvwc2m5eHIdrsPdc50YOFURQNN5bfJuFU8gzUUURWxj1E2sZ7CNBdwbQwtHmKSZIqcvcV3TtkyeNquFGAL1P9E9zJfNaRXPtbISxeVkaWxZHIVwZtvF0XkmpCiXK5MT28T0aHhrwsyrVA+yHS4JwWMbV3P1i7I5gs77eLRMLtCsTjIZo5JAR6uzHjOJohnarDiV3DaifxmTZQi3ykKWuqcBThF4j+JkRcWszgPtNKaYAaAc6q/x6O9Go6XeGt4w1dCW92KArt9OZfIM8BhL4evbxeadhw+hlqhoTeWEI+IAWPUbhKjfiLfy5CJdvdqoDuxT7aJlJ6TSnizaV18pY3M1kHZpxNU+i4DD6kN1g0L5XRftOLlwsb38NGIRpYoQNmzZEFwZc/BKC5lBCW4hjpKODFzf7N0d9w7LMvgCeKUcoP9JcXgCW6bijQit4rIS2Q8n+l//hcPv0XQLIPZZ1GawoAHlC6lqiHEoovHQ1bLnTkG+Zf5McC+ZFw3nDfy0BKCuwPKSwjyJvnC8nv985K4ox/3VobuwFKqa/RIULwXZ4im0GehlE1cDJur8zcbGRqcdHCstqbwHbf6Ps8lUIkvAwNTOoXRKwAOz7SgwFLw8o9TTTpsvRF5aULKjiGfi3bIoL4HDpTcs9t+gbgkUHD3Ost03SxOo5VKPiIvHomTSy7IkIBaVJ4NcHcXLwvJZhLvc5TkPjMBkGk1c4Z0FnGnTUT60O67GboiHQq4hMAW6wkVWJxf0mw1LXTzRW5mJSRBlD8dW5XoBbbXoOZMLYSwUtRXRwH84ryizUo9wOkOYrwnA5v/qIAb1irLJoT1e72itNCZzcLNvOvwgti0WkS2AmBXurtgMheXIwgP/CazRcK1F5KenPWAy3SARMt9TzqiaBklJdAbRfiCpLBQxvQHigVX+HTymCIRF3hjHg+hXV6ohHuYUzK6qZV/ulzGOH7+MwHX5lTi+igXirZHdX3qaV8LNcZsnshVaxiSEazJ2jx8cjo4n058C4bMa5unsVryj5Oi7E3x0Mjkrux7wjv3c5YYIa/hsanH5X1mz/+E/WGWOqAZqme5BwLJeisrtfjeNLeE66DTzi0uW1h8QycdXvXpdvezKlSKrhpWm0mje4RVfH1415nVdetAr5A8USMdedsRPz17tclWGrIRczO1Ca13f5EPMymPZFoMEuu26ksHmXac8BA3rRSFDIJdVSSRKwSBSqgyv2Lwn4k0knOjM3dbC0km6HcMchAQfyYEq08lXJLO8TLph9V7Mwuk6Jrvbj+rIhobsc/WROa9APcQyP1onFsaRNBU2vXrcoppIBEb5IKqlUlSsucBxJj+D22iBv6yD8jmWAtseNzSgmze7DqEvvJD5VMioW9LBLKXKluBlhNXsUD+G6EsZ9+0Vmi45aIfqpquDMtBA/OOcGaQqVmzGI+mRMLmR/+Z7q0r/ZTgLce6NoW25tJ1JHzDRspJd7M8g4S1H9U3LHWh7TUAKQb7X5M+1noStV3YIxBOPX287ERU7EZAU+JJ5dddaKfmubRueXZBxDnMu5nxpMu64hLsYCvVlGr8MWfFoRPfrNiYKtrREAPCOk4uWZFHmijaEJySmt2FP0S7W2K+2svoT0+a+lC2oe3Ojf/9Bz0r1iEOuJKzbMhiBicxEBojlEZpaRQak9izW0Is1vdgHmh2Nxi0sOTEtotcP6oXECdYZLLRQEIIJRBXOGGahJ6YqT9sLiuy1XqQX3Z5DeFGOIOrCxYGYlInTRQQTgzufS2i1S8uk3UW8jtkTgVMkXM1AkvfJYAHJM63mdGWKBgpUSGuchyGTM5fXDp7WUohd5DAHoc1fIkcIHVLU+Oj36ozDbKqLKPpqXQtLOugg+sidWSc9G8yI79qJCaurUNo4KIKhx7uUzjof/1JklDz7TrF5N22/iN3Gowe1tm3sW1I87WLK/2PNUnipqgCAtV9Q5w82Q66sbxRpab14GP4PXvKh3o3wvTQVdnyNbloLDeKdWSBkWQRcBVHwl6TW2NESqb/R8HwNNPC0YN6Pq6T3NRnWXxx3ka2qhQCjAuQKigJk9lEMQVwu6vOBNBCXgRPjFMJtFpMLR1m8pRLaes+4g1u3ECIhI5SgnP1L1gQX4TCxO4VCkWUQrx9lp8hnSl9VMNT96e5JLsMtsZ2GURxXZzIc/MG0BytuGBCeMu7Pru7upPrkLUGxJ/x7Ja5wzNYtYo3SKrNbjKMyv0cvijPFZGcpczSq+VYomeIP83RcwrNxEuFjNEeefGfSkLS57c9j21ncdpKZkQe5pCAjqtYlS5fhidgKa2GMR1aa2HpFwusg3o4g7hSh/Wn/WFc5MjxMaQShIiZDRCBl0MUrjtK1FyrbBFByWAh4vZTc+kEKy1MdXVqbsuf2Jwf7GwehM+jLVWhLNtNG/xGoR15uE3644b33l+bx2jowoWkxWr4M/fo+LLS8unMg6XKuLGssllm8tLngcLLmGQ+l+CRUlv7WjLKsaYNJu/bvqg6vV9zRv5DXbf/AdDNLJBIyzkE/TPqVX5RybOm3/iMyB+8vORfrRbz1797rpR+P/N/WpBsz+2HhRhr3W2dnsk5FruQLG444J+28Us4nKQ3sYfVvFYn6oUUX2das68mz7x4+3epke9MRpNijhy1ce60PTR+zzbSUV7bZ3hlspG3De/c8LD2K8D5dwZ9aAf1Lr7BaQ2mHu8m+U3O7py4A5EUihSXCBV5t6gxFU4F4MyyZPHdNTGpbW4D52wwdLxZrObwx3hu08Tg7lEy6q6kbJNdsd3+MjqcUYDWLXan57PAlc9mhVbzEg2izFmGqWAhQnXnAR78lsXZITuHMuNYT4QO0dAg4Yc17Opu2qg7xLOalZZai05SkCPfEtn6pFPpwtlwGki98pavjZk71BmLMKXrJUaUtMW7hwDxNjuRYZcy4TFijW2liy2TQ3xdrLRMuyp/+AK2Wc68Vcm+SlmSjgElCuRUi5opzPD8+rMqNHv9f/0G3ro/e3Ijnw/25BK+xlX7H3j9q1XO3kB0RzbZaak9UKBzDhz7/KdfKQ+oqlJXA28Y+3L8HRiid6yVDK1VthZMTW1rnt3/axJhmbr244+oz9ronIRiO23ZDVa3gdXJ0qt19h4E0J2dlOFBmUG/zlYV+s817us1XrkVFDOeNpOkXBTkLFyV78kQ10CZ4SU943FXpD3Vwt1/8nPOcPWVLNehQtf59I7yREpnP1RqJVeBJ17oxqYfRnEhvbOTCOQpdmTJaJTPAz10bzgWsyhzFlm7Uc884kl8XsO44rYTDbBA1yYXB264CqLfwzzc6vqjVT52obdziwVTLGmTkSgiRfmFOMeKf4V1p6MSz483PpjGLesFMl2g+UcaoSk8SRmsu8iaP6bLrVMk40NsZ7c3MWNYfZ+usa6FBfCUXCNhI5RCCJGvYAautDrWsNpQhm2w7VzQp8sm13Xif5M4xjUQd8H5Zi4Qcm/NsVL1fqzu9DKiKo6GFFzGOoyTBAmE4vOK/vag7mX1IUlC1NMwXCG+IDptEo5xdvimnS9ie79wxjTu9j3DFsDlh+ylOJswSGpozGgBtqcsFyeyhbEtsQ+LAgkhcS+cljuqkUZTyWu7pue2Zxn5u7k7pR0yq7e7SD94GLTYEtcYq6iE/darlrT49+AL7Oc1eK+ne3HL6vIi7rsVvpT0pZ+G6nS/l/EFqOC4upZfWHfuGG1F8JjiIW+4EuK2fyfiNPTaoRdG0XAnmaOUdhWGHFaBU8yKU2Im00tFuar45atC0cNRE7i0JwqlcGCkSRNIIqBU3aQ8lqkTjG8QByd0C7PzhaHkxEl1lbD0UOJ5XV4EX/QCDNPammVkAdXQ4YpxmxUCE1zg5budeSPgl/Y6jTTzlfsSVsu2OTqT8UckAF7XmcTJTKq+BdhKPE6DEzpwwe1y4NgpPWkaFWZ2a27eYxtorVNY9EBEIrkO1I6e3AAUk7MSbCXKUpfRPEspWVaRs9P39TvQv7RVplNAz7HfQvY4NmfXJKTg4SJoLJlEnrb08iBoMlXjbyD3UE9HQrIIkzpcAdThR+VprxJ9hovpMjliGdamKWddVWB1qqPwoO2xh9wkiankNXLdXa2BvFOyte172Sd0M/0hfcmvIw0b3oyKhXTXWVyBd++ZWunSvv0z38tDAtntVQMVn034kiBpiUy0i54Y0PuYoWvwkpBKHV/EByLvCcQaJryee1ES/aslCEXUAR3BxFG50O2Y6WSC9K6rqebbb/HCyipz5l4FWlmso9ycVlpoXvnwPL2TgPYeZQ/RQmQkBgjg/oSGo55j8+Gngxb+zqNBaPKiJCNFyxNsCtwNoYNtVTg5pUcPXmNOLUcMozRW5m66/w7hCkWeVi0t/5ujQh+m0zvv2Q9oNHGH9EXSUPgTxQrBdCAeG5Xozuyh3+0Q06mMLvF2o7Vv1KeGDbktXrYfTEdx8hDFk2byvDF5PPpz8idO1DK3LuQUhm8NnAipVPDkIB+2KH60ojV6/HttveLi/saVkRg5sHIVs/jH8XYbKUmgiah58fs3fy76xzYL6rbNd0KBMJN+H4h4CKxAxbGfvhusVwZH2W8w8PFI7aVlEWl9/MeWsqyE+tATdVbwDlkl1TSAThJHRVrDO6DOidMFOFUTIW/GJLiP6Gq4EGnH07+JbNAs9tbSuj7etcxPUonBnObuWb4Kc/fqFJ3jq+K2LQDffAewm57JgVTmeTywjk/P7VjTHChSDHMpyIkZ9xs5porSESZAIgzDHVvDp6x+/e/Hq2fdvf3zye/Miu2DCvjcjcVuHthF8klg7l7R8qBO02r2Mwcz0llnFBNqglS4BNM1OD+GhozzGVNN+5m6qM4JMNFLCHo0mJxJeuPjDPHyLpDLIVAoLF9yK30kq6Yv4MT+YhjeoUVx8fb1bWapZ1LC1EeVEp0JrdARwIxdZ0GOZjXsd9sgwBrW+hWa71okuA9HbKpVeLC/iqOo76/r1GMi0fS2auHr7yULI+fm6oHGz5I9Uqr5b2+KLa0EZIAtMU6Zpw5gkZsrQ7qTQFZOgyBV1H9MWbVPkm+fp7cIUSWt2wN7y+b23OGby9+bWWleSj62Bd12LamADUbHwx+Tv6HYuuvdLLdq4Hof11sVhpUChLPLqLZpYuDitt9fFab1txmkdj89a47TeTqZvV8ZpvQWte+vitFZondNQfyu4VhpzqLGJGHZHB/0WTKFNSgku1kKqDJizwyKlrg+CjES/xQglQYjgZ+rS3+HUmHXhN+nKHu4K3lU5gnc3O3GbgVC3J4EFXyzlebhSFgvGX4Wqnk/+y8RCrUI1U/Jk/eK3CFXhO7B4PfHijqh3xr9orplQjKmWstrN5PG9JkgqFHxWfExwXpCRJ0qsdNBl95niNci1BgIkFtzkQ94lTIJEs0p0qsedqIl9dEWNErdifc2Woxc1d9iYfEtftoG6MtGPx3StpWf5aT4KYm84in/85/+ldmG41uR5Kc6Yg0lkR+1XSX4igkvUsms4CI2ITdFTDSGY/5ivQE4H70LQm4sqgpZUzO+CMpqmMOtRHjdXnU3KcEan4XgCwKIeIovut6GspGwOgn3BTJ8y3RKqFqYqiApYvhq02dWSbNjPQXR56KDGquXssNQot3AY799FprsjYaD64ajOn86GoyfLciNBmVoqLsKBfYfQuH3WfgBhGzhpg3pfDNN6TDhByQXLHEw4Wik1bOBNZHjXNWUr6FrrVGdnJ5MjpkG8/V8Ws2mqKQEj/Od+Gdr4p9Bs929v94pO3yrq84tut32z6Qb7FO6NcL7jif943b4FhVixd1/NlnFQ8NkF6VF6olWL7duuNLxYeBcOEWUsj1NCxMQlUjupKRXnGqAbF2vhWiMGs3NCU8kbflkCh5gpWVN+n37xSny8A+uXAK49jQm3qJGghTa0mlTszZ6zRDvJCJfF4sNfnm7sMrvnn0c0VlCJFUeyw4THfg+GB9x9nb8SRYmIOHWCIqqZMzgmhMkd4H+YOI/0LQ6wleTYYVwxp4GZwi6WI8/5gV8LOLrQxGO+e7r7O71FaXKuJYVUiNJV3oKH58ulJr8IBSMyuT62p3m3Xp9RsscWi42xXAsAMxBkI/4yJxElO4LoTO6Hpwc0pU8iULpr/0+kA/kN7nGiWThL3XE2Mlr+lyRCXKjr6Y8igiAYgYyXZar+KCfH+HzNXXoyC5z1biiMyJXj0fLFcnTKT7ZTvw9DofcSguJ98CGVko25LRngJbI0LBUgB8k0BeovinMaQydLQ1VVnXgiQZQlpRPqMRzYm9CbbWGyIkQ+kjjOZv18lDrICI/uxgoRPcGbCwpF28CFO3Zjp7mtreQiTZE12D5Tn2RTns1H4wV3ZfrCbUnt7f9b3Lctt3FlWb73V6QU3QOgDMCkq+x2k0U7ZFnuUliyNKIcahVHwUgRIIkQCKBwkUSrGOHHfu6Ieet57V/omHmcT/GX9Flr731umQkkZddMRZRFnDz3yz777MvaNykOpDxaCyA04I68BkRXsJK5B6pycUEUXGjIDoq78I2P/MCLxXzCsAiL+WqijoXuTimVWbzLSu9Ck82HPuE5BJoWsCHPzQV7AuyHq3ECm+TvA/ThgWu9i9h0U7KgkG2zao2pIVRaBjPkB0OR5w8zH5GH35FuWV0AIivgQ08y6w85UdJkDVXQcdcAq5vxG6HVgXTt/v1WoFjhGyXbRet18/BUKrpfLnibaDFJ5D3ml3wzW20WoN627FpJCCnruumIE0lhhlWdZuGqwnw04wGxNflqqQ03gGnvjnXG9ATotuna/hTPpGTLupl2V5NC6aftjSRYru3NynjksOwakuWKRmVVyaf6pcMXBt59O25XYINX7IjxedtkPytnZ+NpKBLNg83Crt1V28XG3DX9a8ybdU5uj117fvR6ml6GtvrqQxUtvqZEY/Z2PP5cv5tUoKkvx9PFk1k4oYoSjeSOymj3O807wUIqQzoyYKG+VikY8/JdvuShRbVpx4/8Cbyv45Om7nVNHoesCSAN01SZrUq1uzgIH88r4Wl2jiQPBqnDm7PL++ziXmf3MOe3G2JlBYtAjl1LjpI78j3AawHCciUckfAvSCa/2ItFjXLfTWZufl6oY+be8It+YeHTIXF09UvF4Ba+QXXucN/nZ1qnAYLwfE1jt17PwokbPEI3vqIjguUn3P5QNlIDmIMgLx3vsFxfdzsDDmzwrhOIGW/rHXUIUaipRoK66PrJ0BqnVISbuyb1D/Gk7n+xtxdm8KWbGN/RDG9/uHLP4fWWiYWKYFB8sWNa9VJFXfFWbmjsUHLWbXI3YMgL79yxaY4WTArVrw4+DS7z5YlL1KyFFkoWw2a7eYfLGPr53ZhFTLrlPHvr5rCi8YJ+uYf2ohMBwHXd8/8I+Tt+m1/zPgT9meXzPyCG9W+ydO0WxD0Nq6dl+3JIkfRkaKPsdxx2lw+6JBBH0t9702m3c4JnKiNRvOogghrL7KTJGIcrguBVXS2CelzyEFWlGhTQ5JdFjv0+owMXahlenxwP8aoup68s3q3bSdcDt6adXhQciqq6Gfl+iIDlLz6DrUD5vr7A5cQKXE60wE01aAK8ccaIx4Bo8CfaATfN03mHlijSgEu4nACa1bU5GdXdX1msuISF8xEOorhTUo9sgv6WielWEhEhInrpvjsZSzgJHbruwK+11QN9HbKa0cTtqMIim+DSQnI34204rYwv1uaiHrmd5R7KeTclHNRLnN/tTQfPEtk44/Uxa6BFX939RHcga6dHhPkFG3+OCH/FUZTymEAp0e/yffL72J2lwuK6iCSJAfPMn1Pcj6QlnXBr2BX7FnJP/0tDcfFCRlK0KcMnHrym720njMHiPAQ1vA1dzTRZG6B+H8JG4nN9SxRqv/fj/u0MpOUXo1sZVi9sFkpLMKEDFRbt3jap+NIVfWoiTPyAkgoiQ8qMNMq9fJHDEwnUrHlFWP/IHmhp3wkP4x7+1G5YuwhEX6sdPahVbXozRirF7t6jnhKC4b6YDgokzjVxYlR9uVJnZ2hbpwx6EgdqZ1vaX5yUS/XOTektW1N5K9qM0Mjx7KdIWbWewUBQjfpEqMn4Rmp1iJsTu+nHZ4+6kcKhF2DatT2JToBocF4ybbpHMy80daNNa7I/uZ/0e4sFjeeiG69rT1j9uL7bVZTvMmoyW54bl59TnpXND4+avQ7IorSp+0NAnTZHRX+lVGpTAioBVZKnkeZ01TYWU6pjIvzD+IKwEuqUccuOKxI6dRE1VQUZYYf4TlR46V/QfZdTxMxMhiNf+zVNpLVEj8DmwYliaMDt0HZkMfT3UU0lPjRy0m1VQbHjcQ2teupoEJxHHlBbIbHHFpAy0ngZQeXEApFhpYRhM9wzQV9WAsPgaStVZM02CBUu+DKlEB5yqAohsSwerM7KhcA9KhFkcxAl09Cbftpac9eHMqGu3KwNAOIGVc6oKC/KCe1ttG2VPUtn6WB4dqnxXT0tPbCmvQ21OCaor41qNAgV7RV08pQ+hF3V23GUzdyhpkQdvyyn5x5dU604DT3Ccm4WbopjQlzOzi7nyydws1GumFtCGWSZ27W7+MfrcjJdDVW1JHGmmCdmfjCL3xhOQ5czqwj0us+YVAkgGD+zrD8KcN/0qvLeC1Lj/Pzc0awXwT2w0A/yJFErbv8C+3KrmGJgsBD7XxAwemhgKkj2r+XDpAWxt+tkqeowFTTIPv4WC3xVxI3/yYC+PuvVV9KtyS3vRXGo+yJ0Tf178v6xZa0vyZ9YeGIngwhEBgJbHmPpzugX9vsvm9I98a4HZ5eThWc7opcLz7vbrWwtxHNj6MAtzE99dNM14KPWYbtqbwZsxHSHGpdQNY26924yA/qG9kWMUXkfqUOnrzd168xGrOgPfsgMneoBAiRyIn+nXfRVqXpaLTurERb9lLss9p6/s4A9qEX888pmLAodSFa0GrgW8qjmcMv5O4Vkmqxot0lc/ZgczZdvPIBt1DVln5okDVcl38aNJ5qd0qFtO/aRIaI7v59DnBJhr+PYRgF33fj1CFhJD7Iph9h/j4lEJoNj9eomENMO9HioaEXIYoTiSyEW60u35waFz+XPm/bOLsIbuxKobCs9P63313rpbh/YKx6Eq0lN5QV+UgMXkfFWXluszCyWZbR2ciKu8fAWOZgtVsszZ4CGnEo5fPBAcuduZDiCzYfcIoTKBsfLDxU5gsYoZyuJHlp/TE1s19zdW4RAJg/gLrGWPYUKrKZPLaOcNikJJeg0173r1lU6JQtLOgagJvbKwPJoeQPK3rP9ErElv/z8H/1tFqsCP/puXhAMtGtOZMpXmU8XrMfNXfIHt8Qa7khcL8VK29r2rrL+jUQALLUukuaBWw7kU/qNRIbi+vkKMU+Hte59lwC3jfxh1FfUGpcQH4nxOX3z1OvMjeEhRoMBMlqrAvjEhrmrzQUC0yCiWPDygCGuPn5FkiJQc75ZKUMVuJBFnfOVoCrqsY2yCfNnhpdw0ZbQtXBgMQxEsLqNgQasaTFa8KZB7B2DPCEcOl2nzVTLB0BR/lPHuyaJL5cSjJp7abJSbTn8BAi23o/N79Fs4I4DVmQYHfnfYXGXG+auWIW9vg52psvNlE4BCjDiMezROXf7YPAbPs+njme/+9nw87sYkrVNLw0QYsE5pxs/y8X8q1T4lN3ku0USBhoAlaiMTJFNbRJdzt8L+seGMrQu6QSppiS/KEOxqLX85ItQWjK8GoLf8OJlSYxZLBjgviDQpjfbjmJxektJMZV8/vC772ijk0fhdFzdj+Dn77tWGJWr8/SHfxbzqn4ogDChMLhy/w47cSeeIq4ArEL8UZDgbZwbMRLp0wPqG56ho8StyeYscoKUo2sf6AJpNQiwS2I+FC2fiZZtPXybEl99GKwMqEiRMl8dFfv0huGvP3qHpHNzm5IPB8UP5Q+ewYXBy2N3YPA4xEG0ODrcPjyV8u50V6fuWdlxuk2bjphZYvozpuYt6SZFw5gDiOxfl6MaHROF6iEuKX4yUDwTjKVj0YaY8RNH/aYTSBs74Go7pp/hK4Z6maYSmF4NVKqPszuuIUvMbMf2cR7FYD6a8+wZIYn33QFeB0VKamkpG9Afs66q/5EO14J5cqIv5sGUJTuzA8eyxpHPUcFhlMk2aZ4v5MBubf5qxrrx1ztsRPtzMc+miEP7WkwjQ3SjA02Qah33iXIq9Y9LJMa+2616aZZoVo7WwCq+6NRYq+ZijvnCpJ74scpxBljCsZLRiKiCLZ+/EXbQfbODKb/yg8kxjyYrDxI9f3MYnY+MFvfxPeLmw5KcYV9VLR/eAFZ5ulmW6nJFX6qUMnwqPYuEuCSZO/O75/o+TVmxHgQjZbjaqaFHXEXWrUJ+PZUX4tw59A8xQzaD4yCnjlqSv6s8GcGyutrVjtf9WgXiGMUIek00iDuqQfMEl9QWJrNvZGrAcPNQR/Od0lohTKF/6VrVHohsLZMDrt75sirBRT+S5/e8gq/V0naexJgOmf96gEP7nXn+G3ZGFFIhL6CxHmoKLN43LKM82ry4NeMemnWrEUXsJQMWQ6022sqYAfGkSFRsN1UZy47tU71pd+qFzcGQzPFRTQ2qC4zcjbPtSJ9WFB+6i+jKVUt1KD1wuuBIhp4bkRp75FNo8T0yg0/+MEMd1d+mmQ7jyY6FzrW2Bo1sXWRuoJzdLQwOPlTPsxcB3bLFDyoYqtrVzPHglWpxd6oI91Dz17AIZB6pkMxL4XyRa8Au5xuXeq2mmQz7t61a6UOuWErLG5vXtHLRlP2q1z8RuPWGCrLDdXkh5vPI8/CHpz8+75huOmJAPLnrZZzJkGtv6u0mgcPWA/zxgoLAwzbICkyJp4zHVkpTNXPOvG+qus+tbjgRLxEFt/s4pwUf9SZ44lLaQaRbiA3ck5YJEwkVk8bhq7mABRm5siY7ViqvpeEJWbFx1a90Ov18L/OACQKJowR34vBWDAEGVENGSWUzytpIOhuoOix8clLsavXhKwCtljZ6UK37pkh7V0Td8gxQSgCCVqkbrUm/sqkzQ3dBP8aGrX9ggNzlj9givF91rQ6lniDOk58CqBLEe26dFzB5Ky9oVZubJN1WytqwLUmQ4mQ30jVUkIkapFehbY4Lw/qn09WGl1BnHA9zxa1YQa8imiplOB5r3w4h3cQ1kJ+3EFG3I29PQmiIN5TQnQn67IpeegXD06xicQ+KmiQntv9oSar6Te+XqqXJLufC5tnTriSKcnvOKKMjL4p8wy4VWgXI8TzzeB5AM7vvfte8pM6BuK6PE8X8ISj7eVbxal1eQbMSuSZ2aBMonolXdB76tPs/Rh++vBm4/35x01MnxR4RlfXpeGRno8bHvs5ZHkz0KbvHAjLzn6j/PfsU1Dk6L8kALfx89OhPL7Ho6RFPfOzSKU8s85EUmyAD8/AJU3OQCsFgPN+G4pMCqP3unz/aNFu4kcknn/SithrdI7GZ5cDwQSEwPpNXfGZ0CUcwEUALJLgzFXIZzo370Ateq/Hq+vB4FRMo347oONUMKq7hNs6jsgv6CpMhW4FLqWvnW+sNF+XoGHEsuoq5AjMQ+p0eWNC5ekus7a6nSpC8dNlrIBRxx0Jag9aYSZsttdYa1UQSa99xc/pdkIG9JOvQT0L19dMAfcAMDnurRbRCzexVk9IFj3tmn32IQv2uIcBs/aMu1wXzq4/Vl0fqi5vupw2xkjTI138XTTMqxPu5E1gh2ZcNvYkCehken7YkcR7cKfgHSyHujHuYd/U3AWJC7l5WNEK1mRSfNlZy2VAJkGl8nmjHhWEx+qOtrsmd/fH8Cof3w7bjX7cdeaA/5vCV0+mp0mKrs2MHzCb/VkcsGmhdd6odAMLTKX1otSOd4EF+kHUfJOmxfusyoLMU0ds4UORTj1yzdaYZGx5SmcaZcOu+uZqtDuJ9ICBEPkkRiVwtPo76Sd1meOUDzfcDsE4+g5hVEFcFBVAjUTXl08ReCO++5aZQR3qVYWa0R5FqegG/Jp2VIJvMd1rseY9o0ZD2fho89mV1yeYkvvl/a+/8VOr8W3rne1doY+Qa3Ol/O4f6Y/cEWGxzpydPMZm5HQvv4G3cjBk9xIJG7xc/Wdrr1m9H0eI6diWgitQG/vb4V927l/N3dw8SK3KUp9GwRSotVBupcWfQgyi4LF6hGlyWynbFrxB4QBqF2yHnfivPfC8ajcsjpFv/iOh4dbCyx3xAMPq24uXmzwNVa+rrwAzLWz8QIizGILkMWK+9voHpRV+p66x9RQhy+6pIA1/0PfCpBsnofjBTi3MpcerIURT+YqWhMuAXE/9WyhQRli2Put8EZClnr4O71CzlI6k4/PpqPhofiaaQaF477Pd3AQc3mvAHe686U35LfF0ul+G54JMFaO2o8IvgRRphMaKkoATnKsYwxYqgyWRCFLvcNAqbjVer06WAn3GyhvyxCn7idvE2ux30i06KIieERrcrHlYBVI53iHSgIq77jeleCyCRGm3JRRvvrGbJjR301F+iOAgCAmFralVWbXxrfOWRxwyUsDqZq1iQ0QSe1i9ouhFZbMAAYrEWzSGwSADFuxyfj5fjGfBHW+hv0GTQ3/hfjfqbqqm/K2J6FC+gHi0mmFXviJHXq34YBzVNmvNXLO3oCJi7F3f4NrEk2P5I8MzTDuXMFttjIg2Zqc9AbXRbKWiahNyLYEgU1w3irjxQYsL+17/qbzP1eX0d67PaWVA3jGLozv34/ZNzM8v3Bzczor+Djqm1b50NcqWDf6vu1a1sJLtVdwDTUGnv/ibi0x1raabRuWR1sUOmeiczNjddnC/RadiPix0yVok9cYFrcTyauInvJw4pZp/RWXlWyRz8VrH7jFhtQ4Gk0tfYioT6oouxOtSMjByjPTo/tqCL6NN37Ga3I0WEwrIO69HtqvGlIqUiPw2o2tru7oYcUpFSeilJ6KGtJZnDSmbtYnZaDAJOuEIRxEHLzYZeSy38QqOmqp6hxUF7586kouoc3uLOax5P23FYx39Nh2X313b3ar5ZjbereOVgcgTVKnvxfkmwTj7KpNzPl9e11Gl7syabNLxNRLC4Cf4pQileC6z9OiCkgfEGZ8CxqIdQdSQCcYG9JtQh2fyuXtGpEPwR/7UQdCl1KYHhV296FgUKCUCG4JespzRdc6cu2E5a012RMpmAqEf8W2mqBmFyWJy4/vUFwlbwDJnT8aSvoocgRyKdp0blpHM1dpfe2cpveLdDNGnAzKnfptnxQ6ii2U59mmN0WYTS76Z6rJ2z1VtGm5BMz9wg3Bl/ZYAdJ/QGDb3qWAyD5Wn6waoR7EnUoiVRGaqBtrS+nvRLXUWSQ2uiDz6C2q7L7wlgib/ur96+elVFn+BSCMylq5cgFPHMB93I33exq2DLeguaVAP86lrp9mLUVyz8py5V0F4z8E0SKXyMtAihJzspfZ2CsSJBiGvbbQzRCpO0+spuA/b5QWWzjXMQQ5C6BIUgDVCgMcqpVRPEtonwcpv4cjIrmgSLB/UgoYnA8f8VIOj/H2Hj3+zZ3f047PEKkngkq0zd0yroH7dD36hBDomuf34FvnwrI0gCl4h5x0BhTBaEUE3wQnoBscSRrHflsg3JSSr/ZFvlOS0EHdjxWnVZwIe4f3aLIFxjz8vXoC7+6eb+jqZMAeDch7UjHjtsupJ3kX+H7AKmm8096ps+x49qWz60B3eGtiftuL33y8//h7fmLz//747PXGObSI0u8Px8ScLV+ZZQB2Ht6tpWO/qsrEQ/yssnqXUAeAGXbotJbJMQQZw5Yp40Mag/hhr0fE0Xejhvk3SB6WAcCiTjOp+5PXjF5wF+rJfzKS02lvP1/Gw+HUCw2nkVXCdR02TUyx//5m5ue0gEtXdisy1xJqmuvpxmdADfA3RRpbpeU2MGmEAASTVni78Lt2o+KlYaJuNirHS1mLoFOr38qaN1uGMUrOLEI49h1dy7eEq/xrM5ohvYm1na185xzu/LmGoQo2yO28jQksWWtYQKfoqz/i1+dpMpSWpXHI1ehaGYji/KM9ztrHBoCz1cuRJX5anjzz2I/n7fG7e5W7qcjV5wIrtZSZne1QlC0knt0TVUn1fshzLl8402hxvw2Zxa1/t+Rx6vJ1eb6QZ6VRmAO3jf35/yqD0lXvBbaCMtu/l6nq5gHwLBumhWfbq7NmtSfYo1d4rYIus4PanuQaglxtW+UZPx2rFXp94mZeXI6Bh2EvtBzTKZTmWpWUZATrbtr/Tg3nqXvStnYiwZsHLyGs1rkS+8+UZidG3dTJFhVW1lkHZqZYetd7v19EjdgPJdeT4ZT0eMSLN0BJYBWCTHCmEi3HaBXYlPnI3xGH47PtUtRBvX46Fy351n9GkrWCdUNTC0CFVDWDaalBez+Wo9OeOTcWUvaDHgNMMM3bvwVyVitftz/JfNBM5XhZ28zwo/ho7KOteOIw0Ij4FzOvbD5QStTmwaIXQ/Q8zYbFaSGayUlimFwW5SyNfVWMALWbHrP+w4+837/6pcdBE5apK9hPyrhOh/CSl6F+9DbIXk8ZDvE2LBz98Nv90sKfo4xS4s7A6gpxw+p0TCfUeinvxD/ZGSh7o8gaikX4tB9tkKobcT62jcCnoIWpdaP3Gy60cjl1Xy2FGfVpc/eiI00+n8NlDCcJjno1OntO7GsWZfk5PpyLRs+tPNTDO7ux5Ue7MaL0+VN+g0dkSjXeR1PpKDhH0xPtvg5A5kDf3hYRv3NytGErdYkf6jp9QfRWerqrPvHj56cByLIZCO5zp5ra1keNvrn9Kno2p1JDHuwlUNA35WIWFCMAF8l2CjXxVfEt/1sz/oPz0Kv94x6sEDPNC6wfkBAq/VFUKYAeaxnLmyj78ZikWhRKIzFlzkBGyGVngui94JpgI0qQHziDzNit8cMtU9NMcRjL4Xc9IOAVa8sCUHnmPE5hUR7mO8WtFByV/BzCkjzY1mbiLhjU342XRctsFJ/JAP2LwNs5WzCw9emr/9qPOXtL9uJ1eCp3XbjZhtw5oKW2xF2YgMXlkuV6AikiK/Co2OxV/dyk5KAlWUtKDNd6zfU2J2JrKHsiC5Yd0ygK6aDZgIgrbOivSG7e7vO9bXSXeRvL158wuImb9SMV1ZeK+TbuXSw2vrVU+vOJY5eTcEl/sKaj03vu7zh48fnH7/4OWxZHojlgnvTt7A/UVYTFeJ/sYGAmovUg+MWcO3Xo9xu7axQSj0KiL9JvIDTmoq8/PkEiy5/TiwdbvpbefW+sYDHfmJq7AjVVZD+uFTDitMDojJzWFgzuSgHQ83M3edAv6TR1ZcaDKing1ZX55x/9IRhI7Ki5TDkY2uZXv6ptNHtH93BeO6unsE/zseyhFiXoSoeTvBVXvAwy5/96lucNfYQ8v5y8//TiFcNkmy+d3S+pNAC8xffv5fQxHa2UjdzUnAGuFgVU47Kq4BAXlQsO/FZM2wXzOXC6BQl4LqwOuUwaOlrDtlm5lWO4S5b8aq2tG5BflNbBarhCamnDeNzx4DDryluF0I0q1PQsN24dGo0SpkZCJbxRAAcHe4udjFxribQdTJjo8zd+uFiGaVISF/ur3yQjGoyezJ3KQT1RuySbMq7smSB17SNbg6c6+w6cPZes7w9B8KBily86Fhozuc41o3y/oTmT7onmNvu9P0XHaz5nEHRv3u8pdXVqf30qvtwPEQvPB46AZHs+xydJ22rleOa//YXYqbRYi9KUXQD/FR40FcbtIInXW9c8+QzWJbz46H6uxmzENGQOMKk0MhpwCrbK9wS3LPo7BH+ran8lNBRk9m7oBNS9yoX0MforEola2hPbGIVSdnN3Q/RqS5ddXbYiZIWcueYLWLIco0skWFgkfZxr8Xvc/gfDweqWaiRp2EwYKSE5iLzluyJahtgeHfZj1eVqxntxC0RNdkSowWHfEB+Ei61HqoqUPDTvNC72wuD1uX2XdO5xetbYAgZmb+CJiZbbdFQQ/VWJEW8vscy1OFLqQNvWzHptYXiqIzgtiZRvA0qlJ0tKUhEcvtqzEzaZKlk7UyyE0vdEaBAQ05Pu5ijMPI0l6qz+AK7iAqMoA3dZ6W112eSdlR7msVYsXxEZNVoHk02yB8O6oaFqLeNdz+5Rh4nuYqYKkW8/PlfMPBsqz45UOV1Ykbr41p9uuYlYi2jCars7Y6vxY8awO3G1jZNKSB/8FIdZPxiJl6Etcadv2TUspV77DAqd6506rSOoAAlXtsfeY2vVJlo0/LC3CmcBmhv0lkj+hBYrHy1LIRKVaCwoptUFw8hGgslwzdDiPh2Lbnyus0YTsy0LKtUSrsEbjMlG4BDni99EafLpdbA/eP10Khi/nLWG08RMaV5bUAFgiFvuD743OhGv+MhO7neH2RjwAUelciGXt9bkfMrzSIlwH2WmID6KnLycbaYaygd0FLYAXVgpvdlb2kjX/nZlu4oCgh4T7U2dNDkXgudDodnPE/4m23vXuo6T4slmToUm4zW9+iNM2UQnGxqGpf/jHzP8ceC4PwTLoIBNrC5Utukb1Uq1Cb+M0CjLZKxFPlcm4y9jFdfz2ezt8NfhovG1wp6nq9vpxv1t+g4J9dOel/Uk/A+ceOhQypXQicIgAYu8TCkU2DoaOvtnieubY/YkPbUYuVX+m2PiQBYuiBjJDhZYx0vpWLe8U5BDMi5lSA1yUspN/EpA119fVWF5QdVICL3vGQhwHOtITyRfNNGOuP0EIhtsCsGF8t3JTBHkDwmFGT+LuhoV/+/X8yx4SIqMMCdnnSonjRaQG7YwOlxeW8AewbBAFz79onZZVY+248WJ1ZDejIZjaa68q4GZBICqN4vZCpWw2o66qxiAYSU1fGBvSPcukSXTHi0a7EDbIn7ef4lh30EwYlsL4g592V2Xw+n/elF5mtCuVXrtBzLFfXmyL0C18uEt6A9qFIZuDwIflwrCQ/otLuk14Q2gUC5U83SyF63j4h+mo3qdQGW+lLe/lbA0nQroKldgGsybR0NWyNMGy97RVUra5QgwA+kHndUf62yFk9wfHKGtkGX+Wnr9lQuqnwrWY6XbSbHePmi9xPmwCqCqiPZycpwBItZpuXpmNruJsG+/pInbkDfZuC+9HzpRlOsIZAGqpNVM4T8pbhy8B+XbqrjAC4bdsNJdQ0LK7hI3qgbFCLtpMHjTiXShe0Cr2C5ct9wIkHzxu5l58vy9kK0jqRyvxl41hF8edN+6Re/LeNEZRGB0rmozZoUIqzmUXfiUE2Qwt1kXWE5mu1JO0GToxrHuYO88L3kPcFHBlRz8wG/NN8fjUQOr1l6yLX8zm2rex1lip3FPIkisLDeLuzOOND7CqNnAh817VTxqJuTLcpO0gKL8cXIt1s0fNnzNu946Z9KOWiUWBXtqiGm9ctJJSJPsrXYtz05qnj8K2w1KUkOD7Lk6sd1os2INnzgbXSkmZOFF8h5k3b1lTX134stoDwWfc1htram+dWKvwkqtAYntQ/DgyT+s1Btsx4AAfmmnINHCpCZgRHFcfDmOGo1iI8nn/FMpc5yAiygQZcJBwp765VxVUCqN9UQUugzvkS4qp+oX/Q82JRTt29Oqadp/7J9KWuv/4haY5lnynuCn78wL/56RKRiMfvaRhK+1EkPLDfCAJaZ65q2ONjYgl7cmZq89RSlaTMJUHXSetRlVkchAsbjI+Ek9drvMXGtuwgR3+Sua9Vn/uIjS2duduFaIzD7bX3Ztzp8lYNwBg/TgyDoPizI0SqwCslopiGcT8oROditm2JV6d8Wvl9KY8L+mLdXy+neN7cvxoh8DmeHKR/K/cOsh44/v8MkTkcG79eujO9cItfLpcTvAVKSjDHU+HNLNCEVMsG3POEP8DPoQEifdAuhRyNXi8XEw1ZJquN7U+MZS9qWkFp1v3lX/9tNJ6uy5fFpwxzfKCQgKW8iPrJqM/dmCFIEr9ZhLhZQTx5Xi4nlPcL9AmG5moHtiQ9a2UY3k+NvdDYOTGGSACoCPAgfLGsy4scNGRN0a3szkbfu/YMtpnV9grcLmeYXO8t7Pnsbe1sFh/bSh4Oq7kZcNjRma2pQyfEQ3FUquC+qsW81L6iwu/HtAOVBPglfU/4aNOfNb0YdBHNU/8pDpG2k0r3ZIP4IEbaYnDz1iaZ4gfpvT25WR/TsM/N6R4+EIKsfG2O47Kde8Ufi8/3rGWcv3vrLqrDn0Amk258LaVxEgZxcT0NjrDy+2LuuKbh3t7+5/0iyUfbkuF7/Oda9BYw41mt3OE7kHVpsSx6Umj4mK9O03RjKDaTNqw2W0DbCldAm8YIbmsNBvj7+sWOJ1t6+DtztMXJT+erfY8defmYuRFxaiSZulcQFESIJKyYAZYqfvqz4veFYJLxNmCs+n4xGzMWkUT4kvuC3nZ9EHapIKKRh0Lo7iIbeH5hVe86Wt+3Hri74fhycr424k17oXcrQVd2bNHFdOzxBvi+4u2ByiCsemGhiCizJvCrBdzsavhpcVKQNALFjuYbCPiZ0gO6a+kDKQ13zH4L928LX2kUogUlwAphdmnq8/7AaAEjzv1Lv7hOU14C+2c2OigC7x8o1ApTCXrxdXFiOyvgM8lmEDRcTyhD1yVooiMVnwlR5n2NPuEZdqCSj34BQ2L8Wmgce9R7s+v64cyh5Xq1iQ7fZkLVJPizzpfoThhGleLpvDkuH8WH73uGy5hneWlZrnuAtf193lI8V3L9rSeLbidgF6ErrADrsZMQoDpkhCO4L0W3mPALOLc+3t9eFOPPtqT0/IXgY7oFBsVtUUDidAqd6b2qvNhzaxAZx08YB11161uvaaE3/MkNJ/id5WvhPl5nH8MqHGqmvLSmVwrGy+DvMkJg2tPeUvvCH7xHU+/dZfaTJlwj4RoJvbZbWHib7Xs3Os9CcCl6ckzy5vXrxC7Fb3JqA/1+q9cD6t6xbZRsKGlG4jLNeChEMPKN+8WkXi+SGpjAytyQaoU3UXxf9MXHd5yfn1MFEQgyHqE/8R4It4qSW1JjiEtNP3HoyfPc9BKrBd6ugn6mugboSDgt3lpXrpsd1LkRxLw6szqnExXr37Nj2jd16wQRhPYcGVUw6JmBJR9Ek22GlTWyfNW4hpgJmC5eUpwPC5DJ2HWiwpYplCk9c30BvBdC/K17YVrjq0unjKFQEbIsPOKepFM7WWvcvlRDY7em4PXEt+Sle3/R5AA3+847cfR6ukUfrneIXSfxmeUnt1XxrwjVweQe+z+ThXG30h5ILNRkXV+iZ2pwfxdBtcLPuKW46Z/qp25IDpQuCAHbsIuNlxh5LHko1xN/24XJjda0CeWi8TsQ2kLVMyRbkQCsUH0/Oa9+Agzr//1PWmM1ZREsVvy3n7McCZ1tPzVUTGWvsjAtcm7sDs31o7IBBacykjvQmRYHweBeTLhRFgQcVCAXwTXmi3Il0fEoPZvzVJF/pIuef3Zb7dSeeiXlQzUf8QEELJsoJgUT5hwJXxYQYPjvYoE1m14XZ8v5apUywSHGppsi2HNSYSsopBxvqqS8nMzWQYkwwE9MFxwkDOBSLIFVSCbTimwV4EyRqywVmuFxefbXydNL17r7pxx9KhgNs/Lt5KJcz5fuDJVrgAhy89y/Ej8syFYUsfcTlfPYtA6Lb0Gj3S/M7dAbF1Ie0bRFxtTepVuE8iTjsdRGFnsVpJWEY0L3Lncc/vTwh+encIU4DoeII68Ecwf5As4pI32Jf3iaiNmLdLY6mY4eWIa4h9KDTz7JEAaxQFwfjzEYMvfkAcweP3tw71s1v5nVN5IOwgtbvtxzz+7DSIm4dW7rjp/OrU5LXTNtpyoyzLWAHJSyXZWLA/rH4B1J0CS5unF1wUBAH5c8KpXAsmYTcOX6FADzMa/poYD8O5wJykT7cqSP6m4dESa7bE0hvA9reMTK+wHyYomY/Sn/FoD83xlLoCj1lTcFhdQIeC6lFDffF/Nw9l4lj35Wl1Rpy+1em01CADHUREO4EaXm++WCXJ+WkcSHo9ilajNbbRbiXYLL94ZT7m+psC23jqDuwsxHsCr+Gxj3vHq332KZKGx6RPRJqbveSWCub9N4/GxIGPMqHAL5euS5GGOHyLZaSELDvkK1F2A8xm+67jGzZ6/P/IG2v2XDuT8vho+wgS6GL3oMHeCq2mfy+71eL7pF05dC+oBDDF9K1McRf6uOsEf4+rX7jyBoenB5e/BqNrBTnIByfWh/+K/RS8HUuyITj+/O+NkOs/YPoe9+BS5JpRy1Vyb9QxEl2w2XKOtrXSvcYMiR+07vqFfzZ1WbabPrpVUY8WCaJDeiMlgkms8xz6tu70RzvDqxuX3lqMF3k/fjUff3UsxNoYBRmyEtuICa7BqJZdWJGD3ffpXB85/aM3aJufL2g5RfLltW0fbJHVm6hvUMu6lufdt1qo44Nu79w7AJd5HLeDe5l0i8sWreejbxNQ99s3tHcgvCS6rhqsOGiAhsDPAE7SeUtzYFFEjYD3UwcQmwbveJbK7zqtcwm6izTg2SU0ebtSZxxpaZD4tdR3a22HGaGABb5yAhNRfzxPQwJjcNm6ZW4VND6HWjiFTmt9IF7b47YGTr3tx7+y7fr70yUtGRiBNctXohWTvZvQJTP+CEWfuxnidR8/QhWZj5brquXQxBIaUm/LUXiSWrH9HNL1yWpHXXG2k+9B5vnpKvY/Z/EHUcHXAXMjYUBuR+6akt/hha+eAyfHJkv1358pD1aQLoDGVmxVeaRAyQEnhur31v9w/ZjP9hFwS3i5vOaLH82Evb8u/3a2bJ9Tumm5wCq+0r31KYLhisjAf/1IuaDJVHxHSLcds2TVz9YYnFOvae+OjmzdzAZEJi0EJWzts4SG/05+2cLkhYtGRG7oXlqCXWSuQl/HgvWPWU64pZT63rRLikmju9Q+2xpeM5s+QZHu9oAtneUcxFkAPx/Ee/eGuEGPE68O+JDc2LMSNmp47FuT3X4cq8NsP+SjF+KWjjH7L7JYgrRypCarmk8uwvmwkwmeazwiUdFurjW2iICR+i5C3aeRszW3AInhc0BOrsYL1kWmpYrwOBJ20pINuxE1ImStfWDsftIaAN8/HCv4LzOMKI8QoPu/mj+TswHqsx7hlBqjqSkjDGEHM85PZJYoacpkG+xNCKhwYpR1Fn8BwKAQkMmTrXZbhr/UUuhMvtifpu7iGy++Vf/03ld4jzBXECUveaBN8HlLFplE9xKoA7JI6+RD5ETepXQZPU+bpio2TOBQpbnzofUOiVKZ3Yd8ciKlod3fPATw0vy9WTdzMgaY8dWRzitHT//OTJY4J42AYCtAul3l1ZEvx5EuzaaGw3mswjq7uAtZiut2uyF4cH2G5JIinlFKri6quX0Rl8X098Vx1JWags35bOrTQcEHmLupq8Zevr68SygZqmWNbexCvVPes4NXVq4MyiSA35CfeW6sZNUxZPkJpIPpUu6OswV9zV6MwaOn5TeDcX9WKpGD1gz6WKGB3LDs8EcZTzCpLm7sfqpw2VI4293dHyPbhIPnJspfhZRC4FH+NhwdqeQeqVVvdxtRU6IZQQ+79w3h7PHRmSJGrRGG3Yh55z6SBkTDTcZjPh91bQO+Kfb1cQe3h7QXLV7ISZkyLPmOHJa0dS3o4Fddl1M03tauHhXBKqVLZ3WFdMwcpjm9WkDjMABV0x3xL8/3wzEy9yg581Rq+ZnctUxrmdLG6zd5OlQkXQjqvrQaodgXw9HV8FeINtbtiwgyVuXMVJXKuJAADYmM6zI7Znkw0p0R9owOIILjpPtX+/QGf7xTeb835xH7ghfdxW/eL4+ctHD/pi8tA3aiXRr+4x8JWPqNaP4tv2o3BhfTHn/KZcPtJYibr7vkcYIpxgMXAX69/74pUsESKsAvn1/Xixtr/dOkPn0I8iliRRR9x2gx+Jmx4T7q2+Q+wW767WV+2iRRwhStiB8iBRUtfNyFAx/SRE2pu+54r7wbm4b8xwP8Hy6EfB5frGpFlZRdE9FmkyeqLhph6NEWiP0LDjmQa77Wew7/3AIvs/v2W6Tq/EtuwXo/HZ5KqcWuRNBIZeEU1+OTk7lrgJkQtw7M9bcXBN4qEEa223Y548+vHxD6cvnjz71u2q4x8fP7737OXpsycvjmEHJTEC1fJMuYcCera/u+nhPPwXUEsBAhQDFAAAAAgAAAAhXAncLEPsKwAA/bUAAAkAAAAAAAAAAAAAAKSBAAAAAHN0ZXA0LmNzc1BLAQIUAxQAAAAIAAAAIVwXAir1qh0AAMZ3AAAKAAAAAAAAAAAAAACkgRMsAABzdGVwNC5odG1sUEsBAhQDFAAAAAgAAAAhXGwZMCYPiwEA8SYFAAgAAAAAAAAAAAAAAKSB5UkAAHN0ZXA0LmpzUEsFBgAAAAADAAMApQAAABrVAQAAAA=="]
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

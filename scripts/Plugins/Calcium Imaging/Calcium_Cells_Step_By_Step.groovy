// Calcium imaging: cells step by step                                                    version 2026.11.20.22
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
        say("Calcium imaging plugin, version 2026.11.20.22")
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
    static final String VERSION = "2026.11.20.22"
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
    static final String SHA256 = "54e411a74a2374a32c9037afa27c1a5055628f28405b84024a1a90ec11ff124d"
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
        boolean respondersNext = false
        try { respondersNext = respondersWaiting.call() && all.length > 0 && nLikely == 0 && !detecting; addMissedBtn.setVisible(respondersWaiting.call() && measuredNow) } catch (Throwable ignored) { }
        // one blue button: the next thing to press (cosmetic, so a window without these buttons carries on)
        try {
            String blue = detecting ? '' : respondersNext ? 'responders' : measuredNow ? '' : (proposals.isEmpty() && all.length == 0) || (!steps.blind && proposals.any { it.blind }) ? 'suggest' : nLikely > 0 ? 'accept' : all.length > 0 ? 'continue' : 'suggest'
            CalciumTheme.primary(suggestBtn, blue == 'suggest'); CalciumTheme.primary(acceptAllBtn, blue == 'accept'); CalciumTheme.primary(continueBtn, blue == 'continue'); CalciumTheme.primary(addRespondersBtn, blue == 'responders')
        } catch (Throwable ignored) { }
        if (detecting) sayNext("wait while the cell outlines are made (see the bar below).")
        else if (measuredNow && respondersNext) sayNext("these " + cells.length + " cells are measured. Cells that only show when they respond have no outline yet: press <b>Add cells that respond…</b> and measure again, or go on to <b>Analyze responses</b>." + hiddenNote)
        else if (measuredNow) sayNext("these " + cells.length + " cells are measured. Go on to <b>Analyze responses</b>, or change the cells here and measure again." + hiddenNote)
        else if (proposals.isEmpty() && all.length == 0) sayNext("press <b>Suggest cells automatically</b>, or draw the cells yourself." + hiddenNote)
        else if (!steps.blind && proposals.any { it.blind }) sayNext("press <b>Suggest cells again</b>: the outlines on the picture were made before the response, not with it." + hiddenNote)
        else if (nLikely > 0) sayNext("press <b>" + (withResponse ? "Accept responding and likely" : "Accept likely") + " (" + nLikely + ")</b>." + hiddenNote)
        else if (steps.blind && respondersNext) sayNext((nUnsure > 0 ? "click the dashed outlines that are real cells (" + nUnsure + " left). Then press " : "press ") + "<b>Add cells that respond…</b>: cells that only show when they respond have no outline yet. Then <b>Continue to measure</b>." + hiddenNote)
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
    boolean explained = respondersExplained; respondersExplained = false   // the question before measuring has already said it
    if (steps.blind && !testGui && !explained && javax.swing.JOptionPane.showConfirmDialog(frame, "<html><body style='width:320px'>The cells you have now stay recorded as <b>chosen blind</b>, as long as you do not change them.<br><br>The response is shown next and the plugin looks for the cells that only show when they respond. Cells you add from now on are recorded as <b>added after the response was shown</b>. Step 4 gives % responding for both groups.</body></html>", "Add cells that respond", javax.swing.JOptionPane.OK_CANCEL_OPTION) != javax.swing.JOptionPane.OK_OPTION) return
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
    int ans = javax.swing.JOptionPane.showConfirmDialog(frame, "Accept all " + nL + " outlines graded " + (steps.blind ? "" : "\"responds\" or ") + "\"likely a cell\"?\n\nIn tests about 7 in 10 of these were real cells, so about 3 in 10 are not.\nLook through them afterwards and remove the wrong ones.", "Accept all likely ones", javax.swing.JOptionPane.OK_CANCEL_OPTION)
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
                sensitivityHtml = m != null ? "Sensitivity check: the blind detection also found <b>" + m.found + " of " + m.cellLike + "</b> response-assisted outlines that could be cells (" + CalciumSensitivity.pct((int) m.found, (int) m.cellLike) + "); " + (m.cellLike - m.found) + " missed" + ((m.cellLike - m.found) > 0 ? ": <b>Add cells that respond…</b> below adds them" : "") + ". Details in what_was_done.txt.<br>" :
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
if (testGui) { binding.setVariable("hooks", [nextLine: nextLine, addRespondersBtn: addRespondersBtn, addMissedBtn: addMissedBtn, continueBtn: continueBtn, acceptAllBtn: acceptAllBtn, autoSuggestAsked: { -> autoSuggestAsked }, legendClicked: { String k, boolean o -> legendClicked.call(k, o) }, hiddenKinds: hiddenKinds, cellRois: cellRois, sensBox: sensBox, basisFor: basisFor, responseSeen: { -> responseSeen }, reveal: { -> reveal() }, runSensitivity: runSensitivity, trace: trace, respondsLegend: respondsLegend, viewAllBtn: viewAllBtn, viewBaselineBtn: viewBaselineBtn, viewResponseBtn: viewResponseBtn, selectionNote: selectionNote, selectionBasis: { -> selectionBasis }, changeBtn: changeBtn, changeRecording: changeRecording, continueBtn: continueBtn, backToCellsBtn: backToCellsBtn, backToMeasureBtn: backToMeasureBtn, cards: [s1, s2, s3, s4], recordingNote: recordingNote, detectionProgress: detectionProgress, suggestBtn: suggestBtn, edgeBox: edgeBox, addBtn: addBtn, setForceFresh: { boolean v -> forceFresh = v }, setupCellpose: setupCellpose, findCellposePython: findCellposePython, setupBtn: setupBtn, setupRow: setupRow, cpLabel: cpLabel, showCellposeState: showCellposeState, setupText: SETUP_TEXT, saveMovie: saveMovie, makeMovie: makeMovie, rm: rm, refresh: refresh, go: go, status: status, steps: steps, outDir: outDir, showBg: showBg, frame: frame, suggest: suggest, acceptAt: acceptAt, acceptAll: acceptAll, hoverAt: hoverAt, proposals: proposals, result: result, info: info, removeSelected: removeSelected, undoLast: undoLast, root: root]); refresh(true); return }
frame.setVisible(true)
timer.start()
refresh(true)

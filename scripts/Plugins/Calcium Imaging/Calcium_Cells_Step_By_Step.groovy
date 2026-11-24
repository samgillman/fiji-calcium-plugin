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

    // ---- the recording can be closed in Fiji while the plugin's window stays open ----
    /** Said wherever an action needs the recording's frames and they are gone. */
    static final String CLOSED = "The recording is closed in Fiji. Open it again to go on."
    /** Closing a recording's window in Fiji empties the picture: its frames are gone, and reading one gives
     *  "Stack argument out of range". True when that has happened. Looked at without touching the picture. */
    static boolean closed(ImagePlus recording) { return recording == null || !recording.isProcessor() }
    /** What tells this recording from any other: width, height, frame count and bit depth, its frame interval and pixel size,
     *  and the pixels of its first, middle and last frame. Taken while the recording is open; used to know it again when it has
     *  been closed and opened again. A copy with another time base or pixel size is another recording. */
    String recordingPrint = ""
    static String printOf(ImagePlus recording) {
        ImageStack st = recording.getStack(); int frames = st.getSize()
        ij.measure.Calibration cal = recording.getCalibration()
        StringBuilder print = new StringBuilder(sizeOf(recording)).append(':')
            .append([cal.frameInterval, cal.getTimeUnit(), cal.pixelWidth, cal.pixelHeight, cal.getUnit()].join('/').replaceAll('[\\n:]', ' '))
        for (int f : [1, (int) ((frames + 1) / 2), frames]) {
            Object px = st.getProcessor(f).getPixels()
            int hash = px instanceof short[] ? Arrays.hashCode((short[]) px) : px instanceof byte[] ? Arrays.hashCode((byte[]) px) : px instanceof float[] ? Arrays.hashCode((float[]) px) : Arrays.hashCode((int[]) px)
            print.append(':').append(hash)
        }
        return print.toString()
    }
    static String sizeOf(ImagePlus recording) { return recording.getWidth() + "x" + recording.getHeight() + "x" + recording.getStackSize() + "x" + recording.getBitDepth() }
    /** The prints of the open pictures already compared, so that a picture is read once, not every time the window looks. */
    static final Map<ImagePlus, String> PRINTS = Collections.synchronizedMap(new WeakHashMap<ImagePlus, String>())
    /** The recording with this print among the pictures open in Fiji, or null. */
    static ImagePlus openWithPrint(String print) {
        if (print == null || print.isEmpty()) return null
        int[] ids = WindowManager.getIDList()
        if (ids == null) return null
        for (int id : ids) {
            ImagePlus open = WindowManager.getImage(id)
            if (open == null || closed(open) || !print.startsWith(sizeOf(open) + ":")) continue
            try {
                String seen = PRINTS.get(open)
                if (seen == null) { seen = printOf(open); PRINTS.put(open, seen) }
                if (seen == print) return open
            } catch (Throwable ignored) { }
        }
        return null
    }
    /** The recording to read frames from: this one while it holds its frames, else the same recording opened again in Fiji
     *  (the same print), else null. The one check for the plugin window, the Step 4 page and the Fiji Step 4 window. */
    static ImagePlus live(ImagePlus recording, String print) { return closed(recording) ? openWithPrint(print) : recording }
    /** True when the recording holds its frames. When it was closed and the same recording is open in Fiji again (same size,
     *  frame count, time base and pixels), the plugin goes on with that one and says so in the record. */
    boolean recordingReady() {
        ImagePlus open = live(movie, recordingPrint)
        if (open == null) return false
        if (!open.is(movie)) {
            movie = open
            say("Recording: it was closed in Fiji and opened again. The plugin went on with the reopened recording after comparing its size, frame count, time base, pixel size and first, middle and last frame.")
        }
        return true
    }

    /** Steps 1-2: pictures. The movie is only read. */
    void prepare(ImagePlus imp, boolean lookForResponse) {
        movie = imp; n = imp.getStackSize(); w = imp.getWidth(); h = imp.getHeight()
        try { recordingPrint = printOf(imp) } catch (Throwable ignored) { recordingPrint = "" }
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
        String id=random?null:shapeId(roi)
        if(id==null)id='Cell_'+UUID.randomUUID().toString().replace('-','')
        roi.setProperty('calcium.roi.id',id);id
    }
    /** The ID worked out from where the outline is and which pixels it covers (null when it cannot be). Writes nothing. */
    static String shapeId(Roi roi) {
        try{
            java.security.MessageDigest md=java.security.MessageDigest.getInstance('SHA-256')
            java.awt.Rectangle b=roi.getBounds();md.update((b.x+','+b.y+','+b.width+','+b.height+';').getBytes('UTF-8'))
            def mask=roi.getMask();if(mask!=null)md.update((byte[])mask.getPixels())
            return 'Cell_'+md.digest().encodeHex().toString().substring(0,32)
        }catch(Exception ignored){return null}
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
    String recordingPrint=''   // CalciumSteps.printOf of the attached recording
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
        // A matching original recording remains read-only while movie frames are displayed. One closed in Fiji now is attached
        // by its print: its movie views say it is closed, and go on with the same recording once it is open again.
        ImagePlus given=data.movie instanceof ImagePlus?data.movie as ImagePlus:null;String print=data.recordingPrint instanceof String?data.recordingPrint as String:''
        ImagePlus open=CalciumSteps.live(given,print)
        if(open!=null&&snapshot!=null&&open.width==snapshot.width&&open.height==snapshot.height&&open.stackSize>=(data.frames as int)) {
            movie=open;view.addItem('Movie: change from rest');view.addItem('Movie: as recorded')
            if(print)recordingPrint=print else try{recordingPrint=CalciumSteps.printOf(movie)}catch(Throwable ignored){recordingPrint=''}
        } else if(given!=null&&snapshot!=null&&print&&CalciumSteps.closed(given)) {
            movie=given;recordingPrint=print;view.addItem('Movie: change from rest');view.addItem('Movie: as recorded')
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
        // The recording's window was closed in Fiji: its frames are gone. The same recording opened again is taken instead.
        ImagePlus open=movie!=null&&mode>0?CalciumSteps.live(movie,recordingPrint):null
        if(open!=null)movie=open
        if(movie!=null&&mode>0&&open==null) {
            if(player!=null)stopPlaying()
            cellPicture.image=snapshot;cellPicture.onGray=false;cellPicture.note=CalciumSteps.CLOSED
        } else if(movie!=null&&mode>0) {
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
    /** The row of the cells whose outlines were kept for this recording without being chosen in the plugin on it: from another
     *  recording, or drawn on this one before the plugin started. */
    static final String CARRIED_ROW='Carried over'
    /** % responding to the primary stimulus, worked out in this one place, saved beside the tables and shown the same way on
     *  every tab. Rows: all cells; the cells chosen blind and those added after the response was shown, when that was
     *  recorded; and the quiet cells (required rise at most QUIET_RISE), where noisy cells are set aside whether they responded
     *  or not, so the row is not tilted toward responders. Each row is given twice: with every cell, and without the cells whose
     *  background-corrected fluorescence went below zero (dF/F0 under -1 at some frame).
     *  Cells whose outlines were carried over (measurement.carriedOver: from another recording, or drawn on this one before the
     *  plugin started) were not chosen in the plugin on this recording: they are in neither of the two groups and have a row of
     *  their own, which exists only when there are such cells.
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
            Map carried=measurement?.carriedOver instanceof Map&&!(measurement.carriedOver as Map).isEmpty()?measurement.carriedOver as Map:null
            def blank={[counted:[0,0,0],without:[0,0,0]]}   // responded, evaluable, registered
            Map rows=[all:blank(),blind:blank(),later:blank(),carried:blank(),quiet:blank()];List under=[]
            calls.drop(1).findAll{role<0||it[role]=='Stimulus'}.each{row->
                String name=''+row[cell],call=''+row[status];boolean responded=call=='Responded',evaluable=responded||call=='Did not respond'
                if(tags!=null&&!tags.containsKey(name))return   // a cell that is not in this measurement
                boolean below=evaluable&&low.contains(name)     // an excluded cell is out of the percentages already
                if(below)under<<name
                double rise=need<0?Double.NaN:CalciumCsv.number(''+row[need])
                List into=['all']
                if(carried!=null&&carried.containsKey(name))into<<'carried' else if(tags!=null)into<<(tags[name]?'blind':'later')
                if(evaluable&&Double.isFinite(rise)&&rise<=QUIET_RISE)into<<'quiet'
                into.each{String key->List both=[rows[key].counted];if(!below)both<<rows[key].without   // (a "?[" would read as a safe index)
                    both.each{List t->t[2]++;if(evaluable)t[1]++;if(responded)t[0]++}}
            }
            List out=[[key:'all',label:'All cells']+rows.all]
            if(tags!=null)out+=[[key:'blind',label:'Chosen blind']+rows.blind,[key:'later',label:'Added after the response was shown']+rows.later]
            if(carried!=null)out<<([key:'carried',label:CARRIED_ROW]+rows.carried)
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
            def of={String key->(rows.find{it.key==key} as Map)?.counted as List}
            def line={String label,List t->label+': '+t[0]+' of '+t[1]+' responded'+((t[1] as int)>0?' ('+percentOf(t)+' %)':'')}
            List carried=of('carried')
            if(carried!=null){   // outlines carried over: never among the cells chosen blind
                List parts=[line('All cells',of('all'))]
                if(found.tagged){parts<<line('Chosen blind',of('blind'));if((of('later')[2] as int)>0)parts<<line('Added after the response was shown',of('later'))}
                parts<<line(CARRIED_ROW,carried)
                return ' '+parts.join('. ')+'.'
            }
            if(!found.tagged)return ''
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
    /** What tells the attached recording from any other (CalciumSteps.printOf), to know it again when it was closed in Fiji and opened again. */
    String recordingPrint=''
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
        // One closed in Fiji now is attached by its print: its frames are refused in plain words until it is open again.
        ImagePlus given=data.movie instanceof ImagePlus?data.movie as ImagePlus:null;String print=data.recordingPrint instanceof String?data.recordingPrint as String:''
        ImagePlus open=CalciumSteps.live(given,print)
        if(open!=null&&picture!=null&&open.width==width&&open.height==height&&open.stackSize>=(data.frames as int)){movie=open;if(print)recordingPrint=print else try{recordingPrint=CalciumSteps.printOf(movie)}catch(Throwable ignored){recordingPrint=''}}
        else if(given!=null&&picture!=null&&print&&CalciumSteps.closed(given)){movie=given;recordingPrint=print}
        if(data.motion instanceof Map)motion=data.motion as Map
        else if(movie!=null&&!CalciumSteps.closed(movie))try{Map chosen=data.fijiFrames instanceof Map?data.fijiFrames as Map:null;motion=CalciumMotion.summary(CalciumMotion.shifts(movie.getStack(),data.frames as int),chosen?chosen.first as int:0,chosen?chosen.last as int:0)}catch(Throwable ignored){motion=null}
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
        // The recording's window was closed in Fiji: its frames are gone. The same recording opened again is taken instead.
        ImagePlus open=CalciumSteps.live(movie,recordingPrint);if(open==null)throw new FileNotFoundException(CalciumSteps.CLOSED);movie=open
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
                        String version=now.frameVersion();byte[] image
                        // the recording was closed in Fiji: the page says so in plain words instead of a frame
                        try{image=now.frame(frame as int,query.mode,size as int)}catch(FileNotFoundException gone){if(gone.message!=CalciumSteps.CLOSED)throw gone;replyJson(exchange,410,[error:gone.message]);return}
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
    static final String SHA256 = "8f883342b88f94744630bfb83fb04881a5ad0401bd0dcb268d834ea10306a254"
    static final List<String> PARTS = [
        "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"]
    static final String CONTENT = PARTS.join("")
}
// END GENERATED WEB BUNDLE

/** One recording, its own outlines. The ROI Manager is one list shared by everything in Fiji: outlines chosen on one recording
 *  are still in it when the plugin starts on the next. The plugin keeps two notes with that list: for each outline it added,
 *  the recording it was added on; and for each outline kept from elsewhere, the recording it was kept for. An outline is named
 *  by its ID (CalciumMeasurement.idFor). An outline with no note was not added by the plugin: it is from before. Nothing is
 *  written into the outlines, so the files saved with a measurement are what they always were. The notes of where outlines
 *  were added last as long as the list: an outline cleared and later opened again from a saved RoiSet is still known. Closing
 *  the ROI Manager's window ends the list; the next ROI Manager is a new list with no notes.
 *  A recording is known by its print (CalciumSteps.printOf), not by its name: two recordings of the same name are two. For its
 *  outlines the time base and pixel size do not count (same): setting them in Image > Properties keeps the outlines its own.
 *  Every start of the plugin runs in a class loader of its own, so the notes are kept outside its classes (Swing's table of
 *  defaults, by list, dropped with the list) and hold nothing but text. */
class CalciumOutlines {
    static final String MADE='made on',KEPT='kept for'
    /** Kept at Change recording, for whichever recording is chosen next. */
    static final String NEXT='\u0000next'
    /** Put before a recording's mark: drawn by hand on that recording's own window before the plugin started on it. */
    static final String BEFORE='\u0000before\u0000'
    /** The note of an outline cleared at the start whose recording was not known: if Undo brings it back in the window that
     *  cleared it, it is asked about again there. Anywhere else it counts as no note. */
    static final String UNKNOWN='\n\n'
    /** Where a kept outline came from, when it was drawn on this same recording before the plugin started. */
    static final String HERE='\u0000here'
    /** The name under which the open plugin window keeps a note on itself, so that a second start finds it. */
    static final String WINDOW='calcium.steps.window'
    /** The name under which a start of the plugin notes that it is under way, from Step 1 until it ends. */
    static final String STARTING='calcium.steps.starting'
    static final String STARTING_TEXT='The plugin is already starting. Finish or cancel that start first.'
    /** The buttons of the question about n outlines. */
    static List<String> choices(int n){n==1 ? ['Clear it','Keep it for this recording','Cancel'] : ['Clear them','Keep them for this recording','Cancel']}
    static final String KEEP_ONLY_IF='Keep outlines only if this recording shows the same cells in the same place. Kept outlines are recorded as carried over, never as chosen blind.'

    static Map<String,String> notes(Object list,String which) {
        synchronized(javax.swing.UIManager) {
            def all=javax.swing.UIManager.get('calcium.outlines')
            if(!(all instanceof Map)){all=new WeakHashMap();javax.swing.UIManager.put('calcium.outlines',all)}
            Map ofList=(Map)all.get(list);if(ofList==null){ofList=new HashMap();all.put(list,ofList)}
            Map held=(Map)ofList.get(which);if(held==null){held=new java.util.concurrent.ConcurrentHashMap();ofList.put(which,held)}
            return held
        }
    }
    /** One start at a time: false when another start is still under way (making its pictures, or waiting at one of its
     *  questions) and has no window yet. A start that ends, also by a return or an error, lets go with ended(). */
    static boolean starting() {
        synchronized(javax.swing.UIManager) {
            def other=javax.swing.UIManager.get(STARTING)
            if(other instanceof Thread&&!other.is(Thread.currentThread())&&((Thread)other).isAlive())return false
            javax.swing.UIManager.put(STARTING,Thread.currentThread());return true
        }
    }
    static void ended(){synchronized(javax.swing.UIManager){def now=javax.swing.UIManager.get(STARTING);if(now!=null&&now.is(Thread.currentThread()))javax.swing.UIManager.put(STARTING,null)}}

    /** What names a recording: its title, its file when Fiji knows it, and its print. */
    static String mark(ImagePlus recording,String print) { recording.getTitle()+'\n'+fileOf(recording)+'\n'+(print?:'') }
    /** The file a recording was opened from, or '' when Fiji does not know it or it is no longer there. */
    static String fileOf(ImagePlus recording) {
        try{def info=recording.getOriginalFileInfo();if(info!=null&&info.directory&&info.fileName){File file=new File(''+info.directory,''+info.fileName);if(file.isFile())return file.path}}catch(Throwable ignored){}
        ''
    }
    static String nameOf(String mark){mark==null?'':mark.split('\n',-1)[0]}
    static String printPart(String mark){List parts=mark.split('\n',-1) as List;parts.size()==3?parts[2]:''}
    /** A recording's print without its time base and pixel size: size, frame count, bit depth and the first, middle and last frame. */
    static String pixels(String print){List p=print.split(':',-1) as List;p.size()==5?([p[0]]+p.drop(2)).join(':'):print}
    /** The same recording for its outlines: the same size, frame count, bit depth and first, middle and last frame, whatever its
     *  name, file, time base or pixel size. (Taking up a reopened recording also needs the same time base and pixel size.) */
    static boolean same(String a,String b) {
        if(a==null||b==null||a==NEXT||b==NEXT||a.startsWith(BEFORE)||b.startsWith(BEFORE))return false
        String x=printPart(a);x.contains(':')&&pixels(x)==pixels(printPart(b))
    }
    /** Where an outline came from, seen from a recording: the name of the recording it was added on, HERE when it was drawn on
     *  this same recording before the plugin started, or '' when that is not known. */
    static String origin(String made,String recording) {
        if(made==null)return ''
        if(made.startsWith(BEFORE)){String on=made.substring(BEFORE.length());return same(on,recording)?HERE:nameOf(on)}
        nameOf(made)
    }
    /** An outline's ID as the plugin gives it, worked out without writing anything into the outline. */
    static String idOf(Roi outline) {
        String id=outline.getProperty('calcium.roi.id')
        (id&&id==~/Cell_[a-f0-9]{32}/)?id:(CalciumMeasurement.shapeId(outline)?:'')
    }
    /** An outline the plugin added on a recording. addCell gives it an ID that no other outline in the list has. */
    static void made(Object list,Roi outline,String recording){String id=CalciumMeasurement.idFor(outline);notes(list,KEPT).remove(id);notes(list,MADE)[id]=recording}
    /** Outlines drawn on this window's picture and added with the ROI Manager's own Add (the t key) are this recording's too.
     *  One that is new in the list since the last look (seen: the IDs of that look) is noted as made here, also when an outline
     *  of the same shape was noted on another recording before. Outlines that were already in the list keep their notes. */
    static void claim(Object list,Roi[] outlines,String recording,int pictureId,Set<String> seen) {
        Map made=notes(list,MADE),kept=notes(list,KEPT)
        outlines.each{Roi outline->if(outline!=null&&outline.isArea()&&pictureId!=0&&outline.getImageID()==pictureId){String id=CalciumMeasurement.idFor(outline)
            if(!seen.contains(id)){made[id]=recording;kept.remove(id)}else if(!made.containsKey(id)&&!kept.containsKey(id))made[id]=recording}}
        seen.clear();outlines.each{Roi outline->if(outline!=null)seen.add(CalciumMeasurement.idFor(outline))}
    }
    /** At the start on a recording, before anything is changed: what each outline the list holds is. 'own': added by the plugin
     *  on this recording; 'carried': kept for it, at Change recording or by an earlier answer; 'here': drawn on this recording's
     *  own window (recordingId) before the plugin started; 'ask': from another recording, or from before. Writes nothing. */
    static Map sort(Object list,Roi[] held,String recording,int recordingId) {
        Map made=notes(list,MADE),kept=notes(list,KEPT);List kinds=[];Map by=[own:[],carried:[],here:[],ask:[]]
        held.each{Roi outline->
            String kind=null
            if(outline!=null&&outline.isArea()){
                String id=idOf(outline),on=made[id],forWhich=kept[id]
                if(on==UNKNOWN)on=null   // cleared at another start, and not known: as if never noted
                // Only the ROI Manager's own Add takes an outline from the recording's own window (the plugin adds from its colour
                // picture), so such an outline was drawn there, also when an outline of the same shape was noted elsewhere before.
                boolean here=(recordingId!=0&&outline.getImageID()==recordingId)||(on!=null&&origin(on,recording)==HERE)
                kind=same(on,recording)?'own':(forWhich==NEXT||same(forWhich,recording))?'carried':here?'here':'ask'
                (by[kind] as List)<<outline
            }
            kinds<<kind
        }
        by+[kinds:kinds]
    }
    /** The start goes on: the notes say what was decided. IDs are written into the outlines as the plugin window does; a second,
     *  identical outline gets an ID of its own with the same notes. With clear, the outlines asked about keep a note of where
     *  they came from, so that they are asked about again if Undo brings them back, and their IDs are returned for removing.
     *  Notes of outlines kept that are no longer in the list are dropped. */
    static Set<String> settle(RoiManager list,Roi[] held,Map sorted,String recording,boolean clear) {
        Map made=notes(list,MADE),kept=notes(list,KEPT);Set<String> keep=new HashSet<String>(),cleared=new LinkedHashSet<String>()
        Roi[] now=list.getRoisAsArray(),unique=now.collect{it} as Roi[]
        Map<Roi,String> was=new IdentityHashMap<Roi,String>();now.each{Roi r->if(r!=null)was.put(r,idOf(r))}   // the IDs the notes were read by
        CalciumMeasurement.ensureIds(unique)
        for(int i=0;i<now.length;i++)if(!now[i].is(unique[i]))list.setRoi(unique[i],i)
        held.eachWithIndex{Roi outline,int h->
            String kind=(sorted.kinds as List)[h];int i=now.findIndexOf{it.is(outline)}
            if(kind==null||i<0)return
            String id=CalciumMeasurement.idFor(unique[i]),on=made[was.get(outline)]
            if(kind=='own'){made[id]=recording;return}
            if(kind=='ask'&&clear){made[id]=on?:UNKNOWN;cleared<<id;return}
            if(kind=='here')on=BEFORE+recording
            if(on!=null)made[id]=on
            kept[id]=recording;keep<<id
        }
        kept.keySet().retainAll(keep)
        cleared
    }
    static void keep(Object list,Collection<Roi> outlines,String recording){outlines.each{Roi outline->notes(list,KEPT)[CalciumMeasurement.idFor(outline)]=recording}}
    /** Change recording with the outlines kept: they go to whichever recording is chosen next, without a second question. */
    static void keepForNext(Object list,Roi[] outlines,String recording) {
        Map made=notes(list,MADE),kept=notes(list,KEPT)
        outlines.each{Roi outline->if(outline==null||!outline.isArea())return;String id=CalciumMeasurement.idFor(outline);if(!made.containsKey(id)&&!kept.containsKey(id))made[id]=recording;kept[id]=NEXT}
    }
    /** The list was emptied: no outline is kept for any recording. Where outlines were added stays noted. */
    static void forget(Object list){notes(list,KEPT).clear()}
    /** Outlines kept for this recording without being chosen in the plugin on it: outline ID -> where it came from (origin). */
    static Map<String,String> carried(Object list,Roi[] outlines,String recording) {
        Map made=notes(list,MADE),kept=notes(list,KEPT);Map out=[:]
        outlines.each{Roi outline->String id=CalciumMeasurement.idFor(outline);if(same(kept[id],recording)&&!same(made[id],recording))out[id]=origin(made[id],recording)}
        out
    }
    /** Outlines added on another recording, or cleared at this window's start (clearedHere, by ID), that were not kept for this
     *  one. An outline cleared at another window's start, whose recording was not known, is not among them. */
    static List<Roi> foreign(Object list,Roi[] outlines,String recording,Set<String> clearedHere) {
        Map made=notes(list,MADE),kept=notes(list,KEPT)
        outlines.findAll{Roi outline->String id=CalciumMeasurement.idFor(outline),on=made[id];if(on==UNKNOWN&&!clearedHere.contains(id))on=null
            on!=null&&!same(on,recording)&&!same(kept[id],recording)&&origin(on,recording)!=HERE} as List<Roi>
    }
    /** How many of these outlines the plugin added on this recording. */
    static int madeHere(Object list,Roi[] outlines,String recording){Map made=notes(list,MADE);outlines.count{Roi outline->same(made[CalciumMeasurement.idFor(outline)],recording)} as int}
    /** The recordings a set of outlines came from: their names, with '' for outlines from before. Reads the notes only. */
    static Set<String> origins(Object list,Collection<Roi> outlines,String recording){Map made=notes(list,MADE);outlines.collect{Roi outline->origin(made[idOf(outline)],recording)} as LinkedHashSet}
    static String count(int n){n+(n==1?' outline':' outlines')}
    /** "The ROI Manager holds 115 outlines from another recording: P3 Run1." */
    static String held(Set<String> from,int n) {
        from==([''] as Set)?'The ROI Manager already holds '+count(n)+', and the plugin cannot tell which recording '+(n==1?'it is':'they are')+' from.':
            from.size()==1?'The ROI Manager holds '+count(n)+' from another recording: '+from.first()+'.':from.contains('')?'The ROI Manager holds '+count(n)+' '+CalciumOutlines.from(from)+'.':'The ROI Manager holds '+count(n)+' from other recordings.'
    }
    /** The question about outlines that are not this recording's: asked of them; of the others, own were added by the plugin on
     *  this recording, and all of them stay. */
    static String question(Set<String> from,int asked,int own,int others) {
        String rest=others<=0?'':others==1?(own==1?' The other outline was made on this recording and stays.':' The other outline stays.'):' The other '+count(others)+(own==others?' were made on this recording and stay.':' stay.')
        held(from,asked)+rest+' '+KEEP_ONLY_IF
    }
    /** Where outlines came from, for the record: "from P3 Run1", "from an unknown recording", "from other recordings", and with
     *  outlines from before among them "from P3 Run1 or from before" or "from other recordings or from before". */
    static String from(Collection<String> names) {
        Set distinct=names as LinkedHashSet
        if(distinct==([''] as Set))return 'from an unknown recording'
        if(distinct.size()==1)return 'from '+distinct.first()
        Set named=distinct-['']
        (named.size()==1?'from '+named.first():'from other recordings')+(distinct.contains('')?' or from before':'')
    }
    /** How one carried-over cell is named in the registry and in the list of how cells were chosen. */
    static String carriedText(String origin){origin==HERE?'Drawn on this recording before the plugin started':'Carried over '+from([origin?:''])}
    /** One question with buttons; the position of the button pressed, or -1 when the dialog was closed. A test puts its
     *  own answer in place of this method, so no dialog opens. */
    static int ask(java.awt.Component parent,String title,String text,List<String> options) {
        if(java.awt.GraphicsEnvironment.isHeadless())return -1
        Object[] buttons=options as Object[]
        javax.swing.JOptionPane.showOptionDialog(parent,'<html><body style="width:340px">'+CalciumInstaller.escape(text)+'</body></html>',title,
            javax.swing.JOptionPane.DEFAULT_OPTION,javax.swing.JOptionPane.QUESTION_MESSAGE,null,buttons,buttons[0])
    }
}

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

// One plugin window. Every start of the plugin runs in a class loader of its own, so the open window is found by the note it
// keeps on itself (CalciumOutlines.WINDOW, set where the window is made), not by anything static. Starting the plugin while
// that window is open makes no second one: the open window comes forward and offers what its Change recording button offers.
def openAlready = java.awt.Window.getWindows().find { it instanceof javax.swing.JFrame && it.isDisplayable() && ((javax.swing.JFrame) it).getRootPane().getClientProperty(CalciumOutlines.WINDOW) instanceof Map }
if (openAlready != null) {
    Map open = (Map) ((javax.swing.JFrame) openAlready).getRootPane().getClientProperty(CalciumOutlines.WINDOW)
    binding.setVariable("alreadyOpen", "" + open.recording)
    try { ((Closure) open.again).call() } catch (Throwable problem) { IJ.log(CalciumRelease.NAME + ": " + problem) }
    return
}
// One start at a time. Until it ends, a start notes that it is under way: a second start from the menu while the first makes
// its pictures or waits at one of its questions, before its window exists, opens nothing and says why.
if (!CalciumOutlines.starting()) {
    binding.setVariable("alreadyStarting", true)
    if (!testGui) IJ.showMessage(CalciumRelease.NAME, CalciumOutlines.STARTING_TEXT)
    return
}
try {

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

// What the ROI Manager holds now was put there before this recording. (A test gives the plugin a list of its own.) The ROI
// Manager is not opened here when there is none: it opens with the pictures, as it always did.
RoiManager rm = testGui ? new RoiManager(false) : RoiManager.getInstance()
Roi[] heldBefore = rm == null ? new Roi[0] : rm.getRoisAsArray()

ImagePlus imp = null
if (picked == fromDisk) {
    String recordingPath = start.file      // chosen in the Step 1 dialog
    if (recordingPath == null) return
    IJ.run("Bio-Formats Importer", "open=[" + recordingPath + "] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT")
    imp = WindowManager.getCurrentImage()
} else imp = movies[labels.indexOf(picked)]
if (imp == null || imp.getStackSize() < 2) { IJ.error(CalciumRelease.NAME, "No movie was opened."); return }

// One recording, its own outlines. Each outline the ROI Manager already holds is this recording's own (the plugin added it on
// this recording in an earlier window), carried over (kept for it at Change recording, or by an earlier answer), drawn on this
// recording's own window before the plugin started, or none of these: then the plugin asks once, before any cell is shown or
// counted. Nothing is changed until the start goes on, so Cancel here or at the frames leaves the outlines as they were.
// Outlines that are kept are never counted as chosen blind.
List<String> outlinesSaid = []           // what was done about those outlines, for the record of steps
Closure settleOutlines = { -> }          // what the answer does, once the start goes on
Set<String> clearedAtStart = [] as Set   // the outlines cleared, by ID: removed in the window, where Undo brings them back
String clearedFrom = ""
if (rm != null && heldBefore.length == 0) CalciumOutlines.forget(rm)
else if (heldBefore.length > 0) {
    String startPrint = ""; try { startPrint = CalciumSteps.printOf(imp) } catch (Throwable ignored) { }
    String startMark = CalciumOutlines.mark(imp, startPrint)
    Map sorted = CalciumOutlines.sort(rm, heldBefore, startMark, imp.getID())
    List<Roi> asked = sorted.ask as List<Roi>
    int answer = 1
    if (!asked.isEmpty()) {
        int others = (sorted.kinds as List).count { it != null } - asked.size()
        answer = CalciumOutlines.ask(null, "Outlines in the ROI Manager", CalciumOutlines.question(CalciumOutlines.origins(rm, asked, startMark), asked.size(), (sorted.own as List).size(), others), CalciumOutlines.choices(asked.size()))
        if (answer != 0 && answer != 1) { IJ.showStatus(CalciumRelease.NAME + ": cancelled"); return }
    }
    Map madeOn = CalciumOutlines.notes(rm, CalciumOutlines.MADE)
    List<String> keptFrom = ((sorted.carried as List<Roi>) + (answer == 1 ? asked : [])).collect { Roi outline -> CalciumOutlines.origin(madeOn[CalciumOutlines.idOf(outline)], startMark) }
    List<String> away = keptFrom.findAll { it != CalciumOutlines.HERE }
    int here = (sorted.here as List).size() + (keptFrom.count { it == CalciumOutlines.HERE } as int)
    if (answer == 0) {
        clearedFrom = CalciumOutlines.from(CalciumOutlines.origins(rm, asked, startMark))
        outlinesSaid << "Outlines cleared: " + CalciumOutlines.count(asked.size()) + " " + clearedFrom + " " + (asked.size() == 1 ? "was" : "were") + " in the ROI Manager when this recording was chosen and " + (asked.size() == 1 ? "was" : "were") + " cleared before any cell was chosen."
    }
    if (!away.isEmpty()) outlinesSaid << "Outlines carried over: " + CalciumOutlines.count(away.size()) + " " + CalciumOutlines.from(away) + " " + (away.size() == 1 ? "was" : "were") + " in the ROI Manager when this recording was chosen and " + (away.size() == 1 ? "was" : "were") + " kept for it."
    if (here > 0) outlinesSaid << "Outlines drawn before the plugin started: " + CalciumOutlines.count(here) + " drawn on this recording " + (here == 1 ? "was" : "were") + " in the ROI Manager when the plugin started and " + (here == 1 ? "was" : "were") + " kept."
    settleOutlines = { -> clearedAtStart.addAll(CalciumOutlines.settle(rm, heldBefore, sorted, startMark, answer == 0)) }
}

IJ.showStatus("Making the cells visible...")
if (imp.getBitDepth() == 32) IJ.log(CalciumRelease.NAME + ": this movie is 32-bit, so it looks like a processed copy. Background is meant to be taken on the original movie.")
def steps = new CalciumSteps()
steps.prepare(imp, look || blindStart)   // blind: the response is looked for only to offer when the stimulus came; it is never shown while choosing
outlinesSaid.each { steps.say(it) }
if (blindStart) {
    int[] frames = CalciumStartDialog.askBlindFrames(steps)
    if (frames == null) { IJ.showStatus(CalciumRelease.NAME + ": cancelled"); return }
    steps.startBlind(frames[0], frames[1])
}
settleOutlines()   // the start goes on: the answer about the outlines takes effect
steps.colourPic.show()
IJ.run(steps.colourPic, "Set... ", "zoom=75")
if (!testGui) rm = RoiManager.getRoiManager()
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
// The recording's own window can be closed in Fiji while this window stays open: its frames are gone then, and nothing can be
// measured, pictured or traced until the same recording is open again. The card says so, and so does every action that needs
// the frames, in the place where it reports.
String ownMark = CalciumOutlines.mark(imp, steps.recordingPrint)   // names this recording in the notes kept with the ROI Manager
String recordingFile = CalciumOutlines.fileOf(imp)                 // '' when Fiji does not know the file
String reopenProblem = ""        // what happened when Open the recording again was pressed and the plugin could not go on
boolean closedShown = false      // what recordingReady (below) last found: the card says closed only then
String whyShown = ""             // and what the card says after it
Closure whenReopened = null      // set below: Step 4 goes on with the reopened recording too
def describeRecording = { ->
    recordingName.setWhole(imp.getTitle())
    // a recording that did not come from a file cannot be opened again, so the card does not say to open it again
    String closedHead = recordingFile == "" ? "The recording is closed in Fiji." : CalciumSteps.CLOSED
    String closedLine = closedShown ? "<font color='" + CalciumTheme.hex(CalciumTheme.BAD) + "'><b>" + closedHead + "</b>" + (whyShown ? " " + CalciumInstaller.escape(whyShown) : "") + "</font><br>" : ""
    recordingNote.setText("<html>" + closedLine + "<font color='" + MUTED + "'>" + steps.n + " frames · " + (steps.blind ? "cells chosen blind on frames " + steps.blindFirst + "–" + steps.blindLast : steps.found ? "response near frame " + steps.onset + " (magenta)" : (steps.looked ? "no response found (no magenta)" : "response not looked for")) + "</font>" + (movedNote() ? "<br><font color='#9A5B00'>" + movedNote() + "</font>" : "") + "</html>")
    recordingNote.setToolTipText(imp.getTitle() + ". " + (steps.blind ? "Cells are chosen on frames " + steps.blindFirst + "–" + steps.blindLast + ", before the response" + (steps.found ? "; the cells start to rise near frame " + steps.onset : "") + ". The response is not shown while choosing." : steps.found ? "Brightness rises near frame " + steps.onset + "; the response can be shown in magenta." : (steps.looked ? "No response was found, so there is no magenta. Every cell is still shown in green." : "The response has not been looked for. Choose Response in Step 2 to look for it.")))
}
describeRecording()
def changeBtn = btn("‹ Change recording")
changeBtn.setToolTipText("Back to the first screen, to choose another recording or results folder")
def reopenBtn = btn("Open the recording again"); reopenBtn.setVisible(false)   // shown only while the recording is closed and its file is there
reopenBtn.setToolTipText("Reads the recording's file again. The plugin goes on with the same cells when it is the same recording.")
s1.body.add(recordingName); s1.body.add(recordingNote); s1.body.add(row(changeBtn, reopenBtn))
def canReopen = { -> recordingFile != "" && new File(recordingFile).isFile() }
// Every action that reads the recording's frames asks here first. False while the recording is closed in Fiji. When the same
// recording has been opened again (same size, frame count, time base and pixels), the plugin goes on with it.
def recordingReady = { ->
    boolean ready = !CalciumSteps.closed(steps.movie)
    if (!ready) try { ready = steps.recordingReady() } catch (Throwable ignored) { ready = false }
    if (ready && !steps.movie.is(imp)) { imp = steps.movie; reopenProblem = ""; if (whenReopened != null) whenReopened.call(imp) }
    String why = ready ? "" : recordingFile == "" ? "It did not come from a file, so it cannot be opened again." : canReopen() ? reopenProblem : "The file is no longer there: " + recordingFile
    if (closedShown == ready || why != whyShown) {
        closedShown = !ready; whyShown = why
        describeRecording(); reopenBtn.setVisible(!ready && canReopen()); frame.refreshContents()
    }
    ready
}
// The same question for work running on its own thread: changes nothing on the window.
def recordingClosed = { -> CalciumSteps.live(steps.movie, steps.recordingPrint) == null }

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
whenReopened = { again -> if (analysisMeasurement != null) analysisMeasurement.movie = again }
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
// Outlines kept for this recording were not chosen in the plugin on it: outline ID -> where they came from (the recording they
// were made on, '' when not known, or CalciumOutlines.HERE when drawn on this recording before the plugin started).
def carriedIn = { Roi[] all -> CalciumOutlines.carried(rm, all, ownMark) }
def basisFor = { Roi[] all ->
    Map carried = carriedIn(all)
    List away = carried.values().findAll { it != CalciumOutlines.HERE }; int here = carried.size() - away.size()
    String kept = (away.isEmpty() ? '' : '; ' + away.size() + ' of the ' + all.length + ' outlines carried over ' + CalciumOutlines.from(away)) + (here == 0 ? '' : '; ' + here + ' of the ' + all.length + ' outlines drawn before the plugin started')
    if (!steps.blind) {
        // cells chosen blind first and others added afterwards: say so, with the counts (the registry has it cell by cell)
        int blindly = steps.blindLast > 0 && blindCellKeys != null ? (all.count { blindCellKeys.contains(cellKey(it)) && !carried.containsKey(CalciumMeasurement.idFor(it)) } as int) : 0
        return (blindly == 0 ? selectionBasis : 'Blind on frames ' + steps.blindFirst + '-' + steps.blindLast + ', then ' + (selectionBasis.startsWith('Response') ? 'response-assisted' : 'all frames') + ': ' + blindly + ' chosen blind, ' + (all.length - blindly - carried.size()) + ' added afterwards') + kept
    }
    ((!responseSeen || measuredKey(all) == blindKey) ? steps.blindBasis() :
        'Picture of frames ' + steps.blindFirst + '-' + steps.blindLast + ' (before the response), not blind: the response had been shown before these outlines were final') + kept
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
Set<String> seenIds = new HashSet<String>()   // the outline IDs of the last look: an outline not among them is new in the list
def cellRois = { ->
    Roi[] accepted=rm.getRoisAsArray()
    Roi[] unique=accepted.collect{it} as Roi[]
    CalciumMeasurement.ensureIds(unique)
    for(int i=0;i<accepted.length;i++)if(!accepted[i].is(unique[i]))rm.setRoi(unique[i],i)
    CalciumOutlines.claim(rm,unique,ownMark,steps.colourPic.getID(),seenIds)
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
        recordingReady()   // the Recording card says when the recording has been closed in Fiji, and when it is back
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
    if (!recordingReady()) { showView(currentView(), CalciumSteps.CLOSED); return }
    int[] range = CalciumStartDialog.askBlindFrames(steps); if (range == null) return
    busyView(true); selectionNote.setText("<html>Making the picture of frames " + range[0] + "–" + range[1] + "…</html>")
    Thread.start('Calcium pre-stimulus picture') {
        float[] picture = null; String problem = null
        try { picture = steps.framesPicture(range[0], range[1]) } catch (Throwable t) { problem = t.message ?: t.toString() }
        if (CalciumSteps.closed(steps.movie)) picture = null       // closed in Fiji while the picture was made
        javax.swing.SwingUtilities.invokeLater {
            busyView(false)
            if (picture == null) { showView(currentView(), CalciumSteps.closed(steps.movie) ? CalciumSteps.CLOSED : 'Could not make the picture of the frames before the response: ' + problem); return }
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
    if (!steps.looked && !recordingReady()) { showView(currentView(), CalciumSteps.CLOSED); return }   // looking for the response reads the frames
    boolean explained = respondersExplained; respondersExplained = false   // the question before measuring has already said it
    if (steps.blind && !testGui && !explained && javax.swing.JOptionPane.showConfirmDialog(frame, "<html><body style='width:320px'>The cells you have now stay recorded as <b>chosen blind</b>, as long as you do not change them.<br><br>The response is shown next and the plugin looks for the cells that only show when they respond. Cells you add from now on are recorded as <b>added after the response was shown</b>. Step 4 gives % responding for both groups.</body></html>", "Add cells that respond", javax.swing.JOptionPane.OK_CANCEL_OPTION) != javax.swing.JOptionPane.OK_OPTION) return
    boolean first = !steps.looked
    busyView(true); if (first) selectionNote.setText("<html>Looking for the response in all " + steps.n + " frames…</html>")
    Thread.start('Calcium response picture') {
        String problem = null
        try { if (!steps.looked) steps.findResponse() } catch (Throwable t) { problem = t.message ?: t.toString() }
        if (first && CalciumSteps.closed(steps.movie)) { steps.looked = false; steps.found = false; problem = CalciumSteps.CLOSED }   // closed in Fiji meanwhile: looked for again later
        javax.swing.SwingUtilities.invokeLater {
            busyView(false); describeRecording()
            String previous = currentView()
            if (problem != null) { showView(previous, problem == CalciumSteps.CLOSED ? problem : 'Could not look for the response: ' + problem); return }
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
def addCell = { Roi r0 ->
    def r = (Roi) r0.clone(); r.setPosition(0); r.setStrokeColor(null)
    // An outline identical to one the list already holds gets an ID of its own now, as it would at the next look, so that the
    // other one keeps its note of where it came from.
    String id = CalciumOutlines.idOf(r)
    if (rm.getRoisAsArray().any { it != null && CalciumOutlines.idOf(it) == id }) { String name = r.getName(); r = new ij.io.RoiDecoder(RoiEncoder.saveAsByteArray(r), 'cell.roi').getRoi(); r.setName(name); CalciumMeasurement.assignId(r, true) }   // the copy keeps the name, so the ROI Manager labels it as before
    int before = rm.getCount(); rm.addRoi(r)
    if (rm.getCount() > before) CalciumOutlines.made(rm, rm.getRoi(rm.getCount() - 1), ownMark)   // (the ROI Manager turns down an outline equal to its last one)
}
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
    else {
        if (!recordingReady()) trace.show(null, CalciumSteps.CLOSED, 0, 0)   // a trace is read from the recording's frames
        else { reveal(); trace.show(steps.traceOf(target, cells), cap, steps.found ? steps.p0 : 0, steps.found ? steps.p1 : 0) }
        if (!trace.isVisible()) { trace.setVisible(true); frame.refreshContents() }
    }
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
    if (!recordingReady()) { info.setText("<html>" + CalciumSteps.CLOSED + "</html>"); return }
    // "again" asks for a new Cellpose run only when the outlines on screen came from the same kind of picture as now
    forceFresh = proposals.any { (it.blind == true) == (steps.blind == true) }
    detecting = true
    detectionRunner.start(suggestBtn, detectionProgress,
        { task -> suggest(task) },   // reads only the pictures made at the start, so it needs no frame once it runs
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
    if (choice == 0) { rm.reset(); CalciumOutlines.forget(rm) }
    else CalciumOutlines.keepForNext(rm, cellRois(), ownMark)   // the next recording takes them as carried over, and does not ask again
    timer.stop(); detectionProgress.close(); frame.dispose()
    steps.colourPic.changes = false; steps.colourPic.close()
    if (!testGui) relaunch()
}
// Change recording asks about the outlines; lead says first that the plugin was started again from the menu, and then it asks
// also when the ROI Manager is empty, so that the open window is never closed without a question.
def askChange = { String lead ->
    int count = rm.getCount()
    if (count == 0) {
        if (!lead) { changeRecording(1); return }
        if (CalciumOutlines.ask(frame, 'Change recording', lead + 'Go back to choose another recording?', ['Go back', 'Cancel']) == 0) changeRecording(1)
        return
    }
    String unmeasured = analysisMeasurement == null ? (count == 1 ? ' It has not been measured yet.' : ' They have not been measured yet.') : ''
    int answer = CalciumOutlines.ask(frame, 'Change recording', lead + 'The ROI Manager holds ' + CalciumOutlines.count(count) + '.' + unmeasured +
        (count == 1 ? ' Keep it' : ' Keep them') + ' only if the next recording shows the same cells in the same place. On another recording, kept outlines are recorded as carried over, never as chosen blind.',
        ['Clear outlines and go back', 'Keep outlines and go back', 'Cancel'])
    if (answer == 0 || answer == 1) changeRecording(answer)
}
changeBtn.addActionListener({ e -> askChange('') } as java.awt.event.ActionListener)
// The note by which a later start of the plugin finds this window (see "One plugin window" at the top of this part).
frame.getRootPane().putClientProperty(CalciumOutlines.WINDOW, [recording: imp.getTitle(), again: { ->
    Runnable offer = { if (frame.getState() == java.awt.Frame.ICONIFIED) frame.setState(java.awt.Frame.NORMAL); frame.toFront(); askChange('The plugin is already open for ' + imp.getTitle() + '. ') } as Runnable
    if (javax.swing.SwingUtilities.isEventDispatchThread()) offer.run() else javax.swing.SwingUtilities.invokeAndWait(offer)
}])
// Open the recording again: the file is read the way the plugin opens a file at the start. Whatever opens is taken only when it
// is the same recording (recordingReady); the cells, the pictures and everything chosen so far stay as they are.
Closure openAgain = { String path -> IJ.run("Bio-Formats Importer", "open=[" + path + "] color_mode=Grayscale view=Hyperstack stack_order=XYCZT") }
reopenBtn.addActionListener({ e ->
    reopenBtn.setEnabled(false)
    Set<Integer> before = (WindowManager.getIDList() ?: new int[0]) as Set
    Thread.start('Calcium open the recording again') {
        if (canReopen()) try { openAgain.call(recordingFile) } catch (Throwable failure) { IJ.log(CalciumRelease.NAME + ": the recording could not be opened again" + (failure.message ? ": " + failure.message : ".")) }
        boolean opened = ((WindowManager.getIDList() ?: new int[0]) as List).any { !before.contains(it) }
        javax.swing.SwingUtilities.invokeLater {
            reopenBtn.setEnabled(true)
            // the Recording card says what happened: the file is gone, it did not open, or it opened as another recording
            if (!recordingReady() && canReopen()) { reopenProblem = opened ? "The file did not open as the same recording." : "The file could not be opened."; recordingReady() }
            refresh(true)
        }
    }
} as java.awt.event.ActionListener)
openBtn.addActionListener({ e -> try { java.awt.Desktop.getDesktop().open(outDir) } catch (Throwable t) { IJ.log("Results folder: " + outDir.getPath()) } } as java.awt.event.ActionListener)
frame.addWindowListener(new java.awt.event.WindowAdapter() { void windowClosing(java.awt.event.WindowEvent e) { timer.stop(); detectionProgress.close() } })
// Sensitivity check after blind selection: response-assisted detection once, compared with the blind suggestions. Returns the
// comparison and its report, or a message saying why it did not run. Its outlines are never added to the cells or measured.
def runSensitivity = { Roi[] cells, String prefix, Closure tell ->
    if (recordingClosed()) return [text: 'Sensitivity check not run: ' + CalciumSteps.CLOSED.uncapitalize()]
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
    // Nothing is measured, saved or changed when the recording is closed, or when the ROI Manager holds outlines the plugin
    // added on another recording that were not kept for this one.
    if (!recordingReady()) { result.setText("<html>" + CalciumSteps.CLOSED + "</html>"); return }
    def sp = splitCells()
    List<Roi> foreign = CalciumOutlines.foreign(rm, sp.all as Roi[], ownMark, clearedAtStart)
    if (!foreign.isEmpty()) {
        Set<String> origins = CalciumOutlines.origins(rm, foreign, ownMark)
        int others = (sp.all as Roi[]).length - foreign.size()
        boolean one = foreign.size() == 1
        String pressed = go.getText()   // Measure and save, or Measure again once measured
        int answer = CalciumOutlines.ask(frame, pressed, CalciumOutlines.question(origins, foreign.size(), CalciumOutlines.madeHere(rm, sp.all as Roi[], ownMark), others), CalciumOutlines.choices(foreign.size()))
        if (answer == 1) {
            CalciumOutlines.keep(rm, foreign, ownMark)
            steps.say("Outlines carried over: " + CalciumOutlines.count(foreign.size()) + " " + CalciumOutlines.from(origins) + " " + (one ? "was" : "were") + " in the ROI Manager when the cells were measured and " + (one ? "was" : "were") + " kept for this recording.")
        } else {
            if (answer == 0) {
                Set<String> gone = foreign.collect { CalciumMeasurement.idFor(it) } as Set
                Roi[] now = rm.getRoisAsArray()
                cellEdits.remove((0..<now.length).findAll { gone.contains(CalciumMeasurement.idFor(now[it])) })
                refresh(true)
                result.setText("<html>" + CalciumOutlines.count(foreign.size()) + " " + CalciumInstaller.escape(CalciumOutlines.from(origins)) + " " + (one ? "was" : "were") + " removed. Nothing was measured. Check the cells, then press " + CalciumInstaller.escape(pressed) + ".</html>")
            } else result.setText("<html>Not measured. " + CalciumInstaller.escape(CalciumOutlines.held(origins, foreign.size())) + (one ? " Clear it or keep it" : " Clear them or keep them") + " for this recording, then measure.</html>")
            return
        }
    }
    Roi[] cells = (Roi[]) sp.cells; Roi[] edgeCells = (Roi[]) sp.edge
    final Map<String, String> carried = carriedIn(sp.all as Roi[])   // outline ID -> where it came from
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
            }
            result.setText("<html>Measuring " + cells.length + " cells on " + steps.n + " frames...</html>")
            steps.log.removeAll{it.startsWith('Cell selection basis:') || it.startsWith('Sensitivity check') || it.startsWith('Chosen blind or later:') || it.startsWith('Carried over:') || it.startsWith('Drawn before the plugin started:')};steps.say('Cell selection basis: '+measuredSelectionBasis+'. Accepted ROI registry includes premeasurement exclusions.')
            // A cell whose outline was kept from another recording was not chosen on this one: it is in neither group below (no
            // value in chosenBlind) and is named as carried over in the list, the registry and Step 4.
            Map<String, String> measuredCarried = carried.subMap(names.findAll { carried.containsKey(it) })
            List<String> carriedAway = measuredCarried.values().findAll { it != CalciumOutlines.HERE } as List<String>   // from another recording
            int drawnHere = measuredCarried.size() - carriedAway.size()                                            // drawn on this one before the plugin started
            Map chosenBlind = null
            List chosenRows = null   // written once the frames have been read, so that a recording closed meanwhile leaves no file changed
            if (steps.blindLast > 0 && (steps.blind || blindCellKeys != null)) {
                chosenBlind = [:]
                cells.eachWithIndex { Roi c, int i -> chosenBlind[names[i]] = measuredCarried.containsKey(names[i]) ? null : (blindCellKeys != null ? blindCellKeys.contains(cellKey(c)) : !notBlindKeys.contains(cellKey(c))) }
                int nBlind = chosenBlind.values().count { it } as int
                chosenRows = [['cell', 'chosen']] + names.collect { [it, measuredCarried.containsKey(it) ? CalciumOutlines.carriedText(measuredCarried[it]).uncapitalize() : chosenBlind[it] ? 'blind (frames ' + steps.blindFirst + '-' + steps.blindLast + ')' : 'after the response was shown'] }
                steps.say('Chosen blind or later: ' + nBlind + ' of the ' + cells.length + ' cells were chosen blind on frames ' + steps.blindFirst + '-' + steps.blindLast + ' and are unchanged since; ' + (cells.length - nBlind - measuredCarried.size()) + ' were added or reshaped after the response was shown' + (carriedAway.isEmpty() ? '' : '; ' + carriedAway.size() + (carriedAway.size() == 1 ? ' was' : ' were') + ' carried over and not chosen on this recording') + (drawnHere == 0 ? '' : '; ' + drawnHere + (drawnHere == 1 ? ' was' : ' were') + ' drawn before the plugin started') + '. Step 4 gives % responding for each group; only the blind group is free of selection on the response. List: ' + base + '_cells_chosen_blind_or_after_response.csv')
            }
            if (carriedAway.size() == 1) steps.say('Carried over: the outline of 1 of the ' + cells.length + ' cells was carried over ' + CalciumOutlines.from(carriedAway) + '. It was not chosen on this recording and is never counted as chosen blind. Step 4 gives % responding for it separately.')
            else if (!carriedAway.isEmpty()) steps.say('Carried over: the outlines of ' + carriedAway.size() + ' of the ' + cells.length + ' cells were carried over ' + CalciumOutlines.from(carriedAway) + '. They were not chosen on this recording and are never counted as chosen blind. Step 4 gives % responding for them separately.')
            if (drawnHere == 1) steps.say('Drawn before the plugin started: the outline of 1 of the ' + cells.length + ' cells was drawn on this recording before the plugin started. It was not chosen in the plugin and is never counted as chosen blind. Step 4 counts it under Carried over.')
            else if (drawnHere > 0) steps.say('Drawn before the plugin started: the outlines of ' + drawnHere + ' of the ' + cells.length + ' cells were drawn on this recording before the plugin started. They were not chosen in the plugin and are never counted as chosen blind. Step 4 counts them under Carried over.')
            def summary = steps.measureAndSave(cells, names, outDir.getPath() + File.separator + base + "_")
            if (chosenRows != null) CalciumCsv.write(new File(outDir, base + '_cells_chosen_blind_or_after_response.csv'), chosenRows)
            if (corrected != null && saveMovie.isSelected()) IJ.saveAsTiff(corrected, outDir.getPath() + File.separator + base + "_background_corrected_movie.tif")
            if (corrected != null && !testGui) { corrected.show(); corrected.setSlice(Math.min(steps.n, steps.found ? steps.p0 + 5 : 1)); corrected.setDisplayRange(0, CalciumSteps.SHOW_MAX); corrected.updateAndDraw() }
            IJ.log(summary); IJ.log("Saved in: " + outDir.getPath()); IJ.log("Step 4 uses full-precision unlifted fluorescence: " + base + "_gray_values_background_corrected_native.csv")
            File registryFile=new File(outDir,base+'_roi_registry.csv')
            Map selectionByCell=chosenBlind==null?null:chosenBlind.collectEntries{id,blindly->[id,blindly?steps.blindBasis():'Added after the response or all frames were shown (blind frames '+steps.blindFirst+'-'+steps.blindLast+')']}
            if(!carried.isEmpty())selectionByCell=(selectionByCell?:[:])+carried.collectEntries{id,from->[id,CalciumOutlines.carriedText(from)]}
            CalciumCsv.write(registryFile,CalciumMeasurement.registry(sp.all as Roi[],cells,steps.lastBackgroundPixels,steps,measuredSelectionBasis,selectionByCell))
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
                acceptedRois:acceptedRoisFile,acceptedRoisHash:CalciumInstaller.digest(acceptedRoisFile.bytes),recordingPrint:steps.recordingPrint]
            measured.inputHash=CalciumInstaller.digest(measured.csv.bytes)
            if(!measuredCarried.isEmpty())measured.carriedOver=measuredCarried
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
        } catch (Throwable t) {
            boolean closed = recordingClosed()   // the recording was closed in Fiji while it was measured: said in plain words
            IJ.log(CalciumRelease.NAME + ": " + (closed ? CalciumSteps.CLOSED : t))
            result.setText("<html>" + (closed ? CalciumSteps.CLOSED : "Could not measure: " + t) + "</html>"); go.setText("Measure and save"); go.setEnabled(true)
        }
    }
} as java.awt.event.ActionListener)
// Clear them at the start: the outlines go now, as one change to the cells, so that Undo brings them back.
if (!clearedAtStart.isEmpty()) {
    Roi[] held = rm.getRoisAsArray()
    int gone = cellEdits.remove((0..<held.length).findAll { clearedAtStart.contains(CalciumOutlines.idOf(held[it])) })
    if (gone > 0) updateEditHint("Cleared " + CalciumOutlines.count(gone) + " " + CalciumInstaller.escape(clearedFrom) + " when this recording was chosen. <b>Undo</b> brings " + (gone == 1 ? "it" : "them") + " back.")
}
if (testGui) { binding.setVariable("hooks", [nextLine: nextLine, addRespondersBtn: addRespondersBtn, addMissedBtn: addMissedBtn, continueBtn: continueBtn, acceptAllBtn: acceptAllBtn, autoSuggestAsked: { -> autoSuggestAsked }, legendClicked: { String k, boolean o -> legendClicked.call(k, o) }, hiddenKinds: hiddenKinds, cellRois: cellRois, sensBox: sensBox, basisFor: basisFor, responseSeen: { -> responseSeen }, reveal: { -> reveal() }, runSensitivity: runSensitivity, trace: trace, respondsLegend: respondsLegend, viewAllBtn: viewAllBtn, viewBaselineBtn: viewBaselineBtn, viewResponseBtn: viewResponseBtn, selectionNote: selectionNote, selectionBasis: { -> selectionBasis }, changeBtn: changeBtn, changeRecording: changeRecording, continueBtn: continueBtn, backToCellsBtn: backToCellsBtn, backToMeasureBtn: backToMeasureBtn, cards: [s1, s2, s3, s4], recordingNote: recordingNote, detectionProgress: detectionProgress, suggestBtn: suggestBtn, edgeBox: edgeBox, addBtn: addBtn, setForceFresh: { boolean v -> forceFresh = v }, setupCellpose: setupCellpose, findCellposePython: findCellposePython, setupBtn: setupBtn, setupRow: setupRow, cpLabel: cpLabel, showCellposeState: showCellposeState, setupText: SETUP_TEXT, saveMovie: saveMovie, makeMovie: makeMovie, rm: rm, refresh: refresh, go: go, status: status, steps: steps, outDir: outDir, viewSwitch: viewSwitch, hiddenLine: hiddenLine, showSwitches: legendEntries, moreBtn: moreBtn, recordingName: recordingName, frame: frame, suggest: suggest, acceptAt: acceptAt, acceptAll: acceptAll, hoverAt: hoverAt, proposals: proposals, result: result, info: info, removeSelected: removeSelected, undoLast: undoLast, root: root, ownMark: ownMark, reopenBtn: reopenBtn, recordingReady: recordingReady, recordingClosed: recordingClosed, carriedIn: carriedIn, setOpenAgain: { Closure c -> openAgain = c }, analysisMeasurement: { -> analysisMeasurement }, editHint: editHint, undoBtn: undoBtn, askChange: askChange, cellEdits: cellEdits]); refresh(true); return }
frame.setVisible(true)
timer.start()
refresh(true)
} finally { CalciumOutlines.ended() }   // the start is over: a later start may begin

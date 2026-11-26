// Ca²⁺ Imaging (Wang Lab)                                                                  version 2026.11.24.1
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
    /** Whether a recording's unit is micrometres: one test for the scale bar and for the widths of the rings around each cell. */
    static boolean micrometres(String unit) { return unit != null && (unit.startsWith("micro") || unit == "um" || unit == "µm") }

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
        say(CalciumRelease.NAME + ", version 2026.11.24.1")
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
    int[] pixelsOf(Roi roi) { return pixelsOf(roi, w, h) }
    /** The same for a picture of the given size: the one pixel mask of an outline, used wherever an outline is measured. */
    static int[] pixelsOf(Roi roi, int w, int h) {
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
        // Around each cell (CalciumRingAndCore): measured in the same pass and saved beside the tables. It is a record of its
        // own: a failure here leaves the measurement as it is and is said in the record.
        CalciumRingAndCore rings = null; Throwable ringProblem = null
        String ringFile = null   // the file being written when the measurement around the cells stopped
        try {
            ij.measure.Calibration scale = movie.getCalibration()
            rings = new CalciumRingAndCore(cellPx, w, h, scale.pixelWidth, scale.pixelHeight, scale.getUnit())
            ringFile = CalciumCallChecks.VALUES
            rings.open(new File(prefix + CalciumCallChecks.VALUES), names)
        } catch (Throwable t) { ringProblem = t; if (rings != null) rings.abandon(); rings = null }
        try {
            for (int f = 0; f < n; f++) {
                ImageProcessor ip = st.getProcessor(f + 1)
                for (int c = 0; c < nc; c++) {
                    double s = 0; for (int i : cellPx[c]) s += ip.getf(i); raw[f][c] = cellPx[c].length > 0 ? s / cellPx[c].length : Double.NaN
                    s = 0; for (int i : bg[c]) s += ip.getf(i); back[f][c] = bg[c].length > 0 ? s / bg[c].length : Double.NaN
                }
                if (rings != null) try { rings.frame(f + 1, ip) } catch (Throwable t) { ringProblem = t; rings.abandon(); rings = null }
            }
        } catch (Throwable t) { if (rings != null) rings.abandon(); throw t }   // Measure stops as it always did, and leaves no part of the ring table
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
        try {
            new File(prefix + "gray_values_raw.csv").text = a.toString()
            new File(prefix + "background_next_to_each_cell.csv").text = b.toString()
            new File(prefix + "gray_values_background_corrected_for_SimpleCa.csv").text = d.toString()
            // Native analysis preserves the measured doubles. Legacy three-decimal exports remain compatible.
            new File(prefix + "gray_values_background_corrected_native.csv").setText(CalciumMeasurement.fullPrecision(raw, back, names), "UTF-8")
        } catch (Throwable t) { if (rings != null) rings.abandon(); throw t }
        // What was measured around each cell goes beside the tables, with the raw table it belongs to named in it. The field's
        // movement is saved only when Step 1 could measure it (CalciumMotion.measurable).
        String ringRecord
        boolean movement = motionShifts != null && motionShifts.length > 0 && CalciumMotion.measurable(w, h, motionShifts.length)
        try {
            if (rings == null) throw ringProblem != null ? ringProblem : new IOException("not measured")
            rings.close()
            ringFile = CalciumCallChecks.MOVEMENT
            CalciumCallChecks.saveMovement(new File(prefix + CalciumCallChecks.MOVEMENT), movement ? motionShifts : null)
            ringFile = CalciumCallChecks.PIXELS
            CalciumCallChecks.savePixels(prefix, Arrays.asList(names), cells, cellPx, rings)
            ringRecord = CalciumCallChecks.record(rings, new File(prefix).getName(), movement)
        } catch (Throwable t) {
            if (rings != null) rings.abandon()
            CalciumCallChecks.removeBeside(prefix)
            Throwable why = ringProblem != null ? ringProblem : t
            IJ.log(CalciumRelease.NAME + ": around each cell, not saved: " + why)   // the details go to the Log; the record stays plain
            ringRecord = "Around each cell: the rings and the core were not saved" + (why instanceof IOException && ringFile != null ? ": " + new File(prefix).getName() + ringFile + " could not be written" : "") + ". The cell measurements are unchanged."
        }
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
        log.removeAll { String s -> s.startsWith("Around each cell: ") }               // only the latest measurement's line stays in the record
        say(ringRecord)
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
    static final String VERSION = "2026.11.24.1"
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

/** One line of text in a small font that keeps its height whatever it says: cut at its end with an ellipsis when it does not
 *  fit, with the whole line in its tooltip. */
class CalciumLine extends javax.swing.JComponent {
    String text = ""
    CalciumLine(Color colour) { setFont(CalciumTheme.font(java.awt.Font.PLAIN, -1f)); setForeground(colour); setAlignmentX(0f) }
    void setText(String words) { text = words ?: ""; setToolTipText(text ?: null); repaint() }
    java.awt.Dimension getPreferredSize() { def fm = getFontMetrics(getFont()); def i = getInsets(); new java.awt.Dimension(Math.min(CalciumText.NARROW, fm.stringWidth(text)) + i.left + i.right, fm.getHeight() + i.top + i.bottom) }
    java.awt.Dimension getMinimumSize() { new java.awt.Dimension(0, (int) getPreferredSize().height) }
    java.awt.Dimension getMaximumSize() { new java.awt.Dimension(Integer.MAX_VALUE, (int) getPreferredSize().height) }
    protected void paintComponent(java.awt.Graphics g0) {
        def g = (java.awt.Graphics2D) g0.create()
        try { CalciumTheme.smooth(g); g.setFont(getFont()); def fm = g.getFontMetrics(); def i = getInsets(); g.setColor(getForeground()); g.drawString(CalciumTheme.fit(text, fm, getWidth() - i.left - i.right), i.left, i.top + fm.getAscent()) } finally { g.dispose() }
    }
}

/** Empty room as tall as a closure says, worked out each time it is laid out: what the rows above it leave of a height kept
 *  for them, so the room is under them and not between them. */
class CalciumRoom extends javax.swing.JComponent {
    Closure tall
    CalciumRoom(Closure tall) { this.tall = tall; setAlignmentX(0f); setOpaque(false) }
    java.awt.Dimension getPreferredSize() { new java.awt.Dimension(0, Math.max(0, tall.call() as int)) }
    java.awt.Dimension getMinimumSize() { getPreferredSize() }
    java.awt.Dimension getMaximumSize() { new java.awt.Dimension(Integer.MAX_VALUE, (int) getPreferredSize().height) }
}

/** A one-line text box that shows a hint while it is empty, as a page's placeholder does, and takes at most "most" characters
 *  (counted as String.length counts them), as a page's maxlength does. */
class CalciumHintField extends javax.swing.JTextField {
    String hint = ""
    CalciumHintField(String hint, int most, int columns) {
        super(columns); this.hint = hint ?: ""; setFont(CalciumTheme.font())
        ((javax.swing.text.AbstractDocument) getDocument()).setDocumentFilter(new javax.swing.text.DocumentFilter() {
            void insertString(javax.swing.text.DocumentFilter.FilterBypass fb, int offset, String text, javax.swing.text.AttributeSet a) { replace(fb, offset, 0, text, a) }
            void replace(javax.swing.text.DocumentFilter.FilterBypass fb, int offset, int length, String text, javax.swing.text.AttributeSet a) {
                String add = text ?: ""; int room = Math.max(0, most - (fb.getDocument().getLength() - length))
                if (add.length() > room) { int cut = room; if (cut > 0 && Character.isHighSurrogate(add.charAt(cut - 1))) cut--; add = add.substring(0, cut) }
                fb.replace(offset, length, add, a)
            }
        })
    }
    void setHint(String text) { hint = text ?: ""; repaint() }
    protected void paintComponent(java.awt.Graphics g0) {
        super.paintComponent(g0)
        if (!hint || getDocument().getLength() > 0) return
        def g = (java.awt.Graphics2D) g0.create()
        try { CalciumTheme.smooth(g); g.setFont(getFont()); def fm = g.getFontMetrics(); def i = getInsets(); g.setColor(CalciumTheme.FAINT); g.drawString(hint, i.left, i.top + (int) ((getHeight() - i.top - i.bottom - fm.getHeight()) / 2) + fm.getAscent()) } finally { g.dispose() }
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
    /** The outline as the ROI registry identifies it (ROI_sha256): the ROI file of the outline under its ID, on no frame. */
    static String outlineHash(Roi roi,String id){Roi canonical=(Roi)roi.clone();canonical.setName(id);canonical.setPosition(0);CalciumInstaller.digest(RoiEncoder.saveAsByteArray(canonical))}
    static List registry(Roi[] accepted,Roi[] measured,int[] backgroundPixels,CalciumSteps steps,String selection,Map selectionByCell=null) {
        def measuredIds=measured.collect{idFor(it)}
        def rows=[['Cell','ROI_ID','Accepted','Included_in_measurement','Exclusion_reason','ROI_sha256','Area_pixels','Background_pixels','Near_edge','Selection_basis']]
        accepted.each {roi->
            String id=idFor(roi);int index=measuredIds.indexOf(id);boolean included=index>=0
            rows<<[id,id,true,included,included?'':'Edge exclusion policy',outlineHash(roi,id),roi.getContainedPoints().length,
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
    /** Whether shifts() can measure a movie of this size: it gives zeros, not a measurement, when it cannot. */
    static boolean measurable(int width, int height, int frames) {
        int sw = (int) (width / REDUCE), sh = (int) (height / REDUCE)
        return frames >= 3 && sw > 4 * REACH && sh > 4 * REACH
    }
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
    /** How far the field is displaced at every frame, in pixels, counted from the place (rx, ry): the middle value of five
     *  frames around it, so single noisy frames do not count, with a jolt of one or two frames kept. This is the size "summary"
     *  compares with SLIGHT and MOVED. */
    static double[] displacement(double[][] shifts, double rx, double ry) {
        int n = shifts.length
        double[] raw = new double[n], size = new double[n]
        for (int f = 0; f < n; f++) raw[f] = Math.hypot(shifts[f][0] - rx, shifts[f][1] - ry)
        for (int f = 0; f < n; f++) { int a = Math.max(0, f - 2), b = Math.min(n, f + 3); double[] part = new double[b - a]; for (int k = a; k < b; k++) part[k - a] = raw[k]; size[f] = Math.max(median(part), n < 5 ? raw[f] : 0d) }
        // a jolt of one or two frames is kept: the middle value would hide it
        for (int f = 0; f < n; f++) if (raw[f] >= MOVED && size[f] < MOVED) size[f] = raw[f]
        return size
    }
    /** Where the field is over the frames first..last (counted from 1, both in the recording): [x, y, within], the median shift
     *  and how far the field moves inside those frames themselves (the second-largest distance from that median, so one noisy
     *  frame does not count). */
    static double[] place(double[][] shifts, int first, int last) {
        int n = shifts.length, a = first - 1, b = Math.min(n, last); double[] xs = new double[b - a], ys = new double[b - a]
        for (int f = a; f < b; f++) { xs[f - a] = shifts[f][0]; ys[f - a] = shifts[f][1] }
        double cx = median(xs), cy = median(ys)
        double[] far = new double[b - a]; for (int f = a; f < b; f++) far[f - a] = Math.hypot(shifts[f][0] - cx, shifts[f][1] - cy)
        Arrays.sort(far)
        return [cx, cy, far.length >= 2 ? far[far.length - 2] : far[0]] as double[]
    }
    /** Whether the field moved near a rise, for the checks on the calls: at a frame from..to (counted from 1) the field is MOVED
     *  pixels or more from where it was over the F0 frames f0First..f0Last (from where it sits for most of the recording when it
     *  moved within the F0 frames themselves), with the displacement as "displacement" gives it; or, when no frame is, it moved
     *  MOVED pixels or more between two frames from..to, each the middle value of five frames.
     *  @return null when it did not; else first and last (the longest stretch of displaced frames, the first such when two are as
     *          long; or the two closest frames it moved MOVED pixels or more between), runs (every stretch of displaced frames,
     *          [a, b]), largest (pixels: the largest displacement, or the largest distance between two frames) and from: 'F0',
     *          'recording' or 'between' */
    static Map movedNear(double[][] shifts, int f0First, int f0Last, int from, int to) {
        int n = shifts.length
        from = Math.max(1, from); to = Math.min(n, to)
        if (n == 0 || to < from || f0First < 1 || f0Last < f0First || f0First > n) return null
        double[] at = place(shifts, f0First, f0Last)
        boolean steady = at[2] < MOVED
        double[] size = displacement(shifts, steady ? at[0] : 0d, steady ? at[1] : 0d)
        List<List<Integer>> runs = new ArrayList<List<Integer>>()
        int start = -1, best = -1; double largest = 0d
        for (int f = from - 1; f <= to; f++) {
            boolean off = f < to && size[f] >= MOVED
            if (off) { largest = Math.max(largest, size[f]); if (start < 0) start = f }
            else if (start >= 0) { runs.add([start + 1, f] as List<Integer>); start = -1 }
        }
        if (!runs.isEmpty()) {
            List<Integer> longest = runs.get(0)
            for (List<Integer> run : runs) if (run.get(1) - run.get(0) > longest.get(1) - longest.get(0)) longest = run
            return [first: longest.get(0), last: longest.get(1), runs: runs, largest: largest, from: steady ? 'F0' : 'recording'] as Map
        }
        // moved within the frames themselves: the largest distance between two of them
        double[] xs = new double[n], ys = new double[n]
        for (int f = 0; f < n; f++) {
            int a = Math.max(0, f - 2), b = Math.min(n, f + 3); double[] px = new double[b - a], py = new double[b - a]
            for (int k = a; k < b; k++) { px[k - a] = shifts[k][0]; py[k - a] = shifts[k][1] }
            xs[f] = median(px); ys[f] = median(py)
        }
        int bi = -1, bj = -1; double far = 0d
        for (int i = from - 1; i < to; i++) for (int j = i + 1; j < to; j++) {
            double d = Math.hypot(xs[j] - xs[i], ys[j] - ys[i])
            if (d > far) far = d
            if (d >= MOVED && (bi < 0 || j - i < bj - bi)) { bi = i; bj = j }   // the two closest frames it moved that far between
        }
        if (bi < 0) return null
        return [first: bi + 1, last: bj + 1, runs: new ArrayList<List<Integer>>(), largest: far, from: 'between'] as Map
    }
    /** What the shifts mean, counted from the frames refFirst..refLast (the frames the cells were chosen on) when given, else
     *  from where the field sits for most of the recording. "runs": frames displaced by MOVED pixels or more; "slightRuns":
     *  by SLIGHT up to MOVED. Frames are counted from 1, both ends included. A frame's displacement is the middle value of
     *  five frames around it, so single noisy frames do not break a run. */
    static Map summary(double[][] shifts, int refFirst, int refLast) {
        int n = shifts.length; boolean chosen = refFirst >= 1 && refLast >= refFirst && refFirst <= n
        double rx = 0d, ry = 0d, within = 0d
        if (chosen) {
            double[] at = place(shifts, refFirst, refLast)
            within = at[2]
            // a picture made from frames that moved is smeared and has no one place: then count from most of the recording instead
            if (within < MOVED) { rx = at[0]; ry = at[1] }
        }
        boolean steady = !chosen || within < MOVED
        double[] size = displacement(shifts, rx, ry)
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
        return pixel > 0 && CalciumSteps.micrometres(unit) ? pixel : Double.NaN
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
    /** What set a cell's threshold, from the two columns that say it between them. Threshold_source names the term: "Floor"
     *  (the minimum rise) or "Baseline noise" (the SD multiplier times the noise, whichever estimate the noise is).
     *  Noise_source, in a table made with the whole-trace estimator, names that estimate: "Baseline SD", "Whole-trace
     *  noise" or "Noise floor". A table without the column was made on the baseline SD and keeps "Baseline noise". */
    static String thresholdSetBy(Map call) { String term = call.Threshold_source ?: ''; term == 'Baseline noise' && call.Noise_source ? call.Noise_source as String : term }
}

class CalciumAnalysisPlot extends javax.swing.JPanel {
    double[] time=[], mean=[], sem=[], selected=[]
    String caption='Run the analysis to see traces', cell=''
    double baselineStart=0, baselineEnd=0, stimulusStart=0, stimulusEnd=0
    // every cell's trace behind the mean, coloured by its response call; "moment" is the time shown on the picture
    double[][] all=[]; List<String> allCalls=[]; boolean showAll=true; String allLabel='all cells'
    double moment=Double.NaN; Closure seek={double t->}
    Map responseCall=[:]
    /** Under the plot: the selected cell's required rise, what set it, and the run that met it, with the decimals of the call
     *  line above the plot (a span in whole seconds, or to a tenth when it is not whole). As in the call line, a value the
     *  table does not hold is left out with its words; without a required rise there is no line. */
    String ruleLine(){
        def number={String key->CalciumCsv.number(responseCall[key]?.toString())}
        double need=number('Threshold_increase'),span=number('Run_duration_s'),samples=number('Run_samples')
        if(!Double.isFinite(need))return ''
        String setBy=CalciumCalls.thresholdSetBy(responseCall)
        List parts=['Required rise '+String.format(Locale.US,'%.2f',need)+(setBy?' ('+setBy+')':'')]
        if(Double.isFinite(span))parts<<'span '+new java.text.DecimalFormat('0.#',java.text.DecimalFormatSymbols.getInstance(Locale.US)).format(span)+' s'
        if(Double.isFinite(samples))parts<<String.format(Locale.US,'%.0f',samples)+' samples'
        parts.join(' · ')
    }
    /** Under the plot on screen, and never in a saved figure (nothing but the axes and the key is drawn on the plot): the
     *  caption and the rule line, one line each, so they take the same room for every cell. */
    final List<CalciumLine> noteLines=[new CalciumLine(CalciumTheme.MUTED),new CalciumLine(CalciumTheme.MUTED)]
    final javax.swing.JPanel notes=new javax.swing.JPanel()
    /** True only while Save plots paints trace.png: its key then names the dashed lines and the bands too. On screen the key
     *  holds the cell, the mean and the traces behind it, so the plot keeps its height beside a large picture, and pointing at
     *  a band names it. */
    boolean figure=false
    /** The response calls of every cell of the window shown: trace.png's key names the bands the window has for every cell,
     *  and the key keeps the room of the cell whose entries take the most lines, so the plot keeps its place from cell to cell. */
    List<Map> windowCalls=[]
    /** The calls a mark can give the cells of the group shown (a mark moves a cell of the primary stimulus from one call to the
     *  other; the window says which): the key keeps room for a count of each, so a mark never moves the plot. */
    List<String> markCalls=[]
    /** The bands drawn last: [left, right, name], for the tooltip that names the band under the pointer. */
    private List bandPlaces=[]
    private java.awt.image.BufferedImage layer; private String layerKey=''
    private int l=62,t=40,w=1,h=1; private double lo=0,hi=1
    CalciumAnalysisPlot() {
        setPreferredSize(new java.awt.Dimension(520,320)); setMinimumSize(new java.awt.Dimension(240,180)); setBackground(Color.white)
        notes.setLayout(new javax.swing.BoxLayout(notes,javax.swing.BoxLayout.Y_AXIS));notes.setBackground(Color.white);notes.setBorder(javax.swing.BorderFactory.createEmptyBorder(0,l,6,18))
        noteLines.each{notes.add(it)}
        def mouse=new java.awt.event.MouseAdapter(){
            void mousePressed(java.awt.event.MouseEvent e){move(e)}
            void mouseDragged(java.awt.event.MouseEvent e){move(e)}
            void move(java.awt.event.MouseEvent e){if(time.length>1&&w>1&&javax.swing.SwingUtilities.isLeftMouseButton(e))seek.call(time[0]+Math.max(0d,Math.min(1d,(e.x-l)/(double)w))*(time[-1]-time[0]))}
        }
        addMouseListener(mouse);addMouseMotionListener(mouse)
        setToolTipText(HINT)
    }
    static final String HINT='Click or drag to move through time; the picture follows'
    /** The band under the pointer by name, then how to move through time. */
    String getToolTipText(java.awt.event.MouseEvent e){
        boolean inside=e!=null&&(e.y>=t)&&(e.y<=t+h)
        List names=inside?bandPlaces.findAll{(e.x>=(it[0] as int))&&(e.x<(it[1] as int))}.collect{it[2]}:[]
        names?'<html>'+names.join(' · ')+'<br>'+HINT+'</html>':HINT
    }
    /** Call after changing what is drawn (traces, selection, bands); moving the time cursor needs only repaint(). */
    void changed(){layerKey='';showNotes();repaint()}
    /** The caption, and the rule line of a cell with a trace: a cell without one (excluded, or unscorable without a trace) has
     *  only the caption, which says why. */
    void showNotes(){List words=time.length<2 ? ['',''] : [caption,selected.length>0?ruleLine():''];noteLines.eachWithIndex{line,i->line.setText(words[i] as String)}}
    Map referenceValues(Map call=responseCall) {
        double base=CalciumCsv.number(call.Baseline_mean?.toString())
        def values=[baseline:base]
        ['floor':'Floor_threshold','noise':'Noise_threshold','threshold':'Threshold_increase'].each{k,v->
            double fraction=CalciumCsv.number(call[v]?.toString());values[k]=Double.isFinite(base)&&Double.isFinite(fraction)?base+(1d+base)*fraction:Double.NaN
        }
        values
    }
    /** The dashed lines drawn for a call between the start and the end of its response window, in the key's order and words
     *  (the Step 4 page's): [colour, width, name, 'dash']. The required rise is drawn last, over the others; a line at the
     *  value of one named before it is hidden under that one and is not named. */
    List levelsOf(Map call,double start,double end) {
        if(!Double.isFinite(start)||!Double.isFinite(end))return []
        Map values=referenceValues(call);List named=[],at=[]
        [['threshold',CalciumTheme.GOOD,2f,'required rise'],['noise',CalciumTheme.WARN,1f,'noise threshold'],['floor',CalciumTheme.WARN,1f,'minimum rise'],['baseline',CalciumTheme.MUTED,1f,'baseline mean']].each{List level->
            double v=values[level[0]] as double;if(Double.isFinite(v)&&!at.contains(v)){at<<v;named<<[level[1],level[2],level[3],'dash']}}
        named
    }
    /** The bands of a call whose response window runs from start to end, in the order they are painted and in the words of the
     *  Step 4 page's legend: [from, to, colour, name] for the F₀ frames, the baseline just before the stimulus, the response
     *  window and the response (the run that counted). */
    List bandsOf(Map call,double start,double end) {
        double localEnd=CalciumCsv.number(call.Baseline_end_s?.toString()),localStart=CalciumCsv.number(call.Baseline_start_s?.toString())
        if(!Double.isFinite(localEnd))localEnd=start
        if(!Double.isFinite(localStart))localStart=localEnd-CalciumCsv.number(call.Baseline_s?.toString())
        [[baselineStart,baselineEnd,new Color(0xEEF1F4),'F₀ frames'],[localStart,localEnd,new Color(0xDCEBFA),'baseline'],
         [start,end,new Color(0xFFF0DC),'response window'],[CalciumCsv.number(call.Run_start_s?.toString()),CalciumCsv.number(call.Run_end_s?.toString()),new Color(0xC3EBC8),'response']]
    }
    /** Where a call's response window runs, as selectCell puts it on the plot. */
    static double[] windowOf(Map call){[CalciumCsv.number(call.Start_s?.toString()),CalciumCsv.number((call.Effective_end_s?:call.End_s)?.toString())] as double[]}
    /** Whether a band covers some of the time shown. */
    boolean shown(List band){double a=band[0] as double,b=band[1] as double;Double.isFinite(a)&&Double.isFinite(b)&&Math.max(a,time[0])<Math.min(b,time[-1])}
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
        // the key above the plot: the cell, the mean and the traces behind them by call, and in trace.png (figure) the dashed
        // lines and the bands too, each named once and in the order of the Step 4 page's legend; on screen pointing at a band
        // names it. The key takes as many lines as it needs, and the plot starts below it
        def countsOf={int yes,int no,int other->List list=[]
            if(yes>0)list<<[CalciumCalls.RESPONDED,1.4f,'responded ('+yes+')']
            if(no>0)list<<[CalciumCalls.NOT,1.4f,'did not respond ('+no+')']
            if(other>0)list<<[CalciumCalls.OTHER,1.4f,(yes+no>0?'no call':allLabel)+' ('+other+')'];list}
        int yes=many?(allCalls.count{it=='Responded'} as int):0,no=many?(allCalls.count{it=='Did not respond'} as int):0,other=many?all.length-yes-no:0
        List counts=countsOf(yes,no,other)
        // room for a count of each call a mark can bring into the group, as wide as the count of every cell of the window with
        // one of the two calls (a mark changes no other count)
        int most=Math.max(yes+no,windowCalls.count{it.Status=='Responded'||it.Status=='Did not respond'} as int)
        List countRoom=!showAll||!markCalls ? counts : countsOf(yes>0||'Responded' in markCalls?most:0,no>0||'Did not respond' in markCalls?most:0,other)
        // in trace.png a band is named whenever the window has it, for a cell without a run too
        List calls=windowCalls+[responseCall],bands=bandsOf(responseCall,stimulusStart,stimulusEnd)
        List named=!figure ? [] : bands.withIndex().findAll{band,int i->i==0?shown(band):calls.any{Map call->double[] span=windowOf(call);shown(bandsOf(call,span[0],span[1])[i] as List)}}.collect{band,int i->[band[2],0f,band[3],'box']}
        def keyOf={boolean traced,List levels,List countsShown->List key=[];if(traced)key<<[CalciumTheme.INK,1.7f,cell];key<<[CalciumTheme.BLUE,2.4f,'mean ± SEM'];key+countsShown+(figure?levels+named:[])}
        def entries=keyOf(cell&&selected.length>0,levelsOf(responseCall,stimulusStart,stimulusEnd),counts)   // a cell without a trace has no entry: why is said under the plot
        // on screen the key is laid from the frame's left edge across the frame's width, as before the review controls, so the
        // frame starts where it did; trace.png's longer key runs across the whole picture from 8 px in
        int from=figure?8:l,across=figure?Math.max(1,getWidth()-16):w
        def laid={List list->int line=0,run=0;def places=[];list.each{entry->int need=21+fm.stringWidth(entry[2] as String)+14;if(run>0&&run+need>across){line++;run=0};places<<[run,line];run+=need};[places,line]}
        def (places,lines)=laid(entries)
        // the key keeps the room of the cell of this window whose entries take the most lines (its own trace, and in trace.png
        // its dashed lines), with the counts a mark can bring, so the plot does not move from cell to cell or at a mark
        List levelSets=!figure ? [[]] : calls.collect{Map call->double[] span=windowOf(call);levelsOf(call,span[0],span[1])}.unique{List levels->levels.collect{it[2]}.join(' | ')}
        if(cell)levelSets.each{List levels->lines=Math.max(lines as int,laid(keyOf(true,levels,countRoom))[1] as int)}
        t=30+16*((lines as int)+1);h=Math.max(1,getHeight()-t-46)
        def values=[]; mean.eachWithIndex {v,i->if(Double.isFinite(v)) {values<<v;if(i<sem.length&&Double.isFinite(sem[i])) {values<<v-sem[i];values<<v+sem[i]}}}
        selected.each {if(Double.isFinite(it)) values<<it}; values<<0d
        Map references=referenceValues();references.values().each{if(Double.isFinite(it as double))values<<it}
        // With every cell drawn, the axis covers 99 % of all values (plus the mean and the selected cell); rarer extremes run off the edge.
        if(many) {values<<CalciumFrameView.share(all,0.005d,false);values<<CalciumFrameView.share(all,0.995d,false)}
        lo=values.min() as double;hi=values.max() as double
        if(hi-lo<0.001) hi=lo+1
        double pad=(hi-lo)*0.08;lo-=pad;hi+=pad
        double t0=time[0],t1=time[-1]
        // the bands are named in the key, not on the plot; the tooltip names the one under the pointer in the same words
        List drawnBands=[]
        bands.each { band ->
            if(!Double.isFinite(band[0] as double)||!Double.isFinite(band[1] as double))return
            int x0=Math.max(l,X(band[0] as double)),x1=Math.min(l+w,X(band[1] as double))
            if(x1>x0) {g.setColor(band[2] as Color);g.fillRect(x0,t,x1-x0,h);drawnBands<<[x0,x1,band[3]]}
        }
        bandPlaces=drawnBands
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
        // a line for a line, dashed as it is drawn; a filled box for a band; a name longer than the key's width is cut short with …
        entries.eachWithIndex{entry,int i->int kx=from+(places[i][0] as int),ky=20+16*(places[i][1] as int);String text=CalciumTheme.fit(entry[2] as String,fm,across-21),kind=entry.size()>3?entry[3]:'line'
            g.setColor(entry[0] as Color)
            if(kind=='box'){g.fillRect(kx,ky-5,16,10);g.setColor(new Color(0xC5CED8));g.setStroke(new java.awt.BasicStroke(1f));g.drawRect(kx,ky-5,16,10)}
            else{g.setStroke(kind=='dash'?new java.awt.BasicStroke(entry[1] as float,java.awt.BasicStroke.CAP_BUTT,java.awt.BasicStroke.JOIN_MITER,10f,[5f,4f] as float[],0f):new java.awt.BasicStroke(entry[1] as float));g.drawLine(kx,ky,kx+16,ky)}
            g.setStroke(new java.awt.BasicStroke(1f));g.setColor(CalciumTheme.INK);g.drawString(text,kx+21,ky+4)}
        g.setColor(CalciumTheme.MUTED);String xl='Time (s)';g.drawString(xl,l+(int)((w-fm.stringWidth(xl))/2),getHeight()-8)
        def turned=(java.awt.Graphics2D)g.create();try{turned.rotate(-Math.PI/2);String yl='Global ΔF/F₀';turned.drawString(yl,-(t+(int)((h+fm.stringWidth(yl))/2)),14)}finally{turned.dispose()}
    }
}
/** The plot of the Around the cell tab: the cell and the zones around it as the passing-object check used them
 *  (CalciumCallChecks.aroundCell), in gray levels. Each zone the flag rests on has a colour of its own and is named in the key
 *  by its place; every other usable near-ring sector is thin and gray under one entry. The cell's response window is a band,
 *  as on Traces, and the frames the flag names are another. Nothing else is drawn on it; the key is the same on screen and in
 *  around_the_cell.png. Without data the sentence that says why is shown where the plot would be. */
class CalciumAroundPlot extends javax.swing.JPanel {
    /** What aroundCell gave for the cell shown: why is '' when there is something to draw. null: no analysis is shown. */
    Map data
    /** The cell as the key on Traces names it. */
    String cell=''
    /** The cells of the window: the key keeps the room of the widest name, so the plot keeps its place from cell to cell. */
    List<String> names=[]
    /** The cell's response window, as Traces draws it: from Start_s to Effective_end_s of its call. */
    double stimulusStart=Double.NaN,stimulusEnd=Double.NaN
    double moment=Double.NaN; Closure seek={double t->}
    static final String HINT='Click or drag to move through time'
    static final String WAITING='Run the analysis to see traces'
    /** The colour of each part a zone of the flag can be, as on the Step 4 page: before and rim, after and ring, the core. A part
     *  met a second time for one cell takes the next spare colour. */
    static final Map COLOURS=[before:new Color(0x0072B2),rim:new Color(0x0072B2),after:new Color(0xAA3377),ring:new Color(0xAA3377),core:new Color(0x009E73)]
    static final List SPARE=[new Color(0xE69F00),new Color(0x56B4E9),new Color(0xCC79A7),new Color(0x882255)]
    static final Color OTHER=new Color(0xAEB8C2),WINDOW=new Color(0xFFF0DC),FLAGGED=new Color(0xE2D6F3)
    /** The longest names three zones of a flag can have: the key keeps room for them. */
    static final List WIDEST=['far ring, upper right','far ring, upper right','far ring, upper right']
    /** The bands drawn: [left, right, name], for the tooltip that names the band under the pointer. */
    private List bandPlaces=[]
    private java.awt.image.BufferedImage layer; private String layerKey=''
    private int l=62,t=40,w=1,h=1; private double lo=0,hi=1
    CalciumAroundPlot() {
        setPreferredSize(new java.awt.Dimension(520,320)); setMinimumSize(new java.awt.Dimension(240,180)); setBackground(Color.white)
        def mouse=new java.awt.event.MouseAdapter(){
            void mousePressed(java.awt.event.MouseEvent e){move(e)}
            void mouseDragged(java.awt.event.MouseEvent e){move(e)}
            void move(java.awt.event.MouseEvent e){double[] x=times();if(x.length>1&&w>1&&javax.swing.SwingUtilities.isLeftMouseButton(e))seek.call(x[0]+Math.max(0d,Math.min(1d,(e.x-l)/(double)w))*(x[-1]-x[0]))}
        }
        addMouseListener(mouse);addMouseMotionListener(mouse)
        setToolTipText(HINT)
    }
    /** Whether there is a plot to draw: data with times and no why. */
    boolean drawn(){data!=null&&!data.why&&data.time instanceof double[]&&(data.time as double[]).length>1}
    double[] times(){drawn()?data.time as double[]:new double[0]}
    /** What is shown in place of the plot: why, or the sentence before an analysis. */
    String sentence(){data==null?WAITING:(''+(data.why?:''))}
    /** The band under the pointer by name, then how to move through time. */
    String getToolTipText(java.awt.event.MouseEvent e){
        if(!drawn())return null
        boolean inside=e!=null&&(e.y>=t)&&(e.y<=t+h)
        List found=inside?bandPlaces.findAll{(e.x>=(it[0] as int))&&(e.x<(it[1] as int))}.collect{it[2]}:[]
        found?'<html>'+found.join(' · ')+'<br>'+HINT+'</html>':HINT
    }
    /** Call after changing what is drawn; moving the time cursor needs only repaint(). */
    void changed(){layerKey='';repaint()}
    /** The time between two frames, as most frames are apart: the median frame step; of an even number of steps the upper of the
     *  two middle ones, as on the Step 4 page. */
    double step(){double[] x=times();if(x.length<2)return 1d;List d=(1..<x.length).collect{x[it]-x[it-1]}.sort();d[(int)(d.size()/2)] as double}
    /** The lines under the cell, in the order they are drawn, each [colour, width, name, trace]: every other usable near-ring
     *  sector thin and gray, then the zones the flag rests on, each in the colour of its part and named by its place. */
    List lines() {
        if(!drawn())return []
        Map byZone=((data.zones?:[]) as List).collectEntries{[(it.zone as int),it]}
        List out=[],met=[],spares=[]
        ((data.others?:[]) as List).each{def z->double[] y=byZone[z as int]?.trace as double[];if(y!=null)out<<[OTHER,0.9f,'other ring sectors',y]}
        ((data.lines?:[]) as List).each{Map line->double[] y=byZone[line.zone as int]?.trace as double[];if(y==null)return
            Color c=met.contains(line.part)?null:COLOURS[line.part] as Color;met<<line.part
            if(c==null){c=(SPARE.find{!spares.contains(it)}?:CalciumTheme.MUTED) as Color;spares<<c}
            out<<[c,1.5f,''+(line.place?:CalciumCallChecks.place(line.zone as int)),y]}
        out
    }
    /** The bands, [from, to, colour, name] in seconds: the cell's response window as on Traces, then the frames the flag names,
     *  from the first one's time less half the median frame step to the last one's time plus half of it, as on the page. */
    List bands() {
        double[] x=times();List out=[];if(x.length<2)return out
        if(Double.isFinite(stimulusStart)&&Double.isFinite(stimulusEnd)&&Math.max(stimulusStart,x[0])<Math.min(stimulusEnd,x[-1]))out<<[stimulusStart,stimulusEnd,WINDOW,'response window']
        List band=data.band instanceof List?data.band as List:null
        if(band!=null&&band.size()==2&&band.every{it instanceof Number}){int a=band[0] as int,b=band[1] as int;double half=step()/2d
            if(a>=1&&b>=a&&b<=x.length)out<<[x[a-1]-half,x[b-1]+half,FLAGGED,'flagged frames']}
        out
    }
    /** The key, in the order it is drawn: the cell, the zones of the flag, the other ring sectors, then the bands. */
    List entries(List shown) {
        List key=[[CalciumTheme.INK,1.7f,cell,'line']]
        shown.findAll{it[2]!='other ring sectors'}.each{key<<[it[0],it[1],it[2],'line']}
        if(shown.any{it[2]=='other ring sectors'})key<<[OTHER,0.9f,'other ring sectors','line']
        bands().each{key<<[it[2],0f,it[3],'box']}
        key
    }
    int X(double v){double[] x=times();l+(int)Math.round((v-x[0])/(x[-1]-x[0])*w)}
    int Y(double v){t+(int)Math.round((hi-v)/(hi-lo)*h)}
    protected void paintComponent(java.awt.Graphics graphics) {
        super.paintComponent(graphics)
        def g=(java.awt.Graphics2D)graphics.create()
        try {
            if(!drawn()){still(g);return}
            // everything but the time cursor is drawn once into a picture and reused, so dragging through time stays smooth
            double scale=Math.max(1d,g.getTransform().getScaleX())
            String key=getWidth()+'x'+getHeight()+'@'+scale
            if(layer==null||layerKey!=key) {
                layer=new java.awt.image.BufferedImage((int)Math.ceil(getWidth()*scale),(int)Math.ceil(getHeight()*scale),java.awt.image.BufferedImage.TYPE_INT_RGB)
                def lg=layer.createGraphics();try{lg.scale(scale,scale);lg.setColor(Color.white);lg.fillRect(0,0,getWidth(),getHeight());still(lg)}finally{lg.dispose()}
                layerKey=key
            }
            g.drawImage(layer,0,0,getWidth(),getHeight(),null)
            double[] x=times()
            if(Double.isFinite(moment)&&moment>=x[0]&&moment<=x[-1]) {
                CalciumTheme.smooth(g);g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f));def fm=g.getFontMetrics()
                int cx=X(moment);g.setColor(new Color(31,41,51,200));g.setStroke(new java.awt.BasicStroke(1.2f));g.drawLine(cx,t,cx,t+h)
                int i=0;double best=Double.MAX_VALUE;for(int k=0;k<x.length;k++){double d=Math.abs(x[k]-moment);if(d<best){best=d;i=k}}
                String s=String.format(Locale.US,'%.1f s',x[i]);double[] rise=data.rise instanceof double[]?data.rise as double[]:new double[0]
                if(i<rise.length&&Double.isFinite(rise[i])){int y=Math.max(t,Math.min(t+h,Y(rise[i])));g.setColor(Color.white);g.fillOval(cx-5,y-5,10,10);g.setColor(CalciumTheme.INK);g.fillOval(cx-3,y-3,7,7)
                    s+=String.format(Locale.US,' · %s %.1f',cell,rise[i])}
                int sw=fm.stringWidth(s)+10,sx=Math.max(l,Math.min(l+w-sw,cx-(int)(sw/2)))
                g.setColor(CalciumTheme.INK);g.fillRoundRect(sx,t+h-fm.getHeight()-6,sw,fm.getHeight()+4,8,8);g.setColor(Color.white);g.drawString(s,sx+5,t+h-8)
            }
        } finally {g.dispose()}
    }
    /** Words in lines that fit a width. */
    static List wrapped(String text,java.awt.FontMetrics fm,int width) {
        List out=[];String line=''
        (text?:'').split(' ').each{String word->String next=line?line+' '+word:word;if(line&&fm.stringWidth(next)>width){out<<line;line=word}else line=next}
        if(line)out<<line;out
    }
    /** Bands, grid, the other ring sectors, the zones of the flag, the cell, the key and the axis names; without data, why. */
    void still(java.awt.Graphics2D g) {
        CalciumTheme.smooth(g); g.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f)); def fm=g.getFontMetrics()
        if(!drawn()) {
            bandPlaces=[];g.setColor(CalciumTheme.MUTED)
            List said=wrapped(sentence(),fm,Math.max(40,getWidth()-48));int top=(int)((getHeight()-said.size()*fm.getHeight())/2)+fm.getAscent()
            said.eachWithIndex{String s,int i->g.drawString(s,(int)((getWidth()-fm.stringWidth(s))/2),top+i*fm.getHeight())}
            return
        }
        double[] x=times(),rise=data.rise instanceof double[]?data.rise as double[]:new double[0]
        l=62;w=Math.max(1,getWidth()-l-18)
        List shown=lines(),entries=entries(shown)
        // the key above the plot, laid from the frame's left edge across its width; it keeps the room of the most it can hold
        // (the widest cell name, three zones with the longest names, the other ring sectors and both bands), so the plot keeps
        // its place from cell to cell
        def laid={List list->int line=0,run=0;def places=[];list.each{entry->int need=21+fm.stringWidth(entry[2] as String)+14;if(run>0&&run+need>w){line++;run=0};places<<[run,line];run+=need};[places,line]}
        def (places,lines)=laid(entries)
        String widest=((names?:[])+[cell]).max{fm.stringWidth((it?:'') as String)}?:''
        int most=laid([[null,0f,widest]]+WIDEST.collect{[null,0f,it]}+[[null,0f,'other ring sectors'],[null,0f,'response window'],[null,0f,'flagged frames']])[1] as int
        t=30+16*(Math.max(lines as int,most)+1);h=Math.max(1,getHeight()-t-46)
        List values=[0d];([rise]+shown.collect{it[3]}).each{double[] y->for(double v:y)if(Double.isFinite(v))values<<v}
        lo=values.min() as double;hi=values.max() as double
        if(hi-lo<0.001)hi=lo+1
        double pad=(hi-lo)*0.08;lo-=pad;hi+=pad
        double t0=x[0],t1=x[-1]
        // the bands are named in the key, not on the plot; the tooltip names the one under the pointer in the same words
        List drawnBands=[]
        bands().each{band->int x0=Math.max(l,X(band[0] as double)),x1=Math.min(l+w,X(band[1] as double))
            if(x1>x0){g.setColor(band[2] as Color);g.fillRect(x0,t,x1-x0,h);drawnBands<<[x0,x1,band[3]]}}
        bandPlaces=drawnBands
        def yTicks=CalciumTheme.ticks(lo,hi,5); double yStep=yTicks.size()>1?yTicks[1]-yTicks[0]:1d
        yTicks.each { double v -> int y=Y(v);g.setColor(v==0d?new Color(0xAEB8C2):CalciumTheme.LINE);g.drawLine(l,y,l+w,y)
            String s=CalciumTheme.tick(v,yStep);g.setColor(CalciumTheme.MUTED);g.drawString(s,l-8-fm.stringWidth(s),y+(int)(fm.getAscent()/2)-1) }
        def xTicks=CalciumTheme.ticks(t0,t1,Math.max(2,(int)(w/90))); double xStep=xTicks.size()>1?xTicks[1]-xTicks[0]:1d
        xTicks.each { double v -> int px=X(v);g.setColor(CalciumTheme.LINE);g.drawLine(px,t+h,px,t+h+4)
            String s=CalciumTheme.tick(v,xStep);g.setColor(CalciumTheme.MUTED);g.drawString(s,px-(int)(fm.stringWidth(s)/2),t+h+6+fm.getAscent()) }
        g.setColor(new Color(0xC5CED8));g.drawRect(l,t,w,h)
        def clip=g.getClip();g.clipRect(l,t,w+1,h+1)
        // a line breaks where the frames are further apart than one and a half times as most are
        double gap=1.5d*step()
        def path={double[] y->def p=new java.awt.geom.Path2D.Double();boolean pen=false
            for(int i=0;i<Math.min(x.length,y.length);i++){if(i>0&&x[i]-x[i-1]>gap)pen=false;if(Double.isFinite(y[i])){double px=l+(x[i]-t0)/(t1-t0)*w,py=t+(hi-y[i])/(hi-lo)*h;if(pen)p.lineTo(px,py) else p.moveTo(px,py);pen=true} else pen=false};p}
        def round={float width->new java.awt.BasicStroke(width,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND)}
        shown.findAll{it[2]=='other ring sectors'}.each{g.setColor(it[0] as Color);g.setStroke(round(it[1] as float));g.draw(path(it[3] as double[]))}
        shown.findAll{it[2]!='other ring sectors'}.each{g.setColor(Color.white);g.setStroke(round(3.2f));g.draw(path(it[3] as double[]));g.setColor(it[0] as Color);g.setStroke(round(it[1] as float));g.draw(path(it[3] as double[]))}
        if(rise.length>0){g.setColor(Color.white);g.setStroke(round(3.6f));g.draw(path(rise));g.setColor(CalciumTheme.INK);g.setStroke(round(1.7f));g.draw(path(rise))}
        g.setClip(clip);g.setStroke(new java.awt.BasicStroke(1f))
        // the key: a line for a line, a filled box for a band; a name longer than the key's width is cut short with …
        entries.eachWithIndex{entry,int i->int kx=l+(places[i][0] as int),ky=20+16*(places[i][1] as int);String text=CalciumTheme.fit(entry[2] as String,fm,w-21)
            g.setColor(entry[0] as Color)
            if(entry[3]=='box'){g.fillRect(kx,ky-5,16,10);g.setColor(new Color(0xC5CED8));g.drawRect(kx,ky-5,16,10)}
            else{g.setStroke(new java.awt.BasicStroke(entry[1] as float));g.drawLine(kx,ky,kx+16,ky)}
            g.setStroke(new java.awt.BasicStroke(1f));g.setColor(CalciumTheme.INK);g.drawString(text,kx+21,ky+4)}
        g.setColor(CalciumTheme.MUTED);String xl='Time (s)';g.drawString(xl,l+(int)((w-fm.stringWidth(xl))/2),getHeight()-8)
        def turned=(java.awt.Graphics2D)g.create();try{turned.rotate(-Math.PI/2);String yl='ΔF (gray levels)';turned.drawString(yl,-(t+(int)((h+fm.stringWidth(yl))/2)),14)}finally{turned.dispose()}
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
    Set<String> only=null                 // when set, only these cells are drawn, whatever their call (the Marked entry of the Show list)
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
                if(only!=null?!only.contains(i<names.size()?names[i]:''):!calls.isEmpty()&&!CalciumCalls.inGroup(call,group))return
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
    javax.swing.JLabel value=CalciumTheme.label('–',CalciumTheme.INK,java.awt.Font.BOLD,7f)
    /** Wraps at the tile's width, so a long caption ("ATP: 1 of 143 evaluable cells after review; 3 before review") is read whole. */
    CalciumText caption=new CalciumText('&nbsp;',CalciumTheme.MUTED)
    /** The width the tile asks for, when set; its height follows the caption at the width the tile is given. */
    int wide=0
    Closure pressed=null; boolean active=false; private boolean hover=false
    CalciumStat(String text,Color colour=CalciumTheme.INK) {
        setLayout(new javax.swing.BoxLayout(this,javax.swing.BoxLayout.Y_AXIS));setOpaque(false);caption.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f))
        setBorder(javax.swing.BorderFactory.createEmptyBorder(8,12,8,12));value.setForeground(colour);caption.say(CalciumOutlineProgress.escape(text));add(value);add(caption)
        addMouseListener(new java.awt.event.MouseAdapter(){
            void mouseEntered(java.awt.event.MouseEvent e){hover=true;repaint()}
            void mouseExited(java.awt.event.MouseEvent e){hover=false;repaint()}
            void mouseReleased(java.awt.event.MouseEvent e){if(pressed!=null&&javax.swing.SwingUtilities.isLeftMouseButton(e))pressed.call()}
        })
    }
    void show(String number,String text){value.setText(number);caption.say(CalciumOutlineProgress.escape(text))}
    java.awt.Dimension getPreferredSize(){def d=super.getPreferredSize();wide>0?new java.awt.Dimension(wide,(int)d.height):d}
    java.awt.Dimension getMinimumSize(){getPreferredSize()}
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
    /** What a cell's noise is taken from (SD_estimator in a version 2 protocol), each with the earliest engine release that
     *  defines it. A protocol records that release, not the engine that wrote it or runs it, so a rule keeps its name, its
     *  hash and its results when the engine is updated. The engine refuses any other pairing (responder_engine_estimators and
     *  responder_validate_protocol in the pinned R source; tests/test_estimator.groovy holds this table to that source). */
    static final Map<String,String> ESTIMATOR_ENGINE=[sample_sd:'2.0.0',larger_of_baseline_sd_and_trace_noise:'2.1.0']
    /** The same in a few words, for the window, the page and messages. */
    static final Map<String,String> ESTIMATOR_NAME=[sample_sd:'baseline SD',larger_of_baseline_sd_and_trace_noise:'larger of baseline SD and whole-trace noise']
    /** The estimator of a protocol made here. This is the one place that decides it; a protocol that was saved or imported
     *  keeps the estimator it names. */
    static String newEstimator(){'larger_of_baseline_sd_and_trace_noise'}
    /** The legacy rule has one noise only. */
    static final String LEGACY_ESTIMATOR='sample_sd'
    /** What a protocol's calls take noise from: its own estimator for version 2 (a protocol that was read or written here names
     *  it at the top), the baseline SD for the legacy rule. */
    static String noiseOf(Map protocol){protocol.schema_version==2?protocol.sd_estimator as String:LEGACY_ESTIMATOR}
    /** The Step 4 page writes this engine version and SD_estimator into every version 2 protocol it sends, whatever it was
     *  given, and leaves sd_estimator at the top alone (web/step4.js, collectDraft). Three things here exist for that reason
     *  only: the page argument of estimator, forPage and exchange. Drop them when the page sends the copies it was given. */
    static final String PAGE_ENGINE='2.0.0',PAGE_ESTIMATOR='sample_sd'
    /** The estimator a version 2 protocol names: sd_estimator at the top when given, otherwise the one its windows share (a
     *  file from SimpleCa, or one this plugin wrote before it carried engine 2.1.0, names it there only). The copies in the
     *  windows and the engine version must agree with it.
     *  @param page the protocol is a draft from the Step 4 page, or may be one: the name at the top then stands over the
     *  page's fixed copies. A protocol file is read without it, so a file that contradicts itself is refused. */
    static String estimator(Map protocol,boolean page=false) {
        Set copies=(protocol.windows as List).collect{(it as Map).SD_estimator} as Set
        Object top=protocol.sd_estimator
        if(top==null){if(copies.size()!=1)throw new IOException('All windows must share one SD_estimator.');top=copies.first()}
        if(!(top instanceof String)||!ESTIMATOR_ENGINE.containsKey(top))throw new IOException('Unsupported or missing SD_estimator.')
        boolean agree=copies==([top] as Set)&&protocol.engine_version==ESTIMATOR_ENGINE[top]
        boolean pageCopies=page&&protocol.sd_estimator!=null&&copies==([PAGE_ESTIMATOR] as Set)&&protocol.engine_version==PAGE_ENGINE
        if(!agree&&!pageCopies)throw new IOException('The noise estimator and the engine version disagree: '+top+' is recorded with engine version '+ESTIMATOR_ENGINE[top]+', but this protocol has SD_estimator '+copies.join(', ')+' with engine version '+protocol.engine_version+'.')
        top as String
    }
    /** The estimator of a draft from the page, or null when the draft is the legacy rule or cannot be read yet. */
    static String draftEstimator(Object protocol){try{protocol instanceof Map&&protocol.schema_version==2?estimator(protocol as Map,true):null}catch(Exception unfinished){null}}
    /** A protocol in the words the page will send it back in, so that the page can tell an untouched draft from a changed
     *  one: the estimator at the top, and the page's fixed copies. A legacy protocol is given the estimator of a new protocol:
     *  the page carries it over when the rule type is changed, so the version 2 rule made there is a protocol made here.
     *  An unfinished draft goes back as it came. Changes and returns the given map.
     *  @param estimator what a version 2 rule is to name from here on; left out, the rule keeps the one it names */
    static Map forPage(Map protocol,String estimator=null) {
        if(protocol.schema_version==1){protocol.sd_estimator=newEstimator();return protocol}
        String named=estimator?:draftEstimator(protocol);if(named==null)return protocol
        protocol.sd_estimator=named;protocol.engine_version=PAGE_ENGINE
        (protocol.windows as List).each{Map w->w.Engine_version=PAGE_ENGINE;w.SD_estimator=PAGE_ESTIMATOR}
        protocol
    }
    /** The page's one door for both Import protocol JSON (a file comes in, the answer goes into its draft) and Export protocol
     *  (its draft comes in, the answer is the file). Both get the protocol as a file would hold it. The page cannot say which
     *  it is asking for, so: a legacy protocol that names an estimator is the page's own draft (forPage) and the answer is the
     *  file, without the name; one that names none is a file, and the answer carries the name for the page's draft. */
    static Map exchange(Object sent) {
        Map out=canonical(sent,true)
        out.schema_version==1&&!(sent instanceof Map&&sent.containsKey('sd_estimator'))?forPage(out):out
    }
    /** For a status line: says what noise a draft uses when that is not what a new protocol takes. Such a draft came from
     *  remembered settings or an imported file. Empty otherwise, and for the legacy rule, which the page labels itself. */
    static String noiseNote(Object protocol){String used=draftEstimator(protocol);used==null||used==newEstimator()?'':' Noise: '+ESTIMATOR_NAME[used]+', as saved.'}
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
    /** @param estimator the estimator of the protocol being written again (a saved or imported one); leave it out for a
     *  protocol made here, which gets newEstimator(). The legacy rule has no choice of estimator. */
    static Map payload(List<Map> windows,int version,String name,String validation,List overrides=null,String estimator=null) {
        if(!(version in [1,2])||windows.size()<1||windows.size()>2)throw new IllegalArgumentException('Use one primary stimulus and at most one positive control.')
        String noise=estimator?:newEstimator()
        if(version==2&&!ESTIMATOR_ENGINE.containsKey(noise))throw new IllegalArgumentException('Unsupported noise estimator: '+noise+'.')
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
                row.putAll([Protocol_version:2,Protocol_name:name,Engine_version:ESTIMATOR_ENGINE[noise],Validation_status:validation,Signal_reference:'fixed_reference_fractional_change',SD_estimator:noise,Threshold_operator:'>=',Coverage_policy:'complete_windows',Response_direction:'increase',Smoothing:'none'])
                out<<row
            }
        }
        if(windows.size()==2)crossCheck(windows.collect{Map w->[name:w.Stimulus?.toString()?.trim(),baselineStart:numeric(w.Baseline_start_s,'Baseline_start_s'),baselineEnd:numeric(w.Baseline_end_s,'Baseline_end_s'),
            onset:numeric(w.Stimulus_onset_s,'Stimulus_onset_s'),end:numeric(w.End_s,'End_s')]})
        if(!name?.trim()||!(validation in ['provisional','user_defined','validated_by_user']))throw new IllegalArgumentException('Enter a protocol name and validation status.')
        if(version==1&&overrides)throw new IllegalArgumentException('Legacy protocols cannot contain field timing overrides.')
        Map written=[schema:'simpleca.response-protocol',schema_version:version,engine_version:version==2?ESTIMATOR_ENGINE[noise]:'legacy-1']
        // sd_estimator at the top is this plugin's addition to the SimpleCa file: SimpleCa reads the copies in the windows and passes over it.
        if(version==2)written.sd_estimator=noise
        written.putAll([name:version==2?name:'Legacy consecutive-sample protocol',validation_status:version==2?validation:'legacy_unvalidated',signal_reference:version==2?'fixed_reference_fractional_change':'fixed_reference_fractional_dff0',windows:out,field_overrides:overrides?:null])
        written
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
    /** @param page see estimator */
    static Map parse(Object input,boolean page=false) {
        def value
        try { value=new groovy.json.JsonSlurper().parseText(input instanceof String?input:groovy.json.JsonOutput.toJson(input)) }
        catch(groovy.json.JsonException bad) { throw new IOException('This protocol file is not valid JSON. Choose a protocol file.') }
        if(!(value instanceof Map)||value.schema!='simpleca.response-protocol'||!(value.schema_version in [1,2])||!(value.windows instanceof List)||value.windows.size()<1||value.windows.size()>2)throw new IOException('Unsupported protocol schema. Use a version 1 or 2 protocol file with one stimulus and an optional control.')
        if(value.signal_reference!=(value.schema_version==2?'fixed_reference_fractional_change':'fixed_reference_fractional_dff0'))throw new IOException('Unsupported signal reference.')
        if(value.schema_version==2?!(value.engine_version in ESTIMATOR_ENGINE.values()):value.engine_version!='legacy-1')throw new IOException('Unsupported protocol engine version.')
        value.windows.eachWithIndex {row,i->
            if(!(row instanceof Map)||row.Role!=(i==0?'Stimulus':'Control'))throw new IOException('Protocol roles must be Stimulus then optional Control.')
            if(value.schema_version==2) {
                [Signal_reference:'fixed_reference_fractional_change',Threshold_operator:'>=',Coverage_policy:'complete_windows',Response_direction:'increase',Smoothing:'none'].each{k,v->if(row[k]!=v)throw new IOException('Unsupported or missing '+k+'.')}
                if(!ESTIMATOR_ENGINE.containsKey(row.SD_estimator))throw new IOException('Unsupported or missing SD_estimator.')
                if(row.Protocol_version!=2)throw new IOException('Protocol versions disagree.')
                [Engine_version:'engine_version',Protocol_name:'name',Validation_status:'validation_status'].each{k,v->if(row[k]!=value[v])throw new IOException(k=='Protocol_name'?'The protocol name at the top of this file ("'+value[v]+'") differs from the name stored with its '+row.Role+' window ("'+row[k]+'"). Change both to match, or import the file unchanged and rename the protocol on the page.':k=='Validation_status'?'The validation status at the top of this file differs from the one stored with its '+row.Role+' window. Change both to match, or import it unchanged and set the status on the page.':'Protocol metadata disagree: '+k)}
                double baseline=numeric(row.Baseline_s,'Baseline_s'),duration=numeric(row.Baseline_end_s,'Baseline_end_s')-numeric(row.Baseline_start_s,'Baseline_start_s')
                if(Math.abs(baseline-duration)>16d*Math.ulp(1d)*Math.max(1d,Math.abs(baseline)))throw new IOException('Baseline_s must match the explicit baseline start and end.')
            } else {
                row.Baseline_start_s=(row.Start_s as double)-(row.Baseline_s as double);row.Baseline_end_s=row.Start_s;row.Stimulus_onset_s=row.Start_s
                row.Duration_mode='samples';row.Min_duration_s=1d;row.Min_baseline_samples=3;row.Noise_floor_local_fraction=0d;row.Expected_interval_s=null
            }
        }
        // From here on a version 2 protocol names its estimator at the top, wherever the file named it.
        if(value.schema_version==2)value.sd_estimator=estimator(value as Map,page)
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
                payload(value.windows.collect{it.Role==override.Role?combined:it} as List<Map>,2,value.name.toString(),value.validation_status.toString(),null,value.sd_estimator as String)
            }
        }
        // Validate supported fields without executing content. The R import validates the full contract again.
        rewritten(value as Map,(value.name?:'Imported protocol').toString())
        value as Map
    }
    /** A protocol that was read, written again with what it names: its rule type, validation status, overrides and estimator. */
    static Map rewritten(Map value,String name) {
        boolean legacy=value.schema_version==1
        payload(value.windows as List<Map>,value.schema_version as int,name,legacy?'user_defined':value.validation_status.toString(),value.field_overrides as List,legacy?null:value.sd_estimator as String)
    }
    static Map canonical(Object input,boolean page=false) {Map value=parse(input,page);rewritten(value,value.name.toString())}
    static String json(Map payload) {groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(payload))}
}

/** A cell of a results table: when its column cuts its text, the whole text is its tooltip. */
class CalciumCutCell extends javax.swing.table.DefaultTableCellRenderer {
    java.awt.Component getTableCellRendererComponent(javax.swing.JTable table,Object value,boolean selected,boolean focus,int row,int column) {
        def shown=super.getTableCellRendererComponent(table,value,selected,focus,row,column)
        String text=value==null?'':''+value;java.awt.Insets room=getInsets()
        setToolTipText(text&&getFontMetrics(getFont()).stringWidth(text)>table.getCellRect(row,column,false).width-room.left-room.right?text:null)
        shown
    }
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
    /** What noise is taken from in the protocol being edited: the plugin's own for a protocol made here, a loaded protocol's own. */
    String estimator=CalciumResponseProtocol.newEstimator()
    def noiseNote=CalciumTheme.label(' ',CalciumTheme.MUTED)
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
    def preset=new javax.swing.JComboBox<String>(['Provisional slow evoked: 2 samples','Custom version 2','Legacy sample rule (starts at 2)'] as String[])
    def validation=new javax.swing.JComboBox<String>(['provisional','user_defined','validated_by_user'] as String[])
    def settingsButton=CalciumTheme.button('Hide settings')
    def previousCell=CalciumTheme.button('‹'), nextCell=CalciumTheme.button('›')
    def status=new javax.swing.JTextArea('Check acquisition timing and a protocol before analysis.',4,24)
    def appliedNote=new javax.swing.JTextArea('No completed analysis in this window.',3,24)
    def timingNote=new javax.swing.JLabel('Time is reconstructed from confirmed frames per second unless timestamps are loaded.')
    def diagnosticNote=new javax.swing.JLabel('Optional diagnostics are saved separately; they never replace the applied responder counts.')
    def progress=new javax.swing.JProgressBar(0,100)
    def confirmed=new javax.swing.JCheckBox('I checked the timing, baseline windows and cell selection',false)
    def kcl=new javax.swing.JCheckBox('Include a positive control',false)
    /** At the end of each analysis, a flagged responder is marked as not responding (CalciumReview.automatic). Not part of the
     *  protocol: the engine never sees it. */
    def automaticBox=new javax.swing.JCheckBox('Mark possible passing objects as not responders',false)
    static final String AUTOMATIC_TIP='After each analysis, a responder flagged as a possible passing object or as a rise on the edge is marked Did not respond, with the reason Passing object. A mark can be removed.'
    def cells=new javax.swing.JComboBox<String>()
    def stimuli=new javax.swing.JComboBox<String>()
    def group=new javax.swing.JComboBox<String>(['All cells'] as String[])
    def allTraces=new javax.swing.JCheckBox('Every trace',true)
    /** The selected cell's call line; it wraps at the width of the Traces tab, so a marked cell's line is read whole. */
    def cellStatus=new CalciumText(' ',CalciumTheme.MUTED).with{lines=1;it}
    /** The Around the cell tab: the call line and the flags of the cell selected on Traces (showCall, fillFlags), and under them
     *  the plot of what the passing-object check used around it. */
    def aroundStatus=new CalciumText(' ',CalciumTheme.MUTED).with{lines=1;it}
    def aroundFlags=CalciumTheme.row(),aroundBox=new javax.swing.JPanel()
    def aroundPlot=new CalciumAroundPlot()
    /** Under the call line and the flags of the Around the cell tab: what they leave of the room kept for them (aroundKept), so
     *  the plot keeps its place from cell to cell. */
    def aroundRoom=new CalciumRoom({->results==null?0:aroundKept()-([aroundStatus,aroundFlags].findAll{it.visible}.collect{(int)it.preferredSize.height}.sum(0) as int)})
    /** What is shown around the cells of the analysis shown (CalciumCallChecks.aroundContext), read once for it, off the event
     *  thread; aroundFolder is the analysis it was read for. null while it is read. */
    volatile Map aroundContext; volatile String aroundFolder
    /** The reads asked for: one that ends after another was asked for is not used. */
    int aroundAsked=0
    /** Runs a read of aroundContext, off the event thread (a test can run it at once). */
    Closure aroundLater={Closure read->Thread.startDaemon('Around the cell'){read.call()}}
    /** What is shown around one cell (CalciumCallChecks.aroundCell), with the window's own checks. */
    Closure aroundOf={Map context,String cell->CalciumCallChecks.aroundCell(context,checks,cell)}
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
    /** The review marks of the analysis shown (CalciumReview.read): they change the call of the primary stimulus. */
    Map reviewMarks=[:]
    /** Why the review file beside the analysis shown is not applied (CalciumReview.read), or ''. There are no mark controls then. */
    String reviewProblem=''
    /** The tally of the analysis shown with its review marks (CalciumAnalysisSettings.tally), as summarize worked it out. */
    Map reviewTally
    /** The Cell response calls and Responder counts tables as saved (table), before any column is added for the review marks. */
    javax.swing.table.TableModel callsSaved, countsSaved
    /** Asks for the keyboard focus; a window that is never shown (in a test) can see what was asked. */
    Closure takeFocus={java.awt.Component c->c.requestFocusInWindow()}
    /** The component that has the keyboard focus; a window that is never shown (in a test) can say which one would. */
    Closure focused={->java.awt.KeyboardFocusManager.currentKeyboardFocusManager.focusOwner}
    /** What the last analysis run here carried over from the review of an earlier one (CalciumReview.carry). */
    volatile Map carried
    /** What the checks of the calls gave for the last analysis run here: folder and found (CalciumCallChecks.run). Kept, as the
     *  Step 4 page keeps it, because a run that could not check the calls saves no table and its note would be lost. */
    volatile Map checkedRun
    /** The automatic marks the last analysis run here made (CalciumReview.automatic): folder, cells, names and the finish
     *  message's part; null when the option was off. */
    volatile Map markedRun
    /** The checks of the analysis shown, as the Step 4 page is given them (CalciumCallChecks.answer); null when there are none. */
    Map checks
    /** The automatic call of every cell for the primary stimulus (CalciumReview.scored). */
    Map automatic=[:]
    /** The parts of the finish message of the last analysis run here (finished), so that it follows the review marks. */
    Map finish
    // Under the selected cell's call line: its flags, then its mark controls, one button or the form that sets a mark.
    def flagRow=CalciumTheme.row(),markRow=CalciumTheme.row(),reviewBox=new javax.swing.JPanel()
    /** Under the call line, its flags, the mark controls and a refusal, each directly under the one before: what they leave of
     *  the room kept for them while results are shown (headRoom), so the head has one height for every cell and after a mark,
     *  and the picture and the plot under it keep their size. */
    def headRoom=new CalciumRoom({->results==null?0:headKept()-([cellStatus,flagRow,markRow,markProblem].findAll{it.visible}.collect{(int)it.preferredSize.height}.sum(0) as int)})
    def markButton=CalciumTheme.button('Mark as not a responder'),removeMark=CalciumTheme.button('Remove mark'),saveMark=CalciumTheme.button('Save',true),cancelMark=CalciumTheme.button('Cancel')
    def reasonLabel=CalciumTheme.label('Reason'),noteLabel=CalciumTheme.label('Note')
    def reasonBox=new javax.swing.JComboBox<String>(CalciumReview.REASONS as String[])
    def noteField=new CalciumHintField('Optional',CalciumReview.NOTE_MAX,16)
    /** A mark refused, in the words of CalciumReview, shown under the form until the cell changes; or why the review file is not
     *  applied. */
    def markProblem=new CalciumText('&nbsp;',CalciumTheme.WARN)
    /** The cell whose mark form is open, or null; the cell the refusal under the form is about, and the words said. */
    String formCell, problemCell, problemSaid
    /** Enter was pressed on the mark controls and not let go yet: it acts when it is let go, once however long it is held. */
    boolean enterPressed=false
    /** Enter went to the open Reason list, which chooses with it and closes: while the key is held its repeats do not ready
     *  the mark controls. */
    boolean enterTaken=false
    /** Whether the Reason list is open; a window that is never shown (in a test) can say it is. */
    Closure listOpen={->reasonBox.popupVisible}
    /** The reason already chosen when the form opens on a flagged responder, by the code of its flag. */
    static final Map REASON_OF_CODE=[POSSIBLE_PASSING_OBJECT:'Passing object',RISE_ON_EDGE:'Passing object',FIELD_MOVED:'Movement',SHORT_SEPARATE_CROSSINGS:'Noise']
    /** How the engine joins the codes in QC_flags of Cell quality and of the whole-recording metrics. */
    static final String QUALITY_JOIN=';'
    /** The last entry of the Show list while the primary stimulus has marks: the cells with a mark. */
    static final int MARKED=4
    // the movie beside the traces
    ImagePlus movie; java.awt.image.BufferedImage snapshot, frameImage; float[] rest, work, work2; double[] rawRange; CalciumRise riseView; int[] riseShade
    String recordingPrint=''   // CalciumSteps.printOf of the attached recording
    /** The measurement with the recording as this window reads frames from it (CalciumSteps.live): a recording closed in Fiji
     *  and opened again is the one the call checks measure around the cells. As it was when there is none. */
    Map withRecording(){ImagePlus now=CalciumSteps.live(movie,recordingPrint);now!=null?measurement+[movie:now,recordingPrint:recordingPrint]:measurement}
    int frame=1, restFirst=1, restLast=1
    def view=new javax.swing.JComboBox<String>(['Static colour overview'] as String[])
    def outlinesBox=new javax.swing.JCheckBox('All outlines',true)
    def zoomButton=CalciumTheme.button('Zoom to cell'), wholeButton=CalciumTheme.button('Whole picture')
    def playButton=CalciumTheme.button('▶ Play')
    def slider=new javax.swing.JSlider(1,2,1)
    def timeLabel=CalciumTheme.label(' ',CalciumTheme.MUTED)
    javax.swing.Timer clock, player, inputWatch
    private boolean listing=false
    /** The cell just marked or unmarked whose call left the group Show lists: it stays selected, and in the Cell list in its
     *  place, until another cell is selected (letGo), so the arrows go on through the group from it. It is not drawn with the
     *  group. null: none. */
    String kept

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
        identity.body.add(CalciumTheme.label('Imported settings remain editable for each stimulus and control.',CalciumTheme.MUTED));identity.body.add(noiseNote)
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
        def marksCard=new CalciumCard('','Review marks',false);automaticBox.toolTipText=AUTOMATIC_TIP
        marksCard.body.add(row(automaticBox));config.add(marksCard);config.add(javax.swing.Box.createVerticalStrut(8))
        def diagnosticCard=new CalciumCard('','Optional diagnostics',true)
        diagnosticCard.body.add(row(sensitivityButton,diagnosticButton));diagnosticCard.body.add(row(clearDiagnosticsButton,templateButton));diagnosticCard.body.add(diagnosticNote);config.add(diagnosticCard)
        [confirmed,kcl,automaticBox].each{CalciumTheme.plain(it)};confirmed.text='<html>I checked the timing and baselines<br>and reviewed cell selection</html>'
        [status,appliedNote].each{it.editable=false;it.lineWrap=true;it.wrapStyleWord=true;it.opaque=false;it.focusable=false;it.border=null;it.font=CalciumTheme.font();it.foreground=CalciumTheme.INK;it.alignmentX=0f;it.maximumSize=new java.awt.Dimension(Integer.MAX_VALUE,(int)it.preferredSize.height)}
        progress.setStringPainted(true);progress.alignmentX=0f;progress.maximumSize=new java.awt.Dimension(Integer.MAX_VALUE,24)
        def side=new javax.swing.JPanel();side.layout=new javax.swing.BoxLayout(side,javax.swing.BoxLayout.Y_AXIS);side.opaque=false;side.preferredSize=new java.awt.Dimension(305,640)
        def recording=CalciumTheme.label(''+(data.recording?:''),CalciumTheme.MUTED);recording.toolTipText=''+(data.recording?:'')
        def source=CalciumTheme.label("<html><body style='width:245px'>Uses full-precision background-corrected fluorescence. Background is subtracted once; no intensity lift. QC warnings remain reviewable.</body></html>",CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f)
        [CalciumTheme.label('Analyze responses',CalciumTheme.INK,java.awt.Font.BOLD,5f),recording,javax.swing.Box.createVerticalStrut(8),timing,javax.swing.Box.createVerticalStrut(8),confirmed,row(runButton,cancelButton),progress,status,appliedNote,row(setup),row(logButton),source,javax.swing.Box.createVerticalGlue()].each{side.add(it)}
        root.add(side,java.awt.BorderLayout.WEST)
        def stats=new javax.swing.JPanel(new java.awt.GridBagLayout());stats.opaque=false
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->def c=new java.awt.GridBagConstraints();c.gridx=i;c.fill=java.awt.GridBagConstraints.BOTH;c.weightx=i==0?2.2d:1d;c.insets=new java.awt.Insets(0,0,0,i<3?8:0);tile.wide=i==0?220:100;stats.add(tile,c)}
        // the Traces tab: which cells, the picture at the cursor's moment, the traces, and the time slider under both
        [cells,group,view,heatmapOrder,stimuli,reasonBox].each{it.setFont(CalciumTheme.font())}
        // the Cell box is as wide as its longest name and never wider than for a name such as Cell001, so Window and Every
        // trace stay in the row at the window's usual size: a recording's names (Cell_ and 32 letters) are cut in the box, and
        // its tooltip gives the selected cell's whole name
        def across={String s->cells.getFontMetrics(cells.font).stringWidth(s)}
        String widest=(measurement.names as List<String>)?.max{across(it)}
        cells.prototypeDisplayValue=widest!=null&&across(widest)<=across('Cell000')?widest:'Cell000'
        previousCell.setToolTipText('Previous cell (left arrow)');nextCell.setToolTipText('Next cell (right arrow)')
        group.setToolTipText('List and draw only the cells with this response call');allTraces.setToolTipText('Draw every listed cell behind the mean, coloured by its response call')
        def white={javax.swing.JPanel p->p.setOpaque(true);p.setBackground(Color.white);p}
        def pickers=white(row(CalciumTheme.label('Show'),group,CalciumTheme.label('Cell'),cells,previousCell,nextCell,CalciumTheme.label('Window'),stimuli,allTraces))
        // the call line, then the cell's flags, then its mark controls, then a refusal, each directly under the one before (a row
        // with nothing to show is hidden); the room they leave of the room kept for them is under them (headRoom), so the
        // head has one height for every cell
        reviewBox.setLayout(new javax.swing.BoxLayout(reviewBox,javax.swing.BoxLayout.Y_AXIS));white(reviewBox);reviewBox.setBorder(javax.swing.BorderFactory.createEmptyBorder(0,8,6,8))
        [markButton,removeMark,reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{markRow.add(it)}
        [cellStatus,flagRow,markRow,markProblem,headRoom].each{reviewBox.add(it)}
        def head=new javax.swing.JPanel(new java.awt.BorderLayout());head.setBackground(Color.white);head.add(pickers,java.awt.BorderLayout.NORTH);head.add(reviewBox,java.awt.BorderLayout.SOUTH)
        def tracePanel=new javax.swing.JPanel(new java.awt.BorderLayout());tracePanel.setBackground(Color.white);tracePanel.add(head,java.awt.BorderLayout.NORTH)
        // the plot, and under it its caption and rule line (on screen only; Save plots saves the plot)
        def plotPane=new javax.swing.JPanel(new java.awt.BorderLayout());plotPane.setBackground(Color.white);plotPane.add(plot);plotPane.add(plot.notes,java.awt.BorderLayout.SOUTH)
        if(cellPicture.image!=null) {
            view.setToolTipText('What the picture shows. "Change from rest" is each frame minus the average of the baseline frames, so a responding cell lights up.')
            outlinesBox.setToolTipText('Untick to see the cells themselves: only the selected cell keeps its outline')
            def pictureTools=row(view,outlinesBox,zoomButton,wholeButton);pictureTools.setOpaque(true);pictureTools.setBackground(cellPicture.getBackground());outlinesBox.setForeground(new Color(0xC5CED8))
            def picturePane=new javax.swing.JPanel(new java.awt.BorderLayout());picturePane.add(cellPicture);picturePane.add(pictureTools,java.awt.BorderLayout.NORTH)
            def split=new javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT,picturePane,plotPane);split.setResizeWeight(0.5d);split.setBorder(null);split.setContinuousLayout(true);tracePanel.add(split)
            // Give the picture as much width as its shape can use at this height (at most 62 %), until the divider is dragged by hand.
            boolean dragged=false
            if(split.getUI() instanceof javax.swing.plaf.basic.BasicSplitPaneUI)(split.getUI() as javax.swing.plaf.basic.BasicSplitPaneUI).getDivider().addMouseListener(new java.awt.event.MouseAdapter(){void mousePressed(java.awt.event.MouseEvent e){dragged=true}})
            split.addComponentListener(new java.awt.event.ComponentAdapter(){void componentResized(java.awt.event.ComponentEvent e){
                if(dragged||split.getWidth()<=0||split.getHeight()<=0)return
                int ideal=(int)((split.getHeight()-pictureTools.getHeight()-30)*snapshot.getWidth()/(double)snapshot.getHeight())+8
                split.setDividerLocation((int)Math.max(240,Math.min(split.getWidth()*0.62d,ideal)))
            }})
        } else tracePanel.add(plotPane)
        slider.setMaximum(Math.max(2,n));slider.setOpaque(false);slider.setToolTipText('Move through the recording; the picture and the line on the plot follow')
        def timeRow=new javax.swing.JPanel(new java.awt.BorderLayout(8,0));timeRow.setBackground(Color.white);timeRow.setBorder(javax.swing.BorderFactory.createEmptyBorder(4,8,6,8))
        timeLabel.setPreferredSize(new java.awt.Dimension(170,(int)timeLabel.getPreferredSize().height))
        timeRow.add(playButton,java.awt.BorderLayout.WEST);timeRow.add(slider);timeRow.add(timeLabel,java.awt.BorderLayout.EAST)
        tracePanel.add(timeRow,java.awt.BorderLayout.SOUTH)
        // the Around the cell tab: the call line and the flags of the cell selected on Traces, the room kept for them under
        // them, and the plot of what the passing-object check used around the cell
        aroundBox.setLayout(new javax.swing.BoxLayout(aroundBox,javax.swing.BoxLayout.Y_AXIS));white(aroundBox);aroundBox.setBorder(javax.swing.BorderFactory.createEmptyBorder(8,8,6,8))
        [aroundStatus,aroundFlags,aroundRoom].each{aroundBox.add(it)}
        def aroundPanel=new javax.swing.JPanel(new java.awt.BorderLayout());aroundPanel.setBackground(Color.white);aroundPanel.add(aroundBox,java.awt.BorderLayout.NORTH);aroundPanel.add(aroundPlot)
        aroundPlot.names=((measurement.names?:[]) as List).collect{''+it}
        tabs.setFont(CalciumTheme.font())
        def heatmapPanel=new javax.swing.JPanel(new java.awt.BorderLayout());heatmapPanel.setBackground(Color.white)
        heatmapPanel.add(white(row(CalciumTheme.label('Order'),heatmapOrder,CalciumTheme.label('strip on the left: response call · click a row to see that cell',CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f))),java.awt.BorderLayout.NORTH)
        def heatmapScroll=new javax.swing.JScrollPane(heatmap);heatmapScroll.setBorder(null);heatmapScroll.getVerticalScrollBar().setUnitIncrement(CalciumAnalysisHeatmap.ROW*3);heatmapPanel.add(heatmapScroll)
        tabs.addTab('Protocol',new javax.swing.JScrollPane(config));tabs.addTab('Traces',tracePanel);tabs.addTab('Around the cell',aroundPanel);tabs.addTab('Heatmap',heatmapPanel)
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
        automaticBox.addActionListener({changed()} as java.awt.event.ActionListener)
        kcl.addActionListener({if(!loading){if(kcl.selected&&windowModel.rowCount==1){Map primary=windowRows()[0];double baseEnd=(primary.End_s as double)+(primary.Baseline_end_s as double)-(primary.Baseline_start_s as double)+1d;Map control=CalciumResponseProtocol.window('Control','KCl',primary.End_s as double,baseEnd,baseEnd+10d);if(preset.selectedIndex==2)control.Duration_mode='samples';addWindow(control)};if(!kcl.selected&&windowModel.rowCount>1)windowModel.removeRow(1);changed()}} as java.awt.event.ActionListener)
        fields.values().each{field->field.document.addDocumentListener([insertUpdate:{e->changed()},removeUpdate:{e->changed()},changedUpdate:{e->changed()}] as javax.swing.event.DocumentListener)}
        windowModel.addTableModelListener({e->if(!loading&&!busy){if(preset.selectedIndex==2&&e.column==CalciumResponseProtocol.COLUMNS.indexOf('Start_s')&&e.firstRow>=0){loading=true;def start=windowModel.getValueAt(e.firstRow,e.column);windowModel.setValueAt(start,e.firstRow,3);windowModel.setValueAt(start,e.firstRow,4);loading=false};if(preset.selectedIndex==0){loading=true;preset.selectedIndex=1;validation.selectedItem='user_defined';loading=false};windowTable.repaint();changed()}} as javax.swing.event.TableModelListener)
        stimuli.addActionListener({if(!listing){kept=null;selectWindow()}} as java.awt.event.ActionListener)
        cells.addActionListener({if(!listing){letGo();selectCell(cells.selectedItem as String,true)}} as java.awt.event.ActionListener)
        def step={int by->int count=cells.getItemCount();if(count>0)cells.setSelectedIndex((int)(((cells.getSelectedIndex()+by)%count+count)%count))}
        previousCell.addActionListener({step(-1)} as java.awt.event.ActionListener);nextCell.addActionListener({step(1)} as java.awt.event.ActionListener)
        group.addActionListener({if(!listing){kept=null;listCells()}} as java.awt.event.ActionListener)
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
        aroundPlot.seek={double t->plot.seek.call(t)}
        slider.addChangeListener({e->if(slider.getValue()!=frame){stopPlaying();showFrame(slider.getValue())}} as javax.swing.event.ChangeListener)
        player=new javax.swing.Timer(50,{e->if(frame>=slider.getMaximum())stopPlaying() else showFrame(frame+1)} as java.awt.event.ActionListener)
        playButton.toolTipText='Display playback at 20 frames/s; cursor positions use the recording timestamps.'
        playButton.enabled=movie!=null&&view.selectedIndex>0
        playButton.addActionListener({if(player.isRunning())stopPlaying() else {if(frame>=slider.getMaximum())showFrame(1);player.start();playButton.setText('❚❚ Pause')}} as java.awt.event.ActionListener)
        view.addActionListener({if(view.selectedIndex==0)stopPlaying();playButton.enabled=movie!=null&&view.selectedIndex>0;showFrame(frame)} as java.awt.event.ActionListener)
        outlinesBox.addActionListener({cellPicture.outlines=outlinesBox.selected;cellPicture.repaint()} as java.awt.event.ActionListener)
        zoomButton.addActionListener({cellPicture.zoomTo(cellPicture.selected)} as java.awt.event.ActionListener)
        wholeButton.addActionListener({cellPicture.whole()} as java.awt.event.ActionListener)
        // left and right arrows step through the listed cells, unless a text box or list is using them, or a mark form is open
        def keys=getRootPane()
        keys.registerKeyboardAction({e->if(previousCell.enabled&&formCell==null)step(-1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_LEFT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        keys.registerKeyboardAction({e->if(nextCell.enabled&&formCell==null)step(1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_RIGHT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        // the mark controls: one button opens the form in its place; in it Tab moves on, Enter saves (or presses the button that
        // has the focus) and Esc closes it, keeping the cell selected; the focus then goes back to the button. Enter acts when
        // the key is let go, so a key held down (which repeats only its press) opens the form or saves once, as on the page
        markButton.addActionListener({openForm(cells.selectedItem as String)} as java.awt.event.ActionListener)
        removeMark.addActionListener({removeReviewMark(cells.selectedItem as String)} as java.awt.event.ActionListener)
        saveMark.addActionListener({saveReviewMark()} as java.awt.event.ActionListener);cancelMark.addActionListener({closeForm()} as java.awt.event.ActionListener)
        reasonBox.addActionListener({saveMark.enabled=reasonBox.selectedIndex>=0;noteField.hint=reasonBox.selectedItem=='Other'?'':'Optional'} as java.awt.event.ActionListener)
        // (the press only readies it: one taken by the open Reason list does not save when the key is let go, nor do its
        // repeats while the key is held, which come with the list closed; a key listener hears the press before the list
        // takes it, and the key let go before the mark controls do)
        reasonBox.addKeyListener(new java.awt.event.KeyAdapter(){
            void keyPressed(java.awt.event.KeyEvent e){if(e.keyCode==java.awt.event.KeyEvent.VK_ENTER&&listOpen.call())enterTaken=true}
            void keyReleased(java.awt.event.KeyEvent e){if(e.keyCode==java.awt.event.KeyEvent.VK_ENTER)enterTaken=false}})
        markRow.registerKeyboardAction({e->if(!enterTaken)enterPressed=true} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ENTER,0),javax.swing.JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT)
        markRow.registerKeyboardAction({e->enterTaken=false;if(!enterPressed)return;enterPressed=false;def owner=focusOwner()
            if(owner instanceof javax.swing.AbstractButton&&owner.parent==markRow)(owner as javax.swing.AbstractButton).doClick() else if(formCell!=null&&saveMark.enabled)saveReviewMark()} as java.awt.event.ActionListener,
            javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ENTER,0,true),javax.swing.JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT)
        markRow.registerKeyboardAction({e->if(formCell!=null)closeForm()} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ESCAPE,0),javax.swing.JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT)
        flagRow.visible=false;aroundFlags.visible=false;markRow.visible=false;closeForm()
        clock=new javax.swing.Timer(500,{e->if(busy)status.text=stage+' · '+(long)((System.currentTimeMillis()-started)/1000)+' s elapsed'} as java.awt.event.ActionListener)
        addWindowListener(new java.awt.event.WindowAdapter(){void windowClosed(java.awt.event.WindowEvent e){closed=true;job?.cancel();clock.stop();player?.stop();inputWatch?.stop();(measurement.picture as ImagePlus)?.close()}})
        inputWatch=new javax.swing.Timer(1500,{e->if(!busy&&!closed&&results!=null&&(timestampsFile||sensitivityFile||diagnosticFile))changed()} as java.awt.event.ActionListener);inputWatch.start()
        loading=false;[allTraces,outlinesBox].each{CalciumTheme.plain(it)};showEstimator();updatePreviewTiming();if(movie!=null)view.selectedIndex=2;showFrame(1)
    }
    void showEstimator(){noiseNote.text='Noise: '+CalciumResponseProtocol.ESTIMATOR_NAME[preset.selectedIndex==2?CalciumResponseProtocol.LEGACY_ESTIMATOR:estimator]}
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
        plot.moment=t;plot.repaint();aroundPlot.moment=t;aroundPlot.repaint()
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
    Map protocol(){stopEditing();CalciumResponseProtocol.payload(windowRows(),preset.selectedIndex==2?1:2,fields.protocol_name.text.trim(),validation.selectedItem as String,fieldOverrides,estimator)}
    Map fileSignature(File file,String loaded){if(file==null)return null;String current;try{current=file.isFile()?CalciumInstaller.digest(file.bytes):'missing'}catch(Throwable ignored){current='unreadable'};[path:file.absolutePath,loaded_sha256:loaded,current_sha256:current]}
    /** The settings as they are now, to tell whether they changed since the analysis shown. The automatic marks are in it only
     *  when they are on, so it is as before while they are off. */
    String fingerprint(){Map now=[fields:fields.collectEntries{k,v->[k,v.text]},windows:windowRows(),preset:preset.selectedIndex,validation:validation.selectedItem,estimator:estimator,overrides:fieldOverrides,timestamps:fileSignature(timestampsFile,timestampsHash),sensitivity:fileSignature(sensitivityFile,sensitivityHash),diagnostic:fileSignature(diagnosticFile,diagnosticHash)]
        if(automaticBox.selected)now[CalciumAnalysisSettings.AUTOMATIC_MARKS]=true
        groovy.json.JsonOutput.toJson(now)}
    void changed(){if(loading||busy)return;if(results==null){updatePreviewTiming();return};if(results!=null){boolean dirty=fingerprint()!=appliedFingerprint;appliedNote.text=dirty?'Settings changed — rerun to apply. Display and plot export still show the completed run.':'Displayed results match the applied settings.';appliedNote.foreground=dirty?new Color(155,75,0):new Color(20,95,45)}}
    void applyPreset(){if(loading)return;loading=true
        boolean legacy=preset.selectedIndex==2
        if(preset.selectedIndex!=1){estimator=CalciumResponseProtocol.newEstimator();fields.protocol_name.text=legacy?'Legacy consecutive-sample protocol':'Slow evoked response — provisional';validation.selectedItem=legacy?'user_defined':'provisional';for(int r=0;r<windowModel.rowCount;r++){['Floor_local_fraction':0.10d,'SD_multiplier':3d,'Duration_mode':'samples','Min_duration_s':1d,'Consecutive_samples':2].each{k,v->windowModel.setValueAt(v,r,CalciumResponseProtocol.COLUMNS.indexOf(k))};if(legacy){['Min_baseline_samples':3,'Noise_floor_local_fraction':0d,'Expected_interval_s':null].each{k,v->windowModel.setValueAt(v,r,CalciumResponseProtocol.COLUMNS.indexOf(k))};windowModel.setValueAt(windowModel.getValueAt(r,5),r,3);windowModel.setValueAt(windowModel.getValueAt(r,5),r,4);if(r==1)windowModel.setValueAt('KCl',r,1)}}}
        loading=false;windowTable.repaint();showEstimator();updateDiagnosticAvailability();changed()
    }
    void updateDiagnosticAvailability(){boolean available=preset.selectedIndex!=2&&!busy;sensitivityButton.enabled=available;diagnosticButton.enabled=available;templateButton.enabled=available;if(preset.selectedIndex==2){sensitivityFile=null;diagnosticFile=null;sensitivityHash='';diagnosticHash='';diagnosticNote.text='Optional sensitivity and negative-control diagnostics require version 2. Previously selected diagnostic files are detached.'}}
    void chooseFile(String title,boolean save,Closure action){try{def chooser=new javax.swing.JFileChooser(measurement.folder as File);chooser.dialogTitle=title;int answer=save?chooser.showSaveDialog(this):chooser.showOpenDialog(this);if(answer==javax.swing.JFileChooser.APPROVE_OPTION)action.call(chooser.selectedFile)}catch(Throwable t){status.text=t.message?:t.toString()}}
    void saveProtocol(File file){file.setText(CalciumResponseProtocol.json(protocol())+'\n','UTF-8');status.text='Saved draft protocol: '+file.name}
    void loadProtocol(File file){Map p=CalciumResponseProtocol.read(file);loading=true;try{preset.selectedIndex=p.schema_version==1?2:1;estimator=p.schema_version==1?CalciumResponseProtocol.newEstimator():p.sd_estimator as String;fields.protocol_name.text=p.name;validation.selectedItem=p.schema_version==1?'user_defined':p.validation_status;fieldOverrides=p.field_overrides?:[];windowModel.rowCount=0;p.windows.each{addWindow(it as Map)};kcl.selected=p.windows.size()==2;confirmed.selected=false}finally{loading=false};showEstimator();updateOverrideNote();updateDiagnosticAvailability();changed();status.text='Loaded '+file.name+(fieldOverrides?' · inspect the retained field timing overrides.':'')}
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
        closeForm();showReview(cells.selectedItem as String)   // no mark is set or removed while an analysis runs
        fields.values().each{it.enabled=false};windowTable.enabled=false;[confirmed,kcl,automaticBox,preset,validation,loadProtocolButton,saveProtocolButton,loadTimesButton,clearTimesButton,sensitivityButton,diagnosticButton,clearDiagnosticsButton,templateButton,overrideButton,clearOverridesButton].each{it.enabled=false};savePlots.enabled=false
        job=new CalciumLocalJob(logFile:new File(runtime.home,installing?'setup.log':'analysis.log'))
        job.progress={int p,String message->javax.swing.SwingUtilities.invokeLater{if(!closed){progress.value=p;progress.setString(p+'% · '+message);stage=message;logButton.enabled=job.logFile.isFile()}}}
        clock.start()
        Thread.startDaemon('SimpleCa analysis'){
            File output=null;String error=null
            try{output=action.call(job)}catch(java.util.concurrent.CancellationException ignored){error='Cancelled. Completed analyses are kept; partial output is not used.'}catch(Throwable t){error=t.message?:t.toString()}
            final File answer=output;final String failure=error
            javax.swing.SwingUtilities.invokeLater{
                busy=false;clock.stop();if(closed)return
                showReview(cells.selectedItem as String)
                progress.indeterminate=false;cancelButton.enabled=false;setup.enabled=true;runButton.enabled=true;fields.values().each{it.enabled=true};fields.sampling_hz.enabled=timestampsFile==null;windowTable.enabled=true
                [confirmed,kcl,automaticBox,preset,validation,loadProtocolButton,saveProtocolButton,loadTimesButton,clearTimesButton,sensitivityButton,diagnosticButton,clearDiagnosticsButton,templateButton,overrideButton,clearOverridesButton].each{it.enabled=true};logButton.enabled=job.logFile.isFile();savePlots.enabled=results!=null;updateDiagnosticAvailability()
                if(failure){status.text=failure;status.foreground=CalciumTheme.INK;progress.string='Not completed';changed();return}
                if(installing){progress.value=100;status.text='R and protocol support are ready. Check settings, then analyze.';return}
                try{loadResults(answer);appliedFingerprint=submittedFingerprint;changed();String empty=CalciumAnalysisSettings.noEvaluableCells(answer)
                    if(!empty)CalciumAnalysisSettings.percentResponding(answer,measurement)   // saves the percent table
                    finish=[folder:answer.name,empty:empty,saved:'Complete. Protocol, inputs and tables saved in '+answer.name+'.'];status.text=finished()
                    status.foreground=empty?CalciumTheme.WARN:CalciumTheme.INK;progress.value=100;progress.string='Complete';tabs.selectedIndex=1}catch(Throwable t){status.text='Could not display results: '+t.message;progress.string='Display error; see saved outputs'}
            }
        }
    }
    void beginSetup(){launch(true){task->runtime.setup(task)}}
    /** Runs the analysis, then what follows every analysis (CalciumReview.afterAnalysis): the checks of the calls, the carry-over
     *  of the review marks and, with the option on, the automatic marks. The option is read here, as the run starts. */
    void beginAnalysis(){try{Map s=settings();submittedFingerprint=fingerprint();boolean on=automaticBox.selected
        launch(false){task->File made=performAnalysis.call(runtime,measurement.csv,s,measurement.folder,task) as File
            Map after=CalciumReview.afterAnalysis(made,withRecording(),measurement,on,task.logFile);checkedRun=[folder:made.name,found:after.checked];carried=after.carried;markedRun=after.marked;made}}catch(Throwable t){status.text=t.message?:t.toString()}}
    /** The finish message of the last analysis run here, with the figures of the tally as they are now, so that it follows the
     *  review marks as the Step 4 page's does; then the cells marked automatically, or why none could be; then what was carried
     *  over from an earlier analysis's review; both are said until the next run, marks or not (as on the page); then what the
     *  checks of the calls could not do (checks.note). null when the analysis shown was not run here. Writes nothing, and reads
     *  no traces: it takes the window's. */
    String finished() {
        Map parts=finish;if(parts==null||results?.name!=parts.folder)return null
        String carry=CalciumReview.note(carried?.folder==results.name?carried:null),checked=checks?.note?' '+checks.note:''
        String sentence='';if(!parts.empty)try{sentence=CalciumAnalysisSettings.percentSentence(CalciumAnalysisSettings.tally(results,measurement,traces))}catch(Throwable ignored){}   // as percentResponding: '' when the tally cannot be made
        (parts.empty?parts.empty as String:(parts.saved as String)+sentence)+CalciumReview.automaticNote(markedRun?.folder==results.name?markedRun:null)+carry+checked
    }
    List<Map> table(javax.swing.JTable target,File file){
        if(!file.isFile()){target.model=new javax.swing.table.DefaultTableModel();return []}
        def rows=CalciumCsv.read(file);if(!rows)return []
        showTable(target,new javax.swing.table.DefaultTableModel(rows.drop(1).collect{it as Object[]} as Object[][],rows[0] as Object[]){boolean isCellEditable(int r,int c){false}})
        rows.drop(1).collect{row->rows[0].withIndex().collectEntries{name,i->[name,row[i]]}}
    }
    void showTable(javax.swing.JTable target,javax.swing.table.TableModel model){target.model=model;for(int i=0;i<target.columnCount;i++)target.columnModel.getColumn(i).preferredWidth=160
        if(!(target.getDefaultRenderer(Object) instanceof CalciumCutCell))target.setDefaultRenderer(Object,new CalciumCutCell())}
    /** The Cell response calls and Responder counts tabs, as saved; while review marks apply, with columns on their right as on
     *  the Step 4 page: the reviewed call and the review reason (with its note) on the primary row of each cell whose call a
     *  mark changed, and % responding after review with the calls changed on the stimulus row. On screen only: no table is
     *  written, and nothing is read again. */
    void reviewTables() {
        // a saved table with two columns more, filled for each row from the row by its column names
        def extended={javax.swing.table.TableModel saved,List<String> more,Closure fill->
            List names=(0..<saved.columnCount).collect{saved.getColumnName(it)}
            Object[][] data=(0..<saved.rowCount).collect{int r->List row=(0..<saved.columnCount).collect{saved.getValueAt(r,it)}
                (row+(fill.call([names,row].transpose().collectEntries{it}) as List)) as Object[]} as Object[][]
            new javax.swing.table.DefaultTableModel(data,(names+more) as Object[]){boolean isCellEditable(int r,int c){false}}}
        def stimulusRow={Map row->!row.containsKey('Role')||row.Role=='Stimulus'}
        Map marks=(reviewMarks?:[:]) as Map
        if(callsSaved!=null)showTable(calls,marks.isEmpty()?callsSaved:extended(callsSaved,['Reviewed call','Review reason'],{Map row->
            String cell=''+(row.Cell?:''),status=''+(row.Status?:'');Map mark=marks[cell] as Map
            stimulusRow(row)&&mark!=null&&CalciumReview.reviewed(marks,cell,status)!=status ? [mark.reviewed,(mark.reason as String)+(mark.note?': '+mark.note:'')+(mark.by=='automatic'?', automatic':'')] : ['','']}))
        Map all=reviewTally?.reviewed==true?(reviewTally.rows as List).find{it.key=='all'} as Map:null
        Map after=marks.isEmpty()||all==null?null:all.reviewed as Map
        if(countsSaved!=null)showTable(counts,after==null?countsSaved:extended(countsSaved,['After review','Calls changed'],{Map row->
            List r=after.counted as List
            stimulusRow(row) ? [(r[1] as int)>0?CalciumAnalysisSettings.percentOf(r)+'% ('+r[0]+' of '+r[1]+')':'–',''+after.changed] : ['','']}))
    }
    /** The headline numbers and the per-cell calls, read from the tables the analysis just saved. Display only. */
    void summarize(File dir) {
        reviewTally=null
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
        // with review marks the tile gives % responding after review, with the count before review beside it (CalciumReview);
        // the tally takes the traces this window holds, so normalized_traces.csv is not read again for every mark
        Map review=CalciumReview.read(dir,measurement);reviewMarks=review.marks as Map;reviewProblem=(review.problem?:'') as String
        reviewTally=CalciumAnalysisSettings.tally(dir,measurement,traces,review)
        List after=CalciumReview.headline(reviewTally,cell('Stimulus'))
        if(after!=null){responding.show(after[0] as String,after[1] as String);responding.setToolTipText('Primary stimulus after the review marks in '+review.file+'. Selecting another trace window does not change this summary.')}
        else if(review.problem)responding.setToolTipText(review.problem as String)
        int callRole=column(callRows,'Role'),callCell=column(callRows,'Cell'),callStatus=column(callRows,'Status');Map<String,String> found=[:]
        if(callCell>=0&&callStatus>=0)callRows.drop(1).each{if(callRole<0||it[callRole]=='Stimulus')found[''+it[callCell]]=CalciumReview.reviewed(reviewMarks,''+it[callCell],''+it[callStatus])}
        cellPicture.calls=found;cellPicture.repaint()
    }
    /** The call drawn and listed for a cell: for the primary stimulus after review, for a control window the automatic call. */
    String callOf(String name){Map row=callRows[name];row==null?null:row.Role=='Control'?row.Status:CalciumReview.reviewed(reviewMarks,name,row.Status as String)}
    /** Whether a cell is in an entry of the Show list: a group of calls (CalciumCalls.GROUPS), or MARKED. */
    boolean listedIn(String name,int which){which==MARKED?reviewMarks?.containsKey(name)&&callRows[name]!=null&&callRows[name].Role!='Control':CalciumCalls.inGroup(callOf(name),which)}
    /** All linked views use one saved stimulus/control window, never draft settings. */
    void selectWindow(){
        callRows=responseRows.findAll{it.Stimulus==stimuli.selectedItem}.collectEntries{[(it.Cell):it]}
        plot.windowCalls=callRows.values() as List   // the key names the bands of the window for every cell
        cellPicture.calls=callRows.collectEntries{k,v->[k,callOf(k as String)]}
        Map first=callRows.values().find{true}?:[:];plot.stimulusStart=CalciumCsv.number(first.Start_s);plot.stimulusEnd=CalciumCsv.number(first.Effective_end_s?:first.End_s)
        boolean after=reviewMarks&&callRows.values().any{it.Role!='Control'}   // a control window keeps its automatic call
        // Show is as wide as its longest entry can be, in both windows and with or without marks, so the Cell box and the
        // arrows beside it do not move when marks come and go (as on the Step 4 page)
        group.prototypeDisplayValue='Did not respond after review ('+(measurement.names as List).size()+')'
        listing=true;try{int keep=Math.max(0,group.selectedIndex);group.removeAllItems();CalciumCalls.GROUPS.eachWithIndex{String name,int i->group.addItem(name+(after&&(i==1||i==2)?' after review':'')+' ('+(measurement.names as List<String>).count{callRows.isEmpty()||CalciumCalls.inGroup(callOf(it),i)}+')')}
            if(after)group.addItem('Marked ('+(measurement.names as List<String>).count{listedIn(it,MARKED)}+')')
            group.selectedIndex=callRows.isEmpty()||keep>=group.itemCount?0:keep}finally{listing=false}
        listCells();orderHeatmap();plot.changed()
    }
    /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. A cell kept after its
     *  mark (kept) is listed too, in its place, but not drawn with the group. */
    void listCells() {
        int which=Math.max(0,group.getSelectedIndex());String keep=cells.getSelectedItem() as String
        List<String> listed=(measurement.names as List<String>).findAll{callRows.isEmpty()||listedIn(it,which)}
        if(kept!=null&&listed.contains(kept))kept=null   // in the group again: listed as the others are
        List<String> shown=kept==null?listed:(measurement.names as List<String>).findAll{it==kept||listed.contains(it)}
        listing=true
        try {cells.removeAllItems();shown.each{cells.addItem(it)};if(shown.contains(keep))cells.setSelectedItem(keep)} finally {listing=false}
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->tile.active=!callRows.isEmpty()&&((i==0&&which==1)||(i==1&&which==0)||(i>=2&&which==3));tile.repaint()}
        cellPicture.group=which;cellPicture.only=which==MARKED?new LinkedHashSet<String>(listed):null
        List<String> drawn=listed.findAll{traces.containsKey(it)}
        plot.all=drawn.collect{traces[it]} as double[][];plot.allCalls=drawn.collect{callOf(it)};plot.allLabel=which==0?'all cells':which==MARKED?'marked':CalciumCalls.GROUPS[which].toLowerCase()
        // the calls a mark can bring into the group (on the primary stimulus): the plot's key keeps room for a count of each
        boolean marks=callRows.values().any{it.Role!='Control'&&it.Status in CalciumReview.CALLS}
        plot.markCalls=!marks ? [] : which==0||which==MARKED ? CalciumReview.CALLS : which==1||which==2 ? [CalciumReview.CALLS[which-1]] : []
        previousCell.enabled=shown.size()>1;nextCell.enabled=shown.size()>1
        if(cells.getItemCount()>0)selectCell(cells.getSelectedItem() as String,true)
        else {plot.cell='';plot.selected=new double[0];plot.responseCall=[:];plot.caption='';cells.setToolTipText(null);showCall('No cells in this group.');showReview(null);showAround(null);cellPicture.selected=-1;plot.changed();cellPicture.repaint()}
    }
    void showGroup(int which){if(which<group.getItemCount()&&group.getSelectedIndex()!=which)group.setSelectedIndex(which)}
    /** Another cell is selected: the cell kept after its mark leaves the Cell list. */
    void letGo(){
        if(kept==null||cells.selectedItem==kept)return
        String gone=kept;kept=null;listing=true;try{cells.removeItem(gone)}finally{listing=false}
        previousCell.enabled=cells.itemCount>1;nextCell.enabled=cells.itemCount>1
    }
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
        Map review=CalciumReview.read(dir,measurement);reviewMarks=review.marks as Map;reviewProblem=(review.problem?:'') as String   // read again by summarize, which also shows the tile after review
        closeForm();automatic=CalciumReview.scored(dir)
        // the checks of this analysis: what its run here gave (as the Step 4 page takes them), else the table saved beside it
        Map run=checkedRun;checks=CalciumCallChecks.answer(run?.folder==dir.name?run.found as Map:CalciumCallChecks.saved(dir,measurement))
        aroundAsked++;aroundContext=null;aroundFolder=null   // what is shown around the cells is read again for this analysis (readAround)
        results=dir;table(metrics,new File(dir,'metrics.csv'));table(qc,new File(dir,'cell_quality.csv'));table(registry,new File(dir,'roi_registry.csv'));responseRows=table(calls,new File(dir,'responder_calls.csv'));table(counts,new File(dir,'responder_counts.csv'));table(windowMetrics,new File(dir,'window_metrics.csv'));table(replicates,new File(dir,'replicate_counts.csv'));table(sensitivity,new File(dir,'sensitivity_counts.csv'));table(diagnostics,new File(dir,'negative_control_counts.csv'))
        withCodes()
        callsSaved=calls.model;countsSaved=counts.model   // as saved, with the codes of the checks: the columns of the review marks are added to copies (reviewTables)
        responding.caption.lines=2   // the caption after review takes two lines at the window's usual width: the tiles keep that height
        def normalized=CalciumCsv.read(new File(dir,'normalized_traces.csv')),avg=CalciumCsv.read(new File(dir,'mean_trace.csv'))
        plot.time=normalized.drop(1).collect{CalciumCsv.number(it[0])} as double[]
        traces=[:];normalized[0].drop(1).eachWithIndex{name,i->traces[name]=normalized.drop(1).collect{CalciumCsv.number(it[i+1])} as double[]}
        plot.mean=avg.drop(1).collect{CalciumCsv.number(it[1])} as double[];plot.sem=avg.drop(1).collect{CalciumCsv.number(it[2])} as double[]
        appliedManifest=CalciumCsv.read(new File(dir,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
        int first=Math.max(0,(CalciumCsv.number(appliedManifest.baseline_start) as int)-1),last=Math.max(0,(CalciumCsv.number(appliedManifest.baseline_end) as int)-1)
        plot.baselineStart=plot.time[Math.min(first,plot.time.length-1)];plot.baselineEnd=plot.time[Math.min(last,plot.time.length-1)]
        listing=true;try{stimuli.removeAllItems();responseRows.collect{it.Stimulus}.unique().each{stimuli.addItem(it)}}finally{listing=false}
        try{summarize(dir)}catch(Throwable ignored){}
        reviewTables()
        restFirst=first+1;restLast=last+1;rest=null
        kept=null;selectWindow();folderButton.enabled=true;savePlots.enabled=true
        int peak=1;double top=Double.NEGATIVE_INFINITY;plot.mean.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top){top=v;peak=i+1}}
        showFrame(peak)
        // an analysis whose review file holds an automatic mark was run with the option on: the box says so, before the settings
        // it was run with are taken from the window
        if(((reviewMarks?:[:]) as Map).values().any{(it as Map).by=='automatic'})automaticBox.selected=true
        if(!appliedFingerprint)appliedFingerprint=fingerprint();changed()
        readAround(dir)
    }
    /** The check codes of the primary stimulus after the engine's own codes in each tab that shows QC_flags, joined as that
     *  table joins them: as the checks say ("; ") in the calls and stimulus-window tables, ";" in Cell quality and the
     *  whole-recording metrics; without an engine code, the check codes alone. A row of the control window gets none. On screen
     *  only: no table is written. */
    void withCodes() {
        Map byCell=((checks?.cells?:[:]) as Map);if(byCell.isEmpty())return
        String join=(checks.join?:CalciumCallChecks.JOIN) as String
        [[calls,join],[windowMetrics,join],[qc,QUALITY_JOIN],[metrics,QUALITY_JOIN]].each{List pair->
            def model=(pair[0] as javax.swing.JTable).model;String by=pair[1] as String
            List heads=(0..<model.columnCount).collect{model.getColumnName(it)}
            int flagsAt=heads.indexOf('QC_flags'),cellAt=heads.indexOf('Cell'),roleAt=heads.indexOf('Role')
            if(flagsAt<0||cellAt<0)return
            for(int r=0;r<model.rowCount;r++) {
                if(roleAt>=0&&''+model.getValueAt(r,roleAt)!='Stimulus')continue
                List codes=(((byCell[''+model.getValueAt(r,cellAt)] as Map)?.codes?:[]) as List)
                if(codes.isEmpty())continue
                String engine=''+(model.getValueAt(r,flagsAt)?:'')
                model.setValueAt((engine?engine+by:'')+codes.join(by),r,flagsAt)
            }
        }
    }
    /** Reads what is shown around the cells of an analysis (CalciumCallChecks.aroundContext) once, off the event thread, and
     *  shows it for the cell selected when it is read, unless another analysis was loaded meanwhile. */
    void readAround(File dir) {
        int asked=++aroundAsked;Map given=measurement
        aroundLater.call{->
            Map context=CalciumCallChecks.aroundContext(dir,given)
            def use={->if(asked==aroundAsked&&!closed){aroundContext=context;aroundFolder=dir.name;showAround(cells.selectedItem as String)}}
            if(javax.swing.SwingUtilities.isEventDispatchThread())use() else javax.swing.SwingUtilities.invokeLater(use as Runnable)
        }
    }
    /** The Around the cell tab for a cell of the window shown: what the passing-object check used around it (aroundOf), or why
     *  nothing is shown: before an analysis, on the control window, while what is shown is read, and as aroundCell says. */
    void showAround(String name) {
        aroundPlot.cell=name?:''
        Map call=(name==null?null:callRows[name])?:[:]
        double[] span=CalciumAnalysisPlot.windowOf(call);aroundPlot.stimulusStart=span[0];aroundPlot.stimulusEnd=span[1]
        boolean control=callRows.values().any{it.Role=='Control'}
        aroundPlot.data=results==null ? null : name==null ? [why:'No cells in this group.'] : control ? [why:CalciumCallChecks.AROUND_CONTROL]
            : aroundContext==null||aroundFolder!=results.name ? [why:'Reading the values around the cells.'] : aroundOf.call(aroundContext,name) as Map
        aroundPlot.changed()
    }
    /** The call line of the selected cell, on Traces and on Around the cell alike; a tip of null keeps the tooltip. */
    void showCall(String text,String tip=null){[cellStatus,aroundStatus].each{it.setText(text);if(tip!=null)it.setToolTipText(tip)}}
    /** The room kept under the Around the cell tab's own controls: three lines of the call line, as on Traces, and a row of flags. */
    int aroundKept(){int line=aroundStatus.getFontMetrics(aroundStatus.font).height;3*line+line+2*CalciumRowLayout.PAD}
    void selectCell(String name,boolean highlight) {
        if(name==null)return
        cells.setToolTipText(name)   // the whole name: a long one is cut in the box
        plot.cell=name;plot.selected=traces[name]?:new double[0]
        Map call=callRows[name]?:[:];plot.responseCall=call
        plot.stimulusStart=CalciumCsv.number(call.Start_s);plot.stimulusEnd=CalciumCsv.number(call.Effective_end_s?:call.End_s)
        // the engine's reason explains the automatic call: it is left out when a mark changed the call shown
        String callShown=callOf(name)
        plot.caption=name+(traces.containsKey(name)?' · '+(callShown?:'')+(callShown==call.Status?' · '+(call.Reason?:''):''):' · unavailable: inspect Cell quality for the exclusion reason')
        plot.changed()
        int index=measurement.names.findIndexOf{it==name}
        if(highlight&&index>=0&&measurement.picture!=null){Roi roi=(Roi)measurement.rois[index].clone();roi.setStrokeColor(Color.white);roi.setStrokeWidth(4);(measurement.picture as ImagePlus).setRoi(roi)}
        cellPicture.selected=index;cellPicture.reveal(cellPicture.selected);cellPicture.repaint()
        zoomButton.enabled=cellPicture.selected>=0
        heatmap.selected=name;heatmap.repaint()
        // one line saying what the analysis decided for this cell, in the words and numbers of the saved table; with a mark the
        // call after review comes first, then the mark's reason and note, then the automatic call
        Map<String,String> row=callRows[name]
        if(row==null) showCall(callRows.isEmpty()?' ':'<html><b>'+name+'</b> · no response call was saved for this cell</html>')
        else {
            def number={String key,String format->double v=CalciumCsv.number(row[key]);Double.isFinite(v)?String.format(Locale.US,format,v):null}
            def esc={Object text->CalciumOutlineProgress.escape(text==null?'':''+text)}
            def dot={String shown->"<font color='"+CalciumTheme.hex(CalciumCalls.of(shown))+"'>●</font> "}
            Map mark=row.Role=='Control'?null:reviewMarks[name] as Map;String after=callOf(name);boolean marked=mark!=null&&after!=row.Status
            List parts=marked ? [dot(after)+'<b>'+name+': '+esc(after)+' after review</b>',esc((mark.reason as String).toLowerCase()+(mark.note?': '+mark.note:''))]+(mark.by=='automatic' ? ['marked automatically'] : [])+(mark.carriedFrom ? ['from an earlier analysis'] : [])+['automatic call '+esc(row.Status)]
                : [dot(row.Status)+'<b>'+name+': '+esc(row.Status?:'')+'</b>']
            parts<<esc((row.Reason?:'').toLowerCase())
            String rise=number('Peak_increase','%.2f'),at=number('Peak_s','%.0f'),need=number('Threshold_increase','%.2f')
            if(rise!=null)parts<<'largest rise '+rise+(at!=null?' at '+at+' s':'')+(need!=null?' (needed '+need+')':'')
            String tip='Rise and threshold are fractions of the fluorescence just before the stimulus, as in the Cell response calls table.'
            // when the mark was made, to the minute as the Step 4 page says it, and in which analysis when it was carried over
            String when=marked&&mark.time?CalciumReview.minute(mark.time as String):''
            showCall('<html>'+parts.findAll{it}.join(' · ')+'</html>',marked?'<html>'+esc('Marked'+(when?' '+when:'')+(mark.carriedFrom?' in '+mark.carriedFrom:''))+'<br>'+esc(tip)+'</html>':tip)
        }
        showReview(name)
        showAround(name)
    }
    /** Under the selected cell's call line: the flags of its checks, each in the checks' own words, and its mark controls. A flag
     *  is a link to the first frame it names (fillFlags). A control window keeps its automatic call: it has neither.
     *  The Around the cell tab has the same flags under its call line. */
    void showReview(String name) {
        Map row=name==null?null:callRows[name];boolean primary=row!=null&&row.Role!='Control'
        Map cell=primary?((checks?.cells?:[:]) as Map)[name] as Map:null
        String auto=primary?automatic[name] as String:null;boolean evaluable=auto in CalciumReview.CALLS
        [flagRow,aroundFlags].each{fillFlags(it,cell,evaluable)};aroundBox.revalidate();aroundBox.repaint()
        // one button: mark the cell the other way, or remove its mark; none while an analysis runs or while the review file
        // is not applied, which is said under the controls of every cell of the primary stimulus instead
        Map can=markable(name)
        if(formCell!=null&&(formCell!=name||!can.can||can.marked))closeForm()
        String blocked=reviewProblem?CalciumReview.blocked(reviewProblem):null
        if(problemCell!=name||(blocked&&problemSaid==blocked))say(null,null)
        if(blocked&&primary&&problemCell==null)say(blocked,name)
        markButton.text=auto=='Responded'?'Mark as not a responder':'Mark as a responder'
        markButton.visible=can.can&&!can.marked&&formCell==null;removeMark.visible=can.marked;enterPressed=false;enterTaken=false
        showMarkRow()
    }
    /** A row of the selected cell's flags, each in the checks' own words; a flag is a link to the first frame its text names
     *  (its frames), else to the first frame of the event. Then, for an evaluable cell that could not be checked in a recording
     *  whose surroundings were measured, why (else checks.note says it once). Shown while it holds something. */
    void fillFlags(javax.swing.JPanel into,Map cell,boolean evaluable) {
        def owner=this;into.removeAll();boolean around=into.is(aroundFlags)
        ((cell?.flags?:[]) as List).each{Map flag->
            if(into.componentCount>0)into.add(CalciumTheme.label('·',CalciumTheme.MUTED))
            List named=flag.frames instanceof List&&(flag.frames as List).size()==2&&(flag.frames as List)[0] instanceof Number?flag.frames as List:null
            Object at=named!=null?named[0]:cell.firstFrame
            def label=CalciumTheme.label(flag.text as String,at!=null?CalciumTheme.BLUE:CalciumTheme.MUTED)
            String detail=cell.detail?:''
            if(at!=null) {
                int first=at as int
                label.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR))
                label.setToolTipText('<html>Go to frame '+first+(detail?'<br>'+CalciumOutlineProgress.escape(detail):'')+'</html>')
                label.addMouseListener(new java.awt.event.MouseAdapter(){void mouseClicked(java.awt.event.MouseEvent e){if(javax.swing.SwingUtilities.isLeftMouseButton(e))owner.goToFrame(first,around)}})
            } else if(detail)label.setToolTipText(detail)
            into.add(label)
        }
        if(cell!=null&&evaluable&&checks.ringMeasured&&!cell.checked&&cell.why){if(into.componentCount>0)into.add(CalciumTheme.label('·',CalciumTheme.MUTED));into.add(CalciumTheme.label('passing object not checked: '+cell.why,CalciumTheme.MUTED))}
        into.visible=results!=null&&into.componentCount>0
    }
    /** The row of mark controls, while it has one to show: the button, Remove mark or the form. */
    void showMarkRow(){markRow.visible=results!=null&&(markButton.visible||removeMark.visible||formCell!=null);reviewBox.revalidate();reviewBox.repaint()}
    /** The room kept under the pickers while results are shown: three lines of the call line (a mark with the longest note
     *  takes three at the window's usual width), a row of flags, a row of mark controls and two lines of a refusal (one that
     *  names the percent table of a recording). The rows take what they need of it, one under the other; the rest is under them. */
    int headKept(){int line=cellStatus.getFontMetrics(cellStatus.font).height,refusal=markProblem.getFontMetrics(markProblem.font).height
        int controls=(int)[markButton,removeMark,reasonBox,noteField,saveMark,cancelMark].collect{it.preferredSize.height}.max()
        3*line+line+controls+4*CalciumRowLayout.PAD+2*refusal}
    /** What the mark controls can do for a cell: can, when it has a call to review on the primary stimulus, no analysis is
     *  running and the review file is applied; marked, when it has a mark then. */
    Map markable(String name) {
        Map row=name==null?null:callRows[name];String auto=row!=null&&row.Role!='Control'?automatic[name] as String:null
        boolean can=!busy&&results!=null&&auto in CalciumReview.CALLS&&!reviewProblem
        [can:can,marked:can&&(reviewMarks?.get(name) as Map)?.automatic==auto]
    }
    /** A flag's link: the Traces tab, the recording rather than the still colour picture, at the frame. A flag on the Around the
     *  cell tab (around) stays on that tab and moves its time cursor there. */
    void goToFrame(int first,boolean around=false){if(!around){tabs.selectedIndex=1;if(movie!=null&&view.selectedIndex==0)view.selectedIndex=2};stopPlaying();showFrame(first)}
    /** A refusal under the mark controls, as CalciumReview words it; null takes it away. */
    void say(String text,String cell){problemCell=text==null?null:cell;problemSaid=text;markProblem.say(text==null?'&nbsp;':CalciumOutlineProgress.escape(text));markProblem.visible=text!=null;reviewBox.revalidate();reviewBox.repaint()}
    /** The form that sets a mark, in place of the button: the reason (chosen already for a responder with a flag), an optional
     *  note, Save and Cancel. The reviewed call is the other one of the two, so it is not asked for. */
    void openForm(String name) {
        Map can=markable(name);if(!can.can||can.marked)return
        Map cell=((checks?.cells?:[:]) as Map)[name] as Map
        String chosen=automatic[name]=='Responded'?((cell?.codes?:[]) as List).collect{REASON_OF_CODE[it]}.find{it}:null
        formCell=name;enterPressed=false;enterTaken=false;reasonBox.selectedIndex=chosen?CalciumReview.REASONS.indexOf(chosen):-1;noteField.text=''
        [reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{it.visible=true};markButton.visible=false;say(null,null);showMarkRow()
        takeFocus.call(chosen?noteField:reasonBox)
    }
    /** Closes the form; the button comes back in its place, and the focus with it when the focus was in the form (back). */
    void closeForm(boolean back=true){boolean had=back&&formFocused();formCell=null;[reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{it.visible=false}
        Map can=markable(cells.selectedItem as String);markButton.visible=can.can&&!can.marked;say(null,null);showMarkRow();refocus(had)}
    /** The component that has the keyboard focus (focused; a closure, where the window's own "focused" is not this). */
    def focusOwner(){focused.call()}
    /** Whether the keyboard focus is on the mark controls. */
    boolean formFocused(){def owner=focused.call();owner instanceof java.awt.Component&&javax.swing.SwingUtilities.isDescendingFrom(owner as java.awt.Component,markRow)}
    /** The focus back on the button of the selected cell, Remove mark or the mark button, when it was on the mark controls. */
    void refocus(boolean had){if(!had)return;def button=removeMark.visible?removeMark:markButton.visible?markButton:null;if(button!=null)takeFocus.call(button)}
    /** Saves the mark of the open form (CalciumReview.change): checked first with CalciumReview's own rules and words, which
     *  then decide again. A refusal stays under the form, and nothing is written. */
    void saveReviewMark() {
        String name=formCell;if(name==null)return
        if(busy){say('Wait for the analysis to finish.',name);return}
        String auto=automatic[name] as String,call=CalciumReview.CALLS.find{it!=auto} as String,reason=reasonBox.selectedItem as String,note=CalciumReview.tidy(noteField.text?:'')
        String why=CalciumReview.refused(name,call,reason,note,automatic,CalciumReview.shown(CalciumReview.order(results,measurement),name))
        if(why){say(why,name);return}
        try{CalciumReview.change(results,measurement,[cell:name,call:call,reason:reason,note:note],traces)}catch(Throwable t){say(t.message?:t.toString(),name);return}
        boolean had=formFocused();closeForm(false);afterReview(name);refocus(had)
    }
    void removeReviewMark(String name) {
        if(name==null||results==null)return
        if(busy){say('Wait for the analysis to finish.',name);return}
        boolean had=formFocused()
        try{CalciumReview.change(results,measurement,[cell:name,remove:true],traces)}catch(Throwable t){say(t.message?:t.toString(),name);return}
        afterReview(name);refocus(had)
    }
    /** After a mark is set or removed, without loading the analysis again: the tiles, the Show list, the picture, the heatmap,
     *  the plot and the calls and counts tabs follow the marks; Show stays as it was, and the cell stays selected (kept in the
     *  Cell list while it is, when its call left the group); the time cursor stays where it is, and the finish message follows. */
    void afterReview(String name) {
        Map review=CalciumReview.read(results,measurement);reviewMarks=review.marks as Map;reviewProblem=(review.problem?:'') as String
        try{summarize(results)}catch(Throwable ignored){}
        reviewTables()
        kept=name
        selectWindow()
        cells.setSelectedItem(name)
        String done=finished();if(done!=null)status.text=done
    }
    void exportPlots(){
        if(results==null)return
        changed()
        // the saved trace has no time cursor or its reading: only the axes, the traces and the key, which names the dashed
        // lines and the bands too (on screen pointing at a band names it)
        double at=plot.moment;plot.moment=Double.NaN;aroundPlot.moment=Double.NaN
        try{
            // what is shown around the cell, read now when it is still being read, so the picture is of the cell selected
            if(aroundContext==null||aroundFolder!=results.name){aroundAsked++;aroundContext=CalciumCallChecks.aroundContext(results,measurement);aroundFolder=results.name;showAround(cells.selectedItem as String)}
            // Around the cell is saved when it has a plot (a sentence is no figure); the picture of an earlier cell is taken away
            File around=new File(results,'around_the_cell.png');boolean aroundSaved=aroundPlot.drawn()
            if(!aroundSaved)around.delete()
            // each drawn across the whole picture, on white, at the size saved (the heatmap with every row), then given its size back
            Map panels=['trace':plot,'heatmap':heatmap];if(aroundSaved)panels.around_the_cell=aroundPlot
            panels.each{name,panel->
                int w=Math.max(650,panel.width),h=Math.max(panel.is(aroundPlot)?360:260,Math.max(panel.height,panel.is(heatmap)?(int)panel.preferredSize.height:0))
                java.awt.Dimension was=panel.size;def image=new java.awt.image.BufferedImage(w,h,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics()
                boolean trace=panel.is(plot)
                try{g.setColor(Color.white);g.fillRect(0,0,w,h);panel.setSize(w,h);if(trace){plot.figure=true;plot.changed()};if(panel.is(aroundPlot))aroundPlot.changed();panel.paint(g)}
                finally{g.dispose();panel.setSize(was);if(trace){plot.figure=false;plot.changed()};if(panel.is(aroundPlot))aroundPlot.changed()}
                javax.imageio.ImageIO.write(image,'png',new File(results,name+'.png'))}
            Map provenance=[result_folder:results.name,cell:cells.selectedItem,stimulus:stimuli.selectedItem,protocol_hash:appliedManifest.protocol_hash,manifest_sha256:CalciumInstaller.digest(new File(results,'analysis_manifest.csv').bytes),settings_pending:fingerprint()!=appliedFingerprint,source:'Completed analysis snapshot; display is unsmoothed']
            if(aroundSaved)provenance.around_the_cell=around.name
            // the groups, colours and order follow the review marks on the primary stimulus: say so, and which file they come
            // from; a control window keeps its automatic calls, so its plots say nothing of the review
            File review=CalciumReview.file(measurement,results)
            if(reviewMarks&&review?.isFile()&&callRows.values().any{it.Role!='Control'})provenance+=[calls:'after review',review_file:review.name,review_sha256:CalciumInstaller.digest(review.bytes)]
            new File(results,'plot_provenance.json').setText(CalciumResponseProtocol.json(provenance)+'\n','UTF-8')
            status.text='Saved completed-run plots and plot_provenance.json in '+results.name
        }catch(Throwable t){status.text='Could not save plots: '+t.message}
        finally{plot.moment=at;plot.repaint();aroundPlot.moment=at;aroundPlot.repaint()}
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
    /** What the page last sent for this recording. A version 2 rule that names its noise estimator at the top (sd_estimator)
     *  keeps it. One that names none was written by a build that had the baseline SD only, so nobody chose it: like a rule
     *  still on the former default length, it takes what a new protocol takes, and the status line says so.
     *  @param note when given, note.ruleMoved says the remembered rule was on the former default and now counts samples, and
     *  note.noiseMoved that it named no estimator and now takes its noise from what a new protocol does */
    static Map remembered(Map data,Map note=null) {
        File file=rememberedFile(data);if(file==null||!file.isFile())return null
        try{def saved=new groovy.json.JsonSlurper().parseText(file.getText('UTF-8'))
            if(!(saved instanceof Map)||saved.recording!=''+data.recording||(saved.frames as int)!=(data.frames as int)||(saved.baseline_start as int)!=fijiFrame(data,'Start')||(saved.baseline_end as int)!=fijiFrame(data,'End'))return null
            if(!(saved.values instanceof Map)||!(saved.protocol instanceof Map))return null
            boolean moved=!(saved.format instanceof Number&&(saved.format as int)>=REMEMBERED_FORMAT)&&toSampleRule(saved.protocol as Map)
            Map candidate=[values:saved.values,protocol:saved.protocol,timestamps:null,sensitivity:null,negative_control:null]
            if(diagnose(candidate,data).errors)return null
            Map kept=candidate.protocol as Map;String now=CalciumResponseProtocol.newEstimator()
            boolean unnamed=kept.schema_version==2&&kept.sd_estimator==null,noiseMoved=unnamed&&CalciumResponseProtocol.estimator(kept,true)!=now
            if(unnamed)CalciumResponseProtocol.forPage(kept,now)
            if(note!=null){note.ruleMoved=moved;note.noiseMoved=noiseMoved}
            return tidyAutomatic(copy(candidate))
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
     *  The review marks of the analysis (CalciumReview) are applied here and nowhere else. With marks, every row also gives the
     *  figures after review: "reviewed" holds counted and without as above and "changed", the cells of the row whose call a
     *  mark changed; the tally then says reviewed: true and lists in "changed" the cells of all cells whose call a mark changed.
     *  Without marks the tally is as it always was.
     *  @param traces dF/F0 by cell when already read; otherwise normalized_traces.csv is read
     *  @param review the review marks as CalciumReview.read gives them; read from beside the tables when not given
     *  @return null when the analysis has no calls */
    static Map tally(File dir,Map measurement,Map traces=null,Map review=null) {
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
            Map marks=((review!=null?review:CalciumReview.read(dir,measurement))?.marks?:[:]) as Map;boolean reviewing=!marks.isEmpty()
            def blank={reviewing ? [counted:[0,0,0],without:[0,0,0],after:[0,0,0],afterWithout:[0,0,0],changed:0] : [counted:[0,0,0],without:[0,0,0]]}   // responded, evaluable, registered
            Map rows=[all:blank(),blind:blank(),later:blank(),carried:blank(),quiet:blank()];List under=[],changed=[]
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
                if(reviewing){   // the same cell after review: a mark moves it between the two calls and nowhere else
                    String after=CalciumReview.reviewed(marks,name,call);boolean moved=after!=call
                    if(moved)changed<<name
                    into.each{String key->List both=[rows[key].after];if(!below)both<<rows[key].afterWithout
                        both.each{List t->t[2]++;if(evaluable)t[1]++;if(after=='Responded')t[0]++};if(moved)rows[key].changed++}
                }
            }
            def shown={String key,String label->Map r=rows[key] as Map;Map row=[key:key,label:label,counted:r.counted,without:r.without];if(reviewing)row.reviewed=[counted:r.after,without:r.afterWithout,changed:r.changed];row}
            List out=[shown('all','All cells')]
            if(tags!=null)out+=[shown('blind','Chosen blind'),shown('later','Added after the response was shown')]
            if(carried!=null)out<<shown('carried',CARRIED_ROW)
            if(need>=0)out<<shown('quiet','Quiet enough to detect a '+Math.round(100*QUIET_RISE)+' % rise')
            Map found=[quietRise:QUIET_RISE,belowZero:under,tagged:tags!=null,rows:out]
            if(reviewing){found.reviewed=true;found.changed=changed}
            return found
        }catch(Throwable ignored){return null}
    }
    static String percentOf(List t){(t[1] as int)>0?String.format(Locale.US,'%.1f',100d*(t[0] as int)/(t[1] as int)):''}
    static final List PERCENT_COLUMNS=['cells','registered','evaluable','responded','percent_responding','cells_below_zero','evaluable_without_them','responded_without_them','percent_without_them']
    /** Appended on the right while review marks are applied. A mark moves a cell between responded and did not respond only, so
     *  the registered and evaluable cells of every row stay as they are. */
    static final List REVIEW_COLUMNS=['responded_after_review','percent_responding_after_review','responded_without_them_after_review','percent_without_them_after_review','cells_reclassified']
    static File percentFile(File dir,Map measurement){measurement?.folder instanceof File?new File(measurement.folder as File,measurement.recording+'_percent_responding_'+dir.name+'.csv'):null}
    /** The rows of the percent table: with review marks the same rows and columns, then REVIEW_COLUMNS; without them the
     *  table is as it always was. */
    static List percentTable(Map found) {
        boolean reviewed=found.reviewed==true
        [PERCENT_COLUMNS+(reviewed?REVIEW_COLUMNS:[])]+(found.rows as List).collect{Map row->
            List t=row.counted as List,w=row.without as List;List line=[row.label,t[2],t[1],t[0],percentOf(t),(t[2] as int)-(w[2] as int),w[1],w[0],percentOf(w)]
            if(reviewed){Map r=row.reviewed as Map;line+=[(r.counted as List)[0],percentOf(r.counted as List),(r.without as List)[0],percentOf(r.without as List),r.changed]}
            line}
    }
    /** The tally beside the measurement tables, as "<recording>_percent_responding_<analysis>.csv" (percentTable). */
    static void savePercent(File dir,Map measurement,Map found) {
        if(found==null||!(measurement?.folder instanceof File))return
        CalciumCsv.write(percentFile(dir,measurement),percentTable(found))
    }
    /** Brings the percent table in line with the review marks as they are now, when the two disagree: a review file taken
     *  away, or one that no longer matches, leaves after-review columns behind. A table with no marks and no such columns is
     *  not touched, nor is one that is already right. */
    static void syncPercent(File dir,Map measurement,Map found) {
        try {
            File target=percentFile(dir,measurement);if(found==null||target==null)return
            boolean columns=target.isFile()&&(target.withReader('UTF-8'){it.readLine()}?:'').split(',').contains(REVIEW_COLUMNS[0])
            if(found.reviewed!=true&&!columns)return
            String want=percentTable(found).collect{List row->row.collect{CalciumCsv.encode(it)}.join(',')+'\n'}.join('')
            if(target.isFile()&&target.getText('UTF-8')==want)return
            savePercent(dir,measurement,found)
        }catch(Throwable ignored){}
    }
    /** Saves the tally beside the measurement tables and returns the sentence shown when the analysis finishes (percentSentence). */
    static String percentResponding(File dir,Map measurement) {
        try {
            Map found=tally(dir,measurement);if(found==null)return ''
            savePercent(dir,measurement,found)
            return percentSentence(found)
        }catch(Throwable ignored){return ''}
    }
    /** The sentence shown when the analysis finishes: all cells first, then the two groups when cells were added after the
     *  response was shown. '' when it was not recorded how cells were chosen. With review marks every figure is followed by
     *  the figure after review, and all cells are given in any case. Writes nothing. */
    static String percentSentence(Map found) {
        try {
            if(found==null)return ''
            List rows=found.rows as List
            boolean reviewed=found.reviewed==true
            def row={String key->rows.find{it.key==key} as Map}
            def of={String key->row(key)?.counted as List}
            def after={String key->if(!reviewed)return '';List r=(row(key).reviewed as Map).counted as List;', after review '+r[0]+' of '+r[1]+((r[1] as int)>0?' ('+percentOf(r)+' %)':'')}
            def line={String label,String key->List t=of(key);label+': '+t[0]+' of '+t[1]+' responded'+((t[1] as int)>0?' ('+percentOf(t)+' %)':'')+after(key)}
            if(row('carried')!=null){   // outlines carried over: never among the cells chosen blind
                List parts=[line('All cells','all')]
                if(found.tagged){parts<<line('Chosen blind','blind');if((of('later')[2] as int)>0)parts<<line('Added after the response was shown','later')}
                parts<<line(CARRIED_ROW,'carried')
                return ' '+parts.join('. ')+'.'
            }
            if(!found.tagged&&!reviewed)return ''
            if(!found.tagged)return ' '+line('All cells','all')+'.'
            List blind=of('blind'),later=of('later')
            if((later[2] as int)==0)return ' All '+blind[2]+' cells were chosen blind: '+blind[0]+' of '+blind[1]+' responded ('+percentOf(blind)+' %)'+after('blind')+'.'
            return ' '+line('All cells','all')+'. '+line('Chosen blind','blind')+'. '+line('Added after the response was shown','later')+'.'
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
    /** Read-only feedback for incomplete drafts from the page; the same authoritative checks also run on submission. */
    static Map diagnose(Map draft,Map measurement) {
        List errors=[],warnings=[]
        try {
            if(!(draft?.values instanceof Map)||!(draft?.protocol instanceof Map))throw new IllegalArgumentException('Complete the recording settings and response protocol.')
            Map flags=[:],uploads=[:]
            FILE_KEYS.each{k,v->uploads[k]=upload(draft[k],k);if(uploads[k])flags[v]=true}
            check(draft.values as Map,draft.protocol as Map,true,measurement,flags,true)
            double[] times=uploads.timestamps?CalciumMeasurement.readTimesRows(CalciumCsv.parse(uploads.timestamps.content as String),measurement.frames as int):
                ((0..<(measurement.frames as int)).collect{it/CalciumCsv.number(draft.values.sampling_hz.toString())} as double[])
            Map protocol=CalciumResponseProtocol.canonical(draft.protocol,true)
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
    /** A finished run where no cell got a call for the primary stimulus is not a success: says why, from the saved
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
            return 'Finished, but no cell got a call for '+stimulus+': '+why+'. Check the settings; the results were saved for the record.'
        } catch(Throwable ignored) { return null }
    }
    /** A copy is necessary: neither an incoming HTTP map nor later edits can mutate a submitted run. */
    static Map copy(Map draft){new groovy.json.JsonSlurper().parseText(groovy.json.JsonOutput.toJson(draft)) as Map}
    /** A draft as the page is given it: a copy, its protocol in the words the page sends back (CalciumResponseProtocol.forPage). */
    static Map forPage(Map draft){Map out=copy(draft);if(out.protocol instanceof Map)CalciumResponseProtocol.forPage(out.protocol as Map);out}
    /** The option "Mark possible passing objects as not responders" (CalciumReview.automatic): a setting of the analysis, kept
     *  in the draft's values and remembered with them (remember), never given to the engine and not in its fingerprint. On:
     *  true; off: no key at all, so a draft with the option off is what it always was. */
    static final String AUTOMATIC_MARKS='automatic_marks'
    static boolean automaticMarks(Map draft){draft?.values instanceof Map&&(draft.values as Map).get(AUTOMATIC_MARKS)==true}
    /** The draft with the option as it is kept: true, or no key. */
    static Map tidyAutomatic(Map draft){if(draft?.values instanceof Map&&(draft.values as Map).containsKey(AUTOMATIC_MARKS)&&(draft.values as Map).get(AUTOMATIC_MARKS)!=true)(draft.values as Map).remove(AUTOMATIC_MARKS);draft}
    /** What the engine is told the draft was (draft_sha256): without the option of the automatic marks, which it never sees. */
    static String fingerprint(Map draft) {
        def sorted;sorted={v->v instanceof Map?new TreeMap(v.collectEntries{k,x->[''+k,sorted(x)]}):v instanceof Collection?v.collect{sorted(it)}:v}
        Map plain=draft
        if(draft?.values instanceof Map&&(draft.values as Map).containsKey(AUTOMATIC_MARKS)){plain=new LinkedHashMap(draft);plain.values=(draft.values as Map).findAll{it.key!=AUTOMATIC_MARKS}}
        CalciumInstaller.digest(groovy.json.JsonOutput.toJson(sorted(plain)).getBytes('UTF-8'))
    }
    /** Acquisition/source checks are shared with the native window. Protocol response math stays in the pinned R engine.
     *  @param page the protocol is a draft from the Step 4 page (CalciumResponseProtocol.estimator) */
    static Map check(Map typed,Map protocol,boolean confirmed,Map measurement,Map files=[:],boolean page=false) {
        if(!confirmed)throw new IllegalArgumentException('Confirm the frame rate/timestamps, baseline windows and cell selection first.')
        Map values=VALUE_KEYS.collectEntries{k->[k,full(typed[k]).trim()]}
        double first=CalciumResponseProtocol.numeric(values.baseline_start,'normalization first frame'),last=CalciumResponseProtocol.numeric(values.baseline_end,'normalization last frame')
        if(first<1||first!=Math.floor(first)||last!=Math.floor(last)||last<first||last>=(measurement.frames as int))throw new IllegalArgumentException('Check the normalization baseline frame range.')
        if(files.timestamps_file){values.sampling_hz=''}else{double hz=CalciumResponseProtocol.numeric(values.sampling_hz,'frame rate');if(hz<=0)throw new IllegalArgumentException('Use a positive frame rate.')}
        if(protocol==null)throw new IllegalArgumentException('A portable response protocol is required; choose a version 2 or explicit legacy protocol.')
        Map canonical=CalciumResponseProtocol.canonical(protocol,page)
        if(canonical.schema_version==1&&(files.sensitivity_file||files.negative_control_file))throw new IllegalArgumentException('Sensitivity and negative-control diagnostics require a version 2 protocol.')
        values.protocol_json=groovy.json.JsonOutput.toJson(canonical);values.sd_estimator=CalciumResponseProtocol.noiseOf(canonical)   // the runner stops when this is not the protocol's own
        values.recording=measurement.recording;values.input_mode='background_corrected_fluorescence';values.selection_basis=measurement.selectionBasis?:'Not recorded'
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

/** What is measured around and inside every cell for the check that looks for a bright object passing through an outline
 *  (CalciumCallChecks; GitHub issue 29). Each cell gets 25 zones, fixed pixel sets built once from the outlines' own pixel
 *  masks (CalciumSteps.pixelsOf):
 *   - a near ring and a far ring outside the outline, each in 8 sectors of 45 degrees centred on right, lower right, bottom,
 *     lower left, left, upper left, top and upper right, seen from the mean x, y of the outline's pixels. Near ring: more than
 *     NEAR_FROM_UM and up to NEAR_TO_UM from the outline. Far ring: beyond that, up to FAR_TO_UM. Both leave out every pixel
 *     within OTHER_UM of another measured outline and keep bright pixels;
 *   - core: the pixels inside that lie deeper than the median depth. Rim: the other pixels inside, by sector.
 *  Distances are Euclidean, pixel centre to pixel centre: outside, to the nearest pixel of the outline; inside, to the nearest
 *  pixel that is not in the outline, and the picture's border counts as outside. A zone with fewer than MIN_PIXELS pixels is
 *  not used: its column stays empty.
 *  The widths were set in pixels (1, 10, 22 and 2) on recordings at SET_AT_UM micrometres per pixel. They are kept in
 *  micrometres and turned into pixels with the recording's own pixel size, to a hundredth of a pixel. A recording without a
 *  pixel size in micrometres takes the pixel values as they are.
 *  Measured: the mean raw gray value of every zone on every frame, in one pass, written as it goes. Nothing here enters a
 *  cell's own measurement or its background. */
@CompileStatic
class CalciumRingAndCore {
    static final double SET_AT_UM = 0.645d
    static final double NEAR_FROM_UM = 0.645d, NEAR_TO_UM = 6.45d, FAR_TO_UM = 14.19d, OTHER_UM = 1.29d
    static final int MIN_PIXELS = 12, ZONES = 25, CORE = 16
    static final List<String> SIDES = ['right', 'lower right', 'bottom', 'lower left', 'left', 'upper left', 'top', 'upper right']
    /** near_right ... near_upper_right, far_right ... far_upper_right, core, rim_right ... rim_upper_right */
    static final List<String> NAMES = zoneNames()
    static List<String> zoneNames() {
        List<String> out = new ArrayList<String>()
        for (String ring : ['near', 'far']) for (String side : SIDES) out.add(ring + '_' + side.replace(' ', '_'))
        out.add('core')
        for (String side : SIDES) out.add('rim_' + side.replace(' ', '_'))
        return out
    }
    /** How the widths are turned into pixels: 'micrometres' for a positive pixel size in micrometres that is the same both ways,
     *  as the rings are round; 'not square' for micrometres with pixels that are not square; 'no size' otherwise. The last two
     *  take the pixel values as they are. The unit is tested as for the scale bar (CalciumSteps.micrometres). Groovy counts NaN
     *  as larger than any number, so it is tested for first. */
    static String scaleOf(double pixelWidth, double pixelHeight, String unit) {
        if (!CalciumSteps.micrometres(unit) || Double.isNaN(pixelWidth) || Double.isInfinite(pixelWidth) || pixelWidth <= 0d) return 'no size'
        if (Double.isNaN(pixelHeight) || Double.isInfinite(pixelHeight) || Math.abs(pixelHeight - pixelWidth) > 0.01d * pixelWidth) return 'not square'
        return 'micrometres'
    }
    /** Whether a recording's pixel size can be used (scaleOf). */
    static boolean inMicrometres(double pixelWidth, double pixelHeight, String unit) { return scaleOf(pixelWidth, pixelHeight, unit) == 'micrometres' }
    /** The first line of the values table: Frame, then the 25 zones of every cell as "cell|zone". */
    static String header(List<String> names) {
        StringBuilder head = new StringBuilder('Frame')
        for (String name : names) for (String zone : NAMES) head.append(',').append(name).append('|').append(zone)
        return head.append('\n').toString()
    }
    /** Where the near ring starts, where it ends and the far ring starts, where the far ring ends, and the margin around other
     *  outlines, in pixels. */
    static double[] widths(double pixelWidth, double pixelHeight, String unit) {
        double size = inMicrometres(pixelWidth, pixelHeight, unit) ? pixelWidth : SET_AT_UM
        double[] um = [NEAR_FROM_UM, NEAR_TO_UM, FAR_TO_UM, OTHER_UM] as double[], out = new double[4]
        for (int i = 0; i < 4; i++) out[i] = Math.rint(um[i] / size * 100d) / 100d
        return out
    }

    final int cells
    /** Pixel indices of every zone, cell after cell: zone z of cell c is at c * ZONES + z. */
    final int[][] zonePx
    final double[] width
    final boolean micrometres
    /** scaleOf: 'micrometres', 'not square' or 'no size'. */
    final String scale
    final double pixelSize
    private Writer out
    private File temp, target
    private final StringBuilder line = new StringBuilder(1 << 16)

    CalciumRingAndCore(int[][] cellPx, int w, int h, double pixelWidth, double pixelHeight, String unit) {
        cells = cellPx.length
        scale = scaleOf(pixelWidth, pixelHeight, unit); micrometres = scale == 'micrometres'; pixelSize = micrometres ? pixelWidth : Double.NaN
        width = widths(pixelWidth, pixelHeight, unit)
        int[] cover = new int[w * h]                      // how many outlines hold each pixel
        for (int[] px : cellPx) for (int i : px) cover[i]++
        zonePx = new int[cells * ZONES][]
        for (int c = 0; c < cells; c++) { int[][] z = zonesOf(cellPx[c], cover, w, h, width); for (int k = 0; k < ZONES; k++) zonePx[c * ZONES + k] = z[k] }
    }
    int pixels(int cell, int zone) { return zonePx[cell * ZONES + zone].length }
    /** How many of the 8 near-ring sectors of a cell hold enough pixels to be used. */
    int nearUsable(int cell) { int k = 0; for (int z = 0; z < 8; z++) if (pixels(cell, z) >= MIN_PIXELS) k++; return k }

    /** Squared distance from every pixel to the nearest pixel of "feature" (exact, Meijster's two passes), or null when there is none. */
    static int[] distances(boolean[] feature, int cw, int ch) {
        boolean any = false
        for (boolean b : feature) if (b) { any = true; break }
        if (!any) return null
        int inf = cw + ch
        int[] g = new int[cw * ch]
        for (int x = 0; x < cw; x++) {
            g[x] = feature[x] ? 0 : inf
            for (int y = 1; y < ch; y++) { int i = y * cw + x; g[i] = feature[i] ? 0 : Math.min(inf, g[i - cw] + 1) }
            for (int y = ch - 2; y >= 0; y--) { int i = y * cw + x; if (g[i + cw] + 1 < g[i]) g[i] = g[i + cw] + 1 }
        }
        int[] d = new int[cw * ch], s = new int[cw], t = new int[cw]
        for (int y = 0; y < ch; y++) {
            int row = y * cw, q = 0
            s[0] = 0; t[0] = 0
            for (int u = 1; u < cw; u++) {
                long gu = g[row + u]
                while (q >= 0) {
                    long gi = g[row + s[q]], a = t[q] - s[q], b = t[q] - u
                    if (a * a + gi * gi > b * b + gu * gu) q--; else break
                }
                if (q < 0) { q = 0; s[0] = u }
                else {
                    long i = s[q], gi = g[row + s[q]]
                    long at = Math.floorDiv((long) u * u - i * i + gu * gu - gi * gi, 2L * (u - i)) + 1L
                    if (at < cw) { q++; s[q] = u; t[q] = (int) at }
                }
            }
            for (int u = cw - 1; u >= 0; u--) {
                long a = u - s[q], gi = g[row + s[q]]
                d[row + u] = (int) Math.min((long) Integer.MAX_VALUE, a * a + gi * gi)
                if (u == t[q]) q--
            }
        }
        return d
    }
    /** The sector (0 = right, then clockwise on the screen: lower right, bottom ...) of a pixel dy below and dx right of the centre. */
    static int sector(double dy, double dx) { return Math.floorMod((int) Math.rint(Math.toDegrees(Math.atan2(dy, dx)) / 45d), 8) }

    /** The 25 zones of one cell as pixel indices. "cover" counts the outlines that hold each pixel of the picture. */
    static int[][] zonesOf(int[] px, int[] cover, int w, int h, double[] width) {
        int[][] zones = new int[ZONES][]
        if (px.length == 0) { for (int z = 0; z < ZONES; z++) zones[z] = new int[0]; return zones }
        int xmin = w, xmax = -1, ymin = h, ymax = -1
        long sx = 0L, sy = 0L
        for (int i : px) {
            int y = Math.floorDiv(i, w), x = i - y * w
            if (x < xmin) xmin = x
            if (x > xmax) xmax = x
            if (y < ymin) ymin = y
            if (y > ymax) ymax = y
            sx += x; sy += y
        }
        int pad = (int) Math.ceil(width[2] + width[3] + 2d)
        int x0 = Math.max(0, xmin - pad), x1 = Math.min(w, xmax + pad + 1), y0 = Math.max(0, ymin - pad), y1 = Math.min(h, ymax + pad + 1)
        int cw = x1 - x0, ch = y1 - y0, pw = cw + 2
        boolean[] inside = new boolean[cw * ch], other = new boolean[cw * ch], outside = new boolean[pw * (ch + 2)]
        for (int i : px) { int y = Math.floorDiv(i, w); inside[(y - y0) * cw + (i - y * w) - x0] = true }
        Arrays.fill(outside, true)                        // one more pixel all round counts as outside
        for (int yy = 0; yy < ch; yy++) for (int xx = 0; xx < cw; xx++) {
            int k = yy * cw + xx
            if (inside[k]) outside[(yy + 1) * pw + xx + 1] = false
            if (cover[(y0 + yy) * w + x0 + xx] - (inside[k] ? 1 : 0) > 0) other[k] = true
        }
        int[] toOutline = distances(inside, cw, ch), depth = distances(outside, pw, ch + 2), toOther = distances(other, cw, ch)
        double cy = ((double) sy) / px.length - y0, cx = ((double) sx) / px.length - x0
        double[] depths = new double[px.length]
        int n = 0
        for (int yy = 0; yy < ch; yy++) for (int xx = 0; xx < cw; xx++) if (inside[yy * cw + xx]) depths[n++] = Math.sqrt((double) depth[(yy + 1) * pw + xx + 1])
        Arrays.sort(depths)
        double middle = (n & 1) == 1 ? depths[n >> 1] : (depths[(n >> 1) - 1] + depths[n >> 1]) / 2d
        int[] zoneOf = new int[cw * ch], count = new int[ZONES]
        for (int yy = 0; yy < ch; yy++) for (int xx = 0; xx < cw; xx++) {
            int k = yy * cw + xx, z = -1
            if (inside[k]) z = Math.sqrt((double) depth[(yy + 1) * pw + xx + 1]) > middle ? CORE : CORE + 1 + sector(yy - cy, xx - cx)
            else if (toOther == null || Math.sqrt((double) toOther[k]) > width[3]) {
                double far = Math.sqrt((double) toOutline[k])
                if (far > width[0] && far <= width[1]) z = sector(yy - cy, xx - cx)
                else if (far > width[1] && far <= width[2]) z = 8 + sector(yy - cy, xx - cx)
            }
            zoneOf[k] = z
            if (z >= 0) count[z]++
        }
        for (int z = 0; z < ZONES; z++) { zones[z] = new int[count[z]]; count[z] = 0 }
        for (int yy = 0; yy < ch; yy++) for (int xx = 0; xx < cw; xx++) { int z = zoneOf[yy * cw + xx]; if (z >= 0) zones[z][count[z]++] = (y0 + yy) * w + x0 + xx }
        return zones
    }

    /** Starts the table (gzip, one row per frame): Frame, then the 25 zones of every cell as "cell|zone". It is written beside
     *  its final name, under a name of its own so two measurements at once do not write into one file, and moved there by close(). */
    void open(File file, String[] names) {
        target = file; temp = new File(file.getPath() + '.' + UUID.randomUUID().toString().substring(0, 8) + '.part')
        out = new BufferedWriter(new OutputStreamWriter(new java.util.zip.GZIPOutputStream(new FileOutputStream(temp), 1 << 16), 'UTF-8'), 1 << 16)
        out.write(header(Arrays.asList(names)))
    }
    /** One frame: the mean raw gray value of every zone, two decimals; empty where the zone has too few pixels. */
    void frame(int number, ImageProcessor ip) {
        StringBuilder sb = line
        sb.setLength(0); sb.append(number)
        Object pixels = ip.getPixels()
        short[] s16 = pixels instanceof short[] ? (short[]) pixels : null
        float[] f32 = pixels instanceof float[] ? (float[]) pixels : null
        for (int z = 0; z < zonePx.length; z++) {
            int[] px = zonePx[z]
            sb.append(',')
            if (px.length < MIN_PIXELS) continue
            double sum = 0d
            if (s16 != null) { long total = 0L; for (int i : px) total += (s16[i] & 0xffff); sum = (double) total }
            else if (f32 != null) for (int i : px) sum += f32[i]
            else for (int i : px) sum += ip.getf(i)
            twoDecimals(sb, sum / px.length)
        }
        sb.append('\n')
        out.append(sb)
    }
    static void twoDecimals(StringBuilder sb, double v) {
        if (Double.isNaN(v) || Double.isInfinite(v)) return
        long cents = Math.round(Math.abs(v) * 100d)
        if (v < 0d && cents != 0L) sb.append('-')
        int rest = (int) Math.floorMod(cents, 100L)
        sb.append(Math.floorDiv(cents, 100L)).append('.').append((char) (48 + Math.floorDiv(rest, 10))).append((char) (48 + Math.floorMod(rest, 10)))
    }
    void close() {
        out.close(); out = null
        java.nio.file.Files.move(temp.toPath(), target.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    }
    /** Gives up a table that could not be finished; nothing is left behind. */
    void abandon() { try { if (out != null) out.close() } catch (Throwable ignored) { }; out = null; if (temp != null) temp.delete() }
    /** The whole table in one pass over the frames. */
    void measure(ImageStack stack, int frames, File file, String[] names) {
        open(file, names)
        try { for (int f = 1; f <= frames; f++) frame(f, stack.getProcessor(f)); close() } catch (Throwable t) { abandon(); throw t }
    }
    /** The table read back: [column][frame], NaN where a value is empty. The header must be the one open() writes for these names. */
    static float[][] read(File file, List<String> names, int frames) {
        ByteArrayOutputStream all = new ByteArrayOutputStream(1 << 20)
        InputStream input = new java.util.zip.GZIPInputStream(new FileInputStream(file), 1 << 16)
        try { byte[] buffer = new byte[1 << 16]; int got; while ((got = input.read(buffer)) > 0) all.write(buffer, 0, got) } finally { input.close() }
        byte[] data = all.toByteArray()
        byte[] expected = header(names).getBytes('UTF-8')
        if (data.length < expected.length) throw new IOException('The ring and core table is for other cells.')
        for (int i = 0; i < expected.length; i++) if (data[i] != expected[i]) throw new IOException('The ring and core table is for other cells.')
        int columns = names.size() * ZONES, at = expected.length
        float[][] values = new float[columns][frames]
        for (int f = 0; f < frames; f++) {
            long frame = 0L
            while (at < data.length && data[at] != (byte) 44) { frame = frame * 10L + (data[at] - 48); at++ }
            if (frame != f + 1) throw new IOException('The ring and core table does not hold frame ' + (f + 1) + '.')
            for (int c = 0; c < columns; c++) {
                if (at >= data.length || data[at] != (byte) 44) throw new IOException('The ring and core table has too few columns.')
                at++
                boolean negative = false, digits = false
                long whole = 0L; double part = 0d, scale = 0.1d
                if (at < data.length && data[at] == (byte) 45) { negative = true; at++ }
                while (at < data.length && data[at] >= (byte) 48 && data[at] <= (byte) 57) { whole = whole * 10L + (data[at] - 48); at++; digits = true }
                if (at < data.length && data[at] == (byte) 46) { at++; while (at < data.length && data[at] >= (byte) 48 && data[at] <= (byte) 57) { part += (data[at] - 48) * scale; scale *= 0.1d; at++; digits = true } }
                double v = digits ? whole + part : Double.NaN
                values[c][f] = (float) (negative ? -v : v)
            }
            if (at < data.length && data[at] == (byte) 13) at++
            if (at >= data.length || data[at] != (byte) 10) throw new IOException('The ring and core table has too many columns.')
            at++
        }
        if (at != data.length) throw new IOException('The ring and core table holds more frames than the recording.')
        return values
    }
}

/** The check on one cell's call: was the rise a bright object passing through the outline, or a rise on its edge that is not
 *  in the cell? It only flags; no call changes. The rule and its numbers were set and checked on nine recordings (GitHub
 *  issue 29) and are kept as they were tested.
 *  Traces, in gray levels: the cell is raw minus its background minus its mean over the F0 frames; a zone is its mean minus
 *  its own F0 mean minus the change of the cell's background. Both are averaged over AVERAGE frames (centred; the ends use the
 *  frames that exist). The event is the stretch around the peak of the largest rise where the cell stays at or above HALF its
 *  peak; a dip of up to DIP frames does not end it. The cell's rise is its value at the peak.
 *  A zone's noise is the larger of the SD of its averaged trace over the F0 frames and 1.4826 x the median absolute deviation
 *  of the frame-to-frame differences of its unaveraged trace over the whole recording, divided by root 2 and by root AVERAGE.
 *  Passes through: a near- or far-ring sector rises in the "span" frames before the event and a sector at least APART sectors
 *  (90 degrees) away rises in the span frames after it. Each rise is at least CELL_RISE x the cell's own rise and NOISE x that
 *  sector's noise, gone while the cell is up (its mean over the event at most HALF its peak), and at a frame where the cell
 *  is still below HALF its rise.
 *  On the edge: over the event the core rises less than CORE_SHARE x the whole outline, and the near-ring sector outside the
 *  rim sector that rises most rises at least as much as that rim sector and at least NOISE x its own noise.
 *  Crosses the cell (crossing): the centre of the rise over the outline and its near ring moves across the outline, from one
 *  side of its centre to the other, at least CROSS_TRAVEL x the outline's width, steadily and well over its noise. Its shape
 *  and numbers were chosen on eight of the recordings of issue 29.
 *  The span was set as 10 frames at 1 frame per second (SPAN_S seconds). */
@CompileStatic
class CalciumPassingObject {
    static final int AVERAGE = 3, DIP = 3, APART = 2
    static final double SPAN_S = 10d, HALF = 0.5d, CELL_RISE = 1d, NOISE = 5d, CORE_SHARE = 0.25d
    /** Crosses the cell: the window is where the cell stays at or above CROSS_LEVEL of its rise; it flags a travel of at least
     *  CROSS_TRAVEL x the outline's width, a steadiness of at least CROSS_STEADY and a travel of at least CROSS_NOISE x its noise. */
    static final double CROSS_LEVEL = 0.25d, CROSS_TRAVEL = 0.5d, CROSS_STEADY = 0.6d, CROSS_NOISE = 3d
    /** The zones crossing uses: the core, the 8 rim sectors and the 8 near-ring sectors (CalciumRingAndCore.NAMES). */
    static final int[] CROSS_ZONES = [16, 17, 18, 19, 20, 21, 22, 23, 24, 0, 1, 2, 3, 4, 5, 6, 7] as int[]
    static double[] averaged(double[] x) {
        int n = x.length, k = AVERAGE >> 1
        double[] out = new double[n]
        for (int i = 0; i < n; i++) { int a = Math.max(0, i - k), b = Math.min(n, i + k + 1); double s = 0d; for (int j = a; j < b; j++) s += x[j]; out[i] = s / (b - a) }
        return out
    }
    static double mean(double[] x, int first, int last) { double s = 0d; for (int i = first; i <= last; i++) s += x[i]; return s / (last - first + 1) }
    /** Whether every value is a number. A trace with a missing value is not used: Groovy counts NaN as larger than any number,
     *  so a missing value would pass every "at least" in the rules. */
    static boolean finite(double[] x) { for (double v : x) if (Double.isNaN(v) || Double.isInfinite(v)) return false; return true }
    static boolean finite(float[] x) { for (float v : x) if (Float.isNaN(v) || Float.isInfinite(v)) return false; return true }
    /** "frame 7" or "frames 7 to 9", counted from 1. */
    static String frames(int first, int last) { return first == last ? 'frame ' + first : 'frames ' + first + ' to ' + last }
    /** The stretch around the peak where the cell stays at or above "level": [first, last], counted from 0. A dip of up to DIP
     *  frames below it does not end it. */
    static int[] event(double[] cell, int peak, double level) {
        int n = cell.length, a = peak, b = peak
        while (a > 0 && cell[a - 1] >= level) a--
        while (b < n - 1 && cell[b + 1] >= level) b++
        boolean moved = true
        while (moved) {
            moved = false
            for (int d = 2; d <= DIP + 1; d++) if (a - d >= 0 && cell[a - d] >= level) { a -= d; moved = true; while (a > 0 && cell[a - 1] >= level) a--; break }
            for (int d = 2; d <= DIP + 1; d++) if (b + d < n && cell[b + d] >= level) { b += d; moved = true; while (b < n - 1 && cell[b + 1] >= level) b++; break }
        }
        return [a, b] as int[]
    }
    /** 1.4826 x the median absolute deviation of the frame-to-frame differences, divided by root 2: the noise of one frame. */
    static double spread(double[] x) {
        int n = x.length
        double[] d = new double[n - 1]
        for (int f = 1; f < n; f++) d[f - 1] = x[f] - x[f - 1]
        double m = CalciumMotion.median(d)
        for (int f = 0; f < n - 1; f++) d[f] = Math.abs(d[f] - m)
        return 1.4826d * CalciumMotion.median(d) / Math.sqrt(2d)
    }
    /** The largest value among the frames from..to where the cell is below "half": [value, frame], or null when there is no such frame. */
    private static double[] highest(double[] zone, double[] cell, double half, int from, int to) {
        int at = -1
        for (int f = from; f <= to; f++) if (cell[f] < half && (at < 0 || zone[f] > zone[at])) at = f
        return at < 0 ? null : [zone[at], at] as double[]
    }
    /** The traces the check works on, worked out in this one place for the check (check) and for what is shown around a cell
     *  (CalciumCallChecks.aroundCell, which takes them from check): the same numbers, never a second calculation.
     *  @param raw the cell's raw gray values, @param background its background, both one per frame
     *  @param zone the 25 zone traces (CalciumRingAndCore.NAMES); null where a zone is not usable; the whole array null when
     *         nothing was measured around the cell, which still gives the event
     *  @param firstF0, lastF0 the F0 frames and @param peak the frame of the call's peak, all counted from 0
     *  @return null when a gray value is missing or the cell does not rise at the peak; else cell (its rise per frame in gray
     *          levels: raw minus background minus its F0 mean, averaged over AVERAGE frames), a and b (the event's first and
     *          last frame, counted from 0), rise, half, and trace and noise (the zones' traces, null where a zone is not usable,
     *          and their noise), both null when nothing was measured around the cell */
    static Map<String, Object> traces(double[] raw, double[] background, float[][] zone, int firstF0, int lastF0, int peak) {
        int n = raw.length
        if (!finite(raw) || !finite(background)) return null   // the caller says "missing gray values"
        double[] cell = new double[n], change = new double[n]
        for (int f = 0; f < n; f++) cell[f] = raw[f] - background[f]
        double own = mean(cell, firstF0, lastF0), ground = mean(background, firstF0, lastF0)
        for (int f = 0; f < n; f++) { cell[f] = cell[f] - own; change[f] = background[f] - ground }
        cell = averaged(cell)
        double rise = cell[peak]
        if (Double.isNaN(rise) || Double.isInfinite(rise) || rise <= 0d) return null
        double half = HALF * rise
        int[] around = event(cell, peak, half)
        int a = around[0], b = around[1]
        Map<String, Object> out = new LinkedHashMap<String, Object>()
        out.put('cell', cell); out.put('a', a); out.put('b', b); out.put('rise', rise); out.put('half', half); out.put('trace', null); out.put('noise', null)
        if (zone == null) return out
        double[][] trace = new double[CalciumRingAndCore.ZONES][]
        double[] noise = new double[CalciumRingAndCore.ZONES]
        for (int z = 0; z < CalciumRingAndCore.ZONES; z++) {
            if (zone[z] == null || zone[z].length != n || !finite(zone[z])) continue   // not usable, like a zone under MIN_PIXELS
            double[] plain = new double[n]
            double base = 0d
            for (int f = firstF0; f <= lastF0; f++) base += zone[z][f]
            base /= (lastF0 - firstF0 + 1)
            for (int f = 0; f < n; f++) plain[f] = zone[z][f] - base - change[f]
            trace[z] = averaged(plain)
            double sd = Double.NaN
            if (lastF0 > firstF0) { double m = mean(trace[z], firstF0, lastF0), ss = 0d; for (int f = firstF0; f <= lastF0; f++) ss += (trace[z][f] - m) * (trace[z][f] - m); sd = Math.sqrt(ss / (lastF0 - firstF0)) }
            double steps = Double.NaN
            if (n > 1) steps = spread(plain) / Math.sqrt((double) AVERAGE)
            noise[z] = Double.isNaN(sd) ? steps : Double.isNaN(steps) ? sd : Math.max(sd, steps)
        }
        out.put('trace', trace); out.put('noise', noise)
        return out
    }
    /** @param raw, background, zone, firstF0, lastF0, peak as for traces
     *  @param span frames searched before and after the event
     *  @param centre, outline where the zones and the outline's pixels lie, as shape gives them; null when they are not known,
     *         and then the rule "crosses the cell" is not checked
     *  @param asBefore true to choose the two sectors of "passes through" as 2026.11.24.1 did: on equal strength the first pair in
     *         zone order (CalciumCallChecks.aroundCell, for a table of that version)
     *  @return null when the cell does not rise at the peak; else first and last (the event's frames, counted from 1), rule
     *          ('passes through', 'on the edge', 'crosses the cell' or ''), text, detail; traces (what traces gave); where a
     *          passage was seen, before and after ([zone, frame counted from 1] of the sector seen before and after the event);
     *          where the rise is on an edge, edge ([rim sector 0 to 7, the stretch's first frame, its last frame, counted from
     *          1]); cross (what crossing gave; null when it was not checked); and also, when rule A or B flags the cell and the
     *          rise crosses it too, the words added for that at the end of detail ("also crosses the cell from ...") */
    static Map<String, Object> check(double[] raw, double[] background, float[][] zone, int firstF0, int lastF0, int peak, int span, double[][] centre = null, double[][] outline = null, boolean asBefore = false) {
        Map<String, Object> t = traces(raw, background, zone, firstF0, lastF0, peak)
        if (t == null) return null
        int n = raw.length, a = (int) t.get('a'), b = (int) t.get('b')
        double[] cell = (double[]) t.get('cell')
        double rise = (double) t.get('rise'), half = (double) t.get('half')
        Map<String, Object> out = new LinkedHashMap<String, Object>()
        out.put('first', a + 1); out.put('last', b + 1); out.put('rule', ''); out.put('text', ''); out.put('detail', ''); out.put('traces', t)
        if (t.get('trace') == null) return out
        double[][] trace = (double[][]) t.get('trace')
        double[] noise = (double[]) t.get('noise')
        // passes through: seen outside on one side before the event and on another side after it
        List<double[]> before = new ArrayList<double[]>(), after = new ArrayList<double[]>()
        for (int z = 0; z < 16; z++) {
            if (trace[z] == null) continue
            for (int side = 0; side < 2; side++) {
                int from = side == 0 ? Math.max(0, a - span) : b + 1, to = side == 0 ? a - 1 : Math.min(n, b + 1 + span) - 1
                if (to < from) continue
                double[] top = highest(trace[z], cell, half, from, to)
                if (top == null) continue
                double v = top[0]
                if (v < CELL_RISE * rise || v < NOISE * noise[z] || mean(trace[z], a, b) > HALF * v) continue
                (side == 0 ? before : after).add([z, v, top[1]] as double[])
            }
        }
        // the pair whose weaker side is strongest; on equal strength the one whose stronger side is stronger (asBefore: the
        // first in zone order, as 2026.11.24.1 chose)
        double[] entering = null, leaving = null
        double strength = 0d, other = 0d
        for (double[] x : before) for (double[] y : after) {
            int apart = Math.abs(((int) x[0]) % 8 - ((int) y[0]) % 8)
            apart = Math.min(apart, 8 - apart)
            if (apart < APART) continue
            double s = Math.min(x[1], y[1]) / rise, o = Math.max(x[1], y[1]) / rise
            if (entering == null || s > strength || (!asBefore && s == strength && o > other)) { entering = x; leaving = y; strength = s; other = o }
        }
        if (entering != null) {
            out.put('before', [(int) entering[0], ((int) entering[2]) + 1] as List<Integer>); out.put('after', [(int) leaving[0], ((int) leaving[2]) + 1] as List<Integer>)
            out.put('rule', 'passes through')
            out.put('text', 'possible passing object, ' + frames(a + 1, b + 1))
            out.put('detail', 'seen outside at the ' + CalciumRingAndCore.SIDES[((int) entering[0]) % 8] + ' at frame ' + (((int) entering[2]) + 1) + ', ' + String.format(Locale.US, '%.1f', entering[1] / rise) + " x the cell's rise, and at the " +
                CalciumRingAndCore.SIDES[((int) leaving[0]) % 8] + ' at frame ' + (((int) leaving[2]) + 1) + ', ' + String.format(Locale.US, '%.1f', leaving[1] / rise) + ' x')
        }
        // on the edge: the rise is in the ring outside the rim sector that rises most, and not in the core
        double whole = mean(cell, a, b)
        double core = trace[CalciumRingAndCore.CORE] == null ? Double.NaN : mean(trace[CalciumRingAndCore.CORE], a, b)
        int hot = -1
        double most = Double.NEGATIVE_INFINITY
        for (int s = 0; s < 8; s++) { if (trace[CalciumRingAndCore.CORE + 1 + s] == null) continue; double v = mean(trace[CalciumRingAndCore.CORE + 1 + s], a, b); if (hot < 0 || v > most) { hot = s; most = v } }
        if (hot >= 0 && whole > 0d && trace[hot] != null) {
            double ring = mean(trace[hot], a, b)
            if (core < CORE_SHARE * whole && ring >= most && ring >= NOISE * noise[hot]) {
                int first = a, last = b
                while (first > 0 && trace[hot][first - 1] > HALF * ring) first--
                while (last < n - 1 && trace[hot][last + 1] > HALF * ring) last++
                out.put('edge', [hot, first + 1, last + 1] as List<Integer>)
                double over = ring / most   // the rim can be flat: then the two rises are given instead
                String seen = 'core rises ' + String.format(Locale.US, '%.2f', core / whole) + ' x the whole outline, ring ' +
                    (Double.isNaN(over) || Double.isInfinite(over) ? String.format(Locale.US, 'rises %.1f and the rim %.1f gray levels', ring, most) : String.format(Locale.US, '%.1f', over) + ' x the rim')
                if (entering == null) {
                    out.put('rule', 'on the edge')
                    out.put('text', 'rise on the ' + CalciumRingAndCore.SIDES[hot] + ' edge, not in the cell, ' + frames(first + 1, last + 1))
                    out.put('detail', seen)
                } else out.put('detail', ((String) out.get('detail')) + '; also on the ' + CalciumRingAndCore.SIDES[hot] + ' edge: ' + seen)
            }
        }
        // crosses the cell: the centre of the rise moves across the outline; when rule A or B flags the cell, said as well
        Map<String, Object> cross = centre == null || outline == null ? null : crossing(cell, peak, background, zone, centre, outline, firstF0, lastF0)
        out.put('cross', cross)
        if (cross != null && cross.get('flag') == true) {
            String seen = 'from the ' + CalciumRingAndCore.SIDES[(int) cross.get('from')] + ' to the ' + CalciumRingAndCore.SIDES[(int) cross.get('to')] + ', ' + String.format(Locale.US, '%.2f', (double) cross.get('ratio')) + ' x its width'
            if (((String) out.get('rule')).isEmpty()) {
                out.put('rule', 'crosses the cell')
                out.put('text', 'possible passing object, ' + frames((int) cross.get('first'), (int) cross.get('last')))
                out.put('detail', 'the rise crossed the cell ' + seen)
            } else { String also = 'also crosses the cell ' + seen; out.put('also', also); out.put('detail', ((String) out.get('detail')) + '; ' + also) }
        }
        return out
    }
    /** Where the zones of one cell and its outline lie, for crossing: from the pixel indices of the outline and of its 25 zones
     *  (CalciumRingAndCore.zonesOf) in a picture w pixels wide.
     *  @return [centre, outline]: centre is each zone's [pixels, x, y], its pixel count and the mean x and y of its pixels minus
     *          the mean x and y of the outline's pixels, in pixels with y down (x and y 0 for a zone without pixels); outline is
     *          the outline's pixels about that centre, [x values, y values] */
    static double[][][] shape(int[] outline, int[][] zones, int w) {
        long sx = 0L, sy = 0L
        for (int i : outline) { int y = Math.floorDiv(i, w); sx += i - y * w; sy += y }
        double x0 = outline.length == 0 ? 0d : ((double) sx) / outline.length, y0 = outline.length == 0 ? 0d : ((double) sy) / outline.length
        double[][] centre = new double[zones.length][], around = new double[2][outline.length]
        for (int k = 0; k < outline.length; k++) { int y = Math.floorDiv(outline[k], w); around[0][k] = outline[k] - y * w - x0; around[1][k] = y - y0 }
        for (int z = 0; z < zones.length; z++) {
            double mx = 0d, my = 0d
            for (int i : zones[z]) { int y = Math.floorDiv(i, w); mx += i - y * w - x0; my += y - y0 }
            centre[z] = [(double) zones[z].length, zones[z].length == 0 ? 0d : mx / zones[z].length, zones[z].length == 0 ? 0d : my / zones[z].length] as double[]
        }
        double[][][] out = new double[2][][]
        out[0] = centre; out[1] = around
        return out
    }
    /** Crosses the cell: the centre of the rise over the core, the rim sectors and the near-ring sectors (CROSS_ZONES), each
     *  weighted by its pixels at its centre, moves across the outline while the cell is up.
     *  The window is the stretch around the peak where the cell stays at or above CROSS_LEVEL of its rise (event). Each zone
     *  counts with its pixels: R(f) is the sum of pixels x zone mean less all the pixels x the cell's background, X(f) and Y(f)
     *  the same with each zone's pixels at its centre, each less its mean over the F0 frames and averaged over AVERAGE frames.
     *  The start is (X, Y) / R over the first half of the window, the end over its last half (in a window of an odd number of
     *  frames the middle one is in neither). It flags when all of these hold: the travel from the start to the end is at least
     *  CROSS_TRAVEL x the outline's width along it; the centre moves that way steadily, at least CROSS_STEADY (the frames where
     *  R is above 0, each step weighted by the smaller R of its two frames: the steps forward less the steps back over all
     *  steps); the travel is at least CROSS_NOISE x its noise (from 1.4826 x the median absolute deviation of the frame-to-frame
     *  differences of X and Y before averaging, divided by root 2, carried through the averages that give the start and the
     *  end); and the outline's centre lies between the start and the end along the travel. Nothing new is measured.
     *  @param cell the cell's trace (traces), @param peak the frame of the call's peak, counted from 0
     *  @param background the cell's background, @param zone the zone means as given to traces, taken to the two decimals of
     *         the values table  @param centre, outline as shape gives them  @param firstF0, lastF0 the F0 frames, counted from 0
     *  @return null when the cell does not rise at the peak; else first and last (the window, counted from 1) and flag; once
     *          the start and end are found also start and end ([x, y], pixels from the outline's centre, y down), from and to
     *          (their sectors 0 to 7: right, lower right ...; when both are in one sector and there is a travel, the sector
     *          opposite the travel and the sector of the travel), travel (pixels); and with a travel, width (pixels across the
     *          outline along the travel), ratio (travel / width), steadiness (NaN when it cannot be worked out), noise (pixels),
     *          overNoise (travel / noise), along (where the start and the end lie along the travel, in half widths from the
     *          outline's centre: [start, end]) */
    static Map<String, Object> crossing(double[] cell, int peak, double[] background, float[][] zone, double[][] centre, double[][] outline, int firstF0, int lastF0) {
        int n = cell.length
        double rise = cell[peak]
        if (Double.isNaN(rise) || Double.isInfinite(rise) || rise <= 0d) return null
        int[] window = event(cell, peak, CROSS_LEVEL * rise)
        int a = window[0], b = window[1], frames = b - a + 1
        Map<String, Object> out = new LinkedHashMap<String, Object>()
        out.put('first', a + 1); out.put('last', b + 1); out.put('flag', false)
        if (frames < 2) return out
        // the zones summed with their pixels and centres, less the cell's background on as many pixels. The values table holds
        // two decimals and is read as float (CalciumRingAndCore.read): these sums over hundreds of pixels take the saved two
        // decimals back exactly, as a float would move the travel in its 4th digit
        double[] r = new double[n], x = new double[n], y = new double[n]
        double pixels = 0d, kx = 0d, ky = 0d
        for (int z : CROSS_ZONES) {
            if (zone == null || zone[z] == null || zone[z].length != n || !finite(zone[z]) || centre[z] == null || centre[z][0] < CalciumRingAndCore.MIN_PIXELS) continue
            double count = centre[z][0], ox = centre[z][1], oy = centre[z][2]
            for (int f = 0; f < n; f++) { double m = Math.rint(zone[z][f] * 100d) / 100d; r[f] += count * m; x[f] += count * ox * m; y[f] += count * oy * m }
            pixels += count; kx += count * ox; ky += count * oy
        }
        for (int f = 0; f < n; f++) { r[f] -= pixels * background[f]; x[f] -= kx * background[f]; y[f] -= ky * background[f] }
        double r0 = mean(r, firstF0, lastF0), x0 = mean(x, firstF0, lastF0), y0 = mean(y, firstF0, lastF0)
        double[] rs = new double[n], xs = new double[n], ys = new double[n]
        for (int f = 0; f < n; f++) { rs[f] = r[f] - r0; xs[f] = x[f] - x0; ys[f] = y[f] - y0 }
        rs = averaged(rs); xs = averaged(xs); ys = averaged(ys)
        // the start over the first half of the window, the end over its last half
        int k = frames >> 1
        double d1 = 0d, d2 = 0d, x1 = 0d, y1 = 0d, x2 = 0d, y2 = 0d
        for (int f = a; f < a + k; f++) { d1 += rs[f]; x1 += xs[f]; y1 += ys[f] }
        for (int f = b - k + 1; f <= b; f++) { d2 += rs[f]; x2 += xs[f]; y2 += ys[f] }
        if (!(d1 > 0d) || !(d2 > 0d)) return out
        double[] start = [x1 / d1, y1 / d1] as double[], end = [x2 / d2, y2 / d2] as double[]
        double travel = Math.hypot(end[0] - start[0], end[1] - start[1])
        out.put('start', start); out.put('end', end); out.put('travel', travel)
        out.put('from', CalciumRingAndCore.sector(start[1], start[0])); out.put('to', CalciumRingAndCore.sector(end[1], end[0]))
        if (!(travel > 0d) || Double.isInfinite(travel)) return out
        double ux = (end[0] - start[0]) / travel, uy = (end[1] - start[1]) / travel
        // the start and the end in one sector (an outline crossed near one end): the sides are named from the travel instead,
        // the sector opposite it and the sector it points to
        if (out.get('from') == out.get('to')) { out.put('from', CalciumRingAndCore.sector(-uy, -ux)); out.put('to', CalciumRingAndCore.sector(uy, ux)) }
        double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY
        for (int i = 0; i < outline[0].length; i++) { double p = outline[0][i] * ux + outline[1][i] * uy; if (p < lo) lo = p; if (p > hi) hi = p }
        double width = hi - lo
        // steadiness: the frames of the window where R is above 0, in order
        double forward = 0d, all = 0d, before = Double.NaN, weight = 0d
        int used = 0
        for (int f = a; f <= b; f++) {
            if (!(rs[f] > 0d)) continue
            double q = (xs[f] * ux + ys[f] * uy) / rs[f]
            if (used > 0) { double v = Math.min(weight, rs[f]), s = q - before; forward += v * s; all += v * Math.abs(s) }
            before = q; weight = rs[f]; used++
        }
        double steadiness = used >= 2 && all > 0d ? forward / all : Double.NaN
        // the noise of the travel: the noise of one frame of X and Y, carried through the averages the start and the end are made of
        double sx = spread(x), sy = spread(y)
        double[] w0 = new double[n], w1 = new double[n], w2 = new double[n]
        int h = AVERAGE >> 1
        for (int f = firstF0; f <= lastF0; f++) { int lo2 = Math.max(0, f - h), hi2 = Math.min(n, f + h + 1); for (int g = lo2; g < hi2; g++) w0[g] += 1d / (hi2 - lo2) }
        for (int g = 0; g < n; g++) w0[g] /= (lastF0 - firstF0 + 1)
        for (int f = a; f < a + k; f++) { int lo2 = Math.max(0, f - h), hi2 = Math.min(n, f + h + 1); for (int g = lo2; g < hi2; g++) w1[g] += 1d / (hi2 - lo2) }
        for (int f = b - k + 1; f <= b; f++) { int lo2 = Math.max(0, f - h), hi2 = Math.min(n, f + h + 1); for (int g = lo2; g < hi2; g++) w2[g] += 1d / (hi2 - lo2) }
        double sum = 0d
        for (int g = 0; g < n; g++) { double dw = (w2[g] - k * w0[g]) / d2 - (w1[g] - k * w0[g]) / d1; sum += dw * dw }
        double noise = Math.sqrt(sum * (sx * sx * ux * ux + sy * sy * uy * uy))
        double ratio = travel / width, overNoise = travel / noise
        double[] along = [(start[0] * ux + start[1] * uy) / (width / 2d), (end[0] * ux + end[1] * uy) / (width / 2d)] as double[]
        out.put('width', width); out.put('ratio', ratio); out.put('steadiness', steadiness); out.put('noise', noise); out.put('overNoise', overNoise); out.put('along', along)
        // NaN is tested for first: Groovy can count it as larger than any number
        out.put('flag', width > 0d && !Double.isNaN(ratio) && ratio >= CROSS_TRAVEL && !Double.isNaN(steadiness) && steadiness >= CROSS_STEADY &&
            !Double.isNaN(overNoise) && overNoise >= CROSS_NOISE && !Double.isNaN(along[0]) && along[0] <= 0d && !Double.isNaN(along[1]) && along[1] >= 0d)
        return out
    }
    /** Stretches of samples at or above the rise a cell needs, inside its response window, as the engine counts them: a sample
     *  below, a missing value or a gap in time longer than maxGap ends a stretch. @return how many, and the longest in samples */
    static int[] crossings(double[] time, double[] value, double base, double reference, double need, double start, double end, double maxGap) {
        double tol = 8d * Math.ulp(1d) * Math.max(1d, Math.max(Math.abs(start), Math.max(Math.abs(end), Math.abs(time[time.length - 1]))))
        int runs = 0, longest = 0, current = 0, last = -1
        for (int i = 0; i < time.length; i++) {
            if (time[i] < start - tol || time[i] > end + tol) continue
            double rise = (value[i] - base) / reference
            boolean above = !Double.isNaN(rise) && !Double.isInfinite(rise) && rise >= need
            if (!above) { current = 0; last = -1; continue }
            if (current > 0 && last >= 0 && !Double.isNaN(maxGap) && time[i] - time[last] > maxGap + tol) current = 0
            if (current == 0) runs++
            current++; last = i
            if (current > longest) longest = current
        }
        return [runs, longest] as int[]
    }
}

/** A reason, in plain words, why the calls of an analysis cannot be checked. */
class CalciumUnchecked extends RuntimeException { CalciumUnchecked(String message){super(message)} }

/** Checks on the calls of a finished analysis, saved beside its percent table as <recording>_call_checks_<analysis>.csv, one
 *  row per cell. A check only flags: no call, count or percentage changes.
 *   - "possible passing object, frames a to b" and "rise on the <side> edge, not in the cell, frames a to b"
 *     (CalciumPassingObject), from what was measured around each cell (CalciumRingAndCore). When the field moved near the
 *     event (CalciumMotion.movedNear, from the span before the event to the span after it: MOVED pixels or more from where it
 *     was over the F0 frames, or between two of those frames) the row says "field moved, frames a to b" instead: a field that
 *     moves gives the same picture as an object.
 *   - "possible passing object, frames a to b" of the rule "crosses the cell" (CalciumPassingObject.crossing), with a and b
 *     its own window, and the field moved near that window as above. When rule A or B already flags the cell, that flag
 *     stays and its detail says the rise also crosses the cell. Where the zones lie is rebuilt from the outlines (shapes);
 *     without them, or when they do not give the saved counts, the detail says the crossing was not checked.
 *   - "short separate crossings": a responder whose rise reached what it needs CROSSINGS or more separate times in the
 *     response window and never for more than SHORT samples, counted from the engine's own columns.
 *  A stretch of one frame reads "frame a". Only the event of the largest rise is looked at. The F0 frames are the analysis's
 *  own baseline frames.
 *  It is plugin code in Groovy and not in the R adapter because it needs the outlines' pixels when the ring and core could not
 *  be saved at Measure, and because its table lies beside the percent table, which Groovy writes; the analysis folder is a
 *  finished snapshot and is not written to. */
class CalciumCallChecks {
    /** File names after "<recording>_". */
    static final String VALUES='ring_and_core_gray_values.csv.gz',PIXELS='ring_and_core_pixels.csv',MOVEMENT='field_movement.csv'
    static final int CROSSINGS=3,SHORT=2
    /** The columns of a table written by 2026.11.24.1. It still reads: the codes then come from its Rule column (codesOf). */
    static final List HEADER_BEFORE=['Cell','Call','Flag','Rule','Event_first_frame','Event_last_frame','Detail','Near_ring_sectors_usable','Area_pixels','Ring_measured','Field_movement_measured','Passing_object_check','Not_checked_because','Crossings','Longest_crossing_samples']
    /** QC_flags is last, so every other column keeps its place. */
    static final List HEADER=HEADER_BEFORE+['QC_flags']
    /** The code of each rule's flag, defined here only: the QC_flags column of the table holds them, as the engine's own
     *  QC_flags holds its codes. */
    static final Map CODES=['passes through':'POSSIBLE_PASSING_OBJECT','on the edge':'RISE_ON_EDGE','field moved':'FIELD_MOVED','short separate crossings':'SHORT_SEPARATE_CROSSINGS','crosses the cell':'POSSIBLE_PASSING_OBJECT']
    /** Each code as a table on screen names it: the names the filters of the flagged list use. */
    static final Map LABELS=[POSSIBLE_PASSING_OBJECT:'Possible passing object',RISE_ON_EDGE:'Rise on the edge',FIELD_MOVED:'Field moved',SHORT_SEPARATE_CROSSINGS:'Short separate crossings']
    /** Codes are joined as the engine joins its own (responder_qc_flags in analysis/upstream/responders.R). */
    static final String JOIN='; '
    /** The codes of a row's flags from its Rule column, in the flags' order: '' without a flag. */
    static String codesOf(Object rule){String text=rule==null?'':''+rule;text.isEmpty()?'':(text.split('; ') as List).collect{CODES[it]?:''}.findAll{it}.join(JOIN)}
    /** Why nothing is shown around a cell of the control window. */
    static final String AROUND_CONTROL='Checks are made for the primary stimulus.'
    /** Why nothing is shown around the cells of an analysis without a checks table. */
    static final String AROUND_BEFORE='This analysis was made before the calls were checked.'
    /** Why nothing is shown around a cell when what the check read has changed since. */
    static final String AROUND_CHANGED='Not shown: the tables beside the analysis have changed since the check.'
    /** Why one cell was not checked when the reason is its own and not the recording's (Not_checked_because). */
    static final List CELL_REASONS=['no call','no peak time','missing gray values','no rise at the peak','no ring around the cell']
    /** One row per cell: its outline, its pixels and the pixels of each zone, then what the row belongs to (the raw table and
     *  the ring and core values table of the same measurement) and the widths the zones were built with. */
    static final List PIXEL_HEADER=['Cell','ROI_sha256','Area_pixels']+CalciumRingAndCore.NAMES+['Near_ring_sectors_usable','Raw_table_sha256','Ring_and_core_values_sha256','Pixel_size_um','Near_ring_from_px','Near_ring_to_px','Far_ring_to_px','Other_outline_margin_px']
    /** Held while Step 4 measures around the cells, so two analyses at once do not measure the same recording together. */
    static final Object MEASURING=new Object()
    /** The outline as the ROI registry identifies it (ROI_sha256). */
    static String outlineHash(Roi roi,String id){CalciumMeasurement.outlineHash(roi,id)}
    /** The line for what_was_done.txt. @param base the files' names up to and with "<recording>_"
     *  @param movement whether the field's displacement was saved too (it is not when the movement check of Step 1 could not
     *         measure it) */
    static String record(CalciumRingAndCore rings,String base,boolean movement=true) {
        def n={double v->CalciumTheme.shortNumber(v)}
        double[] w=rings.width
        String where=rings.micrometres?'a near ring '+n(CalciumRingAndCore.NEAR_FROM_UM)+' to '+n(CalciumRingAndCore.NEAR_TO_UM)+' µm and a far ring '+n(CalciumRingAndCore.NEAR_TO_UM)+' to '+n(CalciumRingAndCore.FAR_TO_UM)+' µm outside the outline, here '+n(w[0])+' to '+n(w[1])+' and '+n(w[1])+' to '+n(w[2])+' pixels. The rings leave out pixels within '+n(CalciumRingAndCore.OTHER_UM)+' µm of another outline.':
            'a near ring '+n(w[0])+' to '+n(w[1])+' pixels and a far ring '+n(w[1])+' to '+n(w[2])+' pixels outside the outline. '+(rings.scale=='not square'?'The pixels of the recording are not square':'The recording has no pixel size in micrometres')+', so the widths are the pixel values as they are. The rings leave out pixels within '+n(w[3])+' pixels of another outline.'
        'Around each cell: mean gray value on every frame of '+where+' Each ring is in 8 sectors. Inside the outline: the core, which is the deeper half, and the rim in 8 sectors. A zone under '+CalciumRingAndCore.MIN_PIXELS+' pixels is left empty. Files: '+(movement?base+VALUES+', '+base+PIXELS+' and '+base+MOVEMENT+', which holds the displacement of the field at every frame.':base+VALUES+' and '+base+PIXELS+'. The movement of the field was not measured.')+' Step 4 uses them to flag a possible passing object. No measurement and no call changes.'
    }
    /** Beside the tables, once the values table and the raw table are written: the pixels of each zone, with the outline, the
     *  raw table and the values table each row belongs to. */
    static void savePixels(String prefix,List names,Roi[] rois,int[][] cellPx,CalciumRingAndCore rings) {
        String raw=CalciumInstaller.digest(new File(prefix+'gray_values_raw.csv').bytes),values=CalciumInstaller.digest(new File(prefix+VALUES).bytes)
        List rows=[PIXEL_HEADER]
        names.eachWithIndex{name,int c->rows<<([name,outlineHash(rois[c],''+name),cellPx[c].length]+(0..<CalciumRingAndCore.ZONES).collect{int z->rings.pixels(c,z)}+
            [rings.nearUsable(c),raw,values,rings.micrometres?rings.pixelSize:'',rings.width[0],rings.width[1],rings.width[2],rings.width[3]])}
        CalciumCsv.write(new File(prefix+PIXELS),rows)
    }
    /** Takes away what an earlier measurement left when this one could not save around the cells. */
    static void removeBeside(String prefix){[VALUES,VALUES+'.part',PIXELS,MOVEMENT].each{String name->try{new File(prefix+name).delete()}catch(Throwable ignored){}}}
    /** The field's displacement at every frame, counted from where the field sits for most of the recording, with the shifts it
     *  comes from. Without shifts the file is taken away: the movement was not measured. */
    static void saveMovement(File file,double[][] shifts) {
        if(shifts==null||shifts.length==0){file.delete();return}
        double[] size=CalciumMotion.displacement(shifts,0d,0d)
        def text=new StringBuilder('Frame,Displacement_px,Shift_x_px,Shift_y_px\n')
        for(int f=0;f<shifts.length;f++)text.append(f+1).append(',').append(String.format(Locale.US,'%.2f,%.2f,%.2f',size[f],shifts[f][0],shifts[f][1])).append('\n')
        file.setText(text.toString(),'UTF-8')
    }
    /** The saved files, when they belong to these cells, these outlines and this raw table.
     *  @return values, counts, area, near and size (the saved pixel size in micrometres, '' without one), or why not */
    static Map savedAround(String prefix,List names,Roi[] rois,String rawHash,int frames) {
        File pixels=new File(prefix+PIXELS),values=new File(prefix+VALUES)
        if(!pixels.isFile()||!values.isFile())return [why:'ring and core were not measured']
        try {
            def rows=CalciumCsv.read(pixels)
            if(rows.isEmpty()||rows[0]!=PIXEL_HEADER)return [why:'the ring and core files have another layout']
            List body=rows.drop(1)
            if(body.collect{it[0]}!=names)return [why:'the ring and core files are for other cells']
            if(rois!=null&&(0..<names.size()).any{int c->body[c][1]!=outlineHash(rois[c],names[c] as String)})return [why:'the ring and core files are for other outlines']
            int rawAt=PIXEL_HEADER.indexOf('Raw_table_sha256'),valuesAt=PIXEL_HEADER.indexOf('Ring_and_core_values_sha256'),nearAt=PIXEL_HEADER.indexOf('Near_ring_sectors_usable')
            if(body.any{it[rawAt]!=rawHash})return [why:'the ring and core files are for another raw table']
            String valuesHash=CalciumInstaller.digest(values.bytes)
            if(body.any{it[valuesAt]!=valuesHash})return [why:'the ring and core values table has changed since it was measured']
            return [values:CalciumRingAndCore.read(values,names.collect{''+it},frames),counts:body.collect{List row->(3..<3+CalciumRingAndCore.ZONES).collect{row[it] as int}},area:body.collect{it[2] as int},near:body.collect{it[nearAt] as int},
                    size:body?''+body[0][PIXEL_HEADER.indexOf('Pixel_size_um')]:'']
        }catch(Throwable problem){return [why:'the ring and core files cannot be read']}
    }
    /** Whether a recording gives the saved raw gray values: every cell on the first, the middle and the last frame. */
    static boolean sameRecording(ImageStack stack,int[][] cellPx,List raw,int frames) {
        for(int f in [1,(frames+1).intdiv(2),frames].unique()) {
            def ip=stack.getProcessor(f)
            for(int c=0;c<cellPx.length;c++){if(cellPx[c].length==0)continue
                double sum=0d;for(int i:cellPx[c])sum+=ip.getf(i)
                double off=Math.abs(sum/cellPx[c].length-CalciumCsv.number(''+raw[f][c+1]))
                if(Double.isNaN(off)||off>0.00051d)return false}   // NaN first: Groovy counts it as larger than any number
        }
        true
    }
    /** A problem the checks could not get past, for whoever looks into it: in the analysis log when there is one, else in Fiji's Log. */
    static void logProblem(File log,String what,Throwable problem) {
        String line='Checks on the calls, '+what+': '+problem
        try{if(log!=null){log.parentFile?.mkdirs();log.withWriterAppend('UTF-8'){it.write(line+'\n')};return}}catch(Throwable ignored){}
        try{IJ.log(CalciumRelease.NAME+': '+line)}catch(Throwable ignored){}
    }
    /** What was measured around the cells of this measurement: read when the files beside the tables belong to it; otherwise
     *  measured now with the same code when the recording and the outlines are at hand, and saved beside the tables (new files
     *  only). The field's movement is not measured again for files that belong to the measurement: when it is missing, the
     *  record says it was not measured, and the checks say a moving field is not told from an object.
     *  @param measure false to read the saved files only, as what is shown around a cell does
     *  @return values, counts, area and near, or why not */
    static Map around(Map measurement,int frames,File log=null,boolean measure=true) {
        if(!(measurement?.folder instanceof File))return [why:'this measurement has no folder']
        String prefix=new File(measurement.folder as File,measurement.recording+'_').path
        List names=((measurement.names?:[]) as List).collect{''+it}
        File raw=new File(prefix+'gray_values_raw.csv')
        if(!raw.isFile())return [why:'the raw gray values are not beside the analysis']
        Roi[] rois=measurement.rois instanceof Roi[]&&(measurement.rois as Roi[]).length==names.size()?measurement.rois as Roi[]:null
        Map saved=savedAround(prefix,names,rois,CalciumInstaller.digest(raw.bytes),frames)
        if(saved.values!=null||!measure)return saved
        if(rois==null)return [why:saved.why+' and the outlines are not at hand']
        // the recording as every frame reader takes it: this one while it holds its frames, else the same recording opened again
        ImagePlus movie=CalciumSteps.live(measurement.movie instanceof ImagePlus?measurement.movie as ImagePlus:null,measurement.recordingPrint instanceof String?measurement.recordingPrint as String:'')
        if(movie==null||movie.getStackSize()<frames)return [why:saved.why+' and the recording is not open']
        synchronized(MEASURING) {
            Map meanwhile=savedAround(prefix,names,rois,CalciumInstaller.digest(raw.bytes),frames);if(meanwhile.values!=null)return meanwhile
            try {
                int w=movie.getWidth(),h=movie.getHeight()
                int[][] cellPx=rois.collect{Roi roi->CalciumSteps.pixelsOf(roi,w,h)} as int[][]
                def table=CalciumCsv.read(raw)
                if(table.size()!=frames+1||table[0].drop(1).collect{''+it}!=names)return [why:saved.why+' and the raw gray values are for other cells or frames']
                if(!sameRecording(movie.getStack(),cellPx,table,frames))return [why:saved.why+' and the open recording does not give the saved gray values']
                def scale=movie.getCalibration()
                def rings=new CalciumRingAndCore(cellPx,w,h,scale.pixelWidth,scale.pixelHeight,scale.getUnit())
                rings.measure(movie.getStack(),frames,new File(prefix+VALUES),names as String[])
                double[][] shifts=null
                if(CalciumMotion.measurable(w,h,frames))try{shifts=CalciumMotion.shifts(movie.getStack(),frames)}catch(Throwable problem){logProblem(log,'the movement of the field',problem)}
                saveMovement(new File(prefix+MOVEMENT),shifts)
                savePixels(prefix,names,rois,cellPx,rings)
                return savedAround(prefix,names,rois,CalciumInstaller.digest(raw.bytes),frames)
            }catch(Throwable problem){
                logProblem(log,'around each cell',problem)
                [VALUES,PIXELS].each{String name->try{new File(prefix+name).delete()}catch(Throwable ignored){}}   // nothing half saved stays
                return [why:saved.why+(problem instanceof IOException?' and could not be saved now':' and could not be measured now')]
            }
        }
    }
    /** The saved shifts of the field at every frame, [frame][x, y] in pixels; null when they were not saved for these frames. */
    static double[][] movement(String prefix,int frames) {
        try {
            File file=new File(prefix+MOVEMENT);if(!file.isFile())return null
            def rows=CalciumCsv.read(file);int x=rows[0].indexOf('Shift_x_px'),y=rows[0].indexOf('Shift_y_px')
            if(x<0||y<0||rows.size()!=frames+1)return null
            double[][] out=new double[frames][2]
            for(int f=0;f<frames;f++){out[f][0]=CalciumCsv.number(''+rows[f+1][x]);out[f][1]=CalciumCsv.number(''+rows[f+1][y]);if(!Double.isFinite(out[f][0])||!Double.isFinite(out[f][1]))return null}
            return out
        }catch(Throwable ignored){return null}
    }
    /** Where the zones of every cell lie, for the rule "crosses the cell": rebuilt from the outlines with the zone code of Measure
     *  (CalciumRingAndCore) and the saved pixel size, on a picture of the recording's size, and used for a cell only when every
     *  rebuilt count is the saved one. Nothing is measured. @param ring what around gave
     *  @return cells, one per cell of names: centre and outline as CalciumPassingObject.shape gives them, or why not; or why,
     *          for every cell */
    static Map shapes(Map measurement,List names,Map ring) {
        Roi[] rois=measurement?.rois instanceof Roi[]&&(measurement.rois as Roi[]).length==names.size()?measurement.rois as Roi[]:null
        ImagePlus picture=[measurement?.picture,measurement?.movie].find{it instanceof ImagePlus&&(it as ImagePlus).getWidth()>0&&(it as ImagePlus).getHeight()>0} as ImagePlus
        if(rois==null||picture==null||ring?.counts==null)return [why:'the outlines are not at hand']
        try {
            int w=picture.getWidth(),h=picture.getHeight()
            int[][] cellPx=rois.collect{Roi roi->CalciumSteps.pixelsOf(roi,w,h)} as int[][]
            double size=ring.size?CalciumCsv.number(''+ring.size):Double.NaN
            def rebuilt=new CalciumRingAndCore(cellPx,w,h,size,size,'micron')
            List cells=(0..<names.size()).collect{int c->
                List counts=ring.counts[c] as List
                if(cellPx[c].length!=(ring.area[c] as int)||(0..<CalciumRingAndCore.ZONES).any{int z->rebuilt.pixels(c,z)!=(counts[z] as int)})return [why:'the outline does not match its ring file']
                double[][][] at=CalciumPassingObject.shape(cellPx[c],(0..<CalciumRingAndCore.ZONES).collect{int z->rebuilt.zonePx[c*CalciumRingAndCore.ZONES+z]} as int[][],w)
                [centre:at[0],outline:at[1]]
            }
            return [cells:cells]
        }catch(Throwable problem){return [why:'the outlines cannot be placed on the recording']}
    }
    /** "frames 32 to 35 and 40", for stretches [a, b] counted from 1. */
    static String runsText(List runs){
        List parts=runs.collect{List r->r[0]==r[1]?''+r[0]:r[0]+' to '+r[1]}
        (runs.size()==1&&runs[0][0]==runs[0][1]?'frame ':'frames ')+(parts.size()==1?parts[0]:parts.dropRight(1).join(', ')+' and '+parts[-1])
    }
    /** What CalciumMotion.movedNear found, for the detail of a "field moved" row. */
    static String movedText(Map moved) {
        String size=String.format(Locale.US,'%.1f',moved.largest as double)
        if(moved.from=='between')return 'the field moved '+size+' px between frames '+moved.first+' and '+moved.last
        'the field was up to '+size+' px from '+(moved.from=='F0'?'where it was over the F0 frames':'where it sits for most of the recording')+', at '+runsText(moved.runs as List)
    }
    /** The measurement's raw and background tables by cell, when they are the ones the analysis was run on. */
    static Map tables(Map measurement,List names,int frames,Map manifest) {
        if(!(measurement?.folder instanceof File))return [why:'this measurement has no folder']
        File folder=measurement.folder as File;String base=''+measurement.recording
        File raw=new File(folder,base+'_gray_values_raw.csv'),background=new File(folder,base+'_background_next_to_each_cell.csv'),corrected=new File(folder,base+'_gray_values_background_corrected_native.csv')
        if(!raw.isFile()||!background.isFile())return [why:'the raw gray values and backgrounds are not beside the analysis']
        String analysed=manifest.measurement_input_sha256?:manifest.expected_input_sha256?:''
        if(analysed&&corrected.isFile()&&CalciumInstaller.digest(corrected.bytes)!=analysed)return [why:'the tables beside the analysis are from another measurement']
        Map out=[:]
        for(entry in [raw:raw,background:background]) {
            def rows=CalciumCsv.read(entry.value)
            if(rows.size()!=frames+1||rows[0].drop(1).collect{''+it}!=names)return [why:'the tables beside the analysis are for other cells or frames']
            Map columns=[:];List body=rows.drop(1)
            names.eachWithIndex{name,int c->double[] v=new double[frames];for(int f=0;f<frames;f++)v[f]=CalciumCsv.number(''+body[f][c+1]);columns[name]=v}
            out[entry.key]=columns
        }
        out
    }
    static File file(Map measurement,File dir){measurement?.folder instanceof File?new File(measurement.folder as File,measurement.recording+'_call_checks_'+dir.name+'.csv'):null}
    /** Checks the calls of a finished analysis and saves them beside its percent table. Never stops an analysis: when the
     *  checks cannot be made the answer says why in plain words, and the details go to the log (logProblem).
     *  @param log the analysis log, or null for Fiji's Log  @return rows, file, ringMeasured, note */
    static Map run(File dir,Map measurement,File log=null) {
        Map found
        try{found=rows(dir,measurement,log)}
        catch(CalciumUnchecked plain){return [rows:[],ringMeasured:false,note:'The calls were not checked: '+plain.message+'.']}
        catch(Throwable problem){logProblem(log,'the saved analysis',problem);return [rows:[],ringMeasured:false,note:'The calls were not checked: the saved analysis could not be read.']}
        File target=file(measurement,dir)
        if(target!=null)try{CalciumCsv.write(target,[HEADER]+(found.rows as List).collect{Map row->HEADER.collect{row[it]}});found.file=target.name}
            catch(Throwable problem){logProblem(log,'saving '+target.name,problem);found.note=((found.note?:'')+' The checks were not saved: '+target.name+' could not be written.').trim()}
        found
    }
    static Map rows(File dir,Map measurement,File log=null) {
        Map context=context(dir,measurement,log)
        List out=[];(context.names as List).eachWithIndex{String name,int c->out<<cellRow(context,c,name)}
        Map ring=context.ring as Map
        [rows:out,ringMeasured:ring.values!=null,note:note(ring.values!=null,context.shifts!=null,context.unchecked as String)]
    }
    /** A number of a row of the calls table, by its column's name: NaN without that column or row. */
    static double number(List header,List row,String name){int i=header.indexOf(name);i<0||row==null?Double.NaN:CalciumCsv.number(''+row[i])}
    /** What the checks of an analysis are made from, read once for all its cells: its calls, times and F0 frames, the raw and
     *  background tables, and what was measured around the cells. Throws CalciumUnchecked with the reason in plain words.
     *  @param measure false to read the files beside the tables only, never measure (around), as what is shown around a cell does
     *  @param before true to make the checks as 2026.11.24.1 made them, for an analysis whose saved table is one of its
     *         (HEADER_BEFORE): without the rule "crosses the cell", and the sectors of "passes through" chosen as it chose them */
    static Map context(File dir,Map measurement,File log=null,boolean measure=true,boolean before=false) {
        def calls=CalciumCsv.read(new File(dir,'responder_calls.csv'));List h=calls[0]
        int cellAt=h.indexOf('Cell'),roleAt=h.indexOf('Role'),statusAt=h.indexOf('Status')
        if(cellAt<0||statusAt<0)throw new CalciumUnchecked('the calls table has no Cell or Status column')
        Map byCell=[:];calls.drop(1).each{List row->if((roleAt<0||row[roleAt]=='Stimulus')&&!byCell.containsKey(''+row[cellAt]))byCell[''+row[cellAt]]=row}
        def normalized=CalciumCsv.read(new File(dir,'normalized_traces.csv'));List samples=normalized.drop(1)
        double[] time=samples.collect{CalciumCsv.number(''+it[0])} as double[];int frames=time.length
        if(frames<3)throw new CalciumUnchecked('the analysis has fewer than 3 frames')
        if(time.any{!Double.isFinite(it)}||(1..<frames).any{int f->!(time[f]>time[f-1])})throw new CalciumUnchecked('the saved traces have no usable times')
        Map traceAt=[:];normalized[0].drop(1).eachWithIndex{name,int i->traceAt[''+name]=i+1}
        Map manifest=CalciumCsv.read(new File(dir,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
        double from=CalciumCsv.number(''+manifest.baseline_start),to=CalciumCsv.number(''+manifest.baseline_end)
        if(!Double.isFinite(from)||!Double.isFinite(to)||from<1||to<from||to>frames)throw new CalciumUnchecked('the saved baseline frames do not match the recording')
        int first=(from as int)-1,last=(to as int)-1
        List names=measurement?.names?(measurement.names as List).collect{''+it}:(byCell.keySet() as List)
        double[] steps=new double[frames-1];for(int f=1;f<frames;f++)steps[f-1]=time[f]-time[f-1]
        int span=Math.max(1,Math.round(CalciumPassingObject.SPAN_S/CalciumMotion.median(steps)) as int)
        Map measured=tables(measurement,names,frames,manifest)
        Map ring=measured.why ? [why:measured.why] : around(measurement,frames,log,measure)
        double[][] shifts=ring.values!=null?movement(new File(measurement.folder as File,measurement.recording+'_').path,frames):null
        String unchecked=ring.values!=null?'':(''+(ring.why?:'ring and core were not measured'))   // why no cell is checked for a passing object
        Map shapes=ring.values!=null&&!before?shapes(measurement,names,ring):null   // where the zones lie, for the rule "crosses the cell"
        [header:h,byCell:byCell,statusAt:statusAt,samples:samples,time:time,frames:frames,traceAt:traceAt,first:first,last:last,names:names,span:span,measured:measured,ring:ring,shifts:shifts,unchecked:unchecked,shapes:shapes,before:before]
    }
    /** The row of one cell of the table, from what context read. @param c the cell's place in context.names
     *  @param keep when given, what the row was made from is put in it: event (CalciumPassingObject.check, with its traces),
     *         zone (the zone traces given to it; null when nothing was measured around the cell), moved (CalciumMotion.movedNear) */
    static Map cellRow(Map context,int c,String name,Map keep=null) {
        List h=context.header as List;Map byCell=context.byCell as Map,measured=context.measured as Map,ring=context.ring as Map,traceAt=context.traceAt as Map
        int statusAt=context.statusAt as int,frames=context.frames as int,first=context.first as int,last=context.last as int,span=context.span as int
        double[] time=context.time as double[];double[][] shifts=context.shifts as double[][];String unchecked=context.unchecked as String;List samples=context.samples as List
        List call=byCell[name] as List;String status=call==null?'':''+call[statusAt]
        boolean scored=status=='Responded'||status=='Did not respond'
        Map row=[Cell:name,Call:status,Flag:'',Rule:'',Event_first_frame:'',Event_last_frame:'',Detail:'',Near_ring_sectors_usable:'',Area_pixels:'',Ring_measured:ring.values!=null?'yes':'no',
                 Field_movement_measured:shifts!=null?'yes':'no',Passing_object_check:'not checked',Not_checked_because:'',Crossings:'',Longest_crossing_samples:'',QC_flags:'']
        if(ring.values!=null){row.Near_ring_sectors_usable=ring.near[c];row.Area_pixels=ring.area[c]}
        List flags=[]   // a flag is [text, rule, detail]
        String why='',crossNote=''   // crossNote: why the rule "crosses the cell" was not checked on a cell that was checked
        Map event=null,moved=null;float[][] zone=null
        double peakAt=number(h,call,'Peak_s')
        if(!scored)why='no call'
        else if(!Double.isFinite(peakAt))why='no peak time'
        else if(measured.why)why=unchecked
        else {
            double[] raw=measured.raw[name] as double[],background=measured.background[name] as double[]
            if(!CalciumPassingObject.finite(raw)||!CalciumPassingObject.finite(background))why='missing gray values'
            else {
                int peak=0;for(int f=1;f<frames;f++)if(Math.abs(time[f]-peakAt)<Math.abs(time[peak]-peakAt))peak=f
                int outside=0
                if(ring.values!=null){zone=new float[CalciumRingAndCore.ZONES][];for(int z=0;z<CalciumRingAndCore.ZONES;z++)if((ring.counts[c][z] as int)>=CalciumRingAndCore.MIN_PIXELS){zone[z]=(ring.values as float[][])[c*CalciumRingAndCore.ZONES+z];if(z<16)outside++}}
                Map shape=context.shapes==null?null:(context.shapes as Map).why?context.shapes as Map:((context.shapes as Map).cells as List)[c] as Map
                event=CalciumPassingObject.check(raw,background,zone,first,last,peak,span,shape?.centre as double[][],shape?.outline as double[][],context.before==true)
                if(event!=null){row.Event_first_frame=event.first;row.Event_last_frame=event.last}
                if(zone==null)why=unchecked
                else if(event==null)why='no rise at the peak'
                else if(outside==0)why='no ring around the cell'   // no ring sector has enough pixels: neither rule can be tested
                else {
                    row.Passing_object_check='checked'
                    if(shape?.why)crossNote='crossing not checked: '+shape.why
                    if(event.rule) {
                        // the field moved near the frames the flag names: the event, or the window of the rule "crosses the cell"
                        Map frameOf=event.rule=='crosses the cell'?event.cross as Map:event
                        moved=shifts==null?null:CalciumMotion.movedNear(shifts,first+1,last+1,(frameOf.first as int)-span,(frameOf.last as int)+span)
                        // rule A or B flags the cell and the rise crosses it too: that is said only when the field did not move
                        // near the window of the rule "crosses the cell", as for that rule alone
                        String own=event.detail as String,also=event.also as String
                        Map cross=event.cross as Map
                        if(also&&shifts!=null&&CalciumMotion.movedNear(shifts,first+1,last+1,(cross.first as int)-span,(cross.last as int)+span)!=null&&own.endsWith('; '+also))own=own.substring(0,own.length()-also.length()-2)
                        String detail=own+(crossNote?'; '+crossNote:'')
                        if(moved!=null)flags<<['field moved, '+CalciumPassingObject.frames(moved.first as int,moved.last as int),'field moved',movedText(moved)+'; otherwise '+event.text+': '+detail]
                        else flags<<[event.text,event.rule,detail]
                        crossNote=''
                    }
                }
            }
        }
        // separate crossings of the rise the cell needs, from the engine's own columns
        if(scored) {
            double need=number(h,call,'Threshold_increase'),base=number(h,call,'Baseline_mean'),reference=number(h,call,'Local_reference_relative')
            if(!Double.isFinite(reference))reference=number(h,call,'Local_F0_relative')
            if(!Double.isFinite(reference))reference=1d+base
            double start=number(h,call,'Start_s'),end=number(h,call,'Effective_end_s');if(!Double.isFinite(end))end=number(h,call,'End_s')
            if(traceAt[name]!=null&&[need,base,reference,start,end].every{Double.isFinite(it as double)}&&need>0d&&reference>0d) {
                int at=traceAt[name] as int;double[] y=samples.collect{CalciumCsv.number(''+it[at])} as double[]
                int[] crossed=CalciumPassingObject.crossings(time,y,base,reference,need,start,end,number(h,call,'Max_gap_s'))
                row.Crossings=crossed[0];row.Longest_crossing_samples=crossed[1]
                if(status=='Responded'&&crossed[0]>=CROSSINGS&&crossed[1]<=SHORT)flags<<['short separate crossings','short separate crossings',crossed[0]+' crossings, the longest '+crossed[1]+' samples']
            }
        }
        row.Flag=flags.collect{it[0]}.join('; ');row.Rule=flags.collect{it[1]}.join('; ');row.Detail=((crossNote ? [crossNote] : [])+flags.collect{it[2]}).join('; ');row.Not_checked_because=why
        row.QC_flags=codesOf(row.Rule)
        if(keep!=null){keep.event=event;keep.zone=zone;keep.moved=moved}
        row
    }
    /** The place of a zone (CalciumRingAndCore.NAMES) as the legend names it: "ring, lower right", "far ring, upper left",
     *  "rim, lower right", "core". */
    static String place(int zone){zone==CalciumRingAndCore.CORE?'core':(zone<8?'ring, ':zone<16?'far ring, ':'rim, ')+CalciumRingAndCore.SIDES[zone<16?zone%8:zone-CalciumRingAndCore.CORE-1]}
    /** What is shown around the cells of an analysis is read from: what context reads for the checks, from the files beside the
     *  tables only (nothing is measured). Read once per analysis; aroundCell takes one cell from it. For an analysis whose saved
     *  table is one of 2026.11.24.1 the checks are made again as that version made them (context), so they give its rows.
     *  @return the context, or why: the sentence that says why nothing can be shown for any cell */
    static Map aroundContext(File dir,Map measurement) {
        try{return context(dir,measurement,null,false,savedBefore(dir,measurement))}
        catch(CalciumUnchecked plain){return [why:'The calls were not checked: '+plain.message+'.']}
        catch(Throwable problem){return [why:'The calls were not checked: the saved analysis could not be read.']}
    }
    /** What is shown around one cell: what the check of the calls used, from the same code (cellRow and
     *  CalciumPassingObject.check), not a second calculation. It is given only when it gives the cell's row of the checks as
     *  they were made; otherwise why says what changed.
     *  @param context what aroundContext gave for the analysis  @param checks the analysis's checks as answer gives them, or null
     *         when the analysis has none
     *  @return cell, name (as the page shows it) and why: '' when there is something to show, else the sentence that says why
     *          not. With something to show: time (s, per frame), rise (the cell's rise per frame in gray levels: raw minus
     *          background minus its F0 mean, averaged over 3 frames), zones (the 25 zones in CalciumRingAndCore.NAMES order:
     *          zone, key, place, usable, trace in gray levels, null when not usable), eventFirst and eventLast (frames, counted
     *          from 1), rule, text, flags and codes (as answer gives them), lines (the zones the flag rests on, each with zone,
     *          place, part: before, after, rim, ring or core, and frame where the sector was seen, else null; for rule C the
     *          near-ring sectors on the side the rise came from and went to, and the core; for a cell without a flag of rule
     *          A, B or C, the core), others (every other usable near-ring sector), band (the frames the flag's text names,
     *          [first, last], as band gives them; null without one) */
    static Map aroundCell(Map context,Map checks,String cell) {
        List names=(context?.names?:[]) as List
        Map out=[cell:cell,name:CalciumReview.shown(names,cell),why:'']
        def refuse={String why->out.why=why;out}
        // no checks saved: made before the checks, or the run could not check the calls, and then the analysis says why
        if(checks==null)return refuse((context?.why?:AROUND_BEFORE) as String)
        // nothing in this analysis was checked for a passing object: the checks' own note says why, for every cell
        if(checks.ringMeasured!=true)return refuse((checks.note?:'Not checked for a passing object.') as String)
        Map saved=(checks.cells as Map)?.get(cell) as Map
        if(saved==null)return refuse('There is no cell '+cell+' in this analysis.')
        if(!saved.checked)return refuse('Passing object not checked: '+saved.why+'.')
        if(context.why)return refuse(context.why as String)
        int c=names.indexOf(cell)
        if(c<0)return refuse('There is no cell '+cell+' in this analysis.')
        Map ring=context.ring as Map
        if(ring.values==null)return refuse('Not shown: '+context.unchecked+'.')
        Map keep=[:],row=cellRow(context,c,cell,keep),now=(answer([rows:[row]]).cells as Map)[cell] as Map
        Map event=keep.event as Map
        if(now!=saved||event==null||keep.zone==null)return refuse(AROUND_CHANGED)
        Map traces=event.traces as Map;double[][] trace=traces.trace as double[][]
        out.time=context.time;out.rise=traces.cell
        out.zones=(0..<CalciumRingAndCore.ZONES).collect{int z->[zone:z,key:CalciumRingAndCore.NAMES[z],place:place(z),usable:trace[z]!=null,trace:trace[z]]}
        out.eventFirst=event.first;out.eventLast=event.last
        out.rule=now.rule;out.text=now.text;out.flags=now.flags;out.codes=now.codes
        List rules=(now.flags as List).collect{it.rule},lines=[]
        def line={int z,String part,Object frame->if(trace[z]!=null)lines<<[zone:z,place:place(z),part:part,frame:frame]}
        // a field that moved replaces rule A, B or C: the zones are the ones that rule rests on; rule C's only when A and B do
        // not flag the cell
        boolean flagged=rules.any{it in ['passes through','on the edge','field moved','crosses the cell']}
        Map cross=event.cross as Map
        if(flagged&&event.before!=null){List a=event.before as List,b=event.after as List;line(a[0] as int,'before',a[1]);line(b[0] as int,'after',b[1])}
        else if(flagged&&event.edge!=null){int hot=(event.edge as List)[0] as int;line(CalciumRingAndCore.CORE+1+hot,'rim',null);line(hot,'ring',null);line(CalciumRingAndCore.CORE,'core',null)}
        else if(flagged&&cross?.flag==true){line(cross.from as int,'before',null);if(cross.to!=cross.from)line(cross.to as int,'after',null);line(CalciumRingAndCore.CORE,'core',null)}
        else line(CalciumRingAndCore.CORE,'core',null)
        out.lines=lines
        Set named=lines.collect{it.zone} as Set
        out.others=(0..<8).findAll{int z->trace[z]!=null&&!named.contains(z)}
        out.band=band(event,rules,keep.moved as Map)
        out
    }
    /** The frames shown as "flagged frames" around a cell: the frames its flag's text names. Rule A (passes through): the
     *  event; rule B (on the edge): the stretch on the edge; rule C (crosses the cell): its window; field moved: the moved
     *  frames. Short separate crossings alone names no frames: no band.
     *  @param event what CalciumPassingObject.check gave  @param rules the rules of the cell's flags
     *  @param moved what CalciumMotion.movedNear gave for the flag (cellRow)
     *  @return [first, last], counted from 1; null without a flag that names frames */
    static List band(Map event,List rules,Map moved=null){
        String rule=(rules?:[]).find{it!='short separate crossings'}
        if(event==null||rule==null)return null
        if(rule=='field moved')return moved==null?null:[moved.first,moved.last]
        if(rule=='on the edge'){List edge=event.edge as List;return edge==null?null:[edge[1],edge[2]]}
        if(rule=='crosses the cell'){Map cross=event.cross as Map;return cross==null?null:[cross.first,cross.last]}
        [event.first,event.last]
    }
    /** The note for the whole analysis: '' when every evaluable cell could be checked as written. */
    static String note(boolean ring,boolean movement,String unchecked){!ring?'Not checked for a passing object'+(unchecked?': '+unchecked:'')+'.':movement?'':'The movement of the field was not saved, so a moving field is not told from an object.'}
    /** The frames a flag's text names, [first, last] counted from 1: "frames 75 to 120", "frame 7"; null when it names none. */
    static List framesOf(String text){def m=(text?:'')=~/\bframes? (\d+)(?: to (\d+))?$/;m.find() ? [m.group(1) as int,(m.group(2)?:m.group(1)) as int] : null}
    /** The checks as the Step 4 page is given them: for every cell the flag's text, its rule and the first frame of the event.
     *  A cell with two flags has both in text and in rule, joined with "; ", and one by one in flags, each with its code and
     *  frames, the frames its text names (framesOf), which its link goes to; null for a flag that names none.
     *  codes: the cell's codes in the flags' order, qcFlags: them as the table's QC_flags column holds them. labels: each code
     *  as a table on screen names it, join: how codes are joined in a QC_flags column. */
    static Map answer(Map found) {
        if(found==null)return null
        def whole={v->String s=v==null?'':''+v;s==~/-?\d+/?(s as int):null}
        def each={Map row->String text=''+(row.Flag?:''),rule=''+(row.Rule?:'');text.isEmpty() ? [] : [text.split('; ') as List,rule.split('; ') as List].transpose().collect{[text:it[0],rule:it[1],code:CODES[it[1]]?:'',frames:framesOf(it[0] as String)]}}
        // a table of 2026.11.24.1 has no QC_flags: the codes come from its Rule column
        def codes={Map row->String given=row.containsKey('QC_flags')&&row.QC_flags!=null?''+row.QC_flags:codesOf(row.Rule);given.isEmpty() ? [] : (given.split(java.util.regex.Pattern.quote(JOIN)) as List)}
        [file:found.file,ringMeasured:found.ringMeasured==true,note:found.note?:'',labels:LABELS,join:JOIN,cells:((found.rows?:[]) as List).collectEntries{Map row->List mine=codes(row)
            [''+row.Cell,[text:''+(row.Flag?:''),rule:''+(row.Rule?:''),flags:each(row),codes:mine,qcFlags:mine.join(JOIN),firstFrame:whole(row.Event_first_frame),lastFrame:whole(row.Event_last_frame),detail:''+(row.Detail?:''),
                checked:row.Passing_object_check=='checked',why:''+(row.Not_checked_because?:''),ringMeasured:row.Ring_measured=='yes',movementMeasured:row.Field_movement_measured=='yes',
                nearSectors:whole(row.Near_ring_sectors_usable),area:whole(row.Area_pixels),crossings:whole(row.Crossings),longestCrossing:whole(row.Longest_crossing_samples)]]}]
    }
    /** Whether the table saved for an analysis is one written by 2026.11.24.1 (HEADER_BEFORE). */
    static boolean savedBefore(File dir,Map measurement){try{File source=file(measurement,dir);return source!=null&&source.isFile()&&CalciumCsv.read(source)[0]==HEADER_BEFORE}catch(Throwable ignored){return false}}
    /** The checks saved for an analysis, with the note rebuilt from the table; null when there are none. A table written by
     *  2026.11.24.1 (HEADER_BEFORE) reads too: its rows get the codes of their Rule column. */
    static Map saved(File dir,Map measurement) {
        try {
            File source=file(measurement,dir);if(source==null||!source.isFile())return null
            def rows=CalciumCsv.read(source);if(rows.isEmpty()||!(rows[0]==HEADER||rows[0]==HEADER_BEFORE))return null
            List columns=rows[0]
            List body=rows.drop(1).collect{List line->Map row=[columns,line].transpose().collectEntries{it};if(columns==HEADER_BEFORE)row.QC_flags=codesOf(row.Rule);row}
            boolean ring=body.any{it.Ring_measured=='yes'}
            String unchecked=ring?'':(body.find{it.Not_checked_because&&!(it.Not_checked_because in CELL_REASONS)}?.Not_checked_because?:'')
            return [rows:body,file:source.name,ringMeasured:ring,note:note(ring,body.any{it.Field_movement_measured=='yes'},unchecked)]
        }catch(Throwable ignored){return null}
    }
}

/** The manual review of the calls of one analysis. An evaluable cell can be marked the other way: a responder as not
 *  responding, a non-responder as responding, each with a reason (REASONS); Other needs a short note. A mark can be removed.
 *  The marks are kept beside the percent table, in "<recording>_reviewed_calls_<analysis>.csv", one row per marked cell; the
 *  file exists only while there is at least one mark. The engine's responder_calls.csv and the analysis folder are never
 *  written to. Every row names the analysis and the SHA-256 of the responder_calls.csv the mark was made on, and a file that
 *  does not match its analysis is not applied: the answer says so instead.
 *  The marks are applied in one place, CalciumAnalysisSettings.tally (with reviewed() here), and everything that counts
 *  responders reads that tally: the percent table, the finish sentence, the Step 4 page and the Fiji window.
 *  When a new analysis of the same recording finishes, the marks of the most recent earlier analysis made by this version are
 *  carried over for the cells whose automatic call and outline are unchanged (carry); the others are dropped and listed. */
class CalciumReview {
    /** The columns of a review file written by 2026.11.24.1, and of one this version writes while every mark is manual. It
     *  reads: every mark in it is manual. */
    static final List HEADER_BEFORE=['Cell','Automatic_call','Reviewed_call','Reason','Note','Time','Analysis','Responder_calls_sha256','ROI_sha256','Carried_from']
    /** Marked_by is last: manual, or automatic for a mark made at the end of an analysis (automatic). Written only while at
     *  least one mark is automatic. */
    static final List HEADER=HEADER_BEFORE+['Marked_by']
    static final List BY=['manual','automatic']
    /** The rules whose flag on a responder gets an automatic mark: a possible passing object (passes through, crosses the
     *  cell) and a rise on the edge. */
    static final List AUTOMATIC_RULES=['passes through','on the edge','crosses the cell']
    static final String AUTOMATIC_REASON='Passing object'
    /** "Cell 16, Cell 97": the first five names, then "and N more". */
    static String listed(List names){names.size()<=5?names.join(', '):names.take(5).join(', ')+' and '+(names.size()-5)+' more'}
    static final List REASONS=['Noise','Movement','Passing object','Bad outline','Other']
    static final List CALLS=['Responded','Did not respond']
    static final int NOTE_MAX=120
    /** A note is one line: no line break and no other control character. */
    static final java.util.regex.Pattern NOT_ONE_LINE=java.util.regex.Pattern.compile('[\\x00-\\x1F\\x7F\\x85\\u2028\\u2029]')
    /** A spreadsheet reads a cell that starts with one of these as a formula. */
    static final String FORMULA='=+-@'
    static final java.time.format.DateTimeFormatter TIME=java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ssXXX")
    /** A mark's time as the Step 4 page shows it, to the minute: "2026-11-24 10:05". */
    static String minute(String time){(time?:'').replace('T',' ').replaceAll(':\\d\\d(\\.\\d+)?([+-]\\d\\d:?\\d\\d|Z)?$','')}
    /** Why no mark can be set or removed while the review file is not applied (read gives the problem), and what to do. */
    static String blocked(String problem){problem+' Move it out of the folder to review this analysis again.'}
    static File file(Map measurement,File dir){measurement?.folder instanceof File?new File(measurement.folder as File,measurement.recording+'_reviewed_calls_'+dir.name+'.csv'):null}
    /** The call of a cell after review: the mark's when there is one for the call the engine made, else the engine's own. */
    static String reviewed(Map marks,String cell,String call){Map mark=marks?.get(cell) as Map;mark!=null&&mark.automatic==call&&call in CALLS?mark.reviewed as String:call}
    /** The engine's call of every cell for the primary stimulus, in the table's order. */
    static Map scored(File dir) {
        Map out=new LinkedHashMap();File source=new File(dir,'responder_calls.csv');if(!source.isFile())return out
        def calls=CalciumCsv.read(source);List h=calls[0];int cell=h.indexOf('Cell'),role=h.indexOf('Role'),status=h.indexOf('Status')
        if(cell<0||status<0)return out
        calls.drop(1).each{List row->if((role<0||row[role]=='Stimulus')&&!out.containsKey(''+row[cell]))out[''+row[cell]]=''+row[status]}
        out
    }
    /** The outline of every cell measured for an analysis, as its ROI registry names it (ROI_sha256); empty without a registry. */
    static Map outlines(File dir) {
        Map out=new LinkedHashMap();File source=new File(dir,'roi_registry.csv');if(!source.isFile())return out
        try{def rows=CalciumCsv.read(source);List h=rows[0];int cell=h.indexOf('Cell'),roi=h.indexOf('ROI_sha256'),measured=h.indexOf('Included_in_measurement')
            if(cell>=0&&roi>=0)rows.drop(1).each{List row->if(measured<0||row[measured]=='true')out[''+row[cell]]=''+row[roi]}}catch(Throwable ignored){}
        out
    }
    /** The cells in the order the page numbers them: the measurement's own, else the analysis's. */
    static List order(File dir,Map measurement){measurement?.names?(measurement.names as List).collect{''+it}:(outlines(dir).keySet() as List)?:(scored(dir).keySet() as List)}
    /** A cell as the page shows it: "Cell 12" for a long ID, counted from 1 in the measurement's order, else the name itself. */
    static String shown(List names,String cell){int at=names==null?-1:names.indexOf(cell);cell.length()>14&&at>=0?'Cell '+(at+1):cell}
    static String digest(File dir){CalciumInstaller.digest(new File(dir,'responder_calls.csv').bytes)}
    static boolean space(int c){Character.isWhitespace(c)||Character.isSpaceChar(c)}
    /** The note without the spaces at either end, no-break spaces too. */
    static String tidy(String text){int a=0,b=text.length();while(a<b&&space(text.codePointAt(a)))a+=Character.charCount(text.codePointAt(a));while(b>a&&space(text.codePointBefore(b)))b-=Character.charCount(text.codePointBefore(b));text.substring(a,b)}
    /** Whether a note holds a character that shows nothing, or turns the text around (Unicode format characters). */
    static boolean hidden(String text){for(int i=0;i<text.length();){int c=text.codePointAt(i);if(Character.getType(c)==Character.FORMAT)return true;i+=Character.charCount(c)};false}
    /** Why a mark cannot be kept for this cell, or null. @param calls the engine's calls */
    static String refused(String cell,Object call,Object reason,Object note,Map calls,String name) {
        String automatic=calls[cell] as String
        if(!(automatic in CALLS))return name+' has no call to review'+(automatic?': it is '+automatic.toLowerCase():'')+'.'
        if(!(call instanceof String)||!(call in CALLS))return 'Give the reviewed call: Responded or Did not respond.'
        if(call==automatic)return name+"'s automatic call is already "+automatic+'. The reviewed call must be '+CALLS.find{it!=automatic}+'.'
        if(!(reason instanceof String)||!(reason in REASONS))return 'Choose a reason: Noise, Movement, Passing object, Bad outline or Other.'
        if(note!=null&&!(note instanceof String))return 'The note must be text.'
        String text=note==null?'':note as String
        if(NOT_ONE_LINE.matcher(text).find())return 'The note must be one line.'
        if(hidden(text))return 'The note must be plain text, without hidden characters.'
        if(text.length()>NOTE_MAX)return 'Keep the note to '+NOTE_MAX+' characters or fewer.'
        if(!text.isEmpty()&&FORMULA.indexOf((int)text.charAt(0))>=0)return 'The note cannot start with =, +, - or @.'
        if(reason=='Other'&&tidy(text).isEmpty())return 'Other needs a short note.'
        null
    }
    /** The marks saved for an analysis. A byte order mark at the start and blank lines at the end are left out; an empty file
     *  holds no marks. A file without Marked_by (HEADER_BEFORE) reads as manual marks; in a file with it, every row says manual
     *  or automatic. An automatic mark is always Responded to Did not respond for Passing object.
     *  @return file (its name, or null when there is none), marks (cell to automatic, reviewed, reason, note, time, roi,
     *  carriedFrom, by) and problem: '' or why the file is not applied, in which case marks is empty */
    static Map read(File dir,Map measurement) {
        File source=file(measurement,dir);Map none=[file:null,marks:[:],problem:'']
        if(source==null||!source.isFile())return none
        String unused='The review file '+source.name
        try {
            String text=source.getText('UTF-8');if(text.startsWith('﻿'))text=text.substring(1)
            text=text.replaceAll(/[\r\n]+\s*$/,'\n')
            if(text.trim().isEmpty())return [file:source.name,marks:[:],problem:'']
            def rows=CalciumCsv.parse(text)
            if(rows.isEmpty()||!(rows[0]==HEADER||rows[0]==HEADER_BEFORE))return none+[file:source.name,problem:unused+' cannot be read, so its marks are not used.']
            List columns=rows[0]
            Map calls=scored(dir);String hash=digest(dir);Map marks=new LinkedHashMap()
            for(List line in rows.drop(1)) {
                Map row=[columns,line].transpose().collectEntries{it}
                // only a file without the column reads as manual: in one with it, an empty or missing Marked_by does not match
                String by=columns==HEADER_BEFORE?'manual':''+row.Marked_by
                if(row.Analysis!=dir.name)return none+[file:source.name,problem:unused+' is for another analysis, so its marks are not used.']
                if(row.Responder_calls_sha256!=hash||marks.containsKey(row.Cell)||row.Automatic_call!=calls[row.Cell]||refused(row.Cell as String,row.Reviewed_call,row.Reason,row.Note,calls,row.Cell as String)!=null||
                   !(by in BY)||(by=='automatic'&&!(row.Automatic_call=='Responded'&&row.Reason==AUTOMATIC_REASON)))
                    return none+[file:source.name,problem:unused+' does not match the calls of this analysis, so its marks are not used.']
                marks[row.Cell]=[automatic:row.Automatic_call,reviewed:row.Reviewed_call,reason:row.Reason,note:row.Note,time:row.Time,roi:row.ROI_sha256,carriedFrom:row.Carried_from,by:by]
            }
            return [file:source.name,marks:marks,problem:'']
        }catch(Throwable problem){return none+[file:source.name,problem:unused+' cannot be read, so its marks are not used.']}
    }
    /** Writes the marks and the percent table they give, both or neither. The tally is worked out from the marks first; each
     *  file is written beside its place and then moved in, the percent table first; when the review file cannot follow, the
     *  percent table is put back as it was. An empty set of marks removes the review file.
     *  @param traces dF/F0 by cell when already read (CalciumAnalysisSettings.tally)
     *  @return review (as read gives it) and tally */
    static Map write(File dir,Map measurement,Map marks,Map traces=null) {
        File target=file(measurement,dir);if(target==null)throw new IllegalStateException('This measurement has no folder to keep the review in.')
        Map review=[file:marks.isEmpty()?null:target.name,marks:marks,problem:'']
        Map tally=CalciumAnalysisSettings.tally(dir,measurement,traces,review)
        File percent=CalciumAnalysisSettings.percentFile(dir,measurement)
        def locked={File f->new IOException('Close '+f.name+' if it is open in another program, then try again.')}
        if(tally!=null&&percent.exists()&&!percent.canWrite())throw locked(percent)
        if(target.exists()&&!target.canWrite())throw locked(target)
        String tag='.'+UUID.randomUUID().toString().substring(0,8)+'.part'
        File reviewPart=new File(target.path+tag),percentPart=new File(percent.path+tag)
        try {
            if(!marks.isEmpty()) {
                String hash=digest(dir);Map roi=outlines(dir);List cells=scored(dir).keySet() as List
                // Marked_by only while a mark is automatic: a file of manual marks is the one 2026.11.24.1 writes, and reads
                boolean by=marks.values().any{(it as Map).by=='automatic'}
                CalciumCsv.write(reviewPart,[by?HEADER:HEADER_BEFORE]+(marks.keySet() as List).sort{cells.indexOf(it)}.collect{String cell->Map m=marks[cell] as Map;[cell,m.automatic,m.reviewed,m.reason,m.note?:'',m.time,dir.name,hash,roi[cell]?:'',m.carriedFrom?:'']+(by ? [m.by?:'manual'] : [])})
            }
            if(tally!=null)CalciumCsv.write(percentPart,CalciumAnalysisSettings.percentTable(tally))
        }catch(Throwable problem){reviewPart.delete();percentPart.delete();throw new IOException('The review could not be saved in '+target.parentFile.name+'.')}
        byte[] before=percent.isFile()?percent.bytes:null
        def move={File from,File to->java.nio.file.Files.move(from.toPath(),to.toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING)}
        if(tally!=null)try{move(percentPart,percent)}catch(Throwable problem){reviewPart.delete();percentPart.delete();throw locked(percent)}
        try{if(marks.isEmpty()){if(target.exists()&&!target.delete())throw new IOException(target.name)}else move(reviewPart,target)}
        catch(Throwable problem){
            reviewPart.delete()
            if(tally!=null)try{if(before==null)percent.delete() else percent.bytes=before}catch(Throwable ignored){}
            throw locked(target)
        }
        [review:review,tally:tally]
    }
    static String now(){TIME.format(java.time.OffsetDateTime.now())}
    /** Sets or removes one mark. body: cell, and either call, reason and note, or remove: true. Refuses with a plain message.
     *  @return what write gives: review and tally */
    static Map change(File dir,Map measurement,Map body,Map traces=null) {
        if(!(body?.cell instanceof String)||!(body.cell as String))throw new IllegalArgumentException('Say which cell to mark.')
        String cell=body.cell as String;Map calls=scored(dir);List names=order(dir,measurement);String name=shown(names,cell)
        if(!calls.containsKey(cell)||(measurement?.names&&!names.contains(cell)))throw new IllegalArgumentException('There is no cell '+cell+' in this analysis.')
        Map saved=read(dir,measurement)
        if(saved.problem)throw new IllegalStateException(blocked(saved.problem as String))
        Map marks=new LinkedHashMap(saved.marks as Map)
        if(body.remove==true){if(!marks.containsKey(cell))throw new IllegalArgumentException(name+' has no review mark.');marks.remove(cell)}
        else if(body.remove!=null&&body.remove!=false)throw new IllegalArgumentException('Remove must be true or false.')
        else {
            Object note=body.note instanceof String?tidy(body.note as String):body.note
            String why=refused(cell,body.call,body.reason,note,calls,name);if(why)throw new IllegalArgumentException(why)
            marks[cell]=[automatic:calls[cell],reviewed:body.call,reason:body.reason,note:(note?:'') as String,time:now(),carriedFrom:'',by:'manual']
        }
        write(dir,measurement,marks,traces)
    }
    /** The flag an automatic mark rests on, for every flagged cell of the saved checks: "possible passing object" or "rise on
     *  the edge", the names of CalciumCallChecks.LABELS. */
    static Map kinds(File dir,Map measurement) {
        Map out=[:]
        ((CalciumCallChecks.saved(dir,measurement)?.rows?:[]) as List).each{Map row->String rule=((''+(row.Rule?:'')).split('; ') as List).find{it in AUTOMATIC_RULES}
            if(rule)out[''+row.Cell]=(CalciumCallChecks.LABELS[CalciumCallChecks.CODES[rule]] as String).toLowerCase()}
        out
    }
    /** A reason or a kind of flag as the Methods sentence counts it: "1 passing object", "2 passing objects"; noise stays as it is. */
    static final Map PLURALS=['movement':'movements','passing object':'passing objects','bad outline':'bad outlines','other':'others','possible passing object':'possible passing objects','rise on the edge':'rises on the edge']
    static String counted(int n,String what){n+' '+(n>1?(PLURALS[what]?:what):what)}
    /** From the marks that changed a call of all cells; '' when none did. Manual marks only:
     *  "3 of 143 cells were reclassified after manual review: 2 from responder to non-responder (1 noise, 1 passing object) and
     *  1 from non-responder to responder (other)." Automatic marks only: "2 of 143 cells were reclassified automatically from
     *  responder to non-responder (1 possible passing object, 1 rise on the edge)." Both: "5 of 143 cells were reclassified:
     *  2 automatically from responder to non-responder (2 possible passing objects); 3 after manual review, 2 from responder to
     *  non-responder (1 noise, 1 passing object) and 1 from non-responder to responder (other)." A count above 1 takes the
     *  plural (counted).
     *  @param kinds the flag of each automatically marked cell (kinds); without it they are counted as passing object */
    static String methods(Map tally,Map review,Map kinds=null) {
        List changed=(tally?.changed?:[]) as List;if(changed.isEmpty())return ''
        int total=(((tally.rows as List).find{it.key=='all'} as Map).counted as List)[1] as int
        List machine=changed.findAll{((review.marks as Map)[it] as Map)?.by=='automatic'},marks=(changed-machine).collect{(review.marks as Map)[it] as Map}
        def part={List group,String from,String to->
            Map count=group.countBy{it.reason};List reasons=(count.keySet() as List).sort{a,b->(count[b] as int)<=>(count[a] as int)?:REASONS.indexOf(a)<=>REASONS.indexOf(b)}
            group.size()+' from '+from+' to '+to+' ('+(reasons.size()==1?(reasons[0] as String).toLowerCase():reasons.collect{counted(count[it] as int,(it as String).toLowerCase())}.join(', '))+')'}
        List down=marks.findAll{it.automatic=='Responded'},up=marks.findAll{it.automatic!='Responded'},parts=[]
        if(down)parts<<part(down,'responder','non-responder');if(up)parts<<part(up,'non-responder','responder')
        String start=changed.size()+' of '+total+' cells '+(changed.size()==1?'was':'were')+' reclassified'
        if(machine.isEmpty())return start+' after manual review: '+parts.join(' and ')+'.'
        List order=AUTOMATIC_RULES.collect{(CalciumCallChecks.LABELS[CalciumCallChecks.CODES[it]] as String).toLowerCase()}
        Map count=machine.countBy{kinds?.get(it)?:AUTOMATIC_REASON.toLowerCase()}
        List named=(count.keySet() as List).sort{a,b->(count[b] as int)<=>(count[a] as int)?:order.indexOf(a)<=>order.indexOf(b)}
        String which=' from responder to non-responder ('+named.collect{counted(count[it] as int,it as String)}.join(', ')+')'
        if(marks.isEmpty())return start+' automatically'+which+'.'
        start+': '+machine.size()+' automatically'+which+'; '+marks.size()+' after manual review, '+parts.join(' and ')+'.'
    }
    /** At the end of an analysis with the option on (CalciumAnalysisSettings.automaticMarks), after the checks and the carry-over
     *  of marks: every cell of the primary stimulus whose automatic call is Responded, flagged as a possible passing object or
     *  as a rise on the edge (AUTOMATIC_RULES; not a field that moved, not short separate crossings) and without a mark gets
     *  one: Did not respond, reason Passing object, no note, marked automatic. Saved like any mark (write), with its percent
     *  table or not at all. None while the review file does not match its analysis, and none when the checks were made but
     *  their table could not be saved: a mark rests on the saved checks.
     *  @param checked what CalciumCallChecks.run gave  @param on whether the option is on
     *  @param said the carry-over's sentence of this run (carry), which the finish message gives after this one: when it
     *         already ends with what to do about a file that cannot be written, that is not said twice
     *  @return null when the option is off; else folder (the analysis), cells (the cells marked), names (as the page shows
     *          them) and text, the sentence for the finish message: '' when no cell was marked */
    static Map automatic(File dir,Map measurement,Map checked,boolean on,File log=null,String said=null) {
        if(!on)return null
        String none='No automatic mark was made'
        if(checked?.rows&&!checked.file)return [folder:dir.name,cells:[],names:[],text:none+': the checks could not be saved.']
        try {
            Map saved=read(dir,measurement)
            if(saved.problem)return [folder:dir.name,cells:[],names:[],text:none+'. '+blocked(saved.problem as String)]
            Map calls=scored(dir),marks=new LinkedHashMap(saved.marks as Map);List names=order(dir,measurement),made=[]
            ((checked?.rows?:[]) as List).each{Map row->
                String cell=''+row.Cell
                if(calls[cell]!='Responded'||marks.containsKey(cell)||(measurement?.names&&!names.contains(cell)))return
                if(!((''+(row.Rule?:'')).split('; ') as List).any{it in AUTOMATIC_RULES})return
                marks[cell]=[automatic:'Responded',reviewed:'Did not respond',reason:AUTOMATIC_REASON,note:'',time:now(),carriedFrom:'',by:'automatic'];made<<cell
            }
            if(made.isEmpty())return [folder:dir.name,cells:[],names:[],text:'']
            write(dir,measurement,marks)
            List shownNames=made.collect{shown(names,it as String)}
            return [folder:dir.name,cells:made,names:shownNames,text:'Marked automatically: '+listed(shownNames)+'.']
        }catch(Throwable problem){
            CalciumCallChecks.logProblem(log,'the automatic marks',problem)
            boolean locked=problem instanceof IOException&&problem.message?.startsWith('Close ')
            return [folder:dir.name,cells:[],names:[],text:none+(locked?(said?.endsWith(problem.message)?'.':'. '+problem.message):': the marks could not be saved for this analysis.')]
        }
    }
    /** What follows the engine at the end of every analysis, on the Step 4 page and in the Fiji window alike: the checks on the
     *  calls (CalciumCallChecks.run), the carry-over of the review marks (carry), then the automatic marks (automatic).
     *  @param withRecording the measurement with the recording as it is now, which the checks may measure around the cells on
     *  @return checked, carried and marked (null when the option is off) */
    static Map afterAnalysis(File dir,Map withRecording,Map measurement,boolean automaticMarks,File log=null) {
        Map checked=CalciumCallChecks.run(dir,withRecording,log)
        Map carried=carry(dir,measurement,log)
        Map marked=automatic(dir,measurement,checked,automaticMarks,log,carried?.text as String)
        [checked:checked,carried:carried,marked:marked]
    }
    /** The finish message's part about the automatic marks: '' when none was made and nothing was refused. */
    static String automaticNote(Map marked){marked?.text?' '+marked.text:''}
    /** What carry gave for this analysis, until the next run, also once every mark it carried over is removed: null for
     *  another analysis. @param marks the marks of the analysis now, which do not change it */
    static Map current(Map carried,File dir,Map marks){carried!=null&&carried.folder==dir.name?carried.findAll{it.key!='folder'}:null}
    /** What the Step 4 page is given about the review of the analysis shown; problem is the whole sentence the page shows, with
     *  what to do (blocked), or ''. @param carried what carry gave when this analysis was run */
    static Map state(File dir,Map measurement,Map review,Map tally,Map carried=null) {
        List names=order(dir,measurement);Map marks=(review.marks?:[:]) as Map
        File source=review.file?file(measurement,dir):null
        [folder:dir.name,file:review.file,sha256:source?.isFile()?CalciumInstaller.digest(source.bytes):null,problem:review.problem?blocked(review.problem as String):'',
         marks:marks.collectEntries{cell,m->[cell,[name:shown(names,cell as String),automatic:m.automatic,reviewed:m.reviewed,reason:m.reason,note:m.note?:'',time:m.time,carriedFrom:m.carriedFrom?:'',markedBy:m.by?:'manual']]},
         after:scored(dir).collectEntries{cell,call->[cell,reviewed(marks,cell as String,call as String)]},
         changed:((tally?.changed?:[]) as List).size(),methods:methods(tally,review,marks.values().any{it.by=='automatic'}?kinds(dir,measurement):null),carried:current(carried,dir,marks),
         reasons:REASONS,calls:CALLS,noteMax:NOTE_MAX]
    }
    /** When "yyyyMMdd-HHmmss" an analysis was made: from its folder's name, else from when the file beside the tables was written. */
    static String made(String analysis,File source){def m=analysis=~/^SimpleCa-(\d{8}-\d{6})-/;m.find()?m.group(1):new java.text.SimpleDateFormat('yyyyMMdd-HHmmss').format(new Date(source.lastModified()))}
    /** At the end of a new analysis: the marks of the most recent earlier analysis of this recording made by this version (it
     *  has a checks file or a review file beside the tables) are kept for every cell whose automatic call and outline are as
     *  they were, and written for the new analysis; the others are dropped. Nothing is carried when that analysis has no marks,
     *  so marks removed there, or dropped there, do not come back from an older one. Manual marks only: automatic marks are
     *  left behind, neither carried nor counted.
     *  @return null when there is nothing to carry; else folder (the new analysis), from, kept, dropped (cell, name, why, text)
     *          and text, the sentence for the finish message without its leading space */
    static Map carry(File dir,Map measurement,File log=null) {
        String from=''
        try {
            if(!(measurement?.folder instanceof File))return null
            File folder=measurement.folder as File;String base=measurement.recording+'_'
            Map found=new LinkedHashMap()
            (folder.listFiles()?:new File[0]).each{File f->
                if(!f.isFile()||!f.name.endsWith('.csv'))return
                for(String kind in ['reviewed_calls_','call_checks_']){String start=base+kind
                    if(f.name.startsWith(start)){String analysis=f.name.substring(start.length(),f.name.length()-4);if(analysis!=dir.name&&!found.containsKey(analysis))found[analysis]=[analysis:analysis,when:made(analysis,f)]}}
            }
            List earlier=(found.values() as List).findAll{new File(new File(folder,it.analysis as String),'responder_calls.csv').isFile()}.sort{a,b->b.when<=>a.when?:b.analysis<=>a.analysis}
            if(earlier.isEmpty())return null
            from=earlier[0].analysis as String;File other=new File(folder,from)
            Map got=read(other,measurement)
            if(got.problem)return [folder:dir.name,from:from,kept:0,dropped:[],text:'The review marks of '+from+' were not carried over: its review file '+((got.problem as String).contains('cannot be read')?'cannot be read.':'does not match that analysis.')]
            // automatic marks are never carried over: the new analysis's own flags make them again when the option is on
            Map marks=(got.marks as Map).findAll{k,v->(v as Map).by!='automatic'};if(marks.isEmpty())return null
            Map calls=scored(dir),roi=outlines(dir),keep=new LinkedHashMap();List names=order(dir,measurement),before=order(other,null),dropped=[]
            marks.each{String cell,Map mark->
                String now=calls[cell]?:'',why=null
                if(now==mark.automatic&&mark.roi&&roi[cell]&&mark.roi!=roi[cell])why='outline changed'
                else if(now==mark.automatic)keep[cell]=mark.findAll{it.key!='roi'}+[carriedFrom:mark.carriedFrom?:from]
                else why=now in CALLS?'call changed':now?'now '+now.toLowerCase():'no longer measured'
                if(why){String name=calls.containsKey(cell)?shown(names,cell):shown(before,cell);if(!calls.containsKey(cell)&&name!=cell)name+=' of the earlier analysis'
                    dropped<<[cell:cell,name:name,why:why,text:name+': '+why]}
            }
            if(!keep.isEmpty())write(dir,measurement,keep)
            String text='Review marks from '+from+': '+(keep.isEmpty()?'none':keep.size())+' carried over, '+dropped.size()+' dropped.'+dropped.collect{' '+it.text+'.'}.join('')
            return [folder:dir.name,from:from,kept:keep.size(),dropped:dropped,text:text]
        }catch(Throwable problem){
            CalciumCallChecks.logProblem(log,'carrying the review marks over',problem)
            boolean locked=problem instanceof IOException&&problem.message?.startsWith('Close ')
            return [folder:dir.name,from:from,kept:0,dropped:[],text:'The review marks of '+(from?:'earlier analyses')+' were not carried over'+(locked?'. '+problem.message:': they could not be saved for this analysis.')]
        }
    }
    /** The finish message's part about the carried marks: '' when nothing was carried. */
    static String note(Map carried){carried?.text?' '+carried.text:''}
    /** The responding tile of the Fiji window after review: [value, caption], or null without marks. The value has one
     *  decimal, as the tile without marks, the finish message and the percent table (percentOf). */
    static List headline(Map tally,String stimulus) {
        if(tally?.reviewed!=true)return null
        Map all=(tally.rows as List).find{it.key=='all'} as Map;List t=all.counted as List,r=(all.reviewed as Map).counted as List
        [(r[1] as int)>0?CalciumAnalysisSettings.percentOf(r)+' %':'–',stimulus+': '+r[0]+' of '+r[1]+' evaluable cells after review; '+t[0]+' before review']
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
    /** What resultsJson was made from, so a review mark changes only the review and the tally in it. */
    Map resultsMap
    /** What the last analysis run here carried over from the review of an earlier one (CalciumReview.carry). */
    volatile Map carried
    /** The automatic marks the last analysis run here made (CalciumReview.automatic): null when the option was off. */
    volatile Map marked
    /** Whether the analysis shown was run with the option of the automatic marks on: a change of the option in the draft is a
     *  change of the settings (state, settingsPending), although the engine never sees it. */
    volatile boolean appliedAutomatic=false
    /** What is shown around the cells is read from (CalciumCallChecks.aroundContext), once per analysis: folder and context. */
    volatile Map aroundCache
    /** The parts of the finish message of the last analysis run here, so it can follow the review marks (finished). */
    volatile Map finish
    /** While the calls are checked after the engine has finished: the analysis can no longer be cancelled. */
    volatile boolean checking=false
    Map draft;long draftRevision=-1L
    String validationFingerprint='';Map draftValidation=[:]
    double[] appliedTime=new double[0]
    String appliedTimeSource=''
    ImagePlus movie;byte[] picture,backgroundPicture;int width=0,height=0
    /** What tells the attached recording from any other (CalciumSteps.printOf), to know it again when it was closed in Fiji and opened again. */
    String recordingPrint=''
    /** The measurement with the recording as this page reads frames from it (CalciumSteps.live): a recording closed in Fiji
     *  and opened again is the one the call checks measure around the cells. As it was when there is none. */
    Map withRecording(){ImagePlus now=CalciumSteps.live(movie,recordingPrint);now!=null?measurement+[movie:now,recordingPrint:recordingPrint]:measurement}
    float[] rest,work,work2,background;int[] tintBuffer;double[] rawRange,flatRange,meanTrace=new double[0];CalciumRise riseView
    int restFirst=1,restLast=1,frameEpoch=0
    java.awt.image.BufferedImage frameImage
    /** Whether the field moved (CalciumMotion.summary): from Fiji, or checked here when the session was started without it. */
    Map motion
    /** Where the starting stimulus frame came from ('fiji', 'last-recording' or 'last-session'), and Fiji's own values when the last session's were used. */
    String stimulusFrom='fiji';Map fijiDraft
    /** The remembered rule named no noise estimator and was given a new protocol's (CalciumAnalysisSettings.remembered). */
    boolean noiseMoved=false
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
                fijiDraft=CalciumAnalysisSettings.copy(draft);stimulusFrom='last-session';draft=before;restFirst=a as int;restLast=b as int;noiseMoved=carried.noiseMoved==true
                message='Your last settings for this recording are filled in.'+(carried.ruleMoved?' The rule now counts samples: 2 in a row. It was "for 1 s".':'')}
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
         recordedHz:(measurement.hz instanceof Number&&Double.isFinite(measurement.hz as double)&&(measurement.hz as double)>0)?measurement.hz:null,fijiFrames:measurement.fijiFrames instanceof Map?measurement.fijiFrames:null,stimulusFrom:stimulusFrom,fijiDraft:fijiDraft==null?null:CalciumAnalysisSettings.forPage(fijiDraft),motion:motion,motionNote:CalciumMotion.sentence(motion),selectionBasis:measurement.selectionBasis?:'Not recorded',
         rois:((measurement.rois?:new Roi[0]) as List).collect{Roi roi->def p=roi.getPolygon();(0..<p.npoints).collect{int i->[p.xpoints[i],p.ypoints[i]]}},
         draft:CalciumAnalysisSettings.forPage(draft),draftRevision:draftRevision,
         // what noise can be taken from, and what a new protocol takes it from: the page sends draft.protocol.sd_estimator back as given
         estimators:CalciumResponseProtocol.ESTIMATOR_ENGINE.collect{String key,String engine->[key:key,engine:engine,name:CalciumResponseProtocol.ESTIMATOR_NAME[key]]},newEstimator:CalciumResponseProtocol.newEstimator(),
         settings:CalciumAnalysisSettings.suggested(measurement).collect{k,v->[key:k,shown:CalciumAnalysisSettings.shown(v),exact:CalciumAnalysisSettings.full(v)]}]+timing()
    }
    /** What the status line says about the draft's noise before the first analysis. The page shows the estimator nowhere
     *  else. It is worked out from the draft each time, so it follows "Use Fiji's values", an import and a change of rule type. */
    synchronized String noiseLine() {
        String now=CalciumResponseProtocol.newEstimator()
        CalciumResponseProtocol.noiseNote(draft.protocol)?:(noiseMoved&&CalciumResponseProtocol.draftEstimator(draft.protocol)==now?' Noise is now the '+CalciumResponseProtocol.ESTIMATOR_NAME[now]+'.':'')
    }
    synchronized Map state() {
        boolean ready=false;try{ready=runtime.environment()!=null}catch(Throwable ignored){}
        String fingerprint=CalciumAnalysisSettings.fingerprint(draft)
        [key:measurementIdentity,busy:busy,cancellable:busy&&!checking,installing:installing,percent:percent,stage:stage,message:message+(results==null&&!busy&&tone!='bad'?noiseLine():''),tone:tone,elapsed:busy?(long)((System.currentTimeMillis()-started)/1000L):0L,
         version:version,folder:results?.name,engineReady:ready,logAvailable:job?.logFile?.isFile()?:false,rest:[restFirst,restLast],frameVersion:frameVersion(),
         settingsPending:results!=null&&(fingerprint!=appliedFingerprint||CalciumAnalysisSettings.automaticMarks(draft)!=appliedAutomatic),draftFingerprint:fingerprint,appliedFingerprint:appliedFingerprint,draftRevision:draftRevision,validation:validation(fingerprint)]+timing()
    }
    synchronized Map validation(String fingerprint) {
        if(validationFingerprint!=fingerprint){draftValidation=CalciumAnalysisSettings.diagnose(draft,measurement);validationFingerprint=fingerprint}
        draftValidation
    }
    /** Drafts may contain incomplete text. They receive feedback without running analysis or writing uploads. */
    synchronized Map updateDraft(Map value,Object revision=null) {
        if(value==null)throw new IllegalArgumentException('Provide the pending draft settings.')
        if(revision!=null){if(!(revision instanceof Number)||!Double.isFinite(revision as double)||(revision as double)!=(revision as long)||revision<0)throw new IllegalArgumentException('Draft revision must be a nonnegative integer.');if((revision as long)<draftRevision)return state();draftRevision=revision as long}
        draft=CalciumAnalysisSettings.tidyAutomatic(CalciumAnalysisSettings.copy(value))
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
                else {checking=true;stage='Checking the calls'
                    // the checks, the carry-over of the review marks and the automatic marks: the same function as in the Fiji window
                    Map after=CalciumReview.afterAnalysis(output,withRecording(),measurement,CalciumAnalysisSettings.automaticMarks(submitted),mine.logFile)
                    Map checked=after.checked as Map;carried=after.carried as Map;marked=after.marked as Map
                    load(output,fingerprint,submitted,checked);percent=100
                    String empty=CalciumAnalysisSettings.noEvaluableCells(output)
                    if(!empty)try{CalciumAnalysisSettings.savePercent(output,measurement,resultsMap?.tally as Map)}catch(Throwable problem){CalciumCallChecks.logProblem(mine.logFile,'saving the percent table',problem)}
                    synchronized(this){finish=[folder:output.name,empty:empty,saved:'Complete. Protocol, inputs and tables saved in '+output.name+'.',noise:CalciumResponseProtocol.noiseNote(submitted?.protocol),marked:CalciumReview.automaticNote(marked)];message=finished();tone=empty?'warn':'good'}}
            }catch(java.util.concurrent.CancellationException ignored){message='Cancelled. Completed analyses are kept; partial output is not used.';tone='plain';stage='';percent=0}
            catch(Throwable t){message=t.message?:t.toString();tone='bad'}
            finally{busy=false;checking=false}
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
        Map submitted=CalciumAnalysisSettings.tidyAutomatic(CalciumAnalysisSettings.copy(incoming)),uploads=[:]
        Map validation=CalciumAnalysisSettings.diagnose(submitted,measurement)
        if(validation.errors)throw new IllegalArgumentException(validation.errors.join(' '))
        CalciumAnalysisSettings.remember(measurement,submitted)
        CalciumAnalysisSettings.FILE_KEYS.each{k,v->uploads[k]=CalciumAnalysisSettings.upload(submitted[k],k)}
        if(uploads.timestamps)CalciumMeasurement.readTimesRows(CalciumCsv.parse(uploads.timestamps.content as String),measurement.frames as int)
        // Check everything independent of filesystem writes before creating the request snapshot.
        Map flags=[:];CalciumAnalysisSettings.FILE_KEYS.each{k,v->if(uploads[k])flags[v]=true}
        Map checked=CalciumAnalysisSettings.check(submitted.values as Map,submitted.protocol as Map,confirmed,measurement,flags,true)
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
    /** Cancels the engine. Once it has finished, the checks on the calls and the review carry-over run to the end: the
     *  analysis is saved by then, so there is nothing left to cancel (state gives cancellable: false). */
    void cancel(){if(checking)return;job?.cancel();if(busy)stage='Cancelling...'}
    /** The finish message of the last analysis run here, with the figures of the tally as it is now, so it follows the review
     *  marks; the parts about the automatic marks, after the figures, and about carried marks stay until the next run. */
    synchronized String finished() {
        Map parts=finish;if(parts==null||results?.name!=parts.folder)return message
        String carry=CalciumReview.note(CalciumReview.current(carried,results,(resultsMap?.review?.marks?:[:]) as Map))
        parts.empty?(parts.empty as String)+(parts.marked?:'')+carry:(parts.saved as String)+CalciumAnalysisSettings.percentSentence(resultsMap?.tally as Map)+(parts.marked?:'')+(parts.noise?:'')+carry
    }
    static Map table(File dir,String name){File file=new File(dir,name);if(!file.isFile())return [header:[],rows:[]];def rows=CalciumCsv.read(file);[header:rows?rows[0]:[],rows:rows?rows.drop(1):[]]}
    static List<Map> records(Map table){(table.rows as List).collect{line->table.header.withIndex().collectEntries{h,i->[h,i<line.size()?line[i]:'']}}}
    /** Draw only completed files. No response classifications or scientific metrics are recomputed in Java.
     *  @param checked the checks on the calls just made (CalciumCallChecks.run); without them the ones saved for this analysis are read */
    void load(File dir,String fingerprint=null,Map submitted=null,Map checked=null) {
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
        Map review=CalciumReview.read(dir,measurement),tally=CalciumAnalysisSettings.tally(dir,measurement,traces,review)
        CalciumAnalysisSettings.syncPercent(dir,measurement,tally)   // after-review columns follow the marks as they are now
        Map corrected=grayTable(measurement.csv instanceof File?measurement.csv as File:null)
        Map raw=grayTable(tableFolder==null?null:new File(tableFolder,measurement.recording+'_gray_values_raw.csv'))
        Map loaded=[folder:dir.name,time:time,frameTimes:time,timeSource:manifest.time_source?:'constant_fps',timingProvisional:false,mean:mean,sem:numbers(average,2),traces:traces,calls:calls,counts:counts,windows:windows,callsByWindow:callsByWindow,countsByWindow:countsByWindow,
            baseline:[time[first-1],time[last-1]],baselineFrames:[first,last],stimulus:[CalciumCsv.number(firstCall.Start_s?:manifest.stimulus_start),CalciumCsv.number(firstCall.Effective_end_s?:firstCall.End_s?:manifest.stimulus_end)],
            corrected:corrected,raw:raw,chosenBlind:(measurement.chosenBlind instanceof Map?measurement.chosenBlind:null),tally:tally,review:CalciumReview.state(dir,measurement,review,tally,carried),checks:CalciumCallChecks.answer(checked?:CalciumCallChecks.saved(dir,measurement)),appliedDraft:submitted==null?null:CalciumAnalysisSettings.copy(submitted),protocol:protocol('protocol.json'),appliedProtocol:protocol('applied_protocol.json'),manifest:manifest,appliedFingerprint:applied,tables:tables]
        String text=json(loaded)
        synchronized(this){meanTrace=mean as double[];restFirst=first;restLast=last;rest=null;frameEpoch++;appliedTime=time as double[];appliedTimeSource=loaded.timeSource;appliedFingerprint=applied
            if(submitted!=null)appliedAutomatic=CalciumAnalysisSettings.automaticMarks(submitted)
            aroundCache=null;results=dir;resultsMap=loaded;resultsJson=text;version++}
    }
    /** What is shown around one cell of the analysis shown (CalciumCallChecks.aroundCell), for the cell asked for: what the
     *  check used, never every cell at once. The files are read once per analysis, outside the session's lock.
     *  @param folder when given, must be the analysis shown  @return folder, and what aroundCell gives */
    Map around(Object cell,Object folder=null) {
        File dir;Map checks
        synchronized(this){
            if(results==null||resultsMap==null)throw new IllegalStateException('Run an analysis first.')
            if(folder!=null&&folder!=results.name)throw new IllegalArgumentException('This is for another analysis. Reload the page.')
            dir=results;checks=resultsMap.checks as Map
        }
        if(!(cell instanceof String)||!(cell as String))throw new IllegalArgumentException('Say which cell.')
        Map cached=aroundCache,context=cached!=null&&cached.folder==dir.name?cached.context as Map:null
        if(context==null){context=CalciumCallChecks.aroundContext(dir,measurement);synchronized(this){if(results?.name==dir.name)aroundCache=[folder:dir.name,context:context]}}
        [folder:dir.name]+CalciumCallChecks.aroundCell(context,checks,cell as String)
    }
    /** Sets or removes one review mark (CalciumReview.change) on the analysis shown, and answers with the review and the tally
     *  after it. The saved results the page reads change with it: their review and tally, nothing else.
     *  @param body cell, and call, reason and note, or remove: true; folder, when given, must be the analysis shown */
    synchronized Map review(Map body) {
        if(busy)throw new IllegalStateException('Wait for the analysis to finish.')
        if(results==null||resultsMap==null)throw new IllegalStateException('Run an analysis first.')
        if(body?.folder!=null&&body.folder!=results.name)throw new IllegalArgumentException('This mark is for another analysis. Reload the page.')
        Map changed=CalciumReview.change(results,measurement,body,resultsMap.traces as Map)
        Map tally=changed.tally as Map,state=CalciumReview.state(results,measurement,changed.review as Map,tally,carried)
        resultsMap.tally=tally;resultsMap.review=state;resultsJson=json(resultsMap)
        // the finish message of the analysis shown, in its own tone: also after a run that failed or was cancelled since
        if(finish?.folder==results.name&&!busy){message=finished();tone=finish.empty?'warn':'good'}
        version++
        [folder:results.name,version:version,review:state,tally:tally,message:message]
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
                    case 'api/around':replyJson(exchange,200,now.around(query.cell,query.folder));return   // what the check used around one cell
                }
            } else if(method=='POST'&&now!=null) {
                byte[] sent=read(exchange.getRequestBody(),1<<26)
                Object parsed=null
                if(sent.length>0)try{parsed=new groovy.json.JsonSlurper().parseText(new String(sent,'UTF-8'))}catch(Throwable bad){throw new IllegalArgumentException('The request is not valid.')}
                if(parsed!=null&&!(parsed instanceof Map))throw new IllegalArgumentException('The request is not valid.')
                Map body=parsed==null ? [:] : (parsed as Map)
                switch(what) {
                    case 'api/draft':replyJson(exchange,200,now.updateDraft(body.draft as Map,body.revision));return
                    case 'api/protocol':replyJson(exchange,200,[protocol:CalciumResponseProtocol.exchange(body.protocol)]);return
                    case 'api/analyze':now.analyze(body.draft as Map,body.confirmed==true);replyJson(exchange,200,now.state());return
                    case 'api/setup':now.setup();replyJson(exchange,200,now.state());return
                    case 'api/cancel':now.cancel();replyJson(exchange,200,now.state());return
                    case 'api/export':replyJson(exchange,200,now.export(body));return
                    case 'api/review':replyJson(exchange,200,now.review(body));return
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
    static final String HASH = "bcb25fe397fb83f1c691bc82519ca8b466963eed90a2b059fb39221507d8186c"
    static final String DATA = "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"
}
// END GENERATED SIMPLECA BUNDLE

// BEGIN GENERATED WEB BUNDLE
class CalciumWebBundle {
    static final String SHA256 = "20f0d19effd6703e98d9af7dc219d324f0668858bbab2f3fd038961e4c2c36ab"
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

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
    /** The colour of each part a zone of the flag can be; a colour met twice for one cell gives way to a spare one. */
    static final Map COLOURS=[before:new Color(0xD55E00),after:new Color(0x0072B2),rim:new Color(0xAA3377),ring:new Color(0xD55E00),core:new Color(0x009E73)]
    static final List SPARE=[new Color(0xCC79A7),new Color(0xE69F00),new Color(0x56B4E9),new Color(0x882255)]
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
    /** The time between two frames, as most frames are apart. */
    double step(){double[] x=times();if(x.length<2)return 1d;List d=(1..<x.length).collect{x[it]-x[it-1]}.sort();d[(int)(d.size()/2)] as double}
    /** The lines under the cell, in the order they are drawn, each [colour, width, name, trace]: every other usable near-ring
     *  sector thin and gray, then the zones the flag rests on, each in the colour of its part and named by its place. */
    List lines() {
        if(!drawn())return []
        Map byZone=((data.zones?:[]) as List).collectEntries{[(it.zone as int),it]}
        List out=[],used=[]
        ((data.others?:[]) as List).each{def z->double[] y=byZone[z as int]?.trace as double[];if(y!=null)out<<[OTHER,0.9f,'other ring sectors',y]}
        ((data.lines?:[]) as List).each{Map line->double[] y=byZone[line.zone as int]?.trace as double[];if(y==null)return
            Color c=COLOURS[line.part] as Color;if(c==null||used.contains(c))c=(SPARE.find{!used.contains(it)}?:CalciumTheme.MUTED) as Color;used<<c
            out<<[c,1.5f,''+(line.place?:CalciumCallChecks.place(line.zone as int)),y]}
        out
    }
    /** The bands, [from, to, colour, name] in seconds: the cell's response window as on Traces, then the frames the flag names,
     *  from half a frame before the first to half a frame after the last. */
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
    void showTable(javax.swing.JTable target,javax.swing.table.TableModel model){target.model=model;for(int i=0;i<target.columnCount;i++)target.columnModel.getColumn(i).preferredWidth=160}
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
     *  with the first frame of its event is a link to that frame. A control window keeps its automatic call: it has neither.
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
    /** A row of the selected cell's flags, each in the checks' own words; a flag with the first frame of its event is a link to
     *  that frame. Then, for an evaluable cell that could not be checked in a recording whose surroundings were measured, why
     *  (else checks.note says it once). Shown while it holds something. */
    void fillFlags(javax.swing.JPanel into,Map cell,boolean evaluable) {
        def owner=this;into.removeAll()
        ((cell?.flags?:[]) as List).each{Map flag->
            if(into.componentCount>0)into.add(CalciumTheme.label('·',CalciumTheme.MUTED))
            def label=CalciumTheme.label(flag.text as String,cell.firstFrame!=null?CalciumTheme.BLUE:CalciumTheme.MUTED)
            String detail=cell.detail?:''
            if(cell.firstFrame!=null) {
                int first=cell.firstFrame as int
                label.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR))
                label.setToolTipText('<html>Go to frame '+first+(detail?'<br>'+CalciumOutlineProgress.escape(detail):'')+'</html>')
                label.addMouseListener(new java.awt.event.MouseAdapter(){void mouseClicked(java.awt.event.MouseEvent e){if(javax.swing.SwingUtilities.isLeftMouseButton(e))owner.goToFrame(first)}})
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
    /** A flag's link: the Traces tab, the recording rather than the still colour picture, at the first frame of the event. */
    void goToFrame(int first){tabs.selectedIndex=1;if(movie!=null&&view.selectedIndex==0)view.selectedIndex=2;stopPlaying();showFrame(first)}
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
     *  @return null when the cell does not rise at the peak; else first and last (the event's frames, counted from 1), rule
     *          ('passes through', 'on the edge', 'crosses the cell' or ''), text, detail; traces (what traces gave); where a
     *          passage was seen, before and after ([zone, frame counted from 1] of the sector seen before and after the event);
     *          where the rise is on an edge, edge ([rim sector 0 to 7, the stretch's first frame, its last frame, counted from
     *          1]); and cross (what crossing gave; null when it was not checked) */
    static Map<String, Object> check(double[] raw, double[] background, float[][] zone, int firstF0, int lastF0, int peak, int span, double[][] centre = null, double[][] outline = null) {
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
        double[] entering = null, leaving = null
        double strength = 0d
        for (double[] x : before) for (double[] y : after) {
            int apart = Math.abs(((int) x[0]) % 8 - ((int) y[0]) % 8)
            apart = Math.min(apart, 8 - apart)
            if (apart < APART) continue
            double s = Math.min(x[1], y[1]) / rise
            if (entering == null || s > strength) { entering = x; leaving = y; strength = s }
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
            } else out.put('detail', ((String) out.get('detail')) + '; also crosses the cell ' + seen)
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
     *          (their sectors 0 to 7: right, lower right ...), travel (pixels); and with a travel, width (pixels across the
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
     *  @param measure false to read the files beside the tables only, never measure (around), as what is shown around a cell does */
    static Map context(File dir,Map measurement,File log=null,boolean measure=true) {
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
        Map shapes=ring.values!=null?shapes(measurement,names,ring):null   // where the zones lie, for the rule "crosses the cell"
        [header:h,byCell:byCell,statusAt:statusAt,samples:samples,time:time,frames:frames,traceAt:traceAt,first:first,last:last,names:names,span:span,measured:measured,ring:ring,shifts:shifts,unchecked:unchecked,shapes:shapes]
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
                event=CalciumPassingObject.check(raw,background,zone,first,last,peak,span,shape?.centre as double[][],shape?.outline as double[][])
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
                        String detail=event.detail+(crossNote?'; '+crossNote:'')
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
     *  tables only (nothing is measured). Read once per analysis; aroundCell takes one cell from it.
     *  @return the context, or why: the sentence that says why nothing can be shown for any cell */
    static Map aroundContext(File dir,Map measurement) {
        try{return context(dir,measurement,null,false)}
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
    /** The checks as the Step 4 page is given them: for every cell the flag's text, its rule and the first frame of the event.
     *  A cell with two flags has both in text and in rule, joined with "; ", and one by one in flags, each with its code.
     *  codes: the cell's codes in the flags' order, qcFlags: them as the table's QC_flags column holds them. labels: each code
     *  as a table on screen names it, join: how codes are joined in a QC_flags column. */
    static Map answer(Map found) {
        if(found==null)return null
        def whole={v->String s=v==null?'':''+v;s==~/-?\d+/?(s as int):null}
        def each={Map row->String text=''+(row.Flag?:''),rule=''+(row.Rule?:'');text.isEmpty() ? [] : [text.split('; ') as List,rule.split('; ') as List].transpose().collect{[text:it[0],rule:it[1],code:CODES[it[1]]?:'']}}
        // a table of 2026.11.24.1 has no QC_flags: the codes come from its Rule column
        def codes={Map row->String given=row.containsKey('QC_flags')&&row.QC_flags!=null?''+row.QC_flags:codesOf(row.Rule);given.isEmpty() ? [] : (given.split(java.util.regex.Pattern.quote(JOIN)) as List)}
        [file:found.file,ringMeasured:found.ringMeasured==true,note:found.note?:'',labels:LABELS,join:JOIN,cells:((found.rows?:[]) as List).collectEntries{Map row->List mine=codes(row)
            [''+row.Cell,[text:''+(row.Flag?:''),rule:''+(row.Rule?:''),flags:each(row),codes:mine,qcFlags:mine.join(JOIN),firstFrame:whole(row.Event_first_frame),lastFrame:whole(row.Event_last_frame),detail:''+(row.Detail?:''),
                checked:row.Passing_object_check=='checked',why:''+(row.Not_checked_because?:''),ringMeasured:row.Ring_measured=='yes',movementMeasured:row.Field_movement_measured=='yes',
                nearSectors:whole(row.Near_ring_sectors_usable),area:whole(row.Area_pixels),crossings:whole(row.Crossings),longestCrossing:whole(row.Longest_crossing_samples)]]}]
    }
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
    /** The columns of a review file written by 2026.11.24.1. It still reads: every mark in it is manual. */
    static final List HEADER_BEFORE=['Cell','Automatic_call','Reviewed_call','Reason','Note','Time','Analysis','Responder_calls_sha256','ROI_sha256','Carried_from']
    /** Marked_by is last: manual, or automatic for a mark made at the end of an analysis (automatic). */
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
     *  holds no marks. A file without Marked_by (HEADER_BEFORE) reads as manual marks, and gets the column when it is next
     *  written; in a file with it, every row says manual or automatic. An automatic mark is always Responded to Did not respond
     *  for Passing object.
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
                CalciumCsv.write(reviewPart,[HEADER]+(marks.keySet() as List).sort{cells.indexOf(it)}.collect{String cell->Map m=marks[cell] as Map;[cell,m.automatic,m.reviewed,m.reason,m.note?:'',m.time,dir.name,hash,roi[cell]?:'',m.carriedFrom?:'',m.by?:'manual']})
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
    /** From the marks that changed a call of all cells; '' when none did. Manual marks only:
     *  "3 of 143 cells were reclassified after manual review: 2 from responder to non-responder (1 noise, 1 passing object) and
     *  1 from non-responder to responder (other)." Automatic marks only: "2 of 143 cells were reclassified automatically from
     *  responder to non-responder (1 possible passing object, 1 rise on the edge)." Both: "5 of 143 cells were reclassified:
     *  2 automatically from responder to non-responder (2 possible passing object) and 3 after manual review: ..."
     *  @param kinds the flag of each automatically marked cell (kinds); without it they are counted as passing object */
    static String methods(Map tally,Map review,Map kinds=null) {
        List changed=(tally?.changed?:[]) as List;if(changed.isEmpty())return ''
        int total=(((tally.rows as List).find{it.key=='all'} as Map).counted as List)[1] as int
        List machine=changed.findAll{((review.marks as Map)[it] as Map)?.by=='automatic'},marks=(changed-machine).collect{(review.marks as Map)[it] as Map}
        def part={List group,String from,String to->
            Map count=group.countBy{it.reason};List reasons=(count.keySet() as List).sort{a,b->(count[b] as int)<=>(count[a] as int)?:REASONS.indexOf(a)<=>REASONS.indexOf(b)}
            group.size()+' from '+from+' to '+to+' ('+(reasons.size()==1?(reasons[0] as String).toLowerCase():reasons.collect{count[it]+' '+(it as String).toLowerCase()}.join(', '))+')'}
        List down=marks.findAll{it.automatic=='Responded'},up=marks.findAll{it.automatic!='Responded'},parts=[]
        if(down)parts<<part(down,'responder','non-responder');if(up)parts<<part(up,'non-responder','responder')
        String start=changed.size()+' of '+total+' cells '+(changed.size()==1?'was':'were')+' reclassified'
        if(machine.isEmpty())return start+' after manual review: '+parts.join(' and ')+'.'
        List order=AUTOMATIC_RULES.collect{(CalciumCallChecks.LABELS[CalciumCallChecks.CODES[it]] as String).toLowerCase()}
        Map count=machine.countBy{kinds?.get(it)?:AUTOMATIC_REASON.toLowerCase()}
        List named=(count.keySet() as List).sort{a,b->(count[b] as int)<=>(count[a] as int)?:order.indexOf(a)<=>order.indexOf(b)}
        String which=' from responder to non-responder ('+named.collect{count[it]+' '+it}.join(', ')+')'
        if(marks.isEmpty())return start+' automatically'+which+'.'
        start+': '+machine.size()+' automatically'+which+' and '+marks.size()+' after manual review: '+parts.join(' and ')+'.'
    }
    /** At the end of an analysis with the option on (CalciumAnalysisSettings.automaticMarks), after the checks and the carry-over
     *  of marks: every cell of the primary stimulus whose automatic call is Responded, flagged as a possible passing object or
     *  as a rise on the edge (AUTOMATIC_RULES; not a field that moved, not short separate crossings) and without a mark gets
     *  one: Did not respond, reason Passing object, no note, marked automatic. Saved like any mark (write), with its percent
     *  table or not at all. None while the review file does not match its analysis.
     *  @param checked what CalciumCallChecks.run gave  @param on whether the option is on
     *  @param said the carry-over's sentence of this run (carry), which the finish message gives after this one: when it
     *         already ends with what to do about a file that cannot be written, that is not said twice
     *  @return null when the option is off; else folder (the analysis), cells (the cells marked), names (as the page shows
     *          them) and text, the sentence for the finish message: '' when no cell was marked */
    static Map automatic(File dir,Map measurement,Map checked,boolean on,File log=null,String said=null) {
        if(!on)return null
        String none='No automatic mark was made'
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
    static final String SHA256 = "4e37d659d5277e8cd008e32aae1a2319295959e3c510d7aee7bc7be1f44969f5"
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

// Ca²⁺ Imaging (Wang Lab)                                                                  version 2026.11.23.1
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
        say(CalciumRelease.NAME + ", version 2026.11.23.1")
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
    static final String VERSION = "2026.11.23.1"
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
    /** Under the plot: the selected cell's threshold, what set it, and the run that met it. */
    String ruleLine(){'Local threshold '+(responseCall.Threshold_increase?:'unavailable')+' ('+CalciumCalls.thresholdSetBy(responseCall)+') · span '+(responseCall.Run_duration_s?:'unavailable')+' s · '+(responseCall.Run_samples?:'')+' samples'}
    static final String DASHED='Green dashed: local threshold transformed to this global axis'
    /** Under the plot on screen, and never in a saved figure (nothing but the axes and the key is drawn on the plot): the
     *  caption, the rule line and what the green dashed line is, one line each, so they take the same room for every cell. */
    final List<CalciumLine> noteLines=[new CalciumLine(CalciumTheme.MUTED),new CalciumLine(CalciumTheme.MUTED),new CalciumLine(CalciumTheme.GOOD)]
    final javax.swing.JPanel notes=new javax.swing.JPanel()
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
    void showNotes(){List words=time.length<2 ? ['','',''] : [caption,ruleLine(),DASHED];noteLines.eachWithIndex{line,i->line.setText(words[i] as String)}}
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
        if(cell&&selected.length>0)entries<<[CalciumTheme.INK,1.7f,cell]   // a cell without a trace has no entry: why is said under the plot
        if(many) {
            int yes=allCalls.count{it=='Responded'} as int,no=allCalls.count{it=='Did not respond'} as int
            if(yes>0)entries<<[CalciumCalls.RESPONDED,1.4f,'responded ('+yes+')']
            if(no>0)entries<<[CalciumCalls.NOT,1.4f,'did not respond ('+no+')']
            if(all.length-yes-no>0)entries<<[CalciumCalls.OTHER,1.4f,(yes+no>0?'no call':allLabel)+' ('+(all.length-yes-no)+')']
        }
        int line=0,run=0;def places=[]
        entries.each{entry->int need=21+fm.stringWidth(entry[2] as String)+14;if(run>0&&run+need>w){line++;run=0};places<<[run,line];run+=need}
        t=30+16*(line+1);h=Math.max(1,getHeight()-t-46)
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
        // the bands are not named on the plot (only the axes and the key are drawn on it); the tooltip names the one under the pointer
        List drawnBands=[]
        [[baselineStart,baselineEnd,new Color(0xEEF1F4),'normalization'],[localStart,localEnd,new Color(0xDCEBFA),'local baseline'],
         [stimulusStart,stimulusEnd,new Color(0xFFF0DC),'response'],[runStart,runEnd,new Color(0xC3EBC8),'qualifying run']].each { band ->
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
        entries.eachWithIndex{entry,int i->int kx=l+(places[i][0] as int),ky=20+16*(places[i][1] as int);String text=entry[2] as String
            g.setColor(entry[0] as Color);g.setStroke(new java.awt.BasicStroke(entry[1] as float));g.drawLine(kx,ky,kx+16,ky);g.setStroke(new java.awt.BasicStroke(1f));g.setColor(CalciumTheme.INK);g.drawString(text,kx+21,ky+4)}
        g.setColor(CalciumTheme.MUTED);String xl='Time (s)';g.drawString(xl,l+(int)((w-fm.stringWidth(xl))/2),getHeight()-8)
        def turned=(java.awt.Graphics2D)g.create();try{turned.rotate(-Math.PI/2);String yl='Global ΔF/F₀';turned.drawString(yl,-(t+(int)((h+fm.stringWidth(yl))/2)),14)}finally{turned.dispose()}
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
    def cells=new javax.swing.JComboBox<String>()
    def stimuli=new javax.swing.JComboBox<String>()
    def group=new javax.swing.JComboBox<String>(['All cells'] as String[])
    def allTraces=new javax.swing.JCheckBox('Every trace',true)
    /** The selected cell's call line; it wraps at the width of the Traces tab, so a marked cell's line is read whole. While
     *  results are shown it keeps room for two lines, so the picture and the plot under it keep their size from cell to cell. */
    def cellStatus=new CalciumText(' ',CalciumTheme.MUTED).with{lines=1;it}
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
    /** The checks of the analysis shown, as the Step 4 page is given them (CalciumCallChecks.answer); null when there are none. */
    Map checks
    /** The automatic call of every cell for the primary stimulus (CalciumReview.scored). */
    Map automatic=[:]
    /** The parts of the finish message of the last analysis run here (finished), so that it follows the review marks. */
    Map finish
    // Under the selected cell's call line: its flags, then its mark controls, one button or the form that sets a mark.
    def flagRow=CalciumTheme.row(),markRow=CalciumTheme.row(),reviewBox=new javax.swing.JPanel()
    def markButton=CalciumTheme.button('Mark as not a responder'),removeMark=CalciumTheme.button('Remove mark'),saveMark=CalciumTheme.button('Save',true),cancelMark=CalciumTheme.button('Cancel')
    def reasonLabel=CalciumTheme.label('Reason'),noteLabel=CalciumTheme.label('Note')
    def reasonBox=new javax.swing.JComboBox<String>(CalciumReview.REASONS as String[])
    def noteField=new CalciumHintField('Optional',CalciumReview.NOTE_MAX,16)
    /** A mark refused, in the words of CalciumReview, shown under the form until the cell changes; or why the review file is not
     *  applied. Two lines are kept for it while results are shown, so a refusal that names a file of a recording does not
     *  move the picture and the plot. */
    def markProblem=new CalciumText('&nbsp;',CalciumTheme.WARN).with{lines=2;it}
    /** The cell whose mark form is open, or null; the cell the refusal under the form is about, and the words said. */
    String formCell, problemCell, problemSaid
    /** The reason already chosen when the form opens on a flagged responder, by the rule of its flag. */
    static final Map REASON_OF_RULE=['passes through':'Passing object','on the edge':'Passing object','field moved':'Movement','short separate crossings':'Noise']
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
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->def c=new java.awt.GridBagConstraints();c.gridx=i;c.fill=java.awt.GridBagConstraints.BOTH;c.weightx=i==0?2.2d:1d;c.insets=new java.awt.Insets(0,0,0,i<3?8:0);tile.wide=i==0?220:100;stats.add(tile,c)}
        // the Traces tab: which cells, the picture at the cursor's moment, the traces, and the time slider under both
        [cells,group,view,heatmapOrder,stimuli,reasonBox].each{it.setFont(CalciumTheme.font())}
        previousCell.setToolTipText('Previous cell (left arrow)');nextCell.setToolTipText('Next cell (right arrow)')
        group.setToolTipText('List and draw only the cells with this response call');allTraces.setToolTipText('Draw every listed cell behind the mean, coloured by its response call')
        def white={javax.swing.JPanel p->p.setOpaque(true);p.setBackground(Color.white);p}
        def pickers=white(row(CalciumTheme.label('Show'),group,CalciumTheme.label('Cell'),cells,previousCell,nextCell,CalciumTheme.label('Window'),stimuli,allTraces))
        // the call line, then the cell's flags, then its mark controls, then a refusal; while results are shown each keeps its
        // room when it has nothing to show (a strut at the end of a row), so the head has one height for every cell
        reviewBox.setLayout(new javax.swing.BoxLayout(reviewBox,javax.swing.BoxLayout.Y_AXIS));white(reviewBox);reviewBox.setBorder(javax.swing.BorderFactory.createEmptyBorder(0,8,6,8))
        [markButton,removeMark,reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{markRow.add(it)}
        markRow.add(javax.swing.Box.createRigidArea(new java.awt.Dimension(0,(int)[markButton,removeMark,reasonBox,noteField,saveMark,cancelMark].collect{it.preferredSize.height}.max())))
        [cellStatus,flagRow,markRow,markProblem].each{reviewBox.add(it)}
        def head=new javax.swing.JPanel(new java.awt.BorderLayout());head.setBackground(Color.white);head.add(pickers,java.awt.BorderLayout.NORTH);head.add(reviewBox,java.awt.BorderLayout.SOUTH)
        def tracePanel=new javax.swing.JPanel(new java.awt.BorderLayout());tracePanel.setBackground(Color.white);tracePanel.add(head,java.awt.BorderLayout.NORTH)
        // the plot, and under it its caption, rule line and what the green dashed line is (on screen only; Save plots saves the plot)
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
        // left and right arrows step through the listed cells, unless a text box or list is using them, or a mark form is open
        def keys=getRootPane()
        keys.registerKeyboardAction({e->if(previousCell.enabled&&formCell==null)step(-1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_LEFT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        keys.registerKeyboardAction({e->if(nextCell.enabled&&formCell==null)step(1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_RIGHT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        // the mark controls: one button opens the form in its place; in it Tab moves on, Enter saves (or presses the button that
        // has the focus) and Esc closes it, keeping the cell selected; the focus then goes back to the button
        markButton.addActionListener({openForm(cells.selectedItem as String)} as java.awt.event.ActionListener)
        removeMark.addActionListener({removeReviewMark(cells.selectedItem as String)} as java.awt.event.ActionListener)
        saveMark.addActionListener({saveReviewMark()} as java.awt.event.ActionListener);cancelMark.addActionListener({closeForm()} as java.awt.event.ActionListener)
        reasonBox.addActionListener({saveMark.enabled=reasonBox.selectedIndex>=0;noteField.hint=reasonBox.selectedItem=='Other'?'':'Optional'} as java.awt.event.ActionListener)
        markRow.registerKeyboardAction({e->def owner=java.awt.KeyboardFocusManager.currentKeyboardFocusManager.focusOwner
            if(owner instanceof javax.swing.AbstractButton&&owner.parent==markRow)(owner as javax.swing.AbstractButton).doClick() else if(formCell!=null&&saveMark.enabled)saveReviewMark()} as java.awt.event.ActionListener,
            javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ENTER,0),javax.swing.JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT)
        markRow.registerKeyboardAction({e->if(formCell!=null)closeForm()} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ESCAPE,0),javax.swing.JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT)
        flagRow.visible=false;markRow.visible=false;closeForm()
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
    Map protocol(){stopEditing();CalciumResponseProtocol.payload(windowRows(),preset.selectedIndex==2?1:2,fields.protocol_name.text.trim(),validation.selectedItem as String,fieldOverrides,estimator)}
    Map fileSignature(File file,String loaded){if(file==null)return null;String current;try{current=file.isFile()?CalciumInstaller.digest(file.bytes):'missing'}catch(Throwable ignored){current='unreadable'};[path:file.absolutePath,loaded_sha256:loaded,current_sha256:current]}
    String fingerprint(){groovy.json.JsonOutput.toJson([fields:fields.collectEntries{k,v->[k,v.text]},windows:windowRows(),preset:preset.selectedIndex,validation:validation.selectedItem,estimator:estimator,overrides:fieldOverrides,timestamps:fileSignature(timestampsFile,timestampsHash),sensitivity:fileSignature(sensitivityFile,sensitivityHash),diagnostic:fileSignature(diagnosticFile,diagnosticHash)])}
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
                showReview(cells.selectedItem as String)
                progress.indeterminate=false;cancelButton.enabled=false;setup.enabled=true;runButton.enabled=true;fields.values().each{it.enabled=true};fields.sampling_hz.enabled=timestampsFile==null;windowTable.enabled=true
                [confirmed,kcl,preset,validation,loadProtocolButton,saveProtocolButton,loadTimesButton,clearTimesButton,sensitivityButton,diagnosticButton,clearDiagnosticsButton,templateButton,overrideButton,clearOverridesButton].each{it.enabled=true};logButton.enabled=job.logFile.isFile();savePlots.enabled=results!=null;updateDiagnosticAvailability()
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
    void beginAnalysis(){try{Map s=settings();submittedFingerprint=fingerprint();launch(false){task->File made=performAnalysis.call(runtime,measurement.csv,s,measurement.folder,task) as File;checkedRun=[folder:made.name,found:CalciumCallChecks.run(made,withRecording(),task.logFile)];carried=CalciumReview.carry(made,measurement,task.logFile);made}}catch(Throwable t){status.text=t.message?:t.toString()}}
    /** The finish message of the last analysis run here, with the figures of the tally as they are now, so that it follows the
     *  review marks as the Step 4 page's does; then what was carried over from an earlier analysis's review, which is said
     *  until the next run, marks or not (as on the page); then what the checks of the calls could not do (checks.note). null
     *  when the analysis shown was not run here. Writes nothing, and reads no traces: it takes the window's. */
    String finished() {
        Map parts=finish;if(parts==null||results?.name!=parts.folder)return null
        String carry=CalciumReview.note(carried?.folder==results.name?carried:null),checked=checks?.note?' '+checks.note:''
        String sentence='';if(!parts.empty)try{sentence=CalciumAnalysisSettings.percentSentence(CalciumAnalysisSettings.tally(results,measurement,traces))}catch(Throwable ignored){}   // as percentResponding: '' when the tally cannot be made
        (parts.empty?parts.empty as String:(parts.saved as String)+sentence)+carry+checked
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
            stimulusRow(row)&&mark!=null&&CalciumReview.reviewed(marks,cell,status)!=status ? [mark.reviewed,(mark.reason as String)+(mark.note?': '+mark.note:'')] : ['','']}))
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
        cellPicture.calls=callRows.collectEntries{k,v->[k,callOf(k as String)]}
        Map first=callRows.values().find{true}?:[:];plot.stimulusStart=CalciumCsv.number(first.Start_s);plot.stimulusEnd=CalciumCsv.number(first.Effective_end_s?:first.End_s)
        boolean after=reviewMarks&&callRows.values().any{it.Role!='Control'}   // a control window keeps its automatic call
        listing=true;try{int keep=Math.max(0,group.selectedIndex);group.removeAllItems();CalciumCalls.GROUPS.eachWithIndex{String name,int i->group.addItem(name+(after&&(i==1||i==2)?' after review':'')+' ('+(measurement.names as List<String>).count{callRows.isEmpty()||CalciumCalls.inGroup(callOf(it),i)}+')')}
            if(after)group.addItem('Marked ('+(measurement.names as List<String>).count{listedIn(it,MARKED)}+')')
            group.selectedIndex=callRows.isEmpty()||keep>=group.itemCount?0:keep}finally{listing=false}
        listCells();orderHeatmap();plot.changed()
    }
    /** Lists the cells of the chosen group, and draws only those on the picture and behind the mean. */
    void listCells() {
        int which=Math.max(0,group.getSelectedIndex());String keep=cells.getSelectedItem() as String
        List<String> listed=(measurement.names as List<String>).findAll{callRows.isEmpty()||listedIn(it,which)}
        listing=true
        try {cells.removeAllItems();listed.each{cells.addItem(it)};if(listed.contains(keep))cells.setSelectedItem(keep)} finally {listing=false}
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->tile.active=!callRows.isEmpty()&&((i==0&&which==1)||(i==1&&which==0)||(i>=2&&which==3));tile.repaint()}
        cellPicture.group=which;cellPicture.only=which==MARKED?new LinkedHashSet<String>(listed):null
        List<String> drawn=listed.findAll{traces.containsKey(it)}
        plot.all=drawn.collect{traces[it]} as double[][];plot.allCalls=drawn.collect{callOf(it)};plot.allLabel=which==0?'all cells':which==MARKED?'marked':CalciumCalls.GROUPS[which].toLowerCase()
        previousCell.enabled=listed.size()>1;nextCell.enabled=listed.size()>1
        if(cells.getItemCount()>0)selectCell(cells.getSelectedItem() as String,true)
        else {plot.cell='';plot.selected=new double[0];plot.responseCall=[:];cellStatus.setText('No cells in this group.');showReview(null);cellPicture.selected=-1;plot.changed();cellPicture.repaint()}
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
        Map review=CalciumReview.read(dir,measurement);reviewMarks=review.marks as Map;reviewProblem=(review.problem?:'') as String   // read again by summarize, which also shows the tile after review
        closeForm();automatic=CalciumReview.scored(dir)
        // the checks of this analysis: what its run here gave (as the Step 4 page takes them), else the table saved beside it
        Map run=checkedRun;checks=CalciumCallChecks.answer(run?.folder==dir.name?run.found as Map:CalciumCallChecks.saved(dir,measurement))
        results=dir;table(metrics,new File(dir,'metrics.csv'));table(qc,new File(dir,'cell_quality.csv'));table(registry,new File(dir,'roi_registry.csv'));responseRows=table(calls,new File(dir,'responder_calls.csv'));table(counts,new File(dir,'responder_counts.csv'));table(windowMetrics,new File(dir,'window_metrics.csv'));table(replicates,new File(dir,'replicate_counts.csv'));table(sensitivity,new File(dir,'sensitivity_counts.csv'));table(diagnostics,new File(dir,'negative_control_counts.csv'))
        callsSaved=calls.model;countsSaved=counts.model   // as saved: the columns of the review marks are added to copies (reviewTables)
        cellStatus.lines=3   // the call line with a mark and the longest note takes three lines at the window's usual width
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
        if(row==null) cellStatus.setText(callRows.isEmpty()?' ':'<html><b>'+name+'</b> · no response call was saved for this cell</html>')
        else {
            def number={String key,String format->double v=CalciumCsv.number(row[key]);Double.isFinite(v)?String.format(Locale.US,format,v):null}
            def esc={Object text->CalciumOutlineProgress.escape(text==null?'':''+text)}
            def dot={String shown->"<font color='"+CalciumTheme.hex(CalciumCalls.of(shown))+"'>●</font> "}
            Map mark=row.Role=='Control'?null:reviewMarks[name] as Map;String after=callOf(name);boolean marked=mark!=null&&after!=row.Status
            List parts=marked ? [dot(after)+'<b>'+name+': '+esc(after)+' after review</b>',esc((mark.reason as String).toLowerCase()+(mark.note?': '+mark.note:''))]+(mark.carriedFrom ? ['from an earlier analysis'] : [])+['automatic call '+esc(row.Status)]
                : [dot(row.Status)+'<b>'+name+': '+esc(row.Status?:'')+'</b>']
            parts<<esc((row.Reason?:'').toLowerCase())
            String rise=number('Peak_increase','%.2f'),at=number('Peak_s','%.0f'),need=number('Threshold_increase','%.2f')
            if(rise!=null)parts<<'largest rise '+rise+(at!=null?' at '+at+' s':'')+(need!=null?' (needed '+need+')':'')
            cellStatus.setText('<html>'+parts.findAll{it}.join(' · ')+'</html>')
            String tip='Rise and threshold are fractions of the fluorescence just before the stimulus, as in the Cell response calls table.'
            // when the mark was made, to the minute as the Step 4 page says it, and in which analysis when it was carried over
            String when=marked&&mark.time?CalciumReview.minute(mark.time as String):''
            cellStatus.setToolTipText(marked?'<html>'+esc('Marked'+(when?' '+when:'')+(mark.carriedFrom?' in '+mark.carriedFrom:''))+'<br>'+esc(tip)+'</html>':tip)
        }
        showReview(name)
    }
    /** Under the selected cell's call line: the flags of its checks, each in the checks' own words, and its mark controls. A flag
     *  with the first frame of its event is a link to that frame. A control window keeps its automatic call: it has neither. */
    void showReview(String name) {
        Map row=name==null?null:callRows[name];boolean primary=row!=null&&row.Role!='Control'
        Map cell=primary?((checks?.cells?:[:]) as Map)[name] as Map:null
        String auto=primary?automatic[name] as String:null;boolean evaluable=auto in CalciumReview.CALLS
        def owner=this;flagRow.removeAll()
        ((cell?.flags?:[]) as List).each{Map flag->
            if(flagRow.componentCount>0)flagRow.add(CalciumTheme.label('·',CalciumTheme.MUTED))
            def label=CalciumTheme.label(flag.text as String,cell.firstFrame!=null?CalciumTheme.BLUE:CalciumTheme.MUTED)
            String detail=cell.detail?:''
            if(cell.firstFrame!=null) {
                int first=cell.firstFrame as int
                label.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR))
                label.setToolTipText('<html>Go to frame '+first+(detail?'<br>'+CalciumOutlineProgress.escape(detail):'')+'</html>')
                label.addMouseListener(new java.awt.event.MouseAdapter(){void mouseClicked(java.awt.event.MouseEvent e){if(javax.swing.SwingUtilities.isLeftMouseButton(e))owner.goToFrame(first)}})
            } else if(detail)label.setToolTipText(detail)
            flagRow.add(label)
        }
        // an evaluable cell that could not be checked, in a recording whose surroundings were measured (else checks.note says it once)
        if(cell!=null&&evaluable&&checks.ringMeasured&&!cell.checked&&cell.why){if(flagRow.componentCount>0)flagRow.add(CalciumTheme.label('·',CalciumTheme.MUTED));flagRow.add(CalciumTheme.label('passing object not checked: '+cell.why,CalciumTheme.MUTED))}
        // the row keeps the height of a flag when there is none; added after the flags, so the dots above count only flags
        flagRow.add(javax.swing.Box.createRigidArea(new java.awt.Dimension(0,(int)CalciumTheme.label('·').preferredSize.height)))
        flagRow.visible=results!=null
        // one button: mark the cell the other way, or remove its mark; none while an analysis runs or while the review file
        // is not applied, which is said under the controls of every cell of the primary stimulus instead
        Map can=markable(name)
        if(formCell!=null&&(formCell!=name||!can.can||can.marked))closeForm()
        String blocked=reviewProblem?CalciumReview.blocked(reviewProblem):null
        if(problemCell!=name||(blocked&&problemSaid==blocked))say(null,null)
        if(blocked&&primary&&problemCell==null)say(blocked,name)
        markButton.text=auto=='Responded'?'Mark as not a responder':'Mark as a responder'
        markButton.visible=can.can&&!can.marked&&formCell==null;removeMark.visible=can.marked;markRow.visible=results!=null
        reviewBox.revalidate();reviewBox.repaint()
    }
    /** What the mark controls can do for a cell: can, when it has a call to review on the primary stimulus, no analysis is
     *  running and the review file is applied; marked, when it has a mark then. */
    Map markable(String name) {
        Map row=name==null?null:callRows[name];String auto=row!=null&&row.Role!='Control'?automatic[name] as String:null
        boolean can=!busy&&results!=null&&auto in CalciumReview.CALLS&&!reviewProblem
        [can:can,marked:can&&(reviewMarks?.get(name) as Map)?.automatic==auto]
    }
    /** A flag's link: the Traces tab, the recording rather than the still colour picture, at the first frame of the event. */
    void goToFrame(int first){tabs.selectedIndex=1;if(movie!=null&&view.selectedIndex==0)view.selectedIndex=2;stopPlaying();showFrame(first)}
    /** A refusal under the mark controls, as CalciumReview words it; null takes it away. While results are shown its line is
     *  kept, empty. */
    void say(String text,String cell){problemCell=text==null?null:cell;problemSaid=text;markProblem.say(text==null?'&nbsp;':CalciumOutlineProgress.escape(text));markProblem.visible=text!=null||results!=null;reviewBox.revalidate();reviewBox.repaint()}
    /** The form that sets a mark, in place of the button: the reason (chosen already for a responder with a flag), an optional
     *  note, Save and Cancel. The reviewed call is the other one of the two, so it is not asked for. */
    void openForm(String name) {
        Map can=markable(name);if(!can.can||can.marked)return
        Map cell=((checks?.cells?:[:]) as Map)[name] as Map
        String chosen=automatic[name]=='Responded'?((cell?.flags?:[]) as List).collect{REASON_OF_RULE[it.rule]}.find{it}:null
        formCell=name;reasonBox.selectedIndex=chosen?CalciumReview.REASONS.indexOf(chosen):-1;noteField.text=''
        [reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{it.visible=true};markButton.visible=false;say(null,null)
        takeFocus.call(chosen?noteField:reasonBox)
    }
    /** Closes the form; the button comes back in its place, and the focus with it when the focus was in the form (back). */
    void closeForm(boolean back=true){boolean had=back&&formFocused();formCell=null;[reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{it.visible=false}
        Map can=markable(cells.selectedItem as String);markButton.visible=can.can&&!can.marked;say(null,null);refocus(had)}
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
     *  the plot and the calls and counts tabs follow the marks, the cell stays selected, the time cursor stays where it is,
     *  and the finish message follows. */
    void afterReview(String name) {
        Map review=CalciumReview.read(results,measurement);reviewMarks=review.marks as Map;reviewProblem=(review.problem?:'') as String
        try{summarize(results)}catch(Throwable ignored){}
        reviewTables()
        selectWindow()
        if((cells.getModel() as javax.swing.DefaultComboBoxModel).getIndexOf(name)<0)showGroup(0)
        cells.setSelectedItem(name)
        String done=finished();if(done!=null)status.text=done
    }
    void exportPlots(){
        if(results==null)return
        changed()
        // the saved trace has no time cursor or its reading: only the axes, the traces and the key
        double at=plot.moment;plot.moment=Double.NaN
        try{
            ['trace':plot,'heatmap':heatmap].each{name,panel->int w=Math.max(650,panel.width),h=Math.max(260,panel.height);def image=new java.awt.image.BufferedImage(w,h,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics();try{panel.paint(g)}finally{g.dispose()};javax.imageio.ImageIO.write(image,'png',new File(results,name+'.png'))}
            Map provenance=[result_folder:results.name,cell:cells.selectedItem,stimulus:stimuli.selectedItem,protocol_hash:appliedManifest.protocol_hash,manifest_sha256:CalciumInstaller.digest(new File(results,'analysis_manifest.csv').bytes),settings_pending:fingerprint()!=appliedFingerprint,source:'Completed analysis snapshot; display is unsmoothed']
            // the groups, colours and order follow the review marks on the primary stimulus: say so, and which file they come
            // from; a control window keeps its automatic calls, so its plots say nothing of the review
            File review=CalciumReview.file(measurement,results)
            if(reviewMarks&&review?.isFile()&&callRows.values().any{it.Role!='Control'})provenance+=[calls:'after review',review_file:review.name,review_sha256:CalciumInstaller.digest(review.bytes)]
            new File(results,'plot_provenance.json').setText(CalciumResponseProtocol.json(provenance)+'\n','UTF-8')
            status.text='Saved completed-run plots and plot_provenance.json in '+results.name
        }catch(Throwable t){status.text='Could not save plots: '+t.message}
        finally{plot.moment=at;plot.repaint()}
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
    static String fingerprint(Map draft) {
        def sorted;sorted={v->v instanceof Map?new TreeMap(v.collectEntries{k,x->[''+k,sorted(x)]}):v instanceof Collection?v.collect{sorted(it)}:v}
        CalciumInstaller.digest(groovy.json.JsonOutput.toJson(sorted(draft)).getBytes('UTF-8'))
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
 *  The span was set as 10 frames at 1 frame per second (SPAN_S seconds). */
@CompileStatic
class CalciumPassingObject {
    static final int AVERAGE = 3, DIP = 3, APART = 2
    static final double SPAN_S = 10d, HALF = 0.5d, CELL_RISE = 1d, NOISE = 5d, CORE_SHARE = 0.25d
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
    /** The largest value among the frames from..to where the cell is below "half": [value, frame], or null when there is no such frame. */
    private static double[] highest(double[] zone, double[] cell, double half, int from, int to) {
        int at = -1
        for (int f = from; f <= to; f++) if (cell[f] < half && (at < 0 || zone[f] > zone[at])) at = f
        return at < 0 ? null : [zone[at], at] as double[]
    }
    /** @param raw the cell's raw gray values, @param background its background, both one per frame
     *  @param zone the 25 zone traces (CalciumRingAndCore.NAMES); null where a zone is not usable; the whole array null when
     *         nothing was measured around the cell, which still gives the event
     *  @param firstF0, lastF0 the F0 frames and @param peak the frame of the call's peak, all counted from 0
     *  @param span frames searched before and after the event
     *  @return null when the cell does not rise at the peak; else first and last (the event's frames, counted from 1), rule
     *          ('passes through', 'on the edge' or ''), text, detail, and where the rule is 'on the edge' edgeFirst, edgeLast */
    static Map<String, Object> check(double[] raw, double[] background, float[][] zone, int firstF0, int lastF0, int peak, int span) {
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
        int a = peak, b = peak
        while (a > 0 && cell[a - 1] >= half) a--
        while (b < n - 1 && cell[b + 1] >= half) b++
        boolean moved = true
        while (moved) {
            moved = false
            for (int d = 2; d <= DIP + 1; d++) if (a - d >= 0 && cell[a - d] >= half) { a -= d; moved = true; while (a > 0 && cell[a - 1] >= half) a--; break }
            for (int d = 2; d <= DIP + 1; d++) if (b + d < n && cell[b + d] >= half) { b += d; moved = true; while (b < n - 1 && cell[b + 1] >= half) b++; break }
        }
        Map<String, Object> out = new LinkedHashMap<String, Object>()
        out.put('first', a + 1); out.put('last', b + 1); out.put('rule', ''); out.put('text', ''); out.put('detail', '')
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
            if (n > 1) {
                double[] d = new double[n - 1]
                for (int f = 1; f < n; f++) d[f - 1] = plain[f] - plain[f - 1]
                double m = CalciumMotion.median(d)
                for (int f = 0; f < n - 1; f++) d[f] = Math.abs(d[f] - m)
                steps = 1.4826d * CalciumMotion.median(d) / Math.sqrt(2d) / Math.sqrt((double) AVERAGE)
            }
            noise[z] = Double.isNaN(sd) ? steps : Double.isNaN(steps) ? sd : Math.max(sd, steps)
        }
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
    static final List HEADER=['Cell','Call','Flag','Rule','Event_first_frame','Event_last_frame','Detail','Near_ring_sectors_usable','Area_pixels','Ring_measured','Field_movement_measured','Passing_object_check','Not_checked_because','Crossings','Longest_crossing_samples']
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
    /** The saved files, when they belong to these cells, these outlines and this raw table. @return values, counts, area and near, or why not */
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
            return [values:CalciumRingAndCore.read(values,names.collect{''+it},frames),counts:body.collect{List row->(3..<3+CalciumRingAndCore.ZONES).collect{row[it] as int}},area:body.collect{it[2] as int},near:body.collect{it[nearAt] as int}]
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
     *  @return values, counts, area and near, or why not */
    static Map around(Map measurement,int frames,File log=null) {
        if(!(measurement?.folder instanceof File))return [why:'this measurement has no folder']
        String prefix=new File(measurement.folder as File,measurement.recording+'_').path
        List names=((measurement.names?:[]) as List).collect{''+it}
        File raw=new File(prefix+'gray_values_raw.csv')
        if(!raw.isFile())return [why:'the raw gray values are not beside the analysis']
        Roi[] rois=measurement.rois instanceof Roi[]&&(measurement.rois as Roi[]).length==names.size()?measurement.rois as Roi[]:null
        Map saved=savedAround(prefix,names,rois,CalciumInstaller.digest(raw.bytes),frames)
        if(saved.values!=null)return saved
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
        def calls=CalciumCsv.read(new File(dir,'responder_calls.csv'));List h=calls[0]
        int cellAt=h.indexOf('Cell'),roleAt=h.indexOf('Role'),statusAt=h.indexOf('Status')
        if(cellAt<0||statusAt<0)throw new CalciumUnchecked('the calls table has no Cell or Status column')
        Map byCell=[:];calls.drop(1).each{List row->if((roleAt<0||row[roleAt]=='Stimulus')&&!byCell.containsKey(''+row[cellAt]))byCell[''+row[cellAt]]=row}
        def number={List row,String name->int i=h.indexOf(name);i<0||row==null?Double.NaN:CalciumCsv.number(''+row[i])}
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
        Map ring=measured.why ? [why:measured.why] : around(measurement,frames,log)
        double[][] shifts=ring.values!=null?movement(new File(measurement.folder as File,measurement.recording+'_').path,frames):null
        String unchecked=ring.values!=null?'':(''+(ring.why?:'ring and core were not measured'))   // why no cell is checked for a passing object
        List out=[]
        names.eachWithIndex{String name,int c->
            List call=byCell[name] as List;String status=call==null?'':''+call[statusAt]
            boolean scored=status=='Responded'||status=='Did not respond'
            Map row=[Cell:name,Call:status,Flag:'',Rule:'',Event_first_frame:'',Event_last_frame:'',Detail:'',Near_ring_sectors_usable:'',Area_pixels:'',Ring_measured:ring.values!=null?'yes':'no',
                     Field_movement_measured:shifts!=null?'yes':'no',Passing_object_check:'not checked',Not_checked_because:'',Crossings:'',Longest_crossing_samples:'']
            if(ring.values!=null){row.Near_ring_sectors_usable=ring.near[c];row.Area_pixels=ring.area[c]}
            List flags=[]   // a flag is [text, rule, detail]
            String why=''
            double peakAt=number(call,'Peak_s')
            if(!scored)why='no call'
            else if(!Double.isFinite(peakAt))why='no peak time'
            else if(measured.why)why=unchecked
            else {
                double[] raw=measured.raw[name] as double[],background=measured.background[name] as double[]
                if(!CalciumPassingObject.finite(raw)||!CalciumPassingObject.finite(background))why='missing gray values'
                else {
                    int peak=0;for(int f=1;f<frames;f++)if(Math.abs(time[f]-peakAt)<Math.abs(time[peak]-peakAt))peak=f
                    float[][] zone=null;int outside=0
                    if(ring.values!=null){zone=new float[CalciumRingAndCore.ZONES][];for(int z=0;z<CalciumRingAndCore.ZONES;z++)if((ring.counts[c][z] as int)>=CalciumRingAndCore.MIN_PIXELS){zone[z]=(ring.values as float[][])[c*CalciumRingAndCore.ZONES+z];if(z<16)outside++}}
                    Map event=CalciumPassingObject.check(raw,background,zone,first,last,peak,span)
                    if(event!=null){row.Event_first_frame=event.first;row.Event_last_frame=event.last}
                    if(zone==null)why=unchecked
                    else if(event==null)why='no rise at the peak'
                    else if(outside==0)why='no ring around the cell'   // no ring sector has enough pixels: neither rule can be tested
                    else {
                        row.Passing_object_check='checked'
                        if(event.rule) {
                            Map moved=shifts==null?null:CalciumMotion.movedNear(shifts,first+1,last+1,(event.first as int)-span,(event.last as int)+span)
                            if(moved!=null)flags<<['field moved, '+CalciumPassingObject.frames(moved.first as int,moved.last as int),'field moved',movedText(moved)+'; otherwise '+event.text+': '+event.detail]
                            else flags<<[event.text,event.rule,event.detail]
                        }
                    }
                }
            }
            // separate crossings of the rise the cell needs, from the engine's own columns
            if(scored) {
                double need=number(call,'Threshold_increase'),base=number(call,'Baseline_mean'),reference=number(call,'Local_reference_relative')
                if(!Double.isFinite(reference))reference=number(call,'Local_F0_relative')
                if(!Double.isFinite(reference))reference=1d+base
                double start=number(call,'Start_s'),end=number(call,'Effective_end_s');if(!Double.isFinite(end))end=number(call,'End_s')
                if(traceAt[name]!=null&&[need,base,reference,start,end].every{Double.isFinite(it as double)}&&need>0d&&reference>0d) {
                    int at=traceAt[name] as int;double[] y=samples.collect{CalciumCsv.number(''+it[at])} as double[]
                    int[] crossed=CalciumPassingObject.crossings(time,y,base,reference,need,start,end,number(call,'Max_gap_s'))
                    row.Crossings=crossed[0];row.Longest_crossing_samples=crossed[1]
                    if(status=='Responded'&&crossed[0]>=CROSSINGS&&crossed[1]<=SHORT)flags<<['short separate crossings','short separate crossings',crossed[0]+' crossings, the longest '+crossed[1]+' samples']
                }
            }
            row.Flag=flags.collect{it[0]}.join('; ');row.Rule=flags.collect{it[1]}.join('; ');row.Detail=flags.collect{it[2]}.join('; ');row.Not_checked_because=why
            out<<row
        }
        [rows:out,ringMeasured:ring.values!=null,note:note(ring.values!=null,shifts!=null,unchecked)]
    }
    /** The note for the whole analysis: '' when every evaluable cell could be checked as written. */
    static String note(boolean ring,boolean movement,String unchecked){!ring?'Not checked for a passing object'+(unchecked?': '+unchecked:'')+'.':movement?'':'The movement of the field was not saved, so a moving field is not told from an object.'}
    /** The checks as the Step 4 page is given them: for every cell the flag's text, its rule and the first frame of the event.
     *  A cell with two flags has both in text and in rule, joined with "; ", and one by one in flags. */
    static Map answer(Map found) {
        if(found==null)return null
        def whole={v->String s=v==null?'':''+v;s==~/-?\d+/?(s as int):null}
        def each={Map row->String text=''+(row.Flag?:''),rule=''+(row.Rule?:'');text.isEmpty() ? [] : [text.split('; ') as List,rule.split('; ') as List].transpose().collect{[text:it[0],rule:it[1]]}}
        [file:found.file,ringMeasured:found.ringMeasured==true,note:found.note?:'',cells:((found.rows?:[]) as List).collectEntries{Map row->
            [''+row.Cell,[text:''+(row.Flag?:''),rule:''+(row.Rule?:''),flags:each(row),firstFrame:whole(row.Event_first_frame),lastFrame:whole(row.Event_last_frame),detail:''+(row.Detail?:''),
                checked:row.Passing_object_check=='checked',why:''+(row.Not_checked_because?:''),ringMeasured:row.Ring_measured=='yes',movementMeasured:row.Field_movement_measured=='yes',
                nearSectors:whole(row.Near_ring_sectors_usable),area:whole(row.Area_pixels),crossings:whole(row.Crossings),longestCrossing:whole(row.Longest_crossing_samples)]]}]
    }
    /** The checks saved for an analysis, with the note rebuilt from the table; null when there are none. */
    static Map saved(File dir,Map measurement) {
        try {
            File source=file(measurement,dir);if(source==null||!source.isFile())return null
            def rows=CalciumCsv.read(source);if(rows.isEmpty()||rows[0]!=HEADER)return null
            List body=rows.drop(1).collect{List line->[HEADER,line].transpose().collectEntries{it}}
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
    static final List HEADER=['Cell','Automatic_call','Reviewed_call','Reason','Note','Time','Analysis','Responder_calls_sha256','ROI_sha256','Carried_from']
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
     *  holds no marks. @return file (its name, or null when there is none), marks (cell to automatic, reviewed, reason, note,
     *  time, roi, carriedFrom) and problem: '' or why the file is not applied, in which case marks is empty */
    static Map read(File dir,Map measurement) {
        File source=file(measurement,dir);Map none=[file:null,marks:[:],problem:'']
        if(source==null||!source.isFile())return none
        String unused='The review file '+source.name
        try {
            String text=source.getText('UTF-8');if(text.startsWith('﻿'))text=text.substring(1)
            text=text.replaceAll(/[\r\n]+\s*$/,'\n')
            if(text.trim().isEmpty())return [file:source.name,marks:[:],problem:'']
            def rows=CalciumCsv.parse(text)
            if(rows.isEmpty()||rows[0]!=HEADER)return none+[file:source.name,problem:unused+' cannot be read, so its marks are not used.']
            Map calls=scored(dir);String hash=digest(dir);Map marks=new LinkedHashMap()
            for(List line in rows.drop(1)) {
                Map row=[HEADER,line].transpose().collectEntries{it}
                if(row.Analysis!=dir.name)return none+[file:source.name,problem:unused+' is for another analysis, so its marks are not used.']
                if(row.Responder_calls_sha256!=hash||marks.containsKey(row.Cell)||row.Automatic_call!=calls[row.Cell]||refused(row.Cell as String,row.Reviewed_call,row.Reason,row.Note,calls,row.Cell as String)!=null)
                    return none+[file:source.name,problem:unused+' does not match the calls of this analysis, so its marks are not used.']
                marks[row.Cell]=[automatic:row.Automatic_call,reviewed:row.Reviewed_call,reason:row.Reason,note:row.Note,time:row.Time,roi:row.ROI_sha256,carriedFrom:row.Carried_from]
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
                CalciumCsv.write(reviewPart,[HEADER]+(marks.keySet() as List).sort{cells.indexOf(it)}.collect{String cell->Map m=marks[cell] as Map;[cell,m.automatic,m.reviewed,m.reason,m.note?:'',m.time,dir.name,hash,roi[cell]?:'',m.carriedFrom?:'']})
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
            marks[cell]=[automatic:calls[cell],reviewed:body.call,reason:body.reason,note:(note?:'') as String,time:now(),carriedFrom:'']
        }
        write(dir,measurement,marks,traces)
    }
    /** "3 of 143 cells were reclassified after manual review: 2 from responder to non-responder (1 noise, 1 passing object) and
     *  1 from non-responder to responder (other)." From the marks that changed a call of all cells; '' when none did. */
    static String methods(Map tally,Map review) {
        List changed=(tally?.changed?:[]) as List;if(changed.isEmpty())return ''
        int total=(((tally.rows as List).find{it.key=='all'} as Map).counted as List)[1] as int
        List marks=changed.collect{(review.marks as Map)[it] as Map}
        def part={List group,String from,String to->
            Map count=group.countBy{it.reason};List reasons=(count.keySet() as List).sort{a,b->(count[b] as int)<=>(count[a] as int)?:REASONS.indexOf(a)<=>REASONS.indexOf(b)}
            group.size()+' from '+from+' to '+to+' ('+(reasons.size()==1?(reasons[0] as String).toLowerCase():reasons.collect{count[it]+' '+(it as String).toLowerCase()}.join(', '))+')'}
        List down=marks.findAll{it.automatic=='Responded'},up=marks.findAll{it.automatic!='Responded'},parts=[]
        if(down)parts<<part(down,'responder','non-responder');if(up)parts<<part(up,'non-responder','responder')
        changed.size()+' of '+total+' cells '+(changed.size()==1?'was':'were')+' reclassified after manual review: '+parts.join(' and ')+'.'
    }
    /** What carry gave for this analysis, until the next run, also once every mark it carried over is removed: null for
     *  another analysis. @param marks the marks of the analysis now, which do not change it */
    static Map current(Map carried,File dir,Map marks){carried!=null&&carried.folder==dir.name?carried.findAll{it.key!='folder'}:null}
    /** What the Step 4 page is given about the review of the analysis shown; problem is the whole sentence the page shows, with
     *  what to do (blocked), or ''. @param carried what carry gave when this analysis was run */
    static Map state(File dir,Map measurement,Map review,Map tally,Map carried=null) {
        List names=order(dir,measurement);Map marks=(review.marks?:[:]) as Map
        File source=review.file?file(measurement,dir):null
        [folder:dir.name,file:review.file,sha256:source?.isFile()?CalciumInstaller.digest(source.bytes):null,problem:review.problem?blocked(review.problem as String):'',
         marks:marks.collectEntries{cell,m->[cell,[name:shown(names,cell as String),automatic:m.automatic,reviewed:m.reviewed,reason:m.reason,note:m.note?:'',time:m.time,carriedFrom:m.carriedFrom?:'']]},
         after:scored(dir).collectEntries{cell,call->[cell,reviewed(marks,cell as String,call as String)]},
         changed:((tally?.changed?:[]) as List).size(),methods:methods(tally,review),carried:current(carried,dir,marks),
         reasons:REASONS,calls:CALLS,noteMax:NOTE_MAX]
    }
    /** When "yyyyMMdd-HHmmss" an analysis was made: from its folder's name, else from when the file beside the tables was written. */
    static String made(String analysis,File source){def m=analysis=~/^SimpleCa-(\d{8}-\d{6})-/;m.find()?m.group(1):new java.text.SimpleDateFormat('yyyyMMdd-HHmmss').format(new Date(source.lastModified()))}
    /** At the end of a new analysis: the marks of the most recent earlier analysis of this recording made by this version (it
     *  has a checks file or a review file beside the tables) are kept for every cell whose automatic call and outline are as
     *  they were, and written for the new analysis; the others are dropped. Nothing is carried when that analysis has no marks,
     *  so marks removed there, or dropped there, do not come back from an older one.
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
            Map marks=got.marks as Map;if(marks.isEmpty())return null
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
                else {checking=true;stage='Checking the calls'
                    Map checked=CalciumCallChecks.run(output,withRecording(),mine.logFile);carried=CalciumReview.carry(output,measurement,mine.logFile);load(output,fingerprint,submitted,checked);percent=100
                    String empty=CalciumAnalysisSettings.noEvaluableCells(output)
                    if(!empty)try{CalciumAnalysisSettings.savePercent(output,measurement,resultsMap?.tally as Map)}catch(Throwable problem){CalciumCallChecks.logProblem(mine.logFile,'saving the percent table',problem)}
                    synchronized(this){finish=[folder:output.name,empty:empty,saved:'Complete. Protocol, inputs and tables saved in '+output.name+'.',noise:CalciumResponseProtocol.noiseNote(submitted?.protocol)];message=finished();tone=empty?'warn':'good'}}
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
        Map submitted=CalciumAnalysisSettings.copy(incoming),uploads=[:]
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
     *  marks; the part about carried marks stays until the next run. */
    synchronized String finished() {
        Map parts=finish;if(parts==null||results?.name!=parts.folder)return message
        String carry=CalciumReview.note(CalciumReview.current(carried,results,(resultsMap?.review?.marks?:[:]) as Map))
        parts.empty?(parts.empty as String)+carry:(parts.saved as String)+CalciumAnalysisSettings.percentSentence(resultsMap?.tally as Map)+(parts.noise?:'')+carry
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
        synchronized(this){meanTrace=mean as double[];restFirst=first;restLast=last;rest=null;frameEpoch++;appliedTime=time as double[];appliedTimeSource=loaded.timeSource;appliedFingerprint=applied;results=dir;resultsMap=loaded;resultsJson=text;version++}
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
        if(finish?.folder==results.name&&!busy)message=finished()
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
    static final String SHA256 = "4c73e534ef6e1e92273f83bfc51feb329f8b481aa6f86b07d50382717207168a"
    static final List<String> PARTS = [
        "UEsDBBQAAAAIAAAAIVzh/p+pMS4AAJa9AAAJAAAAc3RlcDQuY3NzzX1LjxvXmthev+JEDcNNm6T5bjYb17iS3MJcJL43sD2ZhWEERfKQrKtiFVNV7Fbb0CK7rANkVvkVswqCIIv8FP+SfI/zrlMkW5IxY8OWWHXqPL/363zzlfixlgcxEYdkK/viTZEVx7LqivrpIEW1SjLZFdUhWaX5ViT5WlS7ZF08VmJTZFnxKOqdFD+m+0Mm3yT/719+/6//WySHw5cVPt9Lcf3DN/SX/g+d/gshxI+75CAXQj7I8kksi/fiepWUaxgsxWH2Mj92xfJY10XeFVla1V2xSWW27ohdQl2KVVHmshS9Xpms02N1Bw/TXCyTsqLJrWj2onpM6tVOVvDZg8RxdfsetLwTeYFfbUUKH4lDmmWiKOFvq7RcZbAFP6a/2jk6k6ERaD74Za+3KvK6LLLeTuzS7U48pvVODMeH96KW72ta7vewICGhVSpxR5MlzzJL83eVeCzTupa5gJXxKPCylLTKIs+ehHy/koc6pee478mTqArxuJPc6Ilar+UmzeW6L7765sWiLIpa/PYCl7vMjrCGq8HgZvR6dKce9DKYZw2P7+dvJ2/n5vE6Kd9R48l3N6/u6PvHXVpjB2/pH2y5LZOn3nSAz27e3r59bZ4NB/jwu/v7yf0r83BGD2fTm+HNd+bhLT0cvh3djsc8zLYo1vjo/uZ+PMF2j0mZw4PbV9PXgwFNMMEGryej8XDO35SyOhT5WtKHNzeT2RTbwZnid6PpaDTG33DEsoQn81e3k1cD/pJBtzdciIEYwjnhWZXbZXI9HnYnw+502B30h4POnW06wqZTaDacxNqOOmpKBFwLMTu8v/NgbSFG/MjACjyZwyP6bAMPF6IHCJPJXvVU1XLfFa8ROr5PVj/S77fQpCte/ii3hRT/+JeXXfFDsSzqAp79g8weZJ2uEvFXeZTw5lWZJhkga5JXvUqW6YZH2Rd5sYAu3orv4W/QDv5IVtjDmwTwe50mgPRr7ODlG0CeFNDrr/IRfuKHiPny7sWHF1+J3xBhe1X6K6DOAv5ermXZg0d34sOLXb2HkZfF+gma7SSDGQDGF/hSPd4n5TaFw4Vj5YXjnn4z7E+m4iEpr3k7YPMRh2Hj+JmGG3i+TFbvtmVxzNfeyykd2KNcvktr6qJX7QtG8AUgTg17kiaVXN+JdVodsuRpITaZhFnj/3vrtJQrRLIFjnvc5zjhn3fpei3zX2DS5pu8yKX4d0DoirKGXrFZH2bzCG2wI1gN/Jsc6wLfvOjXxcH9mkeEmWzzHiDWHmBlBXRBAjHaJgf4mCDnkKzXNO3hAKBtNMBnzVUjytptutpsNnd8MgSyupWG9c6dmh8ugGa9LJGc/CZ4r4jWDREk+cGjOrwZoh9RgR7BAHbwWCYHpwt4nut+9GcT/KxAdlHDsvvzO2+YCQ7DcNDL5KZWGAM9wikAQCFRVh1ukn2aPenVICh2/L4IsexIt0DlCiDZmwz3gA/wjihxzz6WWZYeqrRqWdkepvWYrusdASnMinjDb81hzUEhBRnyOTE+aErA0+ZfAewSFRlNp139X38477TvNc6hj3QSEdDp5mp9uxnJYYgv0BA6+yD4OySn4XcAMMP1LPgOG3YYqBHQw0/WciWT4BOgzJ0WmNW0AbjafkGktiqydA0dTVez1cCMswCqfd0v8k6IaYxE8LM41hEMC/AqikIzMx08Vk2TaOykStfSdDuAf8dTankhheDxGAAbIPfhxT5J88ikfegKJvX8kXGDUHgKzkodJwJTR5+EewQO5RwyWT0NuKcIi93skWKRdLI4q/4mO1Y7mJtpM4ht1oXLPrl1wH1GAZJO4+TMcCD4d64JD5HwSwj1349VnW6eiJVL5F+Eq72lrB+l1GcztxSOiJsiJDjKt0QwF5u0rOreapdmeHRtzK6VGiRp+RzGQjuB1AC/S8021U8Z9VvukyzGcWc0iQ8v0vxwrAlLf0Z94JdOV9Cjn+n/exAbftFP8P2f8uN+KUv/WSWTcrX7RY29gFcgmaW1T83HOE91xDMCJH3CPCsjP7lgNxA3ltx4dObN9M39d/PnUGXmo2onzBTVDvQfmWCoCQ7HA4Ye3p5NsTqCeF/JDOBW/2KRnn/1HtIqXWbYA9AzkO8kyYVqssQPBt3hcNId3sy7IA/B1FS7XrHZVLKmtdnZrHZJvpWE+bw+ny6zfBAXHFgFsGf7cwJCYy/NHwCE1n+qyyMeKK8k+goZKugdSeRl23SITTQ4kNwQPNP3CWJ7iAz8XYP17pP3PUfAdGlwD/BBS2Bux8fMlz4V/CgJZHijKIHzBfFNFkZauCRol6t3ZxCRqFkF21UrZKSRVqBdY4+HImU8Bezk3uhAzFR7IEKSgAG9rhCno6eMUylyoCl7j9YKxSL6yYo1yMZMLbWiaSKBWQhDZlCpjAk+bZTC2V1YDn4NFGftTmmOHEIzLq1VM6A1KUMTDfkLA2JxZG+g82n0d/nISSrDEl6DXBGCGljsjyNnW4OgXKXMzOz0+CRhOwb7SkhQTroe5jgv7NIXOwRzIsULOMoEyMm6E+X9HpqfIBCmZ92dhXbQv0EBH57WuWJdO3KF3ou13CTHzDnE/qFMAV6eWid/etqe3tPk8eEoz9w3NIWcGF+/b4zi7KEHhq9vvhvfj8MOzWNvMabT7a6o6qAngqQDkF7UPf3evFfxk/7mK/GKzE6u1SlRTEq8k/LA5jWGZlFsRFpXBOBkWVLzKiWKEZq39f9+3B+6gO1Zsu1t0gwAHn5VyYNEY4wUBmUDec1BrEmAVxbhBWnQwX4EKq+7B02lasLUMeyG4SHY3IhKNm1AgT4kQ7Q+Wu64iOjMojSnSeOQISq5ZDQbaM3gUBbbUlaVY46ZnVFT0V7VbmZhpI6pO3asb4HJPDQMQFZwb0V3l0xSczEeuCTQDNFfHqsnM47qeTz9wrPMJDngJHdWZSi4Dfsj7gsElh5qlClwzBy4NPb953fyaVMme1mp1r+JwReWC7OQ0IMxECpxReE7ZebCaQKzr4+VL2240D/q38o9QTe31Hp9iwqvGi2TRhsWjmwTpefHhRXE/p8AuWEDgRBXffpRoZaQrySbrRH3D7JEMSPZSrK2JwKVFkUJise8L/4xB9CBn10BaI1WAqQPZJtOskwktTLQEzRrGEArOD8+5jXawHeSbenQMRIOFGNXMsvYxl3bWaJ/AYbjyfLXZPDepNtjidZ0oDHabo7WzV3x+BLN8jQYdCiYIhH56qN7ISIDNSSfdj1qbgw9gTUNuz4jBy4B8nCrXCnQ12g/XU1vajJN9stzXYYmu5sBS2301rXpacGvPz0l+ql++8t06xk9GC0aW2woHFnVh0yTTtgY1MSo+6U/s9EwIoUNpy080OlHLTKU36aKUx4JzglSZUI2SQa6BT2rgFS4uLIBgQ2YLyLJxkBf1bUIIpPVzkAhjU9/62VKxg52jKZAAwFmrhA1AbNgIg/SAncllhKdb7/KskBVUQ9tp1XdGYSiCdeAok/M9GsUVBCFEbHhh4MhRq2JIsJpebnBqtSOV3H1RptYaWLoRjBbwahllLr3WqlrE+//KOuWAmy2UcKxdYIZPkN1UFvhdBWAMlusmiMuu40vY8C7SrLVdTgP0UMt0MFQgpZl1OjEqByTptUZXhUHmfeYJocM0PF6oIczrmqS6tdUj2W+PmEvbhxnSHyTpSOHDRo9jU1PMfFZA0UgmXhNzVHeGP214dC4aDunejuT5SUqSbtM73rE3NW2KmI8aD9iYWkqPLGujBKB1sdcZi2m7Y81Zg805BdFBhLos7xmg6j9IgQS1XGWLGnyF3c/O8n2hN+xIZ2h7SLaf8vXURJJuptu/Qg6DPCAXZGuZEXrFcm+ABaFYRYVeoTRElkw10ACf/8e5WIBcuJTpWW0kpQ8QD3NN8gmTfIe8QIeq4+SYU8P9Vu7bcs5AGOpFvFevnUVxICINemoc3jfIhfRM/k8AILgO4g4/fTEETNCdX44GM6Gb6zervU938wVgDT3pG1e/lOrTHrDrEbTkbTD3M/uX78dNBTS8XgymKzbhrqEzFyNZuPReN3oefr6zQDDPj6IaM8Ri9F0PVvezNmQKst1uqojso2/122o9aluXGcKcW53qzF4jfEyFpzSnABymRWrd0Zv5dO0Ku1lWrRikKWjeqOpFmbe+y9HANj6KbI9MfWt6fLWghFDdJI/UUgQE2hYSBpzmHpK6Mj3dW7LFAAA/98DFIJnsKNMomFp8OE+eX896IrhBtflOSBR+9wiHh8Krb2XEr5PH+TdOR4RRStWaLeyT+bLiOPaGXWV5A9J5Q0O0gfMG9X6NCevycA5Rdxg3zBhdkAd+IcXVwcAGpD0EWqU/rQqiwqEwpSI9tUhI4iJvUQ6nYDmDFSVY7QSsS6TrdgXWnLXnV//WhR7QJ4073TFKgOCjp+g9gVnty5IdueuRbJFX7LXH/zvkbRzoBtJvoWvkGbztqmOI8sAuF86+8sN+zA3VG9a2i/hnfNNKbfoGzmzRy/6AEME3afP5aM8rvpFG+VQdveefICGlSMDrJKD8jFdzjv8SI+ZNQBc6hRRVIL8OBfKtFF6FjW66RV9qxWCZ1DNgOUJpzPgsikade0TY8YNxCrV4iI2+er+9fzNyN0Rw2k0qyOfWnqIkOOLXNxxRsDLIb+3pgOeb3l8GTFvYRB6yinsAcU6uOTNOXt09WILxEyzE8522YCXOt3L58rAo9Ao9QfoyT4WjM4bbexSnECAEr3Wv0S40yn/prjCnnpaeHdAd3hD6EhQTvuzYMH2jEL2gJ5r+DM/7kFAXy3QBHLMkhIfVIqLAkHLMlJ0P4KzWYz717JbXO1kUu8TD5kUzDa9z6a1MsLEqM0Lthm5Dv4sOVSSdDn6213MZNgindW7rqg95/A0DKO7wBTgnjvaIbhn78iAe6zePUFTxMKo/8F6M5u2j8ZWHSuYmYryULj2ay/N1wQNPLzxLrV7DkEU6B8Pi0WyqVVLxdxeit//+V9e+huprEfwCZxu3vLR/4l+9KKmuNu6dBi1E3Og3y5yUs/SbH2NXLOhILx9/fbN2++8L2IetFiEifmCbQ9+8OJULjdkcK3XfcA8eB9B5MsQFnaHBGsWWhG23N9RiGDDVSiGNkJpEfwGTdAzpz5qDh/dFwtm5tOxY+lVX5LOzoo6mgHq5J0SHKu6TA+SnSfqCL+sSF3XNt4VmgVyepSl76R2pBhPrlLtWw+9sWWXwcDJz86CxOlBLYSoXfpewjasKvG10BoUUyVyJ+XkwFMWbTRYgti8OWYZrvpqz1/2NBVzXKaa98JcSG41bQxtGp208SmxHwTiHNg7f3ytdp+yFErlvFoohxq6C8gjwKvukMBPzqwK7TR4SLv+Jn3PynXDiSP025DSOLEyyHwPHtjTFy61GrTJyi2wyxpzM/QiGhSNqRLzC2Ki25jYqD3+EaGnnxc9CuDVYbzNyGE4k38v5UHgY1wRQsNDKh/74pU4lBL/ir5CIOgAI2Ip4VDQkZlkTyQfswbDjkhBNs+K8FK7WhDMHpMn8i3uhoEl34sApdjieKwQKrKyrmHHADDT2gc6FcAViNymPQsozdCKUDJpRsY9P87YdWZgSLVW/YO5oOrA0cjBpB+L8h1Hyu6UN+0zhts2pd/oiAffpnpSPvRcWC/Q+b4N4+EjlK2V33tYcXseKxpo1WbiwqyhVK61KdXdg0jQnUda5XpzI20IDTPDcShqcRhCbE/cofsHmatTDXnpZO2P0YvZtxRO75Kqpz3/xop1kU3K/dlnI5U15PQUTHr2r8EgMF+5b4cj9VZ30TVWn7o4rna9RKGKATc0osgesJ8qNIZZu4eRZ+kwVNwI/b3Uw4anRDFJw1l3NOqO5t3+bUh/Y4d8WomLEvzY+Xor6suyLBxGA6OvJ6tlc5cdI0sAsxxl3jDDEwZPnLjQlv60YcC3qQ6DczKIcJHM6KjaNhI2YsWPBceeUm1d/ZThTP1i8utYfGeeXZkjyiPaWftonhW+Ge3RqjcPIkQGkeEgAaWf7SSzqRaXmMSCcYz3w5HIZnN9rop564P1sLQ/pkCqF/3ymEeiB2Jquyu1e5kr0QyC5WZ9I8dnZRiP5A03883KF2wGoqdzxIL00f5g1nkePzYpK2rJlPbiR1OYALQGWvu7N5ni9p1wJZhBnGhdBbx0OFUCVBl28D/vfnWNa5qy2sTIEHlaM2GiNEN9oP2iapwbDbOc1RK8dKh7LET/5BY7CYpMBQZOVzYVpCvUIwZg81NnRzhbwn6viMJhm6ps2JGfDatSHj86yNS12d2hMEvAACeOKcO0tQ5F2aQZka+2WatgrmYyhs1U0BEsG4ptCsmIMj+5H7RE1esuFNY/z1t/MSE6Bfp2eB0NZwGfXxqbeIzOPMN+fTqi6tOPvjV1ocFooo42d60ay2LSjcnCHbR63lzpu8nkWJNbyzpJs4pjnVAE6Fim70DMhxfVca9JUriwC8x4J4L6Ablz1EDlOq1J5rnQhdOcpM4xokeacPu978bBAjziM9Ah4I3vgDRjuQCzB0G2oRexdNueQm5UlwbOnAauyHyY2EcQOdIwTKWbjKJtObNOSzhdEbxWVti2rIvNeDPd3FiB+eZmPp5vdChU5Vi5bXihftfweOFJcJM2VYzD0KJvm6ehYsBapn5pmseHF5ekFnqulDClUG28tpUk2hxC2g5TJTI55hjerSq9YI0Y0IOEUvoB2WtQsdEmssrSwwH1QCW1Vxx5bVK4XVJ3MzI8Wml67luVi92in43MtzbMJCo6WD+G2qLAEXyZjnK6l4jz1Wn2573E4hrXbsriHBfXIbvttSOTzwdDfIy1W8KsdxYYAkNPayzWB+hB28VaHCFMoNxJUfgoULKHneiJm7kKGoWpNA1OGiVivq/hXB2NECoFPnroQpg4wuix47jWAec1GTSmPps+7PRHXuxLc9u8ARRkGaiaRWL7hQhcY7r1ZBRp7Z+3Jm9DPm//YGPWNa2H2eOjY13CLzjisbN5tvhGV1zJfMth5A0jqBBhJJsVk5pv+2gTNWJbMyb9VHubxanlpIi8G3YQD4czUZQRDTXeg2s1RU1LrS96GvOBOYxowZoISLvLP7hjkUbHU3OVM2joljKJnYqH3zFNJEpMDE7ZYGcXWLSo6ENlS0+MKrz8uOXN7aTFCjfclOodblsPeLs1zY2QEnTFFChKp3UWz0GugCa4BGBgCIA1BbiUD6TyVSs5i0ytMTAffESXhXcx4/bpONmTmifCLTJk7eDCDCmd9oTZTFTvDH5gupqG/wq9ETZlpUD5W1UpQwbsZGcGooZy+Fh4ckt1xIXkLlbvmr69eRX1vzck7SCr+4M3G+PAoo7WiDQqaY8WwwlT/ifsR45mKWzKYt/jfWnTOheiNyJ/Lur1JmcBdvv+PbTOE06UXwE/X+rkG4mizUZc4+HCXHnH0TvF5QGqrslWI5IMsMY5P/A7e0Q/LnqNOIMWUI0Vm2oHUudOZoeO6OMfIZUIikzB/L4rk+2WPItyvZWYUMY+EPYNK5g2rmAESQSZpEuP3LZk1zUNUc4LhDSCxJYSOcpFB497j11mpqwNsGmh9BpzMR7GBGAPG0TEUsKOxk1zJpAQNF/u7q7N2q5Gi4RCfPnlXVQ7ZesE0nG2u0+dpFQvGmziue566PilxfTOaUWRzFmjqTtT5nhsnnhX6McA1HYxrRmygevzij/utXlAqTGGi25VFYvIBkfCWXQZtcirj/fTkL+lq5ZzSFcIPpObLzodLqqnW02olRNlPNCDMt44I3fs+tUMGhsQNghchAj7RLKfY2FiBb09o1E7YqwnNjLeZ89CCbpuzUNpzoKbLmBjr311u+N4e3QliEgoghAU+4AsF6MaHncqvYQUSORMqxIpjw1w2KkEEzuPmAPAZLF7qSKxFOlnxMY7eScn8qlfxMoXXMgzrcXmQt7ZwvLOc9FTekZUewhI90iR7nhXjpAc1Vg+CIfix+Tbs8LlSYHSzS9g9vw9BmeURbEnKIvzsS7Q9qoyDdBM8SVlrYvNEVg1+geZ1cUMDkh+SOrsitmIJFYF29y7CiEhHxzGFdEGO/H8X6LQUK7R5KNSryQ90OVcASv2R0zkTUh8xVk0xcZGGMmNZwxUksBPLsPnjeOk5B+OOeEhR48BomHKVxEICEtAD0wxwGkRPrJcAKC93QEev5NPuyJT+2RlaUuEYrqP6ylyqh848nl/5Lma3BpN5kuvCpT79Y0Si69gFceDazRrQB5sz5skX8mMcqNxwXWBMVfHPOfSxRzdw7V1rzm2VVCV2Me0khxylEkVTQd0ETaMmDI+UFo29IX7x1UGMDJgSfFAG9wxGP6fPAKI1K7C0wEiqAsVcAgfkJiSMscp646gUnrKouc4QNIFr4MCSG7yql+pYmgCdLamCJnnaDtl3nUkPVsmye1+EExokSXtXQddhUMF02xGiX5USaCYfGlrYGh82yJTKrD+Rc3JNaVTH0PJ6z71OhnH5ADhK8oBAkiyOlmTXiFWEhSlpRI/FdopmXDnS9MjLwVt5HpwQO82Et0HpwMrUyq6Pp15sq4v6hpJV7sqiE72sa8wpMv3V1E78lZ5RlrtCIb9+Adg+opSMkVEFrvGcMhqj3VD2Ba4lIas/gjvNbgT9cTUc0yXUriCHLo46LJm7UGOI+uDUobDPzgSnzVL3MLQ1T8atVTL/Bj3v7sBbZWKGw2fJXdeFno3OBfo4Q5vBT0nQMYY+r1copjYd3WzkYPl9M5JVmo0NhTSiK6b9Zg9mQiHsKT9MQP+BH/SMWA5LKkKuldYfJ7CnSoXFOld19pcKg7dZV8UKlBV1+rd9NA1xVA5IWb41FEviKQ6U+cwtFnEoaURUOVo3ueO1VZgiMTeBxURufaXdcgyO9AlX3l3ejq2tuGUvNgvf2H0VnXcwlEZgUFt1cRuVXxBFBJFSfPnyvOedOyfQB1vkFiAmOcV11FaGnbi7nXe6rkTLaiiqzmmsNdMaWvLg/RHS1vG0+kyvhH8vP0j7L+vCx+F7lmtU+rkLvNBBRpRJJb0fjB/NbyPDGBw2uyrxynJ30qr6Y2dcBVXUPLTFvlmAkUwevAPMJHRrD8c9keD/nC68AVqps598TdY4U9EVATfZsFXR+hKWKYGD8K/Lp2lSjCstyBFou3yGxBWEcFQKMGaWdQGFYpESfNKwicx4UrreE2XO7/Sodhxid1t0yid3AhCmHgfBeWSn5E7DipNO2lq1k6eMOHwpApbpT666EtDy5hYtpYFUMUSiV/meJREXUEs/EvNxnbFIPCv+oyt7sDiDKgIFOcgc1KutF6Keglm5XDnThU1PPLr4RhpfbFRr78Wc/gPP6OdM2N9bTRM0DxH7K3NE2Atj8R3EpaaNBQqFx1OViuAZtoETNpggM8SmArsCT35RqshWrb8s2nBDVyjwXgyZvODBybPqSngMT1dODs8WyKoXeGftwrBcx6GMSbxGLyw8IfqIKyw/RwQR04SX960FfT9gcM1Uim7IMIxVkvDYrXMImIezcHEgdrPWmZEdr2SYj5A/uabfnQqGlPcBcGREooAkBOd4Yyx2lS/IEsxPvIRMAezb+rHAsmXZD/FS9yol5jlpmqaYYA3WVlQWAI2A/x2/TOIHkmPpbQ/vRy+/KV5NrFLOdTXXLQ87KKDeSr4u1cX220WtZ2rjVT9mEE/GgKUfGG9fIEIEMzntFxyeVX0Bp2MOOxO1pVirUgnrj+/onqY+tgwbH5qnZcLZUVne1tL4nxCGXWQjKmwrpe9ddpob4EUbe9eKXWA0AAe2muqxy3Ym/TvaU8hppWNWYvweTlVZjJa00cTvMspnDvgRTStxa3rl2FwjKRtAXpxdcHZk7i60KrVsKIRj1z0Q1bdvBJ+Eg/qbASDjkM57LzaFPrELwsbDRz/TW9HQ3drXvjTLkvxyTujNuS5jw/w9qQ/NcRiobI5Xcfw7//8v768LOuPjaeK+ap45J/R/vQLShonh/i/X0Ymox2qqqsePH0XdZyIcKMwq9ne0+FyppZpWcDTiqGP7ReX9lQ2QOF9F7ubwBE7rkqkEjpR52IEawbkuJYuLZm6z5RXRmlDQO1ByrBXBqI4oLwQpHee9WR45dqMO8LeqBALao74O/zgTNeBMXadH5bRxRKwgn4XlIXc1ZU3pVwjai50daWmuS/4vk/p77475bbFFWPdXQig5M8qUwx+ecBIBS5pBYJgQq+Fun6kq4zmydpNtr/AqTNs7ri9pqKNg2TF9pz7hwzHONnvisc8K2BeoQke1CAtrJbk0VsIJ5SwwmseUS1/LNB7d02XLGZFgVZCLx4IlntIgDlvS2hQ5Mx2OrR6zzgZdSF/RHG/SKef444Ax0AaG+Kz1A+M93y6fqBnOO4zs291jz07cMIU//OGad5hJBrpZhrEWCxgDQj1Gor0SiqiPqQsgAqvIOpOGws4TILUqV2RrTEgQoXxqRKcJmrUMzdb5xSafZQ8Yp0fRvueKJIdtPA9hzaNKd6woTVG2qnTaEvIIzOnl4kbXDV12hDMQmyLFAfbb0JMVnxLrd48Es/2MiEvUJYRvqobYB1rvkqAsO8QYSmO+Z18agWws2nLwvRxgjj7ucwpplZljvQ79Bx8wwvjzRrOq+Y9FuNppyE5sOXvW5A5uyLyoiGKeqYZsiQ+Wxol9k80CCSnfYtEo4dX1rti06M7iH2hQ0PaUy95n1YfQRZOF7MXuuO2u9Ls0AEOTJWEEjUd3y500RdtL0S8fiqOxF/IPgisqfcJ/6gYhP+ow1XI4seuxOzJYYS6eiMMpyZT6TsSKM/8I+r5NQreOQmLmtlFSm47ntOgbMXMlcq8edmsAffpZyui2xzqc15WFB/hMpuEV5E30skJjhrdlAs3wc3WD3fdLzAQgkY4z/YyfsN57KTbPT2xVFW/ckSjdMBJYJ6/vp2/mt7FFWKmd1NtOW3U2jGE7MQE2tmCCelC+FU3oVA4Gd6ngpeNY/w/hkpQr0Vu40dINl1L0JvSpQrZ4lpazRo87mI09bRtG9WavY8cI8ZzkuwavZ/ijLZxUB7Zm4OVwehCGNwhJPcsLT1i4Jq6OUjfwgE93im6l9ZOWNeXKqviOqG/8P03iRLG6EZ4MufVBXzWUQmQWjNCN40Sp0iNIHHCXpHTNQ2Rou9RT8A3sOCVdNSot+nfU/hVPQIPQC2KDeTwJb1gkO+Le8+JiDfyQBOcgApuVFd18a33Za6vurHbebVnSu5yT68QtQpzi18b1FbeyzmNUzn07gB0UZhnH42ZmRxeEpvPs0NmfenFx91Y7cjmlC9RgHju7favPh5Cj6Dhstib1qCOWGK5Pk47yGcOL0d8U65oc8ETOx6PFeCpln3w+V8pBgeVHFO3kIJ2eX9TRBf2UjIxIyUHZ0zOVQLdYBV9r7KVUzOnYSW2HxLMa6rfUnmn4VYZmoJQ9HWVpOH16E0EmJwLNzddnby+y7Y6eYMXnMJb97Kshb27T5EOHaYI+3un7u6rOHKuElyqUWdvOznbh2O1w6gFKjuON3zhKSj621NFl9oL4BoK7n3x7UkbZgNTL4HTZ97TGSmbIQer0Xz1vGI+m82NnJ+vyd+sOkEW8agOG02radm+Z9iO27to5kixnbpRLdaLjCGSYGp69G+mjTEitueYEGyvy9UfHwofpmzQka1NOVUxOzor/QvxNU6p41RtVYgLgHE95cD9+ezhMXKE/1oRnoOPivBsxse46aOtGfvDW07oiMbYuzv/rfgqEPMcyy8RjC2m8bJIRP59nXVKQfNEwzkaSUVAUmqPFcBgDKI1yprvUhMcPqjCCOTPfXulf2SySTBdxAvcXDB9jnRxrxNk8cy/wNBKg12Ww0gjJxkQTUckEXLUggrMoQh+tXQMmkHK2/FpJN3F9ls80fGTIU8dDq8O14ITyheRc7BHQKfI+19Fpqqr2/ohCCJsheHnkYSvsNnaj/5vq+DlnjLdxd3kC65o6IBIpKZ4dKuDQbIUSEaWBrLvVFnr3JacXUzQ0Ewb9u8JP1MG+EoF3pidnRiFJVHBYAkXKNZBYECQ4VwJ15SDmZ0v1y+/S9eUoKNqCKtEFNXo+vb2tvOyY1KN3uC4KHfpQGQSslTBWrRQPaIGo/w6xZ49JNuCDMpIKg++Oj6a2tB9m/qC5aaDgsZklVaeUJ4aHAiKTqnRjDDgWEWd4fVUKkqOp4I1RJ9YY6HuesapGngxe0rIdSz0lUzK1Q4t9JpHQm/wDNN/NEQ1nyh0QoNbZatkUwYQhyax4PT7//wfKo2KLWoY3aku0drL/Ih3LHKQX1frd476pS2Zkj84Lz5FL9ILevkDRCrLbj5duHp+ZJC6obRFabpAuooIUi1Og7atfI541d7HxfKVI02d7HHHHlf//YXSlskU911uJ4Dpc7olzy3pVOuzC3T8icE46s8eoeY5yXJuAugHNny+KWZ6d4yMniHE/RuRPhXn8XeKLpAIt0vvAROeFqn79J5HfdwRR7TrJNBrda4t1UmB+s7SiC+gabRpljVsISlR4tCIPPOu7wmdE5fswkUe8/v7t5O3c29bLuo7YuHXeRStW68iefnmRPY+O8mxzN7rMpWVNZRev3xVyuRl1+ZSs2G70xf/ZIPctROpEuYuM30PBY60lisQeCoWs6ERpfyxBVTfRgEsNl9X6JPtCnZy59SL42GPMFRvY/r4f3W15IUmr2dktnmX08zYCeBry16c1eWgfQHUztqh9jSsuFvS3yc1am9dEckE1C9b/VrnK9tHYw5PjmQCw7BQV8+PxDK1uy7fSLKfcd+2BFlDVjizYe6EG0G35+p+NvKdB+3zPbUSV00ZN2yNNvTtla48FFZ94npPnJsOAvexvvPKApAOW6FMiyU49kn+ZCNL+k6Cou79QpfAv4EcxQOwKX/yXeE80xmX8W6VYaXJWGIBB24y5Cc5NNj4+1n9GZHp9U+fpfGkkTdR16+AMwGtkqP/7G3zyiNJGifZUkhlIiUTbf1ceoNs/ZrAU6fvpDywxutmuLuVd6xrwxT0Pl0Zxhjq252GBv60fA/wf+f+6BWwOAQBWwpbv6X4kVWW7A9UJD52ExrdBdG47yB2A3K06bem5teJO4gvvASyNVmidehtkL1nndfR5m4edqzup+uDAWj6S9QsZVL/dxSdyVn6rvPSOCaU1arPYKkkDEDbLKOgLgYlmuSj9sla4DU5eVSdgyqtELx9RA1WEVxDEdYKjdyo61YPtS6mDy3RQKPhJ0X6PCMaCHcSc0wvrBjH7m+6bcyxmy6sVYlCEIoD+7npNtottkWZ8rGPvumfTEvoKkxKVBeZKTrCBd6wNJRMiaDgAfLh8wzYCGPO2k4bwIJxoItDKsmTQIjptjGN8Uw4h5MYZuE5y/ri+lWpb1GjMRfeMY36bO2lwkmSrwGyRTm4TGXUHHmhutiWQG7Ha9TWs8B1SntxJ6VtSlcs2/TMNp6pD2y05plvsIreo9o0HDVC4i1DnXhy3+APtik1czzIDuzfxxW/FzPYMXvv4dl4r8AgY4hQhNh4Zf1swS3f9qxyYQGj/8owzQQ2yEfH6HYrEjp1IHUlKEKrLKWoWEMP9r7hmHegIoewANnpeGjQ0UihMBdF4qVoQyQgR7ULqG1VrFV3I356Boy52KsTaFoWWJTzunc7WMttxxFZb0hkpT9uWKig/f2bihGud4AROQo2RDjISv/yFTqWuEH1UkVNuhWJr3SJdY5E8mQX1yHwE0nsW+RvnJylDhK1ZFU9D8sC8R2uimYhIHSFfJ9QPCe5EmBWaPqHoy5TUINyVrDJXVSTX58Ufl3QVNaG/FlpjGdhU9ca1XTsdQ725spPvb07uFqxcemVtWtNpjqn0ZnmyZo3qiUn/FxIhx1v8nDM7mQkqx1L3/DV7cB902bq+wx2vOh132pdKBp3hf61KYowbu9i4/88jCjyvNb+iJbmXXAv7wfhfboMQmlUpoLXBn7kbuLYwJ2AWuQzLlMPrxwyMLF+ijDSgFN9VBrseDzzd829ucCKyW1DfQwSNaHkVO51ODMgCg9JFcCOD6NDZfIKQXQeQH8pD0CpkFfBHNMVDrpGBev66n5yf3t/LwZfiNH0iy6hCf6YAgKRBPUNlQlWtYKVd64CTQ61BaoVamsTOlZNZ3CXhOlg9mY6XjMswhIgBfAKLL37EE+BeHNUSpowORIjmyNxOnt2Fs+enQWY5A6jbpxoRuxETF+6ROtGPqqkNspWU5F8FRNCFNXR2BvbzG8buKnyP2zvHDuSK/MyGZkMc8b0bWBKQiYliB14bUtV4YgowoPcfqxkdNRPqjIzHRiY+vhSM02BtumVjOj8vCN4FzHWPwVuccxT6I0jONgQp8KUdXkWbIAyhU5SSmGTdrA9a6dCUOseXcIJ/U31E7RoR3UZ0SHgY29oKQ5LZWR5c6V3rqSkVsrSZiI4ZY4zch93ZFzUz1Ku/QtHDYxLVwnG8jhlkW91pFE7FHyOej6tPevL21reOxEHaGv5JdaSzi9ey+BcvZmWzpyzVES/ipICNXRLmtKJ2hNhR+Eyhb6RzpQcbsD+5VVIPr9u510ncGYpXmmNP9XlkUqwtVfUCAcd4i25hooqiK6t9G5KWyMlu3YqD2OtBSKqnDIN6EzkwKkoS6KsrlFM6V8VR6F0osjgrIzm3Sgld8GRjEN+NP7cJxS5vvTsQmwkBv3ucXYmUdYDqduuqSCCBxd2aM495oiOFbO7bKR98eunDXN6hzAbAaDiF8eFoW1NpCh7Ie6O+JyskN21xYcgSz9VBT2m7M9UffGmyujUUnC3TYnbl9YA85p5l8kM+TKZyfSLTuOzyCn95otxramJEcXFMdf69b1VdBvVq/1VOzr2GNdmhi04Fgwvi8GBVVAY1u8DVbx20pD3ioBg7WfMcVMlRBO0ZfIISutnioD6lpIFuEsVRceWbuN9Z8OpCsO13vkNGgbWbhK0nnCX6RYycBwZJ7bgrLIMaxVjkAC1J3MvHDTZg7QhQhWAo5iSjkpmJfEvk8mD8hzpWWQAlTDra9zDH2C6HZ6qSpDHsblMule1nTeiwJtUrvVVdLz8TpdLNpgAPqzLbmIM9ajs/vyV808eqL4wG7x+/2//nbb7a9gdGluRct4OnVXRd+znTpV3tXULkxKTOmZtf890GhqXA7z+6W9/+w8/YvOqlodJ/+9Vh6IHWQaksE6CgVo9wYnfqWNxj3NdJrA5ToVZtqezhXspc1DEdt4p950FU+1lUq/ubOQjxmY8GZXgCY65b9PiFcCdCasaKbsxqVBOrKxbq/RjLOKuQdAWT2srbcozRnMcBjzgBvZwr70M8XM2EW2bsF3R2v0BTWS906pZROt8MnqrZq8o7inZymfjs9CKdmks3+TVdD6bxReD9c7XYRLV+Ytl2ypwNfr/TNVFmv2eyoRuDWVsLp+yzGOPrRvgsnxwM9swF5xglKnTb+2H3i7XfQ5Ye4aEF6ZNa5hqlM+6sCwf7S3xxngtMrQ9zc8UwghNUWiqGs66o1F3NO/2bzunZztqOoJMjXWykvEc91Qu2iV3yihOtR7nPE3fK8c/fIfxJfU8XPvbQNngmjU9ZlPgHgM2Xpl75Qdd/Lc/ixVajF2lZSOSZ8phoii8egHbM9hjvfNKxnbEhPWGz3GUB/Lb6AGo4EzYrk9ixG/eVrr7adSuTYIgcszrNCMjwp3rEHHq7no59Ozk0jGNOkKfnbk4j76Sm/pUFqqnSr/DdBTPGk5HzhRWINWUmqEjNDQcP2ERhNV0JddzR3U0BX/Jca6KTyFHz1i2SSrx8j/B1rzEZ9D9/wdQSwMEFAAAAAgAAAAhXB83EthkHgAAxHoAAAoAAABzdGVwNC5odG1s3V1LkxvJcb7vryjjYjIWwJDcXWpFDhGxyyW1tCmSwaF2Q7psFNAFoJaNbqirMRjwRPskRfgk2ZYdId9s3+2DFdZxdfeP4C9xvqq7qh8AZjgrKXSZAbrrkZWVmfVlVlbh9K+SfFbu1kYty1U6+eAU/6lUZ4sHA5MN8IHRCfxbmVKr2VIXzpQPBptyPvp04B9nemUeDM6t2a7zohyoWZ6VJoNiW5uUyweJObczM6IvQ2UzW1qdjtxMp+bBbWyktGVqJg/1d//97u9+r56s9MJmC3Xja6BCPdXTm+q736mz0qzVx6cnXPaD09Rmr1Vh0gcDV+5S45bGQM/LwswfDPR6PZ45N4iKWaDKF0h0qe8N8f2JDG+aJzsZrCnULNXOPRiU+RrKKHWa2HP/bFroLBnsI/bUrXU28fTSl9MTaKHZUmFmeZFA5YGySfSVxvhgMOiutyjybfhquinLPKveLnNXcouzNHdmtNYLA03CFAPxVBSoxzdKZzrdOetOT/g5tYb0+rbWNk25KZPBKA1UXJrZaxxvuazrK3777u1/ynCFraaYfPBBRHkObVB7KRFZ5DhMEAMQmsmrJbSExCp8qWzpUI4yMystjK7M1WP7rR2rvzXAV/yo8rXJxupMn5sEZthtUqhRmJW2GcgYEeifzvMUiLkP37FOTPu8yFfcnp3Tm6DTJDdOZXmpcG7OTTEWpkeDSvUu35QsJ9rZxNAA8cNA6cLqUaqnKH+f+R5BgUpgoaMq8dT6V6AcwJvUl3BCj5Sa6SJhNq6LvMxnOShT3CbUWd6ZvLIreAZzcad6SrRUAkhyRDP2uAAVBvaD7ENneZZ4wbXZelPyiPRqDdq0+Gb5BjrHx6s8gelLzMyudCqDTQyQbqcmme5Af+zKjJYmXaO4UtcVIWtPxNx6mSh0aUbAbfi6tEliUG3W3RUUNUq1gj7OiHIHhBhlHLyABp3aZKVNFcwRzuzcFiua5TkOWGGXKocxb9brdKdYBUGcsFFXwoDdOKQBZwr73NosybduZBJb5kWgi20ez1BhBp6PrIX0bJpfiJaCsYTJHkzUk2yWbkB+tFrnDozkOQkjvjyCe1KSWeEZ+AjscBGMV8/x+3ZpZ0sv61hJbTVwDWok9XBBAFjqvJSCobepa8tgvsZSOpTBU7dZrXSxmzxFqlH/mAjQE+ACtCyvOzkGfEeOPYSptNhyIIIz/2zgC29JzdapnpklKfmDwXOhpy1zXd18bvM0X1hYi0BQbBl0Nq3efINvjunyHpgV1IWhmoHd2RRmqEw5Gx+QfhbmyUOTAnE5iKsDw0NmB2VwmlrHk5VtVlOYvXwO2pcYMGQww6WqyQTxBQWdodCH0ygTd2gaK1MiJYJZ/Cw519nMJO2JW340ednWGbA4H+0d6mfq4dlXamvLJcg6WCnzDdr6dLMiy80GyA0VrBawRGzJKpH8jtWTkkYJvHcNNY4UNZrpuU0NG+EBmkQmkigIZrumfoTl/YLJn/VsZtYAZwBUnA9Lc1GewIf2rAaGnPokW8YGttkHIqbB5FmusCCMOAVdM5XRlRW9UWeWGl0MGut8tK57tX9pVrBcVct6bJ5wzjbQKYE+wDLKzzzREs1d32pR1fcE8wBiOcLVASqweVDnOt0AoXcGk0epXjuRF3UCmuLKfHV6wuWa5QEePjULPduhqQLB2JBVpKXIuLrSCRNwjJ7xUgMD8ItGqQs0W9yjGyuCOEovAEcADIF1o/Bm0vXLmFiTF56VOL+BdFVMoXmPLMlRVuornVrArMgYeFJuXCfHz6tSLb5DmXPrxC6+qL/08X3jTPFNYpBtAHZ/At+UfOurIZ2b5Jvp7husXlENcz3dKXx0eMICFdK0+qAdamuzx7FPVuhw1BL8N2fPn0VLLetvxCW7YielUupvXZ4NwWUg2wldnuCDel5CbawaMRfcSAyrH11E5PQoYNcEP7YmJZVASUSoWYBsOHUDB3TzFE0OgBpNNMyx6Kgqgyh6Cw19NFBuDUsIgQsopVNnujBZVbECZr71o7Wn2US4/EGPhd7hIjXvHhLhbR0bnUivwP58YfUiA/NmZ/Ta9Rul0LafmYxQky13agGdoYkP8Wv9+vsz8WEnR9v4sNI1GfleJj0Da4o2dFRBP8ayDWZlUuwbDyu/N461ejqaba2a18S7LsRSSyS7F468TmfWGrEHuA64kApqmy11Bl4s4WudwqqhvjIFmlt1pxJ7B+Am3TXlvnKBI4HvtIkhRMgXC7+cNYYrgYSzJaAowOBL9E6RLjBfKXjA1Bjg3kScBI/gGRakeQk9Ud2wfM21gAaAoM2+E+sAvJpk8pU1W/DpF13sbsLT4E3D6S022Qh9pj1ulnfwDrlbWAaXNVVLoTwdoRwPxHFWVLEWvQ7BxqrBquupCR9xjINX7cCxDQVu3wT7SExhEXjzECh28aYZ1KkZ/nKTtWM74WzNEM6nffU9kQ+pVGcDUJQ9h07t+oLfiiK4HjULRg1asSiMc9VKLd88tzyrp7pytntU1vOZzKp8jmZg0vTn92kW6MOmpVMVm89MqTbrZhTsGKPiHWeTjeYGlmWEAH1kNv3wjobO2bxUlU5PKA7VjBvis5FJFlKr/iott01K5dE3WBCGtb5EmxIUZIvzpbc0dVDg3dv/DQONzVios28AM0okFL+MOIQmbGFLmxfSeV5YcH5JxWj8JTrAVe9fFHqBEcOVfh1ToRByF2gKM8Ap+dYUY/VFvoHphMXDggWBStC7KcdhgFWCfvAJg4tt4rd58RpMxAxtsKYILons5HR5m8NF1XuibwDLD8bKJUQpUV8w+bcnMLu+wjzNt2LTpRtZiWhslevrUj+iNeqMepHq3Vj9FKNdYNqKuiNY8KxRUzPPYfkCY5qhhfOySwuRKFYYAJ7qSly4qRHIJ0rMC2n5u9/BklerwM6UzYh3Q80QZptEwnxkzKInh6xl2FaJ0LCSmE1ajuSJVOvQciygpnbhQ4cw71B1jWElmrUpSOmvQUrFz63foWAgF8UOK4wublLww1qjbXcYdbawrjQFrj7NzsDWphjA9iUu2bS5oOBhR8P+zSUb3GQORBNNXbvJ+l13oxJkwrZ+voHFEIAtC0wlzP7xbGnrcKWPOEUV5RkqpA8/hVT7YmvcqxH98c8ADyZsFbtqEMMHUQV55MnB6BuH8+nTaJ7qxYJYfEIPvMps0rhbmEEss0knDXM9eQwtuD6giFosk6WIklAtY6DUHWbmRdKg5QBjmgvml+F44/h///j45PG7v38ruO/dL351W+lSuXwlNuUeUYSL0qLIoRC4CkVBIFzN0w2YDzczgAzUlsKP2JvC3sbqlUW4RNWdjIDAcmrmJcZ+0SskJTLFDK03YIPxXqxWDyUCa8EIGa/1ILRgQeO1lsK1h2DPc9wiineOGutWtR7A50zX+qOnFdrguG05aNJSF1XUty71CL5489qkrWqMw8MMzMsCw58xTV09BEOr+4EVqlzpFqgJOvqSS1yx/ZUpC3CS9rT/Yy6hPlSiL1fsiVTmAMMUW3AnCnbVnigov6erl7xO4MYxFb1iP/hpTy8vihzm3UmVVh/9O8QosiCnXmQbjhVojpFdXvroLXXnEpHnKYJw9iKAwSOC5FKyMkr1HtmI/ffYEgm/YFa4kAqjmX7GRvxu0BsrjHpDS7Vu9EKea9g0FapCo5PPcL/FkFgcCEv6rlimwg0pk1ZulDO6mC0HCvVemIOhDNfYK3pGe+B6g2EAjGKXUDWfzzujd2KxX1Gknmj9a1dtAxUUZB4S9EPoakuOraGFRRrAHufKIQ8AoWHde+rdb/9pqBjomtW6BBiDqNTvlkMrvGdVgclHboZ223kwD4YZDTaKH/ZQSwEPFCVOXnVwMAqjgojlm5ZKyXhfyGuiWt2g5YPg+s3Y9YjKBY5Zy9No9p/hqtHd9zN4Jf0WdrHs7rgqFHX6+5ZK9iEtB5C3DNUmeEeeM6yzbik6ybMyIn+7rtKpl2AYitdRoSZOIAwvoUCKy4T6wStQONCX3jvBcYIY71obC4XeDiYvayeme3MAgBOM9/MKUmCeRn7ev5kA/r+p5+RV6CfxtmGMUHxzmFykqC5jGPqIqOPd299WNgeKrAB6AP549/bfqDWEJuQi8ePB5EePX/RRVpXpaK+vDuM8TJ4hvEdaOtXAeAoXdNeRacd98BS3oDCAjvPTV75mRsS2tb3ADfgSHOFMaaeMni29IamrMLdQeAxoNSzLRckbzVoVyPK76P58fMu35k0CpujAZxzGsGJkDk0VAv7AkADZOzVlVcqn34KouUgQcOPhGAPcB3n3wUehDcyMRPImCm2+f9zdR6APUpAT12LFeM6vgG9gTDMKVgnXW29wH1u5GYkiwjfc4abosUyxa+nUrfHdwQTqZ/27olVcfJLlxap/H+/2+FNqava6r8QdLEGz1Cjmp6JrCmSsn+8IV91jsRJ2+eVkJoI7JzmJgNgYwxQWlNRQx77idmlLM1bP97MLViDhseucLe63zVXsOOCb0N6rf0ggYC789x7i6Tn1dWExFylUP1nHMVst1KPDwtIv8NykeLET9Yx7EIzTRXGwJr7J89UoxDOteOfPckQXjCXa62tjuWsteRQ3ilY0CSWdznR2rl34CpEEP23QGzRIKxmsCZmPGoYP4modFXHBdEE9+c5om3FipO80dpwQXGdCtNTqqpFCSWx1pVk3Qsv0HPMG94RWqVMuw1LkH6gb4LWr12Z3MwQgv/hVHXtukZCCqgnKD7/7cC2gxLkuwoEp6xSWNCiEg8ntP/ymYwuOmrLZwUHYLB4DiPWND3kAkw+DiDlnj7YYGrtP29znLFH3ui2y4bYXDmkL4Lse2I1b2DPBXK2SMPqbR0A4ZO7XYQtdwn8MpRjby7uZRXjbJC3kvwQ7kGDEVyMCKgE/ABkxmWN1trSAkD/kgpSzSjFvDA/X8Aamk3sE++k1mQnqG05Ln+nhOtIjWKTKQI34K2tRXyZnb+PeEITtU6iiap6+cet2tYjl7FVznoc+wUxSHxlvEOMA39R2NsJLUQtj9ZDFomDW4o6CV3+kZxwM75C5ovD7CLwkV9nBxqN4/+dprhOfZ11hX8qwbu1+NbvCxKrKGeDoZaEzJ1kvnegeoT3tQnWvD12bkOQN9C0V7/75f2gnolu0QgtCkzPiHbB+I1J5elS87SjK4/0eYMs2U+rwttDrxqJaEFyv0u8GamUzAl0rfcHgC+w5fajxWDdfm3TFuYAS4nBSm2UJZ2tjgphm2ya2+TfPi23HfmzLZ43ZFz7a58dKr3He0tqYpMq8wS8hFqO3jSkUEVP8ro15PxlM4A8uM72pgLf7394ZTO70v4W297V8K8CFt2/t6eQW9LLn/Sfw/pPo/X7AGEuDRK+YafZc4kKcmq9utwWhbQT6Era8ekvyTb9tWffZlJnmgXKAIcibgNUXsAS4v0lgDCUp1Bci2w1L7iAy7fggLtuIXjZa8bBREog6KiJSDnYWK5yJALq9v871G/nQPq+PzrzQ51GNSf3ulIcX0C6tDfk2S4GlvMiGUTQE9qRfzucGdu9hSU8rk208nfS5kfqIz0azZW5ntQ5HwOSe0ujf8i4L5d8NFWAg9LK6od1YfYZphdVaeC9YH+exfMA3XWLNOFqoKGYMmIIaimwAzga69y1dp0U2ULkGuOrWLVmsJzHB/ZoW5XEgKV4cDuwC1ewHj201CCfbxUzMeMJpH/xGG2bKmQcd0TvEIAwNwemduzmso1k+KkEhW3b+NOJUv/Msbi2GuPxeNCZoyQ4V2GCURngTPu3M3QmZcVkufCGiHiBMP9xyCQBrIZE5TElFonx5jF7n610/PaIlV5oUrusFlqYj4HcPg33s4Rm8FYDo2BEH5lc6TKKEBFArlTrzNpzF6DtIgRZuw2cKi9EM3Pfz7wN9YU3aA9WvfQoHVpru+D+0h8ceKqND08ojJAvdYGBo+VgT3IiTd7o3klosfIwZEXOMIZVDhWk3RCF4I3m6IfUDtsYRCCrgrRqn/sakNfarG8tU8+uhJCAKlb9XFtDCCmbnFBiZF5wleILJMGAKj0gD6loUG6F6aJLj9I0UlFZF3DOpj97UX6NVAQaUAM+rHC/64tnrU+KormykDqL2q4d7cmriZiiY5OJW/LOuxMU2Gl2Ry9GZFNgDKlH0qqqIY1dV/8ETyWDIcn/KQFbGyUujHfQdrjtUraDnHZuH8YrKWOIZ/A139agFxhHdJ8uO275rLUVdbhT1hQpWbSNupivcHELF77SXVKMzg7OVsXmC3OvwSYOObRLzm5+EB3P3SQwz6dJSt897JSVaUuJOrzPtzGJPrO5rBCqViqO1D2JvsWxKe1UqBtZAI6lbcRlyysJIAUxPgYeVKRtD3fjuv9TZox/f9KsQ738AVnMgtpxog9GoIRlP4NoSEZUhKfmMG2rNdWMprMnjIO2BsFGFrny4mems1jFYro1Dwcb9CY8PZUxg1jMegvrQP7sEeZyXcYg+z7+IYYQxJesOSPeHsMpmeTy9TmUTS+wMUjAAkOGZ7cwn9sGT7pTcLq82SKMieNCdQQULy9YJfhDuImpgDk8NzG0Ss3NPyJ5aGfDGlE/sOXZvgTzHah67af0aJzyX6UT/wOUYzNoK2KHNuqoNhrErm9nVZqUK66rVMsstsVwK3uTNidLV8AaVZqx4q/bnG1sgasM6gJB1ugWwS91m+xO+6tEAS15VX/awZHJmF2CSY8+DHjWNAp555jdN49DKjAOat96IxBlEeIyNR7nAY058ZpCzLIh1sOwARlphZhxwx6RzAaDAxnCnOdi65q1ZWELGLT8pQfffE9W7ueyT8waTHyFJFNgHpW7vjve1wFv4estD4vp7ykZUoTrCs0sFPESIdyN9YWvmP7UrkkhB8QLmFJYZq6cGYbDGBBTJXPFHynElBuxqZyrFBsb+iOpPqWb7FoPdaIWbEV33F0RrPbYbS9BTTBfnXjyV3AlK9KTMT09soyN9caWOfrJe93bUe9cHcJMa6twIwf1270c2GAaN4pM9ezCHt/WgpcWybOzkgXqFsfo6RB4FPnAuA/DM3+kAS3RJB9DujIlNLlnL8RGIQmJXgzhYJMelmEEc59diwqk3sia0xmB/4B/EOwDi7Fo8lP4lms6HZZGimr/7h38h48g3eGCNN7lkZdGJZ5zI+/G2k/c2FMmv7wga2EwxwX7WQxxnhGGJsO8PGz0fGWHrCaHhpLxv/Awl4djgGZW9WqisFeJA2q8p2DMMLmBwfokH9bzx8OwrXgtxBYAHeJHGSpMf3xGlgZrk5NOhcsogd5cO3JBevU/Uxgdsuqi+chAHb0K5Em1n4OQP+QssWvoC/2JS91BxsAaplNCXN4YSjoDOqi2+tc9ZfZRgGqRdwACvGpPozdqT8IRaRgGKpSHjd0yEYgnf3uTwdc9JJY5RVFEiBClQPuWET8AZmJl1MFjRPDh2XN5vlRq+5/yM5Bki2yUS7fryGNtvBbDhEfkocadgFa1yiiLINVSuLHI8R0h7tmW+5iB1qgt6xiAzi7Ga5BVbJ1VQdCj3htxFTDFCc0v2GJZnF1QmD0KaI8Lw0w4Ud6Ez4TSWGivOx6+PP4CUrpAqV96T+rbg7gBLyqVTqFF0IoKe870IlA8VJ3Z25C3V6YZR0lJ/LSxboiPnB3pEnRKDknxHChkk+HZELb1aDyYvMJxJF0SVm+SYvkAVABQ8h78WETmx+hJ52GIRomswdGrKsnYyP4d+3r39NSWVwX90SAjXw1c874Jzxu44CQbZ5WUQA6br6Rjqs5ck7Tu8gQMkuGrC4fm3HXv0gTggrndRXiWbsb1zjf4JrrV1llprEEfwFvClW+nBhP/3XthhC5sgqJQPezJesTH6d0TnNpubIssHE/nQ67fYc+5ePhw/+S8pizY6M8Ab9j66QB5VlQsfRBZQ9TJjaTr9qThwbYvoFNNX7NNNzS4XK1zwMS1vok21JA0rQy3JjmgbMYwm/WFJWqHe/fLfCZ/98j/2Tn+RTzd4aOg2xjxy9cMfQlU5LEVX8/QkV2/QQJytcJngimIbj2Zpd84iij4m9lfGE8zXu7f/Kmn5AmVseV9g+dxs2cC1zR0qx16/H0r6MwRkWKnfttu/j1pyESmPsu4bCffTTGIwlLvPloi5slxVbixh670kLum4dHXIUtEJ57j1Nr0NmLEHacviewWwTRA6Qtx+IT8GdGOagJS+HtjtB3Jd26x/THDdczjuzwVfX5G8Pojdj62loyvD64PnX8hD9Tt9XmLodo8whlDL8YHYgW/hQPggkPZ23ODScLk66RjeUtE4J+e3qV7RQeJol4prj0o5YL1/nyoI+GLCDVce8ZV5jcBvGINqlx00bduE9iUob4TnXgpWdPQcK9xny/zYrsOWeS4fZ8v4wLZkUMDag7caXYNhCwd0dYvWRdxlTRlWOcaceYo9d97PpHVSvt+W4XVS+73vrptMRFj9zWTxuf3KarQP5tdSEh/Mb9ydckilZevqsEJ3n0HyKkdJcVdWuFCxghuo+J5jYKtTr/EKZImFUuVxpVRHKGtw4a6RMxObrGxl74XXBcGqhCxe5kkrRcffFQRF1I+5CNhavP57ZrqPIzSMBfH8WkwFz95xhmLWPhx+nQajHtTVzUUHiUOZ9VBwht+7LeHBXI8lufyg3tvMHGdBvPC8l/3gA/wHw2cdMbH2JibtqPfvZMrOY22dZdf0L3kLMbqBDNnTdWxAxO6RLlLr7w12MQcb76LTAnW+MfUgkYXoDpcOC0lle5LwvX+qyx5yojfd+fd7L3c45FqyJO0xsEE3x1lakvPjDC3nxLKIXqOBrcZ0dfsaUnZDQhgF3ahQbzPe/N6NKw3kemxrF6vfw3x2OnFdd2M+b17TzTd3y4XeazmByqFUNvJ0sycGWwO2V6E0DGfxOvHki/b1+4ctuQhnZMg7b9I7iAbl2pf3hoPYzh/VAWPKrwdUMROOU/aicQ/OtUKqekxX1/kmgX8iQMVDuR6tv+SQrhNO4cbkblQnb9Vf+y5QoyItlBUbFeEO31QYybN/JK1/LncZGr4X42hnsRLpfVSsMQH1ofwag9yfHWwn0M+JYBrifTZVz3BvCe904x+AUM/k6rbA4pye8O2RYpEiYoP0/UEdjqyeSEofy0NiNd22G13Ciea3ys/39xPzb9KwGXlMKf7tAzqc+t++Lx6vQ3zx7Ed9IK60COFePXn8eC9uq+DtG9PuGrfPWx1Xx3xc1yGfuOzdi48Hk7vqD79RH6v+yzJ+cPHJYPIDLPXJnlKfXtwNtt0+xeJ39xS/feviB3hi8fYtLAkf95W9c/EplLyDJT9VB66TqPLaqwMYLb4la9ti2208egh/FLzso+OjW+E5R/i2r/BdKDy52yzTS3BNHabFY9YvlK9+dAPERAXP20MKXraGlr7ZBmQ//dnXfSS/sQBGfvak9+qeDG/emuD9W/sH1BGXkjNTePtNfctu82pH1rDG8RlURsxVDI1yfXAf9cLy7+N0HK/xu1UU7+GEArfUa9q+pv1FupjmPgboMU8EL/jm6LZsPuNWNh3pob3QIHbvLx2R3yvAFYaB1lG2iYm7hHUKzjC1fk7tUCo/jl8OeGEUH1f/I3P5uZbPsN6Xzv9I+07aJ8DC8128RwKzACV46zc4qiV7JdX5yxVgnviE2TEp9J7otbbF/hz651lw7PPyB82GEs6dmnJrZL+WzvMFHjxe1CRt0oXwlA/STKUP9F8EgzPnCaoyByKoarOUvWhb5VWhjx2g17gZdsCr4GSUYxuWjFJts83KFHbWIYash4i4VzrbVVIpH9oXibaH5/lMI6xmY20y5y/mwCQgv6FGEQ7v9XvgEs/U1dhTk3EEh3zhTiZ1XJ/aHjVOfv+Im4dNSOi6k7Lus6LInYRVzaFvinmyXZoixFvNgyGFkUvFJC3MH7XwbHhRS/3V2MvjPYK1lDN1Obb+2R/dlqH96Q5v9yzDQhdvWFTrcHi9+DE/nyCtLJrZ8uFZ206/TSr2+W3Xkg1Q/YbCEUkBnh3HnQXklb29wDMk8L9lyA3Hj/oXeJB4nfrVtLa1Kf/UkDTDP+jSmCl5V/+iR8dLPgV4Og1p4saqm8g7ZITL1WDtUPp7yFCuS7+b2vmTqZ0/YdI1qmQXn+WVF8EtaNGTKO3Bj8BfCd285UgwHOUqgjGlXCzv8NWZEi1v3Q9OfsssoqJ62D6h2dHEPKeM/jbjKA25LyDPAJh/iEHABxn5qRxOodNkG5uWIyt2g7JeDWW0JWau8ULu7i3Fo88zeKbvC4y1h3Uphfezc13ZP50WrB5G39VnVciLyw6Jvf4Aw/d/bKHx2xkdZqe01Q+4rpvHy08x3rIulStm/GPO3/L80FNsTX6t+YR/s/r/AVBLAwQUAAAACAAAACFcOC/CZX2sAQC2lQUACAAAAHN0ZXA0LmpztL3bctx2mid4X08BaVyTSCsJJeVTmzStoCmqrGmJcot0uXtYbAaYCZIoJRPpBJIUy6WInr3que3diI29332A3ZuO2Mue+32IepL9ft/hfwCQpOSqcYTFBPA/H77z4fHj5LApFsnnySK/KLJk/7pY3jaX5fwiqS+rm3lyWSyLpKyTZZFPk/NldZU0l0XS5Gezouaf+Tyf3dZUos6vi+koqZa+2LKYVMspNTaoEzR2vsyvijr7zePHyRF9RpfJdJnf1NtJ2STTipqcV00yyWeT1Sxv0EC9mjV1lhwWTUPt1MlFlZzlk7dJU3EPi9nqopyPkpvLcnKZLFdzGVRN/STFnD7RAOXV8/KPZXJTzqfVTfabwaqmQs2ynDSD7d+k56v5pCmreZIOk19+kySTal43ySfJTpKW02Gy8y2NbbK6KuZNdlE0+7MCP7+7fTHF521XYW/35Uuq8wuGvajm02K6lQz+0+ZXX33+5ReDEaaG572/+/LL8T49VzSsJd783e7Xn++OB8l739Trg9MfXuwd/fhmv9vgZ9Pp2VnQ4Pl5/uTzPGwwL87+bvIkavB3b17/+MMhNXY82J3Nkkkxm9VUZfDGWsbDs3LKG6Dd4dX+u8lspZ9/nNe0n9h6PB0UxbROfl7ls7K5pSrXZXGD96/y5Vsqf7KdJAlttDxuJdX5OR2laVLNZ7fYrlmhdZIrKlEn+WIxu6XhPv70Uz4ck2pWrZZJRQeSS+XzKe8jdv9iWa3oMadDWeNf2sqSJrWg7VzR2LeTAsdY1iPh2nQ+czl+NJD4bGbJp4/dMklDOzgItOuHGVfe2dlJBtr4IPnzn6P3fjwDv9pv6ZyhFTrJM24JP6S8X/DkaTJY+qetoFB7I1CUnlFowLMK+irnv6PuF9bdKMFoFnQ6aLLc94V8pnbHGHzqnzeT//yfe4c2bBV8EhdsDw/FaTxJWOWzuIo7R+2mP4/LBWesXfILlNTz9rxapjzBbEbXvLlMvqXJ0ef0wWFmRWb5BVro1CjnPJQ6jYoOh32z+BKNPniAI/r6XOoHF/6QL+dqPivnfMbndCgIUqXDUXJNDU9zAJVDwFH6tsLelFeLatlwWXmen1f2u24I5NmDQj57pPNcU1tbycambvBWMh7xNbYis7Lmho9PRsnboljYez71WwnVw6GlA0QAdyAPe5dVXVCj5/msLkY8d7459oYuCo3jsppNa/dqMauaV9WUWhxcFfmcWqrLC8IB9Dw9PwcUWjW0HgXVaJarwj3/VE6bSx6GvtjjC07FBth9gLPV1Vmx3ANkcr0BqezSNadSNzQOAJ6mvAKU00lXdcMrQcfRrd2imlEHF37IeXNZy8KcVe8K/UmI57yxKvywO3lL/QxGDJhkucf66SjslN/8w6pYUY8/EK4r6yKjdapm1wW2XpCM29X8OhxLOS+bMp89CztfLKumInjHL2l0v7wf2amVh6sCqEq/TGjO2Odg03A8uRM7O9jcsNdJvlyW/tTVeel+n9OI6ktfmJCGAG4AyWpRzLn15LxaXiXpL3zikvfDbQKn/J6QZtnUyU1+C5QMJLvNNbVHBe4Mc0FFrObbBPzzRtAxxpHQnAlKo84spzuFCtt2/IXsOCuo90Kg9mpOfdWMFpIUGz10t/EHOpl8IaflEisnx2+W32LnHAKfUJNNoTg8pbM3v87rAW3bu7EtybtN+6Xrf4nu7F3e+HXOp3SY/QWcF27Fr6jKdLfBhS2+5uvsyl0U1ZXbiclydXYWbhXoobmc0LqgJeTTGmLyly8Ofrf/hia6+cV4rHt1VSs1RCuAu18z9QYq7oYJOP4oCIxRoVtyQYmy8Bix6+X7FwdHp0cvXu0Dxn02kuc3+7vP0PHYdUy90OY2dEp4A0FvNdjeOr+t+SMdij9VRA4CfaMHdM+U5Uj2EIC/LmrcthqHysY9IYxK+1/NJwVow0va87OCuiG4zCdpVlGh9KrG9hs5KRPMkhdUUWAdd51uolOMj0GIlRs6kmJRlfOmQxfQVeUXV+V0SrXSyTu6R7fbNFqiN4rm99TBMEt+wmypmAyDW6JZUX/TineVtpPGmdMdp0nK2lzRXTmjXmjhc1kToplzdNzk84tZMdTTQ+eumE1xdi6YNn34XzEbWlB58RC9VHMqfYaJpG4QriHqBfPgEzW09aaiNsqUDgeTSHgZHBQ/iWoeHh1qY9IsZ9xVsodfb4tbFMM46eOCMOtlMBDZ+EaJ81y/E8ZpQJTndBFqPpFnxUWOAV6iX5pJfKjwlu88Hxta8quccJs/vErf11c4SvRNTvN0CaqvyTx0KCcMHM7y2iGKOkTNGC3jp8k7ezW5DUC+A6+8Oe6+Yu5Re88F38qbYrmslh5eXDBS0ZrYNvskAMbwlWyAGwWttP+GNbQvsmCMpXhR8FZABQjviu6AslxEW1ZXfARpfcsr4ruEbn5bLBqAdllBuYZ0Ws8IqiWfj5NX3/HHz8bj5Jya3qjLP1mL1OOKNnOGU01bPC3rBUHaBCWGIxnAlM74lCn+ohayf17kSwVU+e0lQc8s+fuCAAWQxa2OYUkrQeyi7K2SPcbfTS7pXFNjNyVRfMIM1AUDnbOC1mIq/T7c42LChlqJh+6221CXKER7VgHqNDldF14YnOc5qCDaLuISp/4Efbc65xPkaDHZhLNZdVYL7fcqX4AAOCsb2o74HRF6s/Lisole8goJ4vYvz3PijaZxg7cN9nDskf9g4PC0vQrRKq0p3ed5C6DhEklpfFgW/m7sviHWdR/M6+ej5Nn+3utn+89Od4llPAC8/4w6frX7j6ff/dMRo4Ovx8WX8uY5s5SfjXHkXp/9kQBPNi2IoiiILFoUhIVTunbEE86qnMUAI1o+4qC3lL+iFc1sYTI+XN8Kj4IPtjj+w/vhdgDtaQ9ZEkAzon+Z9SeI7Wdb06/5lJCKHFjhE1mm4ZbGjiEw+eRSBBrlkilMOhFXi9pv/h73gu0XehNAhn8Q0b1gsMG1BHFTc/7G3uQlRBeHdHEJLIxDVP49kSLcpicJ8Jm+5/XtfJI4uUS+KFNQsSBgp7cipLBGqgWK1NQOvjHXQswo78KU+EYi2TC3LTCQG3VTgYV9T3zkL6AqLyucpR9eHx4NAMnyaQGC/JdksFcRZJk3G0e3i2JARcCdlxPmZx7/sa7m1IaMZSv5L4evDzLIUuYX5fltKgPkSdoIhU0k0nInybEayXnRTC5TarR8PEgeJTIznccwrEkcVO5qWTMZBpAOMxoOtSInKf3l/VBrlufEBbqy1dshuBjCKLhP+wDHKVrNGDLjrA2OnByJYER9w1IKDMu1AcZsZQNbFoSz5jwyvHjvNhOEdIM5ptc8pIaWjmiKa+FrhbkZgJ0s6+dclMr5owAeqA4ED0YAP03S1D1lXOq725+Yw8DgadrHh5lwHCcimojKahns9/u4szfUAHXnpARcPB0e48WJVNsOxDGRrEZJJS/642uVWs9SNuOyDL2Ya9iW8y2iHyXp5oQbIvnLFUuL3Dr8ErzG6vh1oXX0U9XB9byTQWzbrkk7VE6BFdEvdSrjNEHCUy0k13E7cavwkxJLAlP1RMvKg04iSh70qVKaiyUhE8JpNUGF1WxVj4SdmisXEwm+CC4Rlo2WQauz6KCzGNg3LSAHIR26CT7gU3GDifG5s7MRzGOXd2NQi9BFqScWxWIba2wNbkJnjqnNwE+WASn2ENO+KK+L+XBLmQudotIWuDmfJp1VWSsH7JwlbnQoqyhyXRMpM+1MC9Msq5nJd0MxNkhu7V54CpWf0F8RMdYgvyE8qUUqyucVSwP2DsuhO+OAMZ/2PSrBlwfClhZEXoYndVsB0tC2yAFnD65st9ynIVPWwfYFrdsidg7BSJd+R28RnRqWiFkFOhjWkb7KeOgEYpbROhtgEVCw1R6xzsMVlf1oLYjBRoMf0abx8hIQyUWc0HfcRjiZrFKILobvy4CXrj9wQHovqKSid0BJmfBQn7JDBvvbfZNoTcAo6Pb5FjYlXzXVFaHOid648Lzj4lwR2k3yi7ykCUaTdRXD+fJUO0uets6CzcG2wsP91+dRS/FhHoaz/cBJjuIbL1c6uMtKmXLVYrnBban+yEjyenXFTTKJRoRGQ3ewDx521qF9Fe9bhh5sFsxQht43R+WFBeyxjgvIoYO5Wmvbxl52LbvoCFXwl1/ovaPHEA2kQ0NN7l4GCEo6EiIDY2iKd000hmvAi3xZF8+JFm/kuxtHQJFQJ9fU9EF+EKEMbRl3sroqmNAEfaPYogEPdVYR3zwtL1gcCAYNoo4KBCrhByLmaSozWXGIQ0TkIrCbNl+EMTnLklXep0CbtuMq46nKbaJusBvF1aKhM0LjJIwzV1hfMMFNZHeVPNwsNoXnezh+tzl+OPJSiGV1QRyBgncQITy3WkdAZ5TIS6eMbAN/zHp6wDV0kW2FQckfMhksHzL6fUVk6rIgNmNSpIMRtGPZIHjzl3/9N7zbGPitaIRchMaHyaWt5PE/H288Onma/mH66A/Z0z9MP/3zHzL6PfzkcUZcIW0l9sxGNHR757ZOhF5+A3kGdEL8RMB08lZgbVlEpdI34Zy20aKU6FvoJ9kXtLrKMkMP1gLZ3NFPPIh1B7OzqG455GVW1i+IH7mgz50TqvM8LBqIlESca9ph5hxxQo1Epaml9WW+XOA2XxLbuTEt5jX0l3K1azmH0nmthAUNCChq7/AwWZTvaOrtM3FeNql0HNMB0wXuotIkUwI5k+IHNPAGfBRAE3GNIOakbjaZlcRysZ6EmLH2++8LcMmeZtCPNyiePKAz84q4qIw3IL2hi0OdM/jTcpdcvV3wUgvyKQ4b7G1uu9VYb1Ny9BzGad75iRDnvyfrmQ6eTHHo6XNWF83RkhgvaBhSamEErhp6F/sZM16/oBKtGy0RlxEoZcf9DSuGFxXtKTPFfGrfgSgvJ287G8cv01kFSRodZ4i9bAeZaEsvy+TbZFYNh77363y2Yj5/Vp047r/F7+bMWlHlDVROHss6XeXv0k3rBlK+G4bV/I0e0k2IWfBwPquIS+Wfs+pic5xSg8MhVVkRkKYaaP+x1A+7xVjQL5f6ZifZzL6gu7LJBJy8+YyenwTPX9HzF/S8OR7S7nGDI50gLAVOpHWsYTor5KryqCZFOaNlo0FwnxvJZrHxNZrA4zaVo7Zp9o/k86f8Ga8f7fCbofaRLVb1pcwzP6txtb8Ja9DgxjS46/b+S2VZ+3j3d3lHBwbSR4qHQGLNbw0HMcwHkqoX+QTAbM6GDAL9oTNKrso5oXxGRSmBaFYu0rfL2wXxO8MYuFF/RwAQYPtHOjuCbtdZUz2nuz5NeZTfQtcu05GjUM7TL0fBYm4Euy2N6KIOA3SxAWQBnBHRMrwchDNn1Y3iOFaQAMw5YweVTRimI4DSuQtcPl0SBwiF6bLBQa2r2aopYqAmJFh0NvTUs9QefwgMBu+v9e2Q79S5ERpDtCQnwDqiFQoPg996vowoLjy6u41j+YwvdbVs0jQfJWe8/jmt31l8cKjUsbuHdtWwEb5hLPooBGpYCDqNrSK0KcMTO3cBMFkw8holBJ5aHCEQ607ySTqgQoNwUkZCUIHsspwSNgI6pENoGNDAKQqgsArmUIqetuOKLHEMQcINSwXfZdX5OYFZRi1SpW5uZwVN6dxBcazK58GqKNoq5/NiyRVp5jf0P5V5Rzd783PawUfJYPFuELbYVAyCblHiS1eAl0pEtxt/5X9JMS2ZexAt/Ih5DCialaEbiYZLNepJ8U50Bp4FmoEcdUI6Fl4yXZy25JjXoaHH+VXj6+gRVnoFJ5lI9WtGrYOB/83k21MjCR+54rRGPywJFEGHkH41xDlPISS0CgNY+WgtHYRddVGKLGEZZ/qwLbpqqyUDq/K8JEQLQeQFrKBEOqOCtisiqFXXsRTDutCc75z4HiW3081sTP998dlX4y+//PLrza8/3xwqqcxfMh4Iy5GnCnpE4SwmV0z/RxAS3TkIifIR9YdTayvDH0VcIg3bmhAVYy94ien5gW7C3KPmuE5rHaXtEBzYvsyj/cCpZrmZfJ0P/Q7OsVNzBuW0DsHOdlugcq13ATN1JPaJCzDAIvmAmtxJFxQdCUGuMIdFayisOg1Vf5ZLYH1uCgJs2wNWHtKMhZQn4rS2vRW8LW0+Tq7oKuVzVxy3dXNkAkEdDNggwnibCbPAT9zAqeHiXT5pZrf0TY4lMU46Pr6CnyVO7jbOvqb/vqA2UIZLZyYMVk2MNMEqFJYCAOO+ndP57qApGQCsuXQxVPMSEW56OrQE08H6G6QJq5geSGf0i38YeP8meeJOFPPXATS9yBdKF3maqIQ9BB1aqhi2Q28ePRpyDcFx/PG4PCEIqj+BSU5CXMCFeQvT9CKENBcMYS6gAxv2DI7r9eM/E0DQILmUGct9i51WejH88Fva4qf85vgKEoZUfvJQ6XjoB1C0T0a42cTG7LiVFWowumSewJHC1M6XdCKUq9v/4fDFy9cH9CKkjx0JwFWGjGY2O9iWT9+L+WLVpIR9lGqFPU0+axEstQrRgaTKqcAQPkV0eU8v/wRO22CUAagtaYe+EOBPI/gEIlf62OKPIWRBDwqFpFfBjFAVEfbNbCDBp1OeBr3javpO+t7RMYBN9M03ZcPfeoHd4DAAy1usw5JyhIaz5Md5U60ml/T9vCxmtF+A26rxhHYfJioCr7JBZ7mBAYLVFl0dYIa/dMHsd3a6E3cH1887vFuNoAjXiEpPBB1Yd7gGUpARwdqLapsQlIWMJD+A3CAUM/jL5wQNJgnjfXU9MLHFpYmG2VUZCi6yAFlBkM7yPHfWxoNaSpVq+2PGmJufZeOvNxNYUNYBAFS8SvW1xHj81Vfjr78Yf7b51ci6FyxeyBhkqiJ9mZawsa4Fr+fQNNiOg+Y5uwUPI6Qi9cb7mrkFgKnX6d/v/9OhN48NtxACVC/qspMcXICQ4mdbKtBpxRR2k7XIYY9PlAO4YpkVtiU+GI3RuhFCfxxciljh6obMVu3fqaHHad0Q0X7Kxu3uXTGfyptDFeieQtXQ2DtXYZ/LnXjC7/mL/ZfPTo9eHL3cF2PfdjdbsABXIxR+Rc3E3YYlCjaqbw+CSjhBM04LF7HW22o462Rf225/5y54GZ9DjHA6q2jlTm3DYbVLB/xqdZUsy5p7enZKPRNfQlT0Mhzr4bPEf6CCr8r56XTFQHnOPVtD9pLK7GEkk1VD9MqpYmrY9/q3hr/RXv7ulBCKNJW/46boWQeP3s7cWrumrE+/5q69g4omdHreP+k5Pib8EWv3blFM6HSeGu3Dgyj0raOIIgcKw3KeDcDhZxYnqQdG3D0P6BjlvWtzXRFrPLkXWwHx7Yh6Viywe40zQHH6WbFPYLJdJOIHjn5azXspJWe3XKe+vVFy+acYPc5ZzQg7dJC0/CNTi5ktYrAFOs57QK3wr65lDx2YwdihjviOGwWvj0K+q22EvTurqhlsydtAmsbalQQrgEYvTBdxdzAPeprsLpf5bYaVTX9R2hkk+/tRkp6OkpK3rSSIInWdAsVPx+kFW2J9NgEpxMCAf8oS8eoGonugJ0dj0dRpSE2rGxWiKNMa6IkzfcdG4LEpQ/glC/eyiRrwQ4pVnKgRS0K0d+qkxfTqF4HOrW/hVmOid842AtQBKnDGDOlhNrWZfpJGBBl4ruDwasFo4p4cadcdiRBkaNvr+OefWDRmtmAGLVc1JIlnt13DguVqRog6hu3JX/77/0lz/CYJAbpIHOC1wlZnQTn8a++pYPuCyhAYYtQQvLUYmY+5nncyNH2X1+kLPTxLr6P7et29rtd339brYWS7dTZO0NoVppa10SbMgjb7vjPKhLHyPPjaRpj0vR5H363ROmyUUWQ0pqaC8PLv7mNBsiw7PhtjiOgIrZ5k5+WMsELaAkNDqrJwIux4AepFPheVV20qL3CMwrbvmBGh/GKFG362WMuxsJZz5Sex26L0U3aSiTJleZfoxlovmWuyHvSJe3n0iGCRl6ZzBRmKDeP9tukzMAz1vuBhsrmiSBm51zZQriCPuWsCMljHDH+Dfd7Atgx9PzpY6BipfO4klOAJAzaa+gIzT3WJbyKA8BbTigjIX5K53qlRYJLRM2Y2k5T+eY7vFXgIKZK4696tfDZmYNh5vYk2cQDSVPa+AYikK8GTFfj5Da6ATj7uT4FGt7u6v7u6v7s67m4Ht+NR2F+sPTkSyoNdH2oCWA1bQpZzZ0MOLuKcSKpR8pAJG2fFBX5nClNvsQQRiubgYZY8/8v/8i/OYDwy2DbTDyNfVVSkYPm8mkGdwTRS2YCRKWcqKYNqowA5tsiBxme3KpFU8qwjSIM0qcnfMtNzl1uymRvIcRGGShgp5fbE51RUKCKAC61DjGzzGhfcQ/VQ4FUpG/i3uAVgg2V4s1SsX9LuhfmioY1CDwetX8t15pfbvXaI3KiuWh2NoalcAzDcJcSgfB4buKPrI/a1YPaQveTYKGWqLW98U3678Q1xgd8mcAekzqvaL7qJAtlJ59xEzWJIo4L4c7b6pBnPRSBv2JD4WRySkerf3PLQMmK8LOF2vmFPxobAA5cwb4dUiTuR7MCVEy5SObaxLKAvtQbKuVOWqbcKS0XEuYPLs3IEjhDcDkEFQno0hs+shWVxlZdz70t15bUFKPBmJQYPBwHoz0NlC+GWAwjEvvZiOu57R/7+VMLUAhK9A5WC8OtMQDSAAFWWVzIEApq49J8NuZvgiyNQj9HYSWQxxMYhLxhN0aZHonkvkfGybieeMeuWRBwmr70YP7bf4Ea8ZuR+c45oDb9jtVlaTuE/MS1kfJ/QM7iuDci35PV2SFax0ZlSVVorJKgaJUMj+cToHjprxLhMF0EahR9fE0QIsAHLaOnzrLjIJ7d+YYhStduwcVVNi8EwM30IXLHbTBWTdNRFyrKBDZEmOYFUZBC64nHQ5VQxqJgIegEp9ay3mJfMLR6EOlYlbPAAVityMNLuCOievUWJA9Fq80kEgnfjL4VMHPOOo6yQ+Odjtv2jaZtQSdympfTmwLNcsKJzdyjF8Y8ncxYRk0zzH49P1heDSIQLQd4esHs6Np4iCt87PkXP6nfBIrBjW6ieccEXuOcrD2d4MuJ5xot1M6eJuHXF4yY/6gdf5Bus/9PkWB4Dta3UAdCA/D5aRR29ny4QhZ1CLDm1tB18bqro46aq/PvXD5jo/v31RVt9HwjN579Hnc/DXdN7JWfLNdod8cH2vd1FLHn70AdLBG6l9zM6hFdk+2PQJ9y8POzom6CDpd7yL9nqaY87C+8pk8Yp77FwcILq7fjI0zfKmTPhLK9YmxOwxwlASNonG8UOW92Dvjo9Ute8SXMaZ+u7SWDP3PDOADVc62f9IwoanbhG97UtejkdRpbmz7vIXfE6aFTx/JvA+vGciTnl+1W4tsZjDj7wxTQ2+2Hs+mY1j+QZjDxYcMPY46pidNR+kQmWF4FV/IWpN3iQBZYHXNr059zTQNESILgbRyqMaLqUwSwBI3CclgAgTxN+3pI/BP3/8i//Kzs4ASAOsz9WBDcGTIANgp7L+pXSI8ERCzuEJbD16E7eUuGUO3zciQlhvq9uxBZL5wDb6yw7M39nYwq6xLMjAY24K+d1CbN51E7HgScKNzVsy1paJNXZUCkyx3d7ZjzgwAGgc3rGPM7or7HgujTpuepElKvWFs4j/pq/qo38lvb3VIo+YuO8MfPbHQpmdzZTIgZ9BNoVIrq5FK2EC2rw86pY3h4WM/pULalmOlDHkQ2QOdUShENIF2kbDlIrGeXvUQmvCM+BCEPgGCij+9l0Hr6djk9YL0uT3YUN3pRtUzheE/yNlNcLSfk+j5uQohuJPVlI1n08fQbnXmeeG0rhIsowIg2VImzTStwSnXl2F3ZM4dP4OVOWhUndNTTUiG11aMy1eDfgQ60U/7dm7dahyDL26a2b7qa7/ZVxmzsZj/cpcDKT7+pHzLaSAvpUt6kqiK59xhKaHMQuwOFLEKlruWZoDuG2JBYHbLQu3f/mTqrEFpdWJuUFNTpqaCtsLzIBGzKlNUTBR7TGV7Q1wD5aotOkSm9ci/rcHV6X0PngtlqDU8WLsJyCtEPSL3EBCrCLdL4WebmkC2nWIg8e4AV7nNBfwx25e7EpL84ckNQTxEMS/D03qoJPd8Hhr5Tx6RfFno//Kuo5ErTKRQf1+WLK/mAjhSp1xNGyk0+L5JZK8sY4FJpZyK/QpLBdUhLMOv3o0LjWIceLOswsJlSGnVQvO3/f3vcwS0NMw/c6wnDo3wfYTZgv0c6lAkfPx0zy0x+g70jkTkN6Ex6CmBfjip48HCtzIa/dd2MttBvHWXB5YywwkmGg7Elk9J6LGiX+QdhOP5cHPbcPirwePUFoA9hfixfGJj7qgwjtlzjDtm51zLuMEv+wftTuUuqgD+8Ya1i2Bxpgd5gZ6XyXm4Ym+Uh1RO10D7hRXLVvHT3BJuMhRczxaRDdQM0MYYOiJzUSYSJmhCfAGE4IMqkDu0NI60Fj8TK2cDQBinmzXuzy8fj5XkGNYWO1nvsgUQnAHBO01RX1UEXgIeQAUaDF9VHptcrovF8LgFsPodtZ4K03MsFaKp4C6YIKLWLCfHHMoBQh/3BCB13AiYn0i0K6sHLEcMEqtRHtKOmws7gCUfEA7Y1aqGs4ukdcFA59clkuDFITdCdKo5YgNtE+nK2ahtV/DGgj0pZIHMBV6vniAhECt7VwyzSeW/bf1GaOAwnBXGNvVkoEUSCOsklu4QVezM6ztaAZA09DnCC7SDCRfg2e91yoIEKRk1aXjbPwMHG2UtDMapZN5nYa/bXlUqMYc/GmMtS2J0PcDqDzy033klnBgVmoQewOhMQTQcl1U+kZVSDt0ePylP+0+Eu8YcPR1H4KQKNRJAzABCEOw7O+dj0Da/WuasU8/ryCJItvTd0SlQh0cYxqr8wEFj6P2oIT59jzhG0MBk+tHxwP3u6NhbgXCwIZts6mSnefwtgsjsLGhix6Khuni8uwRgT+Y05hqDafKtMbvAhsw50OBf3zzkI5mOoJGQYbpC+JjNCZ8PX4j3/3aEFKMYZyH53SaWA7y+3r/sft676HK+XtQp6/EKu644HnNrD/Nftkltdlc4vHeXGRs0WZRnkI7fMOfnz13f6bv7kFIA5hjyEZlw3t5TjebWQXhzc9VnASGFcN3qxW276No+qutWKTaLwdwzVZjcAuekFIXP3fmZ0VeWmMjdnrFJ4wKBKYwVpZDY+S1/ArSW9GvYaPSt578wu6H+6dLz7sGkS6Kn22kMHHZ7asVxJwNDAebBkjbo7WWAx+drdh4HiNMaBox2OPkZuOVbRQSKJwUMEEMx+McdhKul4XqwnfJDoBgXXwD5GKBm2FapoQks3yswJszLqIlvwdqJF/ZBMCifWB8H0D2u2IY5zL+3VNwaIALaFY2/sMU7Q+8sWCdnbvspxNfRiJ2PdtTRe6INnksionAhZrxvWMGUqsLGN5td2npYpt6lXvvSNLFn+TECpu9bwhuesu8lM8VhcCzPMEiLlT2AiU6o4lk0oYc+W4ejXur9Y68IXrVzkVFRvAhOFhHOwiQI1q0dbelNNioFYzvFxYOwTqxUd1fdW9724adL+OuOXPPR4AtMJLASx1asfagrGxw/U1USoaiAd+sRaqQQTeZ0VzU6iAvbmpjLvYFjeEZYHo4PMJW+MVkxzx2fNkWjQ5268XsrxwJAedcIaISSrAL8WAATG545DLyU7H3g2jIzaF3a6dUrR24jKJmwxccpjpRF7yZUMc6tYrbsAQUiCC09tHbbkK3xezRdwEvwlb2LikNyzW4CFkdCz3cwQ3U11utiywiulQbAVkEuKb+f3RK8TtiyVBv0r3bIaWVlIvV+1Hw/jEa+JjfShu4VrsM3Ji6/VXZ1pe496YbLp1ukPBpi8UCbAdGIjGhbB20OWs7/jyM/Sr5VqXVBvm4OkEnH5oBYZhMPWDnn4LQxWML7xk2sEwHl4+ZRULoTyirPvwR8AZGej/ELyznUiT0RC4AWdvIpBZoAEIZy5jo/NOPU1+IdvAIGj/5f7e0SAMJtUGSwb2biTcUr9HpnwcKsBqe3PJV5kJa3rpbyZ+Vw5MUUHHrkVWI7XZjMhsI9byLlx0B7qJgOmcKi/LySDGNq7rwebgfhxF73cbWgliUYt0kC/LfEMgh436zmkydY9DI/vd4aCX99+yZQ9h8LciCuRLyHUjzOIAnXYphd6rsvTzXwazZ9+re+RYGsPOgteIL9Y0r4FGwMBmCIDHAgApClwk6i3E36Irup08ROCth9A1c3wjRBnNovUnRFiqajTYglESS6C9YW0ELcMocdjIaBdHJqS+b/VZNhJtIVfUkQ20AB9iGdU6uTUf13hjUHC43Rquym7WjUtkLTwyLDJMJtG4vta30XhDQY4OLhoEVegMQuV39y0OisWd4Y074M1dBFyJBpo2rTZoqoG2K4NMdXthvslb3B4qx6b9kPXi49VZL34bX1FwNm/tU2t4OLHuk948laPE94PjzMVTURKQqg579gEFw8lhHwln7V9Tzy+RcoJoEKJjcKMIeFmkykADQoVZA7L2ANJoeDGJzpisavh8yqMwAnh+HwyAp/hRA5hSOwRl7xwDLkHXcHBb/FAknvQ0HkkLJMMqtviQsykF29RNc7nhtp6BgkOqoU01RABHkZSrzxDW22OLaTRr2SU+jOraWdIY0yY8qhZpggbFBFNgaSjZ9XxqIFPpWoEHMw6705aHbexdXYkD3YfcG1e6e3n8p2iZ8Vrl0lH1+Cr9WKtAbmCj80XvP3c6M0cEHWZqrv1xp6tnzdwo2ru0ukAg/btp6wVT1q5ka13Ywcl9HEQlRTjhX/iNj4q1gtf0IXVf2s2Axeyewh52kKS5s9gBbOk7IKfjEt7RAGycGeK0K4n77eBl3qoRWY4qYDduSENVxR3HNcSofl0FyLk9+FIEnQ4C1wi91rXTF5gE+E/FsoIk4v/7354/ZnH9KOkoiEMyw8Gm37Q6MwGcAZD1VvSavkMHIG6zs+K6mIUhWq9W8Op7NykkL1hL0atAofdtuCWsUYqqyF46+bXbz6hWU0V1eDPbVdxKuBV4E3tKtBbCGQwGk+QkIxXwwzk7Dbc1w9p/71sbsFOc9UyzlY0sHt6w1UbfpO9uYN0d8waJgWniefnHUmCTl3JNODnUfTBFSvXAEzS5oefbFRNg0vupBYU5DZSNgbaHExtU5ovEqofO2FWn7vUP3RKqJ4chapBUyAK9B6MxGuq7GVK5CSkCEAuH0R7VHmwamQDXlGg1o95b9kHiYDZFlGKlB1WvERrIeJxM0GVmqtUeSxPVuFg4gStO4Ap1UTArswU3nFkRRn5S3yzWvkFqFqEWi05yDy/JxdoYQF5G54JfDaJQEZAh0Al2oRyGWI2UCyobJWRb39JwqRYyVNXbfadWi7XJg1B7N/Cl9NT2KPcUDfYMTou1hgeBZslx/tcuqZTAIPVniwtYCbukqJ/jJh/pYV1LekoxIQF8DX/Ij0IAgh4yryPPp38EtAfJ2KY7rjQO+X0kr5bsTmTDRubLyFqHH2OCIxh+TFvwh2DVer5yNbcjOGha+A5136sgWkaS/nbIKUxYEoWwmWNJhRa11FEQGvJtRdRYM5BI5zUIXiBYTL0Sl0weheoltpLj44Hef1ZQ6s+TUXIcaMxCraSLmnFyEtCc8Ti6as1XrXgfMFMddKa/RvnZ1/sdPfepRt06uvoZ37JqVR+WZzOW3YYnGg5/IW3gYpcJpA7jrHEGwvZMYn2tj1Ai887YlcCxfcqB906oT2sbvPYB2+5ouNXk3WriAx/tZM2Z7R2nk+OiCVxTzYBWrRoRjGqtPninrUT6iegKSh1Ag44x9816fqBPG9LfQWEUz/tAn/1grU7DIon+ig7EUo15TUgrzPgjUAZdFpO3Gt3jzu6Tb0330lEFeXT6oK9tqbVaTAne/qCN70khHVOkvFtXMFJOf6TSiEioVY3iDx6kPQFS8HEogUgkBaabqGt0Ljbi07JGQE+slw6BaqH6dlTaq/furONsceD2s4G8akt4HNxfp0uVHN9viTKkrlAaunq4OH9UN7/WkyPWu8GyslcdNlJ1R2oJMu72A70S9QbVSVtISPfdiUP8DMS3isZvA+gOn1UnI1XqD7scvmjpuaFQc+85Zv3kV5WPHfSp3hMwPe619HFWPR2k2ovd/lYGOT6EGq88hlqaFyRG26M4G0qxB8GUvNq/NVhOFIHNYr9qw/zSgKvTh4ZbFQ2LDjtyg2CH1Rjzzi02mdzQLC57Nsu+mLYMJnS2JjFduWRd2S/v7+ETqJjqF+PghKmcZokCFCiONHcEfI/AhXH+CA0s0e/gMt1KGNiKyBzOlkhGypENHhax6d5Dn7uegN8twjr4lucw90Wm2gDHPWhf2NBBB0bdkQsSsHeANMX7Scye72hnmPyMFMWvyglB0Lx+m5qIEvUUSg+5YSS6KVT+6M6ArPsHOkn1gIZyKiDUDItNVMJmxAQ3W/befWdGLoO1HRmYhzsbstV04yUKk/LXsq0T/CM237IvkhNCI1FEE0YeaTHNaIs8hyMnjYRcwadxNjmKf9PB98clHDbAdO90hS822Of0Pb4KPoBjy2vlQZQI7YEE/+EwLxKpqS9oUicPnaJ5UbZ0QiyJ+2hPO/EQ2dQaRIAfqyr6ZbroE0KXDc0uyqnsMVO2uhZJrSUe5WyzyBJt/rdZTzs+eAtaOlRjZt/YusoQAnGV3SQksFuSEDnwKm/KkFVEjYudc3Cf6AZGILuBNNwyBjiODXDF88+4gDhnrciwgMq8nCGS7Cnn7Nv5m6cTw7eaM4Ap12UxT6Xd2LZDto6Nnvt3KpYP8sn24sFV/WFqGirXVdDgZSQNuLiseO9UT8uVumoZbMpf/uX/qDUFxMDPJ+EqH6UWnMyKfAmo17rNZhnbmfp234U1eSYrMWfSRF9VIeLprl0V04hF0GhaPaog2UlVNw8gs1wBpTkI/T44TqEb2WU+5TCArEfEeW6jQTEj0iyxOYNCDfy31YlMrll+3PFXxuXGJRLywseYKp1bgkle+TlMxX0gdSOG4AggPfId6vQOOTj4dA2LM1UjJKN5YHfu44ErPGMCtkUyGefSxluRBxuzDV3JU9uAf/CGIwOyQOgv//3/YlE099krQjLfTqK1f8vuVVGVWFYUlE3+x/+eHD6TvExsZt87+TumqnFJUVUX8GnCOV6XYrfP+6Ofhs4LQ95ys30EpC0jbPyN+mZfDrXyU3eOME6ILalytfFS+iZeyi8JuhgFCivgd6b545yMkhG9hMAyi1M6xH+/NyOQLfHwNYPVtMDtsTzm3s8i03jZSi+ozi7Z8wRfLX5B1n7tg2LDre5HGMOGgZuJN4V//8WSCFXLlCfaM5vJCAgB6c9cvu0gcJUkrK7WZBc9701K15PBpX3x40DodP3U7BT4aENpk4GSkg/4u/Oak+IdRGKFGU8SjcY+i/MwbP4HRChBDpXkcZJqZX78FKIp5DcvJmzb5tDDKDFvPTVwjEzd7kJCbMoTIxr3Gx+iJCh5I05/fFMH/nOUD4tTx2Sfc+gReA9uuGpBPRNCqjWcrGMsSYr8zCLvVOcH6txV+SQPZWG8Lnvwwa6ugeQ0TD0TOpt2nE/jPvWtdRl5sI66XqvUYVvb2unXiKh21ASLkvCNHqzhr91oF1F70LfTcZ95sy5Yg0ToCEO2rSk2zr7oOwIty5jB2tMQ+UWdWyDfpB1d5mmSduKqu9JAaMdqbxJ44w81voMFJhStsU9mqUJeDl7IzmJhnOoQvsRyzziomllbYBwmE+2VItwl6Arsc5jp7fXtdCvhDXVUGyy/tyILnrYP6MBT9Jj0VXVddpaEkaam3uDkigODW2E/rETGC/ESjNwCGfc650ZkXxj4PXbrSfwhZi8kYytzYIsS5WfkjMMliDjMXxlnrlfYGcRRHiU9pjaRHQ2x8KG9S1cW80mPJKavW1AzpcnTzspqVl2UIKGQpG5dL9J2KJZcEAi6LK7y0+tiWYstlKMucIefbFMR0Yb3FRnIrw0J0fAkG2djKIAz+L/RWJx7SlSHhVj+myP6qML0/HzMTd1RRsQPAxu/ml13aKKJJ8Q2LNmP7i966BGoR8vixeankse+Z9q01lqMYFSrzUjsHrVs8QF2wrjrHyPIlpBVKsL6CMeSfp3KMdc/CVxONMDA2z7hd+xh5zikm+yNJCpt+X04SS32dM+5eXjpCFVhz2gv0tXFojKR9NkKfqiIXe+JOroOkx6pJ9qTu3YcKy0/QDodqd6cq2jSjsgNoYP5ZYbT9sKl8Pt7t5reDdQDqE60736/UXUqDcfXydIS+ZkO2/0NBn2TU78CdXHFod2n7SxpPY4H2PxBSxPa5xzsu/koV+F7HYNP5ErIiT1+y44wJyd+k3r8ck2pZ1BtC6mZ3EtAhC0FLsgZEsK/rQ5EHCW/b0OLrT4QMkoOW3BxqwdWcm4RCMSucrrvzn33tIbR4NElcTmQwZ1Wi2KpBb7d4RwiNJr8ojhdVLNycsv5Q1CvKU6NYxklRmCeTstlYUk+6GDTpZCsJldVxamKObIAUb6B0jt05/XqYnXtNH2h8EldJWKUnOg37bO1yLjGadiMFzBoRPOGeLhfeov630/DzITudZCzIpnkzeQySU+H97XVCYvO4LOjIvbSq6kK+kKnQEeLyEedTxRutZWGITKGDb4ZOQryy9vAqmoGcZeviIuVYGR9TO+gZnv2X01NhJSDI6ZC4kJclTq66hBKd4iwFu2xVtO9pvac+dG7kW63Vg9iZxFuD0YnBEZ8y/J0WsA0k1H8fY3doUdfM572EXzaTiByd/mRhih7MgzimnV3mGNGDDmc+EFlXh16uiMH5b4D7qv4WMPlrPCGGnwEtgOdvm5czIJwnaf8R6i2LaiJ5bUoupnV862wqDmy6UBZDbf4OBG5s98CyMGrZVMntDTVsubLgEDkTDCyektTCiGgG/zQYL5UW2gZy+4ivA3dWvFec+tg56LnstNY9jmFJDGTNYHgpIYj3crLuljdL8ozTmOtGZZRSmT3WfJ6DutT1wJ/n1dSlcPf3HDoeG2Kwz+yz8hNfmtjSHen1zkhkClC3iBbcS1WUjyGNyv2iRypwXaFAO83kLOxmSILtS8tyqgL/p4FIF4XFXkcqZF8OefkTPIIibv91ix5BMCLGyiC0ghRqNuuTktnq5Ip6UHyPtqXbc1+TMWH2jR7RskbkfVrFaSiQHDTrSRufSSG+MhUhQ4Goq7dtlU7uqmEideBsyE/tKoXxNRuScx/ga/LtxZslQOmhnstPQebTYyy8CkjD5N1q4PVgGqDkDWT22h/pA0FY5dKd64QanpbrHNNgemXati/Tp1O7looZSlo23ms7YrSoo7SDkd7nB83Bm2lZxQW4SxIdtRNdbQejwYtqUty4pL5BCmIW8mJXS6fB1Eun57wZf1pffow59oYqUy0t+fE1i0PdMiSPzCaGFBwmHMz/oj8rUMxE29LLJ77MErsoXKGGBIWmsLlhCcQVK/gjCWQMIh0moVyytYAW/4/HBCwk80gGnWrRnda7Sa/gRic59Xj5USUUpj1V4gmP0uOv6ZTFH3u5gfNhqMIfsRcUP6OmXBz39w9zbsrI5FYexXUbcuvQeD/yEsQhBCW5QFsO3POkV54NxdhXHSgnfzsTprObB7n6uXfQgr3Rc6oliVh7ny2XsohZk+RnMNXciKNBXFXJYz61ggzRi5GhLx86qQVUK9uec6WfrbMroqfV+WSFa0qa7CUV85+8MFfIV4I4y8EBok6H5eLIvzqw8dKrAgXHdkK2JkOSpxzeMfekJIuU0g/YWnzlwjfD1zcVtjX9CRHlYSoOKfhFGz5H0WJQ4WxgFBYUxLpRc0GreVREaDqyUPjG2FLLFNA1lQvq5tiuUc3JI0lJNwEFeZzytFnOcxCa6ihgMNGPJCYO+CWYRs3Au1AVOZEhpVZyAfBxUGUp3tNxSSr3R2xg70a0AUQDmNKni8/JgOJP8+w80hcXGDRo6z1Dl8TB3ggYlz92so6k6yr1YoeTBPQMXRzY3yE63nL/aZQYzieIherozmGoa1wAXj1OLCVK00Do6vic3x8hKdyPJhSEqNmGfNDMg45LCzP0vvYm7fHnVyf9aaz4TSbx58yF+GtNLdaGn1Px4IFqKvk08dODukiaT3glCFGrdK46rinUbLGxd0uidmwtIe4LpPHoGUpqr7xwruUTduILDIwY9Q2qVwlTSmt6B2UuZh3iP5cEJtbUpE1tWyt89atWp/ixH3tyXCCb2sTnLiPH5rfJKqAaKISvzoMXx6brT1I8yCHCP3Kh+0t9fMyP2Q0br/b/rgOAEZ2JwaqPUGRR9FFIXuFWED1jRzvlIbbioCKN4iavVzNtxQkRTFdo9AW7pjStJLOOe3MI/YJdifUaauR7UsCd3LwdpzG2E9a7JAYMeGYiXIjLiKetj1N+jn6BjR+qM/Y4SPIwvjlW44r/rdZNRppp6E8aOJAKzsomR6HpsgHJydrlh4Ys98XwIXw7phorLHZIJx0va6t0MRTss+p9KYt7xSLWOLRrmkvChbOsPeayYTEd83l26NJc4bEplosium2a9ZZScGKnVgfDlcXfglGExv29t8OSO7lbCxmq2U+S3vmP1KfWT0WZ6vGZ9FDzbWLI6doWajpX5YcVDrpajWbslpPpp1x/IN3Rgqsb/BveOxCe0EDSUrTg9NPNao1Dasf3ljIs78GYvk+giP8vgUlJw5KTvFrIr8Q2rgDLz00Nyt97b03fkF0Itq2aneDzZFkvAx5sR44mnZgzLy1KUH86f5rPOE0Buvn+esmxxwsWyN6+DnphBMIoZ6bG2sRUTOYYB+MCreCS90NqUQjyaDp2j/F+QBoEAxiXKoeFR+LhJCTTi1QPEp91TrhpuDiI37nGom9Zpy0SVBmIjEdrvJbvsJEztkddqu0j/MSru587RoxHJ/z6vhD20NbujXrOQrhKe/Hph83xSTlqNdioJs+woe2XdDwbmj2V65ELzjgqHkqzeXgLhD05nCPThY0tJkmdvJuNBPeJ3aTkvATMAuHPw8rdTWyhCVNYzRBaCWl/RTR5xYEpJeaRvhmycHtTXSsRIVGtYghlkgs7k2gBh2O+sB006p1sq3Vs/LikrMn6eHttgVscyRuzkAbooWUejOIBs6DRWGNAjZKlFVRKJM2nfqB81k78qsiZ3kED5s92jsb0x8RJW4gRnL9RThI/8B0KirkMAs5NrI0CqNc0qGZNKslB8jR+h+7IGEWOXPEbtFDR04DwYW3RP1T+RhJTXvH4HKjdt4cHIUXK3Oxxzk0JC4c01KSRI7qZx3CRxIT3SXCHN0hwCTyoS7/5GMQBEmFU1aZ9WzUDIFaJS5NsiA6Yrpio1yX4A/p/1Lh8JLI58VcLZ6Pk51ubl1OAxEnJgk+ISj/0ygvXpRXiAlsJuM4Qx3kZ99JFpFOKj18eyMpQ8JvE+ILh93RYm92ax5x3spW4fNUdAms52MzSubrwGv8KAmjaIkx5jg+7pL1AgtrhzqgWX2rdHBljwKCFYeaSwik/8an0OgQLEGxVtcKtNeERYR0SXyD1ifH5fiI9+W/9QMYfqD8xue495KrYXfZsec+5fSD1O0ePfnlEyeyD0qi29q/73q4zZCPO1u/d35o63ePy7j9O+vfu6BQ0PE9O9cT3YxjEt69CBLsMpa8OPFbOJJfsYu9vEZnQ3FRe2/SUQ+hY9LueOF9E30Lv6XONEHmdYbM523I7HFGDyG+nq+Jcks5dBHCga78QZTdhsGJBOpHoZZM5m+PElwstU5EsvYO/Qoovi6XG94QDxR0kMiWexGRgvxATOTSBAW2WTEg9aKiNvBd02wfNL5XXPSr+XUveH8fg7bOgeqJduDMyD7YpPOO8BAnokrqaKkC32v/7hsQqOvVScSHqL1kHSRBsHgkgSKp1WPWNxEeQVgomlVrfGtacLEEg5KtZjSET6x4WrOwbp71opiU57d0YRNW7jeA+7RB4kvnVqE1Yx8244Pmfu/svvkVI9cxsI23KhPr0P7gKogkdseGxe6c7eG2v/YNtBPbxA3RRe4ywwEbXDCeUG3ER1DHzZojYjn7bZg3OdrYmhgqT07uugrpWtsDVyw0GgivjA7tVyhi11hOmHOu7VVoLdAPKO5XiMvUXbiS0C4neinmOd5Sw0+0d/LeDObXzNudSprS2Syfv22dyF6jnTtMdvwp7fn4bV9kXxf0DELKx+trt3hIf9NxHVqeFZF7dAsYfOS9cuNry7yuOlHwGBd12kjECYy55MtqKTguD/2Qfdg3cyPzfUrlLDnSCJ8E9SvYgs9uQaRVs+uQ6V579VuBfeLIaKKv6zFFjm0N1BgnNmUeWfa63PkiOSsB18U99sHhRXggHkRlzX/vqWlDalsmc2Rdjkfl+vDXXGy9YJGM8D7UVys9IHwXWKqx02u2HEwrKLrjDId3J2+D4GdhHDRv1uvvABGsrwPLywUtJE32bFZcbQWytUGdXCCKBZ0QDhehgimOAcJKnHO2pmXf+BIyXGpqWkq42XwG96Bb2BPN8nKehVdbjS9N8BJANRsQQbbu8DOpN1xrvtlrwion9E4L1siS4EpSmxxzvGuvbzJO5PEoeey4FHj2SxRqfi/U0R+mjywoxFO8/Wf16Z8QRq7ZM0DD5j4+6YHrd6+AGYQO2xYLGHVWE02balbXRdYgHZQtDpswfIQ56dp1i41Jf+PvduhsW727x9P2P7GR7zGnnSnFUeCEw3dV7zTTVDs5jZayS9OJ8yRmurIun8DeA54WsC3pSXRjbdF9Rahga5MuxZvibFXOVEXJIUcsPomcejVcYVlkbk7F9HUOxWiioQqmgcH6VQWjprNbjiZCNMscpNk5aMjQerma/4A994Zj3jkA4ARpgHCjeUIuLw2CLuErcLaEqHM5wDKOXMK+2WzdnVXn57QOP9D1nTcO7sGIRqbN5vOwlKpLuJjm96S4uWTQhkpYE3acaW8KR7GibeGjNZRAMuVc0wbJdCUCkdo7R6BQgwy1wKCcUjbnEfvoY6s9SorMHVR42dNRxQ/pcMTZyyy/vP1mjwrn2MZg1cVA6mRy84CzL8kbV4zNkXJLNZZGVueRL/oHpWZR5/FOchZ5f1cEo17f9WDursT6mEVamrjbIBsUVwrsr4JNryQU4733XvDpBoPgURK9M/BC1xdMlmWda1me5RykuWUBGw2Ful/eBgfwl5iI4892IDm5fTllW0X5oE4IRM/8GZSRvHSglC0WyjlnmUA9Dhh+f73tFv1pd/uTqPJwpPIIC0m57lK3m9O4cfdESuc7Fx4b6VsSLDD5qoiRvcvdLqknebhD2lYMXpfivnpHo/msWEpuTXHzGnTWRbOT3xfTG8Vapzta7yh6lQS6QOkWOS9V7JJ2asg9jksFLcyKxoxKVWh6cwkfLE7IkSHuog5bw0MDAD3+Z1nHlNfkz7qcw08eC8JeU9PvGZ0/7XJN0WB4KMADY+DZibIZmkR/kg5Cn78RUk+wtkXArlliC+R9YLAVp/WMbiiVU6Jsq7eCA/4MgTsHw5Rwzh3qI1txUMPCGjvoBILCO71Jo1JmW4mNdnShgFDlqQnViWsnj7ZhZtTAWF6OHKPJKyIYr1Uj6J2ZCPkIhbIttHKTvy04aKAFAmFzldZN+OPq6sMynqFgF0nwWwln77d2A2/tYw/+4M+wEAknrwEVB4lbaH/IpO/oHsZVo5bYiYGlu8QvbhrTITVEpEvvWQcTVGK1pJSstSjUg6oPiLZlFBwGEe5y8dZw7w/TZyBawm32Y7T2EnWRG8+oH7+ZO94Syfk03diEwOfsxbypfl8WNyDF4ToJ93O24WUbn8v8uhRXdvYvVzI8zE7WyuaL+d4VXvdYqSV3FEzMJhOg9ZQ4/jbwYIEhVGsLVdxN2BEA0rnB/tbyL4vmZfUcU+hwbKsD+GLecSdWkniaCq2JutVqrpKQKiAlXhFlOdzuJSJsdEMuD5xnZtgxis9qCIxTaMmHrIDgz8Hu9DXOTYrKSkLt2kTvmOaslGnejQC7CO4X1PEBHaWN+xFe8l6WNMrkWgasc/fcoXxb+qNe0kxHGg8feMC+D+PbqamWxS8VxuTJ+MmX2eZm9mScfSGKPdrSqbBXs4rD7iKEb/LiWZaAOhcji7xOHqK55OBhkgosdu1WS2Sd866lGOKnCTXA5h5VBWVCLXH+uOmaFXnFO3GP5mzN+PbiWc09LsVXBQziCsmcap7HeLw5CkQnYr90y27BnQh61ATOqfi3aKzkB/KgKzdw2VdKr88TXwdWTb4+t8SmWoHlYj7e/+esOoeujgPL8MowxC9FW7aVSDiAIBwYMYt6kXRcdJXbQ8WVZuHd00Rd0tslAI5Z0jjX8FVDhstu8j+vCEs1t4ikqbMPhIUu3DCkHdF6HJ+EzOPPvCba1PGcAZpmuRcHNxi6ceRJtDFKzmf5RY87XXBRmR4FTPw5+4e9Uy7Pbbmnk0iEhVcsvREAwbVVninQ4fHx9p9HJ4+HIqtAcXNPgBzyzf7e6zfPXhz87vT5y93fHTJTEpShWtyn+CMH7z0fy4kW/2uxrMzxylcYHL3Z3ds//W7/5eufTl+9OPjx8PT1wX6PTMUmwWK/56jfbWzv+/29v4fcGKOAKtQCa+RSbVDL0sFnblJdXeUiNhHnfRGY8IIQaHlbBLmiqB53EkKLIL4Tf3PBFai7jb/yPzGuE/mmmm9JOINUj5zEyoUDGUcGsVxgXm0s41bAtatB23kJS8CIWwFRLGRZIeqCSU54edKHi6quS4CtBcLqEBSrOMiOJWAjqMe+Nk82xw9Zmm/GgAMLwCl7JvuBPLzofwvmJ7WIHnfkmIlLInZLzIhhJiW7yMsZtHJ6+P3rN0dM1A8wKKZsl9Xq4nIg2cb7xgvSQNevmF4UKPkGKtPw3Ui5SLHTGmyZwFwe4fzEMLQmxMMuVUQRcRc1ih6u+yaYA4t/FGzlOQcXwEKxyzVL+lNodef57LYua9nNwGxaK+K8wgdxGK2ufHx9HsFAl6mOgU4QEN096NlxsDhlK/hJJqYg0BABFvpXArGYyzcR3Z1Tc8aw7Hh2qaI6keIzhpEwGcGpXufnoTJN89Mwg/NLjXXzacKLwta3GgiDb4GUyzRsh8CNTmhHvL6gWwCiYCS5sZksFGyqkXLx0442zWfOVClMPcQ6mO1ewx1hZLonKxJuik6MpXMpTAxt67SI96CwZRUcwCGf5V5plREvM11d+mPuPuof9j4aCjfwguV9E0nggy2NVUqTTCAXBiQ/HW2sgSbYaoYpOgxS0Fpw+J4rpvIz1W7TYCUUTIezQ0QqvYMW/ffKHAmY7KnmLuvkFPEJeLllk5bWOfe1xUfAMlZhmzTEh6h8tIloZd7s7x6+Pjh9/fz0zY8v99dClPsASbdAC4y8Uivh+2AIW8sowPBUl0ZhfsMrlUqoMmS9uiKWeRLSIe6lKMzf6KJNByF1FlAjsva8ax4kaIZApkqGmRWhR6NkehCxbXV8lgN2QXd5p7Xixw47nwg/oOWgvpaOvUe5vBi6qcjzdky3D6Io1LvqSMB20ZYkFgE3Z5bw+S//+m+bW3yaJtVyKenGzmcrgrn1hONZXlSFlUW22ZE/qDRtAtfc4rzIl/w5S54LLFGY9PA7V/PhthLdMiZ2QiuRkX6Z38xH7P01r/y2Jvk7KqKimLkd+iRdzKqGmY8Y/juqKrqBwa42y3xS1J7iPD4ZqkLuWuPIcEyTa3GVSb5JNjblmj7mJA+0J4Fz3dltmDwTiiiml/BOAfOlAv4rRGHD1ZCQDadOgMaK8UWxbG59+G5VlAPSMm0h3eOWjJxad8LpLOCGT4sUhgRXv5wEkcmLfGqnkuNB+aFnwdWPKNmAxj4ePPvxze7RCzqmQooe7r764eX+6YuDo/03v999OTgJCBJUPj3a/0cmR/gkfrd7uP/yxcH+6eGz092j05d04I9Ax54eff/izbPTl6/3dl+ePh9vBca9YrxV1t7/L4ccbjm1tWGrKj5rzpPRB0qVUPvO9AvtTKurco7YZtgqbn0w4sHdPbOtIJ+jKZGdHUfbUKOatw01tl3C1R5zDB2Aj/TKR/b0s8AgCnFgAXtoGC/pDPHVWhYzNuEDiXkWrReIFoZUwUo4kziPr91EOIehDKKXz6Be14KI7XXg4awQGsAsDaPxAqw4QyI38Sx5I+M+yydvLwjJwJMJ8VhcHkZi5UvwhwOAMn/URKDIF1xcHF3miAFzyHiVsMO7E2OKJmPgEO73FvlLCDwzOtaYbHZ5LypzuWImQ3AIiFRFpQR6GNRcFfm8dm61gCD1JbdirrfAazNZ+ZhazWeSySi9g059mqQB/BIq7LtbgSOKnY4PNZLMiSTSiwtrIXhqvDeZxDUHnGSzpkm2D7EWxCZeGOGCO12DErqmuuOAyP3p0lYOpOtqjqjuqO8oSccistvXlh5ARmlXEF/Lgl8hFUFXuuLbU0zfSk9QsJagXxbIdDobZXh5S4B9u6KFbbO6Q1RaFxzxRzeEgZmQyViG3D2LFv3LEfJr0duL6C1zB+MhS2sMO0MTdVawaHQQJzg+XmpT8xNJolbM1XqEy4PRSL5NXDd4CUMMactIAFP+0ltD/54+ORKFnecxle/UF0+5EDNv+mbLA3WmT/hwCaEi6pj08eljiPITf7EQxUh5ajoJN5UEGyQuw/kiK6JiAbHhpsKM7g1YiRBvO8Aqhkgh5KJlgCh1EjMaXNZzxB+OfkBs3jqwwwr1j4LNiFK7BpB6socK3YtzXAZfMUzROB6DEPa5PbpzJwPpgFCW9YyIBBEimOje72bPR9t734Qrzk5m3JGHBnkTkHO2o+ewGlB40LvbWz3SJK1QLgXH2qlpYONT4GyZzMoEMx0Gj5p/Xi0j8q8tsfQcnS3gg2gFJeC4W+ohwdH5JG/StphNrFyULRSOcBhzLESs/YN0nvbLSfuSofiE2Z+kAx27y4ANla/oVoKPeIMvCu8j4AiQk5+fC2LbScYjz0Oi5LjDy9wNQk0Q217TbScm9cqhX1zHjx7Jd4Bex+HvdHkzGxjK92jBIjbLj4nwFPJryewZPMrGAu8R/dHks6143uyquM31Mv7+YePj8vEgtc9a+xxxmZbJnTP8QBZuiycq6HvH0Rx71HOK4DpG/QQkgHMTCl2BWIYcBBYJWJzyqji2T3Ck4YBZ3xq+V1FQJzI7rZb7JCHVXWAAHbvYHnbceVxXkddNOABQZoLyayZEma4xZ7213liRkbAOLWafWiPrLxPYMvcWCDyqRT0j7GBHevhxKguHjNfrI7SjhnsxPYQpIRy7fx8Txnq/oXlvxTu1yyF+azNeEBcGxyJu8e45FJvd3VHk46AGFWY+6MWyI2FNDcubEHfEBprMHJvLoNpqIsDTbTuaAMRUlvlOxFgtmbsNgrvRpF8saiLMYJ7oHDOLU63m4lVtUkOgDiVGS0Y4zJ9AkjkM7uvERJZtIfHIbLhUysNiYf6VsfmORgZ32RtbhzQw8ZEEqwZJAdugWlKQAo1S9JHIvbhwGDI+QAZCn0RRpB38R/s2GrWvsaXUzH25mqiBCWfJYSEe/HSE6vx2pFtdLjTAZdk4mo2rEdc+rR3Kj2xnZ7cHsm7RcOLRKMfGuUrRnoxOcV/bfMVaeepMTPTNCLYY7NoHwKNEAVhAN4intnkZbt0rUzNDzyp7is6ZaYxFmnV4NNm6ZaCeF1yhsxdYjk5Y8I8e+CNBTPDllex8AhzM918b8YhJXWX5YWBJ+/j6xmvAOwb9j/DhYVBwEWBYKlIX7zuk34KY3swi4xyIGXbTIgMzCCBhuZNUi2LuNWOisM+XS9WHQMTLmnqWOIeAJqIbs+QnkYfLe8SZAKBgWn0kNuE8KMegSmDw0JLb2msp5o0u7MD+9VSPQvmQ8BH3rPiuKvEciHFbd3u7lWBCRi9J+MAMbdhpoJKyCj0Z94KTx2XiUPLSQHTvjWbz7yPnE1NHGYwPtE0EuAVjq/BJNQGdsGRiaQlBsSQPvR2pNYhXUF0VDWJSC88X2TxeQoC5o/MV25qj7+ldCvxYF8sG6DYWAEgCL3a5g3h4MOoqmojxe3m69/rlj68OXFj+PRaKwXnTnVp45Nn+c2a/dVY+zSUbul62rrmm9cMcYvPWy44zy1k1vW3P8zt612szbfS57lswSHA9aCpcnJFsoFAnjrRlFDUS+f0O//HmKM3ShSNglm5HDU3oPdtnwRqPqOaLi1mRDti1VcLRskU+akgbLpPaa5o6w6mWoYk3t6LSMmLsWNoGl+1qsXV5q67aepfzebH8/ujVSwzgG9gii23ZzsNp1TxMOM3fzkMveNzC+HAujt/SEWdoO2Sy9eG33zxG9W/ZULye5ItCvgp+x8oNpUs3312nukBBIbI4Yr3TaXzozD07q/zmH6tyng6S//h3tgkLeQy2+1pjJRtKmYRJYRJDQrwysRFd+sKLcCRAG2Aiyqo9hGqg6yucKyzfgmAIdVG8k6RjAS13t75pCA91eN2VdTXLm1aC9XsN68AIc+De+w1vUbBreMtvI9tazHpDJsxYt0Vv7ciWcKY/5L0Oktgm2kdn/0KBGke9CDA2L74woCHmbyH7ToW7Eb3NVo/j+imw5R32sCWCQ1JbCQGOAoxOhXETlBLN90MSacck69rxsOHHFj8hgbYcrJ3kC8Zq/BSEQPkCRdisz3ykIlGLmJayepBBAn07ys8QzE1I49AWt2VZiWlF1uXtTXygQ/ZyMA+Bmfy53ycCxXoOHl4P9Gv3IPn+Ek0U2Ro56gWYpdcSNIAYnoUD/b0tqsRCOADzCo2peTYIqxTy8jFRQaVGQul3jzRtonY0NGtJR4cizpcWZ+NEZpm8786vN7N1BrZ3GcXyN08kiUmyAtCAVPpbGrFpLp5qtrqaMxGkShIQtNNSVHhL8doVha1Yyzp7HSfufPb8OQ6PAVikM9vfwxv4OP+wBzn44Lf083dvdv+JXQeW+a3Yc6PA7pv93bA6UolzyrTdo/3oPdyO8eH5m91X+4d8XCVu2yjZ23/5kt8IbKBeX/zjvryRkIEDExI/KyZEhs1qPw2oZJFta+QD/kzZ7lblFLJAW5Kr5snI45HPRhgSrO3ZhH1zlAQTS55EEuFn+3svXu3yiH5Jjml5Tjjv3TEtHf2ilo5pnegXNXKMZdLPWBv9ieXQorIAJyxMPOap62+Z9YmIFX9jmb+BL2+WZUPH0UXmlkkNaqGZ0oc/FPlbP7WHQ7BKfnKEIJfIOMio9+FvH+KweB3+NreHsFoSTBVWP9HkD2i4z05/PHhxFKn1aRVGODtyMOQgyLbb/p3Esn2mCMJTS5BgKzmWN2KGI0SptyiTjZVQkcIPmi03UZ3VsjjJ+A8xEbOSOIG0fS/UboSqFsAXNQ4LZyUTz2zmHODLzg7wokiAtil0fJ/AbthpvbkHVcsjs8slUXJsuuKKbGlNFjhp/956KtD/zsHEgNsVuUSJUM+V2H1UhiY5xYS6l6YyMTbXsJdmyEET+S0Co+vtGIkdYKbdc6omsWcHqbRiM01bKF5yFrSh1CXb5qeHmnfRJZ48fDZ0phISg5M54qoCy0wY9kzBy1byEAOnI2Z2hMKWET9SX7uwWKeKvPFylDwUfvChN1t0KggWBgi3qO1sJw/z+e1DjQitFglNxrZvMiHkOhNjfR2BKIKwTGqtrxHD2LZfrgEflm13sZKrUsQ1JraSVTY0VYdL54KSplL3JzzE1kTCJJ7+9PrNs0Nn10KzoMOvWNZxmfKXfXDYVnbA2c1V0irvLWko/+U3cfpQ/c1ffqgWq1lu0aSCJ2mJpZfSBf/ityZYNj6Wf+j/v4j4gW7Kkgio9ycjP/4wEaJ/4Ba/a+VFREZT9yqRVIko9zsmz+ir/NAxLogUJDQtw7QHmTcrPRB75CyvSx7nHkc+k88gCyRmFldwQ32T35zSJPPrvJyxbh4t5zcqfdWmJbko2DheVX6UTxwXUdLo0hd+sv6WdPnKPxVTnmYtfdu7RN5FJSVwitMTyxTtIRyx0x6NkuhBsAeadMYXHCwlePBFfhCEcKowSP27fptEz8BjQcevX5y+eIZiryX2aUJP3NjuZFIsGhmR+63rw6IpRCU6VWhttpomIZWC+++oYC0rYAcRJiGz4rxBrFQ9ADQEArNPvvgyHMZlXl+GK7RLTZwqrYAR0aMjI+QMOiY9KPZdYDLUqXFQ5MtTs2PHg7dRDTvWuZ/SOOtwMVy02HgPsBhuTU6Ld7JWqPhSp21G6lJOdsU38Ka4AJNEi3hKpU5NMyIHwz75Gm6c8ZFzETGi14m99t3F1cLhxjWDL74yQ7A6AmExMDF8xRRPdN71ntenRXiu/fu7mkGYrCLHqdLSp0oK6RGcS/6loCk+9Lj1z/orHT5bWwXw+nSBtHbTsIY4tMn76FbpH2DBEPwLwu2sazRNxulbuCMF9Cli0AeEPjcgI9SMl4lHxnrAY4LBWEqpK+VxQLyhS8MIbojPxwPNijYwki/U4Rq4NpsXRrIRREvsVU9lnh9oODiJjAfRaQiIkNNlfmPg2mmtD5/dMaCwcsmw2oD2/Q2Eg2oPon+EQZbtuCGl+btXSuN7xUjgEEIx2GRxIzmCsvAxYtZCzjeIIc2dK8eaqaNwv8ATnE6fy7bhYaDk+qClwLWd0xTXu2b0afXU6T+wbOpvJlwu33m0UMg6edpUpzocfoZOcCHPzFK02w1cFCw/4S8Btd2hRKyPJ190unnyxW9xkf9n9PbFuNPbF+P/ab191Z3bV3/7ucGh6xTt85Vjq2h5WNP6Jk2YRvJ1NO++HqztN4Czm+PTr8d8VgBAJbIhd7YUqk84yr+iN4+rf9xjDM1/hF3tNot4a0yheMD30/ev+EbL3zWTvymnyHgLD8jZ+f3zR2One8W8rgRf84tksoJ927mBfXFpkOvHzaqNAjJUtgZ4yhn6Tr8DMWPNIawHvUzO9OWakd/Vy9rxf08FT3/CnNEbnjZu9GlNNzzB//HvyZP7tikAtUKkhWBaQE8HzB8evPF0OpJrqJW+rGNnLNgbjCVGASEl8cQNSf8IG9qPsxmh9uIDMYGKpmFXD1KkCEMVnE4lLFrYtkVsnplaDQJ2j/YsKum6/SlIJROU2Le+fnI5RuyrP1ZrghK/CiLLCm3jmEMtFZ2UTgxjoqb0hdug1uJvatUo3ibKOdPWUnm8TvBZ8Rst6sadsM68+kPKDjS9ITsGIvhShLSjuKsvNYECvQiWdE002sHz4gaFW4nO2jTB2rDOI3VjS/jrfevtJtkXG3bkQ8a2XE3crW1txJe6Ef5AR+dTBNvBGvygTkeWxpZJCHNECt7dQbq7Jizyp70IL9g+W/uGvcgb1wc39XsXByeglNxw9FXI7sdEGL/zVtwDO80FrBLyprLD7J+5xJF5Dp3CIcvKubcJDucgOo406PyCGI6K2Itbhif6ir1PW+TZtFwW7mQ8cw8yOA4To4y9f4jXlVnncCGUlwaCFLOd00W7rH7wTmVt/vuNuiW1jtvzjzplv6MLZjedfgezP2AHrBmo5Z4ra7l04zvlVqwlX1GWaNknaFGTzsmsRD67AD4CMbMFgYwnNv0MyrloOQxk2R4LYKbogltO4RdcHWaWA0mIY5nFJixa659XBUyRTg3EnUIlGHWgJkSEJ4Iu3IpUZzVtCDWwrp6gymh8DgT0VG4BhE71LmKEr1NUU18IcHN8qwHNUey5xzqa+yDgnWTEX9fPyPHiB5XQSiIhYD/5wKw5skyq13O1F7PqjK/H4Hf8q0WX3Idh+9jP7+5CN+4gBBX0AHSKuwG/ZLxEVLr55DEq9OsEekqlA33w/TMdrDTjfWO6rblv97U4cuxmd2zRWovIdtmVTYggIAIQHzF010xr3Hc141k5hybKORTRtYInTYseEVcffPiN7Q+bNGolarGVeuXX3C8hD52fbJsyTMIva0+w38JealPnJORR1NVB9cF9RI3QbVPAq03QzdUXH0ZY+X2rq9Vy0t01AN3k7DY6BKJxQlQJnuZWi4hG7hBmWKIBv6KNK+fhFtJ9CHBd3xl7Eo31tWcU+KdDwCwH8E7e66QBHs4dWdBCJRz1uFwiLINOQSDLah4yPfQYsh4f1bOeC5tDW9tEo7kKyH6H1yfyqZh+bJfBDBy6dq0W81/VnsfZqzaf4pq216594cRFfhtDgzuG3KWKLkV9KHD8o5tt6Z7o5sC2bj65DY5Toq/ilQl9b+5rnIFVxFPxDcr/ypN6WV5cFu2z2gJ58NvSiDNBVCBmrQMnBs2iPoO+3BxBeGG704jWJ5rOr10mT2guuIMOY6g4WyVNvZD9Tpgb44poStrx7o97rseIUIyEaB1s4oRocUt1Ry5BbzuyiaBCm25G8T5RRVBlAt6p6FbT9/39TIwH61ldVzm/+BBxB5ojgu+yEtUoVdankH2dlRel6WrttymuTDlYxxrqWHfqbuWNciwVrr/ou41lsVct1bcoBKBB+QBqDRXOYkK6q3v54AYCIBW10kejepLOq3x83I4PIrEihWivepTJgrm3hB0wTx0FOIoK7wcK0P0+lWfg8T9K4qe2RnZPYm85OYc+y/7ax6LLQu7FMbt6WMm/3+utSK/vqGQ93nRY21aHnsV1U0GH3XpBfxFbbD2dtzhe60f5XW23XQqthhwx/wlc283K5jV7WnEUwX1qtyzqVF8V+hia5agT9vF5flXCLgsGnCeSd8BemTM023ZKaTYDGqntmtiqtezTpAktp/3DFexinv6SdKpv8b/OrUKb2vIBTPkLqEU0PjwZDjWngYYu2H314uU/sY2ispxbag8+Mpcj/4INrZDiw9CE/HjlKroPYosVvqmLOYdhoQb32s2E3zoVLa6LbnWncvt7qwFvFrkbGxKyBZ7Yt5o9VsY5b4BlA2PYKcJP4XK/ndOphI1Y3bLnUne0eoH4Ts4BlfXmuMtF3Q4zElh9pTAXGkljcbARdMcB4/iwHctGsZPICadP7HnNDujSFDvacJEsn98GH8R7hBt3cbz4KXQBg84OLi6PT2sLO2+1XXTe4CRuyUMQnIOqjdiZuBuvw46t9PGUDYe3Wod3MIhj/R51bJeD3RO7cctsz7RXLVYNGpzPW29GJnYNnRCN1xhtAt2+IJpo7/IjJorFORj60MzOi++GW7nJxMQUtuXStsAA+6pL6gZ0tPvdy/3TAzV1Di/kT1BGbXiS1qRD4R3lML36OID5KVDX8nYrxGKJGUu3L67RBhux5SS3YxYsW6EJ3oZ7r3d9MBKPHBRztnzuG5+t1VzAGO/wj+4hOYsN+2qGHxCcclMcjON8xnfdrQCNPIer2FYo9JzSy5Kn1wNsDv2rLTGuZWi2BvpEpR2uj6ezBjAdWBwqDZRZt7pbB6966nU6joN9xqY6fCPEO5YtVSSzACyg1T+yWLJrZuz0mKSl867IA9CoicxhF9syTPUujHxICbIwDECUEQMO8DkzL6I9dtAN/DUHzulxG1hJPNWYx+BVRXvMl61rsDvzdQFLo27VkMHFD42aCGOVsU8pe43ySD+8u7t7GHAu5ChWq4Iyg2MSnJX2uZwlD52LVAT2HrI/PDsdj5wlOkcoj/bIhSISkC4H09sLmzVwbANsyi5v3dUvbuwTZnUEiJ5JD+RJLZlDv5SgK/VYy1a2pPXOCDnSvjsY9/Hz759VV6Zx+reZTnv8V52RB3cvGnJkttW2l4oscswuRc1a1ASitZDUjg7BYZd7h+BXTUbRY4DXp3Joz9noxthA3TYsNk1vs2qBvXOHQYtNnQOOLGa++s2dA9MIZwOxVk8XcjMfnKK5z7zAMylCvGr+Vw1GptAbTj4A8ltwtOCwgt6fg4tByHBVzgGdJELow7/867+Ns/F4/FA8SeDbYy5C3v2nTbFqq+l16DCixaPwUi6IqY93i6g8QlwiIJQmb7rOmup5+Y6AurWnlOnjf944HmcnnxrdKbGtrBrH3ZKoYZvqU6jd8BejNgcbWFWaqXi16kCBPX4b+M76GLUcDJCdksURNQiZzSvM0UbY74P9itRhTgZuqwhXkGleX4okcF5p4Ftqy9kjt6JTU58WGW9BJ2wk3mdCHLL10VO37pvjcfJpglLJY/k4guMYtpedb22NNUylYKafLNJTHONR4GIccjgO46Qlhhpy2JZYXxOVO596KhYZvSRglbuyQx/RqT4en7RilHPMDYt7KfFVIF68yW9rdffFikbYTf2I1MhMnZ8KbmVanrMKrfHEYZ2whGi1IDjC92SkXj0hLtXIbmhPUMSI/ea4bxeVxTI7eB8xHfBJluwp+u6yHxHLoV2nzAoMJSdJ1YSLII5GW4Tz1VnroUa4QZzKLKb0GP2Db1KiriOgDn3Nzr0zUzsA06AWe7ay1pgYRtS0yqmTlHhNSkoq9i1k7yVP7gjZ6JlmHke06BZZTvjnwCVM6S4pxeR/khqsZJvTWyczRfJOyRahp+Uyn1pilptKFw0x+7Vnau1nAn8anIvdNYiFiEWRQx4BrjOLUbAMZxWMMoH0/bHzR04GioQyc0TWRaxgts+bvI1JsB9fvdp980+nb17/BGbu+HjgWaqPNbOOCJm1hIlpZgbO+4+OV2CA3ZHyK8qNRrbeJrLf/LF/wEGDZj9LQw2Snk5WS9bef5hBbc8K3KdqyO/sb73uobsiv8ayNKj+QR4Q/d19pOPCuqVqK3Rjg+xeXVn/qesncF/eq6haP7S/XiW7ZqS/Th+7fqC/ThXrLQ3uuMbRYvI9vl/113dO22Q/N/WBW93TTOAJ8eSLxJnLw8dncnn/rqwfX+D08MX44xvubzTwbfjqbzLae9wYNjFydiTwnXUG99HuCmta7Y6u37ug15mgA48+yu7/bmeCuOV7rfu9Mf/JidGDcK2j789Q5hWjYmVzHC1OP8ThWrraEid5yNP4M2NnF7nRXApHQNCR/EWYHXqrSd8CL3uN/DhhN+BAjw/S3awSs+SnCiIodrZkD29/K5kgcNI1YcxAy5eTS0fuSVrwqXnXb3hCGCVCYmE7aBTiN6Z3ZAVApVj48lQZjFojHTjZXRR6cqgRLix4B/pOVLbnV06E98QGLVtSYG0yyLEThGLTrefUSZ4SzHVIc43fKAswL+KZK4XKkQy64RAklkInwkERbrkLigCdjRGomt0DMVEvymulVVk67+aO4AwjpQycyGAr2XUH6qZawT6poNvTWIzJlnpXE0Bp0DQ3MT8OiwkBbUEtWdKqioMDYGZOheFSyDPFeZlf6wy3k7lQuYE8ViNgbQeBHDSE6G/jpQ1PDsc7kLMj4Ytyo1/ToBnQ2GGACKTslZR+xEgS/1l7xYmTB8jmHwpR34qSvwwjl9IBlPiFy0wu2HdExkwte5ULkUcFIqYV0eNqjYpvXyAzDpVU7PONYNuSbZBL8wN4USfjQVdpEDz3xNLqhEH3XemwfVFnCPHuzgdBM+3wAfqn34FthciU+gfkS4ibSsdioL9au9zJiQ+f2HAoahmmRHg91mAGCfI8B42iJKaOJlkkhDCn4fh3p4gv4gMhB44wPU2c57M6biMKIMnHeSd58IC3fk2M495xDT+4NA9haLFMdmwVTAFHPEQp8emTb5LPsOUYlIUQjjSKgRJRLkzuJAQsh0BSPgDqdWJ/vRpOrMNC8T6logW6C4JVvS049w5Kap7qiEZMnvaRjaJOOEEGMFZHdnPsQjPGV06Q0TF1o8mdHvAn",
        "ZMjEjwwRIouljwyqg137XaTAw/XJd0UYlMSVJeWoSHGD2ItOk4zSHKC63pKaS5EyyVbSAy8fN/1NMsbo6N0xnuUspKntCZ0e/cliKYk26cVSQ5/222KVWQZUF874fQRikDKWAYCsv0iMyzgwBfhlrNtpztvjIpK8Un8FYRyJvmHJ9KmAkROOxD4tOCh7N1hicHZGSZZhMSTvckgq+I1XwR5LVHbWmhcgBQbiVuzoIW19rOlgYTxiCuHQw9MEwkeOU2WiWlfIvXggwfCocOvDlouDdSxfIGY52Y7GLYKbHR3dU7V9ObHcB6qnj0brTgc+ShYEVM56Ty2XGXLKPGnZhfMJ0xXokDjwHtrysfU0FmlfD3q05VqMOH3O+mJe398BBCBZEFiV6y5Z9mq9I18OzvtEDvtgYBeAX9iyH+QHNL8DibMvX1WgHuTTGeq8LYnMMWaG1eAEznFgWOCh2gNWMV7ihxEH0pxzWtsArFnFd5JoIE4xI3rpd3yEZQxDl4qWWko5eKtEL5Zj4mhBmqyiRCnTTo0wpA4F+b0btqERQnS+s9jIbF8CVPGWlusdepiu6A6kBHrOeCI5UgeMkDfnMZXZwpLScZtyjW85mdOrvLnM6p+XTbqufnpGZBj6GRJF5h+01fQtEiLwSaTGeyChru9bGSt6h3nC1ZYbAQ3ncTCMt9pUANf64vYyDcobOmSQpsSCwpWRu3iCciXU63aiYT1pBJmOq86wdnQYsJD2zQ0to8H1vCRmxb09ia7ZA71nbWRicRc5FcWxZNXcUsAi9lYcYnRL4S6oWksabaMfJaL56cCyAecndRJl1QbBECHuh6XBnbfOqO79yZ3L7MBTbZmwb+ZKykjWlwtkrGC7IFlbTtTozgNnnYivAa5GbaFJcQU8RcpXJPzUphmHktadmDDed5sS9oQoA/lBbLeZJPtJryshkZA1Sub71qXTWybY5FM9LXqR/Es92fG74srfCT9eOXU2JtOISateGSkEIfbcn75B8h//Nx8Nq0Q99NSxKPVpsPVYSMSxDFSG/thxQCAOrowfcZeQZkjMaWSCG4hgg5/1Igx1TWUCwaomW37hoc9bt/Y9S2T1NJUHDcd6oZ/91W3bQIRII2guitOq2E/ZVb6KYu1n1JCcFX3HdFH0xlFIdF0MmRPMkdsQIBS1eLV6CmyixaZ9xDpawg/m1E1DKBF6o55Bd1n9PSuMTeCyWVOxBGwvr2G29P5kGKuSf+FFwW4IAt+yNZBBawLfTPPD67d2IdjnENrnsm0jXggnnGMFF2NLqOHQqGCmDN+PhP50lpC6H/W1aY2FUGDa+Eb79dTyjby4jnS/2MbrgGzg9IR6aq69defj44ejP8xPzPST7a4e8iH2Bp0P2aDz4UNJ8PAQDdWIbmwQmMCO/vzDXMrQ31AP7GLLeh0scq2eAUDvHf5elyLQEmMTeUU4QXEKcdiSXVmWZtMGhk3g6Yj1cxUHvZKOTPbRL8p4Q52lnugCY8YZEFTKtcFtWBhotaFFGcBblyDLmavF0pEMg3Zry4Pg7qgBDQYcmie70doaBNkuvM1xNxRlbFxgnTAjSv2Eh0DZw3Xa/4Bp1ERWdg2QPkp25fgEWawsz1BA7LqL7hi54XZoav9q/+jNi73T7/df/uAM1HT94rsOuR6RqSNTh4wkGq8Gfhqp0xfH+pyKvpZRcb+I8Ac6h03CueL0pKCsX1G9uh3ZK+GHyuX/O3yWDVqWWQ7EeRllR4EOdx6EAxVbV9Z306IF+UpNWspFcqQjLs+bmiMSZwl0ECIkDtRAMANxHmra0bZ60Jagj7hsHOxITECxcgj/dl4hQJg4kGbJ8yqW+mBlitn5SOSWRbLGOtiWo21M7Bal30hBAnjXoklwis0/ruom9Krz4hXRMJbdZTCzS80xhRTQzoNPxtbK4n1pozzCEW+JLe+78u7LxmUxW3RS6wQHm6+O5WKiavlstqFV9aBtzHLifKIcKZLvkCU1TsDXA2F4UJaYxWtf2hAHCIQuXgE53ZoqvReeqQY1jorsv5b5ebMnSox2PkH+FmcTPOQU5hAcI1PdYaZpq7yUZTKjzcLZhqGvNoHH5RDF/SOnqWusoNgx1bfzyTOUoHlOkNw65fL7y2W1JLj/2RdjReerBQGhwmKWqL21C8sP+W+ROkOweKrcWFpwk5xKQGfAbzIX3j1sh1vJMbjgxPmh2llbP/HgLPJrkeQgXK02IXYKxAsvMZT/cvj6IKsZdZfnt7IGHHZaFl6F8PJgu8KF/mFVrCQ/g3/UhVRDMULiMI5KZTZpxOZ76QDayG9y6EgWZTrgxgbsW4dfwUjeR1yfyii4BRnBGy0ozP3SPekAdyfg493Ura2cg6bK6nNjrX11+BbPStbGdF64AJ1j0G1el+DxYwRKr5CDvbhYzQ1m5bppwB7J1UrSLtJ5yetCYS/bk5pbW8JpSxbwbVgsCXnMcySVddfqMJMVBtt314J9YyIfrdD+TkAoqq6xjnoq2hdvuMnckZ0c23BZTR7jd6tzq8WpKbQLBr+/1/eSUpDzaxaNOJmvL+a2DIRxZ0MIDf8gaxh6Bjx4EFEx6YMH4drZ5KjzBoaB2gCvS+v+xFdNkm34A7jmegMTu6steov4UtQEKUrOBMNFCBngfGwsygmO8gahYc6fkk3LGkB4KhgoLBN/xOxYO4LVkAKSd9tevqquy0LycXErah4UN9IqHkKexbI4r31yMnlkjdl2j+JCGRt6qtOX+S3ix3OmraFI97i2EpLhZ34VfnW53pAdXSeuavlZbuSs2CLCKHdQs7OlxFiHkjamjWJr0fPyAqtUTMsGaYFFLy++FqII56bMPHrJeUdyZjG4gdy5x4ivP3sLrmq9q3cvCcRTw5CpWDSszsS8B0g493kGpsqrZ7gEnVv+yydQWDM0JCs1TJ6XFxnhEJiOsweZUBtUNPjtCIfpooQWSe86moW4A39NcMMgk6q+XrBArCdzE2unOK3ZYlax6oNI/YaYS6Qp89usE5LllnmJdkreEA1J6286h7ureAooKqfRsSU72qE9sskJzVnnuK4GneFNvhb6To1TfGKf/oq2TJ3BEBeMgIKyzJuDaEh71axacY7jgRSLOp7I16BrV+p9tx+xh6nbPcnrPU3VxwATfchrl4LPMplaiTVzgcvTJXFIWCeCKhl+c62+NQ0LI838ms1i1lS2qiEe9RfL8MveEId0FOBlQHD5RVNcpVGd1umwhrA6TKUk6Rk18lbyz71fNycm0oM8ekwR9ZSwZd32OX44kV034RxuitD+SPUoMFa+9XMEIm/+HqIcwvSzfJ57kTO7d0VvbXM+Qa4os6ZtN6lQOWD62sVd7q1WQW7bzbBVTCTnGogN+T+Fi/6Pfyf4UyySzzntG31ot8pl428EedQ/YzxKvpRFsgCL7dkM9vL/+H/+8t/+3+TFVX6B5tKfiMEgTHE2DJuUTCtKKmw7VFs2ZT57ZoQyoq0YbhRSuM2JgEeD4uQqf4exOfLHV4nIpjGSO5Yz6aHVckyZ3k1IOMTPBoTwAH4XjmRzZNMUTN0Fv9XCjLkDr2SsqqTalc+1YtxUHg2y7HA2JKEj2F6jQz3wXHtq+QSF7Yo+AYLUfXDsuhhF9U68BjZsf9hqTyiWoc6TtgtZlNMAB/R3bfeaDfVVkYtzxKQ/v8yq+aziLJ7m/wLY5gcISg3ltrUXdvxsT37IkY8IdzHxymyeNV8voe4dEJ33OKhg99jaa9NmgAkcDH7Lb43cEtlPwwlS3w7xgs5LQJYtq1LlvQtItupIvLfQ1fiB6jx5hiFLoYwO1X4OLu9nbx3zNFlkwEpHVfrz8Zj4/5+PN09AIywy7EP8Gk0RYCSiCW0HPNbCp7U7zNiE9P7Bwhns3ZiKvZizr8ztKLmNH99t0uNG8Dl6jujRn9ngVZv/RdqVO1bO03djjH8MFHIbfbgd68S2pS93K99t+hrRh9tNV+O9VxegA7SFalT+vV8NImnNbkiVdYrPRbB9VU6nkjOWE10phSAZo+CxJk5pkFPiyG3bg64wCw7VIwY2luLLAszC+oM8aW4X7FXPBgySXYkTpaLhOzbITofsEbtS7CDZ2Dv5AzHZuHuooLfuZK2XQsfpWwLnm8Pkt1ZLcA2dq8myqtnGDwrRT3ltkw1ee/i70ZPctuTRjhSlncLvlMs/kk2ikvrtVr9tyrfN8JuXQqiVGi+ALOVx6b3ceK/zszpFv0Oo9YqNr2GV8g4a/M+oQf6AdQif2R4rPcvo6D1K6M8mFP5PRkSsZLfy6lZenfihOK1M+SfdDxmObIZYFXhlPxrfpMVBH2JeQE3y8+14GASfQGtiEiN2H/JMQBc/aD7jYVYjqWRouUCtMFrDSTlEsR0ZlBEFT+XxOHr5La0Npv1kHFKWPSjG5llecarFGFbzy3WwukZxLtGF0o5B1tJmCqjFrwBtX4mkLmXpwnuIaG50VArKrXcW7Lnu47pHwlgiry7E/DyCCZsJgy08QwTMnNUMeyYuN/qoqRKX35zbzAbqtGk9B0jEYwNHCZqcYJIvp5G0+APFAL8RNXGYkFTFMHCfVUEFZFn66YVGXU7xjui38djJR1uiD6kUCIIJR9ErTJoHwlhi6qTCif+848Q7iXAG9wkW+dVgaAeg2zK9aAQm0ZLhAZyPI9+VoBh4KUx4TvlesIwIPLMIfYyQZbMJYleULsEue5zXTqR9Q+TJqr7lW9yRP539/7S9zXIcSZImeO+ncKK6OzySjiAAJlmZgURySCZZVdPMJJdkbU0VEk1xIAJEFIEIdHiAADKLIiV7WNm97o7syK7MXPu+l937zH0fIp9k7ftU1UzN3QNkZtdUdxUR7ub2q6ampj+fhjfqgNnVK+aqc/X+0NogKqUqhn2qdB4y9CWvsbEQyK0Jngrpc0x5ByXZGCawevluQHfyULRh1MeSYbP0dRUADNqp6tNQa0mVzUs+HmoGZqxlapxKEnJZ2HxWjAiAogShmJqFXUBG4LOO6NE6HOYRVon73lvORJlHI5GbkUxLCS6WVi49fBX/FhlTBkON349KUmjupbRRkuXZF3t5bXgFBtqcCNHeuuWX8kMavXUYmiMdLef7QmIxJIQm1VQ+nf15BgUWqQlcQSJ9DYYek+TMYPJhTFOrHslu7nHujxnWEJvFTL9H0LRqxODmTrli1kieYDcHGu/A4AOuF0UMrZ/mI9kraQBtkvyxNVGiNgbpS5p0/PMdPPrTrGVqerHSvF4GUoQ63CT/eNk/Xy4OT6dnlrz39m3d8F/vFXeHvTu+nkx0u2uLxwCCP71mVxMb0o5+aCuaXe8yA6GxJt3Uup/VNuC3sBSfSnYBKnTlb3RI/mrnmed0yquXTH8cLgkPbRdrRcyLzGvDU4LqhpUeEzf44jzteG0otpOlej4HJTK7eG+DbxeLCeuHicL1NZ06tzrf+Vta6CA4hULJBjKBynU6oQW7pADyEZX7Xq5y53LllcZFUxFO9slHr+AF63h4eqoWCD1LPsECIC0MzagVLaDuFMicpd7H9BFwtJEfU6lxnV6Aa/dDS8+fNpgcp039Xo/W7hrgaWp3REkmyWjiMP+zjB3+BNcz9W3LgiCdCM8fGt5oHM8RLFm5bdtzTl+iW6O8kgqzoqaxchOTGJH7yhCiJDEn2aG6B0znk5F48iVbTFjDc9/TLn37Yr3rk5XoKg01HCGIfd2F66OBPsbA0XqOMJvDaPg2nDzNT3/9V7pwRSPYg6LMLGUPJBPHtYYgho0JTqLtYTqAtx35B09dcfz1zwcdrdTh4oroG6ZR/JeL6fJakqsulmGnlYNfCf4Y0PIr+bvha/2BgRKR4VeWIGRT7KEwF4baP74bPmLkTwY3JKebv31z8sMn7zHzahhBaRfKnZ03ib6xLzZphvrE6uJesirCV22J/hNYY6vIOvetWvKufCNeGXorSMrSYU8vPnGqzUIolCcnm/3CeWF/r9n58XV3l6DAoGptbvkFT6FajmxseD26DNEHt9k0osMattRmdX06HV2qWSiqbu5a/fYtlbxQYf9DFvYXww7/Xq4eCBxKx1wu+cj7/JTVb3btZb6ju/3H3Tp7hh6JVt0ck9WIzbzsTuWuDOW0PofgjRdBAJBf8RLpRDjcHG7udHFYT9Z1XL7vqXZ2hl0bOiChe/rLqd33vCfKRzpAGWRtD2LleNHpSu4JYOyW9/pITRCXedoNMNRftoR9wpMuq4npaP/OPz/GGVTUR/9yMWsEDlPins2t1n8yjGRgdcB1mRUkSJ+K2H/iN7e4UgSjcMSNhGM7uVAd9376639+1Ssh/vTX/1KU8+l00kimpfl0NRz1EY2bUTeD1L+IvIihIiIRUtFHFpdy5ZpJHTzlFYKh+xerCJ/1VmApmXw2yK0XEa1N8ObkfBtzN3TfJV/l70eGR5u7ftsGG/zisaszC/04/g20LReoZ2Fl0j2JPbmxRifLt7vm5Pl2AVA+3vPf7muuEt7LH4KWdtNOiN5w6lCgJTOBzYhanSL9mLKDAFHnSVaKYh70XeGsWs4Ez1rulkOb/6wmBKbAvbSU8UP650XZPieypzgOn4hDcHhg+Az0MOVFPohY6EsP3xc+h6+oV+DR/m34ZSIHvUzd7zt3pCXZxmH7QnRSeJumMh/Zi2a6VJ2MWggyDfGkxxnw5IdkoZ+MLLIvSTuI9yRH2CsmIxOyRgmWrSKcO4kPpR6ohV7vSEEYC5Ok0fgsUflqmiBgn9VvovIkUM52aDC78mF/SZQ/O8HWwuTh35F575oOC4KVzMkmaKxrK34tkBHq0UDKnzjxjKbIGOoYBo+WwgxBQEg8Iv8EljjF9yjSY5D88dmqvH3yQyDLF4aAUH4+ZIzQulfqRsCoI5py74ilXwGa7jD8mNtVooQiSK1u4pFeZr8OEymEbU++khmjN4FE6EuGEREBrFRoanGxEtHa+4UPjMtF2MOotcvx4EB6UjFQaSqxLtUOX28hGqTV9fl0kqk+Re2eG9RtYn8LMlSdJXwVQtWbhOJHgJIinN/KViWMPa4jQcrDA5RTlV//GlfaT8htUm14FS05ocxmwbrUmvNZEQHPIxq4xNRYNQ8sfxGAeYGPgdnm2sMt3y/+fa7roChjiXK7uMPmess1DFjhTA9HFPKkww+QhUfd2NAs6QLzhWgfmYm95BHG3w9wcGFhx+u798VQD7ihu/sKT8M6MFA5XF/VegLRiIibiKXqvpaXu/IqSfnuoR65j09mp5NANQLfmULW0KxOcTILXdDKE6+SgLVbTZ+cMqMJ7garlSjz4HUH6gM70KfyMDsb354smMn1vJ4t9X2LlfyeEQ76rp5MnrwPL3B4TOeBgMMFL3BpAF+oJWp6VR+tmPWtbF8no/OGI/c+bVNbeYQISk5XKBpu6pytMvRG4rCT/WcR3vY6id7ycsdut2j7Bra2YH5My2PQ1fNQSpxwNKTQyzmjBuneyzv7339/5+DOcHS+OC89eelFUT2e9zLXXXKD2ZHyA72sKkvY1Xc+LCJKYbiXa5WxXI++ZWqwvfqXEp9+kFOC9fABAN/FQzgiNRVP4HfGZKgKMjIzZT6T3tUMjZR0OZMoXrc1TJeLIG0gL8QmJ7rHoSvpbQZiY8lCZuQMyV4Yx45tiJFjk9fDGxso3XgfqzbI5ovXjc5whGm8UJtQuHYGuTzeJK6nq3ygAOnvd4XLhilaIgWrh62oFR80Kn7fODRVHA2S8UcdKWhIFQPVEhk3ZM6DWPRWeCfPsOZ0BkoOn5xTWfHitL4eFX9cXDD8aTnNzbKmSQ53p7kAU9lyepzgP0CMg/lGTgRDTYacaXYaDbSKc0TKqdQspKInuMOuBGAyPvEScVtlugFI7RvfADVrY9gOV4zmlMg+54vLZJWocvil6NUQEfTDuYgPwj9L/bvvWoPHygwYk5E4Q4lXTpzHOEzQS8GMr1uWp565GtsQIwaszI8J6oyoi47d6vEcZ3Y3VomvKhGxsTYgyM8KThPvx6PikeK5Fmp5w15mpPeuQpb5i0EDwFnB4JeYsiveqNur4O9rmZHIasMa2N8P0t/cF45FaKwHbzMxKoV9D4e/fS+2TDxNHnX6jjIA+bQVdmKeAXBzmge+xrUf9QU5BA5Rwv8s3MwXk+sU7dCxRIufrzdwJyt+qkDxt9rROmpF6dj8Yjv67IbgKxLfnWLz3/gfDQCNO0jU1mQLFsFosQJH9ekREOmnjWoputOX2Zvj3Bk7zBwd9OFgaPatfCvbbXLNVtEwUijxjfsQwRDqIDItgcQjeqCwUl+0S/KVNhKKwrmsl8aXamEnO9G/jdR33Q7YKxwtZ1sBwA7CW8aRy1SKZpAV/GA8Jl2qWzx1PUN2A11vFE+MfQ5I4LrJOEzuBmF9AOeqY85m1mBA3LXAx0PawgrBliLR6KkndXFczyDPIRbs7Hx1LdPOCZeBJR7pzoehhWS17ev6HF++knfClF+NgLodf3Y0FkpslCwYMA6h3wMX+IRoJd9iHVOEt+DED3MEgn0ADN2AlnRQ9VXvP2B176bX5o+4z79R7ewg5THDkBXdDxQpoygN+8xeZUNNCPhL04cIBRvCffR6S+/HAKgQQA4DFBuP5Y42kqwVI3veQh2rCF02zvKN2O0fRH/j1xFpRo0TamHT23Qel009nRYIsx0uM799/e0z3Hi69jvCTUWc/xgkccM9TPyz0cYiOilfwk0KD3Jx7zKqd4iq4lIEqH1O1CRBnjyy32P5beZSPxJ/R1oMvavhLVsvgdaxjAQW/fRKF284jH+ic1EBJp2P7STPa93DecxEKyqETibhYrOZuHZ7jmcaNiUlwg1Ak54gjEpc+sIPwlVoPHhlAICb0qn0mEEw+vMAXi9BYEi3JHWKIEsqN2KCnmZjLB6Pob1hFNR7Gmh5FrbSQPTmhkj+BdvdgQPE42bjMHYnFGLhH3+tdSGnsabTqRrjWnH5u/qyl9J/KUHHSPsuUQ9eGUzFWY5IQ9dt6dwgdqqfaG8OFzSeZRKS0Pj+wDKHDTQLmAusoJciWZcL4jf+KUv0y3d37952SeoizgGCfrqJ/j4yF/qyZ+Y7OA3G5oGDsrhYPZqGi+afpstFi/QVDUWc3SWVEHNKmGTThO8R+4kzPt68os6e+dsMHECystmv1/Aesx/fMjVvX6Ck7XcJ39VYSf6wgEn7TQ6gOphWGfcYyvRYD3+4enBQbuLD+AOAIZA28cD2ihWw374MySqWkF/Ze85BKiA/Y4nIhvr3b8+8cNQUrl0QaXoQR9QZQau77d65vrS2rLtl3YKzQj2fnm6incjyhlzTb8KjngVNoVdQ/BGdxUVcWcSSG38WR0Uz8+T4WGyk4agQas5ioVL/fAvCnKcTjcnFUWWNyJVIkWL2CnQ8/WJru1ln4/aKr2Sa17yXrYkDY8QYGvlrO8qO4m1jGM3OnfyVxK/1Y2tAI/APhsejfvH9CfuYWsZw1oNYTZG6TEAlYR8OBdyduVgDMya2O+GOBSpoiZyYs+nKULH/wPQwafcvJPNPBE9fMNnEcskJL45PL5Cw9WgKf9pLUBEQKi+LH8BtoEMg6DQ6F3HpgbRV/IGXEZ88UJuXm0mKQ280VRCa9ekJd1UFNXbqPyQvi5opTWKo6smxQ1vPTZhtTQXnzNCt/MmaqeVaAFCnp+auTjeZWyuV8ztouRG9FGf0jx8q7afoqG+Vq1H8HVowHRSBTAh9W0L0p9TWw9wfUPtjz2GVw09l6529qgiqqxZ6qszxXuH7gRWLvx/kP81iiagMdLKUCgRyNRoy0AbuNyONt/+xoCA1Zl38s0rwU+EiS2trTAOGB9sHlWQAG6dEYKWUw8thmslxmEh2oookwo6hP1LeSEK6+kD+HZlJe1xsMSLNq4TQf785X0dlK/2k65SUUuHZDHDVMpHi7ThtAEEoQ1YG+PdrEjVic8naCb2nrdRRaqZDNqGsiNMm/Zkbazal9NGeLhOwW5XvlXFXQ4Art6gjLUfLmUgDzSJ698mxIjVlgrFR7+7P2kQI/MVfm4eLq6y6lW2wPt9xyUcX9odgDUZ8xzWkQkaIv8Pi40VxR351wSM1M3cRYfPb0JFzq27eBxsZVqWZXRlcmuJsGJLfMdZ+JrRxQxYLTfAxm/flN1jHSalskSVEYoWzyCllpjwnjUGQTIcoLNKpzFwiDkL5+C8FFa0pNp5NEQ8Byv3pr//q6HbDUgNCwiYlIYxl7HFGkC9wFdNfzBoPG0JNlIw5U7x8EstcHCeC0LW+kSySUfpIYiUHqEJ+Mz5yMBxE0ok08wsQST3hGruNWqBevu45ruPv2eOMz0vlcq0JtUs2Ul6JvFiR8jJpcs8s66dgZaS1HDBbmPGn6dkgQhe6s2F/4AmE1ZCglIwEhd8ASzKgNe56wmFJtx1QqFdvCZ5xggotfyzco8pDhYZdOKzcATTWw66rdrMjiHSBNzabtv/9M0TQ9j3e7n+8g/hcPNB1Cw9InKV7NjyIBv+eOS1b5+9++Jh0wzFkjSZqkPmWKNtPoMnqUwiXce77wNT5MGyfhnGDqzdWtD6abN5BNSx7JVqFmVQ917Bq5w1S5V/feQh10U2pVfQAWfYdIL9Qw0lpgglE9tZk++zo1vJUnwoV5Y0FcnDvJT3pI0tMyg7lz/atXgGXsddSrSZILXogTd/zhkIHJPlCcaVM4HEg+DDHyVlH/yCiNhAb39KognrMF7tV66ibIDiTB7jvwkfhBlGfT8seza5Mb9KP5qpi2uikxBqdaWmUNPCWPmW+TFMgfbhtTkGcstIj5vujQN8mRilvo6iqqAXx9H/dQ+Njl47rNGyb3SSlnbmr1pHLTKMmX16CmHlGb24rg0LlFc4lWELbdMiwY7aVk0qkwe7dUbFzwy5Moug0Tx2W3fFGme4Rde4VfkFutRbE367s6iC17en3D6ibC3+NDu0EtM03ht+aWnEy4GYRauL5LZxQ18+vHl6kpbotnNXkfgktAI6/llTpGN9ngnKKMwgSQYqlQG83MSc9SmK+3M2LZXBTiKBpH/x9Vff7eibxC32X+bwtySIk9l/A3bGMYjlQvCsdWzzhc/MWnMLhhtrxE3dKIu/JTsyecJI7AWOYK5u/Ovw6cIyv7hx+/VVzXs+/Rm/d9lvjxK49eWDLq8xCHn/UZ13XKa3ZV3ek7Y4Du9BfYLJCfjX8NPJgexxwmq3DuQUqtXKfxq+5QdNPvnW0l5dLL4bSRl92KmFPIH/f4tACXOo1GakkqPqEGZNbFbDtYRT1jOGkYOq8julUUzzpfZUNMpB6tkqtDPKcFdKSsKc4aFE4rRm0L9Liza7vfJ/6Tt6E+Ayks4N/kO54d2RJSnBjLvweOSVlAv8BLO8Ujk+rhUnIo+K7xayBE40OGFn2APBCZczlyVTA0k+m10l3IcbolZsEWN7BUv0Nzd9/3I3H38qq6OUjF9CMA+Gjw6nEC66KsteiPkyO4h/flZiqSPWJna7dqXn5Xl6bXnuO2z0qscN9VVHMHbT0fXmlH2fOZc9VKlCgsLsHKSrM2KXnhFFt2cM3oqh+W8lWMjOAJWsSBcnPMODMFbSA7m1wDVB+Q8rsbThWxieSFuBr/5iOxa5pSxBgYoxjYJ+8zJbaQyW2CnmqJJVcZ8EhpniG+Xc5Pb2FBzDSVQBKZ9msaLBqW7ZNmkqFd2MVln30EyqICQfd5+lg+oQKfu8Kd3Sh+7MJbzxMt7S/3+pd5fIdMvlfu3UmOYu/DqJxpettq/do6ZPI5BCQtpIzWbzS5JaDXMt3fFq/hcgCpBqJRgFMieAZqYoeRUR/otetiB9C6IpZQ6yqtwvbPQRY1MsDLmA0s0dbn96txFjxP4iHSK7Z03Q6WBzJLxIKS/KQPnNzdEfKuKLq/hc+a+rCkoRWMUuoavF7zQJV5hw5nTjWXWsfxeV0toT6rKNVMuwS7/6p8jtmaG5QCLdScWoxUzah+8P008A+f/Py+e9fvErQX6qx0BsZk47Jt7R8ORAWA034ZfbomZ1BbZs02hfJyg1C+7CNK5b8uTPssGHHD0UYyxAxDVFLVKWzOVP7lViv58fUlIb7xXw4dMiZQYBWMukYvZW85G3yNeHvXTG0TWZLOjEJWGuG2wQ/jt+KxTS3qz0jEHK61ViKHBF+NJcOiAg1aNJZoTNV7SvGlORHnp7MjN6m9Rw5gSWpsQuVE42GmrcIrlSHAovFOzOL4WshmFI8z2DLOIZPinrPacZZPeMxTIbBYZtS9q0y0CBTuIiWtV7S2PR2oR2ALy4ugaFsTANxpSn7fF9Gxe8ikLN5BoSRYFLmItC4orkyxq1DpoyxtIWeZMxTAd+oF3DyfP9Eonpl6dYEiNQiqDGTvGRKXlRirCpNLRBQY6TCJ2FR35XWjeQuIgBvQ/FIxhTJndViJJN+C+MWJ1eDjTi81qwgY/Pc99jYOl7+Yo6e8KUtXEuzhn4iBlJBvCODk5V1CjLVKMgQWmviR5wtCl8Yjjpq2JwLAPquveqEKHXPT/ITJNvUCfzlLjThfFmuEJlEHuk60e8TQ0CFsEkWFy3IdkA7y1XZPYzrq/znq2LbcwYegdxA3Bvmf2o5HK81LcESl9Ew9kMAiOGYHcfZty1POPK0xHrS8qDu01wKRUseJj0P7AiB/iS+SOeKMsF02MYJ1DctF9aMNX5IIybFNUagtDcxPfrYcRrH8/RHYL2VWXXAYCNJqotSa4jcW1zOXIpJO0dAverL+lqoWfvjp1EOIbA8zPAY/OsUWo8m9jTCe+oWU3aZBBmZV84LIjZAs7ccp2ihzsvjzkx2It1vtTgSMeMdU+JvWT3H5JQQc94m4YT8YDfnURlsN3vubaca8NBmXpzynG62emlmy5AsVTRo0UuO35eVtICWXeMerla+f9De0Z2ot0Bzk0Bdsmf1h/gAX1ol2VgfiPxXDvfx7iDzcqdkJs1oTS3B51YWWCVAQDlPCI+/UzWU8WOJl9KnrTR6Ix1LAtWgXtJa97Jv1rRexvSi15qjeBXTCLKohsHIRdMihmKLSoiKPd8bU13BEpwcZMYpuJ61IQNUE518qDo5mZ5OJG/cdL64eHtCYYNHgEkI0ctgo7logB6AgHDzohVkwyBUboxSDx56eZxa7yZlaWNHOiZpUUBBx6d6bsL6XqwWZ0QWxUe7pqVahWttuCbMDpUf5MWsH4rEVCN+gzmCs5z2SeMfFghzAO2QqJnFEkbWvsd/IHloNkMxIpg4wkMra+68XtKcs5/rAyaL1UZB5J29jYSGTQf5xw+fPWNCa7KU4YHoBfQW/nHCEe2MWVTQPckWA4GfQ7hBpgfRDTLlJr4YyiBEP6fDcEqMerUKM3+xmpayxg9X8lVXoYHHI9Xb5sAi2jeJXYrwA+nRWNLiOG2EWXYLpPuWjzX+hVHhD+KEr+2pL699ldNn3omes0bVfB0u/DmFFe1BxtfDSA5uDrWkUCKAGECLol3XxWrpsFs0xW0aTW34+MW0fvcmcH5MLsiyXvW8bxBMPJ1O8lfR/d1/H2SgR9d01GLKuzERGIIYI1KtPAyPjk8XiyVfCjr3fHZ2ccbuDYoP+0CCnwXh7dE1t5GleOvlS5Bu6+XbqY2OwvIpkthlKA3xFZmU+iKl/DX0pkHaw6jtIDiJT6AHRrrKZpNVfYQf7KPQK6YvwlAqP4EvL+ZvgmwappfkwTTY4X8OwipgzXHsRHNqvRrayRO6+xXbhhl401CzB5UoacJzLfS1FKILSiqUFP3svkx5G7AgNovXw3wbZ7PtFYOqJw+PM91g2R4D+ua/q6FLVMgGbjYOPmpj4+egQKmAfwZy7LaOF1nrQo/OkhN3RN9RK+NU4erRYnGKnOKqO6Xat/O1hfUPXlqazHSeUe/PiPa6MUN3jSteOA2ZnYUh02ag9h6xI3UUpOMC2IQ4fjGxj+og4Ekb3Xxjmm0TKzxUJq4aOEg0QOvpYilMPwqRKEA5MpAbkbDkyTAKZ5v6aWaaZLUq/fSVbOEusTSUKvgLLlStee2rQef2xrqjXhA84fiiqU8LzXWPu5h5SVLuf4u9LWgrekWLIQxU4lT9IYbuDhAj/iI0EX6JYAsRXETIbmRgin0gWqbdUVOlUR+oP9pV+lDCDvBRmJ8fFouzTZWjPXqFydo/C8e86+x7Q4DRiqHGqyWTz8Eq3tFdL5AkGq8M4kTuATa2JIEDWx4DYYZvvWkNNT3D83mZaQwtxlHEl93kTZ+r9iwqyzR7DpN+raLvtZ0XdUuATiJv27eoBoVtJCEYM7ijgESaI3Zne7RVNIE3be3gD8eHjQNsqKatXEg+XkCq6uUCfEXT8w5HxROGsMZstuybJixeyDGVpSiOJ5D3fnDHENP/vmlIbXAoOaZYru7b/uVgoC6F9lA8CXGMMR3cG0mqkxIJZrJwpjxCBaE77cPx2NCHW6cksWCmEEPM79EX0EV6ozP0pvGOBgniCG2mpIt24NnPtHiD7NTRbEE4d7CqeHl+erEMJ44zEWuhIN7xL9VzO/mYNOCPK11jf+65gxKDHMp5nZ2U4p1uH/XQkHNWSonMRWrE6othEN6i4pjF38kT4h8ZUfaPYbJ3B+7xV/L4dJU9/VqevsXT3dhWlJSz5rxdsp3Y/B//5WJhdQj4B+GweOPAlRVg6VVMI70Kh/Jqdh5227dyL9zZ2rm/ub29ufN52FrjrXsbGe3bxQK9kRsJkk5A0Ka0F2ZYq3EIPPF16urg9UBC7VLnx99Pvp+U34++n9wePij3b28e4MH4Af73L38aPlCfEhWt/u5vA5iQ+UWS3TcKwFwZYr23jPmgg76v5FSu1I3MbtQwRtAKGBi/KtTENDgFfBMdzGuFshMIGJ6rK5Wk5VEpFyLlaILmY/YPNR9EgJ0qsjLxagDQ/m50I6NHN3VzoJ7Qt3nLE6CICSgVijpmDNHWT2EyqcX0UTsIZhmUJGRgBswzTaddyf3g8uRa5DDlSvCEc0pznXZR3sfHxDvp09X6kz/jhUApTDKOVAt5CNDMazTo6vApau/Lua+6sllYGy9hHQ+M5sz7JIFbHrUdl2Q1DXNuNGuQuOQt3DTFcM6YNBo7/YOIrZc9/CpacHRCH+TvvUIuE0KMGNnzpyDZqKKMd3kNgAxjWm84kME8iJhm2Ju4njPlSfiyR2QVANzYCuZIwu+gbucnbZg0Ps1w0v58cXZuL+J14TcL7IOELyh/wSt1NJmu6tkp+vn9nFEE9kQRQlnRx9DU3i5eLyRgUPrrkkI4paNVPcx7Z4/dwEEscgHEn/EaZHWyMhLWCBsqJqYnaG10VyRFhG2FP/bzCIJvZhO9dvOhD72OKhFEBOjKu1RHNy75Tes7OK+pbiwWDA0Xi6X0VLRI7Otud+w6alWs82200ONny0IvKOtih0rXnfRMbFNWtU4tXqe7ljSidy0nqj50fHzQCBGNHdqWcGYwVjkfEvdiatSUZUl97l3OPkn31+bVbc7WJjOXp4hMXxwh2pn5wqMbUywqzMSCs7jw2RKR0I+LH/jidckoaMktIykA+b4vyZ++kIKPzfPe52qNAbYciCWP+rsWKJk5iBSy30VdBFZV6zWTaGM0haUDDZdNDcYVhJv5Qm/1i+R9YpELeums/FFlPgACmlbxYiJxVJIEJ5mHFeYPOUW4x5l3BmAzNcJ2paNzlQ6qdHvA5RigRe3lRbGH/LPn4NI+rD1o1P8uieL6AcINBMBIkaZ4QY031dxrho8s0ZMgEreA8/M0VS2jjte5xyI/ijfFVMJTKgspDQwJTyVrUyvaKipt94oedsRgjtkRby4RmwkmH8L1eK9nKdbD24Y9A4jdlfulfSxxIhqfxotZqqVKf8dxuc531JMUf4taicdpVXHY2Lvs+Yf2VdltB6EoIDExrXcCf6Y2RVJgTZUOO6LsqHhMn3B4mkgAoyw68ghRfWNYYIcXgV+qK2xopU/YEu1Ix4mEtIq2wjQdae4Yb13rkHuMNJDdLugoII9hJSOMTzbxE481zMAe46dPa5Cpe5CkjN0IdCJ/aRK1YQ/ctfTBHQrpWwQBaLWiVJJXQ/smP/e0JVLHrhXJvB2ie1dPv52LBlVoaahin0lMF8T2HGWEI+t0vJ+2fSuSMxmG0l9bqkEATEdmIgg1sB+dAhScTxRZUfuqCvSBwD3Up9GV36sPKxf7LPpBctewvspUIcDIoZh4sLJf4EdTi+hw4Gbh0vLddMYOLOYdRu8DYYwffoSl9j0zXlopSLSRxDptpanfevwUXYcUq1lQk10VCq5o/TWviE5dyHqrZZKJgfiKbUzjxqpNdOpeZPJ03FeiqVFIZwXZt5xapsg8yyBubZM4Fgb83RS3LLvazNaVue5ljMqQI2h1mzGPCSy1tEq3oqU6l22hgMhOO95WoTM/g3e1+VXMfNjPVHrCF7tn+MgkkhjH096RHmbc7Tg7s3xwjbu4RgzV+U0Y1Mm1CwWT+DbQB3ExRXTT1rwczfy8XaAsEGy4Qupmke/0fP7lTmdt91ih2J5OJcczfRnd1kA7gSI1Qkr8e3XV0qmuGNb5bCDZQ331TBTEUe7A82/rK4xse2fLIxWKh/KPJqeyJREstD1Kgh92+4ANncXhI1p2P8jjMJVN0q1HKndHm7sv3MqOxvxYjRhCXQo+YTLbuGoYy3tbtZ6KYtaK/i+G3RkzaMfOpPR2Mc7Jz5wyzACGggtzFDLisNOT1pyaEMYdaecO/M5w7oyzWEyiYzcwvuphhDxtXle6mk2uP0FJ/M/73zffX2xtbT/e5D9PD27/pfvo76nb9VpduaRcLuEfJCgmxmIgJnIAKYLsWoOm5qK26wIKh4OcWkO7W4WqmfCT6smV4YajcAfYGTP5Qk6lMvoJIdeCQcG0mK8wO8xOyWKBAZ5JItwb9hzW887+91dbW5vfX20//f7q1+G/X9z7/mJna+cL/u+XB5qeiOmikl0iLufZBTMvUgOKAKCRA9m68/35j4+PP9y5+LQqwtLBrslII6MJpdejkxrG4+my8fWTvUTkQQw31f5PcGtcWRPhDkuPG8wIrBKpPigVjpEnt10xiHywd3vz3w2cz2F4DoCEnkEc1XPIVmSVQjl7VXG7QqaJZfHvfOV2ZCXJk/tJ8nFZvXxRSEqmWrwMhZu2DEIRkbwN0xxE2J4zOjGMlix+y97wttqN0fg3nsZyKqw5k0141us3vuhrN8mdfXsDVcWtIXVlsNGOMdrp4q3YBjJsjXzIeWOrw8bfWs6jvNQQwBko2RLeUDpAY9mIhmrc1wk56uIzopn5EQL0Yz2fbdGACFQfl9Q+UUWS/E95BUp6krX08kvmQ458r/q4eczi6D2fNHlAigtFIV+Gucg8LwyIv3MFqizstHVr74FA18HkEOg+yKpzvrKMYqIn7OgyA/xWVBnWmtAqvAKwDxG9h7QFndbR9Z3jafjwzkw4cQs8nUZczkQMiZnPk3NwqzT2hvbM5UzO0yr7mdBCfjt9cIunp6amHkdkQvSB70k/rge8xSwdm4ssCFW0x9OUqCOMRiJ66CHnLkx6ukTTWy0whKmqUFNotrw8mYXZlpzY4c10FSMEQo/hBCK6pIRgoLCClcPB0tAr8/nFLm6GNNlp85Y0C2TASTDfYRERmFpUb2tzhgwmh0BJvS7vXGA790E3o3qWS70tb2Q0yW58AtYMb28sbBkhbrn0EDkziODwe/aJ8RUNWk8vBCxNN1VKyK5vY0Z2zZYgGbdN5bq2uKlQPLr7bGKz7CHge9JzcEF4K2TotChbYsIKzCcVKroGQBoTHfpamPYYgfmg0Olen0Xj3y7nq8jgou7MQ6sHmvMbRgbqoKKSoLMT4TpDfb3bcCRKCbOtJIKgqdSKQ0HrlNGJEmhweK12fAYmMQZQN8RinmJ0477JlBsaMmjPwtyfhi1hETIxcqbWwN2xKkCuLX5MUYBbDnbavG7YMkaJhN36HHGS4oRRCSiEYvJLANTxqtJDBEE+gYc9jBA4pY907kj5rTVxF02LzYsBOvLAOel5gOI1kMTrsJNbuK2ta5p+z4sNkU8jmxUpL4zvaHF+rWA6CVtPIFTpo8CZ2o33HpPmDR4CVGIpFlMMQ5h7Va5jJalHd5R3FmaWwd32peHAksGRCFj/kOxx4/HifDadbBRtn/c4KFgde7XyNnu/wIxUper3LP6YPmoi6ZxJ1UPxFPfSsEczTv6SmsqJG8rvaQJTnl9van15gi/tQJLA7ImC2LQqbaXAEqcxMw/pTCM1MjCW6bnJ+5MssKjMW12JKcSsI0443h/IssCW/h3lDP46yG9Yvk/DmHqug+Usj/NEqEDOpoyeZIyW8Ibe9i+xX7qfr9n+1LWVkzOtiT8iEUkpU1Jk4NkqNVJSnNfvZ29rJCg4Op2dHy7q5WQEkEkBfIz17qaKesTHw9MFHAgwt2pWOES0LTdcvVLVDFXD1kY6CCCPiHIZ3uOB30HYWSKSChdV49KHi6vMRUf8h9apLS37u7kcRMVjpCE8DfLXwxjiMkDoA8LXB6LBkRLM9X2+0PzKgSKOZ1ek0/R2cV4fSXKQwVZ4kSd08OrQ8IlVy1GV0Z9E1iNOcKxiejU9erw4C9xwIttiwNhGnXTJ7jL0X5r8rG5UuPNYK7JeN5E4/XCkJjMZMWO3KWLzDZNK6haUqL+4B3Mngn7wijGV1UxwYPjKEXaawm6VxE8FcVXHPHUCXkoI8qoHNm7QYcbS/vNj4qlJDpq2KqNxkpTTY0h2QCamsT80CcJazGqN+rwhQY3kbVlX5uXilElMY6DEmnI+DwyOE50GfxQnfYfMkPOHvsX7NtA9T5EIdssWHklw8ExgEnMEwYHmpSlj11of2vPux1BRSU4Xr3OKiX5ypFHO4Nd7kpZV+y4vxWRNhOx9lDo4MBk3Ag5k3kKacMtwpDWPQ4YNbCCqoEcFUi3OL1YxjjpTGSNz80vG26I/BFJwWJFiiSTOHkEUtlBgVy8BsNulePrwIX5LhNEgXdCRPVZcDkyEwhWRdrLGjHNhTEtI9/UpsbCCVCsCaXjNGwmciEXyHDOkXBKMcVJz5IUY4Mo2ILE0Jj09NFMexsRsiNi3QfJRoA8xGW88n6vN+Lb1aqMjCE1X0Pt/u5hMy7PwP0nHca6PsfPCPx0bljoh3ByOwURQzGVxMHBecmTZZurXMAyOEoTI1uztmviNBVZOi+SnRb2chQaR+5HShwDu4WQkHyQfHgyzhAc3WUKcZxvm0aDlMgyFUvAFDMeiLuj62IJZwXMuVBUh1slZLfaNgR/DyrKujSJx6YorRPXYbBXmCkCyUEGb2CcmTOuCV+5G26HOUbHxMmUceA8twnwzWoObDaHcNsmESQizAF6ZeWFFeIR4r/XuGpmUxZ7yr0RkUJIDeyZwI0+TKovuJifRNlzAuspUiM0qkwba1yDWgFeMC42blNYVC8Ay9BLRVolWkYMtN/hqe2cDuEAri2kX7aaqKBX6XOZcF5RckR5wc00CB8mpVOzI4a6LfcktS9If3qUJF3EGoaN7w0EtrzGCUuws0ofXi0UrHQQHueeNazF0N7D2shXBa+D1NDx1j1jItSczgZ2I0BBwDxNO7MK8W/VyyaQvoJ51L30kDdqBp4eNCz/u/PP3k9t/r4Ys+WJ4U4cwSZq/PE+zrQ+1ijyvZqivY0MONZ9zRxxeG/TJrRYkQxsFa+Y7ZWKDKgMqQTPaC4XkDA/1KrDrVg/YwyZicMbh5J8BBzD04B86RW63n/QU2vVb297toycH5Jh/97eJGuGxePQvFzMR3jfF1/zoIjDTZfTFiy6+coJyT0naXz/t/PT17MwBZmcswiVXtVtc+sQsGpq3YSacyxwkw0i/u5BwTOzQq/PA1KENnIcry2xy52Iu2duZ4VtitQ9riPiThTmOXzC5DEeDOKpZfVo8ffFq5HqpGs2UAFaSbWc9tP5oIQ07dG7Mqfgwlc+Icx4Eg6mGNZSrdAqz13vZJO5yCx9OaU8OAtLbGhqp382h3lxd78r5jyKg3q1d0ieqiRQ0u307NQDGwrC0+rApUWx/hhDxlbKRSfg4NIDy0s5k15qeSUJ0HTQf3i62d/ugBV8v63k4sZartvsGZqfPvFi/D/cznsDtCBAN9CzbDuLq02hEVVGJPDv6H0VRpJ7lUZmDgbbiXlOTkOrd13KEWcWxCriT31RFqiF+wiFswmO+gyDV7hDp/tiCYLZbVYQlvkSOwp9Ry9etiYw1NufTabem1ugyD78WPe5bE4HFHVQwUL2fNXSWLDIUHHqv9vhlh04HHvLCfQa7hTrqu3YNwJ16/BV3yQMa3+J3QYr86X/9X1wQJSTSp1A6lDsWAEkk1hhk+TRG09ggUOpO4eAGW9OFKY6o2/mtXt1j0zr0SL1UphBZW52twz5jC48oKDdlTmCgFbNzt9O2iw4UYlSKwaCD9qNWSAZQ4bJ4ulHxCCodUwBDAL6axoxbTKEIxUD9bko0N1Ey0OEaU+e2BjLGrwSZpBUd8tNf/6PgBYWr51uHkQGtyqh4zDMsFPrPL+0M+emv/wVyE6tdLRisAdGsOVpeHKbwFHbgVnvvPxDXj3QeMWly2o0Go9SEGwfc1oqXU3qcz1Ies9CkoB2vNCZGUtd7Ok+wTW2ie2E560nJaN0VGDG/PZYo5rDBcMLOwhgBrUS6A8k0lDuOZ8szqyrasCznPUi3tz4N4/C1jWA9EktquE1PxKB7Xl8o6sQssGQcfrV8M5LI06fCo671yD9XjWrxXZggnKrpHOWxa3O8W2z/t//EVWuKnS2Lfzynjv8IANeDtXTdvT6eTE+RR/FWzoR7QoL6lieIYOv30Nq27vzzp86qSq43b9Bh14UxRgtd1nOayc7q+UV9+tEDUWJSpKgPYXKcUalkz5jYbuRmeraf1VfldqV/z+alf5ixuWHlYcylp8Ph0DO1BI8qbTlZSdoE95YeQWrg9R27XkKqA0XAKFoDsWQueQqPToPsg0VcXKxKFockZD6QSxE2wuxZkdJFyjkp0uM2fKiKZ7/77jdPXgY2/vkW9b4GBThP6poZXBoSyqh5PcLs7AxzTnbJwSEQTxBFtr71nr4TSS6ufEvIs4ix5H2RqCwK2UewDoxbLCPuQV5f+/YbnDaxI/M4aVFbtC6haCCXz9UxqS17VZ2jn/5QOJNXcffthe8Avu8VeCNqDfXcHvrX2+51+d6iIN+HTbwKcqgS+WP0cKQt4Z+o6+OzB3ljgYm5SMEfFQ9jHLoVqKJ8Y23Mwim/s9VdNm4x+Lz5qwrmezmEW2E4BR/OZ/CxXswtu6a8xabjX9ELgPSL+qAgjU9lNDAfBQb/aiY2LucCupgnZreG7d3qBkom6TiFOD+6OB6JnWPEPVYalixVbzhewPMwl0jvbjQkIfn8TLexVuf3eFuaHCYi346WmXSDgQbe0xkvLQ1dZBNg7rYmqXdksCkFQAqjBluwrKvikC/r8PJwSDfc2REpA1Xuk3sRKqzEbyPLsNaCCbY12rpn/YOtjckG7/ZsIHWlKL4LX8uYo6ca8L2woovjcK4u6Fg1MYBvo1YbtR7r8hAaiD6WfDSdnZo6w0nkkoY29N3GeKf49uF/eBMuqAmSSluzljJZ3JqlX3Mv600Uq8ysRdvh9H2XmcDSzaAtjUpFkMD/y/8V/r94ATmD4spP/8f/Q0ITA7V+fqNCWBU8uM/J3srUlGs2Bjgtx/7dgk40xfnF2TmHoaim68+L7mX8yiJ7hwz9FP+6DtUDN1oMJ/FkbTHMY2AsxrUfti7HK6pjE9KNygDKbbwQ8DMWKt91jhrDlfhcEPnWkVrV5fmeA4f7/V4a6dYoShSBlAPrxUgC3SWaA+lub21tDYfFZ9J2iyLnSFSkNm91mkB6+nhZAAxCUYN7Rs3lClByqyR+VIVy0qKb4RrlgaKZljNlQn7HU4hb5h03C67ZOszd1lfvHGwEluhwBkMAQNHeVXoPo90k9u0d5fjQq2Trly8dRw5kK3Hj/zS99vVEMn7Kcu5V28G5u/nf6bZPUyL2ksNwTr5zOOLxTNL42It5I9NsVmS9J7qNj95btTiBuqfYsP9k4xonyFpXR18VP647HvOtnfqkcy7gygooEI8iean7whknVDZGPGMHGCCf4/iRdsgUx57Us75+Xezc2xr29VYVWyYL9NCk7gTC4ULZdjhdXU6nc73X2i2k4fIyKocRFOEIC4LrdtHEOZE2sEputoVC5BXpfTOSz9eyOYcdespKywbuNxclCS2K/E7uiSxD66+0+jmrbknVMlkZ726zZ2oCvhUHzVIdNcNxvVwulkmvGTFyWHpTi5lvSn566cvdFoRG9rh7c2SDuKbekqZVjv8baN8FiUM9ocfF4YUA00sSdfq/8QyI99IgguCsCQyDWcxhESygfgj38mR3V15DBVpViB05yFFgSeame7sY/EUicybT+ON4N/cDsd3VlPqZMUdndwb+3jWYCWoXdtnIGdMgSxroh0X0HbcvoPxEvgoE5d8KXMOw55kpy/wY9orokszGTxeHTRJ/fX+yh205mQWvVzxDtz5B+vBi+XDN1TBe8SRJDIMAprypUUUz+8HFsfFayttdLT4TmPnF6QUXgZc6eLyL7RMZFOWap6gecK9kpnACr7YvfKSEP8wmq5PWhY/NxPP/EiUQBDN/XzeqFJZjwYMTeADG0RYTutyST+ASF5aMDUX7yXHL00chgtnyZxIQUHY/t17os99OZ29PVsPRD+GTUkPiJ9P34Q7xAnrMl5CLGFmX2wTZ1texsa3Rr+8hj31LgDsmekMS0e9veSGddQSp5v4WZBr802VPSTyI+8yYkuyLvYwMwwmVBAH7IrmzyCdxS1xyBH4FH2QbI/qpZ10iqalEKy0EuSlffCaPLTodAZBFY10OdHGKqWefkVJKbEQsgXRn4V/p4GVccD47DzfxmYlm0tyfF4eMkLk01JTwwMqBRYmjqKp6svtkBNjcUydMDi7cIc5nd9j9B8d75FvgYP+IQex5nvaPl/zJLGz/+J5/T+dH4eXvX/7u8eIsVA0nSMdR0l1L4gWl9VGT4Js/3w7UgC8EV4bdZ89i4T8zl+WIvoeqvyp//BCEvRGPDgWrTgJa0mTTE0mqVVW168yt2MDi3ZCOn5eE6nmCSsuWTlxNuIZ0hxhP8dUeJYjibBAJ31lP1FOuWWtoeFpmUxQIYyQAB0jhXXYOmB4QtgVRsRO/zukrvB06hry5h/IsOQLT3I2VpQoaqQCppgMlX1bS+Q8ZYw+XmVSHPwQm03BzmGoH4GHhJZEkEkfGglr0kjy0Qo6cs2X/0TcUe9m5mhPWa6RBXJnGs3cfkgxDk8N8m2YDaV2H3dEZy8f+oKrIgdxQOrxuMsXWyTjLJ7KVyyLnY5XEE3R5Y6IDFgh8Bv+2+cwLnWueku+nLKvV5KryxM3Yd0R2/tu5GVG2+0m4SqxubOVC7fKX1o/UJfZ7vJYlrv+4PXh9j45AHdrP/YcjaMPLUkrB0YZu4L+Dve4RJ05exXIymV2G7DZv77Lhfeg09+wacS+9M2Gv/eRnsxjfnxY3GHNmwRSkOx8SrcsT28wf2co37NAObeU7NL7Odqjfk7HEJ+9JSJf/FGMF5M4ApKFTS88mB0BMcih4rJrvRC6RdwsCzZtdSy+rZY2Y/Hj5RHEaESNMAnOuqIUMIuifASGgGmSNyFJPWwG5bhT0a+57o0Enoa9TzSWOYD24HossbZ7tsONqOByizJzdRmJ3pT+uN4oFwMBLCR49BO7FO01qYzVwSKMCwvw8s5mGMrm+PHzbKLiwxrLWK4ahjTI7gmO3XROfYtR1XFjcDaLtEnMmfr6iDOrjnx/X7ok9LxQxNfFB5nOj1/Wky41C8fZOJgXfNX0eNNJRUp5u3q+K0ikc99sXfyr7kz9cFVU2ThvR0WIPTW0PdaLTpYKz3f0kbV9lnlBU+sEdikPtqvkqcY2SaRLP8nfurmOqWbn+qSEEBhf/geu6fqpvI97ncQT7zNjo2pPyxsNFt9meFysi+3Ucs7yl1ZDI9Pj4CkdbGEwuBOAe+lXx4uHLh8+ePXnGi21WAJpLsjS80UAahhm0hZgwB/r2q2JH9MD9R1GbrVrHo/Zt3dGny/kFu9K9dw175ZP+1uzftrwXA3FZU9TB4Y7O7G24bOdKVXj3QwF8pgYNMLpGg1zpyDI9WwQGd3gxgeaSRugaigAs5jIUsDy77Xu771lX5bLPQ4JXwYOW5kUULtHJT66Lx2F1osb06+KbJ4+ff/PkmzcPXz7//XfftFUw+YltStBPVtJkh3F+5MXNleTyr/Zo5nr0x9dPXvUxwjBLgd6djTftC7Q+o46pbS6Moz/c3z7IB7+ZXtbtlz35m1vzHHozdAd/dwhUNXwxNCV86mw2Edn1RicxXU48REMElY2h3KKzF1DeeJ7Ci2tc7EuhSjCQDkbFxsMm6pM36KlEnJ+3S6gMqaCjb9BJ3e8sQCxcabcTk9uIWSZhYybHXbgqwK1CtV/j6FgrrFT8AOn14RI6sEbm7tFHuKvjIZIYBTHvp7/+7/EJMxZZwp5M67PeT25/f/CrnS++vH98jKhFdYETp7faZe396a//pwM5Hhwc0Lk7VWpKMa1xa2saWN1AI3QbyeL9qyCM3ZenFG4Ol+Ck8CObHs3OZ3QxjiU/35Ie2RU7NKn8MG8Ygvp8Va9puPLIzBLNOOk2wZwltIHbEkVk8WZNsxL6ZK0e8z+fUFerluPTetWto6/L7S/fgrx/1qAlk9pu3x4Kkh862rQigbh7KgGyDfuzvY8qRQGe1JBQw656dK1gQUQaOQk9FHGVvn8z0eMmX4rPCo/UPC42ruE6yQFspKTiUrn1kJkV4OMuUjIqTFnhXWAncx8o0gqc9oigIS5A2nQcMo6wxfGxhDGhoIC5qONRo4kJNOPsyey8cSinWZJWSC5M/CAhACmaSkcCuFNtnEem4vrMOiecdg3iT0tmtk7Hm7jHSG5FLmlZZVNKtDUDjezj/effvXnxu8evf//yychiqSZxywm+vCwxzYLYNu6LOXLSDiYt+PmbviCyRfimf/nip/jfg0wmJxG0E1kr1nme2V19pJG4Ys8m8ln4Vc7SrfHWLb4Ps4R/R0Jh6d4ZZ8dvyUSbrREeABgxiFKl9rG9lz9prJLj0G/N37dpzGcR/eiJB+dwglY0lxJyNRQYc6KEiBscyD5iFFnV7b0hdfFv2+YJoj3bQkiuZDHYlOWANb2SEE/DQNFo0Z5kxRjgY+yttilmRlMMwt5F2JhhKG6hRXWyV/z7V8+/C+cn/GVmx9elFVVV1uIqhm7qB3uaACRFOvUUkiI35hPJZaIWd8TkRO6I1YojsO6Dn9wAQSApGGIa1HXFZjDnoqo8Kh9Pou5aWIT0Be1LNk5fHO8Gu/JCUAkOIR4tXy/OwUOWjNfH+OBEJOluhS1hmMOeCrF4rQrTqeQqS0pwDEGzPQjMPieyJ+UDSkqKZhVMoS7Q3AWbgQ9j0jIZDeEFWgY8cPN0GoSGyYBZo5OkljFe2ZJF+wtnGO+84n2FPexoiJ0nRE+mhzVqkFvxriEuseERnqiJv+dKULjryVPxC/qbCaS0n9DDMObOvh3+LP6h2N7ihtphjCN/b/GmtL0DuXVOWAdX7m6n3F2UW7bLbXfKbaNcs6J4uzoZZBrttry8XkIOH5fqzz907vjChfqjTFyAzLHkG/sWMhX84iCZ/9f/Vw94JEGD2EEP5SMLyZ01lttJLHDhSwKjKGwNfg9TME7p1CvOXxIBX9oWrthQw4mGEHXq2pROL8S62exNVWdOO2uUF+vce6w7omMN71op1DIEo5ayw7To3ruY59gtjar3MVuMg8rpPanN90xZ3aJ8KJF203aJfseZ1wwBEzInpTzVYe8uVo+0lp3nU8agMXQ9jo2+n59usNztRkJkw3N2TJhWcyIeuDpR4aiI8LSpcYkPEour6siZ6SMPpoHr9PVo4FK/fHxLtVYr0NIzeEO24oEknVCrrKSGF+8gmR2ZkEFnvXIHFC9Lf1pWz24Iyg0cu50BTEK+bgpRw7jcmTjWB4jN9QQuFSnzb702WvcJwD7KzQatbSANGGDPMwc9qfX89Nd/HXTzgPq5+dA6hqJfaktfD4Sk67hpvZYyzTHLDPvYhrxxjoz5jAPx7dO8StTHsDkDttuSAFzwFwNk62yuOCtk6WFaN6kKJs5imJnGaZNDe73dDM+j0dU5fTruW7YWUNhDdr4PW6zKluU438WyNK5TTkLoGvN75kI/7VURd0ye61m5pU/rY5GZ73yfT+PfKsNjvCvBM1TTQ+ntW5AAJcyJYI5/l1A15AZ1nF/pZ/OkDRgV0cGLCQ6BaAH9q2BuRhe0cTijFUrq8gRwSLGy49lqxSKSFkzvWcsFMVIjvK/mN2mKEj16Gd4G+vnT8+ffvvn2+avXVjXcXaVDR1ccHv8UfzmOYBGaSfiIMVfkjBCWdsW7nGu4a7SHKVzHWZBvTxVfsV6xj+FuR1EDQJLhOcZrcAdnYijUdGkCuRauUAhWHKsFlgG9KKxNq64kjFJbhLzC2yyH+ur1kxevdNj/9OSPryylJFBvoWKR0oqCC6S3rRgklk36sNJn0ykyf2jj50F+OQGwTaAEzQd2FPb9W8MJI0cId5bzqtj4E5ZnJSgnGxkhRR9DIMlM82QAab2CxHq/coOCxXG7KrZH96pipyruVsXnVRFKfBGehd/b9w8qN2z4SgxuD8YwEA729N/N8O8m/nhjf2yFP7ZyoCeSlVfjuZCdseLlzaenxWlgLYyJXoZb6Ovnz5+9sgjokxmQ8+iBTVXe4rxy9/5IpqezqW6TSNCXeCCCacRmReIEPvmBtPQ675lsyZmBsYWuU2lQusNG4DMNVoZwWwAz5FJETaLkcJYbnEQtoX9eW2f9zpZLxh1uJfdtCv8QB2q9ZD0wcOvEhR5cFv/tPxUnwmho2Q9zhNTy89BWeVUV18NIMDqnuiFVS7NaiOFfoMOwNuUPwibDYIiDROVjZfMf9zMmSNiOLSkQxBZn4ZquFN+zgYsyPPwPFV79cRjz3fJcKy9n4OJY8WHE7OLKM9wdvGnsDkSbkmNUYOyMbVTEeJLHrX5AxfkMPrjLynGo5HyA+Pi54w76Zd2Qydbq7nA6O4OrgYJLQQs8eTttsW1NI5R1DEhJQMYVJOp63toa2i59hJepgAyJZEVaOod0mjo3Kl4t8obFv+KQPstv34rKGr2TY+h6ujIoLHCOMEvovbYdvsWHh2SfK4kEjkiGttnQFX0Go2JbiUaf4dDQSa5Bm112XJlnJ+nRCX2IgbxFKe1SCKEVnAcvBW6TYTzIXIGt0ZaPqSY9hQ++gJ/CLFRYskZ7cDIMlfyAy5lU9Fk8WTMfDBLgniRDhi1sGkgU5FppUnvnZySG+j0WQjNAX5lu3k++1gAJvY2q+B7F0ZOd3KuzOYWnjwvSYPFQLHzyubQNOz4e3g5DibrkOA3ajLbyubaEf1Bz5Stmwc/Df/WV/1RalEGGyfmBLpCZXvrH4mosE4RoM9m9YZVD2cqEgj25DwSJPrxAjWN9E+b+2j7misp6kyZSBde+ghNfAVwtfjAiqArlO+NC6wKTGQst3cEJB3ajL8P4SAd8nlhOnhTx1XTVJJU2CH5C3G/wD5waFMACJ6D5hLO2DUaaTtzIcslgcj4k55fsk6GyERwk7gAhv9Rvma1hChcixCY01F/TcexixdNIJAjxfOLVQ9gRlqgpxG1M3aGKP5l8WPU5oGn78OeKDFFZBJ3yQUaj4g/qn1VL/7qcfppQTTqwiZDGSkkNcBTOpqPrHlesPn1ib9CD+Fz1RSectJ//VrmLXNX32tdrHyxSaTqqM0UyyXz7xxZ7EqXwNYAOkRIq0go2j22JItsXniMqgR9d7UbyL7I9kPFKK32d7qt0wDlJWv1aUKsjQ3adKOvVCLsk7PHw1xUYUfj3B9eylvijlLhOJeT21LqCExUjgwGxOfsaE3pTCA6osvEnu7822KlWwzofSOY8dMReUQQSoFDZUxRHKndc55ts7Kp3plP87ubxw7I9XCmxWsv/FmKNoCA/L6pGDin5eBTPqk8kOvYq33dhFhGpL/Vx3fXvH3A26d9KGed8/wM+uk4fXa//CMRyzvc/ZEsMUFUqjyEkCizfEHZISiZzdzPRMwhmwAhxb8yWkaviGGukUBkYaM6ABCGnu56vQhdC83+b9XSb69PXUiEKIzJhGrliLP5AgeKH8DqymM8QybUD5uNKC0QDcrrB0Rp/5d9/lb6/I9/LkJSs54pZWeY94dkVLnhDTJFjPMYsxjpznVd/xKvrDkroOfJ9tRUA8c46m6tMmjajXAFbEmcVL2VWxaCxy0ySiw+n14t40VHtSxBRgqDBc57ngjDOD31eVA+X07pl9P1bnTptNlwVV5ncWhWb4MHKYCEatV6DFauI0+LGw7ZEloY7zqpoAyLojJSXrQMAdhZIZzJVN9dhMnt50u5V+B1+dPJM95NDy0+FOecomQh8Mwlt3FWrFOUG1mxD1DLNQoIzT0/lngf6WphII/LNUISXgp8lsWfDNBHLlafMdMKURh5DUclKH/xV2um2RL/A13anQi6JKa9MoiuCcpitcXxJMfiDu1HjZl53/HE4TxmhKiZmeM5BQXs7kGEN1qK5s+09R/SVSIx7hbiMR+IU2BYUD8R5m4U8cbrX1/batyMyzd460kvu67EF/iGfbUp9JMe4cbt1GQm2Kru2yvRDV1sO1y6NfVUwUFYL669S+7+XX1spYml39lrX12ECIOnlPHIujxNrznaIkHNkdsZxmuLJFeCuirPp/EKB4opXCL0SS3Pv1hjw92Bo52OdkXx4b0TSI6S/UmoqFe/mx8JR2F7hUoNrJTSG659I38ymd71Dyozkvy/5GZjLoioEOrORnwdIYjebeDQza3K3JSgyCVpUliZFupcWy/2rwK+uw3+vgoB0vX2Q7eh4egzd1u+wJTmUZgYmTpMjXDcap8/p0ZjXPWp3QVARBDDJtWe9Xy0WqukNwzoUP8cVORlwHjMzAOrQvqgCctc00lN4UoF1iAGgT/R5vXiExBuH4X8+6ZD7RdIN1IUzGAqBGSUZR83Vu32E4G1VvB9mHONqC6ofVlKWsSj6vL91AF/S8MfOwbB1ZFXZ/sTput2pJrT+M6vx/bpe369tq/CuVXjdqfDE5ud6fcd+Rj3pvleGkQb+yNh3lAys6z5Z1zWeX2fPkxoqofLYzSDpzGry36g2iw2IDOIVaLEN7aapYsqwiLfDEgxFw1Je4+e1/GxlbfS8z0wcy+nbQLIb2NCyCUaacMWOdDOZtLe8KGg8Z3CwYGcpimy9cuIlm3ZYc2DT0h+4cGpOvVt4L/bveQSRc7SiKR64n7DxNqWGwfC/e3qGXo+RVoSCMqnW3DWIouPFysyV4ySRnAaR5rRS2cWMaFAUTxYQZciJRAukcORlzdx3ZkIQJ09ZM7hsNDlbq4uj5aJxN/jOQiYDRrSJ7YoXqdXtzUsQkarIs8PAVUqcMxPbO/DhYjIjgAt1Yt0cZ9wS9ERgpHpLj6021r10iH8Ng9hKLnYOf2MP7mMPPHynvcNV7rb9iIjI20O6Sv23/zQgVvbbacaTN/nII30YfXF9Wq4VkAMlB/hNpTTbvXnJsFwggg72NF/QYTtD15aZ2I0lZvO8gJ+KeLtUxoPB9GwEqXJQ6SwjR9+akrqrKr9JCTluW/OGb4UbxFao/eN3pNHwC3/TG46v7I1WywH8ItmGHlDNqhwYXHVP1+oV8m77GaBTQrTsV0bPk+lRYNSa8cYIX/PCSspI+JBoHlduqm/r84ys27czZ3EWi17X2Jb7GCQ5RxwO68bdmiTJCnKyIavQZHZ2JpYoSkwKcWwab/E66ruFiSbnoaWHWRojoQuBYycC90x2EOR2eoH/jskP57ClSXhDozzOLPD/chEamy4zdeKgkYQbeIxoq9+txOaqDPNwEdjfWSHRQmJkZcgnFYoWT6gKRTmuaCGmxczHfGg9tNemWc6zeEpsNCzbJomKrwXC/1IOOrkp140BFJmpoT7XTEci6bb5XB89tBKIncvmpyczqhtW0cndsz7ZPaxnl/ZAl21RvNat0+b4KDn+Mqf3n28nVLuc465Rnrr/eeUjpb+oiq3Rzv2oNA0POdxhLnueod0Z41hnwBRnA+OuolUa49vPxKolpsSzk/z7j3x4QqMkg6bxKLF1TKE4qct1kz6cRMc5vxoMi/Zr93K39b3eUlnBSU8F8VId32YQGkaMllw12R8k6KWVfqUQN5CurgtgG7EypkBDRYeYjWbfXhwoTpJlz4Dx9HN1JK3Pexgltg6x3g5BgLXoDw5HQV5VkffrAqot1LaJKZJC1yx07QqdpEIneU1bSXhufe7e5PeEH4uj1VXxQYzKpTij1XprdtoWKFjwOLT0Gfom5ulre3qNp6SPQPvHl/ZY1tt9cHxir3Ql42fap9XViPvy1dlCLoxP5nY2Cw4lSoATEGPEeMEW///sEiQ9TDXh6vsKdAO/zeXbw7rcvl/t7FRhf22N7t0bDnZjsZfIB5mqQTht+xVm8vjEmsFFAkbtk76SFTQp3VdYpeNLvj/DSsOqcXwpJWOnmyDuvpvGbls0lFQF/5w/6CbaHt3bdV/EFgBKJr2Vv/yuPkaz4cNh6/GJPh4mU6Gwv5Z6Fr01NetxywAdLbXqtV+rw0qfBtzQ6XAKNDjU5LwSjQxzrYajHRQAzUBpGQRJlye0J8tBsbw4nY7DKRgOL0jAo6F3TpK2Xin/L0WKhaFBToG73uuLif2Sp2PqH26+PKGktrHvZRLlwys5OqctT6Va/GrpomF6MNUdR586vVQ4c6BVVDvVzvbdL74Q7wBU8xY56CxH3KVgjODsf2rpPN/NxOMl+nRhyiYIhTqpTxfObE8fLASHloziHpujzWW0mSNMRxTVmhPa0g4yLnUD+FYWatpIrzMvPxfDp82mSD6ZEziEitOcpTvURGqzxoANT6c1/QQi0qH4cbrVfvHkO3H7k5mHNhMjHQs2Kwa2M9qpihidid9fkNL0DeMr+TdA3jEl7tvt0d3sW7oh8tttQBjrt+FvfIupGBf74Ui/e2AJ2cLvL6uwiAehgFwnRPbt9PJuq5f8zZbuul7eXdPLL1q9/CL28gvXyy9cL8M4Pve9hDPlzs6BhYdhg/w2rE9cGoslnnnxLfr6VtGdNFpGogIXq0f/DN3NJFgxlAjVY9fLris3HkbnRAuUJNLA4vg42k81yTs/UxFSMswLeWLjQ6SeqxpBcstoEtsaQh0xL2erKa3vMUtK3BsSWl1JBObG89gNKjBJ+yPLF9lD6z60WkRf2YHcJdMZtZ+X9fWa8GWNus1zKUiwghdo9mdh6SRB2Z5PXVVxhLTtnZ4GeWcu4MCtQGiq3eRrCMSz+W8wfaXs91cjTmbFZofDtZYZ0bHskeHwU9hhiQcReiSVMzakL6Ka5XKohR/BisdawwMrUoylnUrX0RWw8zEUScHS+yx9YPGj4yKOUruL1rkdBvkZ9jDnQFW6H76dilVMMiCrBcMSPgpPmgrLXZs4Mmy9nd2CP7e2fh1+/nqoOiqjudNmoUlmI4FL5L2mqFTHTEtSCfqbR5CqqQZgzBeOfyKT+GqqJsCMYaaciIxKJIyJCw4+kxBey+OI1yOmoCzvlN9Pbg+/bz77+ztJOXoGCTnLqXi2v33Au0KrFpdY9xuqzdKxu5gnxZhheBm/sS0mMyuqNFpGdQwGfZDycFayu6XWePQlb/bTMDPwv4AhPrw5xlrLMQR9m+VLMC7Do7dUa8FKRIaB+tUi5Y0+leSsy4RPwHP1aPoM6pMmOgnDuQSsjjEpys/WsMmo/IWcnuy8XalKJ0QUE5yzP4dbfYMMo+SvySik+UxU+QApYa7m4ahYgMENR2pYro3A1abnK/JfPdupBKBTS6YDWMynpnR1A6Jy1Bz15lF2C8wPlUN7shFziK/xvXBqX6poyjVOErvOmG7ZVLXlt7P30/kwiTeSpkxccMxFXY+hC3FcTFZueaKg2GFKDXni2bR+zzQDxtcR/972LYZwig3blGEBAdM5r7giVSHzWnHUOauXU3Cv2IeFhSDknJnASOnAukVP1K0oleeqB5PPc8+7D97UdSRgjQjxDbdbcJIrMEc++Moe3Oa/0YB8hPBfLX0dH3xlD6S0NCYc/dXoPFwxmgRahp9rwCZmZMkyIckg3YagkIj5BYBqDXlAQEqImjZLSWNwKZqunoVvvwkLUzIWveduxSvhVsX/c/dBf8kKyzUC29gn2kX4dRCIiRvqs3SosbBvgUNd2x0HnBHGihbwQxLXWNLWMEW6SaW5GFLfHQarM0CAvv7j9y/qf0T16JtQ74YmG4RUTo22rEorH66j77afPHopddDaGZVe7jmMjVgk9htb+BU8xO/I5Ax1UK63YWeWdv9eSOTm/a0thv6xdWqP+PPtdAXQ5ovAwTilZYSKOFxMrof8/Gl9NiMi20pw/x+ezt5iJw3gbB6kiPTmkaIC4KWwl8HHae+LPtLbGgEjEZ3d7WozolrAzSox4LC98y12SKFRlVezg4RTCriVkvqnTSqPzPxJbdMmtUZ8RIQV6EYOoju9y6YJO7ysd+hgmItSm8cqwREfbXaTd6lF2Rcaau6mvPdqpHIZnFsCcJJBfiwkBfbVOGucJHXdegZyumR4Dqb2LFwvw+lCYAHaoJTt3XaLEBjrAjHVvppK5nws/6DnH9Kmced+KTNlAQOBU3ZtMpKaNJxdUnYomjklG/aML0YyRvn7yv7w+qpPKfshbhQAOyza1teHUexhZE8QLRmbGMRPnL8+17ke6ELpKoHwhFexTXKGiqOIeIiJfDLLb+/PHj568uzNy+fPv6Via4f5tZRXLE4nzTPL28ohaMQVyGQnrI4Mi0ugJ5qrzl5TfguE1/5gu/sBt5zOxAssYpNJAebpSW7XLEQynTPEWc0u9Zx3C8oPEghGh2+VHDKnn41ACxsdlx6VWVwyUvW8M6QgjkAsZ3phUVGmBJFqLBuFWiFkweMLP1UasNA5XaWScyGJ7e3afameZkFGw4REGS3VEH6dQeDiNUTIgaGsnRA2rxHiclIuTgs7FMcrE8xjYPlUFU8WPYkAE+KLWWbZKnNnmkKVsFRbS5BQK4GSQ6FTGPToDlWrDWxUPOcIEF8e/aVoa6eJCmA00qItarpWKTyJagxmK7OVIYSS646e7BK31FyHK1err8uCFrtlF3bnwf3kFJnujP6oCxsVctrGOcVHs2VaB4NmRt+gf3DbVM13c1MDptsiLwQj3k6gWm0pb3gf6Y2dW8vlMqlWzCRbo+17kaXqVO5FLDY1qMlO2ZNoUMtyjP2+n1yQVEh125vHV/J16kixm53CV0OvCTdJNgVbZZXlQm6oLZW7Hg6zAfBSmXdf9W8WwLWCqqTtala6/t1mwaEGtUV4UdgMWFFgdZR7su9VtJHKW19exy9xMc96O39mJ3lnwmWFArOcXBxNy7K5OAsDETiDCwTeYaT5MPcXOGsWo+sD/K+LkK1iVl1rVmi4p1WhpPWtlgteGOTMfFBs4bLQ6Ui4+yyAogGWv4jnQ/4bU5H6tf44jsELLUzYxMTKQzuWNv3TWp8O2wgn+xJNdlDoMLSzk6v4IJ5ot7l1uqROMJDWYWZljTZbkZrnixUn3PUwO1KBwKcdU2A9g+jD81DvRF5W6H8Yqf+Jv/DI/7wdf8ZvXZFYgX8Wy93uKcdnB0MdVFGIyLiPdT/I5EYadYWPZNTlpEWEgq5MOqxXlAlltcdxS5Rawr8eVqYDVAruLeRFwZRkmXhia2cf4NRYoVDbSIXrr+ji/iXSO+gbNTsTRofPNIxGwP2a9oeocyTCij76xC9dMduHkezDHcEKF7X9uWk9PQxdrBM6isqfDAkLJ9fVrgmi9uR619gML/K6UPzaUIszL5DoMfhxr+BdHynFu+mngRd9JAlVx4Iu/ijRji5fDvsubttb2/e3H/ebnqOPUxsGKSEaZUhGOfrjbr8roE4BH/X1vxtx2b5AH07fzuYvaoRG7uqlQbqcYnx0/JvukZYNs3ceswa3kSrG6yOoVGhmPP9HfQNuOWNv8llc4y0Afwn+73XuP60uEi3fZX2aLcsvcjTpV06skBocOfnK1C1bgqP6lI9/qLQLUUfw7xdU5A0EDUvvmy1NI9S3I7NVb0fHbOu70+Bkri7DflROTBBxOAUyeUsgk++N7gke9L1Hj7eePg1Pdkb3A3/OLq9G++t0VDb5uT4qxu46tZSqFJMLzvqrbO/I8hM4U4zEKisztCW7ha5iV5u0vfPzdUhZF1YMzGtpIMTYotrcLOzOe4nBk2NFd46dqs+Z6HZxf9iOy3MxebeJrxNro4981ecwFAruDDNdZtetxvRYv4461C7LwDSrjwpRIXDmb9+tMIbbRPEJ/73vOGbZyz2TK0xUdQiGMWFjPFSWwXlp+EjhE6Opt/3hlHm3qWwH0bc9+nN5rR2esxdd6v0qhbm8YqEhtUf3/AqEd6jgOr1T6JAk3cfYBa8ipfSuoQs3L8SX29X2lzvVzr17YTW273c8nKLVxGTyj+sm7/fpJu92PY86Na9TxtPjoL/hm5yddj6hyfWK6o8qvPIIfnOHzSOwOr6pFpwYXr/g90G2m74P2z5XwFoa2SiVAFrvETZEoL/HlC44oJZwyppU+kBoO+6eote4zl/+UV8C6uiD6mijwATW9XBl3frvFUZ1ThuVzYO2JvuiPB9ddeKWruXFdSduCFbOGJNMtjqHFrMc7Ex8JAE+/coLGzZz4ZxNW32NYEEbooPA6hzKgYRe41wGnCv840Q60z+0F8gVdKxwzGg0JhDSREFIKmMmMDXlAGaCX3zFAsgThLFQsGk4cb8TnmnnHDTMZGzD2NRst4eSY5wUCnm6fgSLfDTbqvtAF3igMheGHxT9jH7zbVWP6L2fz0tKN5axhD8womSU67M7oNTBbgzfMhmN5gaGX+HschaILXV7tSCs7BIg0YGQdVp3AOnMz5LK1gzi4/ukZ3wuPVbk+N09hEg0f6RnJpjPwiHYu8m6XzkrDb5qQWck/05iTGnH7vzcBLkSIzf81HUi7f1dFqQmctXjWhbLdD+0b8Nf6OXiMglYSsfPD/8cPhoB6LCEz71toAf8blyUWkDTQTNPRjkc4nIMG/mHNr009dk59KuvLs4RltwiGgCLFg7lm0ij8OgIR0BoASuPJ6oC2lYz1cpM2dAE4j1sUdEGNuGrCVFD2jqjGsuW3SfM2WQvNZqGHB/tuxD2TjnMPlQD29k9hR4ke6rhLfO1GH1bX715W5+/gfuK6GalvIBryN90dQzkpT3MogZJ+TY5+Zzsi0Ef9axCHZiiUucozpYkhf9qj20/WPN+qJBBNkMYZqXuAmw3pYKLbd22uiC7pvq6bbXe97Vl8kMyMcpLzlhlcxB9yHyOcybeeRUOxU6ao8Wliy+QxWD8HDEOdK1CodErPOLyTOeT/NWT+YQvIkKPe2cm6DdN/N6w6p+064llp6ywNVj7bFyI6VZhc6iosiepZmoPpUTWIhWJMrzN/n7a8A+qOHH2CceO5xdIReW+fnkxdwNsvZDRHOTOfYwpbSKW4KAhwNMMGGVIPCTmNOt3kA+QjKjc4L8bw4TaM5+dXZzxi/CWu3FjKBBrWXV0lzq9rK8zxysBStSjV5pArp+GJhzzN/JtOBBehtzTodShupav7fM3TbhdH8HBro7Fkp95ynKJfJaL04uzueaNmwDaR/wX3yLhUnKt9zgmbF1SeqUed0QEm4FH11gKj/iBlWl3VfwwH2UznjBABvJ79yOfP8UKuK+4IplfBJ/kZM+P3sRxgIAWkiLXFfoOj1yhJGYJ6bNevxdYB3aB/BUYqTT9wMbC2Gft7zi60vqARsa8Nt3QuoE4LTOX95+R9W9i6HmVBfr9f//x6Z2nP/1Pfy3qK7gt0/BcF5b4SW2Km5SLzk8Xq1HxVI2cydzsyNeAZ0KjSCh12ojsmIiVhGQJjLCr4NJmYX9CQ+L/2tk1CQyXs2JW3+lKr+euF1XvPkkVZPWGzgzTdmH6HR24bWq47HF+6aeQoHcFwYPJrGjGhx8jOh/jQ9VyKzbvUbEhqizBNd6gpmbDPHGnlqkmhoLwdis7spYFTtC9hx0P8OVU47ylqyVNeR8/Pg6T31Efi8VsmCa4xbmp646ff4XbTCdFltrE2Pe9QpIJi97FFIRLkj6T0hc9r9VTkQlKsWttlQnEYj8oh3KTsCX5WcnqjFMXw/m+Dakz9v8zWf8Pw05nIZ7s868yNVmlFvkhU7+xn5zyl6mcm8nEfGZzpBAKox26/GFcJWlSjBhtTshRKxd4UGiXhBdU8BRLdNzqzlMtcwP3YkSm1CgtVNpU2jOtSr/TYjdwu2Fl5hupOSbkY449XQqQVavqR6nc1jAXKHR+TEamu7DIi0qSCAlRBB9EJGh5kbu3UFbTlpkYEr5UIL+E8JPq+Jsi5HMXjxmrRJVQ4KCR1ZxMa2Y4hY9l1Aq9AIb5l/eDFPnrHVVJaPSbgGHjuu8YPzTEuKgfz94qIjJcpsx9f5bCpJlXSvjMZLYSwBG4hYAzjTW1qOy9YyqcD1JyvtMZTsOF/1oZoDmdi7aLl4tGwcYB4SDsHF2EdRNnjkoC0kQBfmxRPA4HMfw6Eq34UpUuMZYCwTNOhfb0+XevCfG+T+SQU9DUw+yPqvjt9PT9dBWE7yAsBqF6s5kuZ8eSjvLEXqH4b3t/VIVW0/6W2aJR9DX+KL6bXhYvF2dC1hutZ+Ek5BNYUWIFnCF0M2GUyBOc75KKL8IVpKNgHL18xKWJx2nMj3AjiaVYpWtYF7iSjZ5MFiopZORWN858im0kbKW4QNmRpZJgg6PM6GC3ePHwN0/evHr9x2dPSEJ6JE4vU6UNwIcNomR1em1xDqRsoWilVqEDbOPYKKKvFlfhj/Lp737z+5dP3jz93ZNn3yD0ABGH4lLFdLa2K0jxU8qt9UTx+JYLBVOwpn8/x6zUut3eLutzAuCfLxfh77BTAA0sXv9zJPnbVedRsjWeVpa1TwAbVwBzkMsmv+L+DV3mtLx58fx3gYAjBv3ivD6ara5lhregd9v2srTMiO4qOriWtvlWsxVoA+6ECSm0vlgF6gs9kKOwvcoptiOxC3ASGKo0RFQkMe6AYWURc3GnsuTP2a4y6Fe/ffjyyTcMurRuMxz+6rX/dZ39QlPhfLdNjsGPiy8rqeA1f22D08/fhVLJzpgS/RCsLZwVGrjl+oOOqOMtQh5VRVQ3zeztvPzxQ5X1ugqdVpnk1ylk4nM785Dg5go4dYx0wfor/AWM1hWJoRESYFKFC0mp8HY5mxDraUHIW5JS2NU/TJcLrffqW8Rp6jR9W1/FvwHzGmfMFbp2ha59IcvwhHCxMEFoeiw4WRUvAH+QwUkIOPr7jLQ7RuYMOa1/w08Gv5pMp58TAUYe/yn0Fo+/rOt7h9tAdQkPvpFgUaltZ7Rz78DPlMRdDjTlpY/FxLVFIy/FUfjVk2+JGYZLqoRnHdUN5tBiJHMAf7w6NV9bj24X/oPGdJQMb5WGvmXOapDOr3cOd+KoXqao2McPnz1LGVtHq8WzxeV0+bgGIq0V/w4UxIJhr7SLNDhoBxjFgD+ey27HZO/cq3RAnHqr7rE+SrlN0+x1rnpmNwt/xuHFRXtMJ9fBr7aPd768e9fXk03V4fRk5u5BYDsQIebTi9WyPpXjnhAwpbEq4U7C/6MfKZKShcPpaY0UUzGyOAgb9fFKhBCkebBIOolDp6foakkUG+0dK4kEed/G8hqPdZrJ5dsTze/i7Ardx07E+r50D/PCaXKQEVNjNiXh5K5ezYLUfXrR6JUveq+L858korajeDe7SSoRy4VVNFbh8/kEfsY1ZzuUWLaqiAI1eoA0jItzAmdJJ8aa1dZ+Czcc2E9XMtvb9lC36F3bnzrJr2Ltg18dfn7veOt+JJqk3ZOG5YNH8engV9PJ8dbx54NKh2YF5dejWsLl79mXf9BC4bsv6y8+P7qRylWhIZhHrTt/belyxTldrgGKM7ARL2LNBi2rEQoJss1kIocWycm33laxZTeYigvG7DOt2/cwLtl7ugXK+NPt1BYhu8Fyw97TYrYsLHY3Lkt2VcKU/froi19/+fkge897Hl7er+9Ovqzzl7yv4WX9+f1pvZW/tDvsR7id4+E6L8uLuakICQGMxQCpjxWxel4vl6igYK6XUFptRypnUEpSW6JLvENRNgJcUGsLGXrhVwhQGWwe31wuRPZBVg9NTqqyX6NZH3CrhErYk+7LC7L+yb3pdALebN0zMg2/X4cOhnN+dK/4oMNXyfqXCQz3IJosvYTjhjQ7m3rAMYe4rHKdBA/H6G53LZO7lx6dS2R+HjspsNQ+v8SLJCZi1ixIQOqkQPxzRI8wGFcMv2LBfC9r1w7BA4Qmog1ZIzCgOYPGss7AE3dVVE/wImTyorZDiqiLZtoJcY78amnrHf7UlU0H/6yxQ3Z6D//n2Y8LE1lKBiT1DSDdTjX5tilZpS9MHC2R2FSaMoqE0lBUrjaL09mEjJ4O/wSyEwomPosiSkoXZP8YUhu3D7vwbMoLNbOiCH6FBmSzUaV6MeyHv1aXTIBl8lVjbVVRngJD0VZZfzpWwg+dtHBIHJ/WbyO3ule5vjgR8fmKiCHrT2gp9jTU5YWScBN8D1lAq5ZC39izPiHInZMW/p+fz+IlnXHuyLLy81dTXvH8bUbFc1mYdn2aXTysMOK1UkfO69NpEN6FLmRCJTC8ZlA9I5/m6YixhJYmX6w5Rw+/ONo5Ooos+vn8mbAP8bS6u119vl3d29Yg5FTqm1CrFYJLl/1XHO3aOxLk4kLx46kYvdgoDmgBt1yBE0JPZnqKiLxqKgcFN20SJyc0dMnrNU0NFDMxA+TdNJ4nHpEC5mClH8ZQ/1BR2G+811fx4rpw+KbJCxVaY+JVxosN9C3pUoNfXmryv2NxLyclleWrJF0NvBSD31HMSS/jNSI9ClcFwGXyrBjpnEELBwdBeYih8olXBZnyQw6+5bSrW9glPxRHtLea3OksQVH1XMhFCwE12mCl0xDv1fhxw/TZJY/KaLtYcZC8uvAvu5Dw6yjSizrdC+WfNP3Z5Gbq/pY5wF7LhyTmgR3sLC2HuqilfQeM3RHb1LhdfEVWh1/KqwYHnQ3Qq9kSnU6dKWCq9fq6+nT2Ti3Ls/m7iLhl6tQRADDeTjdoYDM7sUCQxbVfLSSAT1BTAQ4xvUzKWvVULF/ApPd2ir4lvrhSHOKFHGYdITxH0Epavj1iSuEO/CNnwTibqmXuTj8/uvflIGpzqO8kYpSrxPYC/GTTQ96kotOl6ENrP1oYpBtBIC02oO3U6FiqJG3TkDurBg6WBQPHMSCUlSrUpWpVUuu86Cy71NFSCHHKpdyCqX3bRToZ0fyFLfTiKqaIhZhFpR01iBItmZKcRyX622kbykVbNpqJRXvopm0UxEQ8Py45X7k1UAE0usKrlO3EtSkvzdSWtGChIuFbt/Y0wCg92iseLpf19WjW8N9UeEjYUf0hTkDvaZd5j+yErwUr1GrBk8wglPcaBSFmn1/JRxq+rvNvz4yhmX1jXJS0HmgsE5ZdTErhDziEwSYk4KHy6i9/QVTb/g5dxUzOTdpFKeqog4Y7Rmmc12ABg4RPIu7oanf8kG+kZ8+f/5OCc4SFW3O34Hao2vsDbvt89NsnD19/vA4lmapn74kJzDDcRCsTdkKmKnJYhydTyymMbXkq6DaQfezmMRatkdqfNCkoQRBB0BRhT8KcLZW9mFqauDm1ik9H4dCcuPRJ1BmKxA0I9uMLyY9GRgBDyaigo4qgQbmeG/zcYRDEiD7GAzOeRDwytX1Ddzq7gNQOtwl21oLqI3YTe2Y3E+A+BzHSAIjQzYxdRs0Y4bNaobyjLR9csQX1JPzo7siz5l+Wq9Ljb84N2deAPHQQfdWzvi8r88uzXoQysth/o2zgrLhR6omuJYbNVF+KKaKBq4vnpQzLh4EG0GaL5TtN5UoiAReazJrz0xr+67hBmRtWqI7txYyQdkZi1s2mYrnV2BFVCmh/DKzP3LD0QkCJE5j/MHjwYtmBZwpNf3N83HbQ896g6jSxZETgcoS+4i/vosqnwxbejUPKkwpQaBImphWFpCPZK5Yju9pIhndyKqS8Rhz1BVwDf/xgQRWOodO0E+ah26HUUkwKwmbCu318dbBLJ3ZEVO/Bs/2dOK5njuzSO/DRzWJbnNr10fYBbPR4cK0jdy7tas+/JqwLQDkRtb1XXBN/5t3t2xa25aYB+Urewe/0Aunf3gX2+l393S5GLp0Nr6/9AaNNvBcHrC31Wi3fA1OXuDXHW1JHjESxJm0lQo2hdh8foitnBXK2flI3r+nEIw4wLotu8mrhZ3BGU7KKLmO3EHeqFIUOxx9yMoYnfv3yt5GydvNU5bJF3bUMF7Wx8wNrss2hIutslTIQStAZ4tTCLnwbDva4q92WPloslxQZ+R0p3fNBhqW+lq6EmVHsH5xR9akk98KUMAWVTmBpr4d00Leh2lO4izpfbNaceOOs+YYrJw3dKj+1Kcde5cl3EuUn9bDWEv0Z2PQN2A2t3FaX7kovAwtIkwUOuEqTNOz9jv1ylVNts5aJhiVBJgdWNfgNWlJnJlGTpJIDI4dv7RCFfS3UTfs6j+CBLrvMo9oNYeeRP8raWDI2XrglXIuanT4KIhxweCdQNU7Cyks9YzmCf/qf/7c4BMPbE6S9Ss13nIUII9JJ3UUlFjrPKL5GTT7NDey4KiASKekjQD2RHxRSc3I45xleiYARRNpwaxLhdS42vM6zcELQO5qxIuakFIExXwsb8hJ1K9pnrgzQuG7GXA3l2nHa9bwcPUmxEu9DaW2d/NMx2Peeu77fja3gVzg03xuzjbwPkwFX+F5eGybF3n1Nf/oeaQUoI9bKJiv4zKphiqV3OCiYYOndUOvtMtgf1ezSwH44H3O85snfycB4NLUEZhsEXORh+NNf/9UpZDeKJohtwEUSmVWzx+uFaiLaPhomF94kgZvlotFbXZAVrsW9pRKAHLuOg46l9pmiO6k3zXKG/ZKUikyCqLfr0n2GfiIKawWRs7MF4GHw/KIdk9KieY1dA78Kh7voThq9yujF55JUejmyw+iVFhtqSMxoCSMEGdJjqS1H3sVeeDVdlfiOEvAjTC+8AcBQw57iBMEtzzrzwJ6ODmNRdkmxCnvXsZazKqx4kxzTnAuU5YjBMRWRUUXqNxZm4ngTlyKZRmgrUa2LLmChmLHyn1XsCYGpROMuKEqgLbVHJP+DnIFVzo/e6olgqGI8j636/3yUcjkXjffO4pUoyLx14LujNBSkB7u2oYB4IxDXvDMNBpNFnRUseAI7HVGq5wSqta500KyRzu84ugjm/+kxa5H+Cfemrs0XDJ0z9/BazF2hlVPJOSOrInOoh4fgissYbcxcoXFBxxIE6YXx/ubl89+/eNV/5HQgwHf/xkcWyFdQFhQgu+/IuvGMMshs6apPlS5PxoXWLFJ/ZBBZZJemNcjiG9edWdy1oTLAm+hLiO4O79swwivD2RbU76xFobw9abrblq9GhLEnSmeUyG5FHlH6uuOBwForbWVsNBCD4LT1sQ68c0jwkBdUc8fgK+Eu0AYIFRG5sDYDxZmXmnLBSAMPBBkt1EgyTcEFNFvYxbTr5i+esGIfoJ/iqNh4KHtZYwaem4Hvtm3yDcFWTHwl0fBusfEy+d68Tyk155vJJ2dDlRmiKBGL5XVEu5dtJxXP5umOQE3qfEE/0sjI1ETD6eneIMzJMTlnittXE9mVXvfVvTQTQPt206splLqfcgfvxCVILmFU8u1ioqcbb/2D9jX7SKHG0vZFYCazPdm1UKveH2RWnmXX5JPe0wH4Gw2C0ZJiVOg8ctaig3j3VfiqZlWvLpo8ZkLgrFpKPy1BK/I4Kw6A5yRNc7Rpj2YbdKgqUTQJvVMsZjBQc0atZrgJa6ZQZzBCpghj3ldGdtBxnNNwVGN0EiIKVqTXL6BKs+e671Eg12soszWNg+ch+zfO1cDigQDX4idubEuDiRtgGqWbD1RKHQcpJ0qr+GFCa9bTD8a10/QPhwc3TG0KSpRMuLXsIhr6vKPkOKllqeoXa+whM+UpN5Ikd6qNrefJeaCMf6m7fqg0v7yjnm8UbLxsuA+rQhC2gcIHhTtuSrckeGU+JSuXchzLdrLlkATITno4p7r9oNNBoqBziPApSKsw/lLiGxNrXmSTEl369XAjz1mgeclDZ8U2IJRNa8ztsLL8iAuMn8N4Pf6DeBfFQ4LcGSCjk3XpS8RHL4uLOJ4xlcf5yrxJLFhjzCg5Heh//b+5dOzJ9vbd4YaxSxMcN6TEhvpD8eYPq2DMp1JzjuIqhtojC9JqW7OCYUWAXFtEgb3Btpc5K23VsEy3C2bwhUZCO0xVA/PlaTGvLzlRNfTJFABVp+cndYZvcoQ799tDvE4KL/FawFocCYThoJIf2/7HTvrBenXZvDeCbI8bBGhBbVWZPl385IYQlWM8HCP8fzzqKOtK2sJDczw668xuw9afUxwrs7G/NoEoiXfK24r+c0nRKpUFxukS8xqOzDLKexbt1CfO9QtvWcqW3ZZPR0wogenINwAYDskbX+7aHU1zJHAeT+P0JwElLYIEroTLdbnxRun8Dd++4RINLbmtoSXXqbUJG3ij6n/GSytGcBRBwoI6oUahi3WoquNcToGKO40qzhyfJwzmn2aMfJfF64oMsENIhmKsoqjw1qxecnZo3rxv3gSumB5wF6XxgKxLO+zCl2/w4MYzMgg8gs9xZ/+f680ftja/PLh9522FT4eyQwEdVTp15Xo15ZvwxxvoKd/IfedGxeSbydOnW28Wx2/4dWhJvhS940DcYEPD2WaQDvQsuDITr6s2v0TnmC5278ai1cfqOUTDUOLU0ScbXBU8q32ylEJFxqXbpGncVw/Ra/W10e1PYhoLCYbBhzV+wwp6BuVCi/t8KDKJVpDqbxRpjdOqPGUycFXYva6PQmk44C7AMezZzr48P8BVqff5uGhnbkKgxdyOUlHbSwIfQcB3lPhw809tYkyv//nN7b+8uf33fIEDxM65veLnEEuWHQKYBHBiQqzdm2ZwsOt1pLzgydwNWVTibTmashlR3BsSPZktVcW6Nzhr3uDo+0iRuZCEYYkyB2xq16Br/DrjGhb6z0LIwVOGE86ChfsUx5nCuqVE9jHf+1Jwf3Zww5Qglpc9UxVxM1KdLxYg/uCkV0UsI6pfKaJqYClhooDdtiQFbquN8FC/lr90/3MmNBw6fBGnYtijEmZRKfH9fMC04t9LbrBkgVtO61MAqlDOCf+KCegW/nQoDPgJcOz2o+2DyJEezum5qCccE2epoy39pkVPkLRpLJA8hcW3msoqOTq5EyTAzxgNgONjdiDRxpWGb2/B3gB39yEiZ6FP4fkw6UUJwt6Oy3UplubiGZarQNHbKpNcGV3p1Nnqd06HpphZzKptYApawMIA2fbl9PiiqU93o1PToKED1B+l0vw2wcrolI4VQo8qDLbCyFpZf2BalrUJBUA44Z9CfNwJpXwycyXw9QNUYiW2UxaTUBAXFUDPILnXySyyutxY6vNlCZoDPXprs8ALmpAfDdWbr2fOShhzkvQBRsUe7XOPYhD8o4VSdHDQtuLOzqYtMXXEUa6JHVAgDtKvEGZEw5jUqxrPooFXFwxwaaNCYCnHJE4TyG4IRMjCdkfhfDs93RKqCn9tR2+4j81r3/HIGOpSrpn+gNxnI5W0ADtUXAIcjJgTMEF6Ol5t8VxUzal0blzoq0r/2G6X2Y5l4FSxLzOChhLlKtny6ouQifgniPAvf2EvWqrLKwUF9b3X7n7ou+rHxYbKrHk3nfgog1lUmuM6OO7XFGq0TZ58XY0YgtBzkxGj68aksRPOfr9G+kCOSZM/bvUKGl0v2uiaYJRFio/OV+UVA6SvtolCRnXGhH7gbcaX0tVK/FyyiYggpVYJpwn29gpRxntzkkiKaquIYumZlKfBqCPrgfmhotQ04L2iQcT8dyzvcSsuYU2IlfpCIbJYTg/RHqs2gNRLNetlwsopktExfgJkoHoymU523UGUEqvCATi6oFpAm3W8M8gYLycT2QmLG+sIPd7Vxfx02jRJU5X8BymG5e63/j8SSaF5Vn14H/TQs6N3SDrymqK7m0ExqoJyTFESIwq7DOcFAs7Uv9TQedvmmtYtPtTKNJtrRPGqkASrLZw3fAdj1zrIPQs3eMNiAnxmCaqjRGNijj5hUgDKL+Gku9qKv5GrL4wkR6o8FTx0632uW/W716tVWzCMEW7xa3NSSFYabrVx68piBrZXkucLSe3aKmKahHm5GfvOpf68djeY161bi0W6cLeNM1w6MucYEmS3IsDxjBI+0IPC/y7yUOLGatFtmtVhRayOuKfzjmksok5/qIEVgkxhOrflZAUXsqL6t6afklW1Z5KE8Wrb2ryYa3NQTtPDs4yRkiB+oltOjfXr2nmdImnVnOgtpDWFHYua1xhKkppJo+ytEHXsqPy0Psovm1SLlNW5eNDBs3o1SvBxyfbAz5LSWdlTOG/5F50abmW4UZIQsWsIlMMhV+uS5BRC5NgypDODQPuEooIBhY8l37kUmhJVgV9pXFgtLqrNJTLimjluOp00MUcaDzp+XesqMBc10XvGXaa7K9VrACBOP8ZXxcOsMgAF6m0vJZPVoQYYhsZtqZFS1Z9vGybWJ8cCnMJBgqemTX1i1PVBhmMBuGT9uzefD9p6KjrXg5B3kA3bcqEXG63jcCNy8/hmo50rLIUtK3pOcrSh0scMkbxBmctwjA3U/cBl7zsUnpFSynMJRu5Ko2ELnPEe72yNLDziHZiv5XL7I/3PAVsmuFdmWHKZOR8zl5kGEMtzC6cyFJROvoEIaNa+6ksv7ML/o7PlhO5Gk05po9r91L7uN9rbg7yrMSCsUpOUGg9oGCjbX5rVSLEYhpE1yDCjtaKMVooPMY1DmmM1F5tPteL3x3286/afJCY2oEPuo3cizdVekkoymahqUwiCCCCDTONyLao535ebXR4oGyTL8DD6H4oIlehIRon6vyp6x8oPbLV6aOm1REz7FUrheJxMbeDjVcX+5tXl8XwujdGnEBKaRAKjFvRffK6JLZpwccT38nzbheOsFPjCsbFB7MOH9o4QUMK4KxQs4pNJXs6XPrLXV5L5RLCg/n/23m05jiy7Enyvr3BSkoV7MhAEwHsgkRwmkyyyi0ymESgxKQiiOeICxiAQAYUHCEAUzVQ5N2nmTf2gH+iZ57GeS0/PvKn/hF8yZ62997m4ewBkVUndYzZpUhHhfvzczz77unY/up0Qj2h9lZI0rVk/w04m1xJxI209O5ift9ON771IoblRLBRPJ8mj8q1sOLAwX9XwGwUzMaSDiEX6vlz49mtU/ZJunK6gmSt78Pp0tmrU1mpoTjlUNpAmdhBLjNwN7usRIhZBOzZUVfHgAWIJRwug28MolfN2v8qbrhvZsqqTydFIgqLFFdyJdYV3khGZVXJaeu3I4hB/K5rXb2BThpOSQiG4jT8YRUE4ci2KrsmiWcQh20DFQm5KaiQTAf8Jx1ypt/ObUXn0sjwRZzJ/Hcq85GwuBu29ph/TcUzepm7JfGaEUdP7tsOda16JXJKdI0BQ4Mzpx86nbrryD5GmzONrZjYG5GXIdZI+BN+wj27XnCGDG+KuNjR3iv5+8KCgDrLv2oNK7hh4Fh/24nb26bzBiMetVCdtzR5asw2f1rBH3sou6Yrb01lT+xH0FcE67f6jxCSMovLqlXoqdFPbl6j5yoa3gPtFE2xQA1I1+J714M1W0vDB6fSIc8rKm6qQbvDb2nRVRwcjznvbwZRXllDiwZ30CEnJyYILUQXn1eS/JxYAJgemMm+1cE67IebLfMA1Vv4KT1a3y0z+9Kbn2n+CxgG5SKQnoLnNp1DGeJ+Zmhtsq8urWKrpemP10AXHiQ7t7Uor7txMjjXMo4qxEoXSSK0IbOzSi5bJHcWUMD5dcOQHo4t5BPBVBc8fC6WcChRDpUuBabXUvW6CZf7gxmL7MvGZmEk8tFgxQLUyuPwvA3S5AD94BMGScJC6jeaG5iJayNgYsvRjpOq6jOGnDRFiCLr2SDwBopTJioLTJo7IqqgLgRLcCKM8zO8su5P9mUCmgIQzb+/p8qxcDG0rz+BKxMmrliOPX0S8v37k2FRpjCqIV1elxsgDSwJOWYWXkNzoR8cny4t46aCwbxNNRBO9UjJpUVfJMjn6KrZ+obD8s5m80ceR8knXiaWk3BeglPOvrAVJxaSWtaiaobqU7e3/kQQYsVv2XD0gG3n+wYCDPzajYaSykdYUWwHbrZzrWyqfky9B/CCul+xGBgyATzXhREqGfD7eFYM5JFucK1yXotrTBLhy3LazpFIYsrCY0pSuif4Iti1YkBdqaWM1da6rLj8lws1IEnfkI/URS9R14nlNx1Qv6SiDYDFQXm/zPchZOdOUKFp2DxkbsCfFDo7G8o770dGrOTzEr05iXsY3kdCnGwN3lttNTt5cjPI9lHFXubvsN6k/5UNWzKcP7hTMI1HSyV7m2K0U/xAlqvy5DlWjRFHn/uWafwm51ba09TBL5xxr5Y+BP1VWQZf9Rk3d0OJaZqnWfXlfg3VCP9wIH9JnTj7EFsBr25Sxg67GNWxrLx/Kv5a6OKhgvWCLv93U6OoG5Ubk6aDHn7sXfv+05/KP8rwla3K7MBZXIrKSIHSHTFq0CMu+V9qjNmDufn0URVNrjWGnsAJ3lPnRBg4utHYs/B6pZNz7b7FlHIMoZdyDG2m2HGYl3fblgSvvPrizhaJr23zN+m/o36k/B+8LRvdFnh23zZsjo/pioTnfCJU5i65D09LBC4EQZGIRFERNgupVRGk7BYOl4QauOaZITU42769tubByMYh3szs81SeyqlKCHuCjkyI7cNvKbSbtPV7aZueM3sijpENInO2e3pTSNzRL8DfyE9lyesv5T445mlTuIss3kLaLU62VDEaTKetwD7WOtS+oo+krYpb+OjeuHISy4rEVvS7L5f6WVY+lphlc1EPX35LVuR4Daqkd0zFV1C7rDZ8bTwrWZ0ps++V8viVOZ8FatSLGSEakxvmvMMu3qDX50WrmwShYO5/R1da2xU5+EUVlM4zQm+Eji7nuYD7pYSDy1/sJb0V9670zVhnbL4Kx/eISY7vt6aZ5/acphWX13FaooD7mjQ5ZXqyRFIpni/LkRDDEhX3UD7o+NwidrpRlBiOBnIbI+lRtYX3c0YNl0lVvrhSzkebVcdxrRX5YQADlmEuSW3e0BayPDCCcj6wMa5eyhnPepxreo6Tjk7N5NwBjBu0/n88UdR+v/JYezyEwKIw1dmUv+1F14gTrFwLzXgN9VPvFP52EU5mFKfWSPpzzNJyZb6WGXS2EAc5lAZCCs1KfSeIiusmcqAhHHBZEjo7cDn7Nla1Eb8uOddSjSjH3ZqMQ6AnIJWXKtfGD6XxwJBycKPTN12DZYy7JQblQnLWRIm69UKA1CJKGOwdyu7G5RvwoBYIvvb9CSIhQj9zlQMsLJ0NIymWqn1SJKHu5axslPYSnAdVlRqLqdv3JuaN/SObrpkpT/m1np/nG/UJ0yPGzOz5lm/txuwhKZtGwkRvDn2JRCwsSuqngSrIt11w1mxtFGvCnuWlFhcwauSdYJYvDTxJPvGyAu1q6kkoI8rye5pjlRL1aMFdIVwBDVzQZZuVGJo9jXREm0c+RDsO7fvHIQaCFq01MBhhWKodBQGkBlS8HN/hVcNuflFV0jnkItYBsbRHZr2Ok18ORnS7IS8u+8UgnEw5SZ4UgYol3mbmPpj6eZ/TxZM4CRxixNkX8BZg+993eIvDP1G357FnisaNEC1bp+Vlv5jgUSnz6WzbDDe4st7aw8MryuJZYHlgCUcFtK+nZreDw+REhUhteoX4mclTs7al+cML/6aCZ+NkLkrqMN73tOnL+5geezNmMUxgQyCEQVGD1B8/v8cTiPcqLLemvp/Uiu88FnTZgBcwCx+nPsef4VFgOvA1WJeRjjJbCDZ2Fv82iIvJMuEMVy3AH6CTgnfLheBxHxC1SGUmnDqW2jJ0zwAc3d1E2Kg9jqAkueM9IyLbEd/BcQJFARRKuiyxvGm7VcAYTMZgabaDox0csXjkQZaA6ubPQXbmUx+UFC4abJZl23aJVcE3ZsEeNPe/2NaZpI9YsDI7qu0syR5tw5GiYuz0RghH9KdOPj1dMP7qA13txGe3Cd+iArg1ebyWi2slUAXP2FIVpIrLDRUAGiQ4/zq9OFGUdV/5cymua8Bi2SWGcjuBYjkNLDJCuySHSsJ1R3AR9uTSYbJiU12xyeKD2Lf2tDtoZJHnkJeahPgchwFxCnULisZVkL8cNgN9C3qF/wf68wFd4vFWHBpEeRtdgfP15gsJVSC8eybPVuHc8weItEyxaF03nEc/0lDEbeQFeQC7Y2Txmma5DDLze90njUlqjJOYA+cpnHrqZTRizJupN5WR8YKBy/43g8+saVHi3uB7MT7O5GA9UU+k5u6kbO7ApAXzJRuu8HqVdcSoxDvaMUdwCUQ1MzVyT/65gfAJDIbwEaoytR/ULX+JK5O74dluLe6dp9zJNNlczm8d3IYaDAbtv3EIvB+/zm27y/nLnxt/+Jedo76+K/W/+svhb9+TmoUA+sv66HmzKLBWe4/AaONnvyUMaihBj2lG2J4moMeUVq7tmTscIZ64zPSzSmIVVDW5rF3UvReZ8WRqZfnRHp79NXBbvnGir/yC+TrbHeeAoCcUcbZQNA8ID0PwH8wWSJLhdl4sQ5R4WzEYkfhghaZD3fIp8Wc3QUPrEOGbuUdkIvDsTLrU4O03EtZk3kndOEniHEEBBP46u+m7pUcAh0YZDAHKEFOHqIngk06Yin8TEPKt87pquhxSnTQ1w79anqErKVnRT0sngK5NSAv68OMmbt4r1E1YM8rS1oxbkGlnqKdemm/n5bxPxL5EuJL89DsOjqZPrES6NqjrhuU/Huo20hMgUry/x8DFljs7B6XLZaRxLuU66tjZyQWCQ0umeEPb6GazOFePEnSv3SYUbEJvtBq4I+zQw/zezzfjIUZQRL6KP7CXu+B1CcYgw1JNELfbq9WiwzCtpB1IPBCknQt2Vf9aLIiAT+urpTESdPj0tIr77VNAAxUHDEAHRFKBWj0bWD6Yelb0UplJca/QlT4C8q0ZLUFa4l7DiHs6EQj1JiYPR4WT2k7vkcn0AM+LunMO6KEILfOTmUMdnr6Rv9m3cHFv4FJTdfgbEV4fMphtyLsusM7jZu4NZ8y+0MN/fCO9NQIi/jyabL21elXtgbaqAcePn+my1L3I8udEySwcK7SurZeJ6FZrxoGgZrjqdFK2LGZzeepPZEQwOt4rmqm7Ec532xv2zzowLcZ/sIXsG2q+sHWTbGL+tPnrrx5beDd5CJntScBwIsMDrR2gIX1l/NgpdqiOQCjluMb9F8hJdK5FnDTTJQuPHAnBOASpSIrGEezWajikyk6Zj985ECU6q79Z+OjkBvAyw7adR4AoQzJvBcKbzQt3zWWC+BD7MvHDGABCjT4TYpy2IRMEmhNOTpIvz4xHxbzT2jwjTCjNWT0BhIXzGgdT4xBqLJe9W6aFSQh/f5384O7X0enxIBCDz5JVM8CS+HcAkvhP8lbaqQ9rVWdGDv8ITdwcxLvPz3/0vHWv4O213traWKGqvrqHmpYVJz+BAAvndrcr1cnmdIiiSO/Ly1E2BG5GBcNHOY9SQD9pdnMr+OoEkG2nlfTjUMA30s1a3s3j6SypnhqMID6W5RrBmBHO6fuEObbksaEOyi3fTh56+kf1X6oglUyN1+mCkCFo/gGLY59HTvJmQj+mm7UQQOpPBvCqxUgBlnQC7cjg5xGwQe1BBRg3oSbbzP//7z7/7v91j3cWPSz6QJAGq3p0vNJPmsWY6mAnGk8Dpi1PTGXDERRdbGv/YF7mnll2UHB5cZQCUqN1TVUDs9lJLMktWb8nkchrtSLOLx4mpYtKjmBHDybB+pGpkT1YMd2TKM1HNwKXSRAGAd7655yZt7fMv//D5l//x8y//0+ff/a///B//+d//8//2+Xf/+9rn3/1HN2uff/f/7N90HFOlJzFmQPwhkvY8KkqqA+BQt+UjptbljNI2W42eTuflMr+Z7/3lsLd/ozg5v9kbnY8Ggm8vtmHZBttCTkK4Pj9w5QHUiQphxbx7n+fu5BySuzsfKwQr5r5ZmGjVlJagZuy4ifn8y+8+//LL51/+m8+//Leff/nvPv/y33/+5X/4/Mvf+/nq9AjW92qcS5WF11R04qn8/Lv/4/Pv/s/Pv/sPn3/3f/lZbfnWEo+jB99tE+h5DbGcrFMf4Pd6CKeSTx5mnfWNzVu379y9d//Bjc9//4+dYLCf0r0VLMA+82ejrZi1xDRFeh783IufUc+jPdvWBosVxQShSL12towVmnk9DD+OAgEAcZeIciDjZ8pbrCeKYUs45WrzW1A2lnvUs2mQvdLnXtniG6utTvvxLrpatjKv7PVfpdu4VEEiESxEqUWXT30voXO4gwE4cu7mRSoG89OPCzF/WFqmn52vFFwaurLfYz6Sc2sTIFcAmLQ8fCko447nu+U+Do+/lT253sNY1mnPdkcPMuWyfeZCr0IP4sFxPhL1FFPOqr9xOWDs0kq78FbN6D2HtJtgbCqGwSJYOICBGbmw7T5//Juddy9f7eyafTu6jBZABYC6/VTSoEmH+lmCGZBiA4T63OjurkfoUssRDLEyPghTauDFczHzRo5e18yVi74LgNaE28C3oADuR/DcUPeC76JmPW8UjMo1rxlvccHW/ZBo9oPXQ+qxsOXKfbtNpxB5/Q1f4/EN87CInV9YY3lQAbb+2/gLhYW68aHu/hAQtT5qTTI1SLaT4AiYY6LxTIyypivwG0s8zzcHTMVB02QDYoLfDJC0K6RX8+lFYg8IWnVnnh02HDRVx4j1xd34D9ataXci6/bVutPk7vOXT95xvYJ2k506m1tWOG+QTlIDQqe0sAkA9wj0X8RIth4U3cJuXoDlVt+k1gVs0rBHMS9+g1q4/5uIM4w3sS8gwAE73MfL5I2/gzYjP7B4/JHTzRugaruJTPU4RWEpQ7D+OxQzRT22lKy6DDpQl2Ro7XUCCLGk8AG0FXm/co0iBZKgogHopgjezszO2GCzRMalK6oyWmzENuvP3ewtjdV0TPGKqroeA2fuRMxfSGK9GiwnxtPQBBSJM1gNPyd2D2UxS0lBV8kR/d8UVMedbflLdK/fbke9jnQqP+da3jEmb+NK9V6PNDKXFJWxIqrF60Yaww93hdfZRDeCShJNMTUKiDRNq1sVIhdG/tRys8YJWt0hKrP8hdvZhZ3OSKCFQ8X3VGMKC07SIrJy7IkCXzX2ghL0RPYkxWVE6cGnvC6KA1JdMywL1RCBg9FU+MY0udpVE6rrm1BmQCILRc0h4f7d7AXP6kWqI12Ko1ws0iVx0JJi6g0YVzhzN/mkWPyGUrLJnyi38wWq1YZGJ01OuBV4hlS48HTkhbt/0M+InODRG2iY+Pi8EA3grxI4Eg0YoQNYshUO5yKZnmX/6Z8yBg+XPqRaLmbvR7Gcn8SJtelsUIk7BYmIUN23GjEA1bxuFXGWSrLGUgkviVnH9LAol5pYkjeTCfjTi6hmjZHKW9U1s6HPGh/5DJlFe01D/g8Q0JbknVdPQuihEY+guaf9KfEklHpr18PDiWJ6iCJsWcUHwSTdNjWQxuZ6OB2ZRE2M6d2+kKqRVUVzaRe+1f7r8kRKS3eu8G9iyUPLxOuVdVsGukf3AGvb56YcIK1I8KDBTtXWH9lyQt2feOApPmep/hJMjr7wyhi21PPem0KDCILrQ6igBdT4EeUaqNID2IgRHdCnDc9rSKB+eaEZyZQzXs5NcTJWTaQk8fOI9RrrrvwBES0jDJX86fPdViz42P8LAqURHvWqVJewsIRf4QgG6ivxGOhq5bFu64ZvD+B3hcdVRLJoC0+gQ2T5ti3ClUsYy2O4a+lcFqykK/3fXP9v2TzgiiuyJvzgi4TBPs3vQkNh44dm+q4A8Q2pYYmmD+9u+8IhjWkynNe12jf4Ra6jfAg9OzXx6rFpnk+n+f2CEpybnzexA8PGOoa15nq9lr1uQ/trTJujHV85aVI7Bm5zR14gCqCfi8cfPAKFrJmpTIjuQ/GKeG9ucpynfm2hSZMEJUbDcIxoxnbra9dabkNHJ0VPeu1aW4RxoxXOizX3UPwWRUuza3Mr4Z8ymBv2FeGJ2dTDWtxxGFRS76YuWtwDCQfcru8S1nCrvtdgBnxWX++ICYcaZxdxrKw0aHbFhNhIZuA4jl0yHc+60tG3/czGeJrfU3/NXdyc8TjSWeijDjcPb/uZn41dtnjbK5fDxPg9Ao7MfbLrKnrm/p89NpwaRc7K/DH4qBYfe+QYhT53Ig/CWnI8HP+AahEsq7stfC0ulu7bFyyFTfnJg9bETMeTEO4qgmJd7LSoy2q0XGJaDkbjOeUTw+3WKNU6exLk2n6o0wNGwdncvAAEsf4sBocw4GeJCyBaVNFtvY+7LTdwpl4KCTflOKj8sB/dEEXMBpAvnl5YBjbJOtk1ERWjAJwezUS7wQNZXZv1Ap7ohZVofgAoR7U4MGjoZw37jL8n9ZZVjY8TrBTm8WlLYiJOxRRo9unVll5g7fB3bvW361B/eiXqixZIrq537/eftwQzyFGL44s8JZOfbZhCIeYI6uIpUxrFD6nDcVyeJKMPwaOqhaD+QKMxPBK2l6QMy89A+2befxRPkFJKoGpquDIrMGgiZl8UhIFOddz/ElO6Zq1wC/l4fnxy6tqhFJEP54NTZPXrwa5ZsKKnLB1cauX0jaH0s23cwqrDL3N2IZaWyoIhxNe5XHpPnhBxTNbrCzUvyW0xofmNnUqcoS8wU4l2xfSDEikSa1d8aAhriRVnF244yWG8il27qHqqVIgyTqIV8lEfWIDavaLwijhR0p9MiU+LEbmbYctDGk+Hdr3IjYJyvcPeM2h3bIW/QThioW6n/OTbzN2YcZVx1liUqF09aQoMfNhdCVDpN53MmmzLT4mIKI5idRJtdv7VxFftebnPhFxkIQ22EqxOFTkDq1GXrDni3o1Wvop0jtQmWCZguUIOkgvEdW6+CHpLgAdY3hkLSSayVskUi8z2J4ZwFNevGxH5AstIbwGoytV6GyMtxP4CNfVHUzOjnVMlhrZjITB2h0yQLe4Mftp5dH90Q5O8z3yTxxqg3n5j+TYlZFWhyGqkHveLBpNfSu0vj1vJUnxztcTVItUVPFegrINiEAWhgXRj3zk9OQGESNGDL/GVd037Tu7Kofc3SetN1lUAT88K41c3YeBYwWHci/OKAZaqEW6og+OiP1N243RBHZMvHUt1zuS0+fmG/v0NvnyLgkJlwLiJIeMDC0bRqI6bHCwm1QlJZSTxCRNdLtW1KITJKjqipQtt85GrFam5ULmxrcuMFanRTxWJ2xIRH3u+uSOer3JUQ4rUXOcW8aNdsqgb+nYwnQCYhls0fFG5O4M1fopnthws5hXV8RcGHuTjruAwt9IXb5UbXuSB13B0bDjKfbFL3ouGRx7Vcl/rjCc+MsQ5n4+Vta0MC8aDpEgco4F8xFnQPI7Hdf5zXWK7coWqYQKhpUZANRXFTElfi2FSqAYilIb5Z+XJ9X3gZSrTSEZKy589Lns3+znAtt/QehA8GItkLyLlZijuRCScD/miqU1tdQwsuQPxnfvrWX1rqYY2xsD4OhOH57GiZBdA+pC70cOQcg5XgQJvZYeLiUU/Ug11VsrzMhvOobP2un3LcqFrlitWeWT/qCOQoeXaw2476Nh6As4QIXlFVSxOZ8nXr/1v94bHq4bw4N5gcEVstf+ADeu5rkLPd05yl78FskfSxq/d9/oAF2LaCgBSmfGd8nTMzQc/QeKxuJsPk9dN18VyzpjrlITh8siflHASmVWJ3yFi/HAcFZSoSMZKzwAkpJ6Yr93ViCHUvyQSjeQBCvKMNyO/TeAU3IZ+m68XMC0jBW9tDtdrc/gXHHo6h/obswIiZC4oYUuLNl6JkBdYRWUfHPEtuW4DBLsr+t8ZNBpbCYPk3zAEYJgmf7HWU/MWQpmj4BcEveiqBrtA71e1zevR8IxE1WnMz7WSIDKkMk1K1ii6sd9Ggxqoek1ydNiDosfTpFS304KO2dDZtXrVy8r/XPumSK6ey8xOTV/1WCf45mtuzOTLaGu1WIJrHu2HenOe1D3bu6Lj+qKLtE5/Yr1msdp4WHnroesF9XiJU5bGNwptT6DUch+OolkRoMudwNyzNIzcAEXapSANdZfaW6k9WC4Y+WMukBeWi91afj8ZDkfAB7JMgSZxrMLqjTMnr8hlxLMzP9liEChPJOQZ88OYn1jbsakwKL7qd5hyFKvvJlvOfzOf0GCqOe3rWyqkuk8PwQqIo5AtxrLmbWdfhtnavQT11deGJPPzFYVCM9d8+STdbcTA0nWzZSCqhIIULJFJbucMJTxLzOOCcyYhVZS6Z9Fmyspjcf957+ZsoctaSaYJ2YIeZMLtp66It3772aYLfb2MTpAfWkEkAiCtFDiczg/K6SPm+IuLKEhwlAWZj58ixqIEWes3Cyfuj3QTcPtxNpcoupzipuJUbX3ZQPwOWDGYFA535YB8seag/KvGwOoftQ7OY/leMrq0Qxv20gJHLodgkzjnOu5aArtcvS9PJJfy6UxjWt2EjpbIFk0760EEul2NyuMq6V58IX4ZbrSCnJ3GQHLh/BiKQuQcKO6YFh5uwm5yt6CID5GZWHyMsQK4M2OXHYSivQ0JtODjZxByoqM7ShyD1NW77gSkjz9F65Wl0dpRx9cQOXREn2b379pam+cRcgQc7YtPkXZOH7nPpYdWxEvW82oU36+MhBdABpvfOPDrCr8rD2qlys1SoBvq6cZakPc0kTO+gL9mNPCVbljfBcGqyKzfVJCWFd2txbNyEgWTZVZu9RlYAVSmLHPKVndbkaJbIKJTcOiEw1l5nVnijkruKRMpayoN996xF9Fetwet4YYR7DyLrdB2uP28aUq8gRYkkTG4Nbc9vlqBUqtg64vID4evW6fbgKKvY9C3XJpRxuoGUL9U7h8kMpDA8jfw+Bcem64m+FrylcoSdficGMlTpsYIB6QSwhIXuYoHZ7KM15rsYuW1JK+vZMQvY69329lqk9+h42xw2J/qUmF5TtQIgkzNEES0XM6Ps5w6/4P5eaEWmuzEMXriFYYYZrM7qS8VZh1QzU6qgVQsGXZUwJRke0exlgaNbVOtysBJbXPbG76hCG3bvisDKFvn2EvGl81vzMM6OuV9/uJ5T3WyknKO7Cfwz4ugHZXA7l0b1htRUFl4bFg7LZeunjyUqWhVO/pXDRj9ugKmpbFIGdNs1jwSgqdF/YP2Js9XNGmyKzU/bSNqKdDmE7TCiTZrumuKn9MV3pqXzFjqn9liK+zSv6fujyIOHbXJau/hpQ6l8Ddq9u/8C/p37vuXzKkdpHYHGu1jLcZOQykjf9rzXcY5BlaYWV0dSS06xgJLkW4mHqsM0D7zBpro4WFPnFsahvLY5uCuhVk1LR3TUXMpYxUyJsZjypGbL1F0jZqdn55Hb1bugCtHfFEfsej4mNnFs/0XOmbfodgoJVYXhT5OrSCejByqs1nRgIiwNz3DivAP4pcXLeYqJzDCYyo1Ih/GFmQzsaX2N+/jeIxrwG+5LhUEzHAXX/diUhfFeiMZTE0AKZJ8ErKU5SG1d8SnPByR5qPjW/LEoKdjJkFfhYTLUW/A5HhNjtqtI0syrcii4nUFKiZHV4ebAGl+LFhZsSGePVJLegLmgZ7U8kIPytkHstB/mnfoeFzUbKlb9maNCOKdoifKH3riLfwtHFIuwlxOB9YId0QsQdCLM0cEErAszbuKBnDpxaQKTknUYrkvRubDK1pZjewt+j5Gx3qQKpOSXJTLOYe/xVltZPmjDskyIEqH4PTFv4L6yqO1SabMmVfRW/tU+5fZ4WJ0AY2JoLlgEB8mo7MqmLrE/0lSaUCBLL5Zhsrtk2G6yaVcEqXIBJlxi4EOrGFmkrXQL0U80R9FvTgIy+P5DIjLoCs/zv0ospA03k3dLpKiu0/niwXntMME2n1EhwM3FOnZLYujGDp2SGai3OsZmg7zt0YNZtLhfCFJ+FZk+1RLRYy7LLsEplrXfzG7WGxzYqiPLO3DkwWN6OOJI5n8PK5w0WaU7vzJeDyOYiOatug0nho1yemfLJZ0ieIvRwBHC/XthaJcOy8P0kJqNY1KyZOgVmzUmNa31VJbrS5zBprCA3A7Lk8Hrhl2Rt7ZHCJz9lS12Lu4ztwZPc7dNNrN4P/0QkrswDGN5j444Bh55AQlUUY8OY/Eva0d+9iH4kk2E4PYOK1Gag8irsf8vC8u9ebC6bPSiC/ZVXmvHbnxwMB6IBN2cQ8mP+b2O3t/sQ/mZm+vc7Hmau10s850jv99qU0cn9IDWP11tYEOcv7iixIZhzrvJ/xCG7cvpFva1c7+fgQiOD8njZ5ABSxug4pl/KPk+nQFvMbgdHYqEr0UjF2r/cVFnHy+YqEUiEYqQOuo1u2ER8vlYuKEDkcTysWkXJu4nTWdQGndgeSETYOSAiuxnUmfepPqx/LHXM2DD0FuQCBlWUAq3FT6+DNpSrKSrGpNW2kUYrOdFKpG9n5L14UQuY5jKrCQoHkJlwSgtfKwnFjwZZ9F/Mz12mBuPaptoTCFWxlSzYfbXb3za06/aqys4ILI+0Gg1wWZmze/6sMjQSI+lxt2KjfsTLIMmnh+7C7iPJzzpMCXHm7ZfIdGMDD8hhNV5Gd02POeRoc99TU67NHbSH5/4/5NPY4Oe/Q5OuxFXkf6w7GLc/nkWTu73e7i4zoBVnMX/6OOPnB1bLj6pMRbeC3Lv+oTS4d3+/74SFEcHcuV+DGsszAcIb0ftQ6YfV1w0bBg34FuwSAUKepSbcHjRy9e7MEKk6vdKLSO9FQNncFmb7OpEjCz1govkpWjNceSVr+4uiHSTOEVIVeS0DhxcpzNCXIE70G3p1EUFvQ33gFSVBL0f4QHr9u/h/EJUePilAZKjQy2/Da0V506CXxWxdMuMdvqtkhYGWQbRnL28kIiTOWoHc/hq5yVYy6SO3rzqbq+KKJg5FsvdWYT5QaDGbbUfPOdEBiXvVpq2mw4zURuBWb7NAXUeA5Y9nivBKcqiUpDsNzch4MJ5x/zjOLpuJKR0jTHD2trXctyjMkXxlPElcHi9OAADwBIc22HZx0TdOaY8WzbSruXJ6MFqEg5G4x6s/lZDkcQVoHVHD4C7vCL5z/++slr2TN8hU0AD/BtTSgjVEqE3Yg5MLRacTlyJEPS7uTRUK+FoT6MBS2V81TFCy8qfWz2AeQjKXquhUG5RG4RpFOBdWuyvy/v4i7YWtH1MWT/gZFCQ3nLZW85fzo5d8L4htwgVUdY8dbOJu5HKho6khANoIiDx2udLzKDdOlk//wf5AKTEoRP9Uih6We+g7cS9bKmYbG0y2crW2hkiLVm5Ju4+gQJLPYjjLEiBOvS8TXK33H9z/uKIUGvKFGsJd6EOLS+AEJ9NB2x7RJheeiLE6gBWqFUV/MkFJceDXFkuPak4SOggRCdjU0NgPjygAfjuZdn2eUh3O5ucpfuwXnsDcRrzDv9yM2KC23NVdclPM3yjCqjYOrRWMoPulUBkGoJlhBnPinUl+uDAopp3iosP1R/2Xec0LXstn/yrYb4bW4WkZZb7mxcqmvuTRfeAdvabq7aUrqR6SMpfh8soP6JNEIsF1upUhls4+nmg1u3VtkQeKtRMDs4R/tdTsrG/W52t2ZkvEy6uzQsfhVMaS0K/gCxMXc4BzeyBybpBAIXwPASP8+lEdoTdx+9h4bnmlfZpc+R2V6ux2+3E/yJh2KptLe0V/Yzx3Er13k4JVyCD7H2NcqVjJ+KYUI6IJSdgvdT2RmSfp2p18G9yY9v+WOj7vLlt664g/+cw2C60g7SWRwelPmtje7tje6dje56737RafIzG8bPrDIhUaVYs0LIM9meDdtRfM8Q/eJSOhjPbJTc7SI9pbvJKdWWYR1vVnHFdmwMs1wMcoGMu0Ou/N6Ve/vyUxPqu9XbbK0xhp8NVzynDXBa5u5HupkrIc/dKYjGKejs8QPdj2EjGpsTM3gGx+zDGhFVOJufHr6P5R/xXloGBuShHJq+9TGSuBRTpsGKhx2/9IcuLDzPVLkUDwL7kRwBe/bFJ8Fq/JoDcbv1QLT4F95yi3lVIMBXn5NLvBhrS69jg5e/LHy5TISDhgEhRrwXSVgkptmoJNaCZqk6np9Wo26MgYDtI4CVozMFW6hHErGmR0vmrFimau6aCFvTcFMZKKzuwkN11TEDTAWjanLozb7XlGqPqZHkfYS0eSjHPqiq8mcEb8/d1NOSeXxRe/1WX8P50W9OV8m3lKddp9zfJluDA8CTC77d1b+/i9ZyVfeXme7hG6wc4vULSNuo8ZtUTE+SQkimH5+MsD3VQy9x0DmKHXR86CFcxqXg3hH9bNSVBoGGkCndF9LWcIttR5lEvKliOBmWS5pRgr3FsmgYu01ppyHIm7hDvsRX1JvMmMS3Mqa2iFoR3x59EYIsDwQWQ6QmbFv3414DCdC850J1afLqtHf0qXNM9FYCNIeUou6L5WQGoKaW6fwy5QnySl5Ek40+ozPadzfdNiTXizDpPoAVj/ssAy3GcOR49ElXNWPYz+fc0xpp/yvS37U/8D8DfpPaGLfv7S6wwXjx3lKlIfa3H4wZpcHnJ4q8Xva9Ew4+/92/dRzQEv/CzgRDEX9bOIt4QA/KqcG9MiMMDEzkp6rRSblw6wn23jtaZyclMTzp5gwLicosOxPcQQp1m5UnJ13xS7a2DNVKXJiZjU3jGYCIqGgy7nmnki74iAZohIIf9dHowjBjNLOBjL6SeBKxJZ5JOLX45LJiWma2EhAkqEB05g7KhQ9ykNYlqWXlgX0Oan23iFhsGA/ppPlDQnYS3OxlhbQ6YyQRfDULEAMGBehh/sYK2YTZj8J3of2Q4/DToxdPdnefADJPyL1ja6vj0knwnT9ZH67fv38Pqt4/uX13/daDMf+8t7m+Ud7nnw8GG/cejPjnwfDWvft3+efw/p17dw/452h478GdW/xzfPBgfKuUP4eDclPKjtfHDzY3OvuCcfFhspgMJxUbv317fePObWn8/ub9e9LirdHtB/cfyJ8bd+/fl8Y377oS8ufG+MHICtw5uHdP/rw7GozuSA0Hd4ajzQPtx+je5h1r/Lg8tHHjP2l64/76+NZQenF7fXzvrk7Bxvj+hkzBaHN8TyZmMLy9fk+ejjfurt8ZaiMP7t6VpsejQXlfnw7Gw4OxNT2ZjUeL2bzR+MH64LZUeLtcH+ic3ru/MbgrtZR3Ngd316Xx8e3bt+/qpN994Ialk36wrjN9b7hxyzc+Lm9b44PJBz/p6+ubm7d1Ijdv3bkndd86uP3g7oB/3rl3Z6iN33ODvSdre7+8f1cXqLzzYHBPuj+4dXDrrsz/aGMwuKNdGo3u37pvjcNuABcD1/at728/fiwNPuV//PP72+u33V7Z/5XKfxodcIiInvz96JyS+B5Bl5/Plnii19mG4+scL8HkCS2vbyGx6urXd8Da8/V+DNNpVMFgOpVqUTn9+tHLn6hUrSQxHdKBmQUtRsQUIdmIRC6xnyCf7rwfNMhrkSkKn1CD4wTxkk1uZ5t3iHbp2TgnOJxoMEPMwjkm8oQ3vx55uTSp6HPz6UarGkXHHv7WkZP7jxaL8iLP2cgNydJnmq8WF+Nt9qbmXExJA/liWYnjktiHOPVfl2yQFwBr7zcTHJ+SYaVjagvdy6PEPRW9GUhvBu6OvuX+QU/cmPYm6DngovZT0YKNOV5qbwDvcPuJsfLRWhYVwOjHXsDTq91VXmfJp+XicGS4tKIVTwMUKdZ51BnNjXNyMpppinr3d6dKgV2v0+nKMfYno/Louvhg/AT0T5oWlo79QtqpYUuSS3zwCH0wLU2ugApivECgkL15uog2TMSmrYWU",
        "a1zNwMiS02J9RQ3lth9/JXyO/TYtTrKDvDZgva1P6sGOhlIX9kjBy7fmt24/HFuPQdAMK4OxN1syFvVZbyQx00EIxyb9d+UbCccSZMuKkugSfh8WR6W0QviPOutxcCGcx0RyGzBvdH4dP6+nsLawYr3b2X374gnu6BAY1c86tQCYDj2E/HPHergnZK/8MxpVOin2sDA0RN0QphG+AVVI3SFZPSRBx2igEKJdvzsFPWrJlMqAITgJoBUBw1qpHaEvLJVUJkyP+09Co7znkjGq5A6xry2pgrrKzBfLroeTJDbNOYUQi5dKWSVxmOm5HSjt+NbtP2nCeK+5Y+MmkuB5RGYOnbHGShW2AaVY+dM9j3JRdhVQh5+JAwN5YmzrRssS7iUV5OA6Pf8NMzfuZQkdc9cJn1LgkKBnd6MQJNjSOLhdCnWOm+Cj1C8qafD663gTKil5NYMzGbb5dccIj07UoYybYu6EdLDnIYKui2QN9CDRR10LiFwW3SgifjpqtC7QQwpOQAaa+0XEACOci0llObrr8cdyWfKk8Jy5QwUy+M7JoAu3d0Z6V0Ll3jL0eJgnp8uqFseZbjg51BPo2CrgWTKZsXcZFFRLJ/fU59NSX4Tm9UToItv5qCSNg2aFhI8mDn7Xok1ngCUSsPdZJCMJje9GFlVZJezyY5iOQrNM0nGRBQECe6j9ONKFr0REGFwD+xYJHw7VDLA+YScT1QahF1UDdEzpB1/W3CRr8DISQLnTc5zvvEdxp6t0RhiQl44pcWXECon3Zlna0x8M6/ZnFXZFR1LVWoOLwu3bvIQpNfw8aEN+RCPN+M98B352oyc2Bciy0wg9tbcddMWxbAOoGEMlseOBRJF2Zf1pit3pyUYTU3MJ+yshBaKn3HqdoNGyl7SCYuiduuqUpz/uQ24NPjS/YBjx83ClFD5dBwegHSUeI+ELR2l1C/XJFdZRYe/BxBleY1SWjwBToM80yUfN2sZBJh/GSh7hnzmLr+dnUqbHg/yuClbSllnDikUjUogGf1FG+uoP6NmHrB8xNp/SLpIifVUPU6rk+1lvcW1lk0zjCw+hg5gTqqJu118BvQChqs3HgtkJThD7PWwI1fvbaSP4qD9u8VlrMpHpNuh6TyZ0u9ivG1ShXIkbDkcVM+uOpTul/KssaKyLFnOJtjvJNxyIHm82izPunx3os37rF+6t8OVraTWJBx88kYEZorxExGZoTldJxIx7+foscsi7bkSAAfSVYPwsJ5plqayQckOV8zP4BrhpMbQRzn0F76BobtxZFgrGvw4IOcLNLU/lT3mMgypP+Zc+xPnTp/xTHrtmtRHUpwRUnqSGB1H7y4UVYngZvSBLKQmSQ0Su7JZBckRmpApSy57800xwPpME5+6fkK0xIUczS19OBZpWo8l38EU/Lq63aD8bdJnFXHKX10Fe5d6Vqrp6XQrGoNxe69zV3OVFT2/ThgCg953PNSoXagSvrtlJvIjAavQ+RtBBP7qgGYeFrUbXtiy3+5cuG+H+LXDd6z1Pz6/RmfDm4hY/915ToKuMFJgANi9SITBxOBF1mL7D7cXF6Hjk800QjxuM7oSeI7ippYPiNurx9a5viLXpwQO3mCejBWJpJtMRMiqPld1RFMbs+g6gu9U6pYze9b5yWJW9I48lL3sZtlSmDp+U1iUx/Fm58MkK5jxfbkJDrmU1gvrMYzjHonVuVa6IljrRSHudjkyqgGB4zbb/qX0QFXTqft2W0ifTBF+iECfQPsHe58ZQCoM3MqVAtKCCQgpn6Is5pkPNg5gAN7gE3JRfKdAIxE5yd/a7LfmH7l86beeyO9twA3dj3k0PRTdk3mGys7u3qS7qYetrVb3xtFy+TG6UXdU5Kccld1lRFEL8yEbNac91zbH3vAUQqNRJ8OPIMSAQSFiNJIUILlfJ1RGpelAW+qe0JBRQzB+2HnvnlMtcuiHszXpvHRyOwHdH78jN9B48SLhK7kfvm+tDiJHyZ1uBTh1l+QDlY09SEClFupEHwPMP7emKinp2oQS6VPivnZ5tX86cqFfJ7ujwNNcV+pmvma9RP0pGj3FKniyWoQ9U3xApdDgnhDTVFKdTAYbyYISLl5OZIC75B+V5uDyUfyljhIEDn9PlwLVdFoQSxnBKn4z+INjVxIoSV5V8QLVULsnl8NS91TpQ5AZ7v7WiKulIaLKlLmlHxx7X5a8Wn5ouvTB8gJdpWYgHDfswPeuGAm8aLDhi0WGupiVAdy+C3kQFJpMtI+OXiN9BE5PHIlmhorA7l4eHrsHwnaWpj+4Z35hEDwZBemmwpbTAyg0mZ8LnGuqaj7ERppWE0bcSJi20kkipeUysiqz+n5R9tKz5iJdCC3wzAqiucSLuP0dtSXm1nI9+Ec8KekEq4eWbKBGpqLLi5pBafAnqTUvi53/4d5iFz//wP/vGRVpPZ1JAgkwBaAoh06pRq+RzaYCRljqoPLGtEa4/fBRUAAKH6sHGYi2019kk6pUyXD8R9u0ZGNkh7t98gpAK3gJqTQzccSmqlW5D7xIPVky3IT9w3HoROs48H+mX3nWc6YgZhEQrNjKiANVGfZUwAmzkLQYsSSSgLOxEGSwmillxD/5QLsu8ef3VVBe16A05N9s13UdXbbjb8r4nyynkeCWZDoC2Njc7+I3aTExWcdNCNR+XcDvrGebbu7mJxlH39+QA71u8k9z2PLGwK+SGK3hOuq1/lueUD6JSaYRgROPC2KVmGXDCTXS9MnxbjFN+/AGF4o1mMov8j4OfsrezkbNeB++30eaM0E89EWSbg30g8+cdpZF7V3k2TX5lPm8sq0mM8fv05ISOdq5pnjnJZTzz3nZkIFVVDkJS9Wn4m17EiXZCayr/sdY4AkFGum0b5GGdgRjERhLOvbeYYREcRxEhpg8sqep6747AEfP6siffuAI3Ze2Y7bLPp2sI0D5T9iAUwzMqbbxnY9/6FhttfG82uu61MAWKf2w8QpEm+PD0usEseXugAiNjalDkGxoRaYX0Q93TfbV3tO/3mFrw0t+b+/s1oFjLv+udz6GtWY6eTEf4lXckAA5BcSzpw0ZjDPfIacuX89GjcUESA18y7oVEk8qn9UDSLsEJ8ZpRV+wfo9JIqKJudZO2AYzthA1qAqP599ktxbgnfAyDRew+TGJHpMuXCvdiC9nOUnVQm2GwYakdi6V2XIemGsfZ/7JIzyaWu/EXOHgB3i8fuNHGvt43sjEm4Hb3i7aXbx/z38P/7M33s+2wo1w3olfcb/FrPqgX2awX2WwUubVPa/6dNJiUVspxhFZkDzxQEUp9cFMp6W24+IVfYslS6DndD44UcFfcsIK2+JrOUNv+5HUm3H0np8uw9dDhljDHShJBfMxMM97POqbFcJfbb2eV4w9huHHPT/2PTvapa2ib21kkRUIzuTNa5roT60r6pn5dFcA0cRY1XSb6theX3sf91pnNyVF1ggyq4afg3TopxSrlot/bixTJ3VirzAjlGa2unZqRTN5J12y0++J7ke8dSRKwfRFmTaV1ZEnPHY/9dVNwVNi2/+SGZV8NCFyoIJq9WQS6lTgqermir6yL6chUahCNcWNdvLYfZg65xroR8IhyC92EAe/6W6Ars9utM619H04AxgYKbmGE+lkD5tcyqtdkLmPvE/j1KXMBeVz1yYygIj7EVJIyRJIaAxkTozyq29EUCflQ+i1q21E5BdeWyxPFFtsjFN6wF9DGCDSGJxuxKT4oFB2v73M2UnXUD3KAMVNerQh/RzoqegNppAzSpDaSV08hmYPMoFw3l1qTVAAVQGSnmhE7mHehenLiThWjRptHKRsyNBBzPtfEFU8sq4WJKQH+RtMjRT5Odfb8N6OLfNimn7KcPLGOWoGC3cyAQEZ+xZSz3OiGPc5ToSneEm0yYb1sdwrmz16wmZGMUFO9v0+nKUhUfZJPSwKvGLo8c5qOXE5dESFJ0p24p8eqx5TjX9gZhZl76j61FjthmTqi8U6nQqAedPNx52Iuh70auvxX9iBgzQvr0a+BX4feacFOExLo8p55yO+v7NkbRSy/ol81IbRTNxSwzRpFacu7apTi6tSrWOdhlIBVMDecZA+NfC972sjMCqk/jyR7n5zVpBmItchzaqp+ZuoUI8Kk0lRRZ6rMhrJBDA3j+NTWUrcm8n+a+yboOjqaPKcbApNIErRmW3HttXcMqekcLMHrVtydXvai2ceu6f+zXNQfokwfTjCsQZTZVMsjiQ6yv5raLJBK1bZUV5I7M86oMG/jnCzilD0hjVAve23ZOybLumqoS+hmKnmqdFg1olz1BT56RFDwQLaDQkkIq6iydXCJhSBywwGPa5r6g9HybGTGBdwMcdJcR3QhijB5q0itPouduHqJdxb2nQjq4lweza1qLDsiGGNf9zSt7e+Vylaz2LYkrW2/HValBhvGScv/5ZKyfWFC2ExjpYURWx3R/Ienj8X3UmMjW6xh3K1AvOtQh9Wpg925PWzfNfJuh3ehFvNl8xVZIUkxC9REqDw+fhKJlDaddQ0/zNtJjTurQDJDUkTLAzxCfuzTE7EG/vWpI43ee0u2ULGCJSgZnc/sozey8ErhQg/lFy0+W603dahhw9eAEr4C/NDvP0W38Q50gYABqNI6toiQSKkdn8iLNvTGTRY19WIoq8l2h/KKooWHHAgInoc9DTqRTyUZgTt943IRUUiPscFWFhIc5AjkFB2kuQLk66gyVTXyM/sRyoaJR1eHctwyNMjmYK1A7ISQJhsul5pguOh6zq+FQ+yq38nXZPBlpuGvT9+Lhh6iuSuS9/orVcfqOR8VcoYhyweehp8Q/U14iJ5+u+0FiJb8wL7YlQmB087kl3OIoI5XcGqax9ZXDcMIO2EUoNHjr841HJShf1jSYd9VPHsQb9P/vCmJm0nlhj3V5UmUov857Kly75lXpTeIp2csEJciPhV2ITvedMz9DQ5xZLg5lgHUTDUtbEQyM8JYmIWbXcLGEJE/VoRGUwDgVHeM3SxjuGDVKEs3xgWlZVd5lUAkvgGR0ISv0j0VaZFZ3nxYGHCBMYDRIBMRJgDKGZ2EaCij83KQ6KRxNQ0nhxOmvtWB0nKPM7mhLgJe/X23GxEqwaWdzg831vP4y8Ls+o4oC5REri18m93FQfHhqTc+eFAYKVHI+n/H790cNKYWXwGDWPsMnymVZBpV9dQnJe+sQWH1+e//sVNPphZge/Z8UIStqghQbr1EBeAufHCmCmDpNixnUpezgIovqoHbw1cgZ4AV/LuWCrgNPu2r0GyIQzah6tBhSi5Ztg9Utoi5nxvLTeoHqlt0c1Ew0Oywgq4nmjf+sK69tXc4qt/qDnTUT1B/18J7e+eGgO4HdGapLV37TzVDSJrr3tPLkOQ+svvIuzcKebOSECknYEhaoggVprRnHKEne+QTjTALUdT71f3NC7FAvscaKbN8968vy2B7mCbj7GZ1VrIfuuE5SCy/v5Th9ZDJ9RrSlVM1dt7Xu7iWsDzoDv/A5O263yDS/jFyuKOGwAdzPEi4TnzsG7LWAVM/ZF9fkdDX1A1fmcz32ZNHu79vMt9L0twmEpgqEX1UUOQr4j0WVHqu6RToUdGa3PaZNHCZcHe1+GUc89DiDSpB5MtactfWwRuvlDHb0s52MQlgAnDTrbQGshQv+4G/7GG9wrdxN9qy0Q4NI3Io4BPDns9J+58/1ayhcpuaQ/MNRV4j1fFcIpX6alxXJ45KRGXVkwgODuHJ4Vmj+YpOytM4Tkgdj6DiuvAJBw+mJVJDiEuAJOMMOhDzIw66gfakfy8NFjoZb1jrFUj2l2bP7frEd5Y1F39PYF/bsTl5MoNxbJgaZektHe8e28S9AyxKZK+V3K7jlgxJIfeqmbRDKldkNRZxV95p1ktAqQip4o8bmb19UdTgUp0wcMzwxrGhqsbdJcDNrs+L8yxyzL5kBX5gINwiSY/AtRM1vrqDa9R3UTQWCudLfNZ79Az3tFsi7BYC4K8ZoO1h1K0afn89j2YtMewqRSd3s+TObDGbpOk0NXagllCza82DhEbJNP3OYB7K6m2UVjIGYI6z/SoAc124a80jC8tlXcZLCoqqnWDLkWXjWmoJa9rp21N3XpmwU9N0Js6i9c2zKm+wTtBVmTo/RTMaJPI01CpJ0RkJ7kXXlmxbUjRDc2DT/kV5rS/L94MKmyl/lvNmkp8YGPBfOgmoDDj3uT55xltSkiI7uBTd1dwZ1B6Yz1P0xavZCypS+unDH9y5+9fOEVrDaOnHcZcebcaifmlLMAtK3aNWgxoUz0FwcsELorzg1JzOYi9OCrfTyZGZCGZz4DuFOwxkQNgsb8lmZ17wT0o6uSY4MZc1KR6nB/e2EvUdIiGy1qfMJymxqhqaakgWpeCZaiLlkcXdWhDGDLUdwErd+1W8Db88/dKXY0tq6p92PW9yaDkQdy7n6b0ouBfIovcFTk9Q0M7HFmqkgXpAmXD/3wxtOgqYf5b0gF3EZ+4uxb+9I3Ha8L9HGgI1KKLfjie8Afw6+vGgKbODKiASw51GsImyYBL0FAlo+LIWaJYfsIW17KDHmnA95qU+K/VZ0QhRPKA6Jmf+8HMxN9TJnBgb1uEhhisY2f6kNie1bcRPkQnHByB4sv+sW5N9w+5mBWuu7kITPfqCMZqkCeqisMgvAN3G70QnkWC/NjmQ1ktE0B9lsKylPUuenWS33y511ymQi30hk/iCoGNKTmoGiW7DBwHThhyRX9x4ZMB4X1Y+Jlf8BjoqcLZ0aEwng9S8kVB3eCHUiWYwmvYV3IEWPp89mbEKECNhDxWUSk82xV2C+1w+/RENBjuo5jVezo2ixeo7frDd1IsF5Qp1hxu9Dbf8ELEKddUe1vxHEjNNdJf+MfPRNTLMtUPZxBw0ntzYlnGmiKo8TAMR7FLGQdYRk2Y5wi5qKeBCidWJ4C5i9FWqtwKgNXsIkBzbNCphAzSPkyvUyhcLXobNzG5+Dxm6aXk4asmNTV224NLrDa27wa5iERtL0Q1PlsnmE8iDsPXcXiLwRQT8Fjv+CfiI6Z4BrfZQF8WC5w3Ru5+tbYQ9pB/CG7t+E4UVYxmAXtJ3EiDI1BrQCPaFG065LJ2emBViVd+J2coY5cauHJp8FZoT0cl0jImcdNH69E1NplLdx1pNpkrWfRWjG5RimzIvLckQ+fqBTpsaTa4odSMutSqnbroFvfoKB7dvBoW2JAlfnSAzyiXPOE4oDeNN15YRE8n+MkEt2lVDCzwsofzx8musjvgXIlXJWn1Rwspf1W6o/8JyR/7/mRmbmRlF+g2m/Lp+K2RqtMsjzYj4L5Oi8fD3Tr2YHupgHTX+JIQM+uC+4HYkUUASwjNkfLkyIQgJBGNiOH9CMLAwLZqZNqBrvkBqx68393YFcNR9d267yYxDDW4iygT+rAbRN6BK3ByxQ0P5RHkJ1/wzhE7XutAGn744PMjhTTWAd/ONrNOVHxvxj03+KOr54TRz8aSbjiQOLQIBxJUyKVoZiD8OoYvuv/OgxwjdeVZsrbrCYvBibCRZW70GUmbtba6Wy5SaKoLsRUrvFupT0mDMLmpJvvOkb6vyL3xJdty2/lP3LPbDlObN1MN5RT+VsHWjQSfqmsNesD7Wzv5h5DN2HlQ0CaGKlXWAVGnJqArPnJBOFb/84zSR6morpqh+Ul22MrlkcEPQxNAHTIhKvl9T0SMqop/JzX2+0dfbu5l91Z0xbyrHwQuWd0ps+BnHUuC3aJWDzTBJXWo2tNXZSy34gLfqefJsjXmu2xKbWoHW3KaqCIbRHLDv5Ww0XfOtaNFaBkpNO4ez2wy5u7W5zr7FGSyBElqIJyc2ePx9MxRvM63gmRSIa4jmJsqDWU+BKa6LmFJLYrvKTuntrt0kyjiYY4NVwZy9m67WIY7DJ/USzvNLMt99QSY/7TcOvXltBmRwr57UVH8n09NFOXWnTZSJlo+pnvSvlz3WMBVeo0z8ExMUdyv0JEeSe+cEsGaB9VCgp6+/InNg/Qzo7K/ObuBXU7et3+dMa3B4dVqDsLNXZja4JHNBN7s868GCzmeR/1R+QcZrV3IQPCNUPZYkYuAWGtHrBoA/mYSDuhf3wNIjrAWXUD6NEiWsHPJZOZOEYpoPYWUocJomoSjqiRIEzVdTVxCxLJe66SPVEp3WN+XCnhvHflek/D6sjaLU7mME6NcLJl7oW1cxVqG/e+N9cCcIbpdX3yWvHK/yxwXiJ1Klmg3E6UKe4CSfzRe4EAzfaDxf+PR+CKD6VQI8fXosiJFlXI0ayQXpVAE5iAtK5TxC5XLP37oKRBi9Dt9NAjtJTVr5dXXdgPu3On2o0r+0zw0ew3uLwAQyE/Cig4tsOKInK44bwExqKBOMvIrMKTHMq+Ov1/j7eORo1MDCCjhTUThDABhZzBXuwdtMQ8gZJ8RQ2tRYQqgyMcAQngMzgIlSp8YoRqGRFAWdeCxTlAtMp8S/cmIFaTUhJu95oQbNK2N835Mgh3QZfBufKffQiYkHdILksmiuiDSfmhw/bBvRULtOvcEv7Rg/LCJsJrYdoP0WkB0QpTg/22Pt+0XDdY7wkh1xmYvS++iLHx91ApdpToIjt2TZdpaCK6nvgD7kZFkjPkOW5t79UNQxAwE0fFzSL5ij7fkHIc1X7UU/++HJ4+cvH73Y2ZM3p66J/ah86+t+awd9/t3njkQcuh5GHVWfTx8l7ojEdc9Pz+aTanRdfY3oQY8HpnfsVNmP+P2uchfjYCQ5KkJqBz3S7xnZlOXq//J08l9PzImAW3ey9LiBbfzqx09dIywfZQt1bY26YV6tUE86CH80/vVc8pOke1aUnJ247+7Gh5faR7MxJX6jiqwQ1YHdZzvc7U9c1LxAnj56+fzFWyZp3a+l3YSVgqEotuGoFlc6uB2dPt7PO89evd7Veh4S1VWz2gBSMryrRYUzyQyapWXl8asXv335I9A3tH8Pk+dqIsn6ChrrBAz12Rr0hJCiJjHs2Kn6fu7oVDmLc6bHZCMMR/8K0doCRDvoldNJKSwIJkR/x1ndta+SL1Aq8VQm6n7hW+tVJ5ya0ORsyGVvG46uQHZDMyTHNdqS2xj2er1e+zhaauxmV5a+FkrvJ1yyvwfg204brRywaXC/8dfCYnQyIlxQZIjiLer42zU9VotTd216/EXFDvT5ZSdLP9eBnHoit/EFq7jXeToZTYkW8Np1Ev+aT0RnP6yWjpshFhZkfwkBH/SUhBcMSWJfEpZJe1NHrUwvZM7FAgKlYdowgFVIURopebJwtMNNoc1tp/JV0HKi93kZbnPWUr9QiZVLshBdp/VUFeXytKrREE+HdvgWW4iNrCgVkw6sXZPUWMwvYo+lQTHV6PQlky0F9i+rzlMu9ooYPMwerGJEs16dT4I4cUGRYCtdwpwFwR6OccrjpXyk3Jtmmqj8oknaZUSgflRjqipnGUkMJhMeAPRpdPc8/q4EuZf3RCb6prVySr8cR+Z7eh06pqkUrzVlz/xl4tHSZcepHtQbPJCWILv++e//cb23vr6uqPN/dl3x80sDDy0PR1vZ7uvfPmGBp+6mfkJ2OLvQfBmzeYTrpX3w6F4QV+FdJGBMCocfwnolfUEKysHwaH6E+LUojjme2Il4zrgC0+XkpL6XyXJjp+U2dJogxcj/Guv8MZgDWaInOBxe/Sp3pf9ANgegSTpbQfRRVAGUxavPf/dvO+ocPxDMiyReblGeiXlxT5sUOhEbC2Ak2FkinDdvlIvOjD5vXge1ngVTsGu7iPq4FDd/91T7Wyv5KRJW0WvhPSGv2i8wnPXm4glAtIDbtn3dez1jRz81RuHpnkdz9/PNmybk0VjVoI6ms4P56FMngX7eCEme85aF9BgzX91jfZMmvlNWFNPX2CF6Iu17Jx3LNn2Ybayvu+MC1PEPXWvYzq9/AJaYs/7T492iZfE+aMd5tD/FGBPW4LQ8rBqz5yogz7HMb+5t/W13/6aG3qC0ZFJ2f/QAOJw3OScWRYEd7BwPz4P5vnSjYWGiGb35Vzl6/bcS1fmnNye9JWJkMIuN/t78KxRNyuDMODLE6KTZfEV7veX8BcDeHjtpIC/a1pO8ttt+fnengoPcHlFZsmXmQKB7gOq05HiHol9wxqPCiRqe1QbZL/xsO4Dq2HrpcYilBXkpX3G2WH3bDlLJNcyDYvm6ivcb3WgWWdUpo5aND2yVmlu2vo/CvSu1KO8JzlKUrhGSAw9SziRNXa9Sva75UYjF4X0FmG+Qgl3fsg+6e6nrY2WE9qKlBIBIuB2YSkLoyKAHgRKnfUAlmVt//cMRph8fvXzyw7vf/vh8d6cnkgQ6yY3d5SzJAwm/26q3o8HaLQ2JuZMXlvt2LzQ46PmEPlEb/hnbcZSnR4ASTUXHy9H2riwhOBv5VGaCNUeErdMVZx/+DtGDJilJiVhfJZXJS2zVRG+lFew3qJCQnV4mmu5eZ8vvBmG6LP9PqjKP89b4VDOCDkS0JbflsEnL9CtjJbzarVZfQ0PFakye66qXb4MHCQuaSCrtYmwQUApVtn3BNxRuEq2Wutaa8hZXH7eOwSVpOtT17NO+RXYmwo47hBabKTWw+w9juoYnxuHYfcAdwF6D2LpC+DuwQeBj5rOle4jcEh1aPfGrsP3j/vfPQSNWdPdIXWyLokY8dIb64gEkRsMVMi6JK/qul6E83jMCFM63e+OpWnQa8cGKSRPAu8CZNphSaLESmiaRhaXJa48SFXMepVjAFi/oLFfVqJOw6EJZtwJRNJ86CsjEoqezpSm/lYuPkHVAhKpujIfDz0xHLmnTtJxQyZiUljpFEuSoIeEiEJG4AsyYfe7ZybUezBdHkmJArGkHCwY4xr0XWdyfMt0YJvlKN1skrk/oy8m0nAgoTiECCjXrWrV2Qep3N+PfcMOhu0g4W0hQG1HqmCGppi508zxH89Crn7mCGJkAdcxnYT5EzHOzh1ljYJwIfa3a7nzildxY+5SAjAQhlTYv+ubpA3e8ZqPFs92XgOuA9EJtFxFMob2KZPSayTeyon3MOuzAGrUOHXcM3b/8sabzjotEioilwBeSn81iQmx9VfjlC2Wf9ibD/ViTT8N1js44wvGn+KPoDd6PBm5vNMAgeZbcB8xZ8nCVWqNfl+GBpEhkGZ/rpKN9JLxHQ8MEmIG2F99uZ7fbFKPJNsXUf4kBozBCwo+ZF8R9ubcv4amJJUFhF1DGEzckVIuSXQpVcKSi8uwdxlyAjCe3FO3qV11faTQuGQ/vK/WBeKgXJyMn6X9Q/E9qwzuaY8no6QfPGbtxfTBB86FcF4FR/mBC2oc6J/2hBh+s4I9Nwpt/VHILohDIf0TkA/vEf/NEmZf1LyfduHbijgi9sAXSX26MG9IDWXfbbVASc/c9vExnFrkxJzyDHXRBIt595p7mSCBWjRZLSLx0HGH+i1rJ791TeEhho4IQUfW/toFUggiLBd3ZNnMvrqOnMgp/bfiNJVSTOHZRFry5RPaO5/L7GHEwStdsBCcMMMo1MkXlHbTxLXw/LDSHeTSgBoa3XN2H7kh86I7cJ/otjBRgPrjSBAU+Es+6yJNqSZ9Cmw/a5wXUd+E/K7Loh7quaC/4j/m2yI/ttM1VvgbiPEPoqWS3TKhuyqeTaqmTgSWr0+5ILryGoj7c2EM/YwX1S7/8+jd3LgLJ6Fa4I8lC4ogMOxIeORsl3WVU/UhzTNahM9Ga5HjkW5xhQmPz3un8OBfqcYYdoEgGjB6hXz9YGK88C56C0aLEJxT38ziT+QjA1i2jkwstOlMmlrvVw2GqRsue6qNNmSsF3b7AAmOALzCXy/nhobtjO07a6NaLkjBpNEHouuKzhsA0/DTWOFmS5aK5COnk4dOeQOzZeESGK1hWBb3wOCrm2KmisVbuYcdi3exkhFIccDkc5tSO68Xbkdg0Pzx3OPM2MyLZFOxW+vUhgAqIXl5FKL88Uc/1GvCFhPvjrXrNW5xRVQ2I/vhE8xz6kDg2jBQ1SKM21fRrD7PcP7dchvCTjQpsGKJ+oGvf+H6VjV7pi4P0xUNXcs099HdSWfSmc+QZeCxdzfX5gSKRKVfx0ay8pp8rUgNdtHsaWwdedquA8JfvOyaAZNs1CQRRDLX9lQguWwpgLtvrvd9e0c5Kds37K3YNSR+zZQ8fv59Mh/nyfarEmgq4X7OiMaCPOkWTLAA0aTh8Ar8rlB45OuiaBf/sTqddGPE9ZRcYYzYfZmvhHdIkhfvtKPGo1QyuY6Xe6PwCQ+bzqLOSZC3vnJ4QRdwxAh2Z5XQ8oVFJ7IDykFr1A083GNZCfsLs4fXwT93/eRmuFALsHCS/i8LOrFZO6njJxKmTXZItwPu3LZHlbSlBOBW22KJjW2UhCbwTilogBMsVJU2rvwtdCmNNUHN2AyS2eUypBPiG4FuL0YfJ6EwxHdzemiJnOfI/msuVhk4LjrBocvAFHKcI5+qdk1gPcn7OxdFK5TqKsnNkL7GQtqYdVa3SVPKgJyPLY8OrMFivZRstoAPa9Yy22EPZkdfSOiTGNlMVTI31zLyaYyJK8Eavq54afnVY2/oHudCKYkLs1CAvefQoQCZinMkbHfnZOITmc54acs3N2m4Bb4sJlh+7LlrbbHrFSC/Ml0q0M037umUatn177MfOVPGWy/VUEr7yg5kXn47dWfQ6nOOe7RTtcXgggGihGFftRpYf97BTcKBFda2/VVnqK/mS4tkneDftcfp0Fhl5WRi0lUhwQTf3Ot7WKzRfs6Bgt3Ol7LmMTLZllGgaBSZUXmDv9Toy8KRFPTKrWhTXE2v1zfsLj8BAbc7gPWJ1hsFJEDPQweDdSMOWu5wglNnjnT/ve9Jgevp28iDaKputdwPi3cjPd7o0chT5hAuiB580pM0douKhqqWdErmbHeJrG0v3Dzmy6XHVzUyx/sqjmOg696IDUDH9QjgEe51kcuj2Ek9P9ICLtV/sJ+qHP/JxjI+hHZ6uHRUaLnhQQFfk/+WY+D0TBP7HktnXBMmU7IgOp9uyKnzT9WoeV/Efy7eYl8l8MKoqd2yFvAmIuWSZ9Z3/QTbMR4aellMwCuNxx2fKRX7oEtns7qwLSQ03wy4z+/gho6KeVCLWQdSj3nq7nI6oBFhiPxt70Yt95nlOE4puRadT0gn60+nOppJKSbnGyFnYWtQ8guA8DtvA5CvJ4/X8h/pZw4geVx8uz02/tCGH8ccJ3eKk9SEli2wITYO0n0gYepD39jr0cMcWQ6zmO/h+6T5lZWZk6np3t7b4Pck75+EAGMvHA8MQ+72JkEsp5U5I2kSMBamgw77PrjAWLP4tF0nUr9TewWbl3V/OxDDm/o1J7dPJeHkhTvqVqOcZTh2sA3lQF44n6vtKFU5J1i4o1kmJhtl18IDZxvXIWef5Dz5vi/roKN/1VBAAcC0QgtDxRF+hRG8yTwzcuWTnGG1N1duxDHAirvFjRItgc63JcegUPXlVFVqmN5xUBv+mT8QbKTlz1/yhi8tEbFmjGYNQeCXNIcTGWgLvFp9ubWer2VfpiX8f83B/jKMDPKslO5BudTbDHooWcLs97CSUSNCGc6kUGX+zm1INSB4iUsKPwrDJFNbHV8qPu3HrN6LP2g77/lZ6esOncozRwr/O6U0DzBKDHfWlSNStYc9jASMWIln05QR1M6FWPElin4aW2NPYaPDiP5RuI2zYRXlGjU3jxWC+WHBDQu3x69eP3sKL/+nTVD/j3Ys/Su/dsIyQKh9pP9G+CAjm2kTTMjhO+1Ipr/8Uv/2XO08ex5/6gAH3fL8JQdxYk9w3M7P6I4yUIu6fuUH01Y7Z0qjEJ8DIrLaDOCAkFW4GTVd/46JDQ5HvhXhreAVRwsdg9TuoT3i/FQYMCzSpGW3MeLnphXKjHoyod8Sjpsl9KlcDXeeS08WkCmKcuSHHEb4aGoIoZzlyu4sDc/mcF0cn7QLAS10PItIaU5Nt3rBW1p1dpJpJi7Mb321L8yl3iI9+Q3ueMkodz5e9w7t3YYSNgxFzT3xCpxz3x7vDRXnxTkzq0O60fd7gv+whKxk+fbr+bj5+h8ouqSc6hklV4TlrY3emow+jafWOscbvAh6saJ9ccyJrJRenq01DLj7Wr0S354jiabrHvz4dLS4EyWa+eDSd5h1w0Z1CNQurtOh4ZVohVOkDMbb4qiX21W7CSz7GkHGIODT6DXa80joaAsOarxoEC2EY/CNEScvPicY6SIS0behwj/swlRBjmsZ0b9WKicDTwaz9NJ0v3SZeLC8yS/GIT/EY3/FvCSYmSkatIlKDImJ/EueeWaw2TkbtqAlNgctFbz4z5d6cJjNWfi3qZfowjNAxNMKWPZ8t53/uSjr6ejCdD47ggykho53WaPf52Ww6L4e5K3wghDiYdqYTKpNWKbtL9BVleu8XI8Sm/vb1Cy0jJNf9Zr1WzBrLZM/oU6pBMZ9uT+GGAWq00AXU58Y9P4rq8+05grq5vm7ZdsW9jHcyp4ceIYwAC36CppD0cZGpl0uUwsqS8R2XS8d2VabiaMY7js5P7M4SZwFHsUezyr2qQQfAMS6L8xmpXYJB27xiQo6Bm3t/9WjtL8q1v1lfe9B7t7Z/46bra+ddXOKv3t3423c3/pQvwKpSXe04u/tpEpUZfS4QHvMDkGUKUdZ63G/tw6xwR264A7TBfNNVuJ66pTlZ4BguGK4yGFCfOk7mrdtROS4dV12uz1+6i+o9H27U3kjbvKN4XqNXzxD0Lu/iqiYzR3yaz3fc0s2GeJ6qVzi7uTLMfn0Lccjp0kUZnN5NcVEWLx1Gm3CxH2YROpuBkwXTIrxVtH7DLZPfjY84Ql5eyQdhdTZuFYVi9nNSV3lOvjPHSWZ20A0Vi4pvLLslAdZAcERterK0mOEyEhR72X91UjrewRHi6dSxHqpJPNN8n2cap2y2YZaJtfRU58HxUhDwVWKEo5YmmrfcczVVfTeCf1X9hKodxf8oq06P+cFwvqUgDqL8aVPVxIHP2gOtRUPEiWQ3U88kiWURR0VzLB3Ooe5sP8EvR8sS7eTJFNkZ9i6L1cCJ22UfSeiBvj4oeyi2JnVgVfn+3QdHM1zdfYonth37WX2DdgEPcTiZhfJaQh7/uTw1DTTUqZPR8OkEMLNuqmfLWC/UfAuqsgMM0eWo5W3X6vtpMV/OHWfaUpm9kprs5Yk+9ZrxcjYZj6qkN/asC5oOd8Pqp9FMJsExTPo3yJFjGB1Fu3CMnywSKjnTbJSCXCiHSCFj+k7uefXbn5yQIcAx+9YJhTbY6Y1FCltScunLc1wqjmjs6VsiznYzrBy+wL9YJvb8nVt/MvJhLPKERnOdTqEdMrPyVqIQkdHCDKVon2x8X5gKIp88zPb8D2LaRL82qFflA//mfGNfh/I38/kx5BJ971HLoYt0UvOSY9c/u5k27MaAP7ruqOqY+If18cTVxc/cvy/nw1HXKzwDmocb6aT6YTzOlUO+6XhW6iOlaDe7eCdwfysGOZ3Hg3xPOVv6LR++g1nfxnWxt+N1njaf7OXSF9EWQpX6GnPn5sqRmndC698RAQY9pS0Q9MONCWlDo1fhmfq3vIPb0DsmTolya+PhK8eUFN4adPPmSoNHTPRW0MVIlxCbg7pZM4gd1k6SC+FlLHNchnjgXGyfOe+ylLQjQp1dgCdn1ETH7B3vcC/Eu1wtUnjsy1Tvy807d1tKyQuLaa9nkRbA6yhliqYeFNfiblAruuk5a1xkVV+znqg68q3iy4Wro5wx5H+4CFOmF4JlaGFkruVmyZ9ODnua1/nYSenQ958YE194JarcQIStq+hsMjRN6Wg4Wc4FVkh05pZ4VZqrXyjyFKghVd50CEavEcpvfTLPXzx2ex7/MlsyT7h8brdP4u7DHhLFkJVpwP6545w8iIOrLn5TntubhHkKjrU4Sh1unZh67bFLhM7kH1eSJHrRWu8CgZB6QA/4hyMF0dnmIzvK5mzS/9qehDlJduQzYkFgO+rC1RBMns4Xjs/nFTCfnnINsdl2nz99SocP97ziw8VIwOjlSpNAEzoaJLV3Gxuy8pvXd0Apg6rsgVu/xWTHhhcjaF7kYmpMlnySGzOmABjc8oVitrz3GcNL+urF6MCCl/L9/HxUaYQVDh2Ch1TNrWlDhTzha+LjnAUucSRO/0mAFb6nyWzMuXRU52R22JH8ze7H/fO7UDOeTPrZrfX1bjyr7u30b846sSVtMnNcVLCioW6J2ecuFVsI9ik3Wz+817DJznlHlIAS+xmFZ/mRp13F4zX5FbrMh/jbOi5P3J+dev/5Inqkg0m0OtEE5y2aHbjhufVVRaX7VUXLVGR/Gn7tuZf7waygEhyqlzdei1fvlQKYRaBznDgZt2gTlhPYKGLJj7OuC5IXPiuvsKluCaJEScq6ysb9xnbtDlzxcko962nXmNW4aig52YJjNX7yRxW5k4Pik++htDw5599UFbKf7hLky2/jXuENEKy+kdG65SskiOn2f/qnujaAtwvSgU8IygLv4x7VWT0JcgrGuj+CNVopQ5/YSnZdggJCrrKMZUBbrOJ7FSCGiLRZziV8zqdych+RNImesgxZGwflYjFBdJLidIhl2JBXZdBvo7TcnlTxkWXFi8/60+e//u3rJ+9evfkRZoYOv4KZVsBm8RduIfm3POe/O8sRPfQu9MWFvriwFwflQl/hL/dyf6ve3stXO0hven99dFeDKmQH+YTxEcPRt/TedHIuswMo3N3FfVxecFJw82ZEKvEwH3PkuwwMBFeCKeLyKIUcXmgOeJzq0+mQKA89Agsqatn7chGZiA6ZrshM7d2In1AKTbYuwvmiPoj7eWmsER9NyPxlUzgDevR/Y0CYz0T4Jjl/Fg+roI6Tmf5lF0+dXxIXZ+tmaNU2gyQCFkgLetbP5jJFbltVJ6MSOeERSqUxZOz6GednOPc4aJj5pbtcfAiXruzrV69eFh6JXm41fnM4Z1eYQF18+c5glZaa1C0x2Zi74rbBxVjv3e7axNwWimyLbaFyAnMxGmDuTRnCexIjw72KXRFnU/DcqrKEql8UXeGZ3w+uf7NR25nBSGPnvN3aik4qG7/ePDxcu19YmEuvBZlYLPLzt13V6epOuKpcOmPQGnpepuLw3WHLJREd9IZ8zBhHRGTgprAAyGS+RADwcd3cWmW19Lw65QKNVjHPhFG5mE4gvvITtuauUHcYjgLmtLvkxu6T2QAIYT4AWKeHN4g/MVPPwDkuCpdCEfL0NE7noJxhkVXYK2cXAt+jTnqnlVvlkOdIvvQZhBGlh+6E7SMCxmzEGHacGvrL5iowlBdFumcmzPiuYnrTQZOpTcg0Fx4udXU5LeCKzk9gdBTh+q9PR6fgnGmw6WYqGMs7P+/WBVTm28EPHEH8BaOxTMPHT+m+4RSk1FCniWKg6iX2KnHP7kq2UjIqXXLzXfpdIlZoxIPPjMTLxYVVyd1mgbi77oOqryuna5CLxThMP9M9eVlSb3tQpQpqkuw6TXXX+xGoFe0Ap8zcnB3SuEEC5yQ39IYexwiBnBjPfR0OndejiKvhaFw66ZmEC9JXVoEp6ocjgM8l7BbbdMMGx5F2CUZCuFo3X8g07drz4umfFZJ/1yhqmGYVgnmgT2fLyZRncaKiLYzR7+eTgWZjsFxEGiJMMUFEh6DrVRB39JVnZQo4vjTvSHXqdnqZJF3t+BwhqMMDX1lSUcYFykSBGAlXI8qBIMSEyEu54nQvoJXJ0OtYJHWrxab3suugy9ex4OwnQi7ZeYF0Oh6JntiT++BlNZyP9IZ579EDEiq+8/wvnoiL2o5Q1A4pL/54Y39ofApcA9SPgeyN3kcbXW69fraJWP5uFmoS4iwww/bXV9VFHhsQasKKRUILIHtUvsUvThjbTsp6Oajz1P+1omgs8sBNryYvhy/3t35VM1vE92fE/AhECFFbhSL0W4g24+LJGfgbsAMZGltHCHrroj15vPv81Y87JKlPqfzqrBDgE2FUP9998vOurDluYwxt7HaaTJP8K+eJq6T+af3s6asfd3cUvmeuER74I1iENhGj6lfN6l7qv7vRtd2yDU6WYRvctm1w+75uKV+Z8RDiSwke/yurvLfZ0sfJ7Mh3kaybo/ScBx73TtKFyPOhCxif+FfbvLlpZh5aHin+weqYRgPxNzLju9HPfbgcRUIDqMwPxAD3agSv564B0F27pnY+/k2SZM54e/LGNX5A4LnVlWkYBYMTy6l4wFzRrurVG1vt6fMnL37gPv1V0NJHpMc7jteW2C+v30HJGpz7dz9H4l572SDKvW0t6wpH58J9qRq7vc4j3QUmAv4MULt0q8nuDAVFFHQF5Y+VBU1W/FkTCPHk1naxovsHShj3yaTPt1f1ycTTt1f1ycuvb3/fPo3Nydl7O684D+Da8PoFpC+cN8nlI+uGhAM4WO4fHoS4j4eLCc/Zr92/ktdZdszgqJOU86k98OMRF92x8EeufHU1qVjvbd4xanE3GegLaxF6e1//S3hDoDNfV/WGXZe+Vk077arz1SoBCqQIFYMZ6DdP5LVweNt7XY2O8c/Ok5eXrM1BKaRs5z3drvmT3YQciBc/Kjcx0+VptgOq4GeHbtu/z+x0dZyB+rVOF4ELrJ1LpiupptlnGtp8p0kzv7jXnmu5262tjVzIrf1mE6GtL1hoX9mK3r8ShjrUOQ8PmkP4sxauawP6bGXo+fcXd8BDXPkpfO1Br/6FptE3GY07NPqvOHZZsIAh2Y0BJa9eUESaCxJs0b5T3IHDjx9SCNPWmvVmj7v5EtZcFBa1l5Bd+2sFDUDmMxDOkKsc1jT2bWGM/Gsq40xg8e9XUYm4Jz76tRvAxvwuCcS8Ud5u9PCNpTHRW785aMtvLqxa+LulFbVR+9XzCdyj8s0GFqez+BNfNiU8ykbV5kEcauXzvz6dwD97MfnKOoLVn+c++UWBZUUl+2r+9zqQ/7KYspThXvjSj0NmnP9Ps3GaIUdWiY4uYXRh7dhY9FBc6bFf7I/WYsEAQYsVNYttw6h9ImOVT5bzqz9I16PtlCD6FM+pGmPl4lkWdduKLOZna1Sk+Lvc/2oWfu9k7LXROcHTSceeQei2BxoCEL5sntpj++MNT3iV5BNZSYUWk5OEPJigpKX3kZyhXWOIjLV0HKA/xtncXDM8Zj1uK/H0E38RbKZ8GVJB8WKiGn82t9KtqoFHPz8RtcA5AlvbjVduXBd4227B2jcFGoo0bVlhhG/UylAqXjYMUgcjcdyUlFRQUQlsH/5ULZ4MvW/4TdnidCpWADPKi0LQFN18OpjDZArJu1VN3u618ufldDJUT2bqYiVILvFgGQOdH14nsbiqXiYCdikJ48dID0Tgj9FFIercbFu/hsMS/ti7tS8qXojLtVe39+lV9zFCBpY6IvyXQvGp1BxeSx8vqm/pGq3jvoR63WhcmObZ8+HRUIZbcNsHghud/xZAizl6anNisCGGPY93PVEvJ+Bp321n4RVdtaQodc9JyW+1JN5EBSVRXVpnFt4UMS5wmCBSeD+kdH5iAKDwhqd1dtjxeWgi7DzvVl5rBIfYN9Ks8UDcmhXyKHynTNglX0Z9uflXf7K3vvagXBvvf7z7yQNNx6sQVaMNCJOGzzlTyqn1oPjNNQCUm1T+hIcAvgt1hpACNQI0DQ5m8Qrae8PY9FiaXk8v51pV43EdBg6qa02alRIp23w4IfH+Q+/DN9nD7Iaegm+yuKrl/KfFaDCBQjbf2CwshU/kFIPefEkLPJjReQFmIBu8GZeK6jeS7hhnzcUVtLXnAYSu1Fyxnoypre449b/z3pFMHii5RuQqODaHhmPN6oVF6JqJex5suubtQJ6gZ6nC6nfFbNT3fQodAtnVq4bjAfCQ9Ecs07lyhWQzinbq+uh0OU9cAgnp5P3vQjCP9zc/70t4P5zaAki6BufAAip5WV+N6eSotkCmx0S4lPsprs+a5DeuloV6MFXxZkPZF7xlcvz5E0y6/Eyro9ud/rVRaB5Qfl73/rwIRvxBOYP0dQo/V5/1VxdV0nJFtoCTU0KBGyaB2/bMuGHWdVdXT7ZWbJvrCj0be4NtbJLrZj+WP8oGgJI6uCHox73sKaMZ2+5as8ryQhQUWkUlgaNIvKncZjFXDnhtKJ1gi/NFzfQZ+UXQLoWBhQTpcnnkT5/vGvY5rJbTUUVPBc6C7NCzuD+nmgodVjrPNZn9vTopZzadoyYvRbOm9KTM1Fc8NVcrEhc2/6LrbZren0Mb5jGJrKwB0sadlS3vjSI8mjnhOL6HmFZiTU2PNh10cLTllPrnh3NjjgBQbq27AXEa5IoZvpaBIJ04/UFdo0U/IspElC8ZWqLOJdASZPCscq9h/kPa2EF5WgVJPlqmJAjmfQyqo54pknt0NMHUBaQd3ebiKT5phKFIra9H49OqnF7BgSnD9FUMmGOnhJsCjslh7DTclbSoHiVWXQpvfEjujrtFiD7rrIG7/fz3/xhnDLp2KQfpaVrEJoVkc9HDKMFpYJPw0D9R4P2INwoVqf7Fl5X3fdGt2A7oRHmJgl+OMIjG+1AzPMYUwanSB0Ty9G/rYj2dLHWgHzOW07sv++Smao7UqB76w331LfyTdJDEWeeztCb5zo2j0S+6HyboEVQBbNc5O+8VG7uSRnJOoTvEJ3WL3u3hYZRhS9AfZXHRml/DTuIkunc8oWf9uRy5/bAzWSm+RGQGrortxkUorwmfvPZ8Rv7C7Rn7C0BJjGDMFdj+yPazZz76skrEXmacwHYWiEAu95urI3edLPSv8jzOCGcImCjgF4w1KZaxMj1oxn1K4IiOE/Cc4Hvs/UKYXBpzgKf0GX2ZvMGGxCljh9yxbGvfdevK9iczbV9a8u3LPj+WXiXtl/5AhPbXU8idPTjlIzRnWxt2qyETp+KMY24m2ZojkEAD2X3++Dc79MZE9NbhxEN/nLi53VhPcET453R+uLGei/jDTLWF+m/m4dIzYZ5RiBPxVg0uZyHp6+lJETP8sgW+1Upw/Jf6tx8wuzAYTabSAzcA6bTrx2jtARBN+Lvm4KvBxMqtBteDSkR1y8MaxOm+ZVRl4iwkU62C2K0OQJXckHBFObvMmTdh/6UTv0Gr25kyjxHA15G8cELyFe7kkShcCKjyteDIm555kCd/Owj5ATWsC9MFG1e5uvaBJ30oshVpeBhbwE77SYWKop8dzcB8ifZHcyNMLEMC55ZeevCBcsyTG+mwPC4P4bPh18DYzYMF/FGNCV11z1buMraYbYnJlhuW4NXsNNzMPtoE/5udVz/2TspFNZIIbgT/D8qlY5Xyg8X8COCkHwNQo8wvP9VgGELtnqifvYKur4Q++sOXL6wMGvS4t81XDQ0JpbzWi9x/6DoAaUP9AxP1Sn3fbIevUinB1GDmwboE1LZjuCr6/kh8jWjFcIjA1Utr7qWGXUQC2RyhXUFGa7oM1V1QkxOmNT+f/bYaJWcsuULrAy5Ul2Boy/WBb2+vmCQotdSVJvhaqYdW7D21X2vAh3jU6ya8jF883idi1o4jUhDBG4V7bPtwD7kKY4LnUXDoNMf5p2OvTqamzGas0sx88F2VD+5ajAHCRSQniQiQwvrHsQXqlp74ogPU1ccliEdmlkMa90SzlHcrhGtJiO7+juArEgkZrwI0s/sVZ2OKhOOzfhxAgnIW0/LgLjxtG691NPrehGaJrpif9LXJh9ldVwRveyfn7hra2ARWD37y/U/n7u1G705Gs8dWGk4Jz/YJvPkd7SsPnTCEyGo5GeGWGnsENxHTMO0ndM0+Q2hNflhIsMkxvDTrcukoinjg3JtHKeMZMOzoS++4KiNTF8scGCZFlIsL8Ag4sTA0NCIQRELOwe8WtQtv+WpMdbBEch96tCj0BCnt37iOHfaeuTmUMr0z90D+eq9ZFVYWi5H19RPHCshfbrEKpvCYH/cDmpmrqIt6wPtIOVvDT0Ws5RojasgmNIoFoFCMoxzIkQqMsr37okd06+vWZwrmh3u9puuaIvBhKZmSEAsB+kg5c8Vh8KKJtFXVVPgLwaOPT07TcTySEK0W+A0vzy+BhHF38IcSUV2HIwnqOncPN4fAYXUny20jf0DmJ8xnwq4kWlxuAFv9pnbMHST980V5AagY16Euq3YL6Wq1BrQQgXiorQKLR93WtJy1f+R63cC+tZCsyKdZMHB1pU0vDHqI0MhEZ9k3miVqLBKzsjJ1ydJ0OBb9IKoNPXnq0AyuSHeB6kI05NzbwQIepISVW3pXrdQzpCACFvnieq2UG3CgINjt++gxl1O30gkuycH7ujWopkxo21tCemsbIKKQF/WtYHmD0PoX7LYtLaqUOk2gsYKY22jCx0rHL/3a0/rwuZfkkj58U6v1uziODeD9uO8AzPNksQAC1G4UzyPU40DEtqRWKEL+0z/Fz7V6RkRyRXvZY+K/g5FCcA8pMBSTiyjAuZci/MiZ1hqbJ/c0IAbNaneYW+/yeDK90Gu2p7/cncpqjk9Ol6PhDkwMuV9E3hM9+EE/ZekttA/YMTr0gkNJp/Im9s06/y8d9E3sJryIFoG9CP7GZAZc7YC7Zzeso6HIli/AZDDSEHaqt6yxe+6YAaMKf9MPeQpQI1fyAuBWaY2ZGRkp5gbNobhfVJJPRcPPFduSKnEC+dbCuqNayykSPF9E66ZgAq1EkkTDMM9icndxOYXEZ+r/3PKVzsBC5SedeI9m4DpTu0cUgql22rsplcc4Cv0nAeAU4qkRZiGiRsMRlSyD4/Whi1Df9n1YCnysKo1vpYigjEuITmE042gIe9CrmYShQCqvhVhSMUF+xklRvWAJlVwdxm/yYiSAAxVU9mfvPFOEgxdb7ZtxMjvSN3NBery7vs5DnvCHTG0kuWGjc4OD6hPjuk+PJ8MhMwHZy0fITY43GKp7jo03kpQjWN8xtuFQ/tZVOnPrfe7+/xRa327mXuCvJru63rt3x2+LT/EdLue0cZdaTDcm9gkxk5ywMTut6kryikG112Xur1uY5cSH9MU6GUJ6wFTQTYP8ID3qVSoRJr0g4ZwR3MCHec7moo7NvYa2SKJMtxoCJEPEZwrKUs1bNDYav94mUOrRyBNFjr/rWi5eH7ouylh/2np1xXLrvbI6gDiNIMXFQfXfhFzMEzdc06GELBRGgZdzkA3g8XUmUMTcBN4DdgGueQvCJl36FEU4kP3yE5NMVC7wIyGUcWU5byZNoopKGSXM7BY+F6+X2ZoNH8yMzf3okEchc9Fe0Q0ilCgoHFqW/NXZzFIMe7VJmolJ3f2gqSRMSV4xiap/f+HfR2Zhg2WR2H7izVFb+ZguW0X8vXf7QwuWSbnjVZuYrB0PiVhGcaq1OSMmzhy726dhIXkxKgur9EUSIjnTqFi19zrZ5rorcn3Lxzf68EhyvYJcVJAkO56F5scRHkoazthS5yrVCqXE9RXyvxtYTj0VLcuN7BCKucMEDKKhFPzx+VnIaJeeT7mlWZk5Ay1leZFRBhPguBz8Kf1SGcb9RqqEGqqFxNf7lx6C41rILaaVEuHasUjnbTiskxnVouK+fkqEVc6efSL5kFZ9lcg2ksR2GQSKfiRNGOSM4DbZHeFWWTU87kLGhfmoOvLeERpIJ4KiRhK7xwOC9llo5WI0hrOtRXbTiO2VCSawCJphxkzXkp4AGy6Y8AdutOo6IE4dUQBhMMGXS510te66gRkuoHmG6F1Ao7WjjEgfa1IRnOjw3G6GsQavosH69tP5e8rq8ti14xp2HAJ3FIqpJ+HWYcOFZwY3q5eoe1gtH80mMKLPZyJE6X2hbmfJp4nHmlCOzNrWjc6nSC3rwbxS9fNKeY5rnSF/bUdmbI1PQKktrWd4dWJYuh4qZoX05a80WJ4O3zcLekErXH5JP8VjJ4a1763fiZDtc9mjA2Z5ZISFY2hubRZUB9V1eklhCbxNS6cqvk1wReQ7Wq/rB1D0Sf++yXLNQj50EzNwl5jbXq+xqMw6WvhBSdYxE19Rt+WklDdeNuUr+WXvJAlmIvq26S6tU4XPkRl/ngq/rcrNRgXJgogMKVXWRUhhRgGKuohkrGjM3WSYWh4E4TmYixyDTiQ9LwDFUGrUjKVb1UCGJGLF4IW+kT2XiORden6BgW9wVVustwVN+6AEKWbNgCtEtm/rXssd4JGLVIBnH4LcTiAj1sUEVYFX40X8YfSFcj3wFTEQVhKTL5BHtEFCAHBZq89qM5elUtqV6gUACcjx6lknptVokwCeIN4Y2Z+FlHoRzXAjSPHhwztDhU7fx9Nq9cVSYadKANOINDLX+yckLCGPH2ljmZwOU2JGBV+fXqVeNY3K0w2Oe0eyXg97AvPz3TZ3CU3W3vNiSAWHbiB9Je64mMoQjW2Dahr+9nwAeCi9X9T5nI/KMgX+z/fw4ZVoMZlkhLNEoglr8Frgf4mNc0I5e2YAPP0mKl5kxca7FfZWzdiV5C9Pd0cCJ6CSDjKKo3BIfggG7gswzPxYvj+dTBM/w476wiU+fHUFL00raTg22Ahcp22+W3iu13/du7PluvU3s2putjxjqob8FRfxe8GhjedN0S5I4Ey1K+lUpLB5bD3WDJp5gIfTGuDy08Bga8oeHYA8mpiEvfPMZkoebdUWc9hMXpH4DJk1d7Vq94TMrVo4kxy6sGR3wju6rlnb+jB6XevFjqkNLjMuU5Gu+D4dP5dxJlX9ICirHhFkNRYvHfveDR6acp7ojqk+Dq4JaE7A6S7FqAehk87kNIajVervgzTktVEDEjT1v6SOsYKgCr9kGIkcJ06/qclYPWcV6BF+JgoGwzsCujJ1JRXuPRIIDdpqIWhFSCIOCPCpYrRQg2KtNzvIeDo2CgSl4Fx5CpQh861kNA+82CInSi9i9mLxzVL0mcpE4JhyQ6cTQzmJBfIRO3oBE2ZeMBHGp5jAqqirHwptAQ3No80WhLPYM0BO9lbqgm2dDj2TrjKZARWkHbjE3TakWLt/u+JPIkEOk/FFHrx8pB1YJWpyhW+aeeF9gt2A5xoNc4/g2CvhlDAFLc6jhRcu4gxcVSZhKRqLI7kpNJyC1yxgM+IHQHsEYZrMTv3FjY/YJxJp7T1WJwyEr434IZ3OJcTh/S2f3zg931pJy4lF4SjVb+T8acgwMmmRs6dIoKt7gdOMfrBcrSPtdbfkZObHdCf6lEy+qCxWNa1YnpKv4wu6CAVBQkRxjW7x61rHWbN8EHdT88Ng78EzP1rWEHRDiyWmL8qfdTxSBNpFgI8qE9iovuQrBqkjtZHsxbyyFyO/Xir0zU8XA+HwNajM1cnzJ296cQTQaHAkuBiGlaqFH+rLkaTh4Vs3F8KTrJjI8OlkdnIqnJwkxfG5mFlCKIqGuoU+hH0XJ7jDMdR+21GDH2y0Awfzk4tLVle+onjlCnpM2Ln8pU/T5Z3HvznoZJnxRZRrHpmrz/dY276fZP2tHzeTWr/XQFsjiclXEk3nq9DhIwsAnOlYVU7t9Qcx8FplH7OD0wOkeYjytSe85ZYR4zxKcu14WuHqa+FnEh0W6F12+epHSx0TWQ5CvAdJUo0I+vSEtozzL1rDuV9AXb55qxKxvmTz+nr5asx9/YtXyj6QRfK76IvnVcMBP14xm/9vc1+33EaWnHk/T1GivQYwAiFS3Wr3gKPukNRqj7alllbUWG7TCkaRBElYIIBBgaLYGkXMpa8dseEb763fYGNj99J+k36SPd+XmeevqgBQ3bPecbhFVNX5z5MnT/58KVuot9eyV0Tglx2KU/jWzx2H1nUjCrEARWVmcvkVsUMU3XPud825YkFKY8gxjntNN6Sr5plwbS90GjZdSqlo3UrWJYre2jnVOVuz8BJs2jpZ6vAS5ipYBWiA/TnTxbb3VJJ7/1RCVhGxKs+q+WgyIQUEPWnTzBYhgiKx5sSmnGG0DlHIS/1LOQZ49fK3Nm+kWXoLzfdqVupIhIpic+wlm/n/g7XWZbWLSA33so54GeL/9DgfFK9zcgjIl3rTWCr8rh++azOqSG4MFpengVK9QfFEqgSucCNtQfR4DPcp7toFUwrJfUg99meAr9n6fGdLbVTW9hR3Itb2ORpmxFhFA+vl2XkCdhwP2AIekUC3Mm8wBiOKzSuKktP05x4cOm7cDOJasYQJxwY3WK7VE+Hnbxy+F4i6oqPZNTVouzWEvonGOyCEiIAtCdlylYw6RseChtVLqKDFbaY819f+YiaJ+RJzHlQnMs34FJGK15Jny3RknY7dnyJx8nwW3yrXsAUG42hwnEcE0LO+H2LFxAj+j5cn6cVQ0QS0+hB0HWrgrpCcBV3flV4ctqsCB8eWbIOGeMGPe/XxYpPUuhRXSzc86bqlSVQ7ZwgzC0FxKVxCxIaUBVlUPlcr6oxd25nF9wN7JbtjmHUNbYf+wDq6V1eCstsWpkYDKpZV0jPiuf66EVNtamSYdCTkg22vsyaX8NAUzUQmi7TU4FqHchxVjPAmqoe/B/5QIYcj3kFzA8bRfxU6nCjOqJkYXccXOxNnEW8ZXUlGhMFcfeHsW3LjlXseVaV7HqU6eyyc31uizmhJEYK73PTsY/3WyvR7e3aO5W+hMIoHnDCUMGJlX8lbmvGVLXEedUKR4VBZVosJ/2RUHbsHo5MjgB5pBQ1dx5tIsm/QaMB1QdX3ov9JbkNBESXPH47OEfiJzdiLtbWqFY9U6rqftUVJKvpqhgSQOw1iRaSMho6+A8fO48uqm1gU5PCPHEzgozwfRVlpjtg9j8UNlSHTglbkajTIs95m84IObnZ62pPELiHd7GKZZJuFacZdCEQdLoFU+MYnlU1nBfuJr5fl2fc+6er+o5dPXrzq9OTV2O1UUujp6V6aWpUTEtwC1prng4Vg1fKEw09tCOJqmkwEl9B7rcDFTnQmdb2+ZC7clt9D/zv4WLF7qKHHevLEtZrBMVp2v+hbL0dwiF/OfPyjwqofXY4ny20nGokbakBhFmzxPFTrVe6g5V3/gAYAfCPvleIOVfMy8ukmBS7aFNC4gkiVktqIcpYqpkxLD2TzJGVS1M65RinxOe1Bfe9yqE1HUWclE8nkyafEdlhezAR//CKn6QXmLXMn2Zx2vMWo7otryQBU3btXj1BTC9S3WXSZPl6t0qZQWFdFq1hyMwW0v17e8opFghZREW2gPl777H9R5Z/rnlOTY6Y3DALxJygNIfUxHXr4CvUpuT9SFYaTaqJcgvjAFEQai1+DI7Ia9jXtNnGwooII4B54mB2br0hVJ+AFi2qpypigfONTOZw+Uen2MTpycsNnGj+d+pKBgzM0ZMw4anpDuH8t6xhvgV1N/fUf/1KIb9DlHJvuc4W0ofGeUdKsF14TexpqVs69J4CPo+6kcZ5xRqqaFFzLdDXGzSDJaQWHHS/U8vuvY++ju7H30eeSRyzNPFXccbXSyb7rH311v4Cv99fFZ5bJTJ9+xqd33dNdM/0kw3jujTNdC5CzhGaRTVxzl1nFWeYyPo5t+R9zn1pZQjOLiyuteGUu47D4iVgT1WP6Wtg8JUW30QkoA3OGj4dnwJXdcBWBSZAH4OmZs0MDVuKydlNJFOPX93SH7UYwEDEgEyXBauA4V3lyzcg2nxAlul9JDhP5dMFMNdI5n4BQgKMmcAyVP+HD80Zg3w8JvshkmbOr7yW/hcErHsrpjZfAYXxsqItRnOfv4jRgumUQKyauTzXPdnN6hGtkFANKssTMaloDjaROvHv3COVDtEELNWvIMNgP/pRct+f+mno+kjDR8VIJomv5MPF4z2fHPCoXCzhYU99RqkOce9hDB7VlO05lLa+gP4qjGhOKFwIAwfscXDWf+jQzJ3NZSj4/nIrpyzchpOo9DrGmL8Qdre8d2ptiW4ZZOtBCIyezunrMlppXIJ0LGQ5dWUyi0no8r/1sXuE7fbmEV7ZtDEWnzXqThwV9CPC1+caSdK79jMN0Q4T4qyYKzIhWHSPp+bsMMR9gzppdqJ7ylamkkjTce8FhXymlnJaT62qsBJhEQPJY6Wk8Zag98WG2sI0o9Sz5WLdMZUtuaxELmWjWw8iFpK6BwjVdLVka0EokS6Ek7SzHi76mKdEehhQz1lcvTGb7PaF+zu4mtB+SNn+dJzXXBMTo/lSSMWIC/eteP6ear+MBZ7uql5LJQJN6fbAHw3TDhh7HkSBx+/EZ9p+SYD1LBZ7nEtd9Rt5uf+dHZTjJQtZnXVXQkSmCSbeeDnKS9ZWzSFd/fje61sSu7gSyPE3QjvQ1I5Z+5876sjofnbxpYp7fQVFk04xeKC5xt2con18XyeOCqOZ/ccr/wRWxDLIXcHsN8Dp0MmIUifRQEwBlK6aHzcoDLs7N5+coHBacxrFP2aarcHTtJ/ZaJ1RDKqE0o7g5XhrQXtRAdBSJojLLMphw3Zh0rWeHtH/FxOJJ6FA75w8CHK2Pp8vF2Amd8UxKziAmAGP1B/zzYPeNx2C7U9xx5N05hLJNXu64VelZHnPrekKzTFsvbygtZQgGSggW0W7ndZqFeOtifLxwE7FFCV6YH+/WFRUtbj/e7Rf3+kxGcdf9/72d4qc//Vvx7//zwofzlEtx7gXrPZUG6z6VbP5huRCDXldbDdKf+nPrcydauzbn7ihZFBvBUQnG7m0BgQoyoDhgO/H7HkK8Wd3XxT2izcnzu9FzEc7twV4S6pjNXrYTBKeR2UV1tiX21x8EYAtTTdh2Nn7HI8OxEDeLe2KO8JaukAI25C0rr3m0iPtvKTclV1N6rOj8goLRE6FiZY86qy9GCwYsFF8VuIxkS8Ji7jH+UVlqmNyYMDXNFYK7eHWBTRmrq8/WGt5Qm7w9c+ZFfRbgC7QJN09T+mieeAQbZRkBx0ZXI0wqkGQjw2Z4QQbSNMtcKI73YrQkeJk3+CFPOTOF+tR2FaKiSkf8r/GNZl5N96ElJasCTQ3TrOBGS+d0rBp5Z0/97ffdlkiOW9p2V6QU1CCSSy2wLhJFRQRV0TP2bxba0e2Ti7d46c7yZH06EHKlIpHPMapD9+Dw8mLoaawrJsaPMaGEIw0zfb+Ftvqq6ukeXRtnBgkMi9scxACBy+7VYDn71n1+0v2sp2Hp+v664b1mVtUvLKQm/8oSq+pnPnIm/c5bY/y1U2dJ6N9kg47cBA6FMIcYVBfAu3yImZpHOxMv5/HG6pPWXFfYE/zte3C312+bc1lLj+crs/23Tx6/Pvz+wTPFyPdCnvrI0NuYMV+4E18NU0HBbd4Jc6xHMAYSpIh3Zxioe/k3375AHI27DkyXmImFZQvATpKn3jsBDUuMB6+lPllJmkUglbNw4ojocdovVHlxCi4q16/S5PNs87vOj5ZO3jwRBtM3xsLRVuqmbVKwbJ2yup4ehw2ELyWC6TTeNCzZ7YXDjCHUpGhxdJZYXp8lz63wi/FxDDAx1Acw1HAbsoaW8O/oCobjAaE2bkzXo2Uc180KUoRRiBqQGsurEqHgmItup5yP73C6vj69TygCXKD/6mJ2MuJPHQGeveOD0fTYvfv9yycAA3HcaLrsPrw8HahALuBbxIY5dgIuVnc626aeuROsuuzIYPa2Pj55Uy3L5WXFaft8d4coaNJnef2P1czdegaE6+uauPax1+sNRqiEfWCIZeEHFDzsAcMC2A2PmGJx8DRocnEfjomxETd5NJkdqTy8SsETy8BYnhh5rW+JzSnvqcImHIoR0RIxT4HnU45ugbZyYlcW19B62FXpwRKDkowYQUUNq5NBtsgztgpMhTyT4OCkZn91hq9+rH3twXum4QBLNcMRml/UiWGkgAIUO49fntHiX6QZjvzF3QR+UxYwhMfC+E0whntN0HfxxiDRYtOm/KA2q4WYF/V/otFYLOvTKv5w3VLUoeOliXd+FCLh7RXHl8llR6VlmRzfJqZHg2YDEl+p2pU9d9gIyttpuVBvK5sjaNLPRsvgWM3qJD8yKnH8uJz3mZ8UzdASxqkk2YhW55SihzudKlnqvoZNeTg/XktLqnt1Hmj9MXUPYOpQfybrn49GSz19YnNXTQff92G/MTl1g79BjNzkSt8pdj/78kuoN0raaDnhiGMAAv4Oge93/Ok+vgpHuDaqA/uYHdjshFTal0W7dUsFpOuhtEvTsPZZLkqs3lU3LFRujmIox1dRxDCLeoSjUBGCkS3HQvSNuY2FhUwACu8iOpPuEVzf5N1n0Tssy/AGoE0p7P+D4mgCC6mimAiv4rISLxA7+j/+hcPv0yAMePeZ14rwQ4PfF1ZVuw7iE42yLpf9aBekJPMzIcNkXjRI2MnlEtbaghAmH3Mnxh/7Y/7n44H5++CfHRAMS6EK35tggyl0F3ehzUA/mbgM7qjzFzs7O51myK2wpPIevPm/zsZTiVeBIJTtQ+mUQBIm5CjgFjw8/e2pmTdfyb2r4g2RV0W7Ji6L7nuge+kJC/ob5vZFQefjLNt5s7SLuRzqHsfxTJRVelh2CbNFJcwwVWvxsLAsGe4sj7KnO0FgPPWGM/fOwti0aX/PtDMuEzfE7yHVNJhaXkEoy8kVvXHdUhcP9FRmuhPE7sNdVqVnAGZV/ciQQ3AMxYJFjPEfLkvefamPuJgheNgu0uZVGwEX6hFlk0Mrv57RWqlPERHNvlkG3PWvqrxYgOuaO7t8M7x0+6sAUKUgGp1sNagO6L8P8M1okAjE76tkVE7djUt0D94qIQkyFIe9Bg2CVf5b+GERXosyNrYHMbWuVdN8knIwO6qWAzlfTrH9WNLD4aVH4um1/8CfGsn5pbu5FcSOZB7YlmsZk+COSd89FjganY2nLxzjsxoWYe+WPKNk60c7+Hgynnd7MYwe+7lPgnBr+Hhq0f63rFknvGtlEVN13DKcgwB7fS+quy96YWwBLUKnmSXe82v9gav96XU/r6ufHLnySduwwlQazzu65uuj69q8bksP+oX8gQ/CtheKePH4+32uygkroRRzp9Bat3f5ELPytZDFMEB5R11JwPhWKSHBw/r+kiFAzqpsEuWiu5qqLkARfyfioySS6Cw6rUWkkyQ+hmQITYBnB6qUp1wRjP0y6YYAfDVzu+uM4u7AqzVrmrZN9ZqprEB9xjLdWhMLDgkaD5te3W5e3SQXTyng1Vsh1tYc6ziTG0gbDaCaOdRfJFKA7HFCAxB6txfh/rkXMp8KRHVbOpgkarkrKBxuNTvUsyGmU8Z9p0VjJhvtSJ1/dVAGRYh/IhcJqYoVmxFKeiRCrpe/+d6q0n8ZJEP0fBNoGw7tyFEA4NOykj3Qp7vhLUc50ZICjdYE+hDse0v+3OpLMHxpm0D8+1h6L7qighIBdIGSzNa71cjJ941suHfBxjnMhTgJSJOe4gKao/toINN4M7zG4xGdupuEKNjkAgPAO04uWpJFWSiGEZ6Qmd6BXUa7mIlfTd/qT0xbVFJIUGlzZ/DFX/ftqz7RzZWF9RoGI+CTyZUB13IPeK1XBiQMLbbQiy092Ieac41GMiw5kTK8LxHqxY0TojNEaOEghCjwKpxTmJcemMo9kBcU4lt9zy96/Qg3RiUCr1MXt2RyJk4XcVEMRH0hAdtRsiftLqKAzC4J9CORaoaSElAGC6Cfabmgg5Q3dKBCWvVicDPZc2ntkGktMdlVCp7g2nztJULokLzGR8uriw9ztFb+6qt1VZbKMAL+o3RmnYzFYMaRZzvGra4CdGOjCDIfz1K6AH34pdgoZfb7xe7ngfw8Ihy3HtTjRti35fNAxbz/n2ruw/eqCgAE/BXVhBAz5Mj6reI3bRdfuv+H772rd8eVl6YcxWd801qoMe/EkiHLIpAtiK1/T24Nipb4/52aP63jgRcFs4lcB/2x3WHjg+Nz5MBqYMCoABmI/AUyKeQDG5dVPh9Qq753khinEM64mFy43+ItldnWe0Yz3L6NwAsZoYT6HLxnTXA8dhN7v1CAswQ49oNQihRT/qoXQ6XP6JzkMtwWG6wbxVk5l+HgDyZTaDlhwHi6nj57St1B9clTgtce9++1ONgxB7hca5RXmf3j1BsF+vTGmCvSO78y96XMR0PZFH+Y/+QS/pJjD0qjmfeknN2GpM29eD827cW96GZm7EEOKdwRVeuSJOGImViL1dFHOStPbDwi4b3gT0cwd16h493+IVc5MuhMeQQBKMYniGtKAJFbttLKA5VtAn7ZLQS8Z7okfXcLSxMovbc2heYOxm8Odt64zqAv164tIaadwW/APdLvduHd697HXtjcXnff2KWpGi2fun594xZaXn32RpLwXFsuWiyz+H5zweG6zT3uvuITV1n4W/PUsqYdpgI7+Fx1eP3iM/0L2eIO3phuZon0RCY5aMGgX3mtnOOulo0LUTh4+55zsV34U//ze/3w4zfx39ZkNGb2w4KYNJo4F2eSTnmp5IYNe/SUZlkplZOUB/ax+reLwP3QYhQvV6/rweOHXz6620loM2JIvkdfNkjtWR/qvmq7YSmvjdjODYzSyPDevRjsHp/wPG2RT+0D/UuPsKyhQOHRZH+WOfNTFwD2IvHHEjcD7zh1qqKpQLwilkzJuyUmtbt3AR646zpeVFspaDLeG2DyabIpmcpXE0JIBtvewSk6HhKLZRa7rmbJQ0lmyEOreIkH3vYtl6miEvg786v3/k9i7ZBMxYlxrS+XD/DSE4AUazbV2bRRdYhnPtstcx9dhNRHOCf2tKRy6KPZculYvsiVUR3rJdU1zJhT9JSjCiRx2iCBxTzZs2O9Y/plwhrdDhPbDY4BB2L1ZRpH+TPeQO333JWX3HW3JSEUCEn4ruWK2bKPF2dHZXenz/8b/HUv10fv7vj9Ef25hKxxN/z2vf9No567ge3I1exuQ+2BC7lt+GWcVZVrFQP1KkCWgOb6PnxxD4JQ2NdLBmyq2go7x7e0zbKfNjGmmdsuPovqM/G6L4EdkbQdDVW1gqvu0aH2qBwGUp+c1iCjxKDe5HML/WaTT3aTz12DihhOIEHTLwpyflx02ZMHqoG2i5f0hNtdlf5QB/cGxe9TmbOvYqmGMqrWf2CM13Mi893aIrNyMulWz6cKMZ7j+Y2NXCRH4StTxsAkBvhF1EbkSlYmDmfLaNSLWHCkvC4Q4H5aCbJZY2qSYYOnXQmocJGf1zrQqNVPnbFt3OIJleUiMnYljEhLmHON+GfELjl0Brofm59NY+b1goku0XyrTFCVngTk1/TKGzyvu71IlYwNvZfw3sSMZf2JbJ25FhrMVzKMQIxUCcHdZA2RoN3qkOXK4R2yLrZzRYMin1Lb2vMkdY6ppf+A98uWZ+QgzvmofLuVO70MqYqjoYUHMbajpNYCYzi65r99rzuZvQu3oHJpSDK4vCHmbOyNcnb4hkwxjjzPo23qKX2AIEhHnLD9FJMxc4+65owHQFsaZZhkTlK2JbYhcWBBfK8lCROHd/Io3vIazumF0UyNnuvUKf3wqbqjs/RdbIMWG4JaYxVLkUUj1fLdAT0Bnfg5TV4r6969G+nzPJq7fn470KTshVWUL9/FG6nmALmUXlp3rAwJUXwmOIjb0Q6ISD+549dobJjF5jQcCeZoFTscww4r8KvmjSgxGGGlvd3UfHPUoGlBroHdW2qFCzkwQkSJJCdQK27QHkp0isZJiANSdAqw80ej5dVIdJW+dffB2aK8drLoOxikQZtmZgGA0tGI0Z8lAxqeY+dElHslQZ30X/Y28ZBREkfKXrR1POf3Sga4qNW3k5lSeQw0s3jsAGV25szZ58I1cXjyMirMcm5uZTGN2StU1nsjVyC4DmVbTk8BXpBAiesZsr9L6Z9klI2qSCH0g4OO91PtF2GU0DMcdNC9jg2Z9ckuePMmaC6Ymp289v0br8HQG28Tu4d6whua9SKJ/SXwH9FVeaU14meYqDaUiGVY71UxG3UVVocM6x/fnjSI+4QmtWwJUbfbNbBrL/bWvfjuE7rp/pG+pNaQL2vd94qEZtXYQOF5rcztcOiuPkxfpQGHTeeqQJXPpgPPEDVUp6y85IbkQOYoWrwQVonNq6gDlF3hOIN02uOY1Xj/bMlt4XUAx3BxFGl0z+dPqZA0FlX1Y7Hb/HCSiiLzLwO2LINR6k8qIjUPfCkPb2agSLuZQxRSN7kECI79mIagfiTk+6JOFn9osaZZlKldEbzliKcFTgfwwKajnBJSlaF2LOjFqMGZ5tLcC8ffkV8hL7PKwaU/U8zpo7BbFwP7Ie06iTB/BB1lHMp4JYgxBBnDcr2cXXX3B8RJGoAEDiu1fas+xRXoNXTVejgdwc1HBEN+m/aVIfHBh5M/sbuWrnXZt2BkC/hMQKWKJ2/cRrtmoZav0evnp/YbnvIvbSmZ5wOEo0DQz9zfXVdZCHFEzcPNa/5G6MaIBfVbZ3vgQcmV/ACKe1xYgbNhlL3vjlcEWdpvMfNwS90PyyK39e0nU866GuJdS9Bd+TNgGVTXhEdBOBptBduMYiP2F+xU7gp52z/RZURf3ZFAI47+XXyFZqGnltb18Z51boxaFEQtFddSIkjFr9fcwdNI3rpyfPMcEDqplAWrytlibHmeIr9vxYgswTEooSzHYtRnDJ6mXwtIB4ExiHBsHz56/uzhk+8ff3P47MHfmRfZFdMAvhyJ2zq0jZCTxNq5pOVDnaDV7mUCZqK3TComfAetdAH2aXZxBA8dlTGmmkw0dVOdEbqilmj2eDSeSJhi9YeFK4tUNch/CgsX3IrPJUH1lS/MAlP3BjWKi29c790kgS1quLvj74mRCq3WEYCYXCXBk91k3NuwR7oxqPXNNduzTvQY3t5UqfRieeVHlVPW6vUYyrT9WjRxefvBQsj5+XVB42aXP8JXObU2xSlnQRlgC0x+psnImHpmyhDxoNAVk6DcK3If0wZtk5ebF+FtZYqkLdtgh3x+7xDbTP7evbvVk5RmW5Bdt7wa2KBZLIwy+DtGlIvuvc6ilvN4rsMonisEHCURXIdooorivQ5XxXsd1uO9zk7njfFeh+PpYWu81yF43WEU79WidQ5D/U7QsjR2UWMcMeyODvoQQqFNShdSrIVmGdxnh590dX0QKCT6LUY6Ce4Ei6lLf4dTY9aFX4Uj+2RfULS6jDhiJ+4woOrO2Ing1VKeuyOlqhjH5ar6dvyPYwvZctVMKZMNiu8QqsJ3EPH64sXtsfRMftEMNu4zJnBKajeTxzeadqlQSFvxMcF+QZ4ff2Olgy67z8Sx7l5r0EJiwQ0+5D3CLUhUrES5xmgW2bWPrqj+xq0IYrPl6EnmDutTeunLJqhYpg+KkWKzpC8vFiN37XVb8ac//ZvaheFak2a7mDOzk9wdtV9dyhMepCLL2REBc3iMi75qCCH8+ywIsjt4FoLfXJUeCqVk1hh8o8kPkx6l8XflfNx1e3TqtieAMPJQW3S/Cbsl5IgQDA3mD5XpllA1N1XuqoDlywDTrpcUw37vri5fRgBm5XJ2pHF32IxffI78ecciQA3cVl08mp2MHiy7OwEg1RJ8EWTsIULjDlj7G1y2gb42zPtiSNmnBCmUDLPM7IStFRLOOtlEhreqKVvBqLVOOZ9PxsdMrngH4YGhpgCw8A+Drmvjj67Z3l/e6RedgVXEgMJOr9dMbEpgHwvGG/od/2EV3YJDtNDu97OlHxR8dsF6lJ9o1WL7tiMNL6rYhUOuMpYdKuBs4hDJdmpI8LkFQMhqyx1rRHaOnND05g2/LAFZTJSsIWvQoPhefLyd6Bdgs2Me405RY0GVNtTOKl7NvuUXzSzDHRbVu1+eb+wzZ+jPYxotXKJlS3aYRjmmQfeA1Nf5M3EUj6yTMxRRzczhmOAmd4j/MB0f+ZsfYCPLsc3YMqdOmAIVa+wufsOvBRKda+Jrvnu0/7d6itLknKWaVODTNm/Bo8vlUlNquA893rk+tqdpt57PebMHifnG+F0DrDNwaT2qMycRX3YEJ5rSD3cPeIoFHYdj/xP5QHqCx+jT/DhJCDIfGS//JZkQF2o1/1FkEQQjUPCy/NcfZOeYnK8ZUSczJ1nvu48RuXI2Wj5Zji5YZC/0+8h99FZCUGIffNxKKcbckbzyElnqlgpAhhSaHPcXxTmNoeOlYbWqTjywIN4lpRPqMezEG9ebPRGyPPA+UkPOZoN0lDpID7oejRVX9ACaLmgWTQMX6TgaO81tTV9WYYqsweaZ+ihEOV+MTitSZSgRkaT29mOKLimX1gJIDzgjrwH1FbxkHoCrnJ0RWxcWsmGxhdj4KA68mM/GTLYwn1VjDSx0Z0qpwuIWK92CJZsXfcJ8COAt4EdeWQj2GBgSF6MEfsmfB+jDY9d6FxnvJhRBodtm1ZqpQ7i0DGbAF4ZNzx/mPiIXv/tKsroARGjAi558rD9kR0mTDVxBx90A127ObwRsB362+/cbAXhFbJSQi9br5uGFVPSonPM00WLykOeYX/LLaXU5B/e2ZddKQqJa103HnMgKMwTs9BOuKtxHMxkQpMlbS2MSA0x7d6QzpjtAyaZr9CmRSQnJupl2R5MC9KftnUgKXqPN2nhks6wbkn0VjcqqklfNS4c3TOf7brRZgUvcYk+Y9XeTz4/L6fFoEopE82CzsI66GrvY+nVD/1q/zTonp8c6mj85mqSHoa2+xlBFi69PojF7Px6/r6/GNcDr89Fk/nwadqhiT+NxR3W0u512SrBEzdCObLNQX6sU5Hp5L2/yhKXatJNHfgfZ18lJE3e7poxD0QTQiOlTma1ateskCJ8lLJFp1o4kTzGpw5uxy7vs4k5n/TBnNxtibQWLwI5dS46TO/a9jdsClOXKOCLlX9BMfrETqxrlvBtP3fy81sDMncEX/cKSskPj6OqXiiEtPER1bnM/4mt6pwHK8HRJZ7dez5KUGzxCNz6iI4blJ9z+UDFS06KDIS+c7LBYXnc72xzY9lUnMDOe1mvqEKbQUI2kitH1k6G1TqkoN9dN6ufxpO5+sbMTZvAHNzG+oxmK/6By1+HliomFiWC7+GLNtOqhirpiUm5pbE++bCJyN2DoC2/dsmmOFkwKNa8OXm2f58sTl2hYCy2ULIbNdjuFyxj6+dmY5WG64Tx77+awovGCfrmD9qIdARh3pfm/hv4dvy2ueReK/szz+b8gM/YvsnSbLYi7GtZ3y+rlkCLpztBG2e84mS8vdEl6j6S/DyaTbucA11Tmt3jTQV42llnLkzEOVwQpsbpaBPW4xwNUlVpQwJN/KHJE+SkDuFDL4Ppgf4BbdTl5Y1l0HSVdb7s17fSilFM01U0p90MFLH/xGmwFyvfNBc7HVuB8rAU+1lMxIBpnhCwPyDF/oB1w0zyZdeiJIg24B+djQLy6NscnTedXloEuEeF83oQom5XUI0TQXzEx3dpD5J2IbrpXByNJUqFDVwr8Wlsd6u2Q1ZyMHUUVli8FhxYedzPZhtPKrGWbHNQnjrLcRTnvpiSZ+gH7d3XTIbJECGe03GcN9OhrOp8YDmTt9IhbP2fjr5A3sLgfPXlGoJTod/k++b3v9lJh2WJEk8Q0fBbPKeFH0pJOuDXsin0Dvaf/pQm+eCDjUUSU4RU3Xtv7TSeMKeg8sDWiDV3NdFnbRv0+MY5k/fqGaNae9uP+rU3P5RejWxtWLxALtSWY0G1VFq0nm1R96Yq+MBUmfsBIBZUhdUbbvGp25I1snkihZs0rbvsn9kBL+054cPjwp3bD2kV6+0br6LDRtOndGGkU23pAOyUUw31xHRRInGvixKj5stJgZ1hbJ0ylEqd/Z1vaX+yUc43OTfktW1N9K9qMUM1x7adKWa2ewUFQnfpEqcmsSep1iJMT1PT7l0+7kcGhF8DftT3JeYAcc14zbbZHcy80c6NNa0KfpCd9v8GCxnPRjde1J6J+XN/NKsqpjJbMDfeN+55TnpXNN4+6vW5TRNmk7g8BvdoCFf2RUqtNGaikaUmuRvqlq7a1mHIdU+HvxQeEldCgjBt2XBHVaYtoqCroCDvEd6LBS/+C7bucIBNnMhx5229oIq0lugS2D04MQ9skh01HFkOI32+oxCdcTrqtJih2PK5ho546HoTgkce0VkhGszm0jHReRqo68UBksioR2Az3TFCclcEwJVulhqzpJRKQC75MKYyHEqpCSCyKx9VxORe4R2WCbA6qZDp6M05ba+76BCm0lZu3AUDcYMo5Kcqzckx/G21bdc/SWQYYHp9r1ljPS4fWtPehlsAEjbVRiwYhp72BTq7Se/CrejeKPrNwqAnRy8/LyalH11QvTkOPsC8v526KY0ZcTo/PZ4vnCLNRqZgkoQKyzO3SHfyjZTmeVAM1LUn2Kn4TCz+YxYeG09DlzCqSvdIZH9XSEsbXLOuPAuW33ap89ILUODs9dTzrdQgPLPSFXEnUi9vfwL5cqabYNliI3S8IPD0wMBU89rflvaQF8bfrZE81YCpYkH1WLxb4qogb/50Bfd3tNVfSbfha7osSUPdF6JrG9+T9Y8taX/J94uEJSgYTiBwEVlzGUsroF/b7D5elu+Jdbx+fj+de7IhuLtzvjlrZWsgSx4SEK4Sf5pypS8BHLQO5am+22YjZDjXboVoalfY+Zg70Le2LGqN2P9KATl9vGtaZjVjRH/yQmZDVAwRIPkb+Trvoq1LztHp21vM2+il3n9h9/tYc/qCWR9Abm7EoDCCp6DVwLexR3eEWsyuFZBpX9NskPn/MjmaLtx7ANuqaik9tmoaLknfj1h3NTunQVm37yBHR7d97UKdEGO7YtlEaXzd+3QJW0oNsyib272MmkengWL2GCcS8Az0eKFoRPjFG8aUwi+W5o7ntwn/l95v2zg7Cj3Yk0NhWenlaz6/lwp0+8FcchqNJXeUFflLTIVHwVllbvMwsQ2a0drIjrnHxFj2YLdaGe84ADTmVsvkQgeT23YnhCLZvcss7KgSOmx8qcgyNudMqyUnavE1Nbdfe3RskVqYM4A6xDXsKE1hDnzbMndpmJJRU1lz3rltX6ZQsLPkYgJrYKwPLo+cNOHvP6CUSS37607/1V3msCvzo1awgGGjXgshUrrKYLniPW7jk926JNW2ShF6Kl7a17UNl/R2JAFjqXSTNA7ccyKeMG4kcxfX1BTKpDhrD+84BbhvFw2isqDUuqUIS53PG5mnUmRvDE4wGA2QOWAXwiR1zq8szJLhBnrIQ5QFHXL38iiZFoOZ8s1KGJnBhizrnlaAq6raNPhPhzxwvEaItCXERwGIYiBB1WxMWWNPitOBdg9g7JotCknWGTpurlk+kovKnjndJFl8uJMU1aWlcqbUccQIEW+/H7vdoNkjHASsyjI7y76DYIsFsiVfY0XXwM11cThgUoAAjHsMenXOnDwZ/yev5xMnsW3cH97YwJGubURpgxIJzzjB+lovlV6nwBbvJe4s82Na0qkRl5BMhatPocv5eMz42lKF3SSdoNeXx6zIUi1rLd74opeWDNwPIG169LA9jEQsOuK8JtOndtqMMn95TUlwlXz359lv66OS5PZ1U93vI849cK8zu1Xnx/d+Ie1U/FEDyUThcuX8HnbgTL5BXAF4hfitIEjjOjTiJ9BkB9ZB76H4S1mRzFgVByta1FwyBtBoE2CVxH4qWz1TLth6+TcnaPgheBjSkSJmv7he7jIbhr9/6gKRTC5uSF8Pi+/J7L+DC4eWZ2zC4HGIjWj4ekg93pdw73dGpNCsUp2TatsXME9PvMXVvSYkUDWMOoLI/Kk8abExUqodsp/jJ9PN8YCIdi7Zkoh877jcZQ9vYgVTbMfsMbzG0y7SVwPRq+lO9nN1yDdnDzHdsF/tRHOajOc+uEfLwkdvAy2BIST0thQD9Nuuq+R/PEVowS3b02Sy4smR7dtuJrHE+dVSwF31kRJp/F74Atba/NWfd+O0tNqL9OZtlU8ShfS2ukSFL0lAfSLVO+kQ51frHJRJn39VevXRLNC9Ha6CKDzp11mo4mGO5MKknvqxynAGWcKRsNGKqEMtnb0UcdO9sY8qvfGNyzCfjyoNEz97uRfsj48V9vI+k+bAkx6CruufDW8AqTy4XpYZcMZYq5Qx3pGeREpcsc+337rq+S1dWrAfBSJkEd2LoEReRd6uwX8/lhTl39vxFzJDNEDjIqaOV5Fe1KyNEVle7+vG6X1VgjlGuoSOiQdxSC5pnuOS2cJl9K1MDgZubOprvlNcKYwr9S9eqcUNka5lscI3Ol1UJIfqRPr/nDXwbLW3neYzpkMWvBzi0X1vkv2FnRCkV8gKa66GhwPx9yzLKpc2rWzPpod22GnHEXjJgcdTaxFoZCyCeFYmJ7WNdx7KGfOon7Vq7sAUYUji+31CD2gKjcOOMHBnTiuIDdxBduGppDmUEThcSycBLI1Jjj3IKPb5PzOGTP8xRR+236Ud78WTHSudGX4NWsS5yN1DJ7gYOBx/q+9mrgG7Y4gdVDNX9ama48Eq1ODtVhbun3zeICBQeaZDMS2F/UWoAlfOOS7tW20wG+t3UrPQhNyyl5U3Ma1u5aMp+1u2fCNx6QgXd4bI8E/d5fPPk+xe/f9Ux23QkgHh218skkwHX3szbbQqHlRv40xUFQYZt0RWYEU8Fj5Wcpu7mnEXf1G2fK8NwIlkiSpL3aUELPutNiMSltoNIt1AbuCstH4wlVUyaz6/hABZk5NqarFmpvJaWK2TNx1XfMuj03k4WARMUEvcT3Im9GwkEGFADGyWXzThrK+ts4erw8MlZsavVp68AtFra6LBe98ci7V0RdcsLQCkDCFalbrQm/RpRZ47ugn4Mgm2+YIDd5ZfYItxfda32pJ6gzpOfAqgS1HtunedweSvP6FWbuyTdVMvaQpZkSPFjN9IlTJCJGaRX421OCsP6p9O1iSyhwTge5oqkWEOvIpoqdTgea982IcPENZGf9xDRsCPvT0JoiLfU0B0L+mzFKL2C6WmqWN2DoqbJif0/NmRV/bb7S93TZF1wYfvsaVcSQ7ldZ1TQkRtFTrALhVYBcjz3PK4HsMzuut8NN6lTIK7r5UQxfwjKfppVXC3LC1hWotDEDn0CJTLxgsFDd7r/cPLhy4/b7r9ffOxpkGKPiMp6dbxve6Mhxr4pWB5C9CG7xwIy87c1/p59CuYcnZdkgJbGPrr0p4dYdPWIJz4O6ZQrlsVIik+QgXn4BxMLkArJYLzchuLjAqj97p/f2jRbupHx7du9qK3W8EgQs2wYXigExmf8hteMLuEIxgJogQduT4WvDOfGveiFqNV4dX16vJoLlG9HbJzqBhXXcJPgUaGCvsJkCClwKXXtfGu9wbw82Ucei65irsANhHGnQ0s61+yJtTr0VBmS1y57C4Qi7lhqbPAac2mzpdZao5rIYu09Tk5PBRnYS7IO/SRVXz9N0AfM4EBbG2Qr1I+9aVK64HHP7LVPUajvNQWYrX/U5aZkfs25+vJMfXHT/bQhVpIm+fpvYmlGhbg/d4IoJHTZ0psooZfh8WlLkufB7YL/Yk+IO+Mu5l39TYCY8HUvKxqh2oyLO62VnLdUAmQa/01EcWFYzP5oq2t6Z789v8Lm/bBq+zeRIzf0p2y+cjI5VF5sdXZsg9nk32iLRQNt6k69A0B4OmQMrXakEyLIh1n3wZKe6bsuE0NLET2NA0c+9Mg1K2eaOeahlWmdCbfulxfTahjTgYAQ+UeKSORq8fnYD5qI4Y1PWN8PwDr5DGJWwVwVFECdRNWVTx/2Qpr4FSeFBtKrDjPjPYpU0wv4NemsBN1kTmlx5D2yRUPbeydE7MvqUsxJYvP/3NH5qdb5l4zO96HQJsi1hNP/cgH1++4KMF8VTk+ZYjx1FIvo4FXSjDk9xIpGHxc/Xtjt1pOjWHGduBJQRRoTf3v8q+7W+exqa5h4kaM8nYYtU2mh1kjNO4MeRMllcQvV5LI0tit+hcAD0incNjnprTz2vWh1Lo+Qbv0louPNwSoe8wLB7NuKl5tfD9SsqbcDcyzf+IIQYTEGzWXAeu31DUwvektbZ+MtQpDbqyJNfNH3wKeaJKP7wVwtTqXEoWNHUfqLSlNlIC4m/q2cKWIsKy51vwjIUi5eh3CpaSpH0nD49cXsZHRfLIVE81rjv78OOLjVhT/4ezW58tvDo3KxCNcF/1iA1u4XfhG8SiMsRvQoGMG5ijFMsSJo8jEhit3XdAqbjqrqcCHgZ5ysAX9UIU7cDt72sIN+0UlR5ITRKLniYhVA5XiGSAdq6rpfmO9tACTSYC052yQ6q11zYxs9jZcohkFBIGJNo8lqk9gaX3kUMQMjrE5mFSsy2sDT+gVdNyKPDThAzJdiOQQWCaB4F6PT0WI0Bf7oBvYbNBnsN/5Xq/2m7urvipgdxSuoT+ZjzKoPxMjr1TiMYUOTFvwVazs6Aubu1R2+TSwJyB8PvPC0xjizwveYSEPm6rOtProbGWjalNzz4EgU1w3mrjJQ4sL+xz/qb3P1ObqO7VmbeVC3jGLg9v3o/fNTc8v3Gzdzor+Fjqm3b5MPcq2Df67uNa1spLvVcACzUGnv/izq0zVraa7RuWZ1vkaneitzNjdbnC/RaaHH+Rodq+SeOMOxODoZu4nvJwEp5p/RqbyoZAF+VRw+I17bMCCp9jX2IqG96GykATUnxo7RHoMfN+CL6NO37Ga3I0WEw7IO69HNqvGlIqMiX23TtLU63A1fSEXK6aUkoYdWluQXVjJrF7OzwSAQhCscQQK03GzosbRBXGjUVD0ytBhuHtyZVFSfwxucee3j2XQc1vGf02Gh/sbuXswuq9FqE69sTI6gXmUvppcE6+STXMr9fHlbS5O1N2uyzcLbxgSLjyE+RTjFkcDaLwNCGgRvSAYci0YI1UciEBegNeEOCfG7esWmQvBH/NdS0KXcpQSGX7PrWZQoJAAZQl6yntJ1ze264DtpTXdFy2QKoh7xb6WpBoTJQXHg+tcXCFvBM+SXTMN2dc7E85IihWNxvZsYirLyRThMV8OY9SIgs9ICPlGKfI2UCaJL8qz3YuROz+NCPf0d+WjOZd7BBT4ziYc+mUFsfhPdUvmRzCzNPQcdqbPyu9GRrz7a5sdpUKkFGUDjo58d+mdOCmcRqubb6rF2jqt3TIUhH710M9UFuuxNOvHGsEcOGNgaxtCxdAyLw/SFNcqnWZugdeJW0AmsubL0jdUmoJyozcqGWqlPYbUEGEDG3mX5HdE58dej6l365Zs6zgaJTgA9j/G5EtbLkeBVQ7URrWowCv1lF9sJTrw3YMYNiLeu0W4vhrsFxd9xTwXmNkMdJXfGy8h8Enqy9ohrsqzWVCdxbeu9QDYCY62rFzZBOf2gSunWOYixV90DxV4NGKgxvKtV4/XVyULnStxVatzxtGhTsA6bwVITxev/K2DU/xyl659N/dD9NAz2GqJ6pLNNw/RqKCg3QyFpQFCJxCC+Bc7+Rs6gBHARN5dthXOZE0o2wU3pBeQWx9CuysUmHCip/PaqynNOCbaw5tbuPoE85v5Zr4pxjb0qj8Bs/BXW/R1NmQLhuRdLx0vW+LYl90N/H1sH0DedefQ7VUvcb2x5zxQPGeqgtONo76c//R8e0D/96X93/McNPpq0bAPX0JckbJ9vCXUQ3q+pbY0nyMpKFqi8fPK0CQgw4POtcA1uU6ZIUEssmyeBBfswB58uCSWAIHayLogWzMeBxxADpo4GL3hNwo/lYjah58pitpw5yWwbCubOmxBCiprGJ71cCWJh90ZDorC+FbuvSVBNffVlN6MDeB8gnGrV9doaM+AIAmmqW1/8XqR2i9Wx0nCdF6eti/nELdDh+Y8drcNto+AdKJGJTC/nBM8J4zuPZ8jyYLoDaV87xzl/JGNqQM6yOd5El5gstqwlXBEm2Ovf4Gc3mZKkdsUT6dXki8norDzGUc8KB7bQg8qVuCgP3T3FJxPY7XsnP3doO/n8NSeym5WU6a0OkJpPao+OoeZvxY8qM8J/1OZwAr6c0fr8yFPk/nJ8cTm5hH1ZBuA23nePJtxqL4ib/A5WWfvcYl4PK/jJwMAgFmb/3B2bDU/9E2vuEDlWlvHzpLrHoZYYX/yjus43jr0+9TYplWOjI/iL7AZz03gykaVmGQF7WUVf6ca9MZVdlVNxGg2YQXmNFr3Jm+7sUnKVrSSmyMGssTJofbWyvY2p3Xp6X8Ohcqo8HY8mJ8zMs3AMlolo5IsK6TIcucC/xj+cjqAUeDc6VBKir+/+QIXxzkvG9hWsEyYrOJyEqnFzPRmXZ9NZtXQ3V7DXyjQJ4shqDipKu4jbJXK3+3P0h8sxgtAK23l3Cz+Gjup8l04iDUiXQXLa98PlBFUHNo0wPhwjd242K8kM1krLlMJxOSnk62ot4JXNoPoPa/Z+O/1flPMuMmiNs4uRv6QQBTFhRVcxHYIUkstDTifExJ9dDb65XFAFdAgqLOwMYMQgXqdMwr3HQ935e/ojZQ9N3wSmkr4ttrPXVgi9HVtH41bQQ/C61AuMk908GjmsksuOxva676MrQjufzk8DZQx7+XcMbpXW3TiW7GuyMx2bFqI/vJzqx+6sB9e+rEaLQ5UNOq0d0awfeZ1PZSOBLkbHl9i527KGfvOwjUeXFTOqW85M/9Jz6k/is3UT4rdPnj7ej7USeI7bO2WtlWx4lTKAWrj79erIYtyBq5YW/KxD44SkCngvSVe/Kr4kzu3dz/WfHpWAV8z+8BgXtG4IAoHir7pAKjeo/MqpK/vs4UA8KyUjn4ngojZgM/RGdJ/omWCmUFMi8BvRK1rxj3t86i6aoyidgFf30h8D3szwqQeuZSTmFRH+Zbxa0UbJb8H8UkaaOw99jHQ5NuHHk1G5CV7kh3zAFnWZrZwdeIhW/eVHnd+k/XE7vhBcsZsSYkaGDRVuQIpCiEziWS4qcBF5Ir8KzRLGX90aJSUJO0p6EucU62lK3O9E91AWZDesWwbQVfcJU0HQ51sR70Du/rxjfZ2UiuTuzZNfwNz8kYrpytKcHXRrhx5uW296esSxzMHVAFLuG5g33fi6r548e3z43eMf9uWjt+KhcXXwFmFAImK6SvQ3CAjoxXg6NGEN73o95i9bJQah0JuI9ZsGEHixqQrQs0uI5PZjaOv2sbdaWuubDHTfT1xNHKmLGtIP/2SvJuSAmXzcC8KZbLT9weXUHaeAQeWWlVCijKlnQ9abZ9y/dASho3Ij5XCE0LVsT+90eon2967gZNh0juB/+wPZQvz2A80hOGqH3Ozyd59mF3eMPbEvf/rTv1IJl02SEL9bWr8T6In605/+x0CUdjZSd3ISuEckWFXbnhTXgMIcFux7MV7S4DJ1XwEc61zQLXicMom2lHW77HKq1Q7g9pyJqrZ1bsB+E9/NOqOJOefH1muPASjeUPsuDOnGO6GFXLg1GowMGZvIVjEkQlyfdi8ONTLpZjvqZMfn27vxQkSzytSYP97clqFY3BT2ZG7SieoN2KR5V/dkyYMs6Rqsjt0tbPJkupz9LewCHwoma3LzoemzO5zjxnDT5h2ZXuhegbbdbnol1KzfuA2j8Yf5zSur00crNnZgfwBZeDRwg6N7enlynbauR45rf98dipfzkINUiqAfEqvHjbi4TDOVNvXOXUMu56t6tj/QoD8THjIGGleYbArZBVhlu4XbI3c9CjTSN5rKdwUFPZm5IZuW/Fk/hz9EY1Eu28B7YhWrTs76FAYYkX6tq74pdoSUtc8TzHpxyJlEPrkw8KjY+Jdi99k+HY1O1DLRYE7CYMHJCVDGIDYhCVpb4AB5uRwtal7EKxhaYmsyI8YGHfGJCMm61IuqrUODTvtCr20uT9+X+blOZmcb+0JBzczvI4Bqtr0pGnyoxopsoL/PMU1V6ULe0MsoNvVCUTShE6idGQxA5zJFiVsYIrOcvpo7lD4VOlmVQY96pTMKbNOh5dMOxjidLv3G+kwy4TaiIiR4l+9Jed3lnhSKcm/rUDNOjhhXgefRfYUw9qhqUIh51/IXLEbANbWQCXtquU9/mF1ysCwr+AQwZXXixhtzu/08YSXiLSfj6nhTm98GMmuLtBtE2TS1g//BjH3j0Qk/6kl+b8Q3jEspVz/DgqR669ZGlTYBJajeY+U1t+2WKoQ+Kc8gmSJ0hnE3kV+mB8vFytPKRsRcSY4rPlJx8ZCqslwwhT2cpWMfpwtv04QrybaW3Ritwy6Bi8zoFmCRlwvv/Oq+cmvg/vFWKHQxvxmry4fouLJvLZEHUsLPef+4J1zjb/Cgew+3L8oRgITvSkZnb8/tiHuVJjMz4GJ72AL+6r5kY5thzaB3wUpgBdWTnd0VWtLGv3WzLVJQ9CCRPjTo1UOyGJG8Sn3LLMxJNOu69kN60xkwPBTyCrYskVEEVrgYFMSmoTtDtQegl0LOa1Hh+7ylrDGCsAdtn0PzD5Eet6ZRuazQohNG0IYYyRiDUIhvFCmxlELS5hShmNKyebqhKpgVKnrkCXf3r/jlwJYF476B5b2UU+M+5+sBfzhCIgEKVekHbn3kLwUv7GmPiFxm9lJ0zB5ElNLapRVOnS1QQW6jZC+KTgSyIRPea8fvSTsFL4yLxp5Vl0cX49yi3lotVqBl3O1NNA8+uhysAUtaM9bYzrIaB+nGrrBkYn686xraAONp//HTx49erexMMsVrmnwAnv50dLpsohe+fCn559qHaA6p15KkVUWKo9l75kTHkVJOBbQLWyU/6Lnwi1HJi/VmST9QZF8w+gKba6Mp+vFtAhEnx4LV3dugpVpbG9xqPD3Evgvz6+2L0fJ8drLaDwgfPpPvImXCZLJ9zP+Iw+/qSUT/H8GBVU4wKedY6w1KkxOH4uJru3n5Z/z+FUSFMAivaxG97qaEIF+LCr1ehYZ47Q9wAvkVvJxDfaJ2zoxeMrfhTxnJ0Wgyu9r+cbRoCRRsGoRb0cvlQxT8e1dOhpPUE7LYQA6BZaDrhRJ3dd8s21sRsPonM4ZHFw9OcYaqEOBEZc3SQpgSCbp2HfsEGcakq9jfIZVkassCEZRJeDJRFlIKnlNb6oSGU6jmxdClUOcLxAq9jYVb1NXXk1/w5lABrnpue+4FYO8S5nf9bsyst2TAQUSZFqOLuZtReIRJZgLUpIKQa+inf/3v/GJMbPBBAY9uaVHiybWAscQga+N6dgkAVKiCZ144krIqrvtuQKAKTBWeXzOFynUzIDmFTmT1Kv9Rt55a3lVjuX0ku7yMDThY5cI9JMMGtroAAvSk/RzpuYN+wqUQ/nfUvXRlNl/NZn3pReatSAuGK/QKy9X1zmj9wpeL1PdgwiiSubh9SF6ozBXL6e6VymDaBaaMmVwuhFH70y96a3cpqQ1RQ+em+7UGkvSVBUutO0dkWrqawE2u7L3VFdT9blGDQB9RfbGm/E0xJHuCaJk1sgrI8WcIRjea6XTRPq4ZN3WyftoEWlzg7bzcQxOG+LFsomt0F1tS0/auqimnbkPfpOBupMAatQLrNvBLw3eLynmmv2EiT1zAz90pSCj4TdsNJdQ5OK7hE3qgF+EN2k5UWgKzIF3QKvzpTWQwJNYIMahyhr9alNMK9hrRy//hclQtBdki7ZPi2dw0W16aJy+Zj8b0eSnidJaHLoabDi005ZgTnq/VkrUbTD9EAji8zQrfQ54XjNZy9UxtwD/OZhfbwqdXkC6+ejV7JLdXK1WuKeRZFM1HMbmz+LhFhk9K40ukgO3aLmNRN6ablN1OCi9GZ+O2+0Pe85f8tnvLTftAykWjAFVuUA2J1y0k3El8vsv5qE3r1aTjscJSl7LgeC+PLza5tyAimDQfJC0taQ6l8RFiuBKbBmv42vfFGxzoLb7GUNvmARq1Cm9HFZrAk0aKQ2BSDQ6si8yMM7QgzWsgMhI8KoRsOhnGQge0FpHxvB6TX1moqGiyNPUwgbl5dlW5RxpTVtMJSVJWzxYwWPQL/YNhevNy4s7VET399U8+X+j66x/yzIn3U0Ugw4/v+TdfnY9PRtuj9wwNYAQBHjy230iH3RSwYFk4RkTV9+zMHKfSWAWyMvcI3i6MH1Ct9TAc2BB8WJsd4xsQtn0OdvQ7mftGByqfu3hDWJPNkhXHiWc3j+tfG/xdT0UcX04Mjaf4e8eI1IWjlNya89kYAtewEKu7eTcn+AbyqvJ0KZcLhvM+Wi4muN48ujgRlacTf8HtKncPsh44+f8YOaqcGL9cuD09d4tfLhbjdxLNe3U+Gk1ENrOUS1ItG3DXE/4wRSwxr+iZSIlGj5ezsV4LZbVB/sw24I0NFdwmuj/90z+fjCbL8ofiTrG7s9MbKjhuKTeifjLqOL7Zifhn7niDgeq0XIxp8RUQMAzN1Q6UZWJMyDB8xDZ7oVnkYjStANUUgLJ4Y1mWZzl81pLGO6HO1ij0zQVsC6zoFThdjjG5HjfDy9mr2rmcf2oreWLI9mYgYUd7tqEOnRAPSlWrgnTViP6sfUWF342ug04RgarfMZGCeVC03Rh0EQ2z5gU2kbaT2neEQHw6P20xAJ5ok3ziB+l1qyTWZ3TtdnO6gxcE4yyPDEJFyLlX/La4t2MtY/89WHZRHf4ERqd042spjZ2wHRfX3eAYK9/PZ05qGuzs7N7rF8l39C4cvMd/rsVyDUfOqnKbbyjrssGy6E6h6/um+ngMxWbShrUJCWhb4QjYpDHCvFuDIRFM82LHky09/LUpobHz0/navMeOvXzK3IiaKtJMPVDTFJkkDU7uv4JYMy0+KwSdk6cBjsOzvtqqNL2znBeMt+6DsUsFEY/cE0a3hc8g84uouuV4fd964M6G/fPx6dKYNz1GryrJM+DEorPJyMM/8H7F0wOVQVn12pLy0WpJCHRLPd0VaUgtcPKMkOkns0uYePmkB5zz0qcUHKyZ/Q2AUCyRs3GIDTgBVgizS2fP90PjBcy9+nf94jp98gNQ8KYnwyLI/oFDVZhK8IuviwOjrIBUKMQguPCeUYauS/pgxyruClPmeY0+iRVV1Bx92kKHqvrsFwgswS+IMOREHwOBtXJvziP60WxG18mweVGzOf5sii29FQZV5386i07mR/HB+57hFeef/GCfXPcA9/5Z3lI8c3IYLsfzbidg+qErrACrs5YtoDp8CAwSX4phkuEX8N99HtydKPetEaj0/LXgRrvlBv/doIDkrxau03tTu7/n3oEyjh8xDiI5NLfe0EJv8KMbTohDztfCvbzOXoZV2NOP8tL6vFYwXgZ/shEa2i769rQv0sJ7NPXeHW0/6oNrPLjGg96mJCySzmrajXa3sF8qopzIfHl0lPgpeiKnd4int2a/EKUdI6OEoKQZyVc45aYQNclD94uPer1Ih2DqKwtLbVTlRHnv0Ref93h2ekqDRGDPuJL+yFMhnDHKfMmboTw1a8WeZ9Yzs1JUc9xkBRVULQ+wmHBafPSGHD5reHVrco/6zOqcjlXJ/8C2ad/cb8bIrLfjmKomSZhaEoFhNNnmaN+g2VcPnJBLCNPFI4vzYYmjmdNVXJpkCmVKj11fAHuJ1LfLXpjW+CDTKWOKcKTyDFe65+nUjpeZM4uc/XaGDiOvFan43N3G6IKGc37tCXlyNFnhH6Unih0u8Z7lK0eq+FdU7BB59/2fycK4M2onDl+nu5gYyMQPiFCfULogx7AcVT1G5leJvUrsWJE3ENUrAvjA7MEzeUe7X2rfQTfF+Ock77zT0nJARvDvNXZYrYb+MVTE2tv7CY0hhMrdvrvdKedxmrVGassbIx6k25ijRbBF+r7DzBh1xxzJ/Ol9y+rHuU428UJfdcPjUGFQom4ibrce+5RRRdHQfFzavk1kgLZtK0ez37OwtqqdJtm8hHKH18Hz0/orALr/+/+iP3PbJ4Lqjv/2c5EtOZk2nxoa9rJbbZgW4TQmdeT2Zdmygngd6W2U6GceOM6UQ2VB6GKFhJMMCbyRV5Jnl9rHGfkQ5W+6xXi1hdVO67M38j5RB0yfisg+E8OuoMvRXP9lAQWQfy8+zPCZO17MqipLwuqzdbspQkQEDd7i6Mfxpkbe8/F0GYww2/iJ6UKIoUFlSyyNKhllWvFZDYJb9FILBTd6Vh7/cfzi3LXu/ilP7gjK0bR8Nz4rl7OF20PlEs5fJJ5HFxLJDN2UYv/fVj2ZTeug+AanmvuFuR1493zqc9pIZETrZ0oi1MeZVKpRJqBVsAey2jEDpN12+N2T718dIphwP2wijtzj8agsQoYPxHTmDBWElfQhZi+yeetkOn5gH8Q9lB7cvp1hFWOBuD4erTh83BMFAnv88vGDb9RXZNrcSDoIr6z6cmdnJ9oqa+a2afvp3Oq0NDWz6VRFoS2W2otayotyPmSEKe7hhF8UYQcnDxws9HLOrVJLUW8+FReuTyH1DuY13RSwH4Q9QZ1yX7b0/aZzWpTx7jOkkXmIhK2urUdkZy8191pdqq7duKBvn4zgcn6Hf0uqnV/bAaf5bmq3MCr5Z3MrpRl4fDGfGMe7NKCf9SVV3nKz23qbEkVCHdAQTkSp+VE5p5ysZeThk5M4KPlyWl3OJT4TB/9HTrk/pQJZrhxB04GZj6Aq/gpXnbx6R2+xThkuU6I6ptVCzyRcR27SeHzRSq4ydUAh3oTwzdkIFCJkNZcHLXSFas8geIzedt31b8fu6/mVdncFwbk/zwZPQUBng9c9JiFyVe3y8fudXi86RdO7VXrlLZfc5E4sjG4ECiVxH2+/dv8RLG6fpsZUBPoZxClOQLncsz/82+huZeZxsSnEZ2es6ID374fQd78C5+RSjtvrteZDET22Ey5xdmgMTnSD4R3Gd3pNvfp9VrUFB7leWoWRDKaP5ERUAYtM8xXm2QmrB/rFmwOb2zeOG3w7fj866X4mxdwUSloLC0WBFNDwueZ0qzqRoOfbrwt4/tXmgl3iB7x6I+WHy4pVNDq5JUvXsp6BmprWd7NONTHHVtrfC0S4jl3G1IRAgoiwGm7HNvENqhGLHMPjDRgvuYarDgQRMdgYIhHWYxi/bQqowrEf6vvsHiA+zD9kc503vZbZRJ1NZqScO9qstSmAVsx8WOwmtrPCTdYUJyCdYcJqzmbJVThmNy1E02gwa2D0Siiix/qlbGnrzw44NJeL0c6u++7nHhmpsk0UMK5aPZCsnexcgaskkDat/dhOlpjJ+tDFTH03XdfOBuCQUhP+2okUufWX6OYX7pOkddcbaT70Hneekrdj9n876jg64A5kEBQG5H7pri1+G1r54D64fd9+u/LlHuvTB+Az1DIWX+kjomiVQEQ98r3d3WMz/ocdECQXN53RYvmxl0by73cbZsn1O+abnAKr7SvfUpguOPyMtn/Ti5oMlUfMdIVz4CpLZvNmiRVhdp/45ObNXcO0aOIQRFHO+4hIb/TnzcIWyVi0ZMbuReRoZNbK5MdTaBB6wSuqXNbcohqDD8Mh1d7pNYaiFR3PhSUv8PhQTWhD78dSBCUQL3/0i3fGiJH5C/8e2NC84jcSdppEnJtLHa7MkQVR1IrxTcF4ivC5X4K4cjxFck73qDz+w+UYqIazaeEe7RWKklFosiqf7Owd2nkXC1uA1JgVdKTqrBG9ZFoaRK+h4H1vqCBbQwmpEKVra5vj5skkDDX5zN+C82g1ZItHjPrs6ewKgkc1wjkjWI/3pSScWcSdEV/7R+LGnT6DfolJmvcMlJWqzhB7G1IbWY6L3PrjjvXXuRIu98fqu7mHyu6nf/pn1d8hYyjUCXi602YqGFLHpvnCJSgDgALY+pJDGTVpXApdemfLmo+XBWdoApw0eINKr8xMx747EVHxXhngDnlqcF5Wz6+miNIbObY4wG7p/v3z588Ig9UPQX6SHagrS4I/D4JfIJ0VT8azyGsxqMnT9XZN9uJEQ6s9ceRJOYGpvX7rZZ4n39cD31XHUuZq/bClcyuNEH6eoq4m7xl8dJ14htA2F+va22Slpmsdp6bJcJ55ZGkgBAFTU98Csy3GE6Qupi+kC3o7zE2dDVbGlo5/LHyYkEYB1ZxG6nHYm4W8Sqi5t860dz822F3SONLa240DUhmnEoVkfEqEShLB+iEO1PiU2gqdEGqI/V/Yb89mjg3JI9odbzEdkSWxdc/ByPjQMh9YCIT3Il/VobUmdZ8oR7DQ9XMCtUqRl/zg+ZFjKe9GkrfAdTN92tXCg5k8qHPZ3l5TMc3+Efv8JnWYAy34SpCwk2DESTXzwCtE1PR2ELHylosTsOJ3NMLSqNGFl9PpZVVOxLWUtkvYEKKfanY/m4EJC3cuAfsUekBXuapHe8IVvqqQaozmFEmxw6TnlXpzSUAzrJrsY+R0rr7lMLPTy0ty86iln7Zp496t8/dhvSibLIxlMsOGwv+fXk4FfsDg8NfXmLks5F7bkA2uxguFruJUdX3SDHfcHE1GFwFuaRUsDLyyiWNbA63RaiJAIjamVOum8Xh8Sb7+Od2p3PGFztPtpF+gs/3i4eVpv3gEHLM+zv5+sf/qh6ePJUZ3imDnvvTdLe9Lxpn3SR0vpGVsmVMEPUmqlH4BcLmXIqL1C5/sx6fFfjqbve2bxZuRwo8k9NkapKOXHjiSCvUBs6D69Lry8JFmkg+5Y/vi0fywXDzVxNnKQL5DTkrsA4nxEAf4RwLNIunCrAL59d1ovrS/3UhgNupH6euSFHRu+AilYh4qhWv+Fon8fMRmXw3Eln6OUKlDFSOjR123DAMFNpZ8uW/7/mLTDwgrfbvP9BNAs36UabhvcraV1VQC+2IQ4IpJ7tGno3fclpPR2Wh68pTJz/tZKpx+uOX4P7/hc51eSXTeL05Gx+MLsBOJZX/tHlZ9zZu1L2mjogD6OBq+Fg+eJMcLAQuOTJ8//f2z7w9fP3/5jSPl/d8/e/bg5Q+HL5+/3ofznySMVudLFQALmEp/9bGHTfh/AVBLAQIUAxQAAAAIAAAAIVzh/p+pMS4AAJa9AAAJAAAAAAAAAAAAAACkgQAAAABzdGVwNC5jc3NQSwECFAMUAAAACAAAACFcHzcS2GQeAADEegAACgAAAAAAAAAAAAAApIFYLgAAc3RlcDQuaHRtbFBLAQIUAxQAAAAIAAAAIVw4L8JlfawBALaVBQAIAAAAAAAAAAAAAACkgeRMAABzdGVwNC5qc1BLBQYAAAAAAwADAKUAAACH+QEAAAA="]
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

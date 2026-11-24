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
        g.drawString(CalciumTheme.fit(ruleLine(),fm,w),l,getHeight()-21)
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
    /** The selected cell's call line; it wraps at the width of the Traces tab, so a marked cell's line is read whole. */
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
    /** A mark refused, in the words of CalciumReview; shown under the form until the cell changes. */
    def markProblem=new CalciumText('&nbsp;',CalciumTheme.WARN)
    /** The cell whose mark form is open, or null; the cell the refusal under the form is about. */
    String formCell, problemCell
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
        // the call line, then the cell's flags, then its mark controls; a hidden row takes no room (CalciumRowLayout, BoxLayout)
        reviewBox.setLayout(new javax.swing.BoxLayout(reviewBox,javax.swing.BoxLayout.Y_AXIS));white(reviewBox);reviewBox.setBorder(javax.swing.BorderFactory.createEmptyBorder(0,8,6,8))
        [markButton,removeMark,reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{markRow.add(it)}
        [cellStatus,flagRow,markRow,markProblem].each{reviewBox.add(it)}
        def head=new javax.swing.JPanel(new java.awt.BorderLayout());head.setBackground(Color.white);head.add(pickers,java.awt.BorderLayout.NORTH);head.add(reviewBox,java.awt.BorderLayout.SOUTH)
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
        // left and right arrows step through the listed cells, unless a text box or list is using them, or a mark form is open
        def keys=getRootPane()
        keys.registerKeyboardAction({e->if(previousCell.enabled&&formCell==null)step(-1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_LEFT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        keys.registerKeyboardAction({e->if(nextCell.enabled&&formCell==null)step(1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_RIGHT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        // the mark controls: one button opens the form in place; in it Tab moves on, Enter saves (or presses the button that has
        // the focus) and Esc closes it, keeping the cell selected
        markButton.addActionListener({if(formCell!=null)closeForm() else openForm(cells.selectedItem as String)} as java.awt.event.ActionListener)
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
     *  review marks as the Step 4 page's does; then what the checks of the calls could not do (checks.note). null when the
     *  analysis shown was not run here. Writes nothing. */
    String finished() {
        Map parts=finish;if(parts==null||results?.name!=parts.folder)return null
        String carry=CalciumReview.note(CalciumReview.current(carried,results,reviewMarks)),checked=checks?.note?' '+checks.note:''
        String sentence='';if(!parts.empty)try{sentence=CalciumAnalysisSettings.percentSentence(CalciumAnalysisSettings.tally(results,measurement))}catch(Throwable ignored){}   // as percentResponding: '' when the tally cannot be made
        (parts.empty?parts.empty as String:(parts.saved as String)+sentence)+carry+checked
    }
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
        // with review marks the tile gives % responding after review, with the count before review beside it (CalciumReview)
        Map review=CalciumReview.read(dir,measurement);reviewMarks=review.marks as Map
        List after=CalciumReview.headline(CalciumAnalysisSettings.tally(dir,measurement,null,review),cell('Stimulus'))
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
        reviewMarks=CalciumReview.read(dir,measurement).marks as Map   // read again by summarize, which also shows the tile after review
        closeForm();automatic=CalciumReview.scored(dir)
        // the checks of this analysis: what its run here gave (as the Step 4 page takes them), else the table saved beside it
        Map run=checkedRun;checks=CalciumCallChecks.answer(run?.folder==dir.name?run.found as Map:CalciumCallChecks.saved(dir,measurement))
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
        plot.caption=name+(traces.containsKey(name)?' · '+(callOf(name)?:'')+' · '+(call.Reason?:''):' · unavailable: inspect Cell quality for the exclusion reason')
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
            // when the mark was made, and in which analysis when it was carried over from an earlier one
            cellStatus.setToolTipText(marked?'<html>'+esc('Marked'+(mark.time?' '+mark.time:'')+(mark.carriedFrom?' in '+mark.carriedFrom:''))+'<br>'+esc(tip)+'</html>':tip)
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
        flagRow.visible=flagRow.componentCount>0
        // one button: mark the cell the other way, or remove its mark; none while an analysis runs
        boolean can=!busy&&results!=null&&evaluable,marked=can&&(reviewMarks?.get(name) as Map)?.automatic==auto
        if(formCell!=null&&(formCell!=name||!can||marked))closeForm()
        if(problemCell!=name)say(null,null)
        markButton.text=auto=='Responded'?'Mark as not a responder':'Mark as a responder'
        markButton.visible=can&&!marked;removeMark.visible=marked;markRow.visible=can
        reviewBox.revalidate();reviewBox.repaint()
    }
    /** A flag's link: the Traces tab, the recording rather than the still colour picture, at the first frame of the event. */
    void goToFrame(int first){tabs.selectedIndex=1;if(movie!=null&&view.selectedIndex==0)view.selectedIndex=2;stopPlaying();showFrame(first)}
    /** A refusal under the mark controls, as CalciumReview words it; null takes it away. */
    void say(String text,String cell){problemCell=text==null?null:cell;markProblem.say(text==null?'&nbsp;':CalciumOutlineProgress.escape(text));markProblem.visible=text!=null;reviewBox.revalidate();reviewBox.repaint()}
    /** The form that sets a mark, in place: the reason (chosen already for a responder with a flag), an optional note, Save and
     *  Cancel. The reviewed call is the other one of the two, so it is not asked for. */
    void openForm(String name) {
        if(name==null||busy||results==null||!((automatic[name] as String) in CalciumReview.CALLS))return
        Map cell=((checks?.cells?:[:]) as Map)[name] as Map
        String chosen=automatic[name]=='Responded'?((cell?.flags?:[]) as List).collect{REASON_OF_RULE[it.rule]}.find{it}:null
        formCell=name;reasonBox.selectedIndex=chosen?CalciumReview.REASONS.indexOf(chosen):-1;noteField.text=''
        [reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{it.visible=true};CalciumTheme.choose(markButton,true);say(null,null)
        (chosen?noteField:reasonBox).requestFocusInWindow()
    }
    void closeForm(){formCell=null;[reasonLabel,reasonBox,noteLabel,noteField,saveMark,cancelMark].each{it.visible=false};CalciumTheme.choose(markButton,false);say(null,null)}
    /** Saves the mark of the open form (CalciumReview.change): checked first with CalciumReview's own rules and words, which
     *  then decide again. A refusal stays under the form, and nothing is written. */
    void saveReviewMark() {
        String name=formCell;if(name==null)return
        if(busy){say('Wait for the analysis to finish.',name);return}
        String auto=automatic[name] as String,call=CalciumReview.CALLS.find{it!=auto} as String,reason=reasonBox.selectedItem as String,note=CalciumReview.tidy(noteField.text?:'')
        String why=CalciumReview.refused(name,call,reason,note,automatic,CalciumReview.shown(CalciumReview.order(results,measurement),name))
        if(why){say(why,name);return}
        try{CalciumReview.change(results,measurement,[cell:name,call:call,reason:reason,note:note],traces)}catch(Throwable t){say(t.message?:t.toString(),name);return}
        closeForm();afterReview(name)
    }
    void removeReviewMark(String name) {
        if(name==null||results==null)return
        if(busy){say('Wait for the analysis to finish.',name);return}
        try{CalciumReview.change(results,measurement,[cell:name,remove:true],traces)}catch(Throwable t){say(t.message?:t.toString(),name);return}
        afterReview(name)
    }
    /** After a mark is set or removed, without loading the analysis again: the tiles, the Show list, the picture, the heatmap and
     *  the plot follow the marks, the cell stays selected, the time cursor stays where it is, and the finish message follows. */
    void afterReview(String name) {
        reviewMarks=CalciumReview.read(results,measurement).marks as Map
        try{summarize(results)}catch(Throwable ignored){}
        selectWindow()
        if((cells.getModel() as javax.swing.DefaultComboBoxModel).getIndexOf(name)<0)showGroup(0)
        cells.setSelectedItem(name)
        String done=finished();if(done!=null)status.text=done
    }
    void exportPlots(){
        if(results==null)return
        changed()
        try{
            ['trace':plot,'heatmap':heatmap].each{name,panel->int w=Math.max(650,panel.width),h=Math.max(260,panel.height);def image=new java.awt.image.BufferedImage(w,h,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics();try{panel.paint(g)}finally{g.dispose()};javax.imageio.ImageIO.write(image,'png',new File(results,name+'.png'))}
            Map provenance=[result_folder:results.name,cell:cells.selectedItem,stimulus:stimuli.selectedItem,protocol_hash:appliedManifest.protocol_hash,manifest_sha256:CalciumInstaller.digest(new File(results,'analysis_manifest.csv').bytes),settings_pending:fingerprint()!=appliedFingerprint,source:'Completed analysis snapshot; display is unsmoothed']
            // the groups, colours and order follow the review marks: say so, and which file they come from
            File review=CalciumReview.file(measurement,results)
            if(reviewMarks&&review?.isFile())provenance+=[calls:'after review',review_file:review.name,review_sha256:CalciumInstaller.digest(review.bytes)]
            new File(results,'plot_provenance.json').setText(CalciumResponseProtocol.json(provenance)+'\n','UTF-8')
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
        if(saved.problem)throw new IllegalStateException((saved.problem as String)+' Move it out of the folder to review this analysis again.')
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
    /** What carry gave for this analysis while it still holds: null for another analysis, and once every mark it carried over
     *  is gone. */
    static Map current(Map carried,File dir,Map marks){carried!=null&&carried.folder==dir.name&&!((carried.kept as int)>0&&(marks==null||marks.isEmpty()))?carried.findAll{it.key!='folder'}:null}
    /** What the Step 4 page is given about the review of the analysis shown. @param carried what carry gave when this analysis was run */
    static Map state(File dir,Map measurement,Map review,Map tally,Map carried=null) {
        List names=order(dir,measurement);Map marks=(review.marks?:[:]) as Map
        File source=review.file?file(measurement,dir):null
        [folder:dir.name,file:review.file,sha256:source?.isFile()?CalciumInstaller.digest(source.bytes):null,problem:review.problem?:'',
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
                    synchronized(this){finish=[folder:output.name,empty:empty,saved:'Complete. Protocol, input snapshots and tables saved in '+output.name+'.',noise:CalciumResponseProtocol.noiseNote(submitted?.protocol)];message=finished();tone=empty?'warn':'good'}}
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
     *  marks; the part about carried marks goes once every carried mark is gone. */
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
    static final String SHA256 = "79917331038f77b575106919102346c1f50908d5f3e14897ae78c0d40d5c9d90"
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

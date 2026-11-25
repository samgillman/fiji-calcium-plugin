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
    static final String DATA = "UEsDBBQAAAAIAAAAIVyskSrinQEAAHICAAAOAAAAcHJvdmVuYW5jZS50c3ZNksuO2zAMRdfxV/gHmugtq9vBLAoUXXS6FyiSigX4BcvBoH9fJWmA2Ukgee7lBWuZt4kR4s7bWsux7n9P43Fs9fvlci3HeEtnXOdLhflapmmG5VL/T8C2da93bD1zOU4pKDLIYSCCpDMyoLZK+uCFRpnNkIUxGESHIyxXpliWeIylxm26Xcty+jNy/zmuE387dkDul7VU7rkeZYZmrZ9gv/Ie1xwTVJ7KwrFShIXioz8++9u/37lu63IfXhqYe3WWZ9F/8t7KREx9Wfq7cv9UfszAflc87qWPx2Jv0P388fb+6+M91hGUdSdvgiet1JBRGRKYh6R88lZqpOxCSloBJWswJ0jOgfeGEgthuUVhJHRPX8R7Pf9+Qa3yTMI65Yh8pgbIRtuUExvTkI5BKJczSCExgNHCW+MDWTGowJm5ux1l+sJjgckCgZaKQkBMCpsXlNaR1lqabKBRyZshoHbWCQIUQrNMGO7Lda/w4ravx4rrVzhoA6BZuCEP3qY0IFpGq7XITc0zUAuFpc1ZBwE6awwmGTvIdhHAKXT/AFBLAwQUAAAACAAAACFcyJJObBgSAADsMgAACAAAAHJ1bm5lci5SrVvrc9w2kv8+fwVUUZbkhZpIuU3d1cq8KkWWtrwny17Jzn7QulgQiZnBmkNSBKlHNv7f79cNgI+ZkR97myo7wwbQaPS7G/B34lz/QwuZy7pVzVxcZ1qVrV7oTGSyyLpCtroqjciqtRKLplqLdqVErctS5eJar+tCncqffhCm6ppMmflMNksjXhzQgrUs8xN8hm0jdaHL5ZuyeBKJeHf1/iya6YUIC1Uu21VIayKxl4g/XkTCtFUdBmePtcpa7IEZulQi1w2+q+YpFrqsu1acXv8aC6PaFniNeEdfVdfSSD81iGaWrHChCzWvpdvq5uhDLIKuNm2j5Dqg360uzPwqiL56RaNMDb6otG6qtsqqf2l5rhq/zB8EnMOcfJ6rQq/t4v/E4hVgqgHvzk8urs9icLeYl3KtDEBZGHxUT4T1XhadCqJ4Jnb/h1WnhTSGlwXZSjYSPG6w9K6rWkXAIGbJQQfmNO5BK5V97DdkGqwAZfn0sqsLnUnIKvSn2Ac9kZdkP95Li+RCB8X3JaEMpZkX2rTDej5INMiXEc7Kbk3LFl2ZkVKGBBT/FPeMq6trMNX8TTYlrSCcmK8anYXm5gZTP3yIomMx0rr7aC85uhC//y72tJkvdKlbBaCju5YGn8GrErToPIj5SMfY7JMX1iYtscjVQnZFmwRBxBtpoqEoBgJoMwLKEU1+lVCFAY8cfAa9WtKBJrvUsWjVY0unBkdDUze6bBdhcP3q9duLs9OTv7ff5/hj/l5CaH4yqF4UnVnNMyhsVagwwhlgJPMMmoYzspL98UMMq+kao+9VQgYam1X14LmZOJE/NODRhKLHmLQiEjwyz8w9IBvKD9Q8J26qB6tDFh2gxKlZJsuqhIYUKf4U0wM/RuLxJgbDc71YhLwYsBgq/zdd5tVDWiv5MS0qWp2fnwexh5+8Px3AKVQuivOmqu3OH2bfiXdwY7CHbg3nVpFfIr/2sAJ/DuCuMiWUafVawo/YsxmxrASMp4X/EepeNU/CktvK20JBVSuhZLZyOEW1EFIYSf4xNbmdhF0/KlUboVsj6gKbzGdlpeFEHCEp4984/z9hzTLP4QoxAHHDZ8C9hcQCojZlalPGQ8ffhKWWCIOxS/vNHgpewnMzwgZgchZusZm37f835t+nmWWu9cXbRgmUbMPOzMlceARwMoy98jfyLqG1c7KKW2kU7etAe4kbGttiHxNMKWuoZ0u/SdlooVhUjbNSbELax2RtqaLF6+ngcfUI52NCmrZh/a+1MSB9iheH98YZ/hwHV3DN5A0WsPQD+KBMG/BBrJU0XUMEFl2FuZkqieuzvO19PFkLk/XTh3jkYXvrmCNesPVB1iS+E/pLXkJ0buDEnEuKdL190plKKFOYt9GLn9i3lbA6/jw6tN8sXAAQlPaS4LzBd0ADU1feT4u2IjIvwfRcdKW+g4wzVRS9LT1ocF62olCky0eHSBoIVWCJ29M5pRdQnbGDztsbUPMBmjZy2uouhWRDT380UMIEGLHusMEtWTHicIYQfk9m2MqGnTNIOAr6IBXu3cu6Lp5oqwNE47GB6WHTx+j332HUofXRMI24qJZM7dFo/8uqPHArJtIV8BXwGXAcMDk68UI31p4xOwxIxZEFwQKJxgBu2Zv7dFiVuSOc1+8lbQNq7Qfoo0UeRr8B4qEXR3bsBX/Z3/+T9OxztLt4JsqqWePHb5zeCb+3aGS5VG53b7cBG3m6rnKMQGVuZfZx2VRdmcNvNQ0rRTrVcrcZ53tdWegF6c2w7qBfN+WfNkJmmUIm6jkQvAM/4b6/1+X3g0qyfyFOE/Pydt9OcjbN0h4COk2LSKwA0+8Xh+6DXR2PvkgOe/5YVFa3oFithLe3mGLwDFnvUnpFg/hbnE+XFEfJTcDTYUlOkrfkOV9LVAbqkWxLt+lIQ2afbNCn4zTSRlbWBfbZwJiufoOe6EVIo0TmyADsCm8EdWU0EdbvTbjArH2ee3AU/UjTjzepIttpZdmmCyZnRqbsc7NRLADTKeZaZxF7oUTOpl0+xWtHRgInt4RjVeQExy4CJNoRhFDsdPn+4sI6LhsTerVrKp36iSk56iByoh8v7z3pKCDtXLvDxT7nRoUolXIRF8c+BfE49asyKzpAsUfqvDvlyRg5e8QIef2UNKFC8hX8MhhJrR8h5aCPOeRgLP5Bqz2tTlNHjtiP7BMVUTS1Ym8tB1dvXvVc4Z36Zc9QTUdrq6J6UA154r4WCL+0MIqIanAFHoh04Qh/ntToeJCfuuvgMieU33wJL0Kg1R97wvGBRF7B3ZdVK5CPIcXqQ+uGVn2aea9mRZfLVs45/ISksQkb+4xo0y3rzlqREbO+33PJGAZUhwTxVKNnSC1CTXkfBSVZVNBOr+xWHQ3l72S8fv/U6GUJdeREbDPU8dob/YEinvPrMXvyOHjOoWKXO6vttNW+OwOznJ0dlxmDEjFZIMwusdCbiTySxFMRTxI7WuZwOgfa7G/pMpR0B5T8k9vY82TfptLOyo/v9q9b2XbMcGs00IMA4Cs2GwZfkBmPYgXih7yXuqAEmnI85OkIE1CPwO10t//XU3BL+uLZFdbpXWah4TABZnly+r9/vnrz/vJl+v7y5NeTVxcnv1ycMX8FqiL7B+cfWDocA2xFMKVN1Lpun1KnPSlF2Ghw5GLQgnws65H4GJ2VkEXoWy2oBDDfY57sHpP2bqgLKXiqDZU4h4ntqtgTsK5ktyiFQhZ2L+u7z/i878RbbKiae/BZrnWhUf4XLrnts7sVfClXMQSlkonSXmHgPmBz3rSgwlWTzz0hd1RceC96fsje0eU61y/TRj4AYhUj4HSaPeg0AN1F344i2qj5hKfvhuTBhAHmeN3DLNfu4tcMp1wNSBOrEt4iqof9t1R3Et+jKDjx+hlYLQjelwYMYAgVvVw4+Irh6PCHRVHBofz883/oH6eOJu4L+kunQsx78nLf51RSfp8Hsd7wTuT37jL2d9WcStKwoSPE7rDR8cjLTWe4gWgoZ/7rMA5+pcgiWy/fXGUFC9d3vITveHGKthW1/eg0Ynvo1ET1uq6atm+hTUP4JqJxsvSdOKOazeZctkU0b809l3sGG3DaRtTDV0DJKYN3xbhoOvylHiXlbvNnKetJMvvUBugQ3mObmbkvn8RvQDl3t/03HujbhEPy73NEu5ITvRZ+MGWt8GBbugNnqzGumgnSUcHTV/dRTB42ppS1qYpkKLPM/sesAF2kjXlso7ZD5ibbo7CKf35VxLvTgDs8y4P3neADE74JGzOtx0X6bNszzAQktysoK0vz+qUYWELVQg0yyC5o0DkdXlkgTc+erKxXkho8o35OzNOdOgClXY89+yooHBTBCoIQRXNx7bu07QrlZWnL4AqLm1G3CGNW93o7wf4PJafFltCunM+G+dTVWYTbqpfeq8Zw088BoiQ5uoiCvqnkvI0f3ocj7LGinN5lnCYfpkACf/jDMyN7Sf9BzsuVRquhhWsP3/ZCGo7DnnsXzjgQt1074YzoqBUdxKPNdHmvjYaejDiC0mniJvyPHY3G3mvM/+FCwIzqfm0Twh1MRkYBB6OKfMA6FCCsUdT/ib6OMLfVLrrcULpN33AVMCWRqUptCmDCIaewefJOKuP+DJxGuAiQ3GXsqb9axWxK64ja75uyw3prXqksofz2fuj/SV/sMNrTJsPuEczy1PeHXay2LSwEJpM1+pZupST3cjAB2l4ubZpoW7R0LQARQfzwIVRRVc18tutwG03ocHNOtLVq/xRfmtsnoyYnOWcLReG+teQXXbmGUtohlZgsvJ2OBds72gpjs1/8BVKRQ7ffSqtd823EcpM8dClFHLgfVJkDuR0c6Ugw/LY10nQmNCIYq4gddHtsnjcebtOs8Caots41mc4Qj93eHPhU19b+52SFyCld6nlVUU6HVNPG+12Ja3BW6KXm3C/o2U4J6wbLRiksh9PUjEGIiakZ7ZRSzG8dyM8+c5M4GXV9KNr2TT8XBKBgWLORMuCqK0fb0ZffiX7nXbMxtb88uPKZDJUpJWy133gM2L6V8dY5vZrx0PEigg+kjYCewhEoq+6pr7QTLJdqvCvBoY2risToK8FnlekmHq5YpgoR+07NhrVNiow4cIsm+j9Dbvew0ZFgxUp2eskT19VJUeX60mOjwJ+UAJwGcKO3n0AJ2VtqrPT9FQQrW2/vRnh44e4fkZTtfbEHFF/1rT2OoBLO98lok0zoIiIux33mlLtWCe1xl7lewF4ydAKQxE6m9yRPViSTFdeqUFaVkL6BgiH7mA7w/VfrskTlsnRAVBJcjfqU7BAECp5pl8kcC5/jH5haZfxAQtFdkrzVXPMiG0XN8ZxrYV0jLaAWz4qSBvpwo9AVSZrx9RpypVx/MHneI8d/bqquTnY4+5gKVHsIPnIfBD5zr9Rf2Xli9xl9FE2/uYWzG/sQIXh+PBwo9TDHj8adzoWrgad+AB6zW69l87TDu3tcjmhKy0b49nsUMNzd8A0pWC4O/dEx88dQ76gnwLdV7V7QTMAsWVg2PM3SobTHyEewM6tdhRrB3qomI1t2xwbvEt9YjJ45MHXR6o3DWtgXD/otR+qZghp4B+GvlSxH4qvtOQbiY+rifGaYwm9aV4i++fbo4Mt3yZMUzSMeG+Xx9iLLlzjg/++y4C3rHK/uSmmo7Us+kQUc25SAJvF9CCKDqhHMqP0+mOsBJ1evXnKuSjVC337+LHLqbgUj6Ib9e7r71s5//4TMQd7z1airiFxtwVe5Y9pyLZdlBX+Xmd09HqNKvm2C35u2eZaNzp+/lNle9k13MgNZU68wRvsvlktEd0xlkhA+Vd1KszcrhIGcPhGYnHA7FR2vcKnoZMFY3z7t4ru/fUx9p2XCfJt5mOf5/8zyf4sQPMqR7nyrEBz9/w45bJ30q4SxvWoqEZSi1Iq9PntN0d7WnzUtmItfNszZxietbCphVC35tpYsraTKFPa2wN9sqfZZpbvLwobzma2IiMHSzJG8N/pxxEx6wBAd2140d4+Gu25eFx3bsFk9XHdrE9qxaCbvl67XxDHCzY0wjZyzYW9Dl1V2IN5zt96XJ3SpWaT0XmjerO01A6eK5Nb8qHdZ9rZ9ZtR61172CYbb4Gj8CMPk4ePNcJDH6EP0o7lr2pAyvzEc/31p8y1X3V/+cRhKwIgYMkRmtI4vUxKpScqISlbEKCZu8J33qpSlvcEaPcvcft0ZDHOpEQ21ti82nT31zzUT/1izf6s5vMlMxi8y+UEmP49Tj6DJVqNXvmGTTK5sl6q98p0caD91lpCSqt2Ta5l9RIX0q5sf+Nn0Qm7UDk8NvMJa7sbxNZ0kxvade77rhwXdi6jc3uE1XXk87LcxURvuFCrZIFSZ1j8DbhS9K1K2AZsrqIUym81WatZSZ3+KcWgqbWwVDi4ofXbJMzOG49Jpx/3OoYEaDy3kEcIBOMLR77KSZrVrd4KP5m82FT+/zs3e2CaVxbKCxazWSfD65c/UDe897NCo/cv1m0vhOvrtU8yX8ZLeu3QNFT98292iWoqCePr2JRn9jiePjpIgPD84P4x+PD88Fgv9SC/LDqnq4seTtnrjV0OuKe/ekcWjzhG7zyQYiv3h2jgJTgqyWGoB3rast7moyJRvn/jJ+zHO4B7s9FMYByWtFJSfUnq3lNBrL8WnokYCNwdguklwYa8YAPyTfTYq+qAWC5uTwdW3/GJz2cgiiPnVTFul9upQsUtjPjhU1DT504BleEpm2e3rT+H6LUzTtNuQuNtzegntZzcKgYYeERyL6pZvd3N398EFK79Z3b7js3iZR+x/6srm/LzpX0+RiiMzfUoC/0DY3bmpnJfIrq2o55QJFLmQ8UFbHdBFS8NnPKi6li9wlH82E8RQrGqtS7acoK+GxI8iHD5+EC/54Vrr047oWLQqW1lV7ZEZ9zTRX8E6jpo+CvMZRrl1ws09G6nhTlhRmFO2tX28lasPwb1Rd52mAO9Sd6Tx0NBRDkStiroAIfTvEcQo23OpgE04DlzC4f4FA3KLyp2TQ/Jmwo6ggHQmnfybAb619pB4iEbxKA4i8CS2icUyjexb2KQr+eGthUXbuMffN3v58DRpPEAv9GNKZC6Qo9kkYfqwxEbl8ZI48G2jdC1LvQDrx8XWBfl2aluQrpLUdtyqnL6hd+/vzuiCZnTDfhgHJw616NdHs/8DUEsDBBQAAAAIAAAAIVx6uMRcdQIAACwEAAAQAAAAdXBzdHJlYW0vTElDRU5TRV1SS2/bMAy+61cQObWA0Q0FdtlNtZVEmF+QlWY5KrYSa3CswFIW9N+PdNJ2GxDAEMnvRaaQGnLX2jFYxlJ/fpvcsY/w0D7C89fnb9CYE6zcMJzMyFhtp5MLwfkRXIDeTnb/BsfJjNF2CRwma8EfoO3NdLQJRA9mfIOznQIC/D4aN7rxCAZa1GE4GXukCf4Qr2ayONyBCcG3ziAfdL69nOwYTSS9gxtsgIfYW1g0d8TicRbprBmYG4F67y24utj7S4TJhji5ljgScGM7XDry8N4e3MndFQg+hw8MSS8BE5DPBE6+cwf62jnW+bIfXOgT6BxR7y8Ri4GK8xYTyvHFTxDsMDBkcOh7zvrpbp4h62daaLyvKFDl2vvTv0lcYIfLNKKknTGdx5XNir9sG6lC4wc/DP5K0Vo/do4She+MaWyZvf9t5yy3244+otWbBTrA+fOq91bozTDA3t4Xhrq4XvNXnInkQ8TDOzPA2U+z3v8xn1B/LaCplnrLlQDZQK2qV5mJDBa8wfciga3U62qjAScUL/UOqiXwcgc/ZJklIH7WSjQNVIrJos6lwJos03yTyXIFL4grK/wDy0JqJNUVkOCdSoqGyAqh0jU++YvMpd4lbCl1SZzLSgGHmist003OFdQbVVeNQPkMaUtZLhWqiEKU+glVsQbiFR/QrHmekxTjG3SvyB+kVb1TcrXWsK7yTGDxRaAz/pKLmxSGSnMuiwQyXvCVmFEVsihGYzd3sF0LKpEex1+qZVVSjLQqtcJngimV/oBuZSMS4Eo2tJClqoqE0ToRUc0kiCvFjYVWDf9cBEfovWnEByFkgufI1RCYIr4PP7E/UEsDBBQAAAAIAAAAIVy7atflJigAAAydAAAVAAAAdXBzdHJlYW0vcmVzcG9uZGVycy5S5T1rc9y4kd/1K6ixkhvmRrKl3KVy65WrtLacqCLLjmTv1t2cbkLNUBJvZ8gxydEjl/z36xeABgiOJNu5p6uSHQENoNFoNBrdjeaz5ENeb0/z+TyZZvN5k6yaPGmv86TJFnlSVvUimxd/zmfJbTGDijqb5k2SNcnHAqpfV6sawLNylvw+z9pFttzZqPNmWZWzvJ4s66qtptU8+X47uVyV07aoymHTFovVfNWMkqbN6naU5OVslFxkTT4vynyUXM6rqh7BwEWDf9WABcBeZcvRRrL+37Qq2xpG20/eHhyfHY5MwYQGguKTT8fHrhTGlbI0+Q/ou1wtLvLaw/UuTYpmByryupjiX7/8ZTLPy6v2Gv/Y3092j7EIYC6LsmhzKIWOzAyxq7YuFreNnTRWF5fJUDqxxckm9fWXv9B4maqAos3yz9PrrHaFKVGiaavlcHBYtoB05gYtgWA7AzPQJizp8CZbLuf3w3nRtENH9HQkUx4l8+qqgMUf7qYpjYik+X6fF6g7GDGHGY5piwyAjYoyaXIg8KyBlasWBFpDQT0ryqtkOq+mP+8kH6EQgRerpk0uq/m8upU+oas1qK/nkf7pCN+Y5jizFzhL6sX+RZ0l30sVdQl/7aq/YI3aGlhjKAO6rjeHxCaAL/ApUXDIuPDfwCLw3+RV8iL11u71dVXh7kmWVVO0xU3ucGQGGSXZYjkv2hXsPMYWCb3I7oD4C9oT0LisyjK/yqj92Rug6rwtoBESAqGz5Pa6mueGvatL3ACwSCtq0OAAecNUvy3KWXVLfDvL2myH5jk8xdb7yeBMlnwwSszPZD9x2/kMl2/CZcRjh+WM/kRekzUtmo+nnw6HsgVT3nhc1VbAB7nictxfg5+n80Eq5PoEtPrD6zkKH48Fs3mFTAHUyWUzzIrLy7yGv5Il7L+svrfAPNEeBvPkhSco1uwVLU6Al7w+UuEPf6viHPq2zaN2B6PPXWx7fK2RAa7zBaAPKqvi0DvLl1mdtXlAXOgJETa80VRJnk2vk2tYhaJtkuq2BCLnFt4MYTBVPFVfwO+hFIxiLPaa8fU5bADDDx4U/uqfY8RgQQ1D6oVFLP8K/4ODrWwu4bRzCP4gU6E2Tvi8xY04AVGWzSeAPp0UACAy6ezNxO1AKGYZ5WH/2u2/iew/bC8n3bvsbgI7mwblLeMkS5qcHEzqHAZO8jnsBizb+OuGOnIZ+Qmf5PoomyFDwyEPI1wW+RyOXAblLdjiYY7bvh2PB3iyD87P6SRDvoFyht0SykINcbstJ8JC6Ty/JHDDb1LtCAkwz5KjcjpHiYaDAuRiifs5a2XX8sGRJRfVqpzh1q3K+X1ycQ9ib3qNvFtjBRwnO4g3Kxe/TX6V7Lzj+q1ZtbqY5zs5dPsrlJXDXRCEF81wOkT0tNrB0tieTxOZVJsV8yEiN4IzGSGqm7wGwY91CLLtNwGA6bxYLnNqbGBfEW6wB23jfUMOt8J/j0BW8G79JIvnepOfHZDbrC7xRAUQ5BEBM+K0AalXtpfDAUMT0qZTou4vrmBg/H8U5qR/eTOiPcEcNhiI5LZDCJFCEoAMIIa4LqbXRLzk1T7zwzZTQnjM8gYU4jizVS0T8ZsqMNuWhBYTjQ6U2R3RaAhjj6QjXtCsqUpDG2w63j2HBSFsqDXKbb3KKeDFK4tdh1Q8tQrMrIKtWlYt1GUlyclL2JpOrKKw5L0IpOMVA00gz0EdciR/qTUi5A5dmQaSjlmWlkgm0r9WSrEEiuCcfn1sJjN4m9+S6kZ413keFdpGHRhEOxUKQ7/CjBFBZsf7WFXJZX7ryGEkHchYopzSQra5LqlX8zwYGvcvygW7e4kKsPLnsoOlhHED8pwnIn5SWPH4hrM4nvGwqJvld1NaJlpTlI/JvFgUrcJG2gzkvDCSD0TWlAiQ0+EiTOekQB8KviBfFOUwvg+l/VvEaeI4m3+hbFq1gb5GoHigsJh/jTe7fSP7uTtS1lo+XT+VDQyZgcjsOWRPedB9GVOdi4s8w3IzlXhzCw26qYIFPQdYAzTT2aQop9h17lV/yLOfozXRQfhMbk2XiZW1scN6lJzgoeyBr+/e4drAhXcq+OClDHqEUxdP7fqq8IiRvAf+ZjVkfec0Uw8sOV2VrLNEynNRYtZ3Wle3OyXdWuy19zqHy5cpoxsy8hdux2HBGujnCarSV0PiFauc3xvFgIrHxTmpBomnnNxbaOLILWYvgMW2g8M7OvBnAwXBXGUgjsqbbF6A8MxuQZFaVUDuaV5O1YVo+PZFatqX+V1LP/9qEZErMrMoIGNBsDIr74eb7pJ+T/IjDVD2EXpXNA0KaEBrBQQD8lhMgGCBjO/BC1tgZ7hNYFAQuOesEws7vsDKXRAICGiRdZVwCH6u22FEsVmL+jF24ND1yFk0yZ/zuoLbjLkxojQGUd0sQGPsmYfcjEmuts133zUzM5vkuZ2Lg1T7SmmOnlr8K4akNh40ivu1OzcYgelZC36AFp8D53LX6SAIMpmvIvN8SK1A0wg6uyzqpnX6CLbYEi74JbXf4uOwwbZrDkI4i6g/e7eHLjfZukNDuIucgXCLOs8YhWa1GOohx7hH4bf0cJ5qJnAyQzhBH4zU47aH/5g6OQf+2z2WLe16IinT048HLFJOQLuIKD5agqCzhN3BpcYl4DlgQxaEkUGxoQyKcCxrBQ67YABYa7uUFjaQ9NIo4CILHcp5AScBE7D2q3UnTJoMzLHHYw1EmybggR3PE5M+K4DWyXe6mWn7BuQj6p5y2XO9eNs/6OVsBUsBaMzUPpvWFQk30/EJbP91YHYg7+CX8azs8gDO3mhKR+hrTnaz2p2l8M7/zlqLEgaA/vX3UtSkBmTNQ/dfY53m8snnFZxA7b09L0mFK2b2BOXd+Sz5Q54vk+uiaasa7UBiNL+4b/NtOEu38QceF6BIFjVI1wJYFaBQx8LzZJm11ztJ8iNo/QXpb8+SPZTK1bLdhlZ4gShzGBmO66y8ylkhbW8rrSKzpWu2rOB6wvc67Lm6TG52d+S2Nvggk5vc8EiD5BdF+QsyDjdDM3N7jAgQK5ah9d504ZpZyWXa7e8ne8dGGNV5u6rLYc+iTG72Hr0WahGsncZAbv0OtOlFNcuZ6ZGThyxdLYTVueGSMDhYtdUCjrwp3HgH77ISRmBLMxrlQL6z1RNvGtoGkuoRbX/j3oGITXd3/hFOODksF/msgOPfDjN2A75KXvgDnIpaAbQi8+WUJjd4a4xMsMQkZ7y7mz7fcata08+cLZvmyChBJVQrP3Jbo5BTNGpEWrNY42KUwL27rpZGoTynZZpWq7L1Nl9N9rRltVzNMz7HnZVXlMy8sQceAG+Z+8m+FoSpKBl9kKGEFHgEx/7/HiqppP+yFDVsf7B4Q7GaRKh4n+ZXIBZg4bBDIjciCJS2M4Dye7TxnVSlpXZDRsJRcohKBl7D8O9O30aFRsQiM7caNozmLnQ9wOrGl3YG+pDDuVe2E8EPJQvbIeFQ233xAhgbSfkcqEqHh7mFmO1Ja09L+S5bDgMGcItOYFI6wBsqqdi0CAMF173bEE86FtwyCxY41ggsFbkkLGZOxQlxM1mOEGo83js/30IcxurvKP+JQicuBdfBrignPPnx+NfnJAd4kgI+lv/SSCyMfYQ6e2mdyRtt4ttWs6T2gwhVALPuwsyqHcR7SOb4kVSkdBSh78HreZaX1QIOMDjs0C/c2MnLqbdqE5BkFcAtc3T6wRYxXlgybWVXOR5JrE1NYL3ZYKeumrRQo/61YwssyNCrsgK+mPLCDQeB4QWlum8+VSViLcUSbUfp+Lx8woj8FI7TPuG78S4oCAqnXinI2JLfg/82FhhQ5S1NulK0n5nRorzJXW05eQMcZ7hc6vQ8odYairnWyhrgdTyGgM6oBaKfx8mlgbZ7aatqbIg9cydYwgUhHyJ/Gk74bjCKXtGlQ6uD5lYAMj+n9urA/MtsJROJ8ooh7V6HoiKsQI0zJg/YJqtFiVJnLTP1so6hN2pUiBTN+wVNfALQxpvJw4jX18gHbHNu12TTGFGcV1R+jcfc/PzcGUQFbxL6AdahUdHaKuGG1sIkiBnt4a35DwBI9QnUw5Fy+qYBRdHbQz5TEgT7Haki+1p2wX6wCQxGKbmvWPBMq8UClZyLBvSjfGY1M1KKUQ2GYxfUHjLUZ+gWBZiLe7hV1ttssTWq8A50eADshdqp7qfOUW/KER9rZm2Aadpr1ClkeCe2nD5ErqjwOsFBKd7lgFbTGMvbhnicdGECwYVYZmx/EK8z9aL6NEzV446LaasIKlYlNFgbpvK9hk/SSK1rkXl2aKJnrP+BA2d8LxmVTpCWE6aKXw26F6FHpxH50imi5tcsKe0QaBMgGo0HVAQb6Tw1A6KNXCqhgKsIxWBkH87WUIPAWzrLL7PITbGxMSXRZUXvkfxpwDZZAH3r5Q6Iz6OhfuICpnDrGT4QAEcyfxGixBfTFh/LX7UM7p/CTUGbUlkIM/RYWWa5hG9RRpZZQH/pQLE1WpilqKgyRFW2J7vTqqvXjt24vTAprkNbtbLjQCierRaNwI9HeHY49XCEx5hS7LHAKuYj39MCNea464jW1MhWGXhfMMAbVDWnOwAHzGDhWHV0blR0U6VUV7SOdRv4SrxdR3HsqIOWyAmShGzZ/EfQNmlmHPojJmIGAnZL9hO/E38hKaIwupzMuHTjJa99a66LV3W1WvJ9EaOjlDopUBYcVceoVscnHXek78L2Qur2BSEAa83A43GD2y88C5/AzXIjhYOXeElswgLMohyp3AjQO/wt1UDNbKdewMAYRqVgx6QRvOgO1eQLohM6EsyI0lnCnokm/be9oF/glWEJbAT/t53s2pi1fAHDfJ/sdUdRF2oKC2V5pbyMeOele23yjlmLRgaAw3dYmy9INhlt4C2Z4MkIdHlPti2jE9SrcicJfDOkHphQEjy34F5So6kZjd2oC7xZ1Xx7Lxq2rDtoEE9oQqHxbOmIXfxscGNjxvMEI6O0TuCwm0hoSMC3NcUmWsMmBSgQGmS1Un82j4xtMgFCC2ttIlbNF8uW9o0L72TvIB7rV+SclEhXXQQaJsx6Yawauq6DjYtR8rowE1C+yDQ4LUk5oRBWtyucu4WiMJz3B/jXnqA0qXTjaQE+NJomLHaQXVTse3EYkNcBnTmhA2hDeX/UVDfiXlIaTkcwbtJY4wLvYtGOXpJnTXkttTuIokYLUME4ppjjHNCRxb9wJx7jZW0RRBKkbqiCOmbOLDiACEfxGI9uDhJsCDeGYBjrEkK6A2kcLUnaC2sgMjwo+o4QjnnTzsui0HVpFfL37XUBu2vIxYnHLzB7piRVknsKyny2JGxV/bYrOUdWCqlkx7dtfGu1t3lobmbbYIPOfilie9btE/G2ORp5u0Vh6pE99W+8tAdYIP6EEbzbdFKKN5Y9DMmqpBjMS3YTYP13oOz/w2/3fpPcURkr+7g94N7Y5tDXLL8pxAgLVF614qlgwJHpy0TQYrD/Rd7e5nkZixseYYfFTTHjuxgdM3vpTpIctSjim1WNET/IHNtttU0/ZAIY1QdwBxQQBL2gsB6p0N1pbu6OFP93ew0IXFRwVTNUzmqM1MENzV6TCgRwDahjZ9kSF1KCnTjAh64fMCDHQenpEdif4CQwrts/JVc4ReiRrOyB4KdoatFarGc8kP33SkqPEtW1RFs+Vajxq4J8aW9997oE/6uUWfwTxZu7HkY2hL2/8xQ8XY87eLXvI25uj9lMAJ6b5e5RItV2UP0G2t/n6eRynl35N7CdnR2rxIsStioJW/o/rB6hObIbzeIaKHXMaVdAFbFCGd0sOLCMFmbEwoCsUGzTgr32eQXnF7/oEISatmYd1Chog5cDdJTckZpO2lSakh1kDponUXrwMhmsidqd3Ow9KnA36rZcG81L6vfEVKqHLVbsemL4lbE3ovaszm8+9qRcMRk+b9j48pjhTqywbWcujhO/gx9UMcF9STCw7rsbGixBreTqtMRLHxHk6qEsUaVswmNVYOiHvq4f4KFQ2G7PstkVsayLciOMj93UA5vgTReZTc8QLoyoRRIBx7QgkcU3jzNA+Ysmyww92C0QhJTs20r0djptwpjSx4TlHv/visuNWFv6AnPNqoCAtawSRNDaaarI3V/M1oftjtZ27CO06Uf12rU/qdx1y4hw0TYiYW+6QycjhP9Ip41H8GL4mL/70oi3L/5PepPNqgKFTExw6tTkJ0f8mrk+EO6nw+NZ0SHPBcZ4+X4OOWAcjL4xyGq96SrosqBKKD8q+rrjq/mqOOyHkJRbBCJJa1JjyMpQL0iays72WN6/oQUIv7d3/WDDwu2dIhWvK9TtqiRrmrzh9bPKtY0g31jzqIK0hQ5A8KTiMbM21yHL99LoNAdtG6YwMZZ1vKsDkaIkcKSJNWO2Pc1B8wAdZOZmCtQgShjRUFGo2YKDjXgZTTcvcY/g1KDNtMpqYOx7pC4NNlCMvD6mXakP5n66S3ctJcCxwhfgr/a1BMefnra48XCEu6WO0GBiGGOC0p7x4yFIIugWAmfkmWrQfdaQruFhQvmF69mqHut7RmHf7ZYiZ1WfcAdbrtpJaeJf8bZo3JsciL1B8dzuNB/X+ZLATNO04/DcWPcyAJuLLql6kIcC2rsqceAcMGHWTYe9o/eQWNRa+rBv39KthlBPCZZD85zA1VNPva8Lel8SqPLJ1by6yOad8PgfgrPQXUpwjUyEThOpF+4YJRRWPXn7Arqdc+h0OAhD2NitKGByVOJkHRDdGzyIsJsOBA32Vc8nvt1zCakDyq+avPuQ4ls9nMA3k9b+Enn+8LVvKKhUyeSwf+osbga1rIPu6XJ6H3RMmEWrqGM5fJYEFT4pdb0I2MGn17Y+0tmz5CO+Yi2wGjbBwYcjuKMVwACg2y7hrCAVJb/Dny7NAxxsdIPcLvNVW0MzulTv9GA3e/u2DzEL0D9NBIutiaruLo2qJLU8j71/CYGyqzxOSYRY5O11hcw6gAkv8z9XcK45Sc7LPCmrSX6H9RWH+8kF4HBeXBUcWydZJt6xnWuBgQWT7AaOv8yrtwxibvkVDtA6COr3j6/FEAJYgdA8zW5B1vgsCkUXPSLQq4szqhJRKMSy204HgEIXfZbuxRUw1OTngiyig1XZLPMpGXgHRqKxfpQMyPShpJdZBCHfic0nMmHbJTaxILaa96ELgzUdSy+Pet3Ed9Xyqjlo3sIAVe1ZijBEhAk1aWY2Kib5Oc+XrFOC+Ef1f2/nxc6LhLbOTkJxKrfKGGsMhNRhgXEluLko6CThSJ0myZYYGwdl2aXJnpHfFU1Ljq2St5oyJ04oxNt7OSOjgMKOCugcJCYsaHWpVnw2gSuuNkqaZ7thx+b2gSH2P8UsmZ5fLwpl+avjJdFPLkTeM4wT+KKbP1NqCV6776fYYQPMmdWWcGT4xWtKCWQiL1dNEW8spRsT674ZaCui98B/xsIHXd3IWM2CaBltUcPH9sBarPZkzY7So/wY9XVc5sV3mY5RS+ZoMi8q342rsJOfY/nvlqhwNmC5O7WHH/RZg6tdliHfMz53QtI8lFFBjGMsF6EQWdYo5e3puDhfH5YXUosGpTDqzym5xUweEOrW0IjsU8b4VFF6gFuMCUuYUmj5QXgv70d0COeqogcUb1/4U2XG2mK5LO9RPm/Bz91z3ezUV9uifURVQNtl0EPQv9Z1QYb34njReYnzeSto29eztT081LUAdvs3FQFpgiOmdwANZF41cd6Yz936XfsIrhvWqY4tjLhxBxWH6vgnkS4Ljp+B5fBEx2B66FNgLjsbPru4TUSPJqaKhFwekKaT0QQ69HH7dkg/R0BzA6zoQJ3wfWzQxbGYgZ5CCWw+m+dnaLtULwjcRnjMW96eJ2gDLunM4bMB3ZXNDx0u58W0aO2OznEEerPkDWFnKm/BOi4kJopbKfuK1bw/dH73zR5Wi9F4cHrwE8J9Ojn48eDo+OCH40MlS/jA8zkJMdzUeg2aQTVgwF4WnjScdU97P2ojrCgodf55VZCbFZtvO53JDA8qfYGc/xLFJKjGLZvM8imctaSpQJ04kskPrV9pP5LoGNx2cvbpw4f3px8P30zOjn53cnAMFDv6eDYwC9L/aDtqh4vy++DNp9ODj0fvTyY/HB6//2lydvDuw/Hh5Ojk4+HpjwfHA/eW9Enc8iUv3Dsc33Nw4cEZ3Qt+SIMRHW57cC9MIV3c3UO2m8FHUkpXpeXpl7C8qKvzmhcL1E7ZG5LdblsTt9l5rL2ueWW/qc3Y3hP7TeVOvNcRsvcUHBvGA3Ve5PtJXB75OJ9VQwkCeNJL/afx9rujs7Ojk99Nfjo6eQNcB5z26bCXsZ/+6v9RT2s75icB8p63O8CuFSoGLwwuq+DpN3Ie6wujfuNt1aBfPdRjn7bjDRBXhtYqS5GBFUeZOmZFlVLhxRNSJzSkPyROpBboZSzt6yjNf1/GWUcnwExHbybH71+DwPzh4Ozw+OjksI+1+PolJs7vI/kYHBRWe+DP15KLQNfuOfue/SarC+uCQJ+rP/+YmOcQE6J9RKb/y+Hpezt10E3V6cotX+0nu89/HW2rmk0OPk6ODw/OPk7ew98ff390asgKqrzr8pHpNSJj/Xh4+s+Ts3cHx8dhv8YKYIluR+u7fieSlNHjjk5kUSea6EH/Lv4LLvTU6xb/GTBBzwXfvMl3Jeugg2sADydl0u5ZYryQKu7qO9/NevYG2Rmf2Ug6R/PqhzLTNtk9ph7cUSquDil12MSX7qffvwc14ePpwevDycn7o7PDiDaHIYpoWPFjRXXX5N1SxH2VaOIQc1EPqccQcTJqI4lSImR5ufYylmjwle2c3uwhHgQ38FUBjc2gEz5onH4/OOoP1KZ9sgLFjEAeAHFiGSQtx/ZPSXfx6NQxNBg1/C/NHxPJJkTRxI9IyFNq5R8NjrPVNJ8pBWah7Mg90vRbZVDxnRLqtHUdqPgtjjaZVOyWCex8j02t4vmJPFi7kv8dGVgeyljS9S16uHNRBPuvz28i4jbmgXkKuCHIg22UbybCQDZFUNCikyeI3P/ciAIAoq2sH0exXqzTTk2QVyjWrfb8SPf4CG13tBaJ58lQxRCkSl/wo29e6YwCHQSUC0yGxvwL7n2eiSjAtzXXeTYb8iORbX4hSaRzJSh89tL4QM4JGKegj0gg2uMOK+kI7WwWss91FcKG5/HaJx79gpJkxDaZ/ac1dFSjpMREw5izzJ57M957LxN2Qt0nWfLvK6ijJw8gTglYCdPtUp2SX6Slk4Z6+Pbt4euPRz+ibnl6eAYaxZs+Pf0pSarWJq+Tl0V9r46Gfuxb58VRLIbpcUF20Uil0ZrIM6uiMkkiSc8A3y2VFDue+OzhRGa2m0enLMMWUPf0rGUM6kce9PSvn4d0EPMVZJqApx4bOEdpBWoKvU0TRhyoeaiZrT/CHasal3ZMAvyvyVCG08natkabA2xZEVgqR8cgVQFUoy9OIUJz8nNBjFSAXhicxpXxODSs64ski+VbWJ/t7OF3A1+U8MyNRvF/wHE2MZiXEEw04zcrtOajlLZOR+W0Y+cDoFO7wHJ+VkFpGDB8nJFld53JnjaBM+DaF4Q2ORlWDd3MXMKwDZdYq6fpZSUkVO1NWJzuBl8Mw9nyaCy8nF4fNAhtUj0pDXkgwyzDFnr8jSAHhDVRs+nVy12WRvKZ8YYy4IipCRr1QxftIGn3BYirc5U6hUSa9qWaSDvBnzr3WRj8yijFMrsd3qFN24+UHZmZqDxIMbiN9TnkNtclkZO8caM+4CheLDFMlV0+ykPHr/Pce2RbmX5F5rl1NF3TAW52PGqwI1fvxyo/X9v5BofVyJtDm3elukxe7Oz+k//vNywAQDSjXIDj8D5pijvvxQh1hzuR/Pocdy1vs/ObHJPIVXMMZb7IpxlGs8E4oORlGFizzckR9acChDt09Mo0L0D2Xg152kC+3/Q+HVqqt0MEnlJo+65epRMZQT2bU4xlwfpDykchkk/46oXtvudZwNclGWRL+N8gveD6dJCRDIP+4fX/Md/g3yDN4N86u+AXJhVcH1v1N0g3uDa9oPlqjMn41ZuGrJuGb+/4XD2zeGSeQpJQkpHwOExReLw+R2EtKe1MgkIAH6vfkSSEsVAyldHwWKc0lO3FGc80WYC6YVbCNJqW8PgpeQm/dfbAtXn+8A7Rc39Yl9VPAX7jzH7H61P7Ha/J7XdskvtJUr6TajBKgvVS6YrSJyX8O/Yy/nX4wGb/S75R+j9C4XWwB/WgG984sx+/u1i30fWiu2guu9m/JF2f23ImoaCZciSpIO9JwpOaPoMjgOIg0W+GF2SHIIedcBoFrAVxuspsikZ05f7h9XwnikBPRkOXwNDVuMXQX29jQur8If9n8xgeXF3V+GUFeVLamGttLkFuqOWu+DqXXBSVfEsOA2/41txghgpumdWYrYOTgTcF8B82hiktKrxFtJVuX76Ea/Qi4w88cfM6XyAn3BQNmnZ0Ags72JocaiNQ4tsM9RZjBaAgs9lEyb+hxDCPvE8DDpz2wrYPI8mcJSSeK61P1VXpyvrzqelYax2oLQnhMAMRLWowDRWQJaDe9xlPDaUSplSRN4YSCSZVkYU0Jp9tVkT5rcAgbhRhSawFquEtezragg9af82RZQf/Wg98FM/yOUZ34QXJoJOBpn7fFJi3BeQhlAtPotYFa5+TCZ245LnRfJ47LVMGZ/wNvcL1toTpI7zjnpD0muYWysSEKUKZOiZSz7IYIBcEiYQQi5ILAMfQRtcIg8PxAYnDffA7TO026EHNYkJQdFk6mM9VClF6A56102vj/oafwRIHs4EGq8Y+GOcrufSAOUpMiWlkcR8L0LnKZCL5SNYB03hl1uAtjnGUXCYayfEmo3Quz00bd/KPuSamIpL272hLcW99iEg357adJahP4hg8rZO77NgELSS7nfhEAe64MDVJ+zrFjhoTkZf7qixVe/fkYOhwdbGLPvEVRGotngeG3go9Iz3sixbCD1TQCzjpmCOpZHJx7+Q0XT8pyYwbZxRQXmFsPmna7cjPYi0qjbDCXcqPLPwPnOptdxA9seCELPEODBMAjRzPJQnWyIX2PCtipIl52cPi5PFnw+MPgZ9zftSDtqfMXQbHMBxvck969R9dqcvxiP4QfJF9twP7nb7paZIvqoXtvVa4hRshcqm/BvYxFnbAacWfjGns3oZOF9Ej0FTjr7QnUkOVWa9T8HLBmEHuXJEOPelc8UU9JBeTWT5z33d+LDQ7PE88oI6VmH0fnCEyvsZ2JbyV/prFdUlG1y6sHbh/eW0q2rv1WWhj1RI/oaRXl0seyRjOfmRQntQdU9Id2ZEsQK6sRV4SV9MZpvNUGuWSJ/DU3LHSGb5LinTVzSmrLV3EkY6XUt/ExbWO/cxQaACZcDrdAOmwheHYcBzMfxaARj6JYZI3sFbhuIWckbZc8tBimTyG655rmtPsqcY/HjzP+HLyfinPQUhis42d9r/+mtEcG13iQPTGkl682q8YsauRNdwNefAqXlm4plxAm8VLDtqE/QK7SV17qmWLRjc4H+hVtQ7QeP7H19CZORswxmSFKipZ/ulb03NQX6f34rv9u8ZcXtGLjLeb+0R/nJYs//q+w60netJf4yd9vEleBhaihYniH/OhJXw64d49HjMZ9ETsbcR+wmnXou8UXPMyxUMoTQMEPc+m51peQyPtpuWjRnW5RbdnT3szsE6DsyWkJWj9IAxb8PzlEqvhlVkTT+cto/bC2zJ2zHUj/BC689zCFXafVrjHhcE7hVgjqeDB+R0Foec/2uwUqZnoNxkk//ufOfIoT3nqqBDzHmtCE/swkTyXsRgyj85BzBgpbhJpYkbxY2uotRcHY0piy2iDX+jY8wJaTff90ZDKDuiH06kKP+YyaKEC7ILRVJhi0Ca0QHpBiLFyHcqI9S53BW1tEQR8o7vJ0UEyHODBgY5E3x23Eb4Lj7jsrDbKodRmX47Nj61TfLGwv69ciRUFDcVdB7m26tkYXpIBGBiIT1HdO6H5LMw6JS9noxV8zjqvuwR+bRjbYjD0vvZAmAivtItgvEP9gSCnhTEtPN1L+fK6Ho/Cql+hBEui6CUKHdeCiO2Lu3hzq1EFIrQHOjlbhwtoVDKSRzL5pk0GeupMbqmxM23PHlzqyNeHmM1aaN3dFKqpcYlFeBNiJnutw1o/FiCQIPu1y61nR6sla8g6Fzp1xd8XcYm7Y/EZ517eiCiEGVgCqkju2M9xrxr+LgkqVeIDcbSayJvfOKm50moRdmXsdwKQ+rM153zcGW/7CVzv9nR3fbOI0FtJRISeRz4b4xNbtvtwO/mSRsUnCOWZGBllIHK//Wso/IIbmq1LY5Y5NsF6Um1kpR15Rj1zqzg6h7ZTaeH+7jQh4TybFVbZNH4JmdJIjIxOGdqwRgD2L5gaNYyxio1c1xHnoo7n1FLcOIvJxBl+tSQWx9ANSvW6c1m0gkKXtjcequr1YjJ9bUusaSdMtftYxWsvab/WNg8yh33vpbgJQnSlNpr4n3PyhvQMAkYe/KLP03IUxyhLkqi7DiR5xBxjMmjaNxmSSLinu05q4WjvLm94X2SzXPfN8J03FLbCwyUubPW7h2gEtuls/EBHRl7wN8b0GytYW/yKGNyd5Gu8yR5/XpfSGWd1fU+OA7bz82vFPKvhXGla6s/mX5hTjr32OmsTFtxw0W0xOg2dRHDNpC/5zrhvzCzf+OXsLxJmPAhlvRHnFAgaVkojPl8n/t194kZwQd3BEPCDo79Y8/TJYxLEMQGgG3v622v2yzWmA0rQb3EYaDMIYKFtHRLX/CzBVCBVjYmuMBcBfaud/Z+LFYpAdlyLl8veckGjuk/kQ14moS26odGkYB+RALPfS4YM+fIEQpqPWdokAbR8hXncjC+d85k2JTQ5Zde7gePvyR+ZvqqL2ddaFICEgdr5dHvCnrUnnLnZOKeh9aQ5k8JexKQQ+PtwcuzrIw2f/jTDfKgr/EAEdGi9HoqORBfu13zUT5yNEZ2Pbrtaw6OraBFeGSOSgI45rQm6meiIZ2jW2m9DkalRkEp1Tgs/zEAA0WAxzUvgGYyysc0Cb60/b5O4CjZbiXEBnDvOBBUoatO6c7yW1qxiK0mg8jRfvndoclr7y6gYOBbLyWvoxXGaOCSt8ps+Y86EB0lqe1LhM7bM3jdxRiTzsb2XSIkGXfdCIdCBv8TS1af+8rLKipv7rWMB7UyW7ZA1Oy7jGpZtGWi+j9KFQYJuTA1FTylPDQU+cWgTqd/I1+yk1YmbeDbmQ6N2GFJF5Xcc8kDEAKd6xJNBbVX9RJzbkYX7gQga4S0h8Cj503j8J0SasFdJ322EzSPbm0m7L4hinj44c6z6Aw3h77K18p2iPmOhOCXFR9xgBIc5O8yVlPfRtmSXgSUA3qJv5MBpn/EGa/mwYN7zzgsTERUEID7qvLAI9H7m48ukvsSTeUGRXTNyKPODSJDHmWFjBtwz1+CQo73iJ4pMXx0qXklfSIkB8sTuuikbdrCZKFw88DahHft4RPMoCWpxebQQHXvSnS4XysAUs7LRkzclogYszuz7pYFWtDoC2lAwKondc7mmK3SdlHaswAHHnXJzTYw1sZfF/14xTvEXVoa/sWwiNuxBlMN8eS7FW53GXdlu2FG+a/AtJHxnVOZJmNVj5bwRgRGZ1ZX7WxI4dCAS1r16NL5UXhsF2nVt/w88R/4TUEsDBBQAAAAIAAAAIVysQTekZRQAAKZPAAAdAAAAdXBzdHJlYW0vcmVzcG9uc2VfcHJvdG9jb2xzLlLdHGtv3Mbxu3/FinETMjhdfUobBEqvgBLLQAvFCWI7X2T1sDru6RjzyAsfkq6N/3tnZt/L5d2pTVG0Bgxb5O7MzvuxQ33Cfqibjt+WYsLuRdMWdSVy1oh2W1etYNum7uplXbZTxn6Sr9mMFS3r1oL11XLNqztYvy7arm6KJS+ffcKWuHPZd8W9OG35ZlsKfNQ1fNl9rXGwM7bhHwTAKTZFdcd4lbMtvmo7US0FE4/bslgW3fSZPEoumoWo7opKLDSEP52yFZygg/+nGUvOprPpi+QZ4L+kdUBDKXgLKB7WNVCitp2eWZqACKBEgmXtsgZUE/ZQdGuirqqLVgA40cIZOZAHS/lyreGyXKxgXx45obPDP2QJ5KXPGB72BRyWzVkiGbRo82QiX8zki2XqvmJJyZs7wFCvFreAvEQ0bb4Ati2QsWJBp02yDBlwYUiEwwJZuRSX4E1ZwNkMBd2ad4qMlhVdy5oe1aCqO1qvlgFAWvgAMAUsA024oJWac0QM42UjeL5jHyrxgNsbsQKGsg9CbCXwNW/XKGeAR7hE25cdMBxFDwwr8bA8R3HBXpCWxKpwgJyWvGkK4DcoIi+LnCNPARbg6VHKvNqxGvEa8gC/pAPO6klJi0fJyxOSeZmxf4BAjBLBmj1yTjNYW/GNaFO9I7u+59ttuTMPJhbJh6p+qDKrWex3RfU7Rk8nrKzv0JDSWZbdXM9unn1Ekb5dG3Zr8qRaI3FA/EDiLsH63T4L0msk2cWKpZH9aqOzeD5nsyuwvlLc8eXudJYwUTpu4/mlh9GQ8z0anKikraOcNqJb16in4rEjjdjwqlgBg9pza4zsVqwLeIcP3rxkG9CfAtyEaFxiaeWC9HOcPst5WPNkYq1hBtS+eWl1AmllrAV3slxbrUIjh6caAHqA10SacqnKXwJ19YoelDXoAmAgCH3Zt0zb/zSRwI7yDIBoy8G3pj46uVfjOoxcoQRHpXGRsHAtqGMpTgmlklY6m/7hq7Mv2SO934i84BXjt21d9h26z/uCbBg8wG3dk8spGrVsok5g8eXFagV+AfQFOCC6ByEqN9Cos7eslp5DnUdZA0SYCcuL+yKHWHW7Y+0vTZeeZdk0yZRAunqbJu+qtt9uIRzCKleUoKpdBzCm6F9Bf2PKcuZpG7oCDCVVXj+g4Zc1QJmzF9PZiwnTAvlCE7fvT943xKXFps4FBQygqMpBEvZVC89nxwBbNeiiYPXZ1YTd8S387/W7q6tjtkKUdvRLMRtouFL0LAyNx0C7Nx580Xa86xFUAsy8L9DkeJkcAwRyBLEEWS0KcCQNgCRGHEsQ2IgVsBeIrQs8KdopHJNPiXGpkmeWKX35Qftc9QI8EkTXW7AJhrsku0FpFLTke1DXbtEm0t3LcKFhsk8/ZSfJG2Vni3rPUnk+Fjy9DuCBszIIb1A5h8ABysdjzqbpfQdqGwIxURaDHYRnZvwU/KuIb8QvfdHAOzgF5DU/gmFiSvPGOBWWfGNUq+MNnsN9JsCXtRGVGFJEUA2AS9qn2X/CyzI1R4mRmTm+YCBbtfOcAWTpSjWsCTiisuRbMmrEm8FfSTn6NSCax9ToesIshLyp0RpfXVy9ubyZwBka8DftRfsKkmbMI9UrSjMQYtXUDylAxweg+b0f7B4hi+C3otSaQuSXEP27dfqYSdqX6QxtVwuX1sO5pQ6vObhUqAQUbBCx/QHSdMUUSSPKd0soCcEUaZ4DYEe70I6gVMAKQDSpa3kZ+/VXtTF4fiKDrZTGt+sakx08hPSeZqE1rxPPoonGvqKMe2/iZrX7sPeXEn3+Q5AiIO/PrsKXqFwkKvhXvbsc5F/jWanPDgXgp4HfBBgDZ6oWvynuwJUuIEGWgZOcwKp4BI9pni1WKBTyuQtZzCVqt7VI3Oc9IXtkwUNluGr3K4wGC8ofDAp11l6kFCtA2S7A5RZdD4GNnhgGA/ttcme3kejRxN3cz2x66UVKs8kLoLBZL2P4s9n8HYQ3J54Od0MUT2AR+JuNCbxm97c2ETGx0RJLQRe2fztMVyx+/riAgOyiRp0Gu6l6cFvwCvT09QUIDhgqs054hmfij3Qm+MkAu4wExijYSAQdoImsOSa8JvoQOqvU2z2WDzIKc8xYvuGIwKSfISNf22RkVP2chAVBQjXVN2ID5cipdC4DZfTKBdcu1ZK3a7BjyDfzRQ3eUS9M/jxPjIKAzfM7sdjWZbHc0dtljSfvxELFA732R9V7WeQQGszJk6JaNlj36WVvNjWUu9g7wbcVuEZ8A2FBBjWMDfgG0yEKwdrPKG8hTO6qwwj8M5bcHlUvBsmSWaD968Wwq3QgYTpJQk/rpScOhkZ0fVNBRKN4aI/bV8UvvT3LwHNrTCoC6acy+GBoQkvh9jk8OdHgnSh6duXlDm4kGdKs9ksyITNAw4AizCDxhTAQ1x4pROPKeD17XGY2yKe8bC2WlcV8P0FyvXriO0Rt4ol0adIXGkXQhM2xcsmCU9skKvFirzyxG3DVYZNBFKXTBdHySVlpNBP1AlIinZcTZPST0MtpkOOejCgbumXFY+Oc8OehZ9LgA5eEi4eeh0BqR0Pi2BRtq5wOEIpVuSOA0C5961Jbh+VT0WqwTnatnuxJrslFUCJszYCMcnZ+hpbKq93LHtvIYD659QGo5GjJSnFPihw8P3XcuJul+usnaB6OWSDuALVSTtTqrqkhzvk+QaWv2waE0+xMHYUeewJFVFnWD7I5wTtITcAt4moFCxdJmrE/tgB0/aZqx402YjPL9De3joF+/2eVT4scxJqq1qpm/rXLlhun0fqYGSm7sn3MlH68vrD/B136Oy5JofbaPLTwIsu8bmy8PiSdl10w7F2qYsH2AjAyb7bdTtd0Vnd9/exq1ACD3fYUlVgH+KUemXIbT7BU6mJruBzqH8gcOnlMibvqN6IploEaHfLzY35d+b/jfHtcMZ9Q/Ecd7bF+da9P3atfAct8FSPlUitcdYKnK0hXwZc9HqVJzlVYpy1nqE8SJArxFhgpJYo9mFSeDfJsZKmNoVGmR8OHaSw9okoY2q8l4Jsb00pIA4KxdYXNFZksod24PKBH7FPmM4T9imjm7IVmhmKHRGYaEUDwtm4LqpmASKxaiFredzX4mWKpuCb58NERoXMFgdUqoHJMitTMnNK8cGteOpojJu2fmMDbCTybdNo2wyPtlfdZpjzJTbUZUOO5AXuCeN0MRxmcNaiR5RLFSW9htDhlMxDA3iWQAncNKPleQL4mG3/veb9K3HGSn/JQg5a9LQKMqFVOPI1lgNgVsm2tI7hHklRHvKJ7Ksb9xkM7Ih6mC4uiws70bSn85Oajn+CMVhvofDNbVSgpjewaxF+902X1RRl0nVuIXDrDcCITMd02iJiE6CvfyThqXeT4vXmW9C1UFuoCHH/WVUq+uN0t8GWQ/fwUnsAy3AKeMBcsmvsArDy4vbYNeednNJl703fUFa4jTLnAKQVPFFqnArW7D9R/gFGfWp3tayRWXbPjFRUsVQuAbj39UO4gMb4DDyJXajIVF9ymRMgINxULCRvpse5RqfF2KOXg8a6ry5wntVgNQ/b3Rwf0heewB7i+lntublwxeScALbaXtBB9oCbBOQitjBhzfAHSDnrilMDJ1+r0nrzkTMlwuybUVDTyYlw+VoMRNQkAZycm1HnHi+1zcJUlKIUddaBRCTlEIedYED8BpItwRQxn8vJQegYkCQdsJPPU6UAj9KmewiojssOUyvmSEeZx33sFcyvBAER0fEUy0qoZ6BBp46AZPj+mFy5PF2nvzam7N2FhY28e6etNWKSlN3c6ehKLbejNVT8vktpJv2PIy8ZCoc3aUKRmvX0ctcmEqm2Vf3k2+TEaL/yWu44V9oY6Md1Zzz96LXgTCdQ29H9q20iOFPbq57H0aF8feDxXivek2ZmbLEXXDLOl2LIxNgTJx8QZeZHX6k4aJf3CZqQPDi918oQjFGAgaCsd6x7qp2SckGa9gJR9JNmch7nxReVkVMsGMvwGoWCebnM6SY8DpjCXqjla/t9FUzuFj/ZZsjEycvjgSspN8sOickzkAYg/u2WCX2qOQVDlBOjJ+FY/TRZWbiR7LClAVMCvWzk1Z6p62j8xQzam1pC+RcVpxK/28xWGBdqlGkZ1CQ5K3f7NvmSfs+l3fAluRjzP6x4y2inO533Othv+mM5wFqeN1UNedTyyhp2yUbae7hFMBkw3x/Q45cAmSiHkLOVwps6QHAOQXKgwj3LUJd6ks04Tyd4j8T+hpQW1pF+w7u/EGIHpOltRgqqNoaPkW11qfMJe1vS8LSBTw9wPM0ReyaFGDf2zdgiTcZwANHNhBI3iRoEhY3Z+pglm7GdUgy+uQB6FejJgQiif6+JmMq7Z1z/fAJecSwCPX74W+JCIl7g9M7slQy9xytZwU7kIl0K8qmqLXCbFh9mjGfxRx7FhMTdy54ISagCRP8drnlpi7C2Yvlo1iyg0pf7dmPMSeydymMM8s/dZ9qTO2/h1m4+P2jEQil8VoqQCjZrEmVu7uPhc5ySbUJZyxX9GoEi3EZbKE5wQfdQE5QmZn4PuL7BjC8oBCq/7X0PU3Jka92sguUW3MI/pIz5haEiFHgXcZ6Xfb1fHCAotfT8S8nuiIbuzJzT8G7JgyGXZszawntOpxjrY4bLMj2AeIH0rIuH4z6Xa+BcirhAvYfVucHZruBW4a+pD37EGwaEScdsYX5FmyWuYVeF2fYIbHLoXCoiaMP+sE8gjcZIqeR/2IL6vwKViMRUeFDINLWoq89D9+Noe7bKSe4GktUuV+CcRy3I7qub5b9ZS9VvMflgaGDLFHYxLHA7vd5G91qr9qkP3nv0SwCr+iDtRNOvBEbeZbMLHNYXyNJSdr2OWPaG2mksRzZXhK/Q1+gj7WWNjjaw5/VSjaJzZ+tpUJX4gsdSPjFAAYxdSz2PX9xNpA44uq0FW/RHC+PCGe5sfjz8Z88IY77rGQZwQZhvpVAdmPBRGw98zM+qoH6H/REzwOIbEGU3RwSwMm2ACYWBUmtVgi8ozqMApYJqnOBqMVtJ+5+LLGvFvHD+0M2hosjRIyVrxywKChExO5QqbonU4k0+8lDaC78keKH2iCSE0DrU6vPqRpGWWq5wAWJcDUFAl+s55JFHSlPBgsUmd/t2xwI8HNPngGNKW78qaj1nQ/7DpPN0CSOoqbla6+TsI97EpIOOI9jgfqVse6Y6pNQDvUBQG9W/rphukPtnQHD/KMUMPm/IMT59n0v3rwcy1JXbfaDUJ2Hy1sddrOZjsPLdDsxkqH1JsZvTMrL4FS73KdrkWG04fJRTY1VnyqS5qTjUh4GrkMjOdN9f8Mi3V4HOveYyL/hrLTAOEesXz2CVy/J7renYjJarv2iKfg1oiNJbYByH7UA5uqixa+Snaoq/M7ZFF0w67wLHb8cONYYVrz8J8tXphEWt1mTsfBfnaBa/GsghIY0Cx9k8AIrxS6O9fRrwhORqVl79GFdzypUiTn9u6KiEdBJ2CVwJbD3P29sd3l8NRrb+++f41A6WE8/xeHsvm+ViNa1AMIH/gdypzdjz3uF+P+XLcjCCVQdsJHo3n/Lyr8UypAjKhyQAQ/m39qIiY0GwpuQE0qUa25CuyL/S4R31rJE0BU2Szi+XFHd46AIMv0J9Dsmz5Fnp1FA5w8wF7o1d4YUEUSKHRvec3u46ERdsV1b4WSK7HtQBXL8izYxXnTuT6LLNr7GcX/vNg8tV7GUyt4Od0/o0zaYd0sFhSkEZuebd27gr/C9pXtAutQ3eN2JZp8rf379vP37//B6D0CHwmB+z0R6FqY+YtkpbvfFeDD6dU4LYBu8Lj0/a8FrKzTVt0b0p/UcRzqx0a0IQ98KbSEQQvzpZ1ru6J3r19dfoVZpvulOT7ypQojuV1ze5bSiit7ayaekPWg1RDQMFoU6x2L8FeXil7IUGMWQhYKd2F2a+nM9PXIe/rq8a5vGCq8gIXfyfaFmQEexxzMbOY6uDPZZCb7A2F3jSnLaN9EGYowFH9sRXBoIY0hT3gTuLvjhwRNzxScR8vwOT/vFFxs0HF7Far6QhyFcz0V+n7PoHyIufpLIn6D43Gzxhi7AxXBOw8ia+T3BpQ+V8YrNCzcd6XU4c5fWzaQAQdTB32mESQxEwiJ45e71pFIwjMrFZR2oksGpfKVuI958EnHr6j9hpO3KRAvGn4zkaE/Tx1/OzRH4OMfu1hvtMbQKSJ6KPB6S/79NHVF34aajBWfQjs8KPAuH04GIYj2oeQRD8c1HhiHxBqVIPB70OYYl8djmmu01U89vpgxAv7mYnnTrGJq65L1eoWx6I55AJCj92QvlJ6t+YYnfWMNOt2W4EzPfjrbaSvyMWy5PjxiV6iu0ngTdaiIXAdWDtevN/TPG9ZnhJovRAwqFlgOBDVoFOGIyKE4n72GbZjIbUEThEw1dyhuhv3Pgic8qmd37nDZC4knWJfdox+n0+LV6I7ZvzHdF8LOxgZjs8HDz/wyGyH1jSrhyMomaMfzmCKzOlNtB5uO+gh7Cio5M6//hXT3ul258+/+BFT4Fw1qfhViZUB0TDk8YT5F0m2RJdUR5sqT2x9HeOEN6Lj9JUecVm1BQYflQx6DvGvTPw/0V8HEfsOJVK/Dz2U0239IHZ2akoT4PRaR7Oc62vYCio4zG7MmyCrGSSKLgh/NEsfRK2QM1rBrvBS/G29PS3FvSidjye1QBwPZz4X0B6Ocm4AaWYM7T2413LTGhB0JfZ3J0fvoLO9l1A4613JPACHlKidr26gpv9nTdSDLWOvuaqLNn2HEu+D0y/SekoT/N8vvBElXd8erraPa9c8sQ0Ua+sQfw+1do7r2cgGsOxUdcTbTmy2WHunpu405eap/DK8wpq/S3HYufqQ4j48ZZ47PRzq9nxTVCky4m39I39IwUGd9d3qK+IIem7aCGv7Cp1U2tV12Z6fb/I/tv1GQsVm4D8BUEsDBBQAAAAIAAAAIVyH4SqcCWMAAK5VAQAQAAAAdXBzdHJlYW0vdXRpbHMuUuW963YbR7Im+p9PUYbpLaANQCBttW3KVB9ZF1v7yLJHsttnhs2GikCRrBZujSrw4rZn7XeYeZHz6/yfR9lPcuKLS96qQFJ275lZa7yWRVRVZuQtMiIybvlhdnjDf9k3xWxVrLPnm8WkLpeL6sbShzsfZq+WdXGQ1edllZ2WsyKjv9Vys54U0+zkOstXq+HrrHt5Xk7Os7yu88l5UWX5bJat8sm7/Kyoelm+4KL1eUHg6qKqs2pT1kWf4GRlnc039GaxrLW6q0jfqmJ2mv37v/33bFMV2cFB9vdNPitPy0mOrg93CNwTqrycZ0saU14v19kp/X9ODc7KxVn26seXL7OLfLYpqp23H/3yy0dvsy8H2akOvZv3s5NeVp5m3Q/KarjYzGbdvEfdzYoZNXfC4M+Lybvs8rygzq9pINny5G/FpM4my0WdlzR7myo/oUmZ5nXObS839WpTV8Od87wa89uwxate9o+djJu0Fq+oxXVRb9aL7vPHL9886/nvqD48XefzIiy1WC8vAehRNkLZWbE4q8/1xc6v6PQ3m3m+GJyuy4Lm4TqbltVqll9ny1NewPtn6+VmlZXTYlHTXBbrKutuFlP6O1muadbXRT7N8orgVLQSWIZFtprR8tRlPcPa0mrmV4QFs/yEJuqhjP1dcU14UVMzm0W93BAWTHvDnRX1ub4ec8l0Is6qzUm3M/640886Gf1z1eMJX86pVcKFTb0MmgUCFtnfluWC0E4GsKB5qWhpgJP08ZoGV2dLWqEF0IzWnxBsiWW7LGk1cwJGC0YLR/0D4tEUXQ8YoTfUTD6tsulycY8gXBTr09nykhtEB3bQlTF+jbkr0Ti4K1WfoF2NJ+c5zeVh9sdRc5WlXC/75RdbMHtDW1LLZ7bCnQ4W9lf6XwdMLa7yqi5cc5PlbJavaFiHmLmO4cwCXehKpV725aHvVs9AMWpXq3W5qE+7nY90MiuCEverJ6j0A2EDLUY2mRX5mqrLDJQLwnHCEWwIAgREWi7yGfX/7xtsb2yERbh8QwL1mgfntuWsqKvs7IzLrIs5TTvPuDRAnTghyDRkwkXCW8LLYn2BuoyTBA1NUOvFfFVfE97ReM6G2WNaqTnRhgnjZFV4nF3TWi5Xm1le0yCoy0tuTMZxcg2yRI/z5XRDjS9P0Bhv9ylaJrKCibtazcoJgbsGLSPSRANCCyebugYtWq5kErZhCr9pIoa+9nihLz44zPbwFuXyoNQHi595kWnE88tKP/QSBMIEGwpxCVnM58s1zQ4R3/k8X18PaqZc84JATWgjlfU5Ea+snK9m15jp+pzKTopFxcOneVgRtyio3LQ+B3XAnBZX2KsyYUyTq2H2jN89/+mbb+9/k89OBz+hQkZEi3ZYUU3W5UlB0CbYZ0CoNWjMQ+zY5WWxHpzQ/pwmpV0npHBWFaucqH0xux7unPKYaL7nBHKzroqxDGhsVC9cA/1Eu6efLWjTFvlizOyBWFEx1588kXf/r17WtOgEzXo5Xshq+MbQh7waYt1oaop10A+s0gLfq82KiGVV/ZTzHqm6VAFb64yKL3oopg3dVliLcRXfpdtqBZ1HRV5XnUlUDcbyUbn4KJt0O1hg0B6s8ZjXmMiQ0SFmToqTtIFKWiLeDZ7wDD89y/7H/5vhb2fLSvS4Om8W5lWTfBatIHdh/BJoM/4KaNNx+8DtBN8gcbnFIdE7G2mPmrWuETb0pLVfXZvRFPzLv4RzGQzuhoYYwMNse5MRyoTNbwPZ1uVfbc5v7i/eGAI9ajAca0mH1CGRTxsChYD8NbofjggP6/LsvB5YKx3bNi3boW8vuUjPCFMHcGlco15HyNNP63xFVIBID8l3xdlyjRUPRQ1i7US455tZXa6oX2DxkEWzivBjcHI9wF+CE4oHwtPBLvHtklpoE0b6QtWIm+59LlNzQXLt7BpfXCniocKEidsAULfy1fhvL+bLf1l0egEx8bt/jwr++ObZ8NXjb5+9oZIq9vEcfF0SHywgBZ+ABhITEumEWHAl9Fr2Qn0t9LmEeDYtViTp0WuwSuJMBOeEuB9mYJh9W07WSxLtVtdgYcs1Ma7TNUnM0/L0lCjwAtx6BumPOj+fLxfE36pz0HaeQEiBG0j2zO5o5omMVRd9FbvenJeL63sVyfRzgkw8hfnuKge9L+YiHy9KkgmGO4QNm3UxltGMSxK3whWAHFY1maO+BtuLxWH/gaVhfRRRSnFa3inAD7od13JHSBgPT0sR3eAfu64Qevfq8dgt2hgEvawqopv6VRhzUo265Dh0RPHTctYxEkK7Hq6jX7R8Jc8QynRsUtt67uhQS3NWT2mLyH5WXnC5pRbA90CG/z4mkaQbTDC/XQG3B52QXiWDO/IDOnai68hWYEBETHeX9oQHG1ayXZMem+S/D7PHWTXHCXNaUK/n5aKsIPKtlrPrxXJeEs2gsxeh7MUSaAakB2DaQbqlAkgecRlpBatPChQ/zdckklWVEI+KujFYrelsgX23ALBrYmkkSg0dOG6Vxvu6mG4mRddv/eDEOZmASExJlsPv7A/Z3id72cf66qOPsv29Tz/79PNP/vjpZ6EcsqlPP/9h+WJRd4vFZB9POBD2aevTkYeQ3hXtuV+eeYw+vxp8RAQ6C9g9ukr15bxw1bNqv1qjIbGKqXbLPpnn74qh7PObkT5BH/4sZO//xumRyDKIP5F4O5FByIdgTmf+AhIe6BsTwjMoH/J6cj6kqgSTCs1xSqeFAcUqT/j4D5GW61eE0o4kgmbKgWNKHXpXrlb4QS3xK5x2UUXFx+EORNX1FJi5Ls4I09bXTcrVz0iwuiiXGxz/WPh2pOyDV6+/+8lt1kQ+59cA16COnnwVV0BvGg3PNA05nVbfeDDVAE6y3rmnVAYHn2g/zrndermcVQcHgDjGDhhX+aIaF1d112ZrK4VgjOBTP+Q6ANTnow+4Zdn61jP6cGQ9sALHVsGXr/OSZLwAnybdEGRfGtIKHp/GHfzezLv8gXtGghCrXgKmwYjsMfcwxeU+ekOv0WJfdF2HbWRVCCRDE1nCqkyWi2kJrECfaMPRGGppCU+ERcTOD5mw5sB0HpoCIiKzxMfNgid9WpzmxGzHLMqM+WMlxbFrx5NiNgOqjV5KbTn/Vo+r59TPJSshWKjYEYpgjOYDnR4TzKmzLFLQNJFAz2OBRO+GgQcdA37yAKQA9Qc/tCfKDDDlulJ9BzxCgwYK+HL9bLperqzjxw7gLrcaoMg8v+oGHxRQj4iocawYEUIwehAzKrjOF+8CWCSnlkU1pBfnS16003Jd1R2FArbDZIP4y4Q4wKzIVWUxy1mjub4oWUPhCEa2qeQcXejyQmJTWHi5KC5NurtcbmbTQFtFMAqSBIW35Be0K2yhqqKG5Na9C6Lw0PidmyZ3vAE+GDoE5fw8nxf5tOuaTrfXr7rFQL6xhGM/6uTAnRObmmI7FbNpX5SxgoDEwenoNv2N2OckuzLcoCyBOVUKnqBHeWkSZLiTqTsoqwxPa+hrX4lfiNhl/fXwdTBWWkHFDWkZp8mRZ8cJMEE8FE/eMVu8IhE1N7E475bN2gyU6adohdzxGbVkCCQb6xTa9odYxMRqi2pgsZkXdASzHjvMwfhOIXkUXQPB43PwaD5OZ8vlOvrsBBPfrCi43PMjkeV5TI0RZrI5qnAL8HCC/QuBKN3xRyRp0tmoK7WpAIlBpzXTYdfwgJavdyyQA5lXOrKTCM7xjAou4nD9wdm6WM26nb9+eHR0cDUtz8r64Pj4H3/8dbdjWN8yqBia7AGGpmcIq+gWmEcGZlp6Enl0xCCOj+klynDhZI+CeqlIM+AzHs8OYfNFWVwOsx+wx2HMoDeiWc1n1xULSjBzgFyVazX7MLML5SKAiHb9+sJpwJgHry92hT5xb2TULBbpWv/DDcWGeNSytm0cwuQPmxbPlP1rflRtx6pgBTDkgg3UuXjqSnkc3Rov+3xWfgkCf5j98PrHZ72ghSODF55y0k8wbnQ7fRF0wrc95ZqrfFbUqh2r85qkMcLvV3zE83RZ6jeQ0YZvE1xlMHLpAZFeTmsaFr9ScYreqb5392tUaHxV3bB+TaYztOZ0rUPRCAiC8CBFv1bDj7IunSkhepGYddXT8TK/osNX3Z0VFwWLPCJpaZ8OM1smE8T4b9916dB+sXoF5wHPqMbzfEHnClraiGPpS6HBisdBQft9ZD92Zf8S0es0jwqdNrRdn5REkXxDkV6DaUE7KJ5soQSJXPq3ivY0EeSDg3r5r2++e6V8F4AZLkOksy/JMmbYeTwleWMy2azXdKQkWWQzPwF5ZdrAByk7AFS6bbAkxYrFGpFZTotiepJP3hEwMR/W+WLKYoqZEmoYj0wEMn28rFRuxmTVZg13lNtpsbE0nZ605LWsin9u6NeDoqY4l+MOupd8rYqCWaAJhWaG4UpM/2U7oVwPBeVRv3tVoYca6Az5qBBwDql1dIT3x8ew0BHrcYZdFBZ5ntryhZgxBS8+pko78Rm/yrqiHkaRflSciUb/burHe9nTYlay6Ss7KxYwqYMhsIVVLGYXJTHN5fpe5fSTsALS4q5hFa9yyKwf3gMk0w1Ol5cLxhc6kr8TSxyb16Zsac2++eGH71mdODm3o/dPxcnjqirmJ7NrADrZlDNYhDY1rXR9Tj/O9FSPFtBdBlgSFl8u1+/UbndKgmvFHgZEOdYMyPWE2yHmRgvnjIwEbMUeB6zYqaRzonu4UEUnYJB0uih4bTEjuSqiJrOS0Jj1z+J8UIiwzxNSiXaC5wlTiQnKfjCthOidqN9V5qwctAWpi9OielcT5XhtM/p/weo1p8FW0EZJn+xpGBRAI6zoABfhMdhaVed8zqjs5MIrSzsnqHy5LiEeOQxm/aEMhj9VviJ8IKA/hqGeSg13tJmxm+hw+2pH+657fW1L9seHTpFDu3R9VtRZcUXQUeVAC1a0Fu+K7OyM+18uTmWFSnO6YCjTgjGBVdxlzfqTOtV0QL9hvcD24DVAmWK+wnv+iHqHYs0SWYxe9HpOmTmELYbeCDsWiz0tAx1wazqbAddl6kjenE4D4YEo+RPmeLyDTaspA+x6ZW2G8S15AoEyeIi/QhDnT3ykrnSd+JiBprqosVtWtDGd/lyFdv5SlT/rqcSeQIxGoZ61IuTrdl4t48W+JDxerZfQcE6H2fd0CqbhQymWn0FQJNlxcr5cwuFiwTJ3JgZaTBcwe+gJjqocTQnMuw4iBO3Gr0pBPBKT1vklyBqWwvU7mAOl1lybT2LB6GQA2GkJLZvwLoDHEZuiCVEny/mK5ARYlLd1VPF3l7jiVJyOvqU3+VnR7VSwmReTfAzEHKMJuFNAbHGT6fbkocd/9/Hkj5/Se8/EwXr/+OmYGLMOTEuakVBqFus1q4zcDiv86jFpJEnyfHn5alk7jynfH9H+d57aTgW5LKYHnUCJZcOj84fvan29YmmCG6fS0806V43X59ZLh8BgFmC4pu6U1eZzafJNmQ8JSgNCIvq2YW8SkgQmxVxsVGsblFGXr9jxglgHsI0o8ES4CwAJDWanMxGfiCilxKnbCx1M8gVzophbNSnaWLw9IsImxOfF1AuhbloxnbO8Mq1wPxsOh7JIOpacQchIUkB0WJ3wxGpZPHU71hNoPA02/xXgyUymU5i7AX6jbErFS17EA2FJhqAABFZqW99GTE3RNBClZVHhOjtbBoy5Oc0J+2Kfm37MxNTxhj/db+Fo5TTT1RaPHXpWBRi7A01h25DDK6v8ZSJ3aXfTmIteK298bj9EdwrC9XOxXg6I92zm4WhVOuMyiyKElc4Ls40e03ITekPSiSPzFlyM8EkmKHP8UrRnxjMVmGCReio9u6A3Uu/oqJwe0/FSCMEpj5D1QiHXcdqkqOei1XdNdZW1peSqjewAuXf83o9MuCICB55Lp9rQbbSqQZxCus0yoZthtnR3bvPcuYl2GdFSUTyhT6Gxs4HgXrCJV0jk7vJsQTj4Qgx1IgTIFv0JLJ825H958T2Rrck5xExRGPMZDDzAi34/BEJc97I4ed3TKTgnxF8ss+oahpHs53LFgo5s23yxYN86aHsWIrzOV5hJ15552rKcjI7YoEROz0UKpwoi38PfisA8g3uCuQKvUGjqtv/ZzyxA9fpKl18/f5LtffFg3xxmIYRu6LxUwrWAZPulWuiIzWdPnz1/+fiHZ9iRBW0vGgMsnCyTE4EisNDQU9knr598si8iKR6p3wAlcga7IFJnq2FGh+A5vEmh18w6ZzQ3HXR4s4DzMQ6CB3A+LuaAn2c/z8oTgNHGmY4QU3gMIjmg5iZwBK428754KmeA5wdFhDcicvSRJ+g7pkQqHsdUqMqe2JEsuyhgtjEH3QoSNXsdZt3TGU7eJ9fuOE7UjIXFMTUxtnUMyYc2LTRDj8oTZVlYGve9c3kSSa2TGS1Yl4omMiuz63o/1eJwLyCoBQaVK3aBoa3AAuVhtk8yMkkTORvESLapSTIS9l9/+lvgfdoKjwBiNcZOhHRgoQAxQnP28/gGOb8zPPtZFXKJFK/1GoI8Q9SJVazXon5qMz8qdKWvVeSTzHf4xkRfD8dLv/4YYB1ysjC7gBHFWZcyfJY8eWVPT+GhSq/gIQD+z1I1NlWAHFlmbupbRe/kDBJI4Wc/y/ht+mXK+RPPjDkV/2xGtG/pEEJnZgJ5VUwH8lk2EwxOtBeKqxXtBpD7eX5WTrK90+zzk372/OXX6ucAsb48pV0IN0N4bGFaj/aOoUIn3KFed0dXe6c9fr0fv/78RF5/Gr1W+j5ZT2QwR10YBT7rHfDfT3vHWxjKh9JxJU4HQpgg+QMxB4KpMsHF6cwUzwS+Kv4OFO/u7WF6B8R90hY+ZLFkle2NBuxnonPzcfa5PGuTMs/K48Fcf1i+pgE52z2vFm8QgHy5nPDEzxy8LslvLIbtD0es8ITh83NacO1wn2UidvsmUavO5ys9i2AD03SOPh198unJg1HvIZOI7r775X583ngzSnYGTTvv8gCyYo1NW0+q+ve5WmaEMtlrkSu2tOJ4c7ovtQ39Jicr2UtHRwpYn2Fj3js+dltMRQY91mkDhEOH+Nd50ijpSofUB7OPvmBQfduzh/ojkKr8btZfH2efjF7Sv+Hw/aNrySjEZDqmJVwHEIwoQPLO3CD/Ea3w/mivfYV/01oXuy2rTS+ZoDy0x42nL+H6FrvJCt/wIxiDhyujbnQpQY2Wef70j8E8azdsWj/Mni34CDnBDNIGE1lpub7e8bP4x9GDcBbjjkZj1FWwAulrg6ldG7hVtVHGz64FE2eV+/duCxSCX0cYIjTMXhA5Lv41e/3dC+dHKodEODds5guSsb4t8sUeTPT4sd/hEyl7D7NDK8leDPUVS5vWVlmrchG+Vh9mT9VGSf8TfYIkfKDNgQqJlRLCLbMsopGbCcnr9A62CRYi+VDoAo/QYHv0EQcf/RUd7R6NBl8cf9yDobbzBD38y1/2XDzSdxDGIUPSVMhpAVZReJsX04fSo0x4e0WdZD/DxdTNEAdEsdu0i1QwAwfDiCNExNFCp5POJzAi43z/qWxK56zz9IeDA4H2GkCSiqA/DRDyo2dQdq92JV6l2pWyT4tTkR3EmLG9UGDcIAoT2jZcDB29313DIXkdnPL8SxvBv75J/c6zrBOJbn2c2Pp07OmzLe4fF/la7F+H3Y+qXluRh3J+46Pe4eHhPZ3se7Cp3//rYNT908Ff/jIcfdz70+79ISIRzT/jTwJ4WM3KSdHd62UH8uLhr8GxMrIxBaNM1HPG8kEh3Ixjw5HwXxVjOLSPT0jWehcLqxr+IG4JF4XTPza8or1LyQexI0rsZEIiUVvEkPuWBg0xT5NgSbAz9f5gL1N6y25sfvQ95wVUL98ViypwPqFdSZNe1q4l2lb9Tl8kPhOgj+j9sXNdqW7zQZFGelrBOWCrRltg9LJ/ISR00/eLgbY4yWh48uOIgR2HQ1RTeHkK5Yep/qX5ow+keK+fxS9ogF0JRAOdvu51nJUV588NH/M7C0KsTvafBxzPQKf1lRxqxBcf1GydL85o/cX6AULyXx9kgibDnYnAGXMo4vWYq4e4Q+QEtVW3uMpZ9taXkIEH2ZxYHcl49m7vOHSg0nlERUYZgRB604+GD1CBKKxrnAvdzx6wgZzdvvZGfxXHn9nybG/UtcK8bLztXfX7qALTKmw16ua+ZnMotZPtiebJvfwEL/eHD5LXn+G1vtwbGbQ/MGie/a9YRZKrYm8xOecgMpJ963m+Up+8imTiYpg9Nu9+hwt9UwkSIJxH13D8A5sRC5hzy81dzOFJcZ5flMv1QwSbGhq6qmf0hHCKTQVLlhpMoUBYL2fe2sidOiE6RwDfcbAjBza+LOdeqSi2euLSQc8q9rEk9OOYzHwyYd5NzdbQr4rTAKOamkKEURG+EhtZbtCvGZ24EdkL2+m6KH621kztxPomj4k6i2Oev8SD4spCudzgCZUkRtb3l9Z6ZP6rN/9HqDs2t4eReWRoG7cRDlewpyjvgd1W0wr22jeLh5zY13xFI7BiiPqmiXdOW5VJXSJSUEe12qB+TYcNR1idVefvN4VG/sy8xf1Q+ejMTY6pvO04V7rHwPKTKuh8P+NnB13Vwe4FZF+N95WD7U2Tactuc0mbYDyRyP2IkjuAPXBsD/6RDEApD/ruIfSCctzFFmLZnWCq3ATImTjGx1sCJKPCN6CELwO0SJpgb49Gs9gHcpKAmrde0sE+h0+OBrDk06nt1is6YPhNxW46RHB5Nw9WS/rAUHI62pyTOISYmI0QBInTdZFz5UJ6XEEVQ/JtfoLwa3Z+gECZr6+HctDghsYSM68ot/cSsX0lAiy6XY9S93l1enCmLAZfuPimBMajLJmjYIOYOaAhEz4JNowbfAkt6VL9LyROhF5Z2Hf3o6knjQ/RaDnfzLOPpr3QTBB3LiFQJtS5g3y4Id2xnI+AhDYi1lbl2aI8TQf9B56afra3jxpe8ku2oE6FgDAV0SDtY/zcS4BH3oEJKCK2t9U2MUkXhYMruAC/4ABt9YkI+j2Q4XPDIykpbw7lb1/NJjzsQ/0hLz3vPMw+8DtaPnpaduh/72DaI9EKXA8xi2wINPlKWxOTp9JWFnyG4nkjJogptMm8byyePZ/ROdSqD7OfsMPoLAemSFAhufYt18S6qjliKGNPt6UIHypfVhAa/ODyqwJhZMbjTzeIXXMin3OXyy/fXeZr3yl1bvQ8F0MlgtZydlCxrm+dbZwd9KDg5MQPDrN9ECgOC/EUzL4H2UaUgo+EZgqDYlynA25I4l3VQAKV0iKA6u+9Y+1RU8Z0JXpCDdvoeC+DxCJEiwW6yota3D6Ru5rFtmH2FbumnBc7ok8lbF6Kgxk3fV7AE21pSh4dmAW77o0Ge3tUqaSlgjMSyeJn+aoaBlxoO5uBGC2Oh0bKvxSyaLHa6qMji+XEBDd+XloWpP2cOAL7sQdr+zesvc+1HYH2bz2Bbq9/25r9unMjS2B1RlfPBE7qD/CDyJAfTd83Y4lHhtgCh2kL8HXsmcn0ST6bcAoPXsX/rNsHuvCFOHuolovV1QNJB5GxJqlKXOjEy3k8Pc3eyE/Zj2rPA83hjATT589HLBEjJwEeQiMeJ6MAiO808wciscuLcrrJZwP+2ACagoDNe6xJLX7SHEOSl8QASSsyCJAgHWEDSE7EJafJCKHYO8TJ31R3XZ6c0M4Iq7559m2mr9OKqmF5bGRV19eOr3AA67563M9ePebUSvACw9DhbmnBRMOdIF2HAh8r/Yh9Bt069f10qy/NHQ4OIk8EcyzKh/erarP4W+rqFKrOQxMLmL7DaVhiohD0l4Vgp1azCeC3HaBSFA/uPhs58S1pDE/lJ3EX1Y3lp43riOPW/VLgvc5Cx+2SqCtB2Zs644vtOjihA4nvkcyjDNy24m0txm26ZdPG3XNbJ7JB+NoavKXOx211TIh0ikGmor5jqoZqqLKMSYb6vTTBAbbZmESB2RibTX6JqgX7yOqIsDTY/l/2jHXVxA+RF666qSgT4Tc5dRSaq1zVCpLEDSGRVX4KZ5OTQnyxygU7ROeRQ5kSnqvsz1wZng4GsOvUi/d1b/RCnwli+gg38BWVEvEx6AoEBnThPivi7nsvj5hooYveCQwsfLngxBlr18B99mayXo2t8NgXTnKHuJ6lsd8+xZtXzV6lalerHe9GxzZRng/sovBHe9CzGjmiY1EpRppJqCzgg6kqL29tEbm9NMuY5GQ7+uvjwX/JBz+PBl+MB8ecoG0sBpGWwmP3XT+1lvrr+JcxG1jicmE/7XXPT6gc5+W9ygAs1PHidU1udo7iPTrvzuCuxJou8Q6DmC6hJer9G6DU+iJ7vFrB7WECAfLPsi1xAkqRRnBhKl5Z3Lo4CGUafwuTGafp6OyI6WeMkmMzzTci4yJ7yfpiV90xiL5x6KF/o6l7OoAUUTxXxBYXpipB7QNnfprBu5V7Ctl2uS7pCEUyhSSGeYPACHUIqA6IUSsYaEFn6li3vFyYuhMQ2XLngeKDJNsRYU39aRUOTcm8rFAaAK/RKgQBDscQmUAiMdN1chwgmQoifnueuDvPxtszGBgIttuy5cNUDdv3uPoLyJq25k8JcuiNXQ69uMsB2Q/rKqLgi+C0puQkksbtc7ZCSX624oONTH4bGSUMfq3Ya1yli0PvlGgHB6M7hA1JKf1bt7ufQdekwqwjklw6qA2HKXjP+n3nt+J0WYcNlYvJbDMtkAczEE2XmQTl0q9pziFC+OyOyzjcBTB4BE7ExtOAu4Xj9Br8pVsMz4YIxVkVE6aDvBr3xDwNsglgWUYzPBBkdU6OvYe2fyuLb9IcNTx/ut9TqTf2G+aQlI7qncbs+3lVj/dH+w/oADkY7SNdUmeHledjJQ0OREwZ+roukZmNg0zgnwYw/XhCjQXwDIW8h18QcP7LAZ4pUQIhApk2dFcpCd+1G3yMd3SPvsIdDoGYYiD2Nfs3baMmc4TfSzGUsL2EYYWDc+qwsAt0iLTHvvHDN9fV8ClVCCySQVjc8CPsSk3b5euGabo4dwgieXQzSppKcVUIcyNq/q2Q47IawvJC2i5VBwQ5bz5ZLghNxU91Vpzlk2s+xJab+fgZTIXj11hIifpkOYZTcWm7zt9BghpMaWXJtyR2s9ZYbYneNlotPWIk8667OpoxzZyiYpqrAZ+b4ov70MzDZZ/C0HZ+QSBAD7601IFQp+oE0LtOcxJAD13HUeR1SbiKL+O90fiLEYvTY9rE46rjhGOGFzBFatLxQ//myP8EE5FqrCaxBhVxqIQiwRtEmOWWnpSXn7Z2Be99mvE18GmNPATQvFdJgMNEfY5BGH948fx5+OIgm/182Yd7MJ+UFxxRwFFDnDWUeW4a8yqcCN8131yfZHUk3+tn01VJKPzJaHTjsVQ9oy0GJewf3Gh/hpfnyVlEQU7nGjTHOt8kB5CAY7qyNaouspCY4/8cqyPGMBDSs826sM5xqgG2hU1brV/5+sz7tvJg2yKqJOPsYThZlpzPJo2edfZEycwzSP/2LXLwEIPXCTkRfXOahJILdFbT004gfT3JyzXOJNVyhphP8ci/f7KRSKTTpZhmp1kxPymmVUYS0eq8MqIHTvL1uijeKTDWOW1ORASoeqK1BhOpEdXz/dPn1lvn5M/Nao7jx+synw3Tk/QkX+Un5axEkptuZ4LudnphAjma4l0FC7l9/ZR/0/Jy2TGNN5SDeGU73FQkmAb1qMZzjLsbNvNh9j1RHnZlcryB8wXnLMW8FkOzz0fL8/pQXLDW2JULB8jZ4GFTR3a2ixJSPvcpmxOj4iRgIKs1a3vLSkJ4pH5bR6XqYeu3DolpF7CqkUDo3FfsVN+OIsQlAhwJN55lexnSrMfbKyzFe4x3KBzhu7pXOwhlQPYlGoxlJ+LFSxoIHrWP0+UQO6srdG3/4EBiaftc3WlZEbUP6jfPJySmFwO4hfN6SVICyImcG3u9nBSVpDhak9AtsWxDlXX81/Y8C5o3zRJaGy2FNDN5dybOcUaxoC0cczxsSKPMd1ktaRoH23nz4tvvXz578nj852ev37z47hVLTbT84gfHijgXElDUdyrfEOA3i3cLEiI7JnBwXhSfiyA63fmjXfMQF5/gbsvIGEotv3SavbLuxDk2t3XmxnyYt/UozOl312411HG0Qce2iIe2nBrbKvJqNTkv5nlQqPOJmm852I3z64/ny2nBadUEg3aTT31nC81n43clo5U+SUmcG7rbavd2UUVgiOUJEuRYAVjA5/vAk9jOFCASXVWcc2wt+6gTF2HX0/G64JStzKU67E8XTQdJR3BCHsN4UC3NVPs+nVMIvayzNxp1NLZ9r7UvyLs5JxEu6tNeMfhj1j0vEb3E2XtVjUW/2KkKeWvPNvl6+hCrjygwDdKMgbPjb0+b5Tsf6uvxakmy2rUkR0HKoTEKjS/2HsIwnfluII3kZfYJn37Y9PHmKVsfsktx69AILhh2ifVWfq5XRf5uTCyHurUeaxskV3AAXJJwxZTTNpGns/xsHNavLGhaEijR+uJ2BcDyCWiq26DeVNfgW6amcTVfqqnlUFmDGt51GsYur52Dn3xKyrMb+djyxzQr8Rc25sX16GR/e619rVXl89WM09xA0j//OawTfZPizAmE3DZmr5EEx7MN5qMSFtVnRQPN6cnyyp2gFzlPGlHITs8nZNzMCctFiri9MdfrtOJvaJsRqSoQ3cfpf4T0ycyp0qw527rcQS5+2vQbiTvt4DSrBZRZmcyMkYX8KyoUcJLDgK1IGUnp5TIxbb+g4hfLg+1Z+5jQeK15VA5dPkKxdlRdX05uOfHCYycAsV6WYza/kKBc87JM+XzIszEWsjI+HXVSodULaVtZBhK7cf8I5Olsg2tRJqAtYYL5Tj6DYHQ9bo6KagnExZKWEQkVXB9BbW7o0PvyJQTpEcyoX/zGd4BwdrEoZi3dBBWnCcMJjHp8cj0mlI0650AyOZEvLStpUQc3bxCFpfojKst7IljqfgQg3SyB21UvwlENtx7Ppw9u7UNTFvI5soOA9CQN24Xr8tYE2/00YY5PoyIp/pJELS4Jmh6gqfPIr6nfePLT7OMy+1Eep16/kY6ZqU0gNYbLGakVm5PbTopEjCf8yOvxpoaRKVC8gcICWA263fnxhycdVu/xI0vOJp1Qo4jIEcmB2G/hCIah2NNA3w1qNeCMb2pTJTnhYSZJb51pR7InreUqEM/EM6S5JK7IQaacxn5ZRynyT8rlbHnGIgmyY3C6v2oo2N3T435LRraFogIbVfsuDZsuY95I+VmFWVMUrdqT9KqJ9o0YTyUdTXUB9RAJKIXaZ+TOmlx6JnZNmDLHXDw8V3n0c4mPRIeTqmp8bKicnmp3s8fVrLrC+fIKmX29guM522gQ738gvZQg0gEnSa+RJgju49kPnKd7nT1fi4UV1EEhvEYtMfXsj/T6BqpYLkDqah/OJaw+SVnhd9RKZDOLWL6aHRzg75inS4OXF+yUekjtMC8aRwrnfsYK6DGQp4SY3JlLVHJHA0BtwHITjej/OhhYByPr8NA6WXcCYwT1HlYQIG7PVZapGaOyRl/LfxxuWSLcMko1qiMK1SQ8UJprgxFhmJY/KvvZcS+oYIzUcpGyCPQLI3JP9OAMz7JSDCdiOZBDcdx2OoYy+sYuiMGbX3f8r2AGX5xml1DxbRgzYCohDNfYY8B9lO1JCnZBBlcTr8buZB905BFMgMHzwGJFRmGzr0VnYXrbteSmI5iuzJS3ht9MBwe56rUlReBWvOLOHmoHb8AjvyqqGJ/WLuPVTZlU0vRbz5gE8GgsC0p7jqYb02s1DuXhyE8B3SUXE5skTSFfc5fPT8qzDfLaywbSwxrhGL8YQAyFd9xMZTl4gXC6y1x3/nCHz30Q7lUwltVL8jRoPhuXf7USixOHYgnlpSJHQXpW97Ln1MwW9+xgJPEXT+2Ddk3xqOIZFWbk6zavWAuoeWeo97Rw+iOFlmt/xeqiE+W1YEN3p5sPxp06o9P3cP/NbKok3t7MskIwTR0KM8RmLvdZsI2KXbrUAvEY4T0bnHmZCmvPSg4OouMxlOtUJqTM//5v/10uUhNGhrWkMZA4ms84qQq3rXCY9eZ6TYzFAz9582fV1hAohPiKxQtjWBLGYQYBFBjLRmSfE92OdxmOd0PQCoiZ2tuKjYc550hPxuKSuVl+KlkG5+ThS8r5GfluKg7JqmVgMmU0EwunQQYoiwFMmlvzTkCwbxA6FebV4twIrGFYeIMSZ9sBiIHpoSfny3KCa5sKuVou8A8X+pGYkYhEPdaQ38CD1PuVRq6oJRtdZ95OHtsAVWOryzGo/QijrDTr/E5guNwAMtVVC5To1E5SjS7yKQfpLCac2vKbn/s+17Sh48LUfGax91n7eElErIo6vMQBwBd7wVwtmFjUuy/wcVeM0CigjNWRhgxNeV4Yy8FbNMwOxd9Gg3qbOMJ6YqpJj6JdQ4hhLrJv0Z+3Ysx567QTb6UTeU3y4QkiGuzeu8rdO8rjR+ZFbwN2+4wNiW4O3EVLDF46ElPavscW06/7hb+bb228xIdg4ulSRAEITW6L66i629kv6PpOlm1hHcpL5cLL5Yy5R8stmlKK5ZG4w7fEqkWFt0R6xWU4ZDZuIuiOcJ83IbWTK29PCh+cosh2xraXtdwxxFEMHZ/hR/ciLy0G4T2DdUF7YSlRwYXF3DLbgHDknHQjuP0EQN/yiSYLDHHxkffHkjE+kcyfkl/emJAc0w/aiI9innsKEZkrB6OnMpXPvixr6/Xnoe+iOC6K22Kn784/ACDnHdeGm0m3Pdx4ir8jkuME6uY9BDREqsg4YoEleMXYX5MIj8asBfkS9EARGF+mdaSFEDu5UHxhzvdQ5Z4ZfYE/0P5i/OISklBMBJqzaAbuxdruaaHXeipHB30YRnd5TeujvtQ+OExmqts7bgiPdNjneFZk6wdtmeiw+u6WBje0vh6hemayJILHdTqOFHZ6aXIdvXKAZDDfFxuTO/HraGLilRKp6NlJxW3UQs7DwaYIV83d0qpbqXFRq72Hs4x71LM1o2CywhyPpiOKZI9LTu1X6/FJc2M7KYcXIVo5mcXplbMtu9b70vJOpHWM97n1KgTy6vFYo1/H6flBIBNrZPrMG5NTeLuTZy1I4Hdo+5mz0WZjil0bogpOW0ZklKnLws421jEaZdCKJ4rJSgYfeCk9k7zjWoZi9nsvprQWrabrwM3LGSJsBOWG5VTQW9fTNGt/FRvEn7q7vwjD2v2FxUB6XO724vXWOs1VVyVus4PNVUkXP329dfUzSQ+O60cdHxQN8g7rJfjsMqjqayKa/gQTXIHpb3cKpV1VXZzMciQHBiQ+PKsAXwlt5a8iMNMR6M97mVw+w4d4HHznuCOd+vz/gN8Nh8O9noB6wQnNS75XW3DHX8ZQsoeU95S2M6fEpXJ0NyIidyRsckHnrY3LPnEGqOj+Q76/HuZDBLyFfIU6Cd8Yd/3Uu2JVD3e8F6DtBloodh7QVWI0VBd1XW6PAFgb0yiCFV9wZihRLPLPPbngiXsQ/Bx5nSOPdXy3DDHEt9gnR1P5+TkYy9L5XIQh1B4yCO777sM5xkt5UcleXArJEaL+DTQEMqn0ZVb9fV13h9+KM83udLlBBFyxaoXIHDPu4ABS1k1QHEVrDNqTp2if6e0RbXvE7c9fQyL6QcDbgAIG1wIKXJwDEOLIiomwcPfVc/VsEV2mn2YUlMlc21auuWZpKNWXGv+MfFWYdIcP2cgyJvLQuRMJnWfF3RPTwfnYMnM8tjsCQm5+76Pqnoj9eppXyb/PdgHevuViQjSC9SN2YtQzbxblNu48UZWmU1f1+cyIhJNwb8BVlpLM5v5iqZkeggyOXGC5KpDc9gJhwFM971f9pJ0Jaq/F3JJP/r4p9cKsSHOTdb/5GfYkyXrGPIjGopcA3PtapTaRbPhsHdW+N4w0azFKBV8s1bnPta/qmrFIGq5g6CvGgqSXIOMaXmqM2uzfLDg2cn3dQSI2xPYysVP+6/G118/85uq1gTQpuo/DtsrRJq8faVJKQZbfIltHc7hlXv4pYvUHKTe5mZS4csdeWrgTLfEVtxETgXUnYqJFW4iJAUmJib0XYhLTEvDPdqn0dinyt9KPSLQ0bZdKXErgtM/EbpxfTy89+IqslR4K77IF3HrcuAcCtPjNm0BH9xt2QTAZsfU+ZpM9EJLwhbnJyYPZXWIPiiw8vxodiiv8R2+/Ww/2/5Rj/R2mvO1Ir/yqE/peRef528bNVnuL3ngheUuvOSGo+aeI3gNRyuJG7tWyEmwfq95/DI9mZs3h2YscFSNNcFtAnbUeqPCpG7g1p87qy2Ws/quGO9onSQCQnMia5rJbrJiqJXP2M3VkwfI0r1a/m6Bmd4VDddgibkGa3mE/FfWxgJfKsV1AGdGyfAL/9HxxzUGhoQYeppRA2kBiT8mwKpn3YfUa7gitsG/OCymcJn0nA3xfXSHXlBg5/n1kmo27H4F14ISP6uAQGyYCk4Tdc/4St108BcaZ1ejHF3fAy8DK+N6GmAiDX9ImNfMD26wsj2zSc75fh2VNiH1qk9JIyuGO+nOY1t6GkpohGlaE32cs0AvEtLX78JRV8yynWAzSj/B0m7XPy8VmYrVDNMwrYquEJBuUg8EUdwGrIXPuLr8I18Ff4BU7sjCBbNhoWickFjxjkmpe7S3CcOJRUOy6akonY6qyjeIo9XBedrxpITPggi19ZG8zvfzFtqLJFXIRl9B60wOo5B1L4j3T66SbOTDyNxzYdjLNZRJsfGNNfgxgXPEAuMeRI5y12vMw4VWq6TkkG3LYUj8lMP7+2hwZa24WTKUKpFK42w3XlsihJ1TTiu71XMragCNWzg8NfJc/2OoFuZ/96EUSWF5qLiUv6QbO1HpUwh99G03rYfzc9P70uHSoBG3cdA6NZ1CbB0c4R5ozyfItd4TOscuUokAldVoXC7luiok/p2eTW2B8KxLQpXXHTngIVcl+WdlsBRkcflKeg4UlArVumFpHPh6xXnSIdLrhazmgJFKfuxZJrj2z+OLJ8myBfCdBpBOyqtuUn1xrYnX2BvgWpINm96uv+9lX7874nyn/C/ffrzwIUOKvoHrsu8trlmbGF30nwHFq9g77g5QkjK+Aypvc2cX5XglO1wSedpB1vvoa+/XkbLxved27X32NyAzAClrf6wyzp46flW4UvLicip6d5K7ZHM4eEGz/r5fsxFHB8ZlYezAjUTqA+NYtmaYnhnmKLbruLvzaMTVJ+HeCJG8kkL3V+m/VVi9+TaVJjHmyKsOdrUi9XdYICiShffI+JUK/xZa5VT7pAiF7f+qenP1y8o7/n+KfBf3rhvAL66h7nPj+D7vtQssb9ViH50yxHoiKMEFZ1ueoF1Tqi/RDtEch9fItAAJHldnoqnmVuBNvUEnHxjii9wRw/WGI+DarcPZdF0jJGCSAERywuwM4MeM9yeDIvsGMiUwdcR/t2ttTc7msIHI3IckE0s9AVPDu9CAiBftGgVskYk83un6uBYvYx9wmoZ0kNGU0l37FRRxoP95O67e8k/PsrQf2NnOpkLKs+wPb+f0MMgljorqQNXJLoQyONYZvG4GZ6rjCBorFsmXd5LJfdVkZ7lgMRLCXGkJha2BBuFU0wUXbZlR/wWjLtZRs7D7n32GXiTbqtKbKn4Kjxn22+FRz9bvbAU0Pv22iTLMr/eA8N9hLnR8TUIkHZGPTdIlgncPQRcu3qRiBYs3RQaQRFifJpI2mp6TLK4XLN1oGFQ2gbcnTBQkBJUN6tQwoD5z15lCNSYyUsBWXx7xl8N4pM0zJxjJf28xvEfzuLvPFQp5zHInvQSCuty6vXJY8lApFNwg65s3OeCg+XG1IG+CWV52KXZegvdnMq9CnyR/is9HLLXjlVKHBGiQJm9ux7KNKEE0BafKLQYxtIPIfTYXEdxFlwSGIqjMRxvHR1CdwFnRsW6iG727WPo5+Mj+S6VS1+79GhIajDZaXkIBcGAhzX6O/2OhQsU0C/AZ6R8irspl6x6ExjggwV8FFXNrlDdA2PCqJhrz5ehD0WIfAB05XtqmYa9XY2VEoJV8B7fpBjpqqem3wE2TP8LFdfKb+3qVwMB0gXyZ0Afv2KWcoY4lxrhIvHwg0DxNDNF9TFx7uAu5rklq78/xdMVQeYIvAg+i17w6HFFu2Thyrkc+WxMzDhmkXMPijo/DtUXkcWmMIIcpjO8oxxwiWImIc0lUVvh6DBrEPHeJ17lXeEdcsAapUUQ2Qi3X2Koib9EZ1nuSBNXDO8XfmpAxxVEXM4du+Pkgf3rL2R2SAt+pd5mbw7e/URLGgFcjEDUdo9o1wgkqsyh3utEliDVnD3psb7D9ZIXWz92rQKfsZRAqjANIapCcILbrr1kRThOC01On55NofdOWj6l34cz9z+hcNGA408EEmGGGsP0qKi1jfx53ieAmOQd1xnqmy+NuEMOt1iiTtopgm1xi78PxIeDAYgRyU8KpmffUeCof3zLrsWdRlodGC7KMVRIU0AW6Ted5buab4G6jXVGeA8URLGC+dX6wmhO0IE+rAOoYucY6mBJ76vqWtJD6rqQHzDRe3gBbnbKpK7gnO3W2IFZop9Wo9TsDf9Gp1F370NcZDXBRClYEqo9gS+GEYCGDmh8prLAn4uTld+6tzOMU8joOsD2F/KA6aVd4Vr5WP4dYA8nBn6fQktdvV5HfQBHuFbpseVhzeRelJ7cmJ2Cll024K4Yh6SCzkhCv8Fi2113F7WLFo0dpn78f4G5B+J1G8BfOsMtAdVcWetjTmb+uhsTGjotPGvtk2296PwBe4ecqkOw64m4it40o11rfp103ehU5WgrVCfZpuvEu+fh2UD2k3vQoiJ0me46kl9JevvDbvwQ/CbAvRqnm78XsOxhv+Le1WpOV1wvNhIwaDDouRoOD0FPGxEbKDnRpjcSuIR/iwTQXF1z0qNUIVE1weZi26j25P1VUVg2NLkvTAOjlZFusJrMwWjsm2rHxa8vWiGlFW8U0W3luT4w6M7Lc06+Y7PtvcxLfvqkTZpkZpMIhg6oLEGJKfzbT9zQN7GM8Qu7QESMC5BgiFw9QSfSehNHlV62lVAyzyWbVElIXVjrUhdgC9pcV2WaGl502dTz/eE3fsvMUWtaK6GNGBnwFLtBCdttHd3qvbB/hBMsJkQ4eIko7wn6HfCSSKyBgYKg0nTRmnYSlMFoPAXS7X79zBfziMoR+zesCiYUHPiGROzrs+KRnoPc2Zz4z4GC+8W8CiQKZIDM6l7lbxQfvIZ2cYrc+Q4ZtI8mzmYGmEEJL5usBNBiH4QGthiRD5PqLfw++9rZvZGOA11An4jy/pZYt4MJFHMd0Pc9ERHDWsuWyyeHcU+iMBZotaxYkhVK+p81BjuMIVFPE2UYLTnAosdRQ62WBnt5jvGbzlxdoGP4zE/E0NhBGY2xq5EU5bKKcZqF22EPbZc7jN/o//QZw22EJ27RdS1MEBg2DNiot8UTvgtWu1708WPuyb3barHeaqxdiOZtM0QFZG4loK4hUlIGEn8je8W0zCjmbDYK1FCFGSAvXtLosiX8+uw0lGd74MEt83/fo+9O5ubh5gKLPcN1AoVJokna1UJ/lUM8sgIqUS4+NijHJANc3QKMnJnMcY7bnBXrLNYoV4jwNNWX2Bp8gDrOfVKIAkKuio0cYWft9J0NUS/MK6691LKWbI7eFIUbFZqMnDrIwZv/dz4Gj/yKkD/6bqwO4+MVLtms847X0x7oIZf3O+yfoUKgwIhmOfsa9dUKIXC7lp/7vuMTQa/e044Ylh7gyrEGlMQGK9rrXFAdRX0+EEPNFXdpeIqZ9m9gIpyvjkL5ZC03GKL4glNiQydgbNpl7Nty4mBW7SleukBy7t4w5isUSOpK2OpGdnbD1hC6+mtcxxPZyku8jXJyVRJpKiSQKn1Rk4DarkYBvuNHKxxcnSJeloBzcaxxni5H7aVuUddGaidGH/MUtq3wSyPQEdrqQ4PR2ppMa/Le1nx2WDQ6H1emR/XYHUqIP+cNTjSx/2KC+h9fEnce1tqiAEgi9ZuOFJyjiH6TDTKPB0UP2sfUiEPjQI+dd62s94IMCs9Yitz+E4PG/gUqxCDbUGt08CUyldQjmaB6rmwxiaTHdzqtuBWwJHfmwDlYJIq2qGWr7kGuPo2aKal3mEDVzFstCGdf6yGX3yxaev77/+y2Z/9LlL5Sqvn99/rq975pq72jRysMqFUFlF+FtoLqXn5d/KYfafJCGr5VM1hx45zOLiS3HRZQ+QIBUs54YlVrWckxgK92hOBBvnaUHCLU5xCn8BtyaWbteZpNw2NCNpkl70Tre0NVMIb9vNtvub6RnTrIxU+Er3vZF6byu8OXzFLhpOhpJqty20BBd4vtQbPF89dq/phTNl8n3D9MbEfMu8Su99GV+BPjwyk9hVw8weYoajlS6NBceVFutCjbd6M6NeACZ+OJq41m9emgHMxdWRXVDrE8r6LrH8r07yX4Y3RVNt1qBrqk1TcIJXSglxRn/Zc+kcp6hzJB+PfZJFvjcNjYwcDFrZ3cAjypXSG7DhBEegmg6c5lOMpBsJEHgchsb205Gs1YgFG9oZjCCz/CwRPWz1+Z5WHVgvKNjln31kHhNDdpyk2SfzhYxfVI4FhPP0ZfbJy1aY1TmbNVshOkhBUhRdCw4RdgsjEkoLdKQ32dbdi3xd5rwBTY/YMp8f2HzT7y29wPxmn+CaU/dWQ4ub/XGixJgzSI8/CWZvqvmgHfI6UiSuzho3BHW2V2U38MgnH+BPLtk2bs/eG42k5pVXZ3evEI41QhzW6YjZHS5nNC/2LVgWKFzxUdPExlRtx2XVdqZNcbrgDj4fOWYStNfDZEabRqjsazdtlhPzVK8j+cqm783TMSJxWoE6PGkBHQJQZHLlkxKK3yjgMVvKvKH9v6m0dbdUzgJorPHPeK8c9T89GQuCWJZSxZFQ0fXQOPbrIq9cGmMHX+p1XsgzLhWwZRcOf+hS9yYMHjOVinyLdvLb9eG9xE6NkHSUtnV0MK9pwmlAHnEOY0TSWYoT5DvU4Y8sqbtE+/jmlVivbk9i/+puSexf3Zpf/g7UGVOgwRFb8q3uSBQEy4CyjdX5og8MZ/R1u+PQfrGgpPwfN4xr6v8voxu8svQCyPCKv/DKAH9fUvphJ2uRMqhQKqyIDbcpr8j1AkFnp2W1muXX41rTwbru4kXfxM1Oi9STmG/SfNM6gJqvzsrCjHhKV72w2UELTiIFQj//stqcPOL3r+WnHmhMMOUv+hvlw+7xN5PxT9Mv2kE7sbPKo1zYrYlBL+G/o5PC7tYo2Y/GgTji5RzOTTJZjKwBv8db1UrxeVPUO2skGoYoFFzNKfX1DhvWeOVnPuNbdV+o24un1bCxbsqn/LoxJfr965afOT9Lu+9EYRPfLxfXQ3oIjHM5Utufl7MpDY/jOTU0Jnzfz1o67/qKP0nSmjI8HnOPwvaOjgOhO8TjYA6aAGVft8MKOn3Hvu5IdXXLcrTOXVIgYQwn1+Hlge0X44T+UFb5FRtX2ddqeRpeqCORG1WSaNDnB3UQUDc8KXEecPc1xBv3Ut1GdXr0VohUV7HtjgglPKKnZik/cjsKLN5ypURSND5VfIftIvl15MDgOFtwPUApWQ9Z0RDYX+92mNp2/4Rad7p7/Wx/1Nt2xrKwfsjeaTqAo70DOp00HKy1/9pKFK8fhPvKIcmCfp1QGYyJhttdzfOr7t5LdxqSFt8/n2HrfRoyA3v/G+U3jP1gTNd4w6UagfuNY7RJQXeHuemlo49J9sLIAFG5MVgAg2SR5Q3nO53y5a376K4MXDbTHTQNO2F8YFOX0XoLUvMumNabXxQPD7NQm3Kz1YiLpKt0mL6505U473ULTt+d5uLrd1rgtNzRs+PD8p5q0lXNYK+W8VjwhkDBV82d5hM2LXmCfVoUU1igo8tvT5dhmtcgw6uLSiLG0bTSRRFDnK+f/aKfCqOKUvFGN1sn3GIuucFVDmYtKp9POV4n97mPG0kHwdIQ9yYTIj3Lp3+j7TBPbz/nsFHfx7vG74S3iwKE3I/ueEcQOux9nmLbD8vZXgssWkG5mEJqxzbNO1yZrtb5sO0dpzbjY4+7VfYgsOtzXicWMnXX4EZ7SXfN4QwVPJs70RwFnnus6C7cNRlROq6kyZpJ51SWvZGsO8kLZxnkcvaGsL6W9R1612ZMdik3004lrpmN3MnaBDrkPXduby7KkJK2qRL61KnbdAlylyNak4y/z9w3UpHEAcQ3hIRlN85I2JPQYyyI7jugz8NOqrJOpqj/HtEuat/bdl35jZ18yJb1JJd6Jb7RYTp1dbRvWUtpXE/0wqGNxn6FmyElg/xqVtSFv11mwNZaaGoLl+l+vqnz2rwDZshmSFxEg0u/BxEu+UrKfHpf5WIm0Xk5I3IqBkVMjKAt8O9dyQ55fJ0KGyc4wb1SMnUOCFIuik8RkYTsnHN3WoAc9xRGRPkq99Rwvvd8NtnMfDdgCqPJBDNC2hm+gpSDNZBksqzj8FQZes4XlHKIqj848Cb2+Sj57g/cyECzodU2i3q5mZwjhFUHkuSSlytk3uBYp2srw7D0y83cNAujHnbF0fCW80s+mRSrWljaluNIzNj4ag2a/uzFAr6NEkVrS1qscVoucfUPjgh8oz3xD+iEqprP01UdxeYSJp1BAPYRLI5H5XPE9yJFHsKEegzNFSdisuQQs4R5LvzIgumRUGeHtejqhPNsyi1Bqxkti0kCyUWoIKK1MarkEtTG0azvZweClt7zcxerlhuYBbXILW1w/W1op5JLQZu3e/sPHH2oj02/WvP3N+SxuKWbbsGw255VV7TQgBRDN6/4R1I1V7h5U2lilA7lIe/gT2dII+siPa09XkeexHwEcgo9XYOkDbc0wSFHq3hg4RF823awH6yr4h3qnfD8ZXs6UFGxwR7jljjEpBTJvakhHJFVpRHZz7SmHWxExvZ3bnNyr81ULmeX/FbdDr50+iF2sCSllDasuhi7EQald3Zar0Qax7fQehE4PV35i86OymNLfBZNUnc0RP52BalcdNTt4CoqltuBJ9DOE+7dx1CsYIcDgHAkWthWkJyUlh2UWBQR8iLK8OP8cHiWDFGSy+JcVz3KBXsOvk624djvSXQQ3g/Ku0F5lDwtLu1QDS0D1a0kwjbygnJ+mO6OpiAILTpFbzt3B6ELyX/ew3dLxUckT0Ud31LyCAGU1s80iZ5LqOJvnIqy5DVuCPCL0ZqbwQGe1p64+rBA/eU/Nc7AraoVLW7t+zQ7N0Tb8lBIKjBXYXWyikcXBxdz8HwaCuSiGhvBdu5Kqnm+2OSzIPCwUTSauIa7ZDC17KG3HVP3m5jaTDQQuHc3ETSKsmBPEcnvoYzl+SirlvjXXwM0XV4uKpyV5iD/aICOsw6aF58H/tbQyNw2zN6u88u37PRW+YsGne+SA2WuJ9SB7vNsQL3o3Ud/zhE1SdTsW2n92RWdVhYsFUp8iKbrNufr9lQW7xnKgf/+OQEcLZsudkYI997NByRRvd3lstWmFqyJOFtCRuwWaXGsXWap9VT0YMNO2vUobr8tOKc9GUkv2RoKZHdaG8aiA63hwm7TNDnUfhuH0mxKiPF2WuXl+jaulGXhtc3Y4+Wi63NHmEb5htucHSL4e5xbwQTt3ABv38PjWwoczEfBxSktix06GFEl2/A+iuFUI1f4Yq+q9v6bNKcrJEKXW0cv+JpP+ARXTSSQyIMm4wBp1SJxWiazH/l6oQ80EbJe1u60FnKaoyOUhFVxEk+im5xeMKGp3tfGIO5KW66qDa/deHt09DZMLAeTzKjjwlkklYEjRtSDuwIOU4ia4HC1axkXAu/v1qZC/xOCvD44OCn5UufLqjvbNoqO1vL3K+qL3SeFxLAELdlEQb58B4osbucfyCAJ6zjpvC9/9E4zNwR3GOCdDClAYtiDQiD9TLkyvTiSD8c9ZNmFuHIcssqANTJ1vJuOu7E9oKIJo6LOc2hVYE6Bw/NZjiNpmAY+i5y2wh1wR1r0JMiGuUw9r0O/a/YAuZ1Gfci+oVBgUMvUyaJSZm3XmwcnRncVpcjqmm1qARW4gwa3Ezo2XZQXpd2RV7HaplhwbJToCSX1xRLAp5tJIdFTgERIZ+xYuuFjoKGq2tRsMEH3xoRNEIFCwTHr+N5Kxc52gpjIhk6ACtuVI/bvQZBFGDXH6BFNhDlL0Vzo+Gz6GwgiR/AoVWLY1d3vYU+ZPn8+6vV+Z4dVLIyJeXCbvYsdCnrtaH2z+/ojzPQNah7cPD3pdr433Yz6Z92k7VWHrhY5pJgedN5LPxsqckJ/fBbMfVMuR8XJtYnqviEpfHfgerjw0L/T04YJhF2nLeHsuNXm9LRkQ0mSEqzn+6Aw79CJiGTGrnGJjc2RFd9Mg7i2NJhEoWb/4/+LrmxlIiQXwdA5tB9sX8nh5l6InyL9619FWXq8yBcExLYUDaxfLSgRVBBErwK3Se7XegkVyljuPLesr+EUePEcGaMa/mrtF2KEmYA9Gfq1f+NVxLyFrEDLzcMKRIk764j4EO/0P0y3D0NDge7DUAmp6pRdboiQxJ65KJ9IFG3eiHYJ6JHUcSqvD7YBdcYQgXF0ZGoLee5lH2ca1pf6xBoSRePwWGMuoFGDgd5guytisAzxwUsVcK6P+hz0URuTyqa3C7a8qxFfjomuukNmOgGPonAyVgly2elyOGFN7Am7k1l5Lah5cuoSN8yfy3XSLDkohF3MmenU7J30hDbuJ8nmfZh1WmekRXmgegM3Be0X9CjmvFrG9nJk5UdGi5PAEAWLlk4RetW8dWbOAhM30OmrRvRQ/vokQGYMAcnhFA0+r6Auikm05qgkcrmjUjwBDopTnfKPoxR6P5MQvDD+mWQKoRkxkF6QtoZtAqymVL5MHPqM9uv62uwS3nbQu1Oj7HsqOQypzu674lqb2V0oJWuZh/AYwQLy1S5EZFxTFU4PZ9Qqi8u4x3gj2TgjOwd3uSexhlXbcWbWuhidac25AnkYrHHPL8cqTd4dCFWKoXgif3cgrk4MKiL2bUenb/NVNzBgAGhfrFO2rYPa8n03BBoW2X3OVi3GO/oRfpIQfP4O/ovTV+NuEPOcbinK/DrFnl/7mZsFHbHMXpuC7M7TmFZOl8Z99+v8vJxB28XjfYXjVGFcytLuv2/jcaM9sa807BQPqM6bzRxGT845pXHimzWdj4ZDUV8gC+HWha+XY3zvA+ujOZQMXzhn4DsxQgfGsznN52iKAuoxn6kj7x0ZPsT9uOTX2HxJUa4+fvE08f9pjuFUplt67o8Z6E8vfCGnjJ3g0MSD+Y2HJTl3fO89uPRUCHk0C+5MG8gCRDcpsYGPF4q3EY/nl0cyNzImoUYn112eGbl3p5eUEQh0+DUJA4FlY51cDjLD7yTKzN0DVMytqI9vaymf3fdCCG4QZAVgP7m+RWbWS+3V9E6wrfjCialtYK3UkOeEuSl2vvDsm7bCEz7kYWlMN9DFabOfPf7xST97jXsq+Y7o4bBnu8MKbtsg3GaUA1KiLrYd812J7Qd62xN+x+1o/tjkbF/9h2CqZtjDTtVmEoN21AN/VYWbpwlmqBvQga17uY0UZO27XDf/Sw1MDvc/nTXmq/razTS61ztiJtAQKcKB3ODGiSltHaYEgDATKxdmINFCuzzO4EzARtMvLT9yVAoLJ2DCQEMrEqhAmnYRQIW+79OXJN0gXERKBNmeUyhHqGJaw78HvI746YIk7MIFPU+6w33aKMPPHvQ0wUPoz4j/qMF1IQp8QD3SZh9lf0e078ccGdnl3wN6tXfssu34ZBcKIlTZvPPTFIoOKi4EcxV9FRQ6dAPGs0tGgf+i0haBd/SOjzqmsdlW5kjTV2773jsONRIdHdTYa3oFoyzIsxHpF2nFfDYpOQWNg3PSLWeq7ceo6KyT+OCkm8sdN6O32w6ZtotinQKnvmLJuVVkVLrYuC5bRcYIlh0w1TU34fZ2OlQ62gYA1BW7vnNfaFdTEOzDijts4LikwRVNxK1D8opom+UtHs42SKctiZ0c5WO7bTDsemSJjY9hu6IP3zKqeb4oT3GVTehEZL93sju6ajhan0JrBFO0iNNbvTPSY+BxKNj+mgQlyFk5OiOrT3UpsbrxrKgnbnqQC16jDnboLjtdVi7qQLCdv6zyGfW9cFyZ9XsWeOCPc4fBQ+NSI1em8S4OlJAi+ls7E2vjwkdz8k4jNbYaxoX7x3G5jbiQrYbCIHw3jDB9TyCyaxnK4vZA38XdAn2bYSjbws22BaZsk8a2hKq0uulo8pQlu/rw6UMqq2B/aL/6kcjkGJhmevGmH9toVAQuzF6Sc19M0I8cQ6MCyYbQvdVvbuI254R2ne8dD6vh3XLu6Ci3o4mkrrOjLM92cMICRdEXlTTm2CirWzTRdMeMYceNxlH5w/h5x9zLMxY/NPkZ+w1oxigj896Zp9pMMPPIJRh68wLtOS1i4rI7PhXI7KlK7XVCwtZRvV3F6dQd7XJPlSRzAoli/9a66qjHi+ouMcZOg9B0Ih0QQASUBM8BXQgeeeEEZLpdG3kE/P6UCsneD3uguzF6JRuO2w53lcDCduJPsoXw00x0sb3UELuz0+MQgSezIl+rK5eIkUjGgQxgbHrOm8ulnu4/4hUbhxMdL3t7F0GoluOqkiTMchEGWAF46w2n3ZXsYJqgUKMWQt98eONzyDudXDezmnlQZZehE5Lw3VVy7Tn1ZlEzwhUVnKIIrMrjFjlgaZaWsykSL2k/I6f5i+yx+vGb7/mfBcbyhK3kNNpicVGulwtIMoNZ+a5gZ3mSPDcnSYjY2/UFa6Avyqo8wSVrE0z+jQ7r6wsRK/gM5ww8EU/3d/Y2XM8TTu5fORbuXwnrVnf4JvcOP7Rx7vi7Z9fh+4RXB5cNRzy38d74qP/QwhqbH5tM0ZdpSgX+W5NbtnxznHH7lcn+i/K8oKjjef6dZ3j+XTu3c1E+6Ym7nY46S9TF0RG/0Yu7gVLujZ38FTkLYJ5GET2R8Jk8CIdviyMiSRyZytn7ZfIOkTGSL5hjGDVM6J++szQ+hvfNa6EIknOzLRQE8TASC3TLnusLRA3YtxuLHFMytXv7dFvuRwGRpgCMgLXck+E8JHRQZYVr5pUSBndkxIBuuB/DhS+1nDoEOUU+2TaYIHaGjRoBFvUiJE1qKoImrvrvg68NjOVCIcaqauAmc/rvatDmrtmm2T15mntbNs5zjo9l8Wm1xMWrU/q1kBtWLOgOlICj+KENyyfrZVUJZY+uE9qSuYL5FxK+QZMZJ6/4ZkMEY4CICPHysfg93xmEWRVXJH1JGO/YfRnrlyRNb9UMZOKXv/xiwQj8GNu+WUHRGflmqyDP8EZT17ikt3r7WmqGDPaPVOJW0uy0PuBy6pv7CPKPVBKPI21V6ncsB1TlL9GLLraRVOXS5t0bre5jql3D5uPLVgxfHsl29j5prQfnX21VfH+tC4pY30imSHYNXEzVSRfhyGzBynE/dlGdL2fsws84xdG7kLmYo8p1ofTNRWNu4BIIxpavxY93tZz5GHw6JW2I0V9LeHO+gAWfkKTC7bWLgq+k4jjnfD05Z8yWtuiky0Hfl9TJ5aW7d1TzvQbkm9+wi6JzfeKAJ1gqvPDWXpkPJ1Scq7O3EFzUONYoAcb34cCAEdZ2MxW27t9KijWd6DgkE4Mdiz9xOb2SW2/t+jzxfLILeNzEYBpagMAP2UAU0vVbAQR37fr1KgSZfPp5DlHH+ruEBI+B4O7SJdr8bK8p5YJJLSm3mke5bHja+26y+36K+o2p6Cfjki2zMK9209YySH8pnIGWiLZ0djl4Rl3pb0w/E9fr2d5LZtoc7O8EzYZhsNLULWmbFtgZf3eT4fKDNgb5KO1mDMhNeSM3gyZk85lddCHVhwVyutV1i6zc/8aVnZaWFAh5vhhJOqZQJNQa09wVa+7pIS/P/ssmNvQsZ25SPh6qSgshpbd+qpBkym57zXl0djTNO7Ddb1rk1Mhn4H3XnHi9wOE0r4M9IW6uEGSZFcoNMgws2mp9NZWKjsNV4u1hDDyglojnFg6eMyzzBhf3X02E7e5ucz63Fjcx9HgauDFZstZbdxnVth6O8xO9gyqGdyRnH4+avPhHcaHjXmR1M8RoKWV43AoFASduSajkcby8cVdNbIhTC7qm48JI8W8SFthfTP8cFvhNcFLUl0WxYJZVTDYcQO4WU93E9fL3aMWq4M5qC4SjRTxZl9MzXE6TuzbkJp3ZdaBZADhWj5wscTn89G85Z6hWXsd3v5c+gQYakgkfCmfW1HVyLSdguRtJJBAH/gpomEi9QyTJ5qwe2oLBVZKYQJt/pYxOH4VBtnJVKymdFUbqbqJi5Z4sKGd8aDLngL3aL7iTY9qim74doXlW8tVwRm0w353TnM+VnThNT0yAnvO+5medIZqZ3N8O6FdglZfRUBPK9DJ/bzjKi1/F4gbVizDT4QqmslLO20KX78x4nfb5xv9CGo5c7UrF/az27gioSfL371izQftjEUBv/XT9ZF5PR8Zhvj7rute995Uh2B3EgmS3sFFlKs6zSgLkk4Gq7AHmFggHSSlxtMobFUXKCOrFZdoEihRy0v8EQJBvPOn3oywtGkkNwaBVbBCsD3hOArGfAsSQl5dSRyrDi5TeXu/5TJJHXATk/3o/eMvl8bbmso7Wu9L1fvTelXeUO0pPfr0Xs63r/YSN7RnDCt5xGWrnEXXCyU2FDwh2sbQOC/uOMjkhDGMNWB1BpHE+it90USh48wuX8S8aDglU/lED6JdNoI9SoF8mQHf8LVBu4v7FBmq4h2gRvtdpy46Q5bzeO0Ihv5LX++6FLqIvoasXlfg4oynPUKwHPxjfc8BG0m1aNfltmkfjphLCDWZAXPC8HhABrpbIT/L8p2++NQf3u3M6PgvbdSdyi8YdeV8rcxuf57PT7Bv6Z0CEAkmU/bcoZw5cXrAzv4t5uRz29NBm16hy6bC2Uz37mF4Ph//1ccv4GvGml0linZgZ9XV+eTZN8cdnfZ7nUi63zlx+l2H2zJenKaRTpUvhxAvk4V9yLicvPA93nIg0Pr08n4e+PHfhejzXfTeVd+JAzZnbalSpx7PilHXslhidXsmQwndMB51UqneYxbnUeW7COngxfmKIiy/67iUo3vgrDXiKoQCpxj+VU76pKbg84b25ocO9WxK4WjnLU9qCdbdAaJntIOfpjRyZFzfld65HntP5wWR7KO2eH9GkxLVbepPwQcYCxwMlUvaSZ1wsoRy2KfbIxUK2pqOkVXHG5o1gZ9UsvfPeletdz/IV68zOoanN2V0J1tb762IAZ2ROaOGFP7ZeLFkGVuCV3o7J4IoKut6yOvft3NMEF9ppSJ7wgoB5SfeibeLhVu7pjnCaEKi3/RDoSxj3MJBHtgzHvcbc0iQEzp6uCksQtOHGmCTkBadiR1LWLSpH12AgzUKPokI6X3JETnwDrRFafMhy7hHBXXqDBAsuBqOQS1pDRzzrAkDQpvPPLPQIjMWODeifo4JRavdb9DAtG/fjTFJRhNOEd9tUNDa5TeVMPMTEH9EtuhJTCIcSqQZKn9TtZ3vik87T+R83a3ba2TZtjoTwJG1TWTnkaE5JMoDItGcsBG7BcbHwSvhgu0l5jvC2uo90Or2xSvY7pCwtMtAiZsyS5Ygjsh2Dkx8he9NfnqW2cjnglyMHTQrRO94StOA4Ine7yQ7FIz0qvI0zxlxRpuF+tq91o9twW2m8uLcHGz1bNJE3JJRCtSKihakJDoeNb6ZC9+carREeYcbsCsOu0KGA7vRdOANEFI620cwr9xqqNAeNsc9yh8H/WQl8S9FYUScc2OYoaC668seNkDHUjzfF0RYMvAH/bhGvXDOGJLdKWJx9eGcrPvl+D1xvtslbila/mtSoB5THC8zu4NVj572BZLyxfkydBe/zhV2S09x7PqkbUx/ASA7g60h12qoMPusZfIQ4RyqA8s3Ba1XrNQ4G1dTdz5O9eSqZN0SLKOJ5WfsAXA3AIXm8ER4SieIhcD8jqrJJYrZdSEW0KogZGtfLMb425Fv7uP+g8d19ezDa/u2zZj2Gi2il8Q+JXC5vX+PosTcafzHivuI+wHEVFXv845MmQNWrjh/DyafeTIsG6rUeAf45En8/vJopcIKmVelnb169Tg8HgXYceGinUg2zMx30qdwMEFwD4S2MqP9GbiJENgzJrRZm3MXx8LzI6zkJZYzarDdfirpbvPI2VWHaa94JTjeNuGcFoTaSgRha1xv2Y0BOEn8ZmJHHG46HTbcsd9XE+xsf42OhJbGg+dUTzrifpSe+zK7iCk1xdtwhAavBXf6pVxN+6W8mvOHejLQLzXRjPl0Z/sbVd9ryiaVAd0ws2zYwd0VG87zhwksu9O7HG6q3Hke8mddrGd2ro4ZBlo6T0RnUyS+oPXixOPV1rRV/hAwLCYKYlDDkVG3SQrxaXC68/lXqHQYL0WY4PuJyxy0LGKFngDXSkENSL6wxII+ruggK32iGeYTATOPUb7geED5kxuXYKSTDxsV1zsVsRofQIHqaOYMSkdjPpM7fwdoV2Eii5IFmDVsgv1agdusH7I69SNnt8HHq4dEOrlW311aXGxRPnMSTRDUAzPDaOK9u5zagah7bV3cYd/26HL/YcXlhVkieV38nj6e2ka0MTqtlJYz2seV9y9juKkboUjSJcseV5WM2g7h3HetOn99/Puolibg9UxeRI/CVDlefFSK4dTbQ5/HKBKGuUaY+XrK70ex+PEaL1Gql5A70jVag5vk0jodFya2qMaM6rqkm3frNBMsJubAmxjqZq0QF46TyzbwbVeMcxanP3ZZBvsdVTndnPne9s+m2rrXn0vRcKU6ZuY1/RSke46yYbkimCDXfivaskKIPTO+YYodyUARnkXcBECyZN5wI/K28DE+94WDxJwLMomG5GPAVG0wOxN9gkkN6mvHlTU5uCm48mRFIBkf9X8ZCXEQ3gkETNeRFvzoK5+G4rUj04qhFm4r3vWM7TIfb1RbbFdUEraP4LW54lQ8+5UEEvXnbb2QgS8Fz/oZbAPzO5m20bdplAnZnNL9h0i6X63dwwPNS1FXU79aJjBZMBuZni5PBnVRxT+mgvlcMvtjeKl+EG1e5Hz3frVd+buMmYgS8ZaXvMidbu7Z9etxqxskW49buTrvUxMJe443TSwKVuOANuWhTaw1+77KQZh9EgLNP/NSCm6HdpGFP4Uq3Dy4+/QdPftS1nsTd6NEnvIRS1w7NuTs0awkdg18Ztb1AAiN69hHrEpgcurwjWff18gRhhkQeK3UIqev9B6xlqx+M5O9nD0TLWhVOBTfJx8WiltQxwSmZ70ZudpC3xic2Mx9mT4tTIJnTkElG0dXeKMSr7ONsNKRXf2gZspSXfkbl6dX28g+a8B/cBP+zJvzPboL/RRP+F1vg60S8kbMVhD+vpmenN/Ovq8V/2SzZkV8a36RDayq6+zb33jI4YEEL7428munQFnury3LbNlshY0fTd9KDFsgy3e8H+cHoLpA/+w19/uyWPuuKePcIt0mk0b3fMJq9u4zmi98A+ItbADv7B6dYyLs1ZFm9RRKPX8gj2oafTpAuMdrl+A7XkpF+5O3tChjDC2lBdzT8vBXjwexcVWdMUPPwq6KG7JZHWC95f4HnPpPuhSrXLolqFD6F3jDDLawMars/KjsXVmxWlhsM2Uo8b/NOVcOw+afaXqMerYqfl+WUdWoqVIrnKs7X2WbhWgz9j9n7KhCP5bDCFk6gQOTb5st+LMbTfDMZS8ejg0x95EqmZuUjB/C46Z6VyAzbgCTFmhDDVoBCASTptFl6IQm4MTjmOK2dWBxAOnIFj5HMJwAafmEA5wUMHAxg65gicB8rCm8bWVi45yxfOBdaX/9gjTakqYABSkqdxdqGrxsu5PLYeKEMYNuohW/ejwpyi64tnDkvz+cJA2b3H+fJxd/srczLCSvI40qwT42d5fN92XnistXgsNsZZuqr1O7F1G1dOaaDWz2YvOgndibXHJPaAPgubAe+M+HMNYo5s4Mvn8xps0pglRC8jaY6Ku+NFIZHN9mCvMDXsAaF0uM2ixD4fps1CFy7zRIEntu0AjmCfosZyDiEmYFor91gAWpii7MFYcLa7UDR8rWahNIFazEN+dW5g3GIdvmOtws9iXPQqzbX36ACAS647iD2q5zWpiSUsqwkzLN7mOZ77j7UxTS6MKFcqE+TGkNDfaZkLZAY/8c+w5am13LKRi7GodX/CzWwsZY0mIAtWtIiJ4lZ1aStGQ+TGNt+OBu3m7VCe8m03mUBkKWBF+y5fHrtL33STK4cyKRyg67ZjiZmoCLJrSyVxeX25b4RBgXm64rAmtEBlE4S3uMgtgYEB8TCHyRDebZ1DnUQOy5VwpjDob9018S4RsNLYpazKAdZmImxVXHd5QtjtNbx8Q0a64aa2uuZ4myG8FnzyxoVsGtUrMHGR6R6/NLy8UXIYVVcMrexps3TSbbkjDuWK/5DbPwTOcwWLpTcr0UlukRn1/Ff5Ci/mdXjKTPNRq7PcEG0rR+rwmtIlSws3BYL0lWGbgoa1RS2txsUxck+ei9NvS7Et8Ih+ZKoI0ET//1TTlnJHl/A0alD9Em3w+sTJFTu6JS7n9xukOrUVY5TL7teIYGutMBWXr5KTubMFeFUzlqoH4IVfeqHmg4Z2eYlK8AlLrLla4LVDB/EjInIpM0gk+WbzbzqiiDnXketwJL4pZlKwtb7GUfVMYPgq68JE+RWVk9Kp8Vpjhwbmj8GKynZZVJ6XkWk3APwVLxqvQC2rc55cSUN0r9T1NNejGU7aF9CQipdCAxGOlx9vxObhSLS5DMIjkLSBGzOpRcnRAKmg9N1SVwD2YUk/52ZROdYJ72N28ZbuQa/PMw+twZ9v18/wc+v1rhUknYV/x0SXE4e/Ek/Q87+zpui3u80xPkPs6/1aulsLlcdM1DEeHPYa9JSRTSvwJUBm2KMBnrdhRtkmhxEaoUWK+VOlZcj+JYiyUEhWSTcnT+c8kbe/U7h4Z8gKzgkC3oVdGO4o4m8b+LGt3DcZ3yNbHZZiEpsVV4s6/jKRX/X4/8Sdqt34/nfAevF2FUYCcvOixmnZJZz+XR4kXPOKOmBJS3KMW775Dkwf4XbAE8wXyZ7aGRWvxERcx/AQTs7mp9au0N06+3B4Vs5T1iW5JAF8gefD7mFR+Jjwh+FyppnorRkCP11wddeMX5AlBNEk52ttLfLt1nPaHHzK1yu5nFTCyB8QH7eq9wEtOAiIsHZ9gq3CNC2bcjMnRjuKKmeGSN0eCofZMkr2lAumXK/cRDr/GUz+uSLT5/ff47kI5+POjecbDr2NvNvu0n9Xqd5zuoECvvKCuiRqMN/49fJyafznamluCx/zOSjrxcdhDpPaZ+gAgdgydsuKt/fp6W3KnKY67yhAzkVhRIPb9IBZfQ8mn7mq8nxqYM/SAdTz3CZcyc+qeqxs2PGC37uNor5A60ruf/go22l/THXlX4w2lr6sybsz5qwbzz3dvZG//5v/+0LaoMXEMUa03PfgYpPnP7+SHrchiSClMExNN1tNaZ3227jj79nv23fYAzabTDpxW/dYDzj/5N32d5oQMvW+999u0Fn/uMCcvqTDWT9O+891qEN6uXg1VIxs1z+n7Mb01n6n7M7s2/UtbfmLJNQkGgCFWc2oHfqK0cMrFyvizM6TK+9oSObzJaTdyS+scxEZ8VxMdWs4W5rSQyBZNljrSTfZSVaugH8Zj/O+JQl5UTTzW6PGkQgnlVHuJibACBTu8YSTLr0TGd3BBwJOILG+bP+f1BLAQIUAxQAAAAIAAAAIVyskSrinQEAAHICAAAOAAAAAAAAAAAAAACkgQAAAABwcm92ZW5hbmNlLnRzdlBLAQIUAxQAAAAIAAAAIVzIkk5sGBIAAOwyAAAIAAAAAAAAAAAAAACkgckBAABydW5uZXIuUlBLAQIUAxQAAAAIAAAAIVx6uMRcdQIAACwEAAAQAAAAAAAAAAAAAACkgQcUAAB1cHN0cmVhbS9MSUNFTlNFUEsBAhQDFAAAAAgAAAAhXLtq1+UmKAAADJ0AABUAAAAAAAAAAAAAAKSBqhYAAHVwc3RyZWFtL3Jlc3BvbmRlcnMuUlBLAQIUAxQAAAAIAAAAIVysQTekZRQAAKZPAAAdAAAAAAAAAAAAAACkgQM/AAB1cHN0cmVhbS9yZXNwb25zZV9wcm90b2NvbHMuUlBLAQIUAxQAAAAIAAAAIVyH4SqcCWMAAK5VAQAQAAAAAAAAAAAAAACkgaNTAAB1cHN0cmVhbS91dGlscy5SUEsFBgAAAAAGAAYAfAEAANq2AAAAAA=="
}
// END GENERATED SIMPLECA BUNDLE

// BEGIN GENERATED WEB BUNDLE
class CalciumWebBundle {
    static final String SHA256 = "d2c1100aa99d6c77fd9da95e4a050c66d03f7d2f3405bfa9fd22564244b3b690"
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

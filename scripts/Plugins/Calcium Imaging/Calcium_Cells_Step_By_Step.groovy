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
    static final String SHA256 = "9053640b8aa233796fca2acc9e3cffe767a1c16bbb2579a2019b5639cb3d2014"
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

// Calcium imaging: cells step by step                                                    version 2026.11.20.2
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
    static double SUGGEST_BRIGHT = 12    // automatic outlines brighter than this on average are shown in orange: probably precipitate
    // -------------------------------------------------------------------

    ImagePlus movie
    int n, w, h
    float[] rest                 // average of all frames (raw)
    float[] flat                 // the all-frames picture (lamp light removed): used to rank pixels by how empty they are
    float[] respDiff             // mean of the response frames minus mean of the frames before (raw gray levels); null without a response
    int lastClear = 0, lastDim = 0, lastLifted = 0; double lastMeanRise = 0
    static double KEEP_LEVEL = 1.0       // "every cell kept" table: a cell whose lowest 20-frame average (before the response) is below this is lifted to it
    ImagePlus allPic, respPic, colourPic
    boolean[] bright             // widened bright objects
    boolean found = false; int onset, p0, p1, b0, b1; double share; String why = ""
    List<String> log = []
    Roi[] alsoAvoid = new Roi[0]         // outlines that are not measured (cells at the edge that are left out) but must never serve as background
    String edgeNote = ""                 // what was done about cells at the picture's edge, for the record of steps

    void say(String s) { log.add(s) }

    /** Steps 1-2: pictures. The movie is only read. */
    void prepare(ImagePlus imp, boolean lookForResponse) {
        movie = imp; n = imp.getStackSize(); w = imp.getWidth(); h = imp.getHeight()
        say("Calcium imaging plugin, version 2026.11.20.2")
        say("Recording: " + imp.getTitle() + " (" + n + " frames, " + w + " x " + h + ", " + imp.getBitDepth() + "-bit)")
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
        if (lookForResponse) findResponseFrames()
        if (found) {
            ImagePlus pk = ZProjector.run(imp, "sum", p0, p1), pre = ZProjector.run(imp, "sum", b0, b1)
            float[] x = (float[]) pk.getProcessor().convertToFloat().getPixels(), y = (float[]) pre.getProcessor().convertToFloat().getPixels()
            float f = (float) ((p1 - p0 + 1) / (double) (b1 - b0 + 1))
            float[] r = new float[x.length]
            for (int i = 0; i < r.length; i++) r[i] = y[i] > 0 ? (float) (x[i] / y[i] / f) : 1f      // Process > Image Calculator, Divide
            respDiff = new float[x.length]
            for (int i = 0; i < r.length; i++) respDiff[i] = (float) (x[i] / (p1 - p0 + 1) - y[i] / (b1 - b0 + 1))
            respPic = new ImagePlus("response_" + imp.getTitle(), new FloatProcessor(w, h, r))
            IJ.run(respPic, "Gaussian Blur...", "sigma=" + BLUR_RESPONSE)
            respPic.setDisplayRange(RESP_MIN, RESP_MAX)
            say("Response: brightness rises near frame " + onset + " in " + Math.round(100 * share) + " % of the field. Response picture = frames " + p0 + "-" + p1 + " divided by frames " + b0 + "-" + b1 + ", Gaussian Blur sigma=" + BLUR_RESPONSE + ", shown " + RESP_MIN + "-" + RESP_MAX)
        } else if (lookForResponse) say("Response: none found (" + why + "). The colour picture has no magenta.")
        ImageStack st = new ImageStack(w, h)                                         // Image > Color > Merge Channels
        st.addSlice("very bright", allPic.getProcessor().duplicate())
        st.addSlice("all frames", allPic.getProcessor().duplicate())
        if (found) st.addSlice("response", respPic.getProcessor().duplicate())
        ImagePlus tmp = new ImagePlus("colour_" + imp.getTitle(), st)
        tmp.setDimensions(st.getSize(), 1, 1)
        CompositeImage ci = new CompositeImage(tmp, IJ.COMPOSITE)
        ci.setChannelLut(LUT.createLutFromColor(Color.red), 1); ci.setC(1); ci.setDisplayRange(BRIGHT_MIN, BRIGHT_MAX)
        ci.setChannelLut(LUT.createLutFromColor(Color.green), 2); ci.setC(2); ci.setDisplayRange(0, SHOW_MAX)
        if (found) { ci.setChannelLut(LUT.createLutFromColor(Color.magenta), 3); ci.setC(3); ci.setDisplayRange(RESP_MIN, RESP_MAX) }
        ci.setC(1); ci.setCalibration(imp.getCalibration().copy())
        colourPic = ci
        say("Colour picture: red = all-frames picture shown " + BRIGHT_MIN + "-" + BRIGHT_MAX + " (precipitate turns yellow), green = all-frames picture 0-" + SHOW_MAX + (found ? ", magenta = response" : ""))
    }

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
     *  On four reviewed recordings (1396 outlines) the share accepted as cells was 83 %, 80 %, 41 %, 12 % and 2 %. */
    double[] grade(Roi roi) {
        int[] px = pixelsOf(roi); int a = px.length
        if (a == 0) return [4d, 0d, 0d, 0d, 0d, 0d, 0d, 0d] as double[]
        boolean[] in = new boolean[w * h]; for (int i : px) in[i] = true
        boolean[] near = new boolean[w * h]; for (int i : pixelsOf(RoiEnlarger.enlarge(roi, 3))) near[i] = true
        double s = 0, sb = 0, sr = 0; int nb = 0; boolean edge = false
        for (int i : px) { s += flat[i]; if (flat[i] > BRIGHT_ABOVE) nb++; if (respDiff != null) sr += respDiff[i]; int x = i % w, y = (int) (i / w); if (x <= 1 || y <= 1 || x >= w - 2 || y >= h - 2) edge = true }
        double mean = s / a, brightShare = nb / (double) a, riseIn = sr / a
        double ringSum = 0, ringRise = 0; int nr = 0
        for (int i : pixelsOf(RoiEnlarger.enlarge(roi, 12))) if (!near[i]) { ringSum += flat[i]; if (respDiff != null) ringRise += respDiff[i]; nr++ }
        double contrast = nr > 0 ? mean - ringSum / nr : 0, riseRing = nr > 0 ? ringRise / nr : 0
        boolean own = respDiff != null && riseIn > 0.5d && riseIn > 1.25d * Math.max(riseRing, 0.05d)
        int g = 2
        boolean precip = mean > SUGGEST_BRIGHT || brightShare >= 0.25d
        if (precip) g = 3
        else {
            if (brightShare < 0.02d && contrast > 0.2d && !edge && mean <= BRIGHT_ABOVE) g = 1
            if (own && brightShare < 0.10d && mean <= BRIGHT_ABOVE) g = 0
        }
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

    ImagePlus outlinePicture() {
        float[] g = (float[]) allPic.getProcessor().getPixels()
        float[] m = found ? (float[]) respPic.getProcessor().getPixels() : null
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
        // The same table with every cell kept: a cell that is not brighter than its background cannot be given a dF/F0 and
        // SimpleCa2+ leaves it out, which shrinks the number of cells and skews percent responding. Here such a cell is lifted
        // by a constant so that its lowest 20-frame average is KEEP_LEVEL. Its rises and falls in gray levels are unchanged,
        // so whether it responds can still be judged; the SIZE of its dF/F0 means nothing.
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
        IJ.saveAsTiff(allPic, prefix + "picture_all_frames.tif"); IJ.saveAsTiff(colourPic, prefix + "picture_colour.tif")
        say("Clear rise = the highest 5-frame average after frames " + s0 + "-" + s1 + " is at least 5 standard deviations of those frames above their mean: " + lastClear + " of " + nc + " cells. Mean rise " + String.format(Locale.US, "%.2f", lastMeanRise) + " gray levels.")
        if (lastDim > 0) say(lastDim + " cells are not brighter than their own background in frames 1-20 (resting value at or below zero). dF/F0 cannot be worked out for them; SimpleCa2+ leaves such cells out. They are marked in the summary table.")
        say("Table for_SimpleCa_every_cell_kept.csv: identical to the corrected table except that " + lastLifted + " cells whose lowest 20-frame average before the response was below " + KEEP_LEVEL + " gray level were lifted by a constant to reach it (see the last column of the summary table). Use it when counting how many cells respond; do not use the size of dF/F0 of the lifted cells.")
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

/** Release lookup sends no recording data. Installation requires an explicit button click. */
class CalciumRelease {
    static final String VERSION = "2026.11.20.2"
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
    synchronized void report(int value, String what, String extra = "") {
        checkCancelled()
        percent = Math.max(percent, Math.min(99, value))
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
    boolean primary = false
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
            if (p.exitValue() != 0) throw new IOException('Analysis engine returned exit code '+p.exitValue()+'. See '+logFile.name+'.')
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
        File marker = new File(home,'ready-'+platform()+'.txt')
        if (!marker.isFile()) return null
        String name = marker.getText('UTF-8').trim()
        if (!(name ==~ /r-4\.5\.2-[a-f0-9-]+/)) return null
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
            script.setText('stopifnot(as.character(getRversion()) == "'+R_VERSION+'")\ncat("R ready\\n")\n','UTF-8')
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
                catch (Exception broken) { new File(home,'ready-'+platform()+'.txt').delete() }
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
            File env = new File(home,'r-'+R_VERSION+'-'+UUID.randomUUID())
            job.progress(30,'Installing R (first run only). Download details are in the setup log.')
            job.run([executable().path,'--no-rc','-r',new File(home,'packages').path,'create','--yes','--override-channels','--strict-channel-priority','-c','https://conda.anaconda.org/conda-forge','-p',env.path,'r-base='+R_VERSION],home,1800)
            job.progress(90,'Checking the installed analysis engine')
            probe(env,job)
            job.check()
            File markerTemp=File.createTempFile('ready-', '.tmp', home); markerTemp.setText(env.name,'UTF-8')
            java.nio.file.Files.move(markerTemp.toPath(),new File(home,'ready-'+platform+'.txt').toPath(),java.nio.file.StandardCopyOption.ATOMIC_MOVE,java.nio.file.StandardCopyOption.REPLACE_EXISTING)
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
    static List<List<String>> read(File file) {
        def rows=[]; def row=[]; def cell=new StringBuilder(); boolean quoted=false
        String data=file.getText('UTF-8')
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

/** Native Fiji trace view; missing samples remain gaps. */
/** Fast pixel work for showing movie frames next to the traces. Display only: nothing here is measured or saved. */
@CompileStatic
class CalciumFrameView {
    /** Mean of frames first..last (counted from 1, both included), pixel by pixel: the resting picture. */
    static float[] mean(ImageStack stack, int first, int last) {
        int n = stack.getWidth() * stack.getHeight(), count = 0; double[] sum = new double[n]
        for (int f = Math.max(1, first); f <= Math.min(stack.getSize(), last); f++) {
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
    static void paint(ImageProcessor ip, float[] rest, double lo, double hi, int[] rgb, float[] a, float[] b, boolean smooth) {
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
        float low = (float) lo, scale = (float) (255d / Math.max(1e-9d, hi - lo))
        for (int i = 0; i < n; i++) {
            int g = (int) ((shown[i] - low) * scale)
            if (g < 0) g = 0 else if (g > 255) g = 255
            rgb[i] = (g << 16) | (g << 8) | g
        }
    }
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
        t=30+16*(line+1);h=Math.max(1,getHeight()-t-46)
        def values=[]; mean.eachWithIndex {v,i->if(Double.isFinite(v)) {values<<v;if(i<sem.length&&Double.isFinite(sem[i])) {values<<v-sem[i];values<<v+sem[i]}}}
        selected.each {if(Double.isFinite(it)) values<<it}; values<<0d
        // With every cell drawn, the axis covers 99 % of all values (plus the mean and the selected cell); rarer extremes run off the edge.
        if(many) {values<<CalciumFrameView.share(all,0.005d,false);values<<CalciumFrameView.share(all,0.995d,false)}
        lo=values.min() as double;hi=values.max() as double
        if(hi-lo<0.001) hi=lo+1
        double pad=(hi-lo)*0.08;lo-=pad;hi+=pad
        double t0=time[0],t1=time[-1]
        [[baselineStart,baselineEnd,new Color(0xEEF1F4),'baseline'],[stimulusStart,stimulusEnd,new Color(0xFFF4E3),'stimulus']].each { band ->
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
        def path={double[] y->def p=new java.awt.geom.Path2D.Double();boolean pen=false
            for(int i=0;i<Math.min(time.length,y.length);i++){if(Double.isFinite(y[i])){double px=l+(time[i]-t0)/(t1-t0)*w,py=t+(hi-y[i])/(hi-lo)*h;if(pen)p.lineTo(px,py) else p.moveTo(px,py);pen=true} else pen=false};p}
        if(many) {
            int alpha=all.length>120?55:all.length>40?80:120
            g.setStroke(new java.awt.BasicStroke(0.9f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND))
            // non-responders first, so the responders are not hidden under them
            [false,true].each{boolean responders->all.eachWithIndex{double[] y,int i->String call=i<allCalls.size()?allCalls[i]:null
                if((call=='Responded')==responders){Color c=CalciumCalls.of(call);g.setColor(new Color(c.getRed(),c.getGreen(),c.getBlue(),alpha));g.draw(path(y))}}}
        }
        g.setColor(new Color(0,114,178,many?70:50))
        for(int i=1;i<time.length;i++) if(i<mean.length&&i<sem.length&&[mean[i-1],mean[i],sem[i-1],sem[i]].every{Double.isFinite((double)it)}) {
            g.fillPolygon([X(time[i-1]),X(time[i]),X(time[i]),X(time[i-1])] as int[],[Y(mean[i-1]-sem[i-1]),Y(mean[i]-sem[i]),Y(mean[i]+sem[i]),Y(mean[i-1]+sem[i-1])] as int[],4)
        }
        g.setColor(Color.white);g.setStroke(new java.awt.BasicStroke(4.2f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));if(many)g.draw(path(mean))
        g.setColor(CalciumTheme.BLUE);g.setStroke(new java.awt.BasicStroke(2.4f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));g.draw(path(mean))
        if(selected.length>0){g.setColor(Color.white);g.setStroke(new java.awt.BasicStroke(3.6f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));g.draw(path(selected))
            g.setColor(CalciumTheme.INK);g.setStroke(new java.awt.BasicStroke(1.7f,java.awt.BasicStroke.CAP_ROUND,java.awt.BasicStroke.JOIN_ROUND));g.draw(path(selected))}
        g.setClip(clip);g.setStroke(new java.awt.BasicStroke(1f))
        // the key, then the axis names
        entries.eachWithIndex{entry,int i->int kx=l+(places[i][0] as int),ky=20+16*(places[i][1] as int);String text=entry[2] as String
            if(entry[0]!=null){g.setColor(entry[0] as Color);g.setStroke(new java.awt.BasicStroke(entry[1] as float));g.drawLine(kx,ky,kx+16,ky);g.setStroke(new java.awt.BasicStroke(1f));g.setColor(CalciumTheme.INK);g.drawString(text,kx+21,ky+4)}
            else {g.setColor(CalciumTheme.WARN);g.drawString(CalciumTheme.fit(text,fm,l+w-kx),kx,ky+4)}}
        g.setColor(CalciumTheme.MUTED);String xl='Time (s)';g.drawString(xl,l+(int)((w-fm.stringWidth(xl))/2),getHeight()-8)
        def turned=(java.awt.Graphics2D)g.create();try{turned.rotate(-Math.PI/2);String yl='ΔF/F₀';turned.drawString(yl,-(t+(int)((h+fm.stringWidth(yl))/2)),14)}finally{turned.dispose()}
    }
}
class CalciumAnalysisHeatmap extends javax.swing.JPanel {
    static final int TOP=46, BOTTOM=26, LEFT=78, ROW=15
    double[][] traces=[]; List<String> names=[]; Closure choose={int i->}; double[] time=[]
    List<String> calls=[]            // the response call of each row, for the coloured strip on the left
    double stimulusStart=Double.NaN, stimulusEnd=Double.NaN, baselineStart=Double.NaN, baselineEnd=Double.NaN
    String selected=''
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
            for(int c=0;c<names.size();c++) {
                int y=TOP+(int)(h*c/names.size()),y2=TOP+(int)(h*(c+1)/names.size())
                String call=c<calls.size()?calls[c]:null
                if(call!=null){g.setColor(CalciumCalls.of(call));g.fillRect(4,y,4,Math.max(1,y2-y-1))}
                if(y2-y>=11){boolean chosen=names[c]==selected;g.setFont(CalciumTheme.font(chosen?java.awt.Font.BOLD:java.awt.Font.PLAIN,-2f));g.setColor(chosen?CalciumTheme.INK:CalciumTheme.MUTED)
                    g.drawString(names[c],14,y+(int)((y2-y+small.getAscent())/2)-1)}
                for(int f=0;f<time.length;f++) {
                    double v=traces[c][f];int fade=Double.isFinite(v)?(int)(255*(1-Math.min(1,Math.abs(v)/scale))):0
                    g.setColor(!Double.isFinite(v)?new Color(0xD9DEE3):v>=0?new Color(255,fade,fade):new Color(fade,fade,255))
                    int x=LEFT+(int)(w*f/time.length),x2=LEFT+(int)(w*(f+1)/time.length)
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

/** Step 4 uses the measured ROI snapshot, never the intensity-lifted export. */
class CalciumAnalysisWindow extends javax.swing.JFrame {
    Map measurement
    CalciumRRuntime runtime = new CalciumRRuntime()
    Closure performAnalysis = {rt,input,settings,out,job->CalciumSimpleCaEngine.analyze(rt,input,settings,out,job)}
    CalciumLocalJob job
    boolean busy=false, closed=false
    File results
    Map<String,javax.swing.JTextField> fields=[:]
    Map<String,List<String>> exact=[:]   // box name -> [what is shown, the exact value]: an untouched box still sends the exact value
    def setup=CalciumTheme.button('Set up analysis engine')
    def runButton=CalciumTheme.button('Analyze responses',true)
    def cancelButton=CalciumTheme.button('Cancel')
    def logButton=CalciumTheme.button('View log')
    def folderButton=CalciumTheme.button('Open analysis folder')
    def savePlots=CalciumTheme.button('Save plots as PNG')
    def settingsButton=CalciumTheme.button('Hide settings')
    def previousCell=CalciumTheme.button('‹'), nextCell=CalciumTheme.button('›')
    def status=new javax.swing.JTextArea('Check the settings, tick the box, then analyze.',3,20)
    def progress=new javax.swing.JProgressBar(0,100)
    def confirmed=new javax.swing.JCheckBox("<html><body style='width:205px'>I checked the frame rate, baseline and stimulus times</body></html>",false)
    def kcl=new javax.swing.JCheckBox('Include a KCl control',false)
    def cells=new javax.swing.JComboBox<String>()
    def group=new javax.swing.JComboBox<String>(['All cells'] as String[])       // which cells are listed and drawn
    def allTraces=new javax.swing.JCheckBox('Every trace',true)
    def cellStatus=CalciumTheme.label(' ',CalciumTheme.MUTED)
    def plot=new CalciumAnalysisPlot()
    def heatmap=new CalciumAnalysisHeatmap()
    def heatmapOrder=new javax.swing.JComboBox<String>(['Responders first, strongest on top','As measured (Cell001, Cell002, …)'] as String[])
    def cellPicture=new CalciumCellPicture()
    def responding=new CalciumStat('responding',CalciumTheme.BLUE), registered=new CalciumStat('cells registered'), excluded=new CalciumStat('excluded'), unscorable=new CalciumStat('unscorable')
    def metrics=new javax.swing.JTable()
    def calls=new javax.swing.JTable()
    def counts=new javax.swing.JTable()
    def tabs=new javax.swing.JTabbedPane()
    Map<String,double[]> traces=[:]
    Map<String,Map<String,String>> callRows=[:]   // cell -> its row of the saved Cell response calls table
    // the movie beside the traces
    ImagePlus movie; java.awt.image.BufferedImage snapshot, frameImage; float[] rest, work, work2; double[] rawRange, changeRange
    int frame=1, restFirst=1, restLast=1
    def view=new javax.swing.JComboBox<String>(['Colour picture'] as String[])
    def outlinesBox=new javax.swing.JCheckBox('All outlines',true)
    def zoomButton=CalciumTheme.button('Zoom to cell'), wholeButton=CalciumTheme.button('Whole picture')
    def playButton=CalciumTheme.button('▶ Play')
    def slider=new javax.swing.JSlider(1,2,1)
    def timeLabel=CalciumTheme.label(' ',CalciumTheme.MUTED)
    javax.swing.Timer clock, player
    private boolean listing=false
    long started
    String stage=''
    CalciumAnalysisWindow(Map data) {
        super('Step 4 · SimpleCa²⁺ analysis in Fiji');measurement=new LinkedHashMap(data)
        if(data.picture!=null) {
            // A private copy of the colour picture: selecting a cell here never touches the Step 2 picture or its outlines.
            snapshot=(data.picture as ImagePlus).getBufferedImage()
            measurement.picture=new ImagePlus('Analysis cells · '+data.recording,new ij.process.ColorProcessor(snapshot))
            cellPicture.image=snapshot;cellPicture.rois=(data.rois?:new Roi[0]) as Roi[];cellPicture.names=(data.names?:[]) as List<String>
        }
        int n=data.frames as int
        // The recording itself is only read, to show the frame the time cursor is on. It must match the picture the outlines were drawn on.
        if(data.movie instanceof ImagePlus&&snapshot!=null&&(data.movie as ImagePlus).getWidth()==snapshot.getWidth()&&(data.movie as ImagePlus).getHeight()==snapshot.getHeight()&&(data.movie as ImagePlus).getStackSize()>=n) {
            movie=data.movie as ImagePlus;view.addItem('Movie: change from rest');view.addItem('Movie: as recorded')
        }
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
        def root=new javax.swing.JPanel(new java.awt.BorderLayout(12,12));root.setBorder(javax.swing.BorderFactory.createEmptyBorder(12,12,12,12));root.setBackground(CalciumTheme.PAGE)
        int end=Math.min(n-1,Math.max(3,(data.baselineEnd?:20) as int))
        restFirst=(data.baselineStart?:1) as int;restLast=end
        double rate=(data.hz?:Double.NaN) as double; boolean timed=Double.isFinite(rate)&&rate>0
        def box={String key,Object value,String tip->
            String full=value==null?'':''+value, shown=(value instanceof Double||value instanceof BigDecimal||value instanceof Float)?CalciumTheme.shortNumber(value as double):full
            def f=new javax.swing.JTextField(shown,5);f.setFont(CalciumTheme.font());f.setToolTipText(shown!=full?tip+' (exact value: '+full+')':tip)
            if(shown!=full)exact[key]=[shown,full]
            fields[key]=f;f}
        def form={->def p=new javax.swing.JPanel(new java.awt.GridBagLayout());p.setOpaque(false);p.setAlignmentX(0f);p}
        def line={javax.swing.JPanel p,String label,java.awt.Component... parts->
            def c=new java.awt.GridBagConstraints();c.gridy=p.getComponentCount();c.gridx=0;c.anchor=java.awt.GridBagConstraints.WEST;c.insets=new java.awt.Insets(3,0,3,8);c.weightx=1d;c.fill=java.awt.GridBagConstraints.HORIZONTAL
            def name=CalciumTheme.label(label);p.add(name,c)
            def holder=new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT,4,0));holder.setOpaque(false);parts.each{holder.add(it)}
            c.gridx=1;c.weightx=0d;c.fill=java.awt.GridBagConstraints.NONE;c.anchor=java.awt.GridBagConstraints.EAST;c.insets=new java.awt.Insets(3,0,3,0);p.add(holder,c)}
        def to={->CalciumTheme.label('to',CalciumTheme.MUTED)}
        def timing=new CalciumCard('','Timing',false);def timingForm=form()
        line(timingForm,'Frames per second',box('sampling_hz',timed?rate:'','From the recording; check it against the acquisition settings'))
        line(timingForm,'Baseline frames',box('baseline_start',data.baselineStart?:1,'First frame of the resting level (frame 1 is the first)'),to(),box('baseline_end',end,'Last frame of the resting level'))
        def stimulusName=box('stimulus','Stimulus','Name used in the tables');stimulusName.setColumns(9)
        line(timingForm,'Stimulus name',stimulusName)
        line(timingForm,'Stimulus window (s)',box('stimulus_start',timed?(end-1)/rate:'','Start of the window that is scored, in seconds'),to(),box('stimulus_end',timed?(n-1)/rate:'','End of the window that is scored, in seconds'))
        line(timingForm,'Baseline before it (s)',box('response_baseline_seconds',timed?(end-1)/rate:'','Seconds just before the stimulus used as the local reference'))
        timing.body.add(timingForm)
        def rule=new CalciumCard('','Response rule',false);def ruleForm=form()
        line(ruleForm,'Minimum rise (fraction)',box('amplitude_floor',0.10,'Smallest rise that counts, as a fraction of the local reference'))
        line(ruleForm,'Baseline SD multiplier',box('noise_multiplier',3,'The rise must also exceed this many standard deviations of the local reference'))
        line(ruleForm,'Consecutive samples',box('consecutive_samples',2,'Samples in a row that must be above the threshold'))
        def kclForm=form();line(kclForm,'KCl window (s)',box('kcl_start','','When the KCl was on: start, in seconds'),to(),box('kcl_end','','When the KCl was on: end, in seconds'))
        kclForm.setVisible(false)
        kcl.setToolTipText('<html>For recordings where KCl was added after the stimulus as a positive control.<br>Tick it and enter when the KCl was on; each cell is then also scored for that window.<br>Leave it unticked if no KCl was used.</html>')
        [kcl,confirmed,allTraces,outlinesBox].each{CalciumTheme.plain(it)}
        kcl.addActionListener({kclForm.setVisible(kcl.selected);root.revalidate();root.repaint()} as java.awt.event.ActionListener)
        rule.body.add(ruleForm);rule.body.add(kcl);rule.body.add(kclForm)
        def row={java.awt.Component... parts->def p=new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT,6,3));p.setOpaque(false);p.setAlignmentX(0f);parts.each{p.add(it)}
            p.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE,(int)p.getPreferredSize().height));p}
        status.setEditable(false);status.setLineWrap(true);status.setWrapStyleWord(true);status.setOpaque(false);status.setFocusable(false);status.setBorder(null)
        status.setFont(CalciumTheme.font());status.setForeground(CalciumTheme.INK);status.setAlignmentX(0f)
        status.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE,(int)status.getPreferredSize().height))
        progress.setAlignmentX(0f);progress.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE,14))
        setup.setToolTipText('One time: downloads a private copy of R for this plugin (needs internet)')
        def side=new javax.swing.JPanel();side.setLayout(new javax.swing.BoxLayout(side,javax.swing.BoxLayout.Y_AXIS));side.setOpaque(false)
        def heading=CalciumTheme.label('Analyze responses',CalciumTheme.INK,java.awt.Font.BOLD,5f)
        def recording=CalciumTheme.label(''+(data.recording?:''),CalciumTheme.MUTED);recording.setToolTipText(''+(data.recording?:''))
        recording.setMaximumSize(new java.awt.Dimension(280,24));recording.setPreferredSize(new java.awt.Dimension(280,(int)recording.getPreferredSize().height))
        def source=CalciumTheme.label("<html><body style='width:205px'>Uses the background-corrected fluorescence as measured: no lift, no second background subtraction.</body></html>",CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f)
        [heading,recording,javax.swing.Box.createVerticalStrut(8),timing,javax.swing.Box.createVerticalStrut(8),rule,javax.swing.Box.createVerticalStrut(8),confirmed,
         row(runButton,cancelButton),javax.swing.Box.createVerticalStrut(4),progress,javax.swing.Box.createVerticalStrut(4),status,row(setup,logButton),source,javax.swing.Box.createVerticalGlue()].each{side.add(it)}
        root.add(side,java.awt.BorderLayout.WEST)
        def stats=new javax.swing.JPanel(new java.awt.GridBagLayout());stats.setOpaque(false)
        [responding,registered,excluded,unscorable].eachWithIndex{tile,i->def c=new java.awt.GridBagConstraints();c.gridx=i;c.fill=java.awt.GridBagConstraints.BOTH;c.weightx=i==0?2.2d:1d;c.insets=new java.awt.Insets(0,0,0,i<3?8:0)
            tile.setPreferredSize(new java.awt.Dimension(i==0?220:100,(int)tile.getPreferredSize().height));stats.add(tile,c)}
        // the Traces tab: which cells, the picture at the cursor's moment, the traces, and the time slider under both
        [cells,group,view,heatmapOrder].each{it.setFont(CalciumTheme.font())}
        previousCell.setToolTipText('Previous cell (left arrow)');nextCell.setToolTipText('Next cell (right arrow)')
        group.setToolTipText('List and draw only the cells with this response call');allTraces.setToolTipText('Draw every listed cell behind the mean, coloured by its response call')
        def white={javax.swing.JPanel p->p.setOpaque(true);p.setBackground(Color.white);p}
        def pickers=white(row(CalciumTheme.label('Show'),group,CalciumTheme.label('Cell'),cells,previousCell,nextCell,allTraces))
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
        tabs.addTab('Traces',tracePanel);tabs.addTab('Heatmap',heatmapPanel)
        [metrics,calls,counts].each{it.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);it.setAutoCreateRowSorter(true);it.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION)
            it.setFont(CalciumTheme.font(java.awt.Font.PLAIN,-1f));it.setRowHeight(22);it.setGridColor(CalciumTheme.LINE);it.getTableHeader().setFont(CalciumTheme.font(java.awt.Font.BOLD,-1f))}
        tabs.addTab('Metrics + quality',new javax.swing.JScrollPane(metrics));tabs.addTab('Cell response calls',new javax.swing.JScrollPane(calls));tabs.addTab('Responder counts',new javax.swing.JScrollPane(counts))
        def bottom=row(folderButton,savePlots,settingsButton,CalciumTheme.label('Counts describe chosen cells; cell N is not animal N.',CalciumTheme.MUTED,java.awt.Font.PLAIN,-1f))
        def main=new javax.swing.JPanel(new java.awt.BorderLayout(0,10));main.setOpaque(false)
        main.add(stats,java.awt.BorderLayout.NORTH);main.add(tabs,java.awt.BorderLayout.CENTER);main.add(bottom,java.awt.BorderLayout.SOUTH)
        root.add(main,java.awt.BorderLayout.CENTER)
        setContentPane(new javax.swing.JScrollPane(root));root.setMinimumSize(new java.awt.Dimension(900,600));root.setPreferredSize(new java.awt.Dimension(1120,720))
        setSize(1320,860);setMinimumSize(new java.awt.Dimension(360,300))
        def bounds=getGraphicsConfiguration().getBounds();setSize(Math.min(getWidth(),(int)bounds.getWidth()-50),Math.min(getHeight(),(int)bounds.getHeight()-80))
        cancelButton.enabled=false;folderButton.enabled=false;savePlots.enabled=false;logButton.enabled=false;previousCell.enabled=false;nextCell.enabled=false
        zoomButton.enabled=false
        setup.addActionListener({beginSetup()} as java.awt.event.ActionListener)
        runButton.addActionListener({beginAnalysis()} as java.awt.event.ActionListener)
        cancelButton.addActionListener({job?.cancel();stage='Cancelling...';status.setText(stage)} as java.awt.event.ActionListener)
        logButton.addActionListener({if(job?.logFile?.isFile())IJ.open(job.logFile.path)} as java.awt.event.ActionListener)
        folderButton.addActionListener({if(results!=null)java.awt.Desktop.getDesktop().open(results)} as java.awt.event.ActionListener)
        savePlots.addActionListener({exportPlots()} as java.awt.event.ActionListener)
        settingsButton.addActionListener({side.setVisible(!side.isVisible());settingsButton.setText(side.isVisible()?'Hide settings':'Show settings');root.revalidate();root.repaint()} as java.awt.event.ActionListener)
        cells.addActionListener({if(!listing)selectCell(cells.selectedItem as String,true)} as java.awt.event.ActionListener)
        def step={int by->int count=cells.getItemCount();if(count>0)cells.setSelectedIndex((int)(((cells.getSelectedIndex()+by)%count+count)%count))}
        previousCell.addActionListener({step(-1)} as java.awt.event.ActionListener);nextCell.addActionListener({step(1)} as java.awt.event.ActionListener)
        group.addActionListener({if(!listing)listCells()} as java.awt.event.ActionListener)
        allTraces.addActionListener({plot.showAll=allTraces.selected;plot.changed()} as java.awt.event.ActionListener)
        heatmapOrder.addActionListener({orderHeatmap()} as java.awt.event.ActionListener)
        responding.pressed={->showGroup(1)};registered.pressed={->showGroup(0)};excluded.pressed={->showGroup(3)};unscorable.pressed={->showGroup(3)}
        [responding,registered,excluded,unscorable].each{it.setToolTipText('Click to list and draw only these cells')}
        heatmap.choose={int i->showGroup(0);cells.setSelectedItem(heatmap.names[i]);tabs.selectedIndex=0}
        cellPicture.choose={int i->if(!busy&&!closed&&cells.getItemCount()>0){String name=measurement.names[i];if((cells.getModel() as javax.swing.DefaultComboBoxModel).getIndexOf(name)<0)showGroup(0);cells.setSelectedItem(name)}}
        [metrics,calls].each{table->table.selectionModel.addListSelectionListener({e->
            if(!e.valueIsAdjusting&&table.selectedRow>=0) {
                int col=(0..<table.columnCount).find{table.getColumnName(it)=='Cell'} ?: 0
                String name=table.getValueAt(table.selectedRow,col).toString()
                if((cells.getModel() as javax.swing.DefaultComboBoxModel).getIndexOf(name)<0)showGroup(0)
                cells.setSelectedItem(name)
            }
        } as javax.swing.event.ListSelectionListener)}
        // time: the slider, the plot and the Play button all move the same time cursor
        plot.seek={double t->if(plot.time.length>1){int best=0;double gap=Double.MAX_VALUE;plot.time.eachWithIndex{double v,int i->double d=Math.abs(v-t);if(d<gap){gap=d;best=i}};showFrame(best+1)}}
        slider.addChangeListener({e->if(slider.getValue()!=frame)showFrame(slider.getValue())} as javax.swing.event.ChangeListener)
        player=new javax.swing.Timer(50,{e->if(frame>=slider.getMaximum())stopPlaying() else showFrame(frame+1)} as java.awt.event.ActionListener)
        playButton.addActionListener({if(player.isRunning())stopPlaying() else {if(frame>=slider.getMaximum())showFrame(1);player.start();playButton.setText('❚❚ Pause')}} as java.awt.event.ActionListener)
        view.addActionListener({showFrame(frame)} as java.awt.event.ActionListener)
        outlinesBox.addActionListener({cellPicture.outlines=outlinesBox.selected;cellPicture.repaint()} as java.awt.event.ActionListener)
        zoomButton.addActionListener({cellPicture.zoomTo(cellPicture.selected)} as java.awt.event.ActionListener)
        wholeButton.addActionListener({cellPicture.whole()} as java.awt.event.ActionListener)
        // left and right arrows step through the listed cells, unless a text box or list is using them
        def keys=getRootPane()
        keys.registerKeyboardAction({e->if(previousCell.enabled)step(-1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_LEFT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        keys.registerKeyboardAction({e->if(nextCell.enabled)step(1)} as java.awt.event.ActionListener,javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_RIGHT,0),javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
        clock=new javax.swing.Timer(500,{e->if(busy)status.setText(stage+' · '+(long)((System.currentTimeMillis()-started)/1000)+' s elapsed')} as java.awt.event.ActionListener)
        addWindowListener(new java.awt.event.WindowAdapter(){void windowClosed(java.awt.event.WindowEvent e){closed=true;job?.cancel();clock.stop();player.stop();(measurement.picture as ImagePlus)?.close()}})
        showFrame(1)
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
        double rise=CalciumFrameView.level(brightest,rest,0.9995d);if(!(rise>0d))rise=1d
        changeRange=[-0.1d*rise,rise] as double[]
        double low=CalciumFrameView.level(resting,null,0.005d),high=Math.max(CalciumFrameView.level(resting,null,0.9995d),CalciumFrameView.level(brightest,null,0.9995d));if(!(high>low))high=low+1d
        rawRange=[low,high] as double[]
        frameImage=new java.awt.image.BufferedImage(stack.getWidth(),stack.getHeight(),java.awt.image.BufferedImage.TYPE_INT_RGB);work=new float[stack.getWidth()*stack.getHeight()];work2=new float[work.length]
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
            double[] range=mode==1?changeRange:rawRange
            CalciumFrameView.paint(movie.getStack().getProcessor(frame),mode==1?rest:null,range[0],range[1],((java.awt.image.DataBufferInt)frameImage.getRaster().getDataBuffer()).getData(),work,work2,mode==1)
            cellPicture.image=frameImage;cellPicture.onGray=true
            cellPicture.note=mode==1?'frame '+frame+' minus the mean of frames '+restFirst+'–'+restLast+', lightly smoothed':'frame '+frame+' as recorded'
        } else if(snapshot!=null) {cellPicture.image=snapshot;cellPicture.onGray=false;cellPicture.note=movie!=null?'colour picture (does not change with time)':''}
        cellPicture.repaint()
    }
    /** Uses the free space to the left of the Step 1-3 window, so the two never cover each other. */
    void besideControls(java.awt.Window panel) {
        try {
            def configuration=panel.getGraphicsConfiguration();def bounds=configuration.getBounds();def insets=java.awt.Toolkit.getDefaultToolkit().getScreenInsets(configuration)
            int left=(int)bounds.getX()+insets.left+8,top=(int)bounds.getY()+insets.top+8,free=panel.getX()-left-8,tall=(int)bounds.getHeight()-insets.top-insets.bottom-16
            if(free>=760)setBounds(left,top,Math.min(free,1500),Math.min(tall,920))
        } catch(Throwable ignored) { }
    }
    Map settings() {
        // A box that still shows its rounded suggestion sends the exact value behind it; an edited box sends what was typed.
        def values=fields.collectEntries{k,v->String typed=v.text.trim();[k,exact[k]!=null&&exact[k][0]==typed?exact[k][1]:typed]}
        return CalciumAnalysisSettings.check(values,confirmed.selected,kcl.selected,measurement)
    }
    void launch(boolean installing,Closure action) {
        if(busy)return
        busy=true;started=System.currentTimeMillis();stage=installing?'Preparing R setup':'Preparing analysis'
        setup.enabled=false;runButton.enabled=false;cancelButton.enabled=true;progress.value=0;progress.indeterminate=installing
        fields.values().each{it.enabled=false};confirmed.enabled=false;kcl.enabled=false;savePlots.enabled=false
        status.setForeground(CalciumTheme.INK);status.setText(stage)
        job=new CalciumLocalJob(logFile:new File(runtime.home,installing?'setup.log':'analysis.log'))
        job.progress={int p,String message->javax.swing.SwingUtilities.invokeLater{if(!closed){progress.value=p;progress.setString(p+'% · '+message);stage=p+' % · '+message;logButton.enabled=job.logFile.isFile()}}}
        clock.start()
        Thread.startDaemon('SimpleCa analysis'){
            File output=null;String error=null
            try {output=action.call(job)}catch(java.util.concurrent.CancellationException ignored){error='Cancelled. Completed analyses are kept; partial output is not used.'}catch(Throwable t){error=t.message?:t.toString()}
            final File answer=output;final String failure=error
            javax.swing.SwingUtilities.invokeLater{
                busy=false;clock.stop();if(closed)return
                progress.indeterminate=false;cancelButton.enabled=false;setup.enabled=true;runButton.enabled=true
                fields.values().each{it.enabled=true};confirmed.enabled=true;kcl.enabled=true;logButton.enabled=job.logFile.isFile()
                savePlots.enabled=results!=null
                if(failure){status.setForeground(failure.startsWith('Cancelled')?CalciumTheme.INK:CalciumTheme.BAD);status.setText(failure);progress.setString('Not completed');return}
                if(installing){progress.value=100;status.setText('R is ready. Check the settings, then Analyze responses.');return}
                try {loadResults(answer);status.setForeground(CalciumTheme.GOOD);status.setText('Complete. Tables and settings saved in '+answer.name);progress.value=100;progress.setString('Complete')}
                catch(Throwable t){status.setForeground(CalciumTheme.BAD);status.setText('Could not display results: '+t.message);progress.setString('Display error; see saved outputs')}
            }
        }
    }
    void beginSetup(){launch(true){task->runtime.setup(task)}}
    void beginAnalysis(){
        try {Map s=settings();launch(false){task->performAnalysis.call(runtime,measurement.csv,s,measurement.folder,task)}}
        catch(Throwable t){status.setForeground(CalciumTheme.WARN);status.setText(t.message?:t.toString())}
    }
    void table(javax.swing.JTable target,File file) {
        def rows=CalciumCsv.read(file)
        def model=new javax.swing.table.DefaultTableModel(rows.drop(1).collect{it as Object[]} as Object[][],rows[0] as Object[]){boolean isCellEditable(int r,int c){false}}
        target.setModel(model);for(int i=0;i<target.columnCount;i++)target.columnModel.getColumn(i).preferredWidth=140
    }
    String callOf(String name){callRows[name]?.Status}
    /** The headline numbers and the per-cell calls, read from the tables the analysis just saved. Display only. */
    void summarize(File dir) {
        def countRows=CalciumCsv.read(new File(dir,'responder_counts.csv')),callTable=CalciumCsv.read(new File(dir,'responder_calls.csv'))
        def column={List rows,String name->(rows[0] as List).findIndexOf{it==name}}
        int role=column(countRows,'Role')
        def row=countRows.drop(1).find{role<0||it[role]=='Stimulus'}
        def cell={String name->int i=column(countRows,name);row!=null&&i>=0?''+row[i]:''}
        def whole={String name->double v=CalciumCsv.number(cell(name));Double.isFinite(v)?''+(long)Math.round(v):'–'}
        double percent=CalciumCsv.number(cell('Percent_responding'))
        responding.show(Double.isFinite(percent)?String.format(Locale.US,'%.0f %%',percent):'–','responding: '+whole('Responded')+' of '+whole('Evaluable')+' evaluable cells')
        registered.show(whole('Registered'),'cells registered');excluded.show(whole('Excluded'),'excluded');unscorable.show(whole('Unscorable'),'unscorable')
        int callRole=column(callTable,'Role'),callCell=column(callTable,'Cell');List header=callTable[0] as List;callRows=[:]
        if(callCell>=0)callTable.drop(1).each{line->if(callRole<0||line[callRole]=='Stimulus'){Map<String,String> one=[:];header.eachWithIndex{h,i->one[''+h]=i<line.size()?''+line[i]:''};callRows[''+line[callCell]]=one}}
        cellPicture.calls=callRows.collectEntries{k,v->[k,v.Status]} as Map<String,String>;cellPicture.repaint()
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
        else {plot.cell='';plot.selected=new double[0];cellStatus.setText('No cells in this group.');cellPicture.selected=-1;plot.changed();cellPicture.repaint()}
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
        heatmap.fitRows()
    }
    void loadResults(File dir) {
        results=dir;table(metrics,new File(dir,'metrics.csv'));table(calls,new File(dir,'responder_calls.csv'));table(counts,new File(dir,'responder_counts.csv'))
        def normalized=CalciumCsv.read(new File(dir,'normalized_traces.csv')),avg=CalciumCsv.read(new File(dir,'mean_trace.csv'))
        plot.time=normalized.drop(1).collect{CalciumCsv.number(it[0])} as double[]
        traces=[:];normalized[0].drop(1).eachWithIndex{name,i->traces[name]=normalized.drop(1).collect{CalciumCsv.number(it[i+1])} as double[]}
        plot.mean=avg.drop(1).collect{CalciumCsv.number(it[1])} as double[];plot.sem=avg.drop(1).collect{CalciumCsv.number(it[2])} as double[]
        Map manifest=CalciumCsv.read(new File(dir,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
        double hz=CalciumCsv.number(manifest.sampling_hz)
        plot.baselineStart=(CalciumCsv.number(manifest.baseline_start)-1)/hz;plot.baselineEnd=(CalciumCsv.number(manifest.baseline_end)-1)/hz
        plot.stimulusStart=CalciumCsv.number(manifest.stimulus_start);plot.stimulusEnd=CalciumCsv.number(manifest.stimulus_end)
        callRows=[:]
        try {summarize(dir)} catch(Throwable ignored) { } // The tables below stay the record; a missing headline must not hide them.
        orderHeatmap()
        // the cell groups, with how many cells each holds
        listing=true
        try {int keep=Math.max(0,group.getSelectedIndex());group.removeAllItems()
            CalciumCalls.GROUPS.eachWithIndex{String name,int i->group.addItem(name+' ('+(measurement.names as List<String>).count{callRows.isEmpty()||CalciumCalls.inGroup(callOf(it),i)}+')')}
            group.setSelectedIndex(callRows.isEmpty()?0:Math.min(keep,group.getItemCount()-1))} finally {listing=false}
        // the movie's resting picture follows the baseline frames of this analysis
        double first=CalciumCsv.number(manifest.baseline_start),last=CalciumCsv.number(manifest.baseline_end)
        if(Double.isFinite(first)&&Double.isFinite(last)){restFirst=(int)first;restLast=(int)last};rest=null
        if(movie!=null&&view.getSelectedIndex()==0&&view.getItemCount()>1)view.setSelectedIndex(1)   // show the coupled movie as soon as there are traces to follow
        listCells()
        folderButton.enabled=true;savePlots.enabled=true
        // start at the moment the cells are brightest on average
        int peak=1;double top=Double.NEGATIVE_INFINITY;plot.mean.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top){top=v;peak=i+1}}
        showFrame(peak)
    }
    void selectCell(String name,boolean highlight) {
        if(name==null)return
        plot.cell=name;plot.selected=traces[name]?:new double[0]
        plot.caption=name+(traces.containsKey(name)?' · gray: baseline; pale orange: stimulus window':' · unavailable: inspect Metrics + quality for the exclusion reason')
        plot.changed()
        int index=measurement.names.findIndexOf{it==name}
        if(highlight&&index>=0&&measurement.picture!=null){Roi roi=(Roi)measurement.rois[index].clone();roi.setStrokeColor(Color.white);roi.setStrokeWidth(4);(measurement.picture as ImagePlus).setRoi(roi)}
        cellPicture.selected=highlight?index:-1;cellPicture.reveal(cellPicture.selected);cellPicture.repaint()
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
    void exportPlots() {
        if(results==null)return
        try {
            ['trace':plot,'heatmap':heatmap].each{name,panel->
                int w=Math.max(650,panel.width),h=Math.max(260,panel.height)
                def image=new java.awt.image.BufferedImage(w,h,java.awt.image.BufferedImage.TYPE_INT_RGB)
                def g=image.createGraphics();try{panel.paint(g)}finally{g.dispose()}
                javax.imageio.ImageIO.write(image,'png',new File(results,name+'.png'))
            }
            status.setText('Saved trace.png and heatmap.png in '+results.name)
        }catch(Throwable t){status.setText('Could not save plots: '+t.message)}
    }
}

/** The Step 4 settings: what is suggested, how a suggestion is shown, and the checks before an analysis may start.
 *  Shared by the Fiji window and the web page, so both send exactly the same values to the analysis. */
class CalciumAnalysisSettings {
    static Map<String,Object> suggested(Map data) {
        int n=data.frames as int, end=Math.min(n-1,Math.max(3,(data.baselineEnd?:20) as int))
        double rate=(data.hz?:Double.NaN) as double; boolean timed=Double.isFinite(rate)&&rate>0
        return [sampling_hz:timed?rate:'', baseline_start:data.baselineStart?:1, baseline_end:end, stimulus:'Stimulus',
            stimulus_start:timed?(end-1)/rate:'', stimulus_end:timed?(n-1)/rate:'', response_baseline_seconds:timed?(end-1)/rate:'',
            amplitude_floor:0.10, noise_multiplier:3, consecutive_samples:2, kcl_start:'', kcl_end:'']
    }
    static String full(Object value){value==null?'':''+value}
    /** A suggested number is shown with at most three decimals; the exact value is what the analysis receives. */
    static String shown(Object value){(value instanceof Double||value instanceof BigDecimal||value instanceof Float)?CalciumTheme.shortNumber(value as double):full(value)}
    /** The checks before an analysis may start. Returns the values to send, or throws with a message for the user. */
    static Map check(Map typed,boolean confirmed,boolean kcl,Map measurement) {
        Map values=new LinkedHashMap(typed)
        if(!confirmed)throw new IllegalArgumentException('Confirm the frame rate, baseline and stimulus times first.')
        def number={key->double v=CalciumCsv.number(values[key]);if(!Double.isFinite(v))throw new IllegalArgumentException('Enter a number for '+key);v}
        double hz=number('sampling_hz'),first=number('baseline_start'),last=number('baseline_end'),start=number('stimulus_start'),end=number('stimulus_end'),pre=number('response_baseline_seconds')
        if(hz<=0||first<1||first!=Math.floor(first)||last!=Math.floor(last)||last<first||last>=measurement.frames)throw new IllegalArgumentException('Check the positive frame rate and baseline frame range.')
        if(start<0||end<=start||pre<=0||start-pre< -1e-8||end>(measurement.frames-1)/hz+1.5/hz)throw new IllegalArgumentException('The stimulus and prestimulus baseline must fit this recording.')
        if(!values.stimulus||number('amplitude_floor')<=0||number('noise_multiplier')<0||number('consecutive_samples')<1||number('consecutive_samples')!=Math.floor(number('consecutive_samples')))throw new IllegalArgumentException('Check the stimulus name and response thresholds.')
        if(kcl){number('kcl_start');number('kcl_end')}
        values.kcl_enabled=''+kcl;values.recording=measurement.recording;values.input_mode='background_corrected_fluorescence'
        if(measurement.inputHash)values.expected_input_sha256=measurement.inputHash
        return values
    }
}

/** One measurement's Step 4 analysis without any window: what the Step 4 web page shows and drives.
 *  It runs the same engine with the same settings as the Fiji window, and reads the tables that engine saves. */
class CalciumAnalysisSession {
    Map measurement
    CalciumRRuntime runtime=new CalciumRRuntime()
    Closure performAnalysis={rt,input,settings,out,job->CalciumSimpleCaEngine.analyze(rt,input,settings,out,job)}
    Closure open={File target->java.awt.Desktop.getDesktop().open(target)}
    CalciumLocalJob job
    volatile boolean busy=false, installing=false
    volatile int percent=0, version=0
    volatile String stage='', message='Check the settings, tick the box, then analyze.', tone='plain'
    volatile long started=0L
    volatile File results
    volatile String resultsJson='{}'
    Map<String,Object> suggested
    ImagePlus movie; byte[] picture; int width=0, height=0
    float[] rest, work, work2; double[] rawRange, changeRange, meanTrace=new double[0]
    int restFirst=1, restLast=1
    java.awt.image.BufferedImage frameImage
    CalciumAnalysisSession(Map data) {
        measurement=new LinkedHashMap(data);suggested=CalciumAnalysisSettings.suggested(data)
        restFirst=(data.baselineStart?:1) as int;restLast=suggested.baseline_end as int
        if(data.picture!=null) {
            def image=(data.picture as ImagePlus).getBufferedImage();width=image.getWidth();height=image.getHeight()
            def out=new ByteArrayOutputStream();javax.imageio.ImageIO.write(image,'png',out);picture=out.toByteArray()
        }
        // The recording is only read, to show the frame the page's time cursor is on. It must match the picture the outlines were drawn on.
        if(data.movie instanceof ImagePlus&&picture!=null&&(data.movie as ImagePlus).getWidth()==width&&(data.movie as ImagePlus).getHeight()==height&&(data.movie as ImagePlus).getStackSize()>=(data.frames as int))movie=data.movie as ImagePlus
    }
    /** JSON text; numbers that are not finite become null, which is what a web page expects for "no value". */
    static String json(Object value) {
        def clean;clean={v->v instanceof Map?v.collectEntries{k,x->[''+k,clean(x)]}:v instanceof Collection?v.collect{clean(it)}:(v instanceof double[]||v instanceof float[]||v instanceof int[])?(v as List).collect{clean(it)}:(v instanceof Double||v instanceof Float)?(Double.isFinite(v as double)?v:null):v instanceof File?v.name:v}
        return groovy.json.JsonOutput.toJson(clean(value))
    }
    /** What never changes for this measurement: the cells, their outlines and the suggested settings. */
    Map described() {
        [recording:''+(measurement.recording?:''),key:''+(measurement.key?:''),frames:measurement.frames as int,names:(measurement.names?:[]) as List,width:width,height:height,
         hasPicture:picture!=null,hasMovie:movie!=null,pluginVersion:CalciumRelease.VERSION,
         rois:((measurement.rois?:new Roi[0]) as List).collect{Roi roi->def p=roi.getPolygon();(0..<p.npoints).collect{int i->[p.xpoints[i],p.ypoints[i]]}},
         settings:suggested.collect{k,v->[key:k,shown:CalciumAnalysisSettings.shown(v),exact:CalciumAnalysisSettings.full(v)]}]
    }
    /** What can change at any moment: whether a job is running, how far it is, and which results are loaded. */
    Map state() {
        boolean ready=false;try{ready=runtime.environment()!=null}catch(Throwable ignored){}
        [key:''+(measurement.key?:''),busy:busy,installing:installing,percent:percent,stage:stage,message:message,tone:tone,elapsed:busy?(long)((System.currentTimeMillis()-started)/1000L):0L,
         version:version,folder:results?.name,engineReady:ready,logAvailable:job?.logFile?.isFile()?:false,rest:[restFirst,restLast]]
    }
    void launch(boolean setup,Closure action) {
        synchronized(this){if(busy)throw new IllegalStateException('A job is already running.');busy=true}
        installing=setup;started=System.currentTimeMillis();percent=0;stage=setup?'Preparing R setup':'Preparing analysis';message=stage;tone='plain'
        final CalciumLocalJob mine=new CalciumLocalJob(logFile:new File(runtime.home,setup?'setup.log':'analysis.log'));job=mine
        mine.progress={int p,String text->percent=p;stage=text}
        Thread.startDaemon('SimpleCa analysis for the Step 4 page') {
            try {
                File output=action.call(mine) as File
                if(setup){percent=100;message='R is ready. Check the settings, then Analyze responses.';tone='good'}
                else {load(output);percent=100;message='Complete. Tables and settings saved in '+output.name;tone='good'}
            } catch(java.util.concurrent.CancellationException ignored){message='Cancelled. Completed analyses are kept; partial output is not used.';tone='plain'}
            catch(Throwable t){message=t.message?:t.toString();tone='bad'}
            finally {busy=false}
        }
    }
    void setup(){launch(true){task->runtime.setup(task)}}
    /** "typed" is what the page's boxes hold. A box that still shows its rounded suggestion sends the exact value behind it. */
    void analyze(Map typed,boolean confirmed,boolean kcl) {
        Map values=[:]
        suggested.each{k,v->String text=(typed[k]==null?'':''+typed[k]).trim();values[k]=text==CalciumAnalysisSettings.shown(v)?CalciumAnalysisSettings.full(v):text}
        Map checked=CalciumAnalysisSettings.check(values,confirmed,kcl,measurement)
        launch(false){task->performAnalysis.call(runtime,measurement.csv,checked,measurement.folder,task)}
    }
    void cancel(){job?.cancel();if(busy)stage='Cancelling...'}
    /** Reads the tables the analysis saved into the form the page draws from. Display only: nothing is recalculated. */
    void load(File dir) {
        def normalized=CalciumCsv.read(new File(dir,'normalized_traces.csv')),average=CalciumCsv.read(new File(dir,'mean_trace.csv'))
        def numbers={List rows,int column->rows.drop(1).collect{CalciumCsv.number(it[column])}}
        Map traces=[:];normalized[0].drop(1).eachWithIndex{name,i->traces[''+name]=numbers(normalized,i+1)}
        Map manifest=CalciumCsv.read(new File(dir,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
        double hz=CalciumCsv.number(manifest.sampling_hz),first=CalciumCsv.number(manifest.baseline_start),last=CalciumCsv.number(manifest.baseline_end)
        def table={String name->def rows=CalciumCsv.read(new File(dir,name));[header:rows[0],rows:rows.drop(1)]}
        def callTable=table('responder_calls.csv'),countTable=table('responder_counts.csv')
        def column={Map t,String name->(t.header as List).findIndexOf{it==name}}
        Map calls=[:];int role=column(callTable,'Role'),cell=column(callTable,'Cell')
        if(cell>=0)callTable.rows.each{line->if(role<0||line[role]=='Stimulus'){Map one=[:];callTable.header.eachWithIndex{h,i->one[''+h]=i<line.size()?''+line[i]:''};calls[''+line[cell]]=one}}
        int countRole=column(countTable,'Role');def countRow=countTable.rows.find{countRole<0||it[countRole]=='Stimulus'};Map counts=[:]
        if(countRow!=null)countTable.header.eachWithIndex{h,i->counts[''+h]=i<countRow.size()?''+countRow[i]:''}
        List mean=numbers(average,1)
        Map loaded=[folder:dir.name,time:numbers(normalized,0),mean:mean,sem:numbers(average,2),traces:traces,calls:calls,counts:counts,
            baseline:[(first-1)/hz,(last-1)/hz],baselineFrames:[first,last],stimulus:[CalciumCsv.number(manifest.stimulus_start),CalciumCsv.number(manifest.stimulus_end)],
            tables:[metrics:table('metrics.csv'),calls:callTable,counts:countTable]]
        String text=json(loaded)
        synchronized(this) {
            meanTrace=mean as double[]
            // the movie's resting picture follows the baseline frames of this analysis
            if(Double.isFinite(first)&&Double.isFinite(last)){restFirst=(int)first;restLast=(int)last};rest=null
            results=dir;resultsJson=text;version++
        }
    }
    /** One frame of the recording as a JPEG: as recorded, or minus the mean of the baseline frames and lightly smoothed. Display only. */
    synchronized byte[] frame(int wanted,String mode) {
        if(movie==null)throw new FileNotFoundException('No recording is attached to this measurement.')
        def stack=movie.getStack();int n=measurement.frames as int,f=Math.max(1,Math.min(n,wanted))
        if(rest==null) {
            rest=CalciumFrameView.mean(stack,restFirst,restLast)
            // the frame where the cells are brightest on average sets the display range, so it does not flicker from frame to frame
            int peak=Math.max(1,Math.min(n,(restLast+n).intdiv(2)));double top=Double.NEGATIVE_INFINITY
            meanTrace.eachWithIndex{double v,int i->if(Double.isFinite(v)&&v>top&&i<n){top=v;peak=i+1}}
            def brightest=stack.getProcessor(peak);def resting=new ij.process.FloatProcessor(width,height,rest)
            double rise=CalciumFrameView.level(brightest,rest,0.9995d);if(!(rise>0d))rise=1d
            changeRange=[-0.1d*rise,rise] as double[]
            double low=CalciumFrameView.level(resting,null,0.005d),high=Math.max(CalciumFrameView.level(resting,null,0.9995d),CalciumFrameView.level(brightest,null,0.9995d));if(!(high>low))high=low+1d
            rawRange=[low,high] as double[]
            frameImage=new java.awt.image.BufferedImage(width,height,java.awt.image.BufferedImage.TYPE_INT_RGB);work=new float[width*height];work2=new float[width*height]
        }
        boolean change=mode!='raw';double[] range=change?changeRange:rawRange
        CalciumFrameView.paint(stack.getProcessor(f),change?rest:null,range[0],range[1],((java.awt.image.DataBufferInt)frameImage.getRaster().getDataBuffer()).getData(),work,work2,change)
        def out=new ByteArrayOutputStream(1<<18);javax.imageio.ImageIO.write(frameImage,'jpg',out);return out.toByteArray()
    }
}

/** Serves the Step 4 page to a browser on this computer only. Nothing is reachable from other computers, and every
 *  address carries a random key, so no other program or web page can read or drive the analysis. */
class CalciumWebServer {
    volatile CalciumAnalysisSession session
    def server; int port=0
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
    String url(){'http://127.0.0.1:'+port+'/'+token+'/'}
    /** Throws when this Java cannot run the small built-in web server; the caller then opens the Fiji window instead. */
    void start() {
        files=page()
        // Looked up by name, so the plugin still loads on a Java that leaves this optional part out.
        def serverClass=Class.forName('com.sun.net.httpserver.HttpServer'),handlerClass=Class.forName('com.sun.net.httpserver.HttpHandler')
        server=serverClass.create(new InetSocketAddress(InetAddress.getLoopbackAddress(),0),0)
        server.createContext('/',{exchange->handle(exchange)}.asType(handlerClass))
        server.setExecutor(java.util.concurrent.Executors.newCachedThreadPool({Runnable r->def t=new Thread(r,'Calcium Step 4 page');t.setDaemon(true);t} as java.util.concurrent.ThreadFactory))
        server.start();port=server.getAddress().getPort()
    }
    void stop(){try{server?.stop(0)}catch(Throwable ignored){}}
    static final Map TYPES=[html:'text/html; charset=utf-8',css:'text/css; charset=utf-8',js:'text/javascript; charset=utf-8',json:'application/json; charset=utf-8',png:'image/png',jpg:'image/jpeg',txt:'text/plain; charset=utf-8']
    void reply(exchange,int status,String type,byte[] body) {
        exchange.getResponseHeaders().set('Content-Type',TYPES[type]?:type);exchange.getResponseHeaders().set('Cache-Control','no-store');exchange.getResponseHeaders().set('X-Content-Type-Options','nosniff')
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
                    case 'api/frame':reply(exchange,200,'jpg',now.frame((query.f?:'1') as int,''+(query.mode?:'change')));return
                }
            } else if(method=='POST'&&now!=null) {
                byte[] sent=read(exchange.getRequestBody(),1<<16)
                Map body=sent.length>0?(new groovy.json.JsonSlurper().parseText(new String(sent,'UTF-8')) as Map):[:]
                switch(what) {
                    case 'api/analyze':now.analyze((body.values?:[:]) as Map,body.confirmed==true,body.kcl==true);replyJson(exchange,200,now.state());return
                    case 'api/setup':now.setup();replyJson(exchange,200,now.state());return
                    case 'api/cancel':now.cancel();replyJson(exchange,200,now.state());return
                    case 'api/open':
                        File target=body.what=='log'?now.job?.logFile:now.results
                        if(target==null||!target.exists()){replyJson(exchange,404,[error:body.what=='log'?'There is no log yet.':'There are no saved results yet.']);return}
                        now.open.call(target);replyJson(exchange,200,[opened:target.name]);return
                }
            }
            reply(exchange,404,'txt','Not found'.bytes)
        } catch(IllegalArgumentException|IllegalStateException refused){try{replyJson(exchange,400,[error:refused.message?:'That is not possible right now.'])}catch(Throwable ignored){}}
        catch(Throwable problem){try{replyJson(exchange,500,[error:problem.message?:problem.toString()])}catch(Throwable ignored){}}
        finally {exchange.close()}
    }
}

// BEGIN GENERATED SIMPLECA BUNDLE
class CalciumSimpleCaBundle {
    static final String HASH = "dcb8a18826fbfdf630e7b0637e577b7418544633b325ff15ff9b46bec833743a"
    static final String DATA = "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"
}
// END GENERATED SIMPLECA BUNDLE

// BEGIN GENERATED WEB BUNDLE
class CalciumWebBundle {
    static final String SHA256 = "9ce57580b923fed98610992b405199d10563897b40e0197020d3f8489b0a9a49"
    static final String CONTENT = "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"
}
// END GENERATED WEB BUNDLE

/** Cell edits keep ROI order/names, and can be undone without overwriting outside edits. */
class CalciumCellEdits {
    RoiManager manager
    ImagePlus picture
    List history = []

    List snapshot() {
        (0..<manager.getCount()).collect { int i -> [roi: (Roi) manager.getRoi(i).clone(), name: manager.getName(i)] }
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

/** Step 1: which recording, whether to look for a response, and where the results go. */
class CalciumStartDialog {
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
            def look = new javax.swing.JCheckBox("Look for a response and show it in magenta", true); CalciumTheme.plain(look)
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
            def showFileRow = { -> fileRow.setVisible(choice.getSelectedItem() == diskLabel); dialog.pack() }
            fileButton.addActionListener({ e -> pickFile() } as java.awt.event.ActionListener)
            choice.addActionListener({ e -> showFileRow(); if (dialog.isVisible() && choice.getSelectedItem() == diskLabel) javax.swing.SwingUtilities.invokeLater { pickFile() } } as java.awt.event.ActionListener)
            fileRow.setVisible(choice.getSelectedItem() == diskLabel)
            [choice, fileRow, javax.swing.Box.createVerticalStrut(4), hint("Pick the original movie, not a processed copy."), javax.swing.Box.createVerticalStrut(8), look].each { recording.body.add(it) }
            def results = new CalciumCard("", "Folder for the results", false)
            def folder = new javax.swing.JTextField(defaultFolder, 30); folder.setFont(CalciumTheme.font()); def browse = CalciumTheme.button("Browse…")
            def line = new javax.swing.JPanel(new java.awt.BorderLayout(6, 0)); line.setOpaque(false); line.setAlignmentX(0f); line.add(folder, java.awt.BorderLayout.CENTER); line.add(browse, java.awt.BorderLayout.EAST)
            line.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) line.getPreferredSize().height))
            [line, javax.swing.Box.createVerticalStrut(4), hint("A subfolder with the recording's name is made inside this folder.")].each { results.body.add(it) }
            def next = CalciumTheme.button("Next", true), cancel = CalciumTheme.button("Cancel")
            def buttons = new javax.swing.JPanel(new java.awt.BorderLayout()); buttons.setOpaque(false); buttons.setAlignmentX(0f)
            def right = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT, 6, 0)); right.setOpaque(false); right.add(cancel); right.add(next)
            buttons.add(CalciumTheme.label("Version " + CalciumRelease.VERSION, CalciumTheme.FAINT, java.awt.Font.PLAIN, -1f), java.awt.BorderLayout.WEST); buttons.add(right, java.awt.BorderLayout.EAST)
            [recording, javax.swing.Box.createVerticalStrut(10), results, javax.swing.Box.createVerticalStrut(14), buttons].each { content.add(it) }
            browse.addActionListener({ e ->
                String chosen = chooseFolder.call(folder.getText())
                if (chosen) folder.setText(chosen)
            } as java.awt.event.ActionListener)
            next.addActionListener({ e ->
                boolean disk = choice.getSelectedItem() == diskLabel
                if (disk && file == null) pickFile()          // Next without a file still asks for one,
                if (disk && file == null) return              // and stays here if that was cancelled
                answer = [recording: choice.getSelectedItem() as String, look: look.isSelected(), folder: folder.getText(), file: disk ? file : null]; dialog.dispose()
            } as java.awt.event.ActionListener)
            cancel.addActionListener({ e -> dialog.dispose() } as java.awt.event.ActionListener)
            dialog.getRootPane().setDefaultButton(next)
            dialog.getRootPane().registerKeyboardAction({ e -> dialog.dispose() } as java.awt.event.ActionListener, javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_ESCAPE, 0), javax.swing.JComponent.WHEN_IN_FOCUSED_WINDOW)
            dialog.setContentPane(content); dialog.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
            dialog.pack(); dialog.setResizable(false); dialog.setLocationRelativeTo(null); dialog.setVisible(true)
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
steps.prepare(imp, look)
steps.colourPic.show()
IJ.run(steps.colourPic, "Set... ", "zoom=75")
def rm = testGui ? new RoiManager(false) : RoiManager.getRoiManager()
if (!testGui) { rm.runCommand(steps.colourPic, "Show None"); IJ.setTool("oval") } // The plugin draws its own contrasted overlay.
def base = imp.getTitle(); int dot = base.lastIndexOf("."); if (dot > 0) base = base.substring(0, dot)
def outDir = new File(outRoot, base)

// ---- the window that stays open for steps 2 and 3 ----
def GRADE_TEXT = ["responds on its own", "likely a cell", "unsure", "very bright: probably precipitate"]
def GRADE_SHARE = ["about 8 in 10 such outlines were real cells in tests", "about 8 in 10 such outlines were real cells in tests", "about 4 in 10 such outlines were real cells in tests", "about 1 in 10 such outlines were real cells in tests"]
def CYAN = new Color(0, 220, 255), GREY = new Color(205, 205, 205), ORANGE = new Color(255, 150, 0), BLUE = new Color(40, 150, 255)
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
def legendRow = { javax.swing.Icon ic, String s -> def l = new javax.swing.JLabel(s, ic, javax.swing.SwingConstants.LEFT); l.setFont(CalciumTheme.font(java.awt.Font.PLAIN, -1f)); l.setForeground(CalciumTheme.MUTED); l.setIconTextGap(6); l.setAlignmentX(0f); l }
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
def recordingNote = note("<b>" + CalciumTheme.middle(imp.getTitle(), 40) + "</b><br><font color='" + MUTED + "'>" + steps.n + " frames · " + (steps.found ? "response near frame " + steps.onset + " (magenta)" : (look ? "no response found (no magenta)" : "response not looked for")) + "</font>")
recordingNote.setToolTipText(imp.getTitle() + ". " + (steps.found ? "Brightness rises near frame " + steps.onset + "; the response is shown in magenta." : (look ? "No response was found, so the picture has no magenta. Every cell is still shown in green." : "The response was not looked for.")))
s1.body.add(recordingNote)

def s2 = section("2", "Choose the cells")
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
def cpLabel = hint("&nbsp;")
def setupBtn = btn("Set up Cellpose (one time)"); def setupRow = row(setupBtn)
def detectionProgress = new CalciumOutlineProgress()
def continueBtn = btn("Continue to measure  ›")
s2.body.add(row(suggestBtn, acceptAllBtn)); s2.body.add(detectionProgress); s2.body.add(cpLabel); s2.body.add(setupRow)
s2.body.add(row(CalciumTheme.label("Or draw:", CalciumTheme.MUTED), ovalBtn, freeBtn, addBtn))
s2.body.add(row(showSug, showBg)); s2.body.add(row(edgeBox))
s2.body.add(javax.swing.Box.createVerticalStrut(4)); s2.body.add(status)
s2.body.add(row(removeBtn, undoBtn, removeMode)); s2.body.add(editHint)
s2.body.add(info); s2.body.add(trace); s2.body.add(javax.swing.Box.createVerticalStrut(6))
def legend = new javax.swing.JPanel(new java.awt.GridLayout(0, 3, 6, 2)); legend.setAlignmentX(0f); legend.setOpaque(false)
def acceptedIcon = new java.awt.image.BufferedImage(30, 14, java.awt.image.BufferedImage.TYPE_INT_ARGB)
def acceptedGraphics = acceptedIcon.createGraphics()
acceptedGraphics.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
acceptedGraphics.setColor(new Color(24, 30, 36)); acceptedGraphics.fillRoundRect(0, 0, 30, 14, 6, 6)
acceptedGraphics.setColor(Color.black); acceptedGraphics.setStroke(new java.awt.BasicStroke(5f)); acceptedGraphics.drawOval(6, 3, 18, 8)
acceptedGraphics.setColor(Color.white); acceptedGraphics.setStroke(new java.awt.BasicStroke(2f)); acceptedGraphics.drawOval(6, 3, 18, 8); acceptedGraphics.dispose()
legend.add(legendRow(new javax.swing.ImageIcon(acceptedIcon), "your cells"))
legend.add(legendRow(lineIcon(CYAN, 2, false, false), "responds"))
legend.add(legendRow(lineIcon(GREY, 1, true, false), "unsure"))
legend.add(legendRow(lineIcon(BLUE, 1, false, true), "background"))
legend.add(legendRow(lineIcon(CYAN, 1, false, false), "likely a cell"))
legend.add(legendRow(lineIcon(ORANGE, 1, false, false), "precipitate?"))
legend.setToolTipText("White: your cells. Blue: the background used. The other four are suggestions, by grade.")
legend.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, (int) legend.getPreferredSize().height))
s2.body.add(legend); s2.body.add(row(continueBtn))

def s3 = section("3", "Measure and save")
def makeMovie = box("Also open a background-corrected movie", true)
def saveMovie = box("Save that movie as a TIFF too (large)", false)
saveMovie.setToolTipText("About 1.7 GB for 300 frames of 1388 x 1040")
def go = CalciumTheme.button("Measure and save", true)
def savedIn = "..." + File.separator + CalciumTheme.middle(outDir.getParentFile().getName(), 24) + File.separator + "<br>" + CalciumTheme.middle(outDir.getName(), 44)
def result = note("Saves the tables in: " + savedIn); result.setToolTipText(outDir.getPath())
def openBtn = btn("Open results folder"); openBtn.setVisible(false)
def closeBtn = btn("Close")
closeBtn.setToolTipText("Close this window; the pictures and saved results stay")
s3.body.add(makeMovie); s3.body.add(saveMovie); s3.body.add(row(go, openBtn)); s3.body.add(result)
def s4 = section("4", "Analyze responses")
def analyzeBtn = CalciumTheme.button("Analyze responses", true); analyzeBtn.setEnabled(false)
analyzeBtn.setToolTipText("Opens the analysis page in your browser. It runs on this computer only; nothing is sent anywhere.")
def classicBtn = btn("Use a Fiji window instead"); classicBtn.setEnabled(false)
classicBtn.setToolTipText("The same analysis in a Fiji window, for when a browser is not wanted")
analyzeBtn.addPropertyChangeListener("enabled", { classicBtn.setEnabled(analyzeBtn.isEnabled()) } as java.beans.PropertyChangeListener)
def analysisNote = hint("Measure and save first. Then review traces, metrics and responder counts here. One-time engine setup needs internet.")
Map analysisMeasurement = null
s4.body.add(row(analyzeBtn)); s4.body.add(analysisNote); s4.body.add(row(classicBtn))
def stepCards = [s2, s3, s4]
def openStep = { CalciumCard card -> stepCards.each { it.setExpanded(it.is(card)) }; frame.refreshContents() }
stepCards.each { it.onToggle = { CalciumCard card -> if (card.expanded) { card.setExpanded(false); frame.refreshContents() } else openStep(card) } }
s3.setExpanded(false); s4.setExpanded(false); s3.setSummary("after choosing cells"); s4.setSummary("SimpleCa²⁺, after measuring")
continueBtn.addActionListener({ e -> openStep(s3) } as java.awt.event.ActionListener)
// Version and update controls stay out of the way until asked for (or until an update is found).
def updatesCard = new CalciumCard("", "Updates", false); updatesCard.setVisible(false); updatesCard.setSummary("Version " + CalciumRelease.VERSION)
def updatesBtn = btn("Updates"); updatesBtn.setToolTipText("Check for a newer version of this plugin")
updatesBtn.addActionListener({ e -> updatesCard.setVisible(!updatesCard.isVisible()); frame.refreshContents() } as java.awt.event.ActionListener)
[s1, s2, s3, s4].each { root.add(it); root.add(javax.swing.Box.createVerticalStrut(8)) }
root.add(row(closeBtn, updatesBtn, CalciumTheme.label("Version " + CalciumRelease.VERSION, CalciumTheme.FAINT, java.awt.Font.PLAIN, -1f)))
root.add(javax.swing.Box.createVerticalStrut(4)); root.add(updatesCard)
def updateBtn = btn("Check for updates")
def installBtn = btn("Install update"); installBtn.setVisible(false)
def releaseBtn = btn("Release notes"); releaseBtn.setVisible(false)
def installBar = new javax.swing.JProgressBar(0,100); installBar.setStringPainted(true); installBar.setVisible(false); installBar.setAlignmentX(0f)
installBar.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE,22))
def autoUpdateBox = box("Check for updates when this plugin opens", Prefs.get("calcium.steps.check_updates", true))
def updateCheck = new CalciumUpdateCheck(enabled: autoUpdateBox.isSelected())
def updateNote = note("New releases offer Install update or Later. Restart Fiji after installation.")
def installBusy = new java.util.concurrent.atomic.AtomicBoolean(false)
Map selectedRelease = null
String installedVersion = null
def installRelease = { Map found, Closure progress ->
    String directory = IJ.getDirectory("imagej")
    new CalciumInstaller(applicationDirectory: directory ? new File(directory) : null).install(found, progress)
}
updatesCard.body.add(row(updateBtn, installBtn, releaseBtn)); updatesCard.body.add(autoUpdateBox); updatesCard.body.add(updateNote); updatesCard.body.add(installBar)
def releasePage = CalciumRelease.REPOSITORY + "/releases/latest"
def openRelease = { ->
    try { java.awt.Desktop.getDesktop().browse(new java.net.URI(releasePage)) }
    catch (Exception problem) { IJ.showMessage("Plugin updates", "Open this page in your browser:\n" + releasePage) }
}
def startInstall = { Map found ->
    if (installedVersion != null || found == null || !found.installable || !installBusy.compareAndSet(false,true)) return
    installBtn.setEnabled(false); updateBtn.setEnabled(false)
    installBar.setValue(0); installBar.setString("Preparing installation..."); installBar.setVisible(true)
    updateNote.setText("Downloading and verifying the update. Your current session stays open.")
    frame.refreshContents()
    Thread.start("Calcium plugin installation") {
        Map installed = null; String problem = null
        try {
            installed = installRelease.call(found, { int percent, String message ->
                javax.swing.SwingUtilities.invokeLater {
                    if (frame.isDisplayable()) { installBar.setValue(percent); installBar.setString(message) }
                }
            })
        } catch (Exception failure) { problem = failure.message ?: failure.toString() }
        final Map completed = installed; final String error = problem
        javax.swing.SwingUtilities.invokeLater {
            installBusy.set(false)
            if (completed != null) { installedVersion = completed.version; updateCheck.close() }
            if (!frame.isDisplayable()) return
            if (error != null) {
                installBar.setString("Update was not installed")
                updateNote.setText("<html><body style='width:300px'>Could not install: " + CalciumInstaller.escape(error) + "<br>You can retry Install update or use the release page.</body></html>")
                installBtn.setEnabled(true); updateBtn.setEnabled(true)
            } else {
                installBar.setValue(100); installBar.setString("Installed " + completed.version)
                installBtn.setText("Installed — restart Fiji"); installBtn.setEnabled(false)
                updateNote.setText("<html><body style='width:300px'><b>Installed " + completed.version + ".</b> Save your work, then restart Fiji to use it. This session continues with the version already loaded." +
                    (completed.backup ? "<br>The previous file was backed up in calcium-plugin-backups inside Fiji's folder." : "") + "</body></html>")
                updateNote.setToolTipText("Installed file: " + completed.target + (completed.backup ? "; backup: " + completed.backup : ""))
            }
            frame.refreshContents()
        }
    }
}
def checkUpdates = { boolean automatic ->
    if (installBusy.get() || installedVersion != null) return
    boolean accepted = updateCheck.check(automatic, { Map found ->
        if (installBusy.get() || installedVersion != null) return
        updateBtn.setEnabled(true); updateBtn.setText("Check for updates")
        if (automatic && found.error) return
        String message
        if (found.error) message = "Could not check for updates: " + CalciumInstaller.escape(found.error) + ". You can open the release page manually."
        else {
            releasePage = found.page; selectedRelease = found
            installBtn.setVisible(found.newer && found.installable)
            updatesCard.setSummary(found.newer ? "Update available: " + found.version : "Up to date · " + CalciumRelease.VERSION)
            if (found.newer) updatesCard.setVisible(true)
            message = found.newer ? "<b>Update available: " + found.version + ".</b> " + (found.installable ? "Click Install update; downloading, backup and replacement are handled for you." : "This release requires installation from its download page.") : "You are up to date (installed " + CalciumRelease.VERSION + "; latest release " + found.version + ")."
        }
        releaseBtn.setVisible(true)
        updateNote.setText("<html><body style='width:300px'>" + message + "</body></html>")
        frame.refreshContents()
        if (found.prompt) {
            def prompt = new javax.swing.JDialog(frame, "Calcium plugin update available", false)
            def content = new javax.swing.JPanel(); content.setBorder(javax.swing.BorderFactory.createEmptyBorder(12, 12, 12, 12))
            content.setLayout(new javax.swing.BoxLayout(content, javax.swing.BoxLayout.Y_AXIS))
            content.add(note("<b>Version " + found.version + " is available.</b><br>" + (found.installable ? "Install update downloads, verifies and replaces the plugin for you. Restart Fiji after saving your work." : "Open the release page to download this version.")))
            def install = btn(found.installable ? "Install update" : "Open release page"), later = btn("Later")
            content.add(row(install, later)); prompt.add(content)
            install.addActionListener({ e -> prompt.dispose(); if (found.installable) startInstall(found) else openRelease() } as java.awt.event.ActionListener)
            later.addActionListener({ e -> prompt.dispose() } as java.awt.event.ActionListener)
            prompt.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
            prompt.pack(); prompt.setLocationRelativeTo(frame); prompt.setVisible(true)
        }
    })
    if (accepted) {
        updateBtn.setEnabled(false); updateBtn.setText("Checking...")
        if (!automatic) updateNote.setText("Checking the latest GitHub release...")
    }
}
autoUpdateBox.addActionListener({ e ->
    updateCheck.enabled = autoUpdateBox.isSelected()
    Prefs.set("calcium.steps.check_updates", updateCheck.enabled)
    updateBtn.setEnabled(!installBusy.get() && installedVersion == null); updateBtn.setText("Check for updates")
} as java.awt.event.ActionListener)
updateBtn.addActionListener({ e -> checkUpdates(false) } as java.awt.event.ActionListener)
installBtn.addActionListener({ e -> startInstall(selectedRelease) } as java.awt.event.ActionListener)
releaseBtn.addActionListener({ e -> openRelease() } as java.awt.event.ActionListener)
frame.addWindowListener(new java.awt.event.WindowAdapter() {
    void windowClosed(java.awt.event.WindowEvent e) { updateCheck.close() }
    void windowClosing(java.awt.event.WindowEvent e) { updateCheck.close() }
})
def scroll = new javax.swing.JScrollPane(root, javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED, javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED)
scroll.setBorder(javax.swing.BorderFactory.createEmptyBorder()); scroll.getViewport().setBackground(CalciumTheme.PAGE)
scroll.setWheelScrollingEnabled(true)
scroll.getVerticalScrollBar().setUnitIncrement(24)
scroll.getHorizontalScrollBar().setUnitIncrement(24)
frame.getContentPane().add(scroll)
frame.pack(); frame.setAlwaysOnTop(true)
if (!testGui) { def scr = java.awt.Toolkit.getDefaultToolkit().getScreenSize(); frame.setLocation((int) Math.max(0, scr.width - frame.getWidth() - 12), 40) }
if (!testGui) CalciumTheme.beside(steps.colourPic, frame)   // the colour picture takes the free space on the left; this window stays on the right

// ---- what the window does ----
def lastKey = ""; int[][] lastBg = new int[0][]
def proposals = []        // automatic outlines: [roi, grade (0 responds, 1 likely, 2 unsure, 3 bright), ev (the measured evidence), cx, cy]
def cellEdits = new CalciumCellEdits(manager: rm, picture: steps.colourPic)
def measuredKey = { Roi[] rois -> rois.collect { CalciumInstaller.digest(RoiEncoder.saveAsByteArray(it)) }.join(':') }
// Step 4 opens as a page in the browser, served from this computer by the plugin. The Fiji window does the same analysis
// and is used when asked for, or when the page cannot be served or opened here.
CalciumWebServer stepFourPage = null
def openFijiWindow = { -> def window = new CalciumAnalysisWindow(analysisMeasurement); window.besideControls(frame); window.setVisible(true) }
analyzeBtn.addActionListener({ e ->
    if (analysisMeasurement == null) return
    if (testGui) { openFijiWindow(); return }
    try {
        if (stepFourPage == null) { def starting = new CalciumWebServer(); starting.start(); stepFourPage = starting }
        if (stepFourPage.session == null || stepFourPage.session.measurement.key != analysisMeasurement.key) { stepFourPage.session?.cancel(); stepFourPage.session = new CalciumAnalysisSession(analysisMeasurement) }
        java.awt.Desktop.getDesktop().browse(new java.net.URI(stepFourPage.url()))
        s4.setSummary("open in your browser")
    } catch (Throwable problem) {
        IJ.log("Calcium imaging plugin: the analysis page could not be opened in a browser (" + problem + "). Using the Fiji window instead.")
        openFijiWindow()
    }
} as java.awt.event.ActionListener)
classicBtn.addActionListener({ e -> if (analysisMeasurement != null) openFijiWindow() } as java.awt.event.ActionListener)
frame.addWindowListener(new java.awt.event.WindowAdapter() {
    void windowClosed(java.awt.event.WindowEvent e) { stepFourPage?.session?.cancel(); stepFourPage?.stop() }
})
def cellRois = { ->
    def all = rm.getRoisAsArray() as List
    all.findAll { it != null && it.isArea() }.collect { def r = (Roi) it.clone(); r.setPosition(0); r } as Roi[]
}
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
    if (steps.respDiff != null) parts << String.format(Locale.US, "rise after the stimulus %.1f gray levels inside, %.1f around it", e[5], e[6])
    parts << String.format(Locale.US, "%.1f gray levels brighter than its surroundings", e[4])
    if (e[3] >= 0.02) parts << String.format(Locale.US, "%.0f %% of it is very bright", 100 * e[3])
    if (e[7] > 0) parts << "touches the edge of the picture"
    if (e[1] > 2300) parts << "large: check it is one cell, not two"
    "<b>Suggestion: " + GRADE_TEXT[p.grade] + "</b> (" + GRADE_SHARE[p.grade] + ")<br>" + parts.join("; ") + "."
}
def refresh = { boolean force ->
    try {
        def sp = splitCells(); Roi[] all = (Roi[]) sp.all; Roi[] cells = (Roi[]) sp.cells; Roi[] edgeCells = (Roi[]) sp.edge
        def key = all.collect { it.getBounds().toString() + it.getContainedPoints().length }.join("|") + "#" + Arrays.toString(rm.getSelectedIndexes()) + "#" + showBg.isSelected() + "#" + proposals.size() + "#" + showSug.isSelected() + "#" + edgeBox.isSelected() + (analysisMeasurement!=null ? '#'+measuredKey(cells) : '')
        if (!force && key == lastKey) return
        lastKey = key
        if (analysisMeasurement != null) {
            boolean current = measuredKey(cells) == analysisMeasurement.key
            analyzeBtn.setEnabled(current && go.isEnabled())
            s3.setDone(current); s4.setSummary(current ? "ready" : "measure again first")
            analysisNote.setText("<html><body style='width:300px'>" + (current ? "Ready. The analysis opens as a page in your browser; it runs on this computer only." : "Cell outlines changed. Measure and save again before analysis.") + "</body></html>")
        }
        steps.alsoAvoid = edgeCells
        lastBg = steps.backgrounds(cells)
        int none = lastBg.count { it.length == 0 } as int
        def ov = new Overlay()
        if (showBg.isSelected() && cells.length > 0) {
            int sel = rm.getSelectedIndex(); int only = -1
            if (sel >= 0 && rm.getSelectedIndexes().length == 1) { def allR = rm.getRoisAsArray(); int k = -1; for (int i = 0; i <= sel && i < allR.length; i++) if (allR[i].isArea()) k++; only = sp.keptIdx.indexOf(k) }
            ov.add(steps.blueLayer(lastBg, only))
        }
        def left = unusedProposals(all)
        if (showSug.isSelected()) left.each { p ->
            def r = (Roi) p.roi.clone()
            if (p.grade == 0) { r.setStrokeColor(CYAN); r.setStrokeWidth(2) }
            else if (p.grade == 1) { r.setStrokeColor(CYAN); r.setStrokeWidth(1) }
            else if (p.grade == 2) { r.setStrokeColor(GREY); r.setStroke(DASH) }
            else { r.setStrokeColor(ORANGE); r.setStrokeWidth(1) }
            ov.add(r)
        }
        def selectedAreas = []; int areaIndex = 0
        rm.getRoisAsArray().eachWithIndex { r, i -> if (r.isArea()) { if (rm.isSelected(i)) selectedAreas << areaIndex; areaIndex++ } }
        all.eachWithIndex { c, i -> CalciumCellAppearance.add(ov, c, !sp.keptIdx.contains(i), selectedAreas.contains(i)) }
        steps.colourPic.setOverlay(ov)
        int nLikely = left.count { it.grade <= 1 } as int
        int edgeSug = edgeBox.isSelected() ? (proposals.count { it.edge } as int) : 0
        status.setText("<html><body style='width:300px'><b>Cells chosen: " + cells.length + "</b>" + (edgeCells.length > 0 ? " &nbsp;(+" + edgeCells.length + " at the edge, left out)" : "") + (none > 0 ? " &nbsp;(" + none + " without background)" : "") + (proposals.size() > 0 ? "<br><font color='" + MUTED + "'>Suggestions left: " + nLikely + " likely, " + left.count { it.grade == 2 } + " unsure, " + left.count { it.grade == 3 } + " bright" + (edgeSug > 0 ? "; " + edgeSug + " at the edge not shown" : "") + "</font>" : "") + "</body></html>")
        s2.setSummary(cells.length + (cells.length == 1 ? " cell chosen" : " cells chosen")); s2.setDone(cells.length > 0 && !s2.expanded)
        acceptAllBtn.setText("Accept likely (" + nLikely + ")"); acceptAllBtn.setEnabled(nLikely > 0)
        undoBtn.setEnabled(!cellEdits.history.isEmpty())
    } catch (Throwable t) { status.setText("" + t) }
}
def addCell = { Roi r0 -> def r = (Roi) r0.clone(); r.setPosition(0); r.setStrokeColor(null); rm.addRoi(r) }
def acceptAt = { int x, int y ->
    def hit = unusedProposals(cellRois()).findAll { it.roi.contains(x, y) }
    if (hit.isEmpty()) return false
    cellEdits.change('add cell') { addCell(hit.min { it.roi.getBounds().width * it.roi.getBounds().height }.roi) }
    rm.deselect(); steps.colourPic.deleteRoi()
    refresh(true); return true
}
def acceptAll = { -> def todo = unusedProposals(cellRois()).findAll { it.grade <= 1 }; cellEdits.change('accept all likely cells') { todo.each { addCell(it.roi) } }; rm.deselect(); steps.colourPic.deleteRoi(); refresh(true); todo.size() }
def lastHover = null
def hoverAt = { int x, int y ->
    Roi[] cells = cellRois()
    def p = showSug.isSelected() ? unusedProposals(cells).findAll { it.roi.contains(x, y) }.min { it.roi.getBounds().width * it.roi.getBounds().height } : null
    def target = null; def text = "&nbsp;"; def cap = ""
    if (p != null) { target = p.roi; text = evidenceText(p); cap = "this suggestion: gray levels above its background, over time" }
    else {
        int ci = -1; cells.eachWithIndex { c, i -> if (ci < 0 && c.contains(x, y)) ci = i }
        if (ci >= 0) { target = cells[ci]; text = "<b>Your cell " + (ci + 1) + "</b>"; cap = "cell " + (ci + 1) + ": gray levels above its background, over time" }
    }
    if (target == null || target.is(lastHover)) { if (target == null && lastHover != null) { lastHover = null; info.setText("<html><body style='width:300px'>&nbsp;</body></html>") }; return text }
    lastHover = target
    info.setText("<html><body style='width:300px'>" + text + "</body></html>")
    trace.show(steps.traceOf(target, cells), cap, steps.found ? steps.p0 : 0, steps.found ? steps.p1 : 0)
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
def suggest = { CalciumOutlineTask task = new CalciumOutlineTask(progress: { p, s, d -> }) ->
    task.report(0, "Checking Cellpose", "Finding the installed Python environment.")
    def py = findCellposePython()
    if (py == null && cellposeBusy) return "Cellpose is still being set up; try again when it says it is found"
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
    def dir = new File(outDir, "automatic_outlines"); dir.mkdirs()
    def tif = new File(dir, base + "_for_automatic_outlines.tif")
    def zip = new File(dir, base + "_for_automatic_outlines_rois.zip")
    def optFile = new File(dir, "cellpose_options_used.txt")
    def pending = new File(dir, "cellpose_run_incomplete.txt")
    def pic = steps.outlinePicture()
    task.imageSize = pic.getWidth() + " x " + pic.getHeight() + " pixels"
    def opts = []; def cmd = []
    boolean reuse = false
    if (!forceFresh && !pending.exists() && tif.isFile() && zip.isFile() && zip.length() > 0) {
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
        def previousLog = new File(dir, "cellpose_output.txt")
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
        if (known == "") {
            task.report(8, "Checking CPU / graphics support", "Checking whether this Cellpose installation can use a graphics processor.")
            // A script file avoids Windows command-line quoting differences for Python -c.
            def probe = new File(dir, "cellpose_device_check.py")
            probe.text = "import torch\nprint('CHIP=GPU' if (torch.cuda.is_available() or (hasattr(torch.backends, 'mps') and torch.backends.mps.is_available())) else 'CHIP=CPU')\n"
            try {
                def probeLog = new File(dir, "cellpose_device_check.txt")
                int probeCode = task.run([py, "-u", probe.getAbsolutePath()], dir, probeLog)
                String answer = probeLog.getText("UTF-8")
                if (probeCode == 0) known = answer.contains("CHIP=GPU") ? "GPU" : answer.contains("CHIP=CPU") ? "CPU" : ""
                if (known != "") Prefs.set(key, known)
            } catch (java.util.concurrent.CancellationException cancelled) { throw cancelled }
            catch (Exception ignored) { } // Keep the uploaded version's GPU-first behaviour if the check is inconclusive.
        }
        task.checkCancelled()
        task.device = known
        boolean gpu = known != "CPU"
        task.settings = "Attempt 1 of 2; " + (gpu ? "full-size detection" : "faster CPU setting (diameter 60)")
        opts = ["--pretrained_model", "cpsam_v2", "--use_gpu"] + (gpu ? [] : ["--diameter", "60"]) + ["--save_rois", "--no_npy", "--cellprob_threshold", "-2", "--flow_threshold", "0.8", "--verbose"]
        cmd = [py, "-u", "-m", "cellpose", "--image_path", tif.getPath()] + opts
        task.report(10, "Starting Cellpose (attempt 1 of 2)", gpu ? "Using full-size detection. First use may download the model." : "CPU detection uses the faster diameter-60 setting from this version; it may still take several minutes.")
        int code = task.run(cmd, dir, new File(dir, "cellpose_output.txt"))
        if (code != 0 || !zip.isFile() || zip.length() == 0) {
            task.checkCancelled(); removeOldZip()
            task.downloading = false; task.device = "CPU"
            task.settings = "Attempt 2 of 2; default model, full-size detection"
            task.report(task.percent, "Retrying Cellpose on the CPU (attempt 2 of 2)", "The first attempt did not produce outlines. Trying this version's default model and full-size image.")
            opts = ["--save_rois", "--no_npy", "--cellprob_threshold", "-2", "--flow_threshold", "0.8", "--verbose"]
            cmd = [py, "-u", "-m", "cellpose", "--image_path", tif.getPath()] + opts
            code = task.run(cmd, dir, new File(dir, "cellpose_output_second_try.txt"))
        }
        task.checkCancelled()
        if (code != 0 || !zip.isFile() || zip.length() == 0) return "Cellpose made no outlines (exit code " + code + "). Open the detection log: " + task.logFile.getPath()
        optFile.text = opts.join(" ")
        how = "Cellpose took " + (long) ((System.nanoTime() - started) / 1000000000L) + " s" + (task.device ? " (" + task.device + ")" : "") +
            (opts.contains("--diameter") ? "; used the faster CPU diameter-60 setting" : "; used full-size detection")
    }
    task.report(85, "Reading the cell outlines", "Opening the outlines and checking their measurements.")
    def got = []; int hidden = 0
    def zf = new java.util.zip.ZipFile(zip)
    try {
        def entries = zf.entries().toList().findAll { !it.isDirectory() && it.getName().endsWith(".roi") }.sort { it.getName() }
        entries.eachWithIndex { e, i ->
            task.checkCancelled()
            def r = new ij.io.RoiDecoder(zf.getInputStream(e).withCloseable { it.bytes }, e.getName()).getRoi()
            if (r != null && r.isArea()) {
                double[] ev = steps.grade(r)
                if (ev[0] >= 4) hidden++
                else { def cen = r.getContourCentroid(); got << [roi:r, grade:(int) ev[0], ev:ev, cx:(int) cen[0], cy:(int) cen[1], edge:steps.nearEdge(r)] }
            }
            task.report(85 + (int) (12 * (i + 1) / entries.size()), "Checking the suggested cells", "Checked " + (i + 1) + " of " + entries.size() + " outlines.")
        }
    } finally { zf.close() }
    task.report(98, "Saving the suggestion measurements", "Saving the evidence behind each outline's grade.")
    def csv = new StringBuilder("suggestion,grade,area_pixels,mean_brightness,share_very_bright,brighter_than_surroundings,rise_inside,rise_around,touches_edge,x,y\n")
    got.eachWithIndex { p, i -> csv.append(String.format(Locale.US, "%d,%s,%.0f,%.2f,%.3f,%.2f,%.2f,%.2f,%s,%d,%d\n", i + 1, GRADE_TEXT[p.grade], p.ev[1], p.ev[2], p.ev[3], p.ev[4], p.ev[5], p.ev[6], p.ev[7] > 0 ? "yes" : "no", p.cx, p.cy)) }
    new File(dir, base + "_suggestions_with_evidence.csv").text = csv.toString()
    steps.say("Automatic outlines (suggestions only): Cellpose, command: " + cmd.join(" ") + " (" + how + "). " + got.size() + " outlines shown: " + got.count { it.grade == 0 } + " respond on their own, " + got.count { it.grade == 1 } + " likely cells, " + got.count { it.grade == 2 } + " unsure, " + got.count { it.grade == 3 } + " very bright; " + hidden + " left out for their size; " + got.count { it.edge } + " lie within " + CalciumSteps.EDGE_MARGIN + " pixels of the picture's edge" + (edgeBox.isSelected() ? " and are not shown (box ticked)" : "") + ". Grades are measured (brightness, share of very bright pixels, contrast, rise after the stimulus inside against around). Only the outlines the user accepted became cells.")
    task.checkCancelled()
    javax.swing.SwingUtilities.invokeAndWait { proposals.clear(); proposals.addAll(got) }
    if (pending.exists() && !pending.delete()) throw new IOException("Could not finish the detection record: " + pending.getPath())
    task.succeeded = true
    return got.size() + " suggestions; " + how
}


def detectionRunner = new CalciumDetectionRunner()
suggestBtn.addActionListener({ e ->
    forceFresh = !proposals.isEmpty()
    detectionRunner.start(suggestBtn, detectionProgress,
        { task -> suggest(task) },
        { String message ->
            suggestBtn.setText(proposals.isEmpty() ? "Suggest cells automatically" : "Suggest cells again")
            info.setText("<html><body style='width:300px'>Suggestions: " + CalciumOutlineProgress.escape(message) + "</body></html>")
            if (!cellposeBusy) showCellposeState(null)
            refresh(true)
        },
        { -> frame.refreshContents() })
} as java.awt.event.ActionListener)
acceptAllBtn.addActionListener({ e ->
    int nL = unusedProposals(cellRois()).count { it.grade <= 1 } as int
    if (nL == 0) return
    int ans = javax.swing.JOptionPane.showConfirmDialog(frame, "Accept all " + nL + " outlines graded \"responds\" or \"likely a cell\"?\n\nIn tests about 8 in 10 of these were real cells, so about 2 in 10 are not.\nLook through them afterwards and remove the wrong ones.", "Accept all likely ones", javax.swing.JOptionPane.OK_CANCEL_OPTION)
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
openBtn.addActionListener({ e -> try { java.awt.Desktop.getDesktop().open(outDir) } catch (Throwable t) { IJ.log("Results folder: " + outDir.getPath()) } } as java.awt.event.ActionListener)
frame.addWindowListener(new java.awt.event.WindowAdapter() { void windowClosing(java.awt.event.WindowEvent e) { timer.stop(); detectionProgress.close() } })
go.addActionListener({ e ->
    def sp = splitCells(); Roi[] cells = (Roi[]) sp.cells; Roi[] edgeCells = (Roi[]) sp.edge
    steps.alsoAvoid = edgeCells
    steps.edgeNote = edgeBox.isSelected() ? "Cells at the picture's edge: outlines within " + CalciumSteps.EDGE_MARGIN + " pixels of the border are left out, because their background can only be taken from the inner, better-lit side (" + edgeCells.length + " outlines in the ROI Manager left out; they are not in the tables or the saved outlines)." : "Cells at the picture's edge were NOT left out (box unticked); the summary table marks them."
    if (cells.length == 0) { result.setText("<html><body style='width:300px'>There are no cells yet. Choose cells in step 2 first.</body></html>"); return }
    go.setEnabled(false); go.setText("Working..."); analyzeBtn.setEnabled(false); analysisMeasurement = null
    Thread.start {
        try {
            outDir.mkdirs()
            def names = (1..cells.length).collect { "Cell" + String.format("%03d", it) } as String[]
            ImagePlus corrected = null
            if (makeMovie.isSelected()) {
                result.setText("<html><body style='width:300px'>Making the background-corrected movie...</body></html>")
                corrected = steps.correctedMovie(cells)
                if (saveMovie.isSelected()) IJ.saveAsTiff(corrected, outDir.getPath() + File.separator + base + "_background_corrected_movie.tif")
            }
            result.setText("<html><body style='width:300px'>Measuring " + cells.length + " cells on " + steps.n + " frames...</body></html>")
            def summary = steps.measureAndSave(cells, names, outDir.getPath() + File.separator + base + "_")
            if (corrected != null && !testGui) { corrected.show(); corrected.setSlice(Math.min(steps.n, steps.found ? steps.p0 + 5 : 1)); corrected.setDisplayRange(0, CalciumSteps.SHOW_MAX); corrected.updateAndDraw() }
            IJ.log(summary); IJ.log("Saved in: " + outDir.getPath()); IJ.log("Step 4 uses unlifted fluorescence: " + base + "_gray_values_background_corrected_for_SimpleCa.csv")
            double interval = imp.getCalibration().frameInterval
            String unit = imp.getCalibration().getTimeUnit().toLowerCase(Locale.ROOT)
            if (unit in ['ms','msec','millisecond']) interval /= 1000d
            else if (unit in ['min','minute']) interval *= 60d
            else if (!(unit in ['sec','s','second'])) interval = Double.NaN
            final Map measured = [csv:new File(outDir,base + "_gray_values_background_corrected_for_SimpleCa.csv"),folder:outDir,
                recording:base,rois:cells.collect{(Roi)it.clone()} as Roi[],names:names,frames:steps.n,picture:steps.colourPic,movie:imp,
                hz:interval > 0 ? 1d/interval : Double.NaN,baselineStart:steps.found?steps.b0:1,baselineEnd:steps.found?steps.b1:Math.min(20,steps.n-1),key:measuredKey(cells)]
            measured.inputHash=CalciumInstaller.digest(measured.csv.bytes)
            def html = "<b>Done.</b> " + cells.length + " cells measured on " + steps.n + " frames.<br>" +
                "<b>" + steps.lastClear + " of " + cells.length + "</b> show a clear rise; mean rise " + String.format(Locale.US, "%.1f", steps.lastMeanRise) + " gray levels.<br>" +
                (edgeCells.length > 0 ? edgeCells.length + " outlines at the picture's edge were left out (dashed).<br>" : "") +
                (steps.lastDim > 0 ? "<font color='#b05000'><b>" + steps.lastDim + " cells are not brighter than their background at rest.</b> Review their quality status in Step 4.</font><br>" : "") +
                "<b>Continue to Step 4 to analyze responses inside Fiji.</b> It uses unlifted fluorescence and reports excluded and unscorable cells explicitly.<br>Saved in: " + savedIn
            javax.swing.SwingUtilities.invokeLater { analysisMeasurement = measured; analyzeBtn.setEnabled(true); result.setText("<html><body style='width:300px'>" + html + "</body></html>"); go.setText("Measure again"); go.setEnabled(true); s3.setSummary(cells.length + " cells measured"); s4.setExpanded(true); refresh(true); openBtn.setVisible(true); frame.refreshContents() }
        } catch (Throwable t) { go.setEnabled(true); go.setText("Measure and save"); result.setText("<html><body style='width:300px'>Could not measure: " + t + "</body></html>"); IJ.log("Calcium imaging plugin: " + t) }
    }
} as java.awt.event.ActionListener)
if (testGui) { binding.setVariable("hooks", [detectionProgress: detectionProgress, suggestBtn: suggestBtn, edgeBox: edgeBox, addBtn: addBtn, setForceFresh: { boolean v -> forceFresh = v }, setupCellpose: setupCellpose, findCellposePython: findCellposePython, setupBtn: setupBtn, setupRow: setupRow, cpLabel: cpLabel, showCellposeState: showCellposeState, setupText: SETUP_TEXT, saveMovie: saveMovie, makeMovie: makeMovie, rm: rm, refresh: refresh, go: go, status: status, steps: steps, outDir: outDir, showBg: showBg, frame: frame, suggest: suggest, acceptAt: acceptAt, acceptAll: acceptAll, hoverAt: hoverAt, proposals: proposals, result: result, info: info, removeSelected: removeSelected, undoLast: undoLast, root: root]); refresh(true); return }
frame.setVisible(true)
timer.start()
refresh(true)
javax.swing.SwingUtilities.invokeLater { checkUpdates(true) }

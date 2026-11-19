// Calcium imaging: cells step by step                                                    version 2026.11.19.8
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
        say("Calcium imaging plugin, version 2026.11.19.8")
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
    static final String VERSION = "2026.11.19.8"
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

/** One compact Step 2 panel, updated only on Swing's event thread. */
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
        bar.setStringPainted(true); bar.setString("0% - Stage progress"); bar.setAlignmentX(0f)
        bar.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, 22))
        bar.setToolTipText("Progress through detection stages, not a time estimate. It can pause while Cellpose runs.")
        text.setAlignmentX(0f); text.setFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 12))
        // Reserve space for a quiet/long-running message without growing beyond the screen mid-run.
        text.setPreferredSize(new java.awt.Dimension(400, 180))
        text.setMinimumSize(new java.awt.Dimension(400, 180))
        text.setMaximumSize(new java.awt.Dimension(Integer.MAX_VALUE, 180))
        def buttons = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 4, 2)); buttons.setAlignmentX(0f)
        [cancelButton, detailsButton, logButton].each { it.setFont(text.getFont()); it.setMargin(new java.awt.Insets(2, 8, 2, 8)); buttons.add(it) }
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
        bar.setValue(0); cancelButton.setText("Cancel detection"); cancelButton.setEnabled(true)
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
        text.setText("<html><body style='width:300px'><b>" + escape(state.stage) + "</b><br>" + escape(elapsed) + "<br>" + escape(context) + "<br>" + escape(message) +
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
class CalciumAnalysisPlot extends javax.swing.JPanel {
    double[] time=[], mean=[], sem=[], selected=[]
    String caption='Run analysis to see traces', cell=''
    double baselineStart=0, baselineEnd=0, stimulusStart=0, stimulusEnd=0
    CalciumAnalysisPlot() {setPreferredSize(new java.awt.Dimension(650,300)); setBackground(Color.white)}
    protected void paintComponent(java.awt.Graphics graphics) {
        super.paintComponent(graphics)
        def g=(java.awt.Graphics2D)graphics.create()
        try {
            g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
            g.setColor(Color.darkGray); g.drawString(caption.take(105),15,20)
            if(time.length<2) return
            int l=62,t=38,w=Math.max(1,getWidth()-82),h=Math.max(1,getHeight()-90)
            def values=[]; mean.eachWithIndex {v,i->if(Double.isFinite(v)) {values<<v;if(i<sem.length&&Double.isFinite(sem[i])) {values<<v-sem[i];values<<v+sem[i]}}}
            selected.each {if(Double.isFinite(it)) values<<it}; values<<0d
            double lo=values.min() as double,hi=values.max() as double
            if(hi-lo<0.001) hi=lo+1
            double pad=(hi-lo)*0.08;lo-=pad;hi+=pad
            def X={double x->l+(int)((x-time[0])/(time[-1]-time[0])*w)}
            def Y={double y->t+(int)((hi-y)/(hi-lo)*h)}
            g.setColor(new Color(235,235,235));g.fillRect(X(baselineStart),t,Math.max(1,X(baselineEnd)-X(baselineStart)),h)
            g.setColor(new Color(255,238,215));g.fillRect(X(stimulusStart),t,Math.max(1,X(stimulusEnd)-X(stimulusStart)),h)
            g.setColor(Color.lightGray);g.drawLine(l,Y(0d),l+w,Y(0d));g.drawRect(l,t,w,h)
            g.setColor(new Color(70,130,210,65))
            for(int i=1;i<time.length;i++) if(i<mean.length&&i<sem.length&&[mean[i-1],mean[i],sem[i-1],sem[i]].every{Double.isFinite((double)it)}) {
                g.fillPolygon([X(time[i-1]),X(time[i]),X(time[i]),X(time[i-1])] as int[],[Y(mean[i-1]-sem[i-1]),Y(mean[i]-sem[i]),Y(mean[i]+sem[i]),Y(mean[i-1]+sem[i-1])] as int[],4)
            }
            [[mean,new Color(35,90,175),2f],[selected,new Color(160,30,130),2f]].each { entry ->
                double[] y=entry[0] as double[];g.setColor((Color)entry[1]);g.setStroke(new java.awt.BasicStroke((float)entry[2]))
                for(int i=1;i<Math.min(time.length,y.length);i++) if(Double.isFinite(y[i-1])&&Double.isFinite(y[i])) g.drawLine(X(time[i-1]),Y(y[i-1]),X(time[i]),Y(y[i]))
            }
            g.setColor(Color.darkGray);g.drawString(String.format(Locale.US,'%.3g',hi),4,t+10);g.drawString(String.format(Locale.US,'%.3g',lo),4,t+h)
            g.drawString('0',l,t+h+18);g.drawString(String.format(Locale.US,'%.3g s',time[-1]),l+w-40,t+h+18)
            g.drawString('ΔF/F₀: blue = cell mean ± SEM; purple = '+cell, l,getHeight()-12)
        } finally {g.dispose()}
    }
}
class CalciumAnalysisHeatmap extends javax.swing.JPanel {
    double[][] traces=[]; List<String> names=[]; Closure choose={int i->}; double[] time=[]
    CalciumAnalysisHeatmap() {
        setBackground(Color.white);setPreferredSize(new java.awt.Dimension(650,260))
        addMouseListener(new java.awt.event.MouseAdapter(){void mouseClicked(java.awt.event.MouseEvent e){
            if(e.y>=30&&e.y<getHeight()-25&&!names.isEmpty()) choose.call(Math.min(names.size()-1,(int)((e.y-30)*names.size()/(getHeight()-55))))
        }})
    }
    protected void paintComponent(java.awt.Graphics graphics) {
        super.paintComponent(graphics);def g=(java.awt.Graphics2D)graphics.create()
        try {
            g.setColor(Color.darkGray);g.drawString('ΔF/F₀ heatmap · blue < 0 · white = 0 · red > 0 · gray = unavailable; click a row',10,18)
            if(names.isEmpty()||time.length<2)return
            double scale=0
            traces.each{it.each{v->if(Double.isFinite(v))scale=Math.max(scale,Math.abs(v))}}
            if(scale==0)scale=1
            int w=Math.max(1,getWidth()-85),h=Math.max(1,getHeight()-55)
            for(int c=0;c<names.size();c++) {
                int y=30+(int)(h*c/names.size()),y2=30+(int)(h*(c+1)/names.size())
                if(y2-y>=11){g.setColor(Color.darkGray);g.drawString(names[c],2,y+10)}
                for(int f=0;f<time.length;f++) {
                    double v=traces[c][f];int fade=Double.isFinite(v)?(int)(255*(1-Math.min(1,Math.abs(v)/scale))):170
                    g.setColor(!Double.isFinite(v)?new Color(170,170,170):v>=0?new Color(255,fade,fade):new Color(fade,fade,255))
                    int x=75+(int)(w*f/time.length),x2=75+(int)(w*(f+1)/time.length)
                    g.fillRect(x,y,Math.max(1,x2-x),Math.max(1,y2-y))
                }
            }
            g.setColor(Color.darkGray);g.drawString('0 s',75,getHeight()-5);g.drawString(String.format(Locale.US,'%.3g s; scale ±%.3g',time[-1],scale),getWidth()-190,getHeight()-5)
        } finally {g.dispose()}
    }
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
    def setup=new javax.swing.JButton('Set up analysis engine (one time)')
    def runButton=new javax.swing.JButton('Analyze responses')
    def cancelButton=new javax.swing.JButton('Cancel')
    def logButton=new javax.swing.JButton('View log')
    def folderButton=new javax.swing.JButton('Open analysis folder')
    def savePlots=new javax.swing.JButton('Save plots as PNG')
    def status=new javax.swing.JLabel('Measure cells, check the settings, then analyze.')
    def progress=new javax.swing.JProgressBar(0,100)
    def confirmed=new javax.swing.JCheckBox('I checked the frame rate, baseline and stimulus times',false)
    def kcl=new javax.swing.JCheckBox('Also score a separate KCl window',false)
    def cells=new javax.swing.JComboBox<String>()
    def plot=new CalciumAnalysisPlot()
    def heatmap=new CalciumAnalysisHeatmap()
    def metrics=new javax.swing.JTable()
    def calls=new javax.swing.JTable()
    def counts=new javax.swing.JTable()
    def tabs=new javax.swing.JTabbedPane()
    Map<String,double[]> traces=[:]
    java.awt.event.MouseListener pictureListener
    javax.swing.Timer clock
    long started
    String stage=''
    CalciumAnalysisWindow(Map data) {
        super('Step 4 · SimpleCa²⁺ analysis in Fiji');measurement=new LinkedHashMap(data)
        if(data.picture!=null) {
            ImagePlus review=new ImagePlus('Analysis cells · '+data.recording,(data.picture as ImagePlus).getProcessor().duplicate())
            def overlay=new Overlay();data.rois.each{CalciumCellAppearance.add(overlay,it,false,false)};review.setOverlay(overlay);review.show()
            measurement.picture=review
        }
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
        def root=new javax.swing.JPanel(new java.awt.BorderLayout(6,6));root.setBorder(javax.swing.BorderFactory.createEmptyBorder(8,8,8,8))
        def settings=new javax.swing.JPanel(new java.awt.GridLayout(0,4,6,4))
        int n=data.frames as int, end=Math.min(n-1,Math.max(3,(data.baselineEnd?:20) as int))
        double rate=(data.hz?:Double.NaN) as double
        def add={String key,String label,Object value->settings.add(new javax.swing.JLabel(label));def field=new javax.swing.JTextField(''+value,9);fields[key]=field;settings.add(field)}
        add('sampling_hz','Frames per second',Double.isFinite(rate)&&rate>0?rate:'')
        add('baseline_start','Baseline first frame',data.baselineStart?:1)
        add('baseline_end','Baseline last frame',end)
        add('stimulus','Stimulus name','Stimulus')
        add('stimulus_start','Stimulus start (s)',Double.isFinite(rate)&&rate>0?(end-1)/rate:'')
        add('stimulus_end','Stimulus end (s)',Double.isFinite(rate)&&rate>0?(n-1)/rate:'')
        add('response_baseline_seconds','Prestimulus baseline (s)',Double.isFinite(rate)&&rate>0?(end-1)/rate:'')
        add('amplitude_floor','Minimum rise (fraction)',0.10)
        add('noise_multiplier','Baseline SD multiplier',3)
        add('consecutive_samples','Consecutive samples',2)
        add('kcl_start','KCl start (s)','');add('kcl_end','KCl end (s)','')
        def top=new javax.swing.JPanel();top.setLayout(new javax.swing.BoxLayout(top,javax.swing.BoxLayout.Y_AXIS))
        top.add(new javax.swing.JLabel('Uses unlifted, background-corrected fluorescence. No second background subtraction.'))
        top.add(settings);top.add(kcl);top.add(confirmed)
        def buttons=new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT))
        [setup,runButton,cancelButton,logButton].each{buttons.add(it)};top.add(buttons)
        progress.setStringPainted(true);top.add(progress);top.add(status)
        root.add(top,java.awt.BorderLayout.NORTH)
        def tracePanel=new javax.swing.JPanel(new java.awt.BorderLayout());tracePanel.add(cells,java.awt.BorderLayout.NORTH);tracePanel.add(plot)
        tabs.addTab('Traces',tracePanel);tabs.addTab('Heatmap',new javax.swing.JScrollPane(heatmap))
        [metrics,calls,counts].each{it.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);it.setAutoCreateRowSorter(true);it.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION)}
        tabs.addTab('Metrics + quality',new javax.swing.JScrollPane(metrics));tabs.addTab('Cell response calls',new javax.swing.JScrollPane(calls));tabs.addTab('Responder counts',new javax.swing.JScrollPane(counts))
        root.add(tabs,java.awt.BorderLayout.CENTER)
        def bottom=new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT));bottom.add(folderButton);bottom.add(savePlots)
        bottom.add(new javax.swing.JLabel('Counts describe chosen cells; cell N is not animal N.'));root.add(bottom,java.awt.BorderLayout.SOUTH)
        setContentPane(new javax.swing.JScrollPane(root));root.setMinimumSize(new java.awt.Dimension(740,650));root.setPreferredSize(new java.awt.Dimension(850,740))
        setSize(920,800);setMinimumSize(new java.awt.Dimension(360,300))
        def bounds=getGraphicsConfiguration().getBounds();setSize(Math.min(getWidth(),(int)bounds.getWidth()-50),Math.min(getHeight(),(int)bounds.getHeight()-80))
        cancelButton.enabled=false;folderButton.enabled=false;savePlots.enabled=false;logButton.enabled=false
        setup.addActionListener({beginSetup()} as java.awt.event.ActionListener)
        runButton.addActionListener({beginAnalysis()} as java.awt.event.ActionListener)
        cancelButton.addActionListener({job?.cancel();stage='Cancelling...';status.setText(stage)} as java.awt.event.ActionListener)
        logButton.addActionListener({if(job?.logFile?.isFile())IJ.open(job.logFile.path)} as java.awt.event.ActionListener)
        folderButton.addActionListener({if(results!=null)java.awt.Desktop.getDesktop().open(results)} as java.awt.event.ActionListener)
        savePlots.addActionListener({exportPlots()} as java.awt.event.ActionListener)
        cells.addActionListener({selectCell(cells.selectedItem as String,true)} as java.awt.event.ActionListener)
        heatmap.choose={int i->cells.setSelectedItem(heatmap.names[i]);tabs.selectedIndex=0}
        [metrics,calls].each{table->table.selectionModel.addListSelectionListener({e->
            if(!e.valueIsAdjusting&&table.selectedRow>=0) {
                int col=(0..<table.columnCount).find{table.getColumnName(it)=='Cell'} ?: 0
                cells.setSelectedItem(table.getValueAt(table.selectedRow,col).toString())
            }
        } as javax.swing.event.ListSelectionListener)}
        def canvas=(measurement.picture as ImagePlus)?.getCanvas()
        if(canvas!=null) {
            pictureListener=new java.awt.event.MouseAdapter(){void mouseReleased(java.awt.event.MouseEvent e){
                if(!javax.swing.SwingUtilities.isLeftMouseButton(e)||busy)return
                int x=canvas.offScreenX(e.x),y=canvas.offScreenY(e.y)
                int i=measurement.rois.findIndexOf{it.contains(x,y)}
                if(i>=0)javax.swing.SwingUtilities.invokeLater{if(!closed)cells.setSelectedItem(measurement.names[i])}
            }}
            canvas.addMouseListener(pictureListener)
        }
        clock=new javax.swing.Timer(500,{e->if(busy)status.setText(stage+' · '+(long)((System.currentTimeMillis()-started)/1000)+' s elapsed')} as java.awt.event.ActionListener)
        addWindowListener(new java.awt.event.WindowAdapter(){void windowClosed(java.awt.event.WindowEvent e){closed=true;job?.cancel();clock.stop();if(canvas!=null)canvas.removeMouseListener(pictureListener);(measurement.picture as ImagePlus)?.close()}})
    }
    Map settings() {
        def values=fields.collectEntries{k,v->[k,v.text.trim()]}
        if(!confirmed.selected)throw new IllegalArgumentException('Confirm the frame rate, baseline and stimulus times first.')
        def number={key->double v=CalciumCsv.number(values[key]);if(!Double.isFinite(v))throw new IllegalArgumentException('Enter a number for '+key);v}
        double hz=number('sampling_hz'),first=number('baseline_start'),last=number('baseline_end'),start=number('stimulus_start'),end=number('stimulus_end'),pre=number('response_baseline_seconds')
        if(hz<=0||first<1||first!=Math.floor(first)||last!=Math.floor(last)||last<first||last>=measurement.frames)throw new IllegalArgumentException('Check the positive frame rate and baseline frame range.')
        if(start<0||end<=start||pre<=0||start-pre< -1e-8||end>(measurement.frames-1)/hz+1.5/hz)throw new IllegalArgumentException('The stimulus and prestimulus baseline must fit this recording.')
        if(!values.stimulus||number('amplitude_floor')<=0||number('noise_multiplier')<0||number('consecutive_samples')<1||number('consecutive_samples')!=Math.floor(number('consecutive_samples')))throw new IllegalArgumentException('Check the stimulus name and response thresholds.')
        if(kcl.selected){number('kcl_start');number('kcl_end')}
        values.kcl_enabled=''+kcl.selected;values.recording=measurement.recording;values.input_mode='background_corrected_fluorescence'
        if(measurement.inputHash)values.expected_input_sha256=measurement.inputHash
        values
    }
    void launch(boolean installing,Closure action) {
        if(busy)return
        busy=true;started=System.currentTimeMillis();stage=installing?'Preparing R setup':'Preparing analysis'
        setup.enabled=false;runButton.enabled=false;cancelButton.enabled=true;progress.value=0;progress.indeterminate=installing
        fields.values().each{it.enabled=false};confirmed.enabled=false;kcl.enabled=false;savePlots.enabled=false
        job=new CalciumLocalJob(logFile:new File(runtime.home,installing?'setup.log':'analysis.log'))
        job.progress={int p,String message->javax.swing.SwingUtilities.invokeLater{if(!closed){progress.value=p;progress.setString(p+'% · '+message);stage=message;logButton.enabled=job.logFile.isFile()}}}
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
                if(failure){status.setText(failure);progress.setString('Not completed');return}
                if(installing){progress.value=100;status.setText('R is ready. Check the settings, then Analyze responses.');return}
                try {loadResults(answer);status.setText('Complete. Tables and settings saved in '+answer.name);progress.value=100;progress.setString('Complete')}
                catch(Throwable t){status.setText('Could not display results: '+t.message);progress.setString('Display error; see saved outputs')}
            }
        }
    }
    void beginSetup(){launch(true){task->runtime.setup(task)}}
    void beginAnalysis(){
        try {Map s=settings();launch(false){task->performAnalysis.call(runtime,measurement.csv,s,measurement.folder,task)}}
        catch(Throwable t){status.setText(t.message?:t.toString())}
    }
    void table(javax.swing.JTable target,File file) {
        def rows=CalciumCsv.read(file)
        def model=new javax.swing.table.DefaultTableModel(rows.drop(1).collect{it as Object[]} as Object[][],rows[0] as Object[]){boolean isCellEditable(int r,int c){false}}
        target.setModel(model);for(int i=0;i<target.columnCount;i++)target.columnModel.getColumn(i).preferredWidth=140
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
        heatmap.names=measurement.names as List;heatmap.time=plot.time
        heatmap.traces=heatmap.names.collect{traces[it]?:new double[plot.time.length].collect{Double.NaN} as double[]} as double[][]
        heatmap.setPreferredSize(new java.awt.Dimension(650,Math.max(260,heatmap.names.size()*12+55)));heatmap.revalidate();heatmap.repaint()
        cells.removeAllItems();measurement.names.each{cells.addItem(it)}
        folderButton.enabled=true;savePlots.enabled=true
    }
    void selectCell(String name,boolean highlight) {
        if(name==null)return
        plot.cell=name;plot.selected=traces[name]?:new double[0]
        plot.caption=name+(traces.containsKey(name)?' · gray: baseline; pale orange: stimulus window':' · unavailable: inspect Metrics + quality for the exclusion reason')
        plot.repaint()
        int index=measurement.names.findIndexOf{it==name}
        if(highlight&&index>=0&&measurement.picture!=null){Roi roi=(Roi)measurement.rois[index].clone();roi.setStrokeColor(Color.white);roi.setStrokeWidth(4);(measurement.picture as ImagePlus).setRoi(roi)}
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

// BEGIN GENERATED SIMPLECA BUNDLE
class CalciumSimpleCaBundle {
    static final String HASH = "d4a92a3143a7c40ba7e6b8f200c9c9076aa1bcc8b275295aa18413ef9c886676"
    static final String DATA = "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"
}
// END GENERATED SIMPLECA BUNDLE

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
    CalciumTracePanel() { setPreferredSize(new java.awt.Dimension(330, 96)); setBackground(Color.white) }
    void show(double[] v, String cap, int b0, int b1) { vals = v; caption = cap; band0 = b0; band1 = b1; repaint() }
    protected void paintComponent(java.awt.Graphics g0) {
        super.paintComponent(g0)
        def g = (java.awt.Graphics2D) g0
        g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
        int W = getWidth(), H = getHeight(), L = 30, T = 16, B = 6
        g.setFont(getFont().deriveFont(10f)); g.setColor(new Color(90, 90, 90)); g.drawString(caption, 6, 11)
        g.setColor(new Color(225, 225, 225)); g.drawRect(L, T, W - L - 6, H - T - B)
        if (vals == null || vals.length < 2) return
        double lo = vals.min(), hi = vals.max(); if (hi - lo < 1) { hi = lo + 1 }
        def X = { int f -> L + (int) ((W - L - 6) * f / (double) (vals.length - 1)) }
        def Y = { double v -> T + (int) ((H - T - B) * (1 - (v - lo) / (hi - lo))) }
        if (band1 > band0) { g.setColor(new Color(255, 225, 255)); g.fillRect(X(band0 - 1), T + 1, Math.max(2, X(band1 - 1) - X(band0 - 1)), H - T - B - 1) }
        if (lo < 0 && hi > 0) { g.setColor(new Color(200, 200, 200)); g.drawLine(L, Y(0d), W - 6, Y(0d)) }
        g.setColor(new Color(40, 40, 40))
        for (int f = 1; f < vals.length; f++) g.drawLine(X(f - 1), Y(vals[f - 1]), X(f), Y(vals[f]))
        g.setColor(new Color(110, 110, 110)); g.drawString(String.format(Locale.US, "%.0f", hi), 4, T + 9); g.drawString(String.format(Locale.US, "%.0f", lo), 4, H - B)
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

def gd = new GenericDialog("Calcium imaging: cells step by step  (version 2026.11.19.8)")
gd.addMessage("STEP 1 of 4:  choose the recording\n \nPick the ORIGINAL movie, not a processed copy.")
gd.addChoice("Recording:", labels as String[], pre)
gd.addCheckbox("Look for a response and show it in magenta", true)
gd.addDirectoryField("Folder for the results:", defaultDir, 40)
gd.addMessage("A subfolder with the recording's name is made inside that folder.\nNext: the cells are made visible and you choose them (step 2), then they are measured (step 3).")
gd.setOKLabel("Next")
gd.showDialog()
if (gd.wasCanceled()) return
def picked = gd.getNextChoice()
boolean look = gd.getNextBoolean()
def outRoot = gd.getNextString()
if (outRoot == null || outRoot.trim().isEmpty()) outRoot = defaultDir
Prefs.set("calcium.steps.folder", outRoot)

ImagePlus imp = null
if (picked == fromDisk) {
    def od = new OpenDialog("Choose the recording (.zvi, .czi, .tif ...)")
    if (od.getFileName() == null) return
    IJ.run("Bio-Formats Importer", "open=[" + od.getPath() + "] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT")
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
def small = new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 12), bold = new java.awt.Font("SansSerif", java.awt.Font.BOLD, 13)
def lineIcon = { Color c, int width, boolean dashed, boolean filled ->
    def bi = new java.awt.image.BufferedImage(30, 14, java.awt.image.BufferedImage.TYPE_INT_ARGB); def g = bi.createGraphics()
    g.setColor(new Color(40, 40, 40)); g.fillRect(0, 0, 30, 14); g.setColor(c)
    if (filled) g.fillRect(4, 3, 22, 8)
    else { g.setStroke(dashed ? DASH : new java.awt.BasicStroke((float) width)); g.drawOval(5, 2, 20, 10) }
    g.dispose(); new javax.swing.ImageIcon(bi)
}
def legendRow = { javax.swing.Icon ic, String s -> def l = new javax.swing.JLabel(s, ic, javax.swing.SwingConstants.LEFT); l.setFont(small); l.setIconTextGap(8); l.setAlignmentX(0f); l }
def section = { String title ->
    def p = new javax.swing.JPanel(); p.setLayout(new javax.swing.BoxLayout(p, javax.swing.BoxLayout.Y_AXIS)); p.setAlignmentX(0f)
    def tb = javax.swing.BorderFactory.createTitledBorder(javax.swing.BorderFactory.createLineBorder(new Color(170, 170, 170)), title); tb.setTitleFont(bold)
    p.setBorder(javax.swing.BorderFactory.createCompoundBorder(tb, javax.swing.BorderFactory.createEmptyBorder(2, 6, 6, 6))); p
}
def row = { java.awt.Component... cs -> def p = new javax.swing.JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 4, 2)); cs.each { p.add(it) }; p.setAlignmentX(0f); p }
def note = { String s -> def l = new javax.swing.JLabel("<html><body style='width:300px'>" + s + "</body></html>"); l.setFont(small); l.setAlignmentX(0f); l }
def btn = { String s -> def b = new javax.swing.JButton(s); b.setFont(small); b.setMargin(new java.awt.Insets(2, 8, 2, 8)); b }
def box = { String s, boolean on -> def b = new javax.swing.JCheckBox(s, on); b.setFont(small); b.setAlignmentX(0f); b }

def frame = new CalciumWindow("Calcium imaging: cells step by step")
frame.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
def root = new CalciumScrollPanel(); root.setBorder(javax.swing.BorderFactory.createEmptyBorder(8, 8, 8, 8))

def s1 = section("Step 1   Recording   ✓")
s1.add(note("<b>" + imp.getTitle() + "</b><br>" + steps.n + " frames. " + (steps.found ? "Response found: brightness rises near frame " + steps.onset + " (shown in magenta)." : (look ? "No response found: the picture has no magenta. Every cell is still shown in green." : "Response not looked for."))))

def s2 = section("Step 2   Choose the cells")
def suggestBtn = btn("Suggest cells automatically")
def acceptAllBtn = btn("Accept all likely ones"); acceptAllBtn.setEnabled(false)
def ovalBtn = btn("Oval"), freeBtn = btn("Freehand"), addBtn = btn("Add drawn outline  (t)")
def removeBtn = btn("Remove selected cell"), undoBtn = btn("Undo last change"); undoBtn.setEnabled(false)
def removeMode = box("Click cells to remove", false)
def editHint = note("Click a white outline to select it, then Remove. Or turn on <b>Click cells to remove</b>.")
def showSug = box("Show suggestions", true), showBg = box("Show background (blue)", true)
def edgeBox = box("Leave out cells at the picture's edge (within " + CalciumSteps.EDGE_MARGIN + " px)", true)
edgeBox.setToolTipText("<html>At the edge a cell's background can only be taken from the inner, better-lit side,<br>so the cell comes out too dim. Untick to suggest, accept and measure edge cells like the others.</html>")
def status = new javax.swing.JLabel("Cells chosen: 0"); status.setFont(bold); status.setAlignmentX(0f)
def info = note("&nbsp;")
def trace = new CalciumTracePanel(); trace.setAlignmentX(0f); trace.setFont(small)
s2.add(note("<b>Click:</b> press the button, then click inside each outline that is a real cell."))
def cpLabel = note("&nbsp;")
def setupBtn = btn("Set up Cellpose (one time)"); def setupRow = row(setupBtn)
def detectionProgress = new CalciumOutlineProgress()
s2.add(row(suggestBtn)); s2.add(detectionProgress); s2.add(cpLabel); s2.add(setupRow); s2.add(row(acceptAllBtn))
s2.add(note("<b>Or draw:</b> pick a tool, draw around a cell, then add it."))
s2.add(row(ovalBtn, freeBtn, addBtn))
def legend = new javax.swing.JPanel(new java.awt.GridLayout(0, 2, 6, 1)); legend.setAlignmentX(0f)
def acceptedIcon = new java.awt.image.BufferedImage(30, 14, java.awt.image.BufferedImage.TYPE_INT_ARGB)
def acceptedGraphics = acceptedIcon.createGraphics()
acceptedGraphics.setColor(new Color(40, 40, 40)); acceptedGraphics.fillRect(0, 0, 30, 14)
acceptedGraphics.setColor(Color.black); acceptedGraphics.setStroke(new java.awt.BasicStroke(6f)); acceptedGraphics.drawOval(5, 3, 20, 8)
acceptedGraphics.setColor(Color.white); acceptedGraphics.setStroke(new java.awt.BasicStroke(2.5f)); acceptedGraphics.drawOval(5, 3, 20, 8); acceptedGraphics.dispose()
legend.add(legendRow(new javax.swing.ImageIcon(acceptedIcon), "your cells: white + black border"))
legend.add(legendRow(lineIcon(BLUE, 1, false, true), "background used"))
legend.add(legendRow(lineIcon(CYAN, 2, false, false), "suggested: responds"))
legend.add(legendRow(lineIcon(CYAN, 1, false, false), "suggested: likely a cell"))
legend.add(legendRow(lineIcon(GREY, 1, true, false), "suggested: unsure"))
legend.add(legendRow(lineIcon(ORANGE, 1, false, false), "probably precipitate"))
s2.add(legend)
s2.add(row(showSug, showBg)); s2.add(row(edgeBox))
s2.add(status); s2.add(row(removeBtn, undoBtn)); s2.add(row(removeMode)); s2.add(editHint)
s2.add(info); s2.add(trace)

def s3 = section("Step 3   Measure and save")
def makeMovie = box("Also open a background-corrected copy of the movie", true)
def saveMovie = box("... and save that copy as a TIFF (about 1.7 GB for 300 frames)", false)
def go = new javax.swing.JButton("Measure and save"); go.setFont(bold)
def result = note("Saves the tables in: ..." + File.separator + outDir.getParentFile().getName() + File.separator + outDir.getName()); result.setToolTipText(outDir.getPath())
def openBtn = btn("Open the results folder"); openBtn.setVisible(false)
def closeBtn = btn("Close")
s3.add(makeMovie); s3.add(saveMovie); s3.add(row(go, openBtn, closeBtn)); s3.add(result)
def s4 = section("Step 4   Analyze responses with SimpleCa²⁺")
def analyzeBtn = btn("Analyze responses in Fiji"); analyzeBtn.setEnabled(false)
def analysisNote = note("Measure and save first. Then review traces, metrics and responder counts here. One-time engine setup needs internet.")
Map analysisMeasurement = null
s4.add(row(analyzeBtn)); s4.add(analysisNote)
[s1, s2, s3, s4].each { root.add(it); root.add(javax.swing.Box.createVerticalStrut(6)) }
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
root.add(row(new javax.swing.JLabel("Version " + CalciumRelease.VERSION), updateBtn))
root.add(row(installBtn, releaseBtn)); root.add(autoUpdateBox); root.add(updateNote); root.add(installBar)
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
scroll.setBorder(javax.swing.BorderFactory.createEmptyBorder())
scroll.setWheelScrollingEnabled(true)
scroll.getVerticalScrollBar().setUnitIncrement(24)
scroll.getHorizontalScrollBar().setUnitIncrement(24)
frame.getContentPane().add(scroll)
frame.pack(); frame.setAlwaysOnTop(true)
if (!testGui) { def scr = java.awt.Toolkit.getDefaultToolkit().getScreenSize(); frame.setLocation((int) Math.max(0, scr.width - frame.getWidth() - 12), 40) }

// ---- what the window does ----
def lastKey = ""; int[][] lastBg = new int[0][]
def proposals = []        // automatic outlines: [roi, grade (0 responds, 1 likely, 2 unsure, 3 bright), ev (the measured evidence), cx, cy]
def cellEdits = new CalciumCellEdits(manager: rm, picture: steps.colourPic)
def measuredKey = { Roi[] rois -> rois.collect { CalciumInstaller.digest(RoiEncoder.saveAsByteArray(it)) }.join(':') }
analyzeBtn.addActionListener({ e ->
    if (analysisMeasurement != null) { def window = new CalciumAnalysisWindow(analysisMeasurement); window.setVisible(true) }
} as java.awt.event.ActionListener)
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
            analysisNote.setText("<html><body style='width:300px'>" + (current ? "Ready: analyze the last measurement in Fiji." : "Cell outlines changed. Measure and save again before analysis.") + "</body></html>")
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
        status.setText("Cells chosen: " + cells.length + (edgeCells.length > 0 ? "   (+" + edgeCells.length + " at the edge, left out)" : "") + (none > 0 ? "   (" + none + " without background)" : "") + (proposals.size() > 0 ? "      suggestions left: " + nLikely + " likely, " + left.count { it.grade == 2 } + " unsure, " + left.count { it.grade == 3 } + " bright" + (edgeSug > 0 ? "; " + edgeSug + " at the edge not shown" : "") : ""))
        acceptAllBtn.setText("Accept all likely ones (" + nLikely + ")"); acceptAllBtn.setEnabled(nLikely > 0)
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
    cpLabel.setText("<html><body style='width:300px'><font color='#666666'>" + (busyText != null ? busyText : found ? "Cellpose: found on this computer." : "Cellpose (the program that suggests outlines) is not installed on this computer yet.") + "</font></body></html>")
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
                recording:base,rois:cells.collect{(Roi)it.clone()} as Roi[],names:names,frames:steps.n,picture:steps.colourPic,
                hz:interval > 0 ? 1d/interval : Double.NaN,baselineStart:steps.found?steps.b0:1,baselineEnd:steps.found?steps.b1:Math.min(20,steps.n-1),key:measuredKey(cells)]
            measured.inputHash=CalciumInstaller.digest(measured.csv.bytes)
            def html = "<b>Done.</b> " + cells.length + " cells measured on " + steps.n + " frames.<br>" +
                "<b>" + steps.lastClear + " of " + cells.length + "</b> show a clear rise; mean rise " + String.format(Locale.US, "%.1f", steps.lastMeanRise) + " gray levels.<br>" +
                (edgeCells.length > 0 ? edgeCells.length + " outlines at the picture's edge were left out (dashed).<br>" : "") +
                (steps.lastDim > 0 ? "<font color='#b05000'><b>" + steps.lastDim + " cells are not brighter than their background at rest.</b> Review their quality status in Step 4.</font><br>" : "") +
                "<b>Continue to Step 4 to analyze responses inside Fiji.</b> It uses unlifted fluorescence and reports excluded and unscorable cells explicitly.<br>Saved in: ..." + File.separator + outDir.getParentFile().getName() + File.separator + outDir.getName()
            javax.swing.SwingUtilities.invokeLater { analysisMeasurement = measured; analyzeBtn.setEnabled(true); result.setText("<html><body style='width:300px'>" + html + "</body></html>"); go.setText("Measure and save again"); go.setEnabled(true); refresh(true); openBtn.setVisible(true); frame.refreshContents() }
        } catch (Throwable t) { go.setEnabled(true); go.setText("Measure and save"); result.setText("<html><body style='width:300px'>Could not measure: " + t + "</body></html>"); IJ.log("Calcium imaging plugin: " + t) }
    }
} as java.awt.event.ActionListener)
if (testGui) { binding.setVariable("hooks", [detectionProgress: detectionProgress, suggestBtn: suggestBtn, edgeBox: edgeBox, addBtn: addBtn, setForceFresh: { boolean v -> forceFresh = v }, setupCellpose: setupCellpose, findCellposePython: findCellposePython, setupBtn: setupBtn, setupRow: setupRow, cpLabel: cpLabel, showCellposeState: showCellposeState, setupText: SETUP_TEXT, saveMovie: saveMovie, makeMovie: makeMovie, rm: rm, refresh: refresh, go: go, status: status, steps: steps, outDir: outDir, showBg: showBg, frame: frame, suggest: suggest, acceptAt: acceptAt, acceptAll: acceptAll, hoverAt: hoverAt, proposals: proposals, result: result, info: info, removeSelected: removeSelected, undoLast: undoLast, root: root]); refresh(true); return }
frame.setVisible(true)
timer.start()
refresh(true)
javax.swing.SwingUtilities.invokeLater { checkUpdates(true) }

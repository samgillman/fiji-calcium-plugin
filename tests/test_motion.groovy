// The movement check: which frames are displaced from where the field sits for most of the recording, and by how much.
import ij.ImageStack
import ij.process.ShortProcessor

def loader = new GroovyClassLoader(getClass().classLoader); loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Motion = loader.loadClass('CalciumMotion'), View = loader.loadClass('CalciumFrameView')
int w = 400, h = 320, n = 30
def rnd = new Random(7), blobs = (1..60).collect { [rnd.nextInt(w - 60) + 30, rnd.nextInt(h - 60) + 30, 6 + rnd.nextInt(8), 200 + rnd.nextInt(600)] }
def frame = { int dx, int dy, int f ->
    short[] px = new short[w * h]; def noise = new Random(100 + f)
    for (int i = 0; i < px.length; i++) px[i] = (short) (300 + noise.nextInt(40) - f)            // camera noise on a field that fades a little
    blobs.each { b -> for (int y = -b[2]; y <= b[2]; y++) for (int x = -b[2]; x <= b[2]; x++) if (x * x + y * y <= b[2] * b[2]) { int xx = b[0] + x + dx, yy = b[1] + y + dy; if (xx >= 0 && xx < w && yy >= 0 && yy < h) px[yy * w + xx] += (short) b[3] } }
    new ShortProcessor(w, h, px, null) }
def stack = new ImageStack(w, h)
(1..n).each { f -> stack.addSlice('' + f, f <= 3 ? frame(8, -16, f) : f <= 6 ? frame(0, -8, f) : f == 20 ? frame(12, 0, f) : frame(0, 0, f)) }
int[][] shifts = Motion.shifts(stack, n)
assert shifts[0] as List == [8, -16] && shifts[2] as List == [8, -16] && shifts[3] as List == [0, -8] && shifts[6] as List == [0, 0] && shifts[19] as List == [12, 0] && shifts[29] as List == [0, 0] : shifts.collect { it as List }
def summary = Motion.summary(shifts)
assert summary.moved == 7 && summary.largest == 16 && summary.runs == [[1, 6], [20, 20]] && summary.frames == 30 : summary
// a recording that did not move, and one too small to check
def still = new ImageStack(w, h); (1..n).each { f -> still.addSlice('' + f, frame(0, 0, f)) }
assert Motion.summary(Motion.shifts(still, n)).moved == 0
def tiny = new ImageStack(40, 30); (1..5).each { tiny.addSlice('' + it, new ShortProcessor(40, 30)) }
assert Motion.summary(Motion.shifts(tiny, 5)) == [moved: 0, largest: 0, runs: [], frames: 5]
// the resting picture leaves moved frames out, unless nothing else is left
float[] all = View.mean(stack, 1, 10), settled = View.mean(stack, 1, 10, [1, 2, 3, 4, 5, 6] as Set), none = View.mean(stack, 1, 3, [1, 2, 3] as Set)
float[] seven = View.mean(stack, 7, 10)
assert (0..<settled.length).every { Math.abs(settled[it] - seven[it]) < 1e-3f } && (0..<all.length).any { Math.abs(all[it] - seven[it]) > 50f }
assert (0..<none.length).every { Math.abs(none[it] - View.mean(stack, 1, 3)[it]) < 1e-3f } : 'with every frame moved, the frames are used as they are'
println 'PASS: movement check (shifts per frame, runs of moved frames, still and tiny recordings) and a resting picture without moved frames'

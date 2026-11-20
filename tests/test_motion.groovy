// The movement check: how far each frame is displaced, counted from where the field sits for most of the recording or from the
// frames the cells are chosen on; which frames moved (8 px or more) and which sit slightly off (2.5 to 8 px); the sentence shown.
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
// frames 1-3 displaced by (8, -16), 4-6 by (0, -12), 7-9 by (0, -4), frame 20 by (12, 0); the rest in place
def placed = { int f -> f <= 3 ? [8, -16] : f <= 6 ? [0, -12] : f <= 9 ? [0, -4] : f == 20 ? [12, 0] : [0, 0] }
def stack = new ImageStack(w, h); (1..n).each { f -> stack.addSlice('' + f, frame(placed(f)[0], placed(f)[1], f)) }
double[][] shifts = Motion.shifts(stack, n)
(1..n).each { f -> assert Math.abs(shifts[f - 1][0] - placed(f)[0]) < 1d && Math.abs(shifts[f - 1][1] - placed(f)[1]) < 1d : 'frame ' + f + ': ' + (shifts[f - 1] as List) }
def summary = Motion.summary(shifts)
assert summary.moved == 7 && summary.runs == [[1, 6], [20, 20]] && summary.largest == 18 : summary
assert summary.slight == 3 && summary.slightRuns == [[7, 9]] && summary.slightLargest == 4 && summary.reference == null && summary.smeared == null : summary
assert Motion.sentence(summary) == 'The field moved (about 18 px) during frames 1\u20136 and 20. The field drifted by about 4 px during frames 7\u20139.' : Motion.sentence(summary)
// counted from the frames the cells are chosen on: steady frames give a reference; frames that moved themselves give "smeared"
def chosenLate = Motion.summary(shifts, 10, 19)
assert chosenLate.reference == [10, 19] && chosenLate.runs == [[1, 6], [20, 20]] && Motion.sentence(chosenLate).startsWith('The field moved (about 18 px from where the cells are chosen (frames 10\u201319)) during frames 1\u20136 and 20.')
def chosenEarly = Motion.summary(shifts, 7, 9)   // the cells were chosen while the field sat 4 px off: later frames are the ones that are off
assert chosenEarly.reference == [7, 9] && chosenEarly.slightRuns.containsAll([[10, 19], [21, 30]]) && chosenEarly.runs[0][0] == 1 && chosenEarly.runs.contains([20, 20]) : chosenEarly
def chosenMoving = Motion.summary(shifts, 1, 9)
assert chosenMoving.reference == null && chosenMoving.smeared == [1, 9] && Motion.sentence(chosenMoving).endsWith('The cells are chosen on frames 1\u20139, which include that movement, so their picture is smeared.') : Motion.sentence(chosenMoving)
// a slow drift that ends before the last frame reads "from frame N on" when counted from the first frames
double[][] drift = new double[100][2]; (0..<100).each { drift[it][1] = it < 40 ? 0d : Math.min(5d, (it - 40) * 0.25d) }
def late = Motion.summary(drift, 1, 20)
assert late.moved == 0 && late.slightRuns.size() == 1 && late.slightRuns[0][1] == 100 && Motion.sentence(late) == 'The field sits about 5 px from where the cells are chosen (frames 1\u201320) from frame ' + late.slightRuns[0][0] + ' on.' : Motion.sentence(late)
assert Motion.framesText([[2, 2], [4, 9], [12, 13], [20, 30]]) == '20 frames between 2 and 30' && Motion.framesText([[5, 5]]) == 'frame 5' && Motion.framesText([]) == ''
// a recording that did not move, and one too small to check
def still = new ImageStack(w, h); (1..n).each { f -> still.addSlice('' + f, frame(0, 0, f)) }
def quiet = Motion.summary(Motion.shifts(still, n)); assert quiet.moved == 0 && quiet.slight == 0 && Motion.sentence(quiet) == ''
def tiny = new ImageStack(40, 30); (1..5).each { tiny.addSlice('' + it, new ShortProcessor(40, 30)) }
assert Motion.summary(Motion.shifts(tiny, 5)).moved == 0 && Motion.sentence(null) == ''
// the resting picture leaves moved frames out, unless nothing else is left
float[] all = View.mean(stack, 1, 15), settled = View.mean(stack, 1, 15, (1..9) as Set), none = View.mean(stack, 1, 3, [1, 2, 3] as Set)
float[] rest = View.mean(stack, 10, 15)
assert (0..<settled.length).every { Math.abs(settled[it] - rest[it]) < 1e-3f } && (0..<all.length).any { Math.abs(all[it] - rest[it]) > 50f }
assert (0..<none.length).every { Math.abs(none[it] - View.mean(stack, 1, 3)[it]) < 1e-3f } : 'with every frame moved, the frames are used as they are'
println 'PASS: movement check (shifts between steps, moved and slightly-off runs, counted from the frames the cells are chosen on, smeared reference, sentences) and a resting picture without moved frames'

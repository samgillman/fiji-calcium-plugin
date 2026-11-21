// How a suggested outline is graded: likely a cell (1), unsure (2), very bright (3). Since 2026.11.20.20 a cell-sized outline that
// stands out and is not bright on average is likely a cell even with bright parts inside.
import ij.gui.Roi

def loader = new GroovyClassLoader(getClass().classLoader); loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Steps = loader.loadClass('CalciumSteps')
int w = 240, h = 240
// A square of the given side at (x, y): `base` gray levels, with the first `brightShare` of its pixels at `bright`; the surroundings at `around`.
def gradeOf = { int x, int y, int side, double base, double brightShare, double bright, double around = 0d ->
    float[] pic = new float[w * h]; java.util.Arrays.fill(pic, (float) around)
    int total = side * side, lit = (int) Math.round(brightShare * total), k = 0
    for (int yy = y; yy < y + side; yy++) for (int xx = x; xx < x + side; xx++) pic[yy * w + xx] = (float) (k++ < lit ? bright : base)
    def s = Steps.newInstance(); s.w = w; s.h = h; s.blind = true; s.blindPic = pic; s.flat = pic
    double[] ev = s.grade(new Roi(x, y, side, side))
    [grade: (int) ev[0], area: (int) ev[1], mean: ev[2], share: ev[3], contrast: ev[4]]
}
assert Steps.LIKELY_MIN_AREA == 500 && Steps.LIKELY_CONTRAST == 0.3d && Steps.LIKELY_BRIGHT_SHARE == 0.6d && Steps.BRIGHT_ABOVE == 8d
// as before: an even, dim cell is likely; one that hardly stands out is unsure
assert gradeOf(100, 100, 30, 3d, 0d, 0d).grade == 1
assert gradeOf(100, 100, 30, 0.1d, 0d, 0d).grade == 2 : 'contrast 0.1 is not enough for likely'
// a bright dot inside (5 % of the pixels) used to make a cell "unsure"
def dot = gradeOf(100, 100, 30, 3d, 0.05d, 20d)
assert dot.grade == 1 && dot.area == 900 && dot.mean < 8d && dot.share > 0.02d && dot.contrast > 0.3d : dot
assert gradeOf(100, 100, 20, 3d, 0.05d, 20d).grade == 2 : 'under 500 pixels the dot still makes it unsure'
def faint = gradeOf(100, 100, 30, 3d, 0.05d, 9d, 3.1d)   // mean 3.3 on surroundings of 3.1
assert faint.grade == 2 && faint.contrast > 0.15d && faint.contrast <= 0.3d : 'with the dot, standing out by 0.3 gray levels or less is not enough: ' + faint
// "very bright" only for its share of bright pixels (30 %), not bright on average: likely
def patchy = gradeOf(100, 100, 30, 3d, 0.30d, 12d)
assert patchy.grade == 1 && patchy.mean <= 8d && patchy.share >= 0.25d : patchy
// bright on average, or mostly bright pixels: stays very bright
assert gradeOf(100, 100, 30, 3d, 0.65d, 9.5d).grade == 3 : '65 % bright pixels'
assert gradeOf(100, 100, 30, 9d, 0.30d, 14d).grade == 3 : 'brighter than 8 gray levels on average'
assert gradeOf(100, 100, 30, 14d, 0d, 0d).grade == 3
// touching the picture's border: never accepted in bulk
assert gradeOf(0, 100, 30, 3d, 0.05d, 20d).grade == 2
// wrong size: not shown at all
assert gradeOf(100, 100, 10, 3d, 0d, 0d).grade == 4
println 'PASS: grading of suggested outlines: even dim cells, low contrast, a bright dot inside, patchy brightness, bright on average, the border, wrong size'

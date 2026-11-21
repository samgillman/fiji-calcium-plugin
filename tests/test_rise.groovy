// The response as the movie views show it (CalciumRise): a cell that rises by a gray level or two shows, a recording in which
// nothing responds stays dark, and neither the fading of the field nor a bright object that fades counts as a response.
import ij.ImageStack
import ij.process.ShortProcessor

def loader = new GroovyClassLoader(getClass().classLoader); loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Rise = loader.loadClass('CalciumRise'), View = loader.loadClass('CalciumFrameView')
int w = 320, h = 256, n = 60
def inBlob = { int x, int y, int cx, int cy, int r -> (x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r }
// camera noise of 4 gray levels on an unevenly lit field that fades by 5 gray levels; three blobs:
//   A at (80,80): a dim cell that rises by `rise` gray levels from frame 30;  B at (230,80): a bright object that fades by a third;
//   C at (160,190): a dim cell that does nothing
def movie = { double rise ->
    def stack = new ImageStack(w, h)
    (1..n).each { int f ->
        def rnd = new Random(1000 + f); short[] px = new short[w * h]
        for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) {
            double v = 100d + 0.05d * x - 5d * f / n
            if (inBlob(x, y, 80, 80, 13)) v += 6d + (f >= 30 ? rise : 0d)
            if (inBlob(x, y, 230, 80, 13)) v += 40d * (1d - 0.33d * f / n)
            if (inBlob(x, y, 160, 190, 13)) v += 8d
            px[y * w + x] = (short) Math.max(0, Math.round(v + 4d * rnd.nextGaussian()))
        }
        stack.addSlice(new ShortProcessor(w, h, px, null))
    }
    stack
}
def shareOn = { int[] shade, Closure where -> int on = 0, all = 0; for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) if (where(x, y)) { all++; if (shade[y * w + x] >= 64) on++ }; on / (double) Math.max(1, all) }
def outsideA = { int x, int y -> !inBlob(x, y, 80, 80, 20) }

def responding = movie(2.5d)
float[] rest = View.mean(responding, 1, 20)
def view = Rise.newInstance(responding, rest, 1, 20, 45)
assert view.noise > 0.1d && view.noise < 0.5d : 'averaging over 9 x 9 pixels and 5 frames leaves a fraction of a gray level of noise: ' + view.noise
assert Math.abs(view.start - 3d * view.noise) < 1e-9 && view.full >= view.start + 3d * view.noise - 1e-9
int[] shade = new int[w * h]
view.paint(responding, 10, shade)
assert shareOn(shade) { x, y -> true } < 0.005d : 'before the response the picture is dark: ' + shareOn(shade) { x, y -> true }
view.paint(responding, 45, shade)
assert shareOn(shade) { int x, int y -> inBlob(x, y, 80, 80, 9) } > 0.8d : 'the cell that rose by 2.5 gray levels shows: ' + shareOn(shade) { int x, int y -> inBlob(x, y, 80, 80, 9) }
assert shareOn(shade, outsideA) < 0.01d : 'nothing else does: ' + shareOn(shade, outsideA)
assert shareOn(shade) { int x, int y -> inBlob(x, y, 230, 80, 13) } < 0.02d : 'a bright object that fades is not a response'
assert shareOn(shade) { int x, int y -> inBlob(x, y, 160, 190, 13) } < 0.02d : 'a cell that does nothing stays dark'
view.paint(responding, 60, shade)
assert shareOn(shade, outsideA) < 0.01d : 'the fading of the field late in the recording is not a response: ' + shareOn(shade, outsideA)
assert shareOn(shade) { int x, int y -> inBlob(x, y, 80, 80, 9) } > 0.8d : 'the cell is still up at the end'

// a recording in which nothing responds: dark throughout, and its scale is the narrowest allowed
def quiet = movie(0d); float[] quietRest = View.mean(quiet, 1, 20); def quietView = Rise.newInstance(quiet, quietRest, 1, 20, 45)
assert Math.abs(quietView.full - (quietView.start + 3d * quietView.noise)) < 1e-6 : 'no response: full strength is 3 SDs above the start'
[10, 30, 45, 60].each { int f -> quietView.paint(quiet, f, shade); assert shareOn(shade) { x, y -> true } < 0.01d : 'frame ' + f + ' of a quiet recording: ' + shareOn(shade) { x, y -> true } }

// the box mean: a flat picture is unchanged and one bright pixel is spread over 81
float[] flat = new float[w * h]; java.util.Arrays.fill(flat, 7f); float[] work = new float[w * h]
Rise.box(flat, work, w, h, 4); assert flat.every { Math.abs(it - 7f) < 1e-4f }
float[] spike = new float[w * h]; spike[100 * w + 100] = 81f; Rise.box(spike, work, w, h, 4)
assert Math.abs(spike[100 * w + 100] - 1f) < 1e-4f && Math.abs(spike[104 * w + 104] - 1f) < 1e-4f && spike[105 * w + 100] == 0f
// the helpers for the other views
int[] one = [0x808080] as int[]; View.green(one); assert one[0] == 0x008000
int[] rgb = new int[1]; View.gray([200] as int[], rgb); assert rgb[0] == 0xC8C8C8
println 'PASS: the response as the views show it: a small rise shows, a quiet recording stays dark, fading and bright fading objects do not count; box mean; green and gray frames'

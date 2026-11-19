import ij.ImagePlus
import ij.gui.Roi
import ij.process.ColorProcessor

def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Steps = loader.loadClass('CalciumSteps')
def steps = Steps.newInstance(w:400,h:400)
assert steps.nearEdge(new Roi(39,100,20,20))
assert !steps.nearEdge(new Roi(40,100,20,20))
assert !steps.nearEdge(new Roi(340,100,20,20))
assert steps.nearEdge(new Roi(341,100,20,20))
steps.alsoAvoid = [new Roi(100,100,20,20)] as Roi[]
def occupied = steps.takenBy(new Roi[0])
assert occupied[110*400+110] && !occupied[200*400+200]
def a = new ImagePlus('a',new ColorProcessor(5,5))
def b = new ImagePlus('b',new ColorProcessor(5,5))
assert Steps.samePicture(a,b)
b.getProcessor().set(0,0x020202)
assert Steps.samePicture(a,b)
b.getProcessor().set(0,0x030000)
assert !Steps.samePicture(a,b)
assert !Steps.samePicture(a,new ImagePlus('size',new ColorProcessor(6,5)))
println 'PASS: edge threshold boundaries, excluded outlines stay out of background, cache image tolerance'

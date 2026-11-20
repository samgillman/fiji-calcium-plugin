// Sensitivity check after blind selection: which response-assisted outlines the blind detection also found, the report
// and the saved files, on hand-made outlines with known answers.
import ij.gui.*
import java.util.zip.ZipFile

def shell = new GroovyShell(getClass().classLoader, new Binding([CALCIUM_TEST: true]))
def Steps = shell.evaluate(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Sens = Steps.classLoader.loadClass('CalciumSensitivity')
// a suggestion as the window keeps it: a 20-pixel circle around (x, y)
def s = { int x, int y, int g, boolean edge = false -> [roi: new OvalRoi(x - 10, y - 10, 20, 20), grade: g, cx: x, cy: y, edge: edge, ev: [g, 314, 0, 0, 0, 0, 0, edge ? 1 : 0] as double[]] }
// response-assisted: two that respond (one found blind), one likely (found), one unsure (missed), one very bright (found), one at the edge
def assisted = [s(50, 50, 0), s(100, 50, 0), s(150, 50, 1), s(200, 50, 2), s(250, 50, 3), s(20, 20, 1, true)]
def blind = [s(52, 50, 1), s(150, 53, 2), s(250, 50, 3), s(400, 400, 1)]
def measured = [new OvalRoi(140, 40, 20, 20)] as Roi[]

def m = Sens.compare(assisted, blind, measured, true)
assert m.cellLike == 4 && m.found == 2 && m.bright == 1 && m.brightFound == 1 && m.edgeLeftOut == 1 : 'Counts over outlines that could be cells; edge and very bright apart'
assert m.byGrade[0] == [2, 1] && m.byGrade[1] == [1, 1] && m.byGrade[2] == [1, 0] && m.byGrade[3] == [1, 1]
assert m.inCells == 1 && m.blindShown == 4 && m.assistedShown == 6
def near = m.rows.find { it.x == 50 }; assert near.found && near.overlap > 0.7 && near.overlap < 1 : 'Two circles 2 px apart overlap strongly'
assert !m.rows.find { it.x == 100 }.found && m.rows.find { it.x == 100 }.overlap == 0
assert Sens.compare(assisted, blind, measured, false).cellLike == 5 : 'Edge outlines count when edge cells are kept'
assert Sens.compare(assisted, [], measured, true).found == 0
println 'PASS: sensitivity comparison counts found and missed by grade, overlap, edge and measured cells'

String text = Sens.report(m, 'Cellpose took 18 s (GPU)')
assert text.contains('found 2 of them (50 %) and missed 2') && text.contains('responds on its own 1 of 2') && text.contains('unsure 0 of 1')
assert text.contains('were not added to the cells and were not measured') && text.contains('How to read it')
def dir = new File('build/results/sensitivity'); dir.mkdirs(); String prefix = new File(dir, 'recording_').path
Sens.save(m, prefix)
def lines = new File(prefix + 'sensitivity_check.csv').readLines('UTF-8')
assert lines.size() == 6 && lines[0] == 'outline,grade_with_response,x,y,area_pixels,found_by_blind_detection,best_overlap_with_a_blind_outline,coincides_with_a_measured_cell'
assert lines.count { it.contains(',no,') && it.split(',')[5] == 'no' } == 2
def zip = new ZipFile(prefix + 'sensitivity_check_missed_RoiSet.zip')
try { assert zip.entries().toList()*.name.sort() == ['missed_002_responds_on_its_own.roi', 'missed_004_unsure.roi'] } finally { zip.close() }
// a missed very bright outline stays in the table but not in the missed set, which holds the ones counted as missed
def m2 = Sens.compare(assisted, blind.findAll { it.cx != 250 }, measured, true); String prefix2 = new File(dir, 'recording2_').path
Sens.save(m2, prefix2)
assert m2.bright == 1 && m2.brightFound == 0 && new File(prefix2 + 'sensitivity_check.csv').readLines('UTF-8').size() == 6
def zip2 = new ZipFile(prefix2 + 'sensitivity_check_missed_RoiSet.zip')
try { assert zip2.entries().toList()*.name.sort() == ['missed_002_responds_on_its_own.roi', 'missed_004_unsure.roi'] } finally { zip2.close() }
println 'PASS: sensitivity report wording and saved files (one row per outline, missed outlines that could be cells as an ROI set)'

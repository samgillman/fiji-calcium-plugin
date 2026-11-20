import ij.ImagePlus
import ij.process.ByteProcessor
import ij.gui.Roi
import ij.io.RoiEncoder
import java.util.zip.ZipOutputStream
import java.util.zip.ZipEntry
import javax.swing.SwingUtilities

def plugin = new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy')
def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(plugin)
def Task = loader.loadClass('CalciumOutlineTask')
def Panel = loader.loadClass('CalciumOutlineProgress')
def events = []
def task = Task.newInstance(progress: { p, s, d -> events << [p,s,d] })
task.report(10, 'Starting Cellpose')
task.logLine('INFO >>>> using CPU')
task.logLine('INFO Downloading: "model" to cache')
task.logLine(' 50%|#####| 500M/1G')
assert task.percent == 20
task.logLine('100%|##########| 1G/1G')
assert task.percent == 25 && !task.succeeded
task.logLine('INFO >>>> loading model cpsam')
task.logLine('  0%|          | 0/1 [00:00<?, ?it/s]')
assert task.stage == 'Finding cells in the image' && task.percent == 40
def explanation = task.detail
long stageAt = task.stageStarted
task.logLine('INFO routine technical diagnostic')
assert task.detail == explanation && task.stageStarted == stageAt
assert task.snapshot().history.contains('routine technical diagnostic')
task.logLine('INFO network run in 42.0s')
task.logLine('INFO masks created in 4.1s')
task.logLine('100%|##########| 1/1 [00:42<00:00]')
assert task.percent == 80
task.logLine('INFO >>>> completed in 46.1 sec')
assert task.percent == 84
assert events.collect{it[0]} == events.collect{it[0]}.sort(false)
println 'PASS: real log stages, CPU, model download, outer 100% never treated as completion'

def work = new File('build/test-output').absoluteFile; work.mkdirs()
def python = System.getenv('CALCIUM_TEST_PYTHON') ?: 'python3'
def output = new File(work, 'stream.log')
events.clear()
def streaming = Task.newInstance(progress: { p,s,d -> events << [p,s,d] })
def streamFixture = new File(work, 'stream-fixture.py')
streamFixture.text = 'import sys; print("Downloading: model"); sys.stderr.write("50%|#####| 5M/10M\\r100%|##########| 10M/10M\\r"); print("network run in 1s"); sys.exit(7)'
int exit = streaming.run([python, '-u', streamFixture.getPath()], work, output)
assert exit == 7 && !streaming.succeeded && events.any{it[0] == 20} && events.any{it[0] == 25} : output.text
assert output.text.contains('network run') && streaming.process == null
println 'PASS: subprocess stdout/stderr, carriage returns, exit failure and log file'

def cancelling = Task.newInstance(progress: { p,s,d -> })
def caught = new java.util.concurrent.atomic.AtomicReference()
def cancelFixture = new File(work, 'cancel-fixture.py')
cancelFixture.text = 'import time; print("started"); time.sleep(30)'
def worker = Thread.start {
    try { cancelling.run([python, '-u', cancelFixture.getPath()], work, new File(work,'cancel.log')) }
    catch (Throwable t) { caught.set(t) }
}
long deadline = System.nanoTime() + 5000000000L
while (cancelling.process == null && System.nanoTime() < deadline) Thread.sleep(10)
def child = cancelling.process; assert child != null
cancelling.cancel(); worker.join(5000)
assert !worker.alive && !child.alive && caught.get() instanceof java.util.concurrent.CancellationException
println 'PASS: cancelling a silent process stops it and exits promptly'

// Exercise the actual suggest closure with a mock Cellpose CLI and real ROI ZIP files.
def runDir = new File(work, 'recording with spaces/automatic_outlines').absoluteFile; runDir.mkdirs()
['fail-first','fail-all'].each { new File(runDir,it).delete() }
def module = new File(runDir, 'cellpose'); module.mkdirs(); new File(module,'__init__.py').text = ''
new File(module, '__main__.py').text = '''import sys, pathlib, shutil
p=pathlib.Path(sys.argv[sys.argv.index('--image_path')+1])
print('>>>> using CPU', flush=True)
if pathlib.Path('fail-all').exists() or ('--use_gpu' in sys.argv and pathlib.Path('fail-first').exists()):
    print('simulated failure', flush=True); sys.exit(7)
print('>>>> loading model cpsam', flush=True)
print('0%|          | 0/1', flush=True)
shutil.copyfile('fixture.zip', p.with_name(p.stem+'_rois.zip'))
print('>>>> completed in 0.1 sec', flush=True)
'''
def fixture = new File(runDir,'fixture.zip')
new ZipOutputStream(new FileOutputStream(fixture)).withCloseable { zip ->
    (1..3).each { n ->
        def bytes = new ByteArrayOutputStream(); new RoiEncoder(bytes).write(new Roi(n*10, 10, 15, 15))
        zip.putNextEntry(new ZipEntry('cell'+n+'.roi')); zip.write(bytes.toByteArray()); zip.closeEntry()
    }
}
def proposals = [[roi:'previous']]
System.setProperty('calcium.steps.pretend_no_gpu', 'true')
def stepStub = new Expando()
stepStub.outlinePicture = { -> new ImagePlus('mock', new ij.process.ColorProcessor(100,100)) }
def gradedWith = []
stepStub.grade = { roi, blindGrade = null -> gradedWith << blindGrade; [1,225,1,0,1,0,0,0] as double[] }
stepStub.say = { s -> }
stepStub.nearEdge = { roi -> roi.bounds.x < 20 }
def source = plugin.text
def closureText = source.substring(source.indexOf('def suggest = {'), source.indexOf('suggestBtn.addActionListener'))
def binding = new Binding([findCellposePython: { -> python }, cellposeBusy:false, testGui:true,
    outDir:runDir.parentFile, base:'recording with spaces', steps:stepStub, proposals:proposals, forceFresh:true, edgeBox:new javax.swing.JCheckBox('',true),
    GRADE_TEXT:['responds','likely','unsure','bright']])
def shell = new GroovyShell(loader,binding)
def suggest = shell.evaluate('import ij.IJ\nimport ij.Prefs\nimport ij.ImagePlus\nimport ij.gui.GenericDialog\n'+closureText+'\nreturn suggest')
events.clear()
def success = Task.newInstance(progress: { p,s,d -> events << [p,s,d] })
def message = suggest(success)
assert message.startsWith('3 suggestions;') && success.succeeded && success.percent == 98 && proposals.size() == 3
assert events.findAll{it[1]=='Checking the suggested cells'}.collect{it[2]} == ['Checked 1 of 3 outlines.','Checked 2 of 3 outlines.','Checked 3 of 3 outlines.']
assert new File(runDir,'recording with spaces_suggestions_with_evidence.csv').readLines().size() == 4
println 'PASS: actual suggest closure imports/grades ROIs, writes CSV, handles paths with spaces'
assert proposals.count { it.edge } == 1
assert new File(runDir,'cellpose_options_used.txt').text.contains('--diameter 60')
assert success.snapshot().context.contains('New detection')
assert success.snapshot().context.contains('100 x 100 pixels')
assert success.snapshot().context.contains('diameter 60')
assert !new File(runDir,'cellpose_run_incomplete.txt').exists()
// The sensitivity check's comparison run: saved outlines of the same picture are reused even after "Suggest cells again",
// the outlines are graded with the response and go to their own list, and the suggestions on screen are untouched.
def comparisonOutlines = [], shownBefore = proposals.collect { it }; gradedWith.clear()
def comparison = Task.newInstance(progress: { p,s,d -> })
def compared = suggest(comparison, [picture: new ImagePlus('response-assisted', new ij.process.ColorProcessor(100,100)), dir: runDir, into: comparisonOutlines])
assert compared.startsWith('3 response-assisted outlines;') && compared.contains('were used again') && comparison.succeeded
assert comparisonOutlines.size() == 3 && proposals == shownBefore && gradedWith == [false, false, false] && !comparisonOutlines.any { it.blind }
println 'PASS: comparison run for the sensitivity check reuses saved outlines, grades with the response, leaves suggestions alone'
new File(runDir,'fail-all').text = ''
binding.setVariable('forceFresh',false)
def reused = Task.newInstance(progress: { p,s,d -> })
assert suggest(reused).contains('were used again') && reused.succeeded
assert reused.snapshot().context.contains('Saved outlines reused (no new Cellpose run)')
assert !reused.snapshot().context.contains('Attempt 1')
new File(runDir,'cellpose_run_incomplete.txt').text = 'incomplete'
def unfinished = Task.newInstance(progress: { p,s,d -> })
assert suggest(unfinished).contains('exit code 7') && !unfinished.succeeded
assert new File(runDir,'cellpose_run_incomplete.txt').exists()
new File(runDir,'fail-all').delete()
binding.setVariable('forceFresh',true)
println 'PASS: preserved edge flags and CPU settings, matched-image cache, interrupted-run invalidation'


new File(runDir,'fail-first').text = ''
def retry = Task.newInstance(progress: { p,s,d -> })
assert suggest(retry).startsWith('3 suggestions;') && retry.succeeded && retry.logFile.name == 'cellpose_output_second_try.txt'
assert retry.snapshot().context.contains('Attempt 2 of 2; default model, full-size detection')
new File(runDir,'fail-all').text = ''
def failure = Task.newInstance(progress: { p,s,d -> })
assert suggest(failure).contains('exit code 7') && !failure.succeeded && proposals.size() == 3
assert !new File(runDir,'recording with spaces_for_automatic_outlines_rois.zip').exists()
println 'PASS: fallback retry; failures retain old accepted suggestions without reusing stale ZIP'
['fail-first','fail-all'].each { new File(runDir,it).delete() }

// The merge rule for blind detection: every outline from the 0-2 picture is kept; an outline from the 0-8 picture is added
// only when it shares less than 10 % of the smaller outline's pixels with every outline from the 0-2 picture.
def Steps = loader.loadClass('CalciumSteps')
assert Steps.BLIND_CAP == 2d && Steps.BLIND_CAP_HIGH == 8d && Steps.BLIND_MERGE_OVERLAP == 0.10d
def sq = { x, y, s -> new Roi(x, y, s, s) }
assert Steps.blindOutlinesToAdd([sq(0,0,10)], [sq(9,0,10)], 100) == [] : 'exactly 10 % shared: the same cell'
assert Steps.blindOutlinesToAdd([sq(0,0,10)], [sq(9,1,10)], 100) == [0] : '9 % shared: a new cell'
assert Steps.blindOutlinesToAdd([sq(0,0,30)], [sq(5,5,3)], 100) == [] : 'a small outline inside a big one'
assert Steps.blindOutlinesToAdd([sq(50,50,3)], [sq(40,40,30)], 100) == [] : 'a big outline over a small one (share of the smaller)'
assert Steps.blindOutlinesToAdd([sq(0,0,10), sq(40,0,10)], [sq(60,60,10), sq(1,1,10), sq(20,0,10)], 100) == [0, 2]
assert Steps.blindOutlinesToAdd([], [sq(0,0,10)], 100) == [0] && Steps.blindOutlinesToAdd([sq(0,0,10)], [], 100) == []
println 'PASS: blind merge rule (10 % of the smaller outline, kept outlines never removed)'

// Blind detection runs Cellpose on both pictures, each in its own folder, combines the outlines and splits the bar.
def blindDir = new File(runDir.parentFile, 'automatic_outlines_blind')
blindDir.deleteDir()
def picDirs = ['picture_1_shown_0_to_2_gray_levels', 'picture_2_shown_0_to_8_gray_levels'].collect { new File(blindDir, it) }
def writeFixture = { File f, List rois -> new ZipOutputStream(new FileOutputStream(f)).withCloseable { zip ->
    rois.eachWithIndex { r, n -> def bytes = new ByteArrayOutputStream(); new RoiEncoder(bytes).write(r); zip.putNextEntry(new ZipEntry('cell' + (n + 1) + '.roi')); zip.write(bytes.toByteArray()); zip.closeEntry() } } }
picDirs.each { d -> def m = new File(d, 'cellpose'); m.mkdirs(); new File(m, '__init__.py').text = ''; new File(m, '__main__.py').text = new File(module, '__main__.py').text }
writeFixture(new File(picDirs[0], 'fixture.zip'), (1..3).collect { n -> new Roi(n*10, 10, 15, 15) })
// picture 2: one on cell 1 (shares most of it), one touching cells 1-3 slightly (under 10 %), one on its own
writeFixture(new File(picDirs[1], 'fixture.zip'), [new Roi(11, 11, 15, 15), new Roi(24, 24, 15, 15), new Roi(60, 60, 15, 15)])
def said = []
def blindStub = new Expando()
blindStub.blind = true; blindStub.blindFirst = 1; blindStub.blindLast = 113
def capPicture = { double cap -> def cp = new ij.process.ColorProcessor(100, 100); cp.setColor(new java.awt.Color(0, (int) (cap * 20), 0)); cp.fill(); new ImagePlus('blind', cp) }
blindStub.blindOutlinePicture = { double cap -> capPicture(cap) }
blindStub.outlinePicture = { -> capPicture(Steps.BLIND_CAP) }
blindStub.blindBasis = { -> 'Pre-stimulus, blind (frames 1-113)' }
blindStub.blindRecipe = { -> 'shown to it twice, from 0 to 2 and from 0 to 8 gray levels' }
def blindGraded = []
blindStub.grade = { roi, blindGrade = null -> blindGraded << blindGrade; [1,225,1,0,1,0,0,0] as double[] }
blindStub.say = { s -> said << s }
blindStub.nearEdge = { roi -> roi.bounds.x < 20 }
binding.setVariable('steps', blindStub); binding.setVariable('forceFresh', true)
events.clear()
def blindTask = Task.newInstance(progress: { p,s,d -> events << [p,s,d] })
def blindMessage = suggest(blindTask)
assert blindMessage.startsWith('5 suggestions from two pictures; Cellpose took ') && blindTask.succeeded && proposals.size() == 5 : blindMessage
assert blindMessage ==~ /5 suggestions from two pictures; Cellpose took \d+ s and \d+ s \(CPU\)/ : 'a short line for the window: ' + blindMessage
assert blindTask.snapshot().context.contains('Two pictures of the frames before the stimulus (0-2 and 0-8 gray levels)')
assert proposals.collect { it.picture } == ['0 to 2', '0 to 2', '0 to 2', '0 to 8', '0 to 8'] && proposals.every { it.blind }
assert proposals.collect { [it.roi.bounds.x, it.roi.bounds.y] } == [[10,10], [20,10], [30,10], [24,24], [60,60]]
assert blindGraded == [null] * 5 : 'blind suggestions are graded the blind way'
assert events.collect { it[0] } == events.collect { it[0] }.sort(false) : 'the bar never moves back'
assert events.findAll { it[1] == 'Finishing Cellpose' }.collect { it[0] }.unique() == [42, 84] : 'picture 1 ends at 42 %, picture 2 at 84 %'
assert events.findAll { it[1] == 'Finding cells in the image' }.collect { it[0] }.unique() == [20, 62]
assert blindTask.percent == 98 && blindTask.partFrom == 0 && blindTask.partShare == 1
picDirs.each { d ->
    assert new File(d, 'recording with spaces_for_automatic_outlines.tif').isFile() && new File(d, 'recording with spaces_for_automatic_outlines_rois.zip').isFile()
    assert new File(d, 'cellpose_options_used.txt').isFile() && !new File(d, 'cellpose_run_incomplete.txt').exists()
}
def evidence = new File(blindDir, 'recording with spaces_suggestions_with_evidence.csv').readLines()
assert evidence.size() == 6 && evidence[0].endsWith(',found_on_picture_gray_levels') && evidence.drop(1).collect { it.split(',')[-1] } == ['0 to 2'] * 3 + ['0 to 8'] * 2
def combined = new java.util.zip.ZipFile(new File(blindDir, 'recording with spaces_suggestions_rois.zip'))
assert combined.entries().toList()*.name == (1..5).collect { String.format('suggestion_%03d.roi', it) }; combined.close()
def note = new File(blindDir, 'detection_picture.txt').text
assert note.contains('Pre-stimulus, blind (frames 1-113)') && note.contains('picture_1_shown_0_to_2_gray_levels, picture_2_shown_0_to_8_gray_levels')
assert note.contains('Combined: all 3 outlines from picture 1 (0-2 gray levels) kept; 2 of 3 from picture 2 (0-8 gray levels) added')
assert said.any { it.startsWith('Automatic outlines are blind: Cellpose got two pictures of frames 1-113 only') && it.contains('2 of 3 from picture 2 (0-8 gray levels) added') }
assert said.any { it.contains('picture 1 (0-2 gray levels): ' + python) && it.contains('picture 2 (0-8 gray levels): ' + python) && it.contains('5 outlines shown') && it.contains('picture 2 (0-8 gray levels): Cellpose took') } : 'the record keeps both commands and both timings'
println 'PASS: blind detection runs Cellpose on both pictures, combines the outlines, records both, splits the bar'

// Both pictures are reused when nothing changed; a failure on picture 2 keeps the earlier suggestions and is named.
binding.setVariable('forceFresh', false)
picDirs.each { new File(it, 'fail-all').text = '' }
def blindReuse = Task.newInstance(progress: { p,s,d -> })
def reusedMessage = suggest(blindReuse)
assert reusedMessage == '5 suggestions from two pictures; the outlines made earlier for this recording were used again' && blindReuse.succeeded : reusedMessage
assert blindReuse.snapshot().context.contains('Two pictures of the frames before the stimulus') && blindReuse.snapshot().context.contains('Saved outlines reused')
new File(picDirs[0], 'fail-all').delete()
binding.setVariable('forceFresh', true)
def shownBeforeFailure = proposals.collect { it }
def blindFailure = Task.newInstance(progress: { p,s,d -> })
def failedMessage = suggest(blindFailure)
assert failedMessage.contains('on blind picture 2 (0-8 gray levels)') && failedMessage.contains('exit code 7') && !blindFailure.succeeded
assert proposals == shownBeforeFailure && blindFailure.partFrom == 0 && blindFailure.partShare == 1
assert new File(picDirs[0], 'cellpose_run_incomplete.txt').exists() && new File(picDirs[1], 'cellpose_run_incomplete.txt').exists() : 'neither half is reused after a failed pair'
new File(picDirs[1], 'fail-all').delete()
binding.setVariable('steps', stepStub)
println 'PASS: blind reuse of both pictures; a failure on picture 2 is named and keeps the earlier suggestions'

def panel
SwingUtilities.invokeAndWait {
    panel = Panel.newInstance(); def pt = panel.begin()
    pt.report(40,'Finding cells','CPU inference')
    panel.started = System.nanoTime() - 601000000000L
    panel.latest.set([percent:40,stage:'Finding cells',detail:'CPU inference',at:System.nanoTime()-61000000000L])
    panel.render()
    assert panel.bar.value == 40 && panel.text.text.contains('10 minutes') && panel.text.text.contains('No new progress')
    assert panel.text.text.contains('CPU inference') // Quiet feedback must not replace the explanation.
    assert panel.detailsText.text.contains('Percentages mark completed stages')
    panel.finish(false,'Failed <example>'); assert panel.bar.value < 100 && !panel.clock.running && panel.text.text.contains('&lt;example&gt;')
    def next = panel.begin(); assert panel.bar.value == 0
    assert !panel.detailsText.text.contains('Failed <example>')
    next.report(98,'Saving',''); panel.finish(true,'3 suggestions')
    assert panel.bar.value == 100 && !panel.clock.running && !panel.cancelButton.enabled
    def again = panel.begin(); panel.close(); assert again.cancelled && !panel.clock.running
}
println 'PASS: panel warnings, safe error text, 100% only on success, retry reset and close cleanup'
def bounded = Task.newInstance(progress: { p,s,d -> })
bounded.report(40, 'Finding cells', 'Persistent explanation')
(1..150).each { bounded.logLine('diagnostic ' + it) }
assert bounded.snapshot().history.readLines().size() == 120
assert bounded.detail == 'Persistent explanation' && bounded.percent == 40
assert !bounded.snapshot().history.contains('diagnostic 1\n')
println 'PASS: persistent explanations, stage timers, detection context, bounded event history'
System.clearProperty('calcium.steps.pretend_no_gpu')
println 'All progress checks passed.'

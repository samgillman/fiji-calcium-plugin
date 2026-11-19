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
def stepStub = new Expando()
stepStub.outlinePicture = { -> new ImagePlus('mock', new ByteProcessor(100,100)) }
stepStub.grade = { roi -> [1,225,1,0,1,0,0,0] as double[] }
stepStub.say = { s -> }
def source = plugin.text
def closureText = source.substring(source.indexOf('def suggest = {'), source.indexOf('suggestBtn.addActionListener'))
def binding = new Binding([findCellposePython: { -> python }, cellposeBusy:false, testGui:true,
    outDir:runDir.parentFile, base:'recording with spaces', steps:stepStub, proposals:proposals,
    GRADE_TEXT:['responds','likely','unsure','bright']])
def shell = new GroovyShell(loader,binding)
def suggest = shell.evaluate('import ij.IJ\nimport ij.Prefs\nimport ij.gui.GenericDialog\n'+closureText+'\nreturn suggest')
events.clear()
def success = Task.newInstance(progress: { p,s,d -> events << [p,s,d] })
def message = suggest(success)
assert message == '3 suggestions' && success.succeeded && success.percent == 98 && proposals.size() == 3
assert events.findAll{it[1]=='Checking the suggested cells'}.collect{it[2]} == ['Checked 1 of 3 outlines.','Checked 2 of 3 outlines.','Checked 3 of 3 outlines.']
assert new File(runDir,'recording with spaces_suggestions_with_evidence.csv').readLines().size() == 4
println 'PASS: actual suggest closure imports/grades ROIs, writes CSV, handles paths with spaces'

new File(runDir,'fail-first').text = ''
def retry = Task.newInstance(progress: { p,s,d -> })
assert suggest(retry) == '3 suggestions' && retry.succeeded && retry.logFile.name == 'cellpose_output_second_try.txt'
new File(runDir,'fail-all').text = ''
def failure = Task.newInstance(progress: { p,s,d -> })
assert suggest(failure).contains('exit code 7') && !failure.succeeded && proposals.size() == 3
assert !new File(runDir,'recording with spaces_for_automatic_outlines_rois.zip').exists()
println 'PASS: fallback retry; failures retain old accepted suggestions without reusing stale ZIP'

def panel
SwingUtilities.invokeAndWait {
    panel = Panel.newInstance(); def pt = panel.begin()
    pt.report(40,'Finding cells','CPU inference')
    panel.started = System.nanoTime() - 601000000000L
    panel.latest.set([percent:40,stage:'Finding cells',detail:'CPU inference',at:System.nanoTime()-61000000000L])
    panel.render()
    assert panel.bar.value == 40 && panel.text.text.contains('10 minutes') && panel.text.text.contains('No new progress')
    panel.finish(false,'Failed <example>'); assert panel.bar.value < 100 && !panel.clock.running && panel.text.text.contains('&lt;example&gt;')
    def next = panel.begin(); assert panel.bar.value == 0
    next.report(98,'Saving',''); panel.finish(true,'3 suggestions')
    assert panel.bar.value == 100 && !panel.clock.running && !panel.cancelButton.enabled
    def again = panel.begin(); panel.close(); assert again.cancelled && !panel.clock.running
}
println 'PASS: panel warnings, safe error text, 100% only on success, retry reset and close cleanup'
println 'All progress checks passed.'

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Runtime=loader.loadClass('CalciumRRuntime'),Job=loader.loadClass('CalciumLocalJob'),Engine=loader.loadClass('CalciumSimpleCaEngine'),Csv=loader.loadClass('CalciumCsv')
assert Runtime.platform('Mac OS X','aarch64')=='osx-arm64'
assert Runtime.platform('Windows 11','amd64')=='win-64'
try {Runtime.platform('Windows 11','aarch64');assert false}catch(IOException expected){}
def root=new File('build/analysis-bridge').absoluteFile;root.mkdirs()
def engine=Engine.unpack(root)
assert new File(engine,'upstream/LICENSE').isFile()
assert new File(engine,'provenance.tsv').bytes==new File('analysis/provenance.tsv').bytes
assert new File(engine,'upstream/response_protocols.R').isFile()
def csv=new File(root,'quoted.csv');csv.text='"Cell","Reason"\n"Cell1","line 1,\nline ""2"""\n'
assert Csv.read(csv)==[['Cell','Reason'],['Cell1','line 1,\nline "2"']]
def python=System.getenv('CALCIUM_TEST_PYTHON')?:'python3'
def sleeper=new File(root,'sleep.py');sleeper.text='import time\ntime.sleep(30)\n'
def job=Job.newInstance(logFile:new File(root,'cancel.log'))
def failed=new java.util.concurrent.atomic.AtomicReference()
def thread=Thread.start{try{job.run([python,sleeper.path],root)}catch(Throwable e){failed.set(e)}}
long deadline=System.nanoTime()+5000000000L
while(job.process==null&&System.nanoTime()<deadline)Thread.sleep(10)
assert job.process!=null
job.cancel();thread.join(5000)
assert !thread.alive && failed.get() instanceof java.util.concurrent.CancellationException
job=Job.newInstance(logFile:new File(root,'fail.log'))
try {job.run([python,'-c','raise SystemExit(7)'],root);assert false}catch(IOException expected){assert expected.message.contains('7')}
// A stimulus and a positive control: the engine's two cross-window rules are checked before a run, with what to change.
def CrossProtocol=loader.loadClass('CalciumResponseProtocol')
def stimulusTo={double end->CrossProtocol.window('Stimulus','ATP',93.12d,112.09d,end)}
def control={double b0,double onset,double end->CrossProtocol.window('Control','KCl',b0,onset,end)}
def refusal={List w->try{CrossProtocol.payload(w,2,'Stress','provisional');null}catch(IllegalArgumentException e){e.message}}
assert refusal([stimulusTo(299.24d),control(187.15d,299.24d,300.24d)]).contains('KCl baseline (187.15–299.24 s) falls inside the ATP response window')
assert refusal([stimulusTo(250d),control(200d,250d,299d)]).contains('End the ATP window at or before 200 s')
assert refusal([stimulusTo(260d),control(255d,256d,299d)]).contains('response windows overlap')
assert refusal([stimulusTo(200d),control(200d,250d,299d)])==null : 'The stress test run that worked still works'
assert refusal([stimulusTo(299.24d)])==null
// Edge cases of the same two rules: a control placed before the stimulus, a baseline that only partly overlaps, and touching windows.
assert refusal([stimulusTo(200d),control(40d,50d,100d)]).contains('ATP baseline (93.12–112.09 s) falls inside the KCl response window (50–100 s)')
assert refusal([stimulusTo(200d),control(195d,210d,260d)]).contains('KCl baseline (195–210 s) falls inside the ATP response window')
assert refusal([stimulusTo(200d),control(40d,50d,93.12d)])==null : 'a control window that ends where the stimulus baseline starts is allowed'
// When the engine does stop, its own message is shown instead of a bare exit code.
def failing=new File(root,'engine-error.py');failing.text="import sys\nprint('SIMPLECA\\t70\\tValidating')\nprint('Error in responder_validate_protocol(protocol) : ')\nprint(\"  Each stimulus requires a baseline outside the other stimulus's response interval.\")\nprint('Calls: responder_import_protocol -> responder_validate_protocol')\nsys.exit(1)\n"
job=Job.newInstance(logFile:new File(root,'engine-error.log'))
try {job.run([python,failing.path],root);assert false}catch(IOException expected){assert expected.message=="The analysis engine stopped: Each stimulus requires a baseline outside the other stimulus's response interval." : expected.message}
def plain=new File(root,'plain-error.log');plain.text='Error: cannot open file\n';assert Job.engineError(plain)=='cannot open file'
plain.text='nothing useful\n';assert Job.engineError(plain)==null
// Windows line endings reach the log as a blank line after every line; the message is still found.
plain.text="Error in responder_validate_protocol(protocol) : \n\n  Each stimulus requires a baseline\n\n  outside the other window.\n\nCalls: a -> b\n\n"
assert Job.engineError(plain)=='Each stimulus requires a baseline outside the other window.'
// A run in which no cell could be scored is reported as such, with the commonest reason, not as a green "Complete".
def Settings=loader.loadClass('CalciumAnalysisSettings')
def empty=new File(root,'no-evaluable');empty.mkdirs()
new File(empty,'responder_counts.csv').text='"Field","Stimulus","Role","Registered","Responded","Evaluable","Unscorable"\n"F","ATP","Stimulus",3,0,0,3\n'
new File(empty,'responder_calls.csv').text='"Field","Cell","Role","Status","Reason"\n"F","A","Stimulus","Unscorable","Fewer baseline samples than required"\n"F","B","Stimulus","Unscorable","Fewer baseline samples than required"\n"F","C","Stimulus","Unscorable","Sample gap exceeds the field limit"\n'
assert Settings.noEvaluableCells(empty)=='Finished, but no cell could be scored for ATP: Fewer baseline samples than required (2 cells). Check the settings; the results were saved for the record.'
assert Settings.noEvaluableCells(new File('tests/fixtures/simpleca'))==null : 'A normal run is still Complete'
println 'PASS: empty-run warning, stimulus/control window rules before a run, readable engine errors, platform selection, embedded source provenance, CSV quoting, cancellation and process failure'

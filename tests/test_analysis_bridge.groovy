def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Runtime=loader.loadClass('CalciumRRuntime'),Job=loader.loadClass('CalciumLocalJob'),Engine=loader.loadClass('CalciumSimpleCaEngine'),Csv=loader.loadClass('CalciumCsv')
assert Runtime.platform('Mac OS X','aarch64')=='osx-arm64'
assert Runtime.platform('Windows 11','amd64')=='win-64'
try {Runtime.platform('Windows 11','aarch64');assert false}catch(IOException expected){}
def root=new File('build/analysis-bridge').absoluteFile;root.mkdirs()
def engine=Engine.unpack(root)
assert new File(engine,'upstream/LICENSE').isFile()
assert new File(engine,'provenance.tsv').text.contains('a472f2460075343d5b97c8bac32ddf3da0c77386')
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
println 'PASS: platform selection, embedded source provenance, CSV quoting, cancellation and process failure'

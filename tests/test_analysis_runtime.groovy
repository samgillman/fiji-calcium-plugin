// Install exactly the same runtime users receive, then exercise the complete Java -> R -> CSV path.
import ij.gui.Roi

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Runtime=loader.loadClass('CalciumRRuntime'),Job=loader.loadClass('CalciumLocalJob'),Engine=loader.loadClass('CalciumSimpleCaEngine'),Csv=loader.loadClass('CalciumCsv'),Protocol=loader.loadClass('CalciumResponseProtocol'),Installer=loader.loadClass('CalciumInstaller')
def root=new File('build/runtime with spaces').absoluteFile
root.mkdirs()
def runtime=Runtime.newInstance(home:root)
def job=Job.newInstance(logFile:new File(root,'setup-test.log'),progress:{p,s->println "${p}% ${s}"})
try {
def env=runtime.setup(job)
assert runtime.environment()==env
runtime.verifyExecutable()
// No download or package solve on an already complete environment.
assert runtime.setup(job)==env
File repository=new File('.').canonicalFile
job.run(runtime.scriptCommand(env,new File(repository,'tests/test_analysis.R')),repository,180)
def settings=[:]
new File('build/analysis-fixtures/settings.tsv').eachLine('UTF-8'){line->def pair=line.split('\t',2);settings[pair[0]]=pair[1]}
File result=Engine.analyze(runtime,new File('build/analysis-fixtures/input.csv').absoluteFile,settings,new File('build/analysis output').absoluteFile,job)
assert new File(result,'COMPLETE').isFile()
def metrics=Csv.read(new File(result,'metrics.csv'))
assert metrics.size()==7 && metrics[1][0]=='Cell001' && metrics[3][3]=='Excluded'
assert new File(result,'input_corrected_fluorescence.csv').bytes==new File('build/analysis-fixtures/input.csv').bytes
assert Csv.read(new File(result,'analysis_manifest.csv')).any{it[0]=='plugin_version'}
// The managed Windows/Mac/Linux path also round-trips the portable protocol
// through the native editor contract and checks the independent shared fixture.
def portable=Protocol.read(new File('tests/fixtures/response-protocols/golden-protocol.json'))
def payload=Protocol.payload(portable.windows,2,portable.name,portable.validation_status,null,portable.sd_estimator)
assert payload.sd_estimator=='sample_sd' && payload.engine_version=='2.0.0' : 'a protocol that was read keeps its estimator'
def shared=Csv.read(new File('build/analysis-fixtures/shared-golden-input.csv'))
def nativeInput=new File(root,'native measured.csv'),timestamps=new File(root,'per frame times.csv')
Csv.write(nativeInput,shared.collect{it.take(1)+it.drop(2)})
Csv.write(timestamps,[['Frame','Time_s']]+shared.drop(1).collect{it.take(2)})
def nativeSettings=[baseline_start:'1',baseline_end:'4',sampling_hz:'',recording:'Golden field',input_mode:'background_corrected_fluorescence',
    protocol_json:groovy.json.JsonOutput.toJson(payload),expected_input_sha256:Installer.digest(nativeInput.bytes),
    timestamps_file:timestamps.path,timestamps_sha256:Installer.digest(timestamps.bytes)]
def nativeResult=Engine.analyze(runtime,nativeInput,nativeSettings,root,job)
def nativeCounts=Csv.read(new File(nativeResult,'responder_counts.csv'))
int evaluable=nativeCounts[0].indexOf('Evaluable'),responded=nativeCounts[0].indexOf('Responded')
assert nativeCounts.drop(1).collect{it[evaluable]}==['4','5','3']
assert nativeCounts.drop(1).collect{it[responded]}==['1','4','1']
assert new File(nativeResult,'timestamps.csv').bytes==timestamps.bytes
assert Protocol.read(new File(nativeResult,'protocol.json')).windows[1].Stimulus=='Ionomycin'
assert new File(nativeResult,'measured_fluorescence.csv').bytes==nativeInput.bytes
try {
    Engine.analyze(runtime,new File('build/analysis-fixtures/input.csv').absoluteFile,settings+[expected_input_sha256:'outdated-measurement'],new File('build/analysis output').absoluteFile,job)
    assert false : 'Changed input must not be linked to old ROI snapshots'
} catch(IOException expected) {assert expected.message.contains('measurement file changed')}
// Bundled engine cache can repair a corrupted file, without fetching or running remote source.
def engine=Engine.unpack(root)
new File(engine,'runner.R').text='corrupted'
assert Engine.unpack(root)==engine
assert new File(engine,'runner.R').getText('UTF-8').startsWith('# Fiji adapter.')
new GroovyShell(loader,new Binding([loader:loader,runtime:runtime,root:root])).evaluate(new File('tests/test_web_analysis_runtime.groovy'))
println 'PASS: managed R/jsonlite setup and cache, legacy and independent golden analysis, native protocol round trip, timestamp snapshots, paths with spaces and embedded source repair'
} catch(Throwable problem) {
    if(job.logFile.isFile())println(job.logFile.getText('UTF-8'))
    throw problem
}

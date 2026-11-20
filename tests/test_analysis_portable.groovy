// Exercise the real Java -> R adapter with an explicitly supplied local Rscript.
// Managed installation is covered separately by test_analysis_runtime on each CI OS.
def rscript=System.getenv('CALCIUM_TEST_RSCRIPT')
assert rscript && new File(rscript).isFile() : 'Set CALCIUM_TEST_RSCRIPT to an Rscript with jsonlite installed'
def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Runtime=loader.loadClass('CalciumRRuntime'),Job=loader.loadClass('CalciumLocalJob'),Engine=loader.loadClass('CalciumSimpleCaEngine'),Csv=loader.loadClass('CalciumCsv'),Protocol=loader.loadClass('CalciumResponseProtocol'),Measurement=loader.loadClass('CalciumMeasurement'),Installer=loader.loadClass('CalciumInstaller')
def root=new File('build/portable bridge with spaces').absoluteFile;root.mkdirs()
def runtime=Runtime.newInstance(home:root)
runtime.metaClass.environment={->root}
runtime.metaClass.verifyExecutable={->}
runtime.metaClass.scriptCommand={File env,File script,List arguments=[]->[rscript,'--vanilla',script.path]+arguments}
def job=Job.newInstance(logFile:new File(root,'portable-test.log'))
job.run([rscript,'--vanilla',new File('tests/test_analysis.R').absolutePath],new File('.').canonicalFile,180)
def fixture=new File('build/analysis-fixtures').absoluteFile
def imported=Protocol.read(new File(fixture,'golden-out/protocol.json'))
def nativePayload=Protocol.payload(imported.windows,2,imported.name,imported.validation_status)
def input=new File(fixture,'golden-input.csv')
def settings=[baseline_start:'1',baseline_end:'5',sampling_hz:'',recording:'Golden',input_mode:'background_corrected_fluorescence',condition:'Vehicle',biological_unit:'Culture1',
    protocol_json:groovy.json.JsonOutput.toJson(nativePayload),roi_registry_file:new File(fixture,'golden-out/roi_registry.csv').path,
    sensitivity_file:new File(fixture,'golden-diagnostics/sensitivity_grid.csv').path,negative_control_file:new File(fixture,'golden-diagnostics/negative_control_windows.csv').path]
settings.expected_input_sha256=Installer.digest(input.bytes)
def result=Engine.analyze(runtime,input,settings,root,job)
assert new File(result,'COMPLETE').isFile()
def counts=Csv.read(new File(result,'responder_counts.csv'))
int yes=counts[0].indexOf('Responded'),evaluable=counts[0].indexOf('Evaluable')
assert counts.drop(1).collect{it[yes]}==['1','3','1']
assert counts.drop(1).collect{it[evaluable]}==['3','4','2']
assert new File(result,'sensitivity_grid.csv').isFile() && new File(result,'negative_control_counts.csv').isFile()
assert new File(result,'input_corrected_fluorescence.csv').bytes==input.bytes
def manifest=Csv.read(new File(result,'analysis_manifest.csv')).drop(1).collectEntries{[it[0],it[1]]}
assert manifest.time_source=='explicit_timestamps' && manifest.protocol_schema_version=='2'
assert Protocol.read(new File(result,'protocol.json')).windows[1].Stimulus=='Ionomycin'

// Timestamp sidecars are snapshotted and aligned without rounding the original fluorescence.
def rows=Csv.read(input),noTime=new File(root,'measured.csv'),times=new File(root,'times.csv')
Csv.write(noTime,rows.collect{it.take(1)+it.drop(2)})
Csv.write(times,[['Frame','Time_s']]+rows.drop(1).collect{it.take(2)})
def timed=settings+[expected_input_sha256:Installer.digest(noTime.bytes),timestamps_file:times.path,timestamps_sha256:Installer.digest(times.bytes)]
def timedResult=Engine.analyze(runtime,noTime,timed,root,job)
assert new File(timedResult,'measured_fluorescence.csv').bytes==noTime.bytes
assert Csv.read(new File(timedResult,'input_corrected_fluorescence.csv'))[0]==rows[0]
try {Engine.analyze(runtime,noTime,timed+[timestamps_sha256:'changed'],root,job);assert false}catch(IOException expected){assert expected.message.contains('timestamps changed')}
try {Engine.analyze(runtime,input,settings+[roi_registry_sha256:'changed'],root,job);assert false}catch(IOException expected){assert expected.message.contains('registry changed')}
try {Engine.analyze(runtime,input,settings+[sensitivity_file_sha256:'changed'],root,job);assert false}catch(IOException expected){assert expected.message.contains('diagnostic file changed')}
println 'PASS: native JSON -> shared R protocol, exact outputs, diagnostic snapshots, explicit timing sidecars and changed-source protection through the real Java bridge'

// Install exactly the same runtime users receive, then exercise the complete Java -> R -> CSV path.
import ij.gui.Roi

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Runtime=loader.loadClass('CalciumRRuntime'),Job=loader.loadClass('CalciumLocalJob'),Engine=loader.loadClass('CalciumSimpleCaEngine'),Csv=loader.loadClass('CalciumCsv')
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
try {
    Engine.analyze(runtime,new File('build/analysis-fixtures/input.csv').absoluteFile,settings+[expected_input_sha256:'outdated-measurement'],new File('build/analysis output').absoluteFile,job)
    assert false : 'Changed input must not be linked to old ROI snapshots'
} catch(IOException expected) {assert expected.message.contains('measurement file changed')}
// Bundled engine cache can repair a corrupted file, without fetching or running remote source.
def engine=Engine.unpack(root)
new File(engine,'runner.R').text='corrupted'
assert Engine.unpack(root)==engine
assert new File(engine,'runner.R').getText('UTF-8').startsWith('# Fiji adapter.')
println 'PASS: clean managed R setup, cached setup, real analysis, paths with spaces, embedded code repair, source snapshot and manifest'
} catch(Throwable problem) {
    if(job.logFile.isFile())println(job.logFile.getText('UTF-8'))
    throw problem
}

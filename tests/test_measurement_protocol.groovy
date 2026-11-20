import ij.gui.Roi
import ij.io.RoiEncoder
import ij.io.RoiDecoder

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Measurement=loader.loadClass('CalciumMeasurement'),Csv=loader.loadClass('CalciumCsv'),Steps=loader.loadClass('CalciumSteps')
def directory=new File('build/measurement-protocol').absoluteFile;directory.mkdirs()
def file=new File(directory,'precise.csv')
double measured=1.000123456789d-1d
file.text=Measurement.fullPrecision([[1.000123456789d,Double.NaN],[1.1d,2d]] as double[][],[[1d,1d],[1d,1d]] as double[][],['CellA','CellB'] as String[])
def rows=Csv.read(file)
assert Double.parseDouble(rows[1][1])==measured
assert rows[1][2]=='' && rows[2][2]=='1.0'

def inside=new Roi(70,70,12,12),edge=new Roi(1,1,12,12)
Measurement.ensureIds([inside,edge] as Roi[])
String id=Measurement.idFor(inside)
assert id==Measurement.idFor(inside.clone())
assert id==Measurement.idFor(new RoiDecoder(RoiEncoder.saveAsByteArray(inside),'cell.roi').getRoi())
def duplicate=inside.clone()
Roi[] duplicateSet=[inside,duplicate] as Roi[]
Measurement.ensureIds(duplicateSet)
assert Measurement.idFor(duplicateSet[1])!=id && Measurement.idFor(inside)==id
def steps=Steps.newInstance();steps.w=200;steps.h=200
def registry=Measurement.registry([edge,inside] as Roi[],[inside] as Roi[],[0] as int[],steps,'Baseline frames 1-20')
assert registry.size()==3 && registry[1][3]==false && registry[1][4]=='Edge exclusion policy'
assert registry[2][0]==id && registry[2][3]==true && registry[2][7]==0

def timeFile=new File(directory,'times.csv');timeFile.text='Frame,Time_s\n1,0\n2,0.73\n'
double[] times=Measurement.readTimes(timeFile,2)
assert times[1]==0.73d
def timed=new File(directory,'timed.csv');Measurement.attachTimes(file,timed,times)
assert Csv.read(timed)[0]==['Frame','Time_s','CellA','CellB']
assert Double.parseDouble(Csv.read(timed)[1][2])==measured
for(String invalid in ['Frame,Time_s\n1,0\n2,0\n','Frame,Time_s\n1,0\n3,1\n','Frame,Time_s\n1,0\n','Frame,Time_s\n1,0\n2,NaN\n']) {
    timeFile.text=invalid
    try {Measurement.readTimes(timeFile,2);assert false:'Invalid timestamps accepted'}catch(IOException expected){}
}
println 'PASS: full-precision corrected values, persistent unique ROI identities, premeasurement exclusions and exact timestamp mapping'

def Protocol=loader.loadClass('CalciumResponseProtocol')
def portable=Protocol.read(new File('tests/fixtures/response-protocols/golden-protocol.json'))
def inconsistent=new groovy.json.JsonSlurper().parseText(groovy.json.JsonOutput.toJson(portable))
inconsistent.windows[0].Baseline_s=99d
def badProtocol=new File(directory,'inconsistent.json');badProtocol.text=groovy.json.JsonOutput.toJson(inconsistent)
try {Protocol.read(badProtocol);assert false:'Conflicting baseline metadata accepted'}catch(IOException expected){assert expected.message.contains('Baseline_s')}
def legacy=Protocol.window('Stimulus','Capsaicin',0d,5d,10d)+[Duration_mode:'samples']
assert Protocol.payload([legacy],1,'Historical','user_defined').schema_version==1
['Duration_mode':'seconds','Min_duration_s':2d,'Min_baseline_samples':5,'Noise_floor_local_fraction':0.1d,'Expected_interval_s':0.5d].each{key,value->
    try {Protocol.payload([legacy+[(key):value]],1,'Historical','user_defined');assert false:'Ignored legacy setting accepted: '+key}
    catch(IllegalArgumentException expected){assert expected.message.contains('version-2-only')}
}
println 'PASS: native protocol rejects inconsistent baseline metadata and unsupported legacy settings'

// Two-window rules mirror the pinned engine (analysis/upstream/response_protocols.R), so these never reach R as exit code 1.
def stimulus=Protocol.window('Stimulus','Capsaicin',0d,20d,250d)
def refused={List windows,String text->try{Protocol.payload(windows,2,'Overlap','user_defined');assert false:'Engine-refused protocol accepted: '+text}catch(IllegalArgumentException expected){assert expected.message.contains(text):expected.message}}
// Stress test: KCl baseline 200-250 s inside a primary window that ends at 250 s.
refused([stimulus,Protocol.window('Control','KCl',200d,250d,299d)],'falls inside the Capsaicin response interval')
refused([stimulus,Protocol.window('Control','KCl',100d,240d,299d)],'response intervals overlap')
refused([Protocol.window('Stimulus','Capsaicin',150d,200d,250d),Protocol.window('Control','KCl',0d,100d,180d)],'Capsaicin baseline (150–200 s) falls inside the KCl response interval')
assert Protocol.payload([stimulus,Protocol.window('Control','KCl',250d,260d,299d)],2,'Adjacent','user_defined').windows.size()==2 : 'A control baseline that starts where the primary window ends is allowed, as in the engine'
assert Protocol.payload([Protocol.window('Stimulus','Capsaicin',0.5d,3d,14.5d)+[Duration_mode:'samples'],Protocol.window('Control','KCl',0.5d,3d,14.5d)+[Duration_mode:'samples']],1,'Legacy','user_defined').windows.size()==2 : 'The engine applies the two-window rules to version 2 only'
// Hand-edited or malformed files explain what to fix.
def edited=new groovy.json.JsonSlurper().parseText(groovy.json.JsonOutput.toJson(portable));edited.name='Renamed by hand'
try{Protocol.parse(edited);assert false}catch(IOException expected){assert expected.message.contains('"Renamed by hand" at the top of the file') && expected.message.contains('Make every copy match'):expected.message}
try{Protocol.parse('{"schema": "simpleca.response-protocol",');assert false}catch(IOException expected){assert expected.message.startsWith('This protocol file is not valid JSON.'):expected.message}
println 'PASS: protocols the engine would refuse for overlapping windows are caught first; edited or malformed files get plain messages'

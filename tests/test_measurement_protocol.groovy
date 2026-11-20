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
// The ID comes from the outline itself, so the same outline has the same ID in another session, and another outline another ID.
assert Measurement.idFor(new Roi(70,70,12,12))==id : 'a new outline with the same place and pixels gets the same ID'
assert Measurement.idFor(new Roi(71,70,12,12))!=id && Measurement.idFor(new Roi(70,70,12,13))!=id
assert Measurement.idFor(new ij.gui.OvalRoi(30,30,20,14))==Measurement.idFor(new ij.gui.OvalRoi(30,30,20,14)) && Measurement.idFor(new ij.gui.OvalRoi(30,30,20,14))==~/Cell_[a-f0-9]{32}/
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

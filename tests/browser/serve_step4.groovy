// A clearly synthetic display fixture for the actual local web server and browser UI.
// Response classification is independently exercised through real R by the adapter gates.
import ij.*
import ij.gui.Roi
import ij.process.*
import groovy.json.JsonOutput

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Csv=loader.loadClass('CalciumCsv'),Protocol=loader.loadClass('CalciumResponseProtocol')
def root=new File('build/web-step4').absoluteFile;root.mkdirs();new File(root,'stop').delete();new File(root,'url.txt').delete()
def source=new File('tests/fixtures/simpleca').absoluteFile,fixture=new File(root,'browser-results');fixture.mkdirs()
source.listFiles().findAll{it.isFile()}.each{new File(fixture,it.name).bytes=it.bytes}
List clock=(0..<30).collect{it/2d+(it>=15?5d:0d)}
['normalized_traces.csv','mean_trace.csv'].each{file->def rows=Csv.read(new File(fixture,file));rows.drop(1).eachWithIndex{row,i->row[0]=clock[i]};Csv.write(new File(fixture,file),rows)}
// One cell's corrected fluorescence goes below zero at the last frame, so the page offers to leave it out of the percentages.
def dipped=Csv.read(new File(fixture,'normalized_traces.csv'));dipped[dipped.size()-1][dipped[0].indexOf('Cell005')]='-1.5';Csv.write(new File(fixture,'normalized_traces.csv'),dipped)
// The control follows the stimulus, as the engine requires: neither baseline may fall inside the other's response window.
def primary=Protocol.window('Stimulus','Capsaicin',0.5d,3d,7d)
def control=Protocol.window('Control','Ionomycin',7d,12.5d,19.5d);control.Floor_local_fraction=0.4d;control.SD_multiplier=5d
Map protocol=Protocol.payload([primary,control],2,'Synthetic browser display protocol','provisional')
['protocol.json','applied_protocol.json'].each{new File(fixture,it).text=Protocol.json(protocol)}
def callRows=Csv.read(new File(fixture,'responder_calls.csv'));List header=callRows[0]
header.add('QC_flags');callRows.drop(1).each{it.add('DURATION_BELOW_SAMPLE_INTERVAL')}
// The stimulus is added when the response window starts, as in a real run, so the plot has a stimulus line to place.
header.add('Stimulus_onset_s');callRows.drop(1).each{it.add(it[header.indexOf('Start_s')])}
def set={List row,String name,Object value->int index=header.indexOf(name);if(index>=0)row[index]=''+value}
callRows.drop(1).each{row->set(row,'End_s',19.5);set(row,'Effective_end_s',19.5)}
def extra=callRows.drop(1).collect{original->def row=new ArrayList(original);set(row,'Role','Control');set(row,'Stimulus','Ionomycin');set(row,'Floor_local_fraction',0.4);set(row,'Floor_threshold',0.4);set(row,'Threshold_increase',0.4);set(row,'SD_multiplier',5);set(row,'Baseline_mean',0.5)
    if(row[header.indexOf('Cell')]=='Cell001'){set(row,'Status','Did not respond');set(row,'Reason','Synthetic control fixture below threshold')}
    if(row[header.indexOf('Cell')]=='Cell002'){set(row,'Status','Responded');set(row,'Reason','Synthetic control fixture sustained crossing')}
    if(row[header.indexOf('Cell')]=='Cell005'){set(row,'Status','Did not respond');set(row,'Reason','Synthetic control fixture below threshold')}
    row}
Csv.write(new File(fixture,'responder_calls.csv'),callRows+extra)
def counts=Csv.read(new File(fixture,'responder_counts.csv'));def controlCount=new ArrayList(counts[1]);['Role':'Control','Stimulus':'Ionomycin','Responded':'1','Did_not_respond':'2','Percent_responding':'33.3333333333333'].each{k,v->int i=counts[0].indexOf(k);if(i>=0)controlCount[i]=v};Csv.write(new File(fixture,'responder_counts.csv'),counts+[controlCount])
def manifest=Csv.read(new File(fixture,'analysis_manifest.csv'));manifest<<['time_source','explicit_timestamps'];manifest<<['protocol_schema_version','2'];manifest<<['protocol_hash','synthetic-browser-fixture'];manifest<<['fixture_note','Synthetic display stand-in; no responder calculation in browser fixture'];Csv.write(new File(fixture,'analysis_manifest.csv'),manifest)
def quality=Csv.read(new File(fixture,'cell_quality.csv'));quality[0].add('QC_flags')
quality.drop(1).each{row->row.add(row[0]=='Cell005'?'reference_below_3_baseline_sd_review':'')}
Csv.write(new File(fixture,'cell_quality.csv'),quality)
// An empty optional table should explain that it has no rows.
Csv.write(new File(fixture,'group_counts.csv'),[['Condition','Percent_responding']])
def frames=new ImageStack(250,150)
(1..30).each{int f->def ip=new ShortProcessor(250,150);ip.set(100);if(f>=8&&f<=17){ip.setRoi(30,50,20,20);ip.set(140);ip.resetRoi()};frames.addSlice(''+f,ip)}
def colour=new ColorProcessor(250,150);colour.setColor(0x103a10);colour.fill()
(1..6).each{colour.setColor(0x30c030);colour.fillOval(it*30+2,52,16,16)}
// The recording says how large its pixels are, as a microscope's file does: 0.8 µm, so the field is 200 µm across.
def recording=new ImagePlus('Recording',frames);recording.getCalibration().pixelWidth=0.8d;recording.getCalibration().pixelHeight=0.8d;recording.getCalibration().setUnit('micron')
def data=[csv:new File(source,'metrics.csv'),folder:root,recording:'Fixture recording',frames:30,hz:2d,baselineStart:1,baselineEnd:5,key:'fixture',measurementId:'synthetic-browser-fixture',
    picture:new ImagePlus('Colour picture',colour),movie:recording,names:(1..6).collect{String.format('Cell%03d',it)},rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]]
def session=loader.loadClass('CalciumAnalysisSession').newInstance(data);session.runtime.home=new File(root,'runtime')
session.runtime.metaClass.environment={->root} // Display test only; never downloads or invokes R.
session.draft.protocol=protocol;session.draft.timestamps=[name:'fixture-clock.csv',content:'Frame,Time_s\n'+clock.withIndex().collect{t,i->(i+1)+','+t}.join('\n')+'\n']
session.open={File target->new File(root,'opened.txt').text=target.name}
session.performAnalysis={rt,input,settings,out,job->new File(root,'sent.json').text=JsonOutput.prettyPrint(JsonOutput.toJson(settings));new File(root,'settings.txt').text=settings.collect{k,v->k+'='+v}.join('\n')
    job.logFile.parentFile.mkdirs();job.logFile.text='Synthetic analysis log';(1..4).each{job.check();job.progress.call(it*20,'Scoring synthetic display fixture');Thread.sleep(150)};fixture}
def server=loader.loadClass('CalciumWebServer').newInstance(session:session);server.start()
new File(root,'url.txt').text=server.url();println 'Step 4 page: '+server.url()
long end=System.currentTimeMillis()+600000L
while(!new File(root,'stop').exists()&&System.currentTimeMillis()<end)Thread.sleep(200)
server.stop();System.exit(0)

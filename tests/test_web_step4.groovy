// The Step 4 page's server: the files it serves, its access key, the measurement it describes, the settings rules,
// the saved results as JSON and the movie frames. Runs without a display and without a browser.
import ij.*
import ij.gui.Roi
import ij.process.*
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Session=loader.loadClass('CalciumAnalysisSession'),Server=loader.loadClass('CalciumWebServer'),Settings=loader.loadClass('CalciumAnalysisSettings')
def root=new File('build/web-step4').absoluteFile;root.deleteDir();root.mkdirs()
def fixture=new File('tests/fixtures/simpleca').absoluteFile
// a small recording that matches the picture: cell 1's square is brighter in frames 10-20
def frames=new ImageStack(250,150)
(1..30).each{int f->def ip=new ShortProcessor(250,150);ip.set(100);if(f>=10&&f<=20){ip.setRoi(30,50,20,20);ip.set(140);ip.resetRoi()};frames.addSlice(''+f,ip)}
def recording=new ImagePlus('Recording',frames),picture=new ImagePlus('Colour picture',new ColorProcessor(250,150))
def names=(1..6).collect{String.format('Cell%03d',it)}
def data=[csv:new File(fixture,'metrics.csv'),folder:root,recording:'Fixture recording',frames:30,hz:2.0000001d,baselineStart:1,baselineEnd:5,key:'first',
    picture:picture,movie:recording,names:names,rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]]
def session=Session.newInstance(data);session.runtime.home=new File(root,'runtime')
def opened=[];session.open={File target->opened<<target.name}
def server=Server.newInstance(session:session);server.start()
def slurp=new JsonSlurper()
def call={String path,Object body=null,String base=null->
    def c=(HttpURLConnection)new URL((base?:server.url())+path).openConnection(java.net.Proxy.NO_PROXY);c.connectTimeout=5000;c.readTimeout=15000
    if(body!=null){c.requestMethod='POST';c.doOutput=true;c.setRequestProperty('Content-Type','application/json');c.outputStream.withCloseable{it.write(JsonOutput.toJson(body).getBytes('UTF-8'))}}
    int status=c.responseCode;def stream=status<400?c.inputStream:c.errorStream;byte[] bytes=stream!=null?stream.bytes:new byte[0]
    [status:status,bytes:bytes,type:c.getHeaderField('Content-Type'),cache:c.getHeaderField('Cache-Control'),text:new String(bytes,'UTF-8')]}
def json={String path,Object body=null->def r=call(path,body);[r.status,r.bytes.length>0?slurp.parseText(r.text):null]}
def settle={->long end=System.nanoTime()+8000000000L;while(session.busy&&System.nanoTime()<end)Thread.sleep(20);assert !session.busy}
def image={byte[] bytes->javax.imageio.ImageIO.read(new ByteArrayInputStream(bytes))}
try {
    // Reachable from this computer only, behind a long random key; the page comes from the copy embedded in the plugin.
    assert server.url() ==~ /http:\/\/127\.0\.0\.1:\d+\/[0-9a-f]{64}\// && server.server.address.address.isLoopbackAddress()
    [['','step4.html','text/html'],['app.css','step4.css','text/css'],['app.js','step4.js','text/javascript']].each{path,file,type->
        def r=call(path);assert r.status==200 && r.type.startsWith(type) && r.cache=='no-store' && r.bytes==new File('web',file).bytes : file+' is not the embedded copy; run tools/embed_analysis.py'}
    def port=server.port
    assert call('api/state',null,'http://127.0.0.1:'+port+'/'+('0'*64)+'/').status==404 && call('api/state',null,'http://127.0.0.1:'+port+'/').status==404
    def socket=new Socket(InetAddress.getLoopbackAddress(),port)
    socket.outputStream.write(('GET /'+server.token+'/api/state HTTP/1.1\r\nHost: example.org\r\nConnection: close\r\n\r\n').getBytes('US-ASCII'));socket.outputStream.flush()
    assert new String(socket.inputStream.bytes,'US-ASCII').startsWith('HTTP/1.1 404') : 'A request that names another host must be refused';socket.close()

    // What the page is told about the measurement, and the settings it should offer.
    def (code,about)=json('api/session')
    assert code==200 && about.recording=='Fixture recording' && about.names==names && about.frames==30 && about.width==250 && about.height==150 && about.hasPicture && about.hasMovie
    assert about.rois.size()==6 && about.rois[0].size()==4 && about.rois[0].collect{it[0]}.min()==30 && about.rois[0].collect{it[1]}.max()==70
    def box=about.settings.collectEntries{[it.key,it]}
    assert about.settings*.key==['sampling_hz','baseline_start','baseline_end','stimulus','stimulus_start','stimulus_end','response_baseline_seconds','amplitude_floor','noise_multiplier','consecutive_samples','kcl_start','kcl_end']
    assert box.sampling_hz.shown=='2' && box.sampling_hz.exact=='2.0000001' && box.baseline_end.shown=='5' && box.amplitude_floor.shown=='0.1' && box.amplitude_floor.exact=='0.10'
    assert box.stimulus_start.exact==''+((5-1)/2.0000001d) && box.stimulus_start.shown=='2' && box.stimulus_end.exact==''+((30-1)/2.0000001d) && box.kcl_start.shown==''
    def state=json('api/state')[1]
    assert !state.busy && state.version==0 && !state.engineReady && state.rest==[1,5] && state.key=='first' && state.message.contains('Check the settings') && call('api/results').text=='{}'

    // The same checks as the Fiji window, before anything runs.
    def typed=about.settings.collectEntries{[it.key,it.shown]}
    def refused=json('api/analyze',[values:typed,confirmed:false,kcl:false])
    assert refused[0]==400 && refused[1].error.contains('Confirm') && !session.busy
    refused=json('api/analyze',[values:typed+[sampling_hz:'0'],confirmed:true,kcl:false]);assert refused[0]==400 && refused[1].error.contains('frame rate')
    refused=json('api/analyze',[values:typed+[kcl_start:'soon'],confirmed:true,kcl:true]);assert refused[0]==400 && refused[1].error.contains('kcl_start')

    // A run: an untouched box sends its exact value, an edited box what was typed; a second run cannot start meanwhile.
    def sent=null;def gate=new CountDownLatch(1)
    session.performAnalysis={rt,input,settings,out,job->sent=settings;job.progress.call(40,'Scoring responses');gate.await(5,TimeUnit.SECONDS);job.check();fixture}
    def started=json('api/analyze',[values:typed+[stimulus:'Capsaicin',stimulus_start:'3',stimulus_end:' 14 ',response_baseline_seconds:'2.5'],confirmed:true,kcl:false])
    assert started[0]==200 && started[1].busy
    long wait=System.nanoTime()+3000000000L;while(sent==null&&System.nanoTime()<wait)Thread.sleep(10)
    state=json('api/state')[1];assert state.busy && state.percent==40 && state.stage=='Scoring responses'
    refused=json('api/analyze',[values:typed,confirmed:true,kcl:false]);assert refused[0]==400 && refused[1].error.contains('already running')
    gate.countDown();settle()
    assert sent.sampling_hz=='2.0000001' && sent.amplitude_floor=='0.10' && sent.stimulus=='Capsaicin' && sent.stimulus_start=='3' && sent.stimulus_end=='14' && sent.response_baseline_seconds=='2.5'
    assert sent.baseline_start=='1' && sent.baseline_end=='5' && sent.kcl_enabled=='false' && sent.recording=='Fixture recording' && sent.input_mode=='background_corrected_fluorescence'
    state=json('api/state')[1];assert state.version==1 && state.tone=='good' && state.message.startsWith('Complete') && state.folder=='simpleca' && state.percent==100

    // The saved tables, as the page reads them. Nothing is recalculated; a missing number is null, never NaN.
    def reply=call('api/results');assert reply.status==200 && !reply.text.contains('NaN')
    def results=slurp.parseText(reply.text)
    assert results.time.size()==30 && results.time[1]==0.5 && results.mean.size()==30 && results.sem.size()==30 && results.folder=='simpleca'
    assert results.traces.keySet()==['Cell001','Cell002','Cell005','Cell006'] as Set && results.traces.Cell001.max()==1
    assert results.calls.Cell001.Status=='Responded' && results.calls.Cell002.Status=='Did not respond' && results.calls.Cell003.Status=='Excluded' && results.calls.Cell001.Peak_increase=='1'
    assert results.counts.Registered=='6' && results.counts.Responded=='2' && results.counts.Evaluable=='3' && results.counts.Percent_responding.startsWith('66.6')
    assert results.stimulus==[3,14.5] && results.baseline[0]==0 && results.tables.metrics.header[0]=='Cell' && results.tables.metrics.rows.size()==6 && results.tables.calls.rows.size()==6 && results.tables.counts.rows.size()==1
    assert Session.json([a:Double.NaN,b:[1.5d,Double.POSITIVE_INFINITY],c:[2d] as double[]])=='{"a":null,"b":[1.5,null],"c":[2.0]}'

    // Movie frames for the page's time cursor: the change from the resting frames, or the frame as recorded. The recording is only read.
    def gray={byte[] bytes,int x,int y->image(bytes).getRGB(x,y)&255}
    reply=call('api/frame?f=15&mode=change');assert reply.status==200 && reply.type=='image/jpeg' && image(reply.bytes).width==250 && image(reply.bytes).height==150
    assert gray(reply.bytes,40,60)>200 && gray(reply.bytes,100,60)<60 : 'Frame 15 minus the resting frames: only the responding square is bright'
    assert gray(call('api/frame?f=3&mode=change').bytes,40,60)<60 && call('api/frame?f=999&mode=change').status==200
    reply=call('api/frame?f=15&mode=raw');assert gray(reply.bytes,40,60)>gray(reply.bytes,100,60)
    assert recording.getRoi()==null && recording.getStackSize()==30 && recording.getStack().getProcessor(15).get(40,60)==140
    reply=call('api/picture');assert reply.type=='image/png' && image(reply.bytes).width==250

    // Opening the folder, cancelling, and a failed run: completed results stay.
    assert json('api/open',[what:'folder'])[0]==200 && opened==['simpleca'] && json('api/open',[what:'log'])[0]==404
    session.performAnalysis={rt,input,settings,out,job->while(!job.cancelled)Thread.sleep(20);job.check()}
    assert json('api/analyze',[values:typed+[stimulus_start:'3',stimulus_end:'14.5',response_baseline_seconds:'2.5'],confirmed:true,kcl:false])[0]==200
    assert json('api/cancel',[:])[0]==200;settle()
    state=json('api/state')[1];assert state.message.startsWith('Cancelled') && state.version==1 && slurp.parseText(call('api/results').text).traces.size()==4
    session.performAnalysis={rt,input,settings,out,job->throw new IOException('Simulated engine error')}
    json('api/analyze',[values:typed+[stimulus_start:'3',stimulus_end:'14.5',response_baseline_seconds:'2.5'],confirmed:true,kcl:false]);settle()
    state=json('api/state')[1];assert state.message.contains('Simulated engine error') && state.tone=='bad' && state.version==1

    // Measuring again gives the page a new measurement; without one it says so; after stopping nothing answers.
    server.session=Session.newInstance(data+[key:'second']);state=json('api/state')[1];assert state.key=='second' && state.version==0
    server.session=null;assert json('api/state')[0]==409 && call('').status==200
    assert Settings.check([sampling_hz:'2',baseline_start:'1',baseline_end:'5',stimulus:'S',stimulus_start:'3',stimulus_end:'14',response_baseline_seconds:'2',amplitude_floor:'0.1',noise_multiplier:'3',consecutive_samples:'2'],true,false,[frames:30,recording:'R']).kcl_enabled=='false'
    server.stop()
    try {call('api/state');assert false : 'The server must stop answering'} catch(IOException expected) { }
    println 'PASS: Step 4 page served from the embedded copy on this computer only, access key and host check, session facts, settings rules with exact untouched values, results JSON, movie frames, cancel, failure and re-measurement'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {server.stop();recording.close();picture.close()}
System.exit(0)

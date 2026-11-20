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

    def Protocol=loader.loadClass('CalciumResponseProtocol'),Csv=loader.loadClass('CalciumCsv'),Installer=loader.loadClass('CalciumInstaller')
    def (code,about)=json('api/session')
    assert code==200 && about.recording=='Fixture recording' && about.names==names && about.frames==30 && about.hasMovie && about.hasPicture
    assert about.key==Settings.identity(data) && about.rois.size()==6 && about.rois[0].size()==4
    assert about.draft.values.sampling_hz=='2.0000001' && about.draft.values.baseline_end=='5'
    assert about.draft.protocol.schema_version==2 && about.draft.protocol.windows[0].Duration_mode=='seconds' && about.draft.protocol.windows[0].Min_duration_s==1
    assert about.timingProvisional && about.timeSource=='constant_fps' && Math.abs(about.frameTimes[1]-1/2.0000001d)<1e-14
    def state=json('api/state')[1]
    assert !state.busy && state.version==0 && !state.engineReady && state.rest==[1,5] && !state.settingsPending && call('api/results').text=='{}'
    // Raw recording frames are usable before any analysis, without altering the source recording.
    def gray={byte[] bytes,int x,int y->image(bytes).getRGB(x,y)&255}
    def raw=call('api/frame?f=15&mode=raw');assert raw.status==200 && gray(raw.bytes,40,60)>gray(raw.bytes,100,60)
    assert call('api/frame?f=3').status==200
    ['0','-5','31','abc','1.5',''].each{f->assert call('api/frame?f='+f+'&mode=raw').status==400}
    assert call('api/frame?f=3&mode=unknown').status==400
    assert about.draft.protocol.windows[0].Start_s>about.frameTimes[(about.draft.values.baseline_end as int)-1] : 'Suggested onset must follow the inclusive F0 baseline'
    ['0','-1','abc',''].each{hz->
        def invalid=Settings.copy(about.draft);invalid.values.sampling_hz=hz
        assert Settings.diagnose(invalid,data).errors
        assert json('api/analyze',[draft:invalid,confirmed:true])[0]==400 && !session.busy
    }
    [[baseline_start:'0'],[baseline_end:'5000'],[baseline_start:'8',baseline_end:'2']].each{changes->
        def invalid=Settings.copy(about.draft);invalid.values.putAll(changes);assert Settings.diagnose(invalid,data).errors
    }
    def badRule=Settings.copy(about.draft);badRule.protocol.windows[0].Floor_local_fraction=-0.1d
    assert Settings.diagnose(badRule,data).errors
    def overlap=Settings.copy(about.draft);overlap.values.baseline_end='10'
    assert Settings.diagnose(overlap,data).warnings.any{it.contains('F₀ baseline includes')}

    def pending=Settings.copy(about.draft)
    pending.values.sampling_hz=''
    assert json('api/draft',[draft:pending,revision:1])[0]==200
    assert json('api/state')[1].frameTimes==null && json('api/state')[1].timeSource=='unknown'
    // Stale HTTP edits cannot overwrite a newer draft revision.
    assert json('api/draft',[draft:about.draft,revision:0])[0]==200 && json('api/session')[1].draft.values.sampling_hz==''
    pending=Settings.copy(about.draft)
    def primary=Protocol.window('Stimulus','Capsaicin',0.5d,3d,8d),control=Protocol.window('Control','Ionomycin',8d,10d,14.5d)   // separate windows, as the engine requires
    control.Floor_local_fraction=0.35d;control.SD_multiplier=5d;control.Duration_mode='samples';control.Consecutive_samples=4
    pending.protocol=Protocol.payload([primary,control],2,'Browser protocol','user_defined')
    def imported=json('api/protocol',[protocol:pending.protocol]);assert imported[0]==200 && imported[1].protocol.windows[1].Consecutive_samples==4
    def inconsistent=Settings.copy(pending.protocol);inconsistent.windows[0].Baseline_s=999
    assert json('api/protocol',[protocol:inconsistent])[0]==400
    def unknown=Settings.copy(pending.protocol);unknown.schema_version=99
    assert json('api/protocol',[protocol:unknown])[0]==400
    assert json('api/analyze',[draft:pending,confirmed:false])[0]==400
    assert json('api/analyze',[values:[:],confirmed:true])[0]==400 : 'The former scalar path must not silently run the legacy engine'
    assert json('api/analyze',[draft:pending+[values:pending.values+[sampling_hz:'0']],confirmed:true])[0]==400
    pending.timestamps=[name:'../../clock.csv',content:'Frame,Time_s\n'+(0..<30).collect{(it+1)+','+(it/2d)}.join('\n')+'\n']
    pending.values.condition='Vehicle';pending.values.biological_unit='Culture1'
    pending.sensitivity=[name:'sensitivity.csv',content:'Floor_local_fraction,SD_multiplier,Min_duration_s\n0.2,3,1\n']
    assert json('api/draft',[draft:pending,revision:2])[0]==200
    state=json('api/state')[1];assert state.timeSource=='explicit_timestamps' && state.frameTimes[1]==0.5 && state.timingProvisional
    // A confirmed run snapshots a full protocol and the exact uploaded files, never user-controlled server paths.
    def sent=null;def gate=new CountDownLatch(1)
    session.performAnalysis={rt,input,settings,out,job->sent=settings
        assert new File(settings.timestamps_file).name=='timestamps.csv' && new File(settings.timestamps_file).text==pending.timestamps.content
        assert settings.timestamps_sha256==Installer.digest(new File(settings.timestamps_file).bytes)
        assert new File(settings.sensitivity_file).name=='sensitivity.csv'
        job.progress.call(40,'Scoring responses');gate.await(5,TimeUnit.SECONDS);job.check();fixture}
    def started=json('api/analyze',[draft:pending,confirmed:true]);assert started[0]==200 && started[1].busy
    long wait=System.nanoTime()+3000000000L;while(sent==null&&System.nanoTime()<wait)Thread.sleep(10)
    assert json('api/analyze',[draft:pending,confirmed:true])[0]==400
    def submitted=slurp.parseText(sent.protocol_json)
    assert submitted.schema_version==2 && submitted.windows[1].Floor_local_fraction==0.35 && submitted.windows[1].Consecutive_samples==4
    assert sent.sampling_hz=='' && sent.condition=='Vehicle' && sent.biological_unit=='Culture1' && sent.input_mode=='background_corrected_fluorescence'
    assert sent.draft_sha256==Settings.fingerprint(pending)
    def edited=Settings.copy(pending);edited.values.condition='Changed while running'
    assert json('api/draft',[draft:edited,revision:3])[0]==200
    gate.countDown();settle()
    state=json('api/state')[1];assert state.version==1 && state.tone=='good' && state.settingsPending && !state.timingProvisional
    assert !new File(sent.timestamps_file).exists() : 'Private request files are cleaned after the engine snapshots them'
    def reply=call('api/results'),results=slurp.parseText(reply.text)
    assert !reply.text.contains('NaN') && results.time.size()==30 && results.frameTimes==results.time && results.time[1]==0.5
    assert results.traces.keySet()==['Cell001','Cell002','Cell005','Cell006'] as Set && results.traces.Cell001.max()==1
    assert results.calls.Cell001.Status=='Responded' && results.calls.Cell003.Status=='Excluded'
    assert results.counts.Responded=='2' && results.counts.Evaluable=='3'
    assert results.windows==[[key:'Stimulus::Capsaicin',role:'Stimulus',stimulus:'Capsaicin']]
    assert results.callsByWindow['Stimulus::Capsaicin']==results.calls && results.countsByWindow['Stimulus::Capsaicin']==results.counts
    assert results.baseline==[0,2] && results.stimulus==[3,14.5]
    assert results.tables.quality.rows.size()==6 && results.tables.containsKey('negativeControlCounts') && results.manifest.input_mode=='background_corrected_fluorescence'
    assert results.appliedFingerprint==Settings.fingerprint(pending) && results.appliedDraft==pending
    assert Session.json([a:Double.NaN,b:[1.5d,Double.POSITIVE_INFINITY],c:[2d] as double[]])=='{"a":null,"b":[1.5,null],"c":[2.0]}'
    // Applied time and normalization remain stable through incomplete or changed draft files.
    edited.timestamps.content='Time_s\n'+(0..<30).collect{it/2d+10d}.join('\n')+'\n';edited.values.baseline_start='2';edited.values.baseline_end='9'
    json('api/draft',[draft:edited,revision:4]);state=json('api/state')[1]
    assert state.frameTimes==results.time && state.rest==[1,5] && state.settingsPending
    reply=call('api/frame?f=15&mode=change');assert reply.status==200 && reply.type=='image/jpeg' && image(reply.bytes).width==250
    assert gray(reply.bytes,40,60)>200 && gray(reply.bytes,100,60)<60
    assert gray(call('api/frame?f=3&mode=change').bytes,40,60)<60 && call('api/frame?f=999&mode=raw').status==400
    assert recording.roi==null && recording.stackSize==30 && recording.getStack().getProcessor(15).get(40,60)==140
    assert call('api/picture').type=='image/png'
    // Failure/cancel/incomplete output cannot replace a completed result or its applied snapshot.
    assert json('api/open',[what:'folder'])[0]==200 && opened==['simpleca'] && json('api/open',[what:'log'])[0]==404
    session.performAnalysis={rt,input,settings,out,job->while(!job.cancelled)Thread.sleep(20);job.check()}
    assert json('api/analyze',[draft:pending,confirmed:true])[0]==200
    assert json('api/cancel',[:])[0]==200;settle();assert json('api/state')[1].version==1
    session.performAnalysis={rt,input,settings,out,job->throw new IOException('Simulated engine error')}
    json('api/analyze',[draft:pending,confirmed:true]);settle();assert json('api/state')[1].tone=='bad' && json('api/state')[1].version==1
    File incomplete=new File(root,'incomplete');incomplete.mkdirs()
    try{session.load(incomplete);assert false}catch(IOException expected){assert expected.message.contains('COMPLETE')}
    assert session.results==fixture && session.version==1
    // A source generation, recording identity or source hash change requires a new browser session.
    assert Settings.identity(data+[recording:'Another recording'])!=session.measurementIdentity
    assert Settings.identity(data+[inputHash:'different'])!=session.measurementIdentity
    assert Settings.identity(data+[measurementId:'new measurement'])!=session.measurementIdentity
    server.session=Session.newInstance(data+[key:'second']);assert json('api/state')[1].key!=about.key && json('api/state')[1].version==0
    server.session=null;assert json('api/state')[0]==409 && call('').status==200
    // Remeasuring gets a distinct URL and read-only input snapshot; the first page remains reachable.
    def Pages=loader.loadClass('CalciumWebPages'),pages=Pages.newInstance(snapshotRoot:new File(root,'snapshots'))
    def original=new File(root,'measurement.csv');original.text='original measured values'
    def oldData=data+[csv:original,measurementId:'first snapshot']
    try {
        def firstPage=pages.page(oldData),secondPage=pages.page(oldData+[measurementId:'second snapshot'])
        assert pages.page(oldData).is(firstPage) && firstPage.url()!=secondPage.url()
        original.text='remeasured values'
        assert firstPage.session.measurement.csv.text=='original measured values' && secondPage.session.measurement.csv.text=='original measured values'
        assert call('api/state',null,firstPage.url()).status==200 && call('api/state',null,secondPage.url()).status==200
        File snapshot=firstPage.session.measurement.csv
        assert call('api/close',[:],firstPage.url()).status==200
        long deadline=System.nanoTime()+3000000000L;while(!firstPage.closed&&System.nanoTime()<deadline)Thread.sleep(10)
        assert firstPage.closed && call('api/state',null,secondPage.url()).status==200
        // stop() also cleans the executor and private measurement files.
        firstPage.stop();assert !snapshot.exists()
    }finally{pages.closeAll()}
    server.stop();try{call('api/state');assert false}catch(IOException expected){}
    println 'PASS: local browser HTTP/token, versioned draft protocol, exact values/uploads, pre-analysis recording transport, applied snapshots, all result tables, failure/cancel and measurement identity'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {server.stop();recording.close();picture.close()}
System.exit(0)

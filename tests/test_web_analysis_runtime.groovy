// Invoked by both portable and managed-runtime gates with loader/runtime/root bindings.
// This exercises the actual local HTTP API -> immutable Session -> pinned R engine,
// and compares saved calls/counts to independent hand-calculated CSV constants.
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import ij.ImagePlus
import ij.ImageStack
import ij.gui.Roi
import ij.process.ColorProcessor
import ij.process.ShortProcessor

def Session=loader.loadClass('CalciumAnalysisSession'),Server=loader.loadClass('CalciumWebServer'),Settings=loader.loadClass('CalciumAnalysisSettings')
def Csv=loader.loadClass('CalciumCsv'),Protocol=loader.loadClass('CalciumResponseProtocol'),Installer=loader.loadClass('CalciumInstaller')
def shared=Csv.read(new File('build/analysis-fixtures/shared-golden-input.csv'))
def output=new File(root,'real web session');output.mkdirs()
def input=new File(output,'measured native.csv')
Csv.write(input,shared.collect{it.take(1)+it.drop(2)})
List times=shared.drop(1).collect{CalciumCsv.number(it[1])};times[1]=0.2d // Uneven early timing leaves the hand-calculated response windows unchanged.
List cells=shared[0].drop(2)
def registry=new File(output,'accepted registry.csv')
Csv.write(registry,[['Cell','ROI_ID','Accepted','Included_in_measurement','Exclusion_reason','Selection_basis','Background_pixels']]+cells.collect{[it,it,true,true,'','Baseline review',100]}+[['Edge','Edge',true,false,'Edge exclusion policy','Baseline review','']])
def stack=new ImageStack(64,48)
(1..times.size()).each{f->def ip=new ShortProcessor(64,48);ip.set(100);if(f>10){ip.setRoi(4,4,10,10);ip.set(150);ip.resetRoi()};stack.addSlice(''+f,ip)}
def recording=new ImagePlus('Original acquisition',stack),picture=new ImagePlus('Measured overview',new ColorProcessor(64,48))
def data=[csv:input,folder:output,recording:'Golden field',frames:times.size(),hz:2d,baselineStart:1,baselineEnd:4,names:cells,rois:cells.withIndex().collect{name,i->new Roi(2+i*8,4,6,8)} as Roi[],
    movie:recording,picture:picture,key:'golden measured ROIs',measurementId:'real-http-session',inputHash:Installer.digest(input.bytes),registry:registry,registryHash:Installer.digest(registry.bytes),selectionBasis:'Baseline review']
def session=Session.newInstance(data);session.runtime=runtime
Map protocol=Protocol.read(new File('tests/fixtures/response-protocols/golden-protocol.json'))
List overrides=[[Field:'Golden field',Role:'Stimulus',Start_s:5d]]
protocol=Protocol.payload(protocol.windows,2,protocol.name,protocol.validation_status,overrides,protocol.sd_estimator)   // read from a file: keeps its estimator
Map draft=Settings.copy(session.draft);draft.protocol=protocol;draft.values.condition='Vehicle';draft.values.biological_unit='Culture1'
draft.timestamps=[name:'actual-clock.csv',content:'Frame,Time_s\n'+times.withIndex().collect{time,i->(i+1)+','+time}.join('\n')+'\n']
draft.sensitivity=[name:'grid.csv',content:new File('build/analysis-fixtures/golden-diagnostics/sensitivity_grid.csv').getText('UTF-8')]
draft.negative_control=[name:'vehicle.csv',content:new File('build/analysis-fixtures/golden-diagnostics/negative_control_windows.csv').getText('UTF-8')]
def server=Server.newInstance(session:session);server.start()
def slurp=new JsonSlurper()
def call={String path,Map body=null->def connection=(HttpURLConnection)new URL(server.url()+path).openConnection(java.net.Proxy.NO_PROXY);connection.connectTimeout=10000;connection.readTimeout=10000
    if(body!=null){connection.requestMethod='POST';connection.doOutput=true;connection.setRequestProperty('Content-Type','application/json');connection.outputStream.withCloseable{it.write(JsonOutput.toJson(body).getBytes('UTF-8'))}}
    int code=connection.responseCode;def stream=code<400?connection.inputStream:connection.errorStream;String text=stream?.getText('UTF-8');[code:code,body:text?slurp.parseText(text):null]}
def waitForRun={->long end=System.nanoTime()+120000000000L;while(session.busy&&System.nanoTime()<end)Thread.sleep(25);assert !session.busy : 'R-backed browser analysis timed out';def state=call('api/state').body;if(state.tone!='good'&&session.job?.logFile?.isFile())println session.job.logFile.getText('UTF-8');assert state.tone=='good' : state.message}
def checkIndependent={Map table,File expected->def rows=Csv.read(expected),actual=table.rows;assert actual.size()==rows.size()-1
    rows.drop(1).eachWithIndex{row,r->rows[0].eachWithIndex{column,c->int index=table.header.indexOf(column);assert index>=0
        String got=actual[r][index],want=row[c];double gn=Csv.number(got),wn=Csv.number(want)
        if(Double.isFinite(wn))assert Double.isFinite(gn)&&Math.abs(gn-wn)<1e-10 : column+' row '+r+' got '+got+' expected '+want
        else assert got==want : column+' row '+r+' got '+got+' expected '+want}}
}
try {
    assert call('api/draft',[draft:draft,revision:1]).code==200
    assert call('api/state').body.frameTimes[1]==0.2 && call('api/state').body.timingProvisional
    assert call('api/analyze',[draft:draft,confirmed:true]).code==200;waitForRun()
    Map result=call('api/results').body,state=call('api/state').body
    assert new File(session.results,'COMPLETE').isFile() && !state.settingsPending && !state.timingProvisional
    checkIndependent(result.tables.calls,new File('tests/fixtures/response-protocols/golden-expected-calls.csv'))
    checkIndependent(result.tables.counts,new File('tests/fixtures/response-protocols/golden-expected-counts.csv'))
    assert result.counts.Percent_responding=='25' && result.countsByWindow['Control::Ionomycin'].Percent_responding=='80'
    assert result.callsByWindow['Stimulus::Capsaicin'].Brief.Status=='Did not respond' && result.callsByWindow['Control::Ionomycin'].Brief.Status=='Responded'
    assert result.callsByWindow['Stimulus::Capsaicin'].InvalidF0.Status=='Excluded' && result.callsByWindow['Stimulus::Capsaicin'].Missing.Status=='Unscorable'
    assert result.frameTimes==times && result.baseline==[0,1.5] && state.frameTimes==times
    assert result.protocol.field_overrides[0].Field=='Golden field' && result.appliedProtocol.field_overrides==null
    assert result.manifest.protocol_schema_version=='2' && result.manifest.simpleca_commit=='b92d4ce98ddab3fceac352179703c1f48f044c90'
    assert result.manifest.changed_in_this_plugin=='The whole-trace noise estimator larger_of_baseline_sd_and_trace_noise and response engine 2.1.0 were added in this plugin and are not in SimpleCa'
    // The golden protocol is a sample_sd rule: engine 2.1.0 scores it as 2.0.0 did, and the record says which noise was used.
    assert result.manifest.sd_estimator=='sample_sd' && result.manifest.protocol_engine_version=='2.0.0' && result.manifest.engine_version=='2.1.0'
    assert result.manifest.noise_rule=='Noise is the sample SD of the local prestimulus baseline.' && result.protocol.engine_version=='2.0.0'
    List goldenHeader=result.tables.calls.header;assert !goldenHeader.contains('Noise_source') && goldenHeader.takeRight(2)==['Condition','Biological_unit']
    Map goldenDiagnostics=[sensitivityCalls:result.tables.sensitivityCalls.header,negativeControlCalls:result.tables.negativeControlCalls.header]
    goldenDiagnostics.each{name,header->assert header.size()>40 && !header.contains('Noise_source')}
    // The page shows the estimator nowhere, so the finish message says when the noise was not what a new protocol takes.
    assert state.message.startsWith('Complete.') && state.message.endsWith(' Noise: baseline SD, as saved.') : state.message
    def Calls=loader.loadClass('CalciumCalls')
    Session.records(result.tables.calls).each{row->assert Calls.thresholdSetBy(row)==row.Threshold_source}   // a table without Noise_source reads as before
    assert result.manifest.time_source=='explicit_timestamps' && result.manifest.measurement_input_sha256==data.inputHash
    assert result.manifest.draft_sha256==Settings.fingerprint(draft) && result.appliedFingerprint==Settings.fingerprint(draft)
    assert new File(session.results,'measured_fluorescence.csv').bytes==input.bytes
    assert new File(session.results,'timestamps.csv').getText('UTF-8')==draft.timestamps.content
    assert result.tables.registry.rows.size()==7 && result.tables.cohort.rows[0][result.tables.cohort.header.indexOf('Premeasurement_excluded')]=='1'
    assert result.tables.quality.rows.size()==6 && result.tables.replicates.rows.size()>0
    assert result.tables.sensitivityCounts.rows.size()>0 && result.tables.negativeControlCounts.rows.size()>0
    assert new File(session.results,'sensitivity_grid.csv').getText('UTF-8')==draft.sensitivity.content
    assert new File(session.results,'negative_control_windows.csv').getText('UTF-8')==draft.negative_control.content
    // Pending values/files must not move an applied clock, normalization picture or displayed percentages.
    Map edited=Settings.copy(draft);edited.values.baseline_start='2';edited.timestamps.content='Time_s\n'+times.collect{it+100}.join('\n')+'\n';edited.protocol.windows[0].Floor_local_fraction=0.9d
    assert call('api/draft',[draft:edited,revision:2]).code==200
    assert call('api/state').body.settingsPending && call('api/state').body.frameTimes==times && call('api/state').body.rest==[1,4]
    assert call('api/results').body.counts.Percent_responding=='25'
    assert session.frame(16,'raw').length>0 && recording.getRoi()==null
    // Explicit legacy uses the true schema-1 classifier, which counts the two-sample Brief response.
    def legacyRows=protocol.windows.collect{new LinkedHashMap(it)}
    legacyRows.each{it.Duration_mode='samples';it.Min_duration_s=1d;it.Consecutive_samples=2;it.Min_baseline_samples=3;it.Noise_floor_local_fraction=0d;it.Expected_interval_s=null}
    legacyRows[1].Stimulus='KCl'
    Map legacy=Settings.copy(draft);legacy.protocol=Protocol.payload(legacyRows,1,'Legacy consecutive-sample protocol','user_defined');legacy.sensitivity=null;legacy.negative_control=null
    assert call('api/analyze',[draft:legacy,confirmed:true]).code==200;waitForRun()
    result=call('api/results').body
    assert result.manifest.protocol_schema_version=='1' && result.protocol.engine_version=='legacy-1'
    assert result.callsByWindow['Stimulus::Capsaicin'].Brief.Status=='Responded' && result.counts.Percent_responding=='50'
    assert result.manifest.sd_estimator=='sample_sd' && result.manifest.protocol_engine_version=='legacy-1' && !result.tables.calls.header.contains('Noise_source')
    assert new File(session.results,'settings.tsv').readLines('UTF-8').contains('sd_estimator\tsample_sd') : 'settings.tsv names the legacy rule\'s one noise too'
    assert !call('api/state').body.message.contains('Noise') : 'the page labels the legacy rule itself'
    // A protocol made by the plugin takes noise from the whole trace when that is larger than the baseline SD. Its three
    // columns come after every column of the 2.0.0 table, which keep their places, and the record names the rule.
    Map made=Settings.copy(draft);made.protocol=Protocol.payload(protocol.windows,2,protocol.name,protocol.validation_status,overrides)
    assert made.protocol.sd_estimator==Protocol.newEstimator() && made.protocol.engine_version=='2.1.0'
    assert call('api/analyze',[draft:made,confirmed:true]).code==200;waitForRun()
    result=call('api/results').body
    assert result.manifest.sd_estimator==Protocol.newEstimator() && result.manifest.protocol_engine_version=='2.1.0' && result.manifest.engine_version=='2.1.0'
    assert result.manifest.noise_rule.startsWith('Noise is the larger of the sample SD of the local prestimulus baseline and the whole-trace noise')
    assert result.protocol.engine_version=='2.1.0' && result.protocol.windows*.SD_estimator.unique()==[Protocol.newEstimator()]
    List head=result.tables.calls.header,added=['Whole_trace_noise','Whole_trace_noise_samples','Noise_source']
    assert head.take(goldenHeader.size())==goldenHeader && head.drop(goldenHeader.size())==added : head.takeRight(6)
    goldenDiagnostics.each{name,header->List now=result.tables[name].header;assert now.take(header.size())==header && now.drop(header.size())==added : name+': '+now.takeRight(6)}
    assert new File(session.results,'settings.tsv').readLines('UTF-8').contains('sd_estimator\t'+Protocol.newEstimator()) && !call('api/state').body.message.contains('Noise')
    List scored=Session.records(result.tables.calls).findAll{it.Status in ['Responded','Did not respond']};assert scored.size()>=8
    scored.each{row->double sd=Csv.number(row.Baseline_SD),whole=Csv.number(row.Whole_trace_noise),used=Csv.number(row.Noise_SD_used)
        assert Double.isFinite(whole) && used==Math.max(sd,whole) && row.Noise_source==(whole>sd?'Whole-trace noise':'Baseline SD') : row.Cell+' '+row.Role
        // Threshold_source says "Baseline noise" for the noise term whichever estimate it used; the two columns together say what set the threshold
        assert row.Threshold_source in ['Floor','Baseline noise'] && Calls.thresholdSetBy(row)==(row.Threshold_source=='Floor'?'Floor':row.Noise_source) : row.Cell+' '+row.Role}
    assert Settings.tally(session.results,data).rows.find{it.key=='all'}.counted[2]==6 : 'the Fiji tally reads the longer table'
    // The measurement hash protects stale session data; a failed source check leaves the previous completed run intact.
    File applied=session.results;int version=session.version
    input.append('\n','UTF-8')
    assert call('api/analyze',[draft:legacy,confirmed:true]).code==200
    long end=System.nanoTime()+30000000000L;while(session.busy&&System.nanoTime()<end)Thread.sleep(25)
    assert !session.busy && session.version==version && session.results==applied && call('api/state').body.message.contains('measurement file changed')
    println 'PASS: real HTTP web-session -> pinned R golden calls/counts, versioned legacy behavior, explicit timing, registry/QC, diagnostics/replicates, immutable snapshots and stale-source protection'
}finally{server.stop();recording.close();picture.close()}

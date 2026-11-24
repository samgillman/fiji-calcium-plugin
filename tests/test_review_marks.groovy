// Manual review marks on the calls of an analysis, without the page: set, remove and refuse over HTTP as the Step 4 page
// would; the tally of every group as scored and after review; the percent table with and without marks; the sentence for a
// Methods section; marks carried over to a new analysis of the same recording and the ones dropped; a review file that does
// not match its analysis. Synthetic analysis folders only; no R, no display, no browser.
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Review=loader.loadClass('CalciumReview'),Settings=loader.loadClass('CalciumAnalysisSettings'),Csv=loader.loadClass('CalciumCsv'),Installer=loader.loadClass('CalciumInstaller')
def Session=loader.loadClass('CalciumAnalysisSession'),Server=loader.loadClass('CalciumWebServer')
def slurp=new JsonSlurper()
def root=new File('build/test-output/review-marks').absoluteFile;root.deleteDir();root.mkdirs()
File folder=new File(root,'rec');folder.mkdirs()
int N=10
List names=(1..8).collect{'Cell_'+String.format('%032x',it)}
def id={int k->names[k-1]}
// cell: call, required rise, chosen blind; cell 5 goes below zero; cell 2 also has a control window
Map cells=[1:['Responded',0.20,true],2:['Responded',0.50,true],3:['Did not respond',0.20,true],4:['Did not respond',0.40,false],
           5:['Responded',0.10,false],6:['Excluded','',true],7:['Unscorable','',false],8:['Did not respond',0.25,false]]
Map chosen=cells.collectEntries{k,v->[id(k),v[2]]}
// the analysis folder a run of the engine would leave, with the columns the tally, the review and the page read
def analysis={String name,Map calls,Map outlines=[:]->
    def dir=new File(folder,name);dir.mkdirs()
    List rows=[['Field','Cell','Role','Stimulus','Status','Reason','Threshold_increase']]
    calls.each{k,v->rows<<['rec',id(k as int),'Stimulus','ATP',v[0],'',v[1]];if(k==2)rows<<['rec',id(2),'Control','KCl','Did not respond','',0.3]}
    Csv.write(new File(dir,'responder_calls.csv'),rows)
    List present=calls.keySet().collect{id(it as int)}
    Csv.write(new File(dir,'normalized_traces.csv'),[['Time']+present]+(0..<N).collect{int f->[f]+present.collect{it==id(5)&&f==6?-1.5:0.1*f}})
    Csv.write(new File(dir,'mean_trace.csv'),[['Time','Mean','SEM','N_Cells']]+(0..<N).collect{int f->[f,0.1*f,0,present.size()]})
    Csv.write(new File(dir,'analysis_manifest.csv'),[['key','value'],['baseline_start','1'],['baseline_end','3'],['sampling_hz','1']])
    Csv.write(new File(dir,'roi_registry.csv'),[['Cell','Included_in_measurement','ROI_sha256']]+present.collect{[it,'true',outlines[it]?:'roi-of-'+it]})
    new File(dir,'COMPLETE').text='complete\n'
    dir
}
def snapshot={File dir->dir.listFiles().sort{it.name}.collectEntries{[it.name,Installer.digest(it.bytes)]}}
def measurement={List cellNames->[folder:folder,recording:'rec',frames:N,hz:1d,baselineStart:1,baselineEnd:3,key:'rec',measurementId:'rec',names:cellNames,chosenBlind:chosen.findAll{it.key in cellNames}]}
def open={Map data->
    def session=Session.newInstance(data);session.runtime.home=new File(root,'runtime')
    def server=Server.newInstance(session:session);server.start()
    def call={String path,Object body=null->def c=(HttpURLConnection)new URL(server.url()+path).openConnection(java.net.Proxy.NO_PROXY);c.connectTimeout=5000;c.readTimeout=30000
        if(body!=null){c.requestMethod='POST';c.doOutput=true;c.setRequestProperty('Content-Type','application/json');c.outputStream.withCloseable{it.write(JsonOutput.toJson(body).getBytes('UTF-8'))}}
        int status=c.responseCode;[status,slurp.parseText(new String((status<400?c.inputStream:c.errorStream).bytes,'UTF-8'))]}
    def run={File made->
        session.performAnalysis={rt,input,settings,out,job->made}
        def started=call('api/analyze',[draft:call('api/session')[1].draft,confirmed:true]);assert started[0]==200 : started
        long end=System.nanoTime()+60000000000L;while(session.busy&&System.nanoTime()<end)Thread.sleep(20)
        def state=call('api/state')[1];assert !state.busy && state.tone=='good' : state.message
        [state:state,results:call('api/results')[1]]}
    // a body sent as it is, for requests the page never makes
    def raw={String path,String text->def c=(HttpURLConnection)new URL(server.url()+path).openConnection(java.net.Proxy.NO_PROXY);c.connectTimeout=5000;c.readTimeout=30000
        c.requestMethod='POST';c.doOutput=true;c.outputStream.withCloseable{it.write(text.getBytes('UTF-8'))}
        int status=c.responseCode;[status,slurp.parseText(new String((status<400?c.inputStream:c.errorStream).bytes,'UTF-8'))]}
    [session:session,server:server,call:call,run:run,raw:raw]
}
def percentName={String made->'rec_percent_responding_'+made+'.csv'}
def reviewName={String made->'rec_reviewed_calls_'+made+'.csv'}

// ---- one analysis with no marks: the review is empty and nothing is written for it
String A='SimpleCa-20261101-100000-aaaaaaaa'
File dirA=analysis(A,cells)
Map folderA=snapshot(dirA)
def page=open(measurement(names))
def first=page.run(dirA)
assert first.results.review.findAll{it.key!='after'}==[folder:A,file:null,sha256:null,problem:'',marks:[:],changed:0,methods:'',carried:null,reasons:['Noise','Movement','Passing object','Bad outline','Other'],calls:['Responded','Did not respond'],noteMax:120] : first.results.review
assert first.results.review.after==cells.collectEntries{k,v->[id(k),v[0]]} : 'the call after review of every cell, here the call as scored'
assert first.state.cancellable==false && first.state.busy==false
assert !new File(folder,reviewName(A)).exists() && !first.results.tally.containsKey('reviewed') && first.results.tally.rows.every{!it.containsKey('reviewed')} : 'no marks: the tally is as it always was'
assert first.state.message.endsWith('. All cells: 3 of 6 responded (50.0 %). Chosen blind: 2 of 3 responded (66.7 %). Added after the response was shown: 1 of 3 responded (33.3 %).') : first.state.message
byte[] percentBefore=new File(folder,percentName(A)).bytes
assert new String(percentBefore,'UTF-8')=='cells,registered,evaluable,responded,percent_responding,cells_below_zero,evaluable_without_them,responded_without_them,percent_without_them\n'+
    'All cells,8,6,3,50.0,1,5,2,40.0\nChosen blind,4,3,2,66.7,0,3,2,66.7\nAdded after the response was shown,4,3,1,33.3,1,2,0,0.0\nQuiet enough to detect a 30 % rise,4,4,2,50.0,1,3,1,33.3\n'
println 'PASS: an analysis with no marks has an empty review, no review file and the tally and percent table as before'

// ---- marks both ways, over HTTP
def mark={String cell,String to,String reason,String note=null->def body=[cell:cell,call:to,reason:reason,folder:A];if(note!=null)body.note=note;page.call('api/review',body)}
def set=mark(id(2),'Did not respond','Noise')
assert set[0]==200 && set[1].folder==A && set[1].review.marks.keySet()==[id(2)] as Set && set[1].review.marks[id(2)].name=='Cell 2' : set
assert set[1].review.marks[id(2)].subMap(['automatic','reviewed','reason','note','carriedFrom'])==[automatic:'Responded',reviewed:'Did not respond',reason:'Noise',note:'',carriedFrom:''] && set[1].review.marks[id(2)].time==~/\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(Z|[+-]\d\d:\d\d)/
assert Review.TIME.format(java.time.OffsetDateTime.of(2026,11,23,18,2,0,0,java.time.ZoneOffset.ofHours(-6)))=='2026-11-23T18:02:00-06:00' : 'the time is always given to the second'
assert set[1].review.after[id(2)]=='Did not respond' && set[1].review.after[id(1)]=='Responded' && set[1].review.after[id(6)]=='Excluded'
assert mark(id(4),'Responded','Other','slow ramp, seen "by eye"')[0]==200
def third=mark(id(5),'Did not respond','Passing object')
assert third[0]==200 && third[1].review.changed==3 && third[1].version>first.state.version
// the status line follows the marks: the figures after review, worked out again from the tally
String afterThree=first.state.message.replace('. All cells: 3 of 6 responded (50.0 %). Chosen blind: 2 of 3 responded (66.7 %). Added after the response was shown: 1 of 3 responded (33.3 %).',
    '. All cells: 3 of 6 responded (50.0 %), after review 2 of 6 (33.3 %). Chosen blind: 2 of 3 responded (66.7 %), after review 1 of 3 (33.3 %). Added after the response was shown: 1 of 3 responded (33.3 %), after review 1 of 3 (33.3 %).')
assert third[1].message==afterThree && page.call('api/state')[1].message==afterThree : page.call('api/state')[1].message
assert third[1].review.after==cells.collectEntries{k,v->[id(k),k in [2,5] ? 'Did not respond' : k==4 ? 'Responded' : v[0]]}
assert third[1].review.sha256==Installer.digest(new File(folder,reviewName(A)).bytes) && third[1].review.file==reviewName(A)
def rowOf={Map tally,String key->tally.rows.find{it.key==key}}
Map tally=third[1].tally
// all: as scored 1, 2, 5 responded of 1, 2, 3, 4, 5, 8; after review 1 and 4. Cell 5 went below zero, so "without" leaves it out.
assert tally.reviewed==true && tally.changed==[id(2),id(4),id(5)]
assert rowOf(tally,'all').counted==[3,6,8] && rowOf(tally,'all').without==[2,5,7] && rowOf(tally,'all').reviewed==[counted:[2,6,8],without:[2,5,7],changed:3]
assert rowOf(tally,'blind').counted==[2,3,4] && rowOf(tally,'blind').reviewed==[counted:[1,3,4],without:[1,3,4],changed:1]
assert rowOf(tally,'later').counted==[1,3,4] && rowOf(tally,'later').reviewed==[counted:[1,3,4],without:[1,2,3],changed:2]
assert rowOf(tally,'quiet').counted==[2,4,4] && rowOf(tally,'quiet').reviewed==[counted:[1,4,4],without:[1,3,3],changed:1] : 'quiet: cells 1, 3, 5, 8; the review moves cell 5'
assert third[1].review.methods=='3 of 6 cells were reclassified after manual review: 2 from responder to non-responder (1 noise, 1 passing object) and 1 from non-responder to responder (other).' : third[1].review.methods
// the saved results the page reads say the same, and nothing else in them changed
def reread=page.call('api/results')[1]
assert reread.review==third[1].review && reread.tally==tally && reread.calls==first.results.calls && reread.traces==first.results.traces && reread.checks==first.results.checks
assert page.call('api/state')[1].version==third[1].version
// the review file: one row per mark in the cells' order, with the analysis, the hash of its calls and the outline
def saved=Csv.read(new File(folder,reviewName(A)))
assert saved[0]==['Cell','Automatic_call','Reviewed_call','Reason','Note','Time','Analysis','Responder_calls_sha256','ROI_sha256','Carried_from'] && saved.drop(1).every{it[9]==''}
assert saved.drop(1).collect{it[0..4]+it[6..8]}==[[id(2),'Responded','Did not respond','Noise',''],[id(4),'Did not respond','Responded','Other','slow ramp, seen "by eye"'],[id(5),'Responded','Did not respond','Passing object','']].collect{it+[A,Installer.digest(new File(dirA,'responder_calls.csv').bytes),'roi-of-'+it[0]]}
assert new File(folder,reviewName(A)).text.contains('"slow ramp, seen ""by eye"""') : 'a note with a comma and quotes is quoted'
assert snapshot(dirA)==folderA : 'the analysis folder, responder_calls.csv with it, is not written to'
// the percent table: the same rows and columns, then the figures after review
List percent=new File(folder,percentName(A)).readLines('UTF-8')
assert percent[0]=='cells,registered,evaluable,responded,percent_responding,cells_below_zero,evaluable_without_them,responded_without_them,percent_without_them,responded_after_review,percent_responding_after_review,responded_without_them_after_review,percent_without_them_after_review,cells_reclassified'
assert percent.drop(1)==['All cells,8,6,3,50.0,1,5,2,40.0,2,33.3,2,40.0,3','Chosen blind,4,3,2,66.7,0,3,2,66.7,1,33.3,1,33.3,1','Added after the response was shown,4,3,1,33.3,1,2,0,0.0,1,33.3,1,50.0,2','Quiet enough to detect a 30 % rise,4,4,2,50.0,1,3,1,33.3,1,25.0,1,33.3,1'] : percent
(0..<5).each{int r->assert percent[r].startsWith(new String(percentBefore,'UTF-8').split('\n')[r]+(r==0?',responded':',')) : 'appended on the right: '+percent[r]}
// the finish sentence gives both figures; so does the one place Settings.tally, read straight from the folder
assert Settings.percentResponding(dirA,page.session.measurement)==' All cells: 3 of 6 responded (50.0 %), after review 2 of 6 (33.3 %). Chosen blind: 2 of 3 responded (66.7 %), after review 1 of 3 (33.3 %). Added after the response was shown: 1 of 3 responded (33.3 %), after review 1 of 3 (33.3 %).'
assert slurp.parseText(JsonOutput.toJson(Settings.tally(dirA,page.session.measurement)))==tally
assert Settings.percentResponding(dirA,page.session.measurement+[chosenBlind:null])==' All cells: 3 of 6 responded (50.0 %), after review 2 of 6 (33.3 %).' : 'with marks all cells are given even when it was not recorded how cells were chosen'
// the Fiji window's tile, from the same tally
assert Review.headline(tally,'ATP')==['33.3 %','ATP: 2 of 6 evaluable cells after review; 3 as scored'] && Review.headline(first.results.tally,'ATP')==null
println 'PASS: marks both ways change the tally of every group after review and leave it as scored; the percent table gains columns on the right; the Methods sentence counts the marks'

// ---- refused, with a plain message, and nothing changes
byte[] reviewBefore=new File(folder,reviewName(A)).bytes
def refused={Map body,String message->def answer=page.call('api/review',body);assert answer[0]==400 && answer[1].error==message : answer;assert new File(folder,reviewName(A)).bytes==reviewBefore}
refused([cell:id(6),call:'Responded',reason:'Noise'],'Cell 6 has no call to review: it is excluded.')
refused([cell:id(7),call:'Responded',reason:'Noise'],'Cell 7 has no call to review: it is unscorable.')
refused([cell:id(1),call:'Responded',reason:'Noise'],'Cell 1 was scored as responded. The reviewed call must be the other one.')
refused([cell:id(1),call:'Maybe',reason:'Noise'],'Give the reviewed call: Responded or Did not respond.')
refused([cell:id(1),call:'Did not respond',reason:'Tired'],'Choose a reason: Noise, Movement, Passing object, Bad outline or Other.')
refused([cell:id(1),call:'Did not respond'],'Choose a reason: Noise, Movement, Passing object, Bad outline or Other.')
refused([cell:id(1),call:'Did not respond',reason:'Other'],'Other needs a short note.')
refused([cell:id(1),call:'Did not respond',reason:'Other',note:'   '],'Other needs a short note.')
refused([cell:id(1),call:'Did not respond',reason:'Other',note:'two\nlines'],'The note must be one line.')
refused([cell:id(1),call:'Did not respond',reason:'Noise',note:'tab\there'],'The note must be one line.')
refused([cell:id(1),call:'Did not respond',reason:'Other',note:'x'*121],'Keep the note to 120 characters or fewer.')
refused([cell:id(1),call:'Did not respond',reason:'Noise',note:42],'The note must be text.')
refused([cell:'Cell_nobody',call:'Did not respond',reason:'Noise'],'There is no cell Cell_nobody in this analysis.')
refused([cell:7,call:'Did not respond',reason:'Noise'],'Say which cell to mark.')
refused([cell:id(1),remove:true],'Cell 1 has no review mark.')
refused([cell:id(1),call:'Did not respond',reason:'Noise',folder:'SimpleCa-another'],'This mark is for another analysis. Reload the page.')
['=SUM(1)','+1 slow','-ve start','@cmd','  =padded'].each{refused([cell:id(1),call:'Did not respond',reason:'Noise',note:it],'The note cannot start with =, +, - or @.')}
['\u200B','abc\u202Edcba','\uFEFFnote'].each{refused([cell:id(1),call:'Did not respond',reason:'Other',note:it],'The note must be plain text, without hidden characters.')}
refused([cell:id(1),call:'Did not respond',reason:'Other',note:'\u00A0\u2003'],'Other needs a short note.')
['{not json','[1,2]','"text"','42'].each{def answer=page.raw('api/review',it);assert answer[0]==400 && answer[1].error=='The request is not valid.' : answer}
assert new File(folder,reviewName(A)).bytes==reviewBefore
assert mark(id(1),'Did not respond','Other','x'*120)[0]==200 && page.call('api/review',[cell:id(1),remove:true])[0]==200 : 'a note of 120 characters is kept'
assert new File(folder,reviewName(A)).bytes==reviewBefore : 'a mark set and removed again leaves the file as it was'
println 'PASS: a mark on a cell without a call, to the same call, without a reason, Other without a note, a note of two lines, too long, hidden or starting like a formula, an unknown cell, another analysis or a body that is not an object is refused with a plain message'

// ---- the percent table cannot be written: nothing is saved, the mark is refused, and the message says what to do
File percentA=new File(folder,percentName(A));byte[] percentNow=percentA.bytes;String resultsNow=page.session.resultsJson
assert percentA.setWritable(false)
def locked=page.call('api/review',[cell:id(1),call:'Did not respond',reason:'Noise'])
percentA.setWritable(true)
assert locked[0]==400 && locked[1].error=='Close '+percentName(A)+' if it is open in another program, then try again.' : locked
assert new File(folder,reviewName(A)).bytes==reviewBefore && percentA.bytes==percentNow && page.session.resultsJson==resultsNow && !folder.list().any{it.endsWith('.part')}
println 'PASS: a mark is saved with its percent table or not at all'

// ---- marks removed: one, then all; the file goes and the percent table is byte for byte as it was without marks
def removed=page.call('api/review',[cell:id(5),remove:true])
assert removed[0]==200 && removed[1].review.marks.keySet()==[id(2),id(4)] as Set && removed[1].review.changed==2 && rowOf(removed[1].tally,'all').reviewed.counted==[3,6,8]
assert mark(id(2),'Did not respond','Movement')[1].review.marks[id(2)].reason=='Movement' : 'a mark set again is replaced'
[id(2),id(4)].each{assert page.call('api/review',[cell:it,remove:true])[0]==200}
assert !new File(folder,reviewName(A)).exists() && new File(folder,percentName(A)).bytes==percentBefore : 'no marks left: no review file, and the percent table as before'
def cleared=page.call('api/results')[1]
assert cleared.tally==first.results.tally && cleared.review.marks==[:] && cleared.review.file==null && cleared.review.methods=='' && cleared.review.sha256==null
assert page.call('api/state')[1].message==first.state.message : 'the status line is as it was without marks'
// the review file cannot take its place: the percent table is put back as it was
File inTheWay=new File(folder,reviewName(A));inTheWay.mkdirs();new File(inTheWay,'keep').text='x'
def blockedReview=page.call('api/review',[cell:id(1),call:'Did not respond',reason:'Noise'])
assert blockedReview[0]==400 && blockedReview[1].error=='Close '+reviewName(A)+' if it is open in another program, then try again.' && new File(folder,percentName(A)).bytes==percentBefore && !folder.list().any{it.endsWith('.part')}
new File(inTheWay,'keep').delete();inTheWay.delete()
// a percent table left with after-review columns by a review file taken away by hand is put right when the results are read
assert mark(id(2),'Did not respond','Noise')[0]==200 && new File(folder,percentName(A)).readLines('UTF-8')[0].endsWith(',cells_reclassified')
assert new File(folder,reviewName(A)).delete()
page.session.load(dirA)
assert new File(folder,percentName(A)).bytes==percentBefore && page.call('api/results')[1].review.marks==[:]
page.session.load(dirA);assert new File(folder,percentName(A)).bytes==percentBefore : 'a table that is right is left alone'
println 'PASS: marks are removed one by one; with none left the review file is gone and the percent table is byte for byte the one without marks'

// ---- carried over to a new analysis of the same recording; the marks whose cell changed are dropped and listed
[[2,'Did not respond','Noise',null],[3,'Responded','Movement',null],[4,'Responded','Other','slow ramp'],[8,'Responded','Bad outline',null]].each{assert mark(id(it[0]),it[1],it[2],it[3])[0]==200}
String timeOf2=Csv.read(new File(folder,reviewName(A))).find{it[0]==id(2)}[5]
page.server.stop()
// measured again without cell 3; cell 8 has a new outline; cell 4 is now scored as responding
String B='SimpleCa-20261101-110000-bbbbbbbb'
Map cellsB=cells.findAll{it.key!=3}.collectEntries{k,v->[k,k==4 ? ['Responded',0.40,false] : v]}
File dirB=analysis(B,cellsB,[(id(8)):'roi-new'])
List namesB=names-[id(3)]
page=open(measurement(namesB))
def second=page.run(dirB)
Map carried=second.results.review.carried
String carriedText='Review marks from '+A+': 1 carried over, 3 dropped. Cell 3 of the earlier analysis: no longer measured. Cell 3: now scored as responded. Cell 7: outline changed.'
assert carried==[from:A,kept:1,dropped:[[cell:id(3),name:'Cell 3 of the earlier analysis',why:'no longer measured',text:'Cell 3 of the earlier analysis: no longer measured'],[cell:id(4),name:'Cell 3',why:'now scored as responded',text:'Cell 3: now scored as responded'],[cell:id(8),name:'Cell 7',why:'outline changed',text:'Cell 7: outline changed']],
    text:carriedText] : carried
assert second.state.message.endsWith(' '+carriedText) && second.state.message.contains('All cells: 4 of 5 responded (80.0 %), after review 3 of 5 (60.0 %). Chosen blind: 2 of 2 responded (100.0 %), after review 1 of 2 (50.0 %).') : second.state.message
def savedB=Csv.read(new File(folder,reviewName(B)))
assert savedB.size()==2 && savedB[1][0..5]==[id(2),'Responded','Did not respond','Noise','',timeOf2] && savedB[1][6..9]==[B,Installer.digest(new File(dirB,'responder_calls.csv').bytes),'roi-of-'+id(2),A] : 'the kept mark, made for the new analysis, keeps its own time and says where it was made'
assert second.results.review.marks[id(2)].carriedFrom==A
assert second.results.review.marks.keySet()==[id(2)] as Set && second.results.review.methods=='1 of 5 cells was reclassified after manual review: 1 from responder to non-responder (noise).'
assert new File(folder,percentName(B)).readLines('UTF-8')[0].endsWith(',cells_reclassified') && new File(folder,reviewName(A)).isFile() : 'the earlier analysis keeps its own review'
// the carried mark is removed: the status line no longer speaks of it
def removedB=page.call('api/review',[cell:id(2),remove:true])
assert removedB[0]==200 && removedB[1].review.carried==null && !removedB[1].message.contains('Review marks') && removedB[1].message.contains('All cells: 4 of 5 responded (80.0 %).') : removedB[1].message
// the next analysis takes the marks of the most recent earlier one, B, which has none now: the marks of A do not come back
String C='SimpleCa-20261101-120000-cccccccc'
File dirC=analysis(C,cellsB)
def third2=page.run(dirC)
assert third2.results.review.carried==null && !third2.state.message.contains('Review marks') && !new File(folder,reviewName(C)).exists() && new File(folder,reviewName(A)).isFile()
// the most recent earlier review file does not match its analysis: said, and no older one is used instead
assert page.call('api/review',[cell:id(2),call:'Did not respond',reason:'Noise'])[0]==200
String Z='SimpleCa-20261101-120001-zzzzzzzz'
File dirZ=analysis(Z,cellsB);new File(folder,reviewName(Z)).text=new File(folder,reviewName(C)).text.replace(C,Z)   // more recent than C, but with C's calls
new File(dirZ,'responder_calls.csv').text=new File(dirZ,'responder_calls.csv').text+'rec,Cell_extra,Stimulus,ATP,Responded,,0.2\n'
String E='SimpleCa-20261101-121000-eeeeeeee'
def fifth=page.run(analysis(E,cellsB))
assert fifth.results.review.carried==[from:Z,kept:0,dropped:[],text:'The review marks of '+Z+' were not carried over: its review file does not match that analysis.'] && fifth.state.message.endsWith(' The review marks of '+Z+' were not carried over: its review file does not match that analysis.') : fifth.results.review.carried
assert !new File(folder,reviewName(E)).exists()
// every mark dropped at a run: the next run does not bring them back either
assert page.call('api/review',[cell:id(2),call:'Did not respond',reason:'Noise'])[0]==200
String F='SimpleCa-20261101-122000-ffffffff',G='SimpleCa-20261101-123000-gggggggg'
def sixthRun=page.run(analysis(F,cellsB.collectEntries{k,v->[k,k==2 ? ['Did not respond',0.50,true] : v]}))
assert sixthRun.results.review.carried.kept==0 && sixthRun.results.review.carried.dropped*.why==['now scored as did not respond'] && !new File(folder,reviewName(F)).exists()
def seventhRun=page.run(analysis(G,cellsB))
assert seventhRun.results.review.carried==null && !new File(folder,reviewName(G)).exists() : 'the marks of E are not carried past F'
// a run with no earlier marks says nothing about them
folder.listFiles().findAll{it.name.startsWith('rec_reviewed_calls_')}.each{it.delete()}
File dirD=analysis('SimpleCa-20261101-130000-dddddddd',cellsB)
def fourth=page.run(dirD)
assert fourth.results.review.carried==null && !fourth.state.message.contains('Review marks') && !new File(folder,reviewName('SimpleCa-20261101-130000-dddddddd')).exists()
println 'PASS: at a new analysis the marks of the most recent earlier analysis are carried over where the call and outline are unchanged; the others are dropped and listed; removed and dropped marks do not come back'

// ---- a review file that does not match its analysis is not applied, and says so
assert page.call('api/review',[cell:id(2),call:'Did not respond',reason:'Noise'])[0]==200
File fileD=new File(folder,reviewName('SimpleCa-20261101-130000-dddddddd'));String good=fileD.text
fileD.text=good.replace(Installer.digest(new File(dirD,'responder_calls.csv').bytes),'0'*64)
page.session.load(dirD);def mismatch=page.call('api/results')[1]
assert mismatch.review.problem=='The review file '+fileD.name+' does not match the calls of this analysis, so its marks are not used.' && mismatch.review.marks==[:] && !mismatch.tally.containsKey('reviewed') : mismatch.review
def blocked=page.call('api/review',[cell:id(1),call:'Did not respond',reason:'Noise'])
assert blocked[0]==400 && blocked[1].error=='The review file '+fileD.name+' does not match the calls of this analysis, so its marks are not used. Move it out of the folder to review this analysis again.'
fileD.text=good.replace('SimpleCa-20261101-130000-dddddddd','SimpleCa-elsewhere')
page.session.load(dirD);assert page.call('api/results')[1].review.problem=='The review file '+fileD.name+' is for another analysis, so its marks are not used.'
fileD.text=good.replace(',Noise,',',Tired,')
page.session.load(dirD);assert page.call('api/results')[1].review.problem.endsWith('does not match the calls of this analysis, so its marks are not used.')
File percentD=new File(folder,percentName('SimpleCa-20261101-130000-dddddddd'))
assert !percentD.readLines('UTF-8')[0].contains('after_review') : 'a review file that is not used leaves no after-review columns in the percent table'
fileD.text='not,a,review\n'
page.session.load(dirD);assert page.call('api/results')[1].review.problem=='The review file '+fileD.name+' cannot be read, so its marks are not used.'
fileD.text=good;page.session.load(dirD);def back=page.call('api/results')[1]
assert back.review.problem=='' && back.review.marks.keySet()==[id(2)] as Set && back.tally.reviewed==true && percentD.readLines('UTF-8')[0].endsWith(',cells_reclassified')
// a file saved by a spreadsheet program, with a byte order mark and blank lines at the end, is read; an empty one holds no marks
fileD.setText('\uFEFF'+good.replace('\n','\r\n')+'\r\n\r\n','UTF-8');page.session.load(dirD);def bom=page.call('api/results')[1]
assert bom.review.problem=='' && bom.review.marks.keySet()==[id(2)] as Set : bom.review
fileD.text='';page.session.load(dirD);def empty=page.call('api/results')[1]
assert empty.review.problem=='' && empty.review.marks==[:] && !empty.tally.containsKey('reviewed')
fileD.text=good;page.session.load(dirD)
println 'PASS: a review file with another hash, another analysis, a mark that does not fit or another layout is not applied and the answer says why'

// ---- once the engine has finished, the checks and the carry-over run to the end: cancel does nothing then, and says so
def job=loader.loadClass('CalciumLocalJob').newInstance(logFile:new File(root,'job.log'))
page.session.job=job;page.session.busy=true;page.session.checking=true
assert page.call('api/state')[1].cancellable==false && page.call('api/cancel',[:])[0]==200 && !job.cancelled && page.session.stage!='Cancelling...'
page.session.checking=false
assert page.call('api/state')[1].cancellable==true;page.session.cancel();assert job.cancelled
page.session.busy=false
println 'PASS: cancel stops the engine only; the checks after it are not cancelled halfway'
page.server.stop()
System.exit(0)

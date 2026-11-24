// The Step 4 window of Fiji with the checks of the calls and the review marks: the flags of the selected cell under its call
// line, each in the checks' own words, and a flag's link to the first frame of its event; the mark button on a cell with a
// call, its form in place with the reason already chosen from a flag, Enter, Esc and the refusals in the form; a mark saved
// and removed, with the tile, the Show list, the picture, the heatmap, the plot, the call line, the finish message and the
// record of saved plots following it, and the time cursor staying where it is; no mark controls on excluded or unscorable
// cells, on a control window or while an analysis runs; a mark carried over to the next analysis run here.
// The window is built and never shown (with -Dapple.awt.UIElement=true no Dock icon appears on a Mac); a window that opens all
// the same fails the test. -Dcalcium.test.pictures=<folder> saves pictures of the Traces tab there.
import javax.swing.*
import ij.*
import ij.gui.*
import ij.process.ColorProcessor
import java.awt.event.KeyEvent

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Window=loader.loadClass('CalciumAnalysisWindow'),Review=loader.loadClass('CalciumReview'),Checks=loader.loadClass('CalciumCallChecks'),Csv=loader.loadClass('CalciumCsv')
def Settings=loader.loadClass('CalciumAnalysisSettings'),Installer=loader.loadClass('CalciumInstaller'),Protocol=loader.loadClass('CalciumResponseProtocol'),Calls=loader.loadClass('CalciumCalls'),Theme=loader.loadClass('CalciumTheme')
def onUI={Closure c->SwingUtilities.invokeAndWait(c as Runnable)}
def strip={String s->s==null?'':s.replaceAll('<br>',' | ').replaceAll('<[^>]+>','').replace('&nbsp;',' ').replace('&lt;','<').replace('&gt;','>').replace('&amp;','&').trim()}

// Nothing may appear on screen. A window that opens all the same is closed at once and fails the test.
def appeared=Collections.synchronizedList([])
java.awt.Toolkit.defaultToolkit.addAWTEventListener({event->if(event.ID==java.awt.event.WindowEvent.WINDOW_OPENED){appeared<<(''+event.window);event.window.dispose()}} as java.awt.event.AWTEventListener,java.awt.AWTEvent.WINDOW_EVENT_MASK)
ij.macro.Interpreter.batchMode=true

// A measurement folder of its own, with a recording name of its own: no file another test leaves is carried over here.
def work=new File('build/analysis-review-gui').absoluteFile;work.deleteDir();work.mkdirs()
def pictures=new File(System.getProperty('calcium.test.pictures')?:new File(work,'pictures').path).absoluteFile;pictures.mkdirs()
String A='SimpleCa-20261101-100000-aaaaaaaa',B='SimpleCa-20261101-110000-bbbbbbbb'
def copy={String name->def dir=new File(work,name);dir.mkdirs();new File('tests/fixtures/simpleca').listFiles().findAll{it.isFile()}.each{new File(dir,it.name).bytes=it.bytes};dir}
File dirA=copy(A)
def snapshot={->work.listFiles().findAll{it.isFile()}.sort{it.name}.collectEntries{[it.name,Installer.digest(it.bytes)]}}
List names=(1..6).collect{String.format('Cell%03d',it)}
String unchecked='Not checked for a passing object: the raw gray values and backgrounds are not beside the analysis.'

// A key pressed in the window, as Swing hands it on: the focused control's own keys, then those of the controls around it.
def press={JComponent focus,int code->
    def ks=KeyStroke.getKeyStroke(code,0),e=new KeyEvent(focus,KeyEvent.KEY_PRESSED,System.currentTimeMillis(),0,code,KeyEvent.CHAR_UNDEFINED)
    def fire={JComponent c,int when->def key=c.getInputMap(when)?.get(ks);def action=key!=null?c.actionMap.get(key):null;action!=null&&SwingUtilities.notifyAction(action,ks,e,c,0)}
    if(fire(focus,JComponent.WHEN_FOCUSED))return true
    for(def c=focus;c instanceof JComponent;c=c.parent)if(fire(c as JComponent,JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT))return true
    false}
def click={java.awt.Component c->c.dispatchEvent(new java.awt.event.MouseEvent(c,java.awt.event.MouseEvent.MOUSE_CLICKED,System.currentTimeMillis(),0,4,4,1,false,java.awt.event.MouseEvent.BUTTON1))}

def ui,original,recording
def finish={->long deadline=System.nanoTime()+30000000000L;while(ui.busy&&System.nanoTime()<deadline)Thread.sleep(20);onUI{};assert !ui.busy : 'the analysis finishes'}
// The window laid out at its usual size without being shown, and its content painted into a picture. A window that is never
// shown is never validated, so every layout forgets the sizes it kept first (a shown window does this when it validates).
def layOut={->ui.setSize(1180,830);def forget,lay
    forget={java.awt.Component c->c.invalidate();if(c instanceof java.awt.Container)(c as java.awt.Container).components.each{forget(it)}}
    lay={java.awt.Component c->if(c instanceof java.awt.Container){(c as java.awt.Container).doLayout();(c as java.awt.Container).components.each{lay(it)}}}
    forget(ui);lay(ui);forget(ui);lay(ui)}   // twice: wrapped text takes its height from the width the first pass gave it
def paint={String name->layOut();def pane=ui.contentPane as JComponent
    def image=new java.awt.image.BufferedImage(pane.width,pane.height,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics();try{pane.paint(g)}finally{g.dispose()}
    javax.imageio.ImageIO.write(image,'png',new File(pictures,name+'.png'))
    def head=ui.reviewBox.parent as JComponent;def top=SwingUtilities.convertPoint(head,0,0,pane)
    javax.imageio.ImageIO.write(image.getSubimage(Math.max(0,top.x as int),Math.max(0,top.y as int),Math.min(head.width,pane.width-(top.x as int)),Math.min(head.height,pane.height-(top.y as int))),'png',new File(pictures,name+'_head.png'))}
try {
    onUI {
        // the look and feel's own list font is made unlike the window's, so a list left in it shows (on a Mac the two are alike)
        UIManager.put('ComboBox.font',new javax.swing.plaf.FontUIResource(Theme.font().family,java.awt.Font.ITALIC,9))
        original=new ImagePlus('Review colour picture',new ColorProcessor(250,150))
        def frames=new ImageStack(250,150)
        (1..30).each{int f->def ip=new ij.process.ShortProcessor(250,150);ip.set(100);if(f>=10&&f<=20){ip.setRoi(30,50,20,20);ip.set(140);ip.resetRoi()};frames.addSlice(''+f,ip)}
        recording=new ImagePlus('Review recording',frames)
        ui=Window.newInstance([csv:new File(dirA,'metrics.csv'),folder:work,recording:'Review',frames:30,hz:2d,baselineStart:1,baselineEnd:5,
            picture:original,movie:recording,names:names,rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]])
        ui.runtime.home=new File(work,'runtime')
        assert !ui.visible && !ui.markRow.visible && !ui.flagRow.visible : 'before an analysis there is nothing to mark'
        // ---- an analysis run here: the checks it gives are kept, and what they could not do is said once, in the finish message
        ui.confirmed.selected=true
        ['Baseline_start_s':0.5d,'Baseline_end_s':3d,'Stimulus_onset_s':3d,'Start_s':3d,'End_s':14.5d].each{k,v->ui.windowModel.setValueAt(v,0,Protocol.COLUMNS.indexOf(k))}
        ui.performAnalysis={a,b,c,d,job->dirA}
        ui.runButton.doClick();assert ui.busy : 'the analysis starts: '+ui.status.text
    }
    finish()
    onUI {
        assert ui.results==dirA : ui.status.text
        assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+A+'. '+unchecked : 'the finish message says once what the checks could not do: '+ui.status.text
        assert ui.checks.note==unchecked && !ui.checks.ringMeasured : 'the checks of the run'
        ui.showGroup(0);ui.cells.selectedItem='Cell001'
        assert !ui.flagRow.visible : 'no flag, and no "not checked" on each cell when the surroundings were not measured'
        assert ui.markRow.visible && ui.markButton.visible && ui.markButton.text=='Mark as not a responder' && !ui.removeMark.visible : 'a responder can be marked'
        // a responder without a flag: no reason is chosen for the user, and Save waits for one
        ui.markButton.doClick()
        assert ui.formCell=='Cell001' && ui.reasonBox.visible && ui.noteField.visible && ui.saveMark.visible && ui.cancelMark.visible : 'the form opens in place'
        assert ui.reasonBox.selectedIndex==-1 && !ui.saveMark.enabled : 'no reason is chosen on a cell without a flag'
        assert (0..<ui.reasonBox.itemCount).collect{ui.reasonBox.getItemAt(it)}==Review.REASONS : 'the reasons are the review file\'s'
        ui.reasonBox.selectedItem='Noise';assert ui.saveMark.enabled : 'Save once a reason is chosen'
        ui.cancelMark.doClick();assert ui.formCell==null && !ui.reasonBox.visible && !ui.saveMark.visible && ui.cells.selectedItem=='Cell001' : 'Cancel closes the form and keeps the cell'
        // the run's own answer is kept: with its table taken away, the checks are still the run's (as the page takes them)
        assert new File(work,'Review_call_checks_'+A+'.csv').delete() : 'the checks table of the run can be taken away'
        ui.loadResults(dirA);assert ui.checks?.note==unchecked : 'the checks of the last run here are kept, not only the saved table'
    }

    // ---- an analysis with flags, and a control window beside its stimulus
    File dirB=copy(B)
    def table=Csv.read(new File(dirB,'responder_calls.csv')),header=table[0]
    Csv.write(new File(dirB,'responder_calls.csv'),table+table.drop(1).collect{List line->def row=new ArrayList(line);['Role':'Control','Stimulus':'Ionomycin','Threshold_increase':'0.4','Floor_threshold':'0.4','Baseline_mean':'0.5'].each{k,v->row[header.indexOf(k)]=v};row})
    File checksB=new File(work,'Review_call_checks_'+B+'.csv'),reviewB=new File(work,'Review_reviewed_calls_'+B+'.csv'),percentB=new File(work,'Review_percent_responding_'+B+'.csv')
    def row={String cell,String call,String flag,String rule,Object first,Object last,String detail,String check,String why->[cell,call,flag,rule,first,last,detail,8,400,'yes','yes',check,why,'','']}
    def none={String cell,String call->row(cell,call,'','','','','','not checked','no call')}
    // first: a field that moved and short separate crossings on one cell, with the first frame of its event
    Csv.write(checksB,[Checks.HEADER,row('Cell001','Responded','','','','','','checked',''),row('Cell002','Did not respond','','','','','','checked',''),none('Cell003','Excluded'),none('Cell004','Excluded'),
        row('Cell005','Responded','field moved, frames 20 to 22; short separate crossings','field moved; short separate crossings',18,23,'the field moved 3.1 px between frames 20 and 22; otherwise possible passing object: seen in the near ring first; 3 crossings, the longest 2 samples','checked',''),none('Cell006','Unscorable')])
    Settings.percentResponding(dirB,ui.measurement)   // the percent table beside the tables, as a run leaves it
    byte[] percentBefore=percentB.bytes
    onUI {
        ui.loadResults(dirB);ui.cells.selectedItem='Cell005'
        assert ui.flagRow.visible && ui.flagRow.components*.text==['field moved, frames 20 to 22','·','short separate crossings'] : ui.flagRow.components*.text
        def moved=ui.flagRow.components[0]
        assert moved.foreground==Theme.BLUE && moved.cursor.type==java.awt.Cursor.HAND_CURSOR && moved.toolTipText.startsWith('<html>Go to frame 18<br>the field moved 3.1 px') : 'a flag with a first frame is a link'
        ui.markButton.doClick();assert ui.reasonBox.selectedItem=='Movement' : 'a field that moved: Movement is chosen'
        ui.cancelMark.doClick()
    }
    // then: a possible passing object on a responder, a rise on the edge of a cell that did not respond, short separate
    // crossings on a cell whose ring could not be measured
    String passing='possible passing object, frames 12 to 14',edge='rise on the lower right edge, not in the cell, frames 5 to 9'
    Csv.write(checksB,[Checks.HEADER,row('Cell001','Responded',passing,'passes through',12,14,'seen in the near ring 2 frames before the cell','checked',''),
        row('Cell002','Did not respond',edge,'on the edge',5,9,'the lower right rim rose first','checked',''),none('Cell003','Excluded'),none('Cell004','Excluded'),
        row('Cell005','Responded','short separate crossings','short separate crossings','','','3 crossings, the longest 2 samples','not checked','no ring around the cell'),none('Cell006','Unscorable')])
    Map before=[:]
    onUI {
        ui.loadResults(dirB);ui.showGroup(0);ui.cells.selectedItem='Cell001'
        // without marks every view is as it was before review marks existed
        before=[status:ui.cellStatus.text,tip:ui.cellStatus.toolTipText,groups:(0..<ui.group.itemCount).collect{ui.group.getItemAt(it)},tile:ui.responding.value.text,caption:ui.responding.caption.text,
                picture:new LinkedHashMap(ui.cellPicture.calls),heatmap:[ui.heatmap.names,ui.heatmap.calls].transpose(),plot:ui.plot.caption]
        assert before.groups==['All cells (6)','Responded (2)','Did not respond (1)','Excluded or unscorable (3)'] && before.tile=='66.7 %' && strip(before.caption)=='Capsaicin: 2 of 3 evaluable cells' : before
        assert strip(before.status)=='● Cell001: Responded · sustained threshold crossing · largest rise 1.00 at 6 s (needed 0.10)' && before.plot=='Cell001 · Responded · Sustained threshold crossing' : strip(before.status)+' | '+before.plot
        // the flag of the selected cell, in the checks' own words, and its link to the first frame of the event
        def flag=ui.flagRow.components[0]
        assert ui.flagRow.visible && ui.flagRow.componentCount==1 && flag.text==passing && flag.foreground==Theme.BLUE && flag.cursor.type==java.awt.Cursor.HAND_CURSOR : 'the flag in its own words, a link'
        assert flag.toolTipText=='<html>Go to frame 12<br>seen in the near ring 2 frames before the cell</html>' : flag.toolTipText
        ui.view.selectedIndex=0;ui.tabs.selectedIndex=0;ui.showFrame(25)
        click(flag)
        assert ui.tabs.selectedIndex==1 && ui.view.selectedIndex==2 && ui.frame==12 && ui.slider.value==12 && ui.cellPicture.image.is(ui.frameImage) : 'a click shows the recording at the first frame of the event'
        ui.view.selectedIndex=1;ui.showFrame(25);ui.playButton.doClick();assert ui.player.running : 'the movie plays'
        click(flag);assert !ui.player.running && ui.view.selectedIndex==1 && ui.frame==12 : 'a click stops playing; a movie view stays as it was chosen'
        // a flag without a first frame is plain text; a cell that could not be checked says so, in the recording whose rings were measured
        ui.cells.selectedItem='Cell005'
        assert ui.flagRow.components*.text==['short separate crossings','·','passing object not checked: no ring around the cell'] : ui.flagRow.components*.text
        def plain=ui.flagRow.components[0]
        assert plain.foreground==Theme.MUTED && plain.cursor.type!=java.awt.Cursor.HAND_CURSOR && plain.toolTipText=='3 crossings, the longest 2 samples' : 'without a first frame a flag is plain text'
        ui.markButton.doClick();assert ui.reasonBox.selectedItem=='Noise' : 'short separate crossings: Noise is chosen';ui.cancelMark.doClick()
        // a cell that did not respond: the other button, and no reason chosen even with a flag
        ui.cells.selectedItem='Cell002'
        assert ui.flagRow.components*.text==[edge] && ui.markButton.visible && ui.markButton.text=='Mark as a responder' : 'a cell that did not respond can be marked as a responder'
        ui.markButton.doClick();assert ui.reasonBox.selectedIndex==-1 && !ui.saveMark.enabled : 'never Passing object chosen on a cell that did not respond'
        def nothing=snapshot();press(ui.reasonBox,KeyEvent.VK_ENTER);assert ui.formCell=='Cell002' && snapshot()==nothing && !ui.markProblem.visible : 'Enter does nothing before a reason is chosen: '+strip(ui.markProblem.text)
        ui.cancelMark.doClick()
        // no call to review: no button, no flags
        ['Cell003','Cell004','Cell006'].each{ui.cells.selectedItem=it;assert !ui.markRow.visible && !ui.markButton.visible && !ui.removeMark.visible && !ui.flagRow.visible : it+' has no call to review'}
        // the new controls have the shape and size of the ones beside them
        layOut()
        [ui.markButton,ui.removeMark,ui.saveMark,ui.cancelMark].each{assert it.preferredSize.height==ui.previousCell.preferredSize.height && it.UI.class.name=='CalciumButtonUI' : 'a mark button is a button of this window: '+it.text}
        assert ui.saveMark.UI.primary && !ui.cancelMark.UI.primary && !ui.markButton.UI.primary : 'Save is the filled one'
        assert ui.reasonBox.preferredSize.height==ui.group.preferredSize.height && ui.reasonBox.font==ui.group.font && ui.reasonBox.font==Theme.font() : 'the Reason list is a list like Show'
        assert ui.noteField.preferredSize.height==ui.fields.sampling_hz.preferredSize.height && ui.noteField.font==ui.fields.sampling_hz.font : 'the Note box is a box like the settings'
        assert ui.reasonLabel.font==Theme.font() && ui.noteLabel.font==Theme.font() && ui.reasonLabel.text=='Reason' && ui.noteLabel.text=='Note' : 'the form\'s labels: '+ui.reasonLabel.text+', '+ui.noteLabel.text

        // ---- the form on a flagged responder: Passing object chosen; the refusals in the form, in CalciumReview's words
        ui.cells.selectedItem='Cell001';ui.markButton.doClick()
        assert ui.formCell=='Cell001' && ui.reasonBox.selectedItem=='Passing object' && ui.saveMark.enabled && ui.noteField.hint=='Optional' && ui.noteField.text=='' : 'a cell that may have a passing object: Passing object is chosen'
        def refused={String reason,String note,String message->def files=snapshot();ui.reasonBox.selectedItem=reason;ui.noteField.text=note;ui.saveMark.doClick()
            assert ui.markProblem.visible && strip(ui.markProblem.text)==message : strip(ui.markProblem.text)
            assert ui.formCell=='Cell001' && snapshot()==files : 'nothing is written: '+message}
        refused('Other','   ','Other needs a short note.')
        assert ui.noteField.hint=='' : 'the note is not optional for Other'
        refused('Noise','=SUM(A1)','The note cannot start with =, +, - or @.')
        refused('Noise','  -1 frame','The note cannot start with =, +, - or @.')   // the spaces at either end go first, as in the review file
        refused('Noise','one\ttwo','The note must be one line.')
        ui.noteField.text='x'*130;assert ui.noteField.text.length()==120 : 'the note box takes at most 120 characters'
        // the review's own refusals: a review file that is not this analysis's, a percent table open in another program, a run going on
        reviewB.text=Review.HEADER.join(',')+'\nCell002,Did not respond,Responded,Noise,,2026-11-01T10:00:00Z,another,x,,\n'
        refused('Passing object','','The review file Review_reviewed_calls_'+B+'.csv is for another analysis, so its marks are not used. Move it out of the folder to review this analysis again.')
        assert reviewB.delete() : 'the other analysis\'s review file is taken away'
        assert percentB.setWritable(false) : 'the percent table is locked'
        try{refused('Passing object','','Close Review_percent_responding_'+B+'.csv if it is open in another program, then try again.')}finally{percentB.setWritable(true)}
        ui.busy=true;try{refused('Passing object','','Wait for the analysis to finish.');ui.showReview('Cell001');assert !ui.markRow.visible && ui.formCell==null : 'no mark controls while an analysis runs'}finally{ui.busy=false}
        ui.showReview('Cell001');assert ui.markRow.visible && !ui.markProblem.visible : 'and they come back'
        // the form closes when the cell changes; the arrows do not change it while the form is open; Esc closes it
        ui.markButton.doClick();ui.noteField.text='=x';ui.saveMark.doClick();assert ui.markProblem.visible : 'a refused note is said in the form'
        ui.cells.selectedItem='Cell005';assert ui.formCell==null && !ui.reasonBox.visible && !ui.markProblem.visible : 'the form and its refusal go when the cell changes'
        ui.cells.selectedItem='Cell001';ui.markButton.doClick()
        ui.rootPane.getActionForKeyStroke(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT,0)).actionPerformed(new java.awt.event.ActionEvent(ui.rootPane,0,''))
        assert ui.cells.selectedItem=='Cell001' && ui.formCell=='Cell001' : 'the arrows do not step to another cell while the form is open'
        assert press(ui.noteField,KeyEvent.VK_ESCAPE) && ui.formCell==null && ui.cells.selectedItem=='Cell001' : 'Esc closes the form and keeps the cell'
        ui.markButton.doClick();assert press(ui.reasonBox,KeyEvent.VK_ESCAPE) && ui.formCell==null : 'Esc closes the form from the Reason list too'
        ui.markButton.doClick();ui.markButton.doClick();assert ui.formCell==null : 'the button closes the form it opened'
        assert !reviewB.exists() && percentB.bytes==percentBefore : 'no mark was saved'

        // ---- a mark saved: every view follows it at once, the cell stays selected and the time cursor stays
        ui.showGroup(1);assert ui.cells.selectedItem=='Cell001' : 'Cell001 is a responder'
        ui.markButton.doClick();ui.noteField.text='  drifts in from the left  '
        paint('fiji_traces_form')
        ui.showFrame(17);double moment=ui.plot.moment
        ui.saveMark.doClick()
        assert !ui.markProblem.visible && ui.formCell==null : strip(ui.markProblem.text)
        Map mark=Review.read(dirB,ui.measurement).marks.Cell001
        assert mark.automatic=='Responded' && mark.reviewed=='Did not respond' && mark.reason=='Passing object' && mark.note=='drifts in from the left' : mark
        def saved=Csv.read(reviewB)
        assert saved[0]==Review.HEADER && saved.size()==2 && saved[1][0..6]==['Cell001','Responded','Did not respond','Passing object','drifts in from the left',mark.time,B] && saved[1][7]==Installer.digest(new File(dirB,'responder_calls.csv').bytes) : saved
        assert Csv.read(percentB)[0]==Settings.PERCENT_COLUMNS+Settings.REVIEW_COLUMNS : 'the percent table gives the figures after review beside the others'
        assert ui.responding.value.text=='33.3 %' && strip(ui.responding.caption.text)=='Capsaicin: 1 of 3 evaluable cells after review; 2 before review' : ui.responding.caption.text
        assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==['All cells (6)','Responded after review (1)','Did not respond after review (2)','Excluded or unscorable (3)','Marked (1)'] : (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}
        assert ui.group.selectedIndex==0 && ui.cells.selectedItem=='Cell001' : 'the cell left Responded: all cells are listed and it stays selected'
        assert ui.frame==17 && ui.plot.moment==moment : 'the time cursor did not jump'
        assert ui.cellPicture.calls.Cell001=='Did not respond' && ui.cellPicture.calls.Cell005=='Responded' : ui.cellPicture.calls
        assert ui.heatmap.calls[ui.heatmap.names.indexOf('Cell001')]=='Did not respond' && ui.heatmap.names.indexOf('Cell001')>ui.heatmap.names.indexOf('Cell005') : 'the heatmap orders and colours by the call after review'
        assert ui.plot.allCalls.count{it=='Responded'}==1 && ui.plot.allCalls.count{it=='Did not respond'}==2 : ui.plot.allCalls
        assert strip(ui.cellStatus.text)=='● Cell001: Did not respond after review · passing object: drifts in from the left · automatic call Responded · sustained threshold crossing · largest rise 1.00 at 6 s (needed 0.10)' : strip(ui.cellStatus.text)
        assert ui.cellStatus.text.contains("<font color='"+Theme.hex(Calls.NOT)+"'>●</font>") : 'the dot has the colour of the call after review'
        assert ui.cellStatus.toolTipText=='<html>Marked '+mark.time+'<br>'+before.tip+'</html>' : ui.cellStatus.toolTipText
        assert ui.plot.caption=='Cell001 · Did not respond · Sustained threshold crossing' : ui.plot.caption
        assert ui.removeMark.visible && !ui.markButton.visible && ui.flagRow.visible : 'a marked cell has one button: Remove mark'
        paint('fiji_traces_marked')
        // the long call line and the tile's caption wrap at the width they are given: nothing is cut off
        [ui.cellStatus,ui.responding.caption].each{assert it.width>=it.preferredSize.width && it.height>=it.preferredSize.height : 'cut off: '+strip(it.text)+' '+it.size+' < '+it.preferredSize}
        assert ui.responding.caption.y+ui.responding.caption.height<=ui.responding.height-ui.responding.insets.bottom : 'the tile holds its caption: '+ui.responding.caption.bounds+' in '+ui.responding.size
        // the Marked entry lists the marked cells, on the picture too
        ui.showGroup(4);assert (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell001'] && ui.cellPicture.only==['Cell001'] as Set && ui.plot.all.length==1 : 'the Marked entry lists and draws only the marked cell'
        ui.showGroup(0);assert ui.cellPicture.only==null : 'All cells draws every cell'
        // a control window keeps its automatic call: no marks, no flags, its own groups
        ui.stimuli.selectedItem='Ionomycin'
        assert !ui.markRow.visible && !ui.flagRow.visible && (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==['All cells (6)','Responded (2)','Did not respond (1)','Excluded or unscorable (3)'] : 'the control window is never reviewed'
        assert strip(ui.cellStatus.text).startsWith('● Cell001: Responded · ') && ui.plot.caption.startsWith('Cell001 · Responded · ') : 'the control window shows its automatic call: '+strip(ui.cellStatus.text)
        ui.stimuli.selectedItem='Capsaicin';ui.cells.selectedItem='Cell001';assert ui.removeMark.visible : 'back on the stimulus, the mark can be removed'
        // saved plots say their groups follow the review marks, and which file the marks are in
        layOut();ui.exportPlots()
        def record=new groovy.json.JsonSlurper().parse(new File(dirB,'plot_provenance.json'))
        assert record.calls=='after review' && record.review_file==reviewB.name && record.review_sha256==Installer.digest(reviewB.bytes) : ui.status.text

        // ---- the mark removed: everything is as it was before
        ui.removeMark.doClick()
        assert !reviewB.exists() && percentB.bytes==percentBefore : 'the review file goes and the percent table is as it was'
        assert ui.responding.value.text==before.tile && ui.responding.caption.text==before.caption : strip(ui.responding.caption.text)
        assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==before.groups && ui.cellPicture.calls==before.picture && [ui.heatmap.names,ui.heatmap.calls].transpose()==before.heatmap : 'the Show list, the picture and the heatmap are as before'
        assert ui.cellStatus.text==before.status && ui.cellStatus.toolTipText==before.tip && ui.plot.caption==before.plot && ui.cells.selectedItem=='Cell001' : strip(ui.cellStatus.text)
        assert ui.markButton.visible && ui.markButton.text=='Mark as not a responder' && !ui.removeMark.visible : 'the cell can be marked again'
        ui.exportPlots();record=new groovy.json.JsonSlurper().parse(new File(dirB,'plot_provenance.json'))
        assert !record.containsKey('calls') && !record.containsKey('review_file') : 'without marks the record is as it always was'

        // ---- Enter saves, from the Note box and from the Reason list; a mark without a note
        ui.markButton.doClick();ui.noteField.text='drifts in from the left'
        assert press(ui.noteField,KeyEvent.VK_ENTER) && ui.formCell==null && Review.read(dirB,ui.measurement).marks.Cell001?.reviewed=='Did not respond' : 'Enter in the Note box saves'
        ui.cells.selectedItem='Cell002';ui.markButton.doClick();ui.reasonBox.selectedItem='Bad outline'
        press(ui.reasonBox,KeyEvent.VK_ENTER)
        assert ui.formCell==null && Review.read(dirB,ui.measurement).marks.Cell002?.reviewed=='Responded' : 'Enter in the Reason list saves: '+strip(ui.markProblem.text)
        assert strip(ui.cellStatus.text)=='● Cell002: Responded after review · bad outline · automatic call Did not respond · no sustained threshold crossing · largest rise 0.00 at 3 s (needed 0.10)' : strip(ui.cellStatus.text)
        assert ui.cellStatus.text.contains("<font color='"+Theme.hex(Calls.RESPONDED)+"'>●</font>") : 'the dot has the colour of Responded'
        assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}[-1]=='Marked (2)' && ui.responding.value.text=='66.7 %' : 'two marks: '+ui.responding.value.text
        ui.removeMark.doClick();assert Review.read(dirB,ui.measurement).marks.keySet()==['Cell001'] as Set : 'removing one mark keeps the other'
    }

    // ---- a run going on: no mark controls, and a mark is refused; then the next analysis takes the mark over
    def latch=new java.util.concurrent.CountDownLatch(1)
    onUI {
        ui.cells.selectedItem='Cell005';ui.markButton.doClick();assert ui.formCell=='Cell005' : 'the form is open when the run starts'
        ui.performAnalysis={a,b,c,d,job->latch.await(20,java.util.concurrent.TimeUnit.SECONDS);dirA}
        ui.runButton.doClick()
        assert ui.busy && ui.formCell==null && !ui.markRow.visible && !ui.markButton.visible && !ui.removeMark.visible : 'no mark controls while an analysis runs'
        ui.cells.selectedItem='Cell001';assert !ui.markRow.visible : 'no Remove mark while an analysis runs'
        def files=snapshot();ui.removeReviewMark('Cell001')
        assert ui.markProblem.visible && strip(ui.markProblem.text)=='Wait for the analysis to finish.' && snapshot()==files : 'a mark is refused while an analysis runs'
        ui.cells.selectedItem='Cell002';assert !ui.markProblem.visible : 'the refusal goes when another cell is chosen'
    }
    latch.countDown();finish()
    onUI {
        String carried='Review marks from '+B+': 1 carried over, 0 dropped.'
        assert ui.results==dirA && ui.status.text=='Complete. Protocol, inputs and tables saved in '+A+'. All cells: 2 of 3 responded (66.7 %), after review 1 of 3 (33.3 %). '+carried+' '+unchecked : ui.status.text
        ui.showGroup(0);ui.cells.selectedItem='Cell001'
        Map mark=Review.read(dirA,ui.measurement).marks.Cell001
        assert mark.carriedFrom==B : mark
        assert strip(ui.cellStatus.text)=='● Cell001: Did not respond after review · passing object: drifts in from the left · from an earlier analysis · automatic call Responded · sustained threshold crossing · largest rise 1.00 at 6 s (needed 0.10)' : strip(ui.cellStatus.text)
        assert ui.cellStatus.toolTipText.startsWith('<html>Marked '+mark.time+' in '+B+'<br>') : 'when and where the mark was made, in the tooltip'
        paint('fiji_traces_carried')
        // the carried text stays until the next run, after the last mark is removed too
        ui.removeMark.doClick()
        assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+A+'. '+carried+' '+unchecked : 'the finish message follows the marks: '+ui.status.text
        assert ui.responding.value.text=='66.7 %' && !new File(work,'Review_reviewed_calls_'+A+'.csv').exists() : 'the tile and the review files are as before'
    }
    assert appeared.isEmpty() : 'a window appeared on screen: '+appeared
    onUI{assert !ui.visible && !ui.isDisplayable() : 'the window was never shown'}
    println 'PASS: Step 4 window flags, links to frames, mark controls, reasons chosen from flags, refusals, Enter and Esc, marks saved and removed with every view following, control window and busy without controls, carried marks, plot record'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();original?.close();recording?.close()}}
System.exit(0)

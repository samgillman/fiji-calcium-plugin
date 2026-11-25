// The Step 4 window of Fiji with the checks of the calls and the review marks: the flags of the selected cell under its call
// line, each in the checks' own words, and a flag's link to the first frame of its event; the mark button on a cell with a
// call, the form in its place with the reason already chosen from a flag, Enter (once when held down, and never from the open
// Reason list), Esc and the refusals in the form, and the focus back on the button after it; a mark saved and removed, with
// the tile, the Show list (which keeps
// the group shown, and the cell in the Cell list until another is selected), the picture, the heatmap, the plot, the call
// line, the calls and counts tabs, the finish message and the record of saved plots following it, and the time cursor
// staying where it is; the head over the picture and the plot one height for every cell, on a control window, with a refusal
// of two lines and with the longest note, with the call line, its flags, the mark controls and a refusal directly under one
// another; the tiles, the tabs and the plot not moving at a mark, nor Show, the Cell box and the arrows, in either window;
// nothing under the plot naming a cell when Show lists none; nothing written on the plot or the saved trace.png but its
// axes and key; on screen the key holds the cell, the mean and the traces behind it, so the plot starts where it did before
// the review controls, and a band is named by pointing at it; trace.png's key also names its bands and dashed lines in the
// Step 4 page's words; the plot keeping its place from cell to cell, for a cell without a trace too, and at a mark that
// empties the group shown or leaves no cell of a call, also as narrow as beside a real recording's picture; under the plot
// the caption and the rule line, which gives the required rise with the call line's decimals, and no rule line for a cell
// without a trace; trace.png and heatmap.png drawn across their whole width; no mark controls on
// excluded or unscorable cells, on a control window, while an analysis runs or while the review file is for another
// analysis; a mark carried over to the next analysis run here; with a recording's long cell names the Cell box as wide as
// for Cell001, the name cut in it and whole in its tooltip, and the pickers row whole in the Traces tab. On a recording whose
// rings the run measures from the open movie: the tick box that marks flagged responders automatically, like the one beside
// it and part of the settings, the cells it marked in the finish message, the call line, the Review reason column and the
// Marked group, an automatic mark removed and made by hand; the codes of the checks after the engine's in QC_flags on every
// tab, none on the control window, the saved tables unchanged; the Around the cell tab with the call line and flags of
// Traces, the zones of the flag each in its colour and named by place, the other ring sectors, the flagged frames and the
// response window, the key and the plot in one place for every cell and width, the time cursor, a sentence where there is
// nothing to draw, rule C's lines without frames, and around_the_cell.png as on screen without the cursor.
// The window is built and never shown (with -Dapple.awt.UIElement=true no Dock icon appears on a Mac); a window that opens all
// the same fails the test. -Dcalcium.test.pictures=<folder> saves pictures of the window there, and a saved trace plot.
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
/** A graphics that draws as the one it holds and keeps every string drawn with it, and where it was drawn. */
class TextSpy extends java.awt.Graphics2D {
    @Delegate java.awt.Graphics2D real
    List drawn=[],at=[]
    void drawString(String s,int x,int y){drawn<<s;at<<[s,x,y];real.drawString(s,x,y)}
    void drawString(String s,float x,float y){drawn<<s;at<<[s,(int)x,(int)y];real.drawString(s,x,y)}
}
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
// A file that stays where it is and cannot be read, so a test can show that nothing reads it: on macOS and Linux it is made
// unreadable; Windows has no unreadable files, so there it is locked, and Windows refuses reads of a locked file, also from
// this process. Returns what gives the file back. Fails at once when the file can still be read.
def lockAway={File f->
    Closure back
    if(f.setReadable(false)&&!f.canRead())back={->f.setReadable(true)}
    else{def channel=java.nio.channels.FileChannel.open(f.toPath(),java.nio.file.StandardOpenOption.READ,java.nio.file.StandardOpenOption.WRITE);def lock=channel.lock();back={->lock.release();channel.close()}}
    boolean read=false
    try{f.bytes;read=true}catch(IOException ignored){}
    if(read){back();assert false : 'the file cannot be read while it is locked: '+f}
    assert f.isFile() : 'the locked file is still there: '+f
    back
}
List names=(1..6).collect{String.format('Cell%03d',it)}
String unchecked='Not checked for a passing object: the raw gray values and backgrounds are not beside the analysis.'

// A key pressed in the window, as Swing hands it on: the focused control's own keys, then those of the controls around it.
def press={JComponent focus,int code->
    def ks=KeyStroke.getKeyStroke(code,0),e=new KeyEvent(focus,KeyEvent.KEY_PRESSED,System.currentTimeMillis(),0,code,KeyEvent.CHAR_UNDEFINED)
    def fire={JComponent c,int when->def key=c.getInputMap(when)?.get(ks);def action=key!=null?c.actionMap.get(key):null;action!=null&&SwingUtilities.notifyAction(action,ks,e,c,0)}
    if(fire(focus,JComponent.WHEN_FOCUSED))return true
    for(def c=focus;c instanceof JComponent;c=c.parent)if(fire(c as JComponent,JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT))return true
    false}
// the key let go: Enter acts then, once however long it was held down (its press repeats, it is let go once)
def release={JComponent focus,int code->
    def ks=KeyStroke.getKeyStroke(code,0,true),e=new KeyEvent(focus,KeyEvent.KEY_RELEASED,System.currentTimeMillis(),0,code,KeyEvent.CHAR_UNDEFINED)
    def fire={JComponent c,int when->def key=c.getInputMap(when)?.get(ks);def action=key!=null?c.actionMap.get(key):null;action!=null&&SwingUtilities.notifyAction(action,ks,e,c,0)}
    if(fire(focus,JComponent.WHEN_FOCUSED))return true
    for(def c=focus;c instanceof JComponent;c=c.parent)if(fire(c as JComponent,JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT))return true
    false}
// a key pressed and let go
def enter={JComponent focus->press(focus,KeyEvent.VK_ENTER);release(focus,KeyEvent.VK_ENTER)}
def click={java.awt.Component c->c.dispatchEvent(new java.awt.event.MouseEvent(c,java.awt.event.MouseEvent.MOUSE_CLICKED,System.currentTimeMillis(),0,4,4,1,false,java.awt.event.MouseEvent.BUTTON1))}

def ui,original,recording,longWindow,aroundWindow,aroundMovie
def finish={->long deadline=System.nanoTime()+30000000000L;while(ui.busy&&System.nanoTime()<deadline)Thread.sleep(20);onUI{};assert !ui.busy : 'the analysis finishes'}
// The window laid out at its usual size without being shown, and its content painted into a picture. A window that is never
// shown is never validated, so every layout forgets the sizes it kept first (a shown window does this when it validates).
def layOut={->ui.setSize(1180,830);def forget,lay
    forget={java.awt.Component c->c.invalidate();if(c instanceof java.awt.Container)(c as java.awt.Container).components.each{forget(it)}}
    lay={java.awt.Component c->if(c instanceof java.awt.Container){(c as java.awt.Container).doLayout();(c as java.awt.Container).components.each{lay(it)}}}
    forget(ui);lay(ui);forget(ui);lay(ui)}   // twice: wrapped text takes its height from the width the first pass gave it
def paint={String name,boolean traces=true->layOut();def pane=ui.contentPane as JComponent
    def image=new java.awt.image.BufferedImage(pane.width,pane.height,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics();try{pane.paint(g)}finally{g.dispose()}
    javax.imageio.ImageIO.write(image,'png',new File(pictures,name+'.png'))
    if(!traces)return
    def head=ui.reviewBox.parent as JComponent;def top=SwingUtilities.convertPoint(head,0,0,pane)
    javax.imageio.ImageIO.write(image.getSubimage(Math.max(0,top.x as int),Math.max(0,top.y as int),Math.min(head.width,pane.width-(top.x as int)),Math.min(head.height,pane.height-(top.y as int))),'png',new File(pictures,name+'_head.png'))}
// the flags under the call line, and the refusal under the controls
def flags={->ui.flagRow.components.findAll{it instanceof JLabel}}
def said={->strip(ui.markProblem.text)}
// the head over the picture and the plot, the plot, and where the tabs start under the tiles, as laid out at the window's usual size
def heights={->layOut();[head:(ui.reviewBox.parent as JComponent).height,plot:ui.plot.height,tabs:ui.tabs.y]}
// the rows of the head that are shown, from the top: each directly under the one before, the room left of the head under them
def stacked={->layOut();List rows=[ui.cellStatus,ui.flagRow,ui.markRow,ui.markProblem].findAll{it.visible}
    (1..<rows.size()).every{int i->rows[i].y==rows[i-1].y+rows[i-1].height} && ui.headRoom.y==rows[-1].y+rows[-1].height}
// the width of Show, and where the Cell box and the next arrow beside it start, as laid out at the window's usual size
def pickers={->layOut();[show:ui.group.width,cell:ui.cells.x,next:ui.nextCell.x]}
// which rows of the head are shown, by name
def shownRows={->[call:ui.cellStatus,flags:ui.flagRow,mark:ui.markRow,refusal:ui.markProblem].findAll{k,v->v.visible}.keySet().toList()}
// a table on its tab: its column names, and a row by the values of two of its columns
def columns={JTable t->(0..<t.model.columnCount).collect{t.model.getColumnName(it)}}
def rowOf={JTable t,Map where->List heads=columns(t);int r=(0..<t.model.rowCount).find{int i->where.every{k,v->t.model.getValueAt(i,heads.indexOf(k))==v}} as int;(0..<heads.size()).collect{t.model.getValueAt(r,it)}}
// the strings the plot draws on itself: only tick labels, the axis title and the key may be there
def drawnOnPlot={->def image=new java.awt.image.BufferedImage(Math.max(1,ui.plot.width),Math.max(1,ui.plot.height),java.awt.image.BufferedImage.TYPE_INT_RGB)
    def spy=new TextSpy(real:image.createGraphics());try{ui.plot.still(spy)}finally{spy.real.dispose()};spy.drawn}
// the plot as Save plots draws it into trace.png: at the width saved, with the key of a saved figure; then as on screen again
def asSaved={Closure c->java.awt.Dimension was=ui.plot.size;ui.plot.setSize(Math.max(650,was.width as int),was.height as int);ui.plot.figure=true
    try{c()}finally{ui.plot.figure=false;ui.plot.setSize(was);ui.plot.changed()}}
def drawnInFigure={->asSaved{drawnOnPlot()}}
// trace.png as the plot of the cell shown would be saved now, on white
def savedPicture={int width,int height->def image=new java.awt.image.BufferedImage(width,height,java.awt.image.BufferedImage.TYPE_INT_RGB);def g=image.createGraphics()
    java.awt.Dimension was=ui.plot.size;ui.plot.setSize(width,height);ui.plot.figure=true
    try{g.color=java.awt.Color.white;g.fillRect(0,0,width,height);ui.plot.still(g)}finally{g.dispose();ui.plot.figure=false;ui.plot.setSize(was);ui.plot.changed()};image}
// the swatch beside an entry of trace.png's key, as the plot draws it on white: the colour of each pixel of the row it is
// centred on, from its left end, and the colour 3 px above its middle (a box is 10 px tall, a line is not)
def swatch={String name->asSaved{def image=new java.awt.image.BufferedImage(Math.max(1,ui.plot.width),Math.max(1,ui.plot.height),java.awt.image.BufferedImage.TYPE_INT_RGB)
    def spy=new TextSpy(real:image.createGraphics());try{spy.real.color=java.awt.Color.white;spy.real.fillRect(0,0,image.width,image.height);ui.plot.still(spy)}finally{spy.real.dispose()}
    def text=spy.at.find{it[0]==name};assert text!=null : 'no entry '+name+' in '+spy.drawn
    int x=(text[1] as int)-21,y=(text[2] as int)-4
    [row:(0..16).collect{image.getRGB(x+it,y)&0xFFFFFF},above:image.getRGB(x+8,y-3)&0xFFFFFF]}}
// the names of the dashed lines and the bands, in trace.png's key only
List keyWords=['required rise','noise threshold','minimum rise','baseline mean','F₀ frames','baseline','response window','response']
def counted={String s->s==~/(responded|did not respond|no call|all cells|marked|responded after review) \(\d+\)/}
def allowedOnScreen={String s->s in ['Time (s)','mean ± SEM',ui.plot.cell] || counted(s) || s==~/-?\d+(\.\d+)?/}
def allowedInFigure={String s->allowedOnScreen(s) || s in keyWords}
// the plot's frame: where it starts under the key
def frameTop={->layOut();drawnOnPlot();ui.plot.@t}
// and at the width the plot has beside the picture of a real recording (245 px beside 1388 x 1040), where most entries of
// the key take a line of their own
def narrowTop={->layOut();ui.plot.setSize(245,ui.plot.height);drawnOnPlot();int top=ui.plot.@t;layOut();top}
// and in trace.png, whose key also names the dashed lines and the bands
def savedTop={->layOut();int top=asSaved{drawnOnPlot();ui.plot.@t};layOut();top}
def placeOf={->[frameTop(),narrowTop(),savedTop()]}
// where the frame of the plot as it is now started at a width before the review controls: the key held the mean, the cell
// and the traces behind them by call, laid from the frame's left edge across the frame's width (the panel's less 80 px), a
// line more for each entry that did not fit, and the frame started 16 px under its last line
def mainTop={int width->def g=new java.awt.image.BufferedImage(1,1,java.awt.image.BufferedImage.TYPE_INT_RGB).createGraphics();def fm
    try{Theme.smooth(g);g.font=Theme.font(java.awt.Font.PLAIN,-1f);fm=g.fontMetrics}finally{g.dispose()}
    def p=ui.plot;int yes=p.allCalls.count{it=='Responded'},no=p.allCalls.count{it=='Did not respond'},other=p.all.length-yes-no
    List words=['mean ± SEM',p.cell];if(p.showAll&&p.all.length>0){if(yes)words<<'responded ('+yes+')';if(no)words<<'did not respond ('+no+')';if(other)words<<(yes+no>0?'no call':p.allLabel)+' ('+other+')'}
    int across=Math.max(1,width-62-18),line=0,run=0
    words.each{String s->int need=21+fm.stringWidth(s)+14;if(run>0&&run+need>across){line++;run=0};run+=need}
    30+16*(line+1)}
// the time of a mark as the Step 4 page shows it, to the minute
def minute={String time->assert time==~/\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d.*/ : time;time.substring(0,16).replace('T',' ')}
// where the window asks for the keyboard focus, and which component has it (a window never shown has no focus of its own)
List asked=[];def focusOn=null
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
        ui.takeFocus={java.awt.Component c->asked<<c;true};ui.focused={->focusOn}
        assert !ui.visible && !ui.markRow.visible && !ui.flagRow.visible && !ui.markProblem.visible : 'before an analysis there is nothing to mark'
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
        assert !ui.flagRow.visible && flags().isEmpty() : 'no flag, and no "not checked" on each cell when the surroundings were not measured'
        assert ui.markRow.visible && ui.markButton.visible && ui.markButton.text=='Mark as not a responder' && !ui.removeMark.visible : 'a responder can be marked'
        assert shownRows()==['call','mark'] && stacked() : 'without flags the button is right under the call line: '+shownRows()
        // a responder without a flag: no reason is chosen for the user, and Save waits for one
        asked.clear();ui.markButton.doClick()
        assert ui.formCell=='Cell001' && ui.reasonBox.visible && ui.noteField.visible && ui.saveMark.visible && ui.cancelMark.visible : 'the form opens in place'
        assert !ui.markButton.visible && !ui.removeMark.visible : 'the form takes the place of the button'
        assert asked.size()==1 && asked[0].is(ui.reasonBox) : 'the focus goes to the Reason list: '+asked
        assert ui.reasonBox.selectedIndex==-1 && !ui.saveMark.enabled : 'no reason is chosen on a cell without a flag'
        assert (0..<ui.reasonBox.itemCount).collect{ui.reasonBox.getItemAt(it)}==Review.REASONS : 'the reasons are the review file\'s'
        ui.reasonBox.selectedItem='Noise';assert ui.saveMark.enabled : 'Save once a reason is chosen'
        focusOn=ui.reasonBox;asked.clear()
        ui.cancelMark.doClick();assert ui.formCell==null && !ui.reasonBox.visible && !ui.saveMark.visible && ui.cells.selectedItem=='Cell001' : 'Cancel closes the form and keeps the cell'
        assert ui.markButton.visible && asked.size()==1 && asked[0].is(ui.markButton) : 'the button is back, and the focus with it: '+asked
        focusOn=null
        // the run's own answer is kept: with its table taken away, the checks are still the run's (as the page takes them)
        assert new File(work,'Review_call_checks_'+A+'.csv').delete() : 'the checks table of the run can be taken away'
        ui.loadResults(dirA);assert ui.checks?.note==unchecked : 'the checks of the last run here are kept, not only the saved table'
    }

    // ---- an analysis with flags, and a control window beside its stimulus
    File dirB=copy(B)
    def table=Csv.read(new File(dirB,'responder_calls.csv')),header=table[0]
    Csv.write(new File(dirB,'responder_calls.csv'),table+table.drop(1).collect{List line->def row=new ArrayList(line);['Role':'Control','Stimulus':'Ionomycin','Threshold_increase':'0.4','Floor_threshold':'0.4','Baseline_mean':'0.5'].each{k,v->row[header.indexOf(k)]=v};row})
    def countsTable=Csv.read(new File(dirB,'responder_counts.csv'))   // and a row of counts for the control window
    Csv.write(new File(dirB,'responder_counts.csv'),countsTable+countsTable.drop(1).collect{List line->def row=new ArrayList(line);['Role':'Control','Stimulus':'Ionomycin'].each{k,v->row[countsTable[0].indexOf(k)]=v};row})
    File checksB=new File(work,'Review_call_checks_'+B+'.csv'),reviewB=new File(work,'Review_reviewed_calls_'+B+'.csv'),percentB=new File(work,'Review_percent_responding_'+B+'.csv')
    def row={String cell,String call,String flag,String rule,Object first,Object last,String detail,String check,String why->[cell,call,flag,rule,first,last,detail,8,400,'yes','yes',check,why,'','',Checks.codesOf(rule)]}
    def none={String cell,String call->row(cell,call,'','','','','','not checked','no call')}
    // first: a field that moved and short separate crossings on one cell, with the first frame of its event
    Csv.write(checksB,[Checks.HEADER,row('Cell001','Responded','','','','','','checked',''),row('Cell002','Did not respond','','','','','','checked',''),none('Cell003','Excluded'),none('Cell004','Excluded'),
        row('Cell005','Responded','field moved, frames 20 to 22; short separate crossings','field moved; short separate crossings',18,23,'the field moved 3.1 px between frames 20 and 22; otherwise possible passing object: seen in the near ring first; 3 crossings, the longest 2 samples','checked',''),none('Cell006','Unscorable')])
    Settings.percentResponding(dirB,ui.measurement)   // the percent table beside the tables, as a run leaves it
    byte[] percentBefore=percentB.bytes
    onUI {
        ui.loadResults(dirB);ui.cells.selectedItem='Cell005'
        assert ui.flagRow.visible && flags()*.text==['field moved, frames 20 to 22','·','short separate crossings'] : flags()*.text
        assert shownRows()==['call','flags','mark'] && stacked() && ui.flagRow.y==ui.cellStatus.y+ui.cellStatus.height : 'the flags right under the call line, the button right under them: '+ui.cellStatus.bounds+' '+ui.flagRow.bounds+' '+ui.markRow.bounds
        def moved=flags()[0]
        assert moved.foreground==Theme.BLUE && moved.cursor.type==java.awt.Cursor.HAND_CURSOR && moved.toolTipText.startsWith('<html>Go to frame 20<br>the field moved 3.1 px') : 'a flag links to the first frame its text names, not the event\'s: '+moved.toolTipText
        def blips=flags()[2]
        assert blips.foreground==Theme.BLUE && blips.toolTipText.startsWith('<html>Go to frame 18<br>') : 'a flag whose text names no frame: the first frame of the event, '+blips.toolTipText
        ui.tabs.selectedIndex=0;ui.showFrame(30);click(moved);assert ui.tabs.selectedIndex==1 && ui.frame==20 : 'the click goes there'
        click(blips);assert ui.frame==18
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
        // without marks the calls and counts tabs are the tables as saved, column for column
        assert columns(ui.calls)==Csv.read(new File(dirB,'responder_calls.csv'))[0] && columns(ui.counts)==Csv.read(new File(dirB,'responder_counts.csv'))[0] : 'the tabs as saved'
        assert ui.calls.model.is(ui.callsSaved) && ui.counts.model.is(ui.countsSaved) : 'without marks the tabs show the models the tables were loaded into'
        before.callColumns=columns(ui.calls);before.countColumns=columns(ui.counts)
        // the caption and the rule line are under the plot, not on it, and the rule line gives the required rise with the
        // decimals of the call line ("needed 0.10"), and leaves out the span and the samples the table does not hold; nothing
        // but the axes and the key is drawn on the plot
        assert ui.plot.noteLines*.text==[before.plot,'Required rise 0.10 (Floor)'] && ui.plot.ruleLine()=='Required rise 0.10 (Floor)' : ui.plot.noteLines*.text
        layOut();List drawn=drawnOnPlot()
        assert ui.plot.notes.parent.is(ui.plot.parent) && ui.plot.notes.y>=ui.plot.y+ui.plot.height && ui.plot.notes.height>=2*ui.plot.noteLines[0].preferredSize.height : 'the two lines are under the plot: '+ui.plot.notes.bounds
        // on screen the key holds what it held before the review controls: the cell, the mean and the traces behind it; no band
        // or dashed line is named there, and the frame of the plot starts where it did, at the window's width and as narrow as
        // beside a real recording's picture
        assert drawn.contains('Time (s)') && drawn.contains('mean ± SEM') && drawn.every{allowedOnScreen(it)} : 'drawn on the plot on screen: '+drawn.findAll{!allowedOnScreen(it)}
        assert drawn.findAll{it in ['Cell001','mean ± SEM'] || counted(it)}==['Cell001','mean ± SEM','responded (2)','did not respond (1)','no call (1)'] : 'the key on screen: '+drawn
        int top=frameTop(),wide=ui.plot.width;assert top==mainTop(wide) : 'the frame starts where it did before the review controls: '+top+' '+mainTop(wide)
        top=narrowTop();assert top==mainTop(245) : 'and as narrow as beside a real recording\'s picture: '+top+' '+mainTop(245)
        // pointing at a band names it, in the words of trace.png's key
        layOut();drawnOnPlot()
        def over=new java.awt.event.MouseEvent(ui.plot,java.awt.event.MouseEvent.MOUSE_MOVED,0L,0,ui.plot.X(ui.plot.stimulusEnd)-2,ui.plot.height-100,0,false)
        assert ui.plot.getToolTipText(over)=='<html>response window<br>'+ui.plot.HINT+'</html>' : 'the band is named in the words of the key: '+ui.plot.getToolTipText(over)
        // trace.png's key names the bands and the dashed lines in the words of the Step 4 page's legend: here the required rise,
        // and the noise threshold at the baseline mean (the minimum rise is the required rise, under it)
        List figured=drawnInFigure()
        assert figured.every{allowedInFigure(it)} && ['required rise','noise threshold','F₀ frames','baseline','response window','response'].every{figured.contains(it)} && !figured.contains('minimum rise') && !figured.contains('baseline mean') : 'trace.png\'s key: '+figured
        // in the order of the page's legend: the cell, the mean, the traces behind it, the dashed lines, then the bands
        List keyed=figured.findAll{it in ['Cell001','mean ± SEM']+keyWords || counted(it)}
        assert keyed==['Cell001','mean ± SEM','responded (2)','did not respond (1)','no call (1)','required rise','noise threshold','F₀ frames','baseline','response window','response'] : 'the order of trace.png\'s key: '+keyed
        assert !ui.plot.figure && drawnOnPlot().every{allowedOnScreen(it)} : 'the plot on screen keeps its own key'
        // a band has a filled box of its colour; a dashed line a dashed line (its colour, then a gap); the mean an unbroken line
        def box=swatch('response window'),dash=swatch('required rise'),solid=swatch('mean ± SEM')
        assert box.row[2..14].every{it==0xFFF0DC} && box.above==0xFFF0DC : 'the response window: a filled box '+box
        assert dash.row[1..3].every{it==(Theme.GOOD.RGB&0xFFFFFF)} && dash.row[6..8].every{it==0xFFFFFF} && dash.above==0xFFFFFF : 'the required rise: a dashed line '+dash
        assert solid.row[1..15].every{it==(Theme.BLUE.RGB&0xFFFFFF)} && solid.above==0xFFFFFF : 'the mean: a line '+solid
        // the flag of the selected cell, in the checks' own words, and its link to the first frame of the event
        def flag=flags()[0]
        assert ui.flagRow.visible && flags().size()==1 && flag.text==passing && flag.foreground==Theme.BLUE && flag.cursor.type==java.awt.Cursor.HAND_CURSOR : 'the flag in its own words, a link'
        assert flag.toolTipText=='<html>Go to frame 12<br>seen in the near ring 2 frames before the cell</html>' : flag.toolTipText
        ui.view.selectedIndex=0;ui.tabs.selectedIndex=0;ui.showFrame(25)
        click(flag)
        assert ui.tabs.selectedIndex==1 && ui.view.selectedIndex==2 && ui.frame==12 && ui.slider.value==12 && ui.cellPicture.image.is(ui.frameImage) : 'a click shows the recording at the first frame of the event'
        ui.view.selectedIndex=1;ui.showFrame(25);ui.playButton.doClick();assert ui.player.running : 'the movie plays'
        click(flag);assert !ui.player.running && ui.view.selectedIndex==1 && ui.frame==12 : 'a click stops playing; a movie view stays as it was chosen'
        // a flag without a first frame is plain text; a cell that could not be checked says so, in the recording whose rings were measured
        ui.cells.selectedItem='Cell005'
        assert flags()*.text==['short separate crossings','·','passing object not checked: no ring around the cell'] : flags()*.text
        def plain=flags()[0]
        assert plain.foreground==Theme.MUTED && plain.cursor.type!=java.awt.Cursor.HAND_CURSOR && plain.toolTipText=='3 crossings, the longest 2 samples' : 'without a first frame a flag is plain text'
        ui.markButton.doClick();assert ui.reasonBox.selectedItem=='Noise' : 'short separate crossings: Noise is chosen';ui.cancelMark.doClick()
        // a cell that did not respond: the other button, and no reason chosen even with a flag
        ui.cells.selectedItem='Cell002'
        assert flags()*.text==[edge] && ui.markButton.visible && ui.markButton.text=='Mark as a responder' : 'a cell that did not respond can be marked as a responder'
        ui.markButton.doClick();assert ui.reasonBox.selectedIndex==-1 && !ui.saveMark.enabled : 'never Passing object chosen on a cell that did not respond'
        def nothing=snapshot();enter(ui.reasonBox);assert ui.formCell=='Cell002' && snapshot()==nothing && said()=='' : 'Enter does nothing before a reason is chosen: '+said()
        ui.cancelMark.doClick()
        // no call to review: no button, no flags
        ['Cell003','Cell004','Cell006'].each{ui.cells.selectedItem=it;assert !ui.markButton.visible && !ui.removeMark.visible && flags().isEmpty() && said()=='' && shownRows()==['call'] : it+' has no call to review: '+shownRows()}
        // the head over the picture and the plot has one height for every cell: flagged or not, with a button or none, a
        // call line of one line or more; so the picture and the plot keep their size, and so does a saved trace.png
        Map sizes=['Cell001','Cell002','Cell003','Cell005','Cell006'].collectEntries{ui.cells.selectedItem=it;[it,heights()]}
        assert sizes.values().toList().unique().size()==1 : 'one height for every cell: '+sizes
        ui.cells.selectedItem='Cell003';paint('fiji_cells_in_a_row_1')
        ui.cells.selectedItem='Cell005';paint('fiji_cells_in_a_row_2')
        // a cell without a trace: nothing is written on the plot for it; why is said under the plot, and there is no rule line
        ui.cells.selectedItem='Cell003';layOut()
        assert drawnOnPlot().every{allowedOnScreen(it)} : 'drawn on the plot for an excluded cell: '+drawnOnPlot().findAll{!allowedOnScreen(it)}
        assert ui.plot.noteLines*.text==['Cell003 · unavailable: inspect Cell quality for the exclusion reason',''] : ui.plot.noteLines*.text
        ui.cells.selectedItem='Cell004';assert ui.plot.noteLines[1].text=='' : 'no rule line for an excluded cell: '+ui.plot.noteLines[1].text
        Map rowWas=ui.plot.responseCall;ui.plot.responseCall=[Threshold_increase:'0.2',Threshold_source:'Floor',Run_duration_s:'1',Run_samples:'2'];ui.plot.changed()
        assert ui.plot.selected.length==0 && ui.plot.noteLines[1].text=='' : 'a cell without a trace has no rule line, whatever its row holds: '+ui.plot.noteLines[1].text
        ui.plot.responseCall=rowWas;ui.plot.changed()
        // an unscorable cell with a trace but no required rise has none either; a cell that did not respond has its required rise
        ui.cells.selectedItem='Cell006';assert ui.plot.selected.length>0 && ui.plot.noteLines[1].text=='' : 'no rule line without a required rise: '+ui.plot.noteLines[1].text
        ui.cells.selectedItem='Cell002';assert ui.plot.noteLines[1].text=='Required rise 0.10 (Floor)' : ui.plot.noteLines[1].text
        // the key keeps its room on a cell without a trace or without a run: the plot keeps its place from cell to cell, on
        // screen and in trace.png, and in trace.png a band of the window is named on every cell
        Map tops=names.collectEntries{ui.cells.selectedItem=it;[it,placeOf()]}
        assert (tops.values() as Set).size()==1 : 'the plot keeps its place: '+tops
        // trace.png is saved at the width of the plot on screen (650 px at least), and its key wraps at that width: at every
        // width the frame of each cell starts where it does for the others, with room for the dashed lines of every cell
        def figureTopAt={int width->java.awt.Dimension was=ui.plot.size;ui.plot.setSize(width,was.height as int);ui.plot.figure=true
            try{drawnOnPlot();ui.plot.@t}finally{ui.plot.figure=false;ui.plot.setSize(was);ui.plot.changed()}}
        List savedWidths=(650..1300).step(5)
        Map figureTops=names.collectEntries{ui.cells.selectedItem=it;layOut();[it,savedWidths.collect{figureTopAt(it as int)}]}
        assert (figureTops.values() as Set).size()==1 : 'the plot keeps its place in trace.png at every width saved: '+savedWidths.withIndex().findAll{w,int i->(figureTops.values()*.getAt(i) as Set).size()>1}.collect{w,int i->[w,figureTops.collectEntries{k,v->[k,v[i]]}]}.take(3)
        assert (figureTops.values()[0] as Set).size()>1 : 'the key of trace.png takes more lines at some widths saved than at others: '+(figureTops.values()[0] as Set)
        ui.cells.selectedItem='Cell002';layOut();assert drawnInFigure().containsAll(['response','required rise','Cell002']) : 'a cell without a run: the band is named all the same: '+drawnInFigure()
        ui.cells.selectedItem='Cell003';layOut();List bare=drawnInFigure()
        assert bare.containsAll(['F₀ frames','baseline','response window','response']) && !bare.contains('Cell003') && !bare.contains('required rise') : 'a cell without a trace: '+bare
        assert drawnOnPlot().every{allowedOnScreen(it)} && !drawnOnPlot().contains('Cell003') : 'on screen no band is named for it either: '+drawnOnPlot()
        // on screen the key is laid as before the review controls, from the frame's left edge (62 px) across the frame's width
        // (to 18 px from the right), here as narrow as beside a real recording's picture; trace.png's runs across the whole
        // picture from 8 px in. Each line goes on to the next only when the next entry would not fit; a name longer than the
        // key's width is cut short
        def keyAt={boolean saved->layOut();if(!saved)ui.plot.setSize(245,ui.plot.height)
            def look={->def image=new java.awt.image.BufferedImage(ui.plot.width,ui.plot.height,java.awt.image.BufferedImage.TYPE_INT_RGB)
                def spy=new TextSpy(real:image.createGraphics()),fm;try{ui.plot.still(spy);fm=spy.real.fontMetrics}finally{spy.real.dispose()}
                [ui.plot.width,spy.at.findAll{it[0] in ['mean ± SEM',ui.plot.cell]+keyWords || counted(it[0]) || it[0].startsWith('Cell_')}.collect{it+[fm.stringWidth(it[0])]}]}
            def found=saved?asSaved(look):look();layOut();found}
        def keyLinesOf={boolean saved->def found=keyAt(saved);int width=found[0];List entries=found[1]
            int start=saved?8+21:62+21,right=saved?width-8:width-18,across=saved?width-16:width-80
            List keyLines=entries.groupBy{it[2]}.sort{it.key}.values().collect{it.sort{a,b->a[1]<=>b[1]}}
            String where=saved?'trace.png':'on screen'
            assert keyLines.every{it[0][1]==start} && entries.every{(it[1] as int)+(it[3] as int)<=right} : where+': each line of the key starts at '+start+' px, and nothing runs past '+right+' px: '+entries
            assert (1..<keyLines.size()).every{int i->def last=keyLines[i-1][-1],next=keyLines[i][0];((last[1] as int)-start)+21+(last[3] as int)+14+21+(next[3] as int)+14>across} : where+': a line of the key goes on to the next only when the next entry would not fit: '+keyLines.collect{it*.getAt(0)}
            [keyLines,start,right]}
        ui.cells.selectedItem='Cell001'
        assert keyLinesOf(false)[0].size()>1 && keyLinesOf(true)[0].size()>1 : 'the key takes more than one line on screen and in trace.png'
        String shortName=ui.plot.cell;ui.plot.cell='Cell_ac975e8b6d7dd113b515892b56003127';ui.plot.changed()
        def (List cutLines,int cutStart,int cutRight)=keyLinesOf(false);def cut=cutLines.collectMany{it}.find{it[0].startsWith('Cell_')}
        assert cut!=null && cut[0].endsWith('…') && cut[1]==cutStart && (cut[1] as int)+(cut[3] as int)<=cutRight : 'a long cell name is cut short in the key: '+cut
        ui.plot.cell=shortName;ui.plot.changed();layOut()
        // a refusal of two lines, as one naming the percent table of a recording with a usual name, takes the room kept for it
        ui.cells.selectedItem='Cell001'
        ui.say('Close Review_percent_responding_2026-11-18_50B11_ATP_100uM_dish2_SimpleCa-20261118-142233-1a2b3c4d.csv if it is open in another program, then try again.','Cell001')
        assert heights()==sizes.Cell001 : 'a refusal of two lines: '+heights()+' '+sizes.Cell001
        ui.say(null,null);ui.cells.selectedItem='Cell005'
        // the new controls have the shape and size of the ones beside them
        layOut()
        [ui.markButton,ui.removeMark,ui.saveMark,ui.cancelMark].each{assert it.preferredSize.height==ui.previousCell.preferredSize.height && it.UI.class.name=='CalciumButtonUI' : 'a mark button is a button of this window: '+it.text}
        assert ui.saveMark.UI.primary && !ui.cancelMark.UI.primary && !ui.markButton.UI.primary : 'Save is the filled one'
        assert ui.reasonBox.preferredSize.height==ui.group.preferredSize.height && ui.reasonBox.font==ui.group.font && ui.reasonBox.font==Theme.font() : 'the Reason list is a list like Show'
        assert ui.noteField.preferredSize.height==ui.fields.sampling_hz.preferredSize.height && ui.noteField.font==ui.fields.sampling_hz.font : 'the Note box is a box like the settings'
        assert ui.reasonLabel.font==Theme.font() && ui.noteLabel.font==Theme.font() && ui.reasonLabel.text=='Reason' && ui.noteLabel.text=='Note' : 'the form\'s labels: '+ui.reasonLabel.text+', '+ui.noteLabel.text

        // ---- the form on a flagged responder: Passing object chosen; the refusals in the form, in CalciumReview's words
        ui.cells.selectedItem='Cell001';asked.clear();ui.markButton.doClick()
        assert ui.formCell=='Cell001' && ui.reasonBox.selectedItem=='Passing object' && ui.saveMark.enabled && ui.noteField.hint=='Optional' && ui.noteField.text=='' : 'a cell that may have a passing object: Passing object is chosen'
        assert asked.size()==1 && asked[0].is(ui.noteField) : 'with the reason chosen, the focus goes to the Note box: '+asked
        Map formOpen=heights()
        def refused={String reason,String note,String message->def files=snapshot();ui.reasonBox.selectedItem=reason;ui.noteField.text=note;ui.saveMark.doClick()
            assert said()==message : said()
            assert ui.formCell=='Cell001' && snapshot()==files : 'nothing is written: '+message}
        refused('Other','   ','Other needs a short note.')
        assert ui.noteField.hint=='' : 'the note is not optional for Other'
        assert shownRows()==['call','flags','mark','refusal'] && stacked() : 'the refusal right under the form: '+ui.markRow.bounds+' '+ui.markProblem.bounds
        assert heights()==formOpen && formOpen==sizes.Cell001 : 'the form and a refusal under it take the room kept for them: '+heights()+' '+formOpen+' '+sizes.Cell001
        refused('Noise','=SUM(A1)','The note cannot start with =, +, - or @.')
        refused('Noise','  -1 frame','The note cannot start with =, +, - or @.')   // the spaces at either end go first, as in the review file
        refused('Noise','one\ttwo','The note must be one line.')
        ui.noteField.text='x'*130;assert ui.noteField.text.length()==120 : 'the note box takes at most 120 characters'
        // the review's own refusals: a review file that is not this analysis's, a percent table open in another program, a run going on
        reviewB.text=Review.HEADER.join(',')+'\nCell002,Did not respond,Responded,Noise,,2026-11-01T10:00:00Z,another,x,,,manual\n'
        refused('Passing object','','The review file Review_reviewed_calls_'+B+'.csv is for another analysis, so its marks are not used. Move it out of the folder to review this analysis again.')
        assert reviewB.delete() : 'the other analysis\'s review file is taken away'
        assert percentB.setWritable(false) : 'the percent table is locked'
        try{refused('Passing object','','Close Review_percent_responding_'+B+'.csv if it is open in another program, then try again.')}finally{percentB.setWritable(true)}
        ui.busy=true;try{refused('Passing object','','Wait for the analysis to finish.');ui.showReview('Cell001');assert !ui.markButton.visible && !ui.removeMark.visible && ui.formCell==null : 'no mark controls while an analysis runs'}finally{ui.busy=false}
        ui.showReview('Cell001');assert ui.markButton.visible && said()=='' : 'and they come back'
        // the form closes when the cell changes; the arrows do not change it while the form is open; Esc closes it, and the
        // focus goes back to the button when it was in the form
        ui.markButton.doClick();ui.noteField.text='=x';ui.saveMark.doClick();assert said()!='' : 'a refused note is said in the form'
        ui.cells.selectedItem='Cell005';assert ui.formCell==null && !ui.reasonBox.visible && said()=='' : 'the form and its refusal go when the cell changes'
        ui.cells.selectedItem='Cell001';ui.markButton.doClick()
        ui.rootPane.getActionForKeyStroke(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT,0)).actionPerformed(new java.awt.event.ActionEvent(ui.rootPane,0,''))
        assert ui.cells.selectedItem=='Cell001' && ui.formCell=='Cell001' : 'the arrows do not step to another cell while the form is open'
        focusOn=ui.noteField;asked.clear()
        assert press(ui.noteField,KeyEvent.VK_ESCAPE) && ui.formCell==null && ui.cells.selectedItem=='Cell001' : 'Esc closes the form and keeps the cell'
        assert ui.markButton.visible && asked.size()==1 && asked[0].is(ui.markButton) : 'Esc: the focus goes back to the button: '+asked
        ui.markButton.doClick();focusOn=ui.reasonBox;asked.clear()
        assert press(ui.reasonBox,KeyEvent.VK_ESCAPE) && ui.formCell==null && asked.size()==1 && asked[0].is(ui.markButton) : 'Esc closes the form from the Reason list too'
        ui.markButton.doClick();focusOn=ui.cells;asked.clear();ui.cancelMark.doClick()
        assert ui.formCell==null && ui.markButton.visible && asked.isEmpty() : 'the focus is left where it was when it was not in the form: '+asked
        focusOn=null
        assert !reviewB.exists() && percentB.bytes==percentBefore : 'no mark was saved'

        // ---- a review file for another analysis is not applied: no mark controls, and every cell of the stimulus says why and
        // what to do, in the words a refused mark is given; the control window does not
        reviewB.text=Review.HEADER.join(',')+'\nCell002,Did not respond,Responded,Noise,,2026-11-01T10:00:00Z,another,x,,,manual\n'
        String other='The review file Review_reviewed_calls_'+B+'.csv is for another analysis, so its marks are not used.'
        ui.loadResults(dirB);ui.showGroup(0)
        List kept=['Cell001','Cell002','Cell003'].collect{ui.cells.selectedItem=it
            assert !ui.markButton.visible && !ui.removeMark.visible && ui.formCell==null && said()==other+' Move it out of the folder to review this analysis again.' : it+': '+said()
            assert !ui.markRow.visible && shownRows()[-1]=='refusal' && stacked() : 'the sentence right under the call line and its flags: '+shownRows()
            heights()}
        assert kept.unique().size()==1 : 'the sentence keeps the head one height from cell to cell: '+kept
        ui.cells.selectedItem='Cell001';ui.markButton.doClick();assert ui.formCell==null : 'no form either'
        assert ui.responding.toolTipText==other && ui.responding.value.text=='66.7 %' : ui.responding.toolTipText
        paint('fiji_review_other_analysis')
        ui.stimuli.selectedItem='Ionomycin';ui.cells.selectedItem='Cell001';assert said()=='' && !ui.markButton.visible : 'the control window is never reviewed'
        assert heights()==kept[0] : 'the head keeps its height on the control window, without the sentence: '+heights()+' '+kept[0]
        assert reviewB.delete();ui.loadResults(dirB);ui.showGroup(0);ui.cells.selectedItem='Cell001'
        assert ui.markButton.visible && said()=='' : 'with the file moved away the cell can be marked again'

        // ---- a mark saved: every view follows it at once, the cell stays selected and the time cursor stays
        ui.showGroup(1);assert ui.cells.selectedItem=='Cell001' : 'Cell001 is a responder'
        paint('fiji_traces_one_line')
        Map pickersBefore=pickers()
        ui.markButton.doClick();ui.noteField.text='  drifts in from the left  '
        paint('fiji_traces_form')
        ui.showFrame(17);double moment=ui.plot.moment
        Map tablesBefore=['responder_calls.csv','responder_counts.csv'].collectEntries{[it,new File(dirB,it).bytes]}
        // the window's traces are used for every mark: the traces table is not read again (here it cannot be read)
        Closure unlockB=lockAway(new File(dirB,'normalized_traces.csv'))
        focusOn=ui.noteField;asked.clear()
        try{ui.saveMark.doClick()}finally{unlockB()}
        assert said()=='' && ui.formCell==null : said()
        assert asked.size()==1 && asked[0].is(ui.removeMark) && ui.removeMark.visible : 'after Save the focus goes to Remove mark: '+asked
        focusOn=null
        Map mark=Review.read(dirB,ui.measurement).marks.Cell001
        assert mark.automatic=='Responded' && mark.reviewed=='Did not respond' && mark.reason=='Passing object' && mark.note=='drifts in from the left' : mark
        def saved=Csv.read(reviewB)
        assert saved[0]==Review.HEADER_BEFORE && saved.size()==2 && saved[1][0..6]==['Cell001','Responded','Did not respond','Passing object','drifts in from the left',mark.time,B] && saved[1][7]==Installer.digest(new File(dirB,'responder_calls.csv').bytes) : saved
        assert Csv.read(percentB)[0]==Settings.PERCENT_COLUMNS+Settings.REVIEW_COLUMNS : 'the percent table gives the figures after review beside the others'
        assert ui.responding.value.text=='33.3 %' && strip(ui.responding.caption.text)=='Capsaicin: 1 of 3 evaluable cells after review; 2 before review' : ui.responding.caption.text
        assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==['All cells (6)','Responded after review (1)','Did not respond after review (2)','Excluded or unscorable (3)','Marked (1)'] : (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}
        // Show stays on the responders: the cell that left them stays selected, and in the Cell list in its place while it is,
        // but the group drawn is the responders after review
        assert ui.group.selectedIndex==1 && ui.group.selectedItem=='Responded after review (1)' && ui.cells.selectedItem=='Cell001' : 'Show is kept: '+ui.group.selectedItem+', '+ui.cells.selectedItem
        assert (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell001','Cell005'] && ui.nextCell.enabled : 'the marked cell stays in the Cell list, in its place'
        assert ui.plot.allCalls==['Responded'] && ui.plot.all.length==1 && ui.cellPicture.group==1 && ui.cellPicture.selected==0 : 'the group drawn is the responders: '+ui.plot.allCalls
        assert ui.frame==17 && ui.plot.moment==moment : 'the time cursor did not jump'
        assert ui.cellPicture.calls.Cell001=='Did not respond' && ui.cellPicture.calls.Cell005=='Responded' : ui.cellPicture.calls
        assert ui.heatmap.calls[ui.heatmap.names.indexOf('Cell001')]=='Did not respond' && ui.heatmap.names.indexOf('Cell001')>ui.heatmap.names.indexOf('Cell005') : 'the heatmap orders and colours by the call after review'
        assert strip(ui.cellStatus.text)=='● Cell001: Did not respond after review · passing object: drifts in from the left · automatic call Responded · sustained threshold crossing · largest rise 1.00 at 6 s (needed 0.10)' : strip(ui.cellStatus.text)
        assert ui.cellStatus.text.contains("<font color='"+Theme.hex(Calls.NOT)+"'>●</font>") : 'the dot has the colour of the call after review'
        assert ui.cellStatus.toolTipText=='<html>Marked '+minute(mark.time)+'<br>'+before.tip+'</html>' : 'when it was marked, to the minute as on the page: '+ui.cellStatus.toolTipText
        // the caption gives the call after review without the engine's reason for the automatic call; it is under the plot
        assert ui.plot.caption=='Cell001 · Did not respond' && ui.plot.noteLines[0].text==ui.plot.caption : ui.plot.caption
        assert ui.removeMark.visible && !ui.markButton.visible && flags().size()==1 : 'a marked cell has one button: Remove mark'
        assert heights()==sizes.Cell001 : 'a mark moves neither the head, the tabs nor the plot: '+heights()+' '+sizes.Cell001
        assert pickers()==pickersBefore : 'Show keeps its width at a mark, and the Cell box and the arrows their place: '+pickers()+' '+pickersBefore
        assert shownRows()==['call','flags','mark'] && stacked() : 'Remove mark right under the flags: '+shownRows()
        // the calls and counts tabs give the reviewed call beside the automatic one, and % responding after review, on screen
        // only: the saved tables are as they were
        assert columns(ui.calls)==Csv.read(new File(dirB,'responder_calls.csv'))[0]+['Reviewed call','Review reason'] : columns(ui.calls)
        assert rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])[-2..-1]==['Did not respond','Passing object: drifts in from the left'] : rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])
        assert rowOf(ui.calls,[Cell:'Cell001',Role:'Control'])[-2..-1]==['',''] && rowOf(ui.calls,[Cell:'Cell005',Role:'Stimulus'])[-2..-1]==['',''] : 'only the marked cell\'s row of the stimulus'
        assert rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])[columns(ui.calls).indexOf('Status')]=='Responded' : 'the automatic call stays in its column'
        assert columns(ui.counts)==Csv.read(new File(dirB,'responder_counts.csv'))[0]+['After review','Calls changed'] && rowOf(ui.counts,[Role:'Stimulus'])[-2..-1]==['33.3% (1 of 3)','1'] : rowOf(ui.counts,[Role:'Stimulus'])
        assert rowOf(ui.counts,[Role:'Control'])[-2..-1]==['',''] : 'After review is on the stimulus row only: '+rowOf(ui.counts,[Role:'Control'])
        assert tablesBefore.every{name,bytes->new File(dirB,name).bytes==bytes} : 'the saved tables are not written'
        paint('fiji_traces_marked')
        [['Cell response calls','fiji_calls_tab',ui.calls],['Responder counts','fiji_counts_tab',ui.counts]].each{String tab,String name,JTable shown->
            // a table puts its header above itself when its window is shown; this one never is, so the picture does it
            def scroll=SwingUtilities.getAncestorOfClass(JScrollPane,shown) as JScrollPane;if(scroll.columnHeader?.view==null)scroll.columnHeaderView=shown.tableHeader
            ui.tabs.selectedIndex=ui.tabs.indexOfTab(tab);layOut();def port=shown.parent as JViewport
            port.viewPosition=new java.awt.Point(Math.max(0,shown.width-port.width),0);paint(name,false)}
        ui.tabs.selectedIndex=1
        // a picture of the Cell list as it drops down, at twice its size: the marked cell in its place among the responders
        def popup=ui.cells.UI.getAccessibleChild(ui.cells,0),drop=popup.list as JList;drop.size=drop.preferredSize
        def dropped=new java.awt.image.BufferedImage(Math.max(1,2*drop.width),Math.max(1,2*drop.height),java.awt.image.BufferedImage.TYPE_INT_RGB);def dg=dropped.createGraphics();try{dg.scale(2d,2d);drop.paint(dg)}finally{dg.dispose()}
        javax.imageio.ImageIO.write(dropped,'png',new File(pictures,'fiji_cells_list_after_mark.png'))
        // the next arrow goes on through the responders from the marked cell's place; the marked cell then leaves the list
        ui.rootPane.getActionForKeyStroke(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT,0)).actionPerformed(new java.awt.event.ActionEvent(ui.rootPane,0,''))
        assert ui.cells.selectedItem=='Cell005' && ui.group.selectedIndex==1 && (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell005'] && !ui.nextCell.enabled : 'the arrow goes to the next responder: '+ui.cells.selectedItem+' '+(0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}
        paint('fiji_traces_marked_next')
        // a mark that empties the group: the cell marked stays shown without the group's traces, and the plot keeps its place,
        // as it does when the mark is removed
        List placed=placeOf();ui.markButton.doClick();ui.reasonBox.selectedItem='Noise';ui.saveMark.doClick()
        assert ui.group.selectedItem=='Responded after review (0)' && ui.cells.selectedItem=='Cell005' && ui.plot.all.length==0 && placeOf()==placed : 'the plot keeps its place when a mark empties the group: '+placed+' '+placeOf()
        paint('fiji_traces_group_emptied')
        ui.removeMark.doClick();assert ui.group.selectedItem=='Responded after review (1)' && placeOf()==placed : 'and when the mark is removed: '+placed+' '+placeOf()
        ui.showGroup(0);ui.cells.selectedItem='Cell001';layOut()
        assert ui.plot.allCalls.count{it=='Responded'}==1 && ui.plot.allCalls.count{it=='Did not respond'}==2 : 'all cells, by the call after review: '+ui.plot.allCalls
        // the long call line and the tile's caption wrap at the width they are given: nothing is cut off
        [ui.cellStatus,ui.responding.caption].each{assert it.width>=it.preferredSize.width && it.height>=it.preferredSize.height : 'cut off: '+strip(it.text)+' '+it.size+' < '+it.preferredSize}
        assert ui.responding.caption.y+ui.responding.caption.height<=ui.responding.height-ui.responding.insets.bottom : 'the tile holds its caption: '+ui.responding.caption.bounds+' in '+ui.responding.size
        // the Marked entry lists the marked cells, on the picture too
        ui.showGroup(4);assert (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell001'] && ui.cellPicture.only==['Cell001'] as Set && ui.plot.all.length==1 : 'the Marked entry lists and draws only the marked cell'
        ui.showGroup(0);assert ui.cellPicture.only==null : 'All cells draws every cell'
        // a control window keeps its automatic call: no marks, no flags, its own groups
        ui.stimuli.selectedItem='Ionomycin'
        assert !ui.markButton.visible && !ui.removeMark.visible && flags().isEmpty() && (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==['All cells (6)','Responded (2)','Did not respond (1)','Excluded or unscorable (3)'] : 'the control window is never reviewed'
        assert strip(ui.cellStatus.text).startsWith('● Cell001: Responded · ') && ui.plot.caption.startsWith('Cell001 · Responded · ') : 'the control window shows its automatic call: '+strip(ui.cellStatus.text)
        assert pickers()==pickersBefore : 'Show has the same width in the control window: '+pickers()+' '+pickersBefore
        // its saved plots follow no review, and their record says nothing of one
        layOut();ui.exportPlots()
        def onControl=new groovy.json.JsonSlurper().parse(new File(dirB,'plot_provenance.json'))
        assert onControl.stimulus=='Ionomycin' && !onControl.containsKey('calls') && !onControl.containsKey('review_file') && !onControl.containsKey('review_sha256') : onControl
        ui.stimuli.selectedItem='Capsaicin';ui.cells.selectedItem='Cell001';assert ui.removeMark.visible : 'back on the stimulus, the mark can be removed'
        // saved plots say their groups follow the review marks, and which file the marks are in
        layOut();ui.exportPlots()
        assert !ui.plot.figure && drawnOnPlot().every{allowedOnScreen(it)} : 'after Save plots the key on screen is as before: '+drawnOnPlot().findAll{!allowedOnScreen(it)}
        def record=new groovy.json.JsonSlurper().parse(new File(dirB,'plot_provenance.json'))
        assert record.calls=='after review' && record.review_file==reviewB.name && record.review_sha256==Installer.digest(reviewB.bytes) : ui.status.text
        // the saved trace.png has no caption, rule line, time cursor or anything else written on it: with other words under
        // the plot and the cursor elsewhere it is the same picture
        File traceB=new File(dirB,'trace.png');def first=javax.imageio.ImageIO.read(traceB)
        java.nio.file.Files.copy(traceB.toPath(),new File(pictures,'fiji_trace_export.png').toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        // drawn across the whole picture: the plot on screen is narrower than the 650 px saved, and no column of the picture is
        // left black; the same for the heatmap; and both keep their size on screen
        def heatPicture=javax.imageio.ImageIO.read(new File(dirB,'heatmap.png'))
        java.nio.file.Files.copy(new File(dirB,'heatmap.png').toPath(),new File(pictures,'fiji_heatmap_export.png').toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        def blackColumn={def image->(0..<image.width).find{int x->(0..<image.height).every{int y->(image.getRGB(x,y)&0xFFFFFF)==0}}}
        assert ui.plot.width<650 && first.width==650 && blackColumn(first)==null && (0..<first.height).every{int y->first.getRGB(first.width-1,y)==(0xFFFFFFFF as int)} : 'trace.png is drawn to its right edge: '+ui.plot.width+' '+blackColumn(first)
        // the plot is drawn at the width saved: its frame runs on past where it ends on screen, to 18 px from the right edge
        int onScreen=ui.plot.width
        assert ((onScreen-15)..<645).any{int x->(0..<first.height).count{int y->(first.getRGB(x,y)&0xFFFFFF)!=0xFFFFFF}>100} : 'the plot is drawn at the width saved'
        assert heatPicture.width>=650 && blackColumn(heatPicture)==null : 'heatmap.png is drawn to its right edge: '+blackColumn(heatPicture)
        // trace.png is the plot with trace.png's key, which names its bands and dashed lines, pixel for pixel
        def asKeyed=savedPicture(first.width,first.height)
        assert (0..<first.width).every{int x->(0..<first.height).every{int y->first.getRGB(x,y)==asKeyed.getRGB(x,y)}} : 'trace.png has the key of a saved figure'
        layOut();int plotWidth=ui.plot.width;ui.exportPlots();assert ui.plot.width==plotWidth : 'the plot keeps its size on screen'
        // the rule line under the plot gives the span and the samples with the call line's decimals, not as the table holds them
        Map call=ui.plot.responseCall;ui.plot.responseCall=call+[Run_samples:'8',Run_duration_s:'7.00007009506226'];ui.plot.caption='Cell001 · other words';ui.plot.changed();ui.showFrame(3)
        assert ui.plot.noteLines[0].text=='Cell001 · other words' && ui.plot.noteLines[1].text=='Required rise 0.10 (Floor) · span 7 s · 8 samples' && Double.isFinite(ui.plot.moment) : ui.plot.noteLines*.text
        ui.exportPlots();def again=javax.imageio.ImageIO.read(traceB)
        assert first.width==again.width && first.height==again.height && (0..<first.width).every{int x->(0..<first.height).every{int y->first.getRGB(x,y)==again.getRGB(x,y)}} : 'nothing but the plot is saved'
        ui.cells.selectedItem='Cell005';ui.cells.selectedItem='Cell001';ui.showFrame(17)
        assert ui.plot.responseCall==call && ui.plot.caption=='Cell001 · Did not respond'

        // ---- the mark removed: everything is as it was before
        asked.clear();focusOn=ui.removeMark;ui.removeMark.doClick()
        assert !reviewB.exists() && percentB.bytes==percentBefore : 'the review file goes and the percent table is as it was'
        assert ui.responding.value.text==before.tile && ui.responding.caption.text==before.caption : strip(ui.responding.caption.text)
        assert heights()==sizes.Cell001 && pickers()==pickersBefore : 'the last mark removed moves nothing: '+heights()+' '+sizes.Cell001+' '+pickers()+' '+pickersBefore
        assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==before.groups && ui.cellPicture.calls==before.picture && [ui.heatmap.names,ui.heatmap.calls].transpose()==before.heatmap : 'the Show list, the picture and the heatmap are as before'
        assert ui.cellStatus.text==before.status && ui.cellStatus.toolTipText==before.tip && ui.plot.caption==before.plot && ui.plot.noteLines[0].text==before.plot && ui.cells.selectedItem=='Cell001' : strip(ui.cellStatus.text)
        assert ui.markButton.visible && ui.markButton.text=='Mark as not a responder' && !ui.removeMark.visible : 'the cell can be marked again'
        assert asked.size()==1 && asked[0].is(ui.markButton) : 'the focus goes from Remove mark to the button: '+asked
        assert ui.calls.model.is(ui.callsSaved) && ui.counts.model.is(ui.countsSaved) && columns(ui.calls)==before.callColumns && columns(ui.counts)==before.countColumns : 'without marks the tabs are the saved tables'
        focusOn=null
        ui.exportPlots();record=new groovy.json.JsonSlurper().parse(new File(dirB,'plot_provenance.json'))
        assert !record.containsKey('calls') && !record.containsKey('review_file') : 'without marks the record is as it always was'

        // ---- under All cells, a mark that leaves no cell of a call, or brings one: the plot keeps its place. Cell002 is the
        // only cell that did not respond; marked, the key has no count of cells that did not respond
        ui.cells.selectedItem='Cell002';placed=placeOf()
        ui.markButton.doClick();ui.reasonBox.selectedItem='Noise';ui.saveMark.doClick()
        assert ui.group.selectedIndex==0 && ui.plot.allCalls.count{it=='Did not respond'}==0 && placeOf()==placed : 'a mark that leaves no cell that did not respond: '+placed+' '+placeOf()
        ui.removeMark.doClick();assert !reviewB.exists() && placeOf()==placed : 'the mark removed: '+placed+' '+placeOf()
        // the two responders marked, one after the other: no cell responded
        ['Cell001','Cell005'].each{String c->ui.cells.selectedItem=c;ui.markButton.doClick();ui.reasonBox.selectedItem='Noise';ui.saveMark.doClick()}
        assert ui.plot.allCalls.count{it=='Responded'}==0 && placeOf()==placed : 'marks that leave no responder: '+placed+' '+placeOf()
        // Show on the responders then lists no cell: the call line says so, and nothing under the plot names a cell
        ui.showGroup(1);paint('fiji_traces_no_cells')
        assert ui.group.selectedItem=='Responded after review (0)' && ui.cells.itemCount==0 && ui.cellStatus.text.contains('No cells in this group.') && ui.plot.noteLines*.text==['',''] : ui.plot.noteLines*.text
        ui.showGroup(0)
        ['Cell005','Cell001'].each{String c->ui.cells.selectedItem=c;ui.removeMark.doClick()}
        assert !reviewB.exists() && placeOf()==placed : 'the marks removed'

        // ---- a note as long as the box takes: the call line of the mark takes the room kept for it
        String longNote='bright debris drifts in from the left edge at frame 12, crosses the cell body and leaves the field on the right at 15.5s'
        assert longNote.length()==Review.NOTE_MAX : longNote.length()
        ui.markButton.doClick();ui.reasonBox.selectedItem='Other';ui.noteField.text=longNote;ui.saveMark.doClick()
        assert said()=='' && Review.read(dirB,ui.measurement).marks.Cell001?.note==longNote : said()
        paint('fiji_traces_long_note')
        assert heights()==sizes.Cell001 : 'a call line with the longest note: '+heights()+' '+sizes.Cell001+' '+strip(ui.cellStatus.text)
        ui.removeMark.doClick();assert !reviewB.exists() : 'the long mark is removed'

        // ---- a key held down acts once: its press repeats and does nothing; letting it go opens the form, saves or removes
        focusOn=ui.markButton;def files=snapshot()
        2.times{press(ui.markButton,KeyEvent.VK_ENTER)}
        assert ui.formCell==null && snapshot()==files : 'Enter held down on the button: nothing yet'
        asked.clear();assert release(ui.markButton,KeyEvent.VK_ENTER) && ui.formCell=='Cell001' && ui.reasonBox.selectedItem=='Passing object' && asked==[ui.noteField] : 'let go: the form opens: '+asked
        // a key let go whose press went elsewhere (as to the open Reason list) does nothing
        focusOn=ui.noteField;release(ui.noteField,KeyEvent.VK_ENTER)
        assert ui.formCell=='Cell001' && snapshot()==files : 'Enter let go without its press: nothing is saved'
        3.times{press(ui.noteField,KeyEvent.VK_ENTER)}
        assert ui.formCell=='Cell001' && snapshot()==files : 'Enter held down in the form: nothing is saved yet'
        asked.clear();release(ui.noteField,KeyEvent.VK_ENTER)
        assert ui.formCell==null && Review.read(dirB,ui.measurement).marks.keySet()==['Cell001'] as Set && asked==[ui.removeMark] : 'let go: saved once, and the focus is on Remove mark: '+asked
        def once=snapshot();focusOn=ui.removeMark;3.times{press(ui.removeMark,KeyEvent.VK_ENTER)}
        assert snapshot()==once && ui.removeMark.visible : 'Enter held down on Remove mark: the mark stays'
        release(ui.removeMark,KeyEvent.VK_ENTER);assert !reviewB.exists() && ui.markButton.visible : 'let go: the mark is removed'
        focusOn=null
        // Enter held down in the open Reason list: the list takes its press (it chooses there and closes), so the press never
        // reaches the mark controls; the repeats then do, with the list closed. Neither they nor the key let go save anything.
        // The keys go as Swing hands them on, through the key listeners first (a window never shown has no list to open: the
        // test says when it is open)
        def key={JComponent c,int id->c.processKeyEvent(new KeyEvent(c,id,System.currentTimeMillis(),0,KeyEvent.VK_ENTER,KeyEvent.CHAR_UNDEFINED))}
        def listWas=ui.listOpen
        ui.markButton.doClick();focusOn=ui.reasonBox;def held=snapshot()
        ui.listOpen={->true};key(ui.reasonBox,KeyEvent.KEY_PRESSED);ui.reasonBox.selectedItem='Noise';ui.listOpen={->false}
        3.times{key(ui.reasonBox,KeyEvent.KEY_PRESSED)};key(ui.reasonBox,KeyEvent.KEY_RELEASED)
        assert ui.formCell=='Cell001' && snapshot()==held : 'Enter held in the open Reason list saves nothing: '+Review.read(dirB,ui.measurement).marks.keySet()
        key(ui.reasonBox,KeyEvent.KEY_PRESSED);key(ui.reasonBox,KeyEvent.KEY_RELEASED)
        assert ui.formCell==null && Review.read(dirB,ui.measurement).marks.Cell001?.reason=='Noise' : 'the next Enter saves the reason chosen'
        focusOn=ui.removeMark;ui.removeMark.doClick();assert !reviewB.exists() : 'the mark is removed'
        // the same, with the focus moved on to the Note box before the key is let go: nothing is saved, and the next Enter saves
        ui.markButton.doClick();focusOn=ui.reasonBox
        ui.listOpen={->true};key(ui.reasonBox,KeyEvent.KEY_PRESSED);ui.reasonBox.selectedItem='Movement';ui.listOpen={->false}
        focusOn=ui.noteField;2.times{key(ui.noteField,KeyEvent.KEY_PRESSED)};key(ui.noteField,KeyEvent.KEY_RELEASED)
        assert ui.formCell=='Cell001' && snapshot()==held : 'Enter held from the open Reason list into the Note box saves nothing: '+Review.read(dirB,ui.measurement).marks.keySet()
        key(ui.noteField,KeyEvent.KEY_PRESSED);key(ui.noteField,KeyEvent.KEY_RELEASED)
        assert ui.formCell==null && Review.read(dirB,ui.measurement).marks.Cell001?.reason=='Movement' : 'the next Enter in the Note box saves'
        focusOn=ui.removeMark;ui.removeMark.doClick();assert !reviewB.exists() : 'the mark is removed'
        ui.listOpen=listWas;focusOn=null

        // ---- Enter saves, from the Note box and from the Reason list; a mark without a note
        ui.markButton.doClick();ui.noteField.text='drifts in from the left'
        assert enter(ui.noteField) && ui.formCell==null && Review.read(dirB,ui.measurement).marks.Cell001?.reviewed=='Did not respond' : 'Enter in the Note box saves'
        ui.cells.selectedItem='Cell002';ui.markButton.doClick();ui.reasonBox.selectedItem='Bad outline'
        enter(ui.reasonBox)
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
        assert ui.busy && ui.formCell==null && !ui.markButton.visible && !ui.removeMark.visible : 'no mark controls while an analysis runs'
        ui.cells.selectedItem='Cell001';assert !ui.markButton.visible && !ui.removeMark.visible : 'no Remove mark while an analysis runs'
        def files=snapshot();ui.removeReviewMark('Cell001')
        assert said()=='Wait for the analysis to finish.' && snapshot()==files : 'a mark is refused while an analysis runs'
        ui.cells.selectedItem='Cell002';assert said()=='' : 'the refusal goes when another cell is chosen'
    }
    latch.countDown();finish()
    onUI {
        String carried='Review marks from '+B+': 1 carried over, 0 dropped.'
        assert ui.results==dirA && ui.status.text=='Complete. Protocol, inputs and tables saved in '+A+'. All cells: 2 of 3 responded (66.7 %), after review 1 of 3 (33.3 %). '+carried+' '+unchecked : ui.status.text
        ui.showGroup(0);ui.cells.selectedItem='Cell001'
        Map mark=Review.read(dirA,ui.measurement).marks.Cell001
        assert mark.carriedFrom==B : mark
        assert strip(ui.cellStatus.text)=='● Cell001: Did not respond after review · passing object: drifts in from the left · from an earlier analysis · automatic call Responded · sustained threshold crossing · largest rise 1.00 at 6 s (needed 0.10)' : strip(ui.cellStatus.text)
        assert ui.cellStatus.toolTipText.startsWith('<html>Marked '+minute(mark.time as String)+' in '+B+'<br>') : 'when and where the mark was made, to the minute as on the page: '+ui.cellStatus.toolTipText
        paint('fiji_traces_carried')
        // the finish message counts the marks from the traces the window holds: normalized_traces.csv is not read again
        Closure unlockA=lockAway(new File(dirA,'normalized_traces.csv'))
        try{
            ui.cells.selectedItem='Cell002';ui.markButton.doClick();ui.reasonBox.selectedItem='Bad outline';ui.saveMark.doClick()
            assert Review.read(dirA,ui.measurement).marks.Cell002?.reviewed=='Responded' : said()
            assert rowOf(ui.calls,[Cell:'Cell002',Role:'Stimulus'])[-2..-1]==['Responded','Bad outline'] : 'a mark without a note: '+rowOf(ui.calls,[Cell:'Cell002',Role:'Stimulus'])
            assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+A+'. All cells: 2 of 3 responded (66.7 %), after review 2 of 3 (66.7 %). '+carried+' '+unchecked : ui.status.text
            ui.removeMark.doClick()
        }finally{unlockA()}
        // the carried text stays until the next run, after the last mark is removed too
        ui.cells.selectedItem='Cell001';ui.removeMark.doClick()
        assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+A+'. '+carried+' '+unchecked : 'the finish message follows the marks: '+ui.status.text
        assert ui.responding.value.text=='66.7 %' && !new File(work,'Review_reviewed_calls_'+A+'.csv').exists() : 'the tile and the review files are as before'
    }

    // ---- a recording's cell names, Cell_ and 32 letters (here Cell001 renamed in a copy of the analysis): the Cell box is as
    // wide as for Cell001, so the whole pickers row (Window and Every trace too) stays in the Traces tab at the window's
    // usual size; the name is cut in the box and its tooltip gives it whole
    String longName='Cell_aecb40839a7897469c9c4138ad4734d0',C='SimpleCa-20261101-120000-cccccccc'
    File longFolder=new File(work,'long'),dirC=new File(longFolder,C);dirC.mkdirs()
    new File('tests/fixtures/simpleca').listFiles().findAll{it.isFile()}.each{new File(dirC,it.name).bytes=new String(it.bytes,'UTF-8').replace('Cell001',longName).getBytes('UTF-8')}
    onUI {
        layOut();int shortBox=ui.cells.width   // the Cell box of the window above, for Cell001
        longWindow=Window.newInstance([csv:new File(dirC,'metrics.csv'),folder:longFolder,recording:'Long',frames:30,hz:2d,baselineStart:1,baselineEnd:5,
            picture:original,movie:recording,names:[longName]+names.drop(1),rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]])
        longWindow.runtime.home=new File(work,'runtime')
        def main=ui;ui=longWindow   // the helpers above look at this window for now
        try {
            ui.loadResults(dirC);ui.tabs.selectedIndex=1;ui.showGroup(0);ui.cells.selectedItem=longName
            paint('fiji_long_cell_name')
            def line=ui.cells.parent as JComponent;java.awt.Rectangle seen=line.visibleRect
            assert seen.width==line.width && seen.height==line.height && line.width==ui.tabs.selectedComponent.width : 'the pickers row is seen whole, across the Traces tab: '+seen+' '+line.bounds
            assert [ui.stimuli,ui.allTraces].every{it.visible && seen.contains(it.bounds)} : 'Window and Every trace in the row: '+[ui.stimuli,ui.allTraces]*.bounds+' in '+seen
            assert line.components.every{seen.contains(it.bounds)} : 'every control of the row in it: '+line.components.findAll{!seen.contains(it.bounds)}*.bounds+' in '+seen
            assert ui.cells.width==shortBox : 'the Cell box is as wide as for Cell001: '+ui.cells.width+' '+shortBox
            assert ui.cells.selectedItem==longName && ui.cells.toolTipText==longName : 'the tooltip gives the whole name: '+ui.cells.toolTipText
            ui.cells.selectedItem='Cell002';assert ui.cells.toolTipText=='Cell002' : 'the tooltip follows the cell: '+ui.cells.toolTipText
            ui.showGroup(1);ui.cells.selectedItem=longName;layOut()
            assert ui.cells.toolTipText==longName && [ui.stimuli,ui.allTraces].every{seen.contains(it.bounds)} : 'under Responded too: '+ui.cells.toolTipText+' '+[ui.stimuli,ui.allTraces]*.bounds
        } finally {ui=main}
    }

    // ---- issue 29 in the window, on a synthetic recording whose rings and core the run measures from the open movie: the
    // option that marks flagged responders automatically, the codes of the checks on the tabs, and the Around the cell tab
    def Steps=loader.loadClass('CalciumSteps'),Card=loader.loadClass('CalciumCard'),Around=loader.loadClass('CalciumAroundPlot')
    int W=420,H=320,N=120
    def disc={int cx,int cy,int r->new OvalRoi(cx-r,cy-r,2*r,2*r)}
    // a bright object crosses Cell001 from the left (a possible passing object), another rests on the lower right edge of
    // Cell002 (a rise on the edge), Cell003 responds in place, Cell004 rises on single frames only, Cell005 is excluded
    List aroundNames=(1..5).collect{String.format('Cell%03d',it)}
    Roi[] aroundRois=[[80,80],[200,80],[300,80],[80,220],[300,220]].collect{disc(it[0],it[1],15)} as Roi[]
    List<int[]> masks=aroundRois.collect{Steps.pixelsOf(it,W,H)}
    def spotRandom=new Random(3);List spots=(1..40).collect{[20+spotRandom.nextInt(W-40),20+spotRandom.nextInt(H-40),3+spotRandom.nextInt(3),15+spotRandom.nextInt(20)]}
    def aroundStack=new ImageStack(W,H)
    (1..N).each{int f->
        short[] px=new short[W*H];def noise=new Random(1000+f)
        for(int i=0;i<px.length;i++)px[i]=(short)(100+noise.nextInt(7)-3)
        def put={int i,int amount->px[i]+=(short)amount}
        def blob={double bx,double by,double r,int amount->for(int y=(int)(by-r-1);y<=by+r+1;y++)for(int x=(int)(bx-r-1);x<=bx+r+1;x++)if((x-bx)*(x-bx)+(y-by)*(y-by)<=r*r&&y>=0&&y<H&&x>=0&&x<W)put(y*W+x,amount)}
        spots.each{blob(it[0],it[1],it[2],it[3])}
        masks.each{int[] m->m.each{put(it,50)}}
        if(f>=50&&f<=54)blob(57,80,5,400);if(f>=55&&f<=62)blob(70+(f-55)*20d/7d,80,5,f in [58,59]?440:400);if(f>=63&&f<=67)blob(103,80,5,400)
        if(f>=60&&f<=90)blob(200+22*Math.cos(Math.PI/4),80+22*Math.sin(Math.PI/4),8,1500)
        if(f>=62)masks[2].each{put(it,40)}
        if(f in [60,70,80])masks[3].each{put(it,30)}
        aroundStack.addSlice(''+f,new ij.process.ShortProcessor(W,H,px,null))
    }
    aroundMovie=new ImagePlus('Around recording',aroundStack);def aroundCal=aroundMovie.calibration;aroundCal.pixelWidth=0.645d;aroundCal.pixelHeight=0.645d;aroundCal.setUnit('micron');aroundCal.frameInterval=1d
    // what Measure leaves beside the analysis: the raw gray values of the outlines and the background next to each cell; no
    // ring files, so the run measures them from the open movie
    File aroundFolder=new File(work,'around');aroundFolder.mkdirs();String rec='Around'
    List rawRows=[['Frame']+aroundNames],backRows=[['Frame']+aroundNames]
    (1..N).each{int f->def ip=aroundStack.getProcessor(f);rawRows<<([f]+masks.collect{int[] m->double s=0d;for(int i:m)s+=ip.getf(i);s/m.length});backRows<<([f]+aroundNames.collect{100})}
    Csv.write(new File(aroundFolder,rec+'_gray_values_raw.csv'),rawRows);Csv.write(new File(aroundFolder,rec+'_background_next_to_each_cell.csv'),backRows)
    // the analysis a run of the engine would leave, with QC_flags in each table that has it, and a control window: F0 = frames
    // 1-40, a rise of 10 % and 3 x the baseline SD for 2 samples in a row
    String D='SimpleCa-20261102-100000-dddddddd',E='SimpleCa-20261102-110000-eeeeeeee'
    String engineCalls='reference_below_3_baseline_sd_review; BASELINE_SD_AT_LEAST_ONE_THIRD_LOCAL_F0',engineQuality='reference_below_3_baseline_sd_review'
    File dirD=new File(aroundFolder,D);dirD.mkdirs()
    Map aroundTraces=[:],aroundStatus=[:]
    List callHead=['Field','Cell','Role','Stimulus','Start_s','End_s','Baseline_s','Max_gap_s','Effective_end_s','Status','Reason','Baseline_mean','Baseline_SD','Threshold_increase','Peak_increase','Floor_threshold','Noise_threshold','Threshold_source','Local_reference_relative','Peak_s','Run_start_s','Run_end_s','QC_flags']
    List callBody=[],qualityRows=[['Cell','F0','Baseline_SD_raw','Status','Reason','QC_flags']],metricRows=[['Cell','F0','Status','Reason','Peak_dFF0','QC_flags']],windowRows=[['Field','Cell','Role','Stimulus','Peak_dFF0','QC_flags']]
    aroundNames.eachWithIndex{String cell,int c->
        List y=rawRows.drop(1).collect{(it[c+1] as double)-100d};double f0=y[0..<40].sum()/40d
        List n=y.collect{(it-f0)/f0};double mean=n[0..<40].sum()/40d,sd=Math.sqrt(n[0..<40].collect{(it-mean)*(it-mean)}.sum()/39d),need=Math.max(0.1d,3d*sd)
        List r=n.collect{(it-mean)/(1d+mean)};int peak=(40..<N).max{r[it]};boolean responded=(40..<N-1).any{r[it]>=need&&r[it+1]>=need}
        String status=c==4?'Excluded':responded?'Responded':'Did not respond';aroundStatus[cell]=status
        String reason=c==4?'Invalid raw fluorescence baseline (F0)':responded?'Sustained threshold crossing':'No sustained threshold crossing'
        String engineFlag=c==0?engineCalls:'',quality=c==0?engineQuality:''
        ['Stimulus':'ATP','Control':'KCl'].each{String role,String stimulus->
            callBody<<(c==4 ? [rec,cell,role,stimulus,40,N-1,10,1.5,N-1,status,reason,'','','','',0.1,'','','','','','',engineFlag]
                : [rec,cell,role,stimulus,40,N-1,10,1.5,N-1,status,reason,mean,sd,need,r[peak],0.1,3d*sd,need>0.1d?'Noise':'Floor',1d+mean,peak,responded?peak:'',responded?peak+1:'',engineFlag])
            windowRows<<[rec,cell,role,stimulus,c==4?'':r[peak],engineFlag]}
        qualityRows<<[cell,f0,sd*f0,c==4?'Excluded':'Included',c==4?reason:'',quality];metricRows<<[cell,f0,c==4?'Excluded':'Included',c==4?reason:'',c==4?'':r[peak],quality]
        if(c!=4)aroundTraces[cell]=n
    }
    assert aroundStatus==[Cell001:'Responded',Cell002:'Responded',Cell003:'Responded',Cell004:'Did not respond',Cell005:'Excluded'] : 'the synthetic cells are called as meant: '+aroundStatus
    // the stimulus rows first, then the control rows, as the engine writes them
    Csv.write(new File(dirD,'responder_calls.csv'),[callHead]+callBody.findAll{it[2]=='Stimulus'}+callBody.findAll{it[2]=='Control'})
    Csv.write(new File(dirD,'window_metrics.csv'),[windowRows[0]]+windowRows.drop(1).findAll{it[2]=='Stimulus'}+windowRows.drop(1).findAll{it[2]=='Control'})
    Csv.write(new File(dirD,'cell_quality.csv'),qualityRows);Csv.write(new File(dirD,'metrics.csv'),metricRows)
    List traced=aroundTraces.keySet() as List
    Csv.write(new File(dirD,'normalized_traces.csv'),[['Time']+traced]+(0..<N).collect{int f->[f]+traced.collect{aroundTraces[it][f]}})
    Csv.write(new File(dirD,'mean_trace.csv'),[['Time','Mean','SEM','N_Cells']]+(0..<N).collect{int f->[f,traced.collect{aroundTraces[it][f]}.sum()/traced.size(),0,traced.size()]})
    Csv.write(new File(dirD,'analysis_manifest.csv'),[['key','value'],['sampling_hz','1'],['baseline_start','1'],['baseline_end','40'],['stimulus','ATP'],['protocol_hash','synthetic']])
    def countsHead=Csv.read(new File('tests/fixtures/simpleca/responder_counts.csv'))[0]
    def countRow={String role,String stimulus->Map v=[Field:rec,Stimulus:stimulus,Population:'Cells in field',Registered:5,Responded:3,Nonresponders:1,Evaluable:4,Excluded:1,Unscorable:0,Percent_responding:75,Effective_end_s:N-1,Window_clipped:'FALSE',Role:role,Start_s:40,End_s:N-1];countsHead.collect{v.containsKey(it)?v[it]:''}}
    Csv.write(new File(dirD,'responder_counts.csv'),[countsHead,countRow('Stimulus','ATP'),countRow('Control','KCl')])
    new File(dirD,'COMPLETE').text='complete\n'
    Map engineFiles=['responder_calls.csv','window_metrics.csv','cell_quality.csv','metrics.csv'].collectEntries{[it,Installer.digest(new File(dirD,it).bytes)]}
    File checksD=new File(aroundFolder,rec+'_call_checks_'+D+'.csv'),reviewD=new File(aroundFolder,rec+'_reviewed_calls_'+D+'.csv'),pixelsFile=new File(aroundFolder,rec+'_ring_and_core_pixels.csv')
    // the strings drawn on the around plot, and the plot as it is drawn at a size, on white
    def aroundDrawn={int width,int height->def image=new java.awt.image.BufferedImage(width,height,java.awt.image.BufferedImage.TYPE_INT_RGB)
        def spy=new TextSpy(real:image.createGraphics());java.awt.Dimension was=ui.aroundPlot.size
        try{spy.real.color=java.awt.Color.white;spy.real.fillRect(0,0,width,height);ui.aroundPlot.setSize(width,height);ui.aroundPlot.still(spy)}finally{spy.real.dispose();ui.aroundPlot.setSize(was);ui.aroundPlot.changed()}
        [image:image,drawn:spy.drawn,at:spy.at,top:ui.aroundPlot.@t,tall:ui.aroundPlot.@h]}
    def aroundNow={->layOut();aroundDrawn(ui.aroundPlot.width,ui.aroundPlot.height)}
    // only the axis titles, the tick labels and the key: the cell, the places of the zones, the other ring sectors and the bands
    def allowedAround={String s->s in ['Time (s)','ΔF (gray levels)','other ring sectors','response window','flagged frames',ui.aroundPlot.cell]+((ui.aroundPlot.data?.lines?:[])*.place) || s==~/-?\d+(\.\d+)?/}
    def aroundTab={->ui.tabs.indexOfTab('Around the cell')}
    def flagWords={JComponent shownRow->shownRow.components.findAll{it instanceof JLabel}*.text}
    def aroundRead={->long deadline=System.nanoTime()+30000000000L;boolean done=false
        while(!done&&System.nanoTime()<deadline){onUI{done=ui.results!=null&&ui.aroundFolder==ui.results.name&&ui.aroundContext!=null};if(!done)Thread.sleep(20)}
        assert done : 'what is shown around the cells is read'}
    // the colour of a band at the middle of its frames, down a column of the frame: most of it is the band's
    def bandShare={Map seen,double at,java.awt.Color colour->int x=ui.aroundPlot.X(at);(int)((seen.top..<(seen.top+seen.tall)).count{int y->(seen.image.getRGB(x,y)&0xFFFFFF)==(colour.RGB&0xFFFFFF)})/(double)seen.tall}
    def colourCount={Map seen,java.awt.Color colour->int count=0;for(int x=62;x<seen.image.width-18;x++)for(int y=seen.top;y<seen.top+seen.tall;y++)if((seen.image.getRGB(x,y)&0xFFFFFF)==(colour.RGB&0xFFFFFF))count++;count}
    // the key: each entry on a line of its own or after the one before, none over another and none cut short
    def keyLaid={Map seen->List keyed=seen.at.findAll{(it[2] as int)<(seen.top as int)};def fm;def g=new java.awt.image.BufferedImage(1,1,java.awt.image.BufferedImage.TYPE_INT_RGB).createGraphics();try{Theme.smooth(g);g.font=Theme.font(java.awt.Font.PLAIN,-1f);fm=g.fontMetrics}finally{g.dispose()}
        List rows=keyed.groupBy{it[2]}.values().collect{it.sort{a,b->a[1]<=>b[1]}}
        rows.every{List keyRow->(1..<keyRow.size()).every{int i->(keyRow[i-1][1] as int)+fm.stringWidth(keyRow[i-1][0] as String)<(keyRow[i][1] as int)-21}} && keyed.every{!(it[0] as String).endsWith('…')}}
    def latchE=new java.util.concurrent.CountDownLatch(1)
    File dirE=new File(aroundFolder,E)
    def box
    def main=ui
    try {
        onUI {
            aroundWindow=Window.newInstance([csv:new File(aroundFolder,rec+'_gray_values_raw.csv'),folder:aroundFolder,recording:rec,frames:N,hz:1d,baselineStart:1,baselineEnd:40,
                picture:new ImagePlus('Around colour picture',new ColorProcessor(W,H)),movie:aroundMovie,names:aroundNames,rois:aroundRois.collect{(Roi)it.clone()} as Roi[]])
            aroundWindow.runtime.home=new File(work,'runtime')
            ui=aroundWindow   // the helpers above look at this window for now
            // the tabs: Around the cell right after Traces
            assert (0..<ui.tabs.tabCount).collect{ui.tabs.getTitleAt(it)}.take(4)==['Protocol','Traces','Around the cell','Heatmap'] : (0..<ui.tabs.tabCount).collect{ui.tabs.getTitleAt(it)}
            assert ui.aroundPlot.sentence()=='Run the analysis to see traces' && !ui.aroundPlot.drawn() && !ui.aroundFlags.visible : 'before an analysis the tab says so'
            // the option: a tick box on the Protocol tab like the one beside it, in its own words, off at first
            box=ui.automaticBox
            assert box instanceof JCheckBox && box.text=='Mark possible passing objects as not responders' && !box.selected : box.text
            assert box.toolTipText=='After each analysis, a responder flagged as a possible passing object or as a rise on the edge is marked Did not respond, with the reason Passing object. A mark can be removed.' : box.toolTipText
            assert box.font==ui.kcl.font && box.font==Theme.font() && !box.opaque && box.alignmentX==ui.kcl.alignmentX && box.UI.class==ui.kcl.UI.class : 'the tick box looks like the one for a positive control: '+box.font
            layOut();assert box.preferredSize.height==ui.kcl.preferredSize.height : 'one height: '+box.preferredSize+' '+ui.kcl.preferredSize
            assert SwingUtilities.isDescendingFrom(box,ui.tabs.getComponentAt(ui.tabs.indexOfTab('Protocol'))) && SwingUtilities.getAncestorOfClass(Card,box)?.title=='Review marks' : 'on the Protocol tab'
            // the settings with the option off are what they were before it existed; on, they say so
            assert !ui.fingerprint().contains('automatic_marks') : ui.fingerprint()
            box.doClick();assert box.selected && ui.fingerprint().contains('"automatic_marks":true') : ui.fingerprint()
            ui.tabs.selectedIndex=0;def port=SwingUtilities.getAncestorOfClass(JViewport,box) as JViewport;layOut()
            port.viewPosition=new java.awt.Point(0,Math.max(0,SwingUtilities.convertPoint(box,0,0,port.view).y as int)-port.height+80);paint('fiji_protocol_tick_box',false)
            // a run with the option on: the checks measure the rings from the open movie, then the flagged responders are marked
            ui.confirmed.selected=true
            ui.performAnalysis={a,b,c,d,job->dirD}
            ui.runButton.doClick();assert ui.busy && !box.enabled : 'the option cannot change while the analysis runs'
        }
        finish()
        onUI {
            assert box.enabled : 'and can once it is done'
            assert ui.results==dirD && pixelsFile.isFile() && checksD.isFile() : 'the run measured around the cells from the movie: '+ui.status.text
            assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+D+'. All cells: 3 of 4 responded (75.0 %), after review 1 of 4 (25.0 %). Marked automatically: Cell001, Cell002.' : 'the finish message names the cells marked automatically, after the figures: '+ui.status.text
            assert ui.checks.cells.Cell001.codes==['POSSIBLE_PASSING_OBJECT'] && ui.checks.cells.Cell002.codes==['RISE_ON_EDGE'] && ui.checks.cells.Cell003.codes==[] : 'the flags of the synthetic cells'
            def saved=Csv.read(reviewD)
            assert saved[0]==Review.HEADER && saved.drop(1).collect{[it[0],it[1],it[2],it[3],it[4],it[-1]]}==[['Cell001','Responded','Did not respond','Passing object','','automatic'],['Cell002','Responded','Did not respond','Passing object','','automatic']] : saved
            assert ui.appliedFingerprint.contains('"automatic_marks":true') && ui.appliedNote.text=='Displayed results match the applied settings.' : ui.appliedNote.text
            // a window opened again on this analysis: its review file holds automatic marks, so the box is ticked and the results
            // match the settings
            def again=Window.newInstance(ui.measurement+[rois:(ui.measurement.rois as Roi[]).collect{(Roi)it.clone()} as Roi[]]);again.runtime.home=new File(work,'runtime');again.aroundLater={Closure read->}
            try{
                assert !again.automaticBox.selected : 'a new window starts with the box off'
                again.loadResults(dirD)
                assert again.automaticBox.selected && again.appliedFingerprint.contains('"automatic_marks":true') && again.appliedNote.text=='Displayed results match the applied settings.' : again.appliedNote.text
                again.automaticBox.doClick();assert again.appliedNote.text.startsWith('Settings changed') : 'turned off, it asks for a new run'
            }finally{again.dispose()}
            // the option is a setting: changing it asks for a new run, changing it back does not
            box.doClick();assert ui.appliedNote.text.startsWith('Settings changed') : ui.appliedNote.text
            box.doClick();assert ui.appliedNote.text=='Displayed results match the applied settings.' : ui.appliedNote.text
            // the mark in the call line, the Review reason column and the Marked group, as for any mark
            ui.showGroup(0);ui.cells.selectedItem='Cell001'
            assert strip(ui.cellStatus.text).startsWith('● Cell001: Did not respond after review · passing object · marked automatically · automatic call Responded · sustained threshold crossing · largest rise ') : strip(ui.cellStatus.text)
            assert ui.removeMark.visible && !ui.markButton.visible : 'an automatic mark can be removed'
            paint('fiji_marked_automatically')
            assert rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])[-2..-1]==['Did not respond','Passing object, automatic'] && rowOf(ui.calls,[Cell:'Cell002',Role:'Stimulus'])[-2..-1]==['Did not respond','Passing object, automatic'] : rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])
            assert rowOf(ui.calls,[Cell:'Cell001',Role:'Control'])[-2..-1]==['',''] : 'not on the control window'
            assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}[-1]=='Marked (2)' : (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}
            ui.showGroup(4);assert (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell001','Cell002'] : 'the Marked group lists the automatic marks'
            ui.showGroup(0)
            // the codes of the checks after the engine's codes in QC_flags, joined as each table joins them; none on a row of the
            // control window or of a cell without a flag; the saved tables as they were
            int callFlags=columns(ui.calls).indexOf('QC_flags'),windowFlags=columns(ui.windowMetrics).indexOf('QC_flags'),qualityFlags=columns(ui.qc).indexOf('QC_flags'),metricFlags=columns(ui.metrics).indexOf('QC_flags')
            assert rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])[callFlags]==engineCalls+'; POSSIBLE_PASSING_OBJECT' && rowOf(ui.calls,[Cell:'Cell002',Role:'Stimulus'])[callFlags]=='RISE_ON_EDGE' : rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])[callFlags]
            assert rowOf(ui.calls,[Cell:'Cell001',Role:'Control'])[callFlags]==engineCalls && rowOf(ui.calls,[Cell:'Cell002',Role:'Control'])[callFlags]=='' && rowOf(ui.calls,[Cell:'Cell003',Role:'Stimulus'])[callFlags]=='' : 'no codes on the control window'
            assert rowOf(ui.windowMetrics,[Cell:'Cell001',Role:'Stimulus'])[windowFlags]==engineCalls+'; POSSIBLE_PASSING_OBJECT' && rowOf(ui.windowMetrics,[Cell:'Cell002',Role:'Stimulus'])[windowFlags]=='RISE_ON_EDGE' && rowOf(ui.windowMetrics,[Cell:'Cell001',Role:'Control'])[windowFlags]==engineCalls : 'the stimulus-window metrics'
            assert rowOf(ui.qc,[Cell:'Cell001'])[qualityFlags]==engineQuality+';POSSIBLE_PASSING_OBJECT' && rowOf(ui.qc,[Cell:'Cell002'])[qualityFlags]=='RISE_ON_EDGE' && rowOf(ui.qc,[Cell:'Cell003'])[qualityFlags]=='' : 'Cell quality: '+rowOf(ui.qc,[Cell:'Cell001'])[qualityFlags]
            assert rowOf(ui.metrics,[Cell:'Cell001'])[metricFlags]==engineQuality+';POSSIBLE_PASSING_OBJECT' && rowOf(ui.metrics,[Cell:'Cell002'])[metricFlags]=='RISE_ON_EDGE' : 'the whole-recording metrics'
            assert engineFiles.every{name,hash->Installer.digest(new File(dirD,name).bytes)==hash} : 'no engine table is written'
            assert ui.calls.model.is(ui.callsSaved)==false && ui.callsSaved.getValueAt((0..<ui.callsSaved.rowCount).find{ui.callsSaved.getValueAt(it,1)=='Cell001'} as int,callFlags)==engineCalls+'; POSSIBLE_PASSING_OBJECT' : 'the reviewed columns are added to the table with its codes'
            [['Cell response calls','fiji_codes_calls_tab',ui.calls],['Cell quality','fiji_codes_quality_tab',ui.qc]].each{String tab,String name,JTable shown->
                def scroll=SwingUtilities.getAncestorOfClass(JScrollPane,shown) as JScrollPane;if(scroll.columnHeader?.view==null)scroll.columnHeaderView=shown.tableHeader
                ui.tabs.selectedIndex=ui.tabs.indexOfTab(tab);layOut();def view=shown.parent as JViewport
                int x=shown.getCellRect(0,columns(shown).indexOf('QC_flags'),true).x as int;view.viewPosition=new java.awt.Point(Math.max(0,Math.min(x-200,shown.width-view.width)),0);paint(name,false)
                // a cell whose text its column cuts shows the whole text as its tooltip; one whose text fits has none; no width changes
                int c=columns(shown).indexOf('QC_flags'),cellAt=columns(shown).indexOf('Cell'),r=(0..<shown.rowCount).find{shown.getValueAt(it,cellAt)=='Cell001'} as int
                String whole=shown.getValueAt(r,c) as String;java.awt.Rectangle place=shown.getCellRect(r,c,false)
                def tipAt={int column->java.awt.Rectangle at=shown.getCellRect(r,column,false);shown.getToolTipText(new java.awt.event.MouseEvent(shown,java.awt.event.MouseEvent.MOUSE_MOVED,0L,0,(at.x as int)+3,(at.y as int)+3,0,false))}
                assert shown.getFontMetrics(shown.font).stringWidth(whole)>place.width && tipAt(c)==whole : 'cut: '+[whole,place.width,tipAt(c)]
                assert tipAt(cellAt)==null && (0..<shown.columnCount).every{shown.columnModel.getColumn(it).preferredWidth==160} : tipAt(cellAt)}
        }
        aroundRead()
        onUI {
            // ---- Around the cell, rule A: the cell, the sectors seen before and after it, the other ring sectors thin and gray, the
            // response window as on Traces and the flagged frames; the call line and the flags of Traces above it
            ui.tabs.selectedIndex=aroundTab();ui.cells.selectedItem='Cell001';layOut()
            Map data=ui.aroundPlot.data,direct=Checks.aroundCell(Checks.aroundContext(dirD,ui.measurement),ui.checks,'Cell001')
            assert data.why=='' && data.rule=='passes through' && data.lines*.part==['before','after'] && data.band==[55,62] : data.subMap(['why','rule','lines','band'])
            assert java.util.Arrays.equals(data.rise as double[],direct.rise as double[]) && (0..<data.zones.size()).every{int z->java.util.Arrays.equals(data.zones[z].trace as double[],direct.zones[z].trace as double[])} : 'the traces the check used'
            List shown=ui.aroundPlot.lines()
            assert shown.findAll{it[2]!='other ring sectors'}.collect{it[2]}==data.lines*.place && (0..<data.lines.size()).every{int i->shown.findAll{it[2]!='other ring sectors'}[i][3].is(data.zones.find{it.zone==data.lines[i].zone}.trace)} : 'each zone of the flag is drawn from its trace'
            assert shown.findAll{it[2]=='other ring sectors'}.size()==data.others.size() && shown.findAll{it[2]=='other ring sectors'}.every{it[0]==Around.OTHER && it[1]<1f} && Around.OTHER==new java.awt.Color(0xAEB8C2) : 'the other ring sectors thin and gray'
            // each part in its colour, as on the Step 4 page: before blue, after purple; the band of the flagged frames lilac
            assert shown.findAll{it[2]!='other ring sectors'}.collect{it[0]}==[new java.awt.Color(0x0072B2),new java.awt.Color(0xAA3377)] && Around.FLAGGED==new java.awt.Color(0xE2D6F3) : shown.collect{it[0]}
            List colours=shown.findAll{it[2]!='other ring sectors'}.collect{it[0]};assert colours.unique(false).size()==colours.size() && !colours.contains(Theme.INK) && !colours.contains(Around.OTHER) : 'each zone of the flag in its own colour, not the cell\'s: '+colours
            assert ui.aroundPlot.entries(shown)*.getAt(2)==['Cell001']+data.lines*.place+['other ring sectors','response window','flagged frames'] : ui.aroundPlot.entries(shown)*.getAt(2)
            assert ui.aroundStatus.text==ui.cellStatus.text && ui.aroundStatus.toolTipText==ui.cellStatus.toolTipText && flagWords(ui.aroundFlags)==flagWords(ui.flagRow) && flagWords(ui.aroundFlags)==['possible passing object, frames 55 to 62'] : 'the call line and the flags of Traces: '+flagWords(ui.aroundFlags)
            def around=ui.aroundFlags.components.find{it instanceof JLabel}
            assert around.foreground==Theme.BLUE && around.cursor.type==java.awt.Cursor.HAND_CURSOR && around.toolTipText.startsWith('<html>Go to frame 55<br>') : 'the flag is a link here too'
            // a click here stays on this tab and moves its time cursor to the frame
            ui.showFrame(25);click(around);assert ui.tabs.selectedIndex==aroundTab() && ui.frame==55 && ui.aroundPlot.moment==54d && ui.plot.moment==54d : [ui.tabs.selectedIndex,ui.frame,ui.aroundPlot.moment]
            Map seen=aroundNow()
            assert seen.drawn.every{allowedAround(it)} && seen.drawn.containsAll(['Time (s)','Cell001','ring, left','ring, right','other ring sectors','response window','flagged frames']) : 'drawn on the plot: '+seen.drawn.findAll{!allowedAround(it)}+' of '+seen.drawn
            assert keyLaid(seen) && seen.tall>=250 : 'the key is readable at the usual size and leaves the plot its room: '+seen.tall
            assert bandShare(seen,58.5d,Around.FLAGGED)>0.6d && bandShare(seen,100d,Around.WINDOW)>0.6d && bandShare(seen,20d,Around.FLAGGED)==0d && bandShare(seen,20d,Around.WINDOW)==0d : 'the flagged frames and the response window: '+[bandShare(seen,58.5d,Around.FLAGGED),bandShare(seen,100d,Around.WINDOW)]
            // the band runs half a frame to each side of the frames named
            assert ui.aroundPlot.bands().find{it[3]=='flagged frames'}[0..1]==[54d-0.5d,61d+0.5d] : ui.aroundPlot.bands()
            // half the median frame step, not the step beside the band: of an even number of steps the upper middle one, as on the page
            def uneven=Around.newInstance();uneven.data=[why:'',time:[0d,1d,2d,4d,6d,8d,10d] as double[],band:[3,5],zones:[],lines:[],others:[]]
            assert uneven.bands()==[[1d,7d,Around.FLAGGED,'flagged frames']] : uneven.bands()
            uneven.data=[why:'',time:[0d,1d,2d,4d,6d] as double[],band:[2,3],zones:[],lines:[],others:[]];assert uneven.bands()*.take(2)==[[0d,3d]] : uneven.bands()
            assert data.lines.every{Map line->colourCount(seen,ui.aroundPlot.lines().find{it[2]==line.place}[0] as java.awt.Color)>20} && colourCount(seen,Theme.INK)>20 : 'each line is drawn in its colour'
            paint('fiji_around_rule_a',false)
            // the time cursor follows the slider and the Traces plot; a click on the plot moves it
            ui.showFrame(70);assert ui.aroundPlot.moment==69d && ui.plot.moment==69d : ui.aroundPlot.moment
            ui.aroundPlot.seek.call(30.2d);assert ui.frame==31 && ui.aroundPlot.moment==30d : ui.frame
            def tip=new java.awt.event.MouseEvent(ui.aroundPlot,java.awt.event.MouseEvent.MOUSE_MOVED,0L,0,ui.aroundPlot.X(58.5d),(seen.top as int)+20,0,false)
            assert ui.aroundPlot.getToolTipText(tip)=='<html>response window · flagged frames<br>'+Around.HINT+'</html>' : ui.aroundPlot.getToolTipText(tip)
            ui.showFrame(58);paint('fiji_around_rule_a_cursor',false)
            Map placeA=[top:seen.top,y:ui.aroundPlot.y]
            // rule B: the rim, the ring outside it and the core
            ui.cells.selectedItem='Cell002';data=ui.aroundPlot.data;seen=aroundNow()
            assert data.rule=='on the edge' && data.lines*.part==['rim','ring','core'] && data.band==[60,90] : data.subMap(['rule','lines','band'])
            assert seen.drawn.every{allowedAround(it)} && ['rim, lower right','ring, lower right','core','flagged frames'].every{seen.drawn.contains(it)} && keyLaid(seen) : seen.drawn
            colours=ui.aroundPlot.lines().findAll{it[2]!='other ring sectors'}.collect{it[0]};assert colours==[new java.awt.Color(0x0072B2),new java.awt.Color(0xAA3377),new java.awt.Color(0x009E73)] : 'rim blue, ring purple, core green: '+colours
            assert [top:seen.top,y:ui.aroundPlot.y]==placeA : 'the plot keeps its place from cell to cell: '+seen.top+' '+placeA
            paint('fiji_around_rule_b',false)
            // a responder without a flag: the cell, the core and the eight near-ring sectors, no band
            ui.cells.selectedItem='Cell003';data=ui.aroundPlot.data;seen=aroundNow()
            assert data.band==null && data.lines*.part==['core'] && data.others==(0..7) && !seen.drawn.contains('flagged frames') && colourCount(seen,Around.FLAGGED)==0 && seen.drawn.contains('response window') : data.subMap(['band','lines','others'])
            assert flagWords(ui.aroundFlags)==[] && !ui.aroundFlags.visible && [top:seen.top,y:ui.aroundPlot.y]==placeA : 'no flag, and the plot in its place'
            paint('fiji_around_no_flag',false)
            ui.cells.selectedItem='Cell004';seen=aroundNow();assert [top:seen.top,y:ui.aroundPlot.y]==placeA && seen.drawn.every{allowedAround(it)} : 'a cell that did not respond'
            // and as narrow as at the smallest window, and at the width saved: one place for every cell
            List widths=[400,650,ui.aroundPlot.width]
            Map tops=['Cell001','Cell002','Cell003','Cell004'].collectEntries{ui.cells.selectedItem=it;layOut();[it,widths.collect{aroundDrawn(it as int,ui.aroundPlot.height).top}]}
            assert (tops.values() as Set).size()==1 : 'the plot keeps its place at every width: '+tops
            // without data: one plain sentence where the plot would be, and nothing else
            ui.cells.selectedItem='Cell005';seen=aroundNow()
            assert !ui.aroundPlot.drawn() && ui.aroundPlot.sentence()=='Passing object not checked: no call.' && seen.drawn==['Passing object not checked: no call.'] : seen.drawn
            assert ui.aroundStatus.text==ui.cellStatus.text && ui.aroundPlot.y==placeA.y : 'the call line of Traces here too'
            paint('fiji_around_no_data',false)
            // the control window: the checks are made for the primary stimulus
            ui.stimuli.selectedItem='KCl';ui.cells.selectedItem='Cell001'
            assert ui.aroundPlot.sentence()=='Checks are made for the primary stimulus.' && !ui.aroundFlags.visible && ui.aroundStatus.text==ui.cellStatus.text : ui.aroundPlot.sentence()
            paint('fiji_around_control',false)
            ui.stimuli.selectedItem='ATP';ui.cells.selectedItem='Cell001'
            // Save plots: around_the_cell.png with the same key and no time cursor, drawn across its whole width on white, and the
            // plot record names it
            ui.showFrame(58);layOut();ui.exportPlots()
            File picture=new File(dirD,'around_the_cell.png');def first=javax.imageio.ImageIO.read(picture)
            def record=new groovy.json.JsonSlurper().parse(new File(dirD,'plot_provenance.json'))
            assert record.around_the_cell=='around_the_cell.png' && record.cell=='Cell001' : record
            assert first.width==Math.max(650,ui.aroundPlot.width) && first.height==Math.max(360,ui.aroundPlot.height) && (0..<first.height).every{int y->first.getRGB(first.width-1,y)==(0xFFFFFFFF as int)} : 'around_the_cell.png is drawn to its right edge: '+first.width+' x '+first.height
            def plain=aroundDrawn(first.width,first.height)
            assert (0..<first.width).every{int x->(0..<first.height).every{int y->first.getRGB(x,y)==plain.image.getRGB(x,y)}} : 'around_the_cell.png is the plot with its key and nothing else'
            assert ui.aroundPlot.moment==57d && ui.frame==58 : 'the time cursor is back on screen'
            java.nio.file.Files.copy(picture.toPath(),new File(pictures,'fiji_around_export.png').toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            ui.cells.selectedItem='Cell002';ui.exportPlots();java.nio.file.Files.copy(picture.toPath(),new File(pictures,'fiji_around_export_rule_b.png').toPath(),java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            // a cell without data saves no picture of a sentence, and takes away the one of the cell before
            ui.cells.selectedItem='Cell005';ui.exportPlots();record=new groovy.json.JsonSlurper().parse(new File(dirD,'plot_provenance.json'))
            assert !picture.exists() && !record.containsKey('around_the_cell') && new File(dirD,'trace.png').isFile() : record
            // ---- rule C as the back end gives it: the near-ring sectors on the start and end sides, without frames, and the core;
            // the window draws any mix of parts and keys nothing on the rule
            Map real=ui.aroundOf.call(ui.aroundContext,'Cell001');Closure given=ui.aroundOf
            Map crossing=real+[rule:'crosses the cell',text:'possible passing object, frames 54 to 63',lines:[[zone:4,place:'ring, left',part:'before',frame:null],[zone:0,place:'ring, right',part:'after',frame:null],[zone:16,place:'core',part:'core',frame:null]],others:[1,2,3,5,6,7],band:[54,63]]
            ui.aroundOf={Map context,String cell->cell=='Cell001'?crossing:given.call(context,cell)}
            ui.cells.selectedItem='Cell001';seen=aroundNow();shown=ui.aroundPlot.lines()
            assert ui.aroundPlot.entries(shown)*.getAt(2)==['Cell001','ring, left','ring, right','core','other ring sectors','response window','flagged frames'] && seen.drawn.every{allowedAround(it)} : seen.drawn
            colours=shown.findAll{it[2]!='other ring sectors'}.collect{it[0]};assert colours==[new java.awt.Color(0x0072B2),new java.awt.Color(0xAA3377),new java.awt.Color(0x009E73)] && colours.every{colourCount(seen,it as java.awt.Color)>20} : 'before blue, after purple, core green: '+colours
            assert bandShare(seen,58.5d,Around.FLAGGED)>0.6d && ui.aroundPlot.bands().find{it[3]=='flagged frames'}[0..1]==[53d-0.5d,62d+0.5d] && keyLaid(seen) && seen.top==placeA.top : 'its band and key, the plot in its place'
            paint('fiji_around_rule_c',false)
            // a sector that is not usable is left out; a part met twice takes another colour
            crossing.lines=[[zone:0,place:'ring, right',part:'after',frame:null],[zone:16,place:'core',part:'core',frame:null]];ui.showAround('Cell001');shown=ui.aroundPlot.lines()
            assert ui.aroundPlot.entries(shown)*.getAt(2)==['Cell001','ring, right','core','other ring sectors','response window','flagged frames'] : ui.aroundPlot.entries(shown)*.getAt(2)
            crossing.lines=[[zone:4,place:'ring, left',part:'after',frame:null],[zone:0,place:'ring, right',part:'after',frame:null]];ui.showAround('Cell001')
            colours=ui.aroundPlot.lines().findAll{it[2]!='other ring sectors'}.collect{it[0]};assert colours==[new java.awt.Color(0xAA3377),new java.awt.Color(0xE69F00)] : colours
            crossing.lines=[[zone:4,place:'ring, left',part:'before',frame:null],[zone:5,place:'ring, upper left',part:'before',frame:null],[zone:0,place:'ring, right',part:'after',frame:null],[zone:1,place:'ring, lower right',part:'after',frame:null],[zone:2,place:'ring, bottom',part:'after',frame:null]];ui.showAround('Cell001')
            colours=ui.aroundPlot.lines().findAll{it[2]!='other ring sectors'}.collect{it[0]}
            assert colours==[0x0072B2,0xE69F00,0xAA3377,0x56B4E9,0xCC79A7].collect{new java.awt.Color(it)} : 'each part met again takes the next spare colour: '+colours
            ui.aroundOf=given
            // a flag of rule C on a responder: its code chooses Passing object in the mark form, as for rule A
            Map was=ui.checks.cells.Cell003
            ui.checks.cells.Cell003=was+[text:'possible passing object, frames 62 to 70',rule:'crosses the cell',flags:[[text:'possible passing object, frames 62 to 70',rule:'crosses the cell',code:'POSSIBLE_PASSING_OBJECT',frames:[62,70]]],codes:['POSSIBLE_PASSING_OBJECT'],qcFlags:'POSSIBLE_PASSING_OBJECT']
            ui.tabs.selectedIndex=1;ui.cells.selectedItem='Cell002';ui.cells.selectedItem='Cell003'
            assert flagWords(ui.flagRow)==['possible passing object, frames 62 to 70'] && flagWords(ui.aroundFlags)==flagWords(ui.flagRow) : flagWords(ui.flagRow)
            ui.markButton.doClick();assert ui.reasonBox.selectedItem=='Passing object' : 'rule C: '+ui.reasonBox.selectedItem
            ui.cancelMark.doClick();ui.checks.cells.Cell003=was
            // ---- an automatic mark removed, and the cell marked again by hand: a manual mark
            ui.cells.selectedItem='Cell002';ui.removeMark.doClick()
            assert Review.read(dirD,ui.measurement).marks.keySet()==['Cell001'] as Set && strip(ui.cellStatus.text).startsWith('● Cell002: Responded · ') : strip(ui.cellStatus.text)
            ui.markButton.doClick();assert ui.reasonBox.selectedItem=='Passing object' : 'a rise on the edge: Passing object is chosen';ui.saveMark.doClick()
            Map mark=Review.read(dirD,ui.measurement).marks.Cell002
            assert mark.by=='manual' && mark.reason=='Passing object' : mark
            assert strip(ui.cellStatus.text).startsWith('● Cell002: Did not respond after review · passing object · automatic call Responded · ') && strip(ui.aroundStatus.text)==strip(ui.cellStatus.text) : strip(ui.cellStatus.text)
            assert rowOf(ui.calls,[Cell:'Cell002',Role:'Stimulus'])[-2..-1]==['Did not respond','Passing object'] && rowOf(ui.calls,[Cell:'Cell001',Role:'Stimulus'])[-2..-1]==['Did not respond','Passing object, automatic'] : rowOf(ui.calls,[Cell:'Cell002',Role:'Stimulus'])
            assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+D+'. All cells: 3 of 4 responded (75.0 %), after review 1 of 4 (25.0 %). Marked automatically: Cell001, Cell002.' : 'the finish message stays until the next run: '+ui.status.text
            // while what is shown around the cells is read, the tab says so; then it shows the cell selected
            List reads=[];ui.aroundLater={Closure read->reads<<read}
            ui.loadResults(dirD);ui.tabs.selectedIndex=aroundTab();ui.cells.selectedItem='Cell001'
            assert ui.aroundPlot.sentence()=='Reading the values around the cells.' && reads.size()==1 : ui.aroundPlot.sentence()
            reads[0].call();assert ui.aroundPlot.drawn() && ui.aroundPlot.data.cell=='Cell001' : 'read: '+ui.aroundPlot.sentence()
            // a read that ends after another analysis was loaded is not used
            reads.clear();ui.loadResults(dirD);ui.loadResults(dirD);assert reads.size()==2
            reads[0].call();assert ui.aroundPlot.sentence()=='Reading the values around the cells.' : 'the first read is not used'
            reads[1].call();assert ui.aroundPlot.drawn() : 'the last one is'
            ui.aroundLater={Closure read->read.call()}
            // ---- the next analysis, with the option off: no automatic mark; the manual mark is carried over, the automatic one is not
            box.doClick();assert !box.selected
            dirE.mkdirs();dirD.listFiles().findAll{it.isFile()&&it.name.endsWith('.csv')||it.name=='COMPLETE'}.each{new File(dirE,it.name).bytes=it.bytes}
            ui.performAnalysis={a,b,c,d,job->latchE.await(20,java.util.concurrent.TimeUnit.SECONDS);dirE}
            ui.runButton.doClick();assert ui.busy && !box.enabled : 'the tick box waits for the run'
        }
        latchE.countDown();finish()
        onUI {
            assert ui.results.name==E && box.enabled
            assert ui.status.text=='Complete. Protocol, inputs and tables saved in '+E+'. All cells: 3 of 4 responded (75.0 %), after review 2 of 4 (50.0 %). Review marks from '+D+': 1 carried over, 0 dropped.' : ui.status.text
            Map marks=Review.read(dirE,ui.measurement).marks
            assert marks.keySet()==['Cell002'] as Set && marks.Cell002.by=='manual' && marks.Cell002.carriedFrom==D : marks
            assert !ui.appliedFingerprint.contains('automatic_marks') && ui.appliedNote.text=='Displayed results match the applied settings.' : 'the run without the option'
            // opened again: its review file holds a manual mark only, so the box stays off
            def againE=Window.newInstance(ui.measurement+[rois:(ui.measurement.rois as Roi[]).collect{(Roi)it.clone()} as Roi[]]);againE.runtime.home=new File(work,'runtime');againE.aroundLater={Closure read->}
            try{againE.loadResults(dirE);assert !againE.automaticBox.selected && !againE.appliedFingerprint.contains('automatic_marks') && againE.appliedNote.text=='Displayed results match the applied settings.' : againE.appliedNote.text}finally{againE.dispose()}
            box.doClick();assert ui.appliedNote.text.startsWith('Settings changed') : 'the option turned on asks for a new run'
            box.doClick();assert ui.appliedNote.text=='Displayed results match the applied settings.'
            ui.tabs.selectedIndex=1
        }
    } finally {ui=main}
    assert appeared.isEmpty() : 'a window appeared on screen: '+appeared
    onUI{assert !ui.visible && !ui.isDisplayable() && !longWindow.visible && !longWindow.isDisplayable() && !aroundWindow.visible && !aroundWindow.isDisplayable() : 'the windows were never shown'}
    println 'PASS: Step 4 window flags, links to frames, mark controls, reasons chosen from flags, refusals, Enter (once when held, never from the open Reason list) and Esc, focus back to the button, one head height with its rows one under another, tiles, plot, Show and the Cell box kept at a mark, nothing under the plot when Show lists no cell, nothing written on the plot or trace.png but its key, the key on screen as before the review controls with the frame where it started, bands and dashed lines named in the key of trace.png, laid across it, the plot in its place from cell to cell and at a mark that empties the group, the rule line with the decimals of the call line and none without a trace, saved figures drawn across their width, reviewed columns on the tabs, a review file for another analysis, marks saved and removed with every view following and Show kept, control window and busy without controls, carried marks, plot record, a long cell name cut in the Cell box with the whole row in the tab, the option that marks flagged responders automatically, the codes of the checks in QC_flags, the Around the cell tab and around_the_cell.png'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();longWindow?.dispose();aroundWindow?.dispose();original?.close();recording?.close();aroundMovie?.close()}}
System.exit(0)

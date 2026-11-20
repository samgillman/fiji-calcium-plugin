import javax.swing.*
import ij.*
import ij.gui.*
import ij.process.ColorProcessor

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Window=loader.loadClass('CalciumAnalysisWindow')
def ui, original, imagej, recording
def onUI={Closure c->SwingUtilities.invokeAndWait(c as Runnable)}
def finish={->long deadline=System.nanoTime()+8000000000L;while(ui.busy&&System.nanoTime()<deadline)Thread.sleep(20);onUI{};assert !ui.busy}
try {
    onUI {
        imagej=new ImageJ(ImageJ.NO_SHOW)
        original=new ImagePlus('Original colour picture',new ColorProcessor(250,150))
        def frames=new ImageStack(250,150)
        (1..30).each{int f->def ip=new ij.process.ShortProcessor(250,150);ip.set(100);if(f>=10&&f<=20){ip.setRoi(30,50,20,20);ip.set(100+40);ip.resetRoi()};frames.addSlice(''+f,ip)}
        recording=new ImagePlus('Recording',frames)
        ui=Window.newInstance([csv:new File('tests/fixtures/simpleca/metrics.csv'),folder:new File('build/analysis-gui'),recording:'Fixture',frames:30,hz:2d,baselineStart:1,baselineEnd:5,
            picture:original,movie:recording,names:(1..6).collect{String.format('Cell%03d',it)},rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]])
        ui.visible=true
        assert !ui.measurement.picture.is(original) : 'Analysis must not edit Step 2 cells'
        assert ui.measurement.picture.getWindow()==null && ui.cellPicture.image!=null : 'The cell picture belongs inside the analysis window'
        assert ui.responding.value.text=='–' && !ui.previousCell.enabled
        assert ui.kcl.text=='Include a KCl control' && ui.kcl.toolTipText.contains('positive control')
        // The movie is offered beside the colour picture; before an analysis the picture is still the colour picture.
        assert (0..<ui.view.itemCount).collect{ui.view.getItemAt(it)}==['Colour picture','Movie: change from rest','Movie: as recorded']
        assert ui.view.selectedIndex==0 && ui.cellPicture.image.is(ui.snapshot) && ui.slider.maximum==30
        // Suggested settings are shown rounded, but an untouched box still sends the exact value.
        def precise=Window.newInstance([csv:new File('tests/fixtures/simpleca/metrics.csv'),folder:new File('build/analysis-gui'),recording:'Precise',frames:300,hz:1.0000537066691941d,
            baselineStart:55,baselineEnd:74,names:[],rois:new Roi[0]])
        try {
            assert precise.fields.sampling_hz.text=='1' && precise.fields.stimulus_start.text=='72.996' && precise.fields.stimulus_end.text=='298.984'
            assert precise.fields.sampling_hz.toolTipText.contains('1.0000537066691941')
            // The Fiji window and the Step 4 web page offer the same settings, shown the same way.
            def Settings=loader.loadClass('CalciumAnalysisSettings')
            def offered=Settings.suggested([frames:300,hz:1.0000537066691941d,baselineStart:55,baselineEnd:74])
            assert offered.keySet()==precise.fields.keySet()
            offered.each{key,value->assert precise.fields[key].text==Settings.shown(value) : key}
            precise.confirmed.selected=true
            def sent=precise.settings()
            assert sent.sampling_hz=='1.0000537066691941' && sent.stimulus_start==''+(73/1.0000537066691941d) && sent.stimulus_end==''+(299/1.0000537066691941d)
            assert sent.amplitude_floor=='0.10' && sent.noise_multiplier=='3' && sent.baseline_start=='55' && sent.baseline_end=='74'
            precise.fields.stimulus_start.text='80';assert precise.settings().stimulus_start=='80' && precise.settings().stimulus_end==''+(299/1.0000537066691941d)
            def kclBoxes=precise.fields.kcl_start.parent.parent
            assert !kclBoxes.visible;precise.kcl.doClick();assert kclBoxes.visible;precise.kcl.doClick();assert !kclBoxes.visible
            assert precise.runButton.font.family==UIManager.getFont('Label.font').family : 'Use the interface font, not a fixed family'
        } finally {precise.dispose()}
        ui.runButton.doClick();assert !ui.busy && ui.status.text.contains('Confirm')
        ui.confirmed.selected=true
        ui.fields.stimulus_start.text='3';ui.fields.stimulus_end.text='14.5';ui.fields.response_baseline_seconds.text='2.5'
        ui.performAnalysis={a,b,c,d,job->throw new IOException('Simulated engine error')}
        ui.runButton.doClick();assert ui.busy && ui.cancelButton.enabled
    }
    finish()
    onUI {
        assert ui.status.text.contains('Simulated engine error') && ui.runButton.enabled && !ui.cancelButton.enabled
        ui.performAnalysis={a,b,c,d,job->new File('tests/fixtures/simpleca').absoluteFile}
        ui.runButton.doClick()
    }
    finish()
    onUI {
        assert ui.metrics.rowCount==6 : ui.status.text
        assert ui.calls.rowCount==6 && ui.counts.rowCount==1
        assert ui.progress.value==100 && ui.folderButton.enabled && ui.savePlots.enabled
        // Headline tiles and outline styles repeat the saved Responder counts and calls; they add no new calculation.
        assert ui.responding.value.text=='67 %' && ui.responding.caption.text.contains('2 of 3 evaluable')
        assert ui.registered.value.text=='6' && ui.excluded.value.text=='2' && ui.unscorable.value.text=='1'
        assert ui.cellPicture.calls.Cell001=='Responded' && ui.cellPicture.calls.Cell002=='Did not respond'
        // Clicking a cell on the picture inside the window selects it; zooming keeps the click on the same cell.
        ui.validate()
        def clickCell={int index->double[] f=ui.cellPicture.frame();def b=ui.measurement.rois[index].bounds
            int x=(int)(f[0]+b.centerX*f[2]),y=(int)(f[1]+b.centerY*f[2])
            [java.awt.event.MouseEvent.MOUSE_PRESSED,java.awt.event.MouseEvent.MOUSE_RELEASED].each{ui.cellPicture.dispatchEvent(new java.awt.event.MouseEvent(ui.cellPicture,it,System.currentTimeMillis(),0,x,y,1,false,java.awt.event.MouseEvent.BUTTON1))}}
        assert ui.cellPicture.width>50 && ui.cellPicture.height>50
        clickCell(3);assert ui.cells.selectedItem=='Cell004' && ui.cellPicture.selected==3
        ui.cellPicture.dispatchEvent(new java.awt.event.MouseWheelEvent(ui.cellPicture,java.awt.event.MouseEvent.MOUSE_WHEEL,System.currentTimeMillis(),0,40,40,0,false,java.awt.event.MouseWheelEvent.WHEEL_UNIT_SCROLL,3,-3))
        assert ui.cellPicture.magnify>1d
        ui.cells.selectedItem='Cell006';clickCell(1);assert ui.cells.selectedItem=='Cell002'
        ui.nextCell.doClick();assert ui.cells.selectedItem=='Cell003';ui.previousCell.doClick();ui.previousCell.doClick();assert ui.cells.selectedItem=='Cell001'
        ui.previousCell.doClick();assert ui.cells.selectedItem=='Cell006'
        // One line says what was decided for the selected cell, in the saved table's own words and numbers.
        ui.cells.selectedItem='Cell001'
        def said=ui.cellStatus.text.replaceAll('<[^>]+>','')
        assert said.contains('Cell001: Responded') && said.contains('sustained threshold crossing') && said.contains('largest rise 1.00 at 6 s (needed 0.10)') : said
        ui.cells.selectedItem='Cell003';assert ui.cellStatus.text.contains('Cell003: Excluded') && ui.cellStatus.text.contains('invalid raw fluorescence baseline')
        // Every cell's trace is drawn behind the mean, coloured by its call; the groups list and draw only their cells.
        assert (0..<ui.group.itemCount).collect{ui.group.getItemAt(it)}==['All cells (6)','Responded (2)','Did not respond (1)','Excluded or unscorable (3)']
        assert ui.plot.showAll && ui.plot.all.length==4 && ui.plot.allCalls.count{it=='Responded'}==2 && ui.registered.active
        ui.group.selectedIndex=1
        assert (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell001','Cell005'] && ui.plot.all.length==2 && ui.cellPicture.group==1 && ui.responding.active && !ui.registered.active
        ui.nextCell.doClick();assert ui.cells.selectedItem=='Cell005';ui.nextCell.doClick();assert ui.cells.selectedItem=='Cell001'
        ui.excluded.pressed.call()
        assert ui.group.selectedIndex==3 && (0..<ui.cells.itemCount).collect{ui.cells.getItemAt(it)}==['Cell003','Cell004','Cell006'] && ui.plot.all.length==1
        clickCell(1);assert ui.group.selectedIndex==0 && ui.cells.selectedItem=='Cell002' : 'Clicking a cell outside the group shows all cells again'
        ui.allTraces.doClick();assert !ui.plot.showAll;ui.allTraces.doClick();assert ui.plot.showAll
        // Time: the slider, the plot and Play move one cursor, and in a movie view the picture shows that frame of the recording.
        int brightest=(0..<ui.plot.mean.length).max{ui.plot.mean[it]}+1
        assert ui.view.selectedIndex==1 && ui.frame==brightest && ui.cellPicture.image.is(ui.frameImage) && ui.plot.moment==ui.plot.time[brightest-1]
        def gray={int x,int y->ui.frameImage.getRGB(x,y)&255}
        ui.showFrame(15);assert ui.slider.value==15 && ui.plot.moment==ui.plot.time[14] && ui.timeLabel.text.contains('7.0 s') && ui.timeLabel.text.contains('frame 15 of 30')
        assert gray(40,60)>200 && gray(100,60)<60 : 'Frame 15 minus the resting frames: only the responding square is bright'
        ui.showFrame(3);assert gray(40,60)<60 && ui.cellPicture.note.contains('frame 3 minus the mean of frames 1–5')
        ui.plot.seek.call(7.2d);assert ui.frame==15 && gray(40,60)>200
        ui.slider.value=25;assert ui.frame==25 && gray(40,60)<60
        ui.view.selectedIndex=2;ui.showFrame(15);assert gray(40,60)>gray(100,60) && ui.cellPicture.note=='frame 15 as recorded'
        ui.view.selectedIndex=0;assert ui.cellPicture.image.is(ui.snapshot) && !ui.cellPicture.onGray
        ui.view.selectedIndex=1
        ui.showFrame(28);ui.playButton.doClick();assert ui.player.running && ui.playButton.text.contains('Pause');ui.playButton.doClick();assert !ui.player.running
        ui.outlinesBox.doClick();assert !ui.cellPicture.outlines;ui.outlinesBox.doClick()
        ui.cellPicture.whole();ui.zoomButton.doClick();assert ui.cellPicture.magnify>1d;ui.wholeButton.doClick();assert ui.cellPicture.magnify==1d
        assert recording.getRoi()==null && recording.getStackSize()==30 && recording.getStack().getProcessor(15).get(40,60)==140 : 'The recording is only read'
        ui.cells.selectedItem='Cell003'
        assert ui.plot.caption.contains('unavailable') && ui.plot.selected.length==0
        ui.cells.selectedItem='Cell001'
        assert ui.plot.selected.max()==1d
        assert ui.measurement.picture.getRoi().bounds==ui.measurement.rois[0].bounds
        assert original.getRoi()==null
        // Heatmap rows: responders first (strongest on top), then non-responders, then cells without a call; or as measured.
        assert ui.heatmap.names==['Cell001','Cell005','Cell002','Cell003','Cell004','Cell006'] && ui.heatmap.calls[0]=='Responded' && ui.heatmap.calls[2]=='Did not respond'
        assert ui.heatmap.limit()==1d && ui.heatmap.preferredSize.height>=6*ui.heatmap.ROW
        ui.heatmapOrder.selectedIndex=1
        assert ui.heatmap.names==(1..6).collect{String.format('Cell%03d',it)}
        ui.tabs.selectedIndex=1;ui.heatmap.choose(4)
        assert ui.cells.selectedItem=='Cell005' && ui.tabs.selectedIndex==0 && ui.heatmap.selected=='Cell005'
        ui.metrics.setRowSelectionInterval(1,1)
        assert ui.cells.selectedItem=='Cell002'
        ui.results=new File('build/analysis-gui').absoluteFile;ui.results.mkdirs()
        ui.exportPlots()
        assert new File(ui.results,'trace.png').length()>0 && new File(ui.results,'heatmap.png').length()>0 : ui.status.text
        // A failed or cancelled rerun must not replace the last successful tables.
        ui.performAnalysis={a,b,c,d,job->while(!job.cancelled)Thread.sleep(20);job.check()}
        ui.runButton.doClick();ui.cancelButton.doClick()
    }
    finish()
    onUI {
        assert ui.status.text.startsWith('Cancelled') && ui.metrics.rowCount==6 && ui.runButton.enabled && ui.savePlots.enabled
        ui.fields.sampling_hz.text='0';ui.runButton.doClick();assert !ui.busy && ui.status.text.contains('frame rate')
        ui.setSize(400,350);ui.validate()
        assert ui.contentPane instanceof JScrollPane
    }
    println 'PASS: analysis settings with exact untouched values, responsive worker, failure/retry/cancel, real result tables, headline tiles, cell groups and status line, every-trace plot, movie frames linked to the time cursor, in-window cell picture, trace/ROI links, exclusions, ordered heatmap and PNG export'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();original?.close();recording?.close();imagej?.dispose()}}
System.exit(0)

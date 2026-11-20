import javax.swing.*
import ij.*
import ij.gui.*
import ij.process.ColorProcessor

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Window=loader.loadClass('CalciumAnalysisWindow')
def Protocol=loader.loadClass('CalciumResponseProtocol')
def ui, original, imagej, recording
def onUI={Closure c->SwingUtilities.invokeAndWait(c as Runnable)}
def finish={->long deadline=System.nanoTime()+8000000000L;while(ui.busy&&System.nanoTime()<deadline)Thread.sleep(20);onUI{};assert !ui.busy}
try {
    onUI {
        imagej=new ImageJ(ImageJ.NO_SHOW)
        def Server=loader.loadClass('CalciumWebServer'),controls=new JFrame('Browser launch controls')
        try{controls.setAlwaysOnTop(true)
            try{Server.openBrowser(controls,new java.net.URI('http://127.0.0.1/'),{address->throw new IOException('Browser launch failed')});assert false}catch(IOException expected){}
            assert controls.alwaysOnTop : 'A failed browser launch preserves the controls for native fallback'
            boolean opened=false;Server.openBrowser(controls,new java.net.URI('http://127.0.0.1/'),{address->opened=true})
            assert opened && !controls.alwaysOnTop : 'Successful browser launch must not leave Fiji pinned over the page'
        }finally{controls.dispose()}

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
        assert ui.kcl.text=='Also score a separate positive control'
        // The original recording is visible and scrubbable even before analysis.
        assert (0..<ui.view.itemCount).collect{ui.view.getItemAt(it)}==['Static colour overview','Movie: change from rest','Movie: as recorded']
        assert ui.view.selectedIndex==2 && ui.cellPicture.image.is(ui.frameImage) && ui.slider.maximum==30
        ui.slider.value=15;assert ui.frame==15 && ui.timeLabel.text.contains('7.0 s')
        // Suggested settings are shown rounded, but an untouched box still sends the exact value.
        def precise=Window.newInstance([csv:new File('tests/fixtures/simpleca/metrics.csv'),folder:new File('build/analysis-gui'),recording:'Precise',frames:300,hz:1.0000537066691941d,
            baselineStart:55,baselineEnd:74,names:[],rois:new Roi[0]])
        try {
            assert precise.fields.sampling_hz.text=='1'
            assert precise.fields.sampling_hz.toolTipText.contains('1.0000537066691941')
            def startColumn=Protocol.COLUMNS.indexOf('Start_s'),endColumn=Protocol.COLUMNS.indexOf('End_s')
            assert precise.windowModel.getValueAt(0,startColumn)==74/1.0000537066691941d
            assert precise.windowModel.getValueAt(0,endColumn)==299/1.0000537066691941d
            def Settings=loader.loadClass('CalciumAnalysisSettings')
            def offered=Settings.suggested([frames:300,hz:1.0000537066691941d,baselineStart:55,baselineEnd:74])
            offered.each{key,value->assert precise.fields[key].text==Settings.shown(value) : key}
            def webDraft=Settings.draft([frames:300,hz:1.0000537066691941d,baselineStart:55,baselineEnd:74])
            assert webDraft.protocol.windows[0].Start_s==precise.windowModel.getValueAt(0,startColumn)
            precise.confirmed.selected=true
            def sent=precise.settings(),p=new groovy.json.JsonSlurper().parseText(sent.protocol_json)
            assert sent.sampling_hz=='1.0000537066691941' && p.windows[0].Start_s==74/1.0000537066691941d && p.windows[0].End_s==299/1.0000537066691941d
            assert p.windows[0].Floor_local_fraction==0.10d && p.windows[0].SD_multiplier==3d && sent.baseline_start=='55' && sent.baseline_end=='74'
            precise.windowModel.setValueAt(80d,0,startColumn)
            assert new groovy.json.JsonSlurper().parseText(precise.settings().protocol_json).windows[0].Start_s==80d
            assert precise.windowModel.rowCount==1;precise.kcl.doClick();assert precise.windowModel.rowCount==2;precise.kcl.doClick();assert precise.windowModel.rowCount==1
            assert precise.windowModel.isCellEditable(0,Protocol.COLUMNS.indexOf('Min_duration_s'))
            assert !precise.windowModel.isCellEditable(0,Protocol.COLUMNS.indexOf('Consecutive_samples'))
            precise.preset.selectedIndex=2;precise.kcl.doClick()
            assert precise.protocol().schema_version==1 && precise.protocol().windows.size()==2
            assert precise.windowModel.getValueAt(1,Protocol.COLUMNS.indexOf('Duration_mode'))=='samples'
            assert precise.runButton.font.family==UIManager.getFont('Label.font').family : 'Use the interface font, not a fixed family'
        } finally {precise.dispose()}
        ui.runButton.doClick();assert !ui.busy && ui.status.text.contains('Confirm')
        ui.confirmed.selected=true
        ['Baseline_start_s':0.5d,'Baseline_end_s':3d,'Stimulus_onset_s':3d,'Start_s':3d,'End_s':14.5d].each{k,v->ui.windowModel.setValueAt(v,0,Protocol.COLUMNS.indexOf(k))}
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
        // Keep render artifacts for visual review of the real merged native UI.
        new File('build/analysis-gui').mkdirs()
        [1:'traces',0:'protocol'].each{index,label->
            ui.tabs.selectedIndex=index;ui.validate()
            def screenshot=new java.awt.image.BufferedImage(ui.width,ui.height,java.awt.image.BufferedImage.TYPE_INT_RGB)
            def g=screenshot.createGraphics();try{ui.paint(g)}finally{g.dispose()}
            javax.imageio.ImageIO.write(screenshot,'png',new File('build/analysis-gui/window-'+label+'.png'))
        }
        ui.tabs.selectedIndex=1
        // Headline tiles and outline styles repeat the saved Responder counts and calls; they add no new calculation.
        assert ui.responding.value.text=='67 %' && ui.responding.caption.text.contains('2 of 3 evaluable')
        assert ui.responding.caption.text.contains('Capsaicin')
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
        assert ui.view.selectedIndex==2 && ui.frame==brightest && ui.cellPicture.image.is(ui.frameImage) && ui.plot.moment==ui.plot.time[brightest-1]
        ui.view.selectedIndex=1
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
        assert ui.heatmap.names==['Cell001','Cell005','Cell002','Cell003','Cell004','Cell006'] && ui.heatmap.calls[0]=='Responded' && ui.heatmap.calls[2]=='Did not respond'
        assert ui.heatmap.limit()==1d && ui.heatmap.preferredSize.height>=6*ui.heatmap.ROW
        ui.heatmapOrder.selectedIndex=1
        assert ui.heatmap.names==(1..6).collect{String.format('Cell%03d',it)}
        ui.tabs.selectedIndex=2;ui.heatmap.choose(4)
        assert ui.cells.selectedItem=='Cell005' && ui.tabs.selectedIndex==1 && ui.heatmap.selected=='Cell005'
        ui.metrics.setRowSelectionInterval(1,1)
        assert ui.cells.selectedItem=='Cell002'
        // Saved irregular timestamps drive source-frame selection and resting-image frames.
        def Csv=loader.loadClass('CalciumCsv'),irregular=new File('build/analysis-gui/irregular-completed').absoluteFile
        irregular.mkdirs();new File('tests/fixtures/simpleca').listFiles().findAll{it.isFile()}.each{new File(irregular,it.name).bytes=it.bytes}
        def savedTimes=(0..<30).collect{it/2d+(it>=15?5d:0d)}
        ['normalized_traces.csv','mean_trace.csv'].each{filename->def rows=Csv.read(new File(irregular,filename));rows.drop(1).eachWithIndex{row,i->row[0]=savedTimes[i]};Csv.write(new File(irregular,filename),rows)}
        def responseTable=Csv.read(new File(irregular,'responder_calls.csv')),header=responseTable[0]
        def controlRows=responseTable.drop(1).collect{originalRow->def row=new ArrayList(originalRow);['Role':'Control','Stimulus':'Ionomycin','Threshold_increase':'0.4','Floor_threshold':'0.4','Baseline_mean':'0.5'].each{k,v->row[header.indexOf(k)]=v};row}
        Csv.write(new File(irregular,'responder_calls.csv'),responseTable+controlRows)
        ui.loadResults(irregular);ui.showFrame(16)
        assert ui.plot.moment==12.5d && ui.frame==16 && ui.slider.value==16 && ui.restFirst==1 && ui.restLast==5
        ui.plot.seek.call(9d);assert ui.frame==15 && ui.plot.moment==7d : 'A gap in acquisition does not manufacture movie frames'
        ui.fields.sampling_hz.text='9';ui.fields.baseline_start.text='2';ui.fields.baseline_end.text='9'
        assert ui.plot.time.toList()==savedTimes && ui.restFirst==1 && ui.restLast==5
        ui.showFrame(16);assert ui.plot.moment==12.5d
        // Cached static layers must refresh when selecting a different saved window at an unchanged canvas size.
        ui.cells.selectedItem='Cell001';ui.stimuli.selectedItem='Capsaicin'
        def render={->def image=new java.awt.image.BufferedImage(ui.plot.width,ui.plot.height,java.awt.image.BufferedImage.TYPE_INT_RGB);def graphics=image.createGraphics();try{ui.plot.paint(graphics)}finally{graphics.dispose()};((java.awt.image.DataBufferInt)image.raster.dataBuffer).data.clone() as int[]}
        int[] primaryPixels=render();assert Math.abs(ui.plot.referenceValues().threshold-0.1d)<1e-12
        ui.stimuli.selectedItem='Ionomycin';int[] controlPixels=render()
        assert ui.plot.responseCall.Role=='Control' && Math.abs(ui.plot.referenceValues().threshold-1.1d)<1e-12
        assert !java.util.Arrays.equals(primaryPixels,controlPixels) : 'Changing the saved window must invalidate the cached trace layer'
        ui.responding.pressed.call();assert ui.stimuli.selectedItem=='Capsaicin' && ui.group.selectedIndex==1 && ui.plot.responseCall.Role=='Stimulus'
        ui.showGroup(0)
        ui.view.selectedIndex=0;assert !ui.playButton.enabled
        ui.view.selectedIndex=2;assert ui.playButton.enabled
        ui.playButton.doClick();assert ui.player.running;ui.slider.value=3;assert !ui.player.running && ui.frame==3
        ui.playButton.doClick();ui.plot.seek.call(6d);assert !ui.player.running && ui.frame==13
        def completeResults=ui.results
        try{ui.loadResults(new File('build/analysis-gui/not-complete'));assert false}catch(IOException expected){assert expected.message.contains('COMPLETE')}
        assert ui.results==completeResults
        ui.fields.sampling_hz.text='2';ui.fields.baseline_start.text='1';ui.fields.baseline_end.text='5'
        ui.results=new File('build/analysis-gui').absoluteFile;ui.results.mkdirs()
        new File(ui.results,'analysis_manifest.csv').bytes=new File('tests/fixtures/simpleca/analysis_manifest.csv').bytes
        ui.exportPlots()
        assert new File(ui.results,'trace.png').length()>0 && new File(ui.results,'heatmap.png').length()>0 : ui.status.text
        assert new File(ui.results,'plot_provenance.json').isFile()
        // Editing a protocol never relabels the last completed response table or its exported plots.
        ui.windowModel.setValueAt(0.2d,0,Protocol.COLUMNS.indexOf('Floor_local_fraction'))
        assert ui.appliedNote.text.contains('Settings changed') && ui.calls.rowCount==12
        def portable=new File(ui.results,'portable.json');ui.saveProtocol(portable)
        ui.loadProtocol(portable)
        assert ui.protocol().schema_version==2 && ui.protocol().windows[0].Floor_local_fraction==0.2d
        assert !ui.confirmed.selected
        // A timing edit must never erase independently imported control rules.
        def sharedProtocol=new File('tests/fixtures/response-protocols/golden-protocol.json')
        ui.loadProtocol(sharedProtocol)
        ui.windowModel.setValueAt(0.35d,1,Protocol.COLUMNS.indexOf('Floor_local_fraction'))
        ui.windowModel.setValueAt(5d,1,Protocol.COLUMNS.indexOf('SD_multiplier'))
        ui.windowModel.setValueAt('samples',1,Protocol.COLUMNS.indexOf('Duration_mode'))
        ui.windowModel.setValueAt(4,1,Protocol.COLUMNS.indexOf('Consecutive_samples'))
        ui.saveProtocol(portable);ui.loadProtocol(portable)
        ui.windowModel.setValueAt(7d,0,Protocol.COLUMNS.indexOf('End_s'))
        def retained=ui.protocol().windows[1]
        assert retained.Floor_local_fraction==0.35d && retained.SD_multiplier==5d
        assert retained.Duration_mode=='samples' && retained.Consecutive_samples==4
        ui.sensitivityFile=sharedProtocol
        ui.preset.selectedIndex=2
        assert ui.protocol().schema_version==1 && ui.protocol().windows[1].Stimulus=='KCl'
        assert !ui.sensitivityButton.enabled && ui.sensitivityFile==null
        ['Duration_mode','Min_duration_s','Min_baseline_samples','Noise_floor_local_fraction','Expected_interval_s','Baseline_end_s','Stimulus_onset_s'].each{key->assert !ui.windowModel.isCellEditable(0,Protocol.COLUMNS.indexOf(key))}
        assert ui.windowModel.isCellEditable(0,Protocol.COLUMNS.indexOf('Consecutive_samples'))
        ui.windowModel.setValueAt(5.5d,0,Protocol.COLUMNS.indexOf('Start_s'))
        assert ui.windowModel.getValueAt(0,Protocol.COLUMNS.indexOf('Baseline_end_s'))==5.5d
        assert ui.windowModel.getValueAt(0,Protocol.COLUMNS.indexOf('Stimulus_onset_s'))==5.5d
        assert ui.protocol().windows[0].Start_s==5.5d
        ui.loadProtocol(portable)
        // Same path does not imply the same timing/diagnostic content. Both a
        // deliberate reload and a later external edit must mark applied state stale.
        def times=new File(ui.results,'clock.csv')
        times.text='Time_s\n'+(0..<30).collect{it/2d}.join('\n')+'\n'
        ui.loadTimes(times);ui.appliedFingerprint=ui.fingerprint();ui.changed()
        assert ui.appliedNote.text.contains('match the applied')
        times.text='Time_s\n'+(0..<30).collect{it/2d+0.01d}.join('\n')+'\n'
        ui.loadTimes(times);assert ui.appliedNote.text.contains('Settings changed')
        ui.appliedFingerprint=ui.fingerprint();times.text+='\n';ui.changed()
        assert ui.appliedNote.text.contains('Settings changed')
        ui.clearTimesButton.doClick()
        [false,true].each{negative->
            def diagnostic=new File(ui.results,negative?'negative.csv':'sensitivity.csv');diagnostic.text='Value\n1\n'
            ui.loadDiagnostic(diagnostic,negative);ui.appliedFingerprint=ui.fingerprint()
            diagnostic.text='Value\n2\n';ui.changed();assert ui.appliedNote.text.contains('Settings changed')
            ui.loadDiagnostic(diagnostic,negative);ui.appliedFingerprint=ui.fingerprint()
            diagnostic.delete();ui.changed();assert ui.appliedNote.text.contains('Settings changed')
            assert ui.fileSignature(diagnostic,'').current_sha256=='missing'
        }
        ui.clearDiagnosticsButton.doClick()
        ui.confirmed.selected=true
        ui.plot.responseCall=[Baseline_mean:'0.5',Floor_threshold:'0.1',Noise_threshold:'0.03',Threshold_increase:'0.1']
        assert Math.abs(ui.plot.referenceValues().threshold-0.65d)<1e-12
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
    println 'PASS: merged portable protocols, exact settings, applied timing and threshold overlays, native movie playback, groups, scientific tables, failure/cancel and plot provenance'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();original?.close();recording?.close();imagej?.dispose()}}
System.exit(0)

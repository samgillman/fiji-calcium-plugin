import javax.swing.*
import ij.*
import ij.gui.*
import ij.process.ColorProcessor

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Window=loader.loadClass('CalciumAnalysisWindow')
def Protocol=loader.loadClass('CalciumResponseProtocol')
def ui, original, imagej
def onUI={Closure c->SwingUtilities.invokeAndWait(c as Runnable)}
def finish={->long deadline=System.nanoTime()+8000000000L;while(ui.busy&&System.nanoTime()<deadline)Thread.sleep(20);onUI{};assert !ui.busy}
try {
    onUI {
        imagej=new ImageJ(ImageJ.NO_SHOW)
        original=new ImagePlus('Original colour picture',new ColorProcessor(250,150))
        ui=Window.newInstance([csv:new File('tests/fixtures/simpleca/metrics.csv'),folder:new File('build/analysis-gui'),recording:'Fixture',frames:30,hz:2d,baselineStart:1,baselineEnd:5,
            picture:original,names:(1..6).collect{String.format('Cell%03d',it)},rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]])
        ui.visible=true
        assert !ui.measurement.picture.is(original) : 'Analysis must not edit Step 2 cells'
        assert ui.measurement.picture.getWindow()==null && ui.cellPicture.image!=null : 'The cell picture belongs inside the analysis window'
        assert ui.responding.value.text=='–' && !ui.previousCell.enabled
        // Suggested settings are shown rounded, but an untouched box still sends the exact value.
        def precise=Window.newInstance([csv:new File('tests/fixtures/simpleca/metrics.csv'),folder:new File('build/analysis-gui'),recording:'Precise',frames:300,hz:1.0000537066691941d,
            baselineStart:55,baselineEnd:74,names:[],rois:new Roi[0]])
        try {
            assert precise.fields.sampling_hz.text=='1'
            assert precise.fields.sampling_hz.toolTipText.contains('1.0000537066691941')
            def startColumn=Protocol.COLUMNS.indexOf('Start_s'),endColumn=Protocol.COLUMNS.indexOf('End_s')
            assert precise.windowModel.getValueAt(0,startColumn)==73/1.0000537066691941d
            assert precise.windowModel.getValueAt(0,endColumn)==299/1.0000537066691941d
            precise.confirmed.selected=true
            def sent=precise.settings(),p=new groovy.json.JsonSlurper().parseText(sent.protocol_json)
            assert sent.sampling_hz=='1.0000537066691941' && p.windows[0].Start_s==73/1.0000537066691941d && p.windows[0].End_s==299/1.0000537066691941d
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
        ui.cells.selectedItem='Cell003'
        assert ui.plot.caption.contains('unavailable') && ui.plot.selected.length==0
        ui.cells.selectedItem='Cell001'
        assert ui.plot.selected.max()==1d
        assert ui.measurement.picture.getRoi().bounds==ui.measurement.rois[0].bounds
        assert original.getRoi()==null
        ui.tabs.selectedIndex=2;ui.heatmap.choose(4)
        assert ui.cells.selectedItem=='Cell005' && ui.tabs.selectedIndex==1
        ui.metrics.setRowSelectionInterval(1,1)
        assert ui.cells.selectedItem=='Cell002'
        ui.results=new File('build/analysis-gui').absoluteFile;ui.results.mkdirs()
        new File(ui.results,'analysis_manifest.csv').bytes=new File('tests/fixtures/simpleca/analysis_manifest.csv').bytes
        ui.exportPlots()
        assert new File(ui.results,'trace.png').length()>0 && new File(ui.results,'heatmap.png').length()>0 : ui.status.text
        assert new File(ui.results,'plot_provenance.json').isFile()
        // Editing a protocol never relabels the last completed response table or its exported plots.
        ui.windowModel.setValueAt(0.2d,0,Protocol.COLUMNS.indexOf('Floor_local_fraction'))
        assert ui.appliedNote.text.contains('Settings changed') && ui.calls.rowCount==6
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
    println 'PASS: portable protocol editing and exact values, pending settings, threshold transform, in-window image/navigation, result counts, failure/retry/cancel and plot provenance'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();original?.close();imagej?.dispose()}}
System.exit(0)

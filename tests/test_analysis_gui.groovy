import javax.swing.*
import ij.*
import ij.gui.*
import ij.process.ColorProcessor

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Window=loader.loadClass('CalciumAnalysisWindow')
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
            assert precise.fields.sampling_hz.text=='1' && precise.fields.stimulus_start.text=='72.996' && precise.fields.stimulus_end.text=='298.984'
            assert precise.fields.sampling_hz.toolTipText.contains('1.0000537066691941')
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
        ui.cells.selectedItem='Cell003'
        assert ui.plot.caption.contains('unavailable') && ui.plot.selected.length==0
        ui.cells.selectedItem='Cell001'
        assert ui.plot.selected.max()==1d
        assert ui.measurement.picture.getRoi().bounds==ui.measurement.rois[0].bounds
        assert original.getRoi()==null
        ui.tabs.selectedIndex=1;ui.heatmap.choose(4)
        assert ui.cells.selectedItem=='Cell005' && ui.tabs.selectedIndex==0
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
    println 'PASS: analysis settings with exact untouched values, responsive worker, failure/retry/cancel, real result tables, headline tiles, in-window cell picture, trace/ROI links, exclusions, heatmap selection and PNG export'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();original?.close();imagej?.dispose()}}
System.exit(0)

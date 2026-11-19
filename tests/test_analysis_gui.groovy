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
    println 'PASS: analysis settings, responsive worker, failure/retry/cancel, real result tables, trace/ROI links, exclusions, heatmap selection and PNG export'
} catch(Throwable t){t.printStackTrace();System.exit(1)}
finally {onUI{ui?.dispose();original?.close();imagej?.dispose()}}
System.exit(0)

// Step 2's picture views, run through the whole plugin on a small synthetic recording: magenta can be turned on after the
// recording is open, the baseline view works, the window keeps its outlines and clicks, and the way back to the first screen.
import ij.*
import ij.gui.*
import ij.plugin.frame.RoiManager
import ij.process.*
import javax.swing.*
import java.awt.Window
import java.util.concurrent.TimeUnit

def onUI = { Closure c -> SwingUtilities.invokeAndWait(c as Runnable) }
def waitFor = { String what, Closure ready -> long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(20); boolean ok = false
    while (!ok && System.nanoTime() < end) { onUI { ok = ready() }; if (!ok) Thread.sleep(25) }; assert ok : 'Timed out waiting for ' + what }
def everything; everything = { c -> [c] + (c instanceof java.awt.Container ? c.components.collectMany { everything(it) } : []) }
def active = { ImagePlus p -> (p as CompositeImage).getActiveChannels()[0..<p.getNChannels()].collect { it ? '1' : '0' }.join('') }

// 30 frames of 144 x 96: round cells at rest, and the left part of the field brightens from frame 15 on
int w = 144, h = 96, n = 30; def rnd = new Random(7); def stack = new ImageStack(w, h)
def cells = [[24, 30], [60, 60], [100, 30], [120, 70]]
(1..n).each { int f -> def ip = new ShortProcessor(w, h)
    for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) {
        double v = 100 + rnd.nextGaussian(); cells.each { c -> if (Math.hypot(x - c[0], y - c[1]) < 8) v += 30 }
        if (f >= 15 && x < 96) v += 25
        ip.set(x, y, (int) Math.round(v)) }
    stack.addSlice('' + f, ip) }
def movie = new ImagePlus('views.tif', stack); def before = movie.getStack().getProcessor(20).getPixelsCopy()
def imagej, hooks, roiManager
try {
    onUI { imagej = new ImageJ(ImageJ.NO_SHOW); movie.show(); roiManager = new RoiManager(false) }
    // the first screen: the open recording is offered; "look for a response" is unticked here, so the views start from all frames
    Thread.start {
        Window dialog = null; long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(20)
        while (dialog == null && System.nanoTime() < end) { Thread.sleep(50); dialog = Window.getWindows().find { it instanceof JDialog && it.visible && it.title?.startsWith('Calcium imaging') } }
        SwingUtilities.invokeLater { def box = everything(dialog).find { it instanceof javax.swing.JCheckBox }; assert box.selected : 'ticked by default'; box.setSelected(false); everything(dialog).find { it instanceof JButton && it.text == 'Next' }?.doClick() }
    }
    def binding = new Binding([CALCIUM_TEST_GUI: true])
    new GroovyShell(getClass().classLoader, binding).evaluate(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
    hooks = binding.getVariable('hooks'); def steps = hooks.steps; def picture = steps.colourPic
    def canvas, listeners
    onUI {
        assert !steps.looked && !steps.found && hooks.selectionBasis() == 'All-frames structural review'
        assert picture instanceof CompositeImage && picture.getNChannels() == 3 && active(picture) == '110' : 'Room for magenta, switched off'
        assert steps.colourComposite(false).getNChannels() == 2 : 'The saved colour picture has no empty magenta channel'
        assert hooks.recordingNote.text.contains('response not looked for') && hooks.selectionNote.text.contains('average of all frames')
        assert hooks.viewAllBtn.getUI().chosen && !hooks.viewResponseBtn.getUI().chosen
        canvas = picture.getCanvas(); listeners = canvas.getMouseListeners() as List
        hooks.viewResponseBtn.doClick()
    }
    // magenta after the recording is open: the response is looked for now, and shown on the same window
    waitFor('the response picture') { hooks.viewResponseBtn.enabled && steps.looked }
    onUI {
        assert steps.found && active(picture) == '111' && hooks.selectionBasis() == 'Response-assisted review'
        float[] magenta = picture.getStack().getPixels(3), response = steps.respPic.getProcessor().getPixels()
        assert magenta == response : 'The magenta channel is the response picture'
        assert picture.getCanvas().is(canvas) && (picture.getCanvas().getMouseListeners() as List).containsAll(listeners) : 'Same window, so clicks on cells still work'
        assert hooks.recordingNote.text.contains('response near frame') && hooks.viewResponseBtn.getUI().chosen && hooks.selectionNote.text.contains('response-assisted')
        assert steps.colourComposite(false).getNChannels() == 3 : 'Once found, the saved picture has the response, as when it was looked for at the start'
        hooks.viewAllBtn.doClick()
        assert active(picture) == '110' && hooks.selectionBasis() == 'All-frames structural review' && hooks.viewAllBtn.getUI().chosen
        float[] green = picture.getStack().getPixels(2), all = steps.allPic.getProcessor().getPixels(); assert green == all
    }
    // the baseline view: frames before the stimulus, offered as the default range
    SwingUtilities.invokeLater { hooks.viewBaselineBtn.doClick() }
    def prompt = null; long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(20)
    while (prompt == null && System.nanoTime() < end) { Thread.sleep(25); prompt = Window.getWindows().find { it instanceof GenericDialog && it.visible } }
    assert prompt != null
    def range = null; onUI { range = prompt.getNumericFields().collect { it.text as int } }
    assert range == [steps.b0, steps.b1] : 'The detected baseline frames are offered'
    onUI { prompt.actionPerformed(new java.awt.event.ActionEvent(prompt.getButtons()[0], java.awt.event.ActionEvent.ACTION_PERFORMED, '')) }
    waitFor('the baseline picture') { hooks.selectionBasis().startsWith('Baseline') && hooks.viewBaselineBtn.enabled }
    onUI {
        assert hooks.selectionBasis() == 'Baseline frames ' + steps.b0 + '-' + steps.b1 && active(picture) == '110' && hooks.viewBaselineBtn.getUI().chosen
        float[] green = picture.getStack().getPixels(2), baseline = steps.baselinePicture(steps.b0, steps.b1)
        assert green == baseline : 'The picture shows the baseline frames'
        assert picture.getCanvas().is(canvas)
        // moving between steps: forward and back
        def (s1, s2, s3, s4) = hooks.cards
        hooks.continueBtn.doClick(); assert s3.expanded && !s2.expanded
        hooks.backToCellsBtn.doClick(); assert s2.expanded && !s3.expanded
        s4.onToggle.call(s4); hooks.backToMeasureBtn.doClick(); assert s3.expanded && !s4.expanded
        // back to the first screen, clearing the outlines: the window and the colour picture close, the recording stays
        hooks.rm.addRoi(new OvalRoi(16, 22, 16, 16)); assert hooks.rm.getCount() == 1
        hooks.changeRecording(0)
        assert hooks.rm.getCount() == 0 && !hooks.frame.isDisplayable() && picture.getWindow() == null
    }
    assert movie.getStack().getProcessor(20).getPixelsCopy() == before && movie.getWindow() != null : 'The recording is only read'
    println 'PASS: Step 2 picture views: magenta turned on after opening (same window and clicks), all-frames and baseline views, saved picture unchanged, back buttons, change recording with outlines cleared'
} catch (Throwable t) { t.printStackTrace(); System.exit(1) }
finally { onUI { hooks?.detectionProgress?.close(); hooks?.frame?.dispose(); movie.close(); imagej?.dispose() } }
System.exit(0)

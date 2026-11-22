// Step 2's picture views, run through the whole plugin on a small synthetic recording: blind selection is the default (frames
// before the stimulus only, nothing of the response in the window or in the picture for Cellpose, no traces), magenta can be
// turned on after the recording is open, the record of how cells were chosen stays honest, the window keeps its outlines and
// clicks, and the way back to the first screen.
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
    // the first screen: the open recording is offered and blind selection is the default; then the frames before the stimulus are confirmed
    def offeredFrames = null
    Thread.start {
        Window dialog = null; long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(20)
        while (dialog == null && System.nanoTime() < end) { Thread.sleep(50); dialog = Window.getWindows().find { it instanceof JDialog && it.visible && it.title?.startsWith('Ca\u00b2\u207a Imaging') } }
        SwingUtilities.invokeLater { everything(dialog).find { it instanceof JButton && it.text == 'Next' }?.doClick() }
        GenericDialog frames = null; end = System.nanoTime() + TimeUnit.SECONDS.toNanos(20)
        while (frames == null && System.nanoTime() < end) { Thread.sleep(25); frames = Window.getWindows().find { it instanceof GenericDialog && it.visible && it.title == 'Frames before the response' } }
        SwingUtilities.invokeAndWait { offeredFrames = frames.getNumericFields().collect { it.text as int }; frames.actionPerformed(new java.awt.event.ActionEvent(frames.getButtons()[0], java.awt.event.ActionEvent.ACTION_PERFORMED, '')) }
    }
    def binding = new Binding([CALCIUM_TEST_GUI: true])
    new GroovyShell(getClass().classLoader, binding).evaluate(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
    hooks = binding.getVariable('hooks'); def steps = hooks.steps; def picture = steps.colourPic
    def canvas, listeners
    onUI {
        // blind: the response was looked for only to offer the frames before it; it is in neither the window nor the picture for Cellpose
        int offered = steps.proposedBlindWindow()[1]
        assert steps.looked && steps.found && !steps.riseIsMovement() && offeredFrames == [1, offered] && offered < steps.onset && offered <= steps.b1 : 'The plugin offers the frames before the rise it found, ending well before it'
        assert steps.blind && steps.blindFirst == 1 && steps.blindLast == offered && hooks.selectionBasis() == 'Before the response, blind (frames 1-' + offered + ')'
        // the Next line names the one thing to do; adding responders is offered only once there are blind cells
        assert hooks.nextLine.text.contains('Next:') && hooks.nextLine.text.contains('Suggest cells automatically') && hooks.addRespondersBtn.visible && !hooks.addRespondersBtn.enabled
        assert hooks.acceptAllBtn.text == 'Accept likely (0)'
        assert picture instanceof CompositeImage && picture.getNChannels() == 3 && active(picture) == '110' : 'Room for magenta, switched off'
        assert (picture.getStack().getPixels(3) as float[]).every { it == 0f } : 'Nothing of the response is held in the window while choosing blind'
        assert (picture.getStack().getPixels(2) as float[]) == steps.framesPicture(1, offered) : 'The picture shows the frames before the response'
        int[] forCellpose = steps.outlinePicture().getProcessor().getPixels()
        assert forCellpose.every { ((it >> 16) & 255) == 0 && (it & 255) == 0 } && forCellpose.any { ((it >> 8) & 255) > 0 } : 'Cellpose gets the picture of the frames before the response only, no magenta'
        assert steps.outlinePicture().getProcessor().getPixels() == steps.blindOutlinePicture().getProcessor().getPixels()
        // the two blind pictures for Cellpose: the same pre-stimulus picture, capped at 2 and at 8 gray levels
        assert steps.blindOutlinePicture(steps.BLIND_CAP).getProcessor().getPixels() == steps.blindOutlinePicture().getProcessor().getPixels()
        float[] high = steps.blindDetectionPicture(steps.BLIND_CAP_HIGH)
        assert (0..<high.length).every { int i -> Math.abs(high[i] - Math.max(0d, Math.min(1d, steps.blindPic[i] / 8d))) < 1e-6 } && high != steps.blindDetectionPicture()
        assert steps.blindRecipe().contains('from 0 to 2 and from 0 to 8 gray levels') && steps.blindRecipe().contains('less than 10 %')
        assert hooks.viewBaselineBtn.getUI().chosen && !hooks.viewResponseBtn.getUI().chosen && hooks.selectionNote.text.contains('Blind')
        assert hooks.recordingNote.text.contains('blind') && !hooks.recordingNote.text.contains('magenta')
        assert hooks.sensBox.visible && !hooks.respondsLegend.visible && hooks.trace.vals == null && hooks.trace.caption.contains('hidden')
        assert !hooks.trace.visible : 'No empty frame for a trace that is not shown'
        // the picture choice is one joined switch with plain names; the four actions are always there, gray when they cannot be pressed
        assert [hooks.viewAllBtn, hooks.viewBaselineBtn, hooks.viewResponseBtn]*.text == ['All frames', 'Before response', 'Response'] && [hooks.viewAllBtn, hooks.viewBaselineBtn, hooks.viewResponseBtn].every { it.parent.is(hooks.viewSwitch) }
        assert [hooks.suggestBtn, hooks.acceptAllBtn, hooks.addRespondersBtn, hooks.continueBtn].every { it.visible } && !hooks.acceptAllBtn.enabled && hooks.suggestBtn.getUI().primary
        assert hooks.recordingName.whole == 'views.tif' && hooks.showSwitches.keySet() == ['cells', 'g0', 'g1', 'g2', 'background', 'g3'] as Set
        // hovering over a cell shows what it is, but no trace
        hooks.rm.addRoi(new OvalRoi(16, 22, 16, 16)); hooks.refresh(true)
        assert hooks.hoverAt(24, 30).contains('Your cell 1') && hooks.trace.vals == null && !hooks.responseSeen() : 'No trace while choosing blind'
        assert hooks.basisFor(hooks.cellRois()) == 'Before the response, blind (frames 1-' + offered + ')'
        assert hooks.nextLine.text.contains('Add cells that respond') && hooks.nextLine.text.contains('Continue to measure') && hooks.addRespondersBtn.enabled : 'With blind cells chosen, the next step is offered: ' + hooks.nextLine.text
        // this recording has a response, so cells that only show when they respond have no outline yet: adding them is the next step
        assert hooks.nextLine.text.contains('have no outline yet') && hooks.nextLine.text.indexOf('Add cells that respond') < hooks.nextLine.text.indexOf('Continue to measure') : hooks.nextLine.text
        assert hooks.addRespondersBtn.getUI().primary && !hooks.continueBtn.getUI().primary : 'Add cells that respond is the one blue button, not Continue to measure'
        assert !hooks.addMissedBtn.visible : 'nothing is measured yet'
        // the legend switches kinds off and on: a click hides or shows one kind, a double-click shows only that kind
        int drawnAll = picture.getOverlay().size()
        hooks.legendClicked('cells', false); int withoutCells = picture.getOverlay().size()
        assert withoutCells < drawnAll && hooks.hiddenKinds == ['cells'] as Set : 'Your cells can be hidden'
        assert !hooks.showSwitches.cells.part.on && hooks.showSwitches.background.part.on && hooks.hiddenLine.text.contains('Hidden on the picture now: your cells.') && !hooks.nextLine.text.contains('Hidden') : 'The switch shows what is hidden, and a line under the switches says so'
        hooks.showSwitches.cells.part.clicked(false)   // a click on the switch itself
        assert picture.getOverlay().size() == drawnAll && hooks.hiddenKinds.isEmpty() && hooks.showSwitches.cells.part.on
        assert !hooks.hiddenLine.text.contains('Hidden') && hooks.hiddenLine.preferredSize.height == 0 : 'With everything shown that line takes no room'
        hooks.legendClicked('background', false); hooks.legendClicked('background', true)
        assert hooks.hiddenKinds == ['cells', 'g0', 'g1', 'g2', 'g3'] as Set && picture.getOverlay().size() <= withoutCells : 'Double-click: only the background'
        hooks.legendClicked('background', false); hooks.legendClicked('background', true)
        assert hooks.hiddenKinds.isEmpty() && picture.getOverlay().size() == drawnAll : 'Double-click again: everything'
        // grades: shape, size and brightness only; never "responds", whatever the rise
        double[] ev = steps.grade(new OvalRoi(4, 22, 40, 16)); assert ev[0] != 0d && ev[5] == 0d && ev[6] == 0d
        assert steps.grade(new OvalRoi(4, 22, 40, 16), false)[5] > 1d : 'With the response the same outline does rise'
        // acceptance check: changing the frames after the stimulus changes nothing about blind detection
        def changed = new ImageStack(w, h); (1..n).each { int f -> def ip = movie.getStack().getProcessor(f).duplicate(); if (f > steps.b1) ip.add(f % 2 == 0 ? 400 : 5); changed.addSlice('' + f, ip) }
        def other = steps.getClass().newInstance(); other.movie = new ImagePlus('changed', changed); other.n = n; other.w = w; other.h = h
        other.startBlind(1, steps.b1)
        assert other.blindDetectionPicture() == steps.blindDetectionPicture() && other.blindOutlinePicture().getProcessor().getPixels() == steps.blindOutlinePicture().getProcessor().getPixels()
        assert other.grade(new OvalRoi(16, 22, 16, 16)) == steps.grade(new OvalRoi(16, 22, 16, 16)) : 'Blind grades do not change with the response frames'
        canvas = picture.getCanvas(); listeners = canvas.getMouseListeners() as List
        hooks.viewAllBtn.doClick()
        assert hooks.addRespondersBtn.visible && !hooks.addRespondersBtn.enabled && hooks.nextLine.text.contains('Continue to measure') && !hooks.nextLine.text.contains('Add cells that respond') && !hooks.autoSuggestAsked()
        assert !steps.blind && hooks.selectionBasis() == 'All-frames structural review' && active(picture) == '110' && hooks.viewAllBtn.getUI().chosen
        assert !hooks.sensBox.visible && hooks.respondsLegend.visible && hooks.recordingNote.text.contains('response near frame')
        assert !hooks.trace.visible && !hooks.trace.caption.contains('hidden') : 'Out of blind choosing nothing says traces are hidden; the frame waits for a trace'
        assert hooks.hoverAt(24, 30).contains('Your cell 1') && hooks.trace.vals != null && hooks.trace.visible : 'The trace takes its room when there is one to show'
        assert (picture.getStack().getPixels(2) as float[]) == (steps.allPic.getProcessor().getPixels() as float[])
        hooks.viewResponseBtn.doClick()
    }
    // magenta after the recording is open: the response is looked for now, and shown on the same window
    waitFor('the response picture') { hooks.viewResponseBtn.enabled && steps.looked }
    onUI {
        assert steps.found && active(picture) == '111' && hooks.selectionBasis() == 'Response-assisted review' && hooks.responseSeen() : 'Showing magenta is recorded as seeing the response'
        float[] magenta = picture.getStack().getPixels(3), response = steps.respPic.getProcessor().getPixels()
        assert magenta == response : 'The magenta channel is the response picture'
        assert picture.getCanvas().is(canvas) && (picture.getCanvas().getMouseListeners() as List).containsAll(listeners) : 'Same window, so clicks on cells still work'
        assert hooks.recordingNote.text.contains('response near frame') && hooks.viewResponseBtn.getUI().chosen && hooks.selectionNote.text.contains('response-assisted')
        assert steps.colourComposite(false).getNChannels() == 3 : 'Once found, the saved picture has the response, as when it was looked for at the start'
        hooks.viewAllBtn.doClick()
        assert active(picture) == '110' && hooks.selectionBasis() == 'All-frames structural review' && hooks.viewAllBtn.getUI().chosen
        float[] green = picture.getStack().getPixels(2), all = steps.allPic.getProcessor().getPixels(); assert green == all
    }
    // back to blind: the frames chosen before are offered again; outlines unchanged since the response was shown stay blind, a new one does not
    SwingUtilities.invokeLater { hooks.viewBaselineBtn.doClick() }
    def prompt = null; long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(20)
    while (prompt == null && System.nanoTime() < end) { Thread.sleep(25); prompt = Window.getWindows().find { it instanceof GenericDialog && it.visible } }
    assert prompt != null
    def range = null; onUI { range = prompt.getNumericFields().collect { it.text as int } }
    int last = steps.proposedBlindWindow()[1]
    assert range == [1, last] : 'The frames before the response are offered'
    onUI { prompt.getNumericFields()[0].text = '2'; prompt.actionPerformed(new java.awt.event.ActionEvent(prompt.getButtons()[0], java.awt.event.ActionEvent.ACTION_PERFORMED, '')) }
    waitFor('the blind picture') { steps.blind && hooks.viewBaselineBtn.enabled }
    onUI {
        assert hooks.selectionBasis() == 'Before the response, blind (frames 2-' + last + ')' && active(picture) == '110' && hooks.viewBaselineBtn.getUI().chosen
        float[] green = picture.getStack().getPixels(2), blindFrames = steps.framesPicture(2, last)
        assert green == blindFrames : 'The picture shows the chosen frames'
        assert (picture.getStack().getPixels(3) as float[]).every { it == 0f } && hooks.trace.vals == null
        assert hooks.basisFor(hooks.cellRois()) == 'Before the response, blind (frames 2-' + last + ')' : 'The outline chosen before the response was shown is still blind'
        hooks.rm.addRoi(new OvalRoi(52, 52, 16, 16))
        assert hooks.basisFor(hooks.cellRois()).startsWith('Picture of frames 2-' + last + ' (before the response), not blind') : 'An outline added after the response was shown is not'
        assert picture.getCanvas().is(canvas)
        // moving between steps: forward and back
        def (s1, s2, s3, s4) = hooks.cards
        hooks.continueBtn.doClick(); assert s3.expanded && !s2.expanded
        hooks.backToCellsBtn.doClick(); assert s2.expanded && !s3.expanded
        s4.onToggle.call(s4); hooks.backToMeasureBtn.doClick(); assert s3.expanded && !s4.expanded
        // back to the first screen, clearing the outlines: the window and the colour picture close, the recording stays
        assert hooks.rm.getCount() == 2
        hooks.changeRecording(0)
        assert hooks.rm.getCount() == 0 && !hooks.frame.isDisplayable() && picture.getWindow() == null
    }
    assert movie.getStack().getProcessor(20).getPixelsCopy() == before && movie.getWindow() != null : 'The recording is only read'
    println 'PASS: Step 2 picture views: blind by default (pre-stimulus frames only in the window and for Cellpose, no traces, no response grade, unchanged by the response frames), magenta turned on after opening (same window and clicks), honest selection basis, all-frames and blind views, saved picture unchanged, back buttons, change recording with outlines cleared'
} catch (Throwable t) { t.printStackTrace(); System.exit(1) }
finally { onUI { hooks?.detectionProgress?.close(); hooks?.frame?.dispose(); movie.close(); imagej?.dispose() } }
System.exit(0)

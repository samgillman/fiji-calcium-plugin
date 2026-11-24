// One recording, its own outlines, driven through the whole plugin on small synthetic recordings. No window is shown: the
// pictures are opened in ImageJ's batch mode, the plugin's window is built and never made visible, and every question is
// answered here instead of in a dialog.
//   replay  The replay of 23 Nov 2026: recording A is measured; the plugin is started again for recording B while A's window
//           is open and the ROI Manager holds A's outlines; then Measure and save is pressed in A's window, first with B's
//           outlines in the ROI Manager and then after A's recording was closed.
//   1 to 5  Each rule on its own, each part with recordings of its own: 1 the recording was closed; 2 outlines already in the
//           ROI Manager; 3 kept outlines in the record and the percent table; 4 one plugin window; 5 the guard at Measure and save.
//   The parts after them each check one case the reviews of 23 Nov found, by the name of that case.
// -Dcalcium.test.plugin=<file> runs the checks on another copy of the plugin (the release, or a copy with one fix taken out,
// to see them fail); -Dcalcium.test.parts=replay,1,... runs only those parts; -Dcalcium.test.work=<folder> writes there.
import ij.*
import ij.gui.*
import ij.io.RoiEncoder
import ij.plugin.frame.RoiManager
import ij.process.*
import javax.swing.*
import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger
import java.util.zip.ZipFile

String pluginPath = System.getProperty('calcium.test.plugin') ?: 'scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'
List<String> ALL_PARTS = ['replay', '1', '2', '3', '4', '5', 'double', 'sameid', 'sametitle', 'handdrawn', 'undoclear', 'step4closed', 'wording',
                          'calibration', 'closemeasure', 'claim', 'reopenmsg', 'cancelids', 'zip', 'emptysecond', 'response', 'workclose', 'noside',
                          'unknownzip', 'claimsame', 'calibchange', 'dupname', 'singular']
List<String> parts = System.getProperty('calcium.test.parts') ? System.getProperty('calcium.test.parts').split(',')*.trim() : ALL_PARTS
String CLOSED = 'The recording is closed in Fiji. Open it again to go on.'
String NAME = 'Ca²⁺ Imaging (Wang Lab)'
String KEEP_ONLY_IF = 'Keep outlines only if this recording shows the same cells in the same place. Kept outlines are recorded as carried over, never as chosen blind.'
String CHANGE_KEEP = 'Keep them only if the next recording shows the same cells in the same place. On another recording, kept outlines are recorded as carried over, never as chosen blind.'
String CLOSED_NO_FILE = 'The recording is closed in Fiji. It did not come from a file, so it cannot be opened again.'
List<String> ONE = ['Clear it', 'Keep it for this recording', 'Cancel'], MANY = ['Clear them', 'Keep them for this recording', 'Cancel']
def onUI = { Closure c -> SwingUtilities.invokeAndWait(c as Runnable) }
def clean = { s -> s == null ? '' : s.toString().replaceAll('<[^>]+>', ' ').replaceAll('&nbsp;', ' ').replaceAll('\\s+', ' ').trim() }
def work = new File(System.getProperty('calcium.test.work') ?: 'build/test-output/recording-guard').absoluteFile; work.deleteDir(); work.mkdirs()

// Nothing may appear on screen. A window that opens all the same is closed at once and fails the test.
def appeared = Collections.synchronizedList([])
java.awt.Toolkit.defaultToolkit.addAWTEventListener({ event ->
    if (event.ID == java.awt.event.WindowEvent.WINDOW_OPENED) { appeared << ('' + event.window); event.window.dispose() }
} as java.awt.event.AWTEventListener, java.awt.AWTEvent.WINDOW_EVENT_MASK)

// 30 or 45 frames of 240 x 200 with round cells that brighten from frame 15 on; "seed" gives each recording its own pixels
def recording = { String title, int frames, List cells, int seed ->
    int w = 240, h = 200; def random = new Random(seed); def stack = new ImageStack(w, h)
    (1..frames).each { int f -> def ip = new ShortProcessor(w, h)
        for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) {
            double v = 100 + random.nextGaussian(); cells.each { c -> if (Math.hypot(x - c[0], y - c[1]) < 9) v += 30 + (f >= 15 ? 25 : 0) }
            ip.set(x, y, (int) Math.round(v)) }
        stack.addSlice('' + f, ip) }
    new ImagePlus(title, stack)
}
// a recording saved to a file and opened from it, so that Fiji knows the file it came from
def fromFile = { File file, int frames, List cells, int seed -> file.parentFile.mkdirs(); IJ.saveAsTiff(recording(file.name, frames, cells, seed), file.path); IJ.openImage(file.path) }
// what Fiji does when a recording's window is closed: the picture leaves the list of open pictures and is emptied
def closeInFiji = { ImagePlus imp -> imp.close(); imp.flush() }
def cellsA = [[70, 70], [120, 110], [170, 80]], cellsB = [[80, 120], [150, 130]]

// No ImageJ window is made at all: without one, ImageJ writes its log and its messages to the console instead of opening windows.
RoiManager manager
ij.macro.Interpreter.batchMode = true
Prefs.set('calcium.steps.check_updates', false)
if (new File('/usr/bin/false').canExecute()) Prefs.set('calcium.steps.cellpose_python', '/usr/bin/false')   // Cellpose is never run here
onUI { manager = new RoiManager(false) }
RoiManager.metaClass.constructor = { boolean hidden -> manager }   // every start of the plugin shares one ROI Manager, as in Fiji

// One start of the plugin. Each start that passes no "script" compiles the plugin in a class loader of its own, as Fiji does.
// answers.reply(title, text) gives the button to press for a question; with none, the question is closed (-1). frames(steps)
// answers the question about the frames before the response (null is Cancel); look and all choose the other two pictures;
// beforeRun(loader) runs once the plugin's classes are loaded and before the plugin starts.
def compile = { -> new GroovyShell(getClass().classLoader, new Binding()).parse(new File(pluginPath)) }
def start = { ImagePlus imp, Map options = [:] ->
    Script script = options.script ? (Script) options.script.getClass().getDeclaredConstructor().newInstance() : compile()
    def loader = script.getClass().classLoader, binding = new Binding([CALCIUM_TEST_GUI: true]); script.setBinding(binding)
    Map answers = options.answers ?: [:]; List asked = []
    def Start = loader.loadClass('CalciumStartDialog')
    Start.metaClass.static.ask = { List labels, String pre, String folder, String fromDisk = null ->
        [recording: labels.find { it.contains(':  ' + imp.title + '   (') }, look: options.look == true, blind: options.look != true && options.all != true, folder: work.path, file: null] }
    Start.metaClass.static.askBlindFrames = { def steps -> options.frames ? options.frames.call(steps) : [1, 10] as int[] }
    Class Outlines = null; try { Outlines = loader.loadClass('CalciumOutlines') } catch (ClassNotFoundException none) { }
    if (Outlines != null) Outlines.metaClass.static.ask = { java.awt.Component parent, String title, String text, List choices ->
        asked << [title: title, text: text, choices: choices, pictureMade: WindowManager.getImage('colour_' + imp.title) != null]
        answers.reply != null ? choices.indexOf(answers.reply.call(title, text)) : -1 }
    if (options.beforeRun) options.beforeRun.call(loader)
    script.run()
    [hooks: binding.hasVariable('hooks') ? binding.getVariable('hooks') : null, alreadyOpen: binding.hasVariable('alreadyOpen') ? binding.getVariable('alreadyOpen') : null,
     alreadyStarting: binding.hasVariable('alreadyStarting'), asked: asked, answers: answers, script: script, loader: loader, Outlines: Outlines, options: options]
}
def draw = { hooks, List cells -> onUI { cells.each { c -> hooks.steps.colourPic.setRoi(new OvalRoi(c[0] - 9, c[1] - 9, 18, 18)); hooks.addBtn.doClick() } } }
def measure = { hooks ->
    onUI { hooks.sensBox.setSelected(false); hooks.go.doClick() }
    long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(120); boolean busy = true
    while (busy && System.nanoTime() < end) { Thread.sleep(50); onUI { busy = !hooks.go.isEnabled() } }
    assert !busy : 'Measure and save did not finish'
    onUI { }   // the result is written by a later event
    String said = null; onUI { said = clean(hooks.result.text) }; said
}
// waits until ok() holds, looked at on the event thread; false when it did not within the time
def waitFor = { int seconds, Closure ok ->
    long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(seconds); boolean done = false
    while (!done && System.nanoTime() < end) { onUI { done = ok.call() }; if (!done) Thread.sleep(50) }
    onUI { }; done
}
// what the plugin writes to ImageJ's log (the console here) while c runs
def capture = { Closure c ->
    PrintStream old = System.out; def seen = new ByteArrayOutputStream()
    def both = new OutputStream() { void write(int b) { seen.write(b); old.write(b) }; void write(byte[] b, int o, int n) { seen.write(b, o, n); old.write(b, o, n) } }
    System.setOut(new PrintStream(both, true, 'UTF-8'))
    try { c.call() } finally { System.out.flush(); System.setOut(old) }
    seen.toString('UTF-8')
}
// the note the plugin keeps of the recording an outline was made on, whichever way this copy of the plugin writes it
def markOf = { Outlines, Steps, ImagePlus imp -> try { Outlines.mark(imp, Steps.printOf(imp)) } catch (MissingMethodException older) { Outlines.mark(imp, imp.stackSize) } }
// A file's content: for a ZIP its entries, since a ZIP also records when it was written.
def contentOf = { File f ->
    if (!f.name.endsWith('.zip')) return f.bytes.encodeHex().toString()
    def zip = new ZipFile(f); try { zip.entries().toList().collect { it.name + '=' + zip.getInputStream(it).bytes.encodeHex() }.join(';') } finally { zip.close() }
}
// what a results folder holds, with when each file was written: a refused action must leave every file as it was
def filesOf = { File dir -> dir.isDirectory() ? dir.listFiles().findAll { it.isFile() }.sort { it.name }.collectEntries { [it.name, it.lastModified() + ':' + contentOf(it).hashCode() + ':' + it.length()] } : [:] }
def changedFiles = { Map before, Map after -> (before.keySet() + after.keySet()).findAll { before[it] != after[it] }.sort() }
def outlinesNow = { -> List out = []; onUI { out = (0..<manager.count).collect { RoiEncoder.saveAsByteArray(manager.getRoi(it)).encodeHex().toString() } }; out }
def closeWindow = { hooks -> onUI { hooks.detectionProgress.close(); hooks.frame.dispose(); hooks.steps.colourPic.close() } }
def text = { File dir, String suffix -> new File(dir, dir.name + suffix).getText('UTF-8') }
def chosenList = { File dir -> text(dir, '_cells_chosen_blind_or_after_response.csv').readLines().drop(1)*.split(',', 2)*.getAt(1) }
// Each part starts as Fiji does after a restart: no plugin window, no picture, an empty ROI Manager, and no note of where
// outlines were made (those notes live as long as the ROI Manager).
def fresh = { ->
    onUI {
        UIManager.put('calcium.outlines', null)
        java.awt.Window.getWindows().each { w -> if (w instanceof JFrame && w.rootPane.getClientProperty('calcium.steps.window') != null) w.dispose() }
        int[] ids = WindowManager.getIDList(); if (ids != null) ids.each { int id -> def im = WindowManager.getImage(id); if (im != null) { im.changes = false; im.close() } }
        manager.reset()
    }
}

try {
    // ------------------------------------------------------------------ the replay of 23 Nov 2026
    if ('replay' in parts) {
        fresh()
        List faults = []
        ImagePlus A = recording('A.tif', 30, cellsA, 1), B = recording('B.tif', 45, cellsB, 2)
        onUI { A.show(); B.show() }
        // A: three cells chosen blind and measured
        def a = start(A)
        assert a.hooks != null && a.asked.isEmpty() : 'an empty ROI Manager: the plugin starts without a question'
        draw(a.hooks, cellsA)
        assert manager.getCount() == 3
        String measuredA = measure(a.hooks)
        assert measuredA.startsWith('Done. 3 cells measured on 30 frames') : measuredA
        assert text(a.hooks.outDir, '_what_was_done.txt').contains('Chosen blind or later: 3 of the 3 cells were chosen blind')
        Map tablesA = filesOf(a.hooks.outDir); assert tablesA.size() > 10
        // the plugin is started again from the menu while A's window is open
        a.answers.reply = { title, said -> 'Cancel' }
        def b = start(B, [answers: [reply: { title, said -> 'Keep them for this recording' }]])
        if (b.hooks != null) faults << 'a second plugin window opened while the first was open'
        else {
            assert b.alreadyOpen == 'A.tif' && a.asked.size() == 1 && a.asked[0].title == 'Change recording' && a.asked[0].choices == ['Clear outlines and go back', 'Keep outlines and go back', 'Cancel'] : 'the open window offers what Change recording offers: ' + a.asked
            assert a.asked[0].text == 'The plugin is already open for A.tif. The ROI Manager holds 3 outlines. ' + CHANGE_KEEP : a.asked[0].text
            assert a.hooks.frame.isDisplayable() && manager.getCount() == 3 : 'Cancel leaves the window and its outlines'
            // What happened to the user needs two windows. A's window is hidden from the next start (its note is taken away), so
            // that B starts beside it as it did on 23 Nov; everything below is then checked rule by rule.
            onUI { a.hooks.frame.rootPane.putClientProperty('calcium.steps.window', null) }
            b = start(B, [answers: [reply: { title, said -> 'Keep them for this recording' }]])
            assert b.hooks != null
        }
        // B starts with A's three outlines in the ROI Manager
        boolean askedAtStart = !b.asked.isEmpty()
        if (!askedAtStart) faults << 'B started with the 3 outlines of A in the ROI Manager and no question was asked'
        else assert b.asked.size() == 1 && b.asked[0].title == 'Outlines in the ROI Manager' && b.asked[0].choices == MANY && !b.asked[0].pictureMade &&
            b.asked[0].text == 'The ROI Manager holds 3 outlines from another recording: A.tif. ' + KEEP_ONLY_IF : b.asked
        draw(b.hooks, cellsB)
        assert manager.getCount() == 5
        String measuredB = measure(b.hooks)
        assert measuredB.startsWith('Done. 5 cells measured on 45 frames') : measuredB
        String recordB = text(b.hooks.outDir, '_what_was_done.txt')
        List chosenB = chosenList(b.hooks.outDir)
        if (recordB.contains('5 of the 5 cells were chosen blind') || chosenB.count { it.startsWith('blind') } == 5)
            faults << 'B\'s record counts the 3 outlines of A as chosen blind on B: "' + recordB.readLines().find { it.startsWith('Chosen blind or later:') }.take(96) + '"'
        else {
            assert chosenB == ['carried over from A.tif'] * 3 + ['blind (frames 1-10)'] * 2 : chosenB
            // with no question at the start, Measure and save asks instead (the answer is the same) and the record says when
            assert recordB.readLines().contains('Outlines carried over: 3 outlines from A.tif were in the ROI Manager ' + (askedAtStart ? 'when this recording was chosen and were kept for it.' : 'when the cells were measured and were kept for this recording.'))
            assert recordB.contains('Chosen blind or later: 2 of the 5 cells were chosen blind on frames 1-10 and are unchanged since; 0 were added or reshaped after the response was shown; 3 were carried over and not chosen on this recording.')
            assert recordB.contains('Carried over: the outlines of 3 of the 5 cells were carried over from A.tif. They were not chosen on this recording and are never counted as chosen blind. Step 4 gives % responding for them separately.')
            assert recordB.contains('Cell selection basis: Before the response, blind (frames 1-10); 3 of the 5 outlines carried over from A.tif.')
            List registry = text(b.hooks.outDir, '_roi_registry.csv').readLines().drop(1)
            assert registry.size() == 5 && registry.take(3).every { it.endsWith(',Carried over from A.tif') } && registry.drop(3).every { it.endsWith(',"Before the response, blind (frames 1-10)"') } : registry
            // what Step 4 is handed: the basis names them, no carried-over cell is chosen blind, and they are a group of their own
            Map handed = b.hooks.analysisMeasurement()
            assert handed.selectionBasis == 'Before the response, blind (frames 1-10); 3 of the 5 outlines carried over from A.tif'
            assert handed.chosenBlind.values().count { it == true } == 2 && handed.chosenBlind.values().count { it == null } == 3 && handed.carriedOver.size() == 3 && handed.carriedOver.values().every { it == 'A.tif' }
            assert handed.carriedOver.keySet().every { handed.chosenBlind.containsKey(it) && handed.chosenBlind[it] == null }
        }
        // A's window: the ROI Manager now also holds the two outlines B's window added. Had A's recording still been open...
        a.answers.reply = { title, said -> 'Cancel' }; int askedBefore = a.asked.size()
        String again = measure(a.hooks)
        if (again.startsWith('Done.')) faults << 'A was measured again with B\'s outlines and its tables were overwritten: "' + again.take(40) + '"'
        else {
            assert again == 'Not measured. The ROI Manager holds 2 outlines from another recording: B.tif. Clear them or keep them for this recording, then measure.' : again
            // A was measured, so its button reads Measure again, and so does the question's title
            assert a.asked.size() == askedBefore + 1 && a.asked[-1].title == 'Measure again' && a.asked[-1].text == 'The ROI Manager holds 2 outlines from another recording: B.tif. The other 3 outlines were made on this recording and stay. ' + KEEP_ONLY_IF : a.asked[-1].title + ': ' + a.asked[-1].text
            assert filesOf(a.hooks.outDir) == tablesA : 'a refused Measure and save writes nothing'
        }
        tablesA = filesOf(a.hooks.outDir)
        // ... and then A's recording is closed in Fiji, and Measure and save is pressed in A's window
        onUI { closeInFiji(A) }
        String afterClose = measure(a.hooks)
        if (afterClose != CLOSED) faults << 'Measure and save on the closed recording said: "' + afterClose + '"'
        else {
            assert !afterClose.contains('java') && !afterClose.contains('Exception')
            assert filesOf(a.hooks.outDir) == tablesA : 'nothing is saved or changed'
            onUI { a.hooks.refresh(false); assert clean(a.hooks.recordingNote.text).startsWith(CLOSED_NO_FILE + ' 30 frames') : 'the Recording card says so too: ' + a.hooks.recordingNote.text }
        }
        assert faults.isEmpty() : 'What happened on 23 Nov 2026 happens with this plugin file:\n  - ' + faults.join('\n  - ')
        println 'PASS: replay of 23 Nov 2026: a second start offers Change recording in the open window; started beside it, B asks about A\'s outlines, records the kept ones as carried over from A and none as chosen blind; A\'s window refuses to measure B\'s outlines; with A\'s recording closed it gives the plain message and writes nothing'
        closeWindow(a.hooks); closeWindow(b.hooks)
    }
    def shared = compile()   // the parts below reuse one compiled copy where a start's own class loader does not matter

    // ------------------------------------------------------------------ 1. the recording was closed
    if ('1' in parts) {
        fresh()
        def saved = new File(work, 'C.tif')
        ImagePlus C = fromFile(saved, 30, cellsA, 3); onUI { C.show() }
        def c = start(C, [script: shared])
        draw(c.hooks, cellsA)
        assert measure(c.hooks).startsWith('Done. 3 cells measured')
        Map tablesC = filesOf(c.hooks.outDir); Map contentC = c.hooks.outDir.listFiles().collectEntries { [it.name, contentOf(it)] }
        // the Step 4 page and the Fiji Step 4 window, opened while the recording is open (a plugin file without the fix has no hook for this)
        def measurementC = c.hooks.analysisMeasurement ? c.hooks.analysisMeasurement() : null, session = null, fijiStep4 = null
        if (measurementC != null) {
            session = c.loader.loadClass('CalciumAnalysisSession').newInstance(measurementC)
            assert session.frame(5, 'raw').length > 100 : 'the Step 4 page shows the recording\'s frames'
            onUI { fijiStep4 = c.loader.loadClass('CalciumAnalysisWindow').newInstance(measurementC); fijiStep4.view.selectedIndex = 2; fijiStep4.showFrame(5)
                assert fijiStep4.cellPicture.note == 'frame 5 as recorded' : 'the Fiji Step 4 window shows the recording\'s frames' }
        }
        onUI { assert !clean(c.hooks.recordingNote.text).contains('closed') && !c.hooks.reopenBtn?.visible; closeInFiji(C) }
        String refused = measure(c.hooks)
        assert refused == CLOSED : 'Measure and save on the closed recording said: "' + refused + '"'
        assert filesOf(c.hooks.outDir) == tablesC && c.hooks.analysisMeasurement().is(measurementC) : 'Measure and save: nothing saved, the measurement kept'
        onUI {
            c.hooks.refresh(false)
            assert clean(c.hooks.recordingNote.text).startsWith(CLOSED + ' 30 frames') && c.hooks.reopenBtn.visible : 'the Recording card says the recording is closed and offers to open its file again: ' + c.hooks.recordingNote.text
        }
        onUI {
            // Suggest cells, the picture of other frames before the response, a cell's trace and the sensitivity check
            c.hooks.suggestBtn.doClick(); assert clean(c.hooks.info.text) == CLOSED && c.hooks.suggestBtn.enabled && c.hooks.proposals.isEmpty()
            c.hooks.viewBaselineBtn.doClick(); assert clean(c.hooks.selectionNote.text) == CLOSED && c.hooks.steps.blind : 'nothing is changed by the refused action'
            c.hooks.viewAllBtn.doClick()   // the all-frames picture is already made: it needs no frame
            assert !c.hooks.steps.blind && c.hooks.hoverAt(70, 70).contains('Your cell 1') && c.hooks.trace.vals == null && c.hooks.trace.caption == CLOSED && c.hooks.trace.visible
            assert c.hooks.runSensitivity(c.hooks.cellRois(), new File(work, 'never_').path, { }).text == 'Sensitivity check not run: the recording is closed in Fiji. Open it again to go on.' && !new File(work, 'never_sensitivity_check.csv').exists()
            // the Fiji Step 4 window: the movie view says so instead of a frame
            fijiStep4.showFrame(6); assert fijiStep4.cellPicture.note == CLOSED && !fijiStep4.player.running
        }
        try { session.frame(5, 'raw'); assert false : 'the Step 4 page must not get a frame of a closed recording' }
        catch (FileNotFoundException refusedFrame) { assert refusedFrame.message == CLOSED }
        // Another recording of the same name and size is not taken for it; the same recording opened again is
        ImagePlus other = recording('C.tif', 30, cellsA, 4); onUI { other.show(); c.hooks.refresh(false) }
        assert measure(c.hooks) == CLOSED : 'a recording with other pixels is another recording'
        onUI { closeInFiji(other); c.hooks.setOpenAgain { String path -> assert path == saved.path; IJ.openImage(path).show() }; c.hooks.reopenBtn.doClick() }
        assert waitFor(20) { c.hooks.reopenBtn.enabled && !c.hooks.reopenBtn.visible } : 'Open the recording again reads the file and the plugin goes on with it'
        onUI { assert !clean(c.hooks.recordingNote.text).contains('closed') && !c.hooks.steps.blind && manager.getCount() == 3 : 'the same cells and the same state' }
        String redone = measure(c.hooks)
        assert redone.startsWith('Done. 3 cells measured on 30 frames') : redone
        assert ['C_gray_values_raw.csv', 'C_gray_values_background_corrected_native.csv', 'C_cells_RoiSet.zip'].every { contentOf(new File(c.hooks.outDir, it)) == contentC[it] } : 'the reopened recording gives the same tables'
        assert text(c.hooks.outDir, '_what_was_done.txt').contains('Recording: it was closed in Fiji and opened again. The plugin went on with the reopened recording after comparing its size, frame count, time base, pixel size and first, middle and last frame.')
        assert session.frame(5, 'raw').length > 100 : 'the Step 4 page goes on with the reopened recording too'
        onUI { fijiStep4.showFrame(5); assert fijiStep4.cellPicture.note == 'frame 5 as recorded'; fijiStep4.dispose() }
        println 'PASS: the recording was closed: Measure and save, Suggest cells, the picture of frames before the response, a trace, the Step 4 page\'s and the Fiji Step 4 window\'s frames say "' + CLOSED + '" and change nothing; the sensitivity check says it was not run; the Recording card says so; another recording of the same name is not taken; Open the recording again goes on with the same cells'
        closeWindow(c.hooks)
    }

    // ------------------------------------------------------------------ 2. outlines already in the ROI Manager
    if ('2' in parts) {
        fresh()
        // another recording: asked once, before the picture is made; Cancel stops the start and leaves everything
        ImagePlus F = recording('F.tif', 30, cellsA, 7), D = recording('D.tif', 30, cellsB, 5); onUI { F.show(); D.show() }
        def f = start(F, [script: shared]); draw(f.hooks, cellsA); closeWindow(f.hooks)
        def cancelled = start(D, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert cancelled.asked.size() == 1 : 'D started with the 3 outlines of F in the ROI Manager and ' + (cancelled.asked ? 'asked ' + cancelled.asked.size() + ' questions' : 'no question was asked')
        assert cancelled.hooks == null && manager.getCount() == 3 && WindowManager.getImage('colour_D.tif') == null : 'Cancel stops the start and leaves the outlines'
        assert cancelled.asked[0].title == 'Outlines in the ROI Manager' && cancelled.asked[0].choices == MANY && !cancelled.asked[0].pictureMade &&
            cancelled.asked[0].text == 'The ROI Manager holds 3 outlines from another recording: F.tif. ' + KEEP_ONLY_IF : cancelled.asked
        // the same recording, started again after its window was closed: its own outlines, no question
        def f2 = start(F, [script: shared])
        assert f2.hooks != null && f2.asked.isEmpty() && manager.getCount() == 3 : 'outlines made on this same recording: no question'
        onUI { assert f2.hooks.carriedIn(f2.hooks.cellRois()).isEmpty() && f2.hooks.basisFor(f2.hooks.cellRois()) == 'Before the response, blind (frames 1-10)' }
        closeWindow(f2.hooks)
        // Clear them removes them
        def d = start(D, [script: shared, answers: [reply: { title, said -> 'Clear them' }]])
        assert d.hooks != null && d.asked.size() == 1 && manager.getCount() == 0 : 'Clear them is the first button and empties the ROI Manager of them'
        assert d.hooks.steps.log.any { it == 'Outlines cleared: 3 outlines from F.tif were in the ROI Manager when this recording was chosen and were cleared before any cell was chosen.' } : d.hooks.steps.log
        // outlines the plugin did not add (from before) next to this recording's own: only the others are asked about and cleared
        draw(d.hooks, cellsB)
        closeWindow(d.hooks)
        onUI { manager.addRoi(new OvalRoi(100, 60, 20, 20)) }
        def d2 = start(D, [script: shared, answers: [reply: { title, said -> 'Clear it' }]])
        assert d2.asked.size() == 1 && d2.asked[0].text == 'The ROI Manager already holds 1 outline, and the plugin cannot tell which recording it is from. The other 2 outlines were made on this recording and stay. ' + KEEP_ONLY_IF : d2.asked
        assert d2.asked[0].choices == ONE && !d2.asked[0].pictureMade
        assert manager.getCount() == 2
        // Change recording keeps its own question, and the next start does not ask again: the outlines are carried over
        d2.answers.reply = { title, said -> 'Keep outlines and go back' }
        onUI { d2.hooks.changeBtn.doClick() }
        assert d2.asked.size() == 2 && d2.asked[1].title == 'Change recording' && d2.asked[1].choices == ['Clear outlines and go back', 'Keep outlines and go back', 'Cancel'] && !d2.hooks.frame.isDisplayable() && manager.getCount() == 2
        assert d2.asked[1].text == 'The ROI Manager holds 2 outlines. They have not been measured yet. ' + CHANGE_KEEP : d2.asked[1].text
        ImagePlus E = recording('E.tif', 30, cellsB, 6); onUI { E.show() }
        def e = start(E, [script: shared])
        assert e.hooks != null && e.asked.isEmpty() : 'kept at Change recording: not asked twice'
        assert e.hooks.steps.log.any { it.startsWith('Outlines carried over: 2 outlines from D.tif were in the ROI Manager when this recording was chosen and were kept for it.') } : e.hooks.steps.log
        // Change recording with the outlines cleared: the next start has nothing to ask about
        e.answers.reply = { title, said -> 'Clear outlines and go back' }
        onUI { e.hooks.changeBtn.doClick() }
        assert manager.getCount() == 0 && !e.hooks.frame.isDisplayable()
        def e2 = start(E, [script: shared]); assert e2.hooks != null && e2.asked.isEmpty()
        closeWindow(e2.hooks)
        println 'PASS: outlines already in the ROI Manager: none or this recording\'s own, no question; another recording\'s or from before, one question before the picture is made, with Clear them first; Cancel stops the start; Change recording keeps its question and the next start does not ask again'
    }

    // What % responding is worked out from: every cell responded, so a row's count is the number of cells in it.
    def Settings = shared.getClass().classLoader.loadClass('CalciumAnalysisSettings'), Csv = shared.getClass().classLoader.loadClass('CalciumCsv')
    def percentRows = { Map handed ->
        File dir = new File(work, 'tally_' + handed.recording); dir.mkdirs()
        new File(dir, 'responder_calls.csv').text = 'Cell,Role,Status,Threshold_increase\n' + handed.names.collect { it + ',Stimulus,Responded,0.1' }.join('\n') + '\n'
        Map found = Settings.tally(dir, handed)
        String sentence = Settings.percentResponding(dir, handed + [folder: dir])
        [rows: found.rows.collectEntries { [it.key, it.label + ':' + it.counted[2]] }, sentence: sentence,
         saved: Csv.read(new File(dir, handed.recording + '_percent_responding_' + dir.name + '.csv')).drop(1).collect { it[0] + ':' + it[1] }]
    }

    // ------------------------------------------------------------------ 3. kept outlines in the record and the percent table
    if ('3' in parts) {
        fresh()
        // blind: 2 outlines of G kept at the start of H, 2 chosen blind on H
        ImagePlus G = recording('G.tif', 30, cellsA, 8), H = recording('H.tif', 30, cellsA + cellsB, 9); onUI { G.show(); H.show() }
        def g = start(G, [script: shared]); draw(g.hooks, cellsB); closeWindow(g.hooks)
        def h = start(H, [script: shared, answers: [reply: { title, said -> 'Keep them for this recording' }]])
        draw(h.hooks, cellsA.take(2))
        assert measure(h.hooks).startsWith('Done. 4 cells measured')
        String recordH = text(h.hooks.outDir, '_what_was_done.txt')
        assert recordH.contains('Chosen blind or later: 2 of the 4 cells were chosen blind on frames 1-10 and are unchanged since; 0 were added or reshaped after the response was shown; 2 were carried over and not chosen on this recording.') :
            'the 2 outlines of G kept for H are counted as chosen on H: "' + recordH.readLines().find { it.startsWith('Chosen blind or later') }?.take(140) + '"'
        assert recordH.readLines().contains('Outlines carried over: 2 outlines from G.tif were in the ROI Manager when this recording was chosen and were kept for it.')
        assert recordH.count('never counted as chosen blind') == 1 : 'the record says more than once that the kept outlines are never counted as chosen blind'
        assert recordH.contains('Carried over: the outlines of 2 of the 4 cells were carried over from G.tif.') && recordH.contains('Cell selection basis: Before the response, blind (frames 1-10); 2 of the 4 outlines carried over from G.tif.')
        assert chosenList(h.hooks.outDir) == ['carried over from G.tif'] * 2 + ['blind (frames 1-10)'] * 2
        Map handedH = h.hooks.analysisMeasurement()
        assert handedH.chosenBlind.values().count { it == true } == 2 && handedH.carriedOver.size() == 2
        Map percentH = percentRows(handedH)
        assert percentH.rows == [all: 'All cells:4', blind: 'Chosen blind:2', later: 'Added after the response was shown:0', carried: 'Carried over:2', quiet: 'Quiet enough to detect a 30 % rise:4'] : 'the percent table: ' + percentH.rows
        assert percentH.saved == ['All cells:4', 'Chosen blind:2', 'Added after the response was shown:0', 'Carried over:2', 'Quiet enough to detect a 30 % rise:4'] : percentH.saved
        assert percentH.sentence == ' All cells: 4 of 4 responded (100.0 %). Chosen blind: 2 of 2 responded (100.0 %). Carried over: 2 of 2 responded (100.0 %).' : percentH.sentence
        // not chosen blind: the same outlines named as carried over, with no blind list at all
        h.answers.reply = { title, said -> 'Keep outlines and go back' }; onUI { h.hooks.changeBtn.doClick() }
        ImagePlus K = recording('K.tif', 30, cellsA + cellsB, 10); onUI { K.show() }
        def k = start(K, [script: shared, look: true]); assert k.asked.isEmpty()
        draw(k.hooks, [cellsA[2]])
        assert measure(k.hooks).startsWith('Done. 5 cells measured')
        String recordK = text(k.hooks.outDir, '_what_was_done.txt')
        assert recordK.contains('Outlines carried over: 4 outlines from other recordings were in the ROI Manager when this recording was chosen and were kept for it.') : recordK
        assert recordK.contains('Carried over: the outlines of 4 of the 5 cells were carried over from other recordings.') && !recordK.contains('Chosen blind or later') && !new File(k.hooks.outDir, 'K_cells_chosen_blind_or_after_response.csv').exists()
        List registryK = text(k.hooks.outDir, '_roi_registry.csv').readLines().drop(1)
        assert registryK.count { it.endsWith(',Carried over from G.tif') } == 2 && registryK.count { it.endsWith(',Carried over from H.tif') } == 2 : registryK
        Map handedK = k.hooks.analysisMeasurement(); assert handedK.chosenBlind == null && handedK.carriedOver.size() == 4
        Map percentK = percentRows(handedK)
        assert percentK.rows == [all: 'All cells:5', carried: 'Carried over:4', quiet: 'Quiet enough to detect a 30 % rise:5'] : percentK.rows
        assert percentK.sentence == ' All cells: 5 of 5 responded (100.0 %). Carried over: 4 of 4 responded (100.0 %).' : percentK.sentence
        // no carried-over outline: the record, the table and the sentence as they always were
        k.answers.reply = { title, said -> 'Clear outlines and go back' }; onUI { k.hooks.changeBtn.doClick() }
        ImagePlus M = recording('M.tif', 30, cellsA, 11); onUI { M.show() }
        def m = start(M, [script: shared]); draw(m.hooks, cellsA)
        assert measure(m.hooks).startsWith('Done. 3 cells measured')
        String recordM = text(m.hooks.outDir, '_what_was_done.txt')
        assert !recordM.contains('arried over') && !recordM.contains('before the plugin started') && recordM.contains('Chosen blind or later: 3 of the 3 cells were chosen blind on frames 1-10 and are unchanged since; 0 were added or reshaped after the response was shown. Step 4')
        Map percentM = percentRows(m.hooks.analysisMeasurement())
        assert !percentM.rows.containsKey('carried') && percentM.sentence == ' All 3 cells were chosen blind: 3 of 3 responded (100.0 %).' : percentM
        closeWindow(m.hooks)
        println 'PASS: kept outlines are named as carried over in what_was_done.txt, the blind list, the registry and what Step 4 is handed, and have a row of their own in the percent table, Carried over, never under Chosen blind; without them the Chosen blind line, the percent table and its sentence are as before'
    }

    // ------------------------------------------------------------------ 4. one plugin window (each start is a compile of its own)
    if ('4' in parts) {
        fresh()
        ImagePlus D = recording('D.tif', 30, cellsB, 5), E = recording('E.tif', 30, cellsB, 6); onUI { D.show(); E.show() }
        def e = start(E); draw(e.hooks, cellsB)
        e.answers.reply = { title, said -> 'Cancel' }
        def second = start(D)
        assert second.hooks == null && second.alreadyOpen == 'E.tif' : 'a second start makes no second window'
        assert e.asked.size() == 1 && e.asked[0].title == 'Change recording' && e.asked[0].choices == ['Clear outlines and go back', 'Keep outlines and go back', 'Cancel'] && e.hooks.frame.isDisplayable() : 'the open window is found from another class loader and offers Change recording: ' + e.asked
        assert e.asked[0].text == 'The plugin is already open for E.tif. The ROI Manager holds 2 outlines. They have not been measured yet. ' + CHANGE_KEEP : e.asked[0].text
        assert WindowManager.getImage('colour_D.tif') == null && manager.getCount() == 2 : 'the second start made no picture and left the outlines'
        e.answers.reply = { title, said -> 'Clear outlines and go back' }
        second = start(D)
        assert second.hooks == null && !e.hooks.frame.isDisplayable() && manager.getCount() == 0 : 'Clear outlines and go back closes the open window'
        def third = start(D)
        assert third.hooks != null && third.asked.isEmpty() : 'with no window open the plugin starts'
        closeWindow(third.hooks)
        println 'PASS: one plugin window: a start while a window is open makes no second window, offers Change recording in the open one, and starts normally once that window is gone'
    }

    // ------------------------------------------------------------------ 5. the guard at Measure and save
    if ('5' in parts) {
        fresh()
        ImagePlus D = recording('D.tif', 30, cellsB, 5), E = recording('E.tif', 30, cellsB, 6); onUI { D.show(); E.show() }
        def d = start(D, [script: shared])
        draw(d.hooks, cellsB)
        def Outlines = d.Outlines, Steps = d.loader.loadClass('CalciumSteps')
        assert Outlines != null : 'this plugin file does not note on which recording it added an outline, so Measure and save cannot tell another recording\'s outlines from this one\'s'
        onUI { manager.addRoi(new OvalRoi(100, 60, 20, 20)); d.hooks.refresh(true)
            Outlines.made(manager, manager.getRoi(2), markOf(Outlines, Steps, E)) }   // an outline the plugin added on E, in the list of D's window
        Map tablesD = filesOf(d.hooks.outDir)
        d.answers.reply = { title, said -> 'Cancel' }
        String refused = measure(d.hooks)
        assert refused == 'Not measured. The ROI Manager holds 1 outline from another recording: E.tif. Clear it or keep it for this recording, then measure.' : 'Measure and save with an outline of E in D\'s window said: "' + refused.take(60) + '"'
        assert d.asked.size() == 1 && d.asked[0].title == 'Measure and save' && d.asked[0].choices == ONE && filesOf(d.hooks.outDir) == tablesD
        d.answers.reply = { title, said -> 'Clear it' }
        assert measure(d.hooks) == '1 outline from E.tif was removed. Nothing was measured. Check the cells, then press Measure and save.' && manager.getCount() == 2 && filesOf(d.hooks.outDir) == tablesD
        onUI { d.hooks.undoLast(); assert manager.getCount() == 3 : 'the removal can be undone' }
        d.answers.reply = { title, said -> 'Keep it for this recording' }
        String kept = measure(d.hooks)
        assert kept.startsWith('Done. 3 cells measured') : kept
        String recordD = text(d.hooks.outDir, '_what_was_done.txt')
        assert recordD.readLines().contains('Outlines carried over: 1 outline from E.tif was in the ROI Manager when the cells were measured and was kept for this recording.') : recordD
        assert recordD.readLines().contains('Carried over: the outline of 1 of the 3 cells was carried over from E.tif. It was not chosen on this recording and is never counted as chosen blind. Step 4 gives % responding for it separately.') : recordD
        assert recordD.contains('; 1 was carried over and not chosen on this recording.')
        int questions = d.asked.size(); assert measure(d.hooks).startsWith('Done. 3 cells measured') && d.asked.size() == questions : 'kept once, not asked again'
        closeWindow(d.hooks)
        println 'PASS: Measure and save refuses outlines the plugin added on another recording until they are cleared or kept for this recording, and writes nothing when it refuses'
    }

    // ------------------------------------------------------------------ double: a start while the first start makes its pictures
    if ('double' in parts) {
        fresh()
        ImagePlus D = recording('DS.tif', 30, cellsB, 181), E = recording('ES.tif', 30, cellsB, 182); onUI { D.show(); E.show() }
        // the first start waits at its question about the frames, after its pictures are made and before its window exists
        def waiting = new CountDownLatch(1), goOn = new CountDownLatch(1); def first = null, failed = null
        Thread starting = Thread.start {
            try { first = start(E, [frames: { s -> waiting.countDown(); goOn.await(60, TimeUnit.SECONDS); [1, 10] as int[] }]) } catch (Throwable problem) { failed = problem }
        }
        assert waiting.await(60, TimeUnit.SECONDS) : 'the first start did not reach the question about the frames'
        def second = start(D)
        goOn.countDown(); starting.join(120000)
        assert failed == null : failed
        int windows = 0; onUI { windows = java.awt.Window.getWindows().count { it instanceof JFrame && it.isDisplayable() && it.rootPane.getClientProperty('calcium.steps.window') != null } }
        assert second.hooks == null && windows == 1 : 'a start while the first start made its pictures opened a second plugin window: ' + windows + ' plugin windows'
        assert second.alreadyStarting && WindowManager.getImage('colour_DS.tif') == null && first.hooks != null
        assert first.Outlines.STARTING_TEXT == 'The plugin is already starting. Finish or cancel that start first.'
        // once the first start has its window, a start finds that window
        first.answers.reply = { title, said -> 'Cancel' }
        def third = start(D)
        assert third.hooks == null && third.alreadyOpen == 'ES.tif' && first.asked.size() == 1
        // a start that ends without a window, here by Cancel at the frames, lets the next start begin
        closeWindow(first.hooks)
        def cancelled = start(D, [frames: { s -> null }])
        assert cancelled.hooks == null && !cancelled.alreadyStarting
        def fourth = start(D); assert fourth.hooks != null && !fourth.alreadyStarting
        closeWindow(fourth.hooks)
        println 'PASS: double: a start while the first start makes its pictures opens nothing; once the first has its window, a start finds it; a start that ends lets the next one begin'
    }

    // ------------------------------------------------------------------ sameid: the same outline added on two recordings
    if ('sameid' in parts) {
        fresh()
        ImagePlus SD = recording('SD.tif', 30, cellsB, 171), SE = recording('SE.tif', 30, cellsB, 172); onUI { SD.show(); SE.show() }
        def sd = start(SD, [script: shared]); draw(sd.hooks, cellsB)
        def Steps = sd.loader.loadClass('CalciumSteps')
        // an outline the plugin added on SE in a window beside this one, then the same outline, point for point, added here on SD
        onUI { manager.addRoi(new OvalRoi(100, 60, 20, 20)); sd.Outlines.made(manager, manager.getRoi(2), markOf(sd.Outlines, Steps, SE)) }
        onUI { sd.hooks.steps.colourPic.setRoi(new OvalRoi(100, 60, 20, 20)); sd.hooks.addBtn.doClick() }
        assert manager.getCount() == 4
        sd.answers.reply = { title, said -> 'Cancel' }
        String said = measure(sd.hooks)
        assert said == 'Not measured. The ROI Manager holds 1 outline from another recording: SE.tif. Clear it or keep it for this recording, then measure.' :
            'with the same outline added on SD, the outline of SE ' + (said.startsWith('Done.') ? 'was measured as SD\'s own: "' + said.take(40) + '"' : 'gave: "' + said + '"')
        sd.answers.reply = { title, said2 -> 'Keep it for this recording' }
        assert measure(sd.hooks).startsWith('Done. 4 cells measured')
        assert text(sd.hooks.outDir, '_what_was_done.txt').contains('Chosen blind or later: 3 of the 4 cells were chosen blind on frames 1-10 and are unchanged since; 0 were added or reshaped after the response was shown; 1 was carried over and not chosen on this recording.')
        assert chosenList(sd.hooks.outDir) == ['blind (frames 1-10)'] * 2 + ['carried over from SE.tif', 'blind (frames 1-10)']
        closeWindow(sd.hooks)
        println 'PASS: sameid: an outline added on one recording keeps its note when the same outline is added on another; it is asked about and recorded as carried over'
    }

    // ------------------------------------------------------------------ sametitle: the same name and size, not from a file
    if ('sametitle' in parts) {
        fresh()
        ImagePlus first = recording('Run1.tif', 30, cellsA, 161); onUI { first.show() }
        ImagePlus firstAgain = first.duplicate(); firstAgain.setTitle('Run1.tif')   // the same pixels, as a picture made in Fiji
        def one = start(first, [script: shared]); draw(one.hooks, cellsA); closeWindow(one.hooks)
        onUI { closeInFiji(first) }
        ImagePlus second = recording('Run1.tif', 30, cellsA, 162); onUI { second.show() }   // another recording of the same name and size
        def two = start(second, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert two.asked.size() == 1 : 'another recording named Run1.tif, of the same size and not from a file, took the 3 outlines of the first as its own'
        assert two.asked[0].text == 'The ROI Manager holds 3 outlines from another recording: Run1.tif. ' + KEEP_ONLY_IF : two.asked[0].text
        onUI { closeInFiji(second); firstAgain.show() }
        def three = start(firstAgain, [script: shared])
        assert three.hooks != null && three.asked.isEmpty() : 'the same pixels under the same name: the same recording, no question'
        closeWindow(three.hooks)
        println 'PASS: sametitle: a recording of the same name and size with other pixels is another recording, also when neither came from a file'
    }

    // ------------------------------------------------------------------ handdrawn: outlines drawn on this recording before the start
    if ('handdrawn' in parts) {
        fresh()
        ImagePlus HD = recording('HD.tif', 30, cellsA + cellsB, 41); onUI { HD.show() }
        onUI { cellsB.each { c -> HD.setRoi(new OvalRoi(c[0] - 9, c[1] - 9, 18, 18)); manager.addRoi(HD.getRoi()) }; HD.deleteRoi() }   // drawn on the recording's own window
        def h = start(HD, [script: shared, answers: [reply: { title, said -> 'Clear them' }]])
        assert h.hooks != null && h.asked.isEmpty() && manager.getCount() == 2 : 'outlines drawn on this recording\'s own window before the start were asked about: ' + h.asked*.text
        draw(h.hooks, [cellsA[0]])
        assert measure(h.hooks).startsWith('Done. 3 cells measured')
        String record = text(h.hooks.outDir, '_what_was_done.txt')
        assert record.readLines().contains('Outlines drawn before the plugin started: 2 outlines drawn on this recording were in the ROI Manager when the plugin started and were kept.') : record
        assert record.contains('Chosen blind or later: 1 of the 3 cells were chosen blind on frames 1-10 and are unchanged since; 0 were added or reshaped after the response was shown; 2 were drawn before the plugin started.')
        assert record.contains('Drawn before the plugin started: the outlines of 2 of the 3 cells were drawn on this recording before the plugin started. They were not chosen in the plugin and are never counted as chosen blind. Step 4 counts them under Carried over.')
        assert !record.contains('arried over from') && !record.readLines().any { it.startsWith('Carried over:') || it.startsWith('Outlines carried over:') }
        assert record.contains('Cell selection basis: Before the response, blind (frames 1-10); 2 of the 3 outlines drawn before the plugin started.')
        assert chosenList(h.hooks.outDir) == ['drawn on this recording before the plugin started'] * 2 + ['blind (frames 1-10)']
        List registry = text(h.hooks.outDir, '_roi_registry.csv').readLines().drop(1)
        assert registry.take(2).every { it.endsWith(',Drawn on this recording before the plugin started') } && registry[2].endsWith(',"Before the response, blind (frames 1-10)"') : registry
        Map handed = h.hooks.analysisMeasurement()
        assert handed.chosenBlind.values().count { it == true } == 1 && handed.carriedOver.size() == 2
        assert percentRows(handed).rows.subMap(['all', 'blind', 'later', 'carried']) == [all: 'All cells:3', blind: 'Chosen blind:1', later: 'Added after the response was shown:0', carried: 'Carried over:2']
        // on another recording they are another recording's outlines, named after the one they were drawn on
        closeWindow(h.hooks)
        ImagePlus HE = recording('HE.tif', 30, cellsA, 42); onUI { HE.show() }
        def other = start(HE, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert other.asked.size() == 1 && other.asked[0].text == 'The ROI Manager holds 3 outlines from another recording: HD.tif. ' + KEEP_ONLY_IF : other.asked*.text
        // one outline drawn before the start: every line about it in the singular
        onUI { manager.reset() }
        ImagePlus HF = recording('HF.tif', 30, cellsA + cellsB, 43); onUI { HF.show(); HF.setRoi(new OvalRoi(cellsB[0][0] - 9, cellsB[0][1] - 9, 18, 18)); manager.addRoi(HF.getRoi()); HF.deleteRoi() }
        def one = start(HF, [script: shared]); assert one.hooks != null && one.asked.isEmpty()
        draw(one.hooks, [cellsA[0]])
        assert measure(one.hooks).startsWith('Done. 2 cells measured')
        List oneRecord = text(one.hooks.outDir, '_what_was_done.txt').readLines()
        assert oneRecord.contains('Outlines drawn before the plugin started: 1 outline drawn on this recording was in the ROI Manager when the plugin started and was kept.') : oneRecord
        assert oneRecord.contains('Drawn before the plugin started: the outline of 1 of the 2 cells was drawn on this recording before the plugin started. It was not chosen in the plugin and is never counted as chosen blind. Step 4 counts it under Carried over.') : oneRecord
        assert oneRecord.any { it.startsWith('Chosen blind or later: 1 of the 2 cells were chosen blind on frames 1-10 and are unchanged since; 0 were added or reshaped after the response was shown; 1 was drawn before the plugin started.') } : oneRecord
        closeWindow(one.hooks)
        println 'PASS: handdrawn: outlines drawn on the recording\'s own window before the start are kept without a question and recorded as drawn before the plugin started, never as chosen blind and never as from another recording'
    }

    // ------------------------------------------------------------------ undoclear: Clear them at the start
    if ('undoclear' in parts) {
        fresh()
        ImagePlus F = recording('UF.tif', 30, cellsA, 191), D = recording('UD.tif', 30, cellsB, 192); onUI { F.show(); D.show() }
        def f = start(F, [script: shared]); draw(f.hooks, cellsA); closeWindow(f.hooks)
        List before = outlinesNow()
        // Clear them, then Cancel at the frames: the start does not go on, and the outlines stay as they were
        def stopped = start(D, [script: shared, answers: [reply: { title, said -> 'Clear them' }], frames: { s -> null }])
        assert stopped.hooks == null && stopped.asked.size() == 1
        assert manager.getCount() == 3 : 'Clear them took effect although the start was cancelled at the frames: the ROI Manager holds ' + manager.getCount()
        assert outlinesNow() == before
        // Clear them, and the start goes on: the outlines go as one change to the cells, which Undo brings back
        def d = start(D, [script: shared, answers: [reply: { title, said -> 'Clear them' }]])
        assert d.hooks != null && d.asked.size() == 1 && manager.getCount() == 0
        String hint = null; onUI { hint = clean(d.hooks.editHint?.text); d.hooks.undoLast() }
        assert manager.getCount() == 3 : 'Undo does not bring back the outlines cleared at the start'
        assert hint == 'Cleared 3 outlines from UF.tif when this recording was chosen. Undo brings them back.' : hint
        // brought back, they are still another recording's outlines
        d.answers.reply = { title, said -> 'Cancel' }
        assert measure(d.hooks) == 'Not measured. The ROI Manager holds 3 outlines from another recording: UF.tif. Clear them or keep them for this recording, then measure.'
        closeWindow(d.hooks)
        println 'PASS: undoclear: Clear them takes effect only when the start goes on, as one change that Undo brings back; outlines brought back are asked about again at Measure and save'
    }

    // ------------------------------------------------------------------ step4closed: Step 4 opened while closed, and after opening again
    if ('step4closed' in parts) {
        fresh()
        File saved = new File(work, 'S4/S4.tif')
        ImagePlus S = fromFile(saved, 30, cellsA, 51); onUI { S.show() }
        def s = start(S, [script: shared]); draw(s.hooks, cellsA)
        assert measure(s.hooks).startsWith('Done. 3 cells measured')
        Map handed = s.hooks.analysisMeasurement()
        onUI { closeInFiji(S) }
        def Session = s.loader.loadClass('CalciumAnalysisSession'), Window = s.loader.loadClass('CalciumAnalysisWindow')
        // Step 4 opened while the recording is closed: its movie views are there and say that the recording is closed
        def page = Session.newInstance(handed)
        String pageSays = null; try { page.frame(5, 'raw') } catch (FileNotFoundException refused) { pageSays = refused.message }
        assert pageSays == CLOSED : 'the Step 4 page opened while the recording is closed says: ' + pageSays
        assert page.described().hasMovie
        def fijiStep4 = null
        onUI { fijiStep4 = Window.newInstance(handed); assert fijiStep4.view.itemCount == 3 : 'the Fiji Step 4 window opened while the recording is closed has no movie views'
            fijiStep4.view.selectedIndex = 2; fijiStep4.showFrame(5); assert fijiStep4.cellPicture.note == CLOSED }
        // Open the recording again: both go on with it, and so does a Step 4 opened afterwards without measuring again
        onUI { s.hooks.refresh(false); assert s.hooks.reopenBtn.visible : 'the Recording card offers to open the closed recording again'
            s.hooks.setOpenAgain { String path -> IJ.openImage(path).show() }; s.hooks.reopenBtn.doClick() }
        assert waitFor(20) { s.hooks.reopenBtn.enabled && !s.hooks.reopenBtn.visible } : 'Open the recording again reads the file and the plugin goes on with it'
        assert page.frame(5, 'raw').length > 100 : 'the Step 4 page opened while the recording was closed does not take up the recording opened again'
        onUI { fijiStep4.showFrame(5); assert fijiStep4.cellPicture.note == 'frame 5 as recorded'; fijiStep4.dispose() }
        Map again = s.hooks.analysisMeasurement()
        assert again.is(handed) && !again.movie.is(S) && again.movie.isProcessor() : 'what Step 4 is handed still holds the closed recording'
        def later = Session.newInstance(again)
        assert later.described().hasMovie && later.frame(5, 'raw').length > 100 : 'a Step 4 page opened after Open the recording again has no movie'
        onUI { def w = Window.newInstance(again); assert w.view.itemCount == 3; w.view.selectedIndex = 2; w.showFrame(5); assert w.cellPicture.note == 'frame 5 as recorded'; w.dispose() }
        closeWindow(s.hooks)
        println 'PASS: step4closed: the Step 4 page and the Fiji Step 4 window opened while the recording is closed say so in their movie views, and they, and a Step 4 opened after Open the recording again, go on with the recording opened again'
    }

    // ------------------------------------------------------------------ wording: the question at Measure and save
    if ('wording' in parts) {
        fresh()
        ImagePlus WG = recording('WG.tif', 30, cellsA, 61), WD = recording('WD.tif', 30, cellsA + cellsB, 62), WE = recording('WE.tif', 30, cellsB, 63); onUI { [WG, WD, WE].each { it.show() } }
        def g = start(WG, [script: shared]); draw(g.hooks, cellsA.take(2))
        g.answers.reply = { title, said -> 'Keep outlines and go back' }; onUI { g.hooks.changeBtn.doClick() }
        def d = start(WD, [script: shared]); assert d.asked.isEmpty() && manager.getCount() == 2
        draw(d.hooks, [cellsB[0]])
        def Steps = d.loader.loadClass('CalciumSteps')
        onUI { manager.addRoi(new OvalRoi(100, 60, 20, 20)); d.Outlines.made(manager, manager.getRoi(3), markOf(d.Outlines, Steps, WE)); d.hooks.refresh(true) }
        d.answers.reply = { title, said -> 'Cancel' }
        assert measure(d.hooks) == 'Not measured. The ROI Manager holds 1 outline from another recording: WE.tif. Clear it or keep it for this recording, then measure.'
        assert d.asked[-1].text == 'The ROI Manager holds 1 outline from another recording: WE.tif. The other 3 outlines stay. ' + KEEP_ONLY_IF :
            'the question at Measure and save says: "' + d.asked[-1].text + '"; 2 of the other 3 were carried over from WG.tif'
        closeWindow(d.hooks)
        println 'PASS: wording: the question at Measure and save says the other outlines were made on this recording only when they were'
    }

    // ------------------------------------------------------------------ calibration: the same pixels with another time base
    if ('calibration' in parts) {
        fresh()
        ImagePlus X = recording('X.tif', 30, cellsA, 71)
        def seconds = new ij.measure.Calibration(); seconds.frameInterval = 2; seconds.setTimeUnit('sec'); X.setCalibration(seconds); onUI { X.show() }
        ImagePlus sameTime = X.duplicate(), otherTime = X.duplicate()
        sameTime.setTitle('X again.tif'); sameTime.setCalibration(seconds.copy())
        def none = seconds.copy(); none.frameInterval = 0; otherTime.setTitle('X other.tif'); otherTime.setCalibration(none)
        def x = start(X, [script: shared]); draw(x.hooks, cellsA)
        assert measure(x.hooks).startsWith('Done. 3 cells measured') && x.hooks.analysisMeasurement().hz == 0.5d
        Map tables = filesOf(x.hooks.outDir)
        onUI { closeInFiji(X); otherTime.show(); x.hooks.refresh(false) }
        String refused = measure(x.hooks)
        assert refused == CLOSED : 'a copy with the same pixels and no frame interval was taken for the recording: "' + refused.take(40) + '"'
        assert filesOf(x.hooks.outDir) == tables
        onUI { closeInFiji(otherTime); sameTime.show(); x.hooks.refresh(false) }
        assert measure(x.hooks).startsWith('Done. 3 cells measured') && x.hooks.analysisMeasurement().hz == 0.5d : 'a copy with the same pixels and time base is the same recording'
        closeWindow(x.hooks)
        println 'PASS: calibration: a copy with the same pixels and another time base is another recording; with the same time base it is taken'
    }

    // ------------------------------------------------------------------ closemeasure: the recording closed while it is measured
    if ('closemeasure' in parts) {
        fresh()
        ImagePlus CM = fromFile(new File(work, 'CM/CM.tif'), 30, cellsA + cellsB, 81); onUI { CM.show() }
        def cm = start(CM, [script: shared]); draw(cm.hooks, cellsA)
        onUI { cm.hooks.makeMovie.selected = true; cm.hooks.saveMovie.selected = true }
        assert measure(cm.hooks).startsWith('Done. 3 cells measured')
        draw(cm.hooks, cellsB)
        Map tables = filesOf(cm.hooks.outDir); assert tables.containsKey('CM_background_corrected_movie.tif') && tables.containsKey('CM_cells_chosen_blind_or_after_response.csv')
        // Fiji closes the recording just as its frames are read
        def realMeasure = cm.loader.loadClass('CalciumSteps').getMethod('measureAndSave', Roi[], String[], String)
        cm.hooks.steps.metaClass.measureAndSave = { Roi[] cells, String[] names, String prefix -> closeInFiji(CM); realMeasure.invoke(delegate, cells, names, prefix) }
        String said = null
        String logged = capture { said = measure(cm.hooks) }
        assert said == CLOSED : 'Measure and save, the recording closed while it was measured: "' + said + '"'
        assert filesOf(cm.hooks.outDir) == tables : 'the recording closed while it was measured, and these files were changed: ' + changedFiles(tables, filesOf(cm.hooks.outDir))
        assert logged.contains(NAME + ': ' + CLOSED) && !logged.contains('Exception') : 'the log says: ' + logged.readLines().findAll { it.startsWith(NAME) }
        closeWindow(cm.hooks)
        println 'PASS: closemeasure: the recording closed while it is measured: "' + CLOSED + '" on the window and in the log, and no file changed'
    }

    // ------------------------------------------------------------------ claim: an outline added with the ROI Manager's own Add
    if ('claim' in parts) {
        fresh()
        ImagePlus A = recording('CA.tif', 30, cellsA, 91), B = recording('CB.tif', 30, cellsB, 92); onUI { A.show(); B.show() }
        def a = start(A, [script: shared])
        onUI { def picture = a.hooks.steps.colourPic; picture.setRoi(new OvalRoi(cellsA[0][0] - 9, cellsA[0][1] - 9, 18, 18)); manager.addRoi(picture.getRoi()); picture.deleteRoi(); a.hooks.refresh(true) }
        closeWindow(a.hooks)
        def b = start(B, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert b.asked.size() == 1 && b.asked[0].text == 'The ROI Manager holds 1 outline from another recording: CA.tif. ' + KEEP_ONLY_IF :
            'an outline added with the ROI Manager\'s own Add on the colour picture of CA is not known as CA\'s: ' + b.asked*.text
        println 'PASS: claim: an outline drawn on the colour picture and added with the ROI Manager\'s own Add is noted as made on that recording'
    }

    // ------------------------------------------------------------------ reopenmsg: what the Recording card says about the file
    if ('reopenmsg' in parts) {
        fresh()
        File file = new File(work, 'reopen/RO.tif'), away = new File(work, 'reopen/RO moved.tif')
        ImagePlus R = fromFile(file, 30, cellsA, 101); onUI { R.show() }
        def r = start(R, [script: shared]); draw(r.hooks, cellsA)
        // the file is gone when the recording is closed: no button, and the card says so
        assert file.renameTo(away)
        onUI { closeInFiji(R); r.hooks.refresh(false) }
        onUI { assert clean(r.hooks.recordingNote.text).startsWith(CLOSED + ' The file is no longer there: ' + file.path) && !r.hooks.reopenBtn.visible : 'the file is gone: ' + clean(r.hooks.recordingNote.text) }
        // it is back: so is the button
        assert away.renameTo(file)
        onUI { r.hooks.refresh(false); assert r.hooks.reopenBtn.visible && clean(r.hooks.recordingNote.text).startsWith(CLOSED + ' 30 frames') : clean(r.hooks.recordingNote.text) }
        // the file does not open
        String logged = capture {
            onUI { r.hooks.setOpenAgain { String path -> throw new IOException('the drive is not connected') }; r.hooks.reopenBtn.doClick() }
            assert waitFor(20) { r.hooks.reopenBtn.enabled }
        }
        onUI { assert clean(r.hooks.recordingNote.text).startsWith(CLOSED + ' The file could not be opened. 30 frames') : 'nothing opened: ' + clean(r.hooks.recordingNote.text) }
        assert logged.contains(NAME + ': the recording could not be opened again: the drive is not connected') && !logged.contains('java.io') : logged
        // it opens with other frames
        onUI { r.hooks.setOpenAgain { String path -> recording('RO.tif', 30, cellsA, 102).show() }; r.hooks.reopenBtn.doClick() }
        assert waitFor(20) { r.hooks.reopenBtn.enabled }
        onUI { assert clean(r.hooks.recordingNote.text).startsWith(CLOSED + ' The file did not open as the same recording. 30 frames') : 'other frames: ' + clean(r.hooks.recordingNote.text) }
        closeWindow(r.hooks)
        // a recording that did not come from a file
        ImagePlus N = recording('RN.tif', 30, cellsA, 103); onUI { N.show(); manager.reset() }
        def n = start(N, [script: shared])
        onUI { closeInFiji(N); n.hooks.refresh(false) }
        onUI { assert clean(n.hooks.recordingNote.text).startsWith(CLOSED_NO_FILE + ' 30 frames') && !clean(n.hooks.recordingNote.text).contains('Open it again') && !n.hooks.reopenBtn.visible : clean(n.hooks.recordingNote.text) }
        closeWindow(n.hooks)
        println 'PASS: reopenmsg: the Recording card says when the file is no longer there, could not be opened, or did not open as the same recording, and when the recording did not come from a file, without saying to open it again; Open the recording again is offered only while the file is there'
    }

    // ------------------------------------------------------------------ cancelids: Cancel at the start changes nothing
    if ('cancelids' in parts) {
        fresh()
        ImagePlus D = recording('CI.tif', 30, cellsB, 111); onUI { D.show() }
        onUI { cellsA.each { c -> manager.addRoi(new OvalRoi(c[0] - 9, c[1] - 9, 18, 18)) } }   // outlines from before, with no ID
        List before = outlinesNow()
        def d = start(D, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert d.hooks == null && d.asked.size() == 1
        List ids = []; onUI { ids = (0..<manager.count).collect { manager.getRoi(it).getProperty('calcium.roi.id') } }
        assert ids == [null] * 3 : 'Cancel at the start wrote IDs into the outlines: ' + ids
        assert outlinesNow() == before
        println 'PASS: cancelids: Cancel at the start leaves the outlines byte for byte as they were'
    }

    // ------------------------------------------------------------------ zip: a RoiSet of another recording opened by hand
    if ('zip' in parts) {
        fresh()
        ImagePlus ZA = recording('ZA.tif', 30, cellsA, 121), ZB = recording('ZB.tif', 30, cellsA + cellsB, 122); onUI { ZA.show(); ZB.show() }
        def za = start(ZA, [script: shared]); draw(za.hooks, cellsA)
        assert measure(za.hooks).startsWith('Done. 3 cells measured')
        File zip = new File(za.hooks.outDir, 'ZA_accepted_rois.zip'); assert zip.isFile()
        za.answers.reply = { title, said -> 'Clear outlines and go back' }; onUI { za.hooks.changeBtn.doClick() }
        def zb = start(ZB, [script: shared]); assert zb.hooks != null && zb.asked.isEmpty()
        onUI { manager.open(zip.path) }   // opened by hand in the ROI Manager while cells are chosen blind
        assert manager.getCount() == 3
        draw(zb.hooks, cellsB)
        zb.answers.reply = { title, said -> 'Cancel' }
        String said = measure(zb.hooks)
        assert said == 'Not measured. The ROI Manager holds 3 outlines from another recording: ZA.tif. Clear them or keep them for this recording, then measure.' :
            'the outlines of ZA opened from its RoiSet in the window of ZB: "' + said.take(40) + '"'
        closeWindow(zb.hooks)
        println 'PASS: zip: a RoiSet of another recording measured while the ROI Manager stayed open, opened by hand, is asked about at Measure and save'
    }

    // ------------------------------------------------------------------ emptysecond: a second start with the ROI Manager empty
    if ('emptysecond' in parts) {
        fresh()
        ImagePlus EA = recording('EA.tif', 30, cellsA, 131), EB = recording('EB.tif', 30, cellsB, 132); onUI { EA.show(); EB.show() }
        def ea = start(EA); draw(ea.hooks, cellsA)
        assert measure(ea.hooks).startsWith('Done. 3 cells measured')
        onUI { manager.reset() }   // emptied by hand
        ea.answers.reply = { title, said -> 'Cancel' }
        def second = start(EB)
        assert second.hooks == null && second.alreadyOpen == 'EA.tif'
        assert ea.asked.size() == 1 && ea.hooks.frame.isDisplayable() : 'a start with the ROI Manager empty closed the open window without a question'
        assert ea.asked[0].title == 'Change recording' && ea.asked[0].text == 'The plugin is already open for EA.tif. Go back to choose another recording?' && ea.asked[0].choices == ['Go back', 'Cancel'] : ea.asked
        ea.answers.reply = { title, said -> 'Go back' }
        start(EB)
        assert !ea.hooks.frame.isDisplayable()
        println 'PASS: emptysecond: a start while the window is open with the ROI Manager empty asks before going back, and Cancel keeps the window'
    }

    // ------------------------------------------------------------------ response: the response view while the recording is closed
    if ('response' in parts) {
        fresh()
        ImagePlus RS = recording('RS.tif', 30, cellsA, 141), RT = recording('RT.tif', 30, cellsA, 142); onUI { RS.show(); RT.show() }
        // blind: the response was looked for at the start, so showing it needs no frame
        def rs = start(RS, [script: shared]); draw(rs.hooks, cellsA)
        onUI { closeInFiji(RS); rs.hooks.refresh(false); rs.hooks.viewResponseBtn.doClick() }
        waitFor(20) { !rs.hooks.steps.blind || clean(rs.hooks.selectionNote.text) == CLOSED }
        onUI { assert clean(rs.hooks.selectionNote.text).startsWith('Magenta: cells that brightened near frame') && rs.hooks.responseSeen() : 'Response with the response already found: ' + clean(rs.hooks.selectionNote.text) }
        closeWindow(rs.hooks); onUI { manager.reset() }
        // all frames: the response has not been looked for, and looking for it reads the frames
        def rt = start(RT, [script: shared, all: true])
        onUI { closeInFiji(RT); rt.hooks.refresh(false); rt.hooks.viewResponseBtn.doClick(); assert clean(rt.hooks.selectionNote.text) == CLOSED && !rt.hooks.steps.looked && !rt.hooks.responseSeen() }
        closeWindow(rt.hooks)
        println 'PASS: response: Response shows the response already found while the recording is closed, and refuses to look for it then'
    }

    // ------------------------------------------------------------------ workclose: the recording closed while the work runs
    if ('workclose' in parts) {
        fresh()
        // the sensitivity check
        ImagePlus WS = recording('WS.tif', 30, cellsA, 151); onUI { WS.show() }
        def ws = start(WS, [script: shared]); draw(ws.hooks, cellsA)
        assert ws.hooks.steps.found
        def realPicture = ws.loader.loadClass('CalciumSteps').getMethod('responseOutlinePicture')
        ws.hooks.steps.metaClass.responseOutlinePicture = { -> def picture = realPicture.invoke(delegate); closeInFiji(WS); picture }
        Roi[] cells = null; onUI { cells = ws.hooks.cellRois() }
        def checked = ws.hooks.runSensitivity(cells, new File(work, 'ws_').path, { })
        // it reads only the pictures made at the start, so it runs on, and what it says is its own result (Cellpose fails here)
        // (Cellpose is /usr/bin/false where that exists, so it makes no outlines; elsewhere, as on Windows, it is not found)
        assert checked.text.startsWith('Sensitivity check not completed: Cellpose ') && !checked.text.contains(CLOSED) : 'the recording closed during the sensitivity check: "' + checked.text + '"'
        closeWindow(ws.hooks); onUI { manager.reset() }
        // the picture of other frames before the response
        ImagePlus WB = recording('WB.tif', 30, cellsA, 152); onUI { WB.show() }
        Map optionsWB = [script: shared]; def wb = start(WB, optionsWB)
        optionsWB.frames = { s -> closeInFiji(WB); [1, 10] as int[] }   // closed in Fiji just after the frames are confirmed
        onUI { wb.hooks.viewBaselineBtn.doClick() }
        waitFor(20) { !clean(wb.hooks.selectionNote.text).startsWith('Making the picture') }
        onUI { assert clean(wb.hooks.selectionNote.text) == CLOSED : 'the recording closed while the picture of frames 1-10 was made: "' + clean(wb.hooks.selectionNote.text) + '"' }
        closeWindow(wb.hooks); onUI { manager.reset() }
        // Suggest cells
        ImagePlus WC = recording('WC.tif', 30, cellsA, 153); onUI { WC.show() }
        def wc = start(WC, [beforeRun: { loader ->
            def Progress = loader.loadClass('CalciumOutlineProgress'); def realBegin = Progress.getMethod('begin')
            Progress.metaClass.begin = { -> def task = realBegin.invoke(delegate); closeInFiji(WC); task } }])   // closed as the suggestions begin
        onUI { wc.hooks.suggestBtn.doClick() }
        assert waitFor(60) { wc.hooks.suggestBtn.enabled }
        onUI { assert clean(wc.hooks.info.text).startsWith('Suggestions: Cellpose ') && !clean(wc.hooks.info.text).contains(CLOSED) : 'the recording closed while the suggestions were made: "' + clean(wc.hooks.info.text) + '"' }
        closeWindow(wc.hooks)
        println 'PASS: workclose: the recording closed while the picture of frames before the response is made: it says "' + CLOSED + '"; closed while Suggest cells or the sensitivity check runs, which read only the pictures made at the start: each reports its own result'
    }

    // ------------------------------------------------------------------ noside: looking whether the recording is closed touches nothing
    if ('noside' in parts) {
        fresh()
        def Steps = shared.getClass().classLoader.loadClass('CalciumSteps')
        AtomicInteger touched = new AtomicInteger()
        ImagePlus watched = new ImagePlus('NS.tif', recording('NS.tif', 30, cellsA, 201).getStack()) { ImageProcessor getProcessor() { touched.incrementAndGet(); super.getProcessor() } }
        onUI { watched.show() }
        assert !Steps.closed(watched) && touched.get() == 0 : 'looking whether the recording is closed calls getProcessor on it'
        def ns = start(watched, [script: shared]); draw(ns.hooks, cellsA)
        touched.set(0); onUI { 5.times { ns.hooks.refresh(true) } }
        assert touched.get() == 0 : 'the window, looking every 600 ms, called getProcessor on the recording ' + touched.get() + ' times in 5 looks'
        onUI { closeInFiji(watched) }; touched.set(0)
        assert Steps.closed(watched) && touched.get() == 0
        closeWindow(ns.hooks)
        println 'PASS: noside: looking whether the recording is closed reads nothing from it and changes nothing on it'
    }

    // ------------------------------------------------------------------ unknownzip: outlines of no known recording cleared at a start
    if ('unknownzip' in parts) {
        fresh()
        // UZ was measured in an earlier Fiji session: its RoiSet is all that is left of it
        ImagePlus UZ0 = recording('UZ.tif', 30, cellsA, 372); onUI { UZ0.show() }
        def first = start(UZ0, [script: shared]); draw(first.hooks, cellsA)
        assert measure(first.hooks).startsWith('Done. 3 cells measured')
        File zip = new File(work, 'UZ_accepted_rois.zip'); zip.bytes = new File(first.hooks.outDir, 'UZ_accepted_rois.zip').bytes
        closeWindow(first.hooks)
        fresh()   // Fiji closed and started again: no notes
        ImagePlus UY = recording('UY.tif', 30, cellsB, 371), UZ = recording('UZ.tif', 30, cellsA, 372); onUI { UY.show(); UZ.show(); manager.open(zip.path) }   // opened to look at it
        def y = start(UY, [script: shared, answers: [reply: { title, said -> 'Clear them' }]])
        assert y.asked.size() == 1 && y.asked[0].text.startsWith('The ROI Manager already holds 3 outlines, and the plugin cannot tell which recording they are from.') && manager.getCount() == 0 : y.asked*.text
        // in the window that cleared them, Undo brings them back and they are asked about again
        onUI { y.hooks.undoLast() }; assert manager.getCount() == 3
        y.answers.reply = { title, said -> 'Cancel' }
        assert measure(y.hooks) == 'Not measured. The ROI Manager already holds 3 outlines, and the plugin cannot tell which recording they are from. Clear them or keep them for this recording, then measure.'
        y.answers.reply = { title, said -> 'Clear outlines and go back' }; onUI { y.hooks.changeBtn.doClick() }
        assert manager.getCount() == 0
        // UZ, started with the ROI Manager empty; its own RoiSet is opened once its window is open, as the guide says
        def z = start(UZ, [script: shared]); assert z.hooks != null && z.asked.isEmpty()
        onUI { manager.open(zip.path); z.hooks.refresh(true) }
        z.answers.reply = { title, said -> 'Cancel' }
        String said = measure(z.hooks)
        assert said.startsWith('Done. 3 cells measured on 30 frames') : 'UZ with its own RoiSet, cleared before at the start of UY: "' + said + '"'
        assert z.asked.isEmpty() && chosenList(z.hooks.outDir) == ['blind (frames 1-10)'] * 3 : chosenList(z.hooks.outDir)
        closeWindow(z.hooks)
        // the same with an outline from before that has no ID, and the window that cleared it closed instead
        ImagePlus UB = recording('UB.tif', 30, cellsB, 301), UC = recording('UC.tif', 30, cellsA, 302)
        onUI { manager.reset(); UB.show(); UC.show(); manager.addRoi(new OvalRoi(cellsA[0][0] - 9, cellsA[0][1] - 9, 18, 18)) }
        def b = start(UB, [script: shared, answers: [reply: { title, said2 -> 'Clear it' }]])
        String hint = null; onUI { hint = clean(b.hooks.editHint?.text) }
        assert b.asked.size() == 1 && b.asked[0].choices == ONE && manager.getCount() == 0 : b.asked
        assert hint == 'Cleared 1 outline from an unknown recording when this recording was chosen. Undo brings it back.' : hint
        closeWindow(b.hooks)
        def c = start(UC, [script: shared]); assert c.hooks != null && c.asked.isEmpty()
        onUI { cellsA.each { p -> manager.addRoi(new OvalRoi(p[0] - 9, p[1] - 9, 18, 18)) }; c.hooks.refresh(true) }
        c.answers.reply = { title, said2 -> 'Cancel' }
        String saidC = measure(c.hooks)
        assert saidC.startsWith('Done. 3 cells measured on 30 frames') && c.asked.isEmpty() && chosenList(c.hooks.outDir) == ['blind (frames 1-10)'] * 3 : saidC
        closeWindow(c.hooks)
        println 'PASS: unknownzip: outlines of no known recording cleared at one start are asked about again only in the window that cleared them; opened later in their own recording\'s window they are its own and chosen blind, as before'
    }

    // ------------------------------------------------------------------ claimsame: the ROI Manager's Add of a shape added on another recording
    if ('claimsame' in parts) {
        fresh()
        ImagePlus QA = recording('QA.tif', 30, cellsA, 311), QB = recording('QB.tif', 30, cellsA + cellsB, 312); onUI { QA.show(); QB.show() }
        def a = start(QA, [script: shared]); draw(a.hooks, cellsA)
        assert measure(a.hooks).startsWith('Done. 3 cells measured')
        a.answers.reply = { title, said -> 'Clear outlines and go back' }; onUI { a.hooks.changeBtn.doClick() }
        def b = start(QB, [script: shared]); assert b.hooks != null && b.asked.isEmpty() && manager.getCount() == 0
        // Edit > Selection > Restore Selection puts QA's first oval on QB's colour picture, and the ROI Manager's own Add takes it
        onUI { def picture = b.hooks.steps.colourPic; picture.setRoi(new OvalRoi(cellsA[0][0] - 9, cellsA[0][1] - 9, 18, 18)); manager.addRoi(picture.getRoi()); picture.deleteRoi(); b.hooks.refresh(true) }
        draw(b.hooks, cellsB)
        b.answers.reply = { title, said -> 'Cancel' }
        String said = measure(b.hooks)
        assert said.startsWith('Done. 3 cells measured on 30 frames') : 'an outline drawn on QB in the place of one of QA\'s and added with the ROI Manager\'s Add: "' + said + '"'
        assert b.asked.isEmpty() && chosenList(b.hooks.outDir) == ['blind (frames 1-10)'] * 3 : chosenList(b.hooks.outDir)
        // an outline already in the list keeps its note: one noted on QA after it was seen here stays QA's
        onUI { b.Outlines.made(manager, manager.getRoi(0), markOf(b.Outlines, b.loader.loadClass('CalciumSteps'), QA)); b.hooks.refresh(true) }
        assert measure(b.hooks).startsWith('Not measured. The ROI Manager holds 1 outline from another recording: QA.tif.')
        closeWindow(b.hooks)
        println 'PASS: claimsame: an outline drawn on this recording\'s colour picture and added with the ROI Manager\'s Add is its own, also when an outline of the same shape was added on another recording before; an outline already in the list keeps its note'
    }

    // ------------------------------------------------------------------ calibchange: the frame interval or pixel size set after the outlines
    if ('calibchange' in parts) {
        fresh()
        ImagePlus K = recording('KX.tif', 30, cellsA, 321); onUI { K.show() }
        def k = start(K, [script: shared]); draw(k.hooks, cellsA)
        assert measure(k.hooks).startsWith('Done. 3 cells measured'); closeWindow(k.hooks)
        onUI { def cal = K.getCalibration(); cal.frameInterval = 2; cal.setTimeUnit('sec') }   // Image > Properties
        def again = start(K, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert again.hooks != null && again.asked.isEmpty() : 'the frame interval set in Image > Properties made the recording another one for its own outlines: ' + again.asked*.text
        assert measure(again.hooks).startsWith('Done. 3 cells measured') && again.hooks.analysisMeasurement().hz == 0.5d
        assert chosenList(again.hooks.outDir) == ['blind (frames 1-10)'] * 3 : chosenList(again.hooks.outDir)
        // Change recording with the outlines kept, then the same recording with another pixel size: its own, not carried over
        again.answers.reply = { title, said -> 'Keep outlines and go back' }; onUI { again.hooks.changeBtn.doClick() }
        onUI { K.getCalibration().pixelWidth = 0.5d; K.getCalibration().pixelHeight = 0.5d }
        def third = start(K, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert third.hooks != null && third.asked.isEmpty()
        onUI { assert third.hooks.carriedIn(third.hooks.cellRois()).isEmpty() }
        assert !third.hooks.steps.log.any { it.startsWith('Outlines carried over') } : third.hooks.steps.log
        closeWindow(third.hooks)
        println 'PASS: calibchange: a recording whose frame interval or pixel size was set after its outlines were chosen keeps them as its own, with no question and nothing carried over'
    }

    // ------------------------------------------------------------------ dupname: the same outline added twice with the plugin's Add
    if ('dupname' in parts) {
        fresh()
        ImagePlus DN = recording('DN.tif', 30, cellsA, 331); onUI { DN.show() }
        def d = start(DN, [script: shared]); draw(d.hooks, cellsA); draw(d.hooks, [cellsA[0]])
        List names = [], ids = []
        onUI { d.hooks.refresh(true); names = (0..<manager.count).collect { manager.getName(it) }; ids = (0..<manager.count).collect { manager.getRoi(it).getProperty('calcium.roi.id') } }
        assert names.size() == 4 && names[3] == names[0] : 'the same outline added twice is labelled in the ROI Manager: ' + names
        assert ids.toSet().size() == 4
        closeWindow(d.hooks)
        println 'PASS: dupname: an outline added twice gets an ID of its own and the ROI Manager label it always had'
    }

    // ------------------------------------------------------------------ singular: one outline, said as one
    if ('singular' in parts) {
        fresh()
        ImagePlus NA = recording('NA.tif', 30, cellsB, 341), NB = recording('NB.tif', 30, cellsA + cellsB, 342); onUI { NA.show(); NB.show() }
        def a = start(NA, [script: shared]); draw(a.hooks, cellsB); closeWindow(a.hooks)
        onUI { manager.addRoi(new OvalRoi(100, 60, 20, 20)) }   // and one outline from before
        def mixed = start(NB, [script: shared, answers: [reply: { title, said -> 'Cancel' }]])
        assert mixed.asked.size() == 1 && mixed.asked[0].text == 'The ROI Manager holds 3 outlines from NA.tif or from before. ' + KEEP_ONLY_IF && mixed.asked[0].choices == MANY : mixed.asked*.text
        // one outline from before, cleared at the start
        onUI { manager.reset(); manager.addRoi(new OvalRoi(100, 60, 20, 20)) }
        def b = start(NB, [script: shared, answers: [reply: { title, said -> 'Clear it' }]])
        assert b.hooks != null && b.asked.size() == 1 && b.asked[0].choices == ONE && b.asked[0].text == 'The ROI Manager already holds 1 outline, and the plugin cannot tell which recording it is from. ' + KEEP_ONLY_IF : b.asked
        // one cell chosen here, and one outline the plugin added on NA in this window's list
        draw(b.hooks, [cellsA[0]])
        def Steps = b.loader.loadClass('CalciumSteps')
        def addFromNA = { -> onUI { manager.addRoi(new OvalRoi(120, 60, 18, 18)); b.Outlines.made(manager, manager.getRoi(manager.count - 1), markOf(b.Outlines, Steps, NA)); b.hooks.refresh(true) } }
        addFromNA()
        b.answers.reply = { title, said -> 'Cancel' }
        assert measure(b.hooks) == 'Not measured. The ROI Manager holds 1 outline from another recording: NA.tif. Clear it or keep it for this recording, then measure.'
        assert b.asked[-1].title == 'Measure and save' && b.asked[-1].choices == ONE && b.asked[-1].text == 'The ROI Manager holds 1 outline from another recording: NA.tif. The other outline was made on this recording and stays. ' + KEEP_ONLY_IF : b.asked[-1].text
        b.answers.reply = { title, said -> 'Clear it' }
        assert measure(b.hooks) == '1 outline from NA.tif was removed. Nothing was measured. Check the cells, then press Measure and save.'
        // Change recording with one outline, cancelled
        b.answers.reply = { title, said -> 'Cancel' }; onUI { b.hooks.changeBtn.doClick() }
        assert b.asked[-1].title == 'Change recording' && b.asked[-1].text == 'The ROI Manager holds 1 outline. It has not been measured yet. Keep it only if the next recording shows the same cells in the same place. On another recording, kept outlines are recorded as carried over, never as chosen blind.' : b.asked[-1].text
        // once measured the button reads Measure again: so do the question and the message
        assert measure(b.hooks).startsWith('Done.')
        addFromNA()
        b.answers.reply = { title, said -> 'Clear it' }
        assert measure(b.hooks) == '1 outline from NA.tif was removed. Nothing was measured. Check the cells, then press Measure again.'
        assert b.asked[-1].title == 'Measure again'
        closeWindow(b.hooks)
        println 'PASS: singular: one outline is said as one in the questions, their buttons, the messages and the record; outlines from a recording and from before are named so; Clear at Measure again names the button the window shows'
    }
    assert appeared.isEmpty() : 'windows appeared on screen: ' + appeared
} catch (Throwable failure) {
    failure.printStackTrace()
    if (!appeared.isEmpty()) System.err.println('windows appeared on screen: ' + appeared)
    System.exit(1)
}
System.exit(0)

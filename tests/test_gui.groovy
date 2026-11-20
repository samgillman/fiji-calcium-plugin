// Optional desktop smoke test: run tools/run_tests.py --gui with a display available.
import javax.swing.*
import java.util.concurrent.TimeUnit

def loader = new GroovyClassLoader(getClass().classLoader)
def source = new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy').getText('UTF-8')
loader.parseClass(source)
def Window = loader.loadClass('CalciumWindow')
def Panel = loader.loadClass('CalciumOutlineProgress')
def Release = loader.loadClass('CalciumRelease')
def Installer = loader.loadClass('CalciumInstaller')
def testDirectory = java.nio.file.Files.createTempDirectory('Fiji GUI install ').toFile()
def target = new File(testDirectory,'scripts/Plugins/Calcium Imaging/'+Release.SCRIPT)
target.parentFile.mkdirs(); target.setText(source,'UTF-8')
byte[] updateBytes = source.replace(Release.VERSION,'2099.1.1.1').getBytes('UTF-8')
def next = [version:'2099.1.1.1',page:Release.REPOSITORY+'/releases/tag/v2099.1.1.1',installable:true,
    script:Release.REPOSITORY+'/releases/download/v2099.1.1.1/'+Release.SCRIPT,
    checksums:Release.REPOSITORY+'/releases/download/v2099.1.1.1/'+Release.CHECKSUMS,size:updateBytes.length]
def buttons = { String label -> new JButton(label) }
def binding = new Binding([
    btn:buttons, note:{ String s -> new JLabel('<html>'+s+'</html>') },
    box:{ String s, boolean selected -> new JCheckBox(s,selected) },
    row:{ java.awt.Component... children -> def p = new JPanel(); children.each{p.add(it)}; p },
])
def hooks, frame, progress
def Updates = loader.loadClass('CalciumUpdatesPanel')
// Step 1: the styled dialog returns the recording, the response choice, the folder and any file chosen from disk, or nothing when cancelled.
def Start = loader.loadClass('CalciumStartDialog')
Start.updatesFactory = { -> def panel = Updates.newInstance(); panel.updateCheck.fetch = { -> throw new IOException('offline fixture') }; panel }
def everything; everything = { c -> [c] + (c instanceof java.awt.Container ? c.components.collectMany { everything(it) } : []) }
def movie = '1:  movie.tif   (300 frames, 16-bit)', disk = 'Choose a file from disk...'
String offered = null; int asked = 0
Start.chooseFile = { -> asked++; offered }                       // stands in for the native file chooser
Start.chooseFolder = { String from -> 'browsed folder' }
def answerStep1 = { String preselected, Closure... acts ->
    def answer = 'not answered'
    def asking = Thread.start { answer = Start.ask([movie, disk], preselected, 'default folder', disk) }
    def dialog = null
    long limit = System.nanoTime() + TimeUnit.SECONDS.toNanos(5)
    while (dialog == null && System.nanoTime() < limit) { Thread.sleep(20); dialog = java.awt.Window.getWindows().find { it instanceof JDialog && it.visible && it.modal } }
    assert dialog != null : 'Step 1 dialog did not open'
    acts.each { act -> SwingUtilities.invokeAndWait { act(everything(dialog)) }; SwingUtilities.invokeAndWait { } } // the second call lets queued choosers run
    asking.join(5000); assert !asking.alive && !dialog.displayable
    answer
}
def button = { parts, String text -> parts.find { it instanceof JButton && it.text == text } }
// A recording that is already open: no file row, no file in the answer.
assert answerStep1(disk, { parts ->
    def choice = parts.find { it instanceof JComboBox }; assert choice.selectedItem == disk
    assert button(parts, 'Choose file…').isShowing() && parts.any { it instanceof JLabel && it.text == 'No file chosen yet' }
    choice.selectedIndex = 0
    assert !button(parts, 'Choose file…').isShowing()
    parts.find { it instanceof JTextField }.text = 'chosen folder'
    def radio = { String start -> parts.find { it instanceof JRadioButton && it.text.startsWith(start) } }
    assert radio('Before the response only').selected && !radio('All frames, with the response').selected && !radio('All frames').selected : 'Blind selection is the default; response-assisted selection starts off'
    radio('All frames, with the response').doClick()
    assert !radio('Before the response only').selected : 'One way of choosing cells at a time'
    assert parts.findAll { it instanceof JButton }*.text.containsAll(['Browse…', 'Cancel', 'Next'])
    button(parts, 'Next').doClick()
}) == [recording: movie, look: true, blind: false, folder: 'chosen folder', file: null]
assert asked == 0
assert answerStep1(disk, { parts -> button(parts, 'Cancel').doClick() }) == null
// Picking "from disk" in the list asks for the file at once and shows its name; a cancelled chooser keeps the dialog open.
def chosen = new File(testDirectory, 'recording 1.zvi').path
assert answerStep1(movie, { parts ->
    assert !button(parts, 'Choose file…').isShowing()
    offered = null; parts.find { it instanceof JComboBox }.selectedItem = disk
}, { parts ->
    assert asked == 1 && button(parts, 'Choose file…').isShowing() : 'Choosing "from disk" must open the file chooser straight away'
    button(parts, 'Next').doClick()                              // still no file: asks again, stays open when cancelled
    assert asked == 2 && button(parts, 'Next').isShowing()
    offered = chosen; button(parts, 'Choose file…').doClick()
    assert asked == 3 && parts.any { it instanceof JLabel && it.text == 'recording 1.zvi' && it.toolTipText == chosen }
    assert button(parts, 'Change file…') != null
    button(parts, 'Browse…').doClick(); assert parts.find { it instanceof JTextField }.text == 'browsed folder'
    button(parts, 'Next').doClick()
}) == [recording: disk, look: false, blind: true, folder: 'browsed folder', file: chosen]
assert asked == 3
// Next alone is enough when "from disk" was already selected: it asks for the file and continues.
assert answerStep1(disk, { parts -> button(parts, 'Next').doClick() }) == [recording: disk, look: false, blind: true, folder: 'default folder', file: chosen]
assert asked == 4
try {
    SwingUtilities.invokeAndWait {
        frame = Window.newInstance('Calcium UI regression test')
        def root = new JPanel(); root.setLayout(new BoxLayout(root,BoxLayout.Y_AXIS))
        def panel = Updates.newInstance(); root.add(panel); frame.contentPane.add(root)
        hooks = [check:{ automatic -> panel.checkUpdates(automatic) }, controller:panel.updateCheck,
            button:panel.updateBtn, note:panel.updateNote, checkbox:panel.autoUpdateBox, install:panel.installBtn,
            later:panel.laterBtn, busy:panel.installBusy, bar:panel.installBar, start:{ data -> panel.startInstall(data) },
            setInstaller:{ Closure c -> panel.installRelease = c }, panel:panel]
        frame.pack()
        assert frame.height <= frame.graphicsConfiguration.bounds.height - 60
        hooks.controller.enabled = true
        hooks.controller.fetch = { -> next }
        hooks.check(true)
        assert !hooks.button.enabled
    }
    long end = System.nanoTime() + TimeUnit.SECONDS.toNanos(5)
    while (hooks.controller.busy.get() && System.nanoTime() < end) Thread.sleep(10)
    SwingUtilities.invokeAndWait {
        assert !frame.ownedWindows.any { it.visible } : 'Updates stay inline on Step 1'
        assert hooks.button.enabled && hooks.note.text.contains('2099.1.1.1')
        assert hooks.install.visible && hooks.later.visible
        hooks.later.doClick()
        assert !hooks.install.visible && !hooks.later.visible
        hooks.install.visible=true // retain the found release while exercising failure/retry below
        progress = Panel.newInstance(); frame.contentPane.add(progress,java.awt.BorderLayout.SOUTH)
        progress.begin().report(40,'Finding cells','Testing the details view')
        progress.detailsButton.doClick()
        assert progress.detailsWindow.visible && !progress.detailsWindow.modal
        assert progress.detailsText.text.contains('Testing the details view')
        progress.close(); assert !progress.detailsWindow.displayable
        hooks.controller.fetch = { -> throw new IOException('offline') }
        hooks.check(true)
    }
    end = System.nanoTime() + TimeUnit.SECONDS.toNanos(5)
    while (hooks.controller.busy.get() && System.nanoTime() < end) Thread.sleep(10)
    SwingUtilities.invokeAndWait {
        assert !frame.ownedWindows.any{it.visible}
        assert hooks.button.enabled
        hooks.setInstaller { data, installProgress ->
            assert !SwingUtilities.isEventDispatchThread()
            throw new IOException('simulated interrupted download')
        }
        hooks.install.doClick()
        assert hooks.busy.get() && !hooks.install.enabled
        hooks.start(next) // duplicate click cannot start another worker
    }
    end = System.nanoTime() + TimeUnit.SECONDS.toNanos(5)
    while (hooks.busy.get() && System.nanoTime() < end) Thread.sleep(10)
    SwingUtilities.invokeAndWait {
        assert !hooks.busy.get() && hooks.install.enabled && hooks.note.text.contains('simulated interrupted download')
        assert target.getText('UTF-8') == source
        hooks.setInstaller { data, installProgress ->
            assert !SwingUtilities.isEventDispatchThread()
            def service = Installer.newInstance(applicationDirectory:testDirectory,download:{url,limit,onBytes ->
                url==next.checksums ? (Installer.digest(updateBytes)+'  '+Release.SCRIPT+'\n').getBytes('UTF-8') : updateBytes
            })
            service.install(data,installProgress)
        }
        hooks.install.doClick()
    }
    end = System.nanoTime() + TimeUnit.SECONDS.toNanos(10)
    while (hooks.busy.get() && System.nanoTime() < end) Thread.sleep(10)
    SwingUtilities.invokeAndWait {
        assert !hooks.busy.get() && !hooks.install.enabled && hooks.bar.value==100
        assert hooks.note.text.contains('restart Fiji') && target.bytes==updateBytes
        assert new File(testDirectory,'calcium-plugin-backups').listFiles().size()==1
        assert !hooks.button.enabled
    }
    println 'PASS: Step 1 dialog answer, cancel and in-dialog file choice, inline Install update button, retry after error, background verified replacement, backup, restart feedback, Later and offline startup'
} finally {
    SwingUtilities.invokeAndWait { progress?.close(); hooks?.controller?.close(); frame?.dispose() }
    testDirectory.deleteDir()
}
System.exit(0) // Native desktop event loops may otherwise keep a successful smoke test alive.

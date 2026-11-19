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
try {
    SwingUtilities.invokeAndWait {
        frame = Window.newInstance('Calcium UI regression test')
        def root = new JPanel(); root.setLayout(new BoxLayout(root,BoxLayout.Y_AXIS))
        binding.setVariable('frame',frame); binding.setVariable('root',root)
        def fragment = source.substring(source.indexOf('def updateBtn ='),source.indexOf('frame.pack(); frame.setAlwaysOnTop'))
        hooks = new GroovyShell(loader,binding).evaluate('import ij.IJ\nimport ij.Prefs\n'+fragment+'\nreturn [check:checkUpdates, controller:updateCheck, button:updateBtn, note:updateNote, checkbox:autoUpdateBox, install:installBtn, busy:installBusy, bar:installBar, start:startInstall, setInstaller:{ Closure c -> installRelease = c }]')
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
        def prompt = frame.ownedWindows.find{it instanceof JDialog && it.visible}
        assert prompt != null && !prompt.modal && hooks.button.enabled
        assert hooks.note.text.contains('2099.1.1.1')
        def descendants; descendants = { c -> [c] + (c instanceof java.awt.Container ? c.components.collectMany { descendants(it) } : []) }
        def controls = descendants(prompt).findAll { it instanceof JButton }
        assert controls*.text.containsAll(['Install update','Later'])
        controls.find{it.text == 'Later'}.doClick()
        assert !prompt.displayable
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
    println 'PASS: actual Install update button, retry after error, background verified replacement, backup, restart feedback, Later and offline startup'
} finally {
    SwingUtilities.invokeAndWait { progress?.close(); hooks?.controller?.close(); frame?.dispose() }
    testDirectory.deleteDir()
}
System.exit(0) // Native desktop event loops may otherwise keep a successful smoke test alive.

// Exercise updates on the actual modal recording chooser, before any image is opened.
import javax.swing.*
import java.awt.event.WindowEvent
import java.util.concurrent.*

def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Start = loader.loadClass('CalciumStartDialog'), Updates = loader.loadClass('CalciumUpdatesPanel'), Release = loader.loadClass('CalciumRelease')
def onUI = { Closure c -> SwingUtilities.invokeAndWait(c as Runnable) }
def await = { Closure condition ->
    boolean done=false; long end=System.nanoTime()+TimeUnit.SECONDS.toNanos(8)
    while (!done && System.nanoTime()<end) { onUI { done=condition() }; if (!done) Thread.sleep(20) }
    assert done : 'Step 1 did not reach the expected state'
}
def all; all = { c -> [c] + (c instanceof java.awt.Container ? c.components.collectMany { all(it) } : []) }
def control = { dialog, text -> all(dialog).find { it instanceof AbstractButton && it.text==text } }
def launched=[]; def panels=[]
def open = { Closure configure ->
    def state=[:]
    Start.updatesFactory = { ->
        def panel=Updates.newInstance(); configure(panel); panels<<panel; state.panel=panel; panel
    }
    state.thread=Thread.start { state.answer=Start.ask(['Choose a file from disk...'],'Choose a file from disk...','chosen folder') }
    await { state.dialog=java.awt.Window.getWindows().find { it instanceof JDialog && it.visible && it.modal };state.dialog!=null }
    launched<<state
    state
}
boolean oldPreference=ij.Prefs.get('calcium.steps.check_updates',true)
try {
    UIManager.setLookAndFeel(UIManager.systemLookAndFeelClassName)
    ij.Prefs.set('calcium.steps.check_updates',true)
    def slow=new CountDownLatch(1), called=new CountDownLatch(1)
    def pending=open { panel -> panel.updateCheck.fetch={ -> called.countDown();slow.await(5,TimeUnit.SECONDS);[version:'2099.1.1.1',page:Release.REPOSITORY+'/releases',installable:true] } }
    assert called.await(5,TimeUnit.SECONDS)
    onUI {
        assert pending.panel.updateBtn.text=='Checking...'
        assert control(pending.dialog,'Next').enabled && control(pending.dialog,'Cancel').enabled
        assert control(pending.dialog,'Release notes').isShowing()
        control(pending.dialog,'Cancel').doClick()
        assert pending.panel.closed && pending.panel.updateCheck.closed
    }
    slow.countDown();pending.thread.join(5000);assert !pending.thread.alive
    await { !pending.panel.updateCheck.busy.get() }
    onUI { assert !pending.panel.installBtn.visible : 'Late update response must not reopen closed Step 1' }

    // Disabled startup performs no fetch; manual checks still work before any file is chosen.
    ij.Prefs.set('calcium.steps.check_updates',false)
    int fetched=0;def installGate=new CountDownLatch(1)
    def ready=open { panel ->
        panel.updateCheck.fetch={ -> fetched++; [version:'2099.1.1.1',page:Release.REPOSITORY+'/releases',installable:true] }
        panel.installRelease={ found,progress ->
            assert !SwingUtilities.isEventDispatchThread()
            progress(40,'Verifying download');installGate.await(5,TimeUnit.SECONDS)
            [version:found.version,target:'test plugin.groovy',backup:'backup.groovy']
        }
    }
    onUI {
        assert fetched==0 && !ready.panel.autoUpdateBox.selected
        assert control(ready.dialog,'Check for updates').isShowing()
        control(ready.dialog,'Check for updates').doClick()
    }
    await { ready.panel.installBtn.visible }
    onUI {
        assert ready.panel.laterBtn.visible && fetched==1
        ready.panel.laterBtn.doClick();assert !ready.panel.installBtn.visible && control(ready.dialog,'Next').enabled
        ready.panel.updateBtn.doClick()
    }
    await { ready.panel.installBtn.visible }
    onUI {
        ready.panel.installBtn.doClick()
        assert ready.panel.installBusy.get() && !control(ready.dialog,'Next').enabled && !control(ready.dialog,'Cancel').enabled
        ready.dialog.dispatchEvent(new WindowEvent(ready.dialog,WindowEvent.WINDOW_CLOSING))
        assert ready.dialog.displayable : 'Keep restart/retry feedback available until install finishes'
        assert all(ready.dialog).find { it instanceof JTextField }.text=='chosen folder'
    }
    installGate.countDown()
    await { !ready.panel.installBusy.get() }
    onUI {
        assert ready.panel.installBar.value==100 && ready.panel.updateNote.text.contains('restart Fiji')
        assert control(ready.dialog,'Next').enabled && control(ready.dialog,'Cancel').enabled
        assert !ready.panel.updateBtn.enabled && !ready.panel.installBtn.enabled
        ready.dialog.setSize(420,360)
    }
    // Native resize events arrive after setSize returns. Check and render the settled layout.
    await { ready.dialog.contentPane.width<=420 && ready.dialog.contentPane.height<=360 }
    onUI {
        ready.dialog.validate()
        def scroll=all(ready.dialog).find { it instanceof JScrollPane }
        assert scroll.verticalScrollBar.visible && scroll.wheelScrollingEnabled
        def next=control(ready.dialog,'Next')
        assert !SwingUtilities.isDescendingFrom(next,scroll) : 'Next must stay outside the scrolling content'
        assert ready.dialog.contentPane.bounds.contains(SwingUtilities.convertRectangle(next.parent,next.bounds,ready.dialog.contentPane))
        ready.dialog.setSize(680,700)
    }
    await { ready.dialog.contentPane.width>=600 && ready.dialog.contentPane.height>=600 }
    onUI {
        ready.dialog.validate()
        // Save a rendering of this real UI for visual review without capturing the user's desktop.
        def rendered=new java.awt.image.BufferedImage(ready.dialog.width,ready.dialog.height,java.awt.image.BufferedImage.TYPE_INT_RGB)
        def g=rendered.createGraphics();try{ready.dialog.paint(g)}finally{g.dispose()}
        def file=new File('build/start-update-screen.png');file.parentFile.mkdirs();javax.imageio.ImageIO.write(rendered,'png',file)
        control(ready.dialog,'Cancel').doClick()
    }
    ready.thread.join(5000);assert !ready.thread.alive

    ij.Prefs.set('calcium.steps.check_updates',true)
    def offline=open { panel -> panel.updateCheck.fetch={ -> throw new IOException('offline fixture') } }
    await { !offline.panel.updateCheck.busy.get() }
    onUI {
        assert control(offline.dialog,'Next').enabled && !offline.panel.updateNote.text.contains('offline fixture')
        offline.panel.updateBtn.doClick()
    }
    await { offline.panel.updateNote.text.contains('offline fixture') }
    onUI { assert offline.panel.updateBtn.enabled; control(offline.dialog,'Cancel').doClick() }
    offline.thread.join(5000)

    // The plugin's file was moved as it started: the note stays above the update messages.
    Start.notice='<b>The plugin is now at Plugins › '+Release.NAME+'.</b> Restart Fiji before starting it from the menu again.'
    def moved=open { panel -> panel.updateCheck.fetch={ -> [version:Release.VERSION,page:Release.REPOSITORY+'/releases',installable:true] } }
    await { moved.panel.updateNote.text.contains('up to date') }
    onUI {
        assert moved.dialog.title==Release.NAME
        assert moved.panel.placeNote.visible && moved.panel.placeNote.text.contains('Plugins › '+Release.NAME)
        moved.dialog.validate()
        def rendered=new java.awt.image.BufferedImage(moved.dialog.width,moved.dialog.height,java.awt.image.BufferedImage.TYPE_INT_RGB)
        def g=rendered.createGraphics();try{moved.dialog.paint(g)}finally{g.dispose()}
        javax.imageio.ImageIO.write(rendered,'png',new File('build/start-moved-note.png'))
        control(moved.dialog,'Cancel').doClick()
    }
    moved.thread.join(5000);Start.notice=null
    println 'PASS: first-screen startup/manual checks, no recording required, preference, Later, offline/retry, close during lookup, install navigation guard, restart feedback, small-screen scrolling and the moved-plugin note'
} catch(Throwable failure) {failure.printStackTrace();System.exit(1)}
finally {
    onUI { panels.each {it.close()};launched.each {it.dialog?.dispose()} }
    ij.Prefs.set('calcium.steps.check_updates',oldPreference)
}
System.exit(0)

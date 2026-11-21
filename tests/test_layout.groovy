// Exercise the real plugin controls in resized native windows, including Windows' look and feel.
import javax.swing.*
import java.awt.*
import java.awt.event.MouseWheelEvent

def loader = new GroovyClassLoader(getClass().classLoader)
def source = new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy').getText('UTF-8')
loader.parseClass(source)
def fragment = source.substring(source.indexOf('def GRADE_SHARE ='), source.indexOf('if (!testGui) { def scr'))
def makeUI = {
    def binding = new Binding([
        imp:new ij.ImagePlus('Resize regression recording', new ij.process.ByteProcessor(50,50)),
        steps:new Expando(n:300,found:true,onset:120), look:true,
        outDir:new File('build/results/recording').absoluteFile,
    ])
    new GroovyShell(loader,binding).evaluate('import ij.IJ\nimport ij.Prefs\nimport java.awt.Color\n' +
        fragment.replace('frame.setAlwaysOnTop(true)', 'frame.setAlwaysOnTop(false)') +
        '\nreturn [frame:frame,root:root,scroll:scroll,progress:detectionProgress,close:closeBtn,go:go,add:removeBtn,trace:trace,info:info,\n' +
        '    openStep:openStep,choose:s2,measure:s3,analyze:s4,next:continueBtn,analyzeButton:analyzeBtn,\n' +
        '    suggest:suggestBtn,accept:acceptAllBtn,responders:addRespondersBtn,views:viewSwitch,status:status,nextLine:nextLine,tick:removeMode,name:recordingName]')
}
def onUI = { Closure action -> SwingUtilities.invokeAndWait(action as Runnable) }
def settle = { ui ->
    // Native resize events arrive asynchronously; wait for the viewport to catch up.
    boolean ready = false
    long end = System.nanoTime() + 5000000000L
    while (!ready && System.nanoTime() < end) {
        onUI {
            ui.frame.validate()
            ready = ui.scroll.width == ui.frame.contentPane.width && ui.scroll.height == ui.frame.contentPane.height
        }
        if (!ready) Thread.sleep(20)
    }
    assert ready : 'Viewport did not follow the resized window'
    onUI { ui.frame.validate() }
}
def controlsFit = { ui ->
    assert ui.root.height >= ui.root.preferredSize.height
    assert ui.root.width >= ui.root.preferredSize.width
    assert !ui.trace.visible || ui.trace.height >= ui.trace.preferredSize.height   // the trace takes room only when there is one
    def walk; walk = { c ->
        if (!c.visible) return
        c.components.each { child ->
            if (child.visible) {
                assert child.x >= 0 && child.y >= 0 && child.x + child.width <= c.width && child.y + child.height <= c.height :
                    "Clipped ${child.class.simpleName}: ${child.bounds} in ${c.size}"
                if (child instanceof AbstractButton) {
                    assert child.width >= child.preferredSize.width && child.height >= child.preferredSize.height : "Squashed button: ${child.text}"
                }
                if (child instanceof Container) walk(child)
            }
        }
    }
    walk(ui.root)
}
def reveal = { ui, button ->
    def rectangle = SwingUtilities.convertRectangle(button.parent, button.bounds, ui.root)
    ui.root.scrollRectToVisible(rectangle)
    assert ui.scroll.viewport.viewRect.contains(rectangle) : "Cannot scroll to ${button.text}"
}

try {
    [UIManager.systemLookAndFeelClassName, UIManager.crossPlatformLookAndFeelClassName].unique().each { look ->
        def ui
        try {
            onUI { UIManager.setLookAndFeel(look); ui = makeUI(); ui.frame.setLocation(20,20); ui.frame.setVisible(true) }
            settle(ui)
            onUI { assert ui.frame.isShowing(); ui.frame.pack() } // Regression: Groovy Double coordinates must not break a visible pack.
            settle(ui)
            onUI { ui.frame.setSize(340,360) }
            settle(ui)
            onUI {
                controlsFit(ui)
                assert ui.scroll.verticalScrollBar.visible && ui.scroll.horizontalScrollBar.visible
                assert ui.scroll.wheelScrollingEnabled
                ui.scroll.viewport.setViewPosition(new Point(0,0))
                ui.scroll.dispatchEvent(new MouseWheelEvent(ui.scroll, MouseWheelEvent.MOUSE_WHEEL, System.currentTimeMillis(), 0,
                    100,100,0,false,MouseWheelEvent.WHEEL_UNIT_SCROLL,3,4))
                assert ui.scroll.verticalScrollBar.value > 0 : 'Mouse wheel did not scroll'
                // Steps 2 to 4 share the space: opening one closes the others, and each one's controls can be scrolled to.
                assert ui.choose.expanded && !ui.measure.expanded && !ui.analyze.expanded
                assert !ui.go.isShowing() && ui.add.isShowing()
                ui.next.doClick(); ui.frame.validate()
                assert ui.measure.expanded && !ui.choose.expanded && !ui.add.isShowing()
                reveal(ui,ui.go); reveal(ui,ui.close)
                ui.openStep(ui.analyze); ui.frame.validate()
                assert ui.analyze.expanded && !ui.measure.expanded
                reveal(ui,ui.analyzeButton)
                ui.openStep(ui.choose); ui.frame.validate()
                assert ui.choose.expanded && !ui.analyze.expanded
                reveal(ui,ui.add)
                // Every label and button uses the look and feel's own interface font, never a fixed family.
                String family = UIManager.getFont('Label.font').family
                def fonts; fonts = { c -> if (c instanceof JLabel || c instanceof AbstractButton) assert c.font.family == family : "Fixed font on ${c.class.simpleName}: ${c.font}"
                    if (c instanceof Container) c.components.each { fonts(it) } }
                fonts(ui.root)
                ui.progress.begin().report(40,'Finding cells','Long-running CPU detection')
                ui.frame.refreshContents()
            }
            settle(ui)
            onUI {
                assert ui.frame.size == new Dimension(340,360) : 'Progress or update status changed the chosen window size'
                controlsFit(ui)
                reveal(ui,ui.progress.detailsButton); reveal(ui,ui.close)
                ui.progress.finish(true,'444 suggestions; Cellpose finished')
                ui.info.setText('<html><body style="width:300px">Suggestions are ready to review.</body></html>')
                ui.frame.refreshContents()
                assert ui.frame.size == new Dimension(340,360)
                ui.frame.setSize(740,620)
            }
            settle(ui)
            onUI {
                controlsFit(ui)
                assert !ui.scroll.horizontalScrollBar.visible : 'Horizontal scrollbar should disappear when the layout fits'
                reveal(ui,ui.close)
                // Text fills its card and follows the window; buttons, ticks and text share one left edge.
                def left = { c -> SwingUtilities.convertPoint(c.parent, c.location, ui.choose.body).x }
                int inner = ui.choose.body.width - ui.choose.body.insets.left - ui.choose.body.insets.right
                assert ui.info.width == inner && ui.status.width == inner && ui.nextLine.width == inner : 'Text is as wide as its card'
                assert [ui.add, ui.suggest, ui.responders, ui.tick, ui.status, ui.nextLine].collect { left(it) }.unique() == [ui.choose.body.insets.left] : 'One left edge'
                ui.nextLine.setText('<html>' + ('several words that have to wrap ' * 12) + '</html>'); ui.frame.validate()
                int wide = ui.nextLine.height; assert wide > ui.status.height : 'A long sentence takes the lines it needs'
                // The four actions keep their places and sizes whatever they say and whichever is the blue one.
                def places = { [ui.suggest, ui.accept, ui.responders, ui.next].collect { SwingUtilities.convertRectangle(it.parent, it.bounds, ui.choose.body) } }
                def before = places()
                ui.accept.setText('Accept responding and likely (888)'); ui.suggest.setText('Suggest cells again'); ui.accept.setEnabled(true)
                [ui.suggest, ui.accept, ui.responders, ui.next].each { b -> loader.loadClass('CalciumTheme').primary(b, b.is(ui.accept)) }
                ui.frame.validate()
                assert places() == before && ui.accept.width >= ui.accept.getFontMetrics(ui.accept.font).stringWidth(ui.accept.text) + 20 : 'Labels and the blue button change; places do not'
                assert before[0].y == before[1].y && before[2].y == before[3].y && before[0].x == before[2].x && before[1].x == before[3].x && before[0].width == before[2].width && before[1].width == before[3].width
                // The picture choice is one joined control: each part starts on the last pixel of the one before.
                def parts = ui.views.parts; assert parts.size() == 3 && (1..2).every { parts[it].x == parts[it - 1].x + parts[it - 1].width - 1 && parts[it].height == parts[0].height }
                // A trace takes its room when it is shown, and nothing else is squeezed.
                ui.trace.setVisible(true); ui.frame.refreshContents(); ui.frame.validate(); controlsFit(ui); assert ui.trace.height >= ui.trace.preferredSize.height; ui.trace.setVisible(false)
                // A long file name is cut in its middle to the width there is; the whole name is one hover away.
                ui.name.setWhole('2026-11-18_' + ('a-very-long-part_' * 12) + 'end.zvi'); ui.frame.validate()
                assert ui.name.width == inner && ui.name.toolTipText.endsWith('end.zvi')
                def fm = ui.name.getFontMetrics(ui.name.font); String shown = loader.loadClass('CalciumTheme').fitMiddle(ui.name.whole, fm, ui.name.width)
                assert shown.startsWith('2026-11-18_') && shown.endsWith('end.zvi') && shown.contains('…') && fm.stringWidth(shown) <= ui.name.width
                ui.frame.setSize(320,260)
            }
            settle(ui)
            onUI { controlsFit(ui); reveal(ui,ui.close); assert ui.nextLine.height > 0 && ui.nextLine.width < 700 : 'Text follows a narrower window too' }
            println 'PASS: step cards, interface font, resize, uncompressed controls, wheel scrolling, both scrollbars, bottom/right controls, stable size, text that follows the width, one left edge, actions that keep their places, joined picture switch (' + look + ')'
        } finally {
            onUI { ui?.progress?.close(); ui?.frame?.dispose() }
        }
    }
} catch (Throwable failure) {
    failure.printStackTrace()
    System.exit(1)
}
System.exit(0)

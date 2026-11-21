// Exercise the real Step 2 controls and canvas listeners against a native ROI Manager.
import ij.*
import ij.gui.*
import ij.plugin.frame.RoiManager
import ij.process.ColorProcessor
import javax.swing.*
import java.awt.Color
import java.awt.event.MouseEvent

def source = new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy').getText('UTF-8')
def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(source)
def imports = source.substring(source.indexOf('import groovy.'), source.indexOf('@CompileStatic'))
def construction = source.substring(source.indexOf('def GRADE_TEXT ='), source.indexOf('if (!testGui) { def scr'))
def logic = source.substring(source.indexOf('def lastKey ='), source.indexOf('def findCellposePython ='))
def listeners = source.substring(source.indexOf('def updateEditHint ='), source.indexOf('showCellposeState(null)\nsetupBtn'))
def onUI = { Closure c -> SwingUtilities.invokeAndWait(c as Runnable) }
def ui, manager, picture, imagej, pluginMouse
def cells = [new OvalRoi(60,60,30,30), new OvalRoi(120,60,30,30), new OvalRoi(180,60,30,30)]
def counts = { -> manager.getCount() }
def click = { int x, int y, int button = MouseEvent.BUTTON1, int drag = 0, int clicks = 1 ->
    onUI {
        def canvas = picture.getCanvas()
        def pressed = new MouseEvent(canvas, MouseEvent.MOUSE_PRESSED, System.currentTimeMillis(), 0, x,y,clicks,false,button)
        def released = new MouseEvent(canvas, MouseEvent.MOUSE_RELEASED, System.currentTimeMillis(), 0, x+drag,y,clicks,false,button)
        if (button == MouseEvent.BUTTON1) { canvas.dispatchEvent(pressed); canvas.dispatchEvent(released) }
        else {
            // Exercise the actual plugin guard without opening ImageJ's unrelated native popup,
            // which waits for human dismissal on Windows and blocks a synthetic event dispatch.
            pluginMouse.mousePressed(pressed); pluginMouse.mouseReleased(released)
        }
    }
    onUI {} // drain the actual listener's invokeLater
}
try {
    onUI {
        imagej = new ImageJ(ImageJ.NO_SHOW)
        picture = new ImagePlus('Cell editing regression',new ColorProcessor(280,180))
        def pixels = picture.getProcessor().getPixels()
        for (int y=0;y<180;y++) for(int x=0;x<280;x++) pixels[y*280+x] = x < 100 ? 0x00a000 : (x < 170 ? 0xff00ff : 0xffff00)
        picture.show(); picture.getCanvas().setMagnification(1d)
        def originalMouseListeners = picture.getCanvas().getMouseListeners() as List
        manager = new RoiManager(false)
        manager.runCommand(picture,'Show None')
        def steps = new Expando(n:300,found:false,onset:0,colourPic:picture,movie:picture,respDiff:null)
        steps.nearEdge = { r -> r.bounds.x < 40 }
        steps.backgrounds = { rois -> rois.collect { [0] as int[] } as int[][] }
        steps.traceOf = { r, rois -> [1d,2d,1d] as double[] }
        def binding = new Binding(imp:picture,steps:steps,look:false,rm:manager,testGui:true,outDir:new File('build/results/cells').absoluteFile)
        ui = new GroovyShell(loader,binding).evaluate(imports + construction.replace('frame.setAlwaysOnTop(true)','frame.setAlwaysOnTop(false)') +
            '\nhiddenKinds.add("background")\n' + logic + listeners +
            '\nreturn [frame:frame, progress:detectionProgress, proposals:proposals, refresh:refresh, click:clickCellAt, remove:removeBtn, undo:undoBtn, mode:removeMode, hint:editHint, hidden:hiddenKinds, edge:edgeBox, add:addBtn, oval:ovalBtn, free:freeBtn, acceptAll:acceptAll, edits:cellEdits, status:status]')
        def installed = picture.getCanvas().getMouseListeners().findAll { !originalMouseListeners.contains(it) }
        assert installed.size() == 1
        pluginMouse = installed.first()
        ui.frame.visible = true
        IJ.setTool('oval')
        cells.eachWithIndex { r,i -> manager.addRoi(r); manager.rename(i,'Cell '+(i+1)) }
        ui.hidden.addAll(['g0', 'g1', 'g2', 'g3'])   // no suggestions on the picture
        ui.refresh(true)
        assert !ui.status.text.contains('Exception') : ui.status.text
    }
    // Regression: a real click on a manual cell selects its ROI even without suggestions.
    click(135,75)
    onUI {
        assert manager.selectedIndex == 1 : 'Picture click did not select the accepted cell'
        assert picture.getRoi() == null : 'Yellow editing handles obscure the new outline'
        ui.remove.doClick()
        assert counts() == 2 && manager.getName(1) == 'Cell 3'
        assert ui.undo.enabled
        ui.undo.doClick()
        assert counts() == 3 && (0..<3).collect { manager.getName(it) } == ['Cell 1','Cell 2','Cell 3']
        assert manager.getRoi(1).bounds == cells[1].bounds
        ui.remove.doClick() // no selection must never delete all cells
        assert counts() == 3
        ui.mode.doClick()
    }
    click(75,75)
    click(195,75)
    onUI {
        assert counts() == 1 && manager.getName(0) == 'Cell 2'
        ui.undo.doClick(); ui.undo.doClick()
        assert (0..<3).collect { manager.getName(it) } == ['Cell 1','Cell 2','Cell 3']
        assert ui.mode.selected
    }
    click(10,150) // empty space
    click(75,75,MouseEvent.BUTTON3) // context click
    click(75,75,MouseEvent.BUTTON1,10) // drawing drag
    click(75,75,MouseEvent.BUTTON1,0,2) // second event of double click
    onUI { assert counts() == 3 }
    println 'PASS: actual canvas selection, direct removals, original-order undo, empty/right/drag/double click guards'

    onUI {
        ui.proposals << [roi:new OvalRoi(120,120,30,30),cx:135,cy:135,grade:1,edge:false,ev:[1,900,0,0,1,0,0,0] as double[]]
        ui.hidden.removeAll(['g0', 'g1', 'g2', 'g3'])
    }
    click(135,135)
    onUI { assert counts() == 3 : 'Removal mode accepted a suggestion'; ui.mode.doClick() }
    click(135,135)
    onUI {
        assert counts() == 4
        assert picture.getOverlay().toArray().count { it.strokeColor == Color.white } == 4
        ui.undo.doClick(); assert counts() == 3
        assert ui.acceptAll() == 1 && counts() == 4
        ui.undo.doClick(); assert counts() == 3
        ui.mode.doClick(); ui.oval.doClick(); assert !ui.mode.selected
        ui.mode.doClick(); ui.free.doClick(); assert !ui.mode.selected
        IJ.setTool('oval')
        picture.setRoi(new OvalRoi(200,120,25,25))
        ui.add.doClick(); assert counts() == 4
        ui.undo.doClick(); assert counts() == 3
        manager.setSelectedIndexes([0,2] as int[])
        ui.remove.doClick(); assert counts() == 1
        ui.undo.doClick(); assert counts() == 3
        // Non-area ROIs must not shift targeting or be deleted by cell controls.
        manager.addRoi(new Line(10,10,30,30)); manager.rename(3,'Reference line')
        manager.setSelectedIndexes([1,3] as int[])
        ui.remove.doClick(); assert counts() == 3 && manager.getName(2) == 'Reference line'
        ui.undo.doClick(); assert counts() == 4 && manager.getName(3) == 'Reference line'
        // An external ROI edit invalidates undo instead of erasing the external edit.
        manager.setSelectedIndexes([0] as int[]); ui.remove.doClick()
        manager.addRoi(new Roi(10,100,10,10))
        ui.undo.doClick()
        assert counts() == 4 && ui.hint.text.contains('outside this window')
        assert !ui.undo.enabled
        // Overlapping accepted cells target the smaller outline deterministically.
        manager.reset(); manager.addRoi(new Roi(60,60,60,60)); manager.addRoi(new Roi(70,70,10,10))
        assert ui.edits.hit(75,75) == 1
        ui.mode.doClick(); ui.click(75,75)
        assert counts() == 1 && manager.getRoi(0).bounds.width == 60
        ui.undo.doClick(); assert counts() == 2
    }
    println 'PASS: removal never accepts, undo additions/Accept all, drawing exits removal mode, bulk deletion, non-area ROIs, outside edits, overlapping cells'

    onUI {
        manager.reset()
        cells.each { manager.addRoi(it) }
        manager.addRoi(new Roi(10,60,20,20)) // excluded edge cell remains dashed
        ui.proposals.clear(); manager.deselect(); ui.refresh(true)
        def overlay = picture.getOverlay().toArray()
        assert overlay.size() == 8
        (0..<4).each { int i ->
            assert overlay[i*2].strokeColor == Color.black && overlay[i*2+1].strokeColor == Color.white
            assert overlay[i*2].strokeWidth > overlay[i*2+1].strokeWidth
            assert overlay[i*2].bounds == manager.getRoi(i).bounds
            assert overlay[i*2+1].bounds == manager.getRoi(i).bounds
            assert (overlay[i*2+1].stroke.dashArray != null) == (i==3)
        }
        assert picture.getCanvas().showAllList == null : 'ROI Manager yellow outlines cover the plugin overlay'
        def before = picture.getProcessor().getPixels().clone()
        def flattened = picture.flatten()
        def rendered = flattened.getProcessor().getPixels()
        int whites = rendered.toList().count { (it & 0xffffff) == 0xffffff }
        int blacks = rendered.toList().count { (it & 0xffffff) == 0 }
        assert whites > 100
        assert blacks > 100
        assert picture.getProcessor().getPixels() == before : 'Outline rendering changed image pixels'
        new File('build').mkdirs()
        IJ.saveAs(flattened,'PNG','build/cell-outline-contrast.png')
        flattened.close()
        ui.mode.selected = false; ui.click(75,75)
        assert picture.getOverlay().get(1).strokeWidth > overlay[1].strokeWidth : 'Selected cell has no thicker highlight'
        manager.setSelectedIndexes([0,1] as int[]); ui.refresh(false)
        assert picture.getOverlay().get(3).strokeWidth > overlay[3].strokeWidth : 'Adding to the same selection did not refresh its highlight'
    }
    println 'PASS: black/white rendered outlines over green/magenta/yellow, edge dashes, selected highlight, original geometry and pixels preserved'
} catch (Throwable failure) {
    failure.printStackTrace()
    System.exit(1)
} finally {
    onUI { ui?.progress?.close(); ui?.frame?.dispose(); manager?.close(); picture?.close(); imagej?.dispose() }
}
System.exit(0)

import javax.swing.*

def loader = new GroovyClassLoader(getClass().classLoader)
def source = new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy').getText('UTF-8')
loader.parseClass(source)
def Panel = loader.loadClass('CalciumOutlineProgress')
def BrokenPanel = loader.parseClass('''class BrokenBeginPanel extends CalciumOutlineProgress {
    CalciumOutlineTask begin() { super.begin(); throw new IOException('simulated panel startup error') }
}''')
def fragment = source.substring(source.indexOf('def detectionRunner ='),source.indexOf('acceptAllBtn.addActionListener'))
['window','launch','begin','worker','completion'].each { failure ->
    def button = new JButton('Suggest cells automatically')
    def panel = (failure == 'begin' ? BrokenPanel : Panel).newInstance()
    def info = new JLabel()
    def calls = new java.util.concurrent.atomic.AtomicInteger()
    def binding = new Binding([suggestBtn:button,detectionProgress:panel,info:info,proposals:[],forceFresh:false,cellposeBusy:false,detecting:false,hiddenKinds:[] as Set,paintLegend:{ -> },steps:new Expando(blind:false),
        showCellposeState:{ state -> },
        frame:new Expando(refreshContents:{ -> if(failure=='window') throw new IOException('simulated window error') }),
        refresh:{ force -> if(failure=='completion') throw new IOException('simulated result display error') },
        suggest:{ task -> calls.incrementAndGet(); if(failure=='worker') throw new IOException('simulated detection error'); task.succeeded=true; '3 suggestions' }])
    def runner = new GroovyShell(loader,binding).evaluate(fragment+'\nreturn detectionRunner')
    if (failure == 'launch') runner.launch = { job -> throw new IOException('simulated thread startup error') }
    SwingUtilities.invokeAndWait { button.doClick() }
    long end=System.nanoTime()+5000000000L
    boolean finished=false
    while(!finished && System.nanoTime()<end) {
        SwingUtilities.invokeAndWait { finished=button.enabled && !panel.running && !panel.clock.running }
        if(!finished) Thread.sleep(10)
    }
    SwingUtilities.invokeAndWait {
        assert finished : 'Detection stayed disabled/running after '+failure
        assert calls.get() == (failure in ['launch','begin'] ? 0 : 1)
        assert button.text == 'Suggest cells automatically'
        if (failure in ['window','completion']) assert panel.bar.value == 100
        else { assert panel.bar.value < 100; assert info.text.contains(failure=='worker' ? 'Failed:' : 'Could not start detection:') }
        panel.close()
    }
}
println 'PASS: actual Suggest listener survives window errors and resets timer/button after panel, worker-launch, detection and completion failures'

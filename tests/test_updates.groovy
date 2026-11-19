def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Release = loader.loadClass('CalciumRelease')
assert Release.newer('v2026.11.19.10', '2026.11.19.2')
assert !Release.newer('v2026.11.19.2', '2026.11.19.2')
assert !Release.newer('2026.10.30.9', '2026.11.19.2')
assert Release.newer('2027.1.1.1', '2026.12.31.9')
def valid = { version ->
    def tag = 'v' + version
    [tag_name:tag, html_url:Release.REPOSITORY+'/releases/tag/'+tag, draft:false, prerelease:false,
     assets:[[name:Release.ASSET,state:'uploaded',size:100,browser_download_url:Release.REPOSITORY+'/releases/download/'+tag+'/'+Release.ASSET]]]
}
def nextVersion = (Release.versionParts(Release.VERSION)[0..2] + [Release.versionParts(Release.VERSION)[3] + 1]).join('.')
def next = valid(nextVersion)
assert Release.describe(next).newer
assert !Release.describe(valid(Release.VERSION)).newer
def bad = []
bad << (valid('2026.11.19.4') + [draft:true])
bad << (valid('2026.11.19.4') + [prerelease:true])
bad << (valid('2026.11.19.4') + [assets:[]])
bad << (valid('2026.11.19.4') + [html_url:'https://example.com/fake'])
bad << (valid('2026.11.19.4') + [tag_name:'v2026.11.19.4-rc1'])
def wrongDownload = valid('2026.11.19.4'); wrongDownload.assets[0].browser_download_url='https://example.com/plugin.zip'; bad << wrongDownload
def incomplete = valid('2026.11.19.4'); incomplete.assets[0].state='new'; bad << incomplete
bad.each { data ->
    boolean rejected = false
    try { Release.describe(data) } catch (IOException expected) { rejected = true }
    assert rejected : data
}
println 'PASS: numeric version comparison, stable complete releases, expected repository and download URLs'

def Check = loader.loadClass('CalciumUpdateCheck')
def callCheck = { checker, automatic ->
    def latch = new java.util.concurrent.CountDownLatch(1)
    def result = new java.util.concurrent.atomic.AtomicReference()
    javax.swing.SwingUtilities.invokeAndWait {
        assert checker.check(automatic, { answer ->
            assert javax.swing.SwingUtilities.isEventDispatchThread()
            result.set(answer); latch.countDown()
        })
    }
    assert latch.await(5, java.util.concurrent.TimeUnit.SECONDS)
    result.get()
}
def checker = Check.newInstance(fetch: { ->
    assert !javax.swing.SwingUtilities.isEventDispatchThread()
    Release.describe(next)
})
assert callCheck(checker, true).prompt
assert !callCheck(checker, true).prompt // Later lasts for this plugin window.
assert callCheck(checker, false).prompt // A manual check may offer it again.
checker.fetch = { -> Release.describe(valid(Release.VERSION)) }
assert !callCheck(checker, true).prompt
checker.fetch = { -> throw new IOException('offline') }
assert callCheck(checker, true).error == 'offline'
checker.enabled = false
assert !checker.check(true, { throw new AssertionError('disabled startup') })
assert callCheck(checker, false).error == 'offline'
checker.close()
assert !checker.check(false, { throw new AssertionError('closed window') })
println 'PASS: background update checks, event-thread callbacks, newer/current/offline, Later and disabled startup'

// A blocked request must not freeze Swing, duplicate work, or update a closed/disabled window.
[true, false].each { closeWindow ->
    def entered = new java.util.concurrent.CountDownLatch(1)
    def releaseRequest = new java.util.concurrent.CountDownLatch(1)
    def slow = Check.newInstance(fetch: { ->
        entered.countDown()
        assert releaseRequest.await(5, java.util.concurrent.TimeUnit.SECONDS)
        Release.describe(next)
    })
    boolean delivered = false
    javax.swing.SwingUtilities.invokeAndWait {
        assert slow.check(true, { delivered = true })
    }
    assert entered.await(5, java.util.concurrent.TimeUnit.SECONDS)
    javax.swing.SwingUtilities.invokeAndWait {
        assert !slow.check(false, { delivered = true })
        if (closeWindow) slow.close() else slow.enabled = false
    }
    releaseRequest.countDown()
    long end = System.nanoTime() + 5000000000L
    while (slow.busy.get() && System.nanoTime() < end) Thread.sleep(10)
    javax.swing.SwingUtilities.invokeAndWait { assert !delivered && !slow.busy.get() }
}
println 'PASS: responsive event thread, concurrent request suppression, close/disable during a request'

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
def next = valid('2026.11.19.3')
assert Release.describe(next).newer
assert !Release.describe(valid(Release.VERSION)).newer
def bad = []
bad << (valid('2026.11.19.3') + [draft:true])
bad << (valid('2026.11.19.3') + [prerelease:true])
bad << (valid('2026.11.19.3') + [assets:[]])
bad << (valid('2026.11.19.3') + [html_url:'https://example.com/fake'])
bad << (valid('2026.11.19.3') + [tag_name:'v2026.11.19.3-rc1'])
def wrongDownload = valid('2026.11.19.3'); wrongDownload.assets[0].browser_download_url='https://example.com/plugin.zip'; bad << wrongDownload
def incomplete = valid('2026.11.19.3'); incomplete.assets[0].state='new'; bad << incomplete
bad.each { data ->
    boolean rejected = false
    try { Release.describe(data) } catch (IOException expected) { rejected = true }
    assert rejected : data
}
println 'PASS: numeric version comparison, stable complete releases, expected repository and download URLs'

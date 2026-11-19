import java.nio.file.Files
import java.nio.file.StandardOpenOption
import java.nio.channels.FileChannel

def loader = new GroovyClassLoader(getClass().classLoader)
def source = new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy').getText('UTF-8')
loader.parseClass(source)
def Installer = loader.loadClass('CalciumInstaller')
def Release = loader.loadClass('CalciumRelease')
String next = (Release.versionParts(Release.VERSION)[0..2] + [Release.versionParts(Release.VERSION)[3]+1]).join('.')
byte[] original = source.getBytes('UTF-8')
byte[] payload = source.replace(Release.VERSION,next).getBytes('UTF-8')
def base = Release.REPOSITORY+'/releases/download/v'+next+'/'
def release = [version:next,installable:true,script:base+Release.SCRIPT,checksums:base+Release.CHECKSUMS,size:payload.length]
byte[] manifest = (Installer.digest(payload)+'  '+Release.SCRIPT+'\n').getBytes('UTF-8')
def download = { url,limit,progress ->
    if (url == release.checksums) return manifest
    assert url == release.script && payload.length < limit
    progress.call(payload.length,payload.length)
    payload
}
def directory = Files.createTempDirectory('Fiji install with spaces ').toFile()
def target = new File(directory,'scripts/Plugins/Calcium Imaging/'+Release.SCRIPT)
def events=[]
def instance = Installer.newInstance(applicationDirectory:directory,download:download)
def reset = { -> target.parentFile.mkdirs(); target.bytes = original }
def unchangedAfterFailure = { operation ->
    boolean failed=false
    try { operation.call() } catch (IOException expected) { failed=true }
    assert failed
    assert target.bytes == original : 'Failure modified the installed plugin'
    assert !target.parentFile.listFiles().any{it.name.startsWith('.calcium-update-')}
}
try {
    def fresh = instance.install(release,{p,message -> events << p})
    assert target.bytes == payload && fresh.backup == null && !fresh.already
    assert events.last() == 100 && events == events.sort(false)
    reset()
    def upgraded = instance.install(release)
    assert target.bytes == payload && upgraded.backup.bytes == original
    assert upgraded.backup.name.endsWith('.groovy.bak') && upgraded.backup.parentFile.name == 'calcium-plugin-backups'
    def again = instance.install(release)
    assert again.already && again.backup == null
    println 'PASS: first install, atomic upgrade with exact backup, spaces in path, already-installed result'

    reset()
    unchangedAfterFailure { Installer.newInstance(applicationDirectory:directory,download:{a,b,c -> throw new IOException('network interrupted')}).install(release) }
    unchangedAfterFailure { Installer.newInstance(applicationDirectory:directory,download:{url,limit,progress -> url==release.checksums ? manifest : 'broken'.bytes}).install(release) }
    unchangedAfterFailure { instance.install(release+[script:'https://example.com/plugin.groovy']) }
    byte[] badVersion = original
    unchangedAfterFailure { Installer.newInstance(applicationDirectory:directory,download:{url,limit,progress -> url==release.checksums ? (Installer.digest(badVersion)+'  '+Release.SCRIPT+'\n').bytes : badVersion}).install(release) }
    unchangedAfterFailure { Installer.newInstance(applicationDirectory:directory,download:download,replaceFile:{a,b -> throw new java.nio.file.AccessDeniedException(b.toString())}).install(release) }
    unchangedAfterFailure { Installer.newInstance(applicationDirectory:directory,download:download,replaceFile:{a,b -> throw new java.nio.file.AtomicMoveNotSupportedException(a.toString(),b.toString(),'test filesystem')}).install(release) }
    println 'PASS: interrupted download, wrong checksum/version/URL, denied replacement and unsupported atomic move preserve original'

    def lockPath = new File(directory,'.calcium-plugin-update.lock').toPath()
    FileChannel.open(lockPath,StandardOpenOption.CREATE,StandardOpenOption.WRITE).withCloseable { channel ->
        def lock = channel.lock()
        try { unchangedAfterFailure { instance.install(release) } } finally { lock.release() }
    }
    def backups = new File(directory,'calcium-plugin-backups'); backups.deleteDir(); backups.text='blocked'
    unchangedAfterFailure { instance.install(release) }
    backups.delete()
    target.bytes = source.replace(Release.VERSION,'2099.12.31.9999').getBytes('UTF-8')
    boolean rejected=false
    try { instance.install(release) } catch (IOException expected) { rejected=true }
    assert rejected && Installer.scriptVersion(target.bytes)=='2099.12.31.9999'
    println 'PASS: concurrent installers, failed backup and downgrade protection'
} finally { directory.deleteDir() }

println 'All installer checks passed.'

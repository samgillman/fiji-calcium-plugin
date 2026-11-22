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
byte[] beforeRename = source.replace(Release.VERSION,'2026.11.21.7').getBytes('UTF-8')
def base = Release.REPOSITORY+'/releases/download/v'+next+'/'
def release = [version:next,installable:true,script:base+Release.SCRIPT,checksums:base+Release.CHECKSUMS,size:payload.length]
byte[] manifest = (Installer.digest(payload)+'  '+Release.SCRIPT+'\n').getBytes('UTF-8')
def download = { url,limit,progress ->
    if (url == release.checksums) return manifest
    assert url == release.script && payload.length < limit
    progress.call(payload.length,payload.length)
    payload
}

// The name: menu entry, file and the place used by versions up to 2026.11.22.1.
assert Release.NAME == 'Ca²⁺ Imaging (Wang Lab)'
assert Release.command(Release.PLACE) == Release.NAME
assert Release.command(Release.PLAIN_PLACE) == 'Ca2+ Imaging (Wang Lab)'
assert Release.command(Release.OLD_PLACE) == 'Calcium Cells Step By Step'
assert Release.OLD_PLACE == 'scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'
assert Release.SCRIPT == 'Calcium_Cells_Step_By_Step.groovy' && Release.ASSET == 'Calcium_Imaging_Fiji_Plugin.zip'
assert !Release.newer(Release.FIRST_PLACED, Release.VERSION)
println 'PASS: menu name Ca²⁺ Imaging (Wang Lab) comes from the file name; download names unchanged'

def directory = Files.createTempDirectory('Fiji install with spaces ').toFile()
def target = new File(directory,Release.PLACE)
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
    assert target.bytes == payload && fresh.backup == null && !fresh.already && !fresh.aside && !fresh.left
    assert events.last() == 100 && events == events.sort(false)
    // The file keeps its name on disk: Fiji lists it as Ca²⁺ Imaging (Wang Lab).
    assert target.parentFile.list().contains('Ca²⁺_Imaging_(Wang_Lab).groovy')
    assert Installer.place(directory.canonicalFile) == new File(directory.canonicalFile, Release.PLACE)
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

// Earlier versions lived at Plugins > Calcium Imaging > Calcium Cells Step By Step.
def fiji = { -> Files.createTempDirectory('Fiji place ').toFile() }
def put = { File app, String place, byte[] bytes -> def f = new File(app, place); f.parentFile.mkdirs(); f.bytes = bytes; f }
def backupsOf = { File app -> (new File(app,'calcium-plugin-backups').listFiles() ?: []) as List }

directory = fiji()
try {   // the new installer, run where only an old version is installed
    def old = put(directory, Release.OLD_PLACE, beforeRename)
    def result = Installer.newInstance(applicationDirectory:directory,download:download).install(release)
    def here = new File(directory, Release.PLACE)
    assert here.bytes == payload && !old.exists() && !old.parentFile.exists() : 'failed: here.bytes == payload && !old.exists() && !old.parentFile.exists()'
    assert result.backup.bytes == beforeRename && result.aside.size() == 1 && !result.left && !result.already : 'failed: result.backup.bytes == beforeRename && result.aside.size() == 1 && !result.left && !result'
    assert backupsOf(directory).size() == 1 : 'failed: backupsOf(directory).size() == 1'
    println 'PASS: installing over a version under Calcium Imaging puts that copy in calcium-plugin-backups and removes the empty folder'
} finally { directory.deleteDir() }

directory = fiji()
try {   // what an earlier updater leaves: the new version under Calcium Imaging, next to the Cellpose path file
    def old = put(directory, Release.OLD_PLACE, original)
    def note = put(directory, 'scripts/Plugins/Calcium Imaging/cellpose_python.txt', '/somewhere/bin/python'.bytes)
    def placed = Installer.newInstance(applicationDirectory:directory).settle()
    def here = new File(directory, Release.PLACE)
    assert placed.moved == old.canonicalFile && placed.target == here.canonicalFile && placed.backup == null && !placed.aside && !placed.left : 'failed: placed.moved == old.canonicalFile && placed.target == here.canonicalFile && placed.backup '
    assert here.bytes == original && !old.exists() && note.text == '/somewhere/bin/python' : 'failed: here.bytes == original && !old.exists() && note.text == /somewhere/bin/python'
    assert !backupsOf(directory) : 'failed: !backupsOf(directory)'
    def twice = Installer.newInstance(applicationDirectory:directory).settle()
    assert twice.moved == null && !twice.aside && here.bytes == original : 'failed: twice.moved == null && !twice.aside && here.bytes == original'
    println 'PASS: a new version under Calcium Imaging is moved to Plugins > Ca²⁺ Imaging (Wang Lab); cellpose_python.txt and its folder stay'
} finally { directory.deleteDir() }

directory = fiji()
try {   // what the ZIP leaves: the plain-letter name, renamed to the one with ² and ⁺
    def plain = put(directory, Release.PLAIN_PLACE, original)
    def placed = Installer.newInstance(applicationDirectory:directory).settle()
    def here = new File(directory, Release.PLACE)
    assert placed.moved == plain.canonicalFile && here.bytes == original && !plain.exists() && here.parentFile.exists() : 'the ZIP copy is renamed'
    println 'PASS: the ZIP copy, Ca2+_Imaging_(Wang_Lab).groovy, is renamed to Ca\u00b2\u207a_Imaging_(Wang_Lab).groovy'
} finally { directory.deleteDir() }

directory = fiji()
try {   // only a version from before the rename: it does not know the new place, so nothing moves
    def old = put(directory, Release.OLD_PLACE, beforeRename)
    def placed = Installer.newInstance(applicationDirectory:directory).settle()
    assert placed.moved == null && !placed.aside && old.bytes == beforeRename && !new File(directory, Release.PLACE).exists() : 'failed: placed.moved == null && !placed.aside && old.bytes == beforeRename && !new File(directory,'
    println 'PASS: a version from before the rename is left where it is'
} finally { directory.deleteDir() }

directory = fiji()
try {   // a newer version under Calcium Imaging replaces an older one at the new place, which is backed up
    def here = put(directory, Release.PLACE, original)
    def old = put(directory, Release.OLD_PLACE, payload)
    new File(old.parentFile, '.DS_Store').bytes = [0] as byte[]
    def placed = Installer.newInstance(applicationDirectory:directory).settle()
    assert placed.moved == old.canonicalFile && placed.backup.bytes == original && here.bytes == payload && !old.parentFile.exists() : 'failed: placed.moved == old.canonicalFile && placed.backup.bytes == original && here.bytes == payl'
    println 'PASS: a newer copy under Calcium Imaging replaces the one at the new place, with a backup; a folder with only .DS_Store is removed'
} finally { directory.deleteDir() }

directory = fiji()
try {   // a copy under Calcium Imaging that is not newer is an extra menu entry: it goes to the backups
    def here = put(directory, Release.PLACE, payload)
    def old = put(directory, Release.OLD_PLACE, beforeRename)
    def placed = Installer.newInstance(applicationDirectory:directory).settle()
    assert placed.moved == null && placed.aside.size() == 1 && placed.aside.values()[0].bytes == beforeRename : 'failed: placed.moved == null && placed.aside.size() == 1 && placed.aside.values()[0].bytes == befo'
    assert here.bytes == payload && !old.exists() : 'failed: here.bytes == payload && !old.exists()'
    println 'PASS: an older copy under Calcium Imaging is put in calcium-plugin-backups when the new place has a newer one'
} finally { directory.deleteDir() }

directory = fiji()
try {   // nothing to do: no lock file, nothing written
    put(directory, Release.PLACE, original)
    def placed = Installer.newInstance(applicationDirectory:directory).settle()
    assert placed.moved == null && !placed.aside && !new File(directory, '.calcium-plugin-update.lock').exists() : 'failed: placed.moved == null && !placed.aside && !new File(directory, .calcium-plugin-update.lock)'
    // another installation holds the lock: nothing changes
    def old = put(directory, Release.OLD_PLACE, payload)
    def lockPath = new File(directory,'.calcium-plugin-update.lock').toPath()
    FileChannel.open(lockPath,StandardOpenOption.CREATE,StandardOpenOption.WRITE).withCloseable { channel ->
        def lock = channel.lock()
        boolean refused = false
        try { Installer.newInstance(applicationDirectory:directory).settle() } catch (IOException expected) { refused = true } finally { lock.release() }
        assert refused && old.bytes == payload && new File(directory, Release.PLACE).bytes == original : 'failed: refused && old.bytes == payload && new File(directory, Release.PLACE).bytes == original'
    }
    // a failed move leaves both copies as they were
    boolean failed = false
    try { Installer.newInstance(applicationDirectory:directory,replaceFile:{a,b -> throw new java.nio.file.AtomicMoveNotSupportedException(a.toString(),b.toString(),'test filesystem')}).settle() }
    catch (IOException expected) { failed = true }
    assert failed && old.bytes == payload && new File(directory, Release.PLACE).bytes == original : 'failed: failed && old.bytes == payload && new File(directory, Release.PLACE).bytes == original'
    println 'PASS: nothing to move writes nothing; a held lock or a failed move changes nothing'
} finally { directory.deleteDir() }

directory = fiji()
def elsewhere = Files.createTempDirectory('Plugin source ').toFile()
try {   // a Calcium Imaging folder that links to another folder (a working copy) is never moved
    def linked = new File(directory, 'scripts/Plugins/Calcium Imaging')
    linked.parentFile.mkdirs()
    boolean made = true
    try { Files.createSymbolicLink(linked.toPath(), elsewhere.toPath()) } catch (Exception unsupported) { made = false }
    if (made) {
        new File(elsewhere, Release.SCRIPT).bytes = original
        def placed = Installer.newInstance(applicationDirectory:directory).settle()
        assert placed.moved == null && !placed.aside : 'a copy in a linked folder was moved: ' + placed
        assert new File(elsewhere, Release.SCRIPT).bytes == original && !new File(directory, Release.PLACE).exists() : 'a copy in a linked folder was changed'
        // A link at the plugin's own place is refused, not written through.
        def place = new File(directory.canonicalFile, Release.PLACE)
        Files.createSymbolicLink(place.toPath(), new File(elsewhere, Release.SCRIPT).toPath())
        String refusal = null
        try { Installer.newInstance(applicationDirectory:directory,download:download).install(release) } catch (IOException expected) { refusal = expected.message }
        assert refusal?.contains('symbolic link') && new File(elsewhere, Release.SCRIPT).bytes == original : 'the installer wrote through a link: ' + refusal
        place.delete()
        println 'PASS: a linked Calcium Imaging folder is left alone; a link at the plugin place is refused'
    } else println 'SKIP: this system does not allow symbolic links here'
} finally { new File(directory, 'scripts/Plugins/Calcium Imaging').delete(); directory.deleteDir(); elsewhere.deleteDir() }

println 'All installer checks passed.'

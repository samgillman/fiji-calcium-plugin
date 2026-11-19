def binding = new Binding([CALCIUM_INSTALL_TEST:true])
def shell = new GroovyShell(getClass().classLoader,binding)
def Installer = shell.evaluate(new File('dist/Install_Calcium_Plugin.groovy'))
assert Installer.name == 'CalciumInstaller'
def Release = shell.classLoader.loadClass('CalciumRelease')
byte[] payload = new File('dist/Calcium_Cells_Step_By_Step.groovy').bytes
byte[] sums = new File('dist/SHA256SUMS.txt').bytes
def base = Release.REPOSITORY+'/releases/download/v'+Release.VERSION+'/'
def release = [version:Release.VERSION,installable:true,script:base+Release.SCRIPT,checksums:base+Release.CHECKSUMS,size:payload.length]
def directory = java.nio.file.Files.createTempDirectory('Fiji bootstrap ').toFile()
try {
    def installer = Installer.newInstance(applicationDirectory:directory,download:{url,limit,progress -> url==release.checksums ? sums : payload})
    def result = installer.install(release)
    assert result.target.bytes == payload && !result.already
    println 'PASS: published standalone installer compiles and installs the packaged plugin in a fresh Fiji folder'
} finally { directory.deleteDir() }

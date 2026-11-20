// Serves the Step 4 page for the small fixture measurement, for the browser check (tools/run_tests.py --browser)
// and for looking at the page while working on it. Stops when build/web-step4/stop appears, or after ten minutes.
import ij.*
import ij.gui.Roi
import ij.process.*

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def root=new File('build/web-step4').absoluteFile;root.mkdirs();new File(root,'stop').delete();new File(root,'url.txt').delete()
def fixture=new File('tests/fixtures/simpleca').absoluteFile
def frames=new ImageStack(250,150)
(1..30).each{int f->def ip=new ShortProcessor(250,150);ip.set(100);if(f>=8&&f<=17){ip.setRoi(30,50,20,20);ip.set(140);ip.resetRoi()};frames.addSlice(''+f,ip)}
def colour=new ColorProcessor(250,150);colour.setColor(0x103a10);colour.fill()
(1..6).each{colour.setColor(0x30c030);colour.fillOval(it*30+2,52,16,16)}
def data=[csv:new File(fixture,'metrics.csv'),folder:root,recording:'Fixture recording',frames:30,hz:2d,baselineStart:1,baselineEnd:5,key:'fixture',
    picture:new ImagePlus('Colour picture',colour),movie:new ImagePlus('Recording',frames),names:(1..6).collect{String.format('Cell%03d',it)},rois:(1..6).collect{new Roi(it*30,50,20,20)} as Roi[]]
def session=loader.loadClass('CalciumAnalysisSession').newInstance(data);session.runtime.home=new File(root,'runtime')
session.open={File target->new File(root,'opened.txt').text=target.name}
session.performAnalysis={rt,input,settings,out,job->new File(root,'settings.txt').text=settings.collect{k,v->k+'='+v}.join('\n')
    (1..4).each{job.check();job.progress.call(it*20,'Scoring responses');Thread.sleep(150)};fixture}
def server=loader.loadClass('CalciumWebServer').newInstance(session:session);server.start()
new File(root,'url.txt').text=server.url();println 'Step 4 page: '+server.url()
long end=System.currentTimeMillis()+600000L
while(!new File(root,'stop').exists()&&System.currentTimeMillis()<end)Thread.sleep(200)
server.stop();System.exit(0)

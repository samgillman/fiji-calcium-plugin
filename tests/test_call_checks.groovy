// What Measure saves around each cell (ring and core, field movement) and the checks Step 4 makes on the calls from it:
// a possible passing object, a rise on an edge, a field that moved, short separate crossings. Synthetic recordings only.
// Runs without a display, without R and without a browser; the checks are read as the Step 4 page reads them, over HTTP.
import ij.*
import ij.gui.OvalRoi
import ij.gui.Roi
import ij.process.*
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def loader=new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Rings=loader.loadClass('CalciumRingAndCore'),Passing=loader.loadClass('CalciumPassingObject'),Checks=loader.loadClass('CalciumCallChecks'),Steps=loader.loadClass('CalciumSteps')
def Motion=loader.loadClass('CalciumMotion'),View=loader.loadClass('CalciumFrameView')
def Csv=loader.loadClass('CalciumCsv'),Measurement=loader.loadClass('CalciumMeasurement'),Installer=loader.loadClass('CalciumInstaller'),Settings=loader.loadClass('CalciumAnalysisSettings')
def root=new File('build/test-output/call-checks').absoluteFile;root.deleteDir();root.mkdirs()
int W=420,H=320,N=120
def disc={int cx,int cy,int r->new OvalRoi(cx-r,cy-r,2*r,2*r)}

// ---- distances and zones, against a plain search over every pixel
def rnd=new Random(3);int cw=37,ch=29;boolean[] dots=new boolean[cw*ch];(0..<40).each{dots[rnd.nextInt(cw*ch)]=true}
int[] found=Rings.distances(dots,cw,ch)
for(int y=0;y<ch;y++)for(int x=0;x<cw;x++){int best=Integer.MAX_VALUE;for(int v=0;v<ch;v++)for(int u=0;u<cw;u++)if(dots[v*cw+u])best=Math.min(best,(x-u)*(x-u)+(y-v)*(y-v));assert found[y*cw+x]==best}
assert Rings.distances(new boolean[12],4,3)==null : 'nothing to measure a distance to'
assert Rings.widths(0.645d,0.645d,'micron') as List==[1d,10d,22d,2d] && Rings.widths(0.3225d,0.3225d,'µm') as List==[2d,20d,44d,4d] : 'widths are kept in micrometres'
assert Rings.widths(1d,1d,'pixel') as List==[1d,10d,22d,2d] && Rings.widths(0.645d,0.645d,null) as List==[1d,10d,22d,2d] && Rings.widths(0.4d,0.8d,'micron') as List==[1d,10d,22d,2d] : 'no usable pixel size: the pixel values as they are'
assert Rings.widths(Double.NaN,0.645d,'micron') as List==[1d,10d,22d,2d] && Rings.widths(0.645d,Double.NaN,'micron') as List==[1d,10d,22d,2d] && Rings.widths(Double.POSITIVE_INFINITY,1d,'micron') as List==[1d,10d,22d,2d] : 'a pixel size that is not a number is no pixel size'
assert Rings.scaleOf(0.645d,0.645d,'micron')=='micrometres' && Rings.scaleOf(0.645d,0.7d,'micron')=='not square' && Rings.scaleOf(0.645d,Double.NaN,'µm')=='not square' && Rings.scaleOf(Double.NaN,0.645d,'micron')=='no size' && Rings.scaleOf(0.645d,0.645d,'pixel')=='no size'
assert Steps.micrometres('micron') && Steps.micrometres('microns') && Steps.micrometres('um') && Steps.micrometres('\u00B5m') && !Steps.micrometres('pixel') && !Steps.micrometres(null) : 'one test of the unit, for the scale bar and the rings'
assert Rings.NAMES.size()==25 && Rings.NAMES[0]=='near_right' && Rings.NAMES[7]=='near_upper_right' && Rings.NAMES[15]=='far_upper_right' && Rings.NAMES[16]=='core' && Rings.NAMES[18]=='rim_lower_right'
int[] alone=Steps.pixelsOf(disc(120,100,15),W,H),beside=Steps.pixelsOf(disc(154,100,15),W,H)
def nearest={int i,int[] to->int x=i%W,y=i.intdiv(W);double best=Double.MAX_VALUE;for(int k:to){int dx=k%W-x,dy=k.intdiv(W)-y;best=Math.min(best,(double)(dx*dx+dy*dy))};Math.sqrt(best)}
def one=Rings.newInstance([alone] as int[][],W,H,0.645d,0.645d,'micron')
assert one.micrometres && one.nearUsable(0)==8
(0..<8).each{int z->one.zonePx[z].each{int i->double d=nearest(i,alone);assert d>1d&&d<=10d};one.zonePx[8+z].each{int i->double d=nearest(i,alone);assert d>10d&&d<=22d}}
int nearCount=0,farCount=0;Set own=alone as Set
for(int i=0;i<W*H;i++){if(own.contains(i))continue;int x=i%W,y=i.intdiv(W);if(Math.abs(x-120)>40||Math.abs(y-100)>40)continue;double d=nearest(i,alone);if(d>1d&&d<=10d)nearCount++;if(d>10d&&d<=22d)farCount++}
assert (0..<8).sum{one.pixels(0,it)}==nearCount && (8..<16).sum{one.pixels(0,it)}==farCount : 'every ring pixel is in exactly one sector'
assert (16..<25).sum{one.pixels(0,it)}==alone.length && one.pixels(0,16)<=alone.length/2 && one.pixels(0,16)>alone.length/3 : 'core and rim together are the outline; the core is its deeper half'
double cx=alone.collect{it%W}.sum()/alone.length,cy=alone.collect{it.intdiv(W)}.sum()/alone.length
def side={int[] px->double mx=px.collect{it%W}.sum()/px.length-cx,my=px.collect{it.intdiv(W)}.sum()/px.length-cy;[Math.signum(Math.rint(mx/4d)) as int,Math.signum(Math.rint(my/4d)) as int]}
assert side(one.zonePx[0])==[1,0] && side(one.zonePx[1])==[1,1] && side(one.zonePx[2])==[0,1] && side(one.zonePx[4])==[-1,0] && side(one.zonePx[6])==[0,-1] && side(one.zonePx[7])==[1,-1] : 'right, lower right, bottom ... as seen on the screen'
assert side(one.zonePx[17])==[1,0] && side(one.zonePx[19])==[0,1] && side(one.zonePx[8+5])==[-1,-1]
def two=Rings.newInstance([alone,beside] as int[][],W,H,0.645d,0.645d,'micron')
(0..<16).each{int z->two.zonePx[z].each{int i->assert nearest(i,beside)>2d : 'a ring leaves out pixels within 2 px of another outline'}}
assert two.pixels(0,0)<one.pixels(0,0) && two.pixels(0,4)==one.pixels(0,4) && two.pixels(1,4)<two.pixels(1,0)
def corner=Rings.newInstance([Steps.pixelsOf(disc(6,6,8),W,H),Steps.pixelsOf(new Roi(300,200,3,3),W,H)] as int[][],W,H,1d,1d,'pixel')
assert !corner.micrometres && corner.nearUsable(0)<8 && corner.pixels(0,5)==0 && corner.pixels(1,16)<12 && corner.nearUsable(1)==8 : 'a ring cut by the border loses sectors; an outline of 9 pixels has no usable core'
println 'PASS: ring and core zones are the pixels at the stated distances, in compass sectors, away from other outlines'

// ---- two synthetic recordings: the same cells and events; in the second the field jumps 12 px for two frames
Map at=[crossed:[80,80],edge:[200,80],left:[300,80],right:[334,80],single:[80,220],blips:[200,220]]
Roi[] rois=at.values().collect{disc(it[0],it[1],15)} as Roi[]
Measurement.ensureIds(rois)
String[] names=rois.collect{Measurement.idFor(it)} as String[]
Map cellOf=[at.keySet() as List,names as List].transpose().collectEntries{it}
List<int[]> masks=rois.collect{Steps.pixelsOf(it,W,H)}
int[] halo=(Steps.pixelsOf(ij.plugin.RoiEnlarger.enlarge(rois[3],3),W,H) as List).minus(masks[3] as List) as int[]
def spots=(1..40).collect{[20+rnd.nextInt(W-40),20+rnd.nextInt(H-40),3+rnd.nextInt(3),15+rnd.nextInt(20)]}   // fixed dim structures without an outline
def movie={boolean jolt->
    def stack=new ImageStack(W,H)
    (1..N).each{int f->
        short[] px=new short[W*H];def noise=new Random(1000+f);int sx=jolt&&(f==46||f==47)?12:0
        for(int i=0;i<px.length;i++)px[i]=(short)(100+noise.nextInt(7)-3)
        def put={int i,int amount->int x=i%W+sx,y=i.intdiv(W);if(x>=0&&x<W)px[y*W+x]+=(short)amount}
        def blob={double bx,double by,double r,int amount->for(int y=(int)(by-r-1);y<=by+r+1;y++)for(int x=(int)(bx-r-1);x<=bx+r+1;x++)if((x-bx)*(x-bx)+(y-by)*(y-by)<=r*r&&y>=0&&y<H&&x>=0&&x<W)put(y*W+x,amount)}
        spots.each{blob(it[0],it[1],it[2],it[3])}
        masks.each{int[] m->m.each{put(it,50)}}                                                   // every cell is 50 gray levels over the glass at rest
        // a bright object comes from the left, crosses the first cell (brightest in the middle of it) and leaves to the right
        if(f>=50&&f<=54)blob(57,80,5,400);if(f>=55&&f<=62)blob(70+(f-55)*20d/7d,80,5,f in [58,59]?440:400);if(f>=63&&f<=67)blob(103,80,5,400)
        // a bright object rests on the lower right edge of the second cell and goes
        if(f>=60&&f<=90)blob(200+22*Math.cos(Math.PI/4),80+22*Math.sin(Math.PI/4),8,1500)
        // two neighbours respond in place, the right one first and with a glow around it
        if(f>=62)masks[2].each{put(it,40)}
        if(f>=58){masks[3].each{put(it,40)};halo.each{put(it,20)}}
        // rises of one frame, and of two frames: three times each
        if(f in [60,70,80])masks[4].each{put(it,30)}
        if(f in [60,61,70,71,80,81])masks[5].each{put(it,30)}
        stack.addSlice(''+f,new ShortProcessor(W,H,px,null))
    }
    def imp=new ImagePlus(jolt?'jolt':'still',stack);def cal=imp.getCalibration();cal.pixelWidth=0.645d;cal.pixelHeight=0.645d;cal.setUnit('micron');cal.frameInterval=1d
    imp
}
// the analysis folder a run of the engine would leave, from the saved corrected table: F0 = frames 1-40, the rule a rise of
// 10 % and 3 x the baseline SD for 2 samples in a row. Only the columns the checks read are written.
def analysis={File folder,String base,String name->
    def dir=new File(folder,name);dir.mkdirs()
    def table=Csv.read(new File(folder,base+'_gray_values_background_corrected_native.csv'));List cells=table[0].drop(1)
    List callRows=[['Field','Cell','Role','Stimulus','Start_s','End_s','Effective_end_s','Max_gap_s','Status','Reason','Baseline_mean','Local_reference_relative','Threshold_increase','Peak_increase','Peak_s']]
    Map traces=[:]
    cells.eachWithIndex{cell,int c->
        List y=table.drop(1).collect{Double.parseDouble(it[c+1])};double f0=y[0..<40].sum()/40d
        List n=y.collect{(it-f0)/f0};double mean=n[0..<40].sum()/40d,sd=Math.sqrt(n[0..<40].collect{(it-mean)*(it-mean)}.sum()/39d),need=Math.max(0.1d,3d*sd)
        List rise=n.collect{(it-mean)/(1d+mean)};int peak=(40..<N).max{rise[it]};boolean responded=(40..<N-1).any{rise[it]>=need&&rise[it+1]>=need}
        callRows<<[base,cell,'Stimulus','Stimulus',40,N-1,N-1,1.5,responded?'Responded':'Did not respond','',mean,1d+mean,need,rise[peak],peak];traces[cell]=n
    }
    Csv.write(new File(dir,'responder_calls.csv'),callRows)
    Csv.write(new File(dir,'normalized_traces.csv'),[['Time']+cells]+(0..<N).collect{int f->[f]+cells.collect{traces[it][f]}})
    Csv.write(new File(dir,'mean_trace.csv'),[['Time','Mean','SEM','N_Cells']]+(0..<N).collect{int f->[f,cells.collect{traces[it][f]}.sum()/cells.size(),0,cells.size()]})
    Csv.write(new File(dir,'analysis_manifest.csv'),[['key','value'],['baseline_start','1'],['baseline_end','40'],['sampling_hz','1'],['measurement_input_sha256',Installer.digest(new File(folder,base+'_gray_values_background_corrected_native.csv').bytes)]])
    new File(dir,'COMPLETE').text='complete\n'
    dir
}
// Measure saves two PNG pictures through Fiji's command table, which is built with the menus and so not without a display.
// The one command is entered by hand; nothing else of Measure needs the menus.
def commands=ij.Menus.getDeclaredField('pluginsTable');commands.accessible=true
if(commands.get(null)==null)commands.set(null,new Hashtable(['PNG...':'ij.plugin.PNG_Writer']))
// What Step 1 leaves for Measure, without its two filters that need Fiji's menus: the all-frames picture here is the mean of
// the frames minus the level of the glass.
def prepared={ImagePlus imp->
    int pw=imp.getWidth(),ph=imp.getHeight()
    def steps=Steps.newInstance();steps.movie=imp;steps.n=N;steps.w=pw;steps.h=ph
    steps.say('Recording: '+imp.getTitle());steps.motionShifts=Motion.shifts(imp.getStack(),N)
    float[] mean=View.mean(imp.getStack(),1,N);steps.rest=mean.clone()
    float[] level=new float[pw*ph];for(int i=0;i<level.length;i++)level[i]=mean[i]-100f
    steps.allPic=new ImagePlus('all_frames',new FloatProcessor(pw,ph,level));steps.flat=level.clone()
    def mask=new ByteProcessor(pw,ph);byte[] m=(byte[])mask.getPixels();for(int i=0;i<m.length;i++)if(level[i]>Steps.BRIGHT_ABOVE)m[i]=(byte)255
    new ij.plugin.filter.RankFilters().rank(mask,(double)Steps.BRIGHT_GROW,ij.plugin.filter.RankFilters.MAX)
    boolean[] bright=new boolean[pw*ph];for(int i=0;i<m.length;i++)bright[i]=m[i]!=0
    steps.bright=bright
    steps
}
def measure={String base,boolean jolt->
    def folder=new File(root,base);folder.mkdirs()
    def imp=movie(jolt);def steps=prepared(imp)
    String prefix=new File(folder,base+'_').path
    steps.measureAndSave(rois,names,prefix)
    [folder:folder,recording:base,frames:N,hz:1d,baselineStart:1,baselineEnd:40,key:base,measurementId:base,csv:new File(prefix+'gray_values_background_corrected_native.csv'),
        picture:new ImagePlus('picture',new ColorProcessor(W,H)),movie:imp,names:names as List,rois:rois.collect{(Roi)it.clone()} as Roi[],steps:steps,prefix:prefix]
}
def still=measure('still',false)
File values=new File(still.prefix+'ring_and_core_gray_values.csv.gz'),pixels=new File(still.prefix+'ring_and_core_pixels.csv'),movement=new File(still.prefix+'field_movement.csv'),done=new File(still.prefix+'what_was_done.txt')
assert values.isFile() && pixels.isFile() && movement.isFile() && !new File(values.path+'.part').exists() : 'Measure saves what it measured around the cells beside the tables'
def pixelRows=Csv.read(pixels)
assert pixelRows[0]==Checks.PIXEL_HEADER && pixelRows[0][0..2]==['Cell','ROI_sha256','Area_pixels'] && pixelRows[0][3..27]==Rings.NAMES && pixelRows.drop(1)*.getAt(0)==names as List
assert pixelRows[0][28..-1]==['Near_ring_sectors_usable','Raw_table_sha256','Ring_and_core_values_sha256','Pixel_size_um','Near_ring_from_px','Near_ring_to_px','Far_ring_to_px','Other_outline_margin_px']
def registered=Measurement.registry(rois,rois,new int[rois.length],still.steps,'test')
pixelRows.drop(1).eachWithIndex{row,int c->
    assert row[1]==registered[c+1][5] && (row[2] as int)==rois[c].getContainedPoints().length : 'the outline is named as the ROI registry names it'
    assert row[pixelRows[0].indexOf('Raw_table_sha256')]==Installer.digest(new File(still.prefix+'gray_values_raw.csv').bytes) && row[pixelRows[0].indexOf('Ring_and_core_values_sha256')]==Installer.digest(values.bytes) : 'the files say which raw table they belong to'
    assert row[pixelRows[0].indexOf('Pixel_size_um')]=='0.645' && row[-4..-1]==['1.0','10.0','22.0','2.0']
}
List head=new java.util.zip.GZIPInputStream(new FileInputStream(values)).withCloseable{new BufferedReader(new InputStreamReader(it,'UTF-8')).readLines()}
assert head.size()==N+1 && head[0].split(',').size()==1+25*names.length && head[0].startsWith('Frame,'+names[0]+'|near_right,'+names[0]+'|near_lower_right,') && head[0].endsWith(names[-1]+'|rim_upper_right')
assert head[1].split(',',-1).size()==1+25*names.length && head[1].split(',',-1)[1]==~/\d+\.\d\d/ && head[N].startsWith(N+',')
float[][] back=Rings.read(values,names as List,N)
def ringsOf=Rings.newInstance(masks as int[][],W,H,0.645d,0.645d,'micron');def firstFrame=still.movie.getStack().getProcessor(1)
(0..<25*names.length).each{int k->int[] px=ringsOf.zonePx[k];if(px.length<12)assert Float.isNaN(back[k][0]) else assert Math.abs(back[k][0]-px.collect{firstFrame.get(it)}.sum()/(double)px.length)<=0.00501d}
def moveRows=Csv.read(movement);assert moveRows[0]==['Frame','Displacement_px','Shift_x_px','Shift_y_px'] && moveRows.size()==N+1 && moveRows.drop(1).every{Double.parseDouble(it[1])<2.5d}
List record=done.readLines('UTF-8');List around=record.findAll{it.startsWith('Around each cell: ')}
assert around.size()==1 && around[0].contains('0.645 to 6.45 µm') && around[0].contains('6.45 to 14.19 µm') && around[0].contains('1 to 10 and 10 to 22 pixels') && around[0].contains('within 1.29 µm of another outline') : around
assert around[0].contains('still_ring_and_core_gray_values.csv.gz') && around[0].contains('still_ring_and_core_pixels.csv') && around[0].contains('still_field_movement.csv') && !around[0].contains('(')
assert record.indexOf(around[0])==record.findIndexOf{it.startsWith('Measured on the original movie')}+1 : 'the line follows what was measured in the cells'
byte[] valuesAtMeasure=values.bytes,pixelsAtMeasure=pixels.bytes,movementAtMeasure=movement.bytes
still.steps.measureAndSave(rois,names,still.prefix)
assert done.readLines('UTF-8').count{it.startsWith('Around each cell: ')}==1 && values.bytes==valuesAtMeasure && pixels.bytes==pixelsAtMeasure : 'measuring again leaves one line and the same files'
// a recording without a pixel size in micrometres: pixel values, and the record says so
def plain=movie(false);plain.getCalibration().setUnit('pixel');plain.getCalibration().pixelWidth=1d;plain.getCalibration().pixelHeight=1d
def plainSteps=prepared(plain);new File(root,'plain').mkdirs();plainSteps.measureAndSave(rois,names,new File(root,'plain/plain_').path)
String plainLine=new File(root,'plain/plain_what_was_done.txt').readLines('UTF-8').find{it.startsWith('Around each cell: ')}
assert plainLine.contains('no pixel size in micrometres') && plainLine.contains('a near ring 1 to 10 pixels and a far ring 10 to 22 pixels') && Csv.read(new File(root,'plain/plain_ring_and_core_pixels.csv'))[1][-5..-1]==['','1.0','10.0','22.0','2.0']
assert new File(root,'plain/plain_ring_and_core_gray_values.csv.gz').bytes==valuesAtMeasure : 'the same pixels are measured when the size in micrometres gives the same widths'
// without the field's movement the record does not name its file
String still2=Checks.record(ringsOf,'still_',false)
assert still2.contains('still_ring_and_core_pixels.csv. The movement of the field was not measured.') && !still2.contains('field_movement') && Checks.record(ringsOf,'still_')==around[0]
// micrometre pixels that are not square: said as such, and the pixel values are used
String oblong=Checks.record(Rings.newInstance(masks as int[][],W,H,0.645d,0.7d,'micron'),'still_')
assert oblong.contains('The pixels of the recording are not square, so the widths are the pixel values as they are.') && !oblong.contains('no pixel size') && !oblong.contains('(')
// a table that is not finished leaves nothing behind
def unfinished=new File(root,'unfinished');unfinished.mkdirs();def half=Rings.newInstance(masks as int[][],W,H,0.645d,0.645d,'micron')
half.open(new File(unfinished,'x_ring_and_core_gray_values.csv.gz'),names);half.frame(1,still.movie.getStack().getProcessor(1))
assert unfinished.list().size()==1 && unfinished.list()[0]==~/x_ring_and_core_gray_values\.csv\.gz\.[0-9a-f]{8}\.part/ : unfinished.list()
half.abandon();assert unfinished.list().size()==0
def leftOver=new File(root,'left');leftOver.mkdirs();[Checks.VALUES,Checks.PIXELS,Checks.MOVEMENT].each{new File(leftOver,'x_'+it).text='old'};Checks.removeBeside(new File(leftOver,'x_').path)
assert leftOver.list().size()==0 : 'what an earlier measurement left is taken away, the field movement with it'
// a picture too small to measure the movement on: no movement file, and the record says it was not measured
assert Motion.measurable(W,H,N) && Motion.measurable(196,320,N) && !Motion.measurable(192,320,N) && !Motion.measurable(420,192,N) && !Motion.measurable(W,H,2)
def small=new ImagePlus('small',still.movie.getStack().crop(0,0,0,160,160,N));small.setCalibration(still.movie.getCalibration().copy())
def smallSteps=prepared(small)
assert smallSteps.motionShifts.every{it[0]==0d&&it[1]==0d} : 'Step 1 gives zeros, not a measurement, on a small picture'
Roi[] smallRois=[rois[0]] as Roi[];String[] smallNames=[names[0]] as String[]
def smallFolder=new File(root,'small');smallFolder.mkdirs();smallSteps.measureAndSave(smallRois,smallNames,new File(smallFolder,'small_').path)
String smallLine=new File(smallFolder,'small_what_was_done.txt').readLines('UTF-8').find{it.startsWith('Around each cell: ')}
assert !new File(smallFolder,'small_field_movement.csv').exists() && new File(smallFolder,'small_ring_and_core_pixels.csv').isFile() && smallLine.endsWith('small_ring_and_core_pixels.csv. The movement of the field was not measured. Step 4 uses them to flag a possible passing object. No measurement and no call changes.') : smallLine
println 'PASS: Measure saves the ring and core table, the pixels of each zone with the raw table they belong to, the field movement and one line in the record'

// ---- a failure around the cells never fails Measure: every other file is as without it, and the record says so plainly
def baseFolder=new File(root,'base');baseFolder.mkdirs()
String baseSummary=prepared(movie(false)).measureAndSave(rois,names,new File(baseFolder,'rec_').path)
def others={File folder->folder.listFiles().findAll{it.isFile()&&!(it.name in ['rec_what_was_done.txt','rec_cells_RoiSet.zip'])&&!it.name.contains('ring_and_core')&&!it.name.contains('field_movement')}.collectEntries{[it.name,Installer.digest(it.bytes)]}}
Map wanted=others(baseFolder);assert wanted.size()>=4
[(Checks.VALUES):'rec_ring_and_core_gray_values.csv.gz',(Checks.MOVEMENT):'rec_field_movement.csv',(Checks.PIXELS):'rec_ring_and_core_pixels.csv'].each{String kind,String blocked->
    def broken=new File(root,'broken-'+kind.replace('.','-'));broken.mkdirs()
    new File(broken,blocked).mkdirs();new File(broken,blocked+'/keep').text='x'   // a folder where the file should go
    String summary=prepared(movie(false)).measureAndSave(rois,names,new File(broken,'rec_').path)
    assert summary==baseSummary && others(broken)==wanted : 'Measure goes on and saves the same: '+kind
    List before=new File(baseFolder,'rec_what_was_done.txt').readLines('UTF-8'),after=new File(broken,'rec_what_was_done.txt').readLines('UTF-8')
    assert after.findAll{!it.startsWith('Around each cell: ')}==before.findAll{!it.startsWith('Around each cell: ')} && after.count{it.startsWith('Around each cell: ')}==1
    assert after.find{it.startsWith('Around each cell: ')}=='Around each cell: the rings and the core were not saved: '+blocked+' could not be written. The cell measurements are unchanged.'
    assert broken.listFiles().findAll{it.isFile()&&(it.name.contains('ring_and_core')||it.name.contains('field_movement')||it.name.endsWith('.part'))}==[] : 'nothing half saved is left: '+broken.list()
}
println 'PASS: when the files around the cells cannot be written, Measure saves everything else as it would without them and the record names the file in plain words'

// ---- the checks, read as the page reads them
def Session=loader.loadClass('CalciumAnalysisSession'),Server=loader.loadClass('CalciumWebServer')
def slurp=new JsonSlurper()
def checked={Map data,String name->
    def session=Session.newInstance(data.findAll{!(it.key in ['steps','prefix'])});session.runtime.home=new File(root,'runtime')
    File made=analysis(data.folder,data.recording,name);session.performAnalysis={rt,input,settings,out,job->made}
    def server=Server.newInstance(session:session);server.start()
    def call={String path,Object body=null->def c=(HttpURLConnection)new URL(server.url()+path).openConnection(java.net.Proxy.NO_PROXY);c.connectTimeout=5000;c.readTimeout=30000
        if(body!=null){c.requestMethod='POST';c.doOutput=true;c.setRequestProperty('Content-Type','application/json');c.outputStream.withCloseable{it.write(JsonOutput.toJson(body).getBytes('UTF-8'))}}
        int status=c.responseCode;[status,slurp.parseText(new String((status<400?c.inputStream:c.errorStream).bytes,'UTF-8'))]}
    try {
        def draft=call('api/session')[1].draft
        def started=call('api/analyze',[draft:draft,confirmed:true]);assert started[0]==200 : started
        long end=System.nanoTime()+60000000000L;while(session.busy&&System.nanoTime()<end)Thread.sleep(20)
        def state=call('api/state')[1];assert !state.busy && state.tone=='good' : state.message
        def results=call('api/results')[1]
        // read again from the saved file, as when an analysis is opened without running it
        session.load(made);def reread=slurp.parseText(session.resultsJson).checks
        assert reread.cells==results.checks.cells && reread.file==results.checks.file && reread.ringMeasured==results.checks.ringMeasured && reread.note==results.checks.note : 'the saved checks read back the same, the note with them'
        return results
    } finally {server.stop()}
}
def results=checked(still,'SimpleCa-first')
File saved=new File(still.folder,'still_call_checks_SimpleCa-first.csv')
assert saved.isFile() && new File(still.folder,'still_percent_responding_SimpleCa-first.csv').isFile() && results.checks.file==saved.name && results.checks.ringMeasured && results.checks.note=='' : results.checks.note
def rows=Csv.read(saved);assert rows[0]==Checks.HEADER && rows.drop(1)*.getAt(0)==names as List
def row={String cell->[rows[0],rows.find{it[0]==cellOf[cell]}].transpose().collectEntries{it}}
def answer={String cell->results.checks.cells[cellOf[cell]]}
assert results.calls[cellOf.crossed].Status=='Responded' && results.calls[cellOf.edge].Status=='Responded' && results.calls[cellOf.left].Status=='Responded' && results.calls[cellOf.right].Status=='Responded'
assert results.calls[cellOf.single].Status=='Did not respond' && results.calls[cellOf.blips].Status=='Responded'
// an object crossing an outline: flagged with the frames in which the cell is up, and the sides it was seen on
assert answer('crossed').text=='possible passing object, frames 55 to 62' && answer('crossed').rule=='passes through' && answer('crossed').firstFrame==55 && answer('crossed').lastFrame==62 && answer('crossed').checked : answer('crossed')
assert answer('crossed').detail==~/seen outside at the left at frame 5[0-4], \d+\.\d x the cell's rise, and at the right at frame 6[3-7], \d+\.\d x/ : answer('crossed').detail
assert row('crossed').Flag==answer('crossed').text && row('crossed').Rule=='passes through' && row('crossed').Event_first_frame=='55' && row('crossed').Call=='Responded' && row('crossed').Ring_measured=='yes' && row('crossed').Passing_object_check=='checked'
assert row('crossed').Near_ring_sectors_usable=='8' && (row('crossed').Area_pixels as int)==masks[0].length && answer('crossed').nearSectors==8 && answer('crossed').area==masks[0].length
// an object resting on an edge: never called an object
assert answer('edge').text=='rise on the lower right edge, not in the cell, frames 60 to 90' && answer('edge').rule=='on the edge' && answer('edge').firstFrame==60 && !answer('edge').text.contains('object') : answer('edge')
assert answer('edge').detail==~/core rises -?0\.\d\d x the whole outline, ring \d+\.\d x the rim/ : answer('edge').detail
// a cell that responds in place beside a responding neighbour, and that neighbour: not flagged
['left','right'].each{assert answer(it).text=='' && answer(it).rule=='' && answer(it).checked && answer(it).firstFrame>=58 : answer(it)}
assert answer('left').nearSectors==ringsOf.nearUsable(2) && ringsOf.pixels(2,8)<ringsOf.pixels(2,12)/2 : 'the ring toward the neighbour is mostly left out'
// rises of one frame are counted and not flagged; a responder with three rises of two frames is
assert answer('single').text=='' && answer('single').crossings==3 && answer('single').longestCrossing==1
assert answer('blips').text=='short separate crossings' && answer('blips').rule=='short separate crossings' && answer('blips').crossings==3 && answer('blips').longestCrossing==2 && answer('blips').detail=='3 crossings, the longest 2 samples'
assert answer('crossed').flags==[[text:'possible passing object, frames 55 to 62',rule:'passes through']] && answer('left').flags==[] && answer('blips').flags==[[text:'short separate crossings',rule:'short separate crossings']]
assert Checks.answer([rows:[[Cell:'a',Flag:'rise on the top edge, not in the cell, frames 3 to 9; short separate crossings',Rule:'on the edge; short separate crossings']]]).cells.a.flags*.rule==['on the edge','short separate crossings'] : 'two flags are given one by one too'
assert row('blips').Crossings=='3' && row('blips').Longest_crossing_samples=='2' && row('left').Crossings=='1' && (row('left').Longest_crossing_samples as int)>50
assert results.tally==slurp.parseText(JsonOutput.toJson(Settings.tally(new File(still.folder,'SimpleCa-first'),still))) && results.tally.rows[0].counted==[5,6,6] : 'a flag changes no call and no count'
println 'PASS: a crossing object, an object on an edge, neighbours responding in place and short separate crossings get the stated flags; calls and counts are unchanged'

// ---- a recording measured before this build: the files are measured when the analysis runs, and are the same files
byte[] recordBefore=done.bytes
[values,pixels,movement].each{assert it.delete()}
def again=checked(still,'SimpleCa-second')
assert values.bytes==valuesAtMeasure && pixels.bytes==pixelsAtMeasure && movement.bytes==movementAtMeasure && done.bytes==recordBefore : 'Step 4 measures the same files with the same code and leaves the record alone'
assert again.checks.cells==results.checks.cells
// files of another measurement are not used: with the recording open they are measured again
def jolt=measure('jolt',true)
['ring_and_core_gray_values.csv.gz','ring_and_core_pixels.csv'].each{new File(still.prefix+it).bytes=new File(jolt.prefix+it).bytes}
assert Checks.savedAround(still.prefix,names as List,still.rois,Installer.digest(new File(still.prefix+'gray_values_raw.csv').bytes),N).why=='the ring and core files are for another raw table'
assert checked(still,'SimpleCa-third').checks.cells==results.checks.cells && values.bytes==valuesAtMeasure && pixels.bytes==pixelsAtMeasure
def moved=new Roi[rois.length];rois.eachWithIndex{r,int i->moved[i]=(Roi)r.clone()};moved[0]=disc(81,80,15)
assert Checks.savedAround(still.prefix,names as List,moved,Installer.digest(new File(still.prefix+'gray_values_raw.csv').bytes),N).why=='the ring and core files are for other outlines'
assert Checks.savedAround(still.prefix,(names as List).reverse(),null,Installer.digest(new File(still.prefix+'gray_values_raw.csv').bytes),N).why=='the ring and core files are for other cells'
// without the files and without the recording the check is not made, and every row says so; the crossings are still counted
[values,pixels].each{assert it.delete()}
def closed=checked(still.findAll{it.key!='movie'},'SimpleCa-fourth')
assert !closed.checks.ringMeasured && closed.checks.note=='Not checked for a passing object: ring and core were not measured and the recording is not open.' && !values.exists()
def closedRows=Csv.read(new File(still.folder,'still_call_checks_SimpleCa-fourth.csv')).drop(1)
assert closedRows.every{it[rows[0].indexOf('Passing_object_check')]=='not checked'&&it[rows[0].indexOf('Ring_measured')]=='no'&&it[rows[0].indexOf('Not_checked_because')]=='ring and core were not measured and the recording is not open'}
assert closed.checks.cells[cellOf.crossed].text=='' && closed.checks.cells[cellOf.crossed].firstFrame==55 && !closed.checks.cells[cellOf.crossed].checked && closed.checks.cells[cellOf.blips].text=='short separate crossings'
assert closed.checks.cells[cellOf.crossed].why=='ring and core were not measured and the recording is not open' && closed.checks.cells[cellOf.crossed].detail=='' && closed.checks.cells[cellOf.blips].detail=='3 crossings, the longest 2 samples' : 'the detail holds the flags only; why apart'
// another recording open under the same name is not measured
def other=movie(false);other.getStack().getProcessor(1).add(7)
def wrong=checked(still+[movie:other],'SimpleCa-fifth')
assert !wrong.checks.ringMeasured && wrong.checks.note.contains('the open recording does not give the saved gray values') && !values.exists()
println 'PASS: missing or foreign ring files are measured again when the recording is open, and the check is not made when it is not'

// ---- what the table says when one cell, or the whole analysis, cannot be checked
File log=new File(root,'analysis.log')
def sixth=checked(still,'SimpleCa-sixth');assert sixth.checks.ringMeasured && values.bytes==valuesAtMeasure && movement.bytes==movementAtMeasure
File dir6=new File(still.folder,'SimpleCa-sixth')
def runRows={Map answered->answered.rows.collectEntries{[it.Cell,it]}}
// the field movement alone is missing: it is not worked out again, and the note says what that means
assert movement.delete()
Map noMove=Checks.run(dir6,still,log);Map noMoveRows=runRows(noMove)
assert noMove.note=='The movement of the field was not saved, so a moving field is not told from an object.' && !movement.exists() && noMoveRows.values().every{it.Field_movement_measured=='no'&&it.Ring_measured=='yes'}
assert noMoveRows[cellOf.crossed].Flag=='possible passing object, frames 55 to 62' && !Checks.answer(noMove).cells[cellOf.crossed].movementMeasured && Checks.saved(dir6,still).note==noMove.note
movement.bytes=movementAtMeasure
assert Checks.movement(still.prefix,N)!=null && Checks.movement(still.prefix,N-1)==null
movement.text=movement.text.replaceFirst(/\n2,[^,]*,[^,]*,[^\n]*/,'\n2,NaN,NaN,0.00');assert Checks.movement(still.prefix,N)==null : 'a shift that is not a number: no movement'
movement.bytes=movementAtMeasure
// a cell without a usable ring sector is not checked, and says why; the analysis note is not changed by it
byte[] pixelsKept=pixels.bytes;def edited=Csv.read(pixels);int singleAt=edited.findIndexOf{it[0]==cellOf.single};(3..18).each{edited[singleAt][it]='0'};Csv.write(pixels,edited)
Map noRing=Checks.run(dir6,still,log);Map noRingRow=runRows(noRing)[cellOf.single]
assert noRingRow.Passing_object_check=='not checked' && noRingRow.Not_checked_because=='no ring around the cell' && noRingRow.Crossings==3 && noRing.note=='' && Checks.saved(dir6,still).note=='' : noRingRow
assert runRows(noRing)[cellOf.crossed].Passing_object_check=='checked'
pixels.bytes=pixelsKept
// a scored cell without a peak time: its own reason, and its crossings still counted
File dir7=new File(still.folder,'SimpleCa-seventh');dir7.mkdirs();dir6.listFiles().each{new File(dir7,it.name).bytes=it.bytes}
def calls7=Csv.read(new File(dir7,'responder_calls.csv'));int peakAt=calls7[0].indexOf('Peak_s');calls7.find{it[1]==cellOf.blips}[peakAt]='';Csv.write(new File(dir7,'responder_calls.csv'),calls7)
Map noPeak=runRows(Checks.run(dir7,still,log))[cellOf.blips]
assert noPeak.Not_checked_because=='no peak time' && noPeak.Passing_object_check=='not checked' && noPeak.Flag=='short separate crossings' && noPeak.Crossings==3 && noPeak.Event_first_frame=='' : noPeak
// an analysis that cannot be checked: a plain reason, never a Java message; the details go to the log
log.delete()
def calls8=Csv.read(new File(dir6,'responder_calls.csv'));Csv.write(new File(dir7,'responder_calls.csv'),calls8.collect{it.withIndex().findAll{v,int i->i!=calls8[0].indexOf('Status')}*.getAt(0)})
assert Checks.run(dir7,still,log).note=='The calls were not checked: the calls table has no Cell or Status column.'
new File(dir7,'responder_calls.csv').bytes=new File(dir6,'responder_calls.csv').bytes
new File(dir7,'analysis_manifest.csv').text=new File(dir6,'analysis_manifest.csv').text.replace('baseline_end,40','baseline_end,NaN')
assert Checks.run(dir7,still,log).note=='The calls were not checked: the saved baseline frames do not match the recording.' : 'a baseline that is not a number is refused, not taken as frames'
new File(dir7,'analysis_manifest.csv').bytes=new File(dir6,'analysis_manifest.csv').bytes
new File(dir7,'normalized_traces.csv').text='Time,a\n0\n'
Map unread=Checks.run(dir7,still,log)
assert unread.note=='The calls were not checked: the saved analysis could not be read.' && unread.rows==[] && log.text.startsWith('Checks on the calls, the saved analysis: java.io.IOException') : log.text
List traceLines=new File(dir6,'normalized_traces.csv').readLines('UTF-8')
new File(dir7,'normalized_traces.csv').text=traceLines.take(3).join('\n')+'\n'
assert Checks.run(dir7,still,log).note=='The calls were not checked: the analysis has fewer than 3 frames.'
new File(dir7,'normalized_traces.csv').text=([traceLines[0]]+traceLines.drop(1).collect{it.replaceFirst(/^[^,]*/,'5')}).join('\n')+'\n'
assert Checks.run(dir7,still,log).note=='The calls were not checked: the saved traces have no usable times.'
// measuring around the cells fails at Step 4: said plainly, the exception in the log, nothing half saved
log.delete();[values,pixels].each{it.delete()};values.mkdirs();new File(values,'keep').text='x'
Map unsaved=Checks.run(dir6,still,log)
assert unsaved.note=='Not checked for a passing object: ring and core were not measured and could not be saved now.' && !pixels.exists() && log.text.startsWith('Checks on the calls, around each cell: java.nio.file') : unsaved.note
assert !still.folder.list().any{it.endsWith('.part')}
new File(values,'keep').delete();values.delete()
println 'PASS: a cell without a ring or a peak time, missing movement, an unreadable analysis and a failed measurement are each said in plain words'

// ---- the field moved near the response: said instead of a passing object
def joltRows=Csv.read(new File(jolt.prefix+'field_movement.csv')).drop(1)
assert joltRows.findAll{Double.parseDouble(it[1])>=8d}*.getAt(0)==['46','47'] : joltRows.findAll{Double.parseDouble(it[1])>=2.5d}
def shaken=checked(jolt,'SimpleCa-jolt').checks.cells
assert shaken[cellOf.crossed].text=='field moved, frames 46 to 47' && shaken[cellOf.crossed].rule=='field moved' && shaken[cellOf.crossed].firstFrame==55 : shaken[cellOf.crossed]
assert shaken[cellOf.crossed].detail==~/the field was up to 1\d\.\d px from where it was over the F0 frames, at frames 46 to 47; otherwise possible passing object, frames 55 to 62: seen outside at the left .*/ : shaken[cellOf.crossed].detail
assert shaken[cellOf.edge].text=='rise on the lower right edge, not in the cell, frames 60 to 90' : 'a movement more than the span before the event leaves the flag: '+shaken[cellOf.edge]
assert shaken[cellOf.left].text=='' && shaken[cellOf.blips].text=='short separate crossings'
println 'PASS: a field that moved within the span of the event is named instead of a passing object'

// ---- the two rules on plain numbers
int n=100;double[] flat=new double[n],ground=new double[n];Arrays.fill(flat,150d);Arrays.fill(ground,100d)
double[] up=flat.clone();(50..59).each{up[it]=190d}
float[][] none=new float[25][];(0..<25).each{int z->none[z]=new float[n];(0..<n).each{int f->none[z][f]=(float)(120d+((f*7+z*3)%5-2)*0.5d)}}   // every zone flat, with a little noise
assert Passing.check(flat,ground,none,0,19,55,10)==null : 'a cell that does not rise at its peak has no event'
Map quiet=Passing.check(up,ground,none,0,19,55,10);assert quiet.first==51 && quiet.last==60 && quiet.rule=='' && quiet.text==''
assert Passing.check(up,ground,null,0,19,55,10).first==51 : 'the event is found without the ring'
def with={Map changes->float[][] z=new float[25][];(0..<25).each{z[it]=none[it].clone()};changes.each{int zone,Map frames->frames.each{f,v->z[zone][f as int]=v as float}};z}
def entering=(44..47).collectEntries{[it,220f]},leaving=(63..66).collectEntries{[it,220f]}
assert Passing.check(up,ground,with([(4):entering,(0):leaving]),0,19,55,10).text=='possible passing object, frames 51 to 60'
assert Passing.check(up,ground,with([(4):entering,(8+1):leaving]),0,19,55,10).detail.startsWith('seen outside at the left at frame 4') : 'near and far rings both count'
assert Passing.check(up,ground,with([(4):entering,(5):leaving]),0,19,55,10).rule=='' : 'the two sides must be 90 degrees apart or more'
assert Passing.check(up,ground,with([(4):entering]),0,19,55,10).rule=='' && Passing.check(up,ground,with([(0):leaving]),0,19,55,10).rule=='' : 'one side alone is not a passage'
assert Passing.check(up,ground,with([(4):(44..47).collectEntries{[it,150f]},(0):leaving]),0,19,55,10).rule=='' : 'a rise outside smaller than the cell\'s own does not count'
assert Passing.check(up,ground,with([(4):(44..60).collectEntries{[it,220f]},(0):leaving]),0,19,55,10).rule=='' : 'it must be gone while the cell is up'
assert Passing.check(up,ground,with([(4):(30..33).collectEntries{[it,220f]},(0):leaving]),0,19,55,10).rule=='' : 'more than the span before the event is not seen'
assert Passing.check(up,ground,with([(4):entering,(0):leaving]),0,19,55,5).rule=='passes through' && Passing.check(up,ground,with([(4):entering,(0):leaving]),0,19,55,3).rule=='' : 'the span is a number of frames'
float[][] gone=with([(4):entering,(0):leaving]);gone[4]=null
assert Passing.check(up,ground,gone,0,19,55,10).rule=='' : 'a sector with too few pixels is not used, and no other sector stands in for it'
def during=(50..59).collectEntries{[it,200f]}
float[][] onEdge=with([(17+1):during,(1):(50..59).collectEntries{[it,260f]}])
assert Passing.check(up,ground,onEdge,0,19,55,10).text=='rise on the lower right edge, not in the cell, frames 51 to 60'
float[][] inCell=with([(17+1):during,(1):(50..59).collectEntries{[it,260f]},(16):during]);assert Passing.check(up,ground,inCell,0,19,55,10).rule=='' : 'the core rises with the cell: the rise is in the cell'
float[][] lowRing=with([(17+1):during,(1):(50..59).collectEntries{[it,150f]}]);assert Passing.check(up,ground,lowRing,0,19,55,10).rule=='' : 'the ring must rise at least as much as the rim'
float[][] noCore=with([(17+1):during,(1):(50..59).collectEntries{[it,260f]}]);noCore[16]=null;assert Passing.check(up,ground,noCore,0,19,55,10).rule=='' : 'no core, no edge rule'
double[] dipped=up.clone();(53..55).each{dipped[it]=150d};dipped[60]=190d
assert Passing.check(dipped,ground,none,0,19,57,10).first==51 && Passing.check(dipped,ground,none,0,19,57,10).last==61 : 'a dip of up to three frames does not end the event'
double[] apart=flat.clone();((50..52)+(60..69)).each{apart[it]=190d};assert Passing.check(apart,ground,none,0,19,65,10).first==61 : 'a longer dip does'
double[] time=(0..<20).collect{it as double} as double[]
double[] y=[0,0,0,0,0,0.5,0.5,0,0.5,0,0.5,0.5,0.5,0,0,0.5,0,0,0,0.5] as double[]
assert Passing.crossings(time,y,0d,1d,0.4d,5d,19d,1.5d) as List==[5,3] && Passing.crossings(time,y,0d,1d,0.4d,6d,12d,1.5d) as List==[3,3] && Passing.crossings(time,y,0d,1d,0.6d,0d,19d,1.5d) as List==[0,0]
assert Passing.crossings(time,y,0.1d,2d,0.2d,5d,19d,1.5d) as List==[5,3] && Passing.crossings(time,y,0.1d,2d,0.21d,5d,19d,1.5d) as List==[0,0] : 'the rise is counted from the baseline in units of the reference, as the engine does'
double[] gaps=time.clone();(11..<20).each{gaps[it]+=5d};assert Passing.crossings(gaps,y,0d,1d,0.4d,5d,30d,1.5d) as List==[6,2] : 'a gap in time ends a crossing'
double[] holes=y.clone();holes[11]=Double.NaN;assert Passing.crossings(time,holes,0d,1d,0.4d,5d,19d,1.5d) as List==[6,2]
double[] holey=up.clone();holey[30]=Double.NaN
assert Passing.check(holey,ground,none,0,19,55,10)==null && Passing.check(up,(ground.clone() as double[]).with{it[3]=Double.NaN;it},none,0,19,55,10)==null : 'a gray value that is not a number: no event, not a wrong one'
float[][] nanZone=with([(4):entering,(0):leaving]);nanZone[4][10]=Float.NaN
assert Passing.check(up,ground,nanZone,0,19,55,10).rule=='' : 'a zone with a value that is not a number is not used'
float[][] flatZero=new float[25][];(0..<25).each{flatZero[it]=new float[n];Arrays.fill(flatZero[it],120f)};flatZero[17+1]=(0..<n).collect{f->f in 50..59 ? 200f : 120f} as float[];flatZero[1]=(0..<n).collect{f->f in 50..59 ? 260f : 120f} as float[]
Map exact=Passing.check(up,ground,flatZero,0,19,55,10);assert exact.rule=='on the edge' && !exact.detail.contains('NaN') && !exact.detail.contains('Infinity') : exact
assert Passing.frames(35,35)=='frame 35' && Passing.frames(35,36)=='frames 35 to 36'
println 'PASS: both rules, the event and the crossings follow their written definitions on plain numbers'

// ---- the field-movement guard on plain shifts
def shifts={Closure place->(0..<N).collect{int f->place(f+1) as double[]} as double[][]}
// the field sat 9 px to the left over the F0 frames and stays put after them, so the event is where the field sits for most of the recording, but not where it was over the F0 frames
double[][] drift=shifts({int f->f<=40 ? [-9d,0d] : [0d,0d]})
Map fromF0=Motion.movedNear(drift,1,40,45,72)
assert fromF0==[first:45,last:72,runs:[[45,72]],largest:9d,from:'F0'] && Motion.displacement(drift,0d,0d)[44..71].every{it<8d} : 'counted from the recording alone it would not be seen'
assert Checks.movedText(fromF0)=='the field was up to 9.0 px from where it was over the F0 frames, at frames 45 to 72'
// the F0 frames themselves moved: counted from where the field sits for most of the recording
double[][] smeared=shifts({int f->f<=40 ? [(f%2==0 ? 10d : -10d),0d] : (f in 50..52 ? [9d,0d] : [0d,0d])})
Map fromMost=Motion.movedNear(smeared,1,40,45,72);assert fromMost.from=='recording' && fromMost.runs==[[50,52]] && Checks.movedText(fromMost)=='the field was up to 9.0 px from where it sits for most of the recording, at frames 50 to 52'
// a slow slide that never goes 8 px from either place but moves 8 px or more within the frames around the event
double[][] slide=shifts({int f->f<45 ? [0d,0d] : f<=72 ? [-4.5d+9d*(f-45)/27d,0d] : [4.5d,0d]})
Map within=Motion.movedNear(slide,1,40,45,72);assert within.from=='between' && within.runs==[] && within.first==47 && within.last==71 && Checks.movedText(within)=='the field moved 8.3 px between frames 47 and 71' : within
assert Motion.movedNear(shifts({[0d,0d]}),1,40,45,72)==null && Motion.movedNear(drift,1,40,1,40)==null && Motion.movedNear(drift,1,40,80,120).runs==[[80,120]]
assert Checks.runsText([[46,46]])=='frame 46' && Checks.runsText([[46,47],[60,60]])=='frames 46 to 47 and 60' && Checks.runsText([[1,2],[5,5],[9,12]])=='frames 1 to 2, 5 and 9 to 12'
println 'PASS: the field is counted as moved from where it was over the F0 frames, from most of the recording when those moved, or between two frames near the event'
System.exit(0)

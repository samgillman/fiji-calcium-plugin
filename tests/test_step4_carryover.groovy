// Step 4 starts from what was confirmed in Fiji, remembers the settings last used for a recording, counts baseline samples
// correctly, and gives % responding apart for cells chosen blind and cells added after the response was shown.
def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Settings = loader.loadClass('CalciumAnalysisSettings')
def work = new File('build/test-output/step4-carryover').absoluteFile; work.deleteDir(); work.mkdirs()

// the frames handed to Step 4
def steps = new Expando(n: 300, blindFirst: 1, blindLast: 113, found: true, b0: 94, b1: 113, onset: 122)
def handed = Settings.framesFromFiji(steps)
assert handed.first == 94 && handed.last == 113 && handed.note.contains('frames you confirmed in Fiji (1–113)') : 'a long range: the last 20 frames before the stimulus'
steps.blindLast = 21; handed = Settings.framesFromFiji(steps)
assert handed.first == 2 && handed.last == 21 && handed.note.contains('after frame 21') && handed.chosen == [first: 1, last: 21] : 'the last 20 frames; the stimulus follows the last frame; the frames used to choose cells are kept apart'
steps.blindLast = 15; handed = Settings.framesFromFiji(steps); assert handed.first == 1 && handed.last == 15 : 'fewer than 20 frames: all of them'
steps.blindFirst = 0; steps.blindLast = 0; handed = Settings.framesFromFiji(steps)
assert handed.first == 94 && handed.last == 113 && handed.note.contains('near frame 122') : 'no frames confirmed: the frames before the rise'
steps.found = false; handed = Settings.framesFromFiji(steps)
assert handed.first == 1 && handed.last == 20 && handed.note.contains('guess')
println 'PASS: Step 4 takes its baseline and stimulus time from the frames confirmed in Fiji'

// the draft built from them, and the baseline sample count
def data = [frames: 300, hz: 1d, baselineStart: 1, baselineEnd: 21, recording: 'rec', folder: work]
def draft = Settings.draft(data)
assert draft.values.baseline_start == '1' && draft.values.baseline_end == '21'
def window = draft.protocol.windows[0]
assert (window.Stimulus_onset_s as double) == 21d && (window.Baseline_start_s as double) == 0d && (window.Baseline_end_s as double) == 21d && (window.End_s as double) == 299d
def checked = Settings.diagnose(draft, data)
assert !checked.errors : checked.errors
assert !checked.warnings.any { it.contains('samples but at least') } : 'A 21-sample baseline is not short: ' + checked.warnings
def shortBase = Settings.copy(draft); shortBase.protocol.windows[0].Baseline_start_s = 19.5d; shortBase.protocol.windows[0].Baseline_s = 1.5d
def shortChecked = Settings.diagnose(shortBase, data)
assert shortChecked.warnings.any { it.contains('holds 1 samples but at least 3 are needed') } : shortChecked.warnings
println 'PASS: baseline samples are counted (a full baseline no longer warns "0 samples"; a real short one still warns)'

// settings are remembered for the same recording and the same frames from Fiji only
assert Settings.remembered(data) == null
def edited = Settings.copy(draft); edited.protocol.windows[0].End_s = 250d; edited.values.condition = 'ATP'
Settings.remember(data, edited)
def again = Settings.remembered(data)
assert again != null && (again.protocol.windows[0].End_s as double) == 250d && again.values.condition == 'ATP' && again.timestamps == null
assert Settings.remembered(data + [baselineEnd: 113]) == null : 'other frames in Fiji: start again from them'
// the stimulus frame of that analysis is offered for the next recording, and does not unfile the settings remembered for this one
int usual = Settings.lastStimulusFrame(); assert usual >= 1 : 'the last analysis sets the starting stimulus frame'
assert Settings.remembered(data + [baselineStart: Math.max(1, usual - 19), baselineEnd: usual, fijiBaselineStart: data.baselineStart, fijiBaselineEnd: data.baselineEnd]) != null : 'filed under the frames from Fiji'
ij.Prefs.set('calcium.step4.stimulus_after_frame', 0d)
assert Settings.remembered(data + [recording: 'other']) == null && Settings.remembered(data + [frames: 299]) == null
new File(work, 'step4_last_settings.json').text = '{broken'
assert Settings.remembered(data) == null : 'an unreadable file is ignored'
println 'PASS: Step 4 settings are remembered per recording and dropped when the frames from Fiji change'

// % responding, worked out once: all cells first, then by how the cell was chosen, then the quiet cells; each row with every
// cell and without the cells whose corrected fluorescence went below zero
def results = new File(work, 'SimpleCa-test'); results.mkdirs()
new File(results, 'responder_calls.csv').text = '"Field","Cell","Role","Stimulus","Status","Threshold_increase"\n' + [
    ['a', 'Stimulus', 'Responded', 0.12], ['b', 'Stimulus', 'Did not respond', 0.45], ['c', 'Stimulus', 'Excluded', ''], ['d', 'Stimulus', 'Responded', 0.40],
    ['e', 'Stimulus', 'Responded', 0.30], ['f', 'Stimulus', 'Did not respond', 0.20], ['a', 'Control', 'Did not respond', 0.12], ['zz', 'Stimulus', 'Responded', 0.10]
].collect { '"rec","' + it[0] + '","' + it[1] + '","Stimulus","' + it[2] + '",' + it[3] }.join('\n') + '\n'
// cell d dips below -1 at one frame; cell c does too but is excluded, so it is not among the cells left out
new File(results, 'normalized_traces.csv').text = '"Time","a","b","c","d","e","f","zz"\n0,0,0,0,0,0,0,0\n1,0.5,0.1,-2,-1.2,0.6,0.1,0.4\n2,0.4,0,0,0.9,0.5,-0.9,0.3\n'
def tagged = data + [chosenBlind: [a: true, b: true, c: true, d: false, e: false, f: false]]
def found = Settings.tally(results, tagged), rowOf = { String key -> found.rows.find { it.key == key } }
assert found.rows*.key == ['all', 'blind', 'later', 'quiet'] && found.belowZero == ['d'] && found.quietRise == 0.30d : found
assert rowOf('all').counted == [3, 5, 6] && rowOf('all').without == [2, 4, 5] : 'responded, evaluable, registered; without cell d'
assert rowOf('blind').counted == [1, 2, 3] && rowOf('blind').without == [1, 2, 3] && rowOf('later').counted == [2, 3, 3] && rowOf('later').without == [1, 2, 2]
// quiet: required rise at most 30 %, whether the cell responded or not: a, e and f; b and d need more and are set aside
assert rowOf('quiet').counted == [2, 3, 3] && rowOf('quiet').without == [2, 3, 3] && rowOf('quiet').label == 'Quiet enough to detect a 30 % rise'
assert Settings.tally(results, tagged, [a: [0d, 0.5d], d: [0d, -1.2d], c: [-2d]]).belowZero == ['d'] : 'the same from traces already read'
def sentence = Settings.percentResponding(results, tagged)
assert sentence == ' All cells: 3 of 5 responded (60.0 %). Chosen blind: 1 of 2 responded (50.0 %). Added after the response was shown: 2 of 3 responded (66.7 %).' : sentence
def table = new File(work, 'rec_percent_responding_SimpleCa-test.csv').readLines()
assert table[0] == 'cells,registered,evaluable,responded,percent_responding,cells_below_zero,evaluable_without_them,responded_without_them,percent_without_them' : table[0]
assert table.size() == 5 && table[1] == 'All cells,6,5,3,60.0,1,4,2,50.0' && table[2] == 'Chosen blind,3,2,1,50.0,0,2,1,50.0' && table[3] == 'Added after the response was shown,3,3,2,66.7,1,2,1,50.0' && table[4] == 'Quiet enough to detect a 30 % rise,3,3,2,66.7,0,3,2,66.7' : table
assert Settings.percentResponding(results, data + [chosenBlind: [a: true, b: true]]).contains('All 2 cells were chosen blind: 1 of 2 responded (50.0 %)')
assert Settings.percentResponding(results, data) == '' : 'not recorded how cells were chosen: no sentence'
assert new File(work, 'rec_percent_responding_SimpleCa-test.csv').readLines()*.split(',')*.getAt(0) == ['cells', 'All cells', 'Quiet enough to detect a 30 % rise'] : 'the table is saved all the same'
assert Settings.tally(new File(work, 'no-such-analysis'), tagged) == null
println 'PASS: % responding worked out once: all cells, by how chosen, quiet cells; with and without the cells below zero'

// A rule remembered from an earlier build, still on the former default "for 1 s", now counts 2 samples in a row; a length set
// by hand is kept; settings remembered by this build are taken as they are
def formerDefault = { Map change -> Map p = new groovy.json.JsonSlurper().parseText(groovy.json.JsonOutput.toJson(Settings.draft(data).protocol)) as Map; (p.windows as List).each { it.Duration_mode = 'seconds'; it.Min_duration_s = 1d; it.putAll(change) }; p }
def asSaved = { Map protocol, Map extra -> new File(work, 'step4_last_settings.json').text = groovy.json.JsonOutput.toJson([recording: '' + data.recording, frames: data.frames as int, baseline_start: Settings.fijiFrame(data, 'Start'), baseline_end: Settings.fijiFrame(data, 'End'), values: Settings.draft(data).values, protocol: protocol] + extra) }
asSaved(formerDefault([:]), [:]); Map note = [:]; def moved = Settings.remembered(data, note)
assert moved != null && note.ruleMoved && moved.protocol.windows[0].Duration_mode == 'samples' && moved.protocol.windows[0].Consecutive_samples == 2 : moved?.protocol?.windows
asSaved(formerDefault([Min_duration_s: 3d]), [:]); note = [:]; def kept = Settings.remembered(data, note)
assert kept != null && !note.ruleMoved && kept.protocol.windows[0].Duration_mode == 'seconds' && (kept.protocol.windows[0].Min_duration_s as double) == 3d : 'a length set by hand is kept'
asSaved(formerDefault([:]), [format: Settings.REMEMBERED_FORMAT]); note = [:]; def chosen = Settings.remembered(data, note)
assert chosen != null && !note.ruleMoved && chosen.protocol.windows[0].Duration_mode == 'seconds' : 'chosen again in this build: kept'
assert Settings.draft(data).protocol.windows[0].Duration_mode == 'samples' && Settings.draft(data).protocol.windows[0].Consecutive_samples == 2 : 'the default rule'
println 'PASS: a remembered rule on the former default moves to 2 samples in a row once; a chosen length is kept'

// the page's arrangement is kept with Fiji's preferences; only known keys and sane values are taken
assert Settings.setPagePref('step4.side', '420px').kept == 'step4.side' && Settings.pagePrefs()['step4.side'] == '420px'
Settings.setPagePref('step4.split', '34.1%'); Settings.setPagePref('step4.help', '1')
assert Settings.pagePrefs() == ['step4.side': '420px', 'step4.split': '34.1%', 'step4.help': '1']
Settings.setPagePref('step4.side', null); assert !Settings.pagePrefs().containsKey('step4.side')
['step4.side': 'javascript:1', 'step4.split': '150%', 'other': '1'].each { k, v -> try { Settings.setPagePref(k, v); assert false : k } catch (IllegalArgumentException expected) { } }
Settings.setPagePref('step4.split', null); Settings.setPagePref('step4.help', null)
// How the cells are drawn on the picture (plain outlines, numbers) and how often the zoom hint was shown are kept the same way.
['step4.plain': '1', 'step4.numbers': '1', 'step4.viewhint': '3'].each { k, v -> assert Settings.setPagePref(k, v).kept == k && Settings.pagePrefs()[k] == v }
['step4.plain': 'yes', 'step4.numbers': '2', 'step4.viewhint': '12'].each { k, v -> try { Settings.setPagePref(k, v); assert false : k } catch (IllegalArgumentException expected) { } }
['step4.plain', 'step4.numbers', 'step4.viewhint'].each { Settings.setPagePref(it, null) }; assert Settings.pagePrefs().isEmpty()
// The figure editor's settings are kept as one line of JSON for the plot and one for the heatmap; anything else is refused.
String figure = '{"width":6,"height":4,"font":"arial","grid":false,"xTitle":"Time (s)","ink":"#000000"}'
assert Settings.setPagePref('step4.figure.plot', figure).kept == 'step4.figure.plot' && Settings.pagePrefs()['step4.figure.plot'] == figure
Settings.setPagePref('step4.figure.heatmap', '{}'); assert Settings.pagePrefs()['step4.figure.heatmap'] == '{}'
[['step4.figure.plot', 'width=6'], ['step4.figure.heatmap', '{"a":1}\n{'], ['step4.figure.other', '{}'], ['step4.figure.plot', '{' + 'x' * 1600 + '}']].each { pair -> try { Settings.setPagePref(pair[0], pair[1]); assert false : pair } catch (IllegalArgumentException expected) { } }
Settings.setPagePref('step4.figure.plot', null); Settings.setPagePref('step4.figure.heatmap', null); assert !Settings.pagePrefs().keySet().any { (it as String).startsWith('step4.figure') }
println 'PASS: pane sizes and the explanations switch are remembered with Fiji\'s preferences'
// saved pictures: PNG or TIFF with the resolution in the file; TIFF compressed LZW or ZIP, or not at all
def Picture = loader.loadClass('CalciumPictureFile')
def drawn = new java.awt.image.BufferedImage(300, 200, java.awt.image.BufferedImage.TYPE_INT_ARGB)
def pen = drawn.createGraphics(); pen.setColor(java.awt.Color.WHITE); pen.fillRect(0, 0, 300, 200); pen.setColor(java.awt.Color.BLUE); pen.drawLine(0, 0, 299, 199); pen.dispose()
def sentPng = new ByteArrayOutputStream(); javax.imageio.ImageIO.write(drawn, 'png', sentPng); String dataUrl = 'data:image/png;base64,' + Base64.encoder.encodeToString(sentPng.toByteArray())
def tiffInfo = { byte[] bytes -> def reader = javax.imageio.ImageIO.getImageReadersByFormatName('tiff').next(); reader.setInput(javax.imageio.ImageIO.createImageInputStream(new ByteArrayInputStream(bytes)))
    def dir = javax.imageio.plugins.tiff.TIFFDirectory.createFromMetadata(reader.getImageMetadata(0)); def x = dir.getTIFFField(282).getAsRational(0)
    [width: reader.getWidth(0), dpi: (x[0] / x[1]) as int, unit: dir.getTIFFField(296).getAsInt(0), compression: dir.getTIFFField(259).getAsInt(0), pixel: reader.read(0).getRGB(150, 100) & 0xffffff] }
def lzw = Picture.fromPage(dataUrl, 'tiff', 600, 'lzw'), plain = Picture.fromPage(dataUrl, 'TIFF', '300', 'none'), zipped = Picture.fromPage(dataUrl, 'tiff', 300, 'zip')
assert lzw[1] == 'tif' && tiffInfo(lzw[0]) == [width: 300, dpi: 600, unit: 2, compression: 5, pixel: drawn.getRGB(150, 100) & 0xffffff] : 'LZW TIFF at 600 dpi, pixels unchanged'
assert tiffInfo(plain[0]).with { it.dpi == 300 && it.compression == 1 } && tiffInfo(zipped[0]).compression in [8, 32946]
assert (lzw[0] as byte[]).length < (plain[0] as byte[]).length : 'compression makes the file smaller'
def png600 = Picture.fromPage(dataUrl, 'png', 600, null)
def pngReader = javax.imageio.ImageIO.getImageReadersByFormatName('png').next(); pngReader.setInput(javax.imageio.ImageIO.createImageInputStream(new ByteArrayInputStream(png600[0] as byte[])))
assert png600[1] == 'png' && (pngReader.getImageMetadata(0).getAsTree('javax_imageio_png_1.0').getElementsByTagName('pHYs').item(0).getAttribute('pixelsPerUnitXAxis') as int) == 23622 : '600 dpi is 23622 pixels per metre'
assert Picture.fromPage(dataUrl, null, null, null)[0] == sentPng.toByteArray() : 'with no options the picture is kept as drawn'
['gif': 300, 'tiff': 5000].each { f, d -> try { Picture.fromPage(dataUrl, f, d, 'lzw'); assert false : f } catch (IllegalArgumentException expected) { } }
try { Picture.fromPage('data:image/png;base64,' + Base64.encoder.encodeToString('not a picture'.bytes), 'tiff', 300, 'lzw'); assert false } catch (IllegalArgumentException expected) { }
println 'PASS: pictures saved as PNG or TIFF with 300/600 dpi in the file and LZW, ZIP or no compression'
println 'All Step 4 carry-over checks passed.'

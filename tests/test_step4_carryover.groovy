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

// % responding by how the cell was chosen
def results = new File(work, 'SimpleCa-test'); results.mkdirs()
new File(results, 'responder_calls.csv').text = '"Field","Cell","Role","Stimulus","Status"\n' + [
    ['a', 'Stimulus', 'Responded'], ['b', 'Stimulus', 'Did not respond'], ['c', 'Stimulus', 'Excluded'], ['d', 'Stimulus', 'Responded'],
    ['e', 'Stimulus', 'Responded'], ['f', 'Stimulus', 'Did not respond'], ['a', 'Control', 'Did not respond'], ['zz', 'Stimulus', 'Responded']
].collect { '"rec","' + it[0] + '","' + it[1] + '","Stimulus","' + it[2] + '"' }.join('\n') + '\n'
def tagged = data + [chosenBlind: [a: true, b: true, c: true, d: false, e: false, f: false]]
def sentence = Settings.byHowChosen(results, tagged)
assert sentence.contains('Chosen blind: 1 of 2 responded (50.0 %)') && sentence.contains('Added after the response was shown: 2 of 3 responded (66.7 %)') && sentence.contains('All cells: 3 of 5 responded (60.0 %)') : sentence
def table = new File(work, 'rec_responders_by_how_chosen_SimpleCa-test.csv').readLines()
assert table.size() == 4 && table[1].startsWith('blind,3,2,1,50.0') && table[2].contains('after the response was shown,3,3,2,66.7') && table[3].startsWith('all cells,6,5,3,60.0') : table
assert Settings.byHowChosen(results, data + [chosenBlind: [a: true, b: true]]).contains('All 2 cells were chosen blind: 1 of 2 responded (50.0 %)')
assert Settings.byHowChosen(results, data) == '' : 'no tags, no sentence'
println 'PASS: % responding apart for cells chosen blind and cells added after the response was shown'
// the page's arrangement is kept with Fiji's preferences; only known keys and sane values are taken
assert Settings.setPagePref('step4.side', '420px').kept == 'step4.side' && Settings.pagePrefs()['step4.side'] == '420px'
Settings.setPagePref('step4.split', '34.1%'); Settings.setPagePref('step4.help', '1')
assert Settings.pagePrefs() == ['step4.side': '420px', 'step4.split': '34.1%', 'step4.help': '1']
Settings.setPagePref('step4.side', null); assert !Settings.pagePrefs().containsKey('step4.side')
['step4.side': 'javascript:1', 'step4.split': '150%', 'other': '1'].each { k, v -> try { Settings.setPagePref(k, v); assert false : k } catch (IllegalArgumentException expected) { } }
Settings.setPagePref('step4.split', null); Settings.setPagePref('step4.help', null)
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

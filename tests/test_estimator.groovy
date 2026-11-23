// What a cell's noise is taken from. A protocol made here uses the larger of the baseline SD and the whole-trace noise; a
// protocol that was saved or imported keeps the estimator it names; settings remembered for a recording keep an estimator
// they name and take the new one when they name none; each estimator is recorded with the earliest engine release that
// defines it; and the plugin's own numbers read the calls table of either kind.
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def loader = new GroovyClassLoader(getClass().classLoader)
loader.parseClass(new File('scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'))
def Protocol = loader.loadClass('CalciumResponseProtocol'), Settings = loader.loadClass('CalciumAnalysisSettings'), Session = loader.loadClass('CalciumAnalysisSession')
def Engine = loader.loadClass('CalciumSimpleCaEngine')
def work = new File('build/test-output/estimator').absoluteFile; work.deleteDir(); work.mkdirs()
def slurp = new JsonSlurper()
def copy = { Map p -> slurp.parseText(JsonOutput.toJson(p)) as Map }
String larger = 'larger_of_baseline_sd_and_trace_noise'

// The plugin's table is the engine's own: every estimator with the earliest release that lists it, read from the pinned source.
String source = new File('analysis/upstream/response_protocols.R').getText('UTF-8')
def listed = (source =~ /(?s)responder_engine_estimators <- function\(\) list\((.*?)\)\n\n/)
assert listed.find() : 'responder_engine_estimators() is not where the test expects it in the pinned source'
Map earliest = [:]
(listed.group(1) =~ /"(\d+\.\d+\.\d+)" = ([^\n]+)/).each { whole, release, names -> (names =~ /"([a-z_]+)"/).each { all, name -> if (!earliest.containsKey(name)) earliest[name] = release } }
assert earliest == Protocol.ESTIMATOR_ENGINE : 'pinned engine: ' + earliest + ', plugin: ' + Protocol.ESTIMATOR_ENGINE
assert Protocol.ESTIMATOR_NAME.keySet() == Protocol.ESTIMATOR_ENGINE.keySet()
String running = (source =~ /responder_engine_version <- function\(\) "([0-9.]+)"/)[0][1]
assert running == Protocol.ESTIMATOR_ENGINE.values().max() : 'the engine in the plugin is ' + running
def unpacked = Engine.unpack(work)
['responders.R', 'response_protocols.R', 'utils.R', 'LICENSE'].each { assert new File(unpacked, 'upstream/' + it).bytes == new File('analysis/upstream/' + it).bytes : it + ': run tools/embed_analysis.py' }
assert new File(unpacked, 'runner.R').bytes == new File('analysis/runner.R').bytes
println 'PASS: the estimators and their engine versions are the pinned engine\'s own, and the plugin carries that engine'

// A protocol made here. One function decides its estimator.
assert Protocol.newEstimator() == larger : 'new analyses take noise from the whole trace when it is larger than the baseline SD'
def window = { -> Protocol.window('Stimulus', 'ATP', 0d, 20d, 299d) }
Map made = Protocol.payload([window()], 2, 'Made here', 'provisional')
assert made.sd_estimator == larger && made.engine_version == '2.1.0' && made.windows*.SD_estimator == [larger] && made.windows*.Engine_version == ['2.1.0']
def data = [frames: 300, hz: 1d, baselineStart: 1, baselineEnd: 20, recording: 'rec', folder: work]
Map startingDraft = Settings.draft(data)
assert startingDraft.protocol.sd_estimator == larger && startingDraft.protocol.engine_version == '2.1.0' && startingDraft.protocol.windows[0].Consecutive_samples == 2
Map legacy = Protocol.payload([window()], 1, 'Historical', 'user_defined')
assert legacy.engine_version == 'legacy-1' && !legacy.containsKey('sd_estimator') && !legacy.windows[0].containsKey('SD_estimator') : 'the legacy rule has one noise only'
try { Protocol.payload([window()], 2, 'Unknown', 'provisional', null, 'median_sd'); assert false } catch (IllegalArgumentException expected) { assert expected.message.contains('median_sd') }
println 'PASS: a protocol made here uses the larger of baseline SD and whole-trace noise, recorded with engine 2.1.0'

// A protocol that is read keeps its estimator, however it is written again.
Map older = Protocol.payload([window()], 2, 'Saved earlier', 'user_defined', null, 'sample_sd')
assert older.sd_estimator == 'sample_sd' && older.engine_version == '2.0.0' && older.windows*.Engine_version == ['2.0.0']
[older, made].each { Map p -> Map again = Protocol.canonical(copy(p))
    assert again.sd_estimator == p.sd_estimator && again.engine_version == p.engine_version && again.windows*.SD_estimator == p.windows*.SD_estimator && again.windows*.Engine_version == p.windows*.Engine_version }
// SimpleCa's own files name the estimator in the windows only
def golden = new File('tests/fixtures/response-protocols/golden-protocol.json')
assert !(slurp.parse(golden, 'UTF-8') as Map).containsKey('sd_estimator') && Protocol.read(golden).sd_estimator == 'sample_sd'
Map goldenAgain = Protocol.canonical(golden.getText('UTF-8'))
assert goldenAgain.sd_estimator == 'sample_sd' && goldenAgain.engine_version == '2.0.0' && goldenAgain.windows*.SD_estimator == ['sample_sd', 'sample_sd']
def asSimpleCa = { Map p, Closure change = { } -> Map q = copy(p); q.remove('sd_estimator'); change.call(q); q }
assert Protocol.canonical(asSimpleCa(made)).sd_estimator == larger && Protocol.canonical(asSimpleCa(made)).engine_version == '2.1.0'
// the engine refuses a rule recorded with another release than the earliest that defines it; so does the plugin, before a run
def refused = { Map p -> try { Protocol.canonical(p); null } catch (IOException e) { e.message } }
assert refused(asSimpleCa(made) { it.engine_version = '2.0.0'; it.windows.each { w -> w.Engine_version = '2.0.0' } }).contains('is recorded with engine version 2.1.0')
assert refused(asSimpleCa(older) { it.engine_version = '2.1.0'; it.windows.each { w -> w.Engine_version = '2.1.0' } }).contains('is recorded with engine version 2.0.0')
assert refused(asSimpleCa(made) { it.windows[0].SD_estimator = 'median_sd' }) == 'Unsupported or missing SD_estimator.'
assert refused(asSimpleCa(made) { it.windows[0].remove('SD_estimator') }) == 'Unsupported or missing SD_estimator.'
assert refused(asSimpleCa(made) { it.engine_version = '3.0.0'; it.windows.each { w -> w.Engine_version = '3.0.0' } }) == 'Unsupported protocol engine version.'
Map two = Protocol.payload([window() + [End_s: 100d], Protocol.window('Control', 'KCl', 100d, 150d, 299d)], 2, 'Two windows', 'provisional')
assert two.windows*.SD_estimator == [larger, larger]
assert refused(asSimpleCa(two) { it.windows[1].SD_estimator = 'sample_sd' }) == 'All windows must share one SD_estimator.'
Map contradiction = copy(made); contradiction.sd_estimator = 'sample_sd'
assert refused(contradiction).contains('disagree') : 'the name at the top and the copies in the windows must agree'
Map unknown = copy(made); unknown.sd_estimator = 'median_sd'
assert refused(unknown) == 'Unsupported or missing SD_estimator.'
println 'PASS: a saved or imported protocol keeps its estimator; an estimator recorded with the wrong engine version is refused'

// The Step 4 page writes "2.0.0" and "sample_sd" into every version 2 protocol it sends and leaves sd_estimator as it was
// given (web/step4.js, collectDraft). Where a draft from the page comes in, the estimator at the top stands. A protocol
// file is read strictly: the same text in a file contradicts itself and is refused.
String pageSource = new File('web/step4.js').getText('UTF-8')
assert pageSource.contains("p.engine_version = legacy ? 'legacy-1' : '" + Protocol.PAGE_ENGINE + "'") && pageSource.contains("SD_estimator: '" + Protocol.PAGE_ESTIMATOR + "'") :
    'the page no longer writes these fixed copies: drop PAGE_ENGINE, PAGE_ESTIMATOR, the page argument of estimator, forPage and exchange'
def fromPage = { Map p -> Map q = copy(p); if (q.schema_version == 2) { q.engine_version = '2.0.0'; q.windows.each { it.Engine_version = '2.0.0'; it.SD_estimator = 'sample_sd' } }; q }
Map sent = Protocol.canonical(fromPage(made), true)
assert sent.sd_estimator == larger && sent.engine_version == '2.1.0' && sent.windows*.SD_estimator == [larger] && sent.windows*.Engine_version == ['2.1.0']
assert Protocol.canonical(fromPage(older), true).sd_estimator == 'sample_sd'
assert refused(fromPage(made)).contains('disagree') : 'a file that names one estimator at the top and another in its windows is refused'
File selfContradicting = new File(work, 'names-two-estimators.json'); selfContradicting.setText(JsonOutput.toJson(fromPage(made)), 'UTF-8')
try { Protocol.read(selfContradicting); assert false : 'Load protocol must refuse it' } catch (IOException expected) { assert expected.message.contains('disagree') }
[false, true].each { boolean page -> try { Protocol.canonical(copy(contradiction), page); assert false } catch (IOException expected) { assert expected.message.contains('disagree') : 'the other contradiction is refused at every door' } }
Map keptByEarlierBuild = fromPage(older); keptByEarlierBuild.remove('sd_estimator')
[false, true].each { boolean page -> assert Protocol.canonical(copy(keptByEarlierBuild), page).sd_estimator == 'sample_sd' && Protocol.canonical(copy(keptByEarlierBuild), page).engine_version == '2.0.0' }
assert Protocol.draftEstimator(fromPage(made)) == larger && Protocol.draftEstimator(legacy) == null && Protocol.draftEstimator([schema_version: 2, windows: 'unfinished']) == null && Protocol.draftEstimator(null) == null
// What the page is handed is in the words it sends back, so an untouched draft equals what the page collects from it.
[made, older, two].each { Map p -> Map handed = Protocol.forPage(copy(p))
    assert JsonOutput.toJson(handed) == JsonOutput.toJson(fromPage(p)) && Protocol.canonical(handed, true) == Protocol.canonical(copy(p)) }
Map unfinished = copy(made); unfinished.windows[0].remove('SD_estimator')
assert Protocol.forPage(copy(unfinished)) == unfinished : 'an unfinished draft goes back as it came'
// A legacy rule has one noise, so its file names none. The page is handed it with the estimator of a new protocol at the
// top and carries that over when the rule type is changed: the version 2 rule made there is a protocol made here.
Map legacyForPage = Protocol.forPage(copy(legacy))
assert legacyForPage.sd_estimator == Protocol.newEstimator() && legacyForPage.findAll { k, v -> k != 'sd_estimator' } == copy(legacy)
assert !Protocol.canonical(legacyForPage, true).containsKey('sd_estimator') && Protocol.noiseOf(Protocol.canonical(legacyForPage, true)) == 'sample_sd' && Protocol.noiseOf(made) == larger
Map madeOnPageFromLegacy = fromPage(made); madeOnPageFromLegacy.sd_estimator = legacyForPage.sd_estimator
assert Protocol.canonical(madeOnPageFromLegacy, true).engine_version == '2.1.0'
// Import protocol JSON and Export protocol are one door (api/protocol): a file comes in or the page's draft does, and the answer is the file.
assert Protocol.exchange(fromPage(made)) == Protocol.canonical(copy(made)) && Protocol.exchange(copy(made)) == Protocol.canonical(copy(made)) : 'the page\'s draft is exported with the estimator it names'
assert Protocol.exchange(asSimpleCa(older)).sd_estimator == 'sample_sd' && Protocol.exchange(asSimpleCa(made)).sd_estimator == larger : 'an imported file keeps its estimator'
assert Protocol.exchange(copy(legacy)).sd_estimator == Protocol.newEstimator() : 'an imported legacy file is answered with the name the page carries'
assert Protocol.exchange(legacyForPage) == Protocol.canonical(copy(legacy)) && !Protocol.exchange(legacyForPage).containsKey('sd_estimator') : 'the page\'s legacy draft is exported as the legacy file, which names no estimator'
try { Protocol.exchange([schema: 'simpleca.response-protocol', schema_version: 99]); assert false } catch (IOException expected) { }
// what is handed to the engine names the estimator, and so does settings.tsv, for the legacy rule too
Map checked = Settings.check(startingDraft.values, fromPage(startingDraft.protocol), true, data, [:], true)
assert checked.sd_estimator == larger && (slurp.parseText(checked.protocol_json) as Map).engine_version == '2.1.0'
try { Settings.check(startingDraft.values, fromPage(startingDraft.protocol), true, data); assert false : 'the Fiji window checks its protocol as a file is checked' } catch (IOException expected) { assert expected.message.contains('disagree') }
assert Settings.check(startingDraft.values, keptByEarlierBuild, true, data).sd_estimator == 'sample_sd'
assert Settings.check(startingDraft.values, legacy, true, data).sd_estimator == 'sample_sd' && Settings.check(startingDraft.values, legacyForPage, true, data, [:], true).sd_estimator == 'sample_sd'
assert !Settings.diagnose([values: startingDraft.values, protocol: fromPage(made)], data).errors && !Settings.diagnose([values: startingDraft.values, protocol: legacyForPage], data).errors : 'drafts are checked as the page sends them'
println 'PASS: the page\'s fixed "sample_sd" and "2.0.0" do not undo the estimator named at the top of its draft; a file that names two estimators is refused; a legacy rule is handed to the page with the estimator of a new protocol'

// A new session starts with the new estimator, and the page is handed its draft in the words it sends back.
def folder = { String name -> File f = new File(work, name); f.mkdirs(); f }
def session = Session.newInstance(data + [folder: folder('new')])
Map about = session.described()
assert about.draft.protocol.sd_estimator == larger && about.newEstimator == larger && session.draft.protocol.engine_version == '2.1.0'
assert JsonOutput.toJson(about.draft) == JsonOutput.toJson(session.draft + [protocol: fromPage(session.draft.protocol as Map)]) : 'an untouched page must find its draft unchanged (Discard changes stays hidden)'
assert about.estimators == [[key: 'sample_sd', engine: '2.0.0', name: 'baseline SD'], [key: larger, engine: '2.1.0', name: 'larger of baseline SD and whole-trace noise']]
assert about.fijiDraft == null && !session.state().message.contains('Noise') && !session.noiseMoved
def remembered = { File into, Map protocol, Integer format = Settings.REMEMBERED_FORMAT -> new File(into, 'step4_last_settings.json').text = JsonOutput.toJson((format == null ? [:] : [format: format]) + [recording: 'rec', frames: 300, baseline_start: 1, baseline_end: 20, values: startingDraft.values, protocol: protocol]) }
String filled = 'Your last settings for this recording are filled in.', isNow = ' Noise is now the larger of baseline SD and whole-trace noise.', asSaved = ' Noise: baseline SD, as saved.'
// Settings remembered by a build that had the baseline SD only name no estimator: nobody chose it. Timing and thresholds are
// kept, noise is taken as a new protocol takes it, and the status line says so.
remembered(folder('earlier'), keptByEarlierBuild)
Map carried = [:]; Map kept = Settings.remembered(data + [folder: folder('earlier')], carried)
assert carried.noiseMoved && !carried.ruleMoved && kept.protocol.sd_estimator == larger && Protocol.canonical(kept.protocol, true).engine_version == '2.1.0'
assert kept.protocol.windows[0].subMap(Protocol.COLUMNS) == keptByEarlierBuild.windows[0].subMap(Protocol.COLUMNS) && kept.values == startingDraft.values
def earlier = Session.newInstance(data + [folder: folder('earlier')])
assert earlier.noiseMoved && earlier.stimulusFrom == 'last-session' && Protocol.draftEstimator(earlier.draft.protocol) == larger
assert earlier.state().message == filled + isNow : earlier.state().message
assert JsonOutput.toJson(earlier.described().draft) == JsonOutput.toJson(earlier.draft) && Protocol.draftEstimator(earlier.described().fijiDraft.protocol) == larger
assert JsonOutput.toJson(earlier.described().fijiDraft.protocol) == JsonOutput.toJson(fromPage(earlier.fijiDraft.protocol as Map))
// as the release wrote it: no format, the former default length and no estimator at the top
Map release = copy(keptByEarlierBuild); release.windows.each { it.Duration_mode = 'seconds'; it.Min_duration_s = 1 }
remembered(folder('release'), release, null)
def fromRelease = Session.newInstance(data + [folder: folder('release')])
assert fromRelease.state().message == filled + ' The rule now counts samples: 2 in a row. It was "for 1 s".' + isNow : fromRelease.state().message
assert Settings.check(fromRelease.draft.values, fromRelease.draft.protocol, true, data, [:], true).sd_estimator == larger
// Settings that name their estimator keep it: the baseline SD of an imported protocol, or this build's own.
remembered(folder('named'), fromPage(older))
def named = Session.newInstance(data + [folder: folder('named')])
assert !named.noiseMoved && named.draft.protocol.sd_estimator == 'sample_sd' && named.state().message == filled + asSaved : named.state().message
remembered(folder('this-build'), fromPage(made))
def again = Session.newInstance(data + [folder: folder('this-build')])
assert !again.noiseMoved && Protocol.canonical(again.draft.protocol, true).sd_estimator == larger && again.state().message == filled
// The status line follows the draft. "Use Fiji's values" on the page posts fijiDraft, a protocol made here.
named.updateDraft(named.described().fijiDraft as Map, 1)
assert Protocol.draftEstimator(named.draft.protocol) == larger && named.state().message == filled : 'the status line still says: ' + named.state().message
Map importedOnPage = copy(named.draft as Map); importedOnPage.protocol = fromPage(Protocol.exchange(asSimpleCa(older)))   // Import protocol JSON, then the page's own collectDraft
named.updateDraft(importedOnPage, 2)
assert named.state().message == filled + asSaved
earlier.updateDraft(importedOnPage, 1)
assert earlier.state().message == filled + asSaved : 'a draft on the baseline SD is not said to be on the new noise: ' + earlier.state().message
earlier.updateDraft(earlier.described().fijiDraft as Map, 2)
assert earlier.state().message == filled + isNow
// An error is not followed by a note, and after the first analysis the finish message carries the note of what was applied.
earlier.tone = 'bad'; assert !earlier.state().message.contains('Noise'); earlier.tone = 'plain'
// Remembered legacy settings: the page is handed the estimator of a new protocol with them.
remembered(folder('legacy'), legacy)
def historical = Session.newInstance(data + [folder: folder('legacy')])
assert historical.draft.protocol.schema_version == 1 && historical.state().message == filled && historical.described().draft.protocol.sd_estimator == Protocol.newEstimator()
Map switched = copy(historical.described().draft as Map); switched.protocol = madeOnPageFromLegacy   // the rule type changed on the page to version 2
historical.updateDraft(switched, 1)
assert !historical.state().validation.errors && Settings.check(historical.draft.values, historical.draft.protocol, true, data, [:], true).sd_estimator == larger : 'a version 2 rule made on the page from a legacy one is a protocol made here'
assert historical.state().message == filled
// a draft the page sends back is remembered with its estimator, and named settings are not moved the next time
session.updateDraft(about.draft as Map, 1)
carried = [:]; assert Protocol.canonical(Settings.remembered(data + [folder: folder('new')], carried).protocol, true).sd_estimator == larger && !carried.noiseMoved
carried = [:]; assert Settings.remembered(data + [folder: folder('named')], carried).protocol.sd_estimator == 'sample_sd' && !carried.noiseMoved
println 'PASS: a new session uses the new estimator; remembered settings that name no estimator take it too and the status line says so; settings that name one keep it; the status line follows the draft'

// The plugin's own numbers on both kinds of calls table: engine 2.0.0's, and the one with the three columns the new
// estimator appends (Whole_trace_noise, Whole_trace_noise_samples, Noise_source).
List cells = [['a', 'Stimulus', 'Responded', 0.12, 0.03, 'Whole-trace noise'], ['b', 'Stimulus', 'Did not respond', 0.45, 0.15, 'Whole-trace noise'], ['c', 'Stimulus', 'Excluded', '', '', ''],
    ['d', 'Stimulus', 'Responded', 0.40, 0.02, 'Baseline SD'], ['e', 'Stimulus', 'Responded', 0.30, 0.1, 'Baseline SD'], ['f', 'Stimulus', 'Did not respond', 0.20, 0.01, 'Noise floor'], ['a', 'Control', 'Did not respond', 0.12, 0.03, 'Whole-trace noise']]
def table = { String name, boolean newKind -> File dir = new File(work, name); dir.mkdirs()
    new File(dir, 'responder_calls.csv').text = '"Field","Cell","Role","Stimulus","SD_estimator","Status","Reason","Threshold_increase"' + (newKind ? ',"Whole_trace_noise","Whole_trace_noise_samples","Noise_source"' : '') + ',"Condition"\n' +
        cells.collect { '"rec","' + it[0] + '","' + it[1] + '","ATP","' + (newKind ? larger : 'sample_sd') + '","' + it[2] + '","' + (it[2] == 'Excluded' ? 'Local background unavailable for this ROI' : '') + '",' + it[3] + (newKind ? ',' + it[4] + ',' + (it[4] == '' ? '' : 299) + ',"' + it[5] + '"' : '') + ',""' }.join('\n') + '\n'
    new File(dir, 'responder_counts.csv').text = '"Field","Stimulus","Role","Registered","Responded","Evaluable"\n"rec","ATP","Stimulus",6,3,5\n'
    new File(dir, 'normalized_traces.csv').text = '"Time","a","b","c","d","e","f"\n0,0,0,0,0,0,0\n1,0.5,0.1,-2,-1.2,0.6,0.1\n2,0.4,0,0,0.9,0.5,-0.9\n'
    dir }
File oldKind = table('SimpleCa-old-kind', false), newKind = table('SimpleCa-new-kind', true)
def tagged = data + [folder: folder('tally'), chosenBlind: [a: true, b: true, c: true, d: false, e: false, f: false]]
Map counted = Settings.tally(oldKind, tagged)
assert counted.rows*.key == ['all', 'blind', 'later', 'quiet'] && counted.rows.find { it.key == 'all' }.counted == [3, 5, 6] && counted.rows.find { it.key == 'quiet' }.counted == [2, 3, 3] && counted.belowZero == ['d']
assert Settings.tally(newKind, tagged) == counted : 'the tally is the same with the three columns in the table'
String finish = Settings.percentResponding(oldKind, tagged)
assert finish == ' All cells: 3 of 5 responded (60.0 %). Chosen blind: 1 of 2 responded (50.0 %). Added after the response was shown: 2 of 3 responded (66.7 %).' && Settings.percentResponding(newKind, tagged) == finish
assert new File(tagged.folder as File, 'rec_percent_responding_SimpleCa-old-kind.csv').text == new File(tagged.folder as File, 'rec_percent_responding_SimpleCa-new-kind.csv').text
assert Settings.noEvaluableCells(oldKind) == null && Settings.noEvaluableCells(newKind) == null
[oldKind, newKind].each { File dir -> Map read = Session.records(Session.table(dir, 'responder_calls.csv')).find { it.Cell == 'b' && it.Role == 'Stimulus' }
    assert read.Status == 'Did not respond' && read.Threshold_increase == '0.45' && read.Condition == '' && (dir == newKind ? read.Noise_source == 'Whole-trace noise' : !read.containsKey('Noise_source')) }
println 'PASS: % responding, the quiet row and the finish message are the same from a calls table with or without the whole-trace noise columns'

// What set a cell's threshold. Threshold_source says "Baseline noise" for the noise term whichever estimate the noise is;
// Noise_source names the estimate. The Fiji window's line under the cell plot reads the two together.
String responders = new File('analysis/upstream/responders.R').getText('UTF-8')
['"Baseline noise" else "Floor"', '"Noise floor"', '"Whole-trace noise"', '"Baseline SD"'].each { assert responders.contains(it) : 'the pinned engine no longer writes ' + it + ': see CalciumCalls.thresholdSetBy' }
def Calls = loader.loadClass('CalciumCalls')
[['Baseline noise', 'Whole-trace noise', 'Whole-trace noise'], ['Baseline noise', 'Baseline SD', 'Baseline SD'], ['Baseline noise', 'Noise floor', 'Noise floor'], ['Floor', 'Whole-trace noise', 'Floor'], ['Floor', 'Baseline SD', 'Floor'],
 ['Baseline noise', null, 'Baseline noise'], ['Floor', null, 'Floor'], ['', '', ''], [null, null, '']].each { String term, String estimate, String expected ->
    Map call = [Threshold_source: term]; if (estimate != null) call.Noise_source = estimate
    assert Calls.thresholdSetBy(call) == expected : call }
def plot = loader.loadClass('CalciumAnalysisPlot').newInstance()
plot.responseCall = [Threshold_increase: '0.405', Threshold_source: 'Baseline noise', Noise_source: 'Whole-trace noise', Run_duration_s: '', Run_samples: '0']
assert plot.ruleLine() == 'Local threshold 0.405 (Whole-trace noise) · span unavailable s · 0 samples' : plot.ruleLine()
plot.responseCall = [Threshold_increase: '0.3034', Threshold_source: 'Baseline noise', Run_duration_s: '1', Run_samples: '2']
assert plot.ruleLine() == 'Local threshold 0.3034 (Baseline noise) · span 1 s · 2 samples' : 'a table scored on the baseline SD reads as before'
println 'PASS: the Fiji window names the estimate that set a threshold: whole-trace noise is not called baseline noise'

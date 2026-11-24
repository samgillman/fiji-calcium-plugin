args <- commandArgs(trailingOnly=TRUE)
engine <- normalizePath('analysis'); root <- if(length(args)) args[1] else 'build/analysis-fixtures'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
source(file.path(engine,'upstream/utils.R'));source(file.path(engine,'upstream/response_protocols.R'));source(file.path(engine,'upstream/responders.R'))
n<-30;time<-(0:(n-1))/2
response<-c(rep(0,8),seq(0,1,length.out=5),seq(1,0,length.out=5),rep(0,12))
x<-data.frame(Frame=1:n,Cell001=10*(1+response),Cell002=rep(10,n),Cell003=rep(0,n),Cell004=rep(NA_real_,n),Cell005=10*(1+c(rep(0,8),rep(1,22))))
x$Cell006<-x$Cell001;x$Cell006[14]<-NA_real_
write.csv(x,file.path(root,'input.csv'),row.names=FALSE,na='')
s<-c(sampling_hz='2',baseline_start='1',baseline_end='5',stimulus='Capsaicin',stimulus_start='3',stimulus_end='14.5',response_baseline_seconds='2.5',amplitude_floor='0.1',noise_multiplier='3',consecutive_samples='2',kcl_enabled='false',recording='Fixture',input_mode='background_corrected_fluorescence')
write.table(data.frame(names(s),unname(s)),file.path(root,'settings.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
script<-if(.Platform$OS.type=='windows') file.path(R.home('bin'),'Rscript.exe') else file.path(R.home('bin'),'exec','R')
run_adapter <- function(input, config, output) {
  system2(script,c('--vanilla',if(.Platform$OS.type=='windows') shQuote(file.path(engine,'runner.R')) else c('--slave',paste0('--file=',shQuote(file.path(engine,'runner.R'))),'--args'),shQuote(engine),shQuote(normalizePath(input)),shQuote(normalizePath(config)),shQuote(file.path(normalizePath(root),output))))
}
result<-run_adapter(file.path(root,'input.csv'),file.path(root,'settings.tsv'),'out')
stopifnot(result==0,file.exists(file.path(root,'out/COMPLETE')))
m<-read.csv(file.path(root,'out/metrics.csv'),check.names=FALSE,na.strings='')
y<-read.csv(file.path(root,'out/normalized_traces.csv'),check.names=FALSE,na.strings='')
for(cell in c('Cell001','Cell002','Cell005','Cell006')) {
  expected<-calculate_cell_metrics((x[[cell]]-mean(x[[cell]][1:5]))/mean(x[[cell]][1:5]),time,c(1,5),TRUE)
  got<-m[m$Cell==cell,names(expected),drop=FALSE];rownames(got)<-NULL
  for (column in names(expected)) { if (is.numeric(expected[[column]])) got[[column]] <- as.numeric(got[[column]]) }; stopifnot(isTRUE(all.equal(got,expected,tolerance=1e-12,check.attributes=FALSE)))
}
stopifnot(m$Peak_dFF0[1]==1,abs(m$AUC[1]-2.5)<1e-12,m$Status[3]=='Excluded',m$Status[4]=='Excluded',isTRUE(m$FWHM_Censored[5]),is.na(m$FWHM[5]))
counts<-read.csv(file.path(root,'out/responder_counts.csv'));calls<-read.csv(file.path(root,'out/responder_calls.csv'))
expected<-responder_field_results(y,names(x)[-1],'Fixture',responder_protocol('Capsaicin',3,14.5,2.5,0.1,3,2,NULL))
stopifnot(identical(calls$Status,expected$calls$Status),counts$Registered==6,counts$Excluded==2,counts$Unscorable==1,counts$Evaluable==3,counts$Responded==2)
stopifnot(!'Cell003'%in%names(y),is.na(y$Cell006[14]))
cat('PASS: runner matches pinned SimpleCa metrics and responder calls; invalid baselines remain registered; missing samples and censored widths preserved\n')

# All invalid F0 values: preserve the selected population, with no invented mean or denominator.
write.csv(x[c('Frame','Cell003','Cell004')],file.path(root,'invalid.csv'),row.names=FALSE,na='')
stopifnot(run_adapter(file.path(root,'invalid.csv'),file.path(root,'settings.tsv'),'invalid-out')==0)
invalid<-read.csv(file.path(root,'invalid-out/responder_counts.csv'),na.strings='')
stopifnot(invalid$Registered==2,invalid$Excluded==2,invalid$Evaluable==0,is.na(invalid$Percent_responding))
stopifnot(all(is.na(read.csv(file.path(root,'invalid-out/mean_trace.csv'),na.strings='')$Mean)))
# Fail invalid inputs explicitly; never write a completion marker.
bad<-x;names(bad)[3]<-names(bad)[2]
write.csv(bad,file.path(root,'duplicate.csv'),row.names=FALSE,na='')
stopifnot(run_adapter(file.path(root,'duplicate.csv'),file.path(root,'settings.tsv'),'duplicate-out')!=0)
bad_settings<-s;bad_settings['sampling_hz']<-'0'
write.table(data.frame(names(bad_settings),unname(bad_settings)),file.path(root,'bad-settings.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
stopifnot(run_adapter(file.path(root,'input.csv'),file.path(root,'bad-settings.tsv'),'bad-settings-out')!=0)
stopifnot(!file.exists(file.path(root,'duplicate-out/COMPLETE')),!file.exists(file.path(root,'bad-settings-out/COMPLETE')))
cat('PASS: all-excluded population, unavailable mean/percentage, duplicate columns and invalid frame rate\n')
# Separate KCl window: one KCl responder that did not respond to the primary stimulus.
kcl_input<-x;kcl_input$Cell002[time>=11]<-20
write.csv(kcl_input,file.path(root,'kcl.csv'),row.names=FALSE,na='')
kcl_settings<-s;kcl_settings['stimulus_end']<-'7.5';kcl_settings['kcl_enabled']<-'true'
kcl_settings<-c(kcl_settings,kcl_start='11',kcl_end='14.5')
write.table(data.frame(names(kcl_settings),unname(kcl_settings)),file.path(root,'kcl-settings.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
stopifnot(run_adapter(file.path(root,'kcl.csv'),file.path(root,'kcl-settings.tsv'),'kcl-out')==0)
kcl_counts<-read.csv(file.path(root,'kcl-out/responder_counts.csv'))
kcl_calls<-read.csv(file.path(root,'kcl-out/responder_calls.csv'))
stopifnot(nrow(kcl_calls)==12,nrow(kcl_counts)==3,kcl_counts$Responded[2]==1,kcl_counts$Registered[3]==1,kcl_counts$Responded[3]==0,kcl_counts$Percent_responding[3]==0)
cat('PASS: independent stimulus/KCl calls and explicit KCl-positive denominator\n')

# Independently specified 2 Hz fixture: one three-sample response, one two-sample
# spike, one nonresponder, one invalid F0 and one missing response sample.
clock2 <- seq(0,20,by=0.5)
gold <- data.frame(Frame=seq_along(clock2),Time_s=clock2,A=10,B=10,C=10,Invalid=0,Missing=10)
gold$A[clock2>=5 & clock2<=6] <- 12
gold$B[clock2>=5 & clock2<=5.5] <- 12
gold$Missing[clock2==5.5] <- NA_real_
for(cell in c('A','B','Missing'))gold[[cell]][clock2>=12 & clock2<=14] <- 13
write.csv(gold,file.path(root,'golden-input.csv'),row.names=FALSE,na='')
windows <- data.frame(Role=c('Stimulus','Control'),Stimulus=c('Capsaicin','Ionomycin'),
  Baseline_start_s=c(2,9),Baseline_end_s=c(5,12),Stimulus_onset_s=c(5,12),Start_s=c(5,12),End_s=c(8,15))
gold_protocol <- responder_protocol_v2('Native independent 2 Hz fixture',windows)
prepare_golden <- function(folder) {
  output <- file.path(root,folder);dir.create(output,showWarnings=FALSE)
  invisible(responder_export_protocol(gold_protocol,file.path(output,'protocol_input.json')))
  write.csv(data.frame(Cell=c('A','B','C','Invalid','Missing','Edge'),Included_in_measurement=c(rep(TRUE,5),FALSE),
    Exclusion_reason=c(rep('',5),'Edge exclusion policy'),Background_pixels=c(rep(100,5),NA)),file.path(output,'roi_registry.csv'),row.names=FALSE)
  output
}
out <- prepare_golden('golden-out')
gold_settings <- c(sampling_hz='',baseline_start='1',baseline_end='5',recording='Golden',input_mode='background_corrected_fluorescence',
  protocol_file='protocol_input.json',roi_registry_file='roi_registry.csv',condition='Vehicle',biological_unit='Culture1')
write.table(data.frame(names(gold_settings),unname(gold_settings)),file.path(root,'golden-settings.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
stopifnot(run_adapter(file.path(root,'golden-input.csv'),file.path(root,'golden-settings.tsv'),'golden-out')==0)
actual <- read.csv(file.path(out,'responder_calls.csv'),na.strings='',check.names=FALSE)
primary <- actual[actual$Role=='Stimulus',]
stopifnot(identical(primary$Status,c('Responded','Did not respond','Did not respond','Excluded','Unscorable')))
stopifnot(primary$Run_samples[1]==3,primary$Run_duration_s[1]==1,abs(primary$Window_peak_local_fraction[1]-0.2)<1e-12,
  abs(primary$Window_AUC_fraction_s[1]-0.25)<1e-12,primary$Peak_latency_s[1]==0)
counts2 <- read.csv(file.path(out,'responder_counts.csv'))
stopifnot(identical(counts2$Evaluable,c(3L,4L,2L)),identical(counts2$Responded,c(1L,3L,1L)),
  counts2$Population[3]=='Ionomycin-positive cells',abs(counts2$Percent_responding[1]-100/3)<1e-12,counts2$Percent_responding[3]==50)
quality2 <- read.csv(file.path(out,'cell_quality.csv'))
stopifnot(quality2$Status[1]!='Excluded',grepl('zero_normalization',quality2$QC_flags[1]))
flow <- read.csv(file.path(out,'cohort_flow.csv'))
stopifnot(flow$Accepted_ROIs==6,flow$Premeasurement_excluded==1,flow$Registered_for_analysis==5)
unit <- read.csv(file.path(out,'replicate_counts.csv'))
stopifnot(all(unit$Replicate=='Culture1'),all(unit$Fields==1))
manifest <- read.csv(file.path(out,'analysis_manifest.csv'))
stopifnot(manifest$value[manifest$key=='time_source']=='explicit_timestamps')

# Diagnostic grid and pseudo windows use shared functions but cannot change the
# golden applied calls/counts. Their expected labels/denominators are independent.
out_diag <- prepare_golden('golden-diagnostics')
write.csv(data.frame(Floor_local_fraction=c(.1,.5),SD_multiplier=3,Min_duration_s=1),file.path(out_diag,'sensitivity_grid.csv'),row.names=FALSE)
write.csv(data.frame(Diagnostic_window='Unstimulated',Baseline_start_s=0,Baseline_end_s=2,Stimulus_onset_s=2,Start_s=2,End_s=4),file.path(out_diag,'negative_control_windows.csv'),row.names=FALSE)
diag_settings <- c(gold_settings,sensitivity_file='sensitivity_grid.csv',negative_control_file='negative_control_windows.csv')
write.table(data.frame(names(diag_settings),unname(diag_settings)),file.path(root,'diagnostic-settings.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
stopifnot(run_adapter(file.path(root,'golden-input.csv'),file.path(root,'diagnostic-settings.tsv'),'golden-diagnostics')==0)
stopifnot(identical(read.csv(file.path(out_diag,'responder_counts.csv')),read.csv(file.path(out,'responder_counts.csv'))))
sen <- read.csv(file.path(out_diag,'sensitivity_counts.csv'))
neg <- read.csv(file.path(out_diag,'negative_control_counts.csv'))
stopifnot(all(sen$Analysis_kind=='sensitivity'),all(neg$Analysis_kind=='negative_control_diagnostic'),all(neg$Apparent_response_percent==0))
cat('PASS: independent 2 Hz duration/count/AUC golden values, generic control denominator, QC warnings, accepted-ROI attrition, timestamps and separate diagnostics\n')

# The app and Fiji both consume this frozen, hand-calculated fixture. Expected
# outcomes are CSV constants, not a second call to the implementation under test.
shared <- file.path('tests','fixtures','response-protocols')
fixed <- read.csv(file.path(shared,'golden-traces.csv'),check.names=FALSE,na.strings='')
raw <- data.frame(Frame=seq_len(nrow(fixed)),Time_s=fixed$Time)
for(cell in setdiff(names(fixed),'Time'))raw[[cell]] <- 100*(1+fixed[[cell]])
raw$InvalidF0 <- 0
write.csv(raw,file.path(root,'shared-golden-input.csv'),row.names=FALSE,na='')
shared_out <- file.path(root,'shared-golden-out');dir.create(shared_out,showWarnings=FALSE)
file.copy(file.path(shared,'golden-protocol.json'),file.path(shared_out,'protocol_input.json'),overwrite=TRUE)
shared_settings <- c(sampling_hz='',baseline_start='1',baseline_end='4',recording='Golden field',
  input_mode='background_corrected_fluorescence',protocol_file='protocol_input.json')
write.table(data.frame(names(shared_settings),unname(shared_settings)),file.path(root,'shared-golden-settings.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
stopifnot(run_adapter(file.path(root,'shared-golden-input.csv'),file.path(root,'shared-golden-settings.tsv'),'shared-golden-out')==0)
same_fixture <- function(actual_file,expected_file) {
  actual <- read.csv(actual_file,check.names=FALSE,na.strings='')
  expected <- read.csv(expected_file,check.names=FALSE,na.strings='')
  actual <- actual[,names(expected),drop=FALSE]
  for(column in names(expected)) {
    if(is.numeric(expected[[column]]))actual[[column]] <- as.numeric(actual[[column]])
  }
  stopifnot(isTRUE(all.equal(actual,expected,tolerance=1e-11,check.attributes=FALSE)))
}
same_fixture(file.path(shared_out,'responder_calls.csv'),file.path(shared,'golden-expected-calls.csv'))
same_fixture(file.path(shared_out,'responder_counts.csv'),file.path(shared,'golden-expected-counts.csv'))
stopifnot(nrow(read.csv(file.path(shared_out,'replicate_counts.csv')))==0,
          nrow(read.csv(file.path(shared_out,'unassigned_biological_units.csv')))==1)
cat('PASS: shared independent golden CSV calls/counts and no fabricated biological unit through the native adapter\n')

# What noise is taken from, through the adapter. A sample_sd protocol keeps the table engine 2.0.0 wrote. A protocol with
# the whole-trace estimate records engine 2.1.0 and gets the engine's three columns after the existing ones. The manifest
# names the estimator and carries the engine's sentence. Expected values below are worked out here, not by the engine.
clock3 <- 0:59
noise_input <- data.frame(Frame=seq_along(clock3),Time_s=clock3,Steady=100,Jumpy=100)
noise_input$Steady[clock3>=20 & clock3<=30] <- 150
# Quiet for the ten baseline samples, then jumping by 15 % from sample to sample: two samples in a row at +15 % and +16 %
# pass a threshold set by the baseline alone and fail one set by the trace's own frame-to-frame noise.
jump <- rep(c(-15,15),25)
jump[c(11,12)] <- c(15,16)
noise_input$Jumpy[11:60] <- 100+jump
noise_input$Jumpy[1:10] <- 100+rep(c(-0.1,0.1),5)
write.csv(noise_input,file.path(root,'noise-input.csv'),row.names=FALSE,na='')
noise_window <- data.frame(Role='Stimulus',Stimulus='ATP',Baseline_start_s=0,Baseline_end_s=10,Stimulus_onset_s=10,Start_s=10,End_s=59)
noise_run <- function(folder,protocol,extra=character(),edit=identity) {
  output <- file.path(root,folder);unlink(output,recursive=TRUE);dir.create(output)
  writeLines(edit(responder_export_protocol(protocol)),file.path(output,'protocol_input.json'),useBytes=TRUE)
  settings <- c(sampling_hz='',baseline_start='1',baseline_end='10',recording='Noise field',input_mode='background_corrected_fluorescence',protocol_file='protocol_input.json',extra)
  config <- file.path(root,paste0(folder,'-settings.tsv'))
  write.table(data.frame(names(settings),unname(settings)),config,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  list(exit=run_adapter(file.path(root,'noise-input.csv'),config,folder),dir=output)
}
read_manifest <- function(dir) { m <- read.csv(file.path(dir,'analysis_manifest.csv'),stringsAsFactors=FALSE);setNames(m$value,m$key) }
baseline_rule <- responder_protocol_v2('Baseline SD fixture',noise_window,duration_mode='samples',frames=2L)
trace_rule <- responder_protocol_v2('Whole-trace noise fixture',noise_window,duration_mode='samples',frames=2L,sd_estimator='larger_of_baseline_sd_and_trace_noise')
stopifnot(unique(baseline_rule$Engine_version)=='2.0.0',unique(trace_rule$Engine_version)=='2.1.0')
old_kind <- noise_run('noise-baseline-sd',baseline_rule,c(sd_estimator='sample_sd'))
new_kind <- noise_run('noise-whole-trace',trace_rule,c(sd_estimator='larger_of_baseline_sd_and_trace_noise'))
stopifnot(old_kind$exit==0,new_kind$exit==0,file.exists(file.path(old_kind$dir,'COMPLETE')),file.exists(file.path(new_kind$dir,'COMPLETE')))
old_calls <- read.csv(file.path(old_kind$dir,'responder_calls.csv'),check.names=FALSE,na.strings='',stringsAsFactors=FALSE)
new_calls <- read.csv(file.path(new_kind$dir,'responder_calls.csv'),check.names=FALSE,na.strings='',stringsAsFactors=FALSE)
added <- c('Whole_trace_noise','Whole_trace_noise_samples','Noise_source')
# every column of the sample_sd table keeps its place (Condition and Biological_unit included); the three come after them
stopifnot(!any(added%in%names(old_calls)),identical(names(new_calls)[seq_along(old_calls)],names(old_calls)),
  identical(names(new_calls)[length(old_calls)+1:3],added),ncol(new_calls)==ncol(old_calls)+3L,
  identical(tail(names(old_calls),2),c('Condition','Biological_unit')))
stopifnot(identical(old_calls$Status,c('Responded','Responded')),identical(new_calls$Status,c('Responded','Did not respond')),
  identical(new_calls$Noise_source,c('Baseline SD','Whole-trace noise')),all(new_calls$Whole_trace_noise_samples==59))
dff <- noise_input$Jumpy/mean(noise_input$Jumpy[1:10])-1
steps <- diff(dff)
expected_noise <- 1.4826*median(abs(steps-median(steps)))/sqrt(2)/(1+mean(dff[1:10]))
expected_sd <- sd(dff[1:10])/(1+mean(dff[1:10]))
stopifnot(expected_noise>0.2,expected_sd<0.002,abs(new_calls$Whole_trace_noise[2]-expected_noise)<1e-12,abs(new_calls$Baseline_SD[2]-expected_sd)<1e-12,
  abs(new_calls$Noise_SD_used[2]-expected_noise)<1e-12,abs(new_calls$Threshold_increase[2]-3*expected_noise)<1e-12,
  abs(old_calls$Threshold_increase[2]-0.1)<1e-12,new_calls$Whole_trace_noise[1]==0,abs(new_calls$Threshold_increase[1]-0.1)<1e-12)
# Threshold_source names the term that set a threshold, and says "Baseline noise" for the noise term whichever estimate it
# used; Noise_source names the estimate. The Fiji window reads the two together (CalciumCalls.thresholdSetBy).
stopifnot(identical(new_calls$Threshold_source,c('Floor','Baseline noise')),identical(old_calls$Threshold_source,c('Floor','Floor')))
# the measurements themselves do not depend on the estimator
for(column in c('Baseline_mean','Baseline_SD','Peak_increase','Peak_s','Window_AUC_fraction_s')) stopifnot(identical(old_calls[[column]],new_calls[[column]]))
old_record <- read_manifest(old_kind$dir);new_record <- read_manifest(new_kind$dir)
stopifnot(old_record[['sd_estimator']]=='sample_sd',old_record[['protocol_engine_version']]=='2.0.0',old_record[['engine_version']]==responder_engine_version(),
  old_record[['noise_rule']]=='Noise is the sample SD of the local prestimulus baseline.',
  new_record[['sd_estimator']]=='larger_of_baseline_sd_and_trace_noise',new_record[['protocol_engine_version']]=='2.1.0',new_record[['engine_version']]==responder_engine_version(),
  new_record[['noise_rule']]==responder_noise_rule(trace_rule),grepl('whole-trace noise',new_record[['noise_rule']],fixed=TRUE))
stopifnot(responder_import_protocol(file.path(old_kind$dir,'protocol.json'))$Engine_version=='2.0.0',responder_import_protocol(file.path(new_kind$dir,'protocol.json'))$SD_estimator=='larger_of_baseline_sd_and_trace_noise')
# Settings written before this build name no estimator and run as before.
unnamed <- noise_run('noise-unnamed',baseline_rule)
stopifnot(unnamed$exit==0,identical(readLines(file.path(unnamed$dir,'responder_calls.csv')),readLines(file.path(old_kind$dir,'responder_calls.csv'))))
# Refused, with no completion marker: settings that name another estimator than the protocol's, and a protocol that
# records its estimator with another engine release than the earliest that defines it.
mismatch <- noise_run('noise-settings-mismatch',trace_rule,c(sd_estimator='sample_sd'))
wrong_release <- noise_run('noise-wrong-release',trace_rule,edit=function(json)gsub('"2.1.0"','"2.0.0"',json,fixed=TRUE))
later_release <- noise_run('noise-later-release',baseline_rule,edit=function(json)gsub('"2.0.0"','"2.1.0"',json,fixed=TRUE))
for(run in list(mismatch,wrong_release,later_release)) stopifnot(run$exit!=0,!file.exists(file.path(run$dir,'COMPLETE')))
# The legacy rule records its one noise, the baseline SD.
legacy_record <- read_manifest(file.path(root,'out'))
stopifnot(legacy_record[['sd_estimator']]=='sample_sd',legacy_record[['protocol_engine_version']]=='legacy-1',legacy_record[['noise_rule']]=='Noise is the sample SD of the local prestimulus baseline.',
  !any(added%in%names(calls)))
cat('PASS: sample_sd keeps the 2.0.0 table; whole-trace noise appends its three columns and changes only what the threshold decides; the manifest names the estimator and the engine\'s sentence; mismatched settings or engine release are refused\n')

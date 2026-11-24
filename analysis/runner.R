# Fiji adapter. Scientific calculations come from the pinned SimpleCa2+ sources.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop('Expected engine directory, input CSV, settings TSV, output directory')
source(file.path(args[1], 'upstream', 'utils.R'))
source(file.path(args[1], 'upstream', 'response_protocols.R'))
source(file.path(args[1], 'upstream', 'responders.R'))
settings <- read.delim(args[3], header = FALSE, col.names = c('key', 'value'),
                       colClasses = 'character', quote = '', comment.char = '', check.names = FALSE)
if (anyDuplicated(settings$key)) stop('Duplicate settings')
s <- setNames(as.list(settings$value), settings$key)
num <- function(key) { v <- suppressWarnings(as.numeric(s[[key]])); if (length(v)!=1L || !is.finite(v)) stop(paste('Invalid',key)); v }
setting <- function(key, default='') if (is.null(s[[key]]) || is.na(s[[key]])) default else s[[key]]
progress <- function(p, text) { cat(sprintf('SIMPLECA\t%d\t%s\n', p, text)); flush.console() }
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
write <- function(x,name) write.csv(x,file.path(args[4],name),row.names=FALSE,na='')
canonical_calls <- function(x) x[,setdiff(names(x),c('Window_peak_local_dFF','Window_AUC_local_dFF_s')),drop=FALSE]
# The columns only the whole-trace estimator writes go last in every calls table, so each column of a sample_sd table
# keeps its place.
noise_columns_last <- function(x) {
  added <- intersect(c('Whole_trace_noise','Whole_trace_noise_samples','Noise_source'),names(x))
  x[,c(setdiff(names(x),added),added),drop=FALSE]
}
local_input <- function(key) {
  value <- setting(key)
  if (!nzchar(value) || basename(value)!=value) stop(paste('Expected snapshotted filename for',key))
  path <- file.path(args[4],value)
  if (!file.exists(path)) stop(paste('Missing',key))
  path
}
progress(5,'Reading full-precision measured fluorescence')
dt <- read.csv(args[2],check.names=FALSE,na.strings=c('','NA','NaN'),stringsAsFactors=FALSE)
if (ncol(dt)<2L || nrow(dt)<10L || names(dt)[1]!='Frame' || anyDuplicated(names(dt))) stop('Expected Frame and unique cell columns with at least 10 frames')
if (!identical(as.numeric(dt[[1]]),as.numeric(seq_len(nrow(dt))))) stop('Frames must be consecutive, starting at 1')
if (any(!vapply(dt[-1],function(x) is.numeric(x)||all(is.na(x)),logical(1)))) stop('Non-numeric fluorescence or timestamps')
first <- num('baseline_start'); last <- num('baseline_end')
if (first!=trunc(first)||last!=trunc(last)||first<1||last<first||last>=nrow(dt)) stop('Invalid normalization baseline range')
if (setting('input_mode')!='background_corrected_fluorescence') stop('Only unlifted background-corrected fluorescence is accepted')
if ('Time_s'%in%names(dt)) {
  time <- dt$Time_s
  if (any(!is.finite(time))||any(time<0)||any(diff(time)<=0)) stop('Time_s must contain finite, nonnegative, strictly increasing seconds')
  time_source <- 'explicit_timestamps'
} else {
  rate <- num('sampling_hz');if(rate<=0)stop('Frame rate must be positive')
  time <- (dt$Frame-1)/rate;time_source <- 'constant_fps'
}
cells <- setdiff(names(dt),c('Frame','Time_s'))
if (!length(cells)) stop('No registered cell columns')
registry <- NULL
if (nzchar(setting('roi_registry_file'))) {
  registry <- read.csv(local_input('roi_registry_file'),check.names=FALSE,stringsAsFactors=FALSE)
  needed <- c('Cell','Included_in_measurement','Exclusion_reason','Background_pixels')
  if (!all(needed%in%names(registry))||anyDuplicated(registry$Cell))stop('Invalid accepted-ROI registry')
  registry$Included_in_measurement <- tolower(as.character(registry$Included_in_measurement))%in%c('true','1','yes')
  if (!setequal(registry$Cell[registry$Included_in_measurement],cells))stop('ROI registry does not match measured cell columns')
}
normalized <- data.frame(Time=time)
quality <- metrics <- vector('list',length(cells))
for(i in seq_along(cells)) {
  result <- normalize_signal_trace(as.numeric(dt[[cells[i]]]),c(first,last),'corrected_fluorescence')
  q <- result$quality
  if(!is.null(registry)) {
    r <- registry[registry$Cell==cells[i],,drop=FALSE]
    if(!is.finite(r$Background_pixels)||r$Background_pixels<=0) {
      result$trace <- NULL;q$Status <- 'Excluded';q$Reason <- 'Local background unavailable for this ROI'
      q$QC_flags <- responder_qc_flags(q$QC_flags,'BACKGROUND_UNAVAILABLE')
    }
  }
  if(is.null(result$trace)) row <- empty_metrics_row() else {
    normalized[[cells[i]]] <- result$trace
    row <- calculate_cell_metrics(result$trace,time,c(first,last),data_is_dFF0=TRUE)
  }
  q <- cbind(Cell=cells[i],q,stringsAsFactors=FALSE)
  # Preserve familiar leading columns while keeping the full shared quality record.
  q <- q[,c('Cell','F0','Baseline_SD_raw','Status','Reason',setdiff(names(q),c('Cell','F0','Baseline_SD_raw','Status','Reason'))),drop=FALSE]
  quality[[i]] <- q
  metrics[[i]] <- cbind(q,Metric_status=if(is.finite(row$Peak_dFF0))'Available' else 'Unscorable',row)
  progress(10+floor(55*i/length(cells)),sprintf('Normalizing cell %d of %d',i,length(cells)))
}
qc <- do.call(rbind,quality);metrics <- do.call(rbind,metrics)
progress(70,'Validating the declared response protocol')
if(nzchar(setting('protocol_file'))) {
  protocol <- responder_import_protocol(local_input('protocol_file'))
} else {
  # Existing settings.tsv files retain the historical sample rule exactly.
  protocol <- responder_protocol(s$stimulus,num('stimulus_start'),num('stimulus_end'),
    num('response_baseline_seconds'),num('amplitude_floor'),num('noise_multiplier'),
    num('consecutive_samples'),NULL,control=identical(s$kcl_enabled,'true'),
    control_start=if(identical(s$kcl_enabled,'true'))num('kcl_start')else NULL,
    control_end=if(identical(s$kcl_enabled,'true'))num('kcl_end')else NULL)
}
# The noise behind the SD multiplier is part of the record. The legacy rule has one estimator, the sample SD of the
# baseline (responder_noise_rule). Settings that name another estimator than the protocol's own stop the run.
estimator <- if(responder_protocol_version(protocol)==1L)'sample_sd' else protocol$SD_estimator[1]
if(nzchar(setting('sd_estimator'))&&setting('sd_estimator')!=estimator)
  stop('The settings name the noise estimator ',setting('sd_estimator'),' but the protocol uses ',estimator)
invisible(responder_export_protocol(protocol,file.path(args[4],'protocol.json')))
applied <- responder_protocol_for_field(protocol,setting('recording'))
invisible(responder_export_protocol(applied,file.path(args[4],'applied_protocol.json')))
responders <- responder_field_results(normalized,cells,setting('recording'),protocol,cell_quality=qc)
if(responder_protocol_version(protocol)==1L)
  responders$calls <- responder_legacy_annotations(normalized,cells,setting('recording'),protocol,cell_quality=qc,legacy_result=responders)
# Canonical metric names describe fractional change for every supported indicator.
responders$calls <- canonical_calls(responders$calls)
responders$calls$Condition <- setting('condition');responders$calls$Biological_unit <- setting('biological_unit')
responders$calls <- noise_columns_last(responders$calls)
responders$counts$Condition <- setting('condition');responders$counts$Biological_unit <- setting('biological_unit')
write(metrics,'metrics.csv');write(normalized,'normalized_traces.csv');write(qc,'cell_quality.csv')
write(responders$calls,'responder_calls.csv');write(responders$counts,'responder_counts.csv')
window_columns <- c('Field','Cell','Role','Stimulus','Status','Reason','Eligible','Condition','Biological_unit','Baseline_start_s','Baseline_end_s','Stimulus_onset_s','Start_s','End_s','Peak_increase','Onset_s','Confirmation_s','Run_start_s','Run_end_s','Run_duration_s','Run_samples','Response_latency_s','Peak_latency_s','Window_peak_local_fraction','Window_AUC_fraction_s','Window_AUC_start_s','Window_AUC_end_s','Window_AUC_covered_s','Window_AUC_coverage_fraction','AUC_method','QC_flags')
write(responders$calls[,intersect(window_columns,names(responders$calls)),drop=FALSE],'window_metrics.csv')
flow <- data.frame(Field=setting('recording'),Accepted_ROIs=if(is.null(registry))length(cells)else nrow(registry),
  Premeasurement_excluded=if(is.null(registry))0L else sum(!registry$Included_in_measurement),Registered_for_analysis=length(cells),
  Normalization_valid=sum(qc$Status!='Excluded'),Normalization_excluded=sum(qc$Status=='Excluded'),Selection_basis=setting('selection_basis','Not recorded'),
  Note='Registered counts are measured cells; stimulus-specific evaluability is in responder_counts.csv')
write(flow,'cohort_flow.csv')
metadata <- data.frame(Field=setting('recording'),Replicate=setting('biological_unit'),Group=setting('condition','Unspecified condition'),stringsAsFactors=FALSE)
if(!nzchar(metadata$Group))metadata$Group <- 'Unspecified condition'
write(metadata,'recording_metadata.csv')
replication <- responder_replicate_summary(responders$counts,metadata)
if(!ncol(replication$replicates))replication$replicates <- data.frame(Group=character(),Replicate=character(),Stimulus=character(),Population=character(),Fields=integer(),Responded=integer(),Evaluable=integer(),Percent_responding=numeric())
if(!ncol(replication$groups))replication$groups <- data.frame(Group=character(),Stimulus=character(),Population=character(),Replicates_evaluable=integer(),Mean_replicate_percent=numeric(),SD_replicate_percent=numeric(),Cell_pooled_percent=numeric())
write(replication$replicates,'replicate_counts.csv');write(replication$groups,'group_counts.csv')
write(data.frame(Field=replication$unassigned_fields,Reason=rep('No independent biological-unit ID supplied',length(replication$unassigned_fields))),'unassigned_biological_units.csv')
progress(82,'Saving applied results and independent diagnostics')
if(nzchar(setting('sensitivity_file'))) {
  grid <- read.csv(local_input('sensitivity_file'),check.names=FALSE,stringsAsFactors=FALSE)
  diagnostic <- responder_sensitivity_results(normalized,cells,setting('recording'),protocol,grid,qc)
  write(noise_columns_last(canonical_calls(diagnostic$calls)),'sensitivity_calls.csv');write(diagnostic$counts,'sensitivity_counts.csv')
}
if(nzchar(setting('negative_control_file'))) {
  windows <- read.csv(local_input('negative_control_file'),check.names=FALSE,stringsAsFactors=FALSE)
  diagnostic <- responder_control_diagnostics(normalized,cells,setting('recording'),protocol,windows,qc)
  write(noise_columns_last(canonical_calls(diagnostic$calls)),'negative_control_calls.csv');write(diagnostic$counts,'negative_control_counts.csv')
}
# Cell SEM is descriptive. Biological-unit summaries are separate and never infer units from cell count.
traces <- as.matrix(normalized[-1]);finite <- is.finite(traces);n <- rowSums(finite)
avg <- if(ncol(traces))rowMeans(replace(traces,!finite,NA_real_),na.rm=TRUE)else rep(NA_real_,length(time))
sem <- if(ncol(traces))apply(traces,1,function(x)sd(x[is.finite(x)])/sqrt(sum(is.finite(x))))else rep(NA_real_,length(time))
write(data.frame(Time=time,Mean=avg,SEM=sem,N_Cells=n),'mean_trace.csv')
provenance <- read.delim(file.path(args[1],'provenance.tsv'),header=FALSE,col.names=c('key','value'),colClasses='character',quote='')
extra <- c(R_version=as.character(getRversion()),jsonlite_version=as.character(packageVersion('jsonlite')),
  protocol_schema_version=as.character(responder_protocol_version(protocol)),
  # engine_version scored this run; protocol_engine_version is the earliest engine release that defines the protocol's rule.
  engine_version=responder_engine_version(),protocol_engine_version=responder_protocol_engine_version(protocol),
  sd_estimator=estimator,noise_rule=responder_noise_rule(protocol),protocol_hash=responder_protocol_hash(protocol),applied_protocol_hash=responder_protocol_hash(applied),protocol_hash_algorithm='MD5 of canonical protocol JSON (identity, not a security signature)',
  time_source=time_source,normalization='(F-F0)/F0; fixed F0 mean of selected baseline frames',normalized_units='fraction',background='Already subtracted once by Fiji; no second subtraction',intensity_lift='None',
  AUC_interval='Legacy AUC: whole recording, signed net integral',Time_to_Peak_definition='Legacy peak: recording timestamp, not stimulus latency',
  window_metrics='Local prestimulus reference; observed samples in each declared response window; no extrapolation',
  QC_policy='Warnings retained; no automatic low-F0-to-SD or peak-outlier exclusion',denominator='Responded / (Responded + Did not respond); technical exclusions and unscorable records separate',
  replication='Cells are nested observations; biological-unit summaries require supplied IDs',diagnostics='Exploratory sensitivity and negative-control outputs do not replace applied results')
all_settings <- rbind(settings,provenance,data.frame(key=names(extra),value=unname(extra)))
all_settings <- all_settings[!duplicated(all_settings$key,fromLast=TRUE),,drop=FALSE]
write(all_settings,'analysis_manifest.csv')
writeLines('complete',file.path(args[4],'COMPLETE'))
progress(100,'Analysis complete')

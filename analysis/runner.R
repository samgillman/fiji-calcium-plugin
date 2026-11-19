# Fiji adapter. Scientific calculations are sourced unchanged from the pinned SimpleCa2+ files.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop('Expected engine directory, input CSV, settings TSV, output directory')
source(file.path(args[1], 'upstream', 'utils.R'))
source(file.path(args[1], 'upstream', 'responders.R'))
settings <- read.delim(args[3], header = FALSE, col.names = c('key', 'value'),
                       colClasses = 'character', quote = '', comment.char = '', check.names = FALSE)
if (anyDuplicated(settings$key)) stop('Duplicate settings')
s <- setNames(as.list(settings$value), settings$key)
num <- function(key) { v <- suppressWarnings(as.numeric(s[[key]])); if (length(v)!=1L || !is.finite(v)) stop(paste('Invalid',key)); v }
progress <- function(p, text) { cat(sprintf('SIMPLECA\t%d\t%s\n', p, text)); flush.console() }
progress(5, 'Reading measured fluorescence')
dt <- read.csv(args[2], check.names = FALSE, na.strings = c('', 'NA', 'NaN'), stringsAsFactors = FALSE)
if (ncol(dt)<2L || nrow(dt)<10L || names(dt)[1]!='Frame' || anyDuplicated(names(dt))) stop('Expected Frame and unique cell columns with at least 10 frames')
if (!identical(as.numeric(dt[[1]]), as.numeric(seq_len(nrow(dt))))) stop('Frames must be consecutive, starting at 1')
if (any(!vapply(dt[-1], function(x) is.numeric(x) || all(is.na(x)), logical(1)))) stop('Non-numeric fluorescence')
rate <- num('sampling_hz'); first <- num('baseline_start'); last <- num('baseline_end')
if (rate<=0 || first!=trunc(first) || last!=trunc(last) || first<1 || last<first || last>=nrow(dt)) stop('Invalid frame timing or baseline range')
if (s$input_mode != 'background_corrected_fluorescence') stop('Only unlifted background-corrected fluorescence is accepted')
time <- (dt[[1]]-1)/rate
cells <- names(dt)[-1]; normalized <- data.frame(Time=time)
qc <- data.frame(Cell=cells, F0=NA_real_, Baseline_SD_raw=NA_real_, Status='Included', Reason='', stringsAsFactors=FALSE)
metrics <- vector('list',length(cells))
for (i in seq_along(cells)) {
  x <- as.numeric(dt[[cells[i]]]); base <- x[first:last]
  # Match SimpleCa2+ build_processed_state raw-fluorescence normalization exactly.
  f0 <- mean(base[is.finite(base)], na.rm=TRUE)
  qc$F0[i] <- f0; qc$Baseline_SD_raw[i] <- sd(base[is.finite(base)])
  if (!is.finite(f0) || f0<=1e-6) {
    qc$Status[i] <- 'Excluded'; qc$Reason[i] <- 'Invalid raw fluorescence baseline (F0)'
    row <- empty_metrics_row()
  } else {
    y <- (x-f0)/f0
    normalized[[cells[i]]] <- y
    row <- calculate_cell_metrics(y, time, c(first,last), data_is_dFF0=TRUE)
    if (!is.finite(row$Peak_dFF0)) { qc$Status[i]<-'Unscorable'; qc$Reason[i]<-'No usable post-baseline metric (or fewer than 10 finite samples)' }
  }
  metrics[[i]] <- cbind(qc[i,,drop=FALSE], row)
  progress(10+floor(60*i/length(cells)),sprintf('Analyzing cell %d of %d',i,length(cells)))
}
metrics <- do.call(rbind,metrics)
progress(75,'Scoring stimulus responses')
protocol <- responder_protocol(s$stimulus, num('stimulus_start'), num('stimulus_end'),
  num('response_baseline_seconds'), num('amplitude_floor'), num('noise_multiplier'),
  num('consecutive_samples'), NULL, control=identical(s$kcl_enabled,'true'),
  control_start=if(identical(s$kcl_enabled,'true')) num('kcl_start') else NULL,
  control_end=if(identical(s$kcl_enabled,'true')) num('kcl_end') else NULL)
responders <- responder_field_results(normalized, cells, s$recording, protocol)
progress(85,'Saving metrics, traces and responder counts')
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
write <- function(x,name) write.csv(x,file.path(args[4],name),row.names=FALSE,na='')
write(metrics,'metrics.csv'); write(normalized,'normalized_traces.csv'); write(qc,'cell_quality.csv')
write(responders$calls,'responder_calls.csv'); write(responders$counts,'responder_counts.csv')
# Summary is cell-level, not biological replication. Excluded cells stay registered in calls/counts.
traces <- as.matrix(normalized[-1]); finite <- is.finite(traces)
n <- rowSums(finite)
avg <- if(ncol(traces)) rowMeans(replace(traces,!finite,NA_real_),na.rm=TRUE) else rep(NA_real_,length(time))
sem <- if(ncol(traces)) apply(traces,1,function(x) sd(x[is.finite(x)])/sqrt(sum(is.finite(x)))) else rep(NA_real_,length(time))
write(data.frame(Time=time,Mean=avg,SEM=sem,N_Cells=n),'mean_trace.csv')
provenance <- read.delim(file.path(args[1],'provenance.tsv'),header=FALSE,col.names=c('key','value'),colClasses='character',quote='')
manifest <- rbind(settings,provenance,data.frame(key=c('R_version','normalization','background','intensity_lift','AUC_interval','Time_to_Peak_definition','replication'),
 value=c(as.character(getRversion()),'(F-F0)/F0; F0 mean of selected baseline frames','Already subtracted once by Fiji; no second subtraction','None','Whole recording, signed net integral','Recording time of peak; not stimulus latency','Cell summaries; cells are not independent animals')))
write(manifest,'analysis_manifest.csv')
writeLines('complete',file.path(args[4],'COMPLETE'))
progress(100,'Analysis complete')

args <- commandArgs(trailingOnly=TRUE)
engine <- normalizePath('analysis'); root <- if(length(args)) args[1] else 'build/analysis-fixtures'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
source(file.path(engine,'upstream/utils.R'));source(file.path(engine,'upstream/responders.R'))
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

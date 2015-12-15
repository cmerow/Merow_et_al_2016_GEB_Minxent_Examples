 library(maptools)
 library(raster)
 species='CEOR'
 model.output='Maxent_output' # folder: all model runs go here
 final.model.output='Final_out' # folder: the directory of files used for the figure in ms
 ###############################################################################
 ## Functions and Setup
 ###############################################################################
 # normalization function; ensures predicted probabilities sum to 1. 
 norm=function(surf){ surf/sum(surf,na.rm=T)}
 # the following must be setup up on your own computer
 wd='/Users/ctg/Minxent/App1'
 # set where yout maxent jar file is
 maxent.location='/Users/ctg/Dropbox/MaxEnt/Program/maxent.jar'
 setwd(wd)
 # automatically create places to put maxent output
 if(!file.exists(model.output)) {dir.create(model.output)}
 # settings that apply to all models below, which will be supplied to the maxent software
 all.models=' nowarnings noprefixes responsecurves jackknife outputformat=raw 
  	removeduplicates noaskoverwrite -a -z threads=3 replicates=5 nothreshold 
  	nohinge noautofeature '
 ###############################################################################
 ## Sampling bias models
 ###############################################################################
 ## In this section we build two models for sampling bias based on anthropogenic factors.
 ## target group sampling bias based on anthropogenic predictors for IPANE data (New England)
 #Cory-- need to retain duplicates here, so "all.models" won't work
 bias_TG_IPANE=paste0('java -jar ',maxent.location, all.models,' -N bio3 -N bio4 
  	-N bio5 -N bio12 ')
 output=paste0(wd,'/',model.output,'/IPANE_TG_bias')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/Bias_IPANE_allPoints.csv')
 system(paste0(bias_TG_IPANE,'outputdirectory=',output,' environmentallayers=',
  	environmental,' samplesfile=',samples))
 ## note that you'll see some errors from points missing evnironmental data
 
 ################################################################################
 ## Species models
 ################################################################################
 # =========================================================================
 ## run MaxEnt for IPANE dataset with no sampling bias
 NE_noBias_C=paste0('java -jar ',maxent.location, all.models, ' -N pop_max -N roads_max ')
 output=paste0(wd,'/',model.output,'/',species,'_NE_noBias')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/',species,'_IPANE_5min_nodup.csv')
 system(paste0(NE_noBias_C,'outputdirectory=',output,' environmentallayers=',
  	environmental,' samplesfile=',samples))
 # =========================================================================
 ## run MaxEnt for IPANE dataset with Target group sampling bias
 NE_TGbias_C=paste0('java -jar ',maxent.location, all.models,' -N pop_max 
  	-N roads_max biastype=3 ')
 output=paste0(wd,'/',model.output,'/',species,'_NE_withTGbias')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/',species,'_IPANE_5min_nodup.csv')
 bias=paste0(wd,'/',model.output,'/IPANE_TG_bias/bias_avg.asc')
 system(paste0(NE_TGbias_C,'outputdirectory=',output,' environmentallayers=',
  	environmental,' samplesfile=',samples,' biasfile=',bias))
 ################################################################################
 ## Show output
 ################################################################################
 #read in all the output
 b=stack(c(paste0(paste0(wd,'/',model.output,'/',species),c(
  	#CEOR model with no sampling bias
  	paste0('_NE_noBias/',species,'_avg.asc'),  
  	## CEOR NE models with sampling bias model included
  	paste0('_NE_withTGbias/',species,'_avg.asc'))),
  	#IPANE TG sampling bias model 
  	paste0(wd,'/',model.output,'/IPANE_TG_bias/bias_avg.asc') 
  	))	
 names(b)=c('maxent.pred','sampling.minxent.pred','sampling.prior')
 values(b)=apply(values(b),2,norm)
 plot(b)

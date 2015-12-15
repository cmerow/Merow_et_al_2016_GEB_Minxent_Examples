 ###  Maxent for App 4: Native/Invasive Range Models #####################
 
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
 wd='/Users/ctg/Minxent/'
 # set where your maxent jar file is
 maxent.location='/Users/ctg/Dropbox/MaxEnt/Program/maxent.jar'
 setwd(wd)
 # automatically create places to put maxent output
 if(!file.exists(model.output)) {dir.create(model.output)}
 # settings that apply to all models below, which will be supplied to the maxent software
 all.models=' nowarnings noprefixes responsecurves jackknife outputformat=raw removeduplicates 
  	noaskoverwrite -a -z threads=3 replicates=5 nothreshold nohinge noautofeature '
 ###############################################################################
 ## Sampling bias models
 ###############################################################################
 ## In this section we build two models for sampling bias based on anthropogenic factors.
 ## target group sampling bias based on anthropogenic predictors for IPANE data (New England)
 bias_TG_IPANE=paste0('java -jar ',maxent.location, all.models,' -N bio3 -N bio4
  	 -N bio5 -N bio12 ')
 output=paste0(wd,'/',model.output,'/IPANE_TG_bias')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/Bias_IPANE_allPoints.csv')
 system(paste0(bias_TG_IPANE,'outputdirectory=',output,' environmentallayers=',environmental,' 
  	samplesfile=',samples))
 ## note that you'll see some errors from points missing evnironmental data
 
 ## target group sampling bias based on anthropogenic predictors for PRDB data (Japan)
 bias_TG_PRDB=paste0('java -jar ',maxent.location, all.models, ' -N bio3 -N bio4 
  	-N bio5 -N bio12 ')
 output=paste0(wd,'/',model.output,'/Japan_TG_bias')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/Japan_ASCII')
 samples=paste0(wd,'/Bias_PRDB_allPoints.csv')
 system(paste0(bias_TG_PRDB,'outputdirectory=',output,' environmentallayers=',environmental,' 
  	samplesfile=',samples))
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
 system(paste0(NE_noBias_C,'outputdirectory=',output,' environmentallayers=',environmental,' 
  	samplesfile=',samples))
 # =========================================================================
 ## run Japan model with TG sampling bias, project to NE
 Japan_TGbias_C=paste0('java -jar ',maxent.location, all.models, ' -N pop_max 
  	-N roads_max biastype=3 ')
 output=paste0(wd,'/',model.output,'/',species,'_Japan_NEproj')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/Japan_ASCII')
 samples=paste0(wd,'/PRDB_',species,'_5min_points.csv')
 bias=paste0(wd,'/',model.output,'/JAPAN_TG_bias/bias_avg.asc')
 projection=paste0(wd,'/NE_ASCII')
 system(paste0(Japan_TGbias_C,'outputdirectory=',output,' environmentallayers=',environmental,' 
  	samplesfile=',samples,' biasfile=',bias, ' projectionlayers=',projection))
 # =========================================================================
 ## run MaxEnt for IPANE dataset with TG sampling bias
 NE_TGbias_C=paste0('java -jar ',maxent.location, all.models,' -N pop_max 
  	-N roads_max biastype=3 ')
 output=paste0(wd,'/',model.output,'/',species,'_NE_withTGbias')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/',species,'_IPANE_5min_nodup.csv')
 bias=paste0(wd,'/',model.output,'/IPANE_TG_bias/bias_avg.asc')
 system(paste0(NE_TGbias_C,'outputdirectory=',output,' environmentallayers=',
  	environmental,' samplesfile=',samples,' biasfile=',bias))
 # =========================================================================
 ## run MaxEnt for IPANE dataset with no sampling bias and native range prior
 NE_noBias_nativePrior_C=paste0('java -jar ',maxent.location, all.models,' -N pop_max 
  	-N roads_max biastype=3 ')
 output=paste0(wd,'/',model.output,'/',species,'_NE_JapanPrior')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/',species,'_IPANE_5min_nodup.csv')
 bias=paste0(wd,'/',model.output,'/',species,'_Japan_NEproj/',
  	species,'_NE_ASCII_avg.asc')
 system(paste0(NE_noBias_nativePrior_C,'outputdirectory=',output,' environmentallayers=',
  	environmental,' samplesfile=',samples,' biasfile=',bias))
 # =========================================================================
 # use both native range and sampling prior
 # Create Combined Prior
 	#CEOR Japan projection to NE
 b=readAsciiGrid(paste0(wd,'/',model.output,'/',species,'_Japan_NEproj/',
  	species,'_NE_ASCII_avg.asc')) 
 	#IPANE TG sampling bias model
 b[[2]]=readAsciiGrid(paste0(wd,'/',model.output,'/IPANE_TG_bias/bias_avg.asc'))[[1]] 
 	#CEOR combined prior
 b[[3]]=norm(b[[1]]*b[[2]])  
 write.asciigrid(b[3], paste0(model.output,'/',species,'_NE_combined_prior.asc'))
 ### run MaxEnt for IPANE dataset with combined prior
 NE_cp_C=paste0('java -jar ',maxent.location, all.models,' -N pop_max -N roads_max biastype=3 ')
 output=paste0(wd,'/',model.output,'/',species,'_NE_CombinedPrior')
 if(!file.exists(output)) {dir.create(output)}
 environmental=paste0(wd,'/NE_ASCII')
 samples=paste0(wd,'/',species,'_IPANE_5min_nodup.csv')
 bias=paste0(wd,'/',model.output,'/',species,'_Japan_NEproj/',species,'_NE_ASCII_avg.asc')
 system(paste0(NE_cp_C,'outputdirectory=',output,' environmentallayers=',environmental,
  	' samplesfile=',samples,' biasfile=', bias))
 ################################################################################
 ### Normalize for Plotting 
 ################################################################################
 #read in all the output
 b=stack(c(paste0(paste0(wd,'/',model.output,'/',species),c(
  		#CEOR model with native range offset, no sampling bias
  	paste0('_NE_JapanPrior/',species,'_avg.asc'),  
  		#CEOR Japan projection to NE
  	paste0('_Japan_NEproj/',species,'_NE_ASCII_avg.asc'), 
  		# CEOR model with no sampling bias 
  	paste0('_NE_noBias/',species,'_avg.asc'), 
  		# CEOR NE models with sampling bias model included
  	paste0('_NE_withTGbias/',species,'_avg.asc'), 
  		# CEOR NE combined offset model	
  	paste0('_NE_CombinedPrior/',species,'_avg.asc'))), 
  		#IPANE TG sampling bias model 
  	paste0(wd,'/',model.output,'/IPANE_TG_bias/bias_avg.asc') 
  	))	
 names(b)=c('native.minxent.raw','native.prior','maxent.pred','sampling.minxent.pred',
  	'native.sampling.minxent.raw','sampling.prior')
 # modify the relevant parts of the output by multiplying by the right offsets
 	#CEOR factor Japan offset back into model without sampling bias
 b$native.minxent.pred=b[['native.minxent.raw']]*b[['native.prior']]  
 	#CEOR Combined native range * TG sampling bias offset
 b$native.sampling.prior=b[['sampling.prior']]*b[['native.prior']] 
 	#CEOR factor Japan offset back into model with sampling bias
 b$native.sampling.minxent.pred=b[['native.sampling.minxent.raw']]*b[['native.prior']] 
 values(b)=apply(values(b),2,norm)
 #=====================================================================
 ## write out ascii grids for further analysis and plotting
 if(!file.exists(final.model.output)) dir.create(final.model.output)
 	## CEOR NE model with no sampling bias, native offset factored back in
 writeRaster(b[['native.minxent.pred']],file=paste0( 'Final_out/',species,
  	'_NE_noBias_withNativePrior_norm.asc'),format='ascii',overwrite=T) 
 	## CEOR NE model with no sampling bias, no native range offset
 writeRaster(b[['maxent.pred']],file=paste0( 'Final_out/',species, 
  	 '_NE_noBias_noNativePrior_norm.asc'),format='ascii',overwrite=T) 
 	## CEOR Japan model projected to NE
 writeRaster(b[['native.prior']],file=paste0( 'Final_out/',species, 
  	'_NE_Japan_prediction_norm.asc'),format='ascii',overwrite=T) 
 	## IPANE TG sampling bias model
 writeRaster(b[['sampling.prior']],file=paste0( 'Final_out/',  
  	'IPANE_TG_sampling_bias_norm.asc'),format='ascii',overwrite=T) 
 	#CEOR NE model that includes IPANE TG sampling bias
 writeRaster(b[['sampling.minxent.pred']],file=paste0( 'Final_out/',species,   
  	'_NE_withSamplingBias_norm.asc'),format='ascii',overwrite=T) 
 	 ## CEOR combined offset 
 writeRaster(b[['native.sampling.prior']],file=paste0( 'Final_out/',species,  
  	'_NE_CombinedPriorSurface_norm.asc'),format='ascii',overwrite=T)
 	## CEOR combined offset prediction, with Japan factored back in
 writeRaster(b[['native.sampling.minxent.pred']],file=paste0( 'Final_out/',species,  
  	'_NE_CombinedPrior_prediction_norm.asc'),format='ascii',overwrite=T) 
 #=====================================================================
 ## plotting
 ## make a figure as in figure 5 in the main text
 
 # set color scheme
 cols1=function(x,bias=1) { 
  	colorRampPalette(c('steelblue4','steelblue1','gold','goldenrod1','red1','red4'),
  	bias=bias)(x) 
  }
 path=paste0(getwd(),"/Final_out/")
 tmp1=stack(paste0(path,c('IPANE_TG_sampling_bias_norm.asc', 'CEOR_NE_withSamplingBias_norm.asc',
  	'CEOR_NE_CombinedPriorSurface_norm.asc','CEOR_NE_CombinedPrior_prediction_norm.asc'))) 
 # normalize
 values(tmp1)=apply(values(tmp1),2,function(x) {x/sum(x,na.rm=T)})
 # make the plot
 pdf(paste0(path,'App6_predictions_CEOR.pdf'),h=9,w=6)
 par(mar=c(.2,.1,1.5,.1),mfrow=c(3,2),oma=c(0,6,3,.2))
 z.max=max(values(tmp1),na.rm=T)
 min.quant=0 # the quantile below which grey scale is used.
 breaks=c(-100,quantile(values(tmp1), seq(min.quant,1,length=100),na.rm=T))
 lab=letters[1:10]
 for(i in 1:dim(tmp1)[3]){
   image(tmp1[[i]] ,col=cols1(100,bias=1),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2,
   	 zlim=c(0,z.max), breaks=breaks)
   text(-73,47,lab[i],cex=1.3)
  }
 lab.ats=c(1-.165,.5)
 mtext('Samping\nPrior\n\n',2,at=lab.ats[1],outer=T,line=0,las=1)
 mtext('Samping\nPrior +\n Native\nRange\nPrior',2,at=lab.ats[2],outer=T,line=0,las=1)
 mtext('Prior',3,at=.25,outer=T,line=0)
 mtext('Prediction',3,at=.75,outer=T,line=0)
 # add legend
 image(tmp1[[2]] ,col='white',xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, zlim=c(0,z.max))
 image.plot(tmp1[[2]],legend.only=TRUE,col=cols1(100,bias=1),xaxt='n', yaxt='n', bty='n',main='',
  	cex.main=2.2, legend.args=list(text="relative\noccurrence rate",cex=1.2),
  	smallplot=c(0,1,.55,.65),zlim=range(log(breaks[-1]),na.rm=T),
  	breaks=c(-100,log(breaks[-1])),
  	horizontal=TRUE,axis.args=list( 
  	at=log(quantile(values(tmp1), c(0,.05,.5,1),na.rm=T)), 
  	labels=sprintf("%.0e",quantile(values(tmp1), c(0,.05,.5,1),na.rm=T))))
 dev.off()
 system(paste0('open ',path,'App6_predictions_CEOR.pdf'))
 

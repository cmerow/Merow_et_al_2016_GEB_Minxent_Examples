# this is the demo version of the full code
 # set your working directory, which should contain all files downloaded for this demo
 setwd('/Users/ctg/Minxent/App2')
 # automatically creat places to put maxent output
 if(!file.exists('Maxent_output')) {dir.create('Maxent_output')}
 if(!file.exists('Maxent_output/Default')) {dir.create('Maxent_output/Default')}
 if(!file.exists('Maxent_output/Prior')) {dir.create('Maxent_output/Prior')}
 # set where yout maxent jar file is
 maxent.location='../../../Program/maxent.jar'
 ######################################################################
 ### Load in results of the dispersal model
 ######################################################################
 # load a raster brick with the predictions from another model (the dispersal model from
 # Merow et al. 2011, Am Nat).
 (load(paste0(getwd(),'/Dispersal_Prior_Raster.rdata')))
 # load a list of the occurrence data split up by observation year. 
 # the first element is all points observed by 1940, 
 # the second element is all points observed by 1960, etc.
 (load(paste0(getwd(),'/data_by_year.rdata')))
 ######################################################################
 ### RUN MODELS
 ######################################################################
 #___________________________________________________
 #Run default maxent model with presences at each time slice
 #___________________________________________________
 run.default=paste0('java -jar ',maxent.location, ' nowarnings noprefixes 
  	-e NE_environmental_grids  
  	nowarnings  threads=3 testsamplesfile=CEOR_pres_All.csv nowriteclampgrid 
  	nowritemess
  	 noaskoverwrite -q -p -h nothreshold noautofeature noresponsecurves  
  	 outputformat=raw -a -z')
 #1980
 system(paste(run.default,'-o Maxent_output/Default -s ceor_pres1980.csv'))
 me.unif=raster('Maxent_output/Default/CEOR.asc')
 #___________________________________________________
 #Run maxent with CA model as bias surface
 #___________________________________________________
 run.prior=paste0('java -jar ',maxent.location,' nowarnings noprefixes -e NE_environmental_grids 
  	testsamplesfile=CEOR_pres_All.csv nowriteclampgrid nowritemess noaskoverwrite -q -p 
  	nothreshold -h noautofeature noresponsecurves  outputformat=raw -a -z ')
 	#1980
 system(paste(run.prior,'-o Maxent_output/Prior -s ceor_pres1980.csv 
  	biasfile=CEOR1980.asc biastype=3'))	
 me.pot=raster('Maxent_output/Prior/CEOR.asc')
 #___________________________________________________
 #make realized distribution (Dispersal factored back in)
 #___________________________________________________
 #maxent bias surface minimizes KL to to dispersal model, so you'd think that 
 # you're predictions are like the realized distrbution. but it actually factors that 
 # prior distribution back out, so you're closer to the potential distribution. 
 # we need to multiply predictions from using the CA model 
 # by the CA model to get back to the realized disribution.	
 
 # multiply prior back in to prediction
 me.real=me.pot
 values(me.real)=values(me.real)*values(disp.prior)[,1]
 values(me.real)=values(me.real)/sum(values(me.real),na.rm=T)
 ######################################################################
 # Plots 
 ######################################################################
 pdf('CEOR_All_Dispersal_Priors_Demo.pdf',h=10,w=1.8)
  titles=c("1980")
  par(mar=c(.2,.1,1.5,.1),mfrow=c(5,1),oma=c(0,3,0,.2))
  tmp1=stack(me.unif,me.pot,me.real,disp.prior)
  z.max=max(values(tmp1),na.rm=T)
  breaks=c(-100,quantile(values(tmp1), seq(.25,1,length=100),na.rm=T))
  image(disp.prior[[3]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main=titles[i],cex.main=2.2, 
   	zlim=c(0,z.max), breaks=breaks)
  tmp=data.by.year[[3]]
  coordinates(tmp)=c(2,3)
  points(tmp,col='black',pch=19,cex=.5)
  image(me.unif[[1]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, 
   	zlim=c(0,z.max), breaks=breaks)
   tmp=read.csv('CEOR_pres_All.csv')
 	coordinates(tmp)=c(2,3)
 	points(tmp,col='black',pch=19,cex=.5)
  image(me.pot[[1]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, 
   	zlim=c(0,z.max), breaks=breaks)
   tmp=read.csv('CEOR_pres_All.csv')
 	coordinates(tmp)=c(2,3)
 	points(tmp,col='black',pch=19,cex=.5)
  image(me.real[[1]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, 
   	zlim=c(0,z.max),breaks=breaks)
  	tmp=data.by.year[[3]]
 	coordinates(tmp)=c(2,3)
 	points(tmp,col='black',pch=19,cex=.5)
  mtext('Minxent\nRealized Distribution',2,at=.3,outer=T,line=0)
  mtext('Minxent\nPotential Distribution',2,at=.5,outer=T,line=0)
  mtext('Dispersal Model\n(Prior)',2,at=.9,outer=T,line=0)
  mtext('Maxent (Uniform Prior)\nPotential Distribution',2,at=.7,outer=T,line=0)
  # add legend
  image(me.real[[i]] ,col='white',xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, zlim=c(0,z.max))
  image.plot(me.real[[i]],legend.only=TRUE,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',
   	cex.main=2.2, legend.args=list(text="relative\noccurrence rate",cex=1.2),
  	smallplot=c(0,.9,.75,.85),zlim=range(log(breaks[-1])),breaks=c(-100,log(breaks[-1])),
  	horizontal=TRUE,axis.args=list( at=log(quantile(values(tmp1),
  	c(.01,.35,.75,1),na.rm=T)), 
  	labels=c('0','1e-5','3e-4','2e-2')))
 dev.off()
 system('open CEOR_All_Dispersal_Priors_Demo.pdf')
 # note that the colors are slightly different from those in the main text because 
 # the quantiles of the predictions are used for the color scale and the quantile 
 # of the predictions aggregated across all time 
 # steps differ from those for models from just 1980.

# this is the demo version of the full code

########################################################################
### Setup
###########################################################################
# set your working directory, which should contain all files downloaded for this demo
setwd("~/Dropbox/My_Papers/Merow_et_al_2016_GEB_Minxent_Examples/Appendix_Data/App2_Dispersal")

# automatically creat places to put maxent output
if(!file.exists('Maxent_output')) {dir.create('Maxent_output')}
if(!file.exists('Maxent_output/Default')) {dir.create('Maxent_output/Default')}
if(!file.exists('Maxent_output/Prior')) {dir.create('Maxent_output/Prior')}

# set where yout maxent jar file is
maxent.location='~/Dropbox/MaxEnt/Program/maxent.jar'

#load necessary package for plotting
library(fields)
library(raster)

######################################################################
### Load in results of the dispersal model
######################################################################
# load a raster ('disp.prior') with the predictions from another model (the dispersal model from
# Merow et al. 2011, Am Nat) for the demo year (1980).
# Note: requires 'raster' package
(load(paste0(getwd(),'/Dispersal_Prior_Raster.rdata')))

# load a list of the occurrence data for the demo year 1980 ('data.by.year'). 
(load(paste0(getwd(),'/data_by_year.rdata')))

######################################################################
### RUN MODELS
######################################################################
#___________________________________________________
#Run default maxent model with presences at each time slice
#___________________________________________________
#sets up Maxent model call 
run.default=paste0('java -jar ',maxent.location, ' nowarnings noprefixes -e NE_environmental_grids nowarnings  threads=3 testsamplesfile=CEOR_pres_All.csv nowriteclampgrid nowritemess noaskoverwrite -q -p -h nothreshold noautofeature noresponsecurves outputformat=raw -a -z')

#runs model for 1980 only for demonstration
system(paste(run.default,'-o Maxent_output/Default -s ceor_pres1980.csv'))

me.unif=raster('Maxent_output/Default/CEOR.asc') #raster of default prediction

#___________________________________________________
#Run maxent with CA model as bias surface
#___________________________________________________
#sets up Maxent model call
run.prior=paste0('java -jar ',maxent.location,' nowarnings noprefixes -e NE_environmental_grids testsamplesfile=CEOR_pres_All.csv nowriteclampgrid nowritemess noaskoverwrite -q -p nothreshold -h noautofeature noresponsecurves  outputformat=raw -a -z ')

#runs model for 1980 only for demonstration
system(paste(run.prior,'-o Maxent_output/Prior -s ceor_pres1980.csv 
	biasfile=CEOR1980.asc biastype=3'))	

me.pot=raster('Maxent_output/Prior/CEOR.asc') #raster of prediction with CA model as bias surface

#___________________________________________________
#make realized distribution (Dispersal factored back in)
#___________________________________________________
# Maxent bias surface minimizes KL to dispersal model, so you'd think that 
# your predictions are like the realized distrbution. However, Maxent actually factors that 
# offset back out, so you're closer to the potential distribution. 
# We need to multiply predictions from using the CA model as an offset
# by the CA model itself to factor it back in and get back to the realized disribution.	

# multiply offset (CA model) back into prediction
me.real=me.pot
values(me.real)=values(me.real)*values(disp.prior)
values(me.real)=values(me.real)/sum(values(me.real),na.rm=T)

######################################################################
# Plots 
######################################################################
# cory's favorite color palette
cols1=function(x,bias=1) { colorRampPalette(c('steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}

pdf('CEOR_All_Dispersal_Priors_Demo.pdf',h=10,w=1.8) 
#set up plot for 5 panels to fit maps and legend
titles=c("1980")
par(mar=c(.2,.1,1.5,.1),mfrow=c(5,1),oma=c(0,3,0,.2))
tmp1=stack(me.unif,me.pot,me.real,disp.prior)
z.max=max(values(tmp1),na.rm=T)
breaks=c(-100,quantile(values(tmp1), seq(.25,1,length=100),na.rm=T))

# plot dispersal model for 1980 with known occurrences to 1980
image(disp.prior[[1]] ,col=cols1(100),main=titles, xaxt='n', yaxt='n', bty='n',cex.main=2.2, zlim=c(0,z.max), breaks=breaks)
tmp=data.by.year
coordinates(tmp)=c(2,3)
points(tmp,col='black',pch=19,cex=.5)

# plot default Maxent model for 1980 with all known occurrences
image(me.unif[[1]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, zlim=c(0,z.max), breaks=breaks)
tmp=read.csv('CEOR_pres_All.csv')
coordinates(tmp)=c(2,3)
points(tmp,col='black',pch=19,cex=.5)

# plot potential distribution from Minxent model with all known occurrences
image(me.pot[[1]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, zlim=c(0,z.max), breaks=breaks)
tmp=read.csv('CEOR_pres_All.csv')
coordinates(tmp)=c(2,3)
points(tmp,col='black',pch=19,cex=.5)

# plot realized distribution  from Minxent model with occurrences to 1980
image(me.real[[1]] ,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, zlim=c(0,z.max),breaks=breaks)
tmp=data.by.year
coordinates(tmp)=c(2,3)
points(tmp,col='black',pch=19,cex=.5)

# Adds model lables
mtext('Dispersal Model\n(Prior)',2,at=.9,outer=T,line=0)
mtext('Maxent (Uniform Prior)\nPotential Distribution',2,at=.7,outer=T,line=0)
mtext('Minxent\nPotential Distribution',2,at=.5,outer=T,line=0)
mtext('Minxent\nRealized Distribution',2,at=.3,outer=T,line=0)

# add legend
image(me.real ,col='white',xaxt='n', yaxt='n', bty='n',main='',cex.main=2.2, zlim=c(0,z.max)) # space holder
image.plot(me.real,legend.only=TRUE,col=cols1(100),xaxt='n', yaxt='n', bty='n',main='', cex.main=2.2, legend.args=list(text="relative\noccurrence rate",cex=1.2),smallplot=c(0,.9,.75,.85),zlim=range(breaks[-1]),horizontal=TRUE)

dev.off()
system('open CEOR_All_Dispersal_Priors_Demo.pdf')

# note that the colors are slightly different from those in the main text because 
# the quantiles of the predictions are used for the color scale and the quantile 
# of the predictions aggregated across all time 
# steps differ from those for models from just 1980.

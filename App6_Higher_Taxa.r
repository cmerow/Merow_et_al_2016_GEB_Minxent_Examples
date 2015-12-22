## ******************************************************************** ##
 ## Purpose:
 ## Demonstrate use of incorporating higher taxon information in 
 ## creating species distribution models with MaxEnt by
 ## using clade level information to model **Protea nana** distribution
 ##
 ## NOTE - the code below is a simple example of how higher order taxon 
 ## information can be used to model a relatively rare species. For the
 ## analysis presented in the manuscript, we carried out k-fold modeling
 ## to examine the effects of minxent on model uncertainty. While we
 ## recommend this approach, we present a simpler example below.
 ##
 ## ******************************************************************** ##
 
 ## Set working directory - user system dependent
 setwd( "/Users/ctg/Dropbox/My_Papers/Merow_et_al_2016_GEB_Minxent_Examples/Appendix_Data/App6_Higher_Taxa" )
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 ## Setup
 ## -------------------------------------------------------------------- ##
 ## * Load packages and define a couple of paths and values used 
 ## throughout the script.
 ## * Define functions used in scripts.
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 
 ## Pacakges used
 library( dismo )
 library( rgeos )
 ## Define maxent.jar file location
 ## **NOTE** This needs to be changed by the user
 ## -------------------------------------------------------------------- ##
 maxent.file <- "~/Dropbox/MaxEnt/Program/maxent.jar"
 ## Set directory for environmental variables
 ## ***
 ## In this application we are using the variables used in Latimer et al. 
 ## 2006, and other projects from this group. This include climate layers
 ## from the Schulze dataset and environmental layers derived from soil
 ## maps. In total there were 24 variables
 ## -------------------------------------------------------------------- ##
 environmental.layers <- 
    "SA_ASCII/"
 ## Compile default maxent command
 ## -------------------------------------------------------------------- ##
 run.default <- paste0( "java -jar ", maxent.file, " nowarnings noprefixes ",
                         "-e ", environmental.layers, " nowarnings threads=1 ", 
                         " nowriteclampgrid nowritemess noresponsecurves outputformat=raw -a -z ",
                         "noaskoverwrite ",
                         "nothreshold nohinge noautofeature noproduct " ) 
 ## ******************************************************************** ##
 ## FUNCTION: normalize_ascii
 ## ensures that a raster of predicted probabilities sums to 1
 ## ******************************************************************** ##
 normalize_ascii <- function( x ){
   values( x ) <-values( x ) / sum( values( x ), na.rm = TRUE )
   return( x)
 }
 ## ******************************************************************** ##
 ## FUNCTION: multiply.offset
 ## Factors the offset/prior in to a prediction. By default, Maxent 
 ## factors out the offset and for some predictions (informative offsets)
 ## it is useful to multiple the offset back in to the prediction.
 ## ******************************************************************** ##
 
 multiply.offset <- function( prior.asc, maxent.asc ){
   ## Multiply the prior and maxent output asciis
   asc_with_bias <- maxent.asc
   values( asc_with_bias ) <-
     values( maxent.asc ) * values( prior.asc )
   ## Normalize the new layer
   asc_with_bias <- normalize_ascii( x = asc_with_bias ) 
   ## Return the ascii layer
   return( asc_with_bias)
 }
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 ## Make a Maxent model for P. nana using PRECIS data
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 if( !file.exists( "Maxent_output/")) dir.create( "Maxent_output/" )
 if( !file.exists( "Maxent_output/PRNANA_PRECIS" ) ) { 
  	dir.create( "Maxent_output/PRNANA_PRECIS" ) }
 ## Run the naieve maxent model using the PRECIS data
 system( paste0( run.default, " -o ", "Maxent_output/PRNANA_PRECIS",
                  " -s ", "prnana_precis_train.csv" ) ) 
 ## Read in the resulting ASCII map, which will be used later
 prnana_precis_asc <- raster( "Maxent_output/PRNANA_PRECIS/prnana.asc" )
 ## Normalize prnana_precis_asc
 prnana_precis_asc <- normalize_ascii( prnana_precis_asc )
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 ## Build a Maxent model for the Rose Protea species using PRECIS data
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 
 if( !file.exists( "Maxent_output/PR_ROSE" ) ) { dir.create( "Maxent_output/PR_ROSE" ) }
 ## Run the naieve maxent model using the PRECIS data
 system( paste0( run.default, " -o ", "Maxent_output/PR_ROSE/",
                  " -s ", "pr_rose_train.csv" ) ) 
 ## Read in the ascii layer resulting from the white protea model
 ## This will serve as the prior (bias layer)
 pr_rose_asc <- raster( "Maxent_output/PR_ROSE/pr_rose.asc" )
 ## Normalize this ascii
 pr_rose_asc <- normalize_ascii( x = pr_rose_asc )
 ## Write this layer to ascii
 writeRaster( x = pr_rose_asc, filename = "pr_rose_prior.asc", overwrite = TRUE )
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 ## Run a Maxent model for PRNANA, using Rose Protea model as a prior layer
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 
 if( !file.exists( "Maxent_output/PRNANA_PRECIS_bias" ) ) { dir.create( 
  	"Maxent_output/PRNANA_PRECIS_bias" ) }
 ## Run the naieve maxent model using the PRECIS data
 system( paste0( run.default, " -o ", "Maxent_output/PRNANA_PRECIS_bias/",
                  " -s ", "prnana_precis_train.csv",
                  " biasfile=", "pr_rose_prior.asc",
                  " biastype=3 ") ) 
 ## Read model with bias layer removed
 prnana_rm_bias_asc <- raster( "Maxent_output/PRNANA_PRECIS_bias/prnana.asc" )
 ## Normalize this ascii layer
 prnana_rm_bias_asc <- normalize_ascii( prnana_rm_bias_asc )
 ## -------------------------------------------------------------------- ##
 ## Create minxent layers
 ## -------------------------------------------------------------------- ##
 
 prnana_minxent <- multiply.offset( prior.asc = pr_rose_asc, maxent.asc = prnana_rm_bias_asc )
 ## Write an ascii file for the minxent layer
 writeRaster( x = prnana_minxent, filename = "prnana_minxent.asc", overwrite = TRUE )

 ## ******************************************************************** ##
 ## ******************************************************************** ##
 ## Plot results
 ## ******************************************************************** ##
 ## ******************************************************************** ##
 to.plot=stack( prnana_precis_asc, pr_rose_asc,  prnana_minxent)
 plot(to.plot)
 
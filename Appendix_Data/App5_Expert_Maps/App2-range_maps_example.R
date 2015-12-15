## ******************************************************************** ##
## App2-range_maps_example.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2014-10-26
##
## Purpose:
## Demonstrate the use of binary expert range maps as prior information
## in Maxent models
##
## ******************************************************************** ##

## Set working directory - user system dependent
setwd( "App2-Expert-Range-Maps/" )

## ******************************************************************** ##
## Setup
## -------------------------------------------------------------------- ##
## * Load packages and define a couple of paths and values used 
## throughout the script.
## * Define functions used in scripts.
## ******************************************************************** ##

## Pacakges used
library( alphahull )
library( dismo )
library( rgeos )

## Define maxent.jar file location
## **NOTE** This needs to be changed by the user
## -------------------------------------------------------------------- ##
maxent.file <- "~/Dropbox/Scripts-Programs/Maxent/maxent.jar"

## Set directory for environmental variables
## ***
## In this application we are using the variables used in Latimer et al. 
## 2006, and other projects from this group. This include climate layers
## from the Schulze dataset and environmental layers derived from soil
## maps. In total there were 24 variables
## -------------------------------------------------------------------- ##
environmental.layers <- 
  "SA_ASCII/"

raster_template <- raster( "SA_ASCII/alt.asc" )

## Compile default maxent command
## -------------------------------------------------------------------- ##
run.default <- paste0( "java -jar ", maxent.file, " nowarnings noprefixes ",
                       "-e ", environmental.layers, " nowarnings threads=1 ", 
                       " nowriteclampgrid nowritemess noresponsecurves  outputformat=raw -a -z ",
                       "noaskoverwrite ",
                       "nothreshold nohinge noautofeature noproduct " ) # Model features to exclude

## ******************************************************************** ##
## FUNCTION: normalize_ascii
## ******************************************************************** ##

normalize_ascii <- function( x ){
  values( x ) <-
    values( x ) / sum( values( x ), na.rm = TRUE )
  
  return( x)
}

## ******************************************************************** ##
## FUNCTION: minxent
## Calculate minxent layer using maxent output (from run with bias
## layer) and prior (i.e., the bias layer)
## ******************************************************************** ##

minxent <- function( prior.asc, maxent.asc ){
  
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
## Make range map
## -------------------------------------------------------------------- ##
## This step can be skipped if the user has a digitized expert range 
## map. Here, we create an "expert" map using alpha convex hulls.
## ******************************************************************** ##
## ******************************************************************** ##

## Read jittered Protea punctata Atlas locations
prpunc <- read.csv( "prpunc_atlas_jittered.csv" )

## Construct an expert range map from the Protea Atlas dataset occurence
## locations using an alpha hull with an alpha value of 0.5
prpunc_atlas_ahull <- 
  ahull( x = unique( prpunc ), alpha = 0.5 )

## Convert from ahull object to spatial polygon. For this step
## I am using a function submitted to the R-SIG-GEO group by
## Anrew Bevan.
source( "alphahull_to_shape.R")
prpunc_atlas_polygon <- ah2sp( x = prpunc_atlas_ahull )

## Rasterize the alpha hull polygon to use it as a range map
prpunc_atlas_raster <- 
  rasterize( prpunc_atlas_polygon, y = raster_template )

## ******************************************************************** ##
## ******************************************************************** ##
## Assign probability values to range map
## ***
## Here we are assigning ROR values to the prior based on our 
## a priori assumptions regarding how much probability should be
## assigned to areas classified as occupied versus those not classified
## as occupied in the range map.
## We have chosen to assing 90% of the probability to areas classified
## as occupied, and 10% to other areas.
## ******************************************************************** ##
## ******************************************************************** ##

## Set all points in the prpunc_atlas_prior raster to 0.1
## Note that we could've chosen any number this (not just 0.1) because 
## we assign the same value to all cells below.
prpunc_atlas_prior <- raster_template
values( prpunc_atlas_prior ) <- 
  ifelse( values( prpunc_atlas_prior ) > 0, yes = 0.1, no = NA )

## Set those points that are > 0 in the prpunc_atlas_raster file to 0.9. 
## Note that we could've chosen any number this (not just 0.9) because 
## we assign the same value to all cells below.
values( prpunc_atlas_prior )[ which( values( prpunc_atlas_raster ) > 0 ) ] <- 0.9

## Set the grid cell values were p = 0.9 such that the sum of all
## of these grid cells is equal to 0.9. Do the same for cells
## with values of 0.1
grids_0.9 <- sum( values( prpunc_atlas_prior ) == 0.9, na.rm = TRUE )
grids_0.1 <- sum( values( prpunc_atlas_prior ) == 0.1, na.rm = TRUE )

prpunc_atlas_prior[ which( values( prpunc_atlas_prior ) == 0.9 ) ] <- ( 0.9 / grids_0.9 )
prpunc_atlas_prior[ which( values( prpunc_atlas_prior ) == 0.1 ) ] <- ( 0.1 / grids_0.1 )

## Write this raster to file
writeRaster( x = prpunc_atlas_prior, 
             filename = "prpunc_atlas_prior.asc", overwrite = TRUE )

## ******************************************************************** ##
## ******************************************************************** ##
## Construct a naive Maxent model using the relatively sparse PRECIS
## dataset
## ***
## Note that this dataset is comporable to what is available via GBIF;
## the herbarium specimens in the PRECIS dataset are served on GBIF as
## well.
## ******************************************************************** ##
## ******************************************************************** ##

if( !file.exists( "Maxent_output/PRPUNC_PRECIS" ) ) { dir.create( "Maxent_output/PRPUNC_PRECIS" ) }

## Run the naieve maxent model using the PRECIS data
system( paste0( run.default, " -o ", "Maxent_output/PRPUNC_PRECIS",
                " -s ", "prpunc_precis_train.csv" ) ) 

## Read in the resulting ASCII map, which will be used later
prpunc_precis_asc <- raster( "Maxent_output/PRPUNC_PRECIS/prpunc.asc" )

## Normalize prpunc_precis_asc
prpunc_precis_asc <- normalize_ascii( prpunc_precis_asc )


## ******************************************************************** ##
## ******************************************************************** ##
## Construct a Maxent model for Protea punctata, using the Protea Atlas
## data as a range map prior
## ******************************************************************** ##
## ******************************************************************** ##

if( !file.exists( "Maxent_output/PRPUNC_PRECIS_bias" ) ) { dir.create( "Maxent_output/PRPUNC_PRECIS_bias" ) }

system( paste0( run.default, " -o ", "Maxent_output/PRPUNC_PRECIS_bias",
                " -s ", "prpunc_precis_train.csv",
                " biasfile=", "prpunc_atlas_prior.asc",
                " biastype=3 ") ) 

## Read model with bias layer removed
prpunc_precis_rm_bias_asc <- 
  raster( "Maxent_output/PRPUNC_PRECIS_bias/prpunc.asc" )

## Normalize this ascii
prpunc_precis_rm_bias_asc <- normalize_ascii( x = prpunc_precis_rm_bias_asc )

## Create Minxent layer
prpunc_precis_minxent <- minxent( prior.asc = prpunc_atlas_prior, 
                                  maxent.asc = prpunc_precis_rm_bias_asc )

## Write an ascii file for this final layer
writeRaster( x = prpunc_precis_minxent, 
             filename = "prpunc_minxent.asc", overwrite = TRUE )

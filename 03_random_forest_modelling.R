#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# RandomForest Modelling of change agents

# This script will read the human interpretted plots, assign a change polygon
# to each plot. Merge with the predictors extracted in script 01 and fit random
# forest model to map chabge agents across the AOI.

# Jonathan Dash

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# To do 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# * Increase sample size - Now doubled to 200
# * Add the factor simplification for the response variables
# * Predict for all polygons
# * OOB error and confusion matrix
# * Independent validation
# * Summarise the outputs 

#### Setup libraries ####
library(sf)
library(raster)
library(tidyverse)
library(here)
library(caret)
library(ranger)


#### Set up file directories ####
print('Setting up file directories')
# Set up the directory structure
# Each ftv folder contains a collection of rasters showing the fitted to vertex values for 
# t-cap spectral properties
tcb_ftv_dir<-here('raster', 'ftvRasters', 'tcb_ftv')
tcg_ftv_dir<-here('raster', 'ftvRasters', 'tcg_ftv')
tcw_ftv_dir<-here('raster', 'ftvRasters', 'tcw_ftv')
chg_map_dir<-here('raster', 'changeMaps')

# List the tif files in the ftv directories 
tcb_ftv_files <- list.files(tcb_ftv_dir, pattern = "\\.tif$", full.names = TRUE)
tcg_ftv_files <- list.files(tcg_ftv_dir, pattern = "\\.tif$", full.names = TRUE)
tcw_ftv_files <- list.files(tcw_ftv_dir, pattern = "\\.tif$", full.names = TRUE)

print('File directories setup')


#### Read data ####
print('Reading Data')

# Read Gain poly
gain_shp<-st_read(here('vector', 'lt_gain_poly_v4.shp'))
#plot(gain_shp)

# Read loss poly
loss_shp<-st_read(here('vector', 'lt_disturbance_poly.shp'))

# Read interpretations
ref<-read_csv(here('AttributionPlots', 'samplePlotsAttribution.csv'))
class(ref)

# Read predictors
gain.preds<-readRDS(here('out', 'lt_gain.predictors.rds'))
loss.preds<-readRDS(here('out', 'lt_loss.predictors.rds'))

print('Finished reading Data')


#### Analysis ####

# Summarise the attribution plots
ref %>% group_by(changeAgent) %>%
  tally()

# Convert attribution table to sf and specify CRS
ref<- st_as_sf(ref, coords = c('X', 'Y'), crs= 2193)
class(ref)


# Simplify the change agent

ref<-ref %>% mutate (simplifiedAgent =  ifelse(changeAgent %in% c(0 ,14), 0, # No change
                                               ifelse(changeAgent %in% c(1, 11, 12),1, # Natural Disturbance
                                                      ifelse(changeAgent %in% c(2,3), 2, # Wilding control
                                                             ifelse(changeAgent %in% c(4,5,6,7), 3, # Forest management
                                                                    ifelse(changeAgent == 8,4,
                                                                           ifelse(changeAgent == 9,5,
                                                                                  ifelse(changeAgent == 10,6,
                                                                                         ifelse(changeAgent == 13,7,
                                                                                                'OTHER')))))))))



# Summarise the attribution plots
ref %>% group_by(simplifiedAgent) %>%
  tally()

plot(ref)

# Separate out loss and gain plots
ref.loss<-subset(ref, changeType == 'Loss')
ref.gain<-subset(ref, changeType == 'Gain')

ref.loss<-st_intersection(ref.loss, loss_shp)
st_geometry(ref.loss) <- NULL
ref.gain<-st_intersection(ref.gain, gain_shp)
st_geometry(ref.gain) <- NULL

# Add predictors into reference  dataset
ref.gain<- left_join(ref.gain, gain.preds, by = 'fid')
ref.loss<- left_join(ref.loss, loss.preds, by = 'fid')



ref.gain<- ref.gain %>% 
  select(-id, -changeType, -Comments, -fid, -cat, -ID, -label, -yod, -geometry, -changeAgent) 

ref.loss<- ref.loss %>% 
  select(-id, -changeType, -Comments, -fid, -cat, -ID, -label, -yod, -geometry, -changeAgent) 







#ref.gain<-as.data.frame(ref.gain)
class(ref.gain)
ref.gain$simplifiedAgent<-as.factor(ref.gain$simplifiedAgent)
ref.loss$simplifiedAgent<-as.factor(ref.loss$simplifiedAgent)


# Fit a random forest model for gain and loss




set.seed(1979)
gain_rf<-ranger(simplifiedAgent~., data=ref.gain, num.trees = 500, mtry =1, importance = "permutation")
print(gain_rf)
gain_rf$confusion.matrix

set.seed(1979)
loss_rf<-ranger(simplifiedAgent~., data=ref.loss, num.trees = 500, mtry =1)
print(loss_rf)
loss_rf$confusion.matrix










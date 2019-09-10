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


# * Increase sample size - Now doubled to 200 other selectively identified added
# * Add the factor simplification for the response variables
# * Predict for all polygons - Done
# * OOB error and confusion matrix - Done
# * Independent validation - TODO
# * Summarise the outputs 

#### Setup libraries ####
library(sf)
library(raster)
library(tidyverse)
library(here)
library(caret)
library(ranger)
library(randomForest)


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
gain_shp<-st_read(here('vector', 'lt_gain_poly_v6.shp'))
#plot(gain_shp)

# Read loss poly
loss_shp<-st_read(here('vector', 'lt_disturbance_poly_v6.shp'))

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
  select(-id, -changeType, -Comments, -fid, -cat, -ID, -label, -yod, -geometry, -changeAgent, -dsnr, -dur, -rate) 

ref.loss<- ref.loss %>% 
  select(-id, -changeType, -Comments, -fid, -cat, -ID, -label, -yod, -geometry, -changeAgent, -dsnr, -dur, -rate) 







#ref.gain<-as.data.frame(ref.gain)
class(ref.gain)
ref.gain$simplifiedAgent<-as.factor(ref.gain$simplifiedAgent)
ref.loss$simplifiedAgent<-as.factor(ref.loss$simplifiedAgent)


# Fit a random forest model for gain and loss


set.seed(1979)
gain_rf<-ranger(simplifiedAgent~., data=ref.gain, num.trees = 500, mtry =1, importance = "permutation")
print(gain_rf)
gain_rf$confusion.matrix
gain_rf$variable.importance

set.seed(1979)
gain_rf2<-randomForest(simplifiedAgent~., data=ref.gain,  norm.votes = TRUE, proximity = TRUE)
print(gain_rf2)
plot(gain_rf2)

gain_rf2$confusion

imp.gain<-as.data.frame(importance(gain_rf2))
imp.gain$Variable<-row.names(imp.gain)


set.seed(1979)
loss_rf<-ranger(simplifiedAgent~., data=ref.loss, num.trees = 500, mtry =1)
print(loss_rf)
loss_rf$confusion.matrix

print(1979)
loss_rf2<-randomForest(simplifiedAgent~., data=ref.loss,  norm.votes = TRUE, proximity = TRUE)
print(loss_rf2)
plot(loss_rf2)

loss_rf2$confusion

imp.loss<-as.data.frame(importance(loss_rf2))
imp.loss$Variable<-row.names(imp.loss)

imps<-left_join(imp.gain, imp.loss, by="Variable")




#### Use RF model to predict for the all polygons ####

gain.out.prob<-as.data.frame(predict(gain_rf2, gain.preds, type = "prob"))
loss.out.prob<-as.data.frame(predict(loss_rf2, loss.preds, type = "prob"))
gain.out.resp<-as.data.frame(predict(gain_rf2, gain.preds, type = "response"))
loss.out.resp<-as.data.frame(predict(loss_rf2, loss.preds, type = "response"))



gain.map<-bind_cols(gain.preds, gain.out.resp)
loss.map<-bind_cols(loss.preds, loss.out.resp)

colnames(gain.map)[24] <- "ChangeClass"
colnames(loss.map)[24] <- "ChangeClass"


# Chart wilding invasion by year

gain.map$ChangeClass<-fct_explicit_na(gain.map$ChangeClass)
loss.map$ChangeClass<-fct_explicit_na(loss.map$ChangeClass)

gain.areas<-gain.map %>% 
  #group_by(yod, ChangeClass) %>% 
  #summarize() %>%
  mutate(area = st_area(.))  

loss.areas<-loss.map %>% 
  #group_by(yod, ChangeClass) %>% 
  #summarize() %>%
  mutate(area = st_area(.)) 


gain.areas$area<-unclass(gain.areas$area)
loss.areas$area<-unclass(loss.areas$area)



# Chart wilding gain only
gain.areas %>%
  filter(ChangeClass == 4) %>% 
  group_by(yod) %>%
  summarise(area2 = sum(area)) %>%
  ggplot(aes(x=yod, y=cumsum(area2)/10000)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(y = "Cumulative area (ha)", x= 'Year')


# Chart wilding control only
loss.areas %>%
  filter(ChangeClass == 2) %>% 
  group_by(yod) %>%
  summarise(area2 = sum(area)) %>%
  ggplot(aes(x=yod, y=cumsum(area2)/10000)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(y = "Cumulative area (ha)", x= 'Year')



# Stacked chart of gain attribution
gain.areas %>%
  group_by(yod, ChangeClass) %>%
  summarise(area2 = sum(area)) %>%
  ggplot(aes(x=yod, y=area2/10000, fill = ChangeClass)) +
  geom_area(alpha=0.6 , size=.5, colour="white") +
  scale_fill_viridis_d() +
  theme_bw()


# Stacked chart of loss attribution
loss.areas %>%
  group_by(yod, ChangeClass) %>%
  summarise(area2 = sum(area)) %>%
  ggplot(aes(x=yod, y=area2/10000, fill = ChangeClass)) +
  geom_area(alpha=0.6 , size=.5, colour="white") +
  scale_fill_viridis_d() +
  theme_bw()
  
  



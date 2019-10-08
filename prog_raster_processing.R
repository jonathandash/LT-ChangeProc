# Prepare raster data for patch level RF modelling to attribute the causes of change
# Jonathan Dash 
# September 2019

#-----
# Step 1

# Produce vrts for the ftv t-cap rasters. No need for this now


# Step 2
#
# Polygonise the change maps. This was done in QGIS


#Step3
# Add the landtrendr segment values into the Change maps

#----
#Step 4
# Get the t-cap spectral values at the start and end of each change


##### Setup libraries ####
library(rgdal)
library(raster)
library(sf)
library(here)
library(pryr)
library(gdalUtils)
library(stars)
library(velox)
library(tidyverse)



#### Define functions ####

# Function to extract raster values using velox library
extractor<-function(rast,shp) {
  
  vx <- velox(rast) #Convert raster to velox data
  x  <- vx$extract(sp = shp, fun = function(x) mean(x, na.rm = TRUE)) #Extract mean raster value per poly
  rm(vx) # remove velox data
  xdf<-as.data.frame(x)
  tmp<-bind_cols(shp, xdf)
  return(tmp) # return data frame
  
  
}

nameChangeLT<-function(df){
  
  colnames(df)[which(names(df) == "V1")] <- "yod"
  colnames(df)[which(names(df) == "V2")] <- "mag"
  colnames(df)[which(names(df) == "V3")] <- "dur"
  colnames(df)[which(names(df) == "V4")] <- "preval"
  colnames(df)[which(names(df) == "V5")] <- "rate"
  colnames(df)[which(names(df) == "V6")] <- "dsnr"
  return(df)
  
}

# This function will go through a directory of rasters and extract the mean
# pixel value within each polygon of a shapefile.
# The function wants a directory with some rasters and a sf.
tileextractor<-function(dir, shp)
{
  out<-data.frame()
  for (i in 1:length(dir))
    #i=6
  {
    ms<-stack(dir[i]) # read in raster
    tc<- extractor(rast = ms, shp = shp) #extract values
    tcdf<-as.data.frame(tc)
    tmp<-bind_cols(shp, tcdf) %>%
      drop_na(c(V1, V8, V15, V15, V19))
    out<-bind_rows(out, tmp)
  }
  return(out)
}


# Simple function to provide the oppsite of %in%
`%notin%` = function(x,y) !(x %in% y)



#### Set up file directories ####
print('Setting up file directories')
# Set up the directory structure
# Each ftv folder contains a collection of rasters showing the fitted to vertex values for 
# t-cap spectral properties
tcb_ftv_dir<-here('raster', 'ftvRasters', 'tcb_ftv')
tcg_ftv_dir<-here('raster', 'ftvRasters', 'tcg_ftv')
tcw_ftv_dir<-here('raster', 'ftvRasters', 'tcw_ftv')
nbr_ftv_dir<-here('raster', 'ftvRasters', 'nbr_ftv')

chg_map_dir<-here('raster', 'changeMaps')

# List the tif files in the ftv directories 
tcb_ftv_files <- list.files(tcb_ftv_dir, pattern = "\\.tif$", full.names = TRUE)
tcg_ftv_files <- list.files(tcg_ftv_dir, pattern = "\\.tif$", full.names = TRUE)
tcw_ftv_files <- list.files(tcw_ftv_dir, pattern = "\\.tif$", full.names = TRUE)
nbr_ftv_files <- list.files(nbr_ftv_dir, pattern = "\\.tif$", full.names = TRUE)

print('File directories setup')


#### Read data ####



print('Reading change polygons')

# Read Gain poly
gain_shp_cons<-st_read(here('vector', 'gain_cons.shp'))
gain_shp_bal<-st_read(here('vector', 'gain_bal.shp'))
gain_shp_unlim<-st_read(here('vector', 'gain_unlim.shp'))
#gain_shp<-st_read(here('vector', 'lt_gain_poly_v7.shp'))
#plot(gain_shp)

# Read loss poly
loss_shp_cons<-st_read(here('vector', 'disturbance_cons.shp'))
loss_shp_bal<-st_read(here('vector', 'disturbance_bal.shp'))
loss_shp_unlim<-st_read(here('vector', 'disturbance_unlim.shp'))



print('Polygons read')



#### Read rasters####

print('Reading LT rasters')
# Read in lt outputs
#lt_gain<-stack(here('raster', 'changeMaps', 'lt-gee_gain_map_HighCountry_v3.tif'))
lt_gain<-stack(here('raster', 'changeMaps', 'lt-gee_gain_map_HighCountry_Unlim.tif'))
NAvalue(lt_gain) <- 0

#lt_gain_bal<-stack(here('raster', 'changeMaps', 'lt-gee_gain_map_HighCountry_Bal.tif'))
#NAvalue(lt_gain_bal) <- 0

#lt_gain_unlim<-stack(here('raster', 'changeMaps', 'lt-gee_gain_map_HighCountry_Unlim.tif'))
#NAvalue(lt_gain_unlim) <- 0

lt_loss<-stack(here('raster', 'changeMaps', 'lt-gee_disturbance_map_HighCountry_Unlim.tif'))
NAvalue(lt_loss) <- 0

#lt_loss_bal<-stack(here('raster', 'changeMaps', 'lt-gee_disturbance_map_HighCountryBal.tif'))
#NAvalue(lt_loss_bal) <- 0

#lt_loss_unlim<-stack(here('raster', 'changeMaps', 'lt-gee_disturbance_map_HighCountry_Unlim.tif'))
#NAvalue(lt_loss_unlim) <- 0



print('LT rasters read')

#### Raster processing ####

print('Extracting LT rasters to polygons')
# Apply extractor function to the change maps

lt_gain.out<-extractor(rast = lt_gain, shp = gain_shp_unlim)
lt_gain.out<-nameChangeLT(df = lt_gain.out)

#lt_gain.out.bal<-extractor(rast = lt_gain_bal, shp = gain_shp_bal)
#lt_gain.out.bal<-nameChangeLT(df = lt_gain.out.bal)

#lt_gain.out.unlim<-extractor(rast = lt_gain_unlim, shp = gain_shp_unlim)
#lt_gain.out.unlim<-nameChangeLT(df = lt_gain.out.unlim)


lt_loss.out <- extractor(rast = lt_loss, shp = loss_shp_unlim)
lt_loss.out<-nameChangeLT(df = lt_loss.out)

#lt_loss.out.bal <- extractor(rast = lt_loss_bal, shp = loss_shp_bal)
#lt_loss.out.bal<-nameChangeLT(df = lt_loss.out.bal)

#lt_loss.out.unlim <- extractor(rast = lt_loss_unlim, shp = loss_shp_unlim)
#lt_loss.out.unlim<-nameChangeLT(df = lt_loss.out.unlim)

print('Finished extracting raster polygons')

# Extract values for the gain polygons
# These are producing a small number of duplicates... I think this is because of
# a small number of cases where a polygon overlaps 2 rasters. 
# I don't think this matters and will deal with it later if needed.
# Turns out I need to remove it. This can be done with the call to dplyr::distinct
# This keeps the first row of the duplicate.

print('Extracting tcb raster tiles to polygons')

#gain.shps<-c('gain_shp_cons', 'gain_shp_bal', 'gain_shp_unlim')
#loss.shps<-c('loss_shp_cons', 'loss_shp_bal', 'loss_shp_unlim')








tcb.gain<-tileextractor(dir = tcb_ftv_files, shp = gain_shp_unlim )
tcb.gain$tCap<-'brightness'
tcb.gain$scenario<-'gain_shp_unlim'
tcb.gain<-distinct(tcb.gain, fid, .keep_all= TRUE)


tcw.gain<-tileextractor(dir = tcw_ftv_files, shp = gain_shp_unlim)
tcw.gain$tCap<-'wetness'
tcw.gain$scenario<-'gain_shp_unlim'
tcw.gain<-distinct(tcw.gain, fid, .keep_all= TRUE)


tcg.gain<-tileextractor(dir = tcg_ftv_files, shp = gain_shp_unlim)
tcg.gain$tCap<-'greeness'
tcg.gain$scenario<-'gain_shp_unlim'
tcg.gain<-distinct(tcg.gain, fid, .keep_all= TRUE)


nbr.gain<-tileextractor(dir = nbr_ftv_files, shp = gain_shp_unlim)
nbr.gain$tCap<-'nbr'
nbr.gain$scenario<-'gain_shp_unlim'
nbr.gain<-distinct(nbr.gain, fid, .keep_all= TRUE)





#Extract values for the loss polygons
tcb.loss<-tileextractor(dir = tcb_ftv_files, shp = loss_shp_unlim)
tcb.loss$tCap<-'brightness'
tcb.loss<-distinct(tcb.loss, fid, .keep_all= TRUE)
tcb.loss$scenario<-'loss_shp_unlim'
tcw.loss<-tileextractor(dir = tcw_ftv_files, shp = loss_shp_unlim)
tcw.loss$tCap<-'wetness'
tcw.loss<-distinct(tcw.loss, fid, .keep_all= TRUE)
tcw.loss$scenario<-'loss_shp_unlim'
tcg.loss<-tileextractor(dir = tcg_ftv_files, shp = loss_shp_unlim)
tcg.loss$tCap<-'greeness'
tcg.loss<-distinct(tcg.loss, fid, .keep_all= TRUE)
tcg.loss$scenario<-'loss_shp_unlim'
nbr.loss<-tileextractor(dir = nbr_ftv_files, shp = loss_shp_unlim)
nbr.loss$tCap<-'nbr'
nbr.loss<-distinct(nbr.loss, fid, .keep_all= TRUE)
nbr.loss$scenario<-'loss_shp_unlim'

print('Finished extracting tcb raster tiles')



#### Preparing predictor dataframe ####


print('Munging predictors dataframe')
# Now make the dataframe as a lookup for the predictor set
# Merge the dataframes
tc.gain<-bind_rows(tcb.gain, tcw.gain, tcg.gain, nbr.gain)
tc.loss<-bind_rows(tcb.loss, tcw.loss, tcg.loss, nbr.loss)

# Specify the names to be changed
oldnames = c("V1","V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11",
             "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20")
newnames = c('2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', 
             '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', 
             '2016', '2017', '2018', '2019')
# Change the column names
tc.gain<-tc.gain %>% rename_at(vars(oldnames), ~ newnames)
tc.loss<-tc.loss %>% rename_at(vars(oldnames), ~ newnames)

# I am writing these out here because they take a long time to generate.
saveRDS(tc.gain, here('out', 'tc_gain_unlim.rds'))
saveRDS(tc.loss, here('out', 'tc_loss_unlim.rds'))

tc.gain<-readRDS(here('out', 'tc_gain_unlim.rds'))
tc.loss<-readRDS(here('out', 'tc_loss_unlim.rds'))


# Make a lookup table to add values into the predictor set
tc.gain.lookup<-tc.gain %>% 
  select(-fid1, -cat1, -ID1, -label1, -geometry1, -scenario) %>%
  gather(key = year, value= value, -fid, -cat, -ID, -label, -geometry, -tCap) %>%
  select(fid, year, tCap, value) %>%
  spread(key=tCap, value = value) %>%
  mutate(year = as.numeric(year))

tc.loss.lookup<-tc.loss %>% 
  select(-fid1, -cat1, -ID1, -label1, -geometry1, -scenario) %>%
  gather(key = year, value= value, -fid, -cat, -ID, -label, -geometry, -tCap) %>%
  select(fid, year, tCap, value) %>%
  spread(key=tCap, value = value) %>%
  mutate(year = as.numeric(year))


# Now make the datasets for prediction

oldnames = c("greeness", "brightness", "wetness", "nbr")
newnames = c("start_greeness", "start_brightness", "start_wetness", "start_nbr")

lt_gain.predictors<-left_join(lt_gain.out, tc.gain.lookup, by = c("fid" = "fid",  "yod" = "year")) 
lt_gain.predictors<-lt_gain.predictors %>% rename_at(vars(oldnames), ~ newnames)

lt_loss.predictors<-left_join(lt_loss.out, tc.loss.lookup, by = c("fid" = "fid",  "yod" = "year"))
lt_loss.predictors<-lt_loss.predictors %>% rename_at(vars(oldnames), ~ newnames)

# Calculate the end date of the segment
lt_gain.predictors$dur<-round(lt_gain.predictors$dur, 0)
lt_loss.predictors$dur<-round(lt_loss.predictors$dur, 0)

lt_gain.predictors$yod_end<-ifelse(lt_gain.predictors$yod + lt_gain.predictors$dur>2019,2019, lt_gain.predictors$yod + lt_gain.predictors$dur)
lt_loss.predictors$yod_end<-ifelse(lt_loss.predictors$yod + lt_loss.predictors$dur>2019,2019, lt_loss.predictors$yod + lt_loss.predictors$dur)

# Get one year before detection
lt_gain.predictors$yod_pre <-lt_gain.predictors$yod - 1
lt_loss.predictors$yod_pre <-lt_loss.predictors$yod - 1

# Get one year after detection
lt_gain.predictors$yod_post <-ifelse(lt_gain.predictors$yod == 2019, 2019, lt_gain.predictors$yod + 1)
lt_loss.predictors$yod_post <-ifelse(lt_loss.predictors$yod == 2019, 2019, lt_loss.predictors$yod + 1)

# bring in the t-cap values for the end of the segment
newnames = c("end_greeness", "end_brightness", "end_wetness", "end_nbr")

lt_gain.predictors<-left_join(lt_gain.predictors, tc.gain.lookup, by = c("fid" = "fid",  "yod_end" = "year")) 
lt_gain.predictors<-lt_gain.predictors %>% rename_at(vars(oldnames), ~ newnames)

lt_loss.predictors<-left_join(lt_loss.predictors, tc.loss.lookup, by = c("fid" = "fid",  "yod_end" = "year")) 
lt_loss.predictors<-lt_loss.predictors %>% rename_at(vars(oldnames), ~ newnames)


# bring in the t-cap values for one year before 
newnames = c("pre_greeness", "pre_brightness", "pre_wetness", "pre_nbr")
lt_gain.predictors<-left_join(lt_gain.predictors, tc.gain.lookup, by = c("fid" = "fid",  "yod_pre" = "year")) 
lt_gain.predictors<-lt_gain.predictors %>% rename_at(vars(oldnames), ~ newnames)

lt_loss.predictors<-left_join(lt_loss.predictors, tc.loss.lookup, by = c("fid" = "fid",  "yod_pre" = "year")) 
lt_loss.predictors<-lt_loss.predictors %>% rename_at(vars(oldnames), ~ newnames)


# bring in the t-cap values for one year after detection
newnames = c("post_greeness", "post_brightness", "post_wetness", "post_nbr")
lt_gain.predictors<-left_join(lt_gain.predictors, tc.gain.lookup, by = c("fid" = "fid",  "yod_post" = "year")) 
lt_gain.predictors<-lt_gain.predictors %>% rename_at(vars(oldnames), ~ newnames)

lt_loss.predictors<-left_join(lt_loss.predictors, tc.loss.lookup, by = c("fid" = "fid",  "yod_post" = "year")) 
lt_loss.predictors<-lt_loss.predictors %>% rename_at(vars(oldnames), ~ newnames)


# Calculate segment delta for t-cap values and remove unwanted columns
lt_loss.predictors<-lt_loss.predictors %>% 
  mutate(delta_brightness = end_brightness - start_brightness,
         delta_greeness = end_greeness - start_greeness,
         delta_wetness = end_wetness - start_wetness,
         delta_nbr = end_nbr - start_nbr) %>%
  select(-cat, -ID, -label, -yod_end, -yod_pre, -yod_post)

lt_gain.predictors<-lt_gain.predictors %>% 
  mutate(delta_brightness = end_brightness - start_brightness,
         delta_greeness = end_greeness - start_greeness,
         delta_wetness = end_wetness - start_wetness,
         delta_nbr = end_nbr - start_nbr) %>%
  select(-cat, -ID, -label, -yod_end, -yod_pre, -yod_post)

# Save the predictors dataframes
saveRDS(lt_gain.predictors, here('out', 'lt_gain.predictors_bal.rds'))
saveRDS(lt_loss.predictors, here('out', 'lt_loss.predictors_bal.rds'))

print('Predictors ready')





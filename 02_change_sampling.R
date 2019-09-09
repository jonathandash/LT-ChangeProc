#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Change sampling

# This script will provide a sample the human interpretted plots within the change polygons
# We use stratified random sampling.
# The strata are just disturbance and gain.

# Jonathan Dash


#### Setup libraries ####
library(sf)
library(raster)
library(tidyverse)
library(here)


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

print('Reading change polygons')

# Read Gain poly
gain_shp<-st_read(here('vector', 'lt_gain_poly_v4.shp'))
#plot(gain_shp)

# Read loss poly
loss_shp<-st_read(here('vector', 'lt_disturbance_poly.shp'))

HCR<-st_read(here('vector', 'High_CountryLCR_bdy.shp'))
HCR<- HCR %>% st_transform(crs = 2193)

print('Polygons read')


# Read the manually identified wilding locations to add to the sample

wld.locs<-st_read(here('vector', 'WildingSamples.kml')) %>% 
  st_transform(crs = 2193) %>%
  select(-Name) %>%
  mutate(id = paste('wld_', row_number()+400, sep=''))



wld.locs %>%
  st_write(here('out', 'WildingLocations.kml'))

xy<-st_coordinates(wld.locs)

wld.locs<-cbind(wld.locs, xy)

wld.locs %>%
  write_csv(here('out', 'wld.locs.csv')) 


#### Start sampling ####

# Choose a sample size per strata and a sampling strategy


#set.seed(1977)
set.seed(1987)

n= 100 # Sample size per stratum
samp.strat = "random" # Sampling type

# Place plots at random
loss.plts<-st_sample(loss_shp, n, samp.strat, exact = TRUE)

gain.plts<-st_sample(gain_shp, n, samp.strat, exact = TRUE)


gain.plts.out<-data.frame(st_coordinates(gain.plts)) %>% 
  mutate(id = paste('gain_', row_number()+100, sep='')) %>% 
  dplyr::select(id, X, Y)


loss.plts.out<-data.frame(st_coordinates(loss.plts)) %>% 
  mutate(id = paste('loss_', row_number()+100, sep='')) %>% 
  dplyr::select(id, X, Y)



plts.out<-bind_rows(gain.plts.out,  loss.plts.out) %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(2193)

class(plts.out)

st_coordinates(plts.out)



plts.out %>%
  st_write(here('out', 'samplePlots2.kml'))

plts.out %>%
  st_write(here('out', 'samplePlots2.shp'))

#plts.out<-st_transform(plts.out, 4326) # This will go to lat long

#plts.out<-st_transform(plts.out, 2193)

xy<-st_coordinates(plts.out)

plts.out<-cbind(plts.out, xy)

plts.out %>%
  write_csv(here('out', 'samplePlots2.csv')) 


ggplot() +
  geom_sf(data = HCR, colour = "red", fill = NA)+
  #scale_colour_manual(values= c25) +
  geom_sf(data = loss.plts, colour="red")+
  geom_sf(data = gain.plts, colour="darkblue")+
  #labs(colour = 'Strata')+
  theme_bw()





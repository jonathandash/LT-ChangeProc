REM Export a single band showing the yod and set the null value.
REM Run from OSGEO4W shell

gdal_translate -a_nodata 0 -co compress=lzw lt-gee_gain_map_HighCountry_v7.tif -b 1 lt_yod_gee_gain_map_HighCountry_v7.tif

gdal_translate -a_nodata 0 -co compress=lzw lt-gee_disturbance_map_HighCountry_v7.tif -b 1 lt_yod_gee_disturbance_map_HighCountry_v7.tif 

gdal_translate -a_nodata 0 -co compress=lzw lt-gee_gain_map_HighCountry_Unlim.tif -b 1 lt_yod_gee_gain_map_HighCountry_Unlim.tif

gdal_translate -a_nodata 0 -co compress=lzw lt-gee_disturbance_map_HighCountry_Unlim.tif -b 1 lt_yod_gee_disturbance_map_HighCountry_Unlim.tif 
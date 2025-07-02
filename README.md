
# downscaleToPoint

This package contains the necessary code and data for reproducing all results in the paper "Simulation and evaluation of local daily temperature and precipitation series derived by stochastic downscaling of ERA5 reanalysis".

The scripts used for creating all relevant results and figures are available in the `exec/`
folder. This folder contains six scripts:

- `plot_station_data.R`:  
  This script is used for creating a map that showcases all available GSOD weather stations that we
  use for performing our downscaling.
- `downscale_precipitation.R`:  
  This script trains the precipitation downscaling models, performs leave-one-out cross-validation
  for all relevant precipitation weather stations, and summarises the results of the cross-validation
  study using a collection of different figures.
- `downscale_temperature.R`:  
  This script trains the temperature downscaling models, performs leave-one-out cross-validation
  for all relevant temperature weather stations, and summarises the results of the cross-validation
  study using a collection of different figures.
  

The folder `raw_data/` contains all relevant weather station data, ERA5 data and DEM data used for
performing our downscaling. We ask that you please reference the paper "Simulation and evaluation of local daily temperature and precipitation series derived by stochastic downscaling of ERA5 reanalysis" if you use this data.


# downscaleToPoint

This package contains the necessary code and data for reproducing all results in the paper "Simulation and evaluation of local daily temperature and precipitation series derived by stochastic downscaling of ERA5 reanalysis", available at [https://arxiv.org/abs/2507.01692](https://arxiv.org/abs/2507.01692).

The scripts used for creating all relevant results and figures are available in the `exec/`
folder. This folder contains nine scripts:

- `01-plot_station_data.R`:  
  This script is used for creating a map that showcases all available GSOD weather stations that we
  use for performing our downscaling.
- `02-downscale_precipitation.R`:  
  This script trains the precipitation downscaling models, performs leave-one-out cross-validation
  for all relevant precipitation weather stations, and summarises the results of the cross-validation
  study using a collection of different figures.
- `03-downscale_temperature.R`:  
  This script trains the temperature downscaling models, performs leave-one-out cross-validation
  for all relevant temperature weather stations, and summarises the results of the cross-validation
  study using a collection of different figures.
- `04-cprcm_precip_comparison.R`:  
  This script repeats the model evaluation from `02-downscale_precipitation.R`, but for a smaller spatial domain, and by comparing our stochastic downscaling model to a convection-permitting regional climate model (CPRCM) instead of ERA5.
- `05-cprcm_temp_comparison.R`:  
  This script repeats the model evaluation from `03-downscale_temperature.R`, but for a smaller spatial domain, and by comparing our stochastic downscaling model to a convection-permitting regional climate model instead of ERA5.
- `06-precip_spatial_consistency.R`:  
  This script is used to examine the spatial properties of the full precipitation downscaling model, when used to create spatiotemporal ensembles of local weather at multiple different locations.
- `07-temperature_spatial_consistency.R`:  
  This script is used to examine the spatial properties of the full temperature downscaling model, when used to create spatiotemporal ensembles of local weather at multiple different locations.
- `08-cprcm_precip_spatial_consistency.R`:  
  This script repeats the model evaluation from `06-precip_spatial_consistency.R`, but for a smaller spatial domain, and by comparing our stochastic downscaling model to a CPRCM instead of ERA5.
- `09-cprcm_temp_spatial_consistency.R`:  
  This script repeats the model evaluation from `06-temperature_spatial_consistency.R`, but for a smaller spatial domain, and by comparing our stochastic downscaling model to a CPRCM instead of ERA5.


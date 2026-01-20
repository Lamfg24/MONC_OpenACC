#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 19:14:45 2025

@author: delbeke
"""

import netCDF4 as nc
import numpy as np

fn = '/home/delbeke/Documents/vn1.0.4_monc_bellouin_delbeke/surface_forcing.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')

time = ds.createDimension('time', 3)

times = ds.createVariable('time', np.float64, ('time',))

surface_temperature = ds.createVariable('surface_temperature', np.float64, ('time'))
surface_temperature.units = 'K'

surface_humidity = ds.createVariable('surface_humidity', np.float64, ('time'))
surface_humidity.units = 'kg.kg-1'

surface_sensible_heat_flux = ds.createVariable('surface_sensible_heat_flux', np.float64, ('time'))
surface_sensible_heat_flux.units = 'W.m-2'

surface_latent_heat_flux = ds.createVariable('surface_latent_heat_flux', np.float64, ('time'))
surface_latent_heat_flux.units = 'W.m-2'


times[:] = [1.0,2.0,5.0]

################ surface_temperature ##############
surface_temperature[0] = 298.0
surface_temperature[1] = 298.1
surface_temperature[2] = 298.25

################ surface_humidity ##############
surface_humidity[0] = 0.000002478
surface_humidity[1] = 0.0000025
surface_humidity[2] = 0.0000026

################ surface_temperature ##############
surface_sensible_heat_flux[0] = 0.0
surface_sensible_heat_flux[1] = 0.25
surface_sensible_heat_flux[2] = 0.5

################ surface_temperature ##############
surface_latent_heat_flux[0] = 0.0
surface_latent_heat_flux[1] = 0.25
surface_latent_heat_flux[2] = 0.5

ds.close()
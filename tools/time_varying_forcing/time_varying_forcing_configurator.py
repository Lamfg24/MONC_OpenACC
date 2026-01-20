# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import netCDF4 as nc
import numpy as np

fn = '/home/delbeke/Documents/vn1.0.4_monc_bellouin_delbeke/tendencies_forcing.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')

time = ds.createDimension('time', 3)
altitude = ds.createDimension('z', 10)

times = ds.createVariable('time', np.float64, ('time',))
altitudes = ds.createVariable('z', np.float64, ('z',))
theta = ds.createVariable('theta', np.float64, ('time', 'z',))
theta.units = 'K'
qv = ds.createVariable('qv', np.float64, ('time', 'z',))
qv.units = 'kg.kg-1'
ql = ds.createVariable('ql', np.float64, ('time', 'z',))
ql.units = 'kg.kg-1'
qr = ds.createVariable('qr', np.float64, ('time', 'z',))
qr.units = 'kg.kg-1'
qi = ds.createVariable('qi', np.float64, ('time', 'z',))
qi.units = 'kg.kg-1'
qs = ds.createVariable('qs', np.float64, ('time', 'z',))
qs.units = 'kg.kg-1'
qg = ds.createVariable('qg', np.float64, ('time', 'z',))
qg.units = 'kg.kg-1'
qAitk = ds.createVariable('qAitk', np.float64, ('time', 'z',))
qAitk.units = 'kg.kg-1'
qAccSol = ds.createVariable('qAccSol', np.float64, ('time', 'z',))
qAccSol.units = 'kg.kg-1'

qAccInsol = ds.createVariable('qAccInsol', np.float64, ('time', 'z',))
qAccInsol.units = 'kg.kg-1'

qCoarse = ds.createVariable('qCoarse', np.float64, ('time', 'z',))
qCoarse.units = 'kg.kg-1'

qDust = ds.createVariable('qDust', np.float64, ('time', 'z',))
qDust.units = 'kg.kg-1'

w_sub = ds.createVariable('w_sub', np.float64, ('time', 'z',))
w_sub.units = 'm.s-1'
ug = ds.createVariable('ug', np.float64, ('time', 'z',))
ug.units = 'm.s-1'
vg = ds.createVariable('vg', np.float64, ('time', 'z',))
vg.units = 'm.s-1'


times[:] = [1.0,2.0,5.0]

altitudes[:] = [0.0,2500,5000,7500,9000,10000,11000,12000,14000,16000]

################ theta ##############
theta[0, : ] = np.linspace(0, 10, num=10)
theta[1, :] = np.linspace(0, 20, num=10)
theta[2, :] = np.linspace(0, 30, num=10)

################ qv ##############
qv[0, :] = np.linspace(0, 0.000002478, num=10)
qv[1, :] = np.linspace(0, 0.0000035, num=10)
qv[2, :] = np.linspace(0, 0.0000040, num=10)

################ ql ##############
ql[0, :] = np.linspace(0, 0.0000002478, num=10)
ql[1, :] = np.linspace(0, 0.00000035, num=10)
ql[2, :] = np.linspace(0, 0.00000040, num=10)

################ qr ##############
qr[0, :] = np.linspace(0, 0.00000002478, num=10)
qr[1, :] = np.linspace(0, 0.000000035, num=10)
qr[2, :] = np.linspace(0, 0.000000040, num=10)

################ qi ##############
qi[0, :] = np.linspace(0, 0.0000002478, num=10)
qi[1, :] = np.linspace(0, 0.00000035, num=10)
qi[2, :] = np.linspace(0, 0.00000040, num=10) 

################ qs ##############
qs[0, :] = np.linspace(0, 0.00000002478, num=10)
qs[1, :] = np.linspace(0, 0.000000035, num=10)
qs[2, :] = np.linspace(0, 0.000000040, num=10)

################ qg ##############
qg[0, :] = np.linspace(0, 0.000000002478, num=10)
qg[1, :] = np.linspace(0, 0.0000000035, num=10) 
qg[2, :] = np.linspace(0, 0.0000000040, num=10)

################ qAitk ##############
qAitk[0, :] = np.linspace(0, 1.0e-14, num=10)
qAitk[1, :] = np.linspace(0, 2.0e-14, num=10) 
qAitk[2, :] = np.linspace(0, 3.0e-14, num=10)

################ qAccSol ##############
qAccSol[0, :] = np.linspace(0, 1.0e-14, num=10)
qAccSol[1, :] = np.linspace(0, 2.0e-14, num=10) 
qAccSol[2, :] = np.linspace(0, 3.0e-14, num=10)

################ qAccInsol ##############
qAccInsol[0, :] = np.linspace(0, 1.0e-14, num=10)
qAccInsol[1, :] = np.linspace(0, 2.0e-14, num=10) 
qAccInsol[2, :] = np.linspace(0, 3.0e-14, num=10)

################ qCoarse ##############
qCoarse[0, :] = np.linspace(0, 1.0e-14, num=10)
qCoarse[1, :] = np.linspace(0, 2.0e-14, num=10) 
qCoarse[2, :] = np.linspace(0, 3.0e-14, num=10)

################ qDust ##############
qDust[0, :] = np.linspace(0, 1.0e-14, num=10)
qDust[1, :] = np.linspace(0, 2.0e-14, num=10) 
qDust[2, :] = np.linspace(0, 3.0e-14, num=10)

################ w_sub ##############
w_sub[0, : ] = np.linspace(0, 1, num=10)
w_sub[1, :] = np.linspace(0, 2, num=10) 
w_sub[2, :] = np.linspace(0, 3, num=10)

################ ug ##############
ug[0, : ] = np.linspace(0, 0.1, num=10)
ug[1, :] = np.linspace(0, 0.2, num=10)
ug[2, :] = np.linspace(0, 0.3, num=10)

################ vg ##############
vg[0, : ] = np.linspace(0, 0.1, num=10)
vg[1, :] = np.linspace(0, 0.2, num=10)
vg[2, :] = np.linspace(0, 0.2, num=10)

ds.close()
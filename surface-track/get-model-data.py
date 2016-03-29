# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:56:56 2015

@author: bling
"""

import sys
import netCDF4
import numpy as np

model = ['massbay','GOM3'] # 'massbay',GOM3, 30yr, roms
basic_data = 'OFF' # ON, OFF

if 'massbay' in model:
    url = """http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?
    lon[0:1:98431],lat[0:1:98431],lonc[0:1:165094],latc[0:1:165094],siglay[0:1:9][0:1:98431],h[0:1:98431],nv[0:1:2][0:1:165094],
    time[0:1:144],Times[0:1:144],zeta[0:1:144][0:1:98431],nbe[0:1:2][0:1:165094],u[0:1:144][0:1:9][0:1:165094],v[0:1:144][0:1:9][0:1:165094]"""
    data = netCDF4.Dataset(url)
    if basic_data == 'ON': 
        lon = data.variables['lon'][:]; lat = data.variables['lat'][:]
        lonc = data.variables['lonc'][:]; latc = data.variables['latc'][:]
        #siglay = data.variables['siglay'][:]; h = data.variables['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        np.savez('FVCOM_massbay_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc)#,siglay=siglay,h=h,nbe=nbe
    time = data.variables['time'][:]; Times = data.variables['Times'][:]
    u = data.variables['u'][:,0,:]; v = data.variables['v'][:,0,:];
    #zeta = data.variables['zeta'][:];
    #print 'Start saving...'
    
    np.savez('FVCOM_massbay_realtime_data.npz',time=time,Times=Times,u=u,v=v)#,zeta=zeta

if 'GOM3' in model:
    url = """http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?lon[0:1:53086],lat[0:1:53086],
    lonc[0:1:99136],latc[0:1:99136],siglay[0:1:39][0:1:53086],h[0:1:53086],time[0:1:144],Times[0:1:144],zeta[0:1:144][0:1:53086],
    nbe[0:1:2][0:1:99136],u[0:1:144][0:1:39][0:1:99136],v[0:1:144][0:1:39][0:1:99136]"""
    data = netCDF4.Dataset(url)

    if basic_data == 'ON': 
        lon = data.variables['lon'][:]; lat = data.variables['lat'][:]
        lonc = data.variables['lonc'][:]; latc = data.variables['latc'][:]
        #siglay = data.variables['siglay'][:]; h = data.variables['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        np.savez('FVCOM_GOM3_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc) #,siglay=siglay,h=h,nbe=nbe
    time = data.variables['time'][:]; Times = data.variables['Times'][:]
    u = data.variables['u'][:,0,:]; v = data.variables['v'][:,0,:];
    #zeta = data.variables['zeta'][:];
    #print 'Start saving...'
    
    np.savez('FVCOM_GOM3_realtime_data.npz',time=time,Times=Times,u=u,v=v)#,zeta=zeta

'''if '30yr' in model:
    url = """http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?lon[0:1:48450],lat[0:1:48450],
    lonc[0:1:90414],latc[0:1:90414],siglay[0:1:44][0:1:48450],h[0:1:48450],nbe[0:1:2][0:1:90414],time[0:1:324779],
    Times[0:1:324779],zeta[0:1:324779][0:1:48450],u[0:1:324779][0:1:44][0:1:90414],v[0:1:324779][0:1:44][0:1:90414]"""
    data = netCDF4.Dataset(url)
    time = data.variables['time'][:]; Times = data.variables['Times'][:]
    u = data.variables['u'][:,0,:]; v = data.variables['v'][:,0,:];
    #zeta = data.variables['zeta'][:];
    print 'Start saving...'
    
    np.savez('FVCOM_30yr_realtime_data.npz',time=time,Times=Times,u=u,v=v)#,zeta=zeta'''

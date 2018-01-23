# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:56:56 2015

@author: bling
"""

import sys
import netCDF4
import numpy as np

model = ['massbay','GOM3','30yr', 'roms','global'] # 'massbay',GOM3, 
basic_data = 'ON' # ON, OFF
#boundary_data = 'ON' # ON, OFF
realtime_data = 'OFF' # ON, OFF

if 'massbay' in model:
    url = "http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc"
    data = netCDF4.Dataset(url).variables
    
    if basic_data == 'ON': 
        lon = data['lon'][:]; lat = data['lat'][:]
        lonc = data['lonc'][:]; latc = data['latc'][:]
        siglay = data['siglay'][:]; h = data['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        
        #if boundary_data == 'ON':
        nbe1=data['nbe'][0];nbe2=data['nbe'][1];
        nbe3=data['nbe'][2]
        pointt = np.vstack((nbe1,nbe2,nbe3)).T
        wl=[]
        for i in pointt:
            if 0 in i: 
                wl.append(1)
            else:
                wl.append(0)
        tf = np.array(wl)
        inde = np.where(tf==True)
        #b_index = inde[0]
        lonb = lonc[inde]; latb = latc[inde]        
        b_points = np.vstack((lonb,latb)).T
        #np.save('Boundary_points_massbay',b_points)
        np.savez('FVCOM_massbay_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc,siglay=siglay,h=h,b_points=b_points)#,nbe=nbe
        
    if realtime_data == 'ON':
        time = data['time'][:]; Times = data['Times'][:]
        u = data['u'][:,0,:][:]; v = data['v'][:,0,:][:];
        #zeta = data.variables['zeta'][:];
        #print 'Start saving...'
        
        np.savez('FVCOM_massbay_realtime_data.npz',time=time,Times=Times,u=u,v=v)#,zeta=zeta

if 'GOM3' in model:
    url = "http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc"
    data = netCDF4.Dataset(url).variables

    if basic_data == 'ON': 
        lon = data['lon'][:]; lat = data['lat'][:]
        lonc = data['lonc'][:]; latc = data['latc'][:]
        siglay = data['siglay'][:]; h = data['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        
        #if boundary_data == 'ON':
        nbe1=data['nbe'][0];nbe2=data['nbe'][1];
        nbe3=data['nbe'][2]
        pointt = np.vstack((nbe1,nbe2,nbe3)).T
        wl=[]
        for i in pointt:
            if 0 in i: 
                wl.append(1)
            else:
                wl.append(0)
        tf = np.array(wl)
        inde = np.where(tf==True)
        #b_index = inde[0]
        lonb = lonc[inde]; latb = latc[inde]        
        b_points = np.vstack((lonb,latb)).T
        #np.save('Boundary_points_GOM3',b_points)
        np.savez('FVCOM_GOM3_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc,siglay=siglay,h=h,b_points=b_points) #,nbe=nbe
            
    if realtime_data == 'ON':
        time = data['time'][:]; Times = data['Times'][:]
        u = data['u'][:,0,:][:]; v = data['v'][:,0,:][:];
        #zeta = data.variables['zeta'][:];
        #print 'Start saving...'
        
        np.savez('FVCOM_GOM3_realtime_data.npz',time=time,Times=Times,u=u,v=v)#,zeta=zeta

if '30yr' in model:
    url = """http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3"""
    data = netCDF4.Dataset(url).variables
    if basic_data == 'ON': 
        lon = data['lon'][:]; lat = data['lat'][:]
        lonc = data['lonc'][:]; latc = data['latc'][:]
        siglay = data['siglay'][:]; h = data['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        
        #if boundary_data == 'ON':
        nbe1=data['nbe'][0];nbe2=data['nbe'][1];
        nbe3=data['nbe'][2]
        pointt = np.vstack((nbe1,nbe2,nbe3)).T
        wl=[]
        for i in pointt:
            if 0 in i: 
                wl.append(1)
            else:
                wl.append(0)
        tf = np.array(wl)
        inde = np.where(tf==True)
        #b_index = inde[0]
        lonb = lonc[inde]; latb = latc[inde]        
        b_points = np.vstack((lonb,latb)).T
        #np.save('Boundary_points_30yr',b_points)
        np.savez('FVCOM_30yr_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc,siglay=siglay,h=h,b_points=b_points) #nbe=nbe
        
    if realtime_data == 'ON':
        time = data['time'][:]; Times = data['Times'][:]
        u = data['u'][:,0,:]; v = data['v'][:,0,:];
        #zeta = data.variables['zeta'][:];
        #print 'Start saving...'
        
        np.savez('FVCOM_30yr_realtime_data.npz',time=time,Times=Times,u=u,v=v)#,zeta=zeta'''
        
if 'roms' in model:
    url = """http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2013_da/his/ESPRESSO_Real-Time_v2_History_fmrc.ncd"""
    data = netCDF4.Dataset(url).variables
    if basic_data == 'ON': 
        lon_rho = data['lon_rho'][:]; lat_rho = data['lat_rho'][:] 
        lon_u,lat_u = data['lon_u'][:], data['lat_u'][:]
        lon_v,lat_v = data['lon_v'][:], data['lat_v'][:]
        h = data['h'][:]; s_rho = data['s_rho'][:]
        mask_u = data['mask_u'][:]; mask_v = data['mask_v'][:]; mask_rho = data['mask_rho'][:]
        
        np.savez('roms_basic_data.npz',lon_rho=lon_rho,lat_rho=lat_rho,lon_u=lon_u,lat_u=lat_u,lon_v=lon_v,lat_v=lat_v,s_rho=s_rho,h=h,mask_u=mask_u,mask_v=mask_v,mask_rho=mask_rho) 
        
    if realtime_data == 'ON':
        time = data['time'][:]; run = data['run'][:]
        u = data['u']; v = data['v']; zeta = data['zeta']
        #zeta = data.variables['zeta'][:];
        #print 'Start saving...'
        
        np.savez('roms_realtime_data.npz',time=time,run=run,u=u,v=v,zeta=zeta)#'''
        
if 'global' in model:
    url = "http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_GLOBAL_FORECAST.nc"
    data = netCDF4.Dataset(url).variables

    if basic_data == 'ON': 
        lon = data['lon'][:]; lat = data['lat'][:]
        lonc = data['lonc'][:]; latc = data['latc'][:]
        siglay = data['siglay'][:]; h = data['h'][:]
        
        for i in range(len(lonc)):
                if lonc[i] > 180:
                    lonc[i] = lonc[i]-360
        for i in range(len(lon)):
                if lon[i] > 180:
                    lon[i] = lon[i]-360
        #nbe = data.variables['nbe'][:]#'''
        
        np.savez('FVCOM_global_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc,siglay=siglay,h=h) #,nbe=nbe

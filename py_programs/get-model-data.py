# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:56:56 2015

@author: bling
"""

import sys
import netCDF4
import numpy as np

model = ['massbay','GOM3'] # 'massbay',GOM3, 30yr, roms
basic_data = 'ON' # ON, OFF
boundary_data = 'ON' # ON, OFF

if 'massbay' in model:
    url = "http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?"
    data = netCDF4.Dataset(url).variables
    
    if basic_data == 'ON': 
        lon = data['lon'][:]; lat = data['lat'][:]
        lonc = data['lonc'][:]; latc = data['latc'][:]
        #siglay = data.variables['siglay'][:]; h = data.variables['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        np.savez('FVCOM_massbay_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc)#,siglay=siglay,h=h,nbe=nbe
        if boundary_data == 'ON':
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
            np.save('Boundary_points_massbay',b_points)
            
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
        #siglay = data.variables['siglay'][:]; h = data.variables['h'][:]
        #nbe = data.variables['nbe'][:]#'''
        np.savez('FVCOM_GOM3_basic_data.npz',lon=lon,lat=lat,lonc=lonc,latc=latc) #,siglay=siglay,h=h,nbe=nbe
        if boundary_data == 'ON':
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
            np.save('Boundary_points_GOM3',b_points)
            
    time = data['time'][:]; Times = data['Times'][:]
    u = data['u'][:,0,:][:]; v = data['v'][:,0,:][:];
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

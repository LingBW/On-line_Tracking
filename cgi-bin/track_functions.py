#!/usr/bin/env /anaconda/bin/python
import cgitb
cgitb.enable()
#import sys
import netCDF4
from datetime import datetime,timedelta
import numpy as np
import pandas as pd
from dateutil.parser import parse
#import pytz
from mpl_toolkits.basemap import Basemap
from matplotlib.path import Path
import math
import sys

def dm2dd(lat,lon):
    """
    convert lat, lon from decimal degrees,minutes to decimal degrees
    """
    (a,b)=divmod(float(lat),100.)   
    aa=int(a)
    bb=float(b)
    lat_value=aa+bb/60.

    if float(lon)<0:
        (c,d)=divmod(abs(float(lon)),100.)
        cc=int(c)
        dd=float(d)
        lon_value=cc+(dd/60.)
        lon_value=-lon_value
    else:
        (c,d)=divmod(float(lon),100.)
        cc=int(c)
        dd=float(d)
        lon_value=cc+(dd/60.)
    return lat_value, -lon_value
def getrawdrift(did,filename):
   '''
   routine to get raw drifter data from ascii files posted on the web
   '''
   url='http://nefsc.noaa.gov/drifter/'+filename
   df=pd.read_csv(url,header=None, delimiter=r"\s+")
   # make a datetime
   dtime=[]
   index = np.where(df[0]==int(did))[0]
   newData = df.ix[index]
   now = datetime.now()
   for k in newData[0].index:
      dt1=datetime(now.year, newData[2][k],newData[3][k],newData[4][k],newData[5][k])
      dtime.append(dt1)
   return newData[8],newData[7],dtime,newData[9]

def getdrift(did):
    """
    routine to get drifter data from archive based on drifter id (did)
    -assumes "import pandas as pd" has been issued above
    -get remotely-stored drifter data via ERDDAP
    -input: deployment id ("did") number where "did" is a string
    -output: time(datetime), lat (decimal degrees), lon (decimal degrees), depth (meters)
    
    note: there is another function below called "data_extracted" that does a similar thing returning a dictionary
    
    Jim Manning June 2014
    """
    url = 'http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.csv?time,latitude,longitude,depth&id="'+did+'"&orderBy("time")'
    df=pd.read_csv(url,skiprows=[1]) #returns a dataframe with all that requested
    # generate this datetime 
    for k in range(len(df)):
       df.time[k]=parse(df.time[k]) # note this "parse" routine magically converts ERDDAP time to Python datetime
    return df.latitude.values,df.longitude.values,df.time.values,df.depth.values  

def get_nc_data(url, *args):
    '''
    get specific dataset from url

    *args: dataset name, composed by strings
    ----------------------------------------
    example:
        url = 'http://www.nefsc.noaa.gov/drifter/drift_tcs_2013_1.dat'
        data = get_url_data(url, 'u', 'v')
    '''
    nc = netCDF4.Dataset(url)
    data = {}
    for arg in args:
        try:
            data[arg] = nc.variables[arg]
        except (IndexError, NameError, KeyError):
            print 'Dataset {0} is not found'.format(arg)
    return data
    
class get_roms():
    '''
    ####(2009.10.11, 2013.05.19):version1(old) 2009-2013
    ####(2013.05.19, present): version2(new) 2013-present
    (2006.01.01 01:00, 2014.1.1 00:00)
    '''
    
    def __init__(self):
        pass
    
    def nearest_point(self, lon, lat, lons, lats, length=0.06):  #0.3/5==0.06
        '''Find the nearest point to (lon,lat) from (lons,lats),
           return the nearest-point (lon,lat)
           author: Bingwei'''
        p = Path.circle((lon,lat),radius=length)
        #numpy.vstack(tup):Stack arrays in sequence vertically
        points = np.vstack((lons.flatten(),lats.flatten())).T  
        
        insidep = []
        #collect the points included in Path.
        for i in xrange(len(points)):
            if p.contains_point(points[i]):# .contains_point return 0 or 1
                insidep.append(points[i])  
        # if insidep is null, there is no point in the path.
        if not insidep:
            #print 'This point out of the model area or hits the land.'
            raise Exception()
        #calculate the distance of every points in insidep to (lon,lat)
        distancelist = []
        for i in insidep:
            ss=math.sqrt((lon-i[0])**2+(lat-i[1])**2)
            distancelist.append(ss)
        # find index of the min-distance
        mindex = np.argmin(distancelist)
        # location the point
        lonp = insidep[mindex][0]; latp = insidep[mindex][1]
        
        return lonp,latp
        
    def get_data(self, starttime, endtime):
        '''
        get url according to starttime and endtime.
        '''
        mstime = starttime.replace(hour=0,minute=0,second=0,microsecond=0)
        self.hours = int((endtime-starttime).total_seconds()/60/60) # get total hours
        # time_r = datetime(year=2006,month=1,day=9,hour=1,minute=0)
        urltime = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/espresso/2013_da/his/ESPRESSO_Real-Time_v2_History_fmrc.ncd'
        self.run = netCDF4.Dataset(urltime).variables['run'][:]
        
        cm = (mstime-datetime(2013,05,18)).total_seconds()/60/60
        self.ind1 = np.argmin(abs(self.run-cm))
        modtime = netCDF4.Dataset(urltime).variables['time'][self.ind1]
        ############ Display modetime.#############
        fmodtime = datetime(2013,05,18) + timedelta(hours=modtime[0])
        emodtime = datetime(2013,05,18) + timedelta(hours=modtime[-1]) #fmodtime + timedelta(days=6)
        # judge if the starttime and endtime in the model time horizon
        if starttime<fmodtime or starttime>emodtime or endtime<fmodtime or endtime>emodtime:
            url = 'error'
            return url,fmodtime,emodtime
        
        self.ind2 = int(round((mstime-fmodtime).total_seconds()/60/60))+1
        self.ind3 = self.ind2+self.hours-1
        ldata = netCDF4.Dataset(urltime).variables
       
        self.u = ldata['u']; self.v = ldata['v']#; self.zeta = data['zeta']
        #print '2' 
        # Basic model data
        data = np.load('/var/www/cgi-bin/ioos/track/roms_basic_data.npz')
        self.lon_rho = data['lon_rho']; self.lat_rho = data['lat_rho'] 
        self.lon_u,self.lat_u = data['lon_u'], data['lat_u']
        self.lon_v,self.lat_v = data['lon_v'], data['lat_v']
        self.h = data['h']; self.s_rho = data['s_rho']
        self.mask_u = data['mask_u']; self.mask_v = data['mask_v']; self.mask_rho = data['mask_rho']
        
        return urltime,fmodtime,emodtime
        
    def get_track(self,lon,lat,depth,track_way):#, depth
        '''
        get the nodes of specific time period
        lon, lat: start point
        depth: 0~35, the 0th is the bottom.
        '''
        
        nodes = dict(lon=[lon], lat=[lat])
        lonrp,latrp = self.nearest_point(lon,lat,self.lon_rho,self.lat_rho)
        indexr = np.where(self.lon_rho==lonrp)
        if not self.mask_rho[indexr]:
            return nodes
        
        waterdepth = self.h[indexr]
        if waterdepth<(abs(depth)): 
            print 'This point is too shallow.Less than %d meter.'%abs(depth)
            raise Exception()
        depth_total = self.s_rho*waterdepth  
        layer = np.argmin(abs(depth_total+depth))
        
        t = abs(self.hours)
        for i in xrange(t):  #Roms points update every 2 hour
            
            try:
                
                lonup,latup = self.nearest_point(lon,lat,self.lon_u,self.lat_u)
                lonvp,latvp = self.nearest_point(lon,lat,self.lon_v,self.lat_v)
                indexu = np.where(self.lon_u==lonup)
                indexv = np.where(self.lon_v==lonvp)
            
                if not self.mask_u[indexu]:
                    #print 'No u velocity.'
                    raise Exception()
                if not self.mask_v[indexv]:
                    #print 'No v velocity'
                    raise Exception() 
            except:
                return nodes
            if track_way=='backward': # backwards case
                #u_t = -1*self.u[t-i,layer][indexu][0]; v_t = -1*self.v[t-i,layer][indexv][0]
                u_t = -1*self.u[self.ind1,self.ind3-i,layer][indexu][0] 
                v_t = -1*self.v[self.ind1,self.ind3-i,layer][indexv][0]
            else:
                #u_t = self.u[i,layer][indexu][0]; v_t = self.v[i,layer][indexv][0]
                u_t = self.u[self.ind1,self.ind2+i,layer][indexu][0]
                v_t = self.v[self.ind1,self.ind2+i,layer][indexv][0]
            #print 'u_t,v_t',u_t,v_t
            if np.isnan(u_t) or np.isnan(v_t): #There is no water
                print 'Sorry, the point on the land or hits the land. Info: u or v is NAN'
                return nodes
            dx = 60*60*u_t#float(u_p)
            dy = 60*60*v_t#float(v_p)
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')                        
            #x,y = mapx(lon,lat)
            #lon,lat = mapx(x+dx,y+dy,inverse=True)            
            lon = lon + dx/(111111*np.cos(lat*np.pi/180))
            lat = lat + dy/111111
            #print '%d,lat,lon,layer'%(i+1),lat,lon,layer
            nodes['lon'].append(lon);nodes['lat'].append(lat)
            
        return nodes
        
    def shrink_data(self,lon,lat,lons,lats):
        # ind = argwhere((lonc >= size[0]) & (lonc <= size[1]) & (latc >= size[2]) & (latc <= size[3]))
        lont = []; latt = []
        p = Path.circle((lon,lat),radius=0.6)
        pints = np.vstack((lons.flatten(),lats.flatten())).T
        for i in range(len(pints)):
            if p.contains_point(pints[i]):
                lont.append(pints[i][0])
                latt.append(pints[i][1])
        lonl=np.array(lont); latl=np.array(latt)#'''
        if not lont:
            print 'point position error! shrink_data'
            sys.exit()
        return lonl,latl

class get_fvcom():
    def __init__(self, mod):
        self.modelname = mod
    
    def get_url(self, starttime, trackdays):
        '''
        get different url according to starttime and endtime.
        urls are monthly.
        '''
        #self.hours = int(round((endtime-starttime).total_seconds()/60/60))
        endtime = starttime + timedelta(trackdays)
        self.starttime = starttime; self.endtime = endtime
        self.trackdays = trackdays
        #int(round((endtime-starttime).total_seconds()/60/60/24))
        #print self.hours
                
        if self.modelname == "global":
            turl = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_GLOBAL_FORECAST.nc'''
            
            try:
                self.tdata = netCDF4.Dataset(turl).variables
                #MTime = tdata['time'] 
                MTime = self.tdata['Times']
            except:
                print '"massbay" database is unavailable!'
                raise Exception()
            #Times = netCDF4.num2date(MTime[:],MTime.units)
            Times = []
            for i in MTime:
                strt = '201'+i[3]+'-'+i[5]+i[6]+'-'+i[8]+i[9]+' '+i[11]+i[12]+':'+i[14]+i[15]
                Times.append(datetime.strptime(strt,'%Y-%m-%d %H:%M'))#'''
            fmodtime = Times[0]; emodtime = Times[-1]+timedelta(days=1)    
            if starttime<fmodtime or starttime>emodtime or endtime<fmodtime or endtime>emodtime:
                #print 'Time: Error! Model(global) only works between %s with %s(UTC).'%(fmodtime,emodtime)
                #raise Exception()
                url = 'error' 
                return url,fmodtime,emodtime
            '''
            npTimes = np.array(Times)
            tm1 = npTimes-starttime; #tm2 = mtime-t2
            index1 = np.argmin(abs(tm1))#
            #index1 = netCDF4.date2index(starttime,MTime,select='nearest')
            index2 = index1 + self.days#
            #print 'index1,index2',index1,index2
            #url = url.format(index1, index2)
            self.mTime = Times[index1:index2] #'''
            self.mTime = np.array(Times) 
            #self.url = turl
            #loncs = self.tdata['lonc'][:]; self.latc = self.tdata['latc'][:]  #quantity:165095
            # Basic model data.
            self.basicdata = np.load('/var/www/cgi-bin/ioos/track/FVCOM_global_basic_data.npz')            
            self.lonc, self.latc = self.basicdata['lonc'], self.basicdata['latc']  #quantity:165095
            self.lons, self.lats = self.basicdata['lon'], self.basicdata['lat']
            self.h = self.basicdata['h']; self.siglay = self.basicdata['siglay']#; print '3'
            
            # Real-time model data.
            self.iy = [i for i in range(len(self.mTime)) if self.mTime[i].day==self.starttime.day]
            
            self.mtime = np.array(Times[self.iy[0]:self.iy[0]+int(trackdays)+1])
            
        
        return turl,fmodtime,emodtime
        
    def get_globaltrack(self,lon,lat,depth,track_way): #,b_index,nvdepth,,bcon 
        '''
        Get forecast points start at lon,lat
        '''
        
        modpts = dict(lon=[lon], lat=[lat])
        # For boundary.
        self.gmap=Basemap(projection='cyl',llcrnrlat=lat-self.trackdays, urcrnrlat=lat+self.trackdays, llcrnrlon=lon-self.trackdays, urcrnrlon=lon+self.trackdays,resolution='l')
        
        nodeindex = np.argwhere((self.lons >= lon-0.3) & (self.lons <= lon+0.3) & (self.lats >= lat-0.3) & (self.lats <= lat+0.3))        
        if len(nodeindex)==0:
            print 'No model data around. Out of Model area.'
            return modpts 
        waterdepth = self.h[[nodeindex]]; #print waterdepth
        if np.mean(waterdepth)<(abs(depth)): 
            print 'This point is too shallow.Less than %d meter.'%abs(depth)
            return modpts #raise Exception()
        depth_total = self.siglay[:,nodeindex[0]]*np.mean(waterdepth); #print depth_total
        layer = np.argmin(abs(depth_total+depth)); #print 'layer',layer
        
        u = self.tdata['u'][self.iy[0]:self.iy[0]+int(self.trackdays)+1,layer,:]; #print '5'
        v = self.tdata['v'][self.iy[0]:self.iy[0]+int(self.trackdays)+1,layer,:]; 
        #modpts = [(lon,lat)]#;st = []
        
        ld = self.trackdays*12
        for j in xrange(int(ld)): 
            if self.gmap.is_land(lat,lon):
                return modpts
            inds = np.argwhere((self.lonc >= lon-0.3) & (self.lonc <= lon+0.3) & (self.latc >= lat-0.3) & (self.latc <= lat+0.3))
            if len(inds) == 0:
                return modpts#print 'inds',inds

            if track_way=='backward' : # backwards case
                tratime = self.endtime-timedelta(hours=j*2)
                iday = [i for i in range(len(self.mtime)) if self.mtime[i].day==tratime.day]; #print iday
                u_t1 = np.mean(u[iday[0]][inds])*(-1); v_t1 = np.mean(v[iday[0]][inds])*(-1)
            else:
                tratime = self.starttime+timedelta(hours=j*2)
                iday = [i for i in range(len(self.mtime)) if self.mtime[i].day==tratime.day]; #print iday
                u_t1 = np.mean(u[iday[0]][inds]); v_t1 = np.mean(v[iday[0]][inds])
            
            dx = 2*60*60*u_t1; dy = 2*60*60*v_t1
            #pspeed = math.sqrt(u_t1**2+v_t1**2)
            #modpts['spd'].append(pspeed)
                     
            lon = lon + (dx/(111111*np.cos(lat*np.pi/180)))
            lat = lat + dy/111111 #'''   
            modpts['lon'].append(lon); modpts['lat'].append(lat)
            
        return modpts
        
    def get_data(self, starttime, endtime):
        '''
        get different url according to starttime and trackdays.
        urls are monthly.
        '''
        #endtime = starttime + timedelta(days=trackdays)
        self.hours = int(round((endtime-starttime).total_seconds()/3600))
        #self.hours = trackdays*24
        #print self.hours
                
        if self.modelname == "GOM3":
            
            turl = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'''  
            
            todata = netCDF4.Dataset(turl).variables
            mTime = todata['Times'][:]
            Times = []
            for i in mTime:
                strt = '201'+i[3]+'-'+i[5]+i[6]+'-'+i[8]+i[9]+' '+i[11]+i[12]+':'+i[14]+i[15]
                Times.append(datetime.strptime(strt,'%Y-%m-%d %H:%M'))
            fmodtime = Times[0]; emodtime = Times[-1]
            if starttime<fmodtime or starttime>emodtime or endtime<fmodtime or endtime>emodtime:
                url = 'error'
                return url,fmodtime,emodtime
            npTimes = np.array(Times)
            tm1 = npTimes-starttime.replace(minute=0,second=0,microsecond=0); #tm2 = mtime-t2
            self.ind1 = np.argmin(abs(tm1))
            self.ind2 = self.ind1 + self.hours#'''
            self.datatime = npTimes[self.ind1:self.ind2]
            
            # Basic model data.
            self.basicdata = np.load('/var/www/cgi-bin/ioos/track/FVCOM_GOM3_basic_data.npz')            
            self.lonc, self.latc = self.basicdata['lonc'], self.basicdata['latc']  #quantity:165095
            self.lons, self.lats = self.basicdata['lon'], self.basicdata['lat']
            self.h = self.basicdata['h']; self.siglay = self.basicdata['siglay']#; print '3'
            self.b_points = self.basicdata['b_points']
            # real-time model data
            self.u = todata['u']; self.v = todata['v']; 
            
        elif self.modelname == "massbay":
            
            turl = "http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc"
            
            todata = netCDF4.Dataset(turl).variables
            mTime = todata['Times'][:]
            Times = []
            for i in mTime:
                strt = '201'+i[3]+'-'+i[5]+i[6]+'-'+i[8]+i[9]+' '+i[11]+i[12]+':'+i[14]+i[15]
                Times.append(datetime.strptime(strt,'%Y-%m-%d %H:%M'))
            fmodtime = Times[0]; emodtime = Times[-1]
            if starttime<fmodtime or starttime>emodtime or endtime<fmodtime or endtime>emodtime:
                url = 'error'
                return url,fmodtime,emodtime
            npTimes = np.array(Times)
            tm1 = npTimes-starttime.replace(minute=0,second=0,microsecond=0); #tm2 = mtime-t2
            self.ind1 = np.argmin(abs(tm1))
            self.ind2 = self.ind1 + self.hours#'''
            self.datatime = npTimes[self.ind1:self.ind2]
            
            # Basic model data.
            self.basicdata = np.load('/var/www/cgi-bin/ioos/track/FVCOM_massbay_basic_data.npz')            
            self.lonc, self.latc = self.basicdata['lonc'], self.basicdata['latc']  #quantity:165095
            self.lons, self.lats = self.basicdata['lon'], self.basicdata['lat']
            self.h = self.basicdata['h']; self.siglay = self.basicdata['siglay']
            self.b_points = self.basicdata['b_points']
            
            # real-time model data
            self.u = todata['u']; self.v = todata['v']; 
            
        elif self.modelname == "30yr": #start at 1977/12/31 23:00, end at 2014/1/1 0:0, time units:hours
            turl = """http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3"""
            #index1 = int(round((starttime-datetime(1977,12,31,22,58,4,0,pytz.UTC)).total_seconds()/3600))
            todata = netCDF4.Dataset(turl).variables
            mtime = todata['time'][:]
            fmodtime = datetime(1858,11,17) + timedelta(float(mtime[0]))
            emodtime = datetime(1858,11,17) + timedelta(float(mtime[-1]))
            # get number of days from 11/17/1858
            t1 = (starttime - datetime(1858,11,17)).total_seconds()/86400 
            t2 = (endtime - datetime(1858,11,17)).total_seconds()/86400
            if not mtime[0]<t1<mtime[-1] or not mtime[0]<t2<mtime[-1]:
                #raise Exception('massbay works from 1977/12/31 23:00 to 2014/1/1 0:0.')
                url = 'error'
                return url,fmodtime,emodtime
            tm1 = mtime-t1; #tm2 = mtime-t2
            self.ind1 = np.argmin(abs(tm1)); #index2 = np.argmin(abs(tm2)); print index1,index2
            self.ind2 = self.ind1 + self.hours
            Times = []
            for i in range(self.hours):
                Times.append(starttime+timedelta(i))
            self.datatime = Times

            # real-time model data
            self.u = todata['u']; self.v = todata['v']; 
            
            # Basic model data.
            self.basicdata = np.load('/var/www/cgi-bin/ioos/track/FVCOM_30yr_basic_data.npz')            
            self.lonc, self.latc = self.basicdata['lonc'], self.basicdata['latc']  #quantity:165095
            self.lons, self.lats = self.basicdata['lon'], self.basicdata['lat']
            self.h = self.basicdata['h']; self.siglay = self.basicdata['siglay']
            self.b_points = self.basicdata['b_points']
        
        
        return self.b_points,fmodtime,emodtime #,nv lons,lats,lonc,latc,,h,siglay
        
    
    def get_track(self,lon,lat,depth,track_way): #,b_index,nvdepth, 
        '''
        Get forecast points start at lon,lat
        '''
        
        modpts = dict(lon=[lon], lat=[lat]) #model forecast points, layer=[]
        
        if lon>90:
            lon, lat = dm2dd(lon, lat)
        # produce layer   
        nodeindex = np.argwhere((self.lons >= lon-0.1) & (self.lons <= lon+0.1) & (self.lats >= lat-0.1) & (self.lats <= lat+0.1))        
        
        if len(nodeindex) == 0:
                print 'No model data around. Out of Model area.'
                return modpts        
        waterdepth = self.h[[nodeindex]]; #print waterdepth
        if np.mean(waterdepth)<(abs(depth)): 
            print 'This point is too shallow.Less than %d meter.'%abs(depth)
            return modpts #raise Exception()
        depth_total = self.siglay[:,nodeindex[0]]*np.mean(waterdepth); #print depth_total
        layer = np.argmin(abs(depth_total+depth)); #print 'layer',layer
        
        u = self.u[self.ind1:self.ind2,layer,:] ; v = self.v[self.ind1:self.ind2,layer,:]        
        
        t = len(self.datatime) #print t        
        for i in xrange(t):            

            elementindex = np.argwhere((self.lonc >= lon-0.1) & (self.lonc <= lon+0.1) & (self.latc >= lat-0.1) & (self.latc <= lat+0.1))
            
            if len(elementindex) == 0:
                print 'No model data around. Out of Model area.'
                return modpts
            ################## boundary 1 ####################
            pa = self.eline_path(lon,lat)

            if track_way=='backward' : # backwards case
                #u_t1 = np.mean(self.u[self.ind2-i,layer,elementindex])*(-1); v_t1 = np.mean(self.v[self.ind2-i,layer,elementindex])*(-1)
                u_t1 = np.mean(u[t-i-1,elementindex])*(-1); v_t1 = np.mean(v[t-i-1,elementindex])*(-1)
            else:
                #u_t1 = np.mean(self.u[self.ind1+i,layer][elementindex]); v_t1 = np.mean(self.v[self.ind1+i,layer][elementindex])
                u_t1 = np.mean(u[i,elementindex]); v_t1 = np.mean(v[i,elementindex])
            dx = 60*60*u_t1; dy = 60*60*v_t1
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')                        
            #x,y = mapx(lon,lat)
            #temlon,temlat = mapx(x+dx,y+dy,inverse=True)            
            temlon = lon + (dx/(111111*np.cos(lat*np.pi/180)))
            temlat = lat + dy/111111 #'''
            
            #print '%d,lat,lon,layer'%(i+1),temlat,temlon,layer
            #########case for boundary 1 #############
            if pa:
                teml = [(lon,lat),(temlon,temlat)]
                tempa = Path(teml)
                if pa.intersects_path(tempa): 
                    print 'One point hits land here. Path Condition'
                    return modpts #'''
            
            lon = temlon; lat = temlat
            modpts['lon'].append(lon); modpts['lat'].append(lat)
            
        return modpts
    
    def eline_path(self,lon,lat):
        '''
        When drifter close to boundary(less than 0.1),find one nearest point to drifter from boundary points, 
        then find two nearest boundary points to previous boundary point, create a boundary path using that 
        three boundary points.
        '''
        def boundary_location(locindex,pointt,wl):
            '''
            Return the index of boundary points nearest to 'locindex'.
            '''
            loca = []
            dx = pointt[locindex]; #print 'func',dx 
            for i in dx: # i is a number.
                #print i  
                if i ==0 :
                    continue
                dx1 = pointt[i-1]; #print dx1
                if 0 in dx1:
                    loca.append(i-1)
                else:
                    for j in dx1:
                        if j != locindex+1:
                            if wl[j-1] == 1:
                                loca.append(j-1)
            return loca
        
        p = Path.circle((lon,lat),radius=0.02) #0.06
        dis = []; bps = []; pa = []
        tlons = []; tlats = []; loca = []
        for i in self.b_points:
            if p.contains_point(i):
                bps.append((i[0],i[1]))
                d = math.sqrt((lon-i[0])**2+(lat-i[1])**2)
                dis.append(d)
        bps = np.array(bps)
        if not dis:
            return None
        else:
            #print "Close to boundary."
            dnp = np.array(dis)
            dmin = np.argmin(dnp)
            lonp = bps[dmin][0]; latp = bps[dmin][1]
            index1 = np.where(self.lonc==lonp)
            index2 = np.where(self.latc==latp)
            elementindex = np.intersect1d(index1,index2)[0] # location 753'''
            #print index1,index2,elementindex  
            loc1 = boundary_location(elementindex,self.pointt,self.wl) ; #print 'loc1',loc1
            loca.extend(loc1)
            loca.insert(1,elementindex)               
            for i in range(len(loc1)):
                loc2 = boundary_location(loc1[i],self.pointt,self.wl); #print 'loc2',loc2
                if len(loc2)==1:
                    continue
                for j in loc2:
                    if j != elementindex:
                        if i ==0:
                            loca.insert(0,j)
                        else:
                            loca.append(j)
            
            for i in loca:
                tlons.append(self.lonc[i]); tlats.append(self.latc[i])
                       
            for i in xrange(len(tlons)):
                pa.append((tlons[i],tlats[i]))
            path = Path(pa)#,codes
            return path

class get_drifter():

    def __init__(self, drifter_id, filename=None):
        self.drifter_id = drifter_id
        self.filename = filename
    def get_track(self, starttime=None, days=None):
        '''
        return drifter nodes
        if starttime is given, return nodes started from starttime
        if both starttime and days are given, return nodes of the specific time period
        '''
        if self.filename:
            temp=getrawdrift(self.drifter_id,self.filename)
        else:
            temp=getdrift(self.drifter_id)
        nodes = {}
        nodes['lon'] = np.array(temp[1])
        nodes['lat'] = np.array(temp[0])
        nodes['time'] = np.array(temp[2])
        #starttime = np.array(temp[2][0])
        if not starttime:
            starttime = np.array(temp[2][0])
        if days:
            endtime = starttime + timedelta(days=days)
            i = self.__cmptime(starttime, nodes['time'])
            j = self.__cmptime(endtime, nodes['time'])
            nodes['lon'] = nodes['lon'][i:j+1]
            nodes['lat'] = nodes['lat'][i:j+1]
            nodes['time'] = nodes['time'][i:j+1]
        else:
            #i = self.__cmptime(starttime, nodes['time'])
            nodes['lon'] = nodes['lon']#[i:-1]
            nodes['lat'] = nodes['lat']#[i:-1]
            nodes['time'] = nodes['time']#[i:-1]
        return nodes
        
    def __cmptime(self, time, times):
        '''
        return indies of specific or nearest time in times.
        '''
        tdelta = []
        
        for t in times:
            tdelta.append(abs((time-t).total_seconds()))
            
        index = tdelta.index(min(tdelta))
        
        return index

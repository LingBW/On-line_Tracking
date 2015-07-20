#!/usr/bin/env python
import cgitb
cgitb.enable()
#import sys
import netCDF4
from datetime import datetime,timedelta
import numpy as np
import pandas as pd
from dateutil.parser import parse
import pytz
#from mpl_toolkits.basemap import Basemap
from matplotlib.path import Path
import math

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
   for k in newData[0].index:
      #dt1=dt.datetime(int(filename[-10:-6]),df[2][k],df[3][k],df[4][k],df[5][k],0,0,pytz.utc)
      dt1=datetime(2015, newData[2][k],newData[3][k],newData[4][k],newData[5][k],0,0,pytz.utc)
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

class track(object):
    def __init__(self, startpoint):
        '''
        gets the start point of the water, and the location of datafile.
        '''
        self.startpoint = startpoint
        
    def get_data(self, url):
        '''
        calls get_data
        '''        
        pass                                 
        
    def bbox2ij(self, lon, lat, lons, lats, length=0.06):  #0.3/5==0.06
        """
        Return tuple of indices of points that are completely covered by the 
        specific boundary box.
        i = bbox2ij(lon,lat,bbox)
        lons,lats = 2D arrays (list) that are the target of the subset, type: np.ndarray
        bbox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]
    
        Example
        -------  
        >>> i0,i1,j0,j1 = bbox2ij(lat_rho,lon_rho,[-71, -63., 39., 46])
        >>> h_subset = nc.variables['h'][j0:j1,i0:i1]
        length: the boundary box.
        """
        '''bbox = [lon-length, lon+length, lat-length, lat+length]
        bbox = np.array(bbox)
        mypath = np.array([bbox[[0,1,1,0]],bbox[[2,2,3,3]]]).T#'''
        p = Path.circle((lon,lat),radius=length)
        points = np.vstack((lons.flatten(),lats.flatten())).T  #numpy.vstack(tup):Stack arrays in sequence vertically
        tshape = np.shape(lons)
        inside = []
        
        for i in range(len(points)):
            inside.append(p.contains_point(points[i]))  # .contains_point return 0 or 1
            
        inside = np.array(inside, dtype=bool).reshape(tshape)
        index = np.where(inside==True)
        
        '''check if there are no points inside the given area'''        
        
        if not index[0].tolist():          # bbox covers no area
            print 'This point is out of the model area.'
            raise Exception()
        else:
            return index
            
    def nearest_point_index(self, lon, lat, lons, lats):  #,num=4
        '''
        Return the index of the nearest rho point.
        lon, lat: the coordinate of start point, float
        lats, lons: the coordinate of points to be calculated.
        '''
        def min_distance(lon,lat,lons,lats):
            '''Find out the nearest distance to (lon,lat),and return lon.distance units: meters'''
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')
            dis_set = []
            #x,y = mapx(lon,lat)
            for i,j in zip(lons,lats):
                #x2,y2 = mapx(i,j)
                ss=math.sqrt((lon-i)**2+(lat-j)**2)
                #ss=math.sqrt((x-x2)**2+(y-y2)**2)
                dis_set.append(ss)
            dis = min(dis_set)
            p = dis_set.index(dis)
            lonp = lons[p]; latp = lats[p]
            return lonp,latp,dis       
        index = self.bbox2ij(lon, lat, lons, lats)        
        lon_covered = lons[index];  lat_covered = lats[index]       
        lonp,latp,distance = min_distance(lon,lat,lon_covered,lat_covered)
        index1 = np.where(lons==lonp)
        index2 = np.where(lats==latp)
        index = np.intersect1d(index1,index2)
        #points = np.vstack((lons.flatten(),lats.flatten())).T         
        #index = [i for i in xrange(len(points)) if ([lonp,latp]==points[i]).all()]
        print 'index',index
        return index,lonp,latp,distance
        
    def get_track(self, timeperiod, data):
        pass
    
class get_roms(track):
    '''
    ####(2009.10.11, 2013.05.19):version1(old) 2009-2013
    ####(2013.05.19, present): version2(new) 2013-present
    (2006.01.01 01:00, 2014.1.1 00:00)
    '''
    
    def __init__(self):
        pass
        
    def get_url(self, starttime, endtime):
        '''
        get url according to starttime and endtime.
        '''
        '''
        self.starttime = starttime
        url_oceantime = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2006_da/his?ocean_time[0:1:69911]'
        data_oceantime = netCDF4.Dataset(url_oceantime)
        t1 = (starttime - datetime(2006,01,01,0,0,0,0,pytz.utc)).total_seconds()
        t2 = (endtime - datetime(2006,01,01,0,0,0,0,pytz.utc)).total_seconds()
        index1 = self.__closest_num(t1,data_oceantime.variables['ocean_time'][:])
        index2 = self.__closest_num(t2,data_oceantime.variables['ocean_time'][:])
        url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2006_da/his?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129]'
        url = url.format(index1, index2)
        return url
        '''
        self.starttime = starttime
        self.hours = int((endtime-starttime).total_seconds()/60/60) # get total hours
        # time_r = datetime(year=2006,month=1,day=9,hour=1,minute=0)
        # url_oceantime = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2006_da/his?ocean_time[0:1:69911]'
        url_oceantime = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd?time'
        self.oceantime = netCDF4.Dataset(url_oceantime).variables['time'][:]
        t1 = (starttime - datetime(2013,05,18, tzinfo=pytz.UTC)).total_seconds()/3600 # for url2006 it's 2006,01,01
        t2 = (endtime - datetime(2013,05,18, tzinfo=pytz.UTC)).total_seconds()/3600
        self.index1 = self.__closest_num(t1, self.oceantime)
        self.index2 = self.__closest_num(t2, self.oceantime)
        #print self.index1, self.index2
        # index1 = (starttime - time_r).total_seconds()/60/60
        # index2 = index1 + self.hours
        # url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2006_da/his?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129]'
        # url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2006_da/his?s_rho[0:1:35],h[0:1:81][0:1:129],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],ocean_time'
        # This one is hourly # url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],time,mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129]'     
        url = 'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd?h[0:1:81][0:1:129],mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129]'
             #'http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2013_da/his_Best/ESPRESSO_Real-Time_v2_History_Best_Available_best.ncd?h[0:1:81][0:1:129],s_rho[0:1:35],lon_rho[0:1:81][0:1:129],lat_rho[0:1:81][0:1:129],temp[{0}:1:{1}][0:1:35][0:1:81][0:1:129],time,mask_rho[0:1:81][0:1:129],u[{0}:1:{1}][0:1:35][0:1:81][0:1:128],v[{0}:1:{1}][0:1:35][0:1:80][0:1:129]'      
        url = url.format(self.index1, self.index2)
        return url
    def __closest_num(self, num, numlist, i=0):
        '''
        Return index of the closest number in the list
        '''
        index1,index2 = 0, len(numlist)  
        indx = int(index2/2)
        if not numlist[0] < num < numlist[-1]:
            raise Exception('{0} is not in {1}'.format(str(num), str(numlist)))
            
        if index2 == 2:
            l1, l2 = num-numlist[0], numlist[-1]-num
            
            if l1 < l2:
                i = i
                
            else:
                i = i+1
                
        elif num == numlist[indx]:
            i = i + indx
            
        elif num > numlist[indx]:
            i = self.__closest_num(num, numlist[indx:],
                              i=i+indx)
                              
        elif num < numlist[indx]:
            i = self.__closest_num(num, numlist[0:indx+1], i=i)
        return i
        
    def get_data(self, url):
        '''
        return the data needed.
        url is from get_roms.get_url(starttime, endtime)
        '''
        data = get_nc_data(url, 'lon_rho', 'lat_rho', 'mask_rho','u', 'v', 'h', 's_rho')
        '''depth_layers = data['h'][index[0]][index[1]]*data['s_rho']  #[0]][index[1][0]]
        layer = np.argmin(abs(depth_layers+depth))
        print layer#'''
        layer = 35
        lon_rho = data['lon_rho'][:]
        lat_rho = data['lat_rho'][:]
        u = data['u'][:,layer]; v = data['v'][:,layer]
        return lon_rho,lat_rho,u,v
        
    def get_track(self,lon,lat,lon_rho,lat_rho,u,v,track_way):#, depth
        '''
        get the nodes of specific time period
        lon, lat: start point
        depth: 0~35, the 36th is the bottom.
        '''
        nodes = dict(lon=[lon], lat=[lat])
        try:
            index,lonp,latp,nearestdistance = self.nearest_point_index(lon,lat,lon_rho,lat_rho)
        except:
            return nodes
        t = abs(self.hours)
        for i in xrange(t):  #Roms points update every 2 hour
            if track_way=='backward': # backwards case
                u_t = -1*u[(t-1)-i][index[0],index[1]] #[0][0]][index[1][0]]
                v_t = -1*v[(t-1)-i][index[0],index[1]]
            else:
                u_t = u[i][index[0],index[1]] #[index[0]][index[1]]#[0][0]][index[1][0]]
                v_t = v[i][index[0],index[1]] #[index[0]][index[1]]#[0][0]][index[1][0]]
            dx = 60*60*u_t#float(u_p)
            dy = 60*60*v_t#float(v_p)
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')                        
            #x,y = mapx(lon,lat)
            #lon,lat = mapx(x+dx,y+dy,inverse=True)            
            lon = lon + dx/(111111*np.cos(lat*np.pi/180))
            lat = lat + dy/111111
            print 'lon,lat,i',lon,lat,i
            try:
                index,lonp,latp,nearestdistance = self.nearest_point_index(lon,lat,lon_rho,lat_rho)
                print 'distance',nearestdistance
            except:
                return nodes
            nodes['lon'].append(lon);  nodes['lat'].append(lat)
        return nodes

class get_fvcom(track):
    def __init__(self, mod):
        self.modelname = mod
    def nearest_point_index(self, lon, lat, lons, lats,rad):  #,num=4
        '''
        Return the nearest point(lonp,latp) and distance to origin point(lon,lat).
        lon, lat: the coordinate of start point, float
        latp, lonp: the coordinate of points to be calculated.
        '''
        def bbox2ij(lon, lat, lons, lats, rad):  
            
            p = Path.circle((lon,lat),radius=rad)
            points = np.vstack((lons,lats)).T  #numpy.vstack(tup):Stack arrays in sequence vertically
            
            inside = []
            
            for i in range(len(points)):
                inside.append(p.contains_point(points[i]))  # .contains_point return 0 or 1
                
            sidex = np.array(inside, dtype=bool)#.reshape(tshape)
            index = np.where(sidex==True)
            
            '''check if there are no points inside the given area'''        
            
            if not index[0].tolist():          # bbox covers no area
                #print 'This point is out of the model area.'
                #print '<h2>"This point is out of the model area."</h2>'
                raise Exception()
                
            else:
                return index        
        
        def min_distance(lon,lat,lons,lats):
            '''Find out the nearest distance to (lon,lat),and return lon.distance units: meters'''
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')
            dis_set = []
            #x,y = mapx(lon,lat)
            for i,j in zip(lons,lats):
                #x2,y2 = mapx(i,j)
                ss=math.sqrt((lon-i)**2+(lat-j)**2)
                #ss=math.sqrt((x-x2)**2+(y-y2)**2)
                dis_set.append(ss)
            dis = min(dis_set)
            p = dis_set.index(dis)
            lonp = lons[p]; latp = lats[p]
            return lonp,latp,dis       
        index = bbox2ij(lon, lat, lons, lats,rad)
        lon_covered = lons[index];  lat_covered = lats[index]       
        lonp,latp,distance = min_distance(lon,lat,lon_covered,lat_covered)
        #index1 = np.where(lons==lonp)
        #index2 = np.where(lats==latp)
        #index = np.intersect1d(index1,index2)        
        #points = np.vstack((lons.flatten(),lats.flatten())).T 
        #index = [i for i in xrange(len(points)) if ([lonp,latp]==points[i]).all()]
       
        return lonp,latp,distance  #,lonp,latp
        
    def get_url(self, starttime, endtime):
        '''
        get different url according to starttime and endtime.
        urls are monthly.
        '''
        self.hours = int(round((endtime-starttime).total_seconds()/60/60))
        #print self.hours
                
        if self.modelname == "GOM3":
            url = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?
            lon[0:1:51215],lat[0:1:51215],lonc[0:1:95721],latc[0:1:95721],siglay[0:1:39][0:1:51215],h[0:1:51215],nbe[0:1:2][0:1:95721]'''
            urll = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?
            u[{0}:1:{1}][0:1:39][0:1:95721],v[{0}:1:{1}][0:1:39][0:1:95721],zeta[{0}:1:{1}][0:1:51215]'''
            current_time = pytz.utc.localize(datetime.now().replace(hour=0,minute=0,second=0,microsecond=0))
            period = starttime-(current_time-timedelta(days=3))
            if abs((starttime-current_time).total_seconds())>259200 or abs((endtime-current_time).total_seconds())>259200: #24*3600*3
                raise IndexError('GOM3 only works between 3days before and 3daysafter.')
            '''if period.total_seconds()<0:
                raise IndexError('GOM3 only works between 3days before and 3daysafter.')#'''
            index1 = int(round(period.total_seconds()/3600))
            index2 = index1 + self.hours
            self.url = urll.format(index1, index2)
            
        elif self.modelname == "massbay":
            url = """http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?
            lon[0:1:98431],lat[0:1:98431],lonc[0:1:165094],latc[0:1:165094],siglay[0:1:9][0:1:98431],h[0:1:98431],
            nbe[0:1:2][0:1:165094]"""
            urll = """http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?
            u[{0}:1:{1}][0:1:9][0:1:165094],v[{0}:1:{1}][0:1:9][0:1:165094],zeta[{0}:1:{1}][0:1:98431]"""
            current_time = pytz.utc.localize(datetime.now().replace(hour=0,minute=0,second=0,microsecond=0))
            #print 'current_time',current_time
            period = starttime-(current_time-timedelta(days=3))
            if abs((starttime-current_time).total_seconds())>259200 or abs((endtime-current_time).total_seconds())>259200: #24*3600*3
                raise IndexError('GOM3 only works between 3days before and 3daysafter.')
            index1 = int(round(period.total_seconds()/3600))
            index2 = index1 + self.hours
            self.url = urll.format(index1, index2)#'''

        elif self.modelname == "massbaya":
            url = '''http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_mb?
            lon[0:1:98431],lat[0:1:98431],lonc[0:1:165094],latc[0:1:165094],siglay[0:1:9][0:1:98431],
            h[0:1:98431],u[{0}:1:{1}][0:1:9][0:1:165094],v[{0}:1:{1}][0:1:9][0:1:165094]'''
            index1 = int((starttime-datetime(2011,1,18,0,0,0,0,pytz.UTC)).total_seconds()/3600)
            index2 = index1 + self.hours
            if index2<index1: #case of backwards run
                url = url.format(index2, index1)
            else:
                url = url.format(index1, index2)
            
        elif self.modelname == "GOM3a":
            url = '''http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_gom3v13?
            lon[0:1:51215],lat[0:1:51215],lonc[0:1:95721],latc[0:1:95721],siglay[0:1:39][0:1:51215],
            h[0:1:51215],u[{0}:1:{1}][0:1:39][0:1:95721],v[{0}:1:{1}][0:1:39][0:1:95721]'''
            index1 = int((starttime-datetime(2013,5,9,0,0,0,0,pytz.UTC)).total_seconds()/3600)
            index2 = index1 + self.hours
            url = url.format(index1, index2)
            
        elif self.modelname == "30yr": #start at 1977/12/31 22:58, end at 2014/1/1 0:0, time units:hours
            url = """http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?h[0:1:48450],
            lat[0:1:48450],latc[0:1:90414],lon[0:1:48450],lonc[0:1:90414],nbe[0:1:2][0:1:90414],
            siglay[0:1:44][0:1:48450],nv[0:1:2][0:1:90414]"""
            urll = """http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?u[{0}:1:{1}][0:1:44][0:1:90414],
            v[{0}:1:{1}][0:1:44][0:1:90414],zeta[{0}:1:{1}][0:1:48450]"""
            if (starttime-datetime(2014,1,1,0,0,0,0,pytz.UTC)).total_seconds()>0 or (endtime-datetime(2014,1,1,0,0,0,0,pytz.UTC)).total_seconds()>0:
                raise IndexError('"30yr" only works between 1988-2013.')
            index1 = int(round((starttime-datetime(1977,12,31,22,58,4,0,pytz.UTC)).total_seconds()/3600))
            index2 = index1 + self.hours
            self.url = urll.format(index1, index2)
        #print url
        return url

    def get_data(self,url):
        '''
        ??? Retrieves data?
        '''
        self.data = get_nc_data(url,'lat','lon','latc','lonc','siglay','h','nbe')#,'nv'
        lonc, latc = self.data['lonc'][:], self.data['latc'][:]  #quantity:165095
        lons, lats = self.data['lon'][:], self.data['lat'][:]
        h = self.data['h'][:]; siglay = self.data['siglay'][:]; #nv = self.data['nv'][:]
        nbe1=self.data['nbe'][0];nbe2=self.data['nbe'][1];
        nbe3=self.data['nbe'][2]
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
        b_points = np.vstack((lonb,latb)).T#'''
        self.b_points = b_points
        self.lons=lons;self.lats=lats;self.lonc=lonc;self.latc=latc;self.h=h;self.siglay=siglay
        return b_points #lons,lats,lonc,latc,h,siglay #,nv 
        
    def shrink_data(self,lon,lat,lons,lats):
        lont = []; latt = []
        p = Path.circle((lon,lat),radius=0.5)
        pints = np.vstack((lons,lats)).T
        for i in range(len(pints)):
            if p.contains_point(pints[i]):
                lont.append(pints[i][0])
                latt.append(pints[i][1])
        lonl=np.array(lont); latl=np.array(latt)#'''
        return lonl,latl
        
    def eline_path(self,lon,lat):
        p = Path.circle((lon,lat),radius=0.1) #0.06
        dis = []; bps = []
        for i in self.b_points:
            if p.contains_point(i):
                bps.append((i[0],i[1]))
                d = math.sqrt((lon-i[0])**2+(lat-i[1])**2)
                dis.append(d)
        if len(dis)<3 :
            return None
        dnp = np.array(dis)
        dis.sort()
        if dis[0]>0.04 :
            return None
        
        else :
            cdis = []; cbps = []
            dis0 = dis[0]
            p = np.where(dnp==dis0)   
            bps0 = bps[p[0]]
            p1 = Path.circle(bps0,radius=0.04)
            for j in bps:
                if p1.contains_point(j):
                    cbps.append((j[0],j[1]))
                    d1 = math.sqrt((lon-j[0])**2+(lat-j[1])**2)
                    cdis.append(d1)
            if len(cdis)<3 :
                return None
            dnp1 = np.array(cdis)
            cdis.sort()            
            cdis1 = cdis[1]; cdis2 = cdis[2]
            p1 = np.where(dnp1==cdis1); p2 = np.where(dnp1==cdis2)
            bps1 = cbps[p1[0]]; bps2 = cbps[p2[0]]
            pa = [bps1,bps0,bps2]; #print 'pa',pa
            #codes = [Path.MOVETO,Path.LINETO,Path.LINETO]
            path = Path(pa)#,codes
            return path
    
    def uvt(self,u1,v1,u2,v2):
        t = 2
        a=0; b=0
        if u1==u2:
            a = u1
        else:
            ut = np.arange(u1,u2,float(u2-u1)/t)
            for i in ut:
                a += i
            a = a/t  
        
        if v1==v2:
            b = v1
        else:
            c = float(v2-v1)/t
            vt = np.arange(v1,v2,c)
            for i in vt:
                b += i
            b = b/t
               
        return a, b
        
    def get_track(self,lon,lat,depth,track_way): #,b_index,nvdepth, 
        '''
        Get forecast nodes start at lon,lat
        '''
        nodes = dict(lon=[lon], lat=[lat])
        uvz = netCDF4.Dataset(self.url)
        u = uvz.variables['u']; v = uvz.variables['v']; zeta = uvz.variables['zeta']
        
        if lon>90:
            lon, lat = dm2dd(lon, lat)
        lonl,latl = self.shrink_data(lon,lat,self.lonc,self.latc)
        lonk,latk = self.shrink_data(lon,lat,self.lons,self.lats)
        try:
            if self.modelname == "GOM3" or self.modelname == "30yr":
                lonp,latp,distance = self.nearest_point_index(lon, lat, lonl, latl,0.2)
                lonn,latn,ds = self.nearest_point_index(lon,lat,lonk,latk,0.3)
                
            if self.modelname == "massbay":
                lonp,latp,distance = self.nearest_point_index(lon, lat, lonl, latl,0.03)
                lonn,latn,ds = self.nearest_point_index(lon,lat,lonk,latk,0.05)
        except:
            #print '<script type="text/javascript">window.alert("Start point is out of the Model(%s) grid.")</script>'%self.modelname            
            return nodes
            
        index1 = np.where(self.lonc==lonp)
        index2 = np.where(self.latc==latp)
        elementindex = np.intersect1d(index1,index2)
        index3 = np.where(self.lons==lonn)
        index4 = np.where(self.lats==latn)
        nodeindex = np.intersect1d(index3,index4)
        ################## boundary 1 ####################            
        
        pa = self.eline_path(lon,lat)
        
        if track_way=='backward' : # backwards case
            waterdepth = self.h[nodeindex]+zeta[-1,nodeindex]
        else:
            waterdepth = self.h[nodeindex]+zeta[0,nodeindex]
        if waterdepth<(abs(depth)): 
            #print 'This point is too shallow.Less than %d meter.'
            print '<script type="text/javascript">window.alert("Start point is too shallow. Less than %d meter.")</script>'%abs(depth)
            #raise Exception()
            return nodes
        depth_total = self.siglay[:,nodeindex]*waterdepth  
        layer = np.argmin(abs(depth_total-depth))
        
            
        t = abs(self.hours)
                 
        for i in range(t):            
            if i!=0 and i%24==0 :
                #print 'layer,lon,lat,i',layer,lon,lat,i
                lonl,latl = self.shrink_data(lon,lat,self.lonc,self.latc)
                lonk,latk = self.shrink_data(lon,lat,self.lons,self.lats)
            if track_way=='backward' : # backwards case
                u_t1 = u[t-i,layer,elementindex][0]*(-1); v_t1 = v[t-i,layer,elementindex][0]*(-1)
                u_t2 = u[t-i-1,layer,elementindex][0]*(-1); v_t2 = v[t-i-1,layer,elementindex][0]*(-1)
            else:
                u_t1 = u[i,layer,elementindex][0]; v_t1 = v[i,layer,elementindex][0]
                u_t2 = u[i+1,layer,elementindex][0]; v_t2 = v[i+1,layer,elementindex][0]
            u_t,v_t = self.uvt(u_t1,v_t1,u_t2,v_t2)          
            '''if u_t==0 and v_t==0: #There is no water
                print 'Sorry, hits the land,u,v==0'
                return nodes,1 #'''
            dx = 60*60*u_t; dy = 60*60*v_t
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')                        
            #x,y = mapx(lon,lat)
            #lon,lat = mapx(x+dx,y+dy,inverse=True)            
            temlon = lon + (dx/(111111*np.cos(lat*np.pi/180)))
            temlat = lat + dy/111111
            nodes['lon'].append(temlon);  nodes['lat'].append(temlat)
            #print ''' pts[%d] = new GLatLng(%f,%f);'''%(i+1,temlat,temlon)
            #########case for boundary 1 #############
            if pa:
                teml = [(lon,lat),(temlon,temlat)]
                tempa = Path(teml)
                if pa.intersects_path(tempa): 
                    #print 'Sorry, point hits land here.path'
                    return nodes #'''
                
            #########################
            lon = temlon; lat = temlat
            
            if i!=(t-1):                
                try:
                    if self.modelname == "GOM3" or self.modelname == "30yr":
                        lonp,latp,distance = self.nearest_point_index(lon, lat, lonl, latl,0.2)
                        lonn,latn,ds = self.nearest_point_index(lon,lat,lonk,latk,0.3)
                    if self.modelname == "massbay":
                        lonp,latp,distance = self.nearest_point_index(lon, lat, lonl, latl,0.03)
                        lonn,latn,ds = self.nearest_point_index(lon,lat,lonk,latk,0.05)
                    index1 = np.where(self.lonc==lonp)
                    index2 = np.where(self.latc==latp)
                    elementindex = np.intersect1d(index1,index2);#print 'elementindex',elementindex
                    index3 = np.where(self.lons==lonn)
                    index4 = np.where(self.lats==latn)
                    nodeindex = np.intersect1d(index3,index4)
                    
                    ################## boundary 1 ####################
            
                    pa = self.eline_path(lon,lat)
                                      
                    if track_way=='backward' : # backwards case
                        waterdepth = self.h[nodeindex]+zeta[t-i-1,nodeindex]
                    else:
                        waterdepth = self.h[nodeindex]+zeta[i+1,nodeindex]
                    #print 'waterdepth',h[nodeindex],zeta[i+1,nodeindex],waterdepth
                    if waterdepth<(abs(depth)): 
                        #print 'This point hits the land here.Less than %d meter.'%abs(depth)
                        raise Exception()
                    depth_total = self.siglay[:,nodeindex]*waterdepth  
                    layer = np.argmin(abs(depth_total-depth))           
                except:
                    #print 'window.alert("Sorry, drifter hits the land")</script><script type="text/javascript">'
                    return nodes
                                
        return nodes
                

class get_drifter(track):

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

def draw_basemap(ax, points, interval_lon=0.3, interval_lat=0.3):
    '''
    draw the basemap?
    '''
    
    lons = points['lons']
    lats = points['lats']
    size = max((max(lons)-min(lons)),(max(lats)-min(lats)))/1
    map_lon = [min(lons)-size,max(lons)+size]
    map_lat = [min(lats)-size,max(lats)+size]
    
    #ax = fig.sca(ax)
    dmap = Basemap(projection='cyl',
                   llcrnrlat=map_lat[0], llcrnrlon=map_lon[0],
                   urcrnrlat=map_lat[1], urcrnrlon=map_lon[1],
                   resolution='h',ax=ax)# resolution: c,l,i,h,f.
    dmap.drawparallels(np.arange(int(map_lat[0])-1,
                                 int(map_lat[1])+1,interval_lat),
                       labels=[1,0,0,0])
    dmap.drawmeridians(np.arange(int(map_lon[0])-1,
                                 int(map_lon[1])+1,interval_lon),
                       labels=[0,0,0,1])
    dmap.drawcoastlines()
    dmap.fillcontinents(color='grey')
    dmap.drawmapboundary()

def uniquecolors(N):
    """
    Generate unique RGB colors
    input: number of RGB tuples to generate
    output: list of RGB tuples
    """
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    colors =  [colorsys.hsv_to_rgb(x*1.0/N, 0.5, 0.5) for x in range(N)]
    return colors

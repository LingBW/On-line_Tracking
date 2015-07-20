#!/usr/bin/env /anaconda/bin/python
import cgitb
cgitb.enable()
import cgi,sys
import pytz
from dateutil.parser import parse
from datetime import datetime, timedelta
from track_functions import get_fvcom,get_roms

print """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
<head>
<script src="http://maps.googleapis.com/maps/api/js" type="text/javascript"></script>
<title>Simulation Track</title>
"""

form = cgi.FieldStorage()
lons = form.getvalue('lons')
lats = form.getvalue('lats')
tds = form.getvalue('tds')
track_way = form.getvalue('tway') # Two option: backward, froward; 
stt = form.getvalue('stt')
dep = form.getvalue('depth')
dset = form.getvalue('dset')

st_lon=[]; st_lat=[]
try:
    lon=lons.split(',')
    lat=lats.split(',')
    for i in zip(lon,lat):
        st_lat.append(float(i[1]))
        st_lon.append(float(i[0]))
except:
    print '<h2>"No start points.%s"</h2>'%lons
    print '</head></html>'
    sys.exit()

utcti = datetime.now(pytz.UTC); utct = utcti.strftime('%H')
locti = datetime.now(); loct = locti.strftime('%H')
ditnu = int(utct)-int(loct)
try:
    if not stt:
        start_time = datetime.now(pytz.UTC)#edtl = (edt[0],edt[1],edt[2],edt[3],edt[4])
    else:
        stime = parse(stt) #datetime(2015,2,10,12,0,0,0,pytz.UTC        
        start_time = pytz.utc.localize(stime)+timedelta(hours=ditnu)
except:
    print '<h2>"Time format error."</h2>'  
    print '</head></html>'
    sys.exit() 
stp_num = len(st_lat)

MODEL = dset      # 'FVCOM', 'ROMS'
GRID = ['massbay','GOM3','30yr']  
forecast_days = float(tds)      #MODEL track time(days) 
depth = float(dep)*(-1)   
end_time = start_time + timedelta(forecast_days)
colors = ['magenta','olive','blue','orange','green','red','yellow','black','purple']
if track_way=='backward':
    end_time = start_time 
    start_time = end_time - timedelta(forecast_days)  #'''

stim = (start_time-timedelta(hours=ditnu)).strftime('%D %H:%M')
etim = (end_time-timedelta(hours=ditnu)).strftime('%D %H:%M')
curtime = pytz.utc.localize(datetime.now().replace(hour=0,minute=0,second=0,microsecond=0))
mst = curtime - timedelta(3); met = curtime + timedelta(3)
mstt = mst.strftime('%D %H:%M'); mett = met.strftime('%D %H:%M')
lon_set = [[]]*stp_num; lat_set = [[]]*stp_num
if MODEL in GRID:
    get_obj = get_fvcom(MODEL)
    try:
        url_fvcom = get_obj.get_url(start_time,end_time)
    except:
        if MODEL == '30yr':
            print '<h2>"Start time: Error! Model(30yr) only works from 01/01/1978 to 01/01/2014."</h2>'
        else: 
            print '<h2>Start time: Error! Model(%s) only works the past and future 3 days(%s - %s). <br>Tips: check step 5 or 6.</h2>'%(MODEL,mstt,mett)
        print '</head></html>'
        sys.exit()
    bpoints = get_obj.get_data(url_fvcom)        
    for i in range(stp_num):
        point = get_obj.get_track(st_lon[i],st_lat[i],depth,track_way)
        if len(point['lon'])==1: 
            print '<script type="text/javascript">window.alert("This point is out of the Model(%s) Zone. Tips: change the Model(step 2), try again.")</script>'%MODEL
        lon_set[i] = point['lon']; lat_set[i] = point['lat']        
        
if MODEL=='ROMS': 
    get_obj = get_roms()
    try:
        url_roms = get_obj.get_url(start_time,end_time)
    except:
        print '<h2>"Time Error! Track time out of Model time horizon."</h2></head></html>' 
        sys.exit()
    lon_rho,lat_rho,u,v = get_obj.get_data(url_roms)
    for i in range(stp_num):
        point = get_obj.get_track(st_lon[i],st_lat[i],lon_rho,lat_rho,u,v,track_way)
        if len(point['lon'])==1: 
            print '<script type="text/javascript">window.alert("This point is out of the Model(%s) Zone. Tips: change the Model(step 2), try again.")</script>'%MODEL
        lon_set[i] = point['lon']; lat_set[i] = point['lat']
        
        

print """
<style type="text/css">
body {
    background-image: url(http://comet.nefsc.noaa.gov/ioos/track/20150412.jpg);
    background-repeat: repeat;
    background-position: top center;
    background-attachment: scroll;
}
#header {
    background-color: blue;
    color:white;
    text-align:center;
    padding: 1px;
    width: 958px;
    height: ;
    margin: 0 auto;
}
#googleMap {
    width: 960px;
    height: 500px;
    margin: 0 auto;
    border: ;
}
#footer {
	background-color: blue;
	color:white;
	clear:both;
	text-align:center;
	padding: 5px;
	width: 950px;
	height: ;
	margin: auto;
} 
</style>

</head>

<body onunload="GUnload()">

<script type="text/javascript">
var myCenter=new google.maps.LatLng(%s,%s);
var stim = "%s";
var etim = "%s";
function initialize() {
    var mapProp = {
        center: myCenter,
        zoom:7,
        mapTypeId: google.maps.MapTypeId.SATELLITE
    };
  
    var map = new google.maps.Map(document.getElementById("googleMap"),mapProp);
    
    function startmarker(tlon,tlat) {
        //Show a point(lon,lat) on the map, and add an event click the marker show the point info.
        var pmarker = new google.maps.LatLng(tlon,tlat);
        var marker = new google.maps.Marker({ 
	        position: pmarker,
	        //animation:google.maps.Animation.BOUNCE
	        title:'Start point',
	        icon: 'http://comet.nefsc.noaa.gov/ioos/track/startmarker.png'
        	});     
        marker.setMap(map);   
        var infowindow = new google.maps.InfoWindow({
            content: 'Latitude: ' + parseFloat(pmarker.lat()).toFixed(6) +'<br>Longitude: ' + parseFloat(pmarker.lng()).toFixed(6)+'<br>StartTime: '+stim+' EDT'});
          
        google.maps.event.addListener(marker, 'click', function() {
            map.setZoom(10);
            map.setCenter(marker.getPosition());
            infowindow.open(map,marker);
        	});
    	}
    function endmarker(tlon,tlat) {
        //Show a point(lon,lat) on the map, and add an event click the marker show the point info.
        var pmarker = new google.maps.LatLng(tlon,tlat);
        var marker = new google.maps.Marker({ 
	        position: pmarker,
	        //animation:google.maps.Animation.BOUNCE
	        title:'End point',
	        icon: 'http://comet.nefsc.noaa.gov/ioos/track/endmarker.png'
        	});     
        marker.setMap(map);   
        var infowindow = new google.maps.InfoWindow({
            content: 'Latitude: ' + parseFloat(pmarker.lat()).toFixed(6) +'<br>Longitude: ' + parseFloat(pmarker.lng()).toFixed(6)+'<br>EndTime: '+etim+' EDT'});
          
        google.maps.event.addListener(marker, 'click', function() {
            map.setZoom(10);
            map.setCenter(marker.getPosition());
            infowindow.open(map,marker);
        	});
    	}"""%(st_lat[0],st_lon[0],stim,etim)

####### polyline #########
for i in range(stp_num):
    # plot start markers
    print 'startmarker(%s,%s)'%(st_lat[i],st_lon[i])
    # plot end markers
    print 'endmarker(%s,%s)'%(lat_set[i][-1],lon_set[i][-1])
    print 'var pts=[]'
    for j in xrange(len(lon_set[i])):
        print 'pts[%d] = new google.maps.LatLng(%s,%s);'%(j,lat_set[i][j],lon_set[i][j])
    # plot line
    print """
    var ptsPath=new google.maps.Polyline({
        path:pts,
        strokeColor:"%s",
        strokeOpacity:0.8,
        strokeWeight:4
        });
    ptsPath.setMap(map);"""%colors[i%9]

print """
}
google.maps.event.addDomListener(window, 'load', initialize);
</script>
   
<div id="header"><h2>Model: %s, %s track %s days.</h2></div>

<div id="googleMap" ></div>

<div id="footer" ><a href="http://comet.nefsc.noaa.gov/ioos/track/index.html"><input type="button" value="Track again" /></a></div>
"""%(MODEL,track_way,tds)

print '</script></body></html>'

#!/usr/bin/env python
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
    print '<h2>"No start points."</h2>'
    print '</head></html>'
    sys.exit()

try:
    if not stt:
        start_time = datetime.now(pytz.UTC)#edtl = (edt[0],edt[1],edt[2],edt[3],edt[4])
    else:
        stime = parse(stt) #datetime(2015,2,10,12,0,0,0,pytz.UTC
        start_time = pytz.utc.localize(stime)
except:
    print '<h2>"Time format error."</h2>'  
    print '</head></html>'
    sys.exit() 
stp_num = len(st_lat)

MODEL = 'FVCOM'      # 'FVCOM', 'ROMS'
GRID = dset     # Apply to FVCOM. '30yr', 'massbaya', 'GOM3a', 'GOM3' or 'massbay'
forecast_days = float(tds)      #MODEL track time(days) 
depth = float(dep)*(-1)   
end_time = start_time + timedelta(forecast_days)
colors = ['magenta','olive','blue','orange','green','red','yellow','black','purple']
if track_way=='backward':
    end_time = start_time 
    start_time = end_time - timedelta(forecast_days)  #'''

lon_set = [[]]*stp_num; lat_set = [[]]*stp_num
if MODEL=='FVCOM':
    get_obj = get_fvcom(GRID)
    try:
        url_fvcom = get_obj.get_url(start_time,end_time)
    except:
        print 'window.alert("Model not works this duration.");'
        print '</head></html>'
        sys.exit()
    bpoints = get_obj.get_data(url_fvcom)        
    for i in range(stp_num):
        point = get_obj.get_track(st_lon[i],st_lat[i],depth,track_way)
        lon_set[i] = point['lon']; lat_set[i] = point['lat']        
        
if MODEL=='ROMS': 
    get_obj = get_roms()
    url_roms = get_obj.get_url(start_time,end_time)
    lon_rho,lat_rho,u,v = get_obj.get_data(url_roms)
    for i in range(stp_num):
        point = get_obj.get_track(st_lon[i],st_lat[i],lon_rho,lat_rho,u,v,track_way)
        lon_set[i] = point['lon']; lat_set[i] = point['lat']       
        

print """
<style type="text/css">
    a:link {}
    a:visited {color: ;}
    a:hover {color: ;}
    a:active {color: #900;}
    body {
        background-image: url(http://127.0.0.1:8000/image/20150412.jpg);
        background-repeat: repeat;
        background-position: top center;
        background-attachment: scroll;
    }
    #header {
        width: auto;
        height: 70px;
        margin: 0 auto;
        border: ;
    }
    #container {
        width: 1024px;
        height: 600px;
        margin: 0 auto;
        border: ;
    }
    #googleMap {
        width: 800px;
        height: 596px;
        float: left;
        border: black 2px solid;
    }
    #sidebar {
        width: 200px;
        height: 580px;
        border: ;
        float: left;
        background: ;
        overflow: auto;
        padding: 10px;
        text-align: center;      
    }
    #ad {   
        text-align: center;   
        font-size: 20px;
        vertical-align: middle;
        text-decoration: none;
        padding: 4px 10px;
        border: 1px solid;
        background: #eeddbb;
        margin: auto;
    }

</style>
</head>

<body onunload="GUnload()">

<script type="text/javascript">
var myCenter=new google.maps.LatLng(41.8,-70.3);

function initialize()
{
  var mapProp = {
  center: myCenter,
  zoom:9,
  mapTypeId: google.maps.MapTypeId.ROADMAP
  };
  
  var map = new google.maps.Map(document.getElementById("googleMap"),mapProp);


function setmarker(tlon,tlat) {
  //Show a point(lon,lat) on the map, and add an event click the marker show the point info.
  var pmarker = new google.maps.LatLng(tlon,tlat);
  var marker = new google.maps.Marker({
  position: pmarker,
  animation:google.maps.Animation.BOUNCE
  //title:'Click to zoom'
  });
  
  marker.setMap(map);
  
  var infowindow = new google.maps.InfoWindow({
  content: 'Latitude: ' + pmarker.lat() +
  '<br>Longitude: ' + pmarker.lng()
  });
  
  google.maps.event.addListener(marker, 'click', function() {
  infowindow.open(map,marker);
  });
}
"""

for i in range(stp_num):
    print 'setmarker(%s,%s)'%(st_lat[i],st_lon[i])
    print 'var pts=[]'
    for j in xrange(len(lon_set[i])):
        print 'pts[%d] = new google.maps.LatLng(%s,%s);'%(j,lat_set[i][j],lon_set[i][j])
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
   
<div id="header"><h2 align="center">%s track %s days.</h2></div>

<div id="container">
    <div id="googleMap" ></div>
    <div id="sidebar" ><a id="ad" href="http://127.0.0.1:8000/index.html">Track again</a></div>
</div>
"""%(track_way,tds)

print '</script></body></html>'

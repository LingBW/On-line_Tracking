#!/usr/bin/env /anaconda/bin/python
import cgitb
cgitb.enable()
import cgi,sys
from datetime import timedelta
from track_functions import get_fvcom,get_drifter

print """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
<head>
<script src="http://maps.googleapis.com/maps/api/js" type="text/javascript"></script>
<title>Drifter Track</title>
"""
form = cgi.FieldStorage()
did = form.getvalue('did')
tds = form.getvalue('tds')
dty = form.getvalue('dty')
track_way = form.getvalue('tway') # Two option: backward, froward; 
dset = form.getvalue('dset')
dep = form.getvalue('depth')

MODEL = 'FVCOM'      # 'FVCOM', 'ROMS'
GRID = dset     # Apply to FVCOM. '30yr', 'massbaya', 'GOM3a', 'GOM3' or 'massbay'
drifter_ID = did
forecast_days = float(tds)      #MODEL track time(days) 
depth = float(dep)*(-1)
INPUT_DATA = dty
start_time = None
drifter_days = None

try:
    drifter = get_drifter(drifter_ID, INPUT_DATA)
    dr_points = drifter.get_track(start_time,drifter_days)
except:
    print '<h2>"Drifter ID or Data Type Error."</h2>'
    print '</head></html>'
    sys.exit()

if MODEL=='FVCOM':
    get_obj =  get_fvcom(GRID)
    try:
        url_fvcom = get_obj.get_url(dr_points['time'][-1],dr_points['time'][-1]+timedelta(forecast_days))
    except:
        print '<h2>"Start time: Error! The last updated time of drifter not works model duration."</h2>'
        print '</head></html>'
        sys.exit()
    else:
        bpoints = get_obj.get_data(url_fvcom) # b_points is model boundary points.
        point1 = get_obj.get_track(dr_points['lon'][-1],dr_points['lat'][-1],depth,track_way)
#'''

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

<body>

<script type="text/javascript">

var myCenter=new google.maps.LatLng(%s,%s);

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
  //animation:google.maps.Animation.BOUNCE
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


"""%(dr_points['lat'][-1],dr_points['lon'][-1])

#add markers
print "setmarker(%s,%s);"%(dr_points['lat'][-1],dr_points['lon'][-1])

#add drifter lines
print "var dts=[]"
for i in xrange(len(dr_points['lon'])):
    print 'dts[%d] = new google.maps.LatLng(%s,%s);'%(i,dr_points['lat'][i],dr_points['lon'][i])
print """
var flightPath=new google.maps.Polyline({
  path:dts,
  strokeColor:"blue",
  strokeOpacity:0.8,
  strokeWeight:4
  });

flightPath.setMap(map);"""
#add model lines
print 'var pts=[]'
for j in xrange(len(point1['lon'])):
    print 'pts[%d] = new google.maps.LatLng(%s,%s);'%(j,point1['lat'][j],point1['lon'][j])
print """
var ptsPath=new google.maps.Polyline({
    path:pts,
    strokeColor:"red",
    strokeOpacity:0.8,
    strokeWeight:4
    });
ptsPath.setMap(map);"""


print """
}
google.maps.event.addDomListener(window, 'load', initialize);

</script>

<div id="header"><h2 align="center">%s track %s days.</h2></div>

<div id="container">
    <div id="googleMap" ></div>
    <div id="sidebar" ><a id="ad" href="http://comet.nefsc.noaa.gov/ioos/track/index2.html">Track again</a></div>
</div>

"""%(track_way,tds)

print '</body></html>'

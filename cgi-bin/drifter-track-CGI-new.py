#!/usr/bin/env python
import cgitb
cgitb.enable()
import cgi,sys
#import pytz
from datetime import datetime,timedelta
from track_functions import get_fvcom,get_roms,get_drifter

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

MODEL = dset      # 'FVCOM', 'ROMS'
GRID = ['massbay','GOM3','30yr'] 
drifter_ID = did
forecast_days = float(tds)      #MODEL track time(days) 
depth = abs(float(dep))
INPUT_DATA = 'drift_X.dat'  #dty
start_time = None
drifter_days = None

try:
    drifter = get_drifter(drifter_ID, INPUT_DATA)
    dr_points = drifter.get_track(start_time,drifter_days)
except:
    print '<h2>"Drifter ID Error!"</h2>'
    print '</head></html>'
    sys.exit()

ltime = dr_points['time'][-1].strftime('%D %H:%M')
#curtime = pytz.utc.localize(datetime.now().replace(hour=0,minute=0,second=0,microsecond=0))
#mst = curtime - timedelta(3); met = curtime + timedelta(3)
#mstt = mst.strftime('%D %H:%M'); mett = met.strftime('%D %H:%M')
start_time = dr_points['time'][-1]
end_time = start_time + timedelta(forecast_days)
if track_way=='backward':
    end_time = start_time 
    start_time = end_time - timedelta(forecast_days)  #'''
    
utcti = datetime.utcnow(); utct = utcti.strftime('%H')
locti = datetime.now(); loct = locti.strftime('%H')
ditnu = int(utct)-int(loct) # the deference between UTC and local time 
if ditnu < 0:
    ditnu = int(utct)+24-int(loct)

# pass to html to show start-time and end-time of marker in ET. 
sloctim = start_time-timedelta(hours=ditnu)
eloctim = end_time-timedelta(hours=ditnu)
stim = sloctim.strftime('%D %H:%M')
etim = eloctim.strftime('%D %H:%M')
#stim = (start_time-timedelta(hours=ditnu)).strftime('%D %H:%M')
#etim = (end_time-timedelta(hours=ditnu)).strftime('%D %H:%M')

if MODEL in GRID:
    get_obj =  get_fvcom(MODEL)
    if (MODEL == 'massbay' or MODEL == 'GOM3') and depth < 2 :
        bpoints = get_obj.sf_get_data(start_time,end_time)
        
        point = get_obj.sf_get_track(dr_points['lon'][-1],dr_points['lat'][-1],track_way)
        if len(point['lon'])==1: 
            print '<script type="text/javascript">window.alert("This point on the land.")</script>'
    else:
        try:
            url_fvcom,fmtime,emtime = get_obj.get_data(start_time,end_time)
        except:
            print '<h2>Model(FVCOM:%s) is unavailable temporarily.Try to use another model.</h2>'%(MODEL)
            print '</head></html>'
            sys.exit()  
        if url_fvcom=='error':
            fmodtime = fmtime - timedelta(hours=ditnu)
            emodtime = emtime - timedelta(hours=ditnu)
            #print fmodtime,emodtime
            mstt = fmodtime.strftime('%m/%d/%Y %H:%M')
            mett = emodtime.strftime('%m/%d/%Y %H:%M')
    
            print '<h2>Time: Error! Model(FVCOM:%s) only works between %s with %s. <br>Tips: check step 5.</h2>'%(MODEL,mstt,mett)
            print '</head></html>'
            sys.exit()    
        
        point = get_obj.get_track(dr_points['lon'][-1],dr_points['lat'][-1],depth,track_way)
        if len(point['lon'])==1: 
            print '<script type="text/javascript">window.alert("The last point of drifter is out of the Model(%s) grid. Tips: change the Model grid(step 2), try again.")</script>'%MODEL

if MODEL=='ROMS': 
    get_obj = get_roms()
    
    try:
        url_roms,fmtime,emtime = get_obj.get_data(start_time,end_time)
    except:
        print '<h2>Model(%s) is unavailable temporarily.Try to use another model.</h2>'%(MODEL)
        print '</head></html>'
        sys.exit()
    if url_roms=='error':
        fmodtime = fmtime - timedelta(hours=ditnu)
        emodtime = emtime - timedelta(hours=ditnu)
        mstt = fmodtime.strftime('%m/%d/%Y %H:%M')
        mett = emodtime.strftime('%m/%d/%Y %H:%M')
        print '<h2>Time: Error! Model(%s) only works between %s with %s. <br>Tips: check step 5 or 6.</h2>'%(MODEL,mstt,mett)
        print '</head></html>'
        sys.exit()
    if depth < 2:
        point = get_obj.sf_get_track(dr_points['lon'][-1],dr_points['lat'][-1],track_way)
    else:
        point = get_obj.get_track(dr_points['lon'][-1],dr_points['lat'][-1],depth,track_way)
    if len(point['lon'])==1: 
        print '<script type="text/javascript">window.alert("Start point is out of the Model(%s) Zone. Tips: change the Model(step 2), try again.")</script>'%MODEL
print """
<style type="text/css">
body {
    background-image: url(http://127.0.0.1:8000/20150412.jpg);
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

<body>

<script type="text/javascript">

var myCenter=new google.maps.LatLng(%s,%s);
var stim = "%s";
var etim = "%s";
var ID = "%s";
function initialize() {
    var mapProp = {
        center: myCenter,
        zoom:11,
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
	        icon: 'http://127.0.0.1:8000/startmarker.png'
        	});     
        marker.setMap(map);   
        var infowindow = new google.maps.InfoWindow({
            content: 'Drifter ID: '+ID+'<br>Latitude: ' + parseFloat(pmarker.lat()).toFixed(6) +'<br>Longitude: ' + parseFloat(pmarker.lng()).toFixed(6)+'<br>StartTime: '+stim+' EDT'});
          
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
	        icon: 'http://127.0.0.1:8000/endmarker.png'
        	});     
        marker.setMap(map);   
        var infowindow = new google.maps.InfoWindow({
            content: 'Drifter ID: '+ID+'<br>Latitude: ' + parseFloat(pmarker.lat()).toFixed(6) +'<br>Longitude: ' + parseFloat(pmarker.lng()).toFixed(6)+'<br>EndTime: '+etim+' EDT'});
          
        google.maps.event.addListener(marker, 'click', function() {
            map.setZoom(10);
            map.setCenter(marker.getPosition());
            infowindow.open(map,marker);
        	});
    	}"""%(float(dr_points['lat'][-1]),float(dr_points['lon'][-1]),stim,etim,drifter_ID)

# plot markers
print 'startmarker(%s,%s)'%(float(dr_points['lat'][-1]),float(dr_points['lon'][-1]))
print 'endmarker(%s,%s)'%(point['lat'][-1],point['lon'][-1])
# plot drifter line
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
#plot model lines
print 'var pts=[]'
for j in xrange(len(point['lon'])):
    print 'pts[%d] = new google.maps.LatLng(%s,%s);'%(j,point['lat'][j],point['lon'][j])
print """
var ptsPath=new google.maps.Polyline({
    path:pts,
    strokeColor:"red",
    strokeOpacity:0.8,
    strokeWeight:4
    });
ptsPath.setMap(map);
}
google.maps.event.addDomListener(window, 'load', initialize);
</script>
"""

print '''
<script >
function show() {
    var w = window.open("","_blank","width=500,height=400,scrollbars=1");
    w.document.open();
    w.document.write("Number, Lattitude, Longitude, Time(EDT)</br>");
'''

for j in xrange(len(point['lon'])):
    showtime = (sloctim+timedelta(hours=j)).strftime('%D %H:%M')
    if track_way=='backward':
        showtime = (eloctim-timedelta(hours=j)).strftime('%D %H:%M')
    print "w.document.write('%d-%02d, %f, %f, %s</br>');"%(1,j+1,point['lat'][j],point['lon'][j],showtime)
    
print """
w.document.close();
}
</script>

<script>
function legendinstruction() {
    var w = window.open("","_blank","width=300,height=230,scrollbars=1");
    w.document.write("<p><img src='http://127.0.0.1:8000/startmarker.png' /> : Start-point.</p>")
    w.document.write("<p><img src='http://127.0.0.1:8000/endmarker.png' /> : End-point.</p>")
    w.document.write("</p><font color='blue' >Blue line</font> is drifter trajectory.</p>");
    w.document.write("</p><font color='red' >Red line</font> is forecast trajectory.</p>");
    }
</script>

<div id="header">
<h2>Model: %s, %s track %s days.</h2>
</div>

<div id="googleMap" ></div>

<div id="footer" ><a href="http://127.0.0.1:8000/index2.html"><input type="button" value="Track again" /></a>
<input type="button" value="Show lat/lon" onclick="show()" />
<input type="button" value="Legend" onclick="legendinstruction()" />
</div>
"""%(MODEL,track_way,tds)

print '</body></html>'

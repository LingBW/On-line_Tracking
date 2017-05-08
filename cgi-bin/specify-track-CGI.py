#!/usr/bin/env python
import cgitb
cgitb.enable()
import cgi,sys
#import pytz
#import pandas as pd
from dateutil.parser import parse
from datetime import datetime, timedelta
from track_functions import get_fvcom,get_roms

print """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
<head>
<script src="http://maps.googleapis.com/maps/api/js" type="text/javascript"></script>
<title>Specify Track</title>
"""

form = cgi.FieldStorage()
lons = form.getvalue('alon')
lats = form.getvalue('alat')
tds = form.getvalue('tds')
track_way = form.getvalue('tway') # Two option: backward, froward; 
stt = form.getvalue('stt')
dep = form.getvalue('depth')
dset = form.getvalue('dset')

st_lon=float(lons); st_lat=float(lats)

utcti = datetime.utcnow(); utct = utcti.strftime('%H')
locti = datetime.now(); loct = locti.strftime('%H')
ditnu = int(utct)-int(loct) # the deference between UTC and local time 
if ditnu < 0:
    ditnu = int(utct)+24-int(loct)
try:
    if not stt:
        start_time = datetime.utcnow()#edtl = (edt[0],edt[1],edt[2],edt[3],edt[4])
    else:
        stime = parse(stt) #datetime(2015,2,10,12,0,0,0,pytz.UTC ,dayfirst=True    
        start_time = stime+timedelta(hours=ditnu);#print start_time
except:
    print '<h2>"Start time: Error! Time format error."</h2>' 
    print '</head></html>'
    sys.exit() 
#stp_num = len(st_lat)   
MODEL = dset      # 'FVCOM', 'ROMS'
GRID = ['massbay','GOM3','30yr']     
forecast_days = float(tds)      #MODEL track time(days) 
depth = abs(float(dep))   
end_time = start_time + timedelta(forecast_days)
colors = ['magenta','olive','blue','orange','green','red','yellow','black','purple']
if track_way=='backward':
    end_time = start_time 
    start_time = end_time - timedelta(forecast_days)  #'''

# pass to html to show start-time and end-time of marker in ET. 
sloctim = start_time-timedelta(hours=ditnu)
eloctim = end_time-timedelta(hours=ditnu)
stim = sloctim.strftime('%D %H:%M')
etim = eloctim.strftime('%D %H:%M')
#stim = (start_time-timedelta(hours=ditnu)).strftime('%D %H:%M')
#etim = (end_time-timedelta(hours=ditnu)).strftime('%D %H:%M')

#curtime = pytz.utc.localize(datetime.now().replace(hour=0,minute=0,second=0,microsecond=0))
#mst = curtime - timedelta(3); met = curtime + timedelta(3)
#mstt = mst.strftime('%D %H:%M'); mett = met.strftime('%D %H:%M')
#lon_set = [[]]*stp_num; lat_set = [[]]*stp_num
if MODEL in GRID:
    get_obj = get_fvcom(MODEL)
    
    if (MODEL == 'massbay' or MODEL == 'GOM3') and depth < 2 :
        bpoints = get_obj.sf_get_data(start_time,end_time)
        
        point = get_obj.sf_get_track(st_lon,st_lat,track_way)
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
            
            print '<h2>Time: Error! Model(FVCOM:%s) only works between %s with %s. <br>Tips: check step 5 or 6.</h2>'%(MODEL,mstt,mett)
            print '</head></html>'
            sys.exit()
       
        point = get_obj.get_track(st_lon,st_lat,depth,track_way)
        if len(point['lon'])==1: 
            print '<script type="text/javascript">window.alert("This point on the land.")</script>'        
        
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
        point = get_obj.sf_get_track(st_lon,st_lat,track_way)
    else:
        point = get_obj.get_track(st_lon,st_lat,depth,track_way)
    if len(point['lon'])==1: 
        print '<script type="text/javascript">window.alert("The point on the land.")</script>'      
        
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

<body onunload="GUnload()">

<script type="text/javascript">
var myCenter=new google.maps.LatLng(%s,%s);
var stim = "%s";
var etim = "%s";
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
	        icon: 'http://127.0.0.1:8000/endmarker.png'
        	});     
        marker.setMap(map);   
        var infowindow = new google.maps.InfoWindow({
            content: 'Latitude: ' + parseFloat(pmarker.lat()).toFixed(6) +'<br>Longitude: ' + parseFloat(pmarker.lng()).toFixed(6)+'<br>EndTime: '+etim+' EDT'});
          
        google.maps.event.addListener(marker, 'click', function() {
            map.setZoom(10);
            map.setCenter(marker.getPosition());
            infowindow.open(map,marker);
        	});
    	}"""%(st_lat,st_lon,stim,etim)

# plot markers
print 'startmarker(%s,%s)'%(st_lat,st_lon)
print 'endmarker(%s,%s)'%(point['lat'][-1],point['lon'][-1])
# plot forecast line
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
    w.document.write("Number, Lattitude, Longitude, Time(EDT)</br>");
'''


for j in xrange(len(point['lon'])):
    showtime = (sloctim+timedelta(hours=j)).strftime('%D %H:%M')
    if track_way=='backward':
        showtime = (eloctim-timedelta(hours=j)).strftime('%D %H:%M')
    print "w.document.write('%d-%02d, %f, %f, %s</br>');"%(1,j+1,point['lat'][j],point['lon'][j],showtime)
        
print """
}
</script>

<script >
function legendinstruction() {
    var w = window.open("","_blank","width=300,height=200,scrollbars=1");
    w.document.write("<p><img src='http://127.0.0.1:8000/startmarker.png' /> : Start-point.</p>")
    w.document.write("<p><img src='http://127.0.0.1:8000/endmarker.png' /> : End-point.</p>")
    w.document.write("</p><font color='red' >Red line</font> is forecast trajectory.</p>");
    }
</script>
   
<div id="header">
<h2>Model:%s, %s track %s days.</h2>
</div>

<div id="googleMap" ></div>

<div id="footer" ><a href="http://127.0.0.1:8000/index1.html"><input type="button" value="Track again" /></a>
<input type="button" value="Show lat/lon" onclick="show()" />
<input type="button" value="Legend" onclick="legendinstruction()" />
</div>
"""%(MODEL,track_way,tds)

print '</script></body></html>'

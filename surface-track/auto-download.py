# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:29:54 2016

@author: bling
"""

import time
from datetime import datetime

print 'Automatic run at 19:00 PM(local) everyday.'

# the deference between UTC and local time
utcti = datetime.utcnow(); utct = utcti.strftime('%H')
locti = datetime.now(); loct = locti.strftime('%H')
ditnu = int(utct)-int(loct)  
if ditnu < 0:
    ditnu = int(utct)+24-int(loct)

q = 0; a = 0
while q<2880:# 24*60
    q+=1
    nowtime = datetime.now()
    etime = nowtime.strftime('%m/%d/%Y %H')
    
    if int(etime[-2:]) == (24-ditnu):      
        while a<100:
            a+=1
            nowtime = datetime.now()
            etime = nowtime.strftime('%m/%d/%Y %H:%M')
            try:
                execfile('get-model-data.py')
            except:
                print "Error: Can't run 'Track.py'."
                raise
            finally:    
                endtime = datetime.now()
                sleepseconds = (endtime-nowtime).seconds
                sleps = 86400-sleepseconds #24*3600=86400 one day
                print "%d days, at %s"%(a,etime)
                print 'Waiting to run next day...'
                time.sleep(sleps) 
    time.sleep(30) # one minute
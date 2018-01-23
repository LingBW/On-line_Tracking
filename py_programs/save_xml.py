# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:25:21 2017

@author: bling
"""

import urllib2
url = 'https://www.nefsc.noaa.gov/drifter/drift_X.xml'
s = urllib2.urlopen(url)
contents = s.read()
file = open("drift_X.xml", 'w')
file.write(contents)
file.close()

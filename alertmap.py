#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import ConfigParser
import argparse
from xml.dom.minidom import parse

#local imports
from neicio.cmdoutput import getCommandOutput
from neicio.shake import ShakeGrid
from travel.travel import TravelTimeCalculator

#third party imports
from obspy.core.util import locations2degrees
from obspy.fdsn import Client
import numpy as np
from neicio import GMTGrid

SHAKEHOME = '/home/shake/ShakeMap'

EVENT_DEFAULT = '''<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="[EVENTID]" lat="[LAT]" lon="[LON]" mag="[MAG]" year="[YEAR]" month="[MONTH]" day="[DAY]" hour="[HOUR]" minute="[MINUTE]" second="[SECOND]" timezone="GMT" depth="[DEPTH]" locstring="[LOCSTR]" created="1407055672" otime="1407054613" type="" network="us" />
'''

GRIND_DEFAULT = '''smVs30default : 686
use_gmpe_sc : true
x_grid_interval : [DX]
y_grid_interval : [DY]
lonspan : [LONSPAN]
latspan : [LATSPAN]
bad_station : 8016 9.9 19990101-
bad_station : 8010 9.9 19990101-
bad_station : 8022 9.9 19990101-
bad_station : 8034 9.9 19990101-
bad_station : 8040 9.9 19990101-

gmpe: [GMPE] 0.0 9.9 0 999

outlier_deviation_level : 3
outlier_max_mag         : 7.0

bias_norm         : l1
bias_max_range    : 120
bias_min_stations : 6
bias_max_mag      : 7.0
bias_max_bias     : 2.0
bias_min_bias     : -2.0
bias_log_amp      : true

direct_patch_size : 1000

fwdata_file		: forward.xml

topobin : <HOME>/bin/topo2grd <EVID> <BOUND> regime=active

pgm2mi: [GMICE]
mi2pgm: [GMICE]
'''

def getEventText(eventfile,lat,lon):
    root = parse(eventfile)
    eq = root.getElementsByTagName('earthquake')[0]
    eventdict = {}
    eventdict['lat'] = lat
    eventdict['lon'] = lon
    eventdict['eventid'] = eq.getAttribute('id')
    eventdict['mag'] = float(eq.getAttribute('mag'))
    eventdict['year'] = int(eq.getAttribute('year'))
    eventdict['month'] = int(eq.getAttribute('month'))
    eventdict['day'] = int(eq.getAttribute('day'))
    eventdict['hour'] = int(eq.getAttribute('hour'))
    eventdict['minute'] = int(eq.getAttribute('minute'))
    eventdict['second'] = int(round(float(eq.getAttribute('second'))))
    eventdict['depth'] = float(eq.getAttribute('depth'))
    eventdict['locstr'] = eq.getAttribute('locstring')
    root.unlink()
    eventtext = EVENT_DEFAULT
    for key,value in eventdict.iteritems():
        eventtext = eventtext.replace('['+key.upper()+']',str(value))
    return eventtext
    
def getSlowestStation(lat,lon,calc):
    client = Client("IRIS")
    inventory = client.get_stations(latitude=lat, longitude=lon,maxradius=1.5)
    lats = []
    lons = []
    for network in inventory.networks:
        for station in network.stations:
            lats.append(station.latitude)
            lons.append(station.longitude)
    lats = np.array(lats)
    lons = np.array(lons)
    distances = []
    times = []
    for i in range(0,len(lats)):
        slat = lats[i]
        slon = lons[i]
        distance = locations2degrees(lat,lon,slat,slon)
        distances.append(distance)
        ptime,stime = calc.getTravelTimes(distance)
        times.append(ptime)
    times = np.array(times)
    distances = np.array(distances)
    sortidx = np.argsort(distances)
    distances = distances[sortidx]
    times = times[sortidx]
    lats = lats[sortidx]
    lons = lons[sortidx]
    distances = distances[0:4]
    times = times[0:4]
    lats = lats[0:4]
    lons = lons[0:4]
    idx = times.argmax()
    return (lats[idx],lons[idx],times[idx])

def writeGrind(config,datadir):
    gmpe = config.get('MAP','gmpe')
    gmice = config.get('MAP','gmice')
    ipe = config.get('MAP','ipe')
    xmin = float(config.get('MAP','xmin'))
    xmax = float(config.get('MAP','xmax'))
    ymin = float(config.get('MAP','ymin'))
    ymax = float(config.get('MAP','ymax'))
    dx = float(config.get('MAP','dx'))
    dy = float(config.get('MAP','dy'))
    lonspan = xmax - xmin
    latspan = ymax - ymin
    
    grindfile = os.path.join(datadir,'config','grind.conf')
    grindstr = GRIND_DEFAULT
    grindstr = grindstr.replace('[GMPE]',gmpe)
    grindstr = grindstr.replace('[GMICE]',gmice)
    grindstr = grindstr.replace('[IPE]',ipe)
    grindstr = grindstr.replace('[DX]',str(dx))
    grindstr = grindstr.replace('[DY]',str(dy))
    grindstr = grindstr.replace('[LONSPAN]',str(lonspan))
    grindstr = grindstr.replace('[LATSPAN]',str(latspan))

    f = open(grindfile,'wt')
    f.write(grindstr)
    f.close()

def main(args):
    datadir = os.path.join(SHAKEHOME,'data',args.event)
    if not os.path.isdir(datadir):
        print 'Cannot find event %s on the system' % args.event
        sys.exit(1)
    #now look for config file in top-level folder
    configfile = os.path.join(datadir,'alert.conf')
    if not os.path.isfile(configfile):
        print 'Cannot find alert config file for %s in the data directory' % args.event
        sys.exit(1)
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))

    #get the array of epicenters
    lats = [float(p) for p in config.get('FAULT','lats').split()]
    lons = [float(p) for p in config.get('FAULT','lons').split()]

    #write out a new grind.conf file
    writeGrind(config,datadir)

    #instantiate our p/s travel time calculator
    calc = TravelTimeCalculator()

    #where is the grind binary?
    grindbin = os.path.join(SHAKEHOME,'bin','grind')

    #specify the event.xml file, get the depth of the event
    eventfile = os.path.join(datadir,'input','event.xml')
    root = parse(eventfile)
    eq = root.getElementsByTagName('earthquake')[0]
    depth = float(eq.getAttribute('depth'))
    root.unlink()
       
    #loop over all the event realizations
    for i in range(0,len(lats)):
        lat = lats[i]
        lon = lons[i]
        if i == 0:
            lonoff = 0
            latoff = 0
        else:
            lonoff = -1* (lons[i] - lons[i-1])
            latoff = lats[i] - lats[i-1]
        #modify the event.xml file to have the new lat/lon epicenter
        sourcetext = getEventText(eventfile,lat,lon)
        f = open(eventfile,'wt')
        f.write(sourcetext)
        f.close()

        grindcmd = '%s -latoff %f -lonoff %f -event %s' % (grindbin,latoff,lonoff,args.event)
        res,stdout,stderr = getCommandOutput(grindcmd)
        if not res:
            print 'Grind command failed: "%s", "%s"' % (stdout,stderr)
            sys.exit(1)

        stationlat,stationlon,ptime = getSlowestStation(lat,lon,calc)
            
        #Get the grid.xml output, do some time calculations
        gridfile = os.path.join(datadir,'output','grid.xml')
        mmigrid = ShakeGrid(gridfile,variable='MMI')
        m,n = mmigrid.griddata.shape
        timegrid = np.zeros((m,n))
        for row in range(0,m):
            for col in range(0,n):
                mmilat,mmilon = mmigrid.getLatLon(row,col)
                distance = locations2degrees(stationlat,stationlon,mmilat,mmilon)
                ptime,stime = calc.getTravelTimes(distance)
                timegrid[row,col] = stime - ptime
        timegmt = GMTGrid()
        timegmt.griddata = timegrid
        timegmt.geodict = mmigrid.geodict
        timefile = os.path.join(datadir,'output','timegrid%03i' % i)
        timegmt.save(timefile)
        
if __name__ == '__main__':
    desc = '''This script does the following:
    1) Find ShakeMap event folder from input ID.
    2) Parse alert.conf file in event data folder (i.e. /home/shake/ShakeMap/data/eventID/alert.conf)
    [FAULTS]
    lats = 32.1 32.2 32.3
    lons = -118.1 -118.2 -118.3
    [MAP]
    gmpe = AkkarBommer07
    gmice = WGRW11
    ipe = AW07_CEUS
    xmin = -119.0
    xmax = -117.0
    ymin = 31.0
    ymax = 33.0
    dx = 0.01
    dy = 0.01
    3) Replace grind.conf file in event data folder with one created using information supplied above
    4) For each epicenter in (lats,lons):
      a) Write new event.xml file
      b) Run ShakeMap grind program
      c) Find 4 nearest stations to epicenter, calculate P arrival times for each, return the slowest.
      d) Calculate S arrival times for each cell in ShakeMap, subtract SlowP+8.5 from each.
      e) Save grid of S arrival times
    5) TODO SCIENCE
    '''
    parser = argparse.ArgumentParser(description=desc,formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument('event', help='Select the event ID')
    pargs = parser.parse_args()
    main(pargs)
    

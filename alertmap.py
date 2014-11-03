#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import ConfigParser
import argparse
from xml.dom.minidom import parse
import csv
from operator import itemgetter

#local imports
from neicio.cmdoutput import getCommandOutput
from neicio.shake import ShakeGrid
from neicio.esri import EsriGrid
from travel.travel import TravelTimeCalculator,saveTimeGrid,readTimeGrid

#third party imports
from obspy.core.util import locations2degrees
from obspy.fdsn import Client
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from neicutil.text import roundToNearest,ceilToNearest,floorToNearest
from neicutil.colormap import GMTColormap
from matplotlib.colors import ListedColormap,LinearSegmentedColormap,Normalize,BoundaryNorm

WATER_COLOR = [.47,.60,.81]
#CITY_COLOR = '#33FF00'
CITY_COLOR = '#FF0000'
NMAPCITIES = 10

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

def getCityList(xmin,xmax,ymin,ymax,cityfile):
    cities = []
    f = open(cityfile,'rt')
    for line in f.readlines():
        city = {}
        parts = line.split('\t')
        city['name'] = parts[2].strip()
        city['lat'] = float(parts[4].strip())
        city['lon'] = float(parts[5].strip())
        city['pop'] = int(parts[14].strip())
        if not city['name']:
            #print 'Found a city with no name'
            continue
        if not all(ord(c) < 128 for c in city['name']):
            continue
        if city['lat'] >= ymin and city['lat'] <= ymax and city['lon'] >= xmin and city['lon'] <= xmax:
            cities.append(city)
    f.close()
    cities.sort(key=itemgetter('pop'),reverse=True)
    return cities

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
    codes = []
    for network in inventory.networks:
        for station in network.stations:
            lats.append(station.latitude)
            lons.append(station.longitude)
            codes.append(station.code)
    lats = np.array(lats)
    lons = np.array(lons)
    codes = np.array(codes)
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
    codes = codes[sortidx]
    distances = distances[0:4]
    times = times[0:4]
    lats = lats[0:4]
    lons = lons[0:4]
    codes = codes[0:4]
    idx = times.argmax()
    sdict = {'lat':lats[idx],'lon':lons[idx],'time':times[idx],'code':codes[idx]}
    return sdict

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

def detectShakeHome():
    isFound = False
    for shakehome in SHAKEHOMELIST:
        grind = os.path.join(shakehome,'bin','grind')
        if os.path.isfile(grind):
            isFound = True
            break
    if isFound:
        return shakehome
    else:
        return None

def getMapLines(dmin,dmax):
    NLINES = 4
    drange = dmax-dmin
    if drange > 4:
        near = 1
    else:
        if drange >= 0.5:
            near = 0.25
        else:
            near = 0.125
    inc = roundToNearest(drange/NLINES,near)
    if inc == 0:
        near = pow(10,round(log10(drange))) #make the increment the closest power of 10
        inc = ceilToNearest(drange/NLINES,near)
        newdmin = floorToNearest(dmin,near)
        newdmax = ceilToNearest(dmax,near)
    else:
        newdmin = ceilToNearest(dmin,near)
        newdmax = floorToNearest(dmax,near)
    darray = np.arange(newdmin,newdmax,inc)
    return darray

def getLatLonGrids(shake):
    dims = shake.getData().shape
    nrows = dims[0]
    ncols = dims[1]
    geodict = shake.getGeoDict()
    xmin = geodict['xmin']
    xmax = geodict['xmax']
    ymin = geodict['ymin']
    ymax = geodict['ymax']
    xdim = geodict['xdim']
    ydim = geodict['ydim']

    if xmax < xmin:
        xmax = xmax + 360

    lonrow = np.arange(xmin,xmax,xdim)
    latcol = np.arange(ymin,ymax,ydim)

    if len(lonrow) < ncols:
        lonrow = concatenate((lonrow,[xmax]))
    elif len(lonrow) > ncols:
        lonrow = lonrow[0:-1]
    if len(latcol) < nrows:
        latcol = concatenate((latcol,[ymax]))
    elif len(latcol) > nrows:
        latcol = latcol[0:-1]

    longrid = np.zeros((len(latcol),len(lonrow)),dtype=float)
    latgrid = np.zeros((len(latcol),len(lonrow)),dtype=float)
    for i in range(0,len(lonrow)):
        longrid[:,i] = lonrow[i]
    for i in range(0,len(latcol)):
        latgrid[i,:] = latcol[i]
    return longrid,latgrid

def makeMap(statgrid,timegrid,metadata,method,datadir,popfile,popcolormap,stationdict,citylist,lats,lons):
    figwidth = 8.0
    bounds = timegrid.getRange()
    bounds = list(bounds)
    if bounds[1] < 0 and bounds[0] > bounds[1]:
        bounds[1] = bounds[1] + 360

    clat = bounds[2] + (bounds[3] - bounds[2])/2
    clon = bounds[0] + (bounds[1] - bounds[0])/2
    dx = (bounds[1] - bounds[0])*111191 * np.cos(np.degrees(clat))
    dy = (bounds[3] - bounds[2])*111191
    aspect = np.abs(dy/dx)
    figheight = aspect * figwidth
    fig = plt.figure(figsize=(figwidth,figheight),edgecolor='g',facecolor='g')
    ax1 = fig.add_axes([0,0,1.0,1.0])
    m = Basemap(llcrnrlon=bounds[0],llcrnrlat=bounds[2],
                urcrnrlon=bounds[1],urcrnrlat=bounds[3],
                resolution='h',projection='merc',lat_ts=clat)
    
    #get the population grid
    popgrid = EsriGrid(popfile)
    popgrid.load(bounds=bounds)
    popdata = np.flipud(popgrid.griddata)

    cmap = GMTColormap(popcolormap)
    clist = cmap.getColorList()
    boundaries = cmap.getZValues()
    palette = ListedColormap(clist,'my_colormap')
 
    i = np.where(np.isnan(popdata))
    popdata[i] = -1
    popdatam = np.ma.masked_values(popdata, -1)
    palette.set_bad(WATER_COLOR,1.0)
    
    ncolors = len(boundaries)
    am = m.imshow(popdatam,cmap=palette,norm=BoundaryNorm(boundaries,ncolors),interpolation='nearest')

    statgrid = np.flipud(statgrid)
    (lons,lats) = getLatLonGrids(timegrid)
    (x,y) = m(lons,lats)
    clevels = np.arange(5,45,5)
    cs = m.contour(x,y,statgrid,clevels)
    #plt.clabel(cs, inline=1, fontsize=10)
    proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_color()[0]) for pc in cs.collections]
    labels = [str(c)+' sec' for c in clevels]

    sx,sy = m(stationdict['lon'],stationdict['lat'])
    m.plot(sx,sy,'rD')
    plt.text(sx,sy,stationdict['code'])

    #plot the various epicenters
    for elat,elon in zip(lats,lons):
        ex,ey = m(elon,elat)
        m.plot(ex,ey,'bx')

    #plot the cities
    for i in range(0,NMAPCITIES):
        city = citylist[i]
        cx,cy = m(city['lon'],city['lat'])
        m.plot(cx,cy,'.',color=CITY_COLOR)
        plt.text(cx,cy,city['name'],color=CITY_COLOR)
    
    m.drawrivers(color=WATER_COLOR)
    m.drawcountries(color='k',linewidth=2.0)
    mer = getMapLines(bounds[0],bounds[1])
    par = getMapLines(bounds[2],bounds[3])
    
    xmap_range = m.xmax-m.xmin
    ymap_range = m.ymax-m.ymin
    xoff = -0.09*(xmap_range)
    yoff = -0.04*(ymap_range)

    m.drawmeridians(mer,labels=[0,0,1,0],fontsize=8,
                             linewidth=0.5,color='black',yoffset=yoff,xoffset=xoff,dashes=[1,0.01])
    m.drawparallels(par,labels=[0,1,0,0],fontsize=8,
                             linewidth=0.5,color='black',yoffset=yoff,xoffset=xoff,dashes=[1,0.01])
    m.drawmapboundary(color='k',linewidth=2.0)

    plt.legend(proxy,labels)
    
    outfile = os.path.join(datadir,method+'.pdf')
    plt.savefig(outfile)

def getGlobalConfig():
    configfile = os.path.join(os.path.expanduser('~'),'.alertmap','config.ini')
    if not os.path.isfile(configfile):
        raise Exception,'Could not find global config file "%s".' % configfile
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))
    gdict = {}
    gdict['shakehome'] = config.get('GLOBAL','shakehome')
    gdict['popfile'] = config.get('GLOBAL','popfile')
    gdict['popcolormap'] = config.get('GLOBAL','popcolormap')
    gdict['cityfile'] = config.get('GLOBAL','cityfile')
    return gdict
    
def main(args):
    globaldict = getGlobalConfig()
    shakehome = globaldict['shakehome']
    popfile = globaldict['popfile']
    if shakehome is None:
        print 'Cannot find ShakeMap home folder on this system.'
        sys.exit(1)
    datadir = os.path.join(shakehome,'data',args.event)
    if not os.path.isdir(datadir):
        print 'Cannot find event %s on the system' % args.event
        sys.exit(1)

    #Make sure the timeoutput folder is there (can't put our time grids in output - that gets
    #wiped out every time shakemap runs
    outfolder = os.path.join(datadir,'timeoutput')
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
        
    #now look for config file in top-level folder
    configfile = os.path.join(datadir,'alert.conf')
    if not os.path.isfile(configfile):
        print 'Cannot find alert config file for %s in the data directory' % args.event
        sys.exit(1)
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))

    #get the bounds of the map so we can find cities
    xmin = float(config.get('MAP','xmin'))
    xmax = float(config.get('MAP','xmax'))
    ymin = float(config.get('MAP','ymin'))
    ymax = float(config.get('MAP','ymax'))
    
    citylist = getCityList(xmin,xmax,ymin,ymax,globaldict['cityfile'])
    
    #Get the MMI threshold below which alert times will NOT be saved
    mmithresh = float(config.get('MAP','mmithresh'))

    #get the array of epicenters
    lats = [float(p) for p in config.get('FAULT','lats').split()]
    lons = [float(p) for p in config.get('FAULT','lons').split()]

    #write out a new grind.conf file
    writeGrind(config,datadir)

    #instantiate our p/s travel time calculator
    calc = TravelTimeCalculator()

    #where is the grind binary?
    grindbin = os.path.join(shakehome,'bin','grind')

    #specify the event.xml file, get the depth of the event
    eventfile = os.path.join(datadir,'input','event.xml')
    root = parse(eventfile)
    eq = root.getElementsByTagName('earthquake')[0]
    depth = float(eq.getAttribute('depth'))
    root.unlink()
       
    #loop over all the event realizations
    timefiles = []
    for i in range(0,len(lats)):
        print 'Calculating arrival times for scenario %i of %i' % (i+1,len(lats))
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

        sdict = getSlowestStation(lat,lon,calc)
        ptime = sdict['time']
        stationlat = sdict['lat']
        stationlon = sdict['lon']
        
        grindcmd = '%s -latoff %f -lonoff %f -event %s' % (grindbin,latoff,lonoff,args.event)
        res,stdout,stderr = getCommandOutput(grindcmd)
        if not res:
            print 'Grind command failed: "%s", "%s"' % (stdout,stderr)
            sys.exit(1)
            
        #Get the grid.xml output, do some time calculations
        gridfile = os.path.join(datadir,'output','grid.xml')
        mmigrid = ShakeGrid(gridfile,variable='MMI')
        m,n = mmigrid.griddata.shape
        timegrid = np.zeros((m,n),dtype=np.float32)
        
        for row in range(0,m):
            for col in range(0,n):
                mmilat,mmilon = mmigrid.getLatLon(row,col)
                if mmigrid.griddata[row,col] < mmithresh:
                    timegrid[row,col] = np.nan
                distance = locations2degrees(stationlat,stationlon,mmilat,mmilon)
                tmp,stime = calc.getTravelTimes(distance)
                timegrid[row,col] = stime - ptime
        
        timefile = os.path.join(outfolder,'timegrid%03i.flt' % (i+1))
        timefiles.append(timefile)
        metadict = {'epilat':lat,'epilon':lon,'eventid':args.event}
        saveTimeGrid(timefile,timegrid,mmigrid.geodict,metadict)
    timestack = np.zeros((m,n,len(lats)),dtype=np.float32)
    for i in range(0,len(timefiles)):
        timefile = timefiles[i]
        timegrid,metadata = readTimeGrid(timefile)
        timestack[:,:,i] = timegrid.griddata

    methods = config.get('MAP','output').split(',')
    for method in methods:
        if method == 'median':
            statgrid = np.median(timestack,axis=2)
        if method == 'mean':
            statgrid = np.mean(timestack,axis=2)
        if method == 'min':
            statgrid = np.min(timestack,axis=2)
        if method == 'max':
            statgrid = np.max(timestack,axis=2)    
        makeMap(statgrid,timegrid,metadata,method,outfolder,popfile,globaldict['popcolormap'],sdict,citylist,lats,lons)
        
if __name__ == '__main__':
    desc = '''This script does the following:
    1) Find ShakeMap event folder from input ID.
    2) Parse alert.conf file in event data folder (i.e. /home/shake/ShakeMap/data/eventID/alert.conf)
    [FAULT]
    lats = 32.1 32.2 32.3
    lons = -118.1 -118.2 -118.3
    [MAP]
    #output can be comma separated list with any of min,mean,median,max
    output = median 
    #mmithresh - MMI value below which alert times will NOT be saved
    mmithresh = 6.0 
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
    

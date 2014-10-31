#!/usr/bin/env python

#stdlib imports
import os.path
import glob

#third party
import numpy as np
from scipy import interpolate
import rasterio
from affine import Affine

#local
from neicio.gmt import GMTGrid

class TravelTimeCalculator(object):
    def __init__(self):
        homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
        ttimefile = os.path.join(homedir,'ttimes.csv')
        data = np.loadtxt(ttimefile,delimiter=',',skiprows=1)
        self.fp = interpolate.interp1d(data[:,0],data[:,1])
        self.fs = interpolate.interp1d(data[:,0],data[:,2])

    def getTravelTimes(self,distance):
        ptime = self.fp(distance)
        stime = self.fs(distance)
        return (ptime,stime)

def readTimeHeader(timehdr):
    lines = open(timehdr,'rt').readlines()
    timedict = {}
    for line in lines:
        parts = line.split()
        key = parts[0]
        value = ' '.join(parts[1:])
        try:
            value = float(value)
        except:
            try:
                value = int(value)
            except:
                pass
        timedict[key] = value
    return timedict

def writeTimeHeader(timedict,timehdrfile):
    f = open(timehdrfile,'wt')
    for key,value in timedict.iteritems():
        key = key.replace(' ','_')
        if isinstance(value,str) and value.find(' ') > -1:
            value = '"'+value+'"'
        if isinstance(value,float):
            fmt = '%s %f\n'
        elif isinstance(value,int):
            fmt = '%s %i\n'
        else:
            fmt = '%s %s\n'
        f.write(fmt % (key,value))
    f.close()
    
def saveTimeGrid(timefile,timegrid,geodict,metadict):
    fbase,fname = os.path.split(timefile)
    timebase,timext = os.path.splitext(fname)
    timetiff = os.path.join(fbase,timebase+'.tiff')
    aff = Affine(geodict['xdim'],0.0,geodict['xmin'],0.0,geodict['xdim'],geodict['ymax'])
    crs = {'no_defs': True, 'ellps': 'WGS84', 'datum': 'WGS84', 'proj': 'longlat'}
    m,n = timegrid.shape
    timeio = rasterio.open(timetiff,mode='w',driver='GTiff',
                           dtype=rasterio.float32,transform=aff,
                           crs=crs,count=1,height=m,width=n)
    timeio.write_band(1,timegrid.astype(np.float32))
    timeio.close()
    with rasterio.drivers():
        rasterio.copy(timetiff,timefile,driver='EHdr')

    os.remove(timetiff)
    timehdr = os.path.join(fbase,timebase+'.hdr')
    timedict = readTimeHeader(timehdr)
    for key,value in metadict.iteritems():
        timedict[key] = value
    writeTimeHeader(timedict,timehdr)
    return True

def readTimeGrid(timefile):
    stkeys = ['TOTALROWBYTES','NBITS','LAYOUT','YDIM','NCOLS',
              'BANDROWBYTES','PIXELTYPE','XDIM','NROWS',
              'NBANDS','ULXMAP','ULYMAP','BYTEORDER']
    src = rasterio.open(timefile,'r',driver='EHdr')
    timedata, = src.read()
    m,n = timedata.shape
    aff = src.affine
    xdim = aff[0]
    xmin = aff[2]
    ydim = -aff[4]
    ymax = aff[5]
    src.close()
    timegrid = GMTGrid()
    timegrid.griddata = timedata
    timegrid.geodict = {'nrows':m,'ncols':n,'nbands':1,'bandnames':['Alert Time'],
                    'xmin':xmin,'xmax':xmin+n*xdim,'ymin':ymax-m*ydim,'ymax':ymax}
    
    timepath,timefile = os.path.split(timefile)
    timebase,timext = os.path.splitext(timefile)
    timehdr = os.path.join(timepath,timebase+'.hdr')
    timedict = readTimeHeader(timehdr)
    for key in stkeys:
        timedict.pop(key)
    for key,value in timedict.iteritems():
        if isinstance(value,str):
            timedict[key] = value.replace('"','')
    return (timegrid,timedict)

if __name__ == '__main__':
    timegrid = np.random.rand(8,13)
    geodict = {'xmin':-118.1,'ymax':32.1,'xdim':0.08,'ydim':0.08}
    timefile = 'temptime.flt'
    saveTimeGrid(timefile,timegrid,geodict,{'foo':'bar','my friend':'barney rubble'})
    timegrid2,metadata = readTimeGrid(timefile)
    timebase,timeext = os.path.splitext(timefile)
    timefiles = glob.glob(timebase+'*')
    for timefile in timefiles:
        os.remove(timefile)
    

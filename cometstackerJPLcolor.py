import numpy as np
import math
import ephem
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io.fits import getheader
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import FK5
from PIL import Image
import cv2
import os
import sys
import fnmatch

def mouse_click(event, x, y, flags, params):
    global xcoord, ycoord
    if event == cv2.EVENT_LBUTTONDOWN:
        xcoord, ycoord = x,y
        cv2.destroyAllWindows()

if len(sys.argv) != 4:
    print('Proper use: python cometstackerJPLcolor.py elements.txt latitude longitude')
    exit()
pattern = '*.new*'
filelist = os.listdir('.')
#minval = int(float(sys.argv[1])*(255*255))
#maxval = float(255/int(sys.argv[2]))
#mag = float(sys.argv[3])
observer = ephem.Observer()
observer.lat = str(str(sys.argv[2]))
observer.lon = str(str(sys.argv[3]))
observer.elevation = 0
observer.pressure = 1013
with open(str(sys.argv[1])) as f:
    lines = [line.rstrip('\n') for line in f]
for idx, line in enumerate(lines):
    if "$$SOE" in line:
        targetname = str('target')
        line1 = lines[idx+1]
        line2 = lines[idx+2]
        line3 = lines[idx+3]
        line4 = lines[idx+4]
        line5 = lines[idx+5] 
        linesplit1 = line1.split(' ')
        dateline = float(linesplit1[0]) - 2415020
        observer.date = dateline
        datesplit = str(observer.date).split('/')
        year = datesplit[0]
        month = datesplit[1]
        day = float(datesplit[2].split(' ')[0])
        fractionday = str('0.'+str(linesplit1[0].split('.')[1]))
        fractionday = float(fractionday) - 0.5
        if fractionday < 0:
            fractionday = 1 + fractionday
        day = day + fractionday
        xephemdate = str(str(month) + '/' + str(day) + '/' + str(year))
        ec = float(line2[4:26])
        qr = float(line2[30:52])
        inc = float(line2[56:78])

        om = float(line3[4:26])
        w = float(line3[30:52])
        tp = float(line3[56:78])

        n = float(line4[4:26])
        ma = float(line4[30:52])
        ta = float(line4[56:78])

        a = float(line5[4:26])
        ad = float(line5[30:52])
        pr = float(line5[56:78])
        if ec<1:
            xephem = str(targetname + ',' + 'e' + ',' + str(inc) + ',' + str(om) + ',' + str(w) + ',' + str(a) + ',' + str(n) + ',' + str(ec) + ',' + str(ma) + ',' + xephemdate + ',' + '2000' + ',' + 'g  6.5,4.0')
        else:
            dateline = float(tp) - 2415020
            observer.date = dateline
            datesplit = str(observer.date).split('/')
            year = datesplit[0]
            month = datesplit[1]
            day = float(datesplit[2].split(' ')[0])
            fractionday = str('0.'+str(linesplit1[0].split('.')[1]))
            fractionday = float(fractionday) - 0.5
            if fractionday < 0:
                fractionday = 1 + fractionday
            day = day + fractionday
            xephemdate = str(str(month) + '/' + str(day) + '/' + str(year))
            xephem = str(targetname + ',' + 'h' + ',' + xephemdate + ',' + str(inc) + ',' + str(om) + ',' + str(w) + ',' + str(ec) + ',' + str(qr) + ',' + '2000' + ',' + 'g  6.5,4.0')
target = ephem.readdb(xephem)
firstimage = True
sx = 0
sy = 0
maxx = 0
maxy = 0
minx = 0
miny = 0

u = math.atan(0.996647*math.tan(math.radians(observer.lat)))
psintheta = 0.996647*math.sin(u)+(observer.elevation/6378140)*math.sin(math.radians(observer.lat))
pcostheta = math.cos(u)+(observer.elevation/6378140)*math.cos(math.radians(observer.lat))


for idx, entry in enumerate(filelist):
    if fnmatch.fnmatch(entry,pattern):
        skip = False
        hdu = fits.open(entry)[0]
        hdu.header['NAXIS'] = 2
        hdr = getheader(entry, 0)
        wcs = WCS(hdu)
        dateobs = hdr['DATE-OBS']
        idate, itime = dateobs.split('T')
        year, month, day = idate.split('-')
        hour, minute, seconds = itime.split(':')
        #second, millisecond = seconds.split('.')
        d = ephem.Date(str(year + '/' + month + '/' + day + ' ' + hour + ':' + minute + ':' + seconds))
        observer.date = d
        target.compute(observer)
        targetra = math.degrees(target.a_ra)
        targetdec = math.degrees(target.a_dec)
        targetdistance = target.earth_distance * 149598000 
        #topocentric parallax
        lst = math.degrees(observer.sidereal_time())
        
        hourangle = lst - targetra
        
        delta = math.atan( (pcostheta*math.sin(math.radians(lst)))/(targetdistance*math.cos(math.radians(targetdec))-pcostheta*math.cos(math.radians(hourangle))))
        hourangleprime = hourangle + delta
        coordinates = [str(str(targetra) + ' ' + str(targetdec))]
        coordinates2 = SkyCoord(coordinates, frame=FK5, unit="deg")
        pixels = skycoord_to_pixel(coordinates2,wcs,origin=0,mode='wcs')
        #print(pixels)
        if firstimage is True:
            fx, fy = pixels
            sx = fx[0] 
            sy = fy[0]
            firstimage = False
        else:
            cx, cy = pixels
            cx = cx[0]
            cy = cy[0]
            dx = round(sx - cx)
            dy = round(sy - cy)
            if dx > maxx:
                maxx = dx
            if dx < minx:
                minx = dx
            if dy > maxy:
                maxy = dy
            if dy < miny:
                miny = dy
print(minx, maxx, miny, maxy)
firstimage = True
            
for idx, entry in enumerate(filelist):
    if fnmatch.fnmatch(entry,pattern):
        skip = False
        hdu = fits.open(entry)[0]
        hdu.header['NAXIS'] = 2
        hdr = getheader(entry, 0)
        wcs = WCS(hdu)
        dateobs = hdr['DATE-OBS']
        idate, itime = dateobs.split('T')
        year, month, day = idate.split('-')
        hour, minute, seconds = itime.split(':')
        d = ephem.Date(str(year + '/' + month + '/' + day + ' ' + hour + ':' + minute + ':' + seconds))
        observer.date = d
        target.compute(observer)
        targetra = math.degrees(target.a_ra)
        targetdec = math.degrees(target.a_dec)
        
        #topocentric parallax
        lst = math.degrees(observer.sidereal_time())
        
        hourangle = lst - targetra
        
        delta = math.atan( (pcostheta*math.sin(math.radians(lst)))/(targetdistance*math.cos(math.radians(targetdec))-pcostheta*math.cos(math.radians(hourangle))))
        hourangleprime = hourangle + delta
        coordinates = [str(str(targetra) + ' ' + str(targetdec))]
        coordinates2 = SkyCoord(coordinates, frame=FK5, unit="deg")
        pixels = skycoord_to_pixel(coordinates2,wcs,origin=0,mode='wcs')
        #print(pixels)
        if firstimage is True:
            fx, fy = pixels
            sx = fx[0] 
            sy = fy[0]
            image_file = entry
            image_data = fits.getdata(image_file)
            imagenew = np.array(image_data,dtype = np.uint16)
            print(imagenew.shape)
            print(coordinates2)
            imagenew = np.pad(imagenew,((0,0),(int(abs(miny)),int(abs(maxy))),(int(abs(minx)),int(abs(maxx)))), mode='constant')[:,:,:]
            hdu = fits.PrimaryHDU(imagenew)
            hdul = fits.HDUList([hdu])
            hdul.writeto(str(str(entry.split('.')[0])+str(idx)+'_aligned.fits'))
        else:
            cx, cy = pixels
            cx = cx[0]
            cy = cy[0]
            dx = round(sx - cx)
            dy = round(sy - cy)
            print(coordinates2, dx, dy)
            image_file = entry
            image_data = fits.getdata(image_file)
            image1 = np.array(image_data,dtype = np.uint16)
            if dy < 0:
                imagenew = np.pad(image1,((0,0),(0,int(abs(maxy+dy))),(0,0)), mode='constant')[:,:]
                imagenew = np.pad(imagenew,((0,0),(int(abs(miny-dy)),0),(0,0)), mode='constant')[:,:]
            elif dy > 0:
                imagenew = np.pad(image1,((0,0),(int(abs(miny+dy)),0),(0,0)), mode='constant')[:,:]
                imagenew = np.pad(imagenew,((0,0),(0,int(abs(maxy-dy))),(0,0)), mode='constant')[:,:]
            else:
                imagenew = np.pad(image1,((0,0),(int(abs(miny)),int(abs(maxy))),(0,0)), mode='constant')[:,:]
            if dx < 0:
                imagenew = np.pad(imagenew,((0,0),(0,0),(0,int(abs(maxx+dx)))), mode='constant')[:,:]
                imagenew = np.pad(imagenew,((0,0),(0,0),(int(abs(minx-dx)),0)), mode='constant')[:,:]
            elif dx > 0:
                imagenew = np.pad(imagenew,((0,0),(0,0),(int(abs(minx+dx)),0)), mode='constant')[:,:]
                imagenew = np.pad(imagenew,((0,0),(0,0),(0,int(abs(maxx-dx)))), mode='constant')[:,:]
            else:
                imagenew = np.pad(imagenew,((0,0),(0,0),(int(abs(minx)),int(abs(maxx)))), mode='constant')[:,:]
            hdu = fits.PrimaryHDU(imagenew)
            hdul = fits.HDUList([hdu])
            hdul.writeto(str(str(entry.split('.')[0])+str(idx)+'_aligned.fits'))

        firstimage = False

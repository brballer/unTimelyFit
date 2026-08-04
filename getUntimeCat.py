#!/Library/Frameworks/Python.framework/Versions/3.12/bin/python3

'''''
Create a catalog of unTimely objects near the ra, dec location

Dec 10, 2025: Add resx,resy
Feb 2, 2026: comment out ntry loops
Feb 4, 2026: require flux > 0
'''''

from astropy.wcs import WCS
from astropy.time import Time
import fitsio
import warnings
import glob

warnings.filterwarnings('ignore',category=Warning)

def getUntimeCat(UntimelyCatalogPath,ra,dec,coaddID,nPixels):
    if UntimelyCatalogPath == None:
        print('getUntimeCat: UntimelyCatalogPath path not defined')
        exit()
    fitsFiles = UntimelyCatalogPath+coaddID+'*.cat.fits'
    files = glob.glob(fitsFiles)
    files.sort()
    catList = []
    matchCut = nPixels * 2.75 / 3600.
    wcs = None
    # the wcs from the first file
    wcs = WCS(fitsio.read_header(files[0]))
    if wcs == None:
        print('wcs failed reading file for',coaddID)
        return catList
    yCutCtr, xCutCtr = wcs.wcs_world2pix(ra, dec, 0)
    raLo = ra-matchCut
    raHi = ra+matchCut
    decLo = dec-matchCut
    decHi = dec+matchCut
    matchString = '(ra > {:.7f}) && (ra < {:.7f}) && (dec > {:.7f}) && (dec < {:.7f})'.format(raLo,raHi,decLo,decHi)
    for file in files:
        fits = fitsio.FITS(file)
        w = fits[1].where(matchString)
        if len(w) == 0:
            continue
        if len(fits[1][w]) == 0:
            continue
        match = fits[1][w]
        for row in match:
            # offset pixel x,y to the cutout
            row['x'] += nPixels - xCutCtr
            row['y'] += nPixels - yCutCtr
            if row['x'] > 2*nPixels or row['y'] > 2*nPixels or row['flux']< 0:
                continue
            # require a valid flux Feb 4, 2026
            if row['flux'] <= 0:
                continue
            obj = {
                'x':row['x'],
                'dx':row['dx'],
                'resx':0., # x residual, used by fitAll
                'y':row['y'],
                'dy':row['dy'],
                'resy':0., # y residual, used by fitAll
                'band':row['band'],
                'forward':row['FORWARD'],
                'epoch':row['EPOCH'],
                'mjd':row['MJDMEAN'],
                'year':Time(row['MJDMEAN'],format='mjd').byear,
                'flux':row['flux'],
                'dflux':row['dflux'],
                'fracflux':row['fracflux'],
                'primary':row['primary'],
                'flags_info':row['flags_info'],
                'nm':-2,
                'flag':0
            }
            catList.append(obj)
#        fits.close()
    return catList

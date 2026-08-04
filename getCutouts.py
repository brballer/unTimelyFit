#!/Library/Frameworks/Python.framework/Versions/Current/bin/python3

# mode = 0 --> extract the W1 band cutout for all epochs
# mode = 1 --> extract the W1 and W2 cutout for the first epoch
import os
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time

def getCutouts(ra,de,coaddID,nPixels,epoch):
    
    home = '/Users/bruceballer/Documents/plan9/'
    files = glob.glob(home+'epochCache/e*.fits')
    for f in files:
        os.remove(f)

    tCache = 'tileCache/'

    for band in range(1,3):
        fitsFile = 'unwise-' + coaddID + '-w{}-img-m.fits'.format(band)
        if not os.path.exists(home+tCache+fitsFile):
            baseURL = 'https://portal.nersc.gov/project/cosmo/temp/ameisner/neo8/'
            epochDir = 'e{:0>3d}/'.format(epoch) + coaddID[0:3] + '/' + coaddID + '/'
            import ssl
            myssl = ssl.create_default_context()
            myssl.check_hostname=False
            myssl.verify_mode=ssl.CERT_NONE
            from urllib.request import urlopen
            srcFile = baseURL + epochDir + fitsFile
            print('srcFile',srcFile)
            dstFile = home+tCache+fitsFile
            # copy the file
            g = None
            try:
                g = urlopen(srcFile,context=myssl)
            except Exception:
                print('no response from nersc portal')
                exit(0)
            with open(dstFile,'b+w') as f:
                f.write(g.read())
            f.close()
            print('Copied',fitsFile)
        # extract the cutout
        with fits.open(home+tCache+fitsFile) as hdul:
            header = hdul[0].header
            wcs = WCS(header)
            mjd = Time(header['MJDMIN'],format='mjd')
            xCutCtr, yCutCtr = wcs.wcs_world2pix(ra,de, 0)
            ixCtr = int(np.rint(xCutCtr))
            iyCtr = int(np.rint(yCutCtr))
            cutoutArray = hdul[0].section[iyCtr-nPixels:iyCtr+nPixels+1, ixCtr-nPixels:ixCtr+nPixels+1]
            # save it in epochCache
            newHDU = fits.PrimaryHDU(cutoutArray, header)
            # paste in as a decimal year
            newHDU.header['YROBS'] = mjd.byear
            outfile = './epochCache/e{:0>2d}.w{}.fits'.format(epoch,band)
            newHDU.writeto(outfile,overwrite=True)

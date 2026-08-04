#!/usr/local/bin/python3.12

import os
from astropy.io import fits
import warnings
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.table import vstack
import numpy as np
#from astropy.time import Time
from getCutouts import getCutouts
from findMinMax import findMinMax
from utfCore import utfCore
from utfUtils import printResult,mag2Flux,loadCoaddID,chkFluxVar
from comoCode import findComoCandidates

'''''
This module has evolved to do much more than simply plot the results of an unTimelyFit result. It reprocesses objects in the vicinity of the provided (ra,dec) using utfCore with debugging enabled.

In the two plots, one for W1 and one for W2, the size of a circle marker reflects the epoch of the unTimely data, small = early and large = late epoch. The PM fit result is displayed as an arrow. The results are overlaid on an image cutout, extracted from epoch 1 of a local neo8 fits file that is stored in ./tileCache. The file is downloaded from the NERSC portal if it doesn't exist locally. 

The disposition of the unTimely objects is described, e.g. how many were included or rejected in the PM fit.

Author: Bruce Baller
Created: November 2023
Modified: April 5, 2024
'''''

debug = True
dump = False
utFitv1 = True
plotFlux = False

nPixels = 2
warnings.filterwarnings('ignore',category=Warning)

plt.rcParams["figure.figsize"] = [10, 6]
plt.rcParams["figure.autolayout"] = True

#rade = '142.9953018 -4.5964981'
#'''''
rade = input('>>>>>>>> Enter RA Dec <optional CoaddID>: ')
if rade == 'q':
    exit()
#'''''
rade = rade.replace('\t',' ')
rade = rade.replace(',',' ')
list = rade.split()
if len(list) < 2:
    print('need RA Dec')
    exit()
elif len(list) != 3:
    coaddID = loadCoaddID('coaddID.txt')
else:
    # assume the input is RA Dec CoaddID
    coaddID = list[2]

ra = np.float64(list[0])
de = np.float64(list[1])

untimelyCatalogPath = 'tileCache/'
epoch = 1
cw1mag = None
cw2mag = None
cTableFile = coaddID+'_c2020.tbl'
if os.path.exists(cTableFile):
    cTbl = Table.read(cTableFile,format='ipac')
    matchCut = 3. / 3600.
    cosde = np.cos(np.pi*de/180.)
    match = cTbl[(abs(cTbl['ra']-ra)*cosde<matchCut) & (abs(cTbl['dec']-de)<matchCut)]
    if len(match) == 1:
        cw1mag = match[0]['w1mpro']
        cw2mag = match[0]['w2mpro']
        print('Use C2020 mags {:.2f} {:.2f}'.format(cw1mag,cw2mag),end=' ')
        print('fluxes {:.0f} {:.0f}'.format(mag2Flux(cw1mag),mag2Flux(cw2mag)))
result, catl = utfCore(untimelyCatalogPath,coaddID,ra,de,cw1mag,cw2mag,nPixels,debug,utFitv1)
result['coaddID'] = coaddID
# populate fluxVar
result = chkFluxVar(result,catl)

if dump:
    fout = open('dump.csv','w')
    print('ipt,nm,fwd,band,epoch,year,x,dx,resx,resx/dx,y,dy,resy,resy/dy,flux,dflux,fracFlux',file=fout)
    for ipt in range(len(catl)):
        obj = catl[ipt]
        print('{:d},{:d},{:d},{:d},{:d}'.format(ipt,obj['nm'],obj['forward'],obj['band'],obj['epoch']),end=',',file=fout)
        print('{:.3f}'.format(obj['year']),end=',',file=fout)
        print('{:.4f},{:.4f}'.format(obj['x'],obj['dx']),end=',',file=fout)
        print('{:.4f}'.format(obj['resx']),end=',',file=fout)
        print('{:.3f}'.format(obj['resx']/obj['dx']),end=',',file=fout)
        print('{:.4f},{:.4f}'.format(obj['y'],obj['dy']),end=',',file=fout)
        print('{:.4f}'.format(obj['resy']),end=',',file=fout)
        print('{:.3f}'.format(obj['resy']/obj['dy']),end=',',file=fout)
        print(int(obj['flux']),end=',',file=fout)
        print(int(obj['dflux']),end=',',file=fout)
        print('{:.3f}'.format(obj['fracflux']),file=fout)
    fout.close()
    print('wrote to dump.csv')

printResult('utFit',result,4)
if result['rchi2'] < 0:
    exit()

recQual = 0
dMov = result['sf'] - result['sfPt']
if dMov > 0.2:
    recQual = 3
elif dMov > 0.1:
    recQual = 2
elif dMov > 0.05:
    recQual = 1
print('recommended quality',recQual,'dMov {:.2f}'.format(dMov))

# comover search
if recQual > 1:
    result['nComoCan'] = findComoCandidates(result,coaddID)

try:
    getCutouts(ra,de,coaddID,nPixels,epoch)
except Exception:
    gotCutouts = False
    print('getCutouts failed')
    exit()

fig, axs = plt.subplots(1,2)

if utFitv1:
    fig.suptitle('{:.7f} {:.7f} \n utFitv1'.format(ra,de), fontsize='x-large')
else:
    fig.suptitle('{:.7f} {:.7f}'.format(ra,de), fontsize='x-large')

# define low/high epoch and year to scale the marker size
firstEpoch = 100000
lastEpoch = 0
firstYear = 0.
lastYear = 0.
# accumulate pixel x,y to offset pixel positions for plotting the trajectory
xx = []
yy = []
yrs = []
for row in catl:
    if row['nm'] < 0:
        continue
#    print('pU: {} {:.2f} {} {:.2f} {:.2f}'.format(row['nm'],row['year'],row['band'],row['x'],row['y']))
    xx.append(row['x'])
    yy.append(row['y'])
    yrs.append(row['year'])
    if row['epoch'] < firstEpoch:
        firstEpoch = row['epoch']
        firstYear = row['year']
    if row['epoch'] > lastEpoch:
        lastEpoch = row['epoch']
        lastYear = row['year']
midx = np.average(xx)
midy = np.average(yy)
midYr = np.average(yrs)
dpp = 128 # eyeball size of a pixel

# LHS two plots (W1 and W2) x,y positions
for ii, ax in enumerate(axs):
    band = ii + 1
    ax.set_title('W{}'.format(band))
    # get a representative cutout image
    fileName = './epochCache/e0{}.w{}.fits'.format(epoch,band)
    with fits.open(fileName) as hdul:
        cutout = hdul[0].data
        if cutout.shape[0] <= 0 or cutout.shape[1] <= 0:
            continue
        valMin, valMax = findMinMax(cutout, 4)
        vmin = valMin + 0.1 * (valMax - valMin)
        vmax = valMin + 1.0 * (valMax - valMin)
        im = ax.imshow(cutout, cmap="gray",origin='lower',vmin=vmin,vmax=vmax)

    # plot the x,y pixel positions
    for ipt,row in enumerate(catl):
        if row['band'] != band:
            continue
        bm1 = row['band'] - 1
        err = max(row['dx'],row['dy'])
        size = 1 + err * dpp
        if row['nm'] == -4:
            color = 'y'
        elif row['nm'] < 0:
            color = 'w'
        else:
            color = 'r'
        ax.plot(row['y'],row['x'],marker='+',color=color,markersize=size,markerfacecolor='none')
        msize = 1 * 0.5 * (1 + (row['epoch']-firstEpoch))
        ax.plot(row['y'],row['x'],marker='o',color=color,markersize=msize,fillstyle='none')

    # plot the trajectory
    if result['fitType'] > 0:
        pix2mas = 1000.*2.75
#        cutOutYearObs = 0.5 * (firstYear + lastYear)
        cutOutYearObs = midYr
        slpx = -result['pmDec']/pix2mas
        x1 = midx - (cutOutYearObs-firstYear)*slpx
        x0 = midx + (lastYear-cutOutYearObs)*slpx
        slpy = result['pmRA']/pix2mas
        y1 = midy - (cutOutYearObs-firstYear)*slpy
        y0 = midy + (lastYear-cutOutYearObs)*slpy
        ax.annotate("",xytext=(y0,x0),xy=(y1,x1),arrowprops=dict(arrowstyle='->'))

plt.show()
print() # a blank line to aid in a screen copy

tygoTblName = 'tygos.tbl'
ans = input('Enter quality. 0 = ignore, 1 (possibly), 2 (probably), 3 (definitely): ')
quality = np.int8(ans)
if quality == 0:
    exit()
elif quality < 0 or quality > 3:
    print('bad quality. Try again')
    exit()
result['quality'] = quality
if not os.path.exists(tygoTblName):
    print('making',tygoTblName)
    oTbl = Table([result],masked=False)
    oTbl.write(tygoTblName,format='ipac',overwrite=True)
else:
    bigTbl = Table.read(tygoTblName,format='ipac')
    # look for duplicates
    matchCut = 4./3600.
    selTbl = bigTbl[(abs(bigTbl['ra']-result['ra']) < matchCut) & (abs(bigTbl['dec']-result['dec']) < matchCut)]
    if len(selTbl) > 0:
        print('Duplicate entry')
        exit()
    bigTbl.add_row(result)
    bigTbl.write(tygoTblName,format='ipac',overwrite=True)
    print('Wrote',len(bigTbl),'rows to',tygoTblName)

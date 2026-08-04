#!/Library/Frameworks/Python.framework/Versions/Current/bin/python3
import os
import numpy as np
from sklearn.cluster import DBSCAN
# local
from getUntimeCat import getUntimeCat
from fitAll import fitAll
from utfUtils import fitResult, bcolors,findCoaddID
from astropy.time import Time

'''''
unTimelyFit core module usually called by processTile.py and plotUntimely.py. Passed variables are:
UntimelyCatalogPath (the path to the unTimely catalogs).
coaddID of the tile.
(ra,dec) position, usually from the CatWise2020 catalog.
catW1mag and catW2mag are the W1 and W2 magnitudes usually from the CatWise2020 catalog.
nPixels defines the size of the search region in the unTimely catalogs
debug (boolean) set False for batch processing and True for visualization (plotUntimely)

Objects within a window of 2*nPixels+1 of (ra,dec) in the unTimely catalogs are collected in the catList list. A selection is made in each band to create a catl list that is passed to the DBSCAN module to find clusters of nearby objects. 

Note that the object position errors (dx,dy) are significantly inflated for objects with very large flux to prevent abnormally high fit chisquare. The assumption is that such bright objects are already known (in Gaia for example) and are not candidates. Bright objects are useful however when comparing the unTimelyFit algorithm PM with other sources.

Note that the unTimely catalog column 'nm' is used to store the cluster ID. Values < 0 flag objects that were removed for some reason. 

checkFlux ensures that the fluxes of objects in each cluster are roughly similar.
makeUnique ensures that all objects in a band have a unique epoch.
checkOverlap compares the average cluster pixel positions in W1 and W2.

Author: Bruce Baller
Created: December 2023
Modified: 
Apr 5, 2024
Sep 23: comment out checkOverlap
Oct 30: added sepx,sepy check assuming nPixels = 2
Sep 29, 2025: Scale unTimely position errors by mag
Dec 23, 2025: Load unTimely uncertainties from a file produced by anaCalibration
Feb 6, 2026: Comment out allpars block
Feb 24, 2026: Change checkFlux pull cut to 2.5
Mar 10, 2026: Scale uncertainties if useMyErrs
Mar 14, 2026: Uncertainties were 1.4 and 1.3
'''''
def checkFlux(catl,debug):
    # remove high-flux outliers
    flx = [[],[]]
    for row in catl:
        if row['nm']< 0:
            continue
        bm1 = row['band']-1
        flx[bm1].append(row['flux'])
    for bm1 in range(2):
        if len(flx[bm1]) < 3:
            continue
        ave = np.average(flx[bm1])
        rms = np.std(flx[bm1])
        for row in catl:
            if row['band'] != bm1+1:
                continue
            pull = abs(row['flux']-ave)/rms
            if pull > 2.5:
                row['nm'] = -6
                if debug:
                    print('CheckFlux: W{} clobber nsig {:.2f} flux {:.0f}'.format(bm1+1,pull,row['flux']))

def checkOverlap(catl,nPixels,debug):
    # Ensure that the w1 cluster and the w2 cluster average positions are close to each other
    xc = [[],[]]
    yc = [[],[]]
    for row in catl:
        if row['nm']< 0:
            continue
        bm1 = row['band']-1
        xc[bm1].append(row['x'])
        yc[bm1].append(row['y'])
    if len(xc[0]) == 0 or len(xc[1]) == 0:
        return
    xc1 = np.average(xc[0])
    yc1 = np.average(yc[0])
    xc2 = np.average(xc[1])
    yc2 = np.average(yc[1])
    dx = xc1 - xc2
    dy = yc1 - yc2
    sep = dx*dx+dy*dy
    if debug:
        print('checkOverlap',end=' ')
        print('w1 {:.2f} {:.2f}'.format(xc1,yc1),end=' ')
        print('w2 {:.2f} {:.2f}'.format(xc2,yc2),end=' ')
        print('sep {:.2f}'.format(sep))
    if sep < 1.0:
        return
    # find the one closest to the center
    if debug:
        print(bcolors.red+'Clobbering cluster'+bcolors.reset)
    dx = xc1 - nPixels
    dy = yc1 - nPixels
    sep1 = dx*dx + dy*dy
    dx = xc2 - nPixels
    dy = yc2 - nPixels
    sep2 = dx*dx + dy*dy
    killBand = 1
    if sep1 < sep2:
        killBand = 2
    for row in catl:
        if row['band'] == killBand:
            row['nm'] = -5
    return catl

def makeUnique(catl,band,clid,nPixels):
    # Ensure that the catalog list only includes unique epochs in this band and cluster. Objects
    # that don't belong in this cluster are removed
    ipts = []
    epos = []
    for ipt in range(len(catl)):
        obj = catl[ipt]
        if obj['nm'] == clid and obj['band'] == band:
            ipts.append(ipt)
            epos.append(obj['epoch'])
    npts = len(ipts)
    uniqueEpochs = set(epos)
    if len(uniqueEpochs) == len(ipts):
        return
    # list of distances^2 from the center
    loDist = []
    hiDist = []
    # and the catalog entry index
    loIndx = []
    hiIndx = []
    for ii in range(npts-1):
        for jj in range(ii+1,npts):
            if catl[ipts[jj]]['epoch'] != catl[ipts[ii]]['epoch']:
                continue # unique
            imlo = ii
            imhi = jj
            dx = catl[ipts[imlo]]['x'] - nPixels
            dy = catl[ipts[imlo]]['y'] - nPixels
            loDist.append(dx*dx+dy*dy)
            loIndx.append(ipts[imlo])
            dx = catl[ipts[imhi]]['x'] - nPixels
            dy = catl[ipts[imhi]]['y'] - nPixels
            hiDist.append(dx*dx+dy*dy)
            hiIndx.append(ipts[imhi])
    lodst = np.average(loDist)
    hidst = np.average(hiDist)
    # use the points that are closest to the center
    killme = loIndx
    if lodst < hidst:
        killme = hiIndx
    for indx in killme:
        catl[indx]['nm'] = -3

def utfCore(UntimelyCatalogPath,coaddID,ra,de,catW1mag,catW2mag,nPixels,debug,useMyErrs=False):
    # The core unTimely Fit algorithm without any printing or display code.
    # Define a bad fit result - all calculated variables = 0 and fitType = -1
    result = fitResult()
    result['ra'] = ra
    result['dec'] = de
    result['fitType'] = -1 # assume no ounTimely objects exist at this position
    result['w1mag'] = catW1mag
    result['w2mag'] = catW2mag

    def distanceMatrix(catl):
        # square matrix of nearest neighbor distances used by DBScan for clustering.
        # Neighbor distances are expanded to prevent clustering objects in the same epoch
        # and to minimize clustering objects with wildly different fluxes
        dm = np.zeros((len(catl),len(catl)),dtype=np.float64)
        for ii in range(len(catl)):
            for jj in range(ii+1,len(catl)):
                if catl[ii]['epoch'] == catl[jj]['epoch']:
                    # move them 20 pixels away from each other
                    dm[ii][jj] = 20.
                else:
                    # flux asymmetry weighting after propagating the flux error of the
                    # two objects
                    dii = catl[ii]['dflux']
                    djj = catl[jj]['dflux']
                    fluwt = 1. + 0.20 * abs(catl[ii]['flux']-catl[jj]['flux'])/np.sqrt(dii*dii+djj*djj)
                    dx = catl[ii]['x'] - catl[jj]['x']
                    dy = catl[ii]['y'] - catl[jj]['y']
                    dm[ii][jj] = fluwt*np.sqrt(dx*dx+dy*dy)
                dm[jj][ii] = dm[ii][jj]
        return dm
    
    if coaddID == None:
        coaddID = findCoaddID(ra,de)
        if coaddID == 'NA':
            return result
    
    # get a list of W1 and W2 objects
    catList = getUntimeCat(UntimelyCatalogPath,ra,de,coaddID,nPixels)
    if len(catList) == 0:
        if debug:
            print('utfCore: No objects found in',UntimelyCatalogPath,coaddID,'*.fits',end=' ')
            print('near {:.7f} {:.7f}'.format(ra,de))
        return result, catList

    result['flag_info'] = catList[0]['flags_info']
    # apply the calibration offsets
    offsetsFile = 'tileCache/'+coaddID+'_offsets.npy'
    if not os.path.exists(offsetsFile):
        print('Calibration files are missing')
        exit(1)
    offsets = np.load(offsetsFile)
    for row in catList:
        bm1 = row['band']-1
        row['x'] -= offsets[0,row['epoch'],bm1]
        row['y'] -= offsets[1,row['epoch'],bm1]
    if useMyErrs:
        result['version'] = 1
        # scale the uncertainties. See studyCal3.py
        for row in catList:
            if row['band'] == 1:
                row['dx'] *= 1.5
                row['dy'] *= 1.4
            else:
                row['dx'] *= 1.1
                row['dy'] *= 1.1
    # DBScan eps in pixel space
    eps = [0.9, 1.5]
    # DBScan minimum number of points to define a cluster
    minSamples = 3
    # fit results in each band that are within 2 pixels of the center of the cutout. 
    # The cluster "position" is the average x,y position of unTimely objects in all epochs
    # in each band. 
    for band in range(1,3):
        # analyze the catlog objects for this band
        bm1 = band-1
        catl = [row for row in catList if (row['band']==band) & (row['flux']>0)]
        if len(catl) < minSamples:
            continue
        # find clusters of objects in different epochs
        dm = distanceMatrix(catl)
        db = DBSCAN(metric='precomputed',eps=eps[band-1],min_samples=minSamples).fit(dm)
        # an object tagged with a + label means it is a real cluster and
        # a label of <0 means that it failed the clustering criteria
        labels = db.labels_
        # Assign objects to clusters
        for cnt in range(len(labels)):
            # assign each object to a cluster, using an unused integer field
            catl[cnt]['nm'] = labels[cnt]
        # find the object that is closest to the center of the cutout and is in a cluster
        bestSep2 = 2.25
        useCLID = -10
        for row in catl:
            if row['nm'] < 0:
                continue
            dx = row['x'] - nPixels
            dy = row['y'] - nPixels
            sep2 = dx*dx+dy*dy
            if sep2 < bestSep2:
                bestSep2 = sep2
                useCLID = row['nm']
        if debug:
            print('band',band,'npts',len(catl),'bestSep {:.3f}'.format(np.sqrt(bestSep2)),'useCLID',useCLID)
        if useCLID < 0:
            # clobber all the clusters
            for row in catl:
                row['nm'] = -2
            continue
        # clobber the other clusters
        for row in catl:
            if row['nm'] != useCLID:
                row['nm'] = -2
        # Ensure the objects are in unique epochs
        makeUnique(catl,band,useCLID,nPixels)
        # update the full catList
        for row in catl:
            indx = row['flag']
            catList[indx]['nm'] = row['nm']
            # clear the flag
            row['flag'] = 0
    return fitAll(result,catList,debug), catList

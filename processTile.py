#!/Library/Frameworks/Python.framework/Versions/3.12/bin/python3

'''''
Runs the utfCore algorithm on all CPU's on this machine for the coaddID specified in utfUtils using the RA,Dec positions in a local CatWise2020 source file. Each CPU processes block of ~1k objects, writes the results to ./tileCache and updates a log.txt file that keeps track of completed blocks.

Author: Bruce Baller
Created: April 5, 2024
Modified:
Jan 31, 2025: 16.7 -> 16.5
Jul 14: require w1mag-w2mag > 0.23
Jul 19: require nLocs > 0 in processBlock
Jul 21: blockSize was 1000
Oct 27, 2025: Added Tahoe stuff
Jan 21, 2026: require w1mag-w2mag > 0.17
Feb 26, 2026: don't process low PMSig targets
Mar 10, 2026: Set useMyErrs True in call to utfCore
Mar 16, 2026: Removed w1, w2
'''''

import os
import glob
import warnings
import numpy as np
import multiprocessing as mp
from astropy.table import Table
from utfCore import utfCore

nPixels = 2
warnings.filterwarnings('ignore',category=Warning)

def logEntry(outPath,blockIndex):
    # Update the log file or create it if it doesn't exist
    logFile = outPath+'log.txt'
    # if it doesn't exist, interpret the block index as the
    # total count of all block indexes
    if not os.path.exists(logFile):
        with open(logFile,'w') as lf:
            lf.write('{}, tbl files expected. Completed ones listed below\n'.format(blockIndex))
        lf.close()
    else:
        with open(logFile,'a') as lf:
            lf.write('{}\n'.format(blockIndex))
        lf.close()

def getBlockList(outPath):
    # get a list of blocks that need to be processed
    logFile = outPath+'log.txt'
    locFile = outPath+'locations.npy'
    nBlocks = 0
    if not os.path.exists(logFile):
        print(logFile,'doesnt exist -> clean start. Checking',locFile,'for nBlocks')
        print('getBlockList',locFile)
        if os.path.exists(locFile):
            locations = np.load(locFile)
            nBlocks = len(locations)
            ''''' # let processBlock deal with this
            if nBlocks > 50:
                print('Too many blocks to process. Quitting')
                return []
            '''''
            logEntry(outPath,nBlocks)
        else:
            return []
    # inspect the log file
    completedBlocks = []
    with open(logFile,'r') as lf:
        # read the header having the format nBlocks,<informative text>
        header = lf.readline()
        if header.find(',') < 0:
            print('Wrong header file format.')
            exit()
        list = header.split(',')
        if nBlocks == 0:
            nBlocks = int(list[0])
        # get the list of completed blocks on the following lines
        while line := lf.readline():
            completedBlocks.append(int(line.rstrip()))
    # form a list of all blocks
    blockList = np.arange(nBlocks)
    # and ignore those that have been completed
    blockList = np.setdiff1d(blockList,completedBlocks,assume_unique=True)
    return blockList

def processBlock(args):
    # do a PM fit of all objects in one block, write the results to a
    # cache file and update the log file to flag that it has been done
    coaddID,blockIndex,nLocs,tileCache = args
    locations = np.load(tileCache+'locations.npy')
    # nLocs is the number of locations that should be processed in
    # this block. It is < 0 for production to process all objects but
    # can be set to a smaller count for debugging
    if nLocs < 0:
        nLocs = len(locations[blockIndex])
#    nLocs = 1 # debug
    if nLocs > 0:
        print('processBlock',blockIndex,'with size',nLocs)
        # make a list of fit result rows that will be appended to the output table
        rows = []
        # process all of the not-done objects
        for indx in range(nLocs):
            loc = locations[blockIndex][indx]
            if loc[0] == 0. and loc[1] == 0.:
                # hit the end of the last (partially filled) block
                break
            result, _ = utfCore(tileCache,coaddID,loc[0],loc[1],loc[2],loc[3],nPixels,debug=False,useMyErrs=True)
            # require a good fit
            if result['fitType'] > 0:
                rows.append(result)
                # temp save
        oTbl = Table(rows,masked=False)
        cacheFile = tileCache+coaddID+'_utfcat_{:03d}.tbl'.format(blockIndex)
        print(blockIndex,'->',cacheFile,'size',len(rows))
        oTbl.write(cacheFile,format='ipac',overwrite=True)
        # declare it done
        logEntry(tileCache,blockIndex)

def loadLocations(catFile,blockSize,outPath):
    # stash the ra,dec locations in a local file in the cache directory
    # as a numpy formatted file. This releases memory used by the large
    # CatWISE table.
    # See if it already exists (which is only relevant for debugging)
    locFile = outPath+'locations.npy'
    if os.path.exists(locFile):
        # open it and return the size
        locations = np.load(locFile)
        print(locFile,'exists')
        return len(locations)
    cTbl = Table.read(catFile,format='ipac')
    selTbl = cTbl[((cTbl['w1mpro']<17.) | (cTbl['w2mpro']<17.))]
    nBlocks = 1 + int(np.rint(len(selTbl)/blockSize))
    locations = np.zeros((nBlocks,blockSize,4))
    print('selTbl size',len(selTbl),'nBlocks',nBlocks,'with blockSize',blockSize)
    blockIndex = 0
    indx = 0
    for row in selTbl:
        locations[blockIndex][indx][0] = row['ra']
        locations[blockIndex][indx][1] = row['dec']
        # transfer the CatWise2020 mags into the UTF catalog
        locations[blockIndex][indx][2] = row['w1mpro']
        locations[blockIndex][indx][3] = row['w2mpro']
        indx += 1
        if indx == blockSize-1:
            blockIndex += 1
            indx = 0
    np.save(locFile,locations)
    return len(locations)

if __name__ == '__main__':
    # Process all CatWISE objects in catFile to a proper motion fit
    # using the utfCore algorithm. This is done in blocks of size
    # set below. Output catalog files are stashed in a cache directory.
    # The cache directory also contains a list of all CatWISE [ra,dec]
    # locations and a log file

    # fix for Tahoe
    mp.set_start_method('spawn',force=True)
    
    try:
        import sys
        tileCache = sys.argv[1]
        coaddID = sys.argv[2]
    except Exception:
        print('Bad args to processTile')
        exit()
    catPat = tileCache+coaddID+'_cat.tbl'
    files = glob.glob(catPat)
    if len(files) == 0:
        print('catFile not in',tileCache)
        exit()
    if len(files) > 1:
        print('Found multiple files matching this pattern in',tileCache)
        exit()
    catFile = files[0]
    if not os.path.isfile(catFile):
        raise RuntimeError('Input CatWise catalog file doesnt exist')
    # input the CatWISE table and process it in blockSize blocks
    blockSize = 1200
    print('processTile: Loading locations in',blockSize,'size blocks',end=' ')
    nBlocks = loadLocations(catFile,blockSize,tileCache)
    print(' -> {} blocks'.format(nBlocks))
    if nBlocks > 50:
        print('processTile: too many blocks')
        exit(1)
    # get a list of blocks that haven't been completed
    blockList = getBlockList(tileCache)
    if len(blockList) == 0:
        print('processTile done.')
        exit(0)
    cpuCount = mp.cpu_count()
    print('processTile: start processing cpu_count',cpuCount)
    nLocPerBlock = -1
    pool = mp.Pool(mp.cpu_count())
    jobs = []
    for blockIndex in blockList:
        args = coaddID,blockIndex,nLocPerBlock,tileCache
        jobs.append(args)
    res = pool.map(processBlock,jobs)
    res = os.system('stackCache.py '+coaddID)
    if res != 0:
        raise Exception('stackCache failed')
    print('stackCache done')

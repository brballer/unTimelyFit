#!/Library/Frameworks/Python.framework/Versions/Current/bin/python3

import sys
import os
import numpy as np

'''''
A somewhat general module for stacking data files in a local cache directory, passed by argv, that were produced by multiprocessing modules like processTile and calibrateTile. The first case is to use the astropy vstack method to merge UTF catalog table files produced by processTile (.tbl). The second case is to analyze the tile calibration pixel offset positions produced by calibrateTile (.npy).

In both cases, the resulting merged data is saved in the current directory with the format <coaddID>_*.tbl or <coaddID>_*.npy.

Note that this module also produces a text file, offdump.csv, that can be easily imported into MS Excel or MacOS Numbers.

Author: Bruce Baller
Created: January 2024
Modified: April 5, 2024
Jan 6, 2026: Comment out npy file stacking
'''''

cacheDir = sys.argv[1]
if not os.path.exists(cacheDir):
    raise RuntimeError('Directory doesnt exist: '+cacheDir)
coaddID = sys.argv[2]

# list of the file names
tblFiles = []
for fileName in os.listdir(cacheDir):
    # ignore the C2020 catalog
    if fileName.endswith('cat.tbl'):
        continue
    if fileName.startswith(coaddID+'_utfcat'):
        tblFiles.append(fileName)

if len(tblFiles) > 0:
    # sort
    tblFiles.sort()
    outName = tblFiles[0]
    if not outName.endswith('_000.tbl'):
        raise RuntimeError('The first table file doesnt have the form "_000.tbl": '+outName)
#    outName = cacheDir+'/'+outName.replace('_000','')
    outName = outName.replace('_000','')
    if os.path.exists(outName):
        print('The stacked table '+outName+' exists. Delete it and try again')
        raise RuntimeError('The stacked table '+outName+' exists')
    # stack em
    from astropy.table import vstack
    from astropy.table import Table
    bigTbl = Table()
    first = True
    for tblName in tblFiles:
        fname = cacheDir+'/'+tblName
        if first:
            first = False
            bigTbl = Table.read(fname,format='ipac')
        else:
            bigTbl = vstack([bigTbl,Table.read(fname,format='ipac')])
    bigTbl.write(outName,format='ipac',overwrite=True)
    print('stackCache wrote',len(bigTbl),'rows to',outName)

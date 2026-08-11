#!/Library/Frameworks/Python.framework/Versions/Current/bin/python3

'''''
The "mother" module that does everything except clean up after processing a tile. That should be done using the clean.py script.

Author: Bruce Baller
Created: April 5, 2024
Modified:
Oct 20: Moved time from processTile
Nov 10: Major change with coaddID
Feb 11, 2025: Added next1
Apr 11: Pass logFile to next1, saveCoaddID
Apr 28: Add makeComoCat
Sep 21: Add loRA, hiRA to allow multiple users to work in separate RA ranges
Nov 15, 2025: Quit if no DESI table
Mar 14, 2026: Use loadCalibration
'''''

import os
import re
import time

from pStar import getDESI, queryLS_DR10
from utfUtils import getUntimeFits,saveCoaddID,getCatFiles,loadCalibration
from utfUtils import CID, saveCID, next1

logFile = './coaddID.txt'

# Only process tiles in a specified RA range. Each user should only process tiles in their
# agreed upon RA range if the intent is to submit TYGOs
# user 1 range
loRA = 0.
hiRA = 359. # degrees

try:
    # see if the user has passed a coaddID
    import sys
    coaddID = sys.argv[1]
except Exception:
    # Nope. Get the last coaddID from coaddID.txt and find the next one in the specified range
    coaddID = next1(logFile,loRA,hiRA)

with open(logFile) as f:
    for line in f:
        if re.search(coaddID,line):
            print(coaddID,'exists in coaddID.txt. Already processed',line)
            exit(0)

# stash the current coaddID in a pickle file for easy access
saveCID(CID(coaddID))

# ensure that the cache directory exists
tileCache = './tileCache/'
if not os.path.exists(tileCache):
    print(tileCache,'doesnt exist. Making it')
    os.system('mkdir '+tileCache)

# See if there are leftover files from the previous job
import glob
cacheFiles = glob.glob(tileCache+'*')
if len(cacheFiles) > 0:
    print('files exist in',tileCache)
    print('chk',cacheFiles)
    ans = input('Continue?')
    if ans != 'y':
        exit(0)

tick = time.perf_counter()
# try to create the ls_dr10 table
#queryLS_DR10(tileCache,coaddID)
# try to make the DESI table for this tile
getDESI(coaddID)
if getCatFiles(tileCache, coaddID) > 0:
    os.system('osascript pingme.scpt 6309455586 "No NERSC response"')
    exit(0)
getUntimeFits(tileCache,coaddID)
rets = os.system('chkUTF.py '+coaddID)
if rets != 0:
    print('chkUTF failed. Try again')
    os.system('osascript pingme.scpt 6309455586 "job is done"')
    exit(0)
# try to create the ls_dr10 table again. This is only useful if the first attempt failed.
queryLS_DR10(tileCache,coaddID)
loadCalibration(tileCache,coaddID)
rets = os.system('caffeinate -d -i /Users/bruceballer/Documents/plan9/processTile.py '+tileCache+' '+coaddID)
if rets != 0:
    print('processTile failed',rets)
    saveCoaddID(logFile,coaddID,comment='failed')
else:
    os.system('scanTile.py '+coaddID)
    # save it for reference
    saveCoaddID(logFile,coaddID,comment='OK')
    tock = time.perf_counter()
    mins = (tock-tick)/60.
    print('Elapsed time {:.0f} minutes {:.1f} hours'.format(mins,mins/60))
    if os.path.exists(coaddID+'_candidates.txt'):
        os.system('open ' + coaddID + '_candidates.txt')
    # send me a text message
    os.system('osascript pingme.scpt 6309455586 "job is done"')

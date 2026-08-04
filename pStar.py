'''''
Functions to create,load and test DESI, DES and PS1 tables

Bruce Baller 
Modified:
Nov 1: Moved code here from getDESI.py. Added cosde
Nov 10: clean up getDESI
Feb 23, 2025: Added try in getDESI
May 23: added time.sleep
Jun 9: Replace _ICRS with _ICRS
Mar 8, 2026: use tileCtrPos
'''''
import os
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import time
from utfUtils import tileCtrPos

def readLSTbl(coaddID):
    tblFile = 'tileCache/'+coaddID+'_ls.tbl'
    if not os.path.exists(tblFile):
        return None
    return Table.read(tblFile,format='ipac')

def getObjectType(lsTbl,ra,dec,matchCut,prt=False,coaddID='',chkRemote=False):
    # queries the noirlab tractor catalog to find the type (PSF, REX, etc) of the closest
    # object within matchCut (degrees).
    # Use the remote table if the local one doesn't exist
    if lsTbl == None:
        if chkRemote:
            # local table doesn't exist so do the remote query
            from astroquery.utils.tap.core import TapPlus
            ## Create TAP service object
            tap_service = TapPlus(url="https://datalab.noirlab.edu/tap")
            query = 'SELECT ra, dec, type FROM ls_dr10.tractor'
            query = query + ' WHERE ra BETWEEN {:.7f} AND {:.7f}'.format(ra - matchCut, ra + matchCut)
            query = query + ' AND dec BETWEEN {:.7f} AND {:.7f}'.format(dec - matchCut, dec + matchCut)
            try:
                job = tap_service.launch_job(query)
                result = job.get_data()
                if len(result) == 0:
                    return None, None, None
                if prt:
                    print('ls_dr10:',end=' ')
                minSep = 100.
                minType = 'NA'
                endStr = ', '
                for cnt,r in enumerate(result):
                    dra = r['ra'] - ra
                    dde = r['dec'] - dec
                    sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                    if sep < minSep:
                        minType = r['type']
                        minSep = sep
                    if prt:
                        if cnt == len(result)-1:
                            endStr = ''
                        print(r['type'],'sep {:.1f}"'.format(sep),end=endStr)
                if prt:
                    print('->',minType)
                return minType, len(result), minSep
            except Exception:
                return None, None, None
        else:
            return None, None, None
    # use the local table
    xm = lsTbl[(lsTbl['ra'] > ra - matchCut) & (lsTbl['ra'] < ra + matchCut) &
                (lsTbl['dec'] > dec - matchCut) & (lsTbl['dec'] < dec + matchCut)]
    if len(xm) == 0:
        return None, None, None
    minSep = 100
    minType = 'NA'
    if prt:
        print('ls_dr10:',end=' ')
#        print('{:.7f} {:.7f}'.format(ra,dec),end=' ')
    cosde = np.cos(np.pi*dec/180.)
    for r in xm:
        dra = r['ra'] - ra
        dde = r['dec'] - dec
        sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
        if prt:
            print(r['type'],'sep {:.1f}"'.format(sep),end=' ')
        if sep < minSep:
            minSep = sep
            minType = r['type']
    if prt:
        print('->',minType)
    return minType, len(xm), minSep

def queryLS_DR10(tileCache,coaddID):
    # extract all ls_dr10 objects in the tile
    tblFile = tileCache+coaddID+'_ls.tbl'
    if os.path.exists(tblFile):
        print(tblFile,'exists')
        return True
    from astroquery.utils.tap.core import TapPlus
#    from astroquery.ipac.irsa import Irsa
    ra,dec = tileCtrPos(coaddID)
    print('queryLS_DR10: Querying datalab.noirlab.edu',coaddID,'ctr {:.7f} {:.7f}'.format(ra,dec))
    # get all ls_dr10 objects in the tile region
    matchCut = 0.8 # degrees 1.6 / 2
    tap_service = TapPlus(url="https://datalab.noirlab.edu/tap")
    query = 'SELECT TOP 20000000 ra, dec, type FROM ls_dr10.tractor\n'
    query = query + ' WHERE (ra BETWEEN {:.7f} AND {:.7f}'.format(ra - matchCut, ra + matchCut)
    query = query + ' AND dec BETWEEN {:.7f} AND {:.7f})\n'.format(dec - matchCut, dec + matchCut)
#    print(query)
    try:
        job = tap_service.launch_job(query)
        tbl = job.get_data()
        print('tbl size',len(tbl))
        if len(tbl) > 10:
            # save it
            tbl.write(tblFile,format='ipac',overwrite=False)
            print('wrote size',len(tbl),'to',tblFile)
        return True
    except Exception:
        return False

def readDESITbl(coaddID):
    tbl = None
    if not os.path.exists(coaddID+'_DESI.tbl'):
        return tbl
    return Table.read(coaddID+'_DESI.tbl',format='ipac')

def readDESTbl(coaddID):
    tbl = None
    if not os.path.exists(coaddID+'_DES.tbl'):
        return tbl
    return Table.read(coaddID+'_DES.tbl',format='ipac')

def readPS1Tbl(coaddID):
    tbl = None
    if not os.path.exists(coaddID+'_PS1.tbl'):
        return tbl
    return Table.read(coaddID+'_PS1.tbl',format='ipac')

def pStar(desiTbl,desTbl,pTbl,ra,dec,matchCut,prt=False):
    # returns pstar from DESI or DES or PS1 in this tile
    cosde = np.cos(np.pi*dec/180.)
    starGalClass = [-9,-1,0,0.99,1,0.90,2,0.1,3,0.001]
    if desiTbl:
        match = desiTbl[(abs(desiTbl['RAJ2000']-ra)*cosde<matchCut) & (abs(desiTbl['DEJ2000']-dec)<matchCut)]
        if prt:
            print('DESI:',len(match),'matches')
        if len(match) > 0:
            close = 1000.
            pst = match[0]['pstar']
            for row in match:
                dra = row['RAJ2000'] - ra
                dde = row['DEJ2000'] - dec
                sep = dra*dra+dde*dde
                if sep < close:
                    close = sep
                    pst = row['pstar']
            return pst

    if desTbl:
        match = desTbl[(abs(desTbl['RA_ICRS']-ra)*cosde<matchCut) & (abs(desTbl['DE_ICRS']-dec)<matchCut)]
        if prt:
            print('DES:',len(match),'matches')
        if len(match) > 0:
            close = 1000.
            clcoad = match[0]['ClCoad']
            for row in match:
                dra = row['RA_ICRS'] - ra
                dde = row['DE_ICRS'] - dec
                sep = dra*dra+dde*dde
                if sep < close:
                    close = sep
                    clcoad = row['ClCoad']
            indx = starGalClass.index(clcoad)
            return starGalClass[indx+1]

    if pTbl:
        match = pTbl[(abs(pTbl['RAJ2000']-ra)*cosde<matchCut) & (abs(pTbl['DEJ2000']-dec)<matchCut)]
        if prt:
            print('PS:',len(match),'matches')
        if len(match) > 0:
            close = 1000.
            zdiff = 1.
            for row in match:
                zd = row['zmag'] - row['zKmag']
                if isinstance(zd,np.float64):
                    dra = row['RAJ2000'] - ra
                    dde = row['DEJ2000'] - dec
                    sep = dra*dra+dde*dde
                    if sep < close:
                        close = sep
                        zdiff = zd
            if zdiff < 0.05:
                return 1.
            else:
                return 0.
    return -1

def getDESI(coaddID):
    from astropy.io import fits
    import warnings
    warnings.filterwarnings('ignore',category=Warning)

    desiTableFile = coaddID+'_DESI.tbl'
    if os.path.exists(desiTableFile):
        print('DESI file exists')
        return True

    desTableFile = coaddID+'_DES.tbl'
    ps1TableFile = coaddID+'_PS1.tbl'
    tileFieldOfView = 1.8*u.deg

    ctrRA, ctrDec = tileCtrPos(coaddID)
    print('pS.getDESI: ra/dec center in tile',coaddID,'is {:.7f} {:.7f}'.format(ctrRA,ctrDec))
    skyCoord = SkyCoord(ctrRA, ctrDec, unit='deg', frame='icrs')

    dv = Vizier(columns=['RAJ2000','DEJ2000','pstar'])
    dv.ROW_LIMIT = -1
    ntries = 0
    tableList = None
    while ntries < 4:
        try:
            tableList = dv.query_region(skyCoord,width=tileFieldOfView,catalog=['VII/292/south','VII/292/north'])
            break
        except Exception as e:
            print('error code',e)
            print('getDESI: query_region failed. Server down?',dv.VIZIER_SERVER,'ntries',ntries)
            ntries += 1
            if ntries == 3:
                print('giving up')
                return False
            time.sleep(4)
        if tableList:
            tbl = tableList[0]
            if len(tbl) < 20000:
                print('DESI table too small',len(tbl))
                return False
            tbl.write(desiTableFile,format='ipac',overwrite=False)
            print('wrote size',len(tbl),'to',desiTableFile)
            return True
        else:
            print('no DESI objects found in Vizier')
            return False

    # make the DES table
    if not os.path.exists(desTableFile):
        pv = Vizier(columns=['RA_ICRS','DE_ICRS','ClCoad'])
        pv.ROW_LIMIT = -1
        tableList = None
        try:
            tableList = pv.query_region(skyCoord,width=tileFieldOfView,catalog=['II/371/des_dr2'])
        except Exception:
            print('DES table query failed')
            return False
        if tableList:
            tbl = tableList[0]
            tbl.write(desTableFile,format='ipac',overwrite=False)
            print('wrote size',len(tbl),'to',desTableFile)
            return True
        else:
            print('no DES objects found in Vizier')
        return False
        
    # make the PS1 table
    if not os.path.exists(ps1TableFile):
        pv = Vizier(columns=['RAJ2000','DEJ2000','zmag','zKmag'])
        pv.ROW_LIMIT = -1
        for ntry in range(4):
            try:
                print('try ps1')
                tableList = pv.query_region(skyCoord,width=tileFieldOfView,catalog=['II/349/ps1'])
                if tableList:
                    tbl = tableList[0]
                    tbl.write(ps1TableFile,format='ipac',overwrite=False)
                    print('wrote size',len(tbl),'to',ps1TableFile)
                    return True
                else:
                    print('no PS1 objects found in Vizier')
                    return False
            except Exception:
                print('No Vizier response. waiting 4 seconds')
                time.sleep(4)
        return False
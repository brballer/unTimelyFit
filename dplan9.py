#!/Library/Frameworks/Python.framework/Versions/Current/bin/python3

'''''
This module queries Vizier and Simbad and prints a compact summary of the results in the vicinity of a provided (ra,dec) position. A small number of the most relevant catalogs are queried to improve response time.  The results from a catalog are usually printed on one line per object.

Author: Bruce Baller
Created: Summer 2023
Modified: 
Apr 6, 2024 Increase matchCut from 3 to 4 arcsec
Sep 18: Print tygo log info line using WISE info if C2020 doesn't exist
Sep 19: Print C2020 - Spitzer PM estimate
Mar 18, 2025: Add getObjectType
Apr 8, 2025: Work on PS1 output
Jun 9: change _ICRS to _ICRS
Nov 3, 2025: Add wise w extended source flag
Jul 11, 2026: fixed logic error using getObjectType. 'NA' was counted as a galaxy
Jul 15, 2025: Added chkRemote
Aug 5, 2026: Added Simbad webbrowser
'''''

import warnings
import clipboard
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from utfUtils import bcolors,spType,flux2Mag,photDistWISE,photDistJK,findCoaddID
from pStar import getObjectType,readLSTbl
import time
#from getCoaddID import getCoaddID

# ******** Main **********

# ignore Simbad warnings
warnings.filterwarnings('ignore', category=Warning, append=True)

# Brown dwarf catalogs
bdCatalogs = ['J/ApJS/271/55', # 20 pc survey
              'J/ApJ/889/74',
              'J/ApJ/842/118',
              'J/ApJ/858/41',
              'J/ApJ/867/109',
              'J/ApJ/883/205',
              'J/ApJ/889/74',
              'J/ApJ/899/123',
              'J/ApJ/934/178'
            ]
# White dwarf catalogs
wdCatalogs = 'J/MNRAS/508/3877/maincat'

allCatalogs = ['II/349/ps1', # PS1
                'II/365/catwise',
                'I/353/gsc242', # Guide Star catalog
                'I/355/gaiadr3', 'V/154/sdss16', 'I/351/gps1_p', 
                'I/324/igsl3', # Initial Gaia source list
                'II/319/las9', # UKIDSS-DR9
                'II/367/vhs_dr5', # VISTA hemisphere survey (VHS)
                'J/ApJ/651/502', # Spitzer IRAC photometry of M, L and T Dwarfs
                'II/371/des_dr2',
                'II/246/out', # 2MASS
                'VII/292/south', 'VII/292/north', # DESI
                'II/349/ps1',
                'I/297/out', # NOMAD
                'I/317/sample', # PPMXL
                'II/379/smssdr4', # SkyMapper
                'II/381/hlsp_ps1_tm' # Probabilistic Classifications of Unresolved Point Sources in PanSTARRS1
               ]

viz = Vizier()
viz.ROW_LIMIT = 20
viz.TIMEOUT = 3600
#viz.VIZIER_SERVER = 'viz.nao.ac.jp'
#viz.VIZIER_SERVER = 'viz.cfa.harvard.edu'

matchCut = 4./3600.
fieldOfView = matchCut*u.degree
print('FOV radius:',fieldOfView.to(u.arcsec))

raList = [] # list of RAs that we have visited

# start an infinite loop that can be interrupted
prev_rafl = -1.
prev_defl = -1.
prevCoaddID = 'na'

w1mag = 0.
w2mag = 0.

lsTbl = None

while True:
    print('******************************')
    ra = input('Enter RA Dec <coaddID>, or q to quit: ')
    if len(ra) == 0:
        exit()
    list = ra.split()
    if len(list) == 1 and not ra.isnumeric():
        exit()
    # remove commas
    ra = ra.replace(',',' ')
    # remove tabs
    ra = ra.replace('\t',' ')
    ra = ra.replace('\n',' ')
    clipboard.copy(ra+'\n')
    # look for an embedded space character which indicates that both
    # RA and DE were entered on one line without an embedded newline
    indx = ra.find(' ')
    if indx > 0:
        # the RA and DE were entered on one line wo a newline. 
        list = ra.split()
        ra = list[0]
        de = list[1]
        # truncate to 7 significant figures
        try:
            rafl = np.float64(ra)
        except Exception:
            exit()
        ra = str(round(rafl,7))
        defl = np.float64(de)
        de = str(round(defl,7))
    else:
        de = input() # no prompt for DE

    # see if the coaddID was passed
    if len(list) == 3:
        coaddID = list[2]
    else:
        coaddID = findCoaddID(rafl,defl)
    print('in tile',coaddID)
    if lsTbl == None:
        lsTbl = readLSTbl(coaddID)
    if lsTbl == None:
        print('no LS table')
    else:
        print('lsTbl size',len(lsTbl))
    prevCoaddID = coaddID

    # try to construct a sky coordinate
    try:
        skyCoord = SkyCoord(rafl, defl, unit='deg', frame='icrs')
    except:
        exit()
#        continue
    # make float versions of ra and de
    try:
        rafl = np.float64(ra)
    except Exception:
        exit()
    defl = np.float64(de)
    cosde = np.cos(np.pi*defl/180.)

    # look for this RA in the already-visited list
    if raList.count(rafl) > 0:
        print('You have already visited this position...')
        continue
    raList.append(rafl)

    # reasons to suggest rejecting this position
    cntGal = 0
    cntStar = 0
    inGaia = False
    inBDCat = False
    artifactFlag = False
    isFaint = False
    lowPMSig = False
    isMType = False
    logInfo = 'NA'
    # check t he more recent DESI catalog than the one in Vizier.
    # lsTbl is None so this will trigger a remote query
    objType, nObjs, sep = getObjectType(lsTbl,rafl,defl,matchCut,prt=False,coaddID=coaddID,chkRemote=True)
    if objType != None:
        print('LS_dr10:',objType,'sep {:.1f}"'.format(sep),'nObjs',nObjs,end=' ')
        if objType == 'PSF':
            cntStar += 1
            print('-> star')
        elif objType != 'NA':
            cntGal += 1
            print('-> galaxy')
        else:
            print()
    # look for objects in all catalogs
    tableList = None
    for ntry in range(1):
        try:
#            print('trying',viz.VIZIER_SERVER)
            tableList = viz.query_region(skyCoord,radius=fieldOfView,catalog=allCatalogs)
            break
        except Exception:
            print('No Vizier response. Waiting 4 seconds to try again')
            time.sleep(4)
    if not tableList:
        print('Empty Vizier tableList')
        continue
    for table in tableList:
        if table.meta['name'] == 'I/322A/out':
            row = table[0]
            print('UCAC4: PM {:.0f} {:.0f}'.format(row['pmRA'],row['pmDE']),end=' ')
            print()
        elif table.meta['name'] == 'II/319/las9':
            row = table[0]
#            print(row.colnames)
            print('UKIDSS:',end=' ')
            dra = row['RAJ2000'] - rafl
            dde = row['DEJ2000'] - defl
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            print('sep {:.1f}"'.format(sep),end=' ')
            print('PM {:.0f} {:.0f}'.format(row['pmRA'],row['pmDE']),end=' ')
            print('Epoch {:.1f}'.format(row['Epoch']),end=' ')
            if row['cl'] == -3:
                print('prob galaxy')
                cntGal += 0.7
            elif row['cl'] == -2:
                print('prob star')
                cntStar += 0.7
            elif row['cl'] == 0:
                print('noise')
            elif row['cl'] == -1:
                print('star')
                cntStar += 1
            elif row['cl'] == 1:
                print('galaxy')
                cntGal += 1
        elif table.meta['name'].find('367/vhs_dr5') > 0:
#            print(table[0].colnames)
            for row in table:
                print('VHS5:',end=' ')
                if row['Mclass'] == -3:
                    print('Probable galaxy',end=' ')
                    cntGal += 0.5
                elif row['Mclass'] == -2:
                    print('Probable star',end=' ')
                    cntStar += 0.5
                elif row['Mclass'] == -1:
                    print('Star',end=' ')
                    cntStar += 1
                elif row['Mclass'] == 1:
                    print('Galaxy',end=' ')
                    cntGal += 1
                dra = row['RAJ2000'] - rafl
                dde = row['DEJ2000'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                print('Jpmag {:.1f}'.format(row['Jpmag']),end=' ')
                print('Hmag {:.1f}'.format(row['Hap3']),end=' ')
                print('Kspmag {:.1f}'.format(row['Kspmag']),end=' ')
                print('sep {:.1f}"'.format(sep))
                if row['Jpmag'] > 0:
                    print(' ',end=' ') # indent the photometric distance line
                    photDistJK(w1mag,w2mag,row['Jpmag'],0.1,row['Ksap3'],0.1,True)
                else:
                    print()
        elif table.meta['name'].startswith('II/328/allwise'):
            row = table[0]
#            print(row.colnames)
            print('AllWise:',end=' ')
            print('W1 {:.2f} W2 {:.2f} W3 {:.2f} W4 {:.2f}'.format(row['W1mag'],row['W2mag'],row['W3mag'],row['W4mag']),end=' ')
            print('J {:.2f} H {:.2f} K {:.2f}'.format(row['Jmag'],row['Hmag'],row['Kmag']),end=' ')
            dra = row['RAJ2000'] - rafl
            dde = row['DEJ2000'] - defl
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            print('sep {:.1f}"'.format(sep))
#            print('   ext src flag',row['ex'],'variability flag',row['var'])
        elif table.meta['name'].startswith('I/351/gps1_p'):
            row = table[0]
            print('GPS1_P:',end=' ')
            print('PM {:.0f} \u00B1 {:.0f} {:.0f} \u00B1 {:.0f}'.format(row['pmRA'],row['e_pmRA'],row['pmDE'],row['e_pmDE']),end=' ')
            print('{} SOURCES,'.format(len(table)))
#            print(row.colnames)
        elif table.meta['name'].startswith('II/381/hlsp_ps1_tm'):
            row = table[0]
#            print(row.colnames)
            print('PS1-PSC: psScore {:.2f}'.format(row['psScore']),end=' ')
            dra = row['RAJ2000'] - rafl
            dde = row['DEJ2000'] - defl
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            print('sep {:.1f}"'.format(sep))
            cntStar += row['psScore']
            cntGal += (1-row['psScore'])
        elif table.meta['name'].find('/mis') > 0:
            row = table[0]
            print('GALEX UV source - white dwarf?')
        elif table.meta['name'].startswith('II/379/smssdr4'):
            print('SkyMapper:',end=' ')
            for row in table:
                dra = row['RAICRS'] - rafl
                dde = row['DEICRS'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                print('sep {:.1f}" Stellarity {:.2f}'.format(sep,row['ClassStar']),end=' ')
                cntStar += row['ClassStar']
                cntGal += (1-row['ClassStar'])
            print()
        elif table.meta['name'].startswith('I/353/gsc242'):
#            print(table[0].colnames)
            # https://arxiv.org/pdf/0807.2522
            # positional precision ~0.2" (table 14)
            minSep = 2.
            codes = ['Star','Galaxy','Blend','Non-Star','Unk','Defect']
            minCls = -1
            print('GSC:',end=' ')
            endStr = ', '
            for cnt, r in enumerate(table):
                dra = r['RA_ICRS'] - rafl
                dde = r['DE_ICRS'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                if sep < minSep:
                    minSep = sep
                    minCls = r['Class']
                if cnt == len(table)-1:
                    endStr = ''
                print(codes[r['Class']],'sep {:0.1f}"'.format(sep),end=endStr)
                if isinstance(r['e'],float):
                    print('eccentrity {:.3f}'.format(r['e']),end=endStr)
            if minCls == 0:
                cntStar += 1
                print(' -> star')
            elif minCls == 1:
                cntGal += 1
                print(' -> galaxy')
            elif minCls == 3:
                cntStar += 0.5
                print(' -> star?')
            else:
                print()
        elif table.meta['name'].find('des_dr2') > 0:
            # https://ui.adsabs.harvard.edu/abs/2021ApJS..255...20A/abstract
            # PSF FWHM ~ 1"
#            print(table[0].colnames)
            starGalClass = [
                'hi-confidence star',
                'Star',
                'Galaxy',
                'hi-confidence galaxy']
            swghts = [2,1,0,0] # star weighting
            gwghts = [0,0,1,2] # galaxy weighting
            print('DES:',end=' ')
            minSep = 2.
            minIndx = -1
            endStr = ','
            for cnt, row in enumerate(table):
                if row['ClCoad'] < 0:
                    # no data
                    continue
                dra = row['RA_ICRS'] - rafl
                dde = row['DE_ICRS'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                indx = row['ClCoad']
                if sep < minSep:
                    minIndx = indx
                    minSep = sep
                print(starGalClass[indx],end=' ')
                if cnt == len(table)-1:
                    endStr = ''
                print('sep {:.1f}"'.format(sep),end=endStr)
            print()
            if minIndx > -1:
                cntStar += swghts[minIndx]
                cntGal += gwghts[minIndx]
        elif table.meta['name'].startswith('I/297/out'):
            print('NOMAD:',end=' ')
            for row in table:
                dra = row['RAJ2000'] - rafl
                dde = row['DEJ2000'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde+dde*dde)*3600
                print('sep {:0.1f}" PM {:.0f} {:.0f}'.format(sep,row['pmRA'],row['pmDE']),end=' ')
            print()
        elif table.meta['name'].startswith('I/317/sample'):
            print('PPMXL:',end=' ')
            for row in table:
                dra = row['RAJ2000'] - rafl
                dde = row['DEJ2000'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde+dde*dde)*3600
                print('sep {:0.1f}" PM {:.0f} {:.0f}'.format(sep,row['pmRA'],row['pmDE']),end=' ')
            print()
        elif table.meta['name'].startswith('VII/292/'):
            # Legacy Survey DR8
            # https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.3662D/abstract
            print('DESI_dr8:',end=' ')
            cnt = 0
            minSep = 2.
            pStar = -1
            type = 'NA'
            for row in table:
                dra = row['RAJ2000'] - rafl
                dde = row['DEJ2000'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde+dde*dde)*3600
                if sep < minSep:
                    minSep = sep
                    pStar = row['pstar']
                cnt += 1
                endStr = ', '
                if cnt == len(table):
                    endStr = ''
                print(row['type'],'pstar {:.2f}'.format(row['pstar']),end=' ')
                print('sep {:.1f}"'.format(sep),end=endStr)
            print()
        elif table.meta['name'].find('/246/out') > 0:
            if len(table) > 1:
                print('multiple 2MASS rows')
            row = table[0]
#            print(row.colnames)
            print('2MASS:',end=' ')
            # C2020 epoch ~2015.5 and 2MASS epoch ~1999.9, so dyr ~15.7 years
            # convert dra,dde to mas/yr = 1000/15.7
            scale = 3600 * 1000/15.7
            dra = (rafl - row['RAJ2000']) * scale
            dde = (defl - row['DEJ2000']) * scale
            print('PM {:.0f} {:.0f} (C2020(2015.5) - 1999.7)'.format(dra,dde))
            if w1mag > 0:
                print(' ',end=' ') # indent the photometric distance line
                photDistJK(w1mag,w2mag,row['Jmag'],0.1,row['Kmag'],0.1,True)
            cntStar += 1
        elif table.meta['name'].startswith('II/349/ps1'):
#            print(table[0].colnames)
            # PanSTARRS
            strng = 'PS1: '
            minSep = 8.
            minSepZdiff = 10.
            for row in table:
                dra = row['RAJ2000'] - rafl
                dde = row['DEJ2000'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                flag = row['f_objID']
                if (flag & 4) or (flag & 8) or (flag & 16):
                    strng = strng + 'QSO '
                if (flag & 32) or (flag & 64):
                    strng = strng + 'RR Lyra '
                if (flag & 128) or (flag & 256):
                    strng = strng + 'variable '
                if (flag & 512) or (flag & 1024):
                    strng = strng + 'asteroid '
                if flag & 2048:
                    strng = strng + 'hiPM star '
                zdiff = row['zmag'] - row['zKmag']
                if not isinstance(zdiff,np.float64):
                    strng = strng + 'unk '
                    continue
                if zdiff < 0.05:
                    strng = strng + 'star '
                else:
                    strng = strng + 'gal '
                strng = strng + 'sep {:.2f}",'.format(sep)
                strng = strng + 'qual {:}'.format(row['Qual'])
                if sep < minSep:
                    minSep = sep
                    minSepZdiff = zdiff
            strng = strng.rstrip(',')
        elif table.meta['name'].find('/wise') > 0:
            row = table[0]
#            print(row.colnames)
            print('WISE:',end=' ')
            if row['ccf'].find('D') >= 0:
                print('Diffraction spike',end=' ')
            elif row['ccf'].find('P') >= 0:
                print('Persistence',end=' ')
            elif row['ccf'].find('H') >= 0:
                print('Halo',end=' ')
            elif row['ccf'].find('O') >= 0:
                print('Optical ghost',end=' ')
            print('Var flag (0-9)',row['var'],end=' ')
            print('w1 {:.2f}'.format(row['W1mag']),end=' ')
            print('w2 {:.2f}'.format(row['W2mag']),end=' ')
            print('w3 {:.2f}'.format(row['W3mag']),end=' ')
            print('w4 {:.2f}'.format(row['W4mag']),end=' ')
            if row['ex'] > 0:
                print(bcolors.red+'ext src Flag {}'.format(row['ex'])+bcolors.reset,end=' ')
            else:
                print('ext src Flag {}'.format(row['ex']),end=' ')
            if logInfo == 'NA':
                logInfo = '{:.7f} {:.7f} {:.2f} {:.2f}'.format(rafl,defl,row['W1mag']-row['W2mag'],row['W1mag'])
            print()
#            '''''
        elif table.meta['name'].find('/catwise') > 0:
#            print(table[0].colnames)
            minSep = 4.
            useRow = 0
            for ii, row in enumerate(table):
                dra = row['RA_ICRS'] - rafl
                dde = row['DE_ICRS'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                if sep < minSep:
                    useRow = ii
                    minSep = sep
            for ii, row in enumerate(table):
                pmRA = row['pmRA']*1000
                epmRA = row['e_pmRA']*1000
                pmDec = row['pmDE']*1000
                epmDec = row['e_pmDE']*1000
                sigRA = abs(pmRA)/epmRA
                sigDec = abs(pmDec)/epmDec
                PM = np.sqrt(pmRA*pmRA + epmDec*epmDec)
                pmSig = np.sqrt(sigRA*sigRA+sigDec*sigDec)
                if useRow == ii:
                    if pmSig < 0.8:
                        lowPMSig = True
                    else:
                        lowPMSig = False
                    isMType = row['W1mproPM']-row['W2mproPM'] < 0.23
                    isFaint = row['W1mproPM'] > 16.7 and row['W2mproPM'] > 16.7
                print('C2020: PM {:.0f} \u00B1 {:.0f}'.format(pmRA,epmRA),end=' ')
                print('{:.0f} \u00B1 {:.0f}'.format(pmDec,epmDec),end=' ')
                dra = row['RA_ICRS'] - rafl
                dde = row['DE_ICRS'] - defl
                sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                print('w1 {:.2f} w2 {:.2f}'.format(row['W1mproPM'],row['W2mproPM']),end=' ')
                print('pmSig {:.2f}'.format(pmSig),end=' ')
#                print('rc2s {:.2f} {:.2f}'.format(row['chi2pmRA'],row['chi2pmDE']),end=' ')
                sptyp, indx = spType(row['W1mproPM'],row['W2mproPM'])
                if sptyp:
                    print('type',sptyp,end=' ')
                print('sep {:.1f}"'.format(sep))
                if not row['abf'][0].isdigit():
                    print(' artifact flag',row['abf'])
                    artifactFlag = True
                w1w2 = row['W1mproPM']-row['W2mproPM']
                # save for use elsewhere
                w1mag = row['W1mproPM']
                w2mag = row['W2mproPM']
                logInfo = '{:.7f} {:.7f} {:.2f}'.format(rafl,defl,w1w2)
                print(' ',end=' ') # indent the photometric distance line
                photDistWISE(row['W1mproPM'],row['W2mproPM'],True)
        elif table.meta['name'].find('sdss16') > 0:
            print('SDSS16:',end=' ')
            if len(table) > 1:
                print(len(table),'Observations',end=' ')
            useRow = table[0]
            minSep = 2.
            if len(table) > 0:
                cnt = 0
                for row in table:
                    dra = row['RA_ICRS'] - rafl
                    dde = row['DE_ICRS'] - defl
                    sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
                    if sep < minSep:
                        minSep = sep
                        useRow = row
                    cnt += 1
            row = useRow
            print('gmag {:.1f} zmag {:.1f}'.format(row['gmag'],row['zmag']),end=' ')
            print('class',row['class'],end=' ')
            if row['class'] == 6:
                print('-> star',end=' ')
                cntStar += 1
            elif row['class'] == 3:
                print('-> galaxy',end=' ')
                cntGal += 1
            elif row['class'] == 0:
                print('-> unk',end=' ')
            print('sep {:.1f}"'.format(sep))
        elif table.meta['name'].startswith('I/324/igsl3'):
            row = table[0]
#            print(row.colnames)
            print('IGSL3:',end=' ')
            dra = row['RAJ2000'] - rafl
            dde = row['DEJ2000'] - defl
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            print('sep {:.1f}"'.format(sep),end=' ')
            print('PM {:.0f} {:.0f}'.format(row['pmRA'],row['pmDE']))
        elif table.meta['name'].find('gaiadr3') > 0:
            row = table[0]
#            print(row.colnames)
            print('GAIA: PM {:.0f} \u00B1 {:.0f} {:.0f} \u00B1 {:.0f}'.format(row['pmRA'],row['e_pmRA'],row['pmDE'],row['e_pmDE']),end=' ')
            print('gmag {:.1f}'.format(row['Gmag']),end=' ')
            print('Plx {:.2f} \u00B1 {:.2f}'.format(row['Plx'],row['e_Plx']),end=' ')
            if row['Plx'] > 0:
                dist = 1000. / row['Plx']
                print('-> {:.0f} pc'.format(dist),end=' ')
            dra = row['RA_ICRS'] - rafl
            dde = row['DE_ICRS'] - defl
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            if sep < 10.:
                inGaia = True
            print('sep {:.1f}"'.format(sep),end=' ')
            if isinstance(row['RUWE'],np.float64):
                print('RUWE {:.2f}'.format(row['RUWE']),end=' ')
                if row['RUWE'] > 1.4:
                    cntGal += 1
                else:
                    cntStar += 1
#            print(row['Source'])
            print()
        elif table.meta['name'].find('/502') > 0:
            print('Spitzer BD catalog')
        elif table.meta['name'].find('368/sstsl2') > 0:
            row = table[0]
            # ch1 (3.6 um), ch2 (4.5 um)
#            print(row.colnames)
            print('Spitzer:',end=' ')
#            print('SMID',row['SMID'],end=' ')
            strng = str(row['SMID'])
            firstChar = strng[:1]
            if firstChar == '2':
                print('Zodiacal light',end=' ')
            elif firstChar == '3':
                print('Galactic plane',end=' ')
            elif firstChar == '4':
                print('Galactic source',end=' ')
                cntStar += 1
            elif firstChar == '5':
                print('Extnded source',end=' ')
                cntGal += 1
            elif firstChar == '6':
                print('Extragalactic source',end=' ')
                cntGal += 1
            # C2020 epoch ~2015.5 and Spitzer epoch ~2003.5, so dyr ~12 yearssp
            # convert dra,dde to mas/yr = 1000/12
            scale = 3600 * 1000/12
            dra = rafl - row['RA_ICRS']
            dde = defl - row['DE_ICRS']
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            pmra = dra * scale
            pmde = dde * scale
#            print('PM {:.0f} {:.0f} (C2020(2015.5) - 2003.5)'.format(pmra,pmde),end=' ')
            print('sep {:.1f}"'.format(sep))
        elif table.meta['name'].find('/apop') > 0:
            row = table[0]
#            print(row.colnames)
            print('APOP:',end=' ')
            print('PM {:.0f}\u00B1{:.0f}'.format(row['pmRA'],row['e_pmRA']),end=' ')
            print('{:.0f}\u00B1{:.0f}'.format(row['pmDE'],row['e_pmDE']),end=' ')
            print()
        elif table.meta['name'].find('363') > 0:
            row = table[0]
#            print(row.colnames)
            print('unWISE',end=' ')
#            print('FW1 {:.0f\u00B1{:.0f}'.format(row['FW1'],row['e_FW1']),end=' ')
            w1mag = flux2Mag(row['FW1'])
            w2mag = flux2Mag(row['FW2'])
            print('W1: mag {:.2f}'.format(w1mag),end=' ')
            print('quality {:.2f}'.format(row['q_W1']),end=' ')
            print('fracFlux {:.2f}'.format(row['fFW1']),end=' ')
            print('W2: mag {:.2f}'.format(w2mag),end=' ')
#            print('FW2 {:.0f}\u00B1{:.0f}'.format(row['FW2'],row['e_FW2']),end=' ')
            print('quality {:.2f}'.format(row['q_W2']),end=' ')
            print('fracFlux {:.2f}'.format(row['fFW2']),end=' ')
            dra = row['RA_ICRS'] - rafl
            dde = row['DE_ICRS'] - defl
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)*3600
            print('sep {:.1f}"'.format(sep))
            if logInfo == 'NA':
                logInfo = '{:.7f} {:.7f} {:.2f} {:.2f}'.format(rafl,defl,w1mag-w2mag,w1mag)
        else:
            print('table not decoded',table.meta['name'])
    #chkVarCat(rafl,defl,True)
    # color-color BD selections from arXiv 2302.15156
    print('cntStar {:.2f} cntGal {:.2f}'.format(cntStar,cntGal))
    sTableList = None
    for ntry in range(2):
        try:
            sTableList = Simbad.query_region(skyCoord,radius=4*fieldOfView)
            break
        except Exception:
            print('No Simbad response. Waiting 4 seconds to try again')
            time.sleep(4)
    if sTableList:
        print('Simbad object:',end=' ')
#            print(sTableList)
        for row in sTableList:
            print(row[0],end=' ')
        print()
        print(bcolors.red+'In Simbad --> CHECK IT OUT'+bcolors.reset)
        import webbrowser
        webbrowser.open('https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={:.7f} {:.7f}+&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=20&Radius.unit=arcsec&submit=submit+query&CoordList='.format(rafl,defl))
    else:
        print('not in Simbad')

    if isMType:
        print(bcolors.red+'M type'+bcolors.reset)        
    if inGaia:
        print(bcolors.red+'In Gaia'+bcolors.reset)
    elif cntGal > cntStar:
        print(bcolors.red+'Looks like a galaxy'+bcolors.reset)
    elif isFaint:
        print(bcolors.red+'Too faint'+bcolors.reset)
    elif inBDCat:
        print(bcolors.red+'Known BD'+bcolors.reset)        
    print('Locations visited',len(raList))

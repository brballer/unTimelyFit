#!/Library/Frameworks/Python.framework/Versions/Current/bin/python3
import numpy as np
import os

'''''
Utility functions used by many modules

Author: Bruce Baller
Created: December 2023
Modification history
Aug 1: Change call to profileHist to use 's'/'' instead of True/False
Aug 7: minObs was 50
Oct 30: added falsePMProb
Nov 23: Added baseURLs in getUntimeFits
Jan 31, 2025: Added saveUntimelyInfo
Feb 11: Added next1, parallaxAlg
Mar 17: Added getObjectType
Apr 11: Pass fileName to next1, loadCoaddID, saveCoaddID
Apr 11: Add isInFile
Apr 18: Add addTimeStamp
May 3: Deleted saveUnTimelyInfo
May 28: Changes to parallaxAlg
Jun 9: change J2000 to ICRS
Sep 21, 2025: Pass loRA, hiRA to next1
Dec 1, 2025: Replace FalsePMProb with rchi2Pt
Feb 27, 2026: Add sf to fitResult and printResult
Mar 5, 2026: Include PM fit rejected points in the Point fit
Mar 8, 2026: Added tileCtrPos, removed getTileCenter
Mar 14, 2026: Added loadCalibration
Mar 23: 2026: Moved getObjectType to pStar
May 5, 2026: Was result['sf'] - result['sfPt'] > 0.05
May 17, 2026: Add fitResult version
Jul 5, 2026: Add makeC2020Table, modifications to makeGaiaTable both using IRSA
'''''

def findCoaddID(ra,dec):
    # replacement for getCoaddID
    from astropy.io import fits
    cosde = np.cos(np.pi*dec/180.)
    with fits.open('allsky-atlas.fits') as hdul:
        table = hdul[1].data
        sel = table[(abs(table['ra']-ra)<4.0) & (abs(table['dec']-dec)<4.0)]
        minsep = 100.
        cid = 'NA'
        for row in sel:
            dra = row['ra'] - ra
            dde = row['dec'] - dec
            sep = np.sqrt(dra*dra*cosde*cosde + dde*dde)
            if sep < minsep:
                minsep = sep
                cid = row['coadd_id']
    return cid

def loadCalibration(tileCache,coaddID):
    # Get the calibration from Merged_Summary_Catalog
    calFile = tileCache+coaddID + '_offsets.npy'
    if os.path.exists(calFile):
        return True
    from astropy.io import fits
    with fits.open('Merged_Summary_Catalog.fits') as hdul:
        nEpochs = 50
        offs = np.zeros((2,nEpochs,2),dtype=np.float64)
        tbl = hdul[1].data
        stbl = tbl[(tbl['COADD_ID'] == coaddID)]
        if len(stbl) == 0:
            return False
        for r in stbl:
            bm1 = r['BAND'] - 1
            offs[0,r['EPOCH'],bm1] = r['RECAL_X_SHIFT']
            offs[1,r['EPOCH'],bm1] = r['RECAL_Y_SHIFT']
    np.save(calFile,offs)
    return True


def tileCtrPos(coaddID):
    # return the ra,dec of the center of the tile
    from astropy.io import fits
    with fits.open('allsky-atlas.fits') as hdul:
        atable = hdul[1].data
        row = atable[(atable['coadd_id'] == coaddID)]
        return row['ra'][0], row['dec'][0]

class CID():
    def __init__(self, coaddID):
        self.coaddID = coaddID

def saveCID(coaddID):
    import pickle
    try:
        with open('cid.pickle','wb') as f:
            pickle.dump(coaddID,f,protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print('Error saving coaddID',ex)

def loadCID():
    import pickle
    try:
        with open('cid.pickle','rb') as f:
            return pickle.load(f)
    except Exception as ex:
        print('Error loading coaddID',ex)

def makeC2020Table(destPath,coaddID, prt=False):
    # make the table in the destPath directory if it doesn't exist
    # see c2020Fields.txt for a column listing
    cTableFile = destPath + coaddID + '_c2020.tbl'
    if prt:
        print('makeC2020Table',cTableFile)
    if os.path.exists(cTableFile):
        if prt:
            print(cTableFile,'exists')
        return True
    from astropy.coordinates import SkyCoord
    from astroquery.ipac.irsa import Irsa
    import astropy.units as u
    ctrRA, ctrDec = tileCtrPos(coaddID)
    skyCoord = SkyCoord(ctrRA, ctrDec, unit='deg', frame='icrs')
    columns = 'source_id,ra,dec,w1mpro,w2mpro,pmra,sigpmra,pmdec,sigpmdec,pmcode,niters_pm'
    tbl = Irsa.query_region(catalog='catwise_2020',coordinates=skyCoord,spatial='Box',width=1.8*u.deg,columns=columns)
    if not tbl:
        print('mC2020T table query failed')
        return False
    from astropy.table import Table
    otbl = Table()
    print('makeC2020Table size',len(tbl))
    for row in tbl:
        if not row['source_id'].startswith(coaddID):
            continue
        if len(otbl) == 0:
            otbl = Table(row)
        else:
            otbl.add_row(row)
    otbl.write(cTableFile,format='ipac',overwrite=False)
    if prt:
        print('mC2020T: wrote table size',len(otbl),'to',cTableFile,'query table size was',len(tbl))
    return True

def makeGaiaTable(destPath,coaddID, prt=False):
    # make a Gaia table in the destPath directory if it doesn't exist
    gTableFile = destPath + coaddID + '_gaia.tbl'
    if prt:
        print('makeGaiaTable',gTableFile)
    if os.path.exists(gTableFile):
        if prt:
            print(gTableFile,'exists')
        return True
    from astropy.coordinates import SkyCoord
    from astroquery.ipac.irsa import Irsa
    import astropy.units as u
    ctrRA, ctrDec = tileCtrPos(coaddID)
    skyCoord = SkyCoord(ctrRA, ctrDec, unit='deg', frame='icrs')
    table = Irsa.query_region(catalog='gaia_dr3_source',coordinates=skyCoord,spatial='Box',width=1.8*u.deg,
                              columns='ra,dec,pmra,pmdec,pm,ruwe,parallax,parallax_error')
    if not table:
        print('mGT table query failed')
        return False
    table.write(gTableFile,format='ipac',overwrite=False)
    if prt:
        print('mGT: wrote table size',len(table),'to',gTableFile)
    return True

def today(time_format='%Y%m%d'):
    from datetime import datetime
    return datetime.now().strftime(time_format)

def isInFile(fileName,string):
    with open(fileName,'r') as f:
        for line in f:
            if line.startswith(string):
                return True
    return False

def next1(fileName,loRA,hiRA):
    # finds the next coaddID after the most recently processed coaddID
    if not os.path.exists(fileName):
        print(fileName,'doesnt exist')
        exit()
    coaddID = loadCoaddID(fileName)
    from astropy.io import fits
    with fits.open('allsky-atlas.fits') as hdul:
        table = hdul[1].data
        cnt = -1
        for row in table:
            # only consider tile center RA positions within the range
            if row['ra'] < loRA or row['ra'] >= hiRA:
                continue
            if cnt >= 0:
                cnt += 1
                print('Next coaddID',row['coadd_id'])
                return row['coadd_id']
            if row['coadd_id'].find(coaddID) == 0:
                cnt = 0
    return 'NA'

def loadCoaddID(fileName):
    # return the first line in coaddID.txt
    with open(fileName,'r') as f:
        fields = f.readline().split()
        return fields[0]

def saveCoaddID(fileName,coaddID,comment=None):
    # read all of the previous entries
    with open(fileName,'r') as fin:
        lines = fin.readlines()
    fin.close()
    # check for duplicate
    for line in lines:
        if line.startswith(coaddID):
            print('This coaddID exists in',fileName)
            return False
    # append the today time stamp
    line = coaddID + ' ' + today()
    # and an optional comment
    if comment != None:
        line = line + ' ' + comment
    with open(fileName,'w') as fout:
        print(line,file=fout)
        # write the previous entries
        for line in lines:
            print(line.strip(),file=fout)
    fout.close()
    print(len(lines),'coaddIDs processed')
    return True

# used to highlight important information on the terminal window
class bcolors:
    blue = '\033[92m'
    red = '\033[91m'
    reset = '\033[0m'

class timeScaling:
    dpy = 365.25
    hpy = 24 * dpy # hours per year

# defines the columns of the PM fit catalog.
def fitResult():
    result = {'ra':0.,'dec':0.,'pmRA':0.,'epmRA':0.,'pmDec':0.,'epmDec':0.,'rchi2':-1.,
            'w1mag':-1.,'w2mag':-1.,'fitType':-1,'npts':0,'fracflux':0.,'nIter':0,
            'w1fluxVar':-1.,'w2fluxVar':-1.,'rchi2Pt':-1.,'flag_info':0,
            'version':0,'nComoCan':-1,'coaddID':'NA',
            'quality':-1, # 1 = maybe moving, 2 = probably moving, 3 = definitely moving
            'sf':0., 'sfPt':0. # survival fractions
            }
    return result

def spType(w1,w2):
    # Approximate spectral type using W1mag - W2mag
    # Carnero Rosell  https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.5301C/abstract
    w1w2 = [-0.01, # NA
        0.0,0.05,0.09,0.13,0.15,0.17,0.18,0.20,0.21,0.22, # M types
        0.23,0.24,0.25,0.26,0.30,0.34,0.36,0.40,0.46,0.52, # L types
        0.59,0.75,0.92,1.13,1.35,1.62,1.90,2.21,2.51,2.83, # T types
        3. # Y types
            ]
    sptyp = ['',
            'M0','M1','M2','M3','M4','M5','M6','M7','M8','M9',
            'L0','L1','L2','L3','L4','L5','L6','L7','L8','L9',
            'T0','T1','T2','T3','T4','T5','T6','T7','T8','T9',
            'Y']
    dw12 = w1 - w2
    indx = (np.abs(w1w2 - dw12)).argmin()
    return sptyp[indx], indx

def photDistWISE(w1, w2, prt=False):
    # "THE HAWAII INFRARED PARALLAX PROGRAM. I. ULTRACOOL BINARIES AND THE L/T TRANSITION"
    # The Astrophysical Journal Supplement Series, 201:19 (84pp), 2012 August 
    # Adapted from code by Tom Bickle
    espt = 2 # assumed error on the spectral type
    samp = 10000 # number of samples for estimating the distance uncertainty
    _, spti = spType(w1,w2) # get the spt index, ignore bdt
    sptn = np.random.normal(spti,espt,samp)
    # estimate distance using W1 mag
    W1_Mabsn = (7.14765 * (sptn**0.0)) + (3.55395e-1 * (sptn**1.0)) - (4.38105e-3 * (sptn**2.0)) - (3.33944e-4 * (sptn**3.0)) + (1.58040e-5 * (sptn**4.0))
    W1_Mabs = np.average(W1_Mabsn)
    W1_eMabs = np.std(W1_Mabsn)
    W1_mod = w1 - W1_Mabs
    ew1 = 0.1 # assumed error on w1mag
    W1_emod = np.sqrt(ew1*ew1 + W1_eMabs*W1_eMabs)
    W1_d = 10.0*(10.0**(W1_mod/5.0))
    W1_ed = 4.6052*W1_emod*np.sqrt(10.0**(0.4*W1_mod))
    # estimate distance using W2 mag
    W2_Mabsn = (7.46564 * (sptn**0.0)) + (1.92354e-1 * (sptn**1.0)) + (1.14325e-2 * (sptn**2.0)) - (8.81973e-4 * (sptn**3.0)) + (1.78555e-5 * (sptn**4.0))
    W2_Mabs = np.average(W2_Mabsn)
    W2_eMabs = np.std(W2_Mabsn)
    W2_mod = w2 - W2_Mabs
    ew2 = 0.1 # assumed error on w2mag
    W2_emod = np.sqrt(ew1*ew2 + W2_eMabs*W2_eMabs)
    W2_d = 10.0*(10.0**(W2_mod/5.0))
    W2_ed = 4.6052*W2_emod*np.sqrt(10.0**(0.4*W2_mod))
    if prt:
        print('photDistWISE:',end=' ')
        print('W1 {:.0f} \u00B1 {:.0f} pc'.format(W1_d,W1_ed),end=' ')
        print('W2 {:.0f} \u00B1 {:.0f} pc'.format(W2_d,W2_ed))
    return W1_d, W1_ed, W2_d, W2_ed

def photDistJK(w1,w2,jmag,ejmag,kmag,ekmag,prt=False):
    # ala photDistWISE
    espt = 2 # assumed error on the spectral type
    samp = 10000 # number of samples for estimating the distance uncertainty
    _, spti = spType(w1,w2) # get the spt index, ignore bdt
    sptn = np.random.normal(spti,espt,samp)
    j_d = 0.
    j_ed = 1.
    k_d = 0.
    k_ed = 1.
    if jmag > 0:
        j_Mabsn = (-9.67994 * (sptn**0.0)) + (8.16362 * (sptn**1.0)) - (1.33053 * (sptn**2.0)) + (1.11715e-1 * (sptn**3.0)) - (4.82973e-3 * (sptn**4.0)) + (1.00820e-4 * (sptn**5.0)) - (7.84614e-7 * (sptn**6.0))
        j_Mabs = np.average(j_Mabsn)
        j_eMabs = np.std(j_Mabsn)
        j_mod = jmag - j_Mabs
        j_emod = np.sqrt(ejmag*ejmag + j_eMabs*j_eMabs)
        j_d = 10.0*(10.0**(j_mod/5.0))
        j_ed = 4.6052*j_emod*np.sqrt(10.0**(0.4*j_mod))
    if kmag > 0:
        k_Mabsn = (1.10114e1 * (sptn**0.0)) - (8.67471e-1 * (sptn**1.0)) + (1.34163e-1 * (sptn**2.0)) - (6.42118e-3 * (sptn**3.0)) + (1.06693e-4 * (sptn**4.0))
        k_Mabs = np.average(k_Mabsn)
        k_eMabs = np.std(k_Mabsn)
        k_mod = kmag - k_Mabs
        k_emod = np.sqrt(ekmag*ekmag + k_eMabs*k_eMabs)
        k_d = 10.0*(10.0**(k_mod/5.0))
        k_ed = 4.6052*k_emod*np.sqrt(10.0**(0.4*k_mod))
    if prt:
        print('photDistJK:',end=' ')
        if j_d > 0:
            print('J {:.0f} \u00B1 {:.0f} pc'.format(j_d,j_ed),end=' ')
        if k_d > 0:
            print('K {:.0f} \u00B1 {:.0f} pc'.format(k_d,k_ed),end=' ')
        print()
    return j_d, j_ed, k_d, k_ed

def chkFluxVar(result,catl):
    # finds the weighted average flux in each band, then counts the number of catalog entries
    # that lie outside 1.5 sigma from the average. The fraction should be ~13% if the flux has
    # a Gaussian distribution. Larger values indicate some level of long term variability.
    result['w1fluxVar'] = -1
    result['w2fluxVar'] = -1
    for band in range(1,3):
        fluxs = []
        wghts = []
        for row in catl:
            if row['band'] == band and row['nm'] >= 0:
                fluxs.append(row['flux'])
                wghts.append(1./row['dflux'])
        if len(fluxs) < 4:
            continue
        ave = np.average(fluxs,weights=wghts)
        cnts = 0
        for ipt in range(len(fluxs)):
            pull = wghts[ipt] * abs(fluxs[ipt] - ave)
            if pull > 1.5:
                cnts += 1
        frac = cnts / len(fluxs)
        if band == 1:
            result['w1fluxVar'] = frac
        else:
            result['w2fluxVar'] = frac
        # expect 13% of points > 1.5 sigma for a Gaussian distribution
    return result

def fitQuality(result,catl,ptFit=False):
    # return the x+y rchi2 and return the indx of the highest residual catalog item
    pix2mas = 1000.*2.75
    yr = []
    xx = []
    yy = []
    for row in catl:
        # ignore points not in the original cluster
        if row['nm'] == -2:
            continue
        # ignore points rejected during the PM fit
        if not ptFit and row['nm'] == -4:
            continue
        yr.append(row['year'])
        xx.append(row['x'])
        yy.append(row['y'])
    midYear = np.average(yr)
    midx = np.average(xx)
    midy = np.average(yy)
    # project result to the expected positions and calculate rchi2
    cnt = 0
    chi2s = 0.
    if ptFit:
        npar = 2
        slpx = 0.
        slpy = 0.
    else:
        npar = 4
        slpx = result['pmDec']/pix2mas
        slpy = -result['pmRA']/pix2mas
    bigChi = 0.
    imbig = -1
    for ipt, row in enumerate(catl):
        if row['nm'] < 0 or row['flux'] < 0.:
            continue
        dyr = row['year'] - midYear
        dxx = (midx + slpx * dyr - row['x'])/row['dx']
        dyy = (midy + slpy * dyr - row['y'])/row['dy']
        # Dec 9, 2025 added 0.5
        chi2 = 0.5 * (dxx*dxx + dyy*dyy)
        if chi2 > bigChi:
            bigChi = chi2
            imbig = ipt
        chi2s += chi2
        cnt += 1
    ndof = cnt - npar
    if ndof < 3:
        return 999.,0.,0
    rchi2 = chi2s / ndof
#    print('fQ:',ptFit,rchi2,ndof)
    return rchi2, ndof, imbig

def getCatFiles(tileCache,coaddID):
    # copy the C2020 tables to a local subdirectory of the current path.
    # First see if they already exist
    catFiles = tileCache+coaddID+'_cat.tbl'
    import glob
    files = glob.glob(catFiles)
    if len(files) > 0:
        print('gCF: files already exist in',catFiles)
        return 0
    baseURL = 'https://portal.nersc.gov/project/cosmo/data/CatWISE/2020/'
    srcDir = baseURL + coaddID[0:3] + '/'
    print('gCF: srcDir',srcDir)
    import ssl
    myssl = ssl.create_default_context()
    myssl.check_hostname=False
    myssl.verify_mode=ssl.CERT_NONE
    import os
    from urllib.request import urlopen
    response = None
    try:
        response = urlopen(srcDir,context=myssl)
    except Exception:
        print('gCF: no response from nersc')
        return 1
    string = response.read().decode('utf-8')
    # insert new line characters
    string.replace('/td>','/td>\n')
    # and split into lines
    splits = string.splitlines()
    cnt = 0
    for split in splits:
        # extract the file name with pattern <coaddID>....gz
        icid = split.find(coaddID)
        if icid < 0:
            continue
        igz = split.find('.gz')
        if igz < 0:
            continue
        igz += 3
        if split.find('rej') > 0:
            continue
        fileName = split[icid:igz]
        srcFile = srcDir+fileName
        dstFile = tileCache+coaddID+'_cat.tbl.gz'
        # copy the file
        g = urlopen(srcFile,context=myssl)
        with open(dstFile,'b+w') as f:
            f.write(g.read())
        f.close()
        print('gCF: copied',fileName)
        cnt += 1
    os.system('gunzip '+tileCache+'*.gz')
    return 0

def getUntimeFits(tileCache,coaddID):
    # copy all the unTimely catalog fits files to a local subdirectory.
    # See if the files are already slored locally
    fitsFiles = tileCache+coaddID+'*.fits'
    import glob
    files = glob.glob(fitsFiles)
    if len(files) > 4:
        return
    import ssl
    import os
    # get a list of zipped unTimely files
    baseURLs = ['https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/',
                'https://portal.nersc.gov/project/cosmo/temp/ameisner/crowdsource_tr_neo8/']
    for baseURL in baseURLs:
        srcDir = baseURL + coaddID[0:3] + '/' + coaddID + '/'
        print('srcDir',srcDir)
        myssl = ssl.create_default_context()
        myssl.check_hostname=False
        myssl.verify_mode=ssl.CERT_NONE
        from urllib.request import urlopen
        response = urlopen(srcDir,context=myssl)
        string = response.read().decode('utf-8')
        # insert new line characters
        string.replace('/td>','/td>\n')
        # and split into lines
        splits = string.splitlines()
        cnt = 0
        for split in splits:
            # extract the file name with pattern <coaddID>....gz
            icid = split.find(coaddID)
            if icid < 0:
                continue
            igz = split.find('.gz')
            if igz < 0:
                continue
            igz += 3
            fileName = split[icid:igz]
            srcFile = srcDir+fileName
            dstFile = tileCache+fileName
            if os.path.exists(dstFile):
                print(dstFile,'already exists')
                continue
            # copy the file
            g = urlopen(srcFile,context=myssl)
            with open(dstFile,'b+w') as f:
                f.write(g.read())
            f.close()
            cnt += 1
    os.system('gunzip '+tileCache+'*.gz')

def flux2Mag(flux):
    if flux > 0.:
        return 22.5-2.5*np.log10(flux)
    return -1.

def mag2Flux(mag):
    if mag < 22.5:
        exp = (22.5-mag)/2.5
        return np.power(10,exp)
    return -1.

def printResult(caller,result,printLevel):
    # print the fit result to the terminal window. The caller variable is a string passed
    # by the calling module. More detailed information is displayed when printLevel > 0
    if result['rchi2'] < 0:
        print(caller,'No fit')
        return
    pmSig = -1.
    if result['epmRA'] > 0 and result['epmDec']> 0:
        sigRA = abs(result['pmRA'])/result['epmRA']
        sigDec = abs(result['pmDec'])/result['epmDec']
        pmSig = np.sqrt(sigRA*sigRA+sigDec*sigDec)
    caller = caller + 'v{}:'.format(result['version'])
    print(caller,end=' ')
    if printLevel > 1:
        print('{:.6f} {:.6f}'.format(result['ra'],result['dec']),end=' ')
    print('PM {:.0f} \u00B1 {:.0f} {:.0f} \u00B1 {:.0f}'.format(result['pmRA'],result['epmRA'],
                                            result['pmDec'],result['epmDec']),end=' ')
    fitStr = 'NA'
    if result['fitType'] == 1:
        fitStr = 'w1'
    elif result['fitType'] == 2:
        fitStr = 'w2'
    elif result['fitType'] == 3:
        fitStr = 'w1w2'
    fitStr += ' fit'
    if result['rchi2'] > 1.5 or pmSig < 3.:
        print(fitStr,'pmSig {:.2f}'.format(pmSig),end=' ')
        print('rchi2 {:.2f}'.format(result['rchi2']),end=' ')
    else:
        print(fitStr,bcolors.red+'pmSig {:.2f}'.format(pmSig),end=' ')
        print('rchi2 {:.2f}'.format(result['rchi2'])+bcolors.reset,end=' ')       
    print('rchi2Pt {:.2f}'.format(result['rchi2Pt']),end=' ')
    if result['npts'] < 7:
        print(bcolors.red+'nEpochs {}/{}'.format(result['npts'],result['npts']+result['nIter'])+bcolors.reset,end=' ')
    else:
        print('nEpochs {}/{}'.format(result['npts'],result['npts']+result['nIter']),end=' ')
    print()
    print('\tmags {:.1f} {:.1f}'.format(result['w1mag'],result['w2mag']),end=' ')
    print('SF {:.2f} SFPt {:.2f}'.format(result['sf'],result['sfPt']),end=' ')
    if printLevel > 0:
        if result['fracflux'] < 0.8:
            print(bcolors.red+'fracflux {:.2f}'.format(result['fracflux'])+bcolors.reset,end=' ')
        else:
            print('fracflux {:.2f}'.format(result['fracflux']),end=' ')
        print(result['coaddID'],end=' ')
        if result['w1fluxVar'] > 0:
            print('w1fluxVar {:.2f}'.format(result['w1fluxVar']),end=' ')
        if result['w2fluxVar'] > 0:
            print('w2fluxVar {:.2f}'.format(result['w2fluxVar']),end=' ')
        if result['flag_info'] > 0:
            flag = result['flag_info']
            print(caller,'flag_info bits:')
            # parse the flags_info bits using https://catalog.unwise.me/catalogs.html#flags_info
            info = ['In PSF bright star','HyperLeda galaxy','Big object','Very bright star','Saturated','Nebulosity']
            # add w1Var and w2Var
            for bit in range(6):
                shft = np.right_shift(flag,bit)
                mask = np.bitwise_and(shft,1)
                if mask == 1:
                    print(bcolors.red+info[bit]+bcolors.reset,end=' ')
        print('quality',result['quality'])
        print()

def profileHist(x,y,nbins,xmin,xmax,ax=None,label=None,spread=False):
    # displays a profile histogram of (x,y) with nbins between xmin and xmax on the
    # specified axis
    if len(x) != len(y):
        print('x and y sizes are different')
        return 0
    from scipy import stats
    means, binEdges, _ = stats.binned_statistic(x, y, bins=nbins, range=(xmin,xmax), statistic='mean')
    binCtrs = (binEdges[:-1] + binEdges[1:])/2.
    stdevs, binEdges, _ = stats.binned_statistic(x, y, bins=nbins, range=(xmin,xmax), statistic='std')
    cnts, binEdges, _ = stats.binned_statistic(x, y, bins=nbins, range=(xmin,xmax), statistic='count')
    if not spread:
        stdevs /= np.sqrt(cnts)
    if ax != None:
        ax.errorbar(x=binCtrs, y=means, yerr=stdevs, linestyle='none', marker='.',label=label)
    return binCtrs, means, stdevs

def showStats(vals,ax):
    # displays the mean and std deviation on the specified axis
    mean = np.average(vals)
    stdev = np.std(vals)
    # standard error on the mean
    sem = stdev / np.sqrt(len(vals))
    stats = '{:.0f} entries\nmean:  {:.3f}\nstdev: {:.3f}\nSEM:   {:.3f}'.format(len(vals),mean,stdev,sem)
    ax.text(0.01,0.99,stats,
                  fontsize='small',transform=ax.transAxes,verticalalignment='top')

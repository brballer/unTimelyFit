import numpy as np

'''''
Dec 1, 2025: comment out ptfit (use fitQuality in fitAll.py)
'''''

def linfit(x,y,wgt):
    result = [-1,-1,-1.,-1.,-1]
    if len(x) < 3:
        return result
    sum = np.double(0.)
    sumx = np.double(0.)
    sumy = np.double(0.)
    sumxy = np.double(0.)
    sumx2 = np.double(0.)
    sumy2 = np.double(0.)
    resids = np.zeros(len(x),np.float64)
    cnt = 0.
    for ipt in range(len(x)):
        if wgt[ipt] <= 0.:
            continue
        sum += wgt[ipt]
        sumx += wgt[ipt]*x[ipt]
        sumy += wgt[ipt]*y[ipt]
        sumx2 += wgt[ipt]*x[ipt]*x[ipt]
        sumy2 += wgt[ipt]*y[ipt]*y[ipt]
        sumxy += wgt[ipt]*x[ipt]*y[ipt]
        cnt += 1
    if cnt < 3:
        result[0] = 999.
        return result
    delta = sum*sumx2-sumx*sumx
    if delta == 0.:
        result[0] = 998.
        return result
    # A is the intercept that we don't care about but need it to find the variance
    A = (sumx2*sumy-sumx*sumxy)/delta
    # B is the slope
    B = (sumxy*sum-sumx*sumy)/delta
    ndof = cnt-2
    varnce = (sumy2+A*A*sum+B*B*sumx2-2*(A*sumy+B*sumxy-A*B*sumx))/ndof
    chisum = np.double(0.)
    for ipt in range(len(x)):
        if wgt[ipt] == 0.:
            continue
        resid = y[ipt]-A-B*x[ipt]
        chisum += resid*resid*wgt[ipt]
        resids[ipt] = resid
    BErr = 1000
    rchi2 = -1.
    if varnce > 0:
        arg = varnce*sum/delta
        if arg < 0.:
            result[0] = 997.
            return result
        BErr = np.sqrt(arg)
        rchi2 = chisum/ndof
    return rchi2, ndof, B, BErr, resids

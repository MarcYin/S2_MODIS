import sys
sys.path.insert(0,'python')
import numpy.ma as ma
from collections import Counter
import cPickle as pkl
from scipy import optimize
from functools import partial
import scipy.ndimage as ndimage
import numpy as np
import scipy
from scipy import signal
import scipy.stats
from cloud import *
from fastRWpkl import *
import numpy 
from numpy import clip, where
from scipy.ndimage.morphology import *
import xml.etree.cElementTree as ET
import multiprocessing
from get_r import *
from sklearn import linear_model

keys = 'B02', 'B03','B04','B08','B8A','B11','B12'
bands = [2,3,4,8,13,11,12]

def gaussian(xwin, ywin, xstd, ystd, angle, norm = True):
    win = max(xwin, ywin)
    winx = win*2**0.5
    winy = win*2**0.5
        
    xgaus = signal.gaussian(winx, xstd)
    ygaus = signal.gaussian(winy, ystd)
    gaus  = np.outer(xgaus, ygaus)
    r_gaus = scipy.ndimage.interpolation.rotate(gaus, angle, reshape=True)
    center = np.array(r_gaus.shape)/2
    cgaus = r_gaus[center[0]-xwin/2: center[0]+xwin/2, center[1]-ywin/2:center[1]+ywin/2]
    if norm:
        return cgaus/cgaus.sum()
    else:
        return cgaus

def _psf(sent, sinds, mod, minds,band, psf):    
    xstd,ystd,angle, xs, ys = psf
    xwin,ywin = 150, 150
    
    to_regression =[]          
    cx = sinds[0]
    cy = sinds[1]
    mx = minds[0]
    my = minds[1]
    
    gaus = gaussian(xwin,ywin,xstd,ystd,angle,False)                              
    ker = gaus/(gaus.sum())

    s = signal.fftconvolve(sent, ker, mode='same')
    
    vld_x = ((cx+xs)>xwin/2)&((cx+xs)<10000-xwin/2)
    vld_y = ((cy+ys)>ywin/2)&((cy+ys)<10000-ywin/2)
    vld = vld_x&vld_y
    
    
    indx,indy = np.round((cx+xs)[vld]).astype(int), np.round((cy+ys)[vld]).astype(int)
    vals = s[indx,indy]
    brdf = mod[mx[vld], my[vld]]
    mask = (brdf>0)&(brdf<1)&(vals>0)&(vals<1)
    if sum(mask) ==0:
        print 'Too much cloud again to affect the convolve results'
        return 10000
    else:
        dif = vals[mask] - brdf[mask]
        inliers = (dif>(np.nanmean(dif)-3*np.nanstd(dif)))&(dif<(np.nanmean(dif)+3*np.nanstd(dif)))&(dif<0.2)
        
        #global vals; global mask; global brdf; global inliers
        #x,y = ransaclin(vals[mask][inliers], brdf[mask][inliers])

        vx = indx[mask][inliers]
        vy = indy[mask][inliers]
        m = vals[mask][inliers]#y.ravel()
        s = brdf[mask][inliers]#x.ravel()

    
    return [band,vx,vy, m,s]

def ransaclin(x,y):
    y, x = y.reshape((len(y),1)), x.reshape((len(x),1))
    model = linear_model.LinearRegression()
    model.fit(x, y)

    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
    model_ransac.fit(x, y)
    inlier_mask = model_ransac.inlier_mask_
    return x[inlier_mask], y[inlier_mask]

def ScaleExtent(data, shape): # used for unifine different array,

    re = int(shape[0]/(data.shape[0]))

    a = np.repeat(np.repeat(data, re, axis = 1), re, axis =0)
    
    if (re*(data.shape[0])-shape[0]) != 0:
        extended = np.zeros(shape)
        extended[:re*(data.shape[0]),:re*(data.shape[0])] = a
        extended[re*(data.shape[0]):,re*(data.shape[0]):] = a[re*(data.shape[0])-shape[0]:, re*(data.shape[0])-shape[0]]
        return extended
    else:
        return a

def op(ind,  args=None ):
    fpath, sentm, brdfs, sinds, minds, psfs = args
    Sent = gdal_read(bands[ind], fpath)[keys[ind]]
    sent = ScaleExtent(Sent, (10980,10980)) 
    sent[sentm]= np.nanmean(sent[~sentm])
    sent[np.isnan(sent)] = np.nanmean(sent[~sentm])
    if ind<4:
        brdfs[ind][brdfs[ind].mask] = np.nan                
        to_regs = _psf(sent,sinds, brdfs[ind]*0.001, minds, bands[ind], psfs[ind])

    else:
        brdfs[ind-1][brdfs[ind-1].mask] = np.nan
        to_regs = _psf(sent,sinds, brdfs[ind-1]*0.001, minds, bands[ind], psfs[ind])
    return to_regs
        
def applied(lat, lon,fpath, mfile, psfs, pr=1):
    
    sentm = get_cloud_mask(fpath)
    doy = int(mfile[0].split('.')[1][-3:])
    pos = fpath.split('/')[6]+fpath.split('/')[7]+fpath.split('/')[8]
    if pr:
        print 'DOY: ', doy,'\n', 'Location: ', pos
        print 'Cloud proportion: ', sentm.sum()/(10980.*10980.)
    
    minds, sinds = get_coords(lat,lon) 
    modis_filenames = gdal.Open(mfile[0]).GetSubDatasets()
    modisQA = gdal.Open(mfile[1]).GetSubDatasets()

    brdfs = get_rs(modisQA, modis_filenames, fpath)
    args = fpath, sentm, brdfs, sinds, minds, psfs
    par = partial(op, args=args)
    pool = multiprocessing.Pool(processes = 7)
    retval = pool.map(par, range(7))
    pool.close()
    pool.join()
    return retval
    

   
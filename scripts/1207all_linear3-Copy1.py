import multiprocessing
import sys
sys.path.insert(0,'python')
from fastRWpkl import *
import numpy.ma as ma
from readSent import *
from collections import Counter
import cPickle as pkl
from scipy import optimize
from functools import partial
import scipy.ndimage as ndimage
import numpy as np
import scipy
from scipy import signal
import scipy.stats


fhead = 'data/50SMG20165100'

def gaussian(xwin, ywin, xstd, ystd, angle, norm = True):
    win = max(xwin, ywin)
    winx = win*2**0.5
    winy = win*2**0.5
    
    xstd = xstd*2**0.5
    ystd = ystd*2**0.5
        
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

def cost(p, process):    
    xstd,ystd,angle = p
    sx, sy = -22.544921, 57.609055
    xwin,ywin = 100,100
    
    i, j = process
    print 'patch %s%s'%(i,j)
    
    s1 = slice((i*1000),(i+1)*1000)
    s2 = slice((j*1000),(j+1)*1000)
    
    #modis_cut = modis_sent[s1,s2]*0.001
    sen_cut = Sent[s1,s2]
    
    in_patch = (centers[0]>=i*1000)&(centers[0]<i*1000+1000)&(centers[1]>=j*1000)&(centers[1]<j*1000+1000)
    patchx, patchy = (centers[0])[in_patch], (centers[1])[in_patch]
    to_regression =[]
    for ii,cx in enumerate(patchx):
        cy = patchy[ii]
        minx = cx+sx - 50
        maxx = cx+sx + 50
        miny = cy+sy - 50
        maxy = cy+sy + 50

        to_conv = Sent[max(0+i*1000, minx): min(1000+i*1000, maxx), max(0+j*1000, miny): min(1000+j*1000, maxy)]
        brdf = modis_sent[cx,cy]*0.001

        if (to_conv.shape[0]==100) & (to_conv.shape[1]==100) & (brdf!=np.nan):
            To_conv = to_conv
            Brdf = brdf
            nanval = np.where(~((To_conv < 1)&(To_conv > 0)))
            To_conv[nanval[0], nanval[1]] = np.nanmean(To_conv)

            
            gaus = gaussian(xwin,ywin,xstd,ystd,angle,False)                              
            ker = gaus/(gaus.sum())

            s = signal.fftconvolve(To_conv, ker, mode='valid')
            to_regression.append([s[0][0], Brdf])
    x,y = np.array(to_regression).T
    nanm = (np.isnan(x))|(np.isnan(y))
    r = scipy.stats.linregress(x[~nanm],y[~nanm])
    costs = abs(1-r.rvalue) 
    print 'costs:', costs, 'rvalue: ', r.rvalue, 'slop: ', r.slope, '\n', 'parameters: ', p,'\n'

    return costs
                
def solve(process):
    solved = []
    for i,j in process:
        p = np.array([140, 240,45])
        bound = np.array([(14., 400.),(8.,400.),(20.,160.)])
        psolve = optimize.fmin_l_bfgs_b(cost,p,approx_grad=True,iprint=-1,args=(([i,j],)),bounds=bound)
        solved.append([i,j,psolve])
        print 'solved one: ', psolve, '\n'
    return solved

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

Sindex = parallel_rw_pkl(None, 'Sindex', 'r')
centers = Sindex
data = parallel_rw_pkl(None, 'inter_sent%i'%2, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%2, 'r')
sent = readfile([4,],fhead)['B04']
cm = parallel_rw_pkl(None, '0510diacm', 'r')
sent[cm] = np.nan

print 'finshed reading data'
data[mask]=np.nan
modis_sent = np.array(data)
Sent = sent

patches = np.array(zip(np.mgrid[0:10,0:10][0].ravel(), np.mgrid[0:10,0:10][1].ravel()))
pros = np.array(np.array_split(patches, 16))

par = partial(solve)
pool = multiprocessing.Pool(processes=16)
data = pool.map(par, pros)
pool.close()
pool.join()
parallel_rw_pkl(data, 'bgopb4_modis', 'w')

print 'lol finished bgop b4!!!!!'                                    

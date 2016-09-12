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

def applied(process):
    store = []
    for i, j in process:
        print i,j

        s1 = slice((i*1000),(i+1)*1000)
        s2 = slice((j*1000),(j+1)*1000)
        modis_cut = modis_sent[s1,s2]*0.001
        Stm_cut = Stm[s1,s2]

        sen_cut = Sent[s1,s2]
        ulist = np.unique(Stm_cut)
        
        
        for ii,u in enumerate(ulist):
            mask = (Stm_cut == u)
            if (mask.sum() >= 1800) and (mask.sum() < 4000):
                xmin = np.where(mask)[0].min()
                xmax = np.where(mask)[0].max()
                ymin = np.where(mask)[1].min()
                ymax = np.where(mask)[1].max()
                center0 = np.array([(xmax+xmin)/2, (ymax+ymin)/2])
                
                xs,ys = -21.8540896067 , 60.2630764605

                minx = center0[0]+xs - 50
                maxx = center0[0]+xs + 50
                miny = center0[1]+ys - 50
                maxy = center0[1]+ys + 50

                #offset = -0.040317007446953879
                to_conv = sen_cut[max(0, minx): min(999, maxx), max(0, miny): min(999, maxy)]
                brdf = modis_cut[center0[0],center0[1]]
                
                
                if (to_conv.shape[0]>=100) & (to_conv.shape[1]>=100) & (brdf!=np.nan):
                    To_conv = to_conv
                    Brdf = brdf
                    nanval = np.where(~((To_conv < 1)&(To_conv > 0)))
                    To_conv[nanval[0], nanval[1]] = np.nanmean(To_conv)
                   
                    xstd, ystd, angle = 29.696923, 339.96191637, 45.8861022468

                    xwin,ywin = 100,100
                    gaus = gaussian(xwin,ywin,xstd,ystd,angle,False)
                    if gaus.sum() <= 0:
                        print 'invalid gauss: %s'%([xstd, ystd, angle])
                        pass
                    else:
                        kernel = gaus/(gaus.sum())    
                        s = signal.fftconvolve(to_conv, kernel, mode='valid')
                        
                        store.append([s[0][0], brdf])
                        print store[-1]
                else:
                    pass

            else:
                pass
    
    return store

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


data = parallel_rw_pkl(None, 'inter_sent%i'%3, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%3, 'r')
sent = readfile([8,],fhead)['B08']
cm = parallel_rw_pkl(None, '0510diacm', 'r')
sent[cm] = np.nan
stm = parallel_rw_pkl(None, 'std_m', 'r')
print 'finshed reading data'
data[mask]=np.nan
modis_sent = np.array(data)
Sent = sent
Stm = stm
patches = [(0, 0), (0, 1), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (3, 1), (3, 2), (3, 3), (3, 4), (3, 6), (3, 8), (3, 9), (4, 0), (4, 1), (4, 2), (5, 0), (5, 1), (5, 2), (5, 8), (6, 0), (6, 1), (6, 3), (7, 0), (7, 1), (7, 3), (8, 1), (8, 4), (8, 5), (9, 0), (9, 1), (9, 4)]

pros = np.array(np.array_split(patches, 16))

par = partial(applied)
pool = multiprocessing.Pool(processes=16)
data = pool.map(par, pros)
pool.close()
pool.join()
parallel_rw_pkl(data, 'eropb8_modis', 'w')

print 'lol finished erop b8!!!!!'    

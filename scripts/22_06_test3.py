from scipy import signal
import numpy as np
import scipy
from scipy import signal
from scipy.ndimage.morphology import binary_dilation as bd
import scipy.ndimage as ndimage

from functools import partial
import multiprocessing
import sys
sys.path.insert(0,'python')
from fast_rw import *
import numpy.ma as ma



def gaussian(xwin, ywin, xstd, ystd, angle, norm = True):
    win = max(xwin, ywin)
    winx = win*2**0.5
    winy = win*2**0.5
    xgaus = signal.gaussian(winx, xstd)
    ygaus = signal.gaussian(winy, ystd)
    gaus  = np.outer(xgaus/(xgaus.sum()), ygaus/(ygaus.sum()))
    r_gaus = scipy.ndimage.interpolation.rotate(gaus, angle, reshape=True)
    center = np.array(r_gaus.shape)/2
    cgaus = r_gaus[center[0]-xwin/2.: center[0]+xwin/2., center[1]-ywin/2.:center[1]+ywin/2.]
    if norm:
        return cgaus/cgaus.sum()
    else:
        return cgaus
#Gaussian one
def shift(process, modis_sent = None, Sent = None, star = None, end =None, mode='gaus'): 
    wcosts = []
    i, j = process
    s1 = slice((i*1000),(i+1)*1000)
    s2 = slice((j*1000),(j+1)*1000)
    modis_cut = modis_sent[s1,s2]*0.001
    modis_cut.mask[np.isnan(modis_cut)] = True
    modis_mask = np.zeros_like(modis_cut).astype(int)

    sen_cut = Sent[s1,s2]

    ulist = np.sort(np.unique(modis_cut[~modis_cut.mask]))   
    struct = ndimage.generate_binary_structure(2, 2)

    for ii,u in enumerate(ulist):
        mask = (modis_cut == u)
        if (mask.sum() >= 2000) and (mask.sum() < 3000):
            print u
            dia_mask = bd(mask, structure=struct, iterations=200)
            data = sen_cut.copy()
            data[~dia_mask] = 0
            to_conv = data[min(np.where(dia_mask)[0]):max(np.where(dia_mask)[0])+1,\
                           min(np.where(dia_mask)[1]):max(np.where(dia_mask)[1])+1]
            if mode == 'mean':
                w = 1./(np.nansum(mask))
                k = np.zeros(mask.shape).astype('float')
                k[mask] = w
                conved = signal.fftconvolve(to_conv, k, mode='valid')
                dif = abs(conved-u)
                minm = np.nanmin(dif)
                x = np.where(dif == minm)[0][0]-np.ceil((conved.shape[0])/2.)
                y = np.where(dif == minm)[1][0]-np.ceil((conved.shape[1])/2.)
                wcost.append([x, y , u, np.nanmin(dif), np.where(mask)[0]+1000*i, np.where(mask)[1]+j*1000])

            elif mode == 'gaus':
                xmin = np.where(mask)[0].min()
                xmax = np.where(mask)[0].max()
                ymin = np.where(mask)[1].min()
                ymax = np.where(mask)[1].max()
                xwin = xmax - xmin
                ywin = ymax - ymin
                print 'wind:',xwin,ywin
                if (xwin <= 0) or (ywin<=0):
                    pass
                else:
                    cost = []
                    start = 100
                    if star == 0:
                        star +=1
                    if end == 0:
                        end +=1
                    for xstd in np.arange(star,end,5):
                        for ystd in np.arange(star,end,5):
                            if xstd <=ystd:
                                for angle in xrange(0,180, 5):
                                    gaus = gaussian(xwin,ywin,xstd,ystd,angle, False)
                                    smask = mask[xmin:xmax,ymin:ymax].copy()
                                    gaus[~smask]=0
                                    ker = gaus/(gaus.sum())
                                    conved = signal.fftconvolve(to_conv, ker, mode='valid')
                                    dif = abs(conved-u)
                                    if np.nanmin(dif)<start:
                                        minm = np.nanmin(dif)
                                        x = np.where(dif == minm)[0][0]-np.ceil((conved.shape[0])/2.)
                                        y = np.where(dif == minm)[1][0]-np.ceil((conved.shape[1])/2.)
                                        cost.append([xstd, ystd, angle, x, y , u, np.nanmin(dif), 
                                                     np.where(mask)[0]+1000*i, np.where(mask)[1]+j*1000])
                                        start = np.nanmin(dif)
                                        print 'Find One!!', start
                            else:
                                pass


                wcosts.append(cost[-1])
        else:
            pass
    return wcosts


data = parallel_rw_pkl(None, 'inter_sent%i'%2, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%2, 'r')
grid_z0 = ma.array(data, mask=mask)
b4 = parallel_rw_pkl(None, 'b4', 'r')

patches = np.array(zip(np.mgrid[0:10,0:10][0].ravel(), np.mgrid[0:10,0:10][1].ravel()))
ite = tuple(patches[40:60])

par = partial(shift, modis_sent=grid_z0, Sent = b4, star=1, end=100)
pool = multiprocessing.Pool(processes=20)
data = pool.map(par, ite)
pool.close()
pool.join()
data = np.array(data, dtype=object)
parallel_rw_pkl(data, '22_06_test3', 'w')

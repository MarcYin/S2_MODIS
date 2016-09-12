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
from collections import Counter


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
def shift(pat, mode='gaus'):
    mask,to_conv, brdf, indx, indy = pat

    if (mask.sum() >= 2000) and (mask.sum() < 3000):

        if mode == 'mean':
            w = 1./(np.nansum(mask))
            k = np.zeros(mask.shape).astype('float')
            k[mask] = w
            conved = signal.fftconvolve(to_conv, k, mode='valid')
            dif = abs(conved-u)
            minm = np.nanmin(dif)
            x = np.where(dif == minm)[0][0]-np.ceil((conved.shape[0])/2.)
            y = np.where(dif == minm)[1][0]-np.ceil((conved.shape[1])/2.)
            vals = conved[np.where(dif == minm)[0][0], np.where(dif == minm)[1][0]]
            return [x, y , brdf, vals, indx,indy]
        elif mode == 'gaus':

            xwin, ywin = mask.shape

            if (xwin <= 0) or (ywin<=0):
                pass
            else:
                cost = []
                start = 100
                #if star == 0:
                 #   star +=0.0001
                #if end == 0:
                 #   end +=0.0002
                for xstd in np.arange(1,50,1):
                    for ystd in np.arange(1,50,1):
                        if xstd <=ystd:
                            for angle in xrange(30,160, 2):
                                gaus = gaussian(xwin,ywin,xstd,ystd,angle, False)
                                gaus[~mask]=0
                                ker = gaus/(gaus.sum())
                                conved = signal.fftconvolve(to_conv, ker, mode='valid')
                                dif = abs(conved-brdf)
                                minm = np.nanmin(dif)
                                if minm<start:
                                    x = np.where(dif == minm)[0][0]-np.ceil((conved.shape[0])/2.)
                                    y = np.where(dif == minm)[1][0]-np.ceil((conved.shape[1])/2.)
                                    vals = conved[np.where(dif == minm)[0][0], np.where(dif == minm)[1][0]]
                                    cost.append([xstd, ystd, angle, x, y , brdf, vals, indx, indy])
                                    start = minm
                                    print 'Find One!!', start
                        else:
                            pass
                return cost[-1]

        else:
            pass


    else:
        pass

data = parallel_rw_pkl(None, 'inter_sent%i'%2, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%2, 'r')
stm= parallel_rw_pkl(None, 'std_m', 'r')

grid_z0 = ma.array(data, mask=mask)
b4 = parallel_rw_pkl(None, 'b4', 'r')
print 'finished reading data'
#inds = np.array_split(np.unique(stm),5)
#ite = inds[0]

Sent = b4
modis_sent = grid_z0
struct = ndimage.generate_binary_structure(2, 2)

def get_pixels(i):
    mask = (stm==i)
    brdf = Counter(modis_sent[mask]).most_common(1)[0][0]
    xmin = np.where(mask)[0].min()
    xmax = np.where(mask)[0].max()
    ymin = np.where(mask)[1].min()
    ymax = np.where(mask)[1].max()
    indx = np.where(mask)[0]
    indy = np.where(mask)[1]
    dia_mask = bd(mask, structure=struct, iterations=200)
    Sent[~dia_mask] = 0
    to_conv = Sent[min(np.where(dia_mask)[0]):max(np.where(dia_mask)[0])+1,\
                       min(np.where(dia_mask)[1]):max(np.where(dia_mask)[1])+1]
    mask = mask[xmin:xmax, ymin:ymax]

    return np.array([mask,to_conv, brdf, indx, indy], dtype = object)

randind = np.random.choice(np.unique(stm),1000)

for i,j in enumerate(np.array_split(randind, 25)):
    patches = (get_pixels(k) for k in j)
    par = partial(shift, mode='gaus')
    pool = multiprocessing.Pool(processes=45)
    data = pool.map(par, patches)
    pool.close()
    pool.join()
    data = np.array(data, dtype=object)
    parallel_rw_pkl(data, '2406test4_%i'%i, 'w')

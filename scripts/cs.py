from functools import partial
import multiprocessing
import sys
sys.path.insert(0,'python')
from fast_rw import *
import numpy.ma as ma
from scipy import signal
from scipy.ndimage.morphology import binary_dilation as bd
import scipy.ndimage as ndimage

def shift(process, modis_sent = None, Sent = None): 
    wcosts = []
    for i, j in process:

        s1 = slice((i*1000),(i+1)*1000)
        s2 = slice((j*1000),(j+1)*1000)
        modis_cut = modis_sent[s1,s2]*0.001
        #modis_cut = grid_z0*0.001
        modis_cut.mask[np.isnan(modis_cut)] = True
        #modis_mask = np.zeros_like(modis_cut).astype(int)

        #sen_cut = b4
        sen_cut = Sent[s1,s2]
        #sen_modis = np.zeros_like(modis_cut)

        ulist = np.sort(np.unique(modis_cut[~modis_cut.mask]))   
        struct1 = ndimage.generate_binary_structure(2, 2)

        for ii,u in enumerate(ulist):
            mask = (modis_cut == u)
            if (mask.sum() >= 1800) and (mask.sum() < 3000):
                dia_mask = bd(mask, structure=struct1, iterations=200)
                sen = sen_cut.copy()
                sen[~dia_mask] = 0
                to_conv = sen[min(np.where(dia_mask)[0]):max(np.where(dia_mask)[0])+1,\
                               min(np.where(dia_mask)[1]):max(np.where(dia_mask)[1])+1]

                w = 1./(np.nansum(mask))
                k = np.zeros(mask.shape).astype('float')
                k[mask] = w

                ker = k[min(np.where(mask)[0]):max(np.where(mask)[0])+1, min(np.where(mask)[1]):max(np.where(mask)[1])+1]
                print np.unique(ker)
                conved = signal.fftconvolve(to_conv, ker, mode='valid')
                #x, y = np.where((conved-u == np.nanmin(conved-u)))[0][0]-np.ceil((to_conv.shape()[0])/2.),\
                #np.where((conved-u == np.nanmin(conved-u)))[1][0]-np.ceil((to_conv.shape()[1])/2.)
                wcosts.append([u, conved, np.where(mask)[0]+i*1000, np.where(mask)[1]+j*1000])
                print i,j
            else:
                pass
    return wcosts
    
    
data = parallel_rw_pkl(None, 'inter_sent%i'%2, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%2, 'r')
grid_z0 = ma.array(data, mask=mask)
b4 = parallel_rw_pkl(None, 'b4', 'r')

patches = np.array(zip(np.mgrid[0:10,0:10][0].ravel(), np.mgrid[0:10,0:10][1].ravel()))
pros = np.array(np.array_split(patches, 16))

par = partial(shift, modis_sent=grid_z0, Sent = b4)
pool = multiprocessing.Pool(processes=16)
data = pool.map(par, pros)
pool.close()
pool.join()
parallel_rw_pkl(data, 'avcosts', 'w')

print 'lol finished!!!!!'


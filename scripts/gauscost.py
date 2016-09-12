#ussian one
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
        modis_mask = np.zeros_like(modis_cut).astype(int)

        #sen_cut = b4
        sen_cut = Sent[s1,s2]
        #sen_modis = np.zeros_like(modis_cut)

        ulist = np.sort(np.unique(modis_cut[~modis_cut.mask]))   
        struct1 = ndimage.generate_binary_structure(2, 2)

        for ii,u in enumerate(ulist):
            mask = (modis_cut == u)
            if (mask.sum() >= 1800) and (mask.sum() < 4000):
                dia_mask = bd(mask, structure=struct1, iterations=100)
               
		data = sen_cut.copy()
                data[~dia_mask] = 0
                to_conv = data[np.nanmin(np.where(dia_mask)[0]):np.nanmax(np.where(dia_mask)[0]),\
                               np.nanmin(np.where(dia_mask)[1]):np.nanmax(np.where(dia_mask)[1])]
		print to_conv.shape
                xwin = np.nanmax(np.where(mask)[0]) - np.nanmin(np.where(mask)[0])
                ywin = np.nanmax(np.where(mask)[1]) - np.nanmin(np.where(mask)[1])
                cost = []
                for xstd in np.arange(1,100,1):
                    for ystd in np.arange(1,100,1):
                        print 'window size: ', xwin, ywin
                        xgaus = signal.gaussian(xwin, xstd)
                        ygaus = signal.gaussian(ywin, ystd)
                        gaus  = np.outer(xgaus, ygaus)
                        #x0 = slice(min(np.where(mask))[0], max(np.where(mask))[0]+1)
                        #y0 = slice(min(np.where(mask))[1], max(np.where(mask))[1]+1)
                        m = mask[np.nanmin(np.where(mask)[0]): np.nanmax(np.where(mask)[0]), np.nanmin(np.where(mask)[1]): np.nanmax(np.where(mask)[1])]
                        gaus[~m] = 0
                        print 'Gasussian:', np.unique(gaus)
			print np.unique(np.isnan(gaus))
			if np.nansum(gaus)==np.nan:
			    pass
			else:
                            gaus = gaus/(np.nansum(gaus))
                            conved = signal.fftconvolve(to_conv, gaus, mode='valid')
                            dif = abs(conved-u)
                            x = np.where((dif == np.nanmin(dif)))[0][0]-np.ceil((to_conv.shape[0])/2.)
                            y = np.where((dif == np.nanmin(dif)))[1][0]-np.ceil((to_conv.shape[1])/2.)
                            cost.append([xstd, ystd, x, y , np.nanmin(dif)])
                            print xstd,ystd
                
                cost = np.array(cost)
                c_ind = np.where(cost[:,4] == np.nanmin(cost[:,4]))
                c_min = cost[:,4][c_ind[0]]
                mxstd = cost[:,0][c_ind[0]]; mystd = cost[:,1][c_ind[0]]
                mx = cost[:,2][c_ind[0]]; my = cost[:,3][c_ind[0]]
                
                wcosts.append([mxstd, mystd, mx, my,c_min, np.where(mask)[0]+i*1000, np.where(mask)[1]+j*1000])
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
parallel_rw_pkl(data, 'gauscosts', 'w')
print 'finished!!!! lol started on 20/07!!!'

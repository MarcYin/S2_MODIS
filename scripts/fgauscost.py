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
	ulist = np.sort(np.unique(modis_cut[~(modis_sent.mask)]))
        struct1 = ndimage.generate_binary_structure(2, 2)
    for ii,u in enumerate(ulist):
        mask = (modis_cut == u)
        if (mask.sum() >= 1800) & (mask.sum() < 4000):
            dia_mask = bd(mask, structure=struct1, iterations=100)
            data = sen_cut.copy()
            data[~dia_mask] = 0
            xmin = np.nanmin(np.where(dia_mask)[0])
            xmax = np.nanmax(np.where(dia_mask)[0])
            ymin = np.nanmin(np.where(dia_mask)[1])
            ymax = np.nanmax(np.where(dia_mask)[1])
            to_conv = data[xmin:xmax, ymin:ymax]
            
       
            cost = []
            for xstd in np.arange(1,100,1):
                for ystd in np.arange(1,100,1):
                    print 'x_ystd:',xstd,ystd, 'ii:',ii,'u:',u   
                    mxmin = np.nanmin(np.where(mask)[0])
                    mxmax = np.nanmax(np.where(mask)[0])
                    mymin = np.nanmin(np.where(mask)[1])
                    mymax = np.nanmax(np.where(mask)[1])
                    m = mask[mxmin: mxmax, mymin:mymax]
                    
                    xwin = (mxmax - mxmin)
                    ywin = (mymax - mymin)
                    xgaus = signal.gaussian(xwin, xstd)
                    ygaus = signal.gaussian(ywin, ystd)
                    gaus  = np.outer(xgaus, ygaus)
                    gaus[~m] = 0
                    print 'g_shape:', gaus.shape
                    if gaus.sum()==np.nan:
                        pass
                    else:
                        gaus = gaus/(np.nansum(gaus))
                        conved = signal.fftconvolve(to_conv, gaus, mode='valid')
                        dif = abs(conved-u)
                        x = np.where((dif == np.nanmin(dif)))[0][0]-np.ceil((to_conv.shape[0])/2.)
                        y = np.where((dif == np.nanmin(dif)))[1][0]-np.ceil((to_conv.shape[1])/2.)
                        cost.append([xstd, ystd, x, y , np.nanmin(dif)])

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

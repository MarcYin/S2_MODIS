import multiprocessing
import sys
sys.path.insert(0,'python')
from fastRWpkl import *
import numpy.ma as ma
from readSent import *
from collections import Counter
import cPickle as pkl

fhead = 'data/50SMG20165100'

def shift(process):
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
                center = np.array([(xmax+xmin)/2, (ymax+ymin)/2])
                center0 = np.array([(xmax+xmin)/2, (ymax+ymin)/2])
                center1 = np.array([49, 112])
                offset = center0 - center1
                x = inds[1] + offset[0]
                y = inds[0] + offset[1]
                cmask = (x>=0)&(y>=0)&(x<1000)&(y<1000)
                sx,sy = x[cmask], y[cmask]
                if (len(sx)>0) & (len(sy)>0):
                    mean = np.nanmean(sen_cut[sx, sy])
                    brdf = modis_cut[center[0],center[1]]
                    store.append([mean, brdf])
                    print ii,u,mean, brdf
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

ind = pkl.load(open('pkls/gaus_trans.pkl','rb'))
inds = ind
data = parallel_rw_pkl(None, 'inter_sent%i'%5, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%5, 'r')
sent = readfile([12,],fhead)['B12']
sent = ScaleExtent(sent, (10980,10980))
cm = parallel_rw_pkl(None, '0510diacm', 'r')
sent[cm] = np.nan
stm = parallel_rw_pkl(None, 'std_m', 'r')
print 'finshed reading data'
data[mask] = np.nan
modis_sent = np.array(data)
Sent = sent
Stm = stm
patches = np.array(zip(np.mgrid[0:10,0:10][0].ravel(), np.mgrid[0:10,0:10][1].ravel()))
pros = np.array(np.array_split(patches, 16))

par = partial(shift)
pool = multiprocessing.Pool(processes=16)
data = pool.map(par, pros)
pool.close()
pool.join()
parallel_rw_pkl(data, 'cm12_modis', 'w')

print 'lol finished cb12!!!!!'                                    
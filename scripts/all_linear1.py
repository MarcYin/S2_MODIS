import multiprocessing
import sys
sys.path.insert(0,'python')
from fast_rw import *
import numpy.ma as ma
from readSent import *

fhead = 'data/50SMG20165100'


def shift(process, modis_sent = None, Sent = None, stm=None):
    data = []
    for i, j in process:

        s1 = slice((i*1000),(i+1)*1000)
        s2 = slice((j*1000),(j+1)*1000)
        modis_cut = modis_sent[s1,s2]*0.001
        stm = stm[s1,s2]

        modis_cut.mask[np.isnan(modis_cut)] = True
        sen_cut = Sent[s1,s2]
        ulist = np.unique(stm)
        for ii,u in enumerate(ulist):
            mask = (modis_cut == u)
            if (mask.sum() >= 1800) and (mask.sum() < 4000):
                sx,sy = np.where(mask)
                mean = np.nanmean(sen_cut[sx, sy])
                data.append([mean, u])
                print ii,u,mean
            else:
                pass
    return data


data = parallel_rw_pkl(None, 'inter_sent%i'%2, 'r')
mask = parallel_rw_pkl(None, 'inter_sentm%i'%2, 'r')
grid_z0 = ma.array(data, mask=mask)
sent = readfile([4,],fhead)['B04']
stm = parallel_rw_pkl(None, 'stm', 'r')
print 'finshed reading data'
patches = np.array(zip(np.mgrid[0:10,0:10][0].ravel(), np.mgrid[0:10,0:10][1].ravel()))
pros = np.array(np.array_split(patches, 16))

par = partial(shift, modis_sent=grid_z0, Sent = sent)
pool = multiprocessing.Pool(processes=16)
data = pool.map(par, pros)
pool.close()
pool.join()
parallel_rw_pkl(data, 'b4_modis', 'w')

print 'lol finished!!!!!'                                    

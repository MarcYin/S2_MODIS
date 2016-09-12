import sys 
sys.path.insert(0,'python')
from fast_rw import *


for i in range(1,5):
    print 'costs%i'%(i+1),i
    cost = parallel_rw_pkl(None, 'costs%i'%(i+1), 'r')
    xs = np.zeros((10000,10000))
    xs[:] = np.nan
    ys = xs.copy()
    xstd = xs.copy()
    ystd = xs.copy()
    dif = xs.copy()
    for j in cost:
        for k in j:
            for l in k:
                print i,j,k,l
                xs[l[5], l[6]] = l[2]
                ys[l[5], l[6]] = l[3]
                xstd[l[5], l[6]] = l[0]
                ystd[l[5], l[6]] = l[1]
                dif[l[5], l[6]] = l[4]
    parallel_rw_pkl(xs, 'xs%i'%i, 'w')
    parallel_rw_pkl(ys, 'ys%i'%i, 'w')
    parallel_rw_pkl(xstd, 'xstd%i'%i, 'w')
    parallel_rw_pkl(ystd, 'ystds%i'%i, 'w')
    parallel_rw_pkl(dif, 'dif%i'%i, 'w')

print 'finished!!!'

import sys
sys.path.insert(0,'python')
from readSent import *
import scipy
import scipy.stats
from scipy import interpolate, ndimage,optimize
from classification import *
from scipy.ndimage.morphology import binary_dilation as bd
import numpy as np

def get_cloud_mask(sfpath, do_post=False):


    cl = classification(fhead = sfpath, bands = (2,3,4,8,11,12,13), bounds = None)
    cl.Get_cm_p()

    struct1 = ndimage.generate_binary_structure(2, 1)
    mask = bd(cl.cm, structure=struct1, iterations=4).astype(cl.cm.dtype)

    m = cl.cp>0# a combination of all posible cloud mask to maximumly remove those cloud pixels
    cm = (m|mask)
    rr,cc = np.where(~cm)    
    cf = rr,cc
    del cl
    
    if np.nansum(cm)/(10980.*10980.)>0.6:
        print 'more than 60% are cloud, cloud proportion: ',np.nansum(cm)/(10980.*10980.)
        return cm
    elif do_post:
        bands = readfile([2, 4,8], sfpath)
        b2 = bands['B02']
        b4 = bands['B04']
        b8 = bands['B08']
        cm = post_cloud(b2,b4,b8,cf)
        return cm
    else:
        return cm

def cost(p, bands=None, cf=None):
    rr,cc= cf
    p1, p2, p3  = p
    tc = p1*b2[rr,cc]-p2*b4[rr,cc]
    cl = tc<tc.mean()-p3*tc.std()
    clpix = [rr[cl], cc[cl]]
    r = scipy.stats.linregress(b4[clpix[0], clpix[1]],b2[clpix[0], clpix[1]])
    return 1-r.rvalue

def hot(b2,b4,b8,cf, op=False):
    rr,cc=cf
    p=[ 1.41283864e+04, 5.41603017e+03, 1.06467703e+00]
    if op:
        psolve = optimize.fmin(cost, p ,args=((b2,b4,b8),cf), full_output=1)
        p = psolve[0]
    else:
        pass
    
    tc = p[0]*b2[rr,cc]- p[1]*b4[rr,cc]
    
    cl = tc<(tc.mean()-p[2]*tc.std())
    clpix = [rr[cl], cc[cl]]
    r = scipy.stats.linregress(b4[clpix[0], clpix[1]],b2[clpix[0], clpix[1]])
    arc = np.arctan(r.slope)
    
    arc = np.arctan(r.slope)
    #hzpix = [rr[~cl], cc[~cl]]
    
    hot = np.sin(arc)*b2 - np.cos(arc)*b4
    
    #hzhot = np.sin(arc)*b2[hzpix[0], hzpix[1]] - np.cos(arc)*b4[hzpix[0], hzpix[1]]
    clmean = np.nanmean(np.sin(arc)*b2[clpix[0], clpix[1]] - np.cos(arc)*b4[clpix[0], clpix[1]])
    
    return hot, clmean

def post_cloud(b2,b4,b8,cf, do_it_anyway=False):
    
    all_hot, clmean = hot(b2,b4,b8, cf)
    ints = [clmean+i*0.001 for i in range(100)]
    lbs = np.array([np.histogram(b8[all_hot>i])[1][0] for i in ints])
    
    adjust = lbs[lbs>0] - lbs[0]
    y,x = np.array(ints)[adjust>0], adjust[adjust>0]
    r = scipy.stats.linregress(y,x)
    if r.rvalue<0.95:
        print 'Rvalue is low: ', r.rvalue, 'Maybe alterlative cloud classification is needed.'
        if do_it_anyway:
            print r
            test = np.zeros_like(b8)
            test[:] = np.nan
            test[all_hot>clmean] = all_hot[all_hot>clmean]*r.slope + r.intercept
            test[all_hot<clmean] = 0
            test[test<0] = 0

            cloud = test > 0.15
            struct1 = ndimage.generate_binary_structure(2, 1)
            cloud = bd(cloud, structure=struct1, iterations=3).astype(cloud.dtype)
        else:
            return 0
    else:   
        print r
        test = np.zeros_like(b8)
        test[:] = np.nan
        test[all_hot>clmean] = all_hot[all_hot>clmean]*r.slope + r.intercept
        test[all_hot<clmean] = 0
        test[test<0] = 0
        
        cloud = test > 0.15
        struct1 = ndimage.generate_binary_structure(2, 1)
        cloud = bd(cloud, structure=struct1, iterations=3).astype(cloud.dtype)
        return cloud
import gdal
import numpy as np
import numpy.ma as ma
import kernels
from geo_trans import *
from readSent import *

def r_modis(fname):
    g = gdal.Open(fname)
    if g is None:
        raise IOError
    else:
        return g.ReadAsArray()


def ScaleExtent(data, shape): # used for unifine different array,

    re = int(shape[0]/(data.shape[0]))

    a = np.repeat(np.repeat(data, re, axis = 1), re, axis =0)
    
    if re*(data.shape[0]-shape[0]) != 0:
        extended = np.zeros(shape)
        extended[:re*(data.shape[0]),:re*(data.shape[0])] = a
        extended[re*(data.shape[0]):,re*(data.shape[0]):] = a[re*(data.shape[0])-shape[0]:, re*(data.shape[0])-shape[0]]
        return extended
    else:
        return a
#bands = [2,3,4,8,13,11,12]


def get_kk(fhead):
    mete = readxml('%smetadata.xml'%fhead)
    vZa = np.array([i+ np.zeros((23,23)) for i in mete['mVz']])
    vAa = np.array([i+ np.zeros((23,23)) for i in mete['mVa']])
    sza = mete['mSz'][0]
    vza = mete['mVz'][3]
    rel_a = (mete['mSa']-mete['mVa'])[3]
    kk = kernels.Kernels(vza ,sza,rel_a,\
                         RossHS=False,MODISSPARSE=True,\
                         RecipFlag=True,normalise=1,\
                         doIntegrals=False,LiType='Dense',RossType='Thick')
    return kk



def get_rs(modisQA, modis_filenames, fhead):
    
    kk = get_kk(fhead)
    k_vol = kk.Ross
    k_geo   = kk.Li
    iso =  kk.Isotropic
    
    brdfs = []
    QA = r_modis(modisQA[0][0])


    for i in [2,3,0,1,5,6]:
        br = r_modis(modis_filenames[i][0])
    
        mask = (br[0] > 32766) | (br[1] > 32766) |(br[2] > 32766)| (QA>1)
        brdf = br[0] + br[1]*k_vol + br[2]*k_geo
        brdf = ma.array(brdf, mask = mask)
        brdfs.append(brdf)
        
    return brdfs

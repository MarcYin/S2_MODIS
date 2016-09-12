import gdal
import numpy as np
import numpy.ma as ma
import sys
sys.path.insert(0,'python')
from readSent import *
import kernels
from geo_trans import *
from scipy.interpolate import griddata
from fast_rw import *


modis_filenames = gdal.Open('m_data/MCD43A1.A2016105.h27v05.005.2016122100738.hdf').GetSubDatasets()
modisQA = gdal.Open("m_data/MCD43A2.A2016105.h27v05.005.2016122100739.hdf").GetSubDatasets()
h=27
v=5
fhead = 'data/50SMG20164100'

def get_coords(h,v):
    mgrss = get_lon_lat(27, 5).ravel()

    mgrss = np.array([(i[:5],i[-8:-4],i[-4:]) for i in mgrss]).reshape(2400,2400,3)

    index = np.where(mgrss[:,:,0]=='50SMG')
    Scoords = [9999-mgrss[index[0], index[1],2].astype('int'), mgrss[index[0], index[1],1].astype('int')]
    return index, Scoords

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

    ks = []
    for i in [2,3,4,8,11,12]:
        if i ==8:
            i -= 1
            kk = kernels.Kernels(vZa[i] ,mete['SAG_Z'],mete['SAG_A']-vAa[i],\
                                 RossHS=False,MODISSPARSE=True,\
                                 RecipFlag=True,normalise=1,\
                                 doIntegrals=False,LiType='Dense',RossType='Thick')
            ks.append(kk)
            kk = kernels.Kernels(vZa[i+1] ,mete['SAG_Z'],mete['SAG_A']-vAa[i+1],\
                                 RossHS=False,MODISSPARSE=True,\
                                 RecipFlag=True,normalise=1,\
                                 doIntegrals=False,LiType='Dense',RossType='Thick')
            ks.append(kk)

        else:
            i -= 1
            kk = kernels.Kernels(vZa[i] ,mete['SAG_Z'],mete['SAG_A']-vAa[i],\
                                 RossHS=False,MODISSPARSE=True,\
                                 RecipFlag=True,normalise=1,\
                                 doIntegrals=False,LiType='Dense',RossType='Thick')
            ks.append(kk)
    return ks



def get_rs(modisQA, modis_filenames, fhead):
    
    ks = get_kk(fhead)
    
    brdfs = []
    brdfms = []

    QA = r_modis(modisQA[0][0])

    for i in [2,3,0,1,5,6]:
        br = r_modis(modis_filenames[i][0])
        brdfs.append(br)
        brdfms.append((br[0] > 32766) | (br[1] > 32766) |(br[2] > 32766)| (QA>1))


    Rs = []
    for i,j in enumerate(brdfs):
        if i == 3:
            kk = ks[i]
            ross = kk.Ross
            li   = kk.Li
            iso =  kk.Isotropic
            k_vol = ScaleExtent(ross, (2400, 2400)); k_geo = ScaleExtent(li, (2400,2400))
            Rs.append(j[0] + j[1]*k_vol + j[2]*k_geo)

            kk = ks[i+1]
            ross = kk.Ross
            li   = kk.Li
            iso =  kk.Isotropic
            k_vol = ScaleExtent(ross, (2400, 2400)); k_geo = ScaleExtent(li, (2400,2400))
            Rs.append(j[0] + j[1]*k_vol + j[2]*k_geo)
        elif i>3:
            i +=1
            kk = ks[i]
            ross = kk.Ross
            li   = kk.Li
            iso =  kk.Isotropic
            k_vol = ScaleExtent(ross, (2400, 2400)); k_geo = ScaleExtent(li, (2400,2400))
            Rs.append(j[0] + j[1]*k_vol + j[2]*k_geo)
        else:
            kk = ks[i]
            ross = kk.Ross
            li   = kk.Li
            iso =  kk.Isotropic
            k_vol = ScaleExtent(ross, (2400, 2400)); k_geo = ScaleExtent(li, (2400,2400))
            Rs.append(j[0] + j[1]*k_vol + j[2]*k_geo)
    return Rs, brdfms

def inter_modis(modisQA, modis_filenames, fhead, h, v):
    
    Rs, brdfms = get_rs(modisQA, modis_filenames, fhead)

    grid_x, grid_y = np.mgrid[0:10980, 0:10980]
    index, Scoords = get_coords(h,v)
    for i in xrange(6):
        mask = brdfms[i]
        sentm = np.zeros((10980,10980))

        values = Rs[i][index[0],index[1]]
        inter_sent = griddata(np.array(zip(Scoords[0],Scoords[1])), values, (grid_x, grid_y), method='nearest')
        sentm[Scoords[0], Scoords[1]] = mask[index[0],index[1]]
        inter_sent = ma.array(inter_sent, mask=sentm)
        parallel_rw_pkl(inter_sent.data, 'inter_sent%i'%i, 'w')
        parallel_rw_pkl(inter_sent.mask, 'inter_sentm%i'%i, 'w')
inter_modis(modisQA, modis_filenames, fhead, h, v)

import sys
sys.path.insert(0,'python')
from geo_trans import *
import numpy as np
from scipy.interpolate import griddata
from fast_rw import *

def get_coords(h,v):
    mgrss = get_lon_lat(27, 5).ravel()

    mgrss = np.array([(i[:5],i[-8:-4],i[-4:]) for i in mgrss]).reshape(2400,2400,3)

    index = np.where(mgrss[:,:,0]=='50SMG')
    Scoords = [9999-mgrss[index[0], index[1],2].astype('int'), mgrss[index[0], index[1],1].astype('int')]
    return index, Scoords


h=27; v=5
Rs =  np.arange(2400*2400).reshape(2400,2400)
grid_x, grid_y = np.mgrid[0:10980, 0:10980]
index, Scoords = get_coords(h,v)

values = Rs[index[0],index[1]]
std_int_sent = griddata(np.array(zip(Scoords[0],Scoords[1])), values, (grid_x, grid_y), method='nearest')
parallel_rw_pkl(std_int_sent, 'std_int_sent%i', 'w')

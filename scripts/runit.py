import sys
sys.path.insert(0, 'python')
from get_grid import *
from fast_rw import *

fname = '50SMG20164100'
x_inds, y_inds = get_inds(fname)

parallel_rw_pkl(x_inds,fname+'X', o = 'w')
parallel_rw_pkl(y_inds,fname+'Y', o = 'w')



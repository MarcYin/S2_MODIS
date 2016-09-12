import cPickle as pkl
import multiprocessing
from functools import partial
import numpy as np
def w_r_pkl(data = None, o = 'none'):

    if o =='w':
        pkl.dump(data.values(), open('%s.pkl'%(data.keys()[0]), 'wb'))
    elif o == 'r':
        return {data.keys()[0]: pkl.load(open('%s.pkl'%(data.keys()[0]), 'rb'))}
    else:
        print 'Please specify operation'

def parallel_rw_pkl(data, fname, o = 'w'):

    if o =='w':
        subname = []
        for i in np.arange(16):
            subname.append('pkls/%s%i'%(fname, i))
      
        data = np.array(np.array_split(data, 16))
        dict_data = [{subname[i]: data[i]} for i in range(16)]

        par =  partial(w_r_pkl, o = 'w')
        pool = multiprocessing.Pool(processes = 16)
        pool.map(par, dict_data)
        pool.close()
        pool.join()
    if o == 'r':
        subname = []
        for i in np.arange(16):
            subname.append('pkls/%s%i'%(fname, i))
      

        dict_data = [{subname[i]: []} for i in range(16)] 
        par =  partial(w_r_pkl, o = 'r')
        pool = multiprocessing.Pool(processes = 16)
        dict_data = pool.map(par, dict_data)
        #data = np.array([dict_data[i]['pkls/%s%i'%(fname, i)] for i in range(16)])
        data = np.vstack(tuple([dict_data[i]['pkls/%s%i'%(fname, i)][0] for i in range(16)]))
        pool.close()
        pool.join()
        return data
    else:
        pass


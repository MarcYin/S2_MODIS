import cPickle as pkl
import multiprocessing
from functools import partial
import numpy as np
import glob

cores = multiprocessing.cpu_count()
core_to_use = max(cores-1, 1)
core_to_use = 40

def w_r_pkl(data = None, o = 'none'):

    if o =='w':
        pkl.dump(data.values(), open('%s.pkl'%(data.keys()[0]), 'wb'))
    elif o == 'r':
        return {data.keys()[0]: pkl.load(open('%s'%(data.keys()[0]), 'rb'))}
    else:
        print 'Please specify operation'
        
def chunks(l, n):
    '''split to many chunks'''
    if type(l) is np.ndarray:
        l = list(np.ndarray.tolist(l.ravel()))
    elif type(l) is list or tuple:
        print 'Length of data: %s'%len(l)
    else:
        print 'Date type is %s, but expected data type is ndarray or list or tuple!'%type(l)
    size = len(l)/n
    deci = len(l)-size*n
    inte = l[:size*n]
    data = [inte[i:i+size] for i in range(0,len(inte), size)]
    rest = l[size*n:]
    print rest
    if deci != 0:
        data[-1] = data[-1]+rest
    return data

def parallel_rw_pkl(data, fname, o = 'w'):
    core_to_use = 40
    if o =='w':
        if type(data) is np.ndarray:
            data = list(data.ravel())
        elif type(data) is list or tuple:
            print 'Length of data: %s'%len(data)
        else:
            print 'Date type is %s, but expected data type is ndarray or list or tuple!'%type(data)
        subname = []
        for i in np.arange(core_to_use):
            subname.append('pkls/%s%i'%(fname, i))
        data = chunks(data, core_to_use)  
        dict_data = [{subname[i]: data[i]} for i in xrange(core_to_use)]
        par =  partial(w_r_pkl, o = 'w')
        pool = multiprocessing.Pool(processes = core_to_use)
        pool.map(par, dict_data)
        pool.close()
        pool.join()
    if o == 'r':
        subname = glob.glob('pkls/%s*'%fname)
        dict_data = [{i: []} for i in subname] 
        par =  partial(w_r_pkl, o = 'r')
        pool = multiprocessing.Pool(processes = core_to_use)
        dict_data = pool.map(par, dict_data)
        #data = np.array([dict_data[i]['pkls/%s%i'%(fname, i)] for i in range(16)])
        pool.close()
        pool.join()
        temp = {}
        for i in dict_data:
            temp.update(i)
        data = [temp[i][0] for i in subname]
        
        if data[0] is np.ndarray:
            return np.vstack((data))
        elif data[0] is list or tuple:
            a = []
            for i in data:
                a+=i
                return a
    else:
        pass

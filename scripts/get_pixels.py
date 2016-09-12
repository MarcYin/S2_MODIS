# build a function to read in a time serials of tile pixels and random sample the 
#1000 pixels for each day with 13 bands value stored. 
import sys
sys.path.insert(0,'python')
from readSent import *
from classification import *
import cPickle as pkl
import numpy as np
import glob
import multiprocessing
from functools import partial



def ScaleExtent( data): # used for unifine different array, 
    shape = (10980,10980)
    re = int(shape[0]/(data.shape[0]))
    return np.repeat(np.repeat(data, re, axis = 1), re, axis =0)


def pkl_f(f_wr, o = None):
    if o == 'r':
        return pkl.load(open(f_wr, 'rb'))
    elif o == 'w':
        pkl.dump(f_wr.values()[0], open(f_wr.keys()[0], 'wb'))
    else:
        print "please specify operation as 'r'or 'w' for reading and writing"  

def get_veg(fhead):
    
    img = readfile(range(1,14),fhead)
    b2 = img['B02']
    b4 = img['B04']
    b3 = img['B03']
    b8 = img['B08']
    b11 = img['B11']
    b12 = img['B12']
    b8a = img['B8A']
    b5 = img['B05']
    b6 = img['B06']
    b10 = img['B10']
    b9 = img['B09']
    b7 = img['B07']
    b1 = img['B01']

    e_b = []
    for i in [b1,b5,b6,b7,b9,b10,b8a,b11,b12]:
        e_b.append(ScaleExtent(i))
    b1,b5,b6,b7,b9,b10,b8a,b11,b12 = e_b
    bands = [b1,b2,b3,b4,b5,b6,b7,b8,b8a,b9,b10,b11,b12]
    
    l = classification(img = img)
    l._p8()
    #l._p9()
    vm = l.vm
    
    w_p = np.where(vm)
    a = np.arange(len(w_p[0]))
    inds=[]
    p_rs = []
   
    for i in range(1000):
        ind = np.random.choice(a)
        inds.append(ind)
        p_x = w_p[0][ind]
        p_y = w_p[1][ind]
        p_r = []
        for m in bands:
            p_r.append(m[p_x,p_y])
        p_rs.append(p_r)
    
    
    
    pkl_f({'pkls/%spixels.pkl'%fhead[5:]: p_rs}, 'w')
    pkl_f({'pkls/%sinds.pkl'%fhead[5:]: inds}, 'w')
    pkl_f({'pkls/%svm.pkl'%fhead[5:]: vm}, 'w')

fheads= [i.split('B01')[0] for i in glob.glob('data/50SMG*B01.jp2')]
for i in fheads[4:]:
    get_veg(i)
    print 'finished %s'%i

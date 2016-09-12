import multiprocessing
from concurrent.futures import *
from functools import partial
import gdal



def readimg(band,fhead='data/50SMG20164100', bounds = None):
   
    if band ==8:
        filename = fhead+'B0'+'%s'%band +'.jp2'
  
    elif band ==13:
        filename = fhead+'B'+'%sA'%8 +'.jp2'
   
    elif band<10:
        filename = fhead+'B0'+'%s'%band +'.jp2'
    
    else:
        filename = fhead+'B'+'%s'%band +'.jp2'
     
    

    imgdata = {filename.split('.')[0][-3:]: []}

    g = gdal.Open(filename)
    if g is None:
        raise IOError
    if bounds is None:
        imgdata[filename.split('.')[0][-3:]]=g.ReadAsArray()/10000.
    else:
        imgdata[filename.split('.')[0][-3:]]= g.ReadAsArray(bounds[0],bounds[1],bounds[2],bounds[3])/10000.
    
    return imgdata
    


rim = partial(readimg, fhead = 'data/50SMG20164100', bounds=None)


if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=13)
    img = (pool.map(rim, [1,2,3,4,5,6,7,8,9,10,11,12,13], 1))
    pool.close()
    pool.join()






       

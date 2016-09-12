import multiprocessing
from functools import partial
import numpy as np
import numpy as np
from osgeo import ogr
from osgeo import osr
import mgrs
from scipy.interpolate import griddata
import cPickle as pkl

def s_pc(fhead):
    #fhead = 'data/50SMG20164100'
    # generate a list of the forth and fifth letter
    a = map(chr, range(65, 91))
    a.remove('I')
    a.remove('O')
    
    #get the corner coordinates for each tile
    # the mgrs gives the lower left coordinate

    t = fhead.split('201')[0][-5:]
    m1 = t
    m2 = t.replace(t[-1], a[a.index(t[-1])+1])
    m3 = t.replace(t[-2], a[a.index(t[-2])+1])
    m4 = t.replace(t[-2], a[a.index(t[-2])+1]).replace(t[-1], a[a.index(t[-1])+1])
    
    m = mgrs.MGRS()
    
    b_l = list(m.toLatLon(m1))
    u_l = list(m.toLatLon(m2))
    b_r = list(m.toLatLon(m3))
    u_r = list(m.toLatLon(m4))
    
    return [u_l, u_r, b_l, b_r]



def trans(ps, (grid_x, grid_y), coor = None):
    print ps
    coord = coor
    return griddata(ps, coord, (grid_x, grid_y), method='linear')

def coordarray(tx,fhead):
 
    c = np.array(s_pc(fhead))
    x = [0, 0, 10980, 10980]
    y = [0,10980, 0, 10980]
    ps = np.array([x,y]).T

    lats = c[:,0]
    lons = c[:,1]
    grid_x, grid_y = np.mgrid[0:10980, 0:10980]
    
    p_lat = griddata(ps, lats, (grid_x, grid_y), method='linear')
    p_lon = griddata(ps, lons, (grid_x, grid_y), method='linear')
    #par = partial(trans,ps, (grid_x, grid_y))

    #pool = multiprocessing.Pool(processes=2)
    #p_coords = pool.map(par,[c[:,0], c[:,1]],1)
    #pool.close()
    #pool.join()
    
    f_lat = p_lat.ravel()
    f_lon = p_lon.ravel()
    
    cds = zip(f_lon,f_lat)
    mc = tx.TransformPoints(cds)

    m_x = np.array([i[0] for i in mc])
    m_y = np.array([i[1] for i in mc])
    
    return (m_x, m_y)

def cal_coords(fhead):
    
    tx, htw, hcs, vtw, vcs, modis_x1, modis_y0 = geo_transform()
    
    m_x, m_y = coordarray(tx, fhead)
    
    h = ((m_x - modis_x1)/htw).astype(int) 
    v = ((modis_y0 - m_y)/vtw).astype(int)

    p_x = np.ceil((m_x - modis_x1)/hcs) - (h)*2400 # start from 0 so h not h-1 
    p_y = np.ceil((modis_y0 - m_y)/vcs) - (v)*2400
    
    return np.array([h,v,p_x, p_y, modis_x1, modis_y0])

def geo_transform():
    # specified for the MODIS and Sentinel 2 geo_transform
    
    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    modis_sinu = osr.SpatialReference() # define the SpatialReference object
    # In this case, we get the projection from a Proj4 string
    modis_sinu.ImportFromProj4 ( \
                "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    tx = osr.CoordinateTransformation( wgs84, modis_sinu )

    modis_x0, modis_y, modis_z = tx.TransformPoint (180, 0)
    modis_x1, modis_y, modis_z = tx.TransformPoint (-180, 0)
    modis_x, modis_y0, modis_z = tx.TransformPoint (0, 90)
    modis_x, modis_y1, modis_z = tx.TransformPoint (0, -90)

    #horizontal tilewidth
    htw = (modis_x0 - modis_x1)/36.
    #cell size
    hcs = htw/2400.
    # vertical tilewidth
    vtw = (modis_y0 - modis_y1)/18.
    #vertical cell size 
    vcs = vtw/2400.
    
    return (tx, htw, hcs, vtw, vcs, modis_x1, modis_y0)

coords = cal_coords('50SMG20164100')
pkl.dump(coords, open('pkls/50SMGcoords.pkl', 'wb'))

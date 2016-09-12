import numpy as np
from osgeo import ogr
from osgeo import osr
import mgrs
from scipy.interpolate import griddata

def get_lon_lat(fhead):
    #fhead = 'data/50SMG20164100'
    # generate a list of the forth and fifth letter
    a = map(chr, range(65, 91))
    a.remove('I')
    a.remove('O')
    
    #get the corner coordinates for each tile
    # the mgrs gives the lower left coordinate

    t = fhead.split('201')[0][-5:]
    bl = t
    ul = t.replace(t[-1], a[a.index(t[-1])+1])
    br = t.replace(t[-2], a[a.index(t[-2])+1])
    ur = t.replace(t[-2], a[a.index(t[-2])+1]).replace(t[-1], a[a.index(t[-1])+1])
    
    m = mgrs.MGRS()
    b_l = list(m.toLatLon(bl))
    u_l = list(m.toLatLon(ul))
    b_r = list(m.toLatLon(br))
    u_r = list(m.toLatLon(ur))
    
    return np.array([u_l, u_r, b_l, b_r])

def transform():
    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    modis_sinu = osr.SpatialReference() # define the SpatialReference object
    # In this case, we get the projection from a Proj4 string
    modis_sinu.ImportFromProj4 ( \
                    "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    tx = osr.CoordinateTransformation( wgs84, modis_sinu )
    
    return tx

def get_mextend():
    tx = transform()
    
    m_lon1, modis_lat, modis_z = tx.TransformPoint (180, 0)
    m_lon0, modis_lat, modis_z = tx.TransformPoint (-180, 0)
    modis_lon, m_lat0, modis_z = tx.TransformPoint (0, 90)
    modis_lon, m_lat1, modis_z = tx.TransformPoint (0, -90)
    
    return (m_lon0, m_lon1, m_lat0, m_lat1)

def get_steps():
    m_lon0, m_lon1, m_lat0,m_lat1 = get_mextend()
    
    lon_step = (m_lon1-m_lon0)/36.
    lat_step = (m_lat0-m_lat1)/18.
    lon_cstep = lon_step/2400.
    lat_cstep = lon_step/2400.
    
    return (lon_step, lat_step, lon_cstep, lat_cstep)


def get_coords(m_lon, m_lat):
    
    #tx = transform()
    
    #m_lon, m_lat, m_z = tx.TransformPoint(lon, lat)
    
    lon_step, lat_step, lon_cstep, lat_cstep = get_steps()
    
    m_lon0, m_lon1, m_lat0, m_lat1 = get_mextend()
    
    v = ((m_lat0 - m_lat)/lon_step).astype('int')
    h = ((m_lon - m_lon0)/lat_step).astype('int')
    
    p_x = np.ceil((m_lat0 - m_lat)/lat_cstep) - (v)*2400 
    p_y = np.ceil((m_lon - m_lon0)/lon_cstep) - (h)*2400
    
    return p_x.astype('int'), p_y.astype('int')

def interp(m_lats, m_lons):
    
    x = [0, 0, 10980, 10980]
    y = [0,10980, 0, 10980]
    corinds = np.array([x,y]).T
    
    grid_x, grid_y = np.mgrid[0:10980, 0:10980]
    
    Alats = griddata(corinds, m_lats, (grid_x, grid_y), method='linear')
    Alons = griddata(corinds, m_lons, (grid_x, grid_y), method='linear')
    
    return (Alats, Alons)

def get_inds(fhead):
    
    tx = transform()
    cors = get_lon_lat(fhead)
    #indexes = []
    #for i in cors:
     #   indexes.append(get_coords(cors[1], cors[0]))
    
    cor_lats = cors[:,0]; cor_lons = cors[:,1]
    cds = zip(cor_lons,cor_lats)
    
    mc = tx.TransformPoints(cds)
    
    m_lats = np.array([i[0] for i in mc])
    m_lons = np.array([i[1] for i in mc])
    
    m_lats, m_lons = interp(m_lats, m_lons)
    x_inds, y_inds = get_coords(m_lons, m_lats)
    
    return x_inds, y_inds
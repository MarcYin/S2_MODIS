import numpy as np
from osgeo import ogr
from osgeo import osr
import mgrs
from scipy.interpolate import griddata

def transform(a=True):
    # from prof. lewis
    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    modis_sinu = osr.SpatialReference() # define the SpatialReference object
    # In this case, we get the projection from a Proj4 string
    modis_sinu.ImportFromProj4 ( \
                    "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    if a:
        # from modis to wgs 84
        tx = osr.CoordinateTransformation(modis_sinu, wgs84)
    else:
        # from wgs 84 to modis
        tx = osr.CoordinateTransformation(wgs84,modis_sinu)
    
    return tx

def get_mextend():
    
    # get the begining and ending of the MODIS grid
    tx = transform(False)
    
    m_lon1, modis_lat, modis_z = tx.TransformPoint (180, 0)
    m_lon0, modis_lat, modis_z = tx.TransformPoint (-180, 0)
    modis_lon, m_lat0, modis_z = tx.TransformPoint (0, 90)
    modis_lon, m_lat1, modis_z = tx.TransformPoint (0, -90)
    
    return (m_lon0, m_lon1, m_lat0, m_lat1)

def get_steps():
    
    m_lon0, m_lon1, m_lat0,m_lat1 = get_mextend()
    # find out the modis tile step
    lon_step = (m_lon1-m_lon0)/36.
    lat_step = (m_lat0-m_lat1)/18.
    # find out the modis pixel step 
    lon_cstep = lon_step/2400.
    lat_cstep = lon_step/2400.
    
    return (lon_step, lat_step, lon_cstep, lat_cstep)

def get_lon_lat(h, v):
    # calculate the lat and lon for one MODIS tile
    tx = transform(True)
    
    #m_lon, m_lat, m_z = tx.TransformPoint(lon, lat)
    
    lon_step, lat_step, lon_cstep, lat_cstep = get_steps()
    
    m_lon0, m_lon1, m_lat0, m_lat1 = get_mextend()
    
    h_0 = (h)*lon_step + m_lon0 ; v_0 = m_lat0 - (v)*lat_step 
    h_e = (h+1)*lon_step + m_lon0; v_e = m_lat0 - (v+1)*lat_step 

    hs = np.arange(h_0, h_e, lon_cstep)[0:2400]
    vs = np.arange(v_e, v_0, lat_cstep)[0:2400][::-1]
    
    h_array = np.tile(hs, 2400).reshape(2400,2400).ravel()
    v_array = (np.tile(vs, 2400).reshape(2400,2400).T).ravel()
    
   
    #transformation to the wgs
    wgs = tx.TransformPoints(zip(h_array,v_array))
    m = mgrs.MGRS()
    # transformation to MGRS
    mgr = []
    for i in wgs:
        mgr.append(m.toMGRS(i[1], i[0],MGRSPrecision=4))

    return np.array(mgr).reshape(2400,2400)

def mtile_cal(lat, lon):
    # a function calculate the tile number for MODIS, based on the lat and lon
    tx = transform(False)
    ho,vt,z = tx.TransformPoint(lon, lat)
    lon_step, lat_step, lon_cstep, lat_cstep = get_steps()
    m_lon0, m_lon1, m_lat0, m_lat1 = get_mextend()
    h = int((ho - m_lon0)/lon_step)
    v = -int((vt - m_lat0)/lat_step)
    return h,v

def get_coords(lat,lon):
    # find out pixels in the MODIS tile within the MGRS tile, as MODIS tile is larger than the MGRS tile
    # here only the indexes are returned, which realise the pixel to pixel transformation.
    m = mgrs.MGRS()
    mg_coor = m.toMGRS(lat, lon, MGRSPrecision=4)
    s_area = mg_coor[:5]
    h, v = mtile_cal(lat, lon)
    mgrss = get_lon_lat(h, v).ravel()
    mgrss = np.array([(i[:5],i[-8:-4],i[-4:]) for i in mgrss]).reshape(2400,2400,3)
    index = np.where(mgrss[:,:,0] == s_area)
    d = mgrss[index[0],index[1],:]
    Scoords = [9999-d[:,2].astype('int'), d[:,1].astype('int')]
    return index, Scoords


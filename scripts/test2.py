import numpy as np
from osgeo import ogr
from osgeo import osr
import mgrs
from scipy.interpolate import griddata

def transform(a=True):
    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    modis_sinu = osr.SpatialReference() # define the SpatialReference object
    # In this case, we get the projection from a Proj4 string
    modis_sinu.ImportFromProj4 ( \
                    "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    if a:
        tx = osr.CoordinateTransformation(modis_sinu, wgs84)
    else:
        tx = osr.CoordinateTransformation(wgs84,modis_sinu)
    
    return tx

def get_mextend():
    
    tx = transform(False)
    
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

def get_lon_lat(h, v):
    
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
    
   
    
    wgs = tx.TransformPoints(zip(h_array,v_array))
    m = mgrs.MGRS()
    
    mgr = []
    for i in wgs:
        mgr.append(m.toMGRS(i[1], i[0],MGRSPrecision=0))

    return np.array(mgr).reshape(2400,2400)

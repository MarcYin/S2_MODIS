import numpy as np
from osgeo import ogr
from osgeo import osr
import mgrs
from scipy.interpolate import griddata

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

def m_pc(lat,lon):
    
    #lat = cd[0]
    #lon = cd[1]
    
    m_x, m_y, m_z = tx.TransformPoint (lon, lat)
    h = int((m_x - modis_x1)/htw) # start from 0 so int would be enough
    v = int((modis_y0 - m_y)/vtw) 
    
    
    p_x = np.ceil((m_x - modis_x1)/hcs) - (h)*2400 # start from 0 so h not h-1 
    p_y = np.ceil((modis_y0 - m_y)/vcs) - (v)*2400 # but the pixel start from 1 so use ceil value

    return (h,v,p_x,p_y)

def s_pc(fhead):
    #fhead = 'data/50SMG20164100'
    # generate a list of the forth and fifth letter
    a = map(chr, range(65, 91))
    a.remove('I')
    a.remove('O')
    
    #get the corner coordinates for each tile
    # the mgrs gives the lower left coordinate

    t = fhead[5:10]
    m1 = t
    m2 = t.replace(t[-1], a[a.index(t[-1])+1])
    m3 = t.replace(t[-2], a[a.index(t[-2])+1])
    m4 = t.replace(t[-2], a[a.index(t[-2])+1]).replace(t[-1], a[a.index(t[-1])+1])
    
    m = mgrs.MGRS()
    b_l = m.toLatLon(m1)
    u_l = m.toLatLon(m2)
    b_r = m.toLatLon(m3)
    u_r = m.toLatLon(m4)
    
    '''
    m_bl = m_pc(b_l)
    m_ul = m_pc(u_l)
    m_br = m_pc(b_r)
    m_ur = m_pc(u_r)
    '''
    
    return u_l, u_r, b_l, b_r

c = s_pc('data/50SMG20164100')
x = [0, 0, 10980, 10980]
y = [0,10980, 0, 10980]
ps = np.array([x,y]).T
lats= np.array([lat for lat, lon in c])
lons = np.array([lon for lat, lon in c])

grid_x, grid_y = np.mgrid[0:10980, 0:10980]
p_lat = griddata(ps, lats, (grid_x, grid_y), method='cubic')
p_lon = griddata(ps, lons, (grid_x, grid_y), method='cubic')
f_lat = p_lat.ravel()
f_lon = p_lon.ravel()
cds = zip(f_lon,f_lat)
mc = tx.TransformPoints(cds)

m_x = np.array([i[0] for i in mc])
m_y = np.array([i[1] for i in mc])

h = ((m_x - modis_x1)/htw).astype(int) 
v = ((modis_y0 - m_y)/vtw).astype(int)

p_x = np.ceil((m_x - modis_x1)/hcs) - (h)*2400 # start from 0 so h not h-1 
p_y = np.ceil((modis_y0 - m_y)/vcs) - (v)*2400

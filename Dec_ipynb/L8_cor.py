import sys
sys.path.insert(0, 'python')
from PSF_optimization import *
from read_net import *

lat, lon, year, month, day = 37.474469, 117.346456, 2016, 3, 10
h,v = mtile_cal(lat, lon)
#pr=get_wrs(lat, lon)
#path, row = pr[0]['path'],pr[0]['row']
path, row = 122, 34
Hfiles = glob.glob(directory +'l_data/LC8%03d%03d%d*LGN00_sr_band1.tif'%(path, row, year))
Hfiles

doy = 70
Hfile = directory +'l_data/LC8%03d%03d%d%03dLGN00_toa_'%(path, row, year, doy)
Lfile = glob.glob('m_data/MCD43A1.A%d%03d.h%02dv%02d.006.*.hdf'%(year,doy,h,v))[0]

psf = 1.14719796e+01,   1.51676255e+01,   7.69298694e+00, 3.25671574e+01,   4.54695201e+00,   1.80884279e-02
psf

def read_meta(Hfile):
    
    with open(Hfile[:-4]+'_MTL.txt', 'r') as inF:
        for line in inF:
            if 'CLOUD_COVER ' in line:
                cloud_cover =  float(line.split('= ')[1])
    if cloud_cover<20:
        #print 'Less than 20% cloud.'
        b1 = gdal.Open(Hfile+'band1.tif').ReadAsArray()
        corners = b1.shape
        dic = {}
        with open(Hfile[:-4]+'_MTL.txt', 'r') as inF:
            for line in inF:
                if ('CORNER_' in line)&('LAT_PRODUCT' in line):
                    dic[line.split(' = ')[0].strip()[7:13]] = float(line.split(' = ')[1])
                elif ('CORNER_' in line)&('LON_PRODUCT' in line):
                    dic[line.split(' = ')[0].strip()[7:13]] = float(line.split(' = ')[1])
                elif 'ROLL_ANGLE' in line:
                    vza = float(line.split(' = ')[1])
                elif 'SUN_AZIMUTH' in line:
                    saa = float(line.split(' = ')[1])
                elif 'SUN_ELEVATION' in line:
                    sza = float(line.split(' = ')[1])
        with open('Landsat_azm.pkl', 'r') as savefile:
            Landsat_azm = pkl.load(savefile)

        vaa = np.nanmean(Landsat_azm[(Landsat_azm[:,2]==path)&(Landsat_azm[:,3]==row)].squeeze()[:2])
        
        return sza, saa, vza, vaa, dic, corners
    else:
        print 'To much cloud: ', cloud_cover
        return None  
sza, saa, vza, vaa, dic, corners = read_meta('l_data/LC81220342016070LGN00_sr_')
L_inds, H_inds = ML_geo_trans(lat, lon, dic, corners)
Lx, Ly = L_inds
Hx, Hy = H_inds

tems = np.zeros((3,6))
tems[0,:] = sza
tems[1,:] = vza
tems[2,:] = vaa-saa
brdf, qa = get_brdf_six(Lfile, (tems[0], tems[1], tems[2]), bands=[3,4,1,2,6,7], flag=None, Linds= L_inds)

cloud = gdal.Open(Hfile[:-5]+'_cfmask.tif').ReadAsArray()
cl_mask = cloud==4 # cloud pixels; strictest way is to set the clear pixels with cloud==0
struct = ndimage.generate_binary_structure(2, 2)
dia_cloud = ndimage.binary_dilation(cl_mask, structure=struct, iterations=20).astype(cl_mask.dtype)


def gaussian( xstd, ystd, angle, norm = True):
    win = int(round(max(2*1.69*xstd, 3*ystd)))
    winx = win*2**0.5
    winy = win*2**0.5

    xgaus = signal.gaussian(winx, xstd)
    ygaus = signal.gaussian(winy, ystd)
    gaus  = np.outer(xgaus, ygaus)
    r_gaus = ndimage.interpolation.rotate(gaus, angle, reshape=True)
    center = np.array(r_gaus.shape)/2
    cgaus = r_gaus[center[0]-win/2: center[0]+win/2, center[1]-win/2:center[1]+win/2]
    if norm:
        return cgaus/cgaus.sum()
    else:
        return cgaus 
    
shape =  dia_cloud.shape
xstd,ystd, angle, xs, ys = psf[:5]
shx, shy = (Hx+xs).astype(int), (Hy+ys).astype(int)
val = (Hx+xs<shape[0])&(Hy+ys<shape[1])&(Hx+xs>0)&(Hy+ys>0)
ker = gaussian(xstd,ystd,angle,True)

def L8_get_to_cor(band, Hfile=None, ker = None, dia_cloud=None):
    fname = Hfile + 'band%s.tif'%band
    data = gdal.Open(fname).ReadAsArray()*0.0001
    mask = ~(data<=0).astype('bool')
    struct = ndimage.generate_binary_structure(2, 2)
    small_mask = ndimage.binary_erosion(mask, structure=struct, iterations=20).astype(mask.dtype)
    val_mask = (~dia_cloud)&small_mask
    used_mask = val_mask[shx[val], shy[val]]
    used_data = signal.fftconvolve(data, ker, mode='same')[shx[val], shy[val]]
    
    return used_data, used_mask 
bands = [2,3,4,5,6,7]
par = partial(L8_get_to_cor, Hfile=Hfile, ker = ker, dia_cloud=dia_cloud)
retval = parmap(par, bands, nprocs=len(bands))


L8_mask = np.array(retval)[:,1,:].astype(bool)
L8_data = np.array(retval)[:,0,:]
l8 = L8_data.copy()
Mcomb_mask = np.all(qa==0, axis=0)
Scomb_mask = np.all(L8_mask, axis = 0)
l8[:,(~Scomb_mask)|(~Mcomb_mask[val])]=np.nan
l8[np.isnan(l8)], brdf[np.isnan(brdf)] = -9999999, -9999999
mas = np.all((brdf[:,val]>0)&(brdf[:,val]<1)&(l8>0)&(l8<1), axis=0)
to_cor = shx[val][mas], shy[val][mas],l8[:,mas], brdf[:,val][:,mas]



L8_emus = parallel_rw_pkl(None, '6S_emulation_L8_', 'r')
aot = read_net(year, month,day, np.arange(lat,lat+1, 0.125), np.arange(lon,lon+2, 0.125), dataset='aot')
twv = read_net(year, month,day, np.arange(lat,lat+1, 0.125), np.arange(lon,lon+2, 0.125), dataset='wv')
tco = read_net(year, month,day, np.arange(lat,lat+1, 0.125), np.arange(lon,lon+2, 0.125), dataset='tco')


import gdal
def elevation(lat, lon, north=True, east=True):
    lats = range( int(lat.min()), int(lat.max())+1)
    lons = range( int(lon.min()), int(this_lon.max())+1)
    eles=np.zeros_like(lat)
    for lat0 in lats:
        for lon0 in lons:
            if lat0>=0:
                lat_name = 'N%d'%int(lat0)
            elif lat0<0:
                lat_name = 'S%d'%int(lat0)
            else:
                'Wrong lat given, and float is expected!'
            if lon0>=0:
                lon_name = 'E%d'%int(lon0)
            elif lon0<0:
                lon_name = 'W%d'%int(lon0)
            else:
                'Wrong lon given, and float is expected!'
            fname = 'SRTM/'+lat_name+lon_name+ '.hgt'
            mask = (lat>lat0)&(lat<lat0+1)&(lon>lon0)&(lon<lon0+1)
            g = gdal.Open(fname)
            geo = g.GetGeoTransform()
            l_lon, lon_size,l_lat, lat_size, = geo[0], geo[1],geo[3],geo[5]
            x, y = ((lon-l_lon)/lon_size).astype(int), ((lat-l_lat)/lat_size).astype(int)
            ele = g.ReadAsArray()
            eles[mask] = ele[x[mask],y[mask]]    
    return eles

def cost(p, args = None):
    if any(p>10) or any(p<0):
        return 10000
    else:
        #'aot550', 'water', 'ozone'
        aot550, water = p
        TOA_refs, M_refs, angles, ele, ozone = args        
        sz, sa, vz, va = angles        
        Sur_refs = [L8_emus[ind][0].predict(np.array([[toa_ref, aot550, water, ozone, \
                                                            np.sin(sz), np.sin(vz), np.cos((sa-va)), \
                                                            ele],]))[0][0] for ind, toa_ref in enumerate(TOA_refs)]
        Sur_refs = np.array(Sur_refs)
        M_refs = np.array(M_refs)
        cost = sum(abs(Sur_refs-M_refs)*w)    
        return cost


def opt(ind):
    
    #sent_refs, modis_refs = np.array([refs[ii][tuple(aoi[ind])] for ii in range(7)]).T
    
    TOA_refs = to_cor[2][:,ind]
    M_refs = to_cor[3][:,ind]
    m = mgrs.MGRS()
    pix_lat, pix_lon = cor_inter(np.array([[to_cor[0][ind], to_cor[1][ind]],]).T, dic, corners)
    ele = eles[ind]
    inx_lat, inx_lon = (np.abs(lats-pix_lat)).argmin(),(np.abs(lons-pix_lon)).argmin()
    aot0, tcw0, tco0 = aot[inx_lat, inx_lon], twv[inx_lat, inx_lon]/10., tco[inx_lat, inx_lon]
                                                          
    angles =[i*np.pi/180 for i in [sza, saa, vza, vaa]]
    
    ozone = tco0*46.698
    
    args = TOA_refs, M_refs , angles, ele, ozone
    #print args
    p = aot0, tcw0 
    psolve = optimize.fmin_l_bfgs_b(cost,p, iprint=-1, approx_grad=1, args=(args,))
    return [to_cor[0][ind],to_cor[1][ind],psolve]

wl = np.array([482.04,561.41,654.59,864.67,1608.86,2200.73])/1000
alpha = 1.42 #angstrom exponent for continental type aerosols
w = (np.array(wl)/wl[0])**(-alpha)
lats, lons = np.arange(lat,lat+1, 0.125), np.arange(lon,lon+2, 0.125)
this_lat, this_lon = cor_inter(np.array([to_cor[0], to_cor[1]]), dic, corners)
eles = elevation(this_lat, this_lon)/1000.


pool = multiprocessing.Pool(processes=16)
retval = pool.map(opt, range(len(to_cor[0])))
pool.close()
pool.join()

parallel_rw_pkl(retval, 'L8_122034_aot', 'w')



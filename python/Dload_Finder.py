import sys
sys.path.insert(0, 'python')
from sentinel_downloader import *
from get_modis import *
from geo_trans import *
import glob

directory = '/home/ucfafyi/DATA/Sentinel/'

def dload_sent_mod(lat, lon, start, end, directory = '/home/ucfafyi/DATA/Sentinel/', flist=True):
    try:
        download_sentinel_amazon(lat, lon, start, directory+'s_data/', end_date=end)
    except:
        pass
    h,v = mtile_cal(lat,lon)
    get_modisfiles( 'MOTA', 'MCD43A1.005',start.year , 'h%02dv%02d'%(h,v), None, 
                       doy_start=start.timetuple().tm_yday, doy_end=end.timetuple().tm_yday, out_dir=directory+'m_data/')
    get_modisfiles( 'MOTA', 'MCD43A2.005',start.year , 'h%02dv%02d'%(h,v), None, 
                       doy_start=start.timetuple().tm_yday, doy_end=end.timetuple().tm_yday, out_dir=directory+'m_data/')
    
def file_finder(directory):
    fnames = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            fnames.append( os.path.join(path, name))
    return np.array(fnames)


def file_finder(directory):
    fnames = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            fnames.append( os.path.join(path, name))
    return np.array(fnames)

def get_closet(lat, lon, year, p=True):
    
    m = mgrs.MGRS()
    mg_coor = m.toMGRS(lat, lon, MGRSPrecision=4)
    
    h,v = mtile_cal(lat, lon)
    
    sfnames = np.sort(file_finder('s_data'))
    #print sfnames[0][7:9] ==mg_coor[:2], sfnames[0][10]==mg_coor[2], sfnames[0][12:14]==mg_coor[3:5]
    dates = np.array([i.split('/')[-5:-2] for i in sfnames if ((i[-7:]=='B01.jp2')&(i[7:9]==mg_coor[:2])&(i[10]==mg_coor[2])&(i[12:14]==mg_coor[3:5]))&(i[15:19]=='%s'%year)]).astype(int)
    sdates = np.array([datetime .datetime(i[0], i[1], i[2]).timetuple().tm_yday for i in dates])
    sdates.sort()
    mfnames = file_finder('m_data')
    mdates = np.array([i.split('.')[1][5:] for i in mfnames if (i[len('m_data')+1:len('m_data')+8]=='MCD43A1')&(i.split('.')[1][1:5]=='%s'%year)&(i.split('.')[2]=='h%02dv%02d'%(h,v))]).astype(int)
    mdates.sort()
    if p:
        print 'sentinel dates: ',sdates, '\nmodis dates, ', mdates, '\n'
    
    try:
        dif = abs(sdates.reshape(len(sdates),1)-mdates)
        si,mi = np.where(dif<4)
        sdi = []
        mdi = []
        for i,j in enumerate(sdates[si]):
            sd = datetime.datetime(year, 1, 1) + datetime.timedelta(j - 1)
            y,m,d = sd.year, sd.month, sd.day
            sdi.append(directory + 's_data'+'/%s/%s/%s/%s/%s/%s/0/'%(mg_coor[:2], mg_coor[2], mg_coor[3:5], y,m,d))
            mdi.append([directory+ii for ii in mfnames if ((ii.split('.')[1])[1:]=='%s%03d'%(year,mdates[mi][i]))&(ii.split('.')[2]=='h%02dv%02d'%(h,v))])
        if p:
            print 'modis files: ', np.array(mdi),'\n\n', 'sentinel files directory: ', np.array(sdi), '\n'
        return np.array(sdi), np.array(mdi)
    except:
        print 'One/both of two kinds of file dose not exist, see above if print is True!!'
        return 0
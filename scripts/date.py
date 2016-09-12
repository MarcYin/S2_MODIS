fhead = 'data/50SMG20164100'
import glob
import datetime
import calendar
import sys
sys.path.insert(0,'python')
from get_modis import *
fh = fhead.split('20')[0]
files = glob.glob('%s*B01.jp2'%fh)
dates = []
ds = []
de = []
sys = []
eys = []
for i in files:
    date = datetime.datetime.strptime((i.split('%s'%fh)[1]).split('0B01.jp2')[0], "%Y%m%d")
    dates.append(date)
    doy = date.timetuple().tm_yday
    if doy - 8 <0:
        sys.append(date.year-1)
        days = (datetime.date(date.year,1,1)-datetime.date(date.year-1,1,1)).days
        ds.append(doy-8+days)
    else:
        sys.append(date.year)
        ds.append(doy-8)

    if doy + 8 > 366:
        eys.append(date.year+1)
        days = (datetime.date(date.year+1,1,1)-datetime.date(date.year,1,1)).days
        ds.append(doy+8-days)
    else:
        eys.append(date.year)
        de.append(doy+8)
for i in range(len(sys)-1):
        
    if sys[i] == eys[i]:
        get_modisfiles( 'MOTA', 'MCD43A2.005', sys[i], 'h27v05', None, doy_start=ds[i], doy_end=de[i], out_dir='m_data/' )
        get_modisfiles( 'MOTA', 'MCD43A3.005', sys[i], 'h27v05', None, doy_start=ds[i], doy_end=de[i], out_dir='m_data/' )
    else :
        days = (datetime.date(sys[i]+1,1,1)-datetime.date(sys[i],1,1)).days
        get_modisfiles( 'MOTA', 'MCD43A2.005', 
                       sys[i], 'h27v05', None, doy_start=ds[i], 
                       doy_end=days, out_dir='m_data/' )
	get_modisfiles( 'MOTA', 'MCD43A2.005',
                       eys[i], 'h27v05', None, doy_start=1,
                       doy_end=de[i], out_dir='m_data/' )
	get_modisfiles( 'MOTA', 'MCD43A3.005',
                       sys[i], 'h27v05', None, doy_start=ds[i],
                       doy_end=days, out_dir='m_data/' )
        get_modisfiles( 'MOTA', 'MCD43A3.005',
                       eys[i], 'h27v05', None, doy_start=1,
                       doy_end=de[i], out_dir='m_data/' )
        
       
      
     
       

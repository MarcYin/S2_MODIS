import numpy.ma as ma
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
from matplotlib import gridspec
import pylab as plt


import sys
sys.path.insert(0,'./python')
from readSent import *
img = readfile((2,3,4,8,11,12,13), 'data/50SMG20164100')
from classification import *
l = classification(img = img)
l.Get_cm_p()
vm = l.vm
wm = l.wm
wm2 = l.wm2
wm3 = l.wm3
bsm = l.bsm
Dssm = l.Dssm
sm = l.sm
NDVI = l.NDVI
fcm = l.fcm
cm = l.cm
cp = l.cp
b2 = l.b2
b3 = l.b3
b4 = l.b4
b8 = l.b8
b8a = l.eb8a
b11 = l.eb11
b12 = l.eb12
del l

files = readfile([1,5,6,7,9,10],'data/50SMG20164100')
b5 = files['B05']
b6 = files['B06']
b10 = files['B10']
b9 = files['B09']
b7 = files['B07']
b1 = files['B01']

def ScaleExtent( data): # used for unifine different array, 
    shape = (10980,10980)
    re = int(shape[0]/(data.shape[0]))
    return np.repeat(np.repeat(data, re, axis = 1), re, axis =0)

a = []
for i in [b1,b5,b6,b7,b9,b10]:
    a.append(ScaleExtent(i))
b1,b5,b6,b7,b9,b10 = a

plt.ioff()

def zoomindraw(m_name,mask,rgb,bands):
    maxm = 1.
    phi = 1.
    theta = 1.3

    for k in range(11):
        if k ==0:
            
            w_p = np.where(mask)
            a = np.arange(len(w_p[0]))
        else:
            ind = np.random.choice(a)
            p_x = w_p[0][ind]
            p_y = w_p[1][ind]

          
            cmap = plt.get_cmap('gray')

            fig = plt.figure(figsize=(10,10))
            ax1 = plt.subplot2grid((2, 2), (0, 0))
            ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=2)
            ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2, rowspan=2)

            patch = patches.Rectangle((max(0,p_y/10.-25.), max(0,p_x/10.-25.)), 50+min(0,(p_y/10.-25.)),50+min(0,(p_x/10.-25.)),facecolor='none', linewidth=0.5, edgecolor = 'r')
            ax1.add_patch(patch)
            ax1.imshow(rgb[::10,::10,:]*1.7, interpolation =None)
            #plot1 = (((maxm)/phi)*(rgb[::10,::10,:]/(maxm/theta)))**0.9
            #ax1.imshow(plot1, interpolation = None)
            ax1.set_xlim([0, 1098])
            ax1.set_ylim([1098, 0])
            ax1.plot(p_y/10.,p_x/10.,'ro',markersize=2.8)
            
            
            p2 = rgb[max(0,p_x-50):min(10980,p_x+50), max(0,p_y-50):min(10980,p_y+50),:]
            plot2 = (p2-p2.min())/(p2.max()-p2.min())
            ax2.imshow(plot2, interpolation=None)
            ax2.set_xlim([0, 100])
            ax2.set_ylim([100, 0])
            ax2.plot(50,50,'ro',markersize=5)

            p_r = []
            for m in bands:
                p_r.append(m[p_x,p_y])
            plt.plot(wl,p_r)
            plt.ylim([0,1.])
            plt.ylabel('TOA reflectance')
            plt.xlabel('wavelength (nm)')
            for i in wl:
                line = ax3.plot([i,i],[0,1],'--')
                plt.setp(line, linewidth=0.5, color='gray')

            plt.tight_layout()
            plt.savefig('%s%i.pdf'%(m_name,k), format='pdf')
	    plt.close(fig)
mc_m = (cp < 0.11) & (cp>0.03) & cm
sb = """443	20
490	65
560	35
665	30
705	15
740	15
783	20
842	115
865	20
945	20
1380	30
1610	90
2190	180"""
sb= [i.split('\t') for i in sb.split('\n')]
wl = [float(i) for i,j in sb]
rgb = np.dstack((b4,b3,b2))
bands = [b1,b2,b3,b4,b5,b6,b7,b8,b8a,b9,b10,b11,b12]
for i,j in zip(['cloud','water','bare_soil','Dss','veg'],[mc_m,wm,bsm,Dssm,vm]):
    zoomindraw(i,j,rgb,bands)

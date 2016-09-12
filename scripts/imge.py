import matplotlib.patches as patches
from matplotlib import gridspec
import pylab as plt
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
            for m in bs:
                p_r.append(m[p_x,p_y])
            plt.plot(wl,p_r)
            plt.ylim([0,1.])
            plt.ylabel('TOA reflectance')
            plt.xlabel('wavelength (nm)')
            for i in wl:
                line = ax3.plot([i,i],[0,1],'--')
                plt.setp(line, linewidth=0.5, color='gray')

            plt.tight_layout()
            plt.savefig('image/%s%i.pdf'%(m_name[:-1],k), format='pdf')

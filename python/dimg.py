import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
'''
plt.figure()
ax = plt.gca()
im = ax.imshow(np.arange(100).reshape((10,10)))

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)
'''

def dimg(data, b=False, step=1, cbar = True):
    if b:
        plt.figure(figsize=(15,15))
        ax = plt.gca()
        im = ax.imshow(data[::step, ::step])
        
    else:
        plt.figure()
        ax = plt.gca()
        im = ax.imshow(data[::step,::step])
        
    if cbar:
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)
    else:
        pass

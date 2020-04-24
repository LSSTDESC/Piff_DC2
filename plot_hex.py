import numpy as np
import fitsio
import treecorr
import pickle

import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

plt.style.use('/global/u2/m/mjarvis/SVA1StyleSheet.mplstyle')

def make_fig(u,v,dT,de1,de2,title,lo=-2.3,hi=2.3, ngrid=300):
    vmax = 0.01
    fig, axs = plt.subplots(ncols=3, sharey=True, figsize=(15, 4))
    fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
    ax = axs[0]
    hb = ax.hexbin(u, v, dT, gridsize=ngrid, cmap='inferno', vmin=-vmax, vmax=vmax)
    ax.set_title(r'$\Delta T/T$')
    ax.set_aspect('equal')
    ax.set_xlim(lo,hi)
    ax.set_ylim(lo,hi)

    ax = axs[1]
    hb = ax.hexbin(u, v, de1, gridsize=ngrid, cmap='inferno', vmin=-vmax, vmax=vmax)
    ax.set_title(r'$\Delta e1$')
    ax.set_aspect('equal')
    ax.set_xlim(lo,hi)
    ax.set_ylim(lo,hi)
    ax.set_xlabel(title)

    ax = axs[2]
    hb = ax.hexbin(u, v, de2, gridsize=ngrid, cmap='inferno', vmin=-vmax, vmax=vmax)
    ax.set_title(r'$\Delta e2$')
    ax.set_aspect('equal')
    ax.set_xlim(lo,hi)
    ax.set_ylim(lo,hi)
    cb = fig.colorbar(hb, ax=ax)

    return fig

def plot_hex(data, band):
    if psfex:
        psf = 'psfex'
    else:
        psf = 'psf'

    print('band ',band)

    x = data['x']
    y = data['y']

    e1 = data['g1_data']
    e2 = data['g2_data']
    T = data['T_data']
    de1 = data['g1_data'] - data['g1_model']
    de2 = data['g2_data'] - data['g2_model']
    dT = (data['T_data'] - data['T_model']) / data['T_data']

    fig = make_fig(x,y,dT,de1,de2,'Chip coordinates',0,4096, 60)
    filename = psf + '_resid_chip_%s.pdf'%band
    plt.savefig(filename, bbox_inches='tight')
    print('wrote',filename)
    plt.close(fig)

    if 'u' not in data.dtype.names:
        return
    u = data['u']
    v = data['v']
    u = u / 3600.  # Don't use /=, since that will modify inputs.
    v = v / 3600.
    print('u range = ',np.min(u), np.max(u))
    print('v range = ',np.min(v), np.max(v))

    fig = make_fig(u,v,dT,de1,de2,'Sky coordinates')
    filename = psf + '_resid_sky_%s.pdf'%band
    plt.savefig(filename, bbox_inches='tight')
    print('wrote',filename)
    plt.close(fig)

    if 'u_cam' not in data.dtype.names:
        return
    mask = data['u_cam'] < 1.e5  # 1.e100 is used when no camera coords are available
    u2 = data['u_cam']
    v2 = data['v_cam']
    u2 = u2 / 3600.
    v2 = v2 / 3600.
    print('u2 range = ',np.min(u2[mask]), np.max(u2[mask]))
    print('v2 range = ',np.min(v2[mask]), np.max(v2[mask]))

    fig = make_fig(u2[mask],v2[mask],dT[mask],de1[mask],de2[mask],'Focal plane coordinates')
    filename = psf + '_resid_focal_%s.pdf'%band
    plt.savefig(filename, bbox_inches='tight')
    print('wrote',filename)
    plt.close(fig)


bands = ['u','g','r','i','z','y']
#bands = ['r','i','z']
all_data = {}
psfex = False
psfex = True

# Do each band separately.
for band in bands:
    if psfex:
        data = fitsio.read('psfex_%s.fits'%band)
    else:
        data = fitsio.read('psf_%s.fits'%band)
    all_data[band] = data
    plot_hex(data, band)

# Now do riz combined.
if 'r' in bands and 'i' in bands and 'z' in bands:
    data = np.concatenate([all_data[b] for b in 'riz'])
    plot_hex(data, 'riz')

if len(bands) == 6:
    data = np.concatenate([all_data[b] for b in bands])
    plot_hex(data, 'all')

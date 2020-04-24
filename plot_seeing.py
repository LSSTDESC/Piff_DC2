import numpy as np
import fitsio
import treecorr
import pickle

import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

plt.style.use('/global/u2/m/mjarvis/SVA1StyleSheet.mplstyle')

def plot_seeing(data, band):

    print('band ',band)
    T = data['T_data']
    # The version of piff that ran this labeled this wrong.  T is really sigma.
    #seeing = np.sqrt(T/2) * 2.3548
    seeing = T * 2.3548
    fig, ax = plt.subplots(ncols=1, figsize=(8, 4))
    ax.hist(seeing, bins=50, range=[0.5,1.5] )
    ax.set_title(r'%s band seeing'%band)
    ax.set_xlabel('seeing fwhm (arcsec)')
    ax.set_xlim(0.5,1.5)

    ax.axvline(x=0.75, ymin=0, ymax=1.e6, color='blue')

    filename = 'seeing_hist_%s.pdf'%band
    plt.savefig(filename, bbox_inches='tight')
    print('wrote',filename)


bands = ['u','g','r','i','z','y']
all_data = {}

# Do each band separately.
for band in bands:
    data = fitsio.read('psf_%s.fits'%band)
    all_data[band] = data
    plot_seeing(data, band)

# Now do riz combined.
if 'r' in bands and 'i' in bands and 'z' in bands:
    data = np.concatenate([all_data[b] for b in 'riz'])
    plot_seeing(data, 'riz')

if len(bands) > 1:
    data = np.concatenate([all_data[b] for b in bands])
    plot_seeing(data, 'all')

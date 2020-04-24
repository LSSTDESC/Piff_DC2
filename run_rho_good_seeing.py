import numpy as np
import fitsio
import treecorr
import pickle

def run_rho(data, results, band):
    dt = (data['T_data'] - data['T_model']) / data['T_data']
    ecat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], ra_units='deg', dec_units='deg',
                            g1=data['g1_data'], g2=data['g2_data'])
    qcat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], ra_units='deg', dec_units='deg',
                            g1=data['g1_data'] - data['g1_model'],
                            g2=data['g2_data'] - data['g2_model'])
    wcat = treecorr.Catalog(ra=data['ra'], dec=data['dec'], ra_units='deg', dec_units='deg',
                            g1=data['g1_data']*dt, g2=data['g2_data']*dt)
    ecat.name = 'ecat'
    qcat.name = 'qcat'
    wcat.name = 'wcat'

    bin_config = dict(
        sep_units = 'arcmin',
        bin_slop = 0.1,
        min_sep = 0.5,
        max_sep = 250.,
        bin_size = 0.2,
    )
    pairs = [ (qcat, qcat),
              (ecat, qcat),
              (wcat, wcat),
              (qcat, wcat),
              (ecat, wcat),
              (ecat, ecat) ]
    print('data has %s stars'%ecat.ntot)
    for (cat1, cat2) in pairs:
        print('Doing correlation of %s band %s vs %s'%(band, cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results[(band, cat1.name, cat2.name)] = rho

#bands = ['u','g','r','i','z','y']
bands = ['r','i','z']
all_data = {}
results = {}

# Do each band separately.
for band in bands:
    data = fitsio.read('psf_%s.fits'%band)
    print('len(data) for ',band,' = ',len(data))
    tot = len(data)
    seeing = data['T_data'] * 2.3548
    data = data[seeing < 0.75]
    print('good seeing = ',len(data))
    print('frac = ',len(data) / tot)
    all_data[band] = data
    run_rho(data, results, band)

# Now do riz combined.
if 'r' in bands and 'i' in bands and 'z' in bands:
    data = np.concatenate([all_data[b] for b in 'riz'])
    run_rho(data, results, 'riz')

if len(bands) == 6:
    data = np.concatenate([all_data[b] for b in bands])
    run_rho(data, results, 'all')

file_name = 'run_rho_good_seeing.out'
print('writing ',file_name)
with open(file_name, 'wb') as f:
    pickle.dump(results, f)


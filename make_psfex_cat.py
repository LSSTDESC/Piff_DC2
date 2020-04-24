import fitsio
import numpy as np
import glob
import os
import concurrent.futures

base_dir = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1/src/'

def make_visit_cat(t):
    i, n, d = t
    print(i,'/',n)
    visit_cat = 'psfex_visit_cats/' + d[len(base_dir):-1] + '.fits'
    if os.path.exists(visit_cat):
        return

    #print('Making visit catalog ',visit_cat)

    files = glob.glob(d + '*/src*.fits')
    #print('found %d LSSTDM cat files for this visit'%(len(files)))

    all_data = []

    for f in files:
        #print(f)
        try:
            data = fitsio.read(f)
        except OSError as e:
            print('Caught ',e,'for file',f)
            continue

        use = data['base_ClassificationExtendedness_value'] == 0
        data = data[use]
        flux = data['base_PsfFlux_instFlux']

        max_flux = np.max(flux)
        min_flux = np.min(flux)
        #print('range of flux = ',2.5*np.log10(max_flux/min_flux),' mag')
        # Remove brightest 3 mag to avoid B/F.
        # And only keep 3 mag worth to avoid noisy faint stars.
        use = (flux < max_flux / 15.8) & (flux > max_flux / 251.2)
        data = data[use]
        flux = flux[use]

        ixx_data = data['ext_shapeHSM_HsmSourceMoments_xx']
        ixy_data = data['ext_shapeHSM_HsmSourceMoments_xy']
        iyy_data = data['ext_shapeHSM_HsmSourceMoments_yy']
        ixx_model = data['ext_shapeHSM_HsmPsfMoments_xx']
        ixy_model = data['ext_shapeHSM_HsmPsfMoments_xy']
        iyy_model = data['ext_shapeHSM_HsmPsfMoments_yy']

        use = ~np.isnan(ixx_data) & ~np.isnan(ixy_data) & ~np.isnan(iyy_data)
        use &= ~np.isnan(ixx_model) & ~np.isnan(ixy_model) & ~np.isnan(iyy_model)
        #use &= (ixx_data*iyy_data > ixy_data**2)
        #use &= (ixx_model*iyy_model > ixy_model**2)
        data = data[use]
        flux = flux[use]
        ixx_data = ixx_data[use]
        iyy_data = iyy_data[use]
        ixy_data = ixy_data[use]
        ixx_model = ixx_model[use]
        iyy_model = iyy_model[use]
        ixy_model = ixy_model[use]

        T_data = ixx_data + iyy_data
        g1_data = (ixx_data - iyy_data) / (T_data - (ixx_data*iyy_data - ixy_data**2)**0.5)
        g2_data = 2.*ixy_data / (T_data - (ixx_data*iyy_data - ixy_data**2)**0.5)

        T_model = ixx_model + iyy_model
        g1_model = (ixx_model - iyy_model) / (T_model - (ixx_model*iyy_model - ixy_model**2)**0.5)
        g2_model = 2.*ixy_model / (T_model - (ixx_model*iyy_model - ixy_model**2)**0.5)

        # Avoid nans
        bad = (np.isnan(T_model) | np.isnan(g1_model) | np.isnan(g2_model) | 
               np.isnan(T_data) | np.isnan(g1_data) | np.isnan(g2_data))

        ra = data['coord_ra'] * 180./np.pi
        dec = data['coord_dec'] * 180./np.pi
        x = data['base_SdssCentroid_x']
        y = data['base_SdssCentroid_y']

        # Also avoid bad outliers.
        med_T = np.median(T_data[~bad])
        med_g1 = np.median(g1_data[~bad])
        med_g2 = np.median(g2_data[~bad])
        bad |= abs(T_data-med_T) / med_T > 0.1
        bad |= abs(g1_data-med_g1) > 0.1
        bad |= abs(g2_data-med_g2) > 0.1
        bad |= abs(T_model-med_T) / med_T > 0.1
        bad |= abs(g1_model-med_g1) > 0.1
        bad |= abs(g2_model-med_g2) > 0.1
        use = ~bad

        dtype = np.dtype([('ra', '>f8'), ('dec', '>f8'), ('flux', '>f8'),
                          ('x', '>f8'), ('y', '>f8'),
                          ('T_data', '>f8'), ('g1_data', '>f8'), ('g2_data', '>f8'),
                          ('T_model', '>f8'), ('g1_model', '>f8'), ('g2_model', '>f8')])
        
        data = np.empty(len(ra[use]), dtype=dtype)
        data['ra'] = ra[use]
        data['dec'] = dec[use]
        data['flux'] = flux[use]
        data['x'] = x[use]
        data['y'] = y[use]
        data['T_data'] = T_data[use]
        data['g1_data'] = g1_data[use]
        data['g2_data'] = g2_data[use]
        data['T_model'] = T_model[use]
        data['g1_model'] = g1_model[use]
        data['g2_model'] = g2_model[use]
        all_data.append(data)

    #print('concatenating...')
    data = np.concatenate(all_data)

    #print('writing to',visit_cat)
    fitsio.write(visit_cat, data, clobber=True)

    return i, visit_cat, len(data)


for band in ['u','g','r','i','z','y']:
    dirs = glob.glob(base_dir + '*-%s/'%band)
    print('found %d visit directories for %s band'%(len(dirs),band))

    if True:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for t in executor.map(make_visit_cat,
                                  [(i,len(dirs),dirs[i]) for i in range(len(dirs))]):
                if t is not None:
                    i, visit_cat, n = t
                    print(i, visit_cat, n)
    else: # Single proc version
        for i in range(len(dirs)):
            t = make_visit_cat((i, len(dirs), dirs[i]))
            if t is not None:
                i, visit_cat, n = t
                print(i, visit_cat, n)

    print('Combining all visit data for band %s'%band)
    all_data = []
    for d in dirs:
        visit_cat = 'psfex_visit_cats/' + d[len(base_dir):-1] + '.fits'
        data = fitsio.read(visit_cat)
        all_data.append(data)

    print('concatenating...')
    data = np.concatenate(all_data)
    print('total len = ',len(data))

    file_name = 'psfex_%s.fits'%band
    print('writing to',file_name)
    fitsio.write(file_name, data, clobber=True)


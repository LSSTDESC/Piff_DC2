import fitsio
import numpy as np
import glob

def compute_uv_cam(data):
    """Calculate u_cam, v_cam that are the same as u,v, but rotated so that
    +u is parallel to the +x direction.
    """
    x = data['x']
    y = data['y']
    u = data['u']
    v = data['v']

    #print('u range = ',np.min(u),np.max(u))
    #print('v range = ',np.min(v),np.max(v))

    # Pick a random point on some chip.  It doesn't matter too much, but picking something
    # on the edge is safer.  Otherwise you can get spurious matches to our distance criterion.
    # Using the mimium u works well.
    indx = np.argmin(u)
    u0 = u[indx]
    v0 = v[indx]
    x0 = x[indx]
    y0 = y[indx]

    # Now find the other points on the same chip.
    # The heuristic we use for this is that the distance in u,v coords is about
    # 0.2 * distance in x,y coords
    uvdistsq = (u - u0)**2 + (v - v0)**2
    xydistsq = (x - x0)**2 + (y - y0)**2
    select = (uvdistsq > 0.18**2 * xydistsq) & (uvdistsq < 0.22**2 * xydistsq)
    if np.sum(select) == 0:
        print('No points with the right distance.')
        return None, None
    u1 = u[select]
    v1 = v[select]
    x1 = x[select]
    y1 = y[select]

    # Calculate the rotation required to put the +u direction parallel to +x
    uv_theta = np.arctan2(v1-v0, u1-u0)
    xy_theta = np.arctan2(y1-y0, x1-x0)
    rot = xy_theta - uv_theta
    rot = rot % (2.*np.pi)
    # Median just in case there are a couple spurious matches to our distance criterion.
    # Shouldn't happen I don't think from an edge point, but be safe anyway.
    theta = np.median(rot)
    #print('rotation angle = ',theta)
    if np.isnan(theta):
        return None, None

    # Rotate (u,v) by theta -> (u_cam,v_cam)
    u_cam = u * np.cos(theta) - v * np.sin(theta)
    v_cam = u * np.sin(theta) + v * np.cos(theta)

    # Also rotate shears
    g1 = data['g1_model'] * np.cos(2*theta) - data['g2_model'] * np.sin(2*theta)
    g2 = data['g1_model'] * np.sin(2*theta) + data['g2_model'] * np.cos(2*theta)
    data['g1_model'] = g1
    data['g2_model'] = g2
    g1 = data['g1_data'] * np.cos(2*theta) - data['g2_data'] * np.sin(2*theta)
    g2 = data['g1_data'] * np.sin(2*theta) + data['g2_data'] * np.cos(2*theta)
    data['g1_data'] = g1
    data['g2_data'] = g2

    #print('u_cam range = ',np.min(u_cam),np.max(u_cam))
    #print('v_cam range = ',np.min(v_cam),np.max(v_cam))

    # Recenter the points now that we are oriented right.
    umin = np.min(u_cam)
    umax = np.max(u_cam)
    vmin = np.min(v_cam)
    vmax = np.max(v_cam)
    u_cam -= (umax+umin)/2
    v_cam -= (vmax+vmin)/2
    #print('u_cam range => ',np.min(u_cam),np.max(u_cam))
    #print('v_cam range => ',np.min(v_cam),np.max(v_cam))

    if (np.max(np.abs(u_cam)) > 6400
        or np.max(np.abs(v_cam)) > 6400
        or np.any((np.abs(u_cam) > 3900) & (np.abs(v_cam) > 3900))):
        print('u_cam range = ',np.min(u_cam),np.max(u_cam))
        print('v_cam range = ',np.min(v_cam),np.max(v_cam))

        # Try to correct for the possibility that the exposure is only partially done.
        # This happens when the exposure falls off the edge of the DC2 region.
        # The result is that range of u or v (or both) will be shifted.
        # Once the coordinates are rotated to camera orientation, then either the top or
        # bottom should have 3 rafts (9 ccds) along the edge.  The other side will typically
        # be cut at a diagonal, so less (or maybe cut along the long direction, so more).
        # Find the one that seems to be about the right size along the edge and shift it down 
        umin = np.min(u_cam)  # (recalculate these, since they changed.)
        umax = np.max(u_cam)
        vmin = np.min(v_cam)
        vmax = np.max(v_cam)
        left = u_cam < umin + 2400  # One raft worth
        right = u_cam > umax - 2400
        vrange_left = np.max(v_cam[left]) - np.min(v_cam[left])
        vrange_right = np.max(v_cam[right]) - np.min(v_cam[right])
        bottom = v_cam < vmin + 2400
        top = v_cam > vmax - 2400
        urange_bottom = np.max(u_cam[bottom]) - np.min(u_cam[bottom])
        urange_top = np.max(u_cam[top]) - np.min(u_cam[top])

        print('vrange on left/right = ',vrange_left,vrange_right)
        print('urange on top/bottom = ',urange_top,urange_bottom)
        # expected value is about 7580.  See who is farther off.

        if abs(vrange_left-7580) > max(abs(vrange_right-7580), 200):
            print('Left side seems wrong.  Shift right to expected value.')
            u_cam += 6320 - umax
        elif abs(vrange_right-7580) > max(abs(vrange_left-7580), 200):
            print('Right side seems wrong.  Shift left to expected value.')
            u_cam += -6320 - umin

        if abs(urange_bottom-7580) > max(abs(urange_top-7580), 200):
            print('Bottom side seems wrong.  Shift up to expected value.')
            v_cam += 6320 - vmax
        elif abs(urange_top-7580) > max(abs(urange_bottom-7580), 200):
            print('Top side seems wrong.  Shift down to expected value.')
            v_cam += -6320 - vmin
        print('u_cam range => ',np.min(u_cam),np.max(u_cam))
        print('v_cam range => ',np.min(v_cam),np.max(v_cam))

        if (np.max(np.abs(u_cam)) > 6400
            or np.max(np.abs(v_cam)) > 6400
            or np.any((np.abs(u_cam) > 3900) & (np.abs(v_cam) > 3900))):

            print()
            print(" *************** ")
            print(" Still a problem!")
            print(" *************** ")
            print()
            print(' U check: ',np.max(u_cam) > 6400, np.min(u_cam) < -6400)
            print(' V check: ',np.max(v_cam) > 6400, np.min(v_cam) < -6400)
            print(' Upper right? ',np.sum((u_cam > 3900) & (v_cam > 3900)))
            print(' Upper left? ',np.sum((u_cam < -3900) & (v_cam > 3900)))
            print(' Lower right? ',np.sum((u_cam > 3900) & (v_cam < -3900)))
            print(' Lower left? ',np.sum((u_cam < -3900) & (v_cam < -3900)))
            #ll = (u_cam < -3900) & (v_cam < -3900)
            #print('    u = ',u_cam[ll])
            #print('    v = ',v_cam[ll])
            print()
            return None, None

    return u_cam, v_cam

for band in ['u','g','r','i','z','y']:
    files = glob.glob('piff-run-0001/*/*-%s_cat.fits'%band)
    print('found %d cat files for %s band'%(len(files),band))

    all_data = []

    for f in files:
        print(f)
        try:
            data = fitsio.read(f)
        except OSError as e:
            print('Caught ',e)
            continue
        u_cam, v_cam = compute_uv_cam(data)
        if u_cam is None:
            print('Unable to figure out camera coordinates.  Using 1.e100 for u_cam, v_cam')
            u_cam = np.ones_like(data['u']) * 1.e100
            v_cam = np.ones_like(data['u']) * 1.e100
        dtype2 = np.dtype(data.dtype.descr + [('u_cam','>f8'), ('v_cam','>f8')])
        data2 = np.empty(len(data), dtype=dtype2)
        for key in data.dtype.names:
            data2[key] = data[key]
        data2['u_cam'] = u_cam
        data2['v_cam'] = v_cam
        all_data.append(data2)

    print('concatenating...')
    data = np.concatenate(all_data)
    print('dtype = ',data.dtype)

    file_name = 'psf_%s.fits'%band
    print('writing to',file_name)
    fitsio.write(file_name, data, clobber=True)


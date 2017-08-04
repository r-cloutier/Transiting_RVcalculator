from imports import *


def WM14(rps):
    '''
    Compute the planetary mass from input planetary radii based on the mean 
    mass-radius relation for small planets (radius <= 4 Earth radii) from 
    Weiss & Marcy 2014 (Eqs. 1-3).

    Parameters
    ----------
    `rps': float or array-like
        Value(s) of the planetary radii that will be converted to a planetary 
        mass. Units of the input radii are assumed to be in Earth radii.

    Returns
    -------
    `mps': numpy-array
        Planetary masses in Earth masses corresponding to the planetary radii 
        given in `rps'.
    '''
    rps = np.ascontiguousarray(rps)
    maxE, maxN = 1.5, 4.

    # Check
    if np.any(rps >= maxN):
	raise ValueError('WM14 mass-radius relation is only valid for ' + \
			 'planets with radii < %i Earth radii.'%maxN)

    # Compute mean mass for each planet radius
    mps = np.zeros(rps.size)
    for i in range(rps.size):

        if rps[i] < maxE:
            mps[i] = .441 * rps[i]**3 + .615 * rps[i]**4

        else:
            mps[i] = 2.69 * rps[i]**(.93)

    return mps



def WM14_upperlimit(rps):
    '''
    Compute the 1 sigma upper limit on the planetary mass from input planetary 
    radii based on the mean mass-radius relation plus rms for small planets 
    (radius <= 4 Earth radii) from Weiss & Marcy 2014 (Eqs. 1-3).

    Parameters
    ----------
    `rps': float or array-like
        Value(s) of the planetary radii that will be converted to a planetary
        mass. Units of the input radii are assumed to be in Earth radii.

    Returns
    -------
    `mps': numpy-array
        Planetary masses in Earth masses corresponding to the planetary radii
        given in `rps'.
    '''
    # Compute mean mass
    rps = np.ascontiguousarray(rps)
    mps = WM14(rps)

    # Add 1 sigma to mean masses to get upper limits
    maxE, rmsE, rmsN = 1.5, 2.7, 4.7
    for i in range(mps.size):

        if rps[i] < maxE:
            mps[i] += rmsE

        else:
            mps[i] += rmsN

    return mps



def WM14_lowerlimit(rps):
    '''
    Compute the 1 sigma lower limit on the planetary mass from input planetary 
    radii based on the mean mass-radius relation minus rms for small planets 
    (radius <= 4 Earth radii) from Weiss & Marcy 2014 (Eqs. 1-3). Negative  
    masses are not permitted. 

    Parameters
    ----------
    `rps': float or array-like
        Value(s) of the planetary radii that will be converted to a planetary
        mass. Units of the input radii are assumed to be in Earth radii.

    Returns
    -------
    `mps': numpy-array
        Planetary masses in Earth masses corresponding to the planetary radii
        given in `rps'.
    '''
    # Compute mean mass
    rps = np.ascontiguousarray(rps)
    mps = WM14(rps)

    
    # Subtract 1 sigma from mean masses to get lower limits
    maxE, rmsE, rmsN = 1.5, 2.7, 4.7
    for i in range(mps.size):

        if rps[i] < maxE:
            mps[i] -= rmsE

        else:
            mps[i] -= rmsN

    # Assign small masses to negative mass planets
    if np.any(mps < 0):
        low_mp = 1e-10
        warnings.warn('\nNegative mass planets are being set to ' + \
                      '%.0e Earth masses.'%low_mp, Warning)
        mps[mps < 0] = low_mp
            
    return mps

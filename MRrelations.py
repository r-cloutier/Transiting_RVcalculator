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
    maxE, maxN = 1.5, 4.
    
    rps = np.array(rps)
    assert not np.any(rps > maxN)

    mps = np.zeros(rps.size)
    for i in range(rps.size):

        if rps[i] < maxE:
            mps[i] = .44 * rps[i]**3 + .614 * rps[i]**4

        else:
            mps[i] = 2.69 * rps[i]**(.93)

    return mps

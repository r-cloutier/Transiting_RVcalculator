from imports import *

# Conversion functions
def days2sec(t):
    return t * 24. * 60 * 60
def Msun2kg(m):
    return m * 1.989e30
def Mearth2kg(m):
    return m * 5.972e24


def RV_K(P, Ms, Mp, ecc=0., inc=90.):
    '''
    Compute the RV semi-amplitude of an input planet.

    Parameters
    ----------
    `P': scalar or uncertainty array
        Orbital period [in days]
    `Ms': scalar or uncertainty array
        Stellar mass [in Solar masses]
    `Mp': scalar or uncertainty array
        Planet mass [in Earth masses]
    `ecc': scalar or uncertainty array
        Orbital eccentricity
    `inc': scalar or uncertainty array
        Orbital inclination relative to the plane of the sky (e.g 90 deg for 
        edge-on orbits) [in degrees]

    Returns
    -------
    `K': uncertainty array
        RV semi-amplitude [in m/s]
    '''
    P, Ms, Mp, inc = days2sec(P), Msun2kg(Ms), Mearth2kg(Mp), \
                     unumpy.radians(inc)
    G = 6.67408e-11
    return (2*np.pi*G/(P*Ms*Ms))**(1./3) * \
        Mp*unumpy.sin(inc) / unumpy.sqrt(1-ecc**2)

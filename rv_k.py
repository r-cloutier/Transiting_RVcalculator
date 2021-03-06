from imports import *

global G
G = 6.67408e-11

# Conversion functions
def days2sec(t):
    return t * 24. * 60 * 60
def Msun2kg(m):
    return m * 1.98849925145e30
def Mearth2kg(m):
    return m * 6.04589804468e24
def kg2Mearth(m):
    return m / 6.04589804468e24


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
                     unp.radians(inc)
    return (2*np.pi*G/(P*Ms*Ms))**(1./3) * \
        Mp*unp.sin(inc) / unp.sqrt(1-ecc**2)


def RV_mp(P, Ms, K, ecc=0., inc=90.):
    '''
    Compute the planet mass an input planet given its RV semi-amplitude.

    Parameters
    ----------
    `P': scalar or uncertainty array
        Orbital period [in days]
    `Ms': scalar or uncertainty array
        Stellar mass [in Solar masses]
    `K': scalar or uncertainty array
        RV semi-apmplitude [in m/s]
    `ecc': scalar or uncertainty array
        Orbital eccentricity
    `inc': scalar or uncertainty array
        Orbital inclination relative to the plane of the sky (e.g 90 deg for 
        edge-on orbits) [in degrees]

    Returns
    -------
    `mp': uncertainty array
        Planet mass [in Earth masses]
    '''
    P, Ms, inc = days2sec(P), Msun2kg(Ms), unp.radians(inc)
    mp = K * (P*Ms*Ms/(2*np.pi*G))**(1./3) * \
         unp.sqrt(1-ecc**2)/unp.sin(inc)    
    return kg2Mearth(mp)

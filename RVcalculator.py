from imports import *
from MRrelations import *

global c
c = 3e8

class RVcalculator:

    def __init__(self, startheta, planettheta, MRfunc,
                 detsigs=[3,5], additiveRVjitter=0, texp=0,
                 texpmin=5, texpmax=60, bands=['Y','J','H','K'], R=75e3,
                 aperture=3.58, efficiency=.15, SNRtarget=150):
        '''
        Given the known stellar parameters, planet transit parameters, and 
        various parameters of the velocimeter, estimate the number of RV 
        measurements required to detect the transiting planet at a given 
        detection significance.

        Parameters
        ----------
        `startheta': list containing [list, scalar, scalar, scalar, scalar] (5,)
            List of stellar parameters with the following 4 entries:
            `mags': list of broadband apparent magnitudes for each band used 
                    by the velocimeter to measure stellar RVs
            `SpT': numerical M dwarf spectral type (e.g. for an M4.5 dwarf, 
                   SpT = 4.5)  
            `vsini': the projected stellar rotation [in km/s]
            `Ms': the most-likely host stellar mass [in Solar masses]
            `sigMs': the 1sigma uncertainty on the host stellar mass 
                     [in Solar masses]
        `planettheta': list containing [scalar, scalar, scalar] (3,)
            List of transit planet parameters with the following 3 entries:
            `P': orbital period [in days]
            `sigP': 1 sigma uncertainty on the orbital period [in days]
            `rp': planet radius [in Earth radii]
        `MRfunc': callable
            Mass-radius relation function which takes input planet radii and 
            returns the corresponding planet masses. Needed to estimate the 
            RV semi-amplitude of transiting planets.
        `detsigs': scalar or array-like
            The desired planet mass detection signficances (e.g. 3 for 3 sigma)
        `additiveRVjitter': scalar
            An extra RV jitter term to be added to the effective RV uncertainty 
            in quadrature [in m/s]
        `texp': scalar
            The spectroscopic exposure time [in minutes]. If texp=0 then the 
            exposure time is calculated based on a target signal-to-noise 
            ratio. If texp!=0 then it is fixed to that value.
        `texpmin': scalar
            The minimum permitted exposure time [in minutes].
        `texpmax': scalar
            The maximum permitted exposure time [in minutes].
        `bands': list of strings (Nmags,)
            Names of the velocimeter's broadband photometric bands for which 
            the star's apparent magnitudes have been specified. Must have the 
            same number of entries as `mags' in `startheta'.
        `R': scalar
            Spectral resolution of the velocimeter.
        `aperture': scalar
            Diameter of the telescope's primary [in meters].
        `efficiency': scalar
            Fractional detector efficiency. Must be between 0-1.
        `SNRtarget': scalar
            The target signal-to-noise ratio per resolution element used to 
            estimate the exposure time unless `texp' is specified by the user.
        '''
        # Check
        mags, self.SpT, self.vsini, Ms, sigMs = startheta
        self.Ms = unumpy.uarray(Ms, sigMs)
        self.mags, self.bands = np.asarray(mags), np.asarray(bands)
        if mags.size != bands.size:
            raise ValueError('Do not have the same number of magnitudes ' + \
                             'as bands.')

        # Compute exposure time if not specified
        self.SNRtarget, self.R = SNRtarget, R
        self.texpmin, self.texpmax = texpmin, texpmax
        self.aperture, self.efficiency = aperture, efficiency
        self.element_res = c / self.R
        self.texp = self._compute_texp() if texp == 0 else texp



        
    def _compute_texp(self):
            

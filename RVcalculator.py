from imports import *
from MRrelations import *

global c
c = 299792458.


class RVcalculator:

    def __init__(self, startheta, planettheta, MRfunc,
                 detsigs=[3,5], additiveRVjitter=0, texp=0,
                 texpmin=1, texpmax=60, bands=['Y','J','H','K'], R=75e3,
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

        Example
        -------
        >>> from MRrelations import WM14
        >>> mags, SpT, vsini, Ms, sigMs = [10.,9.2,8.7,8.3], 4.5, .1, .181, .019
        >>> P, sigP, rp = 1.629, 3e-5, 1.16
        >>> self = RVcalculator((mags,SpT,vsini,Ms,sigMs), (P,sigP,rp), WM14)
        '''
        
        # Check
        mags, self.SpT, self.vsini, Ms, sigMs = startheta
        self.Ms = unumpy.uarray(Ms, sigMs)
        self.mags, self.bands = np.asarray(mags), np.asarray(bands)
        if self.mags.size != self.bands.size:
            raise ValueError('Do not have the same number of magnitudes ' + \
                             'as bands.')

        # Compute exposure time if not specified
        self.SNRtarget, self.R = SNRtarget, R
        self.texpmin, self.texpmax = texpmin, texpmax
        if self.texpmin <= 0:
            raise ValueError('texpmin must be greater than 0.')
        self.aperture, self.efficiency = aperture, efficiency
        self.element_res = c / self.R
        if texp == 0:
            self._compute_texp()
        else:
            self.texp = texp

        # Compute expected RV measurement uncertainty
        #self._compute_sigmaRV()


        
    def _compute_texp(self):
        '''
        Compute the exposure time required to achieve a desired SNR within an 
        allowable range of exposure times.
        '''
        # Compute sigmaRV resulting from different assumed exposure times
        texps = np.arange(self.texpmin, self.texpmax+.25, .25)
        SNRs = np.zeros(texps.size)
        for i in range(texps.size):
            snrs_tmp = np.zeros(self.mags.size)
            for j in range(self.mags.size):
                snrs_tmp[j] = self._get_snr(self.bands[j], self.mags[j],
                                            texps[i])
            SNRs[i] = np.median(snrs_tmp)

        # Return texp that gets us closest to the desired SNRtarget
        g = abs(SNRs-self.SNRtarget) == abs(SNRs-self.SNRtarget).min()
        self.texp = texps[g]



    def _get_snr(self, band, mag, texp):
        '''
        Compute the SNR of the spectrum in a certain band (e.g. 'J').
        
        Parameters
        ----------
        `band': string
            Name of the spectral band. Must be one of 'Y','J','H', or 'K'. 
        `mag': scalar
            Stellar apparent magnitude in `band'.
        `texp': scalar
            Exposure time of the spectroscopic observation [in minutes].

        Returns
        -------
        `SNR': float
            The signal-to-noise of the spectroscopic observation of the star. 
        '''
        # Define constants
        area = np.pi * (self.aperture*1e2/2)**2   # in cm^2
        resele = float(self.element_res)
        texp_s = texp*60.     # in sec
        c_As  = c * 1e10      # angstroms/s
        c_kms = c_As / 1e13   # km/s
        h = 6.62607004e-27    # ergs s
        
        # Get flux density zeropoint (for m=0) in ergs s^-1 A^-1 cm^-2,
        # wavelength in microns and bandwidth in microns
        if band == 'Y':
            Fl = 6.063e-10
            l  = 10174.
        elif band == 'J':
            Fl = 3.143e-10
            l  = 12350.
        elif band == 'H':
            Fl = 1.144e-10
            l  = 16620.
        elif band == 'K':
            Fl = 4.306e-11
            l  = 21590.
        else:
            raise ValueError('The band %s has an undefined zeropoint'%band)

        # Get Fl of the star in ergs s^-1 A^-1 cm^-2
        fl = Fl * 10**(-.4 * mag)

        # Get # of photons per angstrom
        Ephot = h*c_As/l
        Nphot_A = fl * texp_s * area * self.efficiency / Ephot

        # Get # photons per resolution element (dl/l)
        dl_l = c_kms / resele  # dimensionless (75000)
        Nphot_res = Nphot_A * (l / dl_l)

        # Add photon noise and dark current to noise budget
        darkcurrent, footprint = 1e-2, 12     # electrons/s, pixels
        SNR = Nphot_res / np.sqrt(Nphot_res + darkcurrent*footprint*texp_s)

        return SNR



    def _compute_sigmaRV(self):
        '''
        Compute the expected RV measurement uncertainty when observing a star 
        spectroscopically in a particular set of bands.
        '''
        # Get SNR in each passband of interest
        SNRs = np.zeros(self.bands.size)
        for i in range(self.bands.size):
            print self.bands[i], self.mags[i], self.texp
            SNRs[i] = self._get_snr(self.bands[i], self.mags[i], self.texp)

        # Get sigmaRV from Figueira+2016 tables
        self.SNRs, self.SNR = SNRs, np.median(SNRs)
        self.sigmaRV = self._rvcontent(self.SNR)


    def _rvcontent(self, SNR):
        '''
        Estimate the RV measurement uncertainty given a spectral type (velocity 
        information), vsini (broadening), spectral resolution, bands, and SNR 
        per resolution element (see get_snr()).

        Parameters
        ----------
        `SNR' = scalar
            The signal-to-noise per resolution element obtained throughout the 
            observations (see self._get_snr()).
        '''
        # Check parameter values for interpolation
        if self.SpT < 0:
            self.SpT = 0
        elif self.SpT > 9:
            self.SpT = 9.
        if self.vsini < 1:
            self.vsini = 1.
        elif self.vsini > 10:
            self.vsini = 10.
        if self.R < 6e4:
            self.R = 6e4
        elif self.R > 1e5:
            self.R = 1e5

        # Get grid of data
        d = np.genfromtxt('input_data/rvcontent_Figueira16.dat', dtype=str)
        spt_t, sigmarv_t = d[:,0].astype(float), d[:,4].astype(float)
        vsini_t, R_t = d[:,2].astype(float), d[:,3].astype(float)
        band_t = d[:,1]
        
        # Get range of resolutions
        if self.R not in R_t:
            Rs = np.array([np.floor(self.R/2e4)*2e4, np.ceil(self.R/2e4)*2e4])
        else:
            Rs = np.ascontiguousarray(self.R)

        # Get SNR in each band
        nbands, nR = self.bands.size, Rs.size
        sigmarvs = np.zeros((nbands, nR))
        rescoeffs = np.zeros(nR)  # for weighted average
        for i in range(nbands):
            for j in range(nR):
                thisband = np.logical_and(band_t==self.bands[i], R_t==Rs[j])
                spt_band     = spt_t[thisband]
                vsini_band   = vsini_t[thisband]
                sigmarv_band = sigmarv_t[thisband]

                # Interpolate along spt and vsini
                lintsigmarv = interp2d(spt_band, vsini_band, sigmarv_band)
                sigmarvs[i,j] = lintsigmarv(self.SpT, self.vsini)
                
                # Get resolution coefficients for weighted mean
                if nR > 1:
                    rescoeffs[j] = abs(1 - abs(self.R - Rs[j]) / np.diff(Rs))

        # Scale RVs because the RV accuracy in Figueira+2016 was set at SNR=100
        sigmarvs *= np.sqrt(1e2/SNR)

        # Combine results of different resolutions with weighted mean
        if nres > 1:
            sigmarvs2 = np.zeros(nbands)
            for i in range(nbands):
                sigmarvs2[i] = sigmarvs[i][0] * rescoeffs[0] + \
                               sigmarvs[i][1] * rescoeffs[1]
        else:
            sigmarvs2 = sigmarvs

        # Combine all bands and floor at a 1 m/s stability
        return 1./np.sqrt(np.sum(1./sigmarvs2**2))


# TEMP
if __name__ == '__main__':
    mags, SpT, vsini, Ms, sigMs = [10.,9.2,8.7,8.3], 4.5, .1, .181, .019
    P, sigP, rp = 1.629, 3e-5, 1.16
    self = RVcalculator((mags,SpT,vsini,Ms,sigMs), (P,sigP,rp), WM14)

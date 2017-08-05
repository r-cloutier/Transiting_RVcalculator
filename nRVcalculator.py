from imports import *
from rv_k import *

global c
c = 299792458.


class nRVcalculator:

    def __init__(self, startheta, planettheta, MRfunc,
                 detsigs=[3,5], sigmaRV=0, bands=['Y','J','H','K'], texp=0,
                 texpmin=1, texpmax=60, additiveRVjitter=0, R=75e3,
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
        `sigmaRV': scalar
            The (constant) RV measurement uncertainty [in m/s]. If sigmaRV=0 
            then sigmaRV is estimated from the stellar and instrument 
            parameters. If sigmaRV!=0 then it is fixed to that value. 
        `bands': list of strings (Nmags,)
            Names of the velocimeter's broadband photometric bands for which 
            the star's apparent magnitudes have been specified. Must have the 
            same number of entries as `mags' in `startheta'.
        `texp': scalar
            The spectroscopic exposure time [in minutes]. If texp=0 then the 
            exposure time is calculated based on a target signal-to-noise 
            ratio. If texp!=0 then it is fixed to that value.
        `texpmin': scalar
            The minimum permitted exposure time [in minutes].
        `texpmax': scalar
            The maximum permitted exposure time [in minutes].
        `additiveRVjitter': scalar
            An extra RV jitter term to be added to the effective RV uncertainty 
            in quadrature [in m/s]
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
        # Define stellar parameters
        mags, self.SpT, self.vsini, Ms, sigMs = startheta
        self.Ms = unp.uarray(Ms, sigMs)
        
        # Check bands and mags
        self.mags, self.bands = np.ascontiguousarray(mags), np.ascontiguousarray(bands)
        if sigmaRV == 0 and texp == 0:
            if self.mags.size != self.bands.size:
                raise ValueError('Do not have the same number of ' + \
                                 'magnitudes as bands.')

        # Compute exposure time if not specified
        self.SNRtarget, self.R = SNRtarget, R
        self.texpmin, self.texpmax = texpmin, texpmax
        if self.texpmin <= 0:
            raise ValueError('texpmin must be greater than 0.')
        self.aperture, self.efficiency = aperture, efficiency
        self.element_res = c*1e-3 / self.R
        if texp == 0 and sigmaRV == 0:
            self._compute_texp()
        else:
            self.texp = texp

        # Compute expected RV measurement uncertainty if not specified
        self.sigmaRV, self.additiveRVjitter = sigmaRV, additiveRVjitter
        if self.sigmaRV == 0:
            self._compute_sigmaRV()
        else:
            self.sigmaRV, self.sigmaRV_eff = np.repeat(sigmaRV, 2)
        self.sigmaRV_eff = np.sqrt(self.sigmaRV_eff**2+self.additiveRVjitter**2)
        
        # Get planet parameters
        self.planettheta, self.MRfunc = planettheta, MRfunc
        self._define_planet()

        # Estimate the number of RV visits required to measure the planet
        # mass at a desired detection significance
        self.detsigs = np.ascontiguousarray(detsigs)
        self._estimate_nRVs()
        

    def _define_planet(self):
        '''
        Define planet parameters including period, radius, mass, and RV 
        semi-amplitude.
        '''
        P, sigP, self.rp = self.planettheta
        self.P = unp.uarray(P, sigP)
        self.mp = float(self.MRfunc(self.rp))
        K = unp.nominal_values(RV_K(unp.nominal_values(self.P),
                                       unp.nominal_values(self.Ms),
                                       self.mp))
        self.K = float(K)
        del self.planettheta

        Ktmp = unp.nominal_values(self.K) 
        if Ktmp <= self.sigmaRV_eff:
            message = '\nExpected K (%.3f m/s) is less than the '%Ktmp + \
                      'effective RV measurement uncertainty ' + \
                      '(%.3f m/s). Planet will be '%self.sigmaRV_eff + \
                      'difficult to detect.'
            warnings.warn(message, Warning)

        
        
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
        self.texp = float(texps[g])



    def _get_snr(self, band, mag, texp):  # WORKING
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
        Ephot = h * c_As / l
        Nphot_A = fl * texp_s * area * self.efficiency / Ephot

        # Get # photons per resolution element (dl/l)
        Nphot_res = Nphot_A * (l / self.R)

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
            SNRs[i] = self._get_snr(self.bands[i], self.mags[i], self.texp)

        # Get sigmaRV from Figueira+2016 tables
        self.SNRs, self.SNR = SNRs, np.median(SNRs)
        self.sigmaRV, self.sigmaRV_eff = self._rvcontent(self.SNR)
        

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
            SpT = 0
        elif self.SpT > 9:
            SpT = 9.
        if self.vsini < 1:
            vsini = 1.
        elif self.vsini > 10:
            vsini = 10.
        if self.R < 6e4:
            R = 6e4
        elif self.R > 1e5:
            R = 1e5
        else:
            R = self.R

        # Get grid of data
        d = np.genfromtxt('input_data/rvcontent_Figueira16.dat', dtype=str)
        spt_t, sigmarv_t = d[:,0].astype(float), d[:,4].astype(float)
        vsini_t, R_t = d[:,2].astype(float), d[:,3].astype(float)
        band_t = d[:,1]
        
        # Get range of resolutions
        if R not in R_t:
            Rs = np.array([np.floor(R/2e4)*2e4, np.ceil(R/2e4)*2e4])
        else:
            Rs = np.ascontiguousarray(R)

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
                    rescoeffs[j] = abs(1 - abs(R - Rs[j]) / np.diff(Rs))

        # Scale RVs because the RV accuracy in Figueira+2016 was set at SNR=100
        sigmarvs *= np.sqrt(1e2/SNR)

        # Combine results of different resolutions with weighted mean
        if nR > 1:
            sigmarvs2 = np.zeros(nbands)
            for i in range(nbands):
                sigmarvs2[i] = sigmarvs[i][0] * rescoeffs[0] + \
                               sigmarvs[i][1] * rescoeffs[1]
        else:
            sigmarvs2 = sigmarvs

        # Combine all bands and floor at a 1 m/s stability
        sigmaRV = 1./np.sqrt(np.sum(1./sigmarvs2**2))
        sigmaRV_eff = sigmaRV if sigmaRV >= 1 else 1.
        return sigmaRV, sigmaRV_eff


    
    def _estimate_nRVs(self):
        '''
        Estimate the number of RV measurements required to detect the planet 
        mass at a given detection significance.
        '''
        self.Ntest, self.mpdetsigs = np.zeros((self.detsigs.size, 100)), \
                                     np.zeros((self.detsigs.size, 100))
        self.nRVs = np.zeros(self.detsigs.size, dtype=np.int)
        for i in range(self.detsigs.size):

            # Compute detection significance for increasing nRVs until
            # the desired value is obtained
            gotdetsig, moreN = False, 0
            while not gotdetsig:
                
                # Compute sigma_K as a function of N
                self.Ntest[i] = np.arange(1+moreN, 101+moreN)
                sigKs = self._compute_sigK(self.Ntest[i])
                Ks = unp.uarray(self.K, sigKs)
                
                # Compute corresponding sigma_mp as a function of N
                mps = RV_mp(self.P, self.Ms, Ks)
                self.mpdetsigs[i] = unp.nominal_values(mps) / \
                                    unp.std_devs(mps)
                
                # Did we get the desired detsig?
                if self.mpdetsigs[i].min() <= self.detsigs[i] <= \
                   self.mpdetsigs[i].max():
                    g = abs(self.mpdetsigs[i]-self.detsigs[i]) == \
                        abs(self.mpdetsigs[i]-self.detsigs[i]).min()
                    self.nRVs[i] = self.Ntest[i][g][0]
                    gotdetsig = True

                elif self.detsigs[i] < self.mpdetsigs[i].min():
                    warnings.warn("It is too easy to measure this planet's " + \
                                  'mass at a detection significance of ' + \
                                  '%.1f sigma.'%self.detsig, Warning)

                else:
                    moreN += 25



    def _compute_sigK(self, nRV_arr):
        '''
        Estimate the measurement uncertainty on K analytically from the Fisher 
        information assuming a single transiting planet on a circular orbit.

        Parameters
        ----------
        `nRV_arr': array-like
            Array of numbers of RV measurements.

        Returns
        -------
        `sigKs': numpy-array (nRV,)
            Array of K measurement uncertainties as a function of `nRV_arr'. 
        '''
        return self.sigmaRV_eff * np.sqrt(2./nRV_arr)


    def report_results(self):
        '''
        Print the results to screen.
        '''
        Ms, sigMs = unp.nominal_values(self.Ms), unp.std_devs(self.Ms)
        P, sigP = unp.nominal_values(self.P), unp.std_devs(self.P)
        print '#'*83, '\n## Stellar Parameters:'
        for i in range(self.mags.size):
            print '## %s\t=\t%.2f'%(self.bands[i], self.mags[i])
        print '## Spectral type\t=\tM%.1f'%self.SpT
        print '## Stellar mass\t=\t%.3f +- %.3f Solar masses'%(Ms, sigMs)
        print '## Projected rotation velocity\t=\t%.2f km/s'%self.vsini
        print '##\n## Planet Parameters:'
        print '## Orbital period\t=\t%.5f +- %.1e days'%(P, sigP)
        print '## Planet radius\t=\t%.2f Earth radii'%self.rp
        print '## Planet mass\t=\t%.2f Earth masses'%self.mp
        print '## RV semi-amplitude\t=\t%.2f m/s\n##'%self.K
        if self.texp != 0:
            print '## Exposure time\t=\t%.3f minutes'%self.texp
        print '## Effective RV uncertainty\t=\t %.3f m/s'%self.sigmaRV_eff
        print '##\n## Results:'
        for i in range(self.detsigs.size):
            print '## %.2f sigma detection significance '%self.detsigs[i] + \
                'requires %i uniformly sampled RV measurements.'%self.nRVs[i]
        print '#'*83
        

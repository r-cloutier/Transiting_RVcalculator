'''
Given a transiting planet candidate, compute the number of visits and total on-sky time 
required to recover the planet mass.
'''
import numpy as np
import rvs
from uncertainties import unumpy as unp
from scipy.interpolate import interp2d


def compute_Nvisits_1circularplanet(startheta, planettheta, MRfunc, detsigs, texp=0,
                                    addjitter2sigmaRV=0, bands=['Y','J','H','K'], res=75e3,
                                    tmin=5., tmax=60., aperture=3.58, efficiency=.15,
                                    element_res=4., SNRtarget=150):
    '''Estimate the number of visits and total on-sky time required to recover the 
    mass of a single transiting planet to a certain precision (i.e. 3 sigma) assuming 
    that the planet is 
    i) on a nearly circular orbit, 
    ii) dominates the RV signal, and
    iii) has its orbit well-sampled by the RV measurements (prbly need at least 3 
    measurements at T0 and T0 +- P/2).

    startheta = mags, SpT, vsini, unMs
    mags = list of near-IR magnitudes (must be the same length as bands)
    SpT = scalar spectral type
    vsini = projected stellar rotation velocity (km/s)
    unMs = stellar mass and uncertainty (Msun); to estimate mp from K

    planettheta = P, rp
    P = orbtial period (days)
    rp = planet radius (Rearth)

    To compute the planet mass measurement precision we need to estimate 
    the planet mass using an assumed mass-radius relation given by the 
    function specified in MRfunc (e.g. MRfunc=mr_WM14).

    detsigs is a list of desired detection significances e.g. [3,5]

    Estimating sigmaRV:
    Using the near-IR magnitudes, spectral type, and vsini we can estimate the 
    RV measurement uncertainty which affects how well the planet's 
    mass can be measured. Then set addjitter2sigmaRV at the resiudal jitter from additional 
    planets and/or jitter to be added to the RV measurement uncertainty. This gives the 
    effective RV uncertainty.

    Setting texp to a non-zero value in minutes fixes the exposure time to the value of 
    texp. Otherwise we try to calculate the minimum value of texp to achieve 1 m/s stability 
    but never exceeding 1 hour exposures.

    bands specifies which near-IR bands are accessed by the spectrograph ('Y','J','H','K' for SPIRou)
    res is the spectral resolution of the spectrograph (75,000 for SPIRou)
    tmin, tmax = are the minimum and maximum permitted integration times in minutes
    aperture = the diameter of the telescope in meters
    efficiency = detector efficiency in [0,1]
    element_res = resolution element size in km/s
    SNRtarget = desired SNR of the spectrum (used to determine the optimal integration time)
    '''
    mags, SpT, vsini, Ms = startheta
    if len(mags) != len(bands):
        raise ValueError('Do not have the same number of magnitudes as bands.')

    # Compute exposure time if not specified
    if texp == 0:
        texp = compute_texp(mags, SpT, vsini, SNRtarget=SNRtarget, bands=bands, res=res,
                            tmin=tmin, tmax=tmax, aperture=aperture, efficiency=efficiency,
                            element_res=element_res)
    
    # Compute the expected RV measurement uncertainty
    sigmaRV, SNR = compute_sigmaRV(mags, SpT, vsini, bands=bands, res=float(res), texp=texp)
    print '\n', '#'*90
    print 'Target SNR = %.3f\nSNR achieved (with %.3f minute integration) = %.3f'%(SNRtarget, texp, SNR)
    effsigmaRV = np.sqrt(sigmaRV**2 + addjitter2sigmaRV**2)
    
    # Compute the planet mass from rp and K
    P, rp = planettheta
    mp = MRfunc(rp)
    K = rvs.RV_K(P, unp.nominal_values(Ms), mp)

    # Check
    if K <= sigmaRV:
        print "\n***WARNING*** K (%.3f m/s) is less than sigmaRV (%.3f m/s). Difficult to detect this planet."%(K,sigmaRV)
    
    # Assuming a single planet (circular) system, analytically compute the number
    # of visits required to achieve a mass detection at detsig sigma
    N_sigs = [Nvisits_to_get_sigma_1circularplanet(detsig, effsigmaRV, unp.nominal_values(Ms),
                                                   unp.std_devs(Ms), P, K) for detsig in detsigs]
    print_conclusions_1circularplanet(effsigmaRV, sigmaRV, K, detsigs, N_sigs, texp)
    print '#'*90


def compute_sigmaRV(mags, spt, vsini, bands=['Y','J','H','K'], res=75e3, texp=15.):
    '''Compute the expected sigmaRV when observing a star in particular set of bands.'''
    if len(mags) != len(bands):
        raise ValueError('Do not have the same number of magnitudes as bands.')
    # Get SNR in each passband of interest
    SNRs = np.zeros(len(bands))
    for i in range(len(bands)):
        SNRs[i] = get_snr(mags[i], bands[i], texp)
    # Get sigmaRV from Figueira+2016
    return rvcontent(spt, vsini, res, bands, np.median(SNRs)), np.median(SNRs)


def compute_texp(mags, spt, vsini, SNRtarget=150, bands=['Y','J','H','K'],
                 res=75e3, tmin=5., tmax=60., aperture=3.58, efficiency=.15, element_res=4):
    '''Compute the integration time (in minutes) required to achieve a desired SNR within 
    a specified range of times nominally 5 min < texp < 60 min.'''
    # Compute sigmaRV resulting from different exposure times
    texps = np.arange(tmin, tmax+.25, .25)
    SNRs = np.zeros(texps.size)
    for i in range(texps.size):
        snrs_tmp = np.zeros(len(mags))
        for j in range(len(mags)):
            snrs_tmp[j] = get_snr(mags[j], bands[j], texps[i], aperture=aperture,
                                  efficiency=efficiency, element_res=element_res)
        SNRs[i] = np.median(snrs_tmp)
    # Get texp that gets us closest to SNRtarget
    return texps[abs(SNRs - SNRtarget) == abs(SNRs - SNRtarget).min()]


def get_snr(mag, band, texp, aperture=3.58, efficiency=.15,
            element_res=4):
    '''Compute the SNR of the spectrum from the apparent magnitude of 
    the star in a certain band (e.g. 'J'), the exposure time in 
    minutes, the aperture of the telescope in meters, detector efficiency, 
    the element resolution in km/s.'''
    # Define constants
    area = np.pi * (aperture*1e2/2)**2   # CFHT aperature area in cm^2
    resele = float(element_res)          # element resolution is 4 km/s
    texp_s = texp*60.      # exposure time in sec
    c_As  = 299792458e10   # angstroms/s
    c_kms = c_As / 1e13    # km/s
    h = 6.62607004e-27     # ergs s
    # Get flux density zeropoint (for m=0) in ergs s^-1 A^-1 cm^-2,
    # wavelength in microns and bandwidth in microns
    if band == 'Y':
        Fl = 6.063e-10  # ergs s^-1 cm^-2 A^-1
        l  = 10174.     # angstroms
    elif band == 'J':
        Fl = 3.143e-10  # ergs s^-1 cm^-2 A^-1
        l  = 12350.     # angstroms
    elif band == 'H':
        Fl = 1.144e-10  # ergs s^-1 cm^-2 A^-1
        l  = 16620.     # agstroms
    elif band == 'K':
        Fl = 4.306e-11  # ergs s^-1 cm^-2 A^-1
        l  = 21590.     # agstroms
    else:
        raise ValueError('Not sure which band to use')
    # Get Fl of the star in ergs s^-1 A^-1 cm^-2
    fl = Fl * 10**(-.4 * mag)
    # Get # of photons per angstrom
    Ephot = h*c_As/l
    Nphot_A = fl * texp_s * area * efficiency / Ephot
    # Get # photons per resolution element (dl/l)
    dl_l = c_kms / resele  # dimensionless (75000)
    Nphot_res = Nphot_A * (l / dl_l)
    # Add photon noise and dark current to noise budget
    darkcurrent, footprint = 1e-2, 12     # electrons/s, pixels
    SNR = Nphot_res / np.sqrt(Nphot_res + darkcurrent*footprint*texp_s)
    return SNR


def rvcontent(spt, vsini, res, bands, SNR):
    '''Get the RV precision given a spectral type (lines), vsini (broadening),
    resolution, bands, and SNR per resolution element (see get_snr())
    spt = numerical spectral type (0-9)
    vsini = must be in [1, 10] km/s
    res = spectral resolution in [6e4, 1e5]
    bands = array of band names (e.g. ['Y','J','H','K'])
    SNR = from stellar magnitude (see get_SNR function)'''
    # Check parameter values
    if spt < 0:
        spt = 0
    elif spt > 9:
        spt = 9.
    if vsini < 1:
        vsini = 1.
    elif vsini > 10:
        vsini = 10.
    if res < 6e4 or res > 1e5:
        raise ValueError('Cannot extrapolate spectral resolution')

    # Get grid of data
    d = np.genfromtxt('input_data/rvcontent.dat', dtype=str)
    spt_t, sigmarv_t = d[:,0].astype(float), d[:,4].astype(float)
    vsini_t, res_t = d[:,2].astype(float), d[:,3].astype(float)
    band_t = d[:,1]

        # Get range of resolutions
    if res not in res_t:
        ress = np.array([np.floor(res/2e4)*2e4,
                         np.ceil(res/2e4)*2e4])
    else:
        ress = np.array([res])

    # Get SNR in each band
    nbands, nres = len(bands), ress.size
    sigmarvs = np.zeros((nbands, nres))
    rescoeffs = np.zeros(nres)  # for weighted average
    for i in range(nbands):
        for j in range(nres):
            thisband = np.where(np.logical_and(band_t == bands[i],
                                               res_t == ress[j]))[0]
            spt_band     = spt_t[thisband]
            vsini_band   = vsini_t[thisband]
            sigmarv_band = sigmarv_t[thisband]
            # Interpolate along spt and vsini
            lintsigmarv = interp2d(spt_band, vsini_band, sigmarv_band)
            sigmarvs[i,j] = lintsigmarv(spt, vsini)
            # Get resolution coefficients for weighted mean
            if nres > 1:
                rescoeffs[j] = abs(1 - abs(res - ress[j]) / np.diff(ress))

    # Scale RVs because the RV accuracy in Figueira+2016 was set at SNR=100
    sigmarvs *= np.sqrt(1e2/SNR)

    # Combine results of different resolutions with weighted mean
    if nres > 1:
        sigmarvs2 = np.zeros(nbands)
        for i in range(nbands):
            sigmarvs2[i] = sigmarvs[i][0]*rescoeffs[0] + \
                           sigmarvs[i][1]*rescoeffs[1]
    else:
        sigmarvs2 = sigmarvs

    # Combine all bands and floor at a 1 m/s stability
    sigmaRV = 1./np.sqrt(np.sum(1./sigmarvs2**2))
    return sigmaRV


def mr_WM14(rp):
    '''Convert the input planetary radius in Rearth to mass in Mearth.'''
    try:
        rp = float(rp)
    except TypeError:
        raise ValueError('Input radius must be a floating point number.')
    # Compute mass
    return .44*rp**3 + .614*rp**4 if rp < 1.5 else 2.69*rp**(.93)


def print_conclusions_1circularplanet(effsigmaRV, sigmaRV, K, detsigs, Nvisits, texp):
    '''Print the resulting number of visits.'''
    print '\nExpected RV semi-amplitude = %.3f m/s'%K
    print 'With %.3f minute integration times:'%texp
    print 'RV measurement uncertainty = %.3f m/s'%sigmaRV
    print 'Effective RV uncertainty = %.3f m/s\n'%effsigmaRV
    for j in range(len(detsigs)):
        print '%i visits are required to obtain a %i sigma mass detection (total on-sky time = %.3f hrs)' %(Nvisits[j],
                                                                                                            detsigs[j],
                                                                                                            (Nvisits[j]*texp)/60)
        if Nvisits[j] < 10:
            print '\n***WARNING*** %i observations may be too few to adequately sample the full planetary orbit.\n'%Nvisits[j]


def Nvisits_to_get_sigma_1circularplanet(detsig, sigRV, Ms, sigMs, P, K):
    '''Estimate the number of visits required to measure the planet mass to 
    a specified detection significance (e.g. detsig = 3 sigma).
    sigRV = RV measurement uncertainty in m/s
    Ms = stellar mass in MSun
    sigMs = stellar mass uncertainty in MSun
    P = planetary orbital period from transits in days
    K = expected planetary RV semiamplitude in m/s'''
    # Search N until we get to the desired detection significance
    atdetsig, moreN = False, 0
    while not atdetsig:
        # Compute sigmaK as a function of N
        Ntest = np.arange(1, 101+moreN)
        sigK = compute_sigK_1circularplanet(sigRV, Ntest)

        # Compute resulting sigmamp as a function of N 
        uKs  = unp.uarray(K, sigK)
        uMs = unp.uarray(Ms, sigMs)
        umps = rvs.RV_Mp(P, uMs, uKs)
        detsigs = unp.nominal_values(umps) / unp.std_devs(umps)

        # Did we get the desired detsig?
        if detsigs.min() <= detsig <= detsigs.max():
            Nvisits = Ntest[abs(detsigs-detsig) == abs(detsigs-detsig).min()][0]
            atdetsig = True
        elif detsig < detsigs.min():
            raise ValueError('It is too easy to get measure the planet at this detection significance.')
        else:        
            moreN += 25
    return Nvisits


def compute_sigK_1circularplanet(sigRV, N):
    '''Evaulate the model for sigmaK given the number of observations 
    and the effective RV uncertainty for a single planet system.'''
    return sigRV * np.sqrt(2./N)


def _compute_sigK_multiplanet(Ps, T0s, bjd, sigRV):
    '''Compute sigmaK from the Fisher information matrix evaluated at the 
    RV epochs of observation bjd. Still assumes that all planets are close to 
    circular.

    Can solve the problem for systems with either 
    2 or 3 known transiting planets. If only one planet is known, use the function
    compute_sigK_1circularplanet instead of this function.'''
    # Compute the Fisher information matrix with partials wrt the RV semiamplitudes of each planet
    nplanets = len(Ps)
    if nplanets <= 1:
        raise ValueError('This is not a multi-planet system.')
    B = np.zeros((nplanets, nplanets))
    for i in range(nplanets):
        for j in range(nplanets):
            B[i,j] = sin_phases(Ps[np.array([i,j])], T0s[np.array([i,j])], bjd) / sigRV**2
    C = abs(np.linalg.inv(B))
    return np.sqrt(np.diag(C))


def _sin_phases(Ps, T0s, bjd):
    '''Compute the factor for the Fisher information matrix sum_k^N sin(phase_i,k)*sin(phase_j,k).'''
    phase0 = foldAt(bjd, Ps[0], T0=T0s[0])
    phase1 = foldAt(bjd, Ps[1], T0=T0s[1])
    return  np.sum(np.sin(2*np.pi*phase0) * np.sin(2*np.pi*phase1))


if __name__ == '__main__':
    spirou = 1
    if spirou:
        mags, bands, res = (10.41, 9.41, 8.67, 8.32), ['Y','J','H','K'], 75e3
    else:
        mags, bands, res = (10.41, 9.41, 8.67), ['Y','J','H'], 1e5
    SpT, vsini, Ms, P, rp = 4.5, .1, unp.uarray(.181,.019), 1.62893, 1.14
    texp = 0.     # zero means it will solve for texp at a SNR of 150
    jitter = 2.4   # add this to sigmaRV in quadrature
    
    # GJ 273
    mags=(6.5,5.714,5.219,4.857)
    bands=('Y','J','H','K')
    SpT, vsini, Ms, P, rp, jitter = 3.5, 2*np.pi*rvs.Rsun2m(.293)*1e-3/rvs.days2sec(99), .29, 18.65, 1.32, 4.2
    compute_Nvisits_1circularplanet((mags,SpT,vsini,Ms), (P,rp), mr_WM14, [5,10], texp=texp,
                                    addjitter2sigmaRV=jitter, res=res, bands=bands)

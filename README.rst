Transiting_RVcalculator
=====================================================================

A tool used to estimate the number of radial velocity measurements required to detect the mass of a known transiting planet at a desired detection significance based on the Fisher information.

Given the stellar apparent magnitudes in the various spectral bins of the velocimeter (e.g. Y,J,H,K), the numerical spectral type, projected rotation velocity, stellar mass, and its uncertainty:

.. code:: python

   mags = [10., 8.9, 8.3, 7.9]
   SpT = 4.5
   vsini = .1  # km/s
   Ms, sigMs = .2, .02  # solar masses

plus the planet's radius, orbital period, and its uncertainty from the transit light curve analysis:

.. code:: python

   P, sigP, rp = 20, 1e-5, 3  # days, days, and Earth radii

one can estimate the number of RV measurements required to obtain a detection of the planet's mass at a specified 
detection significance by assumming a mass-radius relation to estimate the planet's mass from its radius e.g.

.. code:: python

   def simple_MR(rp):
      return float(rp)**2  # planet mass in Earth masses assuming a constant surface gravity

Now we can calculate the number of RV measurements (N) required to detect the planet's mass at 5 sigma using

.. code:: python

    from nRVcalculator import *
    detsig = 5
    self = nRVcalculator((mags,SpT,vsini,Ms,sigMs), (P,sigP,rp), simple_MR, detsigs=detsig)

The estimation of N is accessed via

.. code:: python

    print self.nRVs

and is 62 in this example for a planet with an RV semi-amplitude of 2.26 m/s in an RV time-series with an average measurement uncertainty of 2.38 m/s.

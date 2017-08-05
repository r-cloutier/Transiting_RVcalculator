Transiting_RVcalculator
=====================================================================

A tool used to estimate the number of radial velocity measurements required to detect the mass of a known transiting planet at a desired detection significance based on the Fisher information.

.. code:: python

   import rebound
   sim = rebound.Simulation()
   sim.add(m=1.0)
   sim.add(m=1.0e-3, a=1.0)
   sim.integrate(1000.)
   sim.status()

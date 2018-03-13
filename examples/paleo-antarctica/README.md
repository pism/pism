Antarctica glacial-cycle spin-up example
=========

The scripts in this directory apply the PISM model to the Antarctic
Ice Sheet using climate forcing over the last two glacial cycles.

These scripts shows one selected simulation out of a parameter ensemble (https://github.com/pism/pism-ensemble), 
which has been used as reference simulation in

* J. Kingslake, R. Scherer, T. Albrecht, J. Coenen, R. Powell, R. Reese, N. Stansell, S. Tulaczyk, M. Wearing, P. L. Whitehouse. 
_Extensive Holocene ice sheet grounding line retreat and uplift-driven re-advance in West Antarctica_, **final review**



Input data
---------

The processed input data are provided on request (contact albrecht@pik-potsdam.de). They include

  * Topography: ice thickness, bed topography and ice surface elevation from [Bedmap2](https://www.bas.ac.uk/project/bedmap-2/).
  
  
  * Geothermal heatflux (Shapiro & Ritzwoller 2004) from [Albmap](http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica).
  
  
  * Climate: precipitation and air temperature from [RACMOv2.1](https://www.projects.science.uu.nl/iceclimate/models/racmo.php), forced with HadCM3, 1980-1999.
  
 
 * [Ocean](http://science.sciencemag.org/content/346/6214/1227): salinity and temperatures observations of Antarctic Continental Shelf Bottom Water (contact S. Schmidtko)
  
  
  
  * Temperature forcing: surface temperature reconstruction at [WDC](http://www.pnas.org/content/113/50/14249) (WAIS Divide Ice Core, contact C. Buizert), combined with [EDC](ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt) at EPICA Dome C before 67kyr
  
  
  * Sea-level forcing: eustatic sea-level reconstructions from [ICE-6G](http://www.atmosp.physics.utoronto.ca/~peltier/data.php) (contact D. Peltier), combined with [SPECMAP](https://doi.pangaea.de/10.1594/PANGAEA.734145) before 120kyr
  


Preprocessing of present-day input data for Antarctic simulations with PISM are described in https://github.com/pism/pism-ais 



Running the 15km-grid spinup
---------

First set your environment paths to locate input data, output directory and pismr executable in:

    $ set_environment.sh

Next look at what the run script will do (i.e. a dry-run):

    $ PISM_DO=echo ./run-paleo.sh 8

Then actually do the run; this will take a number of processor-hours:

    $ ./run_paleo.sh 8

This is a 15 km grid (the default) run with 8 processes. The first stage essentially smooths the surface, the second stage improves the enthalpy field (a "no mass continuity" run), and then the third stage optimizes the till friction angle for a best fit to present-day ice thickness. Grounding line position and ice shelves are prescribed in this stage but the 'full' physics are applied.

In a fourth stage the actual paleo-climate forced reference simulation is performed, with applied sea-level forcing, surface temperature and ocean temperature forcing. Grounding line and calving front can freely evolve.


Reproduce reference experiment
---------

To perform the actual reference simulation, see scripts in PISM release `pik_0.7-holocene-gl-rebound`.
Look in `examples/paleo-antarctica`.  These scripts will require modifications
to run under more recent versions of PISM.


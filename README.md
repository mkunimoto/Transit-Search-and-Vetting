# Transit-Search-and-Vetting

Transiting exoplanet search and vetting pipeline described in Kunimoto et al. (2019), written for use on light curve data from the original Kepler mission to find multi-planet systems.

The pipeline also uses: 

- the [Kepler Transit Model Codebase](https://github.com/jasonfrowe/Kepler) for data reduction and the transit search (Rowe, J. 2016, Zenodo. http://doi.org/10.5281/zenodo.60297)
- [Model-Shift](https://github.com/JeffLCoughlin/Model-Shift) for part of the vetting pipeline (Coughlin, J. L. 2017, KSCI-19105-002)
- [isochrones](https://github.com/timothydmorton/isochrones) for fitting of stellar properties (Morton, T. 2015, Astrophysics Source Code Library, ascl:1503.010)
- [vespa](https://github.com/timothydmorton/VESPA) for false positive probability calculation (Morton, T. 2015, Astrophysics Source Code Library, ascl:1503.011)
- [emcee](https://github.com/dfm/emcee) for general MCMC fitting (Foreman-Mackey, D. et al. 2012, PASP, 125, 925)

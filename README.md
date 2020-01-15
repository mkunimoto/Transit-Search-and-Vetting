# Transit-Search-and-Vetting

Transiting exoplanet search and vetting codes described in Kunimoto et al. (2020), for use on light curve data to find multi-planet systems.

The pipeline also uses: 

- the [Kepler Transit Model Codebase](https://github.com/jasonfrowe/Kepler) for data reduction and the transit search (Rowe, J. 2016, Zenodo. http://doi.org/10.5281/zenodo.60297)
- [Model-Shift](https://github.com/JeffLCoughlin/Model-Shift) for part of the vetting pipeline (Coughlin, J. L. 2017, KSCI-19105-002)

If you use the following codes, please cite the following:

Kunimoto, M., Matthews, J., and Ngo, H. 2020, Searching the Entirety of Kepler Data. I. 17 New Planet Candidates Including 1 Habitable Zone World, submitted to AAS journals.

## Codes Provided

**Setup**

`sigma_clip.py`: iteratively sigma-clip outliers in the positive flux direction, using user-provided sigma value. Calculates sigma using the Median Absolute Deviation.

`remove_gaps.py`: removes data near gaps of user-defined width in a light curve.

Note: full pipeline uses `kfitsread` from the Kepler Transit Model Codebase to combine downloaded light curves from MAST into a single, median-subtracted light curve, and `detrend5` to detrend the data with a given detrending width.

**Searching**

`remove_transit.py`: cuts a transit from a light curve by removing all data within 2 transit durations of the centres of each transit.

`first_checkpoint.py`: returns a transit's signal-to-noise ratio (SNR), robust statistic (RS), chi-square statistic (CHI), and the number of transits for use as a first test against false alarms (see Sec. 2.3 of the paper).

Note: full pipeline uses `transitfind2` from the Kepler Transit Model Codebase to search light curves for transit signals.

**Vetting**

`get_full_template.py`: takes a list of periods, epochs, and depths corresponding to an arbitrary number of planets and produces a multi-planet model file for use in the Kepler Transit Model Codebase.

`transit_tests.py`: outputs results for the depth mean-to-median test (Sec. 3.2.3), chases tests (Sec. 3.2.4), and signal event domination test (Sec. 3.2.7).

`model_tests.py`: fits transit model and provides best-fit parameters. Requires stellar properties and limb-darkening parameters as inputs. Also outputs results for the transit model fit test (Sec. 3.2.1), transit model SNR test (Sec. 3.2.2), first odd-even depth test (Sec. 3.3.3), and V-shape test (Sec. 3.3.4).

Note: full pipeline uses Model-Shift for the second uniqueness test (Sec. 3.2.5), transit shape test (Sec. 3.2.6), significant secondary test (Sec. 3.3.1), and second odd-even depth test (Sec. 3.3.3).

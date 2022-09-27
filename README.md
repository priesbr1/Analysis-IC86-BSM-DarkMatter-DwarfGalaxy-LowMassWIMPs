# Analysis-IC86-BSM-DarkMatter-DwarfGalaxy-LowMassWIMPs

Repository containing code for analysis searching for neutrino signals from low-mass WIMPs in dwarf galaxies.  
[Analysis Wiki](https://wiki.icecube.wisc.edu/index.php/Dark_matter_decay_in_dwarf_galaxies_with_DRAGON)

## Requirements
Software requirements:
* `I3Tray` (if data processing from `i3` files is required)
* `skylab` (see [skylab repository](https://github.com/icecube/skylab) for additional instructions)

Data:
* DRAGON dataset
* PPPC spectra
* dwarf galaxy data (RA/dec, angular sizes, J-factors)

## Files
Analysis scripts are located in the `scripts` directory:
* `background_pdf_trials.py` - calculates trials from background and signal PDFs under background hypothesis
* `DM_FluxComputation.py` - contains functions for neutrino oscillations
* `DM_neutrino_oscillations.py` - applies neutrino oscillations to (and plots) PPPC annihilation spectra
* `energy_range_optimization.py` - calculates optimal energy ranges based on PPPC annihilaiton spectra
* `event_cuts_l5.py` - includes cuts applied to data sample for energy range optimization
* `generate_background_pdfs.py` - generates background PDFs from DRAGON exp data
* `generate_signal_pdfs.py` - generates background PDFs from DRAGON MC data and source locations
* `i3_to_npy.py` - converts DRAGON data from `i3` format to `npy` format (requires `I3Tray`)
* `sensitivity_cross_section.py` - calculates sensitivity flux and cross section from trials results (requires `skylab`)
* `signal_pdf_trials.py` - calculates trials from background and signal PDFs under signal hypothesis

Plotting scripts are located in the `plotting_scripts` directory:
* `plot_background_trials.py` - plots TS distributions from background trials
* `plot_cross_sections.py` - plots expected upper limits of annihilation cross sections
* `plot_energy_ranges.py` - plots signal-tobackground ratio in energy range phase space
* `plot_PDFs.py` - plots signal and background PDFs
* `plot_signal_trials.py` - plots TS averages and signal recovery from signal trials
* `plot_source_info.py` - plots information on dwarf galaxy sources (skymap of locations, J-factors, angular sizes)
* `plot_spectra.py` - plots PPPC annihilation spectra
* `plot_variables.py` - plots distributions of variables in DRAGON data

## Recommended Procedure
1. `i3_to_npy.py` (if DRAGON data is not already in `npy` format)
2. `DM_neutrino_oscillations.py`
3. `energy_range_optimization.py`
3. `generate_background_pdfs.py`
4. `generate_signal_pdfs.py`
5. `background_pdf_trials.py`
6. `signal_pdf_trials.py`
7. `sensitivity_cross_section.py`

## Notes
* All scripts support argument parsing through python's `argparse` package for increased modularity.
* Longer-running scripts (e.g., `background_pdf_trials.py`, `energy_range_optimization.py`, `generate_background_pdfs.py`, `generate_signal_pdfs.py`, `i3_to_npy.py`, `signal_pdf_trials.py`) are recommended to be submitted as jobs. Other scripts finish on timescales such that running them in the terminal is sensible if desired.
* All scripts that generate intermediate datasets (`background_pdf_trials.py`, `generate_background_pdfs.py`, `generate_signal_pdfs.py`, `signal_pdf_trials.py`) and plot randomized/scrambled data (`plot_variables.py`) support RNG seeding for reproducibility.
  * Background trials take ~200,000-300,000 compute hours per channel/mass combination. Because of this, background trials have been *heavily* parallelized. The current parallelization technique is: 1 job array of 500 jobs, 2 cores per job, 10 trials per core. At full parallelization, this takes ~4-6 hours per channel/mass combination. However, due to heavy parallelization, RNG seeding is *not* recommended.
  * Though the computation time for signal trials is not nearly as long, signal trials are submitted as a job array of 300 jobs. At full parallelization, signal trials for a channel/mass combination can be calculated ~1-2 minutes.

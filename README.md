# Analysis-IC86-BSM-DarkMatter-DwarfGalaxy-LowMassWIMPs

Repository containing code for analysis searching for neutrino signals from low-mass WIMPs in dwarf galaxies.  
[Analysis Wiki](https://wiki.icecube.wisc.edu/index.php/Dark_matter_decay_in_dwarf_galaxies_with_DRAGON)

## Requirements
Software requirements:
* `I3Tray` (if data processing from `i3` files is required)
* `skylab` (see [skylab repository](https://github.com/icecube/skylab) for additional instructions)

Data:
* DRAGON dataset
* [PPPC](http://www.marcocirelli.net/PPPC4DMID.html) spectra
* dwarf galaxy data (RA/dec, angular sizes, J-factors)
* MC data for energy optimization
* comparison cross section results from other analyses

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

### MSU Cluster
1. `python scripts/i3_to_npy.py --files "Level5p_IC86.2011_data.??????.i3.bz2" --inpath "/gfps/home/binfalse-002/neergarr/icecube/data_symlink" --outpath "<name of output folder for npy DRAGON data>" --MC 0`  
  a. This only needs to be run if the DRAGON data is not already in `npy` format.  
  b. This should be run separately for each year (2011-2017), changing the `--files` argument for each.  
2. `python scripts/i3_to_npy.py --files "Level5p_IC862013_genie_numu.014674.??????.i3.bz2" --inpath "/mnt/research/IceCube/jpandre/Matt/Level5p/numu/14674/" --outpath "<name of output folder for npy DRAGON data>" --MC 1`  
  a. This only needs to be run if the DRAGON data is not already in `npy` format.  
  b. This only processes MC data from 2013. This is taken to be the same MC data for all years (2011-2017).  
3. `python scripts/DM_neutrino oscillations.py --oscillated 0 --output "<name of output folder for plots>" --flux_units 0`  
4. `python scripts/energy_range_optimization.py --datafolder "/mnt/research/IceCube/jpandre/Matt/level5p/" --spectra "<name of npy file where spectra are stored from (3)>" --channel "b" --mass 10 --part 0 --outfolder "<name of output folder for text files>"`  
  a. This should be run as a job array of 10 jobs (10 hours, 20GB).  
  b. The `part` argument should be passed as the job array ID (one of 0-9).  
5. `python plotting_scripts/plot_energy_ranges.py --datafolder "<name of output folder from (4)>" --channel "b" -mass 10 --outfolder "<name of output folder for plots">`  
6. `python scripts/generate_background_pdfs.py --data "/mnt/research/IceCube/datasets/ps_DRAGON/version-001-p00/IC86_201?_exp.npy" --sources "/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt" --channel "b" --mass 10 --band_size 360 --seed 0 --output "<name of output npy file to store background PDFs>"`  
  a. This should be run as a job (12 hours, 20GB).  
  b. `--band_size 360` will generate an all-sky PDF.  
7. `python scripts/generate_signal_pdfs.py --data "/mnt/research/IceCube/datasets/ps_DRAGON/version-001-p00/IC86_201?_MC.npy" --spectra "<name of npy file where spectra are stored from (3)>" --channel "b" --output "<name of output npy file to store signal PDFs>"`  
  a. This should be run as a job (12 hours, 20GB).  
8. `python scripts/background_pdf_trials.py --background_pdfs "<name of npy file where background PDFs are stored from (6)>" --signal_pdfs "<name of npy file where signal PDFs are stored from (7)>" --sources "/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt" --channel "b" --mass 10 --num_trials 10 --seed 0 --outfolder "<name of output folder to store background trials results>"`  
  a. This should be run as a job array of 500 jobs (12 hours, 1 node, 2 cores, 5GB).  
  b. Instead of passing a single seed, you could pass the job array ID as the seed (similar to `energy_range_optimization.py`).  
9. `python scripts/signal_pdf_trials.py --background_pdfs "<name of npy file where background PDFs are stored from (6)>" --signal_pdfs "<name of npy file where signal PDFs are stored from (7)>" --sources "/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt" --channel "b" --mass 10 --seed 0 --outfolder "<name of output folder to store signal trials results>"`  
  a. This should be run as a job array of 300 jobs (10 minutes, 20 GB).  
  b. Instead of passing a single seed, you could pass the job array ID as the seed (as in `background_pdf_trials.py`).  
10. `python scripts/sensitivitiy_cross_section.py --spectra "<name of npy file where spectra are stored from (3)>" --mass 10 --channel "b" --sources "/mnt/home/priesbr1/DM_Search/analysis_sources_ra_dec_jfactors.txt" --background_trials "<name of output folder with background trials from (8)>" --signal_trials "<name of output folder with signal trials from (9)>" --outfile "<name of output npy file to store cross sections>"`  
  a. This ***must*** be run from the main `skylab` directory, wherever this has been cloned. Please make sure to update any paths above accordingly.  
10. `python plotting_scripts/plot_cross_sections.py --file "<name of output npy file with cross sections from (10)>" --comparison "/mnt/home/priesbr1/DM_Search/data/comparison_data/" --output "<name of output folder for plots>"`  

### Madison Cluster
1. `python scripts/DM_neutrino oscillations.py --oscillated 0 --output "<name of output folder for plots>" --flux_units 0`  
2. `python scripts/energy_range_optimization.py --datafolder "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/energy_optimization_MC/" --spectra "<name of npy file where spectra are stored from (1)>" --channel "b" --mass 10 --part 0 --outfolder "<name of output folder for text files>"`  
  a. This should be run as a job array of 10 jobs (10 hours, 20GB).  
  b. The `part` argument should be passed as the job array ID (one of 0-9).  
3. `python plotting_scripts/plot_energy_ranges.py --datafolder "<name of output folder from (1)>" --channel "b" -mass 10 --outfolder "<name of output folder for plots">`  
4. `python scripts/generate_background_pdfs.py --data "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/ps_DRAGON/version-001-p00/IC86_201?_exp.npy" --sources "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/analysis_sources_ra_dec_jfactors.txt" --channel "b" --mass 10 --band_size 360 --seed 0 --output "<name of output npy file to store background PDFs>"`  
  a. This should be run as a job (12 hours, 20GB).  
  b. `--band_size 360` will generate an all-sky PDF.  
5. `python scripts/generate_signal_pdfs.py --data "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/ps_DRAGON/version-001-p00/IC86_201?_MC.npy" --spectra "<name of npy file where spectra are stored from (1)>" --channel "b" --output "<name of output npy file to store signal PDFs>"`  
  a. This should be run as a job (12 hours, 20GB).  
6. `python scripts/background_pdf_trials.py --background_pdfs "<name of npy file where background PDFs are stored from (4)>" --signal_pdfs "<name of npy file where signal PDFs are stored from (5)>" --sources "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/analysis_sources_ra_dec_jfactors.txt" --channel "b" --mass 10 --num_trials 10 --seed 0 --outfolder "<name of output folder to store background trials results>"`  
  a. This should be run as a job array of 500 jobs (12 hours, 1 node, 2 cores, 5GB).  
  b. Instead of passing a single seed, you could pass the job array ID as the seed (similar to `energy_range_optimization.py`).  
7. `python scripts/signal_pdf_trials.py --background_pdfs "<name of npy file where background PDFs are stored from (4)>" --signal_pdfs "<name of npy file where signal PDFs are stored from (5)>" --sources "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/analysis_sources_ra_dec_jfactors.txt" --channel "b" --mass 10 --seed 0 --outfolder "<name of output folder to store signal trials results>"`  
  a. This should be run as a job array of 300 jobs (10 minutes, 20 GB).  
  b. Instead of passing a single seed, you could pass the job array ID as the seed (as in `background_pdf_trials.py`).  
8. `python scripts/sensitivitiy_cross_section.py --spectra "<name of npy file where spectra are stored from (1)>" --mass 10 --channel "b" --sources "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/analysis_sources_ra_dec_jfactors.txt" --background_trials "<name of output folder with background trials from (6)>" --signal_trials "<name of output folder with signal trials from (7)>" --outfile "<name of output npy file to store cross sections>"`  
  a. This ***must*** be run from the main `skylab` directory, wherever this has been cloned. Please make sure to update any paths above accordingly.  
9. `python plotting_scripts/plot_cross_sections.py --file "<name of output npy file with cross sections from (8)>" --comparison "/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/comparison_data/" --output "<name of output folder for plots>"`  

## Notes
* All scripts support argument parsing through python's `argparse` package for increased modularity.
* Longer-running scripts (e.g., `background_pdf_trials.py`, `energy_range_optimization.py`, `generate_background_pdfs.py`, `generate_signal_pdfs.py`, `i3_to_npy.py`, `signal_pdf_trials.py`) are recommended to be submitted as jobs. Other scripts finish on timescales such that running them in the terminal is sensible if desired.
* All scripts that generate intermediate datasets (`background_pdf_trials.py`, `generate_background_pdfs.py`, `generate_signal_pdfs.py`, `signal_pdf_trials.py`) and plot randomized/scrambled data (`plot_variables.py`) support RNG seeding for reproducibility.
  * Background trials take ~200,000-300,000 compute hours per channel/mass combination. Because of this, background trials have been *heavily* parallelized. The current parallelization technique is: 1 job array of 500 jobs, 2 cores per job, 10 trials per core. At full parallelization, this takes ~4-6 hours per channel/mass combination. However, due to heavy parallelization, RNG seeding is *not* recommended.
  * Though the computation time for signal trials is not nearly as long, signal trials are submitted as a job array of 300 jobs. At full parallelization, signal trials for a channel/mass combination can be calculated ~1-2 minutes.

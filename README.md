# GW200105 Targeted Eccentric Search: Data and Figure Scripts

This repository contains plotting scripts and description of data products associated with the targeted eccentric search for GW200105. The scripts were used to produce figures for the paper below.

**Paper:** K. S. Phukon, S. Schmidt, J. Morras, and A. Pratten, [“Detection of GW200105 with a targeted eccentric search,” *Physical Review D* 113 (2026)](https://doi.org/10.1103/cmtb-5q4b).

**Complete search analysis data:** [Zenodo record 20119453](https://zenodo.org/records/20119453).

## Repository contents

This GitHub repository contains the figure-generation scripts and supporting analysis code. **Data products are not stored in this repository.** Download the HDF5 files from the [Zenodo data release](https://zenodo.org/records/20119453).

The archive includes the eccentric and quasi-circular template banks, injection/recovery summaries, coincident search triggers and background data, and statistical maps. The eccentric template bank uses the TaylorF2Ecc waveform approximant. Its datasets `mass1`, `mass2`, `spin1z`, `spin2z`, and `eccentricity` give the component masses, aligned spin components, and eccentricity. The quasi-circular bank's mass and spin datasets have the same meaning.

The repository includes the figure scripts (`fig1.py`–`fig7.py`), supporting analysis scripts, and `physrev.mplstyle`.

## Figure scripts

| Script | Figure content |
| --- | --- |
| `fig1.py` | Eccentric template-bank parameters and GW200105 parameter estimates |
| `fig2.py` | Signal-to-noise ratio (SNR) and reduced chi-squared $(\chi_r^2)$ comparisons for eccentric injections |
| `fig3.py` | SNR and $\chi_r^2$ comparisons for quasi-circular injections |
| `fig4.py` | Sensitive volume-time $\langle VT \rangle$ ratio and injection recovery fractions versus IFAR threshold |
| `fig5.py` | Foreground and background ranking statistics and significance |
| `fig6.py` | Corner plots for offsets between parameters of eccentric injections and those of recovered templates |
| `fig7.py` | L1 trigger properties around GW200105 |

## Data files

The scripts expect the data products to be downloaded from Zenodo and available under the filenames specified in each script. Place the required files in the working directory before running a script. Main input files include:

- Template banks: `Eccentric_template_bank_GW200105_search.hdf` and `Quasicircular_template_bank_GW200105_search.hdf`.
- Injection/recovery summaries: `H1L1V1_INJ_FIND_<INJECTION_TYPE>_INJECTIONS_<SEARCH_TYPE>_SEARCH.hdf`, where each type is `ECCENTRIC` or `QUASICIRCULAR`.
- Coincident search data: `<DETECTOR_COMBINATION>-EXCLUDE_ZEROLAG_FULL_DATA_<N>DET_<SEARCH_TYPE>_SEARCH.hdf`.
- Combined statistical maps: `H1L1V1-COMBINE_STATMAP_<SEARCH_TYPE>_SEARCH.hdf`.
- Single-detector trigger files used by some figures: detector-specific injection-trigger or background-trigger HDF5 files, and `L1_triggers_eccentric_search.hdf` for `fig7.py`.

Here, detector combinations include H1, L1, V1, H1L1, H1V1, L1V1, and H1L1V1. Not every file is required for every figure; see the input filenames in the corresponding script. The Zenodo record is the source for these data products.

## Reproducing figures

Use a Python environment with PyCBC, NumPy, Matplotlib, and the additional packages imported by the selected script (including packages such as `h5py`, `scipy`, `corner`, `mpl-scatter-density`, and `cycler`, as applicable). Some plots also require a LaTeX installation with the AMS math packages.

From the repository directory, run a figure script with:

```bash
python fig1.py
```

The scripts use fixed input filenames and write figure outputs to the current directory. Obtain the required HDF5 inputs from the [Zenodo record](https://zenodo.org/records/20119453) and place them in the working directory before running the relevant scripts. The scripts do not provide an automated data-download workflow.

## Citation

If you use data from this research, please cite the paper and the archived data release:

> K. S. Phukon, S. Schmidt, J. Morras, and A. Pratten, “Detection of GW200105 with a targeted eccentric search,” *Physical Review D* 113 (2026). [https://doi.org/10.1103/cmtb-5q4b](https://doi.org/10.1103/cmtb-5q4b)

> Khun Sang Phukon et al., GW200105 targeted eccentric-search data release, Zenodo (2026). [https://zenodo.org/records/20119453](https://zenodo.org/records/20119453)

## License

The code and repository contents are distributed under the [MIT License](LICENSE).

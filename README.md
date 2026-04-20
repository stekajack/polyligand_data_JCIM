# Polyligand Data for JCIM Submission

This repository contains data and code associated with the manuscript:

**“Rational Polyligand Design for G-Quadruplex Multimer Stabilization”**

## Overview

The purpose of this repository is to provide the data underlying the figures reported in the manuscript, along with the scripts and versioning information used in their generation.
The datasets provided here correspond directly to the quantities reported in the paper and are sufficient to reproduce the numerical results and figures.

---

## Repository Contents

### Data (CSV files)

- `K2_data.csv` — Relative shape anisotropy parameter data
- `Lp_data.csv` — Local persistence length data
- `Rg_data.csv` — Radius of gyration data
- `melting_temp_data.csv` — Melting temperature data
- `stacking_fraction_data.csv` — Fraction of intercalated ligands data

These files contain the processed, machine-readable data from the figures and analyses presented in the manuscript.

---

### Code

- `simulation_script.py` — Script used to perform simulations
- `missing_kernels.py` — Analysis kernels functions external to a project referenced in versioning_and_info.txt
- `example_to_launch_analysis.ipynb` — Notebook example to launch a analysis kernels

---

### Versioning and Dependencies

- `versioning_and_info.txt` — :
  - External code repositories used
  - Specific commits/tags
  - Software versions (e.g., ESPResSo)
  - Example simulation command

This file provides the necessary information to reconstruct the computational environment used in this work.

---

## Reproducibility

All data required to reproduce the analyses and figures in the manuscript are provided in processed form (CSV files).

The repository references:
- The **data underlying the reported results**
- The **analysis procedures applied to these data**
- The **software and versioning context** in which the data were generated

---

## Raw Simulation Data

Raw simulation trajectories/output files are not included due to their large size. 
All derived observables used in the manuscript are provided in the repository.
The simulation protocol, parameters, and computational setup are described in the manuscript and in `versioning_and_info.txt`.

---

## Contact

For questions regarding the data or code, please contact the corresponding author.

# Polyligand Data and Analysis for JCIM

This repository contains the simulation and analysis code associated with the manuscript:

**"Rational Polyligand Design for G-Quadruplex Multimer Stabilization"**

where:

- `G4_polybraco_analysis.ipynb`: analysis and figure-generation notebook.
- `simulation_script.py`: script used to generate the simulations.
- `versioning_and_info.txt`: software versions, source revisions, and an example simulation command.

The simulation data and processed analysis results are archived on Zenodo:

**Zenodo archive: [Zenodo dataset](https://doi.org/10.5281/zenodo.21905850)**

To reproduce the paper plots:

1. Download and extract the Zenodo archive.
2. In `G4_polybraco_analysis.ipynb`, replace `{path_to_data}` with the path to the extracted `G4_polyligand_data` directory.
3. Run the notebook from top to bottom.

The required Pressomancy, pmtools, and ESPResSo revisions are recorded in `versioning_and_info.txt`.

Note that in the supplied `simulation_script.py`, the ligand bending-potential call is commented out; this is the `ligands_nobending` configuration. The `ligands_bending` data were generated with that call enabled.

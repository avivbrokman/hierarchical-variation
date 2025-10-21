# Extinction probability of a lineage under hierarchical variation

This repository contains code used to study how hierarchical variation affects the extinction probability in a metapopulation model modeled by a discrete-time branching process. The implementation is written primarily in Julia, with a Python script for generating plots from experiment output.

Contents
--------
- `src/` — Julia source files:
  - `HierarchicalVariation.jl` — top-level module and entry points.
  - `extinction_probability_multipatch_functions.jl` — core functions to compute extinction probabilities for multipatch systems.
  - `minimize_extinction_probability_multipatch.jl` — routines to search/minimize extinction probability across parameter spaces.
  - `parameter_structs.jl` — parameter structures and helper constructors.
- `experiments/` — experiment drivers (Julia scripts used to run study experiments).
  - `study_experiments.jl` — experiment configuration and orchestration for the study.
- `scripts/` — convenience scripts and plotting utilities:
  - `run_study_experiments.jl` — a wrapper to run experiments (Julia).
  - `generate_plots.py` — Python script to generate figures from experiment outputs.
- `output/` — expected location for experiment results and generated figures (gitignored in typical workflows).

Quick start
-----------
Requirements
- Julia 1.6+ (recommended: the version specified in `Project.toml` / `Manifest.toml`).
- Python 3.8+ with `plotnine` for plotting (only needed to run `generate_plots.py`).

Setup
1. Clone the repository and change to its directory:

   git clone <https://github.com/avivbrokman/hierarchical-variation.git>
   cd hierarchical-variation

2. Activate and instantiate the Julia environment (from the repository root):

   julia --project=. -e 'using Pkg; Pkg.instantiate();'

Running experiments (Julia)
--------------------------
Run the study experiments from the `scripts` folder using Julia's project environment:

   julia --project=scripts run_study_experiments.jl

Generated data and figures will be placed under `output/` by default.

Generating plots (Python)
-------------------------
The repository includes a small plotting helper. From the project root run:

   python3 scripts/generate_plots.py

Ensure you have the required Python packages installed.

# BOIN-ETB Simulation Code

## Overview
This repository provides R code for conducting dose-finding trial simulations
based on the BOIN-ETB (A Bayesian optimal interval design considering efficacy
and toxicity in early phase basket trials) design. 
The code identifies the Optimal Dose (OD) by maximizing a utility function 
that jointly accounts for toxicity and efficacy outcomes. 
Simulations are run across multiple baskets (disease subtypes) and dose
levels, and results are summarized across user-defined scenarios. This
code is intended for use by biostatisticians involved in the design and
evaluation of early-phase oncology basket trials.

## Requirements
- R version: 4.6.0 (2026-04-24 ucrt)
- R packages:
  - `Iso`
  - `UniIsoRegression`
  - `mfp`
  - `openxlsx`
  - `tidyverse`
  - `boinet`
  - `grid`
  - `gridExtra`
  - `readxl`

## Usage
1. Set the working directory path in `FOLDER` within `SIM.BOINETB_example.R`.
2. Prepare the scenario input file (`SIM.BOINETB_scenario.xlsx`) and
   place it in the directory specified by `FOLDER`.
3. Open `SIM.BOINETB_example.R` and configure the simulation parameters
   (see key parameters below).
4. Run `SIM.BOINETB_example.R`. The script automatically loads the required
   function files (`FUNC_CALC.SUM.R`, `FUNC_SIM.BOINETB.R`) via `source()`.
5. Output files are saved to the directory specified by `OUTFILE`.

Key parameters to configure:
- `K` : Number of baskets
- `J` : Number of dose levels
- `nsim` : Number of simulation iterations
- `cohortsize` : Cohort size per basket
- `tgt.t` : Target toxicity probability (by basket)
- `tgt.e` : Target efficacy probability (by basket)
- `dev` : Allowable deviation from target toxicity/efficacy
- `u` : Utility scores for the four outcome categories
  (efficacy without toxicity, neither, both, toxicity without efficacy)

## Input / Output
### Input
- SIM.BOINETB_scenario.xlsx:
  Excel file containing simulation scenarios. Each scenario defines the
  true toxicity probabilities (`pt`), true efficacy probabilities (`pe`),
  target dose indices based on utility (`udose`) and efficacy (`edose`),
  and overly toxic dose indices (`odose`) for each basket and dose level.

### Output
- summary_[FILENAME]_ALL.xlsx: Excel file with three sheets:
  - `Scr` : True toxicity, efficacy, and utility values for each scenario
    and basket
  - `Sum1` : Per-basket simulation summary, including ODC selection
    percentages by dose level and dose category (utility- and
    efficacy-based, 3- and 4-category classifications), mean number of
    patients assigned per dose level, and early stopping rates
  - `Sum2` : Trial-level summary of total sample size statistics
    (mean, SD, min, median, max)
- [FILENAME]_sum.txt: Text file containing detailed simulation
  results per scenario, including escalation/de-escalation thresholds,
  observed toxicity/efficacy rates, ODC selection rates, and patient
  allocation summaries.

## Code Description
- SIM.BOINETB_example.R:
  Main script for running BOIN-ETB simulations. Users configure
  simulation parameters in this file, which then calls functions defined
  in the `_FUNC/` subdirectory to execute simulations across all scenarios
  and compile results into output files.

- FUNC_SIM.BOINETB.R:
  Defines the core simulation function `SIM.BOINETB.ALG1()`, which
  implements the BOIN-ETB ALG1 dose-finding algorithm and runs the
  specified number of simulated trials.

- FUNC_CALC.SUM.R:
  Defines the `CALC.SUM()` function, which processes raw simulation
  outputs to compute summary statistics including ODC selection rates,
  patient allocation by dose category, and total sample size summaries,
  and exports results to Excel and text files.

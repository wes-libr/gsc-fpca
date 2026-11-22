# GSC-FPCA: generalized synthetic control for sparse functional data

This repository holds the code for the paper

> Shao L, Pohl KM, Thompson WK. A generalized synthetic control algorithm for sparse functional data. bioRxiv, 2026. doi:10.64898/2026.02.23.707582

The code is for the revised version of the paper (November 2026). Code for the first version, posted on bioRxiv on 25 February 2026, is in the earlier history of this repository.

**The paper cites Lucy Shao's repository, https://github.com/l66shao/GSC-FPCA, as the reference version of this code.** This repository is a mirror of it, kept by Wesley Thompson. If the two ever differ, use that one.

## The method in brief

GSC-FPCA estimates the average treatment effect on the treated (ATT) when units are observed at irregular, sparse times. Each unit's untreated outcome is modeled as a smooth trajectory: a mean function, covariate effects, and a few functional principal components (FPCs) with unit-specific scores, plus measurement error. The model is fitted with a Gibbs sampler.

- The mean function, covariate effects, FPCs, score variances and error variance are learned from the controls only (a cut posterior).
- The scores of each treated unit are estimated from its pre-treatment observations. Treated scores have their own prior mean, with a flat prior, so treatment may depend on the scores (latent ignorability).
- Each posterior draw gives a counterfactual trajectory for every treated unit. Averaging the observed minus counterfactual outcomes over post-treatment observations, with a draw of measurement error, gives draws of the ATT and its 95% posterior predictive interval.

The sampler also draws the scores with prior mean zero for treated units, which the paper uses as a comparison.

## Contents

| Folder | Contents | Results in the paper |
|---|---|---|
| `method/` | `sampler_tm.cpp`: Gibbs sampler (Rcpp and RcppArmadillo). `mfpca.R`: R function `mfpca_1d()` that runs it. `loo_cut.R`: leave-one-out for the cut posterior, used to choose the number of FPCs. | Sections 2 to 4 |
| `simulation/` | Data generation, the competing methods, the simulation driver, and the scripts that summarize results and make the tables and Figure 3. | Section 3 and Appendix D: Tables 1, 2 and 13 to 17, Figure 3 |
| `application/` | Data preparation, all model fits, diagnostics, figures, and table scripts for the NCANDA example. The data are not included. | Section 4 and Appendix B: Tables 3 to 8, Figures 4 to 10 |
| `appendix_c/` | Numerical checks of the convergence and coverage results, and Figures 1 and 2. `prior_mean_zero/` holds an earlier sampler with prior mean zero for treated scores and a variant that also updates the mean and covariate coefficients with treated pre-treatment data. They are used only here. | Appendix C: Tables 9 to 12, Figures 1 and 2 |

Run each script from its own folder; paths are relative to it. Outputs are written to the same folder. Samplers are compiled into `cpp_cache/` on first use; when starting several workers in the same folder, start one first and let it finish compiling, so the others load the compiled code instead of compiling it at the same time.

## Requirements

R 4.3.3 with a C++ compiler for Rcpp. These package versions were used:

| Package | Version | Package | Version |
|---|---|---|---|
| Rcpp | 1.0.12 | loo | 2.6.0 |
| RcppArmadillo | 0.12.8.1.0 | coda | 0.19-4.1 |
| MASS | 7.3-60.0.1 | dplyr | 1.1.4 |
| Matrix | 1.6-5 | ggplot2 | 3.4.4 |
| data.table | 1.14.10 | gridExtra | 2.3 |
| lme4 | 1.1-35.1 | jsonlite | 1.8.8 |
| nlme | 3.1-164 | pracma | 2.4.4 |
| fect | 2.4.5 | readr | 2.1.5 |

`sandwich` and `lmtest` are needed only for two scripts in `application/scanner_check/`. `fect` is needed only for the interactive fixed effects and matrix completion comparisons in the simulation. Version 2.4.5 is from GitHub (`remotes::install_github("xuyiqing/fect")`); newer versions may give different results.

## Fitting the model to other data

`simulation/sim2_functions.R` shows the whole workflow: `fit_gsc()` prepares the data and calls `mfpca_1d()`, and `gsc_att()` turns the posterior draws into ATT draws and intervals. `application/app2_core.R` does the same for the NCANDA data, with more options. In brief:

- Scale time to [0, 1]. The basis is evaluated on a grid with step 0.001.
- Give `mfpca_1d()` one list element per unit: the outcomes (`Y_sparse`), times (`time_sparse`), covariates whose effects vary over time (`X_sparse`) and covariates with constant effects (`X_c_sparse`), each as a one-element list holding a matrix with one row per observation. Which slot a covariate goes in is a choice about its effect, not about whether its value changes: the paper puts baseline covariates (sex, race, baseline ICV) in `X_sparse` so their effects can vary with age, and site, scanner and study-year indicators in `X_c_sparse` as constant shifts. Covariates that the treatment could affect should not be included. For treated units, pass only the pre-treatment observations.
- Set `b_control` to 1 for controls and 0 for treated units, `Q` to the number of FPCs, `sim = FALSE`, and `PSLOO = 1` to store the pointwise log likelihood for `loo_cut.R`. With the default `sim = TRUE`, element 10 of the result is overwritten by simulation parameters.
- The result is a list of posterior draws. Element 1 holds the scores with prior mean zero for treated units and element 10 the scores with the separate treated mean; element 11 holds the treated mean. Elements 2, 3, 5, 6, 7 and 8 hold the mean-function coefficients, score variances, error variance, FPC loadings, coefficients of the time-varying covariate effects, and coefficients of the constant covariate effects. Element 9 holds the basis, and element 12 the pointwise log likelihood.

Choose the number of FPCs with the leave-one-out comparison and the rule of Section 2.3.

## Reproducing the results

Each simulation task and each model fit sets its own random seed, so results do not depend on the number of workers or the order of runs. With the package versions above, rerunning the scripts reproduced the saved results exactly in every check we made: individual simulation tasks, Table 9, the first replicates behind Table 12, the analysis files of the empirical example, the final application run (Figures 4 to 10), and all tables and quoted numbers made from saved results.

### Simulation (Section 3 and Appendix D)

```sh
cd simulation
Rscript run_sim2.R 1 &      # start as many workers as you have cores;
Rscript run_sim2.R 2 &      # they share the tasks through claims/
wait
Rscript rule_check.R 1      # Section 2.3 rule for k, first 100 datasets of each LOOIC scenario
Rscript summ_sim2.R         # summary2.csv, looic_choice.csv, shape100.csv
Rscript make_tables.R       # rows of Tables 1, 2 and 13 to 17 (tab_*_rows.tex)
Rscript make_numbers.R      # numbers quoted in Section 3 (tokens_sim.json)
Rscript make_figure3.R      # Figure 3
Rscript test_loo_cut.R      # optional: checks the leave-one-out calculation
```

`run_sim2.R` describes the scenarios and methods at the top. The full study has 7,700 tasks and takes about 12 CPU hours. `SIM_MAX=1 Rscript run_sim2.R 1` runs only the first replicate of each cell, and `SIM_CELLS="c(1,10)"` restricts the cells. Results go to `results2/w<worker>.csv`. A task is finished when `claims/<cell>_<replicate>.done` exists; after a crash, with no workers running, delete the claim folders that lack a `.done` file.

### Empirical example (Section 4 and Appendix B)

The example uses data from the National Consortium on Alcohol and Neurodevelopment in Adolescence (NCANDA; NDA collection 4513, release NCANDA_NIAAADA_09Y_IMAGING_V01). The data are available to approved researchers through the NIMH Data Archive and are not included here. The scripts expect two input files in `application/data/`, or at the paths given by environment variables:

- `FullData_Y9.csv` (`NCANDA_FULLDATA`): the merged release data through the year 9 follow-up, one row per participant visit, with the columns used in `extract_ncanda_data.R` and `build_ext.R`.
- `brain-abnormalities.csv` (`NCANDA_ABNORMAL`): participants flagged for brain abnormalities (column `mri_xnat_sid`), who are excluded.

`extract_ncanda_data.R` makes the analysis file `ncanda_dat.csv` (`NCANDA_DAT`) from them. `build_ext.R` then takes the outcome, the average of the left and right superior frontal gray-matter volumes, from `FullData_Y9.csv`, checks that the other variables match, and adds regional outcomes and further covariates. In the version of the extraction script used for the first posting, a column-pairing error averaged the left superior frontal volume with the left frontal pole vertex count. The script now pairs the columns by name; because the analyses take the outcome from `build_ext.R`, the results are the same with either version.

```sh
cd application
Rscript extract_ncanda_data.R   # data/ncanda_dat.csv
Rscript build_ext.R         # ext_all.rds and rois.rds (participant level: keep them private)
Rscript queue_k2.R 1        # k = 2: four chains, placebo, sensitivity analyses, negative controls,
                            #   34 regions, 50 pseudo-treatment analyses; also k = 1 placebo and negative control
Rscript queue_k3.R 1        # the same analyses with k = 3
Rscript queue_loo.R         # fits with k = 1 to 5 and leave-one-out for the cut posterior
Rscript queue_reloo.R       # exact refits for controls with Pareto k above 0.7; negative control with site
Rscript diag_k2.R           # R-hat and effective sample size across the four chains
Rscript diag_k3.R
APP_Q=2 APP_NLOOP=5500 APP_BURN=500 APP_SEED=20261104 APP_OUT=out_final_k2 Rscript application_final.R
Rscript make_tables.R       # rows of Tables 4 and 8, pseudo-treatment summary
Rscript make_numbers.R      # rows of Tables 5 and 6, numbers quoted in Section 4 and Appendix B
```

Several workers can run each queue at once (`Rscript queue_k2.R 2`, and so on); they share the analyses through claim folders, and finished analyses are skipped on restart. Fits are saved in `res2/`. Each fit with k = 2 and 5,500 iterations takes about two minutes.

`application_final.R` makes Figures 4 to 10, Table 3 (`demo.csv`), the mixed-model results behind Table 5 (`lme.csv`), and `att_table.csv`, which holds the primary ATT estimates. Its output files map to the figures as follows: Figure 4, `age_distribution_heavy_binge_drinking.png`; Figure 5, `diagnostic plot.png`; Figure 6, `FPC and mean.png`; Figure 7, `ATT.png`; Figure 8, `individual fit.png`; Figure 9, `covariate plots.png`; Figure 10, `race plot.png`.

### Scanner changes and the negative control (Section 4 and Appendix B, Table 7)

`application/scanner_check/` checks whether the recorded scanner changes over follow-up, or any change shared by the scans at a site in a given study year, explain the intracranial volume (ICV) negative control, and splits the ICV estimate into gray matter, white matter and CSF. Run the scripts from `application/`, after `build_ext.R` and `queue_k2.R`, whose fits `k2_negctl_icv`, `k2_negctl_icv_site`, `k2_s1` and `k2_site` they reuse:

```sh
cd application
Rscript scanner_check/desc.R                # scanner models by site and study year, ICV jumps at each change, scanner mix of treated vs control scans
Rscript scanner_check/fe_check.R            # within-person regression check
Rscript scanner_check/refit.R 0             # GSC-FPCA with scanner indicators and study year; ICV components (two workers: 0 and 1)
Rscript scanner_check/refit.R 1
Rscript scanner_check/placebo.R 0           # held-out last pre-treatment visit, with and without scanner indicators
Rscript scanner_check/placebo.R 1
Rscript scanner_check/refit_siteyear.R 0    # a level for each site and study year, estimated from controls
Rscript scanner_check/refit_siteyear.R 1
Rscript scanner_check/figures.R             # estimates_by_adjustment.csv and icv_components.csv (rows of Table 7), and five diagnostic figures
```

Fits go to `res2/`, aggregate tables to `scanner_check/out/` and figures to `scanner_check/fig/`. `desc.R` also saves `scanner_check/d_icv.rds`, which is participant level: keep it private. `desc.R` and `fe_check.R` need the `sandwich` and `lmtest` packages.

### Appendix C

```sh
cd appendix_c
Rscript oracle_check.R              # Table 9 (the rows with N_tr = 80)
Rscript oracle_tm_check2.R          # Table 10
Rscript nongauss_check.R            # Table 11 and the comparisons quoted with it
Rscript general_oracle.R            # three designs with treated score covariance different from D
Rscript latent_check.R controls     # Table 12
Rscript latent_check.R pooled
Rscript latent_summary.R
Rscript fig_run.R 1 2 1 300 &       # Figures 1 and 2: worker 1 of 2, replicates 1 to 300
Rscript fig_run.R 2 2 1 300 &
wait
Rscript make_figures.R
```

## Implementation notes

- Normal draws use the lower Cholesky factor (`arma::chol(Sigma, "lower")`). Armadillo's default `chol` returns the upper factor, and draws built from it have the wrong covariance when there is more than one FPC.
- The sampler leaves the FPCs unconstrained. Only the covariance function they define is identified, and the ATT does not depend on their rotation or scale (Appendix C).

No participant data are included in this repository.

## License

MIT; see `LICENSE`.

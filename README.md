# Pistacia satellite array simulation

Forward simulation and approximate Bayesian computation (ABC) for the terminal
satellite arrays of *Pistacia vera* 'Kerman' hap2, as reported in Fig. 5e of
Lee et al., "Massive satellite arrays drive 20 million years of chromosome-end
restructuring across the *Pistacia* super-pangenome".

Arrays are simulated with [reCYCLing](https://github.com/monroe-lab/reCYCLing),
a forward-time simulator of repeat-array evolution. Each simulated genome
starts as identical copies of one 178 bp monomer and evolves under local
(tandem) duplication, long-range duplication within a chromosome, duplicative
exchange between chromosomes, hotspot feedback, single-monomer gene conversion,
deletion and point substitution. Real and simulated genomes are reduced to the
same 89 summary statistics (`R/abc_summary_stats.R`), and nine model
parameters are estimated for Kerman hap2 with ABC random forests (abcrf).

## Install

```r
install.packages(c("data.table", "abcrf", "ranger", "Matrix", "ggplot2", "patchwork", "remotes"))
remotes::install_github("monroe-lab/reCYCLing", ref = "v1.9.2")
# Biostrings (Bioconductor) is needed for step 1:
# BiocManager::install("Biostrings")
```

Analyses used R 4.5.1 locally and R 4.4.2 on the UC Davis Farm cluster, abcrf
2.0, ranger 0.18.0 and data.table 1.18.2.1. The reference table was generated
with reCYCLing v1.9.0 and the other steps with v1.9.1 and v1.9.2; the simulator
code is identical across these versions.

## Quick start: reproduce Fig. 5e from the frozen inputs

Run everything from the repository root.

```bash
mkdir -p results
cp inputs/* results/
Rscript R/07_fig5e.R                      # figures/fig_g10_main_compact_hap2.{pdf,png}
Rscript R/08_summary_statistics_table.R   # tables/summary_statistics.csv
```

Steps 1, 6 and 7 also read the Kerman hap2 monomer FASTAs in `data/knobs/`.

## Pipeline

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `R/01_score_kerman_hap2.R` | 89 summary statistics for Kerman hap2 | `results/real_genome_stats_g9.csv` |
| 2 | `R/02_reference_table.R` via `farm/02_reference_table.sh` | 45,000 prior draws, one simulated genome each (Farm, 150 tasks × 300) | `results/reftable_g10_chunk_*.rds` |
| 3 | `R/03_assemble_reference_table.R` | merge the chunks and check design invariants | `results/reftable_g10.rds` |
| 4 | `R/04_fit_abcrf.R` | ABC-RF model choice and posteriors for nine parameters | `results/fit_g10_posteriors.csv`, `results/fit_g10_diagnostics.csv` |
| 5 | `R/05_posterior_densities.R` | full posterior samples for the violins in Fig. 5e | `results/fig_g10_main_compact_hap2_postdensity.csv` |
| 6 | `R/06_display_simulations.R` | 15 simulations at the posterior medians (replicate 4 is shown in Fig. 5e) | `results/sim_g10pick_kermanhap2_rep*.rds` |
| 7 | `R/07_fig5e.R` | Fig. 5e | `figures/fig_g10_main_compact_hap2.{pdf,png}` |
| 8 | `R/08_summary_statistics_table.R` | the 89 statistics with Kerman hap2's values and the reference-table range | `tables/summary_statistics.csv`, `figures/fig_g10_statset.{pdf,png}` |
| – | `R/supp_topology_check.R` | checks that exchange and conversion are distinguishable only when the between-chromosome statistics are fitted | `results/diag_g10_topology.csv` |

Step 2 is the slow one (about 70 minutes per 300-simulation task on 8 cores);
the rest run locally in minutes. `R/calibrate_realtime.R` converts clock steps
to years by anchoring the fitted array age to the *Pistacia* crown age.

## Data

`data/knobs/Kerman_hap2_Chr*_*.fasta.gz`: the rotated, aligned ~178 bp monomers
of the 11 satellite arrays of *P. vera* 'Kerman' hap2, one gzipped FASTA per
chromosome, in array order. These are the observed data the model is fitted to.

## Tables

| File | Contents |
|---|---|
| `tables/model_parameters.csv` | every model quantity: role, whether it is fitted, fixed, derived or calibrated, and its prior or fixed value (transcribed from `R/02_reference_table.R` and the reCYCLing source) |
| `tables/summary_statistics.csv` | the 89 summary statistics: definition, fitted or unused, Kerman hap2's value and where it falls in the reference table (step 8) |

## Inputs

`inputs/` freezes the outputs behind the published numbers, so Fig. 5e
reproduces without rerunning the simulations:

| File | From step | Notes |
|---|---|---|
| `reftable_g10.rds` | 3 | 44,998 simulations, 42,751 used for training |
| `real_genome_stats_g9.csv`, `real_genome_stats_g9_abcset.txt` | 1 | Kerman hap2 |
| `fit_g10_posteriors.csv`, `fit_g10_diagnostics.csv` | 4 | posteriors and out-of-bag errors |
| `fig_g10_main_compact_hap2_postdensity.csv` | 5 | |
| `sim_g10pick_kermanhap2_rep{1,4,9}.rds` | 6 | replicate 4 is the simulated genome shown in Fig. 5e; replicates 4, 9 and 1 give the simulated histograms |

All simulations use fixed seeds, recorded in each script.

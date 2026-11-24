#!/usr/bin/env Rscript
# ============================================================================
# Step 1. Score the observed genome: P. vera 'Kerman' hap2.
#
# Reads the 11 Kerman hap2 monomer FASTAs (data/knobs/Kerman_hap2_Chr*_*.fasta)
# and reduces the genome to the 89-statistic catalog in R/abc_summary_stats.R,
# the same functions that score every simulated genome.
#
# The pair-sampling steps are seeded with 20260817 + 8. The 8 is Kerman hap2's
# position in the genotype-haplotype ordering used when the statistics were
# first computed, so this run reproduces inputs/real_genome_stats_g9.csv exactly.
#
# Out: results/real_genome_stats_g9.csv         (one row, 89 statistics)
#      results/real_genome_stats_g9_abcset.txt  (core fitted statistic names)
# Run from the repository root:  Rscript R/01_score_kerman_hap2.R
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(reCYCLing)})
source({
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  file.path(if (length(.f)) dirname(sub("^--file=", "", .f[1])) else "R",
            "abc_summary_stats.R")
})
stopifnot(file.exists("config.R"))   # run from the repository root
source("config.R")
RES <- file.path(SIM_DATA_DIR, "results")

files <- Sys.glob(file.path(DATA_DIR, "knobs", "Kerman_hap2_Chr*_*.fasta"))
if (length(files) != 11L)
  stop("expected 11 Kerman hap2 monomer FASTAs under ", DATA_DIR, "/knobs/, found ",
       length(files))

dt <- assemble_genome_dt(files)
set.seed(20260817L + 8L)                     # see header
st <- repeat_genome_stats(dt, max_per_group = 2000L, max_total = 500000L)
v  <- repeat_genome_predictors(st)
cat(sprintf("Kerman hap2: chrom=%d nmono=%d div=%.2f xpairs=%.3f overdisp=%.2f young=%.3f\n",
            length(files), nrow(dt), v[["mean_div_bp_m"]], v[["frac_pairs_xchrom"]],
            10^v[["log_windiv_overdisp"]], v[["frac_win_young"]]))

out <- data.table(genotype = "Kerman", hap = "hap2", group = "vera",
                  n_chrom = length(files), n_mono = nrow(dt),
                  as.data.table(as.list(v)))
fwrite(out, file.path(RES, "real_genome_stats_g9.csv"))

ncat <- ncol(out) - 5L
abc_names <- names(repeat_abc_predictors(setNames(rep(0, ncat),
  names(out)[6:ncol(out)])))
writeLines(abc_names, file.path(RES, "real_genome_stats_g9_abcset.txt"))
cat(sprintf("wrote real_genome_stats_g9.csv (%d statistics) and the %d-name ABC set\n",
            ncat, length(abc_names)))

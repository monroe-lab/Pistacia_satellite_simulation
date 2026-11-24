#!/usr/bin/env Rscript
# ============================================================================
# WHY did g10 fail the identifiability gate, and is it fixable?
#
# 04_fit_abcrf.R: exchange inside the conversion-on subset scored OOB 0.333, which
# is exactly chance under the 1/3 spike prior -- the same collapse as g8, even
# though g10's conversion tract was hard-fixed at ONE monomer. So tract length
# was never the cause.
#
# The likely cause is visible in the variable importances: the exchange forest
# and the conversion forest chose the SAME top four statistics --
#   frac_genome_xchrom, log_xchrom_jacc_mean, frac_pairs_xchrom, log_n_xchrom_fam
# Every one of those measures HOW MUCH sequence is shared between chromosomes.
# None measures HOW it is shared. Rare-large-pairwise exchange and frequent-
# single-monomer-many-partner conversion can produce the same total, so with
# this stat set the two mechanisms are genuinely exchangeable and the forest
# has nothing to separate them with.
#
# The statistics that DO encode the topology are the held-out ones (spr, obc,
# xblk, the Jaccard quantiles) -- exactly the ones the g9 validation flagged as
# systematically mispredicted. This script tests, and does not assume, whether
# promoting them to the fitted set rescues exchange:
#
#   arm A  g9/g10 fitted set (56 stats)                     -- the failure
#   arm B  + the 4 topology families (spr, obc, xblk, xjacc) -- the candidate fix
#   arm C  topology stats ONLY, as an upper bound on how much they carry
#
# Reported for each arm: OOB error for "is exchange present?" and "is
# conversion present?", both overall and within the conversion-on subset.
# Below 0.333 = identifiable; at or above = not.
#
# In:  results/reftable_g10.rds
# Out: results/diag_g10_topology.csv  (printed too)
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(ranger)})

# ---- data root ---------------------------------------------------------------
# Code lives in this repo, DATA (results/) lives with the analysis. Point the
# scripts at it with RECYCLING_ABC_DATA; otherwise the current working directory
# is used, which is what the Farm scripts assume (they submit from the data
# root, i.e. the repository root).
DATA_ROOT <- Sys.getenv("RECYCLING_ABC_DATA", unset = getwd())
if (!dir.exists(DATA_ROOT))
  stop("RECYCLING_ABC_DATA points at a directory that does not exist: ", DATA_ROOT)
setwd(DATA_ROOT)

NTREE <- 500L; NCORES <- 6L; SPIKE <- 1/3
set.seed(20260826L)

FAMILY_STATS  <- c("frac_redundant_m", "log_group_m", "log_group_iqr",
                   sprintf("gsz%02d", 1:7), "log_max_grp", "log_pairs_per_mono",
                   "occ_mean_within", "occ_mean_between", "occ_frac_interchrom")
HOTSPOT_STATS <- c("log_windiv_overdisp", "frac_win_young")
FITTED <- c(readLines("results/real_genome_stats_g9_abcset.txt"),
            FAMILY_STATS, HOTSPOT_STATS)
# the held-out topology families (definitions in R/abc_summary_stats.R)
TOPO <- c(sprintf("spr%02d", 1:6), sprintf("obc%02d", 1:6),
          sprintf("xblk%02d", 1:6), sprintf("log_xjacc_q%02d", 1:3))

ref <- readRDS("results/reftable_g10.rds")
ref <- ref[ok == TRUE & oob == FALSE]
have <- function(v) v[paste0("pred_", v) %in% names(ref)]
FITTED <- have(FITTED); TOPO <- have(TOPO)
cc <- paste0("pred_", c(FITTED, TOPO))
ref <- ref[complete.cases(ref[, ..cc])]
cat(sprintf("sims %d | fitted %d | topology %d\n", nrow(ref), length(FITTED),
            length(TOPO)))

ARMS <- list(A_fitted_only = FITTED,
             B_fitted_plus_topology = c(FITTED, TOPO),
             C_topology_only = TOPO)

has_conv <- ref$p_conversion > 0
oob_err <- function(stats, y, rows) {
  X <- as.data.frame(ref[rows, paste0("pred_", stats), with = FALSE])
  rf <- ranger(x = X, y = factor(ifelse(y[rows], "present", "absent")),
               num.trees = NTREE, num.threads = NCORES, seed = 1L)
  rf$prediction.error
}

out <- rbindlist(lapply(names(ARMS), function(a) {
  s <- ARMS[[a]]
  data.table(arm = a, n_stats = length(s),
    xch_overall  = oob_err(s, ref$p_translocation > 0, rep(TRUE, nrow(ref))),
    xch_conv_on  = oob_err(s, ref$p_translocation > 0, has_conv),
    conv_overall = oob_err(s, ref$p_conversion > 0,    rep(TRUE, nrow(ref))))
}))
out[, `:=`(xch_overall = round(xch_overall, 3), xch_conv_on = round(xch_conv_on, 3),
           conv_overall = round(conv_overall, 3))]
cat(sprintf("\nchance under the %.2f spike prior = %.3f; below that = identifiable\n",
            SPIKE, SPIKE))
print(out)
fwrite(out, "results/diag_g10_topology.csv")

verdict <- out[arm == "B_fitted_plus_topology", xch_conv_on]
cat(sprintf("\nVERDICT: adding the topology statistics takes exchange-given-conversion\nfrom %.3f to %.3f (chance %.3f) -- %s\n",
            out[arm == "A_fitted_only", xch_conv_on], verdict, SPIKE,
            if (verdict < SPIKE - 0.02) "RESCUED, the fix is the stat set, not the model."
            else "STILL DEAD; conversion and exchange are not separable by these data."))

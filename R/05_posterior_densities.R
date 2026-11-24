#!/usr/bin/env Rscript
# Extract the FULL ABC-RF posterior (not just quantiles) for Kerman hap2, as a
# weighted sample of the reference-table parameter values, so the main figure
# can show each posterior as a violin/density instead of a CI bar.
#
# regAbcrf posterior = the reference sims reweighted by how often the observed
# genome shares a forest leaf with each sim (predict()$weights). We resample the
# reference parameter values with those weights to get an effective posterior
# sample, rescale to the 0-1 prior axis (same PRIORS as the figure), and write
# a long CSV. Rebuilds the SAME forests as 04_fit_abcrf.R (same NTREE, same subset
# rows per conditional parameter, same log transforms).
#
# Out: results/fig_g10_main_compact_hap2_postdensity.csv  (par, draw [0-1])
suppressPackageStartupMessages({library(data.table); library(abcrf)})
# Data root: the directory containing results/ (default: the repository root).
DATA_ROOT <- Sys.getenv("RECYCLING_ABC_DATA", unset = getwd())
if (!dir.exists(DATA_ROOT))
  stop("RECYCLING_ABC_DATA points at a directory that does not exist: ", DATA_ROOT)
setwd(DATA_ROOT)
set.seed(20260903L)
NTREE <- 500L; NCORES <- 5L; NDRAW <- 4000L

FAMILY_STATS  <- c("frac_redundant_m", "log_group_m", "log_group_iqr",
                   sprintf("gsz%02d", 1:7), "log_max_grp", "log_pairs_per_mono",
                   "occ_mean_within", "occ_mean_between", "occ_frac_interchrom")
HOTSPOT_STATS <- c("log_windiv_overdisp", "frac_win_young")
TOPOLOGY_STATS <- c(sprintf("spr%02d", 1:6), sprintf("obc%02d", 1:6),
                    sprintf("xblk%02d", 1:6), sprintf("log_xjacc_q%02d", 1:3))
abc_names <- c(readLines("results/real_genome_stats_g9_abcset.txt"),
               FAMILY_STATS, HOTSPOT_STATS, TOPOLOGY_STATS)
pred_cols <- paste0("pred_", abc_names)

ref <- readRDS("results/reftable_g10.rds")
ref <- ref[ok == TRUE & oob == FALSE]
ref <- ref[complete.cases(ref[, ..pred_cols])]
real <- fread("results/real_genome_stats_g9.csv")

# same degenerate-stat guard as 04_fit_abcrf.R
splits <- list(ref$p_distal_dup > 0, ref$p_translocation > 0,
               ref$hotspot_alpha > 0, ref$p_conversion > 0)
degenerate <- vapply(pred_cols, function(p) {
  v <- ref[[p]]
  any(vapply(splits, function(g)
    sd(v[g], na.rm = TRUE) == 0 || sd(v[!g], na.rm = TRUE) == 0, TRUE))
}, TRUE)
pred_cols <- pred_cols[!degenerate]; abc_names <- abc_names[!degenerate]

obs_all <- as.data.frame(real[, ..abc_names]); names(obs_all) <- pred_cols
idx <- which(real$genotype == "Kerman" & real$hap == "hap2")
stopifnot(length(idx) == 1L)
obs <- obs_all[idx, , drop = FALSE]
sumsta <- as.data.frame(ref[, ..pred_cols])

# figure priors (must match 07_fig5e.R)
PRIORS <- data.table(
  par = c("age", "p_local", "local_mean", "alpha", "p_distal", "distal_mean",
          "p_xch", "xch_mean", "p_conv"),
  lo    = c(150,  2.5e-5, 1,   0.1, 1e-5, 1,     1e-3, 5,    1e-5),
  hi    = c(1500, 2.5e-3, 300, 5,   0.2,  20000, 0.95, 3000, 2e-3),
  logsc = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
resc <- function(x, lo, hi, logsc)
  if (logsc) (log10(x) - log10(lo)) / (log10(hi) - log10(lo)) else (x - lo) / (hi - lo)

# param -> (reference value vector, subset rows, log-transform-for-forest)
SPEC <- list(
  age         = list(v = ref$n_generations,       rows = rep(TRUE, nrow(ref)), log = FALSE),
  p_local     = list(v = ref$p_local_dup,         rows = rep(TRUE, nrow(ref)), log = TRUE),
  local_mean  = list(v = ref$local_mean,          rows = rep(TRUE, nrow(ref)), log = TRUE),
  p_distal    = list(v = ref$p_distal_dup,        rows = ref$p_distal_dup > 0, log = TRUE),
  distal_mean = list(v = ref$distal_mean,         rows = ref$p_distal_dup > 0, log = TRUE),
  p_xch       = list(v = ref$p_translocation,     rows = ref$p_translocation > 0, log = TRUE),
  xch_mean    = list(v = ref$transloc_chunk_mean, rows = ref$p_translocation > 0, log = TRUE),
  alpha       = list(v = ref$hotspot_alpha,       rows = ref$hotspot_alpha > 0, log = TRUE),
  p_conv      = list(v = ref$p_conversion,        rows = ref$p_conversion > 0, log = TRUE))

out <- list()
for (nm in names(SPEC)) {
  sp <- SPEC[[nm]]
  yt  <- if (sp$log) log10(sp$v[sp$rows]) else sp$v[sp$rows]
  dtr <- data.frame(r = yt, sumsta[sp$rows, ])
  rf  <- regAbcrf(r ~ ., data = dtr, ntree = NTREE, paral = TRUE, ncores = NCORES)
  pr  <- predict(rf, obs, dtr, rf.weights = TRUE, paral = TRUE, ncores = NCORES)
  w   <- as.numeric(pr$weights)                 # weight per training (subset) sim
  stopifnot(length(w) == sum(sp$rows))
  vals <- sp$v[sp$rows]                          # real-unit reference values
  draws <- sample(vals, NDRAW, replace = TRUE, prob = w / sum(w))
  pr_row <- PRIORS[par == nm]
  out[[nm]] <- data.table(par = nm,
    draw = resc(draws, pr_row$lo, pr_row$hi, pr_row$logsc))
  cat(sprintf("[%s] OOB MSE %.4f | eff sample %.0f\n", nm,
              rf$model.rf$prediction.error, sum(w)^2 / sum(w^2)))
}
fwrite(rbindlist(out), "results/fig_g10_main_compact_hap2_postdensity.csv")
cat("done: results/fig_g10_main_compact_hap2_postdensity.csv\n")

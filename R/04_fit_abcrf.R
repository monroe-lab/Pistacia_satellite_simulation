#!/usr/bin/env Rscript
# ============================================================================
# Step 4. ABC random-forest inference for Kerman hap2 (abcrf, 500 trees).
#
# 1. Model choice: one forest per optional mechanism (long-range duplication,
#    exchange, hotspot feedback, gene conversion), present vs absent. Under the
#    1/3 point-mass prior an out-of-bag (OOB) error of 0.333 is chance.
# 2. Identifiability check: within simulations that have conversion switched
#    on, exchange must still be classified clearly better than chance. This
#    holds only when the between-chromosome (topology) statistics are fitted
#    (see R/supp_topology_check.R).
# 3. Posteriors for nine parameters (age, local rate and size, long-range rate
#    and size, exchange rate and size, hotspot exponent, conversion rate). Each
#    mechanism's parameters are trained only on simulations where it is present.
#
# Fitted statistics: the core ABC set (results/real_genome_stats_g9_abcset.txt)
# plus family-size, hotspot and topology statistics; three topology bins that
# are constant within a mechanism class are dropped, leaving 74.
#
# In:  results/reftable_g10.rds, results/real_genome_stats_g9.csv
# Out: results/fit_g10_posteriors.csv, results/fit_g10_diagnostics.csv
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(abcrf)})

# ---- data root ---------------------------------------------------------------
# Code lives in this repo, DATA (results/) lives with the analysis. Point the
# scripts at it with RECYCLING_ABC_DATA; otherwise the current working directory
# is used, which is what the Farm scripts assume (they submit from the data
# root, i.e. the repository root).
DATA_ROOT <- Sys.getenv("RECYCLING_ABC_DATA", unset = getwd())
if (!dir.exists(DATA_ROOT))
  stop("RECYCLING_ABC_DATA points at a directory that does not exist: ", DATA_ROOT)
setwd(DATA_ROOT)

NCORES <- 6L
NTREE  <- 500L
SPIKE  <- 1/3
set.seed(20260826L)

# The error a classifier achieves by ignoring the data entirely and always
# calling the majority class. Under the 1/3 absent / 2/3 present spike prior
# that is 1/3. Any subset OOB at or above this is UNINFORMATIVE.
CHANCE <- SPIKE

FAMILY_STATS  <- c("frac_redundant_m", "log_group_m", "log_group_iqr",
                   sprintf("gsz%02d", 1:7), "log_max_grp", "log_pairs_per_mono",
                   "occ_mean_within", "occ_mean_between", "occ_frac_interchrom")
HOTSPOT_STATS <- c("log_windiv_overdisp", "frac_win_young")   # g9
# g10: promoted from the held-out set. These encode the TOPOLOGY of cross-
# chromosome sharing -- how many chromosomes a family spans (spr), how many
# identical copies sit elsewhere (obc), how long a shared block runs (xblk),
# and the chromosome-pair Jaccard quantiles. Definitions in abc_summary_stats.R.
TOPOLOGY_STATS <- c(sprintf("spr%02d", 1:6), sprintf("obc%02d", 1:6),
                    sprintf("xblk%02d", 1:6), sprintf("log_xjacc_q%02d", 1:3))
abc_names <- c(readLines("results/real_genome_stats_g9_abcset.txt"),
               FAMILY_STATS, HOTSPOT_STATS, TOPOLOGY_STATS)
pred_cols <- paste0("pred_", abc_names)
stopifnot(!anyDuplicated(abc_names))

ref <- readRDS("results/reftable_g10.rds")
ref <- ref[ok == TRUE & oob == FALSE]
ref <- ref[complete.cases(ref[, ..pred_cols])]
cat(sprintf("training sims: %d | fitted stats: %d (g9 set %d + %d promoted topology)\n",
            nrow(ref), length(pred_cols),
            length(pred_cols) - length(TOPOLOGY_STATS), length(TOPOLOGY_STATS)))
cat(sprintf("mechanism presence: long-range %.3f | exchange %.3f | hotspot %.3f | conversion %.3f\n",
            mean(ref$p_distal_dup > 0), mean(ref$p_translocation > 0),
            mean(ref$hotspot_alpha > 0), mean(ref$p_conversion > 0)))

real <- fread("results/real_genome_stats_g9.csv")[genotype == "Kerman" & hap == "hap2"]
stopifnot(nrow(real) == 1L)

# ---- degenerate-statistic guard ---------------------------------------------
# abcrf fits an LDA axis over the fitted statistics, and LDA errors out on any
# predictor that is constant WITHIN a class. Three of the promoted topology
# stats are: the long cross-chromosome block bins are essentially never
# populated. In the 42,751 usable sims, xblk06 is exactly zero in 100% (sd = 0),
# xblk05 in 99.1%, xblk04 in 91.3% -- and each is constant within the
# no-exchange class. They are dropped here.
#
# This costs nothing. The REAL genomes are just as empty in those bins:
# xblk05 and xblk06 are exactly zero in all 38 haplotypes and xblk04 is
# non-zero in 3 of 38 at a median of 0. A statistic that is zero in both the
# simulations and the data cannot discriminate anything -- and note it was
# never a meaningful held-out check either, since 0 vs 0 "validates" trivially.
splits <- list(ref$p_distal_dup > 0, ref$p_translocation > 0,
               ref$hotspot_alpha > 0, ref$p_conversion > 0)
degenerate <- vapply(pred_cols, function(p) {
  v <- ref[[p]]
  any(vapply(splits, function(g)
    sd(v[g], na.rm = TRUE) == 0 || sd(v[!g], na.rm = TRUE) == 0, TRUE))
}, TRUE)
if (any(degenerate)) {
  cat(sprintf("dropping %d degenerate stat(s), constant within a class: %s\n",
              sum(degenerate),
              paste(sub("^pred_", "", pred_cols[degenerate]), collapse = " ")))
  pred_cols <- pred_cols[!degenerate]
  abc_names <- abc_names[!degenerate]
}
cat(sprintf("fitting on %d statistics\n", length(pred_cols)))

obs  <- as.data.frame(real[, ..abc_names]); names(obs) <- pred_cols
stopifnot(!anyNA(obs))

sumsta <- as.data.frame(ref[, ..pred_cols])
diag_rows <- list()

# ---- 1. model choice --------------------------------------------------------
fit_choice <- function(label, present) {
  dtr <- data.frame(mod = factor(ifelse(present, "present", "absent")), sumsta)
  rf  <- abcrf(mod ~ ., data = dtr, ntree = NTREE, paral = TRUE, ncores = NCORES)
  pr  <- predict(rf, obs, dtr, ntree = NTREE, paral = TRUE, ncores = NCORES)
  votes_present <- pr$vote[, "present"] / rowSums(pr$vote)
  imp <- sort(rf$model.rf$variable.importance, decreasing = TRUE)[1:8]
  diag_rows[[label]] <<- data.table(target = label, kind = "model-choice",
    oob = rf$prior.err, top_stats = paste(sub("^pred_", "", names(imp)), collapse = " "))
  cat(sprintf("[%s] OOB error %.3f (chance %.3f) | top stats: %s\n", label,
              rf$prior.err, CHANCE,
              paste(sub("^pred_", "", names(imp))[1:5], collapse = ", ")))
  list(votes = votes_present,
       post_prob = ifelse(pr$allocation == "present", pr$post.prob, 1 - pr$post.prob))
}
mc_distal <- fit_choice("distal_present",  ref$p_distal_dup > 0)
mc_xch    <- fit_choice("xch_present",     ref$p_translocation > 0)
mc_hot    <- fit_choice("hotspot_present", ref$hotspot_alpha > 0)
mc_conv   <- fit_choice("conv_present",    ref$p_conversion > 0)

# ---- 2. IDENTIFIABILITY GATE (the g8 post-mortem, run as a test) ------------
cat("\n=== identifiability gate: does conversion destroy the other mechanisms? ===\n")
cat(sprintf("chance-level OOB error under the %.2f spike prior = %.3f; a subset\n",
            SPIKE, CHANCE))
cat("OOB at or above that means the mechanism is UNIDENTIFIABLE (the g8 failure).\n")
overall <- c(distal = diag_rows$distal_present$oob,
             xch    = diag_rows$xch_present$oob,
             hot    = diag_rows$hotspot_present$oob)
has_conv <- ref$p_conversion > 0
cat(sprintf("sims with conversion on: %d of %d\n", sum(has_conv), nrow(ref)))
gate <- list(distal = ref$p_distal_dup > 0, xch = ref$p_translocation > 0,
             hot = ref$hotspot_alpha > 0)
gate_pass <- TRUE
for (m in names(gate)) {
  dtr_s <- data.frame(mod = factor(ifelse(gate[[m]][has_conv], "present", "absent")),
                      sumsta[has_conv, ])
  rf_s <- abcrf(mod ~ ., data = dtr_s, ntree = NTREE, paral = TRUE, ncores = NCORES)
  e <- rf_s$prior.err
  ok <- e < CHANCE - 0.02          # must be clearly better than always-guessing
  if (m == "xch" && !ok) gate_pass <- FALSE
  cat(sprintf("  [%s | conversion on] OOB %.3f vs %.3f overall vs %.3f chance  %s\n",
              m, e, overall[[m]], CHANCE, if (ok) "OK" else "** UNIDENTIFIABLE **"))
  diag_rows[[paste0(m, "_given_conv")]] <- data.table(
    target = paste0(m, "_given_conv"), kind = "identifiability", oob = e,
    top_stats = sprintf("overall=%.3f chance=%.3f", overall[[m]], CHANCE))
}

# (b) + (c) conversion recovery and confounding, on the conversion-on sims
dtr_c <- data.frame(r = log10(ref$p_conversion[has_conv]), sumsta[has_conv, ])
rf_c  <- regAbcrf(r ~ ., data = dtr_c, ntree = NTREE, paral = TRUE, ncores = NCORES)
oob_c <- rf_c$model.rf$predictions              # ranger OOB predictions
tru_c <- log10(ref$p_conversion[has_conv])
res_c <- oob_c - tru_c
# correlate the residual with the EXCHANGE parameters: if conversion is soaking
# up exchange signal (the g8 trade-off), this is where it shows.
xch_on <- has_conv & ref$p_translocation > 0
r_px  <- cor(res_c[ref$p_translocation[has_conv] > 0],
             log10(ref$p_translocation[xch_on]))
r_mx  <- cor(res_c[ref$p_translocation[has_conv] > 0],
             log10(ref$transloc_chunk_mean[xch_on]))
cat(sprintf("  [conversion recovery] OOB r = %.3f | MAE = %.3f (log10 units over LU[1e-5,2e-3])\n",
            cor(oob_c, tru_c), mean(abs(res_c))))
cat(sprintf("  [conversion/exchange confounding] cor(conv residual, log10 p_xch) = %+.3f | cor(., log10 xch_mean) = %+.3f\n",
            r_px, r_mx))
diag_rows[["conv_recovery"]] <- data.table(target = "conv_recovery",
  kind = "identifiability", oob = rf_c$model.rf$prediction.error,
  top_stats = sprintf("r=%.3f mae=%.3f cor_pxch=%+.3f cor_xchmean=%+.3f",
                      cor(oob_c, tru_c), mean(abs(res_c)), r_px, r_mx))
cat(sprintf("\nGATE: %s\n", if (gate_pass)
  "PASS -- exchange survives with conversion on; g10 is interpretable." else
  "FAIL -- exchange is unidentifiable with conversion on, exactly as in g8. Do NOT interpret the g10 exchange posterior; report g9."))

# ---- 3. posteriors ----------------------------------------------------------
fit_reg <- function(label, y, rows, log10tr) {
  yt  <- if (log10tr) log10(y[rows]) else y[rows]
  dtr <- data.frame(r = yt, sumsta[rows, ])
  rf  <- regAbcrf(r ~ ., data = dtr, ntree = NTREE, paral = TRUE, ncores = NCORES)
  pr  <- predict(rf, obs, dtr, quantiles = c(0.025, 0.5, 0.975),
                 paral = TRUE, ncores = NCORES)
  bt  <- function(v) if (log10tr) 10^v else v
  diag_rows[[label]] <<- data.table(target = label, kind = "regression",
    oob = rf$model.rf$prediction.error,
    top_stats = paste(sub("^pred_", "", names(sort(
      rf$model.rf$variable.importance, decreasing = TRUE)[1:8])), collapse = " "))
  cat(sprintf("[%s] n_train=%d OOB MSE %.4f\n", label, sum(rows),
              rf$model.rf$prediction.error))
  data.table(med = bt(pr$quantiles[, 2]),
             q025 = bt(pr$quantiles[, 1]), q975 = bt(pr$quantiles[, 3]))
}
all_rows  <- rep(TRUE, nrow(ref))
po_age    <- fit_reg("age",         ref$n_generations,       all_rows, FALSE)
po_plocal <- fit_reg("p_local",     ref$p_local_dup,         all_rows, TRUE)
po_local  <- fit_reg("local_mean",  ref$local_mean,          all_rows, TRUE)
po_pdist  <- fit_reg("p_distal",    ref$p_distal_dup,        ref$p_distal_dup > 0, TRUE)
po_mdist  <- fit_reg("distal_mean", ref$distal_mean,         ref$p_distal_dup > 0, TRUE)
po_pxch   <- fit_reg("p_xch",       ref$p_translocation,     ref$p_translocation > 0, TRUE)
po_mxch   <- fit_reg("xch_mean",    ref$transloc_chunk_mean, ref$p_translocation > 0, TRUE)
po_alpha  <- fit_reg("alpha",       ref$hotspot_alpha,       ref$hotspot_alpha > 0, TRUE)
po_conv   <- fit_reg("p_conv",      ref$p_conversion,        has_conv, TRUE)

out <- data.table(real[, .(genotype, hap, group, n_chrom, n_mono)],
  P_distal_votes  = mc_distal$votes, P_distal_adj  = mc_distal$post_prob,
  P_xch_votes     = mc_xch$votes,    P_xch_adj     = mc_xch$post_prob,
  P_hotspot_votes = mc_hot$votes,    P_hotspot_adj = mc_hot$post_prob,
  P_conv_votes    = mc_conv$votes,   P_conv_adj    = mc_conv$post_prob)
add <- function(nm, po) for (cc in names(po)) out[[paste(nm, cc, sep = "_")]] <<- po[[cc]]
add("age", po_age); add("p_local", po_plocal); add("local_mean", po_local)
add("p_distal", po_pdist); add("distal_mean", po_mdist)
add("p_xch", po_pxch); add("xch_mean", po_mxch)
add("alpha", po_alpha); add("p_conv", po_conv)
fwrite(out, "results/fit_g10_posteriors.csv")
fwrite(rbindlist(diag_rows), "results/fit_g10_diagnostics.csv")

cat("\n=== per-genome summary ===\n")
out[, {
  cat(sprintf("%-12s %s P(LR)=%.2f P(X)=%.2f P(HOT)=%.2f P(CONV)=%.2f age=%4.0f pL=%.2g pX=%.3f alpha=%.2f conv=%.2g [%.2g,%.2g]\n",
      genotype, hap, P_distal_votes, P_xch_votes, P_hotspot_votes, P_conv_votes,
      age_med, p_local_med, p_xch_med, alpha_med, p_conv_med,
      p_conv_q025, p_conv_q975))
}, by = seq_len(nrow(out))]
invisible(NULL)

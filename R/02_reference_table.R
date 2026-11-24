#!/usr/bin/env Rscript
# ============================================================================
# Step 2. Reference table: draw parameters from the prior and simulate one
# genome per draw with reCYCLing, then score it with R/abc_summary_stats.R.
#
# Model (schema "g10"): local duplication, long-range duplication, exchange
# between chromosomes, hotspot feedback, single-monomer gene conversion,
# governed deletion and a fixed substitution clock (5e-5 per bp per step).
# Long-range duplication, exchange, hotspot feedback and conversion each have a
# 1/3 prior probability of being absent. The conversion tract is fixed at one
# monomer and the conversion donor is drawn uniformly over the genome, so it
# falls on another chromosome with probability (K - 1) / K.
# Priors are listed in tables/model_parameters.csv.
#
# Genomes with more than 5e7 identical-monomer pairs are flagged (oob = TRUE)
# and not scored.
#
# Usage:  Rscript R/02_reference_table.R <N> <cores> <seed> <outfile> [nbatch]
# On Farm this runs as farm/02_reference_table.sh (150 tasks x 300 = 45,000).
# Checkpoints in sub-batches, so a preempted task resumes where it stopped.
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(reCYCLing); library(parallel)})

# ABC summary stats are PROJECT code since reCYCLing 1.9.0 became a pure
# simulator (2026-08-24 split); locate abc_summary_stats.R next to THIS script
# so both a local run and a Farm run (where the working directory is the DATA
# root, not this repo) work. Sourced before the setwd() below.
source(local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  home <- Sys.getenv("RECYCLING_ABC_HOME", unset = "")
  cand <- c(if (length(f)) file.path(dirname(sub("^--file=", "", f[1])),
                                     "abc_summary_stats.R"),
            "R/abc_summary_stats.R",
            if (nzchar(home)) file.path(home, "R", "abc_summary_stats.R"))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("abc_summary_stats.R not found; looked in: ",
                         paste(cand, collapse = ", "))
  hit[1]
}))
setDTthreads(1L)

# ---- data root ---------------------------------------------------------------
# The outfile argument below is resolved relative to this root. Point it with
# RECYCLING_ABC_DATA; otherwise the current working directory is used, which is
# what the Farm array assumes (it submits from the data root and passes an
# explicit results/... outfile).
DATA_ROOT <- Sys.getenv("RECYCLING_ABC_DATA", unset = getwd())
if (!dir.exists(DATA_ROOT))
  stop("RECYCLING_ABC_DATA points at a directory that does not exist: ", DATA_ROOT)
setwd(DATA_ROOT)

args   <- commandArgs(trailingOnly = TRUE)
N      <- if (length(args) >= 1) as.integer(args[1]) else 300L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 8L
SEED   <- if (length(args) >= 3) as.integer(args[3]) else 1L
OUT    <- if (length(args) >= 4) args[4] else "results/reftable_g10.rds"
NBATCH <- if (length(args) >= 5) as.integer(args[5]) else 10L

# ---- FIXED constants ---------------------------------------------------------
MU_FIX       <- 5e-5      # molecular clock (the time unit; NEVER fitted)
LOCAL_SHAPE  <- 2
DISTAL_SHAPE <- 2
SIZE_BAND    <- 0.10
SPIKE        <- 1/3       # P(mechanism absent): long-range, exchange, hotspot, conversion

# g10: the two fixed conversion settings (see header)
CONV_TRACT     <- list(type = "fixed", value = 1)   # exactly one monomer
# Same-chromosome donors (the 1/K share) are drawn at this separation. Kept
# local on purpose: within-chromosome homogenization is the standard gene-
# conversion picture, and it is only 1 event in K.
CONV_DIST_MEAN <- 20

# ---- stat-computation guards (g7) --------------------------------------------
MAX_PER_GROUP <- 2000L    # per-family pair cap (exact for families <= 63 copies)
MAX_TOTAL     <- 500000L  # uniform per-chromosome pair cap (shape-preserving)
PAIR_OOB      <- 5e7      # exact genome-wide pair count above which sim = OOB

# ---- prior ranges ------------------------------------------------------------
K_LO  <- 5L;   K_HI  <- 12L
TS_LO <- 1000; TS_HI <- 75000
GEN_LO <- 150; GEN_HI <- 1500
P_LOCAL_LO <- 2.5e-5; P_LOCAL_HI <- 2.5e-3       # g7: FREED, log10-uniform
LOCAL_MEAN_LO <- 1;  LOCAL_MEAN_HI <- 300
DISTAL_P_LO <- 1e-5; DISTAL_P_HI <- 0.2
DISTAL_MEAN_LO <- 1; DISTAL_MEAN_HI <- 20000
XCH_P_LO <- 1e-3; XCH_P_HI <- 0.95
XCH_MEAN_LO <- 5; XCH_MEAN_HI <- 3000
ALPHA_LO <- 0.1; ALPHA_HI <- 5                   # g9: hotspot feedback
# --- g10: single-monomer conversion -------------------------------------------
# 1e-5 is ~1% of monomers converted over a 1000-generation array lifetime and is
# indistinguishable from off; 2e-3 is past the point where the fitted stats
# degrade. The prior spans the whole informative range and a little beyond.
CONV_P_LO <- 1e-5; CONV_P_HI <- 2e-3

rlunif <- function(n, lo, hi) 10^runif(n, log10(lo), log10(hi))

draw_prior <- function() {
  K            <- sample.int(K_HI - K_LO + 1L, 1L) + (K_LO - 1L)
  target_sizes <- round(rlunif(K, TS_LO, TS_HI))
  n_gen        <- round(runif(1, GEN_LO, GEN_HI))
  p_local      <- rlunif(1, P_LOCAL_LO, P_LOCAL_HI)
  local_mean   <- rlunif(1, LOCAL_MEAN_LO, LOCAL_MEAN_HI)
  has_distal   <- runif(1) >= SPIKE
  p_distal     <- if (has_distal) rlunif(1, DISTAL_P_LO, DISTAL_P_HI) else 0
  distal_mean  <- if (has_distal) rlunif(1, DISTAL_MEAN_LO, DISTAL_MEAN_HI) else NA_real_
  has_xch      <- runif(1) >= SPIKE
  p_xch        <- if (has_xch) rlunif(1, XCH_P_LO, XCH_P_HI) else 0
  xch_mean     <- if (has_xch) rlunif(1, XCH_MEAN_LO, XCH_MEAN_HI) else NA_real_
  has_alpha    <- runif(1) >= SPIKE
  alpha        <- if (has_alpha) rlunif(1, ALPHA_LO, ALPHA_HI) else 0
  # g10: conversion rate, same 1/3 spike at zero as the other mechanisms
  has_conv     <- runif(1) >= SPIKE
  p_conv       <- if (has_conv) rlunif(1, CONV_P_LO, CONV_P_HI) else 0
  list(K = K, target_sizes = target_sizes, n_generations = n_gen,
       p_local_dup = p_local, local_mean = local_mean,
       p_distal_dup = p_distal, distal_mean = distal_mean,
       p_translocation = p_xch, transloc_chunk_mean = xch_mean,
       hotspot_alpha = alpha, p_conversion = p_conv)
}

# fixed-length predictor template (names + length contract)
.na_stat <- repeat_genome_stats(data.table(chrom = 1L, num = 1, bponly = "A", dir = "+"))
PRED_TEMPLATE <- repeat_genome_predictors(.na_stat) * NA_real_
PRED_COLS     <- paste0("pred_", names(PRED_TEMPLATE))

run_one <- function(i) {
  pr <- draw_prior()
  # donor drawn uniformly over the genome -> lands off-chromosome (K-1)/K of
  # the time. This is a DERIVED quantity, not a fitted parameter.
  p_conv_inter <- (pr$K - 1) / pr$K
  sim <- tryCatch(
    sim_genome(
      K             = pr$K,
      target_sizes  = pr$target_sizes,
      n_generations = pr$n_generations,
      mu_total      = MU_FIX,
      p_local_dup   = pr$p_local_dup,
      local_dist    = list(type = "gamma", shape = LOCAL_SHAPE,
                           scale = pr$local_mean / LOCAL_SHAPE),
      p_distal_dup  = pr$p_distal_dup,
      distal_dist   = list(type = "gamma", shape = DISTAL_SHAPE,
                           scale = (if (is.na(pr$distal_mean)) 500 else pr$distal_mean) / DISTAL_SHAPE),
      p_exchange          = pr$p_translocation,
      exchange_chunk_mean = if (is.na(pr$transloc_chunk_mean)) 500 else pr$transloc_chunk_mean,
      hotspot_alpha = pr$hotspot_alpha,
      # ---- g10 ----
      p_conversion      = pr$p_conversion,
      conv_tract        = CONV_TRACT,
      conv_dist         = list(type = "gamma", shape = 1, scale = CONV_DIST_MEAN),
      p_conv_interchrom = p_conv_inter,
      size_band     = SIZE_BAND,
      verbose       = FALSE),
    error = function(e) NULL)
  if (is.null(sim) || !nrow(sim)) return(NULL)

  # OOB guard: exact pair count from the group-size table BEFORE any scoring
  gs <- sim[, .N, by = .(chrom, bponly)]$N
  n_pairs_exact <- sum(gs * (gs - 1) / 2)
  oob <- n_pairs_exact > PAIR_OOB

  pred <- PRED_TEMPLATE
  stat_ok <- FALSE
  if (!oob) {
    stat <- tryCatch(repeat_genome_stats(sim, max_per_group = MAX_PER_GROUP,
                                         max_total = MAX_TOTAL),
                     error = function(e) NULL)
    p2 <- if (is.null(stat)) NULL else
      tryCatch(repeat_genome_predictors(stat), error = function(e) NULL)
    if (!is.null(p2) && length(p2) == length(PRED_TEMPLATE)) pred <- p2
    stat_ok <- !is.null(stat) && isTRUE(stat$ok)
  }

  realized <- attr(sim, "array_sizes")
  cbind(
    data.table(K = pr$K, n_generations = pr$n_generations,
               p_local_dup = pr$p_local_dup, local_mean = pr$local_mean,
               p_distal_dup = pr$p_distal_dup, distal_mean = pr$distal_mean,
               p_translocation = pr$p_translocation,
               transloc_chunk_mean = pr$transloc_chunk_mean,
               hotspot_alpha = pr$hotspot_alpha,
               p_conversion = pr$p_conversion,
               conv_tract_fixed = 1L, p_conv_interchrom = p_conv_inter,
               mu_total = MU_FIX, seed = SEED,
               target_size_gm  = 10^mean(log10(pr$target_sizes)),
               target_size_min = min(pr$target_sizes),
               target_size_max = max(pr$target_sizes),
               realized_size_m = if (length(realized)) mean(as.numeric(realized)) else NA_real_,
               n_transloc = as.numeric(attr(sim, "n_exchange")),   # attr renamed in reCYCLing 1.9.0; column name kept for reftable compat
               n_pairs_exact = n_pairs_exact, oob = oob, ok = stat_ok),
    as.data.table(as.list(setNames(as.numeric(pred), PRED_COLS))))
}

cat(sprintf(paste0("reftable_g10: N=%d cores=%d seed=%d nbatch=%d | p_local~LU[%g,%g] ",
                   "cap=%d/%d oob>%g | K~U[%d,%d] ts~LU[%g,%g] age~U[%g,%g] ",
                   "spike=%.2f mu=%g | alpha~LU[%g,%g] | p_conv~LU[%g,%g] ",
                   "tract=1 inter=(K-1)/K -> %s\n"),
            N, CORES, SEED, NBATCH, P_LOCAL_LO, P_LOCAL_HI,
            MAX_PER_GROUP, MAX_TOTAL, PAIR_OOB, K_LO, K_HI, TS_LO, TS_HI,
            GEN_LO, GEN_HI, SPIKE, MU_FIX, ALPHA_LO, ALPHA_HI,
            CONV_P_LO, CONV_P_HI, OUT))

done_ok <- function(f) file.exists(f) && tryCatch({ readRDS(f); TRUE }, error = function(e) FALSE)
if (done_ok(OUT)) { cat("OUT exists and loads -- nothing to do.\n"); quit(save = "no") }

outdir <- dirname(OUT); base <- sub("\\.rds$", "", basename(OUT))
parts_dir <- file.path(outdir, "parts")
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)
part_path <- function(b) file.path(parts_dir, sprintf("%s_b%03d.rds", base, b))
save_atomic <- function(obj, path) {
  tmp <- sprintf("%s.tmp.%d", path, Sys.getpid())
  saveRDS(obj, tmp); invisible(file.rename(tmp, path))
}

batch_sizes <- rep(N %/% NBATCH, NBATCH); rem <- N %% NBATCH
if (rem > 0) batch_sizes[seq_len(rem)] <- batch_sizes[seq_len(rem)] + 1L
batch_sizes <- batch_sizes[batch_sizes > 0]

t0 <- Sys.time()
for (b in seq_along(batch_sizes)) {
  pf <- part_path(b)
  if (done_ok(pf)) { cat(sprintf("[seed %d] batch %d/%d: skip (exists)\n",
                                 SEED, b, length(batch_sizes))); next }
  set.seed(SEED * 100000L + b)
  rows <- mclapply(seq_len(batch_sizes[b]), run_one,
                   mc.cores = CORES, mc.preschedule = FALSE)
  bt <- rbindlist(Filter(Negate(is.null), rows), use.names = TRUE, fill = TRUE)
  save_atomic(bt, pf)
  cat(sprintf("[seed %d] batch %d/%d: %d/%d sims (%d oob) -> %s (%.1f min elapsed)\n",
              SEED, b, length(batch_sizes), nrow(bt), batch_sizes[b],
              sum(bt$oob %in% TRUE), basename(pf),
              as.numeric(Sys.time() - t0, units = "mins")))
}

pfs <- vapply(seq_along(batch_sizes), part_path, character(1))
ref <- rbindlist(lapply(pfs, readRDS), use.names = TRUE, fill = TRUE)
save_atomic(ref, OUT)
invisible(file.remove(pfs))
cat(sprintf("Done: %d/%d sims in %.1f min -> %s | ok=%d oob=%d\n",
            nrow(ref), N, as.numeric(Sys.time() - t0, units = "mins"), OUT,
            sum(ref$ok %in% TRUE), sum(ref$oob %in% TRUE)))

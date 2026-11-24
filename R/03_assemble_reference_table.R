#!/usr/bin/env Rscript
# ============================================================================
# Assemble the Farm g10 reference-table chunks into results/reftable_g10.rds.
#
# The Farm array (farm/02_reference_table.sh, 150 tasks x 300 sims) writes one
# results/reftable_g10_chunk_<task>.rds per task. Tasks run on the preemptible
# `low` partition, so at any moment some chunks are finished, some are partial
# (their sub-batch checkpoints live in results/parts/ and are NOT read here),
# and some have not started. This script takes whatever complete chunks exist
# and reports what is missing rather than silently assembling a short table --
# a reference table that is quietly 40% of its intended size would bias every
# posterior downstream without any visible symptom.
#
# Usage:  Rscript 03_assemble_reference_table.R [chunk_dir] [outfile]
# Safe to re-run: rebuilds from scratch each time.
# ============================================================================
suppressPackageStartupMessages(library(data.table))

# ---- data root ---------------------------------------------------------------
# Code lives in this repo, DATA (results/) lives with the analysis. Point the
# scripts at it with RECYCLING_ABC_DATA; otherwise the current working directory
# is used, which is what the Farm scripts assume (they submit from the data
# root, i.e. the repository root).
DATA_ROOT <- Sys.getenv("RECYCLING_ABC_DATA", unset = getwd())
if (!dir.exists(DATA_ROOT))
  stop("RECYCLING_ABC_DATA points at a directory that does not exist: ", DATA_ROOT)
setwd(DATA_ROOT)

args <- commandArgs(trailingOnly = TRUE)
DIR  <- if (length(args) >= 1) args[1] else "results"
OUT  <- if (length(args) >= 2) args[2] else "results/reftable_g10.rds"
N_TASKS_EXPECTED <- 150L
N_PER_TASK       <- 300L

fs <- list.files(DIR, pattern = "^reftable_g10_chunk_[0-9]+\\.rds$",
                 full.names = TRUE)
if (!length(fs)) stop("no reftable_g10_chunk_*.rds found in ", DIR)
got <- as.integer(sub(".*_chunk_([0-9]+)\\.rds$", "\\1", fs))
miss <- setdiff(seq_len(N_TASKS_EXPECTED), got)

load1 <- function(f) tryCatch(as.data.table(readRDS(f)), error = function(e) {
  cat(sprintf("  UNREADABLE (still being written?): %s\n", basename(f))); NULL })
parts <- Filter(Negate(is.null), lapply(fs[order(got)], load1))
ref <- rbindlist(parts, use.names = TRUE, fill = TRUE)

cat(sprintf("chunks: %d of %d expected%s\n", length(parts), N_TASKS_EXPECTED,
            if (length(miss)) sprintf("  MISSING: %s",
              paste(head(miss, 20), collapse = ",")) else ""))
cat(sprintf("sims: %d of %d intended (%.1f%%)\n", nrow(ref),
            N_TASKS_EXPECTED * N_PER_TASK,
            100 * nrow(ref) / (N_TASKS_EXPECTED * N_PER_TASK)))
cat(sprintf("usable (ok & !oob): %d | oob: %d\n",
            sum(ref$ok %in% TRUE & ref$oob %in% FALSE), sum(ref$oob %in% TRUE)))
cat(sprintf("mechanism presence: long-range %.3f | exchange %.3f | hotspot %.3f | conversion %.3f  (prior 2/3 = %.3f)\n",
            mean(ref$p_distal_dup > 0), mean(ref$p_translocation > 0),
            mean(ref$hotspot_alpha > 0), mean(ref$p_conversion > 0), 2/3))

# the two g10 design invariants -- cheap to check, catastrophic if wrong
stopifnot(all(ref$conv_tract_fixed == 1L))
stopifnot(isTRUE(all.equal(ref$p_conv_interchrom, (ref$K - 1) / ref$K)))
cat("design invariants OK: tract fixed at 1 monomer, p_conv_interchrom = (K-1)/K\n")

saveRDS(ref, OUT)
cat(sprintf("-> %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6))

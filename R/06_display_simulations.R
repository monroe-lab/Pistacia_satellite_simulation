#!/usr/bin/env Rscript
# Pick matrix for KERMAN HAP2: real vs 15 simulation replicates, fingerprints
# for all knob-bearing chromosomes. Columns = real + rep1..15, rows = chroms.
# Also the GENERATOR of the display-sim caches sim_g10pick_kermanhap2_rep*.rds
# that the main figure reads (deterministic seeds 202000 + rep).
#
# The parameter point here is Kerman
# hap2's OWN ABC-RF posterior medians (results/fit_g10_posteriors.csv), with
# hap2's array sizes as simulation targets. Reps cached as
# results/sim_g10pick_kermanhap2_rep%d.rds (seeds 202000 + rep).
suppressPackageStartupMessages({
  library(data.table); library(reCYCLing); library(ggplot2); library(parallel)
})
source({
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  file.path(if (length(.f)) dirname(sub("^--file=", "", .f[1])) else "R",
            "abc_summary_stats.R")
})
## Run from the repository root: Rscript R/06_display_simulations.R
stopifnot(file.exists("config.R"))   # run from the repository root
source("config.R")
RES <- file.path(SIM_DATA_DIR, "results")
TS_CAP <- 75000L; MAXPAIR_CHROM <- 8000L; NREP <- 15L

files <- Sys.glob(file.path(DATA_DIR, "knobs", "Kerman_hap2_Chr*_*.fasta*"))
if (!length(files)) stop("no Kerman hap2 monomer FASTAs under ", DATA_DIR,
                         "/knobs/ -- see data/README.md")
real <- assemble_genome_dt(files)
chrom_labs <- sort(unique(sub(".*_(Chr\\d+)_.*", "\\1", basename(files))))
sizes_real <- real[, .N, by = chrom][order(chrom)]
show_idx <- sizes_real$chrom   # all knob-bearing chromosomes
ts  <- pmin(sizes_real$N, TS_CAP)

post <- fread(file.path(RES, "fit_g10_posteriors.csv"))[genotype == "Kerman" & hap == "hap2"]
stopifnot(nrow(post) == 1L)
ex <- list(age = round(post$age_med), p_local = post$p_local_med,
           local_mean = post$local_mean_med, p_distal = post$p_distal_med,
           distal_mean = post$distal_mean_med, p_xch = post$p_xch_med,
           xch_mean = post$xch_mean_med, alpha = post$alpha_med,
           p_conv = post$p_conv_med)
cat(sprintf("Kerman hap2 posterior medians: age %d p_local %.3g local_mean %.3g p_distal %.3g distal_mean %.0f p_xch %.3g xch_mean %.0f alpha %.3g p_conv %.3g\n",
            ex$age, ex$p_local, ex$local_mean, ex$p_distal, ex$distal_mean,
            ex$p_xch, ex$xch_mean, ex$alpha, ex$p_conv))

sim_rep <- function(rep) {
  cache <- file.path(RES, sprintf("sim_g10pick_kermanhap2_rep%d.rds", rep))
  if (file.exists(cache)) return(readRDS(cache))
  set.seed(202000L + rep)
  s <- sim_genome(K = length(ts), target_sizes = as.integer(ts),
    n_generations = ex$age, mu_total = 5e-5, p_local_dup = ex$p_local,
    local_dist  = list(type = "gamma", shape = 2, scale = ex$local_mean / 2),
    p_distal_dup = ex$p_distal,
    distal_dist  = list(type = "gamma", shape = 2, scale = ex$distal_mean / 2),
    p_exchange = ex$p_xch, exchange_chunk_mean = ex$xch_mean,
    hotspot_alpha = ex$alpha,
    p_conversion = ex$p_conv, conv_tract = list(type = "fixed", value = 1),
    conv_dist = list(type = "gamma", shape = 1, scale = 20),
    p_conv_interchrom = (length(ts) - 1) / length(ts), verbose = FALSE)
  saveRDS(s, cache); s
}
cat("simulating/loading", NREP, "replicates...\n")
sims <- mclapply(seq_len(NREP), sim_rep, mc.cores = 5L)

pair_positions <- function(dt) {
  out <- lapply(split(dt, by = "chrom"), function(g) {
    g  <- as.data.table(g)
    ps <- g[, { n <- .N
      if (n < 2L) NULL else {
        i1 <- rep.int(seq_len(n), (n - 1L):0)
        i2 <- unlist(lapply(seq_len(n - 1L), function(j) (j + 1L):n), use.names = FALSE)
        .(n1 = num[i1], n2 = num[i2], inverted = dir[i1] != dir[i2]) }
    }, by = bponly][, bponly := NULL]
    if (!nrow(ps)) return(NULL)
    ps <- ps[n1 != n2]
    ps[, `:=`(lo = pmin(n1, n2), hi = pmax(n1, n2))][, chrom := g$chrom[1]][]
  })
  rbindlist(Filter(Negate(is.null), out))
}
set.seed(1)
fp_of <- function(dt, who) {
  d  <- dt[chrom %in% show_idx]
  pp <- pair_positions(d)
  sz <- d[, .(n = max(num)), by = chrom]
  pp <- merge(pp, sz, by = "chrom")
  pp[, `:=`(x = lo / n, y = hi / n, source = who)][]
}
srcs <- c("REAL\nKerman h2", sprintf("rep%d", seq_len(NREP)))
fp <- rbindlist(c(list(fp_of(real, srcs[1])),
  mclapply(seq_len(NREP), function(r) fp_of(sims[[r]], srcs[r + 1]),
           mc.cores = 5L)))
fp[, chrom_lab := factor(chrom, levels = show_idx,
  labels = sprintf("%s %sk", chrom_labs[show_idx],
                   round(sizes_real$N[show_idx] / 1000)))]
fp[, source := factor(source, levels = srcs)]
fp_plot <- fp[, if (.N > MAXPAIR_CHROM) .SD[sample.int(.N, MAXPAIR_CHROM)] else .SD,
              by = .(source, chrom)]

p <- ggplot(fp_plot, aes(x, y, colour = inverted)) +
  geom_point(size = 0.02, alpha = 0.5, shape = 16, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#c23b22")) +
  facet_grid(chrom_lab ~ source, switch = "y") +
  scale_x_continuous(expand = c(0.02, 0), breaks = NULL) +
  scale_y_continuous(expand = c(0.02, 0), breaks = NULL, position = "right") +
  coord_fixed() + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.25),
        panel.spacing = unit(1.5, "pt"),
        strip.text.x = element_text(size = 6, face = "bold"),
        strip.text.y.left = element_text(size = 6, angle = 90, face = "bold"))
n_chr <- length(show_idx)
ggsave(file.path(SIM_OUT_DIR, "sim_pick_hap2.png"), p,
       width = 0.92 * (NREP + 1) + 0.7, height = 0.92 * n_chr + 0.7,
       dpi = 200, bg = "white", limitsize = FALSE)
ggsave(file.path(SIM_OUT_DIR, "sim_pick_hap2.pdf"), p,
       width = 0.92 * (NREP + 1) + 0.7, height = 0.92 * n_chr + 0.7,
       bg = "white", limitsize = FALSE)
cat("done:", file.path(SIM_OUT_DIR, "sim_pick_hap2.{png,pdf}"), "\n")

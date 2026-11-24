#!/usr/bin/env Rscript
# ============================================================================
# The ABC training statistic set, shown against the reference table.
#   Exemplar genome: Kerman (P. vera) hap2.
#
# PURPOSE. Every other g10 figure shows the FIT. This one shows the INPUT: what
# the 89 summary statistics actually are, what Kerman hap2's values are, and
# where each of those values sits inside the distribution the random forests
# were trained on. Two views of the same thing:
#
#   a) The seven histogram families plotted as the distributions they are --
#      the raw calculation, real genome against the reference-table envelope.
#      These seven blocks account for 62 of the 89 statistics.
#   b) All 89 statistics as the percentile at which Kerman hap2's observed
#      value falls within the 42,751 training simulations, split by whether the
#      statistic was FITTED (74) or HELD OUT (15).
#
# The percentile is the honest common currency here: the statistics have wildly
# different units (proportions, log10 rates, counts), so a raw overlay is
# meaningless across blocks, and a z-score is unstable for the several
# statistics whose reference distribution is near-degenerate.
#
# The stat set is reconstructed by the SAME expressions used in R/04_fit_abcrf.R,
# not copied from a note, so it cannot drift from what was actually fitted.
#
# In:  results/reftable_g10.rds
#      results/real_genome_stats_g9.csv
#      results/real_genome_stats_g9_abcset.txt
# Out: figures/fig_g10_statset.{pdf,png}
#      tables/summary_statistics.csv (all 89 statistics, Kerman hap2 values)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
## Run from the repository root: Rscript R/08_summary_statistics_table.R
stopifnot(file.exists("config.R"))   # run from the repository root
source("config.R")
RES <- file.path(SIM_DATA_DIR, "results")

GENOTYPE <- "Kerman"; HAP <- "hap2"
BASE <- 6

# ---- 1. the stat set, rebuilt exactly as 04_fit_abcrf.R builds it -----------------
FAMILY_STATS  <- c("frac_redundant_m", "log_group_m", "log_group_iqr",
                   sprintf("gsz%02d", 1:7), "log_max_grp", "log_pairs_per_mono",
                   "occ_mean_within", "occ_mean_between", "occ_frac_interchrom")
HOTSPOT_STATS <- c("log_windiv_overdisp", "frac_win_young")
TOPOLOGY_STATS <- c(sprintf("spr%02d", 1:6), sprintf("obc%02d", 1:6),
                    sprintf("xblk%02d", 1:6), sprintf("log_xjacc_q%02d", 1:3))
abc_names <- c(readLines(file.path(RES, "real_genome_stats_g9_abcset.txt")),
               FAMILY_STATS, HOTSPOT_STATS, TOPOLOGY_STATS)
stopifnot(!anyDuplicated(abc_names), length(abc_names) == 77L)

real <- fread(file.path(RES, "real_genome_stats_g9.csv"))
ID_COLS <- c("genotype", "hap", "group", "n_chrom", "n_mono")   # metadata, not statistics
all_names <- setdiff(names(real), ID_COLS)
stopifnot(length(all_names) == 89L)

ref <- readRDS(file.path(RES, "reftable_g10.rds"))
ref <- ref[ok == TRUE & oob == FALSE]
ref <- ref[complete.cases(ref[, paste0("pred_", abc_names), with = FALSE])]
cat(sprintf("training sims: %d\n", nrow(ref)))

# The same degenerate-statistic guard 04_fit_abcrf.R applies: constant within any
# mechanism class => dropped by abcrf's LDA step. Recomputed, not assumed.
splits <- list(ref$p_distal_dup > 0, ref$p_translocation > 0,
               ref$hotspot_alpha > 0, ref$p_conversion > 0)
degenerate <- vapply(paste0("pred_", abc_names), function(p) {
  v <- ref[[p]]
  any(vapply(splits, function(g)
    sd(v[g], na.rm = TRUE) == 0 || sd(v[!g], na.rm = TRUE) == 0, TRUE))
}, TRUE)
dropped  <- abc_names[degenerate]
fitted_names <- abc_names[!degenerate]
heldout_names <- setdiff(all_names, fitted_names)
cat(sprintf("fitted: %d | dropped as degenerate: %s | held out: %d\n",
            length(fitted_names), paste(dropped, collapse = " "),
            length(heldout_names)))
stopifnot(length(fitted_names) == 74L, length(heldout_names) == 15L)

# ---- 2. blocks and human-readable definitions --------------------------------
# Every description below is taken from the definition in R/abc_summary_stats.R
# (bin counts and ranges from the GX_* / DIST_* / DIV_* constants there).
blk <- function(n) {
  fifelse(grepl("^wxd", n),  "within-chromosome pair distance",
  fifelse(grepl("^gdv", n),  "monomer divergence (clock)",
  fifelse(grepl("^gsz", n),  "identical-family size",
  fifelse(grepl("^owc", n),  "occupancy, own chromosome",
  fifelse(grepl("^obc", n),  "occupancy, other chromosomes",
  fifelse(grepl("^spr", n),  "chromosome spread",
  fifelse(grepl("^xblk", n), "cross-chromosome block length",
  fifelse(grepl("windiv|win_young", n), "spatial heterogeneity (hotspot)",
  fifelse(grepl("xchrom|xjacc|family_spread|blk_runlen", n), "cross-chromosome, scalar",
          "other scalar")))))))))
}

FLOOR <- 1e-5   # display floor for the log-y panels, so exact zeros stay visible
HIST <- list(
  # wxd ticks sit on BIN EDGES (i - 0.5 = lower edge of bin i) and are labelled
  # with the exact edge value 10^((i-1)/20 * log10(75000)) -- the previous
  # centre-ish labels ("10^3" at bin 11, "3x10^4" at bin 16) were off by up to
  # 5x against the true bin spans (273.9-480 and 4,532-7,944)
  list(pre = "wxd",  n = 20, lab = "within-chromosome pair distance", log = FALSE,
       xlab = "distance between identical monomers (monomers)",
       ticks = c(0.5, 5.5, 10.5, 15.5, 20.5),
       tlabs = c("1", "17", "274", "4.5k", "75k"), xlim = c(0.5, 20.5)),
  list(pre = "gdv",  n = 11, lab = "monomer divergence (the clock)", log = FALSE,
       xlab = "substitutions from array consensus",
       ticks = c(1, 4, 7, 11), tlabs = c("0", "6", "12", "22")),
  list(pre = "gsz",  n = 7,  lab = "identical-family size", log = TRUE,
       xlab = "monomers per identical family",
       ticks = c(1, 3, 5, 7), tlabs = c("1", "20", "400", "5000")),
  list(pre = "owc",  n = 6,  lab = "occupancy, own chromosome", log = TRUE,
       xlab = "identical copies on own chromosome",
       ticks = c(1, 3, 6), tlabs = c("1", "20", "2000")),
  # obc conditions on >= 1 between-chromosome copy and bins from 1 (see HDEF),
  # so its first tick is "1" like owc, not "0"
  list(pre = "obc",  n = 6,  lab = "occupancy, other chromosomes", log = TRUE,
       xlab = "identical copies elsewhere",
       ticks = c(1, 3, 6), tlabs = c("1", "20", "2000")),
  list(pre = "spr",  n = 6,  lab = "chromosome spread", log = TRUE,
       xlab = "chromosomes carrying the sequence",
       ticks = 1:6, tlabs = c("1", "2", "3", "4", "5", "6+")),
  list(pre = "xblk", n = 6,  lab = "cross-chromosome block length", log = TRUE,
       xlab = "shared block length (monomers)",
       ticks = c(1, 3, 6), tlabs = c("1", "40", "5000"))
)

# ---- 3. Kerman hap2 against the reference table ------------------------------
obs <- real[genotype == GENOTYPE & hap == HAP]
stopifnot(nrow(obs) == 1L)

pctile <- function(nm) {
  v <- ref[[paste0("pred_", nm)]]
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  mean(v < obs[[nm]], na.rm = TRUE) + 0.5 * mean(v == obs[[nm]], na.rm = TRUE)
}
qs <- function(nm, p) {
  v <- ref[[paste0("pred_", nm)]]
  as.numeric(quantile(v[is.finite(v)], p, names = FALSE))
}

# Per-statistic definitions, transcribed from R/abc_summary_stats.R. Histogram
# bins are described by their family + bin edges; scalars individually.
DEFN <- c(
  frac_win_young        = "fraction of 500-monomer windows whose mean divergence is < 7 substitutions (recent near-homogeneous patches)",
  log_windiv_overdisp   = "log10 of the between-window sd of mean divergence divided by the sd expected under exchangeable windows (1-3 = homogeneous turnover, higher = spatially clustered young patches)",
  blk_runlen            = "mean length of runs of consecutive identical monomers within an array, averaged over all runs including length-1 runs",
  family_spread         = "mean number of chromosomes a cross-chromosome-shared sequence occupies (families on >= 2 chromosomes only)",
  frac_ancestral_xchrom = "fraction of retained monomers belonging to the ancestral (genus-wide) sequence family, excluded before all derived cross-chromosome statistics",
  frac_genome_xchrom    = "fraction of retained monomers whose sequence is derived and present on >= 2 chromosomes (absolute cross-chromosome abundance)",
  frac_pairs_xchrom     = "fraction of identical-monomer pairs that are inter-chromosomal",
  frac_xchrom_inv       = "fraction of inter-chromosomal identical pairs in opposite orientation",
  log_blk_xchrom_mean   = "log10 mean length of contiguous runs of cross-chromosome-shared monomers (transferred block size)",
  log_n_xchrom_fam      = "log10 number of distinct derived sequences present on >= 2 chromosomes",
  log_xchrom_jacc_mean  = "log10 mean pairwise Jaccard similarity of distinct-sequence sets between chromosomes",
  log_xjacc_q01         = "log10 50th percentile of the pairwise between-chromosome Jaccard similarities",
  log_xjacc_q02         = "log10 90th percentile of the pairwise between-chromosome Jaccard similarities",
  log_xjacc_q03         = "log10 99th percentile of the pairwise between-chromosome Jaccard similarities",
  frac_redundant_m      = "mean across chromosomes of the fraction of monomers with >= 1 identical copy on the same chromosome",
  frac_redundant_iqr    = "interquartile range across chromosomes of that redundant fraction",
  log_array_size_m      = "log10 mean number of retained monomers per chromosome",
  log_group_m           = "log10 mean across chromosomes of the mean identical-family size (families of size >= 2 only)",
  log_group_iqr         = "log10 interquartile range across chromosomes of the mean identical-family size (families of size >= 2 only)",
  log_max_grp           = "log10 size of the largest genome-wide identical family",
  log_n_chrom           = "log10 number of knob-bearing chromosomes analysed",
  log_pairs_per_mono    = "log10 mean number of identical pairs per retained monomer",
  mean_div_bp_m         = "mean across chromosomes of the mean per-monomer divergence from the array consensus (substitutions, < 22 cut-off)",
  mean_div_bp_iqr       = "interquartile range across chromosomes of that mean divergence",
  occ_frac_interchrom   = "fraction of monomer sequences with at least one identical copy on a different chromosome",
  occ_mean_between      = "log10 mean number of identical copies a monomer has on other chromosomes",
  occ_mean_within       = "log10 mean number of identical copies a monomer has on its own chromosome"
)
HDEF <- c(
  wxd  = "proportion of within-chromosome identical-monomer pairs in log-spaced distance bin %d of 20 (1 to 75,000 monomers)",
  gdv  = "proportion of retained monomers in divergence bin %d of 11 (0 to 22 substitutions from the array consensus)",
  gsz  = "monomer-weighted proportion of identical families in log-spaced size bin %d of 7 (1 to 5,000 monomers)",
  owc  = "monomer-weighted proportion in log-spaced bin %d of 6 of identical-copy count on the monomer's own chromosome (1 to 2,000)",
  obc  = "monomer-weighted proportion in log-spaced bin %d of 6 of identical-copy count on other chromosomes among monomers with at least one identical copy on another chromosome (bins span 1 to 2,000)",
  spr  = "monomer-weighted proportion of sequences carried on %s chromosome(s)",
  xblk = "proportion of cross-chromosome shared blocks in log-spaced length bin %d of 6 (1 to 5,000 monomers)"
)
defn_of <- function(n) {
  pre <- sub("[0-9]+$", "", n); i <- as.integer(sub("^[^0-9]+", "", n))
  if (pre %in% names(HDEF) && !is.na(i)) {
    return(sprintf(HDEF[[pre]], if (pre == "spr") if (i == 6L) "6 or more" else as.character(i) else i))
  }
  if (n %in% names(DEFN)) return(unname(DEFN[[n]]))
  NA_character_
}

cat_dt <- data.table(
  statistic = all_names,
  block     = blk(all_names),
  definition = vapply(all_names, defn_of, ""),
  role      = fifelse(all_names %in% fitted_names, "fitted",
              fifelse(all_names %in% dropped, "held out (degenerate)", "held out")),
  kerman_hap2 = vapply(all_names, function(n) as.numeric(obs[[n]]), 0),
  ref_q05   = vapply(all_names, qs, 0, p = 0.05),
  ref_med   = vapply(all_names, qs, 0, p = 0.50),
  ref_q95   = vapply(all_names, qs, 0, p = 0.95),
  percentile = vapply(all_names, pctile, 0)
)
cat_dt[, in_ref_90 := percentile > 0.05 & percentile < 0.95]

# ---- panel a: the histogram families as distributions ------------------------
mk_hist <- function(h) {
  cl <- if (h$log) function(x) pmax(x, FLOOR) else identity
  nm <- sprintf("%s%02d", h$pre, seq_len(h$n))
  env <- rbindlist(lapply(seq_along(nm), function(i) {
    v <- ref[[paste0("pred_", nm[i])]]; v <- v[is.finite(v)]
    q <- cl(quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE))
    data.table(bin = i, lo = q[1], q1 = q[2], med = q[3], q3 = q[4], hi = q[5])
  }))
  rl <- data.table(bin = seq_along(nm), nm = nm,
                   val = cl(vapply(nm, function(n) as.numeric(obs[[n]]), 0)),
                   role = fifelse(nm %in% fitted_names, "fitted", "held out"))
  # statistic-name labels, hung vertically from the top of the panel so they
  # never collide with the data or the prior-predictive ribbons.
  if (h$log) {
    ytop <- 10; ylim <- c(FLOOR, 10)   # headroom above 1 for the hanging labels
  } else {
    ytop <- max(c(env$hi, rl$val), na.rm = TRUE) * 1.62
    ylim <- c(0, ytop)
  }
  rl[, ylab := ytop]
  p <- ggplot(env, aes(bin)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
    geom_ribbon(aes(ymin = q1, ymax = q3), fill = "grey70") +
    geom_line(aes(y = med), colour = "grey35", linewidth = 0.35) +
    geom_line(data = rl, aes(y = val), colour = "#0B3D91", linewidth = 0.6) +
    geom_point(data = rl, aes(y = val, shape = role), colour = "#0B3D91",
               size = 1.1, fill = "white", stroke = 0.4) +
    geom_text(data = rl, aes(y = ylab, label = nm, colour = role),
              angle = 90, hjust = 1, vjust = 0.4,
              size = if (h$n >= 20) 1.5 else 1.85,
              family = "mono", show.legend = FALSE) +
    scale_shape_manual(values = c(fitted = 16, `held out` = 21), guide = "none") +
    scale_colour_manual(values = c(fitted = "black", `held out` = "#B22222"),
                        guide = "none") +
    scale_x_continuous(breaks = h$ticks, labels = h$tlabs, limits = h$xlim,
                       expand = expansion(mult = 0.02)) +
    labs(title = h$lab, x = h$xlab,
         y = if (h$log) "proportion (log)" else "proportion") +
    theme_classic(base_size = BASE) +
    theme(plot.title = element_text(size = BASE, face = "bold"),
          axis.title = element_text(size = BASE - 0.5),
          legend.position = "none",
          plot.margin = margin(2, 4, 2, 2))
  p <- if (h$log) {
    p + scale_y_log10(limits = ylim, breaks = c(1e-5, 1e-3, 1e-1),
                      labels = c("0", "-3", "-1"))
  } else {
    p + scale_y_continuous(limits = ylim, expand = expansion(mult = c(0, 0.01)))
  }
  p
}
pa <- wrap_plots(lapply(HIST, mk_hist), nrow = 2)

# ---- panel b: all 89 statistics as reference-table percentiles ---------------
ord <- cat_dt[order(block, statistic)]
ord[, idx := .I]
ord[, col := factor(fifelse(role == "fitted", "fitted (74)",
                    fifelse(role == "held out", "held out, informative (12)",
                            "held out, degenerate (3)")),
                    levels = c("fitted (74)", "held out, informative (12)",
                               "held out, degenerate (3)"))]
half <- ceiling(nrow(ord) / 2)
ord[, side := fifelse(idx <= half, "1", "2")]
ord[, y := fifelse(side == "1", half - idx + 1, nrow(ord) - idx + 1)]

pb <- ggplot(ord, aes(percentile, y)) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  geom_vline(xintercept = 0.5, colour = "grey60", linewidth = 0.25) +
  geom_point(aes(colour = col, shape = col), size = 1.15, stroke = 0.45,
             fill = "white") +
  geom_text(aes(x = -0.04, label = statistic), hjust = 1, size = 1.5) +
  scale_shape_manual(values = c(`fitted (74)` = 16,
                                `held out, informative (12)` = 21,
                                `held out, degenerate (3)` = 4), name = NULL) +
  scale_colour_manual(values = c(`fitted (74)` = "#0B3D91",
                                 `held out, informative (12)` = "#C0392B",
                                 `held out, degenerate (3)` = "#C0392B"),
                      name = NULL) +
  scale_x_continuous(limits = c(-0.42, 1.02),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", ".25", ".5", ".75", "1")) +
  facet_wrap(~ side, nrow = 1, scales = "free_y") +
  labs(x = "percentile of Kerman hap2 within the 42,751 training simulations",
       y = NULL) +
  theme_classic(base_size = BASE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_blank(), strip.background = element_blank(),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(10, "pt"))

fig <- (pa / pb) +
  plot_layout(heights = c(1, 1.55)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = BASE + 1, face = "bold"))

ggsave(file.path(SIM_OUT_DIR, "fig_g10_statset.pdf"), fig, width = 7.2, height = 8.4)
ggsave(file.path(SIM_OUT_DIR, "fig_g10_statset.png"), fig, width = 7.2, height = 8.4, dpi = 400)

# ---- 4. supplementary stat-catalog table -------------------------------------
fwrite(cat_dt[order(match(block, c(
  "within-chromosome pair distance", "monomer divergence (clock)",
  "identical-family size", "occupancy, own chromosome",
  "occupancy, other chromosomes", "chromosome spread",
  "cross-chromosome block length", "spatial heterogeneity (hotspot)",
  "cross-chromosome, scalar", "other scalar")), statistic)],
  file.path("tables", "summary_statistics.csv"))

cat(sprintf("\nKerman hap2 inside the reference-table central 90%%: %d of %d\n",
            sum(cat_dt$in_ref_90), nrow(cat_dt)))
cat("outside:\n")
print(cat_dt[in_ref_90 == FALSE, .(statistic, role, kerman_hap2,
                                   ref_med, percentile)])
cat("\nwrote", file.path(SIM_OUT_DIR, "fig_g10_statset.{pdf,png}"),
    "and", file.path("tables", "summary_statistics.csv"), "\n")

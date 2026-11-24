#!/usr/bin/env Rscript
# ============================================================================
# g10 main figure, COMPACT variant -- fingerprints + posteriors only.
#
# The paper applies the model to one genome, P. vera 'Kerman' hap2 (Fig. 5e).
#   a) fingerprints, five display chromosomes
#      (Chr02/03/04/05/07): real hap2 | simulated, one row. The simulation is
#      drawn at hap2's ABC-RF POSTERIOR MEDIANS (rep4 of 15, cache
#      sim_g10pick_kermanhap2_rep*.rds built by R/06_display_simulations.R,
#      seeds 202000+rep, best by the
#      pair-count / inverted-fraction / far-pair ranking) -- the same point as
#      the filled dots in panel b.
#   b) hap2's nine posterior estimates (violins, medians) with the mechanism
#      icons and real-time calibration.
#
# Out: figures/fig_g10_main_compact_hap2.{pdf,png}
#      <SIM_DATA_DIR>/results/fig_g10_main_compact_hap2_posteriors.csv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(reCYCLing); library(ggplot2); library(patchwork)
})
SRC_DIR <- local({
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f)) dirname(sub("^--file=", "", .f[1])) else "R"
})
source(file.path(SRC_DIR, "abc_summary_stats.R"))
## Run from the repository root:
##   Rscript R/07_fig5e.R
stopifnot(file.exists("config.R"))   # run from the repository root
source("config.R")
RES <- file.path(SIM_DATA_DIR, "results")
TS_CAP <- 75000L; MAXPAIR_CHROM <- 10000L; BASE <- 6
NSHOW <- 5L

files <- Sys.glob(file.path(DATA_DIR, "knobs", "Kerman_hap2_Chr*_*.fasta"))
if (!length(files)) stop("no Kerman hap2 monomer FASTAs under ", DATA_DIR,
                         "/knobs/ -- see data/README.md")
real <- assemble_genome_dt(files)
chrom_labs <- sort(unique(sub(".*_(Chr\\d+)_.*", "\\1", basename(files))))
sizes_real <- real[, .N, by = chrom][order(chrom)]
post <- fread(file.path(RES, "fit_g10_posteriors.csv"))[genotype == "Kerman" & hap == "hap2"]

sim <- readRDS(file.path(RES, "sim_g10pick_kermanhap2_rep4.rds"))

# ---- displayed arrays ---------------------------------------------------------
SHOW_CHR <- c("Chr02", "Chr03", "Chr04", "Chr05", "Chr07")
show_idx <- match(SHOW_CHR, chrom_labs)
stopifnot(!anyNA(show_idx), length(show_idx) == NSHOW)

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

L_REAL <- "real (Kerman hap2)"; L_SIM <- "simulated (fitted model)"
set.seed(1)
fp_of <- function(dt, who) {
  d  <- dt[chrom %in% show_idx]
  pp <- pair_positions(d)
  sz <- d[, .(n = max(num)), by = chrom]
  pp <- merge(pp, sz, by = "chrom")
  pp[, `:=`(x = lo / n, y = hi / n, source = who)][]
}
fp <- rbind(fp_of(real, L_REAL), fp_of(sim, L_SIM))
fp_plot <- fp[, if (.N > MAXPAIR_CHROM) .SD[sample.int(.N, MAXPAIR_CHROM)] else .SD,
              by = .(source, chrom)]

# real panels are labelled by chromosome + size (in-triangle); sim panels by
# size ONLY -- simulated chromosomes have no identity, only a realized size
sim_sizes <- sim[chrom %in% show_idx, .(n = max(num)), by = chrom][
  order(match(chrom, show_idx))]
REAL_LABS <- sprintf("%s\n%sk", chrom_labs[show_idx],
                     round(sizes_real$N[show_idx] / 1000))
SIM_LABS  <- sprintf("%sk", round(sim_sizes$n / 1000))

# Chrom/size labels go INSIDE the empty lower-right triangle of each panel
# (points sit above the diagonal); real/simulated identity is a TITLE above the
# block (uses the band that would otherwise be white above the squares).
fp_panel <- function(who, title, title_col, labs_vec, subtitle = NULL) {
  d <- copy(fp_plot[source == who])
  d[, chrom_lab := factor(chrom, levels = show_idx)]
  labdt <- data.table(chrom_lab = factor(show_idx, levels = show_idx),
                      lab = labs_vec)
  ggplot(d, aes(x, y, colour = inverted)) +
    geom_point(size = 0.03, alpha = 0.5, shape = 16, show.legend = FALSE) +
    geom_text(data = labdt, inherit.aes = FALSE, aes(x = 0.97, y = 0.05,
              label = lab), hjust = 1, vjust = 0, lineheight = 0.85,
              size = 4.6 / .pt, family = "Helvetica", colour = "grey20") +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#c23b22")) +
    facet_grid(. ~ chrom_lab) +
    scale_x_continuous(expand = c(0.02, 0), breaks = NULL) +
    scale_y_continuous(expand = c(0.02, 0), breaks = NULL) +
    coord_fixed(clip = "off") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = BASE, base_family = "Helvetica") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.25),
          panel.spacing = unit(1, "pt"),
          strip.text = element_blank(),
          plot.title = element_text(size = 6.5, face = "bold", family = "Helvetica",
                                    hjust = 0.5, colour = title_col,
                                    margin = margin(b = 0.5)),
          plot.subtitle = element_text(size = 5, family = "Helvetica",
                                       hjust = 0.5, colour = "grey35",
                                       margin = margin(b = 0.5)),
          plot.margin = margin(0, 2, 0, 2))
}
# tag theme must be set on the INNER plots -- the outer `&` cannot reach
# inside wrap_elements(), which is how the two letters ended up different sizes
TAG_THEME <- theme(plot.tag = element_text(size = 8, face = "bold",
                                           family = "Helvetica"),
                   plot.tag.position = c(0, 1))
p_real <- fp_panel(L_REAL, "real (Kerman hap2)", "black", REAL_LABS) +
  theme(plot.margin = margin(0, 4, 0, 1))
p_sim  <- fp_panel(L_SIM, "simulated (example)", "#008B8B", SIM_LABS) +
  theme(plot.margin = margin(0, 1, 0, 2))

# note about chromosome coverage -- a REAL-genome fact, so it sits to the LEFT
# of the real block (right-aligned against it), NOT next to the simulation
p_note <- ggplot() +
  annotate("text", x = 0.95, y = 0.5, hjust = 1, vjust = 0.5, lineheight = 0.95,
           size = 4.8 / .pt, family = "Helvetica", colour = "grey35",
           label = "knob arrays\non 11 of 15\nchromosomes\n(5 shown)") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void()

# ---- b) posterior estimates inside their priors -----------------------------
PRIORS <- data.table(
  par = c("age", "p_local", "local_mean", "alpha", "p_distal", "distal_mean",
          "p_xch", "xch_mean", "p_conv"),
  lab = c("array age (clock steps)",
          "local duplication rate (per monomer)",
          "local duplication size (monomers)",
          "hotspot feedback exponent",
          "long-range duplication rate (per chromosome)",
          "long-range duplication size (monomers)",
          "between-chromosome exchange rate (per genome)",
          "between-chromosome chunk size (monomers)",
          "gene conversion rate (per monomer)"),
  lo    = c(150,  2.5e-5, 1,   0.1, 1e-5, 1,     1e-3, 5,    1e-5),
  hi    = c(1500, 2.5e-3, 300, 5,   0.2,  20000, 0.95, 3000, 2e-3),
  logsc = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  grp   = c("age", "loc", "loc", "hot", "lr", "lr", "xch", "xch", "conv"))
PRIORS[, ypos := .N:1]

resc <- function(x, lo, hi, logsc)
  if (logsc) (log10(x) - log10(lo)) / (log10(hi) - log10(lo)) else (x - lo) / (hi - lo)

long_of <- function(d, cols) rbindlist(lapply(seq_len(nrow(PRIORS)), function(i) {
  pr <- PRIORS[i]
  data.table(par = pr$par, ypos = pr$ypos, grp = pr$grp,
             v = resc(d[[paste0(pr$par, cols)]], pr$lo, pr$hi, pr$logsc))
}))
pf <- cbind(PRIORS,
  med  = long_of(post, "_med")$v,
  lo95 = long_of(post, "_q025")$v,
  hi95 = long_of(post, "_q975")$v)
# round, never truncate (formatC format="d" truncates)
fmt <- function(x) trimws(
  ifelse(x >= 100, formatC(round(x), format = "d", big.mark = ","),
  ifelse(x >= 1, formatC(x, format = "fg", digits = 3),
         formatC(x, format = "g", digits = 2))))
pf[, raw := vapply(seq_len(.N), function(i) post[[paste0(par[i], "_med")]], 0)]
pf[, vlab := fmt(raw)]

source(file.path(SRC_DIR, "calibrate_realtime.R"))
ANCHORS <- c(GENUS_AGE$crown, GENUS_AGE$crown_hpd_low, GENUS_AGE$crown_hpd_high)
YPG <- ANCHORS / post$age_med
DENOM <- c(age = NA, p_local = post$n_mono, local_mean = NA, alpha = NA,
           p_distal = post$n_chrom, distal_mean = NA,
           p_xch = 1, xch_mean = NA, p_conv = post$n_mono)
# compact form: unit chosen by the median, range numbers unitless
fmt_n <- function(x) ifelse(x < 10, sprintf("%.1f", x), sprintf("%.0f", x))
fmt_rng <- function(y) {  # y sorted c(lo, med, hi) in years
  d <- if (y[2] >= 1e6) 1e6 else if (y[2] >= 1e3) 1e3 else 1
  u <- if (y[2] >= 1e6) "Myr" else if (y[2] >= 1e3) "kyr" else "yr"
  sprintf("1 per %s %s (%s-%s)", fmt_n(y[2] / d), u,
          fmt_n(y[1] / d), fmt_n(y[3] / d))
}
pf[, rtlab := vapply(seq_len(.N), function(i) {
  d <- DENOM[[par[i]]]
  if (par[i] == "age")
    return(sprintf("%.2f Myr (%.1f-%.2f); Han et al. 2022",
                   ANCHORS[1] / 1e6, ANCHORS[2] / 1e6, ANCHORS[3] / 1e6))
  if (is.na(d)) return("")
  fmt_rng(sort(YPG / (raw[i] * d)))
}, character(1))]

out_cols <- pf[, .(par, lab, prior_lo = lo, prior_hi = hi,
                   kerman_hap2_med = vlab,
                   scaled_med = round(med, 3), scaled_q025 = round(lo95, 3),
                   scaled_q975 = round(hi95, 3), real_time = rtlab)]
fwrite(out_cols, file.path(RES, "fig_g10_main_compact_hap2_posteriors.csv"))

GCOLS <- c(age = "grey25", loc = "#1a9850", hot = "#7b3294", lr = "#d95f02",
           xch = "#2166ac", conv = "#e7298a")
YLIM <- c(1 - 0.45, 9 + 1.05)   # headroom so the column titles clear the top violin

# ABC-RF posterior DENSITY (weighted reference sample) -> violin per parameter.
# Precomputed by R/05_posterior_densities.R; fall back to CI
# bars if absent.
DENS_F <- file.path(RES, "fig_g10_main_compact_hap2_postdensity.csv")
HAS_DENS <- file.exists(DENS_F)
if (HAS_DENS) {
  dens <- fread(DENS_F)
  dens <- merge(dens, PRIORS[, .(par, ypos, grp)], by = "par")
  dens <- dens[draw >= -0.02 & draw <= 1.02]   # clip to the plotted prior span
}

p_i <- ggplot(pf, aes(y = ypos)) +
  geom_segment(aes(x = 0, xend = 1, yend = ypos), colour = "grey90",
               linewidth = 1.4, lineend = "butt") +
  { if (HAS_DENS) geom_violin(data = dens,
      aes(x = draw, y = ypos, group = ypos, fill = grp, colour = grp),
      inherit.aes = FALSE, orientation = "y", width = 0.72, linewidth = 0.12,
      alpha = 0.5, scale = "width", trim = TRUE, bw = 0.03) } +
  geom_point(aes(x = med, colour = grp), size = 0.9, shape = 16) +
  annotate("segment", x = -0.005, xend = 1.005, y = -Inf, yend = -Inf,
           linewidth = 0.3) +
  geom_text(aes(x = 1.04, label = vlab), hjust = 0, size = 5 / .pt,
            family = "Helvetica") +
  geom_text(aes(x = 1.24, label = rtlab), hjust = 0, size = 5 / .pt,
            family = "Helvetica", colour = "grey20") +
  annotate("text", x = 0.5, y = 9 + 0.72, hjust = 0.5, size = 5 / .pt,
           family = "Helvetica", fontface = "bold", colour = "grey20",
           label = "random-forest ABC posterior") +
  annotate("text", x = c(1.04, 1.24), y = 9 + 0.72,
           hjust = 0, size = 5 / .pt, family = "Helvetica", fontface = "bold",
           label = c("model", "real time (approx.)")) +
  scale_colour_manual(values = GCOLS, guide = "none") +
  scale_fill_manual(values = GCOLS, guide = "none") +
  scale_x_continuous(limits = c(-0.005, 2.16), breaks = c(0, 0.5, 1),
                     labels = c("prior min", "", "prior max"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = PRIORS$ypos, labels = PRIORS$lab,
                     limits = YLIM, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = BASE, base_family = "Helvetica") +
  theme(axis.text = element_text(size = 5, family = "Helvetica", colour = "black"),
        axis.text.y = element_text(hjust = 1),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(2, 2, 0, 0))

# ---- b) mechanism cartoon legend: the reCYCLing model's mechanisms, one row --
# Layout: far-left blank gutter (for Grey's hand-drawn "example" arrow), then a
# "reCYCLing model parameters" label + arrow pointing right into 7 labelled
# icons. Colours match the posterior-row groups (deletion is calibrated, not a
# fitted posterior, so it has no row -- shown here for completeness).
LT <- function(hex) { # lighten toward white for the "copy" cell
  v <- grDevices::col2rgb(hex); grDevices::rgb(t(255 - (255 - v) * 0.45), maxColorValue = 255)
}
ARW   <- arrow(length = unit(2.2, "pt"), type = "closed")
ARW_S <- arrow(length = unit(1.7, "pt"), type = "closed")
DEL_COL <- "#8c510a"
IC_DX <- 0.088; IC_W <- 0.072
hstrip <- function(xc, yc, n = 6L, h = 0.26)
  data.table(cell = seq_len(n),
             xmin = xc - (n * IC_DX - (IC_DX - IC_W)) / 2 + (seq_len(n) - 1) * IC_DX,
             ymin = yc - h / 2, ymax = yc + h / 2)[, xmax := xmin + IC_W][]
cx <- function(s, i) (s$xmin[i] + s$xmax[i]) / 2
rects <- list(); mark <- function(s, cells, cols) {
  s[, fill := "grey90"][, col := "grey55"]
  if (length(cells)) s[cells, `:=`(fill = cols, col = cols)]
  rects[[length(rects) + 1]] <<- s
  s
}
YS <- 0.55   # single-strip icon centre
IC_SP <- 0.76     # icon slot spacing (tightened; blank space falls on the right)
XC <- 0.5 + IC_SP * (0:6)   # seven icon slot centres
s_age <- mark(hstrip(XC[1], YS), 3L, "grey25")
# deletion: middle monomer excised (white, dashed), neighbours close the gap
s_del <- mark(hstrip(XC[2], YS), integer(0), character(0))
s_del[3L, `:=`(fill = "white", col = DEL_COL)]
s_loc <- mark(hstrip(XC[3], YS), c(2L, 3L, 4L, 5L),
              c("#1a9850", "#1a9850", LT("#1a9850"), LT("#1a9850")))
s_hot <- mark(hstrip(XC[4], YS), 3L, "#7b3294")
s_lr  <- mark(hstrip(XC[5], YS), c(1L, 2L, 5L, 6L),
              c("#d95f02", "#d95f02", LT("#d95f02"), LT("#d95f02")))
s_xt  <- mark(hstrip(XC[6], 0.72, h = 0.20), c(2L, 3L), rep("#2166ac", 2))
s_xb  <- mark(hstrip(XC[6], 0.38, h = 0.20), c(4L, 5L), rep(LT("#2166ac"), 2))
s_ct  <- mark(hstrip(XC[7], 0.72, h = 0.20), 4L, "#e7298a")
s_cb  <- mark(hstrip(XC[7], 0.38, h = 0.20), 4L, LT("#e7298a"))
icon_rects <- rbindlist(rects)
IC_LABS <- data.table(
  x = XC,
  lab = c("mutation", "deletion", "local dup.", "hotspot",
          "long-range dup.", "exchange", "conversion"),
  col = c("grey25", DEL_COL, "#1a9850", "#7b3294", "#d95f02",
          "#2166ac", "#e7298a"))
XLO <- -2.15  # left edge: small blank gutter for Grey's hand-drawn arrow

p_ic <- ggplot() +
  geom_rect(data = icon_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                   ymax = ymax, fill = fill, colour = col),
            linewidth = 0.25) +
  scale_fill_identity() + scale_colour_identity() +
  geom_text(data = IC_LABS, aes(x, y = 0.10, label = lab, colour = col),
            vjust = 1, size = 4.6 / .pt, family = "Helvetica") +
  # left label + arrow into the icons (kept close to the icons; only a small
  # blank gutter remains on the far left for Grey's hand-drawn arrow)
  annotate("text", x = -0.55, y = YS, hjust = 1, vjust = 0.5, lineheight = 0.9,
           size = 5.2 / .pt, fontface = "bold", family = "Helvetica",
           colour = "grey20", label = "reCYCLing model\nparameters:") +
  annotate("segment", x = -0.45, xend = -0.02, y = YS, yend = YS,
           linewidth = 0.35, colour = "grey20", arrow = ARW) +
  # mutation: substitution star over the hit monomer
  annotate("text", x = cx(s_age, 3), y = YS + 0.16, label = "*",
           size = 11 / .pt, family = "Helvetica", colour = "grey25",
           vjust = 0.3) +
  # deletion: excise cell 3 (drawn X, font lacks U+2715), neighbours close in
  annotate("segment", x = cx(s_del, 3) - 0.028, xend = cx(s_del, 3) + 0.028,
           y = YS - 0.085, yend = YS + 0.085, linewidth = 0.3, colour = DEL_COL) +
  annotate("segment", x = cx(s_del, 3) - 0.028, xend = cx(s_del, 3) + 0.028,
           y = YS + 0.085, yend = YS - 0.085, linewidth = 0.3, colour = DEL_COL) +
  annotate("segment", x = cx(s_del, 1), xend = cx(s_del, 3) - 0.03, y = YS + 0.22,
           yend = YS + 0.22, linewidth = 0.3, colour = DEL_COL, arrow = ARW_S) +
  annotate("segment", x = cx(s_del, 6), xend = cx(s_del, 3) + 0.03, y = YS + 0.22,
           yend = YS + 0.22, linewidth = 0.3, colour = DEL_COL, arrow = ARW_S) +
  # local duplication: short arc, chunk centre to chunk centre
  annotate("curve", x = (cx(s_loc, 2) + cx(s_loc, 3)) / 2, y = YS + 0.16,
           xend = (cx(s_loc, 4) + cx(s_loc, 5)) / 2, yend = YS + 0.16,
           curvature = -0.6, linewidth = 0.3, colour = "#1a9850", arrow = ARW) +
  # hotspot: two arrows converging on the same target
  annotate("curve", x = cx(s_hot, 1), y = YS + 0.40,
           xend = cx(s_hot, 3) - 0.02, yend = YS + 0.20,
           curvature = -0.3, linewidth = 0.3, colour = "#7b3294", arrow = ARW) +
  annotate("curve", x = cx(s_hot, 5), y = YS + 0.40,
           xend = cx(s_hot, 3) + 0.02, yend = YS + 0.20,
           curvature = 0.3, linewidth = 0.3, colour = "#7b3294", arrow = ARW) +
  # long-range duplication: long arc, chunk centre to chunk centre
  annotate("curve", x = (cx(s_lr, 1) + cx(s_lr, 2)) / 2, y = YS + 0.16,
           xend = (cx(s_lr, 5) + cx(s_lr, 6)) / 2, yend = YS + 0.16,
           curvature = -0.35, linewidth = 0.3, colour = "#d95f02", arrow = ARW) +
  # exchange: chunks swap between the two chromosomes (double-headed arrow)
  annotate("segment", x = cx(s_xt, 2) + 0.02, y = 0.72 - 0.12,
           xend = cx(s_xb, 5) - 0.02, yend = 0.38 + 0.12,
           linewidth = 0.3, colour = "#2166ac",
           arrow = arrow(length = unit(1.7, "pt"), type = "closed",
                         ends = "both")) +
  # conversion: donor overwrites the recipient in place
  annotate("segment", x = cx(s_ct, 4), y = 0.72 - 0.11,
           xend = cx(s_cb, 4), yend = 0.38 + 0.13,
           linewidth = 0.3, colour = "#e7298a", arrow = ARW_S) +
  # right gap (icons end at x~6.9; extend to 8.2) leaves blank space for arrows
  scale_x_continuous(limits = c(XLO, 8.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Helvetica") +
  theme(plot.margin = margin(1, 2, 1, 2))

# ---- c) two summary-statistic fits beside the posteriors --------------------
# Observed (grey bars) vs simulated (pink step, mean of the top-3 RF-median
# replicates by the similarity ranking). Same binning as fig_g10_main_single.R.
DBRKS <- seq(0, log10(TS_CAP), length.out = 41)
hist_dist <- function(pp) {
  ld <- pmin(pmax(log10(pp$hi - pp$lo), 0), log10(TS_CAP) - 1e-9)
  ct <- tabulate(findInterval(ld, DBRKS, rightmost.closed = TRUE), 40L)
  data.table(bin = (DBRKS[-1] + DBRKS[-41]) / 2, p = ct / sum(ct))
}
DIV_CAP <- 22L
per_chrom_div <- function(dt) {
  rbindlist(lapply(split(dt, by = "chrom"), function(g) {
    g  <- as.data.table(g)
    al <- if ("aligned" %in% names(g)) g$aligned else NULL
    g[, div := per_monomer_div(bponly, al)][]
  }))[is.finite(div) & div <= DIV_CAP]
}
hist_div <- function(dv) {
  ct <- tabulate(pmin(floor(dv$div), DIV_CAP - 1L) + 1L, DIV_CAP)
  data.table(bin = 0:(DIV_CAP - 1L), p = ct / sum(ct))
}
HREPS <- c(4L, 9L, 1L)
hsims <- lapply(HREPS, function(r)
  readRDS(file.path(RES, sprintf("sim_g10pick_kermanhap2_rep%d.rds", r))))
hr_d <- hist_dist(pair_positions(real))
hs_d <- rbindlist(lapply(hsims, function(s)
  hist_dist(pair_positions(s))))[, .(p = mean(p)), by = bin]
hr_v <- hist_div(per_chrom_div(real))
hs_v <- rbindlist(lapply(hsims, function(s)
  hist_div(per_chrom_div(s))))[, .(p = mean(p)), by = bin]
fwrite(rbind(cbind(hr_d, stat = "dist", who = "real"),
             cbind(hs_d, stat = "dist", who = "sim"),
             cbind(hr_v, stat = "div", who = "real"),
             cbind(hs_v, stat = "div", who = "sim")),
       file.path(RES, "fig_g10_main_compact_hap2_hists.csv"))

SIMCOL <- "#008B8B"   # teal: distinct from the pink conversion mechanism colour
mini_theme <- theme_classic(base_size = BASE, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 5, face = "bold", family = "Helvetica",
                                  margin = margin(b = 1)),
        axis.text = element_text(size = 5, family = "Helvetica"),
        axis.title = element_text(size = 5, family = "Helvetica"),
        legend.position = "none", plot.margin = margin(1, 3, 0, 3))
p_d <- ggplot() +
  geom_col(data = hr_d, aes(bin, p), fill = "grey78",
           width = diff(hr_d$bin[1:2])) +
  geom_step(data = hs_d, aes(bin - diff(hr_d$bin[1:2]) / 2, p),
            colour = SIMCOL, linewidth = 0.35, direction = "hv") +
  scale_x_continuous(breaks = 0:4, labels = c("1", "10", "100", "1k", "10k"),
                     expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  labs(x = "pair distance (monomers)", y = "proportion",
       title = "within-chrom pair distance") +
  mini_theme
p_e <- ggplot() +
  geom_col(data = hr_v, aes(bin, p), fill = "grey78", width = 1) +
  geom_step(data = hs_v, aes(bin - 0.5, p),
            colour = SIMCOL, linewidth = 0.35, direction = "hv") +
  scale_x_continuous(breaks = seq(0, 20, 5), expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  labs(x = "substitutions", y = "proportion",
       title = "distance to consensus") +
  mini_theme

row1 <- p_note + p_real + p_sim + plot_layout(nrow = 1, widths = c(0.20, 1, 1))
row3 <- p_i + p_d + p_e + plot_layout(nrow = 1, widths = c(3.95, 0.85, 0.85))
# 6.65 x 2.4 in (Grey's spec). Rows: a fingerprints / b icon legend /
# c posteriors + two stat-fit histograms.
fig <- wrap_elements(full = row1) /
       wrap_elements(full = p_ic) /
       wrap_elements(full = row3) +
  plot_layout(heights = c(0.80, 0.30, 1.30)) &
  theme(plot.margin = margin(0, 0, 0, 0))
ggsave(file.path(SIM_OUT_DIR, "fig_g10_main_compact_hap2.pdf"), fig, width = 6.65, height = 2.4,
       bg = "transparent")
ggsave(file.path(SIM_OUT_DIR, "fig_g10_main_compact_hap2.png"), fig, width = 6.65, height = 2.4,
       bg = "white", dpi = 400)
cat("done:", file.path(SIM_OUT_DIR, "fig_g10_main_compact_hap2.{pdf,png}"), "\n")

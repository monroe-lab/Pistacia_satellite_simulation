# =============================================================================
# calibrate_realtime.R
#
# Put the ABC posteriors on an absolute time axis.
#
# The simulator's molecular clock is FIXED (MU_FIX = 5e-5 substitutions per bp
# per model-generation; see reftable_g8.R). Nothing else in the model carries
# a time unit, so a "model generation" is *defined* as the interval in which
# 5e-5 subs/bp accumulate. Given any real substitution rate r (per bp per year),
#
#     years_per_generation = MU_FIX / r
#
# and every fitted rate converts by the reciprocal. Equivalently, the quantities
# the data actually identify are the dimensionless ratios age*mu and p/mu; r is
# the single external number that turns them into years.
#
# TWO DIRECTIONS ARE IMPLEMENTED, and they are NOT equally defensible:
#
#   forward_from_rate()  -- assume a literature r, read off an age.
#       Every published bulk-genome plant/tree rate gives an age OLDER than the
#       Pistacia genus itself, so this direction FALSIFIES the bulk clock rather
#       than dating the knob. Keep it in the paper as the falsification, not as
#       the estimate.
#
#   invert_from_age()    -- assume the knob is no older than the genus, and
#       back-derive the satellite's own substitution rate. This is the
#       defensible direction and the one to report.
#
# ALL CONSTANTS BELOW ARE VERIFIED AGAINST PRIMARY SOURCES (see provenance
# comments). Do not edit a number here without re-checking it -- these are the
# load-bearing external facts in the whole real-time argument.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ---- the model's fixed clock ------------------------------------------------
MU_FIX <- 5e-5    # subs / bp / model-generation. MUST match reftable_g8.R:57.

# ---- external calibration constants -----------------------------------------
# Bulk-genome somatic substitution rates in long-lived trees, per bp per YEAR.
# These bracket the published range: Satake is explicitly reported as faster
# than Quercus / Populus / Eucalyptus / Picea, so it is the fast end.
#
# Xian et al. 2025, G3 15(8):jkaf143. PMID 40543050. Hub key: Xian2025-gj.
#   Quercus robur ("Napoleon Oak"), 234 yr old. Verbatim from the PMC full text
#   (NOT the abstract, which only reports "a factor of 5"):
#   "Scaled to the age of the Napoleon Oak, 234 years, the annual mutation rate
#    would be approximately 2.5 x 10-10 bp-1"
#
# Satake et al. 2024, eLife 12:RP88456. PMID 39441734. Hub key: Satake2023-ci.
#   Shorea laevis (slow) and S. leprosula (fast), Borneo. Verbatim from the PMC
#   full text: "with 7.71x10-10 and 8.05x10-10 nucleotide-1 year-1 for the slow-
#   and fast-growing species, respectively".
RATE_BULK <- list(
  oak_quercus   = 2.50e-10,   # Xian2025-gj
  shorea_slow   = 7.71e-10,   # Satake2023-ci
  shorea_fast   = 8.05e-10    # Satake2023-ci  <- fastest published bulk tree rate
)

# Pistacia genus divergence times, in years.
#
# Han et al. 2022, Front Plant Sci 13:1030647. PMID 36438104. Hub key: Han2022-jq.
#   Verbatim: "the stem and crown nodes of Pistacia were 37.74 Ma (95% highest
#   posterior density (HPD): 35.73-39.87) in the later of Eocene and 15.68 Ma
#   (95% HPD: 6.7-26.65) in the middle Miocene, respectively"
#   This is the ONLY source located for the CROWN age. (The vault previously
#   attributed these figures to Yi 2008 / Xie 2014 -- that was wrong.)
#
# Xie et al. 2014, Mol Phylogenet Evol 77:136-46. PMID 24747126. Hub key: Xie2014-mr.
#   Independent check on the stem/origin only. Verbatim: "Pistacia was estimated
#   to have originated at 37.60 mya (with 95% highest posterior density interval
#   (HPD): 25.42-48.51 mya)". Agrees with Han's 37.74 Ma stem.
GENUS_AGE <- list(
  stem_han       = 37.74e6,   # Han2022-jq, point
  stem_xie       = 37.60e6,   # Xie2014-mr, point (independent, concordant)
  crown          = 15.68e6,   # Han2022-jq, point
  crown_hpd_low  =  6.70e6,   # Han2022-jq
  crown_hpd_high = 26.65e6    # Han2022-jq
)

# =============================================================================
# Direction 1: assume a rate, get an age. (The falsification.)
# =============================================================================
#' @param age_gens fitted age in model generations
#' @param r        real substitution rate, subs / bp / year
forward_from_rate <- function(age_gens, r) {
  yr_per_gen <- MU_FIX / r
  list(yr_per_gen = yr_per_gen,
       age_years  = age_gens * yr_per_gen)
}

# =============================================================================
# Direction 2: assume an age, get the satellite's own rate. (The estimate.)
# =============================================================================
#' @param age_gens  fitted age in model generations
#' @param age_years assumed absolute age (e.g. GENUS_AGE$crown)
invert_from_age <- function(age_gens, age_years) {
  divergence <- MU_FIX * age_gens          # expected subs/site over the run
  r_sat      <- divergence / age_years     # subs / bp / year
  list(yr_per_gen    = age_years / age_gens,
       divergence    = divergence,
       rate_satellite = r_sat,
       fold_vs_fastest_bulk = r_sat / RATE_BULK$shorea_fast)
}

# =============================================================================
# Convert a whole posterior table to real time under one anchor.
#
# Rate-unit semantics were read off reCYCLing/src/sim_engine.cpp and are NOT
# uniform -- do not "simplify" this:
#   p_local_dup      per MONOMER     per generation  (sim_engine.cpp:151)
#   p_distal_dup     per CHROMOSOME  per generation  (sim_engine.cpp:207)
#   p_exchange       per GENOME      per generation  (sim_engine.cpp:528, :768; named p_translocation before reCYCLing 1.9.0)
#   p_conversion     per MONOMER     per generation  (g8; sim_engine.cpp:~207)
# so genome-wide event counts per generation need n_mono / n_chrom / 1
# respectively before dividing by years-per-generation.
# =============================================================================
add_realtime_cols <- function(post, anchor_years = GENUS_AGE$crown) {
  post <- as.data.table(copy(post))
  ypg  <- anchor_years / post$age_med         # years per model generation, per row

  post[, yr_per_gen    := ypg]
  post[, rate_satellite := MU_FIX / ypg]      # subs / bp / yr implied for this row
  post[, fold_vs_bulk   := rate_satellite / RATE_BULK$shorea_fast]

  # genome-wide events per year
  post[, ev_local_per_yr  := p_local_med  * n_mono  / ypg]
  post[, ev_distal_per_yr := p_distal_med * n_chrom / ypg]
  post[, ev_xch_per_yr    := p_xch_med              / ypg]

  # readable inverses: years between events
  post[, yr_per_local  := 1 / ev_local_per_yr]
  post[, yr_per_distal := 1 / ev_distal_per_yr]
  post[, yr_per_xch    := 1 / ev_xch_per_yr]

  # g8 only -- present once conversion is fitted
  if ("p_conv_med" %in% names(post)) {
    post[, ev_conv_per_yr := p_conv_med * n_mono / ypg]
    post[, yr_per_conv    := 1 / ev_conv_per_yr]
  }
  post[]
}

# =============================================================================
# Report. Run standalone:  Rscript R/calibrate_realtime.R [posteriors.csv]
# =============================================================================
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  csv  <- if (length(args) >= 1) args[1] else "results/fit_g10_posteriors.csv"
  post <- fread(csv)

  cat("== source:", csv, "|", nrow(post), "genotype-haplotypes ==\n\n")

  ref <- post[genotype == "Kerman" & hap == "hap2"]
  if (!nrow(ref)) ref <- post[1]
  cat(sprintf("reference row: %s %s | age_med = %.0f model gens | n_mono = %d | n_chrom = %d\n\n",
              ref$genotype, ref$hap, ref$age_med, ref$n_mono, ref$n_chrom))

  cat("-- direction 1: assume a published bulk tree rate, read off an age --\n")
  cat("   (Pistacia stem age is ~37.7 Ma; anything older than that is impossible)\n")
  for (nm in names(RATE_BULK)) {
    f <- forward_from_rate(ref$age_med, RATE_BULK[[nm]])
    cat(sprintf("   %-13s r = %.2e /bp/yr -> %8.0f yr/gen -> age %6.1f Myr  %s\n",
                nm, RATE_BULK[[nm]], f$yr_per_gen, f$age_years / 1e6,
                if (f$age_years > GENUS_AGE$stem_han) "IMPOSSIBLE" else "ok"))
  }

  cat("\n-- direction 2: anchor to the genus, back-derive the satellite's rate --\n")
  anchors <- list(`stem 37.74 Ma`      = GENUS_AGE$stem_han,
                  `crown HPD hi 26.65` = GENUS_AGE$crown_hpd_high,
                  `crown 15.68 Ma`     = GENUS_AGE$crown,
                  `crown HPD lo 6.70`  = GENUS_AGE$crown_hpd_low)
  for (nm in names(anchors)) {
    v <- invert_from_age(ref$age_med, anchors[[nm]])
    cat(sprintf("   %-19s -> %7.0f yr/gen | r_sat = %.2e /bp/yr | %5.2fx fastest bulk\n",
                nm, v$yr_per_gen, v$rate_satellite, v$fold_vs_fastest_bulk))
  }

  cat("\n-- genome-wide event rates at the crown anchor (15.68 Ma) --\n")
  rt <- add_realtime_cols(post, GENUS_AGE$crown)
  r1 <- rt[genotype == ref$genotype & hap == ref$hap]
  cat(sprintf("   local duplication      1 per %8.0f yr\n", r1$yr_per_local))
  cat(sprintf("   long-range within-chrom 1 per %8.0f yr\n", r1$yr_per_distal))
  cat(sprintf("   inter-chrom exchange    1 per %8.0f yr\n", r1$yr_per_xch))

  cat(sprintf("\n   across all %d haplotypes: inter-chrom exchange 1 per %.0f-%.0f yr (median %.0f)\n",
              nrow(rt), min(rt$yr_per_xch), max(rt$yr_per_xch), median(rt$yr_per_xch)))

  cat("\nCAVEAT: age prior is U[150, 1500] gens (reftable_g8.R:72). Rows whose\n",
      "age_q975 sits at the ceiling have prior-limited, not data-limited, upper\n",
      "bounds -- their year-scale CIs are artifacts. Widen the prior before\n",
      "publishing absolute ages with intervals.\n", sep = "")
}

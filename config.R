## Paths for the figure and scoring scripts. Run every script from the
## repository root, e.g.  Rscript R/09_fig5e.R
##   results/  simulation and inference outputs (seed it from inputs/, see README)
##   data/     Kerman hap2 monomer FASTAs, data/knobs/Kerman_hap2_Chr*_*.fasta.gz
##   figures/  figure output
SIM_DATA_DIR <- "."
DATA_DIR     <- Sys.getenv("PISTACIA_DATA", unset = "data")
SIM_OUT_DIR  <- "figures"
dir.create("results", showWarnings = FALSE)
dir.create(SIM_OUT_DIR, showWarnings = FALSE)

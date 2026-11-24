#!/bin/bash
# Install reCYCLing v1.9.2 into the user R library on Farm. It compiles C++
# code, so it runs as a job rather than on the login node:
#   sbatch farm/install_recycling.sh
#SBATCH --job-name=recycling_install
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=farm/logs/install_recycling_%j.out
#SBATCH --error=farm/logs/install_recycling_%j.out
source /etc/profile.d/modules.sh
module load R/4.4.2
set -euo pipefail
mkdir -p farm/logs
Rscript -e 'remotes::install_github("monroe-lab/reCYCLing", ref = "v1.9.2", upgrade = "never")'
Rscript -e 'library(reCYCLing); cat("installed reCYCLing", as.character(packageVersion("reCYCLing")), "\n")'

# 00_install_packages.R
# Run once. Installs every package used in the project.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("tidyverse", "ape", "pheatmap", "RColorBrewer", "rmarkdown", "knitr")
bioc_pkgs <- c("Biostrings", "msa")

for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)

message("All packages installed.")

# MD7115 package installer
# Run this script once before class.

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing_cran <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }
}

install_if_missing_bioc <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  }
}

install_if_missing_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
  }
}

cat("Step 1/4: installing package managers\n")
install_if_missing_cran(c("BiocManager", "remotes"))

cat("Step 2/4: installing CRAN packages\n")
cran_packages <- c(
  "tidyverse",
  "readr",
  "tibble",
  "shadowtext",
  "ggsci",
  "igraph",
  "ggraph",
  "tidygraph",
  "ggnetwork",
  "extrafont",
  "ggnewscale",
  "RColorBrewer"
)
install_if_missing_cran(cran_packages)

cat("Step 3/4: installing Bioconductor packages\n")
bioc_packages <- c(
  "phyloseq",
  "ggtree",
  "treeio",
  "ggtreeExtra",
  "simplifyEnrichment"
)
install_if_missing_bioc(bioc_packages)

cat("Step 4/4: installing GitHub packages\n")
install_if_missing_github("microbiomeViz", "lch14forever/microbiomeViz")
install_if_missing_github("r4projects", "jaspershen/r4projects")

required_packages <- c(
  "BiocManager",
  "remotes",
  cran_packages,
  bioc_packages,
  "microbiomeViz",
  "r4projects"
)

missing_after_install <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

cat("\nInstallation check:\n")
if (length(missing_after_install) == 0) {
  cat("All required packages are installed.\n")
} else {
  cat("These packages still need attention:\n")
  print(missing_after_install)
}

cat("\nYou can now open the project and run the class examples.\n")

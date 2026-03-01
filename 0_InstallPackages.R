required_packages <- c(
  "mgcv",
  "rstan",
  "rstanarm",
  "ggplot2",
  "stringr",
  "BMisc",
  "Matrix",
  "loo",
  "bayesplot",
  "dplyr",
  "RColorBrewer",
  "MASS",
  "lubridate",
  "arm",
  "gridExtra",
  "scales",
  "extRemes",
  "extrafont",
  "xtable",
  "texmex",
  "reshape",
  "ggpubr",
  "ggh4x",
  "patchwork",
  "parallel",
  "heatwaveR", 
  "reshape2"
)

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg)
  }
}

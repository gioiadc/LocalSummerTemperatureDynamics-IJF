## Load packages ----
library(mgcv)
library(rstan)
library(rstanarm)
library(ggplot2)
library(stringr)
library(BMisc)
library(Matrix)
library(loo)
library(bayesplot)
library(dplyr)
library(RColorBrewer)
library(MASS)
library(lubridate)
library(arm)
library(gridExtra)
library(scales)
library(extRemes)
library(extrafont)
library(xtable)
library(texmex)
library(reshape)
library(ggpubr)
library(ggh4x)
library(patchwork)
library(parallel)
library(heatwaveR)
library(cmdstanr)

theme_set(bayesplot::theme_default())

colors <- c(rgb(0,114,178, maxColorValue = 255),  
            rgb(230,159,0, maxColorValue = 255)) # Colors for the regimes

rdylgn_palette <- brewer.pal(11, "RdYlGn")
darker_red <- rdylgn_palette[1]
darker_green <- rdylgn_palette[11]

# reshape2 is required but do not load it 

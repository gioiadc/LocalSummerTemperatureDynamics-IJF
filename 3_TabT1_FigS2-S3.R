# This file reproduce results in the Supplementary Material:
# - the posterior parameter summaries and trace plots for Capriva (K = 5, last.regime = "N")

source("functions/load_Packages.R")

ttheme_default(base.size = 8)

extracted_chains <- readRDS("ExtractedResults/extracted_samples.rds")
summary_table <- rstan::monitor(extracted_chains, print = FALSE, warmup = 0)

TABLE <- xtable(round(as.matrix(summary_table)[,1:10],2))
TABLE

labels <- c(
  'betaL1_new[1]' = "alpha[11]",
  'betaL1_new[2]' = "alpha[12]",
  'betaL1_new[3]' = "alpha[13]",
  'betaL1_new[4]' = "alpha[14]",
  'betaL1_new[5]' = "alpha[15]",
  'betaL1_new[6]' = "alpha[16]",
  'betaS1_new[1]' = "gamma[1]",
  'ltauL1' = "tau[1]",
  'betaTp01_new[1]' = "nu[12]",
  'betaTp10_new[1]' = "nu[21]",
  'betaTp12_new[1]' = "nu[23]",
  'betaTp21_new[1]' = "nu[32]",
  'betaTp23_new[1]' = "nu[34]",
  'betaTp32_new[1]' = "nu[43]",
  'betaTp34_new[1]' = "nu[45]",  
  'betaTp43_new[1]' = "nu[54]"
)

R1.1 <- mcmc_trace(extracted_chains, regex_pars = c("betaL1_new", "betaS1_new"),
                   facet_args = list(ncol = 1, labeller = as_labeller(labels[1:7], label_parsed))) +
  guides(colour = "none")

R1.2 <- mcmc_trace(extracted_chains, regex_pars = c("ltauL1", "betaTp01_new", "betaTp10_new", "betaTp12_new",
                                                    "betaTp21_new", "betaTp23_new", "betaTp32_new",
                                                    "betaTp34_new", "betaTp43_new"),
                   facet_args = list(ncol = 1, labeller = as_labeller(labels[8:16], label_parsed))) +
  guides(colour = "none")

as_matrix_chains <- function(extracted_chains, pars) {
    all_pars <- dimnames(extracted_chains)[[3]]
    idx <- unlist(lapply(pars, function(p) grep(paste0("^", p), all_pars)))
    arr <- extracted_chains[, , idx, drop = FALSE]
    matrix(arr,
           nrow = dim(arr)[1] * dim(arr)[2],
           dimnames = list(NULL, all_pars[idx]))
  }
  
R1acf.1 <- mcmc_acf(as_matrix_chains(extracted_chains, c("betaL1_new", "betaS1_new")),
                      facet_args = list(ncol = 1, labeller = as_labeller(labels[1:7], label_parsed))) +
    guides(colour = "none")
  
R1acf.2 <- mcmc_acf(as_matrix_chains(extracted_chains, c("ltauL1", "betaTp01_new", "betaTp10_new", "betaTp12_new",
                                                           "betaTp21_new", "betaTp23_new", "betaTp32_new",
                                                           "betaTp34_new", "betaTp43_new")),
                      facet_args = list(ncol = 1, labeller = as_labeller(labels[8:16], label_parsed))) +
    guides(colour = "none")



g <- arrangeGrob(R1.1, R1acf.1, ncol =2)
ggsave("Plots/DiagnosticL1.pdf", g,  dpi = 600, width = 174, height = 234, units = "mm")
g <- arrangeGrob(R1.2, R1acf.2, ncol =2)
ggsave("Plots/DiagnosticL1_2.pdf", g,  dpi = 600, width = 174, height = 234, units = "mm")

# ==============================================================================
# 02. corHMM simmaps

# Uses the best-fitting corHMM hidden-rate model (selected via corHMMDredge) to
# generate stochastic character maps (SIMMAPs) for the evolution of eusociality
# across the bee phylogeny. 
#
# This script:
#   1. Loads the trait dataset and pruned ML phylogeny used for the eusociality
#      corHMM analyses.
#   2. Aligns the tree and trait data, removes species with unknown sociality
#      states, and retains only species scored as eusocial, non-eusocial, or
#      polymorphic.
#   3. Loads the best-supported corHMM model (2-rate hidden Markov model),
#      and uses it as the transition matrix for stochastic character mapping
#      with makeSimmap (rate.cat = 2; root state fixed as non-eusocial).
#   4. Summarizes simulated transition histories across all SIMMAPs to obtain
#      mean numbers of transitions between sociality states and their variance.
#
# Note: taxa with missing sociality scores were removed prior to SIMMAP,
# but were inferred in the original corHMM model fit.
# ==============================================================================

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
rm(list=ls())
wd <- "/Users/lenarh/Desktop/eusociality_evolution"
wd
setwd(wd)
results_wd <- file.path(wd, "results")
results_wd

library(corHMM)
library(dplyr)
# devtools::install_github("thej022214/corHMM")

#-------------------------------------------------------------------------------
# Organizing dataset
#-------------------------------------------------------------------------------
# Load traits and tree
traits <- read.csv("raw_data/sociality_scored_lasioglossum_inferred.csv")
# traits <- read.csv("raw_data/sociality_scored_lasioglossum_unknown.csv")
phy <- read.tree("raw_data/ML_beetree_4586tips.tre") 

# Drop tips in the tree that don’t match any in the trait data
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

# Find species shared by both datasets
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

# Remove unknown/NA species
dat$ASR_sociality[dat$ASR_sociality == "unknown"] <- NA
dat <- na.omit(dat)
head(dat)

# Reorder the data and tree to align species
dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]

# Keep relevant trait columns
dat <- dat[,c("tips", "ASR_sociality")]

nrow(dat) # 4586 species
Ntip(phy) # 4586 species

unique(dat$ASR_sociality)
# [1] "non-eusocial" "polymorphic"  NA     "eusocial"   

head(dat)

#-------------------------------------------------------------------------------
# Run simmaps
#-------------------------------------------------------------------------------
# Load corHMM results
load(file.path(results_wd, "corhmm_dredge_run1.Rsave"))
# load(file.path(results_wd, "corhmm_dredge_run2.Rsave"))
model <- dredge_sociality[[30]]$solution # 30th = best model
# model <- dredge_sociality[[21]]$solution

# Get simmap from corHMM solution
simmaps <- makeSimmap(tree = phy, 
                      data = dat, 
                      model = model, 
                      rate.cat = 2, 
                      nSim = 100, 
                      root.p = c(0, 1, 0),
                      nCores = 1)

# Load and save the summary
simmap_summaries <- lapply(simmaps, summarize_single_simmap)
summary_df <- summarize_transition_stats(simmap_summaries)

# Save summary data
print(summary_df)

write.csv(
  summary_df,
  file = file.path(results_wd, "transitions_summary_corhmm_run1.csv"),
  row.names = FALSE
)

# write.csv(
#   summary_df,
#   file = file.path(results_wd, "transitions_summary_corhmm_run2.csv"),
#   row.names = FALSE
# )

# Save summary figure
pdf(file = file.path(results_wd, "transitions_summary_plot_run1.pdf"),
    width = 15, height = 8)
# pdf(file = file.path(results_wd, "transitions_summary_plot_run2.pdf"),
#     width = 15, height = 8)
plot_transition_summary(simmap_summaries)
dev.off()









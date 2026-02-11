# ==============================================================================
# 01. Ancestral state reconstruction of eusociality in bees
# ==============================================================================
# Uses corHMM to model the evolution of eusociality in bees under hidden rate models, 
# selects the best-fitting model via AIC, 
# reconstructs ancestral states across the phylogeny,
# and plots a circular phylogeny with the ASR and trait states mapped onto it.

# Species with unknown sociality state will be inferred by the model,
# and polymorphic species are treated as a separate state
# (i.e., the trait state options are eusocial, non-eusocial, or polymorphic)

# run1: unknown Lasioglossum manually inferred as primitively eusocial
# run2: unknown Lasioglossum left as unknown and allowed to be inferred by the model

# ==============================================================================

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
rm(list=ls())
wd <- "/Users/lenarh/Desktop/eusociality_evolution"
setwd(wd)
results_wd <- file.path(wd, "results")

library(corHMM)
library(phytools)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(scales)
#devtools::install_github("thej022214/corHMM")

# Load utility functions
source("00_utility_functions.R")

# Load traits and tree (load different trait dataset for different runs)
# traits <- read.csv("raw_data/sociality_scored_lasioglossum_inferred.csv") # Run 1
traits <- read.csv("raw_data/sociality_scored_lasioglossum_unknown.csv") # Run 2
phy <- read.tree("raw_data/ML_beetree_4586tips.tre") 

# head(traits)
head(traits)
phy

#-------------------------------------------------------------------------------
# Organizing and preparing dataset
#-------------------------------------------------------------------------------
# Drop tips in the tree that don’t match any in the trait data
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

# Find species shared by both datasets
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

# Reorder the data and tree to align species
dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]

# Keep relevant trait columns
dat <- dat[,c("tips","broad_sociality","ASR_sociality")]
head(dat)

# Check that data and phylogeny contain same number of species
nrow(dat) # 4586 species
Ntip(phy) # 4586 species

unique(dat$broad_sociality)
# [1] "solitary"        "social"          "polymorphic"     "unknown"        
# [5] "eusocial"        "social parasite" "kleptoparasite" 

unique(dat$ASR_sociality)
# [1] "non-eusocial" "polymorphic"  "unknown"      "eusocial"   

# Recode "unknown" to "?" so that corHMM recognizes this as missing trait data
# to be inferred, rather than a separate trait state
dat$ASR_sociality[dat$ASR_sociality == "unknown"] <- "?"

unique(dat$ASR_sociality)
# [1] "non-eusocial" "polymorphic"  "?"            "eusocial"  

# Lets put it into the format corHMM expects:
corHMM_dat <- dat %>%
  dplyr::select(Species = tips, state = ASR_sociality)

head(corHMM_dat)

# Now the data is ready to be run through corHMM.
# Use corHMM_dat from here onwards for analyses.

#-------------------------------------------------------------------------------
# Run model selection with corHMMdredge
#-------------------------------------------------------------------------------
dredge_sociality <- corHMM::corHMMDredge(phy, corHMM_dat, 
                                         max.rate.cat = 2, # allows for 2 rate classes
                                         root.p = c(0, 1, 0), # codes the root state as non-eusocial
                                         n.cores = 60, 
                                         nstarts = 5)

# save(dredge_sociality, file = file.path(results_wd, "corhmm_dredge_run1.Rsave"))
save(dredge_sociality, file = file.path(results_wd, "corhmm_dredge_run2.Rsave"))


# Make model comparison table
corhmm_tbl_sociality <- corHMM:::getModelTable(dredge_sociality)

# write.csv(corhmm_tbl_sociality, file = file.path(results_wd, "corhmm_tbl_dredge_run1.csv"))
write.csv(corhmm_tbl_sociality, file = file.path(results_wd, "corhmm_tbl_dredge_run2.csv"))

# Load model comparison table
# corhmm_tbl_sociality <- read.csv(file.path(results_wd, "corhmm_tbl_dredge_run1.csv"), header = TRUE)
corhmm_tbl_sociality <- read.csv(file.path(results_wd, "corhmm_tbl_dredge_run2.csv"), header = TRUE)
corhmm_tbl_sociality

# Extract the transition rate matrix from the best model
# rates_mat <- dredge_sociality[[30]]$solution # for run 1, this was the 30th
rates_mat <- dredge_sociality[[21]]$solution # for run 2, this was the 21st but several were similar

# write.csv(rates_mat, file = file.path(results_wd, "transition_rates_corHMM_run1.csv"))
write.csv(rates_mat, file = file.path(results_wd, "transition_rates_corHMM_run2.csv"))

# Best model includes 2 rate classes (aka there is rate heterogeneity across the tree)

# Load transition matrix
# transition_matrix <- read.csv(file.path(results_wd, "transition_rates_corHMM_run1.csv"))
transition_matrix <- read.csv(file.path(results_wd, "transition_rates_corHMM_run2.csv"))
transition_matrix


#===============================================================================
# Ancestral State Reconstruction & Plotting
#===============================================================================
#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
WtQ <- function(Q, Weights){ # Weighted average of transition rates using AIC weights
  Weights <- Weights[!is.na(Q)]/sum(Weights[!is.na(Q)])
  AvgQ <- sum(Q[!is.na(Q)] * Weights)
  return(AvgQ)
}

getModelAvgRate <- function(file){ # Loads model set, computes AIC weights, and creates a model-averaged rate matrix
  load(file)
  rate.mat <- obj$res_ARD.ARD.2$index.mat
  AICc <- unlist(lapply(obj, function(x) x$AICc))
  AICwt <- exp(-0.5 * AICc - min(AICc))/sum(exp(-0.5 * AICc - min(AICc)))
  # Solutions <- lapply(obj, function(x) x$solution)
  Solutions <- lapply(obj, function(x) c(na.omit(c(x$solution))))
  Solutions[[1]] <-  c(Solutions[[1]][1], NA, Solutions[[1]][2], NA,
                       NA, Solutions[[1]][1], NA, Solutions[[1]][2])
  Solutions[[2]] <-  c(Solutions[[2]][1], NA, Solutions[[2]][2], NA,
                       NA, Solutions[[2]][1], NA, Solutions[[2]][2])
  Rates <- do.call(rbind, Solutions)
  p.wt <- apply(Rates, 2, function(x) WtQ(x, AICwt))
  rate.mat[!is.na(rate.mat)] <- p.wt
  return(rate.mat)
}

getTipRecon <- function(file){ # Reconstructs states using model-averaged transition matrix, returns marginal likelihoods
  load(file)
  phy <- obj$res_ER$phy
  data <- obj$res_ER$data
  root.p <- obj$res_ARD$root.p
  index.mat <- obj$res_ARD.ARD.2$index.mat
  p <- getModelAvgRate(file)[sapply(1:max(index.mat, na.rm = TRUE), function(x) match(x, index.mat))]
  res <- corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = index.mat, node.states = "marginal", p = p, root.p = root.p, get.tip.states = TRUE)
  return(res)
}

# Plots ancestral state reconstructions with colored pie charts at nodes
plotRECON <- function(phy, likelihoods, piecolors=NULL, cex=0.5, pie.cex=0.25, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
  if(is.null(piecolors)){
    piecolors=c("red","navyblue","darkgreen","pink","skyblue","lightgreen")
    # piecolors=c("#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF",
    #             "#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF")
    # piecolors=rev(c("#00204DFF", "#575C6DFF", "#A69D75FF", "#FFEA46FF",
    #                 "#00204DFF", "#575C6DFF", "#A69D75FF", "#FFEA46FF"))
  }
  if(!is.null(file)){
    pdf(file, height=height, width=width,useDingbats=FALSE)
  }
  plot(phy, cex=cex, show.tip.label=show.tip.label, ...)
  
  if(!is.null(title)){
    title(main=title)
  }
  nodelabels(pie=likelihoods,piecol=piecolors, cex=pie.cex)
  states <- colnames(likelihoods)
  legend(x="topleft", states, cex=0.8, pt.bg=piecolors,col="black",pch=21);
  
  if(!is.null(file)){
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# Plot ancestral state reconstruction
#-------------------------------------------------------------------------------
# Extract likelihoods for internal nodes from best model
# anc_recon <- dredge_sociality[[30]]$states
anc_recon <- dredge_sociality[[21]]$states
anc_recon[1,] # checking root state

# pdf(file.path(results_wd, "corhmm_dredge_recon_final_run1.pdf"), height = 45, width = 10)
pdf(file.path(results_wd, "corhmm_dredge_recon_final_run2.pdf"), height = 45, width = 10)

plotRECON(
  phy = phy,
  likelihoods = anc_recon,
  pie.cex = 0.3,
  show.tip.label = T,
  cex = 0.1
)

axisPhylo()
dev.off()

# ==============================================================================
# Circular phylogeny: ASR branch coloring + eusociality tip ring (+ families)
# ==============================================================================

# Make sure we're in the project directory
wd <- "/Users/lenarh/Desktop/eusociality_evolution"
setwd(wd)
results_wd <- file.path(wd, "results")

library(ape)
library(phytools)
library(viridis)
library(dplyr)

# ----------------------------------------------------------------------
# Load tree, trait data, and corHMM results (if not already in memory)
# ----------------------------------------------------------------------
phy <- read.tree("raw_data/ML_beetree_4586tips.tre")
# traits <- read.csv("raw_data/sociality_scored_lasioglossum_inferred.csv")
traits <- read.csv("raw_data/sociality_scored_lasioglossum_unknown.csv")

# Keep only species with trait data and reorder to match tree
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)
dat <- dat[match(shared_species, dat$tips), ]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips), ]

# Eusociality states for ASR
dat <- dat[, c("tips", "ASR_sociality")]
dat$ASR_sociality[dat$ASR_sociality == "unknown"] <- "?"
unique(dat$ASR_sociality)

# For the tip ring, treat unknowns as a fourth category
ring_state <- dat$ASR_sociality
ring_state[ring_state == "?"] <- "unknown"

# Named vector of tip states
sociality_vec <- setNames(ring_state, dat$tips)

# Colors for eusociality states
eusocial_cols <- c(
  "non-eusocial" = "#1b9e77",
  "polymorphic"  = "#d95f02",
  "eusocial"     = "#7570b3",
  "unknown"      = "grey80"
)

# ----------------------------------------------------------------------
# Load best corHMM model and extract ASR
# ----------------------------------------------------------------------
# load(file.path(results_wd, "corhmm_dredge_run1.Rsave"))
load(file.path(results_wd, "corhmm_dredge_run2.Rsave")) 
# anc_recon <- dredge_sociality[[30]]$states
anc_recon <- dredge_sociality[[21]]$states

# Column names in anc_recon should be like:
# "R1 eusocial", "R1 non-eusocial", "R1 polymorphic",
# "R2 eusocial", "R2 non-eusocial", "R2 polymorphic"

colnames(anc_recon)

# Map hidden-rate states to eusociality states
asr_state_key <- c(
  "R1 eusocial"      = "eusocial",
  "R1 non-eusocial"  = "non-eusocial",
  "R1 polymorphic"   = "polymorphic",
  "R2 eusocial"      = "eusocial",
  "R2 non-eusocial"  = "non-eusocial",
  "R2 polymorphic"   = "polymorphic"
)

# Vector of node states: pick max-likelihood hidden state if max > 0.5
n_tips <- length(phy$tip.label)
node_states <- apply(anc_recon, 1, function(x) {
  if (max(x) > 0.5) names(which.max(x)) else NA
})

# Initialize branch colors
branch_colors <- rep("gray80", nrow(phy$edge))

# Recolor internal branches according to ASR state
for (i in 1:nrow(phy$edge)) {
  parent_node <- phy$edge[i, 1]
  if (parent_node > n_tips) {
    st <- node_states[parent_node - n_tips]
    if (!is.na(st)) {
      trait_label <- asr_state_key[st]
      branch_colors[i] <- eusocial_cols[trait_label]
    }
  }
}

# ----------------------------------------------------------------------
# Helper functions: polar/cartesian + radial time axis
# ----------------------------------------------------------------------
toCart <- function(r, th, deg = FALSE) {
  if (deg) th <- th * pi / 180
  list(x = r * cos(th), y = r * sin(th))
}

toPolar <- function(x, y) {
  list(r = sqrt(x^2 + y^2), th = atan2(y, x))
}

add_time_axis_radial <- function(tree, ring_times,
                                 angle_deg = 0,
                                 tick_length = 0.015,
                                 label_cex = 0.4) {
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tips <- length(tree$tip.label)
  root_radius <- max(sqrt(obj$xx^2 + obj$yy^2)[(n_tips + 1):length(obj$xx)])
  tree_height <- max(nodeHeights(tree))
  # Convert times to radii assuming ring_times are from present (0) to root
  reversed_radii <- (tree_height - ring_times) / tree_height * root_radius
  
  theta <- angle_deg * pi / 180
  dx <- cos(theta); dy <- sin(theta)
  segments(0, 0, root_radius * dx, root_radius * dy, col = "black", lwd = 1)
  
  for (i in seq_along(reversed_radii)) {
    r <- reversed_radii[i]
    x_tick <- r * dx; y_tick <- r * dy
    perp_dx <- -dy; perp_dy <- dx
    
    # Tick marks
    tick_x0 <- x_tick - tick_length * root_radius * perp_dx
    tick_y0 <- y_tick - tick_length * root_radius * perp_dy
    tick_x1 <- x_tick + tick_length * root_radius * perp_dx
    tick_y1 <- y_tick + tick_length * root_radius * perp_dy
    segments(tick_x0, tick_y0, tick_x1, tick_y1, col = "black", lwd = 1)
    
    label_offset <- 2 * tick_length * root_radius
    label_x <- x_tick + label_offset * perp_dx
    label_y <- y_tick + label_offset * perp_dy
    
    label_text <- if (ring_times[i] == 0) "0" else as.character(ring_times[i])
    text(label_x, label_y, labels = label_text, cex = label_cex, adj = c(0.5, 0))
  }
}

# ----------------------------------------------------------------------
# Plot circular phylogeny with ASR branches + eusociality tip ring
# ----------------------------------------------------------------------
# pdf(file.path(results_wd, "eusociality_ASR_circular_phylogeny_run1.pdf"),
#     width = 8, height = 8)
pdf(file.path(results_wd, "eusociality_ASR_circular_phylogeny_run2.pdf"),
    width = 8, height = 8)

par(mar = c(2, 2, 2, 2), xpd = TRUE)

plot(phy, type = "fan", edge.color = branch_colors, edge.width = 0.6,
     show.tip.label = FALSE, no.margin = FALSE,
     open.angle = 15, rotate.tree = 7.5, add = TRUE)

# Optional: time axis (here using relative units, e.g. 0–100)
ring_times <- c(0, 20, 40, 60, 80, 100)
add_time_axis_radial(phy, ring_times, angle_deg = 0)

# ----------------------------------------------------------------------
# Eusociality ring at tips
# ----------------------------------------------------------------------
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips_xx <- obj$xx[1:n_tips]
tips_yy <- obj$yy[1:n_tips]

ri    <- 0.2 # radial increment
len   <- 15 # thickness of colored band
space <- 4 # inner + outer white space

for (i in 0:(len + 2 * space)) {
  p <- toPolar(tips_xx, tips_yy)
  p$r <- p$r + ri * i
  c <- toCart(p$r, p$th)
  
  if (i < space || i >= space + len) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.15)
  } else {
    tip_cols <- eusocial_cols[sociality_vec[phy$tip.label]]
    points(c$x, c$y, col = tip_cols, pch = 16, cex = 0.25)
  }
}


# ------------------------------------------------------------------------------
# Tip labels at intervals
# ------------------------------------------------------------------------------

obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)

n_tips <- length(phy$tip.label)
tips_xx <- obj$xx[1:n_tips]
tips_yy <- obj$yy[1:n_tips]

label_every <- 5
tip_idx <- seq(1, n_tips, by = label_every)

r_tip <- sqrt(tips_xx^2 + tips_yy^2)

# Push labels outside the trait ring; increase if still crowded
r_out <- max(r_tip) + ri * (space + len) + 0.03 * max(r_tip)

theta <- atan2(tips_yy[tip_idx], tips_xx[tip_idx])
x_lab <- r_out * cos(theta)
y_lab <- r_out * sin(theta)

rot <- theta * 180 / pi
rot_flip <- ifelse(rot > 90 | rot < -90, rot + 180, rot)

pos_vec <- ifelse(cos(theta) >= 0, 4, 2)  # right vs left

tip_cex <- 0.20

for (k in seq_along(tip_idx)) {
  text(
    x = x_lab[k],
    y = y_lab[k],
    labels = phy$tip.label[tip_idx[k]],
    cex = tip_cex,
    srt = rot_flip[k],
    pos = pos_vec[k],
    offset = 0
  )
}


# ----------------------------------------------------------------------
# Family labels + MRCAs (only if 'family' column exists)
# ----------------------------------------------------------------------
if ("family" %in% names(traits)) {
  
  # Re-align traits to current tip order
  traits_fam <- traits[match(phy$tip.label, traits$tips), ]
  
  family_list <- unique(traits_fam$family)
  
  family_nodes <- sapply(family_list, function(fam) {
    spp <- traits_fam %>% dplyr::filter(family == fam) %>% pull(tips)
    tryCatch(findMRCA(phy, spp, type = "node"), error = function(e) NA)
  })
  
  family_nodes <- family_nodes[!is.na(family_nodes)]
  
  r_ring  <- max(sqrt(tips_xx^2 + tips_yy^2)) + ri * (space + len)
  r_label <- r_ring + 0.02 * r_ring
  
  for (i in seq_along(family_nodes)) {
    fam <- names(family_nodes)[i]
    fam_tips <- which(traits_fam$family == fam)
    fam_tip_labels  <- traits_fam$tips[fam_tips]
    fam_tip_indices <- match(fam_tip_labels, phy$tip.label)
    x_vals <- obj$xx[fam_tip_indices]
    y_vals <- obj$yy[fam_tip_indices]
    thetas <- atan2(y_vals, x_vals)
    mean_theta <- atan2(mean(sin(thetas)), mean(cos(thetas)))
    x_out <- r_label * cos(mean_theta)
    y_out <- r_label * sin(mean_theta)
    hjust <- ifelse(cos(mean_theta) >= 0, 0, 1)
    text(x_out, y_out, labels = fam, cex = 0.6, srt = 0, adj = c(hjust, 0.5))
  }
  
  mrca_color <- "darkgrey"
  for (node in family_nodes) {
    points(obj$xx[node], obj$yy[node], pch = 16, cex = 1, col = mrca_color)
  }
}

# ----------------------------------------------------------------------
# Legend
# ----------------------------------------------------------------------
legend("topleft",
       legend = c("non-eusocial", "polymorphic", "eusocial", "unknown"),
       col    = eusocial_cols[c("non-eusocial", "polymorphic", "eusocial", "unknown")],
       pch = 16, pt.cex = 0.9,
       title = expression(bold("Eusociality state")),
       bty = "n", cex = 0.7,
       x.intersp = 1.3)

dev.off()









# Ancestral State Reconstruction of Eusociality in Bees*

## Ancestral state reconstruction ('01_eusociality_asr.R')
This repository contains code to reconstruct the evolutionary history of eusociality across bees using corHMM. 

First, we fit models of trait evolution, compare them using AIC to identify the best-supported model, and then map ancestral sociality states across a phylogeny of ~4,586 species. Next, we extract transition rate estimates (i.e. number of transitions to/from eusociality per million year), and visualize the distribution of eusociality across the tree.

Two alternative runs are done to assess how uncertainty in trait scoring for one particularly variable and poorly studied clade, _Lasioglossum_, impacts the results:
- **Run 1:** Unknown *Lasioglossum* states manually scored as primitively eusocial  
- **Run 2:** Unknown *Lasioglossum* left as missing and inferred during model fitting

## Stochastic character mapping ('02_corHMM_simmaps.R')

This script uses the best-fitting hidden-rate `corHMM` model to simulate the evolutionary history of eusociality across the bee phylogeny using stochastic character mapping (SIMMAP). That is, this step simulates many possible histories across the tree given our estimated transition rates and model for trait evolution, and then summarizes how often transitions occur on average across these simulated evolutionary histories in order to understand how often and when gains and losses of eusociality occurred across the bee phylogeny.

## Data Inputs

Located in `raw_data/`:

- `ML_beetree_4586tips.tre`  
  Time-scaled bee phylogeny with 4,586 species.

- `sociality_scored_lasioglossum_inferred.csv`  
  Bee species from the phylogeny scored for sociality. Here, unknown _Lasioglossum_ states are manually inferred as primitively eusocial (Run 1).

- `sociality_scored_lasioglossum_unknown.csv`  
  Bee species from the phylogeny scored for sociality. Here, unknown _Lasioglossum_ states remain unknowns to be inferred by the model (Run 2).


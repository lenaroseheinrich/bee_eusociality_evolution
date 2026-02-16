# Evolution of eusociality in bees
This repository contains code to reconstruct the evolutionary history of eusociality across bees using corHMM. 

## Ancestral state reconstruction ('01_eusociality_asr.R')
In this script, we fit models of trait evolution, compare them using AIC to identify the best-supported model, and extract transition rate estimates (i.e. number of transitions to/from eusociality per million year) from the best model. In this case, the best model was a model with two rate-classes, meaning the rate at which eusociality evolved across the tree was variable, with some fast and some slow clades. Next, we reconstruct ancestral eusociality states across the phylogeny and plot the ancestral state reconstruction across the tree.

Two alternative runs were done to assess how uncertainty in trait scoring for one particularly variable clade, _Lasioglossum_, impacts the results:
- **Run 1:** Unknown *Lasioglossum* states manually scored as (primitively) eusocial
- **Run 2:** Unknown *Lasioglossum* left as unknowns and inferred during model fitting

## Stochastic character mapping ('02_corHMM_simmaps.R')
This script uses the best-fitting model identified above and its estimated transition rates to perform many simulations of possible evolutionary histories of eusociality across the tree. It then summarizes how many transitions occur on average across these simulated evolutionary histories in order to understand how often and when gains and losses of eusociality are likely to have occurred while still reflecting the uncertainty inherent to these inferences.

## Data inputs

Located in `raw_data/`:

- `ML_beetree_4586tips.tre`  
  Time-scaled bee phylogeny with 4,586 species.

- `sociality_scored_lasioglossum_inferred.csv`  
  Bee species from the phylogeny scored for sociality. Here, unknown _Lasioglossum_ states are manually inferred as primitively eusocial (Run 1).

- `sociality_scored_lasioglossum_unknown.csv`  
  Bee species from the phylogeny scored for sociality. Here, unknown _Lasioglossum_ states remain unknowns to be inferred by the model (Run 2).

## Outputs
Located in `results/`:

### Model fitting and comparison
- `corhmm_tbl_dredge_run1.csv` / `corhmm_tbl_dredge_run2.csv`  
  Model comparison tables summarizing all tested models, their AICc values, and relative support.

- `transition_rates_corHMM_run1.csv` / `transition_rates_corHMM_run2.csv`  
  Estimated transition rate matrices from the best-supported model (rate = transitions among eusociality states per million years).

### Ancestral state reconstructions
- `corhmm_dredge_recon_final_run1.pdf` / `corhmm_dredge_recon_final_run2.pdf`  
  Phylogenies with ancestral state likelihoods at nodes.
  
- `eusociality_ASR_circular_phylogeny_run1.pdf` / `eusociality_ASR_circular_phylogeny_run2.pdf`  
  Circular phylogenies with eusociality states mapped along the tips and branches colored based on the ancestral state reconstruction. 

### Stochastic character mapping (SIMMAP)
- `transitions_summary_corhmm_run1.csv` / `transitions_summary_corhmm_run2.csv`  
  Summary tables reporting the mean and variance in the number of transitions among eusociality states across simulated evolutionary histories.

- `transitions_summary_plot_run1.pdf` / `transitions_summary_plot_run2.pdf`  
  Figures visualizing the distribution of transition counts (i.e., how many gains and losses of eusociality are estimated to have occurred and how long ago they happened).



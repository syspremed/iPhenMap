# iPhenMap from the Sadanandam Lab at the ICR, London
Integrated PhenMap

This tool integrates multiomics with phenomics data

Integrative Phenotype Mapping (iPhenMap) maximizes information shared between and within datasets and attempts to explain this variation using additional patient phenotypic information 

#install iPhenMap

install.packages("devtools")

devtools::install_github("syspremed/iPhenMap)

Run iPhenMAP Analysis
Description
This function reads clinical, CNA, proteomics, metabolomics, and RNA-seq data, preprocesses and matches samples, and then runs the iPhenMAP model.

Usage

iphen(
  D,
  X,
  q = 10,
  mcmc = 10000,
  g = 1000,
  c = 1,
  d = 1,
  nu0 = 0.01,
  Sigo = 0.01,
  s2o = 0.01,
  nu0k = 0.01,
  noiseMod = "unique",
  scaleType = "unit",
  iniType = "EVD",
  time = NULL,
  status = NULL
)


Arguments
q	
Number of latent variables (default 10)

mcmc	
Number of MCMC iterations (default 100000)

data_dir	
Path to the parent data directory containing subfolders: Clinical, CNA, Proteomics, Metabolomics, RNAseq.

...	
Additional parameters passed to iphen (g, c, d, nu0, Sigo, s2o, nu0k, noiseMod, scaleType, iniType, time, status)

Value
Saves results in "iPhenMAPresults" and returns the final MCMC output list.

## Roadmap: deep learning enhancements

iPhenMap cohorts typically have tens to hundreds of samples with thousands of features. At that size, deep models overfit easily and give up the calibrated uncertainty and interpretable loadings of the Bayesian model. The ideas below therefore keep iPhenMap's structure (shared factors, factors switched on or off per view, phenotypes regressed on factors) and use neural networks only where they add something.

1. **Amortized variational inference.** Fit the same model with an encoder network trained by stochastic variational inference (for example Pyro or NumPyro) instead of the Gibbs sampler. This scales to large feature counts, runs in minutes on a GPU, and gives `predict()` for new samples for free. Gibbs stays as the reference, because variational inference tends to underestimate uncertainty.
2. **Multi-view VAE that handles missing views.** One encoder per view, combined by a product-of-experts posterior, so samples missing a whole view still get factors and the missing view can be imputed. A linear decoder keeps loadings interpretable, and per-view group sparsity on the decoder weights plays the role of the `z` activation matrix. Closest prior art: MOVE (Allesøe et al., *Nature Biotechnology*, 2023).
3. **Survival modelled jointly with the factors.** A Cox partial-likelihood head (DeepSurv-style) on the latent factors, trained together with the rest of the model. This replaces the post-hoc Cox fit on posterior medians and should be checked with a cross-validated C-index.
4. **Pathway-structured networks.** Sparse layers that follow Reactome (gene → pathway → higher-level pathway), as in P-NET (Elmarakeby et al., *Nature*, 2021), so pathways are part of the model instead of a post-hoc enrichment step.
5. **Clustering inside the model.** A Gaussian-mixture prior on the latent space (VaDE), so subtypes are learned together with the factors instead of by separate ConsensusClusterPlus and mclust steps. Most useful for larger cohorts.

Every candidate should be benchmarked against the current Gibbs sampler on:
- simulated data with a linear structure, where it must at least match on factor recovery, view activity and phenotype effects;
- simulated data with a nonlinear structure, where a deep model should win;
- simulated data with missing views;
- a public cohort such as TCGA, using cross-validated survival C-index and comparing against MOFA2 and iClusterBayes.

Suggested first step: item 1 with the linear decoder and group sparsity from item 2.

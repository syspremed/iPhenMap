# iPhenMap
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

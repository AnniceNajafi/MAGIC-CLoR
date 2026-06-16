<h1>MAGIC-CLoR</h1>

<i>Null-calibrated association scoring improves transcription-factor target recovery after single-cell data diffusion.</i>

This repository holds the R package (<code>cvMAGIC</code>) and the reproduction notebooks for the MAGIC-CLoR study.

MAGIC imputes single-cell RNA-seq data by diffusing expression over a cell-cell similarity graph, after which gene-gene relationships are usually read out from correlations of the imputed matrix. The study asks a single question: is it the diffusion (how the data are smoothed) or the readout (how a relationship is scored from the smoothed data) that limits transcription-factor (TF) to target recovery? Across four UMI datasets we find that tuning the smoothing (diffusion-time selection, per-gene/per-module adaptive t, a learned multi-scale mixture, and a better graph) does not give a consistent improvement, while replacing the raw-correlation readout with a Context Likelihood of Relatedness (CLR) calibration does. MAGIC-CLoR keeps the diffusion as it is and changes only the readout.

<h2>🖥️ Installation</h2>

Install the package directly from this repository:

```
install.packages("devtools")
devtools::install_github("AnniceNajafi/MAGIC-CLoR", subdir = "cvMAGIC")
```

or from a local clone:

```
install.packages("cvMAGIC", repos = NULL, type = "source")
```

The package depends on Matrix, RANN, irlba, and stats. Rmagic (and its Python magic-impute backend) is only needed for the notebook cells that compare against the upstream auto-selector.

<h2>The method in one paragraph</h2>

After diffusion, two genes can correlate simply because they ride the same smoothed cell-state manifold (lineage, library size, module activity), so raw correlation mixes direct association with manifold-wide co-expression and hub genes correlate with almost everything. The CLR transform (Faith et al., PLoS Biology 2007) null-calibrates each pair against the background distribution of each gene's correlations, so a pair scores highly only when the two genes are more correlated with each other than each is typically correlated with all genes. CLR has no learned parameters and uses no relationship labels, so it cannot overfit the benchmark.

<h2>Usage</h2>

The exported functions are:

| Function | What it does |
|---|---|
| <code>magic_graph()</code> | build the cell-cell graph (PCA, kNN, adaptive alpha-decay kernel, Markov matrix) |
| <code>magic_impute()</code>, <code>magic()</code> | diffuse the data for t steps (t = 0 returns the data unchanged) |
| <code>knee_select_t()</code> | the original knee-point heuristic for choosing t |
| <code>mcv_select_t()</code> | molecular cross-validation choice of t |
| <code>mcv_select_t_tolerant()</code>, <code>magic_tolerant()</code> | tolerance-regularised MCV |
| <code>magic_per_gene()</code> | per-gene adaptive t |
| <code>magic_modules()</code> | per-module adaptive t |
| <code>magic_bootstrap_se()</code> | bootstrap per-(cell, gene) standard errors |
| <code>clr_from_cor()</code> | CLR calibration of any correlation matrix |
| <code>magic_clr()</code> | run MAGIC imputation and return the CLR readout |

A minimal example:

```
library(cvMAGIC)

# X is a cells x genes matrix of raw integer UMI counts
res <- magic_clr(X, t = 0)   # t = 0 scores CLR on raw counts; t > 0 diffuses first
clr <- res$assoc             # gene-gene CLR association matrix
```

<h2>Reproducing the paper</h2>

Every notebook runs the real cvMAGIC pipeline (and Rmagic for the upstream comparator); there are no precomputed or hard-coded figures. Each notebook reproduces one part of the paper and nothing beyond it.

| Paper element | Notebook |
|---|---|
| Four-dataset benchmark, cross-dataset summary, and the prediction-vs-recovery divergence | <code>paper_benchmark.Rmd</code> |
| Per-gene Pearson MCV is dataset-specific | <code>per_gene_comparison.Rmd</code>, <code>pbmc_realdata.Rmd</code> |
| Tolerance regularisation does not generalise (sparsity sweep and Paul seeds) | <code>tolerance_stress_test.Rmd</code> |
| Tolerance is inert in its ideal regime | <code>tusi_trajectory.Rmd</code> |
| The learned multi-scale mixture (msMAGIC) overfits | <code>msmagic_prototype.Rmd</code> |
| A better graph (AnchorMAGIC) does not raise the oracle | <code>anchor_magic.Rmd</code> |
| Bootstrap uncertainty calibration | <code>bootstrap_calibration.Rmd</code> |
| MAGIC-CLoR positive result (random negatives) | <code>residual_magic.Rmd</code> |
| MAGIC-CLoR matched-negative control | <code>residual_magic_matched.Rmd</code> |
| MAGIC-CLoR control checks at usable t, resampling, and detection-stratified AUC | <code>clr_controls.Rmd</code> |
| External validation on the DoRothEA (A/B) regulons | <code>dorothea_validation.Rmd</code> |
| Single-dataset fidelity checks | <code>paul_realdata.Rmd</code>, <code>tabula_muris.Rmd</code> |

Run order: install cvMAGIC, then from the <code>notebooks/</code> folder render <code>paper_benchmark.Rmd</code> first (it processes all four datasets identically), then the rest. Figures are written to <code>../figures/</code>.

```
install.packages("cvMAGIC", repos = NULL, type = "source")
setwd("notebooks")
rmarkdown::render("paper_benchmark.Rmd")
```

<h2>Data</h2>

The notebooks read raw UMI count matrices from <code>../data/</code>. None of the data is redistributed here; see <code>data/README.md</code> for the accessions and the expected file layout.

| Dataset | Source |
|---|---|
| Paul 2015 (MARS-seq) | GEO GSE72857 |
| PBMC 3k (10x v1) | 10x Genomics public dataset |
| Tabula Muris marrow (10x droplet) | figshare droplet v2 / GEO GSE109774 |
| Tusi 2018 (inDrop) | GEO GSE89754 |

MCV relies on the Poisson nature of UMI counts, so Smart-seq2 and read-count datasets are deliberately excluded.

<h3>Open source disclaimer</h3>

This software is released under the MIT License and is provided for research use without warranty of any kind. If you use it, please cite the MAGIC-CLoR paper and the original MAGIC method (van Dijk et al., Cell 2018).

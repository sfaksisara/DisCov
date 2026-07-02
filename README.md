# DiscCov: Discriminability-Coverage Adaptive Gene Selection

Implementation and full evaluation pipeline for **DiscCov**, an adaptive
filter-based gene selection framework for cancer classification on
high-dimensional microarray data. DiscCov combines a rank-aggregated
discriminability score (SNR, *t*-statistic, ANOVA *F*-statistic) with a
class-conditioned coverage score, balanced by a dataset-size-aware adaptive
weighting mechanism, and uses a dual-index strategy (greedy order for small
subsets, SVM-reranked order for larger subsets).

This repository contains the code used to produce the results in the
accompanying paper, *"DiscCov: A Discriminability-Covariance Adaptive Gene
Selection Framework for Cancer Classification on Microarray Data."*

## Overview

The pipeline:

1. Preprocesses raw microarray expression data (imputation, log₂ transform,
   quantile normalisation, variance pre-filtering).
2. Selects genes using **12 feature-selection methods**:
   - **DiscCov family:** DiscCov (adaptive), AL-fixed (α=β=0.5), disc-only
     (β=1), cov-only (α=1)
   - **Classical baselines:** mRMR, *t*-test, LASSO, ReliefF, RF-importance,
     SVM-RFE, Elastic Net, ANOVA F-test
3. Evaluates each selected gene subset across **7 subset sizes**
   (k = 5, 10, 25, 50, 100, 150, 200) and **9 classifiers** (Logistic
   Regression, Random Forest, linear SVM, LDA, KNN, Gradient Boosting,
   AdaBoost, Voting Ensemble, Stacking Ensemble).
4. Reports **6 metrics** (accuracy, balanced accuracy, macro-F1, MCC,
   Cohen's κ, ROC-AUC) under stratified k-fold cross-validation.
5. Runs statistical significance testing (one-sided Wilcoxon signed-rank,
   Cohen's *d*, Friedman test) and generates 10 summary plots.

## Repository Structure

```
.
├── disccov_pipeline.py     # Main script (this file)
├── prostate.csv             # Input dataset (genes as columns, last column = label)
├── al_gene_outputs_v4/      # Auto-created output directory
│   ├── results_per_classifier.csv
│   ├── results_averaged.csv
│   ├── statistical_tests.csv
│   └── plot01_...10_*.png
└── README.md
```

> If `prostate.csv` is not found, the script generates a synthetic dataset
> automatically so the pipeline can be run end-to-end without real data.

## Requirements

- Python ≥ 3.9
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- seaborn

Install everything with:

```bash
pip install numpy pandas scipy scikit-learn matplotlib seaborn
```

## Usage

Place your gene expression dataset as `prostate.csv` in the working
directory, with genes/probes as columns and the class label as the final
column, then run:

```bash
python disccov_pipeline.py
```

The script will:
- Preprocess the data and cap dimensionality at 2,000 genes via a variance
  pre-filter if needed.
- Run 10-fold stratified cross-validation (5-fold for very small datasets —
  edit `n_folds` in `run_pipeline(...)` inside `main()` as needed).
- Print per-method / per-classifier summary tables and Wilcoxon /
  Friedman statistical test results to the console.
- Save CSV result tables and 10 PNG plots to `al_gene_outputs_v4/`
  (or to Google Drive automatically if run inside Google Colab).

### Using DiscCov as a standalone selector

```python
from disccov_pipeline import ALGeneSelector

selector = ALGeneSelector(n_genes=100, random_state=42)
selector.fit(X_train, y_train)

X_train_k = selector.transform(X_train, k=50)
X_test_k  = selector.transform(X_test,  k=50)
```

- `alpha_fixed` / `beta_fixed`: set both to fix the coverage/discriminability
  weights (e.g. `alpha_fixed=0.5, beta_fixed=0.5` reproduces the AL-fixed
  ablation). Leave both `None` for the fully adaptive DiscCov behaviour.
- `k <= 25` uses the frozen greedy order (`al_indices_`); `k > 25` uses the
  SVM-reranked order (`selected_indices_`).

## Key Parameters

| Constant     | Value | Meaning                                          |
|--------------|-------|---------------------------------------------------|
| `BASE_BETA`  | 0.65  | Base discriminability coefficient (β₀)             |
| `BETA_MAX`   | 0.80  | Late-stage discriminability floor (after k > 25)   |
| `ALPHA_MIN`  | 0.20  | Minimum coverage weight                            |
| `K_VALUES`   | [5, 10, 25, 50, 100, 150, 200] | Evaluated gene subset sizes |

## Outputs

| File                                 | Description                                   |
|---------------------------------------|-----------------------------------------------|
| `results_per_classifier.csv`          | Metrics per method × k × classifier            |
| `results_averaged.csv`                | Metrics per method × k, averaged over classifiers |
| `statistical_tests.csv`               | Wilcoxon p-values, Cohen's d vs DiscCov        |
| `plot01_accuracy_f1_vs_k.png`         | Accuracy & F1-macro vs. gene subset size       |
| `plot02_auc_vs_k.png`                 | ROC-AUC vs. gene subset size                   |
| `plot03–05_heatmap_*.png`             | Method × classifier heatmaps (Acc / F1 / AUC)  |
| `plot06_best_metrics_bar.png`         | Best score per method across all k             |
| `plot07_boxplot_three_metrics.png`    | Per-fold distributions at best k               |
| `plot08_clf_comparison_three_metrics.png` | DiscCov performance per classifier         |
| `plot09_effect_size.png`              | Cohen's d effect sizes vs. DiscCov             |
| `plot10_multi_metric_grouped.png`     | Grouped bar chart of all methods at best k     |

## Reproducibility

All stochastic components (cross-validation splits, classifiers, LinearSVC
re-ranking) are seeded with `random_state=42` by default for full
reproducibility of the reported results.



Sara Sfaksi — sara.sfaksi@univ-biskra.dz
Leila Djerou — l.djerou@univ-biskra.dz
Laboratory of LESIA, University of Biskra, Algeria

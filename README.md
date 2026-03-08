# Somatic Variant Detection Hyperparameter Optimization

Bayesian optimization of filtering thresholds for matching single-cell somatic variants to bulk tumor variants in disseminated tumor cells (DTCs).

## Technologies
- Python, Optuna, cyvcf2, pysam, pandas, NumPy, anndata

## Approach
Uses Optuna's TPE sampler to jointly optimize 7 filtering parameters across cell-level (read depth, alt reads, allele frequency) and tumor-level (blood contamination, tumor support, population frequency) thresholds. F1 score is computed by treating same-patient variant matches as true positives and cross-patient matches as false positives, after excluding germline variants.

## Results
- **Best F1 Score:** 0.843
- **Precision:** 0.979 | **Recall:** 0.740
- Optimized over 2,604 single-cell VCFs across 9 patients (100 trials)

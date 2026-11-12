# qSIP calculations

This folder has the code and input data used to calculate the ¹⁸O excess atom fraction (EAF) of each ASV with quantitative stable isotope probing (qSIP). 

---

## Folder structure

```
qSIP calculations/
├── Data/                                   # input data
│   ├── 18O.lab.data.density.corrected_subset.csv
│   ├── DRIP16S_table_2x2_subset.txt
│   ├── DRIP16S_2x2_RDP50.txt
│   ├── DRIPITS_table_2x2_subset.txt
│   └── DRIPITS_table_2x2_UNITE.csv
├── 16S qSIP calculations/
│   ├── calcs/                              # EAF 18O calculations, one script per treatment
│   │   ├── det normal.R
│   │   ├── rhizo det normal.R
│   │   └── rhizo normal.R
│   ├── remove outliers/
│   │   └── remove outliers.R
│   └── results/                            # original outputs reported in the paper
│       ├── detritus_normal.csv
│       ├── rhizo_det_normal.csv
│       └── rhizo_normal.csv
└── ITS/
    ├── calcs/                              # same three scripts, for ITS
    ├── remove outliers/
    └── results/
        ├── detritus_normal.csv
        ├── rhizo_detritus_normal.csv       # equivalent of 16S rhizo_det_normal.csv
        └── rhizo_normal.csv
```

---

## Software

| Software | Version used for the results in `results/` | Version tested for reproduction |
|---|---|---|
| R | 4.4.237 | 4.4.3 |
| [qSIP2](https://github.com/jeffkimbrel/qSIP2) | version installed approx. Oct 2024 | 0.23.8 (commit `0026bab`) |


### Differences between qSIP2 versions (please read)

The files in `*/results/` are the **original outputs reported in the paper**. They were made with an older version of qSIP2. Running the scripts with qSIP2 0.23.8 gives slightly different output:

1. **Tube IDs must be text.** In newer qSIP2, `unlabeled_source_mat_ids` and `labeled_source_mat_ids` must be character vectors. Wrap the IDs in `as.character()`.
2. **Extra output columns.** `summarize_EAF_values()` now also returns `pval`, `labeled_resamples`, `unlabeled_resamples`, `labeled_sources` and `unlabeled_sources`.
3. **`observed_EAF` is the same** in both versions for every ASV.

qSIP2 settings used in all comparisons:

- feature filter: `min_unlabeled_sources = 2`, `min_labeled_sources = 1`, `min_unlabeled_fractions = 2`, `min_labeled_fractions = 2`
- resampling: `resamples = 1000`, `with_seed = 17`, `allow_failures = TRUE`
- confidence intervals: `confidence = 0.9`

---

# Skin Tumor Analysis

R script for analysis and visualization of longitudinal skin tumor measurements in mice across six genotypes differing in chromosomal instability (CIN) levels.

## Contents

```
├── data/
│   └── tumor_sizes_raw.csv.xlsx   # Input data (not tracked by git)
├── output/                         # Generated plots and CSVs (not tracked by git)
├── tumor_analysis_cleaned.R        # Main analysis script
└── README.md
```

## Input format

The input is an Excel file (`.xlsx`) with the following structure:

- **Column 1:** Week label (e.g. "Week 1")
- **Column 2:** Timepoint in days (numeric)
- **Remaining columns:** Individual tumor sizes in mm³, one column per tumor. Column headers encode mouse ID and genotype (e.g. `7080001_WT/WT`).

## How to run

1. Clone the repository
2. Place the input file at `data/tumor_sizes_raw.csv.xlsx`
3. Open `skin_tumor_analysis.Rproj` in RStudio — this sets the working directory automatically
4. Run `tumor_analysis_cleaned.R`

All plots and CSV outputs are saved to the `output/` folder.

## Dependencies

R packages (version numbers reflect those used during development — update after running `sessionInfo()`):

| Package    | Version  |
|------------|----------|
| stringr    | 1.5.0    |
| ggplot2    | 3.4.4    |
| ggbeeswarm | 0.7.2    |
| openxlsx   | 4.2.5.2  |
| readr      | 2.1.4    |
| tidyr      | 1.3.0    |
| dplyr      | 1.1.3    |
| emmeans    | 1.8.9    |
| ggh4x      | 0.2.6    |
| ggpubr     | 0.6.0    |
| rstatix    | 0.7.2    |
| ggsignif   | 0.6.4    |
| MASS       | 7.3-60   |
| dunn.test  | 1.3.5    |
| nlme       | 3.1-163  |

## Outputs

### Plots
| File | Description |
|------|-------------|
| `tumor_onset_per_genotype.png` | Tumor onset day per genotype |
| `tumor_size_crosssection_day145.png` | Mean tumor size per mouse at day 145 |
| `tumor_growth_curves_per_tumor_by_genotype_days_from_onset_upper.png` | Growth curves per tumor by genotype |
| `tumor_growth_curves_by_genotype_and_max_sizeclass_days_from_onset_lower.png` | Growth curves by genotype and max size class |
| `average_tumor_counts_per_day_genotype_sizeclass_barplot_lower.png` | Average tumor counts over absolute time |
| `average_tumor_counts_by_genotype_sizeclass_weeks_from_mouse_onset_lower.png` | Average tumor counts from mouse onset |
| `tumor_growth_rates_boxplot.png` | Exponential growth rates per genotype |
| `incidence_rate_overall_stats.png` | Incidence rate from experiment start |
| `incidence_rate_from_onset_stats.png` | Incidence rate from first tumor onset |

### CSVs
Key intermediate data files saved for downstream use or verification, including long-format tumor data, growth rates, incidence rates, and statistical test results.

## Statistical methods

Tumor growth rates were estimated by fitting log-linear models to each individual tumor (minimum 3 timepoints). Differences between genotypes were assessed using Kruskal-Wallis tests followed by Dunn post-hoc pairwise comparisons against WT/WT with Benjamini-Hochberg correction. A linear mixed model with mouse as a random intercept was additionally fitted to growth rates to account for multiple tumors per mouse. Tumor incidence rates were calculated as tumors per mouse per week, either from experiment start or from first tumor onset.

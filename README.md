# MASLD_project


1. Overview

This pipeline automates:

Data loading and preprocessing

Statistical testing with a volcano plot

Group-wise correlation analysis and Fisher z-based correlation-difference testing

Machine learning with ROC evaluation and overfitting checks

External validation of key metabolite correlations

Session information capture for reproducibility

[Data Load] -> [Z-score Normalize] -> [Stats + Volcano]
                -> [Group-wise Correlation + Fisher Z-test]
                -> [ML Models + ROC + Overfitting]
                -> [External Validation] -> [Reports]

2. Project Structure and Required Files
project_root/
├─ data/
│  ├─ 250721_MASLD.xlsx
│  ├─ 250726_val_01_plasma.xlsx
│  └─ 250726_val_03_plasma.xlsx
├─ outputs/               # auto-created
│  ├─ significant_plots/  # auto-created
│  └─ ...                 # result files are written here
└─ script.R               # the provided script

2.1 250721_MASLD.xlsx sheet layout

Sheet df

Column 1: Sample

Column 2: Group with values NC or MASLD

Columns 3+: metabolite features

Sheet mapping

Feature ID to human-readable name mapping

Notes

The script uses normalized values for statistical tests and raw values for fold change calculations.

2.2 External validation files

250726_val_01_plasma.xlsx sheet test_all

250726_val_03_plasma.xlsx sheet test_all

Both are standardized to two variables: M_180 and M_263.

3. Environment Setup
3.1 Required R packages

The script auto-installs and loads packages listed in required_packages.
xlsx depends on Java. If Java is problematic, you can switch to readxl as below.

Alternative without Java

Replace the xlsx::read.xlsx calls with readxl::read_excel:

# Original
raw_df     <- xlsx::read.xlsx("data/250721_MASLD.xlsx", sheetName = "df")
mapping_df <- xlsx::read.xlsx("data/250721_MASLD.xlsx", sheetName = "mapping")

# Alternative
raw_df     <- readxl::read_excel("data/250721_MASLD.xlsx", sheet = "df")
mapping_df <- readxl::read_excel("data/250721_MASLD.xlsx", sheet = "mapping")

3.2 Working directory

Set your project root:

setwd("~/projects/masld-pipeline")

4. Quick Start
4.1 Run in R console or RStudio
source("script.R")
results <- main_analysis_pipeline()
print_session_info()

4.2 Run from command line
Rscript script.R


If the script prints “Script loaded. Run main_analysis_pipeline() to execute the complete analysis.”, change the ending block to run automatically:

results <- main_analysis_pipeline()
print_session_info()

5. What Each Step Does
5.1 Data loading and preprocessing

Function: load_and_preprocess_data()

Z-score normalization is applied to metabolite columns (columns 3+).

Prints sample count, feature count, and group levels.

5.2 Statistical analysis and volcano plot

Function: perform_statistical_analysis(df_normalized, df_raw)

Performs LMSstat::All_stats and generates boxplots.

Manual t tests produce log2FC, pvalue, FDR, and volcano plot.

Output

outputs/volcano_plot.svg

Manual results are returned in the manual_results object.

5.3 Correlation analysis and between-group differences

Function: perform_correlation_analysis(df_normalized)

Pearson correlations via Hmisc::rcorr for NC and MASLD separately.

MASLD heatmap uses NC clustering order for easier visual comparison.

Fisher z tests evaluate correlation differences between groups.

Coefficients of variation (CV) per group are computed and merged with fold change and univariate p values.

Outputs

outputs/correlation_heatmap_NC.svg

outputs/correlation_heatmap_MASLD.svg

outputs/significant_correlation_pairs.xlsx

outputs/significant_plots/ top 10 scatterplots

5.4 Machine learning and ROC analysis

Function: perform_ml_analysis(df_normalized)

Feature sets

mean_diff_features as provided examples

interaction_features are placeholders and should be replaced with actual significant metabolites

Derived features

Local correlation

Group-centered interaction term

Models

Logistic regression and Random Forest with 5-fold CV using twoClassSummary and ROC metric

Evaluation

Train and test AUC, sensitivity, specificity, F1

ROC curves saved in base R and ggplot formats

Overfitting summary table

Outputs

outputs/ROC_curves_test.svg

outputs/ROC_curves_ggplot.svg

outputs/model_performance_comparison.xlsx

5.5 External validation

Function: perform_validation_analysis()

Pearson correlations between M_180 and M_263 in two validation sets plus scatterplots.

Outputs

outputs/validation_set_01_scatter.svg

outputs/validation_set_03_scatter.svg

outputs/validation_comparison.svg

5.6 Reproducibility information

Function: print_session_info()

Output file: outputs/session_info.txt

6. Function Summary

load_and_preprocess_data()
Load inputs and normalize metabolite columns.

perform_statistical_analysis(df, raw_df)
LMSstat summary, manual t tests, volcano plot.

perform_correlation_analysis(df)
Group-wise correlations, Fisher z difference tests, heatmaps, and scatterplots.

perform_ml_analysis(df)
Derived features, model training, ROC evaluation, and overfitting checks.

perform_validation_analysis()
External correlation validation and comparison plots.

main_analysis_pipeline()
Orchestration that runs the full workflow and returns a results list.

7. Outputs Summary

Statistics

outputs/volcano_plot.svg

Correlations

outputs/correlation_heatmap_NC.svg

outputs/correlation_heatmap_MASLD.svg

outputs/significant_correlation_pairs.xlsx

outputs/significant_plots/*.svg

Machine learning

outputs/model_performance_comparison.xlsx

outputs/ROC_curves_test.svg

outputs/ROC_curves_ggplot.svg

External validation

outputs/validation_set_01_scatter.svg

outputs/validation_set_03_scatter.svg

outputs/validation_comparison.svg

Reproducibility

outputs/session_info.txt

8. Customization Guide
8.1 Replace feature sets
# Mean difference based features
mean_diff_features <- c("M_107", "M_109", "M_44", "M_185", "M_211", "M_111")

# Interaction candidate features (2 variables expected by the current logic)
interaction_features <- c("M_000", "M_999")  # Replace with significant metabolite names


Update based on volcano results or correlation-difference findings.

Keep interaction_features length at 2 for the current derived-feature logic.

8.2 Significance thresholds

Adjust FDR cutoff and multiple-testing method inside the relevant functions.

Example: change p.adjust(..., method = "BH") to "BY" if needed.

8.3 Cross-validation and tuning

Modify trainControl parameters number and repeats.

Tune Random Forest ntree and tuneLength.

8.4 Output formats

Change ggsave(..., device = "svg") to png or pdf and set width, height, dpi to match journal requirements.

9. Statistical Notes

t tests are performed on normalized features with FDR correction.

Pearson correlations and Fisher z tests compare correlation structures across groups.

If a group has fewer than 4 samples, Fisher z variance is undefined and such pairs are skipped or set to NA.

CV is unstable when the mean approaches zero and is set to NA in that case.

10. Troubleshooting

xlsx installation errors

Install Java or switch to readxl as shown above.

Empty heatmaps

Too many missing values or zero-variance features. Consider imputation or feature filtering.

ROC errors

Ensure binary factor order is correct:

df$Group <- factor(df$Group, levels = c("NC", "MASLD"))


Font issues in plots

SVG is usually robust. If exporting PNG, check system fonts and consider ggtext alternatives if needed.

11. Reproducibility Tips

set.seed(42) is used in several places. For complete consistency, set it once at the top as well.

Use renv to lock package versions.

install.packages("renv")
renv::init()

12. Citations and Licensing

To cite R packages, use:

citation("pROC")
citation("caret")
citation("Hmisc")


Follow your lab or grant policy for data and code licensing.

13. Changelog Template

2025-08-21

Initial README in English

Documented correlation-difference tests and ML outputs

14. Contact

Maintainer: [Your Name]

Email: [your.email@domain]

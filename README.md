# MASLD Metabolomics Analysis Pipeline

## Overview

This R pipeline provides a comprehensive analysis framework for metabolomics data comparing Normal Control (NC) vs MASLD (Metabolic Dysfunction-Associated Steatotic Liver Disease) groups. The pipeline includes statistical testing, correlation analysis, machine learning models, and external validation.

## Features

- **Statistical Analysis**: Differential metabolite analysis with FDR correction
- **Visualization**: Volcano plots, correlation heatmaps, and scatter plots
- **Machine Learning**: Multiple classification models with ROC analysis
- **Correlation Analysis**: Group-specific correlation patterns and differences
- **External Validation**: Independent dataset validation
- **Reproducibility**: Comprehensive session information and seed management

## Requirements

### R Version
- R >= 4.0.0

### Required Packages
```r
required_packages <- c(
  "xlsx", "dplyr", "ggplot2", "tidyverse", "pheatmap", "Hmisc", 
  "svglite", "psych", "LMSstat", "caret", "pROC", "ROCR", "readxl"
)
```

## Installation

1. Clone or download the analysis script
2. Install required packages:
```r
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
```

## Data Structure

### Input Files Required
```
data/
├── 250721_MASLD.xlsx
│   ├── df (sheet): Main metabolomics data
│   └── mapping (sheet): Sample mapping information
├── 250726_val_01_plasma.xlsx
│   └── test_all (sheet): Validation dataset 1
└── 250726_val_03_plasma.xlsx
    └── test_all (sheet): Validation dataset 2
```

### Expected Data Format
- **Sample**: Sample identifiers
- **Group**: NC or MASLD classification
- **Metabolites**: M_001, M_002, ... M_xxx columns with numerical values

## Pipeline Workflow

### Step 1: Data Import and Preprocessing
```r
# Load and normalize data
data_list <- load_and_preprocess_data()
```
- Imports raw metabolomics data
- Applies z-score normalization to metabolite columns
- Creates normalized dataset for analysis

### Step 2: Statistical Analysis
```r
# Perform differential analysis
stat_results <- perform_statistical_analysis(data_list$normalized, data_list$raw)
```
- **LMSstat Analysis**: Comprehensive statistical testing with FDR correction
- **Manual t-tests**: Individual metabolite comparisons
- **Volcano Plot**: Visualization of fold changes vs significance
- **Boxplots**: Group comparison visualization

**Key Functions:**
- `All_stats()`: LMSstat comprehensive analysis
- `perform_manual_ttest()`: Individual t-test calculations
- `create_volcano_plot()`: Volcano plot generation

### Step 3: Correlation Analysis
```r
# Analyze correlation patterns
correlation_results <- perform_correlation_analysis(data_list$normalized)
```
- **Group-specific Correlations**: Separate correlation matrices for NC and MASLD
- **Correlation Heatmaps**: Clustered visualization of metabolite relationships
- **Differential Correlations**: Fisher's Z-test for correlation differences
- **Significant Pairs**: Top correlation differences between groups

**Key Functions:**
- `calculate_group_correlations()`: Pearson correlations with p-values
- `create_correlation_heatmaps()`: Clustered heatmap visualization
- `find_significant_correlation_differences()`: Fisher's Z-test analysis

### Step 4: Machine Learning Analysis
```r
# Build and evaluate classification models
ml_results <- perform_ml_analysis(data_list$normalized)
```
- **Feature Engineering**: Local correlations and group-centered interactions
- **Model Training**: GLM and Random Forest approaches
- **Cross-validation**: 5-fold CV with performance metrics
- **ROC Analysis**: Comprehensive model evaluation
- **Overfitting Detection**: Train vs test performance comparison

**Model Types:**
- **Mean Difference GLM**: Based on significantly different metabolites
- **Interaction GLM**: Using derived correlation features
- **Random Forest**: Ensemble methods for both feature sets

### Step 5: Validation Analysis
```r
# External dataset validation
validation_results <- perform_validation_analysis()
```
- **External Datasets**: Independent validation using plasma samples
- **Correlation Validation**: M_263 vs M_180 relationships
- **Comparative Analysis**: Cross-dataset correlation patterns

## Output Structure

```
outputs/
├── volcano_plot.svg                    # Differential analysis
├── correlation_heatmap_NC.svg          # NC group correlations
├── correlation_heatmap_MASLD.svg       # MASLD group correlations
├── ROC_curves_test.svg                 # Model performance
├── ROC_curves_ggplot.svg              # Enhanced ROC visualization
├── validation_comparison.svg           # External validation
├── significant_correlation_pairs.xlsx  # Statistical results
├── model_performance_comparison.xlsx   # ML performance metrics
├── session_info.txt                   # Reproducibility info
└── significant_plots/                 # Individual scatter plots
    ├── metabolite1_vs_metabolite2.svg
    └── ...
```

## Key Analysis Features

### Statistical Testing
- **Multiple Testing Correction**: Benjamini-Hochberg FDR control
- **Effect Size Calculation**: Log2 fold changes using raw data
- **Significance Thresholds**: FDR < 0.15 for significance

### Advanced Correlation Analysis
- **Fisher's Z-test**: Statistical comparison of correlations between groups
- **Clustering**: Hierarchical clustering for metabolite relationships
- **Coefficient of Variation**: Metabolite variability assessment

### Machine Learning
- **Feature Engineering**:
  - Local correlation calculation
  - Group-centered interaction terms
  - Polynomial features
- **Model Evaluation**:
  - AUC-ROC analysis
  - Sensitivity/Specificity metrics
  - F1-score calculation
  - Cross-validation performance

### Visualization
- **Volcano Plots**: Publication-ready differential analysis
- **Heatmaps**: Clustered correlation matrices
- **ROC Curves**: Multiple model comparison
- **Scatter Plots**: Significant correlation pairs

## Usage

### Complete Pipeline Execution
```r
# Run entire analysis
results <- main_analysis_pipeline()
```

### Individual Step Execution
```r
# Step-by-step analysis
data_list <- load_and_preprocess_data()
stat_results <- perform_statistical_analysis(data_list$normalized, data_list$raw)
correlation_results <- perform_correlation_analysis(data_list$normalized)
ml_results <- perform_ml_analysis(data_list$normalized)
validation_results <- perform_validation_analysis()
```

## Results Interpretation

### Statistical Results
- **Significant Metabolites**: FDR < 0.15
- **Fold Changes**: Log2(MASLD/NC) ratios
- **Effect Sizes**: Clinical relevance assessment

### Correlation Differences
- **Fisher's Z p-values**: Statistical significance of correlation differences
- **Biological Relevance**: Metabolic pathway disruption indicators

### Machine Learning Performance
- **AUC Values**: Model discrimination ability
- **Cross-validation**: Generalization assessment
- **Feature Importance**: Key metabolite identification

## Reproducibility

The pipeline ensures reproducibility through:
- **Seed Management**: Fixed random seeds (set.seed(42))
- **Session Information**: Complete package versions
- **Parameter Documentation**: All analysis parameters recorded

## Troubleshooting

### Common Issues
1. **Missing Packages**: Install all required packages before execution
2. **Data Format**: Ensure proper Excel sheet structure
3. **Memory Usage**: Large correlation matrices may require sufficient RAM
4. **File Paths**: Verify data directory structure

### Performance Optimization
- **Parallel Processing**: Can be enabled in LMSstat analysis
- **Memory Management**: Clear workspace between major steps
- **Output Control**: Selective plot generation for large datasets

## Citation and References

If you use this pipeline, please consider citing relevant packages:
- **LMSstat**: Statistical analysis framework
- **caret**: Machine learning infrastructure  
- **pROC**: ROC analysis
- **pheatmap**: Heatmap visualization

## Author Information

- **Script**: MASLD Metabolomics Analysis Pipeline
- **Language**: R
- **License**: Open source
- **Maintenance**: Active development

## Version History

- **v1.0**: Initial comprehensive pipeline
- **Features**: Complete statistical, correlation, and ML analysis
- **Validation**: External dataset integration

---

*This pipeline provides a comprehensive framework for metabolomics analysis with focus on reproducibility and clinical relevance.*

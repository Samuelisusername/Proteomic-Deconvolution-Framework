# Proteomic Deconvolution Framework

A comprehensive computational framework for benchmarking and evaluating proteomic-based cellular deconvolution methods, with a focus on lymphoma research and personalized treatment strategies.

## Overview

Lymphoma is a blood cancer affecting lymphocytes, a crucial component of our immune system. While RNA-based algorithms have been developed to estimate cell composition in tumor microenvironments, 11%-38% of patients still develop drug resistance, with poor survival outcomes (18 months for progressive CLL, 4 months for Richter's transformation). This highlights the critical need for better diagnostic tools and personalized treatment approaches.

This framework addresses the gap between transcriptomic and proteomic deconvolution by providing one of the first comprehensive benchmarking platform specifically designed for protein expression data analysis.

## What This Framework Does

### Core Functionality
- **Synthetic Data Generation**: Creates realistic bulk proteomic samples with customizable cell-type distributions (healthy and cancerous)
- **Method Comparison**: Systematically evaluates three deconvolution methods:
  - CIBERSORT
  - BayesPrism  
  - NNLS (Non-Negative Least Squares)
- **Signature Matrix Construction**: Supports flexible reference matrix generation from Dirichlet-distributed combinations or fixed single-patient matrices
- **Comprehensive Evaluation**: Standardized metrics for comparing method performance across different conditions

### Key Features
- **Biological Foundation**: Built on extensive data exploration including PCA visualization of immune cell subsets across ~10,000-dimensional proteomic space
- **Preprocessing Strategies**: Multiple normalization approaches optimized for protein expression data
- **Bias Analysis**: Systematic characterization of method-specific biases and limitations
- **Clinical Relevance**: Demonstrates 100% classification accuracy for distinguishing healthy vs. cancerous samples using B-cell fraction thresholds

## Research Findings

### Method Performance
- **NNLS**: Most robust method with outlogged normalization.
- **CIBERSORT**: Strong performance but sensitive to preprocessing (1.16% vs 2.04% error for inlogged normalizations and outlogged respectively)
- **BayesPrism**: Underperformed but showed potential for signature matrix improvement

## Research Findings

### Method Performance
Our comprehensive benchmarking revealed distinct performance characteristics:

- **NNLS**: Most robust method with mean absolute errors as low as 1.26% under optimal conditions
  - Deterministic results with no hyperparameter tuning required
  - Near-perfect correlations (Pearson > 0.99 for B and T cells)
  - Fastest runtime compared to other methods

- **CIBERSORT**: Best performance with inlogged preprocessing but preprocessing-sensitive
  - Significant difference between inlogged (1.00% error) vs outlogged (2.04% error) normalization
  - Requires careful data preprocessing optimization
  - Good balance of accuracy and interpretability

- **BayesPrism**: Underperformed but showed improvement potential
  - Systematic overshooting in signature matrix updates
  - Correct update direction suggests potential for optimization
  - May benefit from modified scaling approaches

### Key Insights
- **Signature Matrix Construction**: Randomized multi-individual matrices reduced prediction errors by ~50% vs single-patient references
- **Systematic Biases**: All methods showed consistent overestimation of T cells and underestimation of B/myeloid cells
- **Clinical Validation**: 100% classification accuracy for healthy vs. cancerous samples using B-cell fraction thresholds
- **Distribution Sensitivity**: Performance degraded with highly skewed cell compositions typical in pathological conditions

## Why This Matters

### Cancer Cell Fraction Estimation
- **Improved Diagnosis**: Better cell composition estimates for lymphoma patients
- **Personalized Treatment**: Foundation for developing patient-specific therapeutic strategies
- **Drug Resistance Monitoring**: Tools to better understand and predict treatment resistance

### Research Impact
- **Reproducible Research**: Open-source framework enables standardized method comparison
- **Future Development**: Provides baseline for developing improved proteomic deconvolution algorithms
- **Cross-Platform Validation**: Framework supports testing across different proteomic technologies

## Getting Started

### Prerequisites
- Python 3.x
- Required packages:
  - `pandas`
  - `numpy`
  - `anndata`
  - `scanpy`
  - `matplotlib`
  - `rpy2` (for BayesPrism integration)
  - `openpyxl`

### Installation
```bash
git clone https://github.com/[username]/proteomic-deconvolution-framework
cd proteomic-deconvolution-framework
pip install -r requirements.txt
```

### Quick Start

#### Basic Usage
```bash
# Generate samples for CIBERSORT and NNLS
python3 choose_input_frac.py 10 5 ciber_and_nnls 1

# When prompted, enter fractions:
# Healthy fractions (NK, B, Myeloid, T): 0.04227402 0.09125808 0.22551629 0.64095161
# Cancerous fractions (B, T, Myeloid): 0.82520507 0.12311171 0.03607365 0.01560957
# Signature matrix type: rand

# Generate samples for BayesPrism
python3 choose_input_frac.py 10 5 bayes 1
```

#### Command Line Parameters
- `<number_of_healthy_bulk_samples>`: Number of healthy samples to generate
- `<number_of_cancerous_bulk_samples>`: Number of cancerous samples to generate  
- `<method>`: Either "ciber_and_nnls" or "bayes"
- `<job_id>`: Unique identifier for this run

#### Input Requirements
You'll need the proteomic reference data file:
- `41590_2017_BFni3693_MOESM10_ESM.xlsx` - Contains the proteomic profiles from healthy patients

## Framework Components

### 🔍 Run NNLS Deconvolution

The `nnls.py` script runs **Non-Negative Least Squares (NNLS)** deconvolution to estimate cell-type fractions from synthetic bulk proteomic samples using a predefined signature matrix.

#### ✅ Usage

```bash
python3 nnls.py <normalization>
```

**Arguments:**
- `<normalization>`: Normalization type used when generating the data. Must be one of:
  - `inlogged`
  - `outlogged`
  - `nonlogged`

#### 💡 Example

```bash
python3 nnls.py outlogged
```

This command will:
- Load `imputed_sig_matrix_outlogged.txt` as the signature matrix
- Load `sample_outlogged_imputed.txt` as the bulk sample file
- Perform NNLS-based estimation using `LinearRegression(positive=True)`
- Save the predicted cell fractions to `NNLS-Results_outlogged.txt`

#### 📥 Input Files

Make sure the following files exist in the working directory:
- `imputed_sig_matrix_<normalization>.txt` – Signature matrix (rows: proteins, columns: cell types)
- `sample_<normalization>_imputed.txt` – Bulk proteomic samples (rows: proteins, columns: samples)
- `real_fracs.tsv` – Ground truth cell fractions (used internally for column mapping)

#### 📤 Output

- `NNLS-Results_<normalization>.txt`: Predicted cell-type fractions (tab-delimited, one row per sample)

#### ⚙️ Notes
- Uses `scikit-learn`'s `LinearRegression(positive=True)` for NNLS approximation.
- Estimated coefficients are normalized to sum to 1 per sample.
- R² scores for model fits are printed to the console.
- Column names are internally mapped to match the proteomic dataset structure.

#### 📦 Dependencies

These packages must be installed (already listed in `requirements.txt`):

```txt
numpy
pandas
scikit-learn
matplotlib
scipy
```


## Framework Components

### Data Generation (`choose_input_frac.py`)
The core script that generates synthetic bulk proteomic samples with realistic cell-type distributions:

#### Key Features:
- **Three Normalization Strategies**:
  - `inlogged`: Log transformation applied before combining individual profiles
  - `outlogged`: Log transformation applied after combining profiles into bulk samples  
  - `nonlogged`: No log transformation applied
- **Flexible Signature Matrix Construction**:
  - `fixed`: Uses single patient profiles (patient 04)
  - `rand`: Creates randomized matrices using Dirichlet-distributed combinations
- **Realistic Cell Compositions**: Based on healthy vs. lymphoma patient distributions

#### Generated Outputs:
- `imputed_sig_matrix_<normalization>.txt` - Signature matrices for CIBERSORT/NNLS
- `sample_<normalization>_imputed.txt` - Bulk samples for deconvolution
- `real_fracs.tsv` - Ground truth cell fractions (26 cell types)
- `real_coarse_fracs.tsv` - Ground truth coarse fractions (4 main cell types)
- `myinput2.gbm.rdata` - BayesPrism-formatted data (when using "bayes" option)

### Cell Type Hierarchy
The framework models 26 fine-grained cell types organized into 4 coarse categories:

**B cells**: B.memory, B.naive, B.plasma  
**T cells**: T4.CM, T4.EM, T4.EMRA, T4.naive, T8.CM, T8.EM, T8.EMRA, T8.naive, Th1, Th17, Th2, mTregs, nTregs  
**Myeloid cells**: Basophil, Eosinophil, MO.classical, MO.intermediate, MO.nonclassical, Neutrophil, mDC, pDC  
**NK cells**: NK.bright, NK.dim

### Evaluation Framework (`evaluate_results.py`)
Comprehensive analysis script for comparing predicted vs. ground truth cell fractions:

#### Key Features:
- **Correlation Analysis**: Computes Pearson and Spearman correlations for each cell type
- **Error Quantification**: Calculates mean absolute errors and visualizes prediction errors
- **Comparative Visualization**: Generates stacked bar charts, scatter plots, and error distributions
- **Multi-level Analysis**: Evaluates both fine-grained (26 cell types) and coarse-grained (4 main types) predictions

#### Generated Outputs:
- `{method}_correlations.csv` - Correlation coefficients for each cell type
- `{method}_fractions_compare.csv` - Predicted vs. real fractions for all samples
- `{method}_errors.png` - Error visualization plots
- `{method}_all_celltypes_percents_correlation.png` - Scatter plots of predicted vs. real values

#### Usage:
```bash
python3 evaluate_results.py <results_file> <method_name>
```

Where:
- `results_file`: Tab-separated file containing deconvolution results
- `method_name`: Identifier for the method (e.g., "cibersort", "nnls", "bayesprism")

### Deconvolution Methods
- **NNLS (Non-Negative Least Squares)**: Most robust method with deterministic results
- **CIBERSORT**: Strong performance but sensitive to preprocessing choices
- **BayesPrism**: Joint estimation of cell fractions and signature matrices

### Evaluation Framework
- **Standardized Metrics**: Mean Absolute Error (MAE), Pearson correlations
- **Bias Analysis**: Systematic characterization of method-specific biases
- **Clinical Validation**: Classification accuracy for healthy vs. cancerous samples

## Limitations and Future Directions

### Current Limitations
- Reference signatures from healthy patients only
- Performance degradation with highly skewed cell distributions
- Systematic biases toward specific cell types

### Future Development Priorities
1. **Expanded Reference Datasets**: Disease-specific proteomic profiles
2. **Bias Correction**: Post-processing algorithms for systematic error correction
3. **Uncertainty Quantification**: Confidence intervals for predictions
4. **Cross-Platform Validation**: Testing across different proteomic technologies

## Contributing

We welcome contributions to improve the framework! Please see [CONTRIBUTING.md] for guidelines on:
- Adding new deconvolution methods
- Implementing additional normalization strategies
- Expanding evaluation metrics
- Improving documentation


## Acknowledgments

This work contributes to the broader effort to develop better diagnostic and treatment tools for lymphoma patients, ultimately advancing personalized medicine in oncology.

## Contact

For questions, suggestions, or collaborations, please contact:
- Samuel Gair, sagair@ethz.ch

---

**Note**: This framework represents ongoing research in proteomic deconvolution. Results should be validated with additional datasets before clinical application.

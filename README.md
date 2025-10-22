# Testing for Missingness Patterns in Incomplete Functional Data

This repository provides an R package and simulation framework designed to test the **Missing Completely at Random (MCAR)** assumption in partially observed functional data.\
It was developed as part of the Bachelor's Thesis *"Testing for Missingness Patterns in Incomplete Functional Data"* at **Ludwig-Maximilians-Universität München (LMU)**.

The project implements and evaluates the mean-based MCAR tests proposed by [**Ofner et al. (2025)**](https://arxiv.org/abs/2505.08721), offering both **asymptotic** and **bootstrap-based** inference methods, along with simulation studies and visualization for empirical validation.

The repository contains:

-   Full R source code of the `mcartest` package
-   Simulation and bootstrap implementations
-   Replication scripts for all data examples
-   Generated figures and type I error tables

## Motivation

In modern statistics, observations have evolved from individual numerical values to entire **trajectories over time or space**, such as medical monitoring curves (e.g., ECG and EEG), intraday electricity demand, satellite temperature profiles, or yield curves. When each observation is a function rather than a single value, traditional methods must be adapted to handle estimation, inference, and diagnostic tasks across continuous domains. These challenges are addressed by **Functional Data Analysis (FDA)**, a framework that offers flexible methods for estimating, interpreting, and comparing entire functions.

In practice, however, these functional observations are often **incomplete**. Sensors may fail, participants may remove monitoring devices, or recordings may stop prematurely. Such interruptions create missing segments that can distort estimators of the mean and covariance functions and lead to **biased** conclusions if ignored. Furthermore, it is crucial to determine whether the missing parts occur **completely at random** or are systematically related to the underlying process. Examining this relationship is a crucial step before applying standard FDA techniques, as violations can lead to biased estimation and misleading conclusions.

This thesis addresses this problem by implementing and empirically evaluating the statistical tests for the **Missing Completely at Random (MCAR)** assumption proposed by [**Ofner et al. (2025)**](https://arxiv.org/abs/2505.08721).\
Rather than developing a new methodology, the focus lies on reproducing and validating their results through a dedicated R implementation within the **tidyfun** framework (Scheipl et al., 2025).\
The objective is to provide a practical, reproducible framework that confirms the findings of the original study and demonstrates that the implemented methods perform as intended.\
Further details, theoretical background, and the full documentation of the replication study can be found in the **Bachelor’s Thesis**, located in the folder **`/Thesis/`**.

## Implemented Tests

The package implements the **two mean-based MCAR tests** proposed by [**Ofner et al. (2025)**](https://arxiv.org/abs/2505.08721):

-   $T_{\mu,L2}$ — compares group means using the $L^2$-norm\
-   $T_{\mu,D}$ — compares group means using the supremum (sup-) norm

In addition, the package provides **simultaneous confidence bands** for the mean difference $\hat{\mu}_A - \hat{\mu}_B$, constructed either asymptotically or via bootstrapped quantiles, which serve as a graphical tool to visualize potential deviations between the two groups.

Both tests are available in two variants:

-   **Asymptotic version:** Approximates the null distribution of the test statistics based on the estimated covariance function and its eigenvalues.

-   **Bootstrap version:** Approximates the distribution of the test statistic by **resampling the data groupwise**, i.e., separately within each group $A$ and $B$ under the null hypothesis of equal means. This approach flexibly adapts to finite-sample structures, avoids the need for eigen-decomposition, and improves numerical stability and small-sample performance. In practice, the bootstrap procedure is **parallelized** to efficiently handle a large number of resamples while ensuring reproducibility through controlled random seeds.

> 🔹 **Note:** The distributional test $T_F$ from Ofner et al. (2025) is not included, as the focus of this thesis lies on the implementation and evaluation of mean-based tests and their bootstrap approximations.

## Evaluation

All methods were evaluated on both **simulated** and **real-world datasets** (heart rate, electricity prices, and temperature series) from [Ofner et al. (2025)](https://arxiv.org/abs/2505.08721). The results confirm that the implemented methods work as intended and align with the findings of the original paper, providing an accessible and reliable R tool for practitioners.

## Project Structure

```         
Bachelor_FDA/
├── cran_package/
│   └── mcartest/
│     ├── DESCRIPTION          
│     ├── NAMESPACE            # Export/import declarations
│     ├── LICENSE              
│     ├── LICENSE.md           
│     ├── .gitignore           
│     ├── .Rbuildignore        
│     ├── .Rhistory            
│     ├── mcartest.Rproj       
│     ├── R/                   
│     │ ├── asym_tests.R         # Asymptotic tests 
│     │ ├── boot_test.R          # Bootstrap test 
│     │ ├── utils_backend.R      # Parallel backend handling
│     │ ├── utils_prepare.R      # Data preparation and setup
│     │ ├── utils_stats.R        # Statistical computations
│     │ └── mcartest_package.R   
│     ├── man/                 # roxygen2-generated documentation
│     └── tests/               # Unit tests 
│
├── Code/
│ ├── brownian_motion.R
│ ├── temperature.R
│ ├── heart_rate.R
│ ├── electricity_market.R
│ ├── grouping.R
│ └── x_o.R
│
├── Plots/
│ ├── brownian_motion_plot.png
│ ├── brownian_motion_rej_probs_plot.png
│ ├── temperature_plot.png
│ ├── temperature_rej_probs_plot.png
│ ├── heart_rate_plot.png
│ ├── electricity_market_plot.png
│ ├── grouping_partition_1_plot.png
│ ├── grouping_partition_2_plot.png
│ └── x_o_plot.png
│
├── Data/
│ ├── brownian_motion_rej_probs_df.rds
│ ├── heart_rate.RData
│ ├── logbidcurves.csv
│ ├── temp_graz.rda
│ ├── type_1_error_n100.rds
│ ├── type_1_error_n250.rds
│ └── type_1_error_n500.rds
│
├── Thesis/
│ └── Bachelor_Thesis.pdf   
│
├── requirements.R      # Installs and loads all required packages and fonts
├── .gitignore
├── Bachelor_FDA.Rproj       
└── README.md
```

## Requirements/Getting Started

-   **R** (≥ 4.3.0)
-   Before using the package, make sure that all required dependencies and fonts are installed and loaded by running:

``` r
# --- Install and load all required packages and fonts ---
source("requirements.R")
```

-   If the package is not loaded with `library(mcartest)`, all functions can be made available manually using:

``` r
# --- Load all package functions ---
devtools::load_all("cran_package/mcartest")
```

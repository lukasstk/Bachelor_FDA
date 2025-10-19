\section*{MCARtest -- Testing the Missing Completely at Random (MCAR) Assumption for Functional Data}

This repository contains the R source code of the \textbf{MCARtest} package, developed as part of the Bachelor's Thesis 
\emph{``Testing for Missingness Patterns in Incomplete Functional Data''} at Ludwig-Maximilians-Universität München (LMU). 
The package implements statistical tests to verify whether partially observed functional data satisfy the 
\textbf{Missing Completely at Random (MCAR)} assumption -- a key prerequisite for many functional data analysis (FDA) methods.

The theoretical foundation is based on:  
\begin{quote}
Ofner, M., Hörmann, S., Kraus, D., \& Liebl, D. (2025).  
\emph{Testing the Missing Completely at Random Assumption for Functional Data.}  
\texttt{arXiv:2505.08721v1}
\end{quote}

\subsection*{Conceptual Overview}

Functional data analysis (FDA) studies samples of curves or trajectories over a continuous domain (e.g., time). 
In many real-world settings, these functions are only \textbf{partially observed}, i.e., available only on subsets of their domain.  
Before applying standard FDA methods, one must check that the missingness mechanism is independent of the data itself.

The \texttt{MCARtest} package provides:
\begin{itemize}
  \item Asymptotic and bootstrap-based tests for the MCAR assumption
  \item Simultaneous confidence bands for group mean comparisons
  \item A framework for detecting dependence between the data process $X$ and its observation process $O$
\end{itemize}

Implemented test statistics include:
\begin{itemize}
  \item $T_{\mu,L2}$ and $T_{\mu,D}$ for mean-function comparison
  \item $T_F$ for conditional distribution comparison
  \item Asymptotic and bootstrap variants for all tests
\end{itemize}

\subsection*{Folder Structure}

\begin{verbatim}
Bachelor_FDA/
├── cran_package/
│   └── mcartest/
│       ├── DESCRIPTION                 # Package metadata
│       ├── NAMESPACE                   # Export/import declarations
│       ├── R/                          # Core implementation
│       │   ├── asym_tests.R            # Asymptotic test statistics (Tµ,L2, Tµ,D, TF)
│       │   ├── boot_test.R             # Bootstrap test variants
│       │   ├── utils_backend.R         # Parallel backend handling
│       │   ├── utils_prepare.R         # Subdomain selection, data preparation
│       │   ├── utils_stats.R           # Covariance estimation and helper functions
│       │   └── mcartest_package.R      # Main package entry point
│       ├── man/                        # Function documentation (roxygen2)
│       └── tests/                      # Unit tests for internal and exported functions
│
├── simulations/
│   ├── brownian_motion.R
│   ├── temperature.R
│   ├── heart_rate.R
│   ├── electricity_market.R
│   ├── grouping.R
│   └── x_o.R
│
├── plots/
│   ├── brownian_motion_plot.png
│   ├── brownian_motion_rej_probs_plot.png
│   ├── temperature_plot.png
│   ├── temperature_rej_probs_plot.png
│   ├── heart_rate_plot.png
│   ├── electricity_market_plot.png
│   ├── grouping_partition_1_plot.png
│   ├── grouping_partition_2_plot.png
│   └── x_o_plot.png
│
├── thesis/
│   └── [LaTeX source files and final written report]
│
├── .gitignore
├── Bachelor_FDA.Rproj
└── README.tex
\end{verbatim}

\subsection*{Usage and Dependencies}

\textbf{Important:} The package functions are interdependent. You must load the entire package (not individual files) 
for it to work properly.

From the main project folder (\texttt{Bachelor\_FDA/}), run:

\begin{verbatim}
devtools::load_all("cran_package/mcartest")
\end{verbatim}

or set the working directory directly to the package:

\begin{verbatim}
setwd("cran_package/mcartest")
devtools::load_all(".")
\end{verbatim}

If \texttt{devtools} is not installed, you can install it using:

\begin{verbatim}
install.packages("devtools")
\end{verbatim}

Alternatively, you can load all R scripts manually:
\begin{verbatim}
pkg_path <- "cran_package/mcartest/R"
files <- list.files(pkg_path, pattern = "\\.R$", full.names = TRUE)
sapply(files, source)
\end{verbatim}

\textbf{Dependencies:}
\begin{verbatim}
checkmate, tidyfun, stats, foreach, doParallel, doRNG, matrixStats
\end{verbatim}

\subsection*{Unit Tests}

Unit tests are located in:
\begin{verbatim}
cran_package/mcartest/tests/testthat/
\end{verbatim}

They verify:
\begin{itemize}
  \item Correct covariance estimation under partial observation
  \item Subdomain coverage and selection logic
  \item Type I error control of bootstrap tests
  \item Parallel backend setup and shutdown
\end{itemize}

Run all tests with:
\begin{verbatim}
devtools::test("cran_package/mcartest")
\end{verbatim}



# Shrinkage Approaches for Ridge-Type Estimators under Multicollinearity - Codes

This repository contains the R scripts, simulation files, and empirical datasets supporting the paper:

> **Marwan Al-Momani\*, Bahadır Yüzbaşı, M.S. Bataineh, Rihab Abdallah, and Athifa Moideenkutty.**  
> *Shrinkage Approaches for Ridge-Type Estimators under Multicollinearity.*  
> University of Sharjah & Inonu University, 2025.

---

## Abstract

Multicollinearity is a common issue in regression analysis that occurs when some predictor variables are highly correlated, leading to unstable least squares estimates 
of model parameters. Various estimation strategies have been proposed to address this problem. In this study, we enhance a ridge-type estimator by incorporating pretest and 
shrinkage techniques. We conducted an analytical comparison to evaluate the performance of the proposed estimators in terms of bias, quadratic risk, and numerical performance 
using both simulated and real data. Additionally, we assessed several penalization methods and three machine-learning algorithms to facilitate a comprehensive comparison. 
Our results demonstrate that the proposed estimators outperform the standard ridge-type estimator with respect to mean squared error in simulated data and mean squared 
prediction error in real data applications.

Both **Monte Carlo simulation experiments** and **real-data applications** were conducted to evaluate the performance of the proposed estimators.

This repository provides full replication materials for the study:
- Parameter settings and simulation scripts.
- Empirical datasets used in the real-data analysis.
- Output summaries (bias, MSE, prediction error metrics).

---


## Simulation Study

### Purpose
The Monte Carlo simulation experiments were conducted to examine:
- The effect of multicollinearity (correlation among predictors) on estimator performance.
- The bias and efficiency of the proposed ridge-type shrinkage estimators compared with classical estimators (OLS, Ridge, Liu, etc.).

### Main Settings

| Parameter | Symbol | Value / Description |
|------------|----------|----------------------|
| Sample size | \(n\) | 100 |
| Number of predictors | \(p\) | 15 |
| Number of restrictions | \(q\) | 10 |
| Ridge correlation parameter | \(\rho\) | 0.95 |
| Noise variance | \(\tau^2\) | 1 |
| Confidence level | \(\alpha\) | 0.05 |
| Shrinkage step | \(d\) | 0.1 |
| Liu-type | \(k\) | optimally selected |
| Baseline parameter | \(\delta_0\) | 0 |
| Number of iterations | |`iter = 1000` |
| Random seed | | 1973 |


- Sample sizes: *n = 25, 50, 100*
- Number of regressors: *p = 4, 6, 8*
- Correlation levels: *ρ = 0.80, 0.90, 0.99*
- Replications: *10,000*
- Performance metrics: **Bias**, **MSE**, **RMSE**, and **Predictive Error (PE)**

### Running the Simulation
```R
setwd("simulation/")
source("sim_main.R")
```

## Real Data Applications

### Running the applications
```R
setwd("R/")
source("ridge_estimators.R")




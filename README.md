# discrete-hull-white-caplets-floorlets

Discrete-time Hull-White pricing of backward-looking caplets and floorlets — R code for the thesis figures and sensitivity analysis.

## Overview

This repository contains the numerical implementation accompanying my Bachelor’s thesis,

*Caplet and Floorlet Pricing under Backward-Looking Compounded Rates in a Discrete Hull-White Framework*  
University of Padova, 2025/2026.

The thesis prices backward-looking caplets and floorlets in a discrete Hull-White model, where the short rate follows a first-order autoregressive recursion, obtained as the Euler-Maruyama discretization of the Hull-White equation. Under the single-period convention `p(k, k+1) = exp(-r_k)`, the compounded rate telescopes exactly into the exponential of the sum of short rates, and the caplet/floorlet prices follow in closed form from truncated Gaussian moments.

The code reproduces the numerical analysis of Chapter 4: the variance decomposition `nu_X = nu_W + nu_{m,N}`, the sensitivity of the price to volatility, mean reversion, strike, the pre-accrual and accrual windows, and the backward-versus-forward price gap.

## Contents

- `bachelor-thesis-numerical-implementation.R` — single self-contained script. It defines the variance components and the closed-form caplet/floorlet prices, prints the base-case diagnostics corresponding to Table 4.1, and draws the figures of Sections 4.2 and 4.4.

## Requirements

- Base R only — no external packages.

## Author

Sebastiano Pinotti  
Bachelor’s thesis in Statistics for Economics and Business  
University of Padova, Department of Statistical Sciences  

Supervisor: Prof. Massimiliano Caporin  
Co-supervisor: Prof. Claudio Fontana

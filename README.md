# discrete-hull-white-caplet-floorlet

Discrete-time Hull–White pricing of backward-looking caplets and floorlets — R code for the thesis and the accompanying working paper.

## Overview

This repository contains the numerical implementation accompanying:

- My Bachelor's thesis, *Caplet and Floorlet Pricing under Backward-Looking Compounded Rates in a Discrete Hull-White Framework*, University of Padova, 2025/2026.
- The working paper *Backward-Looking Caplet Pricing in a Discrete Hull–White Model: Closed Forms, Variance Decomposition, and the Backward–Forward Price Gap* (2026), which grew out of the thesis.

Both works study backward-looking caplets and floorlets in a discrete Hull-White model, where the short rate follows a first-order autoregressive recursion, obtained as the Euler-Maruyama discretization of the Hull-White equation. Under the single-period convention `p(k, k+1) = exp(-r_k)`, the compounded rate telescopes exactly into the exponential of the sum of short rates, and the caplet/floorlet prices follow in closed form from truncated Gaussian moments. The two conventions share a single Black-type formula and differ only through the variance decomposition `nu_X = nu_W + nu_{m,N}`, which governs the backward–forward price gap.

## Contents

- `thesis-caplet-floorlet-sensitivity-analysis.R` — self-contained script reproducing the numerical analysis of Chapter 4 of the thesis: base-case diagnostics (Table 4.1), the variance decomposition, the sensitivity of the price to volatility, mean reversion, strike, the pre-accrual and accrual windows, and the backward-versus-forward price gap (Sections 4.2 and 4.4).
- `paper-caplet-sensitivity-analysis.R` — self-contained script reproducing all figures and numbers of Section 4 of the working paper. It saves the figures as PNG files to `figures/` and recomputes, with the paper's quoted values alongside, every number cited in the text.

The scripts are independent, each runs on its own.

## Requirements

- Base R only — no external packages.

## Author

Sebastiano Pinotti  
Bachelor's degree in Statistics for Economics and Business  
University of Padova, Department of Statistical Sciences  

Supervisor: Prof. Massimiliano Caporin  
Co-supervisor: Prof. Claudio Fontana

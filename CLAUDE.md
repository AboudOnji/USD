# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Financial quantitative model for simulating and forecasting the USD/MXN exchange rate using stochastic processes. Implements a Merton (1976) jump-diffusion model over a 3-year horizon (March 2026–March 2029).

**Author**: Prof. D.Sc. BARSEKH-ONJI Aboud, Facultad de Ingeniería, Universidad Anáhuac México

## Running the Simulation

This is a single-file MATLAB project. Run from the MATLAB command window:

```matlab
>> trial1
```

Outputs 3 PNG files to `/mnt/user-data/outputs/`:
- `usdmxn_trayectorias.png` — trajectory bands with confidence intervals
- `usdmxn_distribuciones.png` — distribution histograms at years 1, 2, 3
- `usdmxn_fanchart.png` — fan chart showing percentile evolution

## Architecture

All logic lives in `trial1.m` (370 lines), organized into 5 sequential sections:

| Section | Lines | Purpose |
|---------|-------|---------|
| 1 | 23–62 | Model parameters (hardcoded) |
| 2 | 65–115 | Monte Carlo simulation engine |
| 3 | 120–146 | Statistical analysis |
| 4 | 149–355 | Visualization (3 figures) |
| 5 | 357–370 | Export PNG files |

## Model Parameters (Section 1)

Key calibrated values — edit these to recalibrate the model:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `S0` | 17.50 | Initial MXN/USD rate (Mar 2026) |
| `mu` | 0.025 | Annual drift (2.5% depreciation) |
| `sigma` | 0.115 | Annual volatility (11.5%) |
| `lambda` | 3 | Jump frequency (events/year, Poisson) |
| `mu_j` | -0.01 | Mean jump size |
| `sigma_j` | 0.03 | Jump size volatility |
| `dt` | 1/252 | Daily time step (trading days) |
| `T` | 3 | Horizon in years |
| `N_sim` | 1000 | Number of Monte Carlo paths |
| `seed` | 42 | Random seed for reproducibility |

## Stochastic Model

The simulation uses Geometric Brownian Motion with Poisson jumps (Merton 1976):

```
dS = S·[μ - λ·(e^(μ_j + σ_j²/2) - 1)]·dt  +  S·σ·dW  +  S·J·dN
```

where `dW` is Brownian motion, `dN` is a Poisson process, and `J` is log-normal jump size.

## Notes

- Comments and output messages are in Spanish
- Historical calibration references: 2023–2024 Banxico data embedded in Section 1 comments
- The model corrects for jump-diffusion drift using the Merton (1976) formula (visible in Section 2)

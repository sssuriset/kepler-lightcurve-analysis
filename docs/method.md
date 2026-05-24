# Method

This project analyzes a Kepler light curve stored in FITS format.

## Data Loading

The script loads the first table extension in the FITS file. It uses PDCSAP_FLUX when available because this column contains corrected flux values. If PDCSAP_FLUX is unavailable, the script falls back to SAP_FLUX.

## Data Cleaning

The script removes points that have invalid time or flux values, nonpositive flux, or nonzero Kepler QUALITY flags. Removed points are saved separately so the cleaning process is traceable.

## Flux Normalization

The cleaned flux is divided by its median value. This produces a normalized light curve centered near 1.

## Period Search

The script uses a Lomb-Scargle periodogram to search for periodic structure. Lomb-Scargle is commonly used in astronomy because it works with unevenly sampled time-series data.

## Period Uncertainty

The period uncertainty is estimated from the width of the periodogram peak. This is a simplified diagnostic uncertainty, not a full posterior uncertainty from a physical model.

## Sinusoidal Model

A sinusoidal model is fit near the Lomb-Scargle period. The model provides a simple periodic baseline for visualization and residual analysis.

## Phase Folding

The light curve is folded on the fitted period. If the period estimate is meaningful, repeated structure should align in phase.

## Residuals

Residuals are calculated by subtracting the sinusoidal model from the normalized flux. The residual plot and RMS value help evaluate how much variation remains after the simple periodic model.

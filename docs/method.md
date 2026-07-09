# Methods

Both pipelines share the period-search implementations in `src/period_search.py`. This document covers the method choices in each.

## Kepler pipeline

### Data loading and cleaning

The script loads the first table extension of the FITS file and uses PDCSAP_FLUX when available, since that column contains systematics-corrected flux, falling back to SAP_FLUX otherwise. Points with invalid time or flux, nonpositive flux, or nonzero Kepler QUALITY flags are removed and saved separately so the cleaning is traceable. Cleaned flux is normalized by its median.

### Period search and uncertainty

A Lomb-Scargle periodogram searches 0.5 to 40 days at 25,000 frequency samples. Lomb-Scargle handles the uneven sampling left after quality filtering. The period uncertainty is estimated from the half-height width of the periodogram peak. This is a diagnostic width, not a posterior uncertainty from a physical model, and it is reported alongside the much smaller formal uncertainty from the sine fit covariance so the difference stays visible.

### Sine fit, folding, residuals

A four-parameter sinusoid (offset, amplitude, period, phase) is fit with bounds keeping the period within 25 percent of the periodogram value. The light curve is folded on the fitted period as a visual check that repeated structure aligns in phase, and residuals from the fit give the RMS scatter that remains after the periodic model.

## Injection-recovery pipeline

### Synthetic light curve

The model starts from unit flux and injects a box transit with period 3.72 days, depth 0.018, and duration 0.18 days, centering each transit at phase zero and marking points within half the duration on each side. A slow sinusoidal trend and Gaussian noise (sigma 0.004) make the data less idealized. The random seed is fixed so runs reproduce.

### Two search methods

Lomb-Scargle is run as a deliberate mismatch: it targets smooth periodic signals, and on the box-shaped transit it returns a strongest period near 1.24 days. That value is kept in the outputs to show that method choice drives the result. Box Least Squares is the matched method. It tests trial periods against a grid of transit durations and recovers 3.72097 days from the injected 3.72, an error of about 1.4 minutes, with a recovered duration of 0.175 days against the injected 0.18.

### Signal-to-noise

The transit SNR is the injected depth divided by the residual scatter after removing the noiseless model, giving 4.36 for the current settings. This is a controlled-simulation estimate; a survey search would use more careful noise modeling and candidate vetting.

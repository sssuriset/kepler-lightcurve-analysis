from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

filename = "kplr000757450-2009350155506_llc.fits"

# Open FITS file
hdul = fits.open(filename)

# Exploratory information
hdul.info()
data = hdul[1].data
print(data.columns)

# Extract data
time = data["TIME"]
flux = data["PDCSAP_FLUX"]

# Clean data
mask = np.isfinite(time) & np.isfinite(flux)
time = time[mask]
flux = flux[mask]

# Normalize flux
flux = flux / np.mean(flux)

# Simple raw light curve plot
plt.figure(figsize=(10, 5))
plt.plot(time, flux, ".", markersize=1)
plt.xlabel("Time (BKJD)")
plt.ylabel("Normalized Flux")
plt.title("Kepler Light Curve")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("lightcurve.png", dpi=300)
plt.show()

# Remove mean for period search
flux_detrended = flux - np.mean(flux)

# Trial periods in days
periods = np.linspace(5, 40, 2000)
scores = []

for P in periods:
    omega = 2 * np.pi / P
    s = np.sin(omega * time)
    c = np.cos(omega * time)
    model = np.column_stack((s, c, np.ones_like(time)))
    coeffs, _, _, _ = np.linalg.lstsq(model, flux_detrended, rcond=None)
    fit_trial = model @ coeffs
    score = np.sum((fit_trial - np.mean(flux_detrended))**2)
    scores.append(score)

scores = np.array(scores)
best_period = periods[np.argmax(scores)]

print("Estimated period =", best_period, "days")

# Best-fit sinusoid
omega = 2 * np.pi / best_period
s = np.sin(omega * time)
c = np.cos(omega * time)
model = np.column_stack((s, c, np.ones_like(time)))
coeffs, _, _, _ = np.linalg.lstsq(model, flux_detrended, rcond=None)
fit = model @ coeffs + np.mean(flux)

# Plot 1: light curve with fit
plt.figure(figsize=(10, 5))
plt.plot(time, flux, ".", markersize=1, label="Observed flux")
plt.plot(time, fit, "-", linewidth=2, label="Best-fit periodic model")
plt.xlabel("Time (days)")
plt.ylabel("Normalized Flux")
plt.title("Kepler Light Curve with Periodic Fit")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("lightcurve_fit.png", dpi=300)
plt.show()

# Plot 2: score vs period
plt.figure(figsize=(10, 5))
plt.plot(periods, scores)
plt.xlabel("Trial Period (days)")
plt.ylabel("Fit Score")
plt.title("Period Search")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("period_search.png", dpi=300)
plt.show()

hdul.close()
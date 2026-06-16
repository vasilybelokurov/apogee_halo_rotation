# APOGEE Halo Rotation Analysis Notes

This note records the public, reproducible workflow implemented in `scripts/measure_apogee_halo_rotation.py`.

## External Data

The analysis uses three public files cached under `data/external/`:

- APOGEE DR17 allStarLite:
  `https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/allStarLite-dr17-synspec_rev1.fits`
- astroNN APOGEE DR17:
  `https://data.sdss.org/sas/dr17/apogee/vac/apogee-astronn/apogee_astroNN-DR17.fits`
- Vasiliev & Baumgardt 2021 Gaia EDR3 globular-cluster members:
  `https://zenodo.org/records/4891252/files/clusters.zip?download=1`

Use `python scripts/download_data.py` to download these files.

## Selection

The clean sample applies red-giant, APOGEE flag, abundance uncertainty, astroNN velocity uncertainty, distance, Magellanic Cloud program, duplicate/telluric, Gaia source-id, and Vasiliev & Baumgardt globular-cluster member-removal cuts. No orbital-energy window is applied. The GC member removal matches APOGEE `GAIAEDR3_SOURCE_ID` to the Zenodo catalogue `source_id` values with `memberprob > 0.8`.

The Belokurov & Kravtsov chemical split is:

- In-situ: `[Fe/H] > -0.4`.
- At `[Fe/H] <= -0.4`, in-situ: `[Al/Fe] >= -0.075`.
- At `[Fe/H] <= -0.4`, accreted: `[Al/Fe] < -0.075`.

The default rotation measurement uses `[Fe/H] < -1.0`. The z-slice comparison overlay uses `[Fe/H] < -1.3`.

## Measurement

The measured bulk rotation is the median astroNN `galvt` in bins of spherical Galactocentric radius,

```text
r_GC = sqrt(galr^2 + galz^2).
```

Bootstrap 16th/84th percentile intervals are computed for the median velocity. Velocity estimates are reported only for radial bins containing at least three stars. Velocity plots use fixed y-axis limits of `-50` to `125 km/s`. The radius histograms use the same radial bin edges and radial axis limits as the velocity plots, and use a base-10 logarithmic count axis. By default, plots are PNG only; PDF outputs are optional command-line arguments.

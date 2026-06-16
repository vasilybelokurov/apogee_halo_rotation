# APOGEE Halo Rotation

This repository measures the bulk azimuthal rotation of chemically selected in-situ and accreted Milky Way halo stars from APOGEE DR17 as a function of Galactocentric radius.

The sample definition follows the APOGEE/astroNN selection used by Belokurov & Kravtsov (2022), with an additional public Gaia EDR3 globular-cluster member removal from Vasiliev & Baumgardt (2021).

## Data

The external catalogues are not committed to git. Download them into `data/external/` with:

```bash
python scripts/download_data.py
```

The script fetches:

| Dataset | Repository cache path | Public source |
|---|---|---|
| APOGEE DR17 allStarLite chemistry and flags | `data/external/allStarLite-dr17-synspec_rev1.fits` | `https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/allStarLite-dr17-synspec_rev1.fits` |
| astroNN APOGEE DR17 distances, kinematics, actions, and energies | `data/external/apogee_astroNN-DR17.fits` | `https://data.sdss.org/sas/dr17/apogee/vac/apogee-astronn/apogee_astroNN-DR17.fits` |
| Vasiliev & Baumgardt Gaia EDR3 globular-cluster members | `data/external/vasiliev_baumgardt_2021_clusters.zip` | `https://zenodo.org/records/4891252/files/clusters.zip?download=1` |

The Vasiliev & Baumgardt catalogue is the Zenodo dataset "Catalogue of stars in Milky Way globular clusters from Gaia EDR3", DOI `10.5281/zenodo.4891252`.

## Sample Selection

The script joins allStarLite and astroNN on `APOGEE_ID`, then applies the following clean-sample cuts:

- Finite required APOGEE abundances, uncertainties, astroNN distances, velocities, velocity uncertainties, energy, and Galactocentric coordinates.
- Red giants: `LOGG < 3`.
- astroNN velocity uncertainties: `galvr_err`, `galvt_err`, `galvz_err < 50 km/s`.
- APOGEE abundance uncertainties: `C_FE_ERR`, `N_FE_ERR`, `O_FE_ERR`, `MG_FE_ERR`, `AL_FE_ERR`, `SI_FE_ERR`, `MN_FE_ERR`, `FE_H_ERR`, `NI_FE_ERR < 0.2 dex`.
- Remove APOGEE `ASPCAPFLAGS` containing `STAR_BAD`, `TEFF_BAD`, `LOGG_BAD`, `VERY_BRIGHT_NEIGHBOR`, `LOW_SNR`, `PERSIST_HIGH`, `PERSIST_JUMP_POS`, `PERSIST_JUMP_NEG`, or `SUSPECT_RV_COMBINATION`.
- Remove `ASPCAPFLAG` bit 23.
- Remove `EXTRATARG` duplicate bit `2**4` and telluric bit `2**2`.
- Remove `PROGRAMNAME` containing `magcloud`.
- Require non-zero Gaia EDR3 source id.
- Remove Vasiliev & Baumgardt GC candidates matched by Gaia source id with `memberprob > 0.8`.
- Require astroNN heliocentric distance `< 15 kpc`.
- Require the Belokurov & Kravtsov energy window: `-0.75 < Energy * 1e-5 < -0.4`.

The in-situ/accreted chemical split is:

- In-situ: `[Fe/H] > -0.4`.
- At `[Fe/H] <= -0.4`, in-situ: `[Al/Fe] >= -0.075`.
- At `[Fe/H] <= -0.4`, accreted: `[Al/Fe] < -0.075`.

The default rotation measurement further restricts the sample to `[Fe/H] < -1.0`. The z-slice comparison plot overlays dashed curves for `[Fe/H] < -1.3`.

## Measurement

The reported rotation is the median astroNN Galactocentric azimuthal velocity `galvt`; positive velocity is prograde. Radius is spherical Galactocentric radius,

```text
r_GC = sqrt(galr^2 + galz^2)
```

where `galr` is the astroNN cylindrical Galactocentric radius. Bootstrap 16th/84th percentile intervals are shown for the median velocity.

Run the analysis with:

```bash
python scripts/measure_apogee_halo_rotation.py
```

By default, plots are written as PNG files only. PDF files are made only when the corresponding `--out-*-pdf` arguments are passed.

## Results

Current run summary:

| Quantity | Value |
|---|---:|
| Joined APOGEE-astroNN rows | 657,134 |
| Vasiliev & Baumgardt GC source ids with `memberprob > 0.8` | 1,008,359 |
| Clean BK22-like sample | 89,993 |
| Measurement sample, `[Fe/H] < -1.0` | 2,608 |
| In-situ measurement stars | 1,098 |
| Accreted measurement stars | 1,510 |
| `\|z\| < 3 kpc` measurement stars | 1,583 |
| `\|z\| > 3 kpc` measurement stars | 1,025 |

Median rotation by spherical Galactocentric radius:

| Population | r range [kpc] | N | median Vphi [km/s] |
|---|---:|---:|---:|
| In-situ | 0-3 | 10 | 134.1 |
| In-situ | 3-5 | 123 | 158.3 |
| In-situ | 5-7 | 241 | 95.4 |
| In-situ | 7-9 | 540 | 93.1 |
| In-situ | 9-11 | 158 | 54.6 |
| In-situ | 11-13 | 26 | 9.3 |
| Accreted | 0-3 | 12 | 154.6 |
| Accreted | 3-5 | 87 | 119.0 |
| Accreted | 5-7 | 213 | 44.1 |
| Accreted | 7-9 | 587 | 19.5 |
| Accreted | 9-11 | 393 | 9.3 |
| Accreted | 11-13 | 200 | 7.0 |
| Accreted | 13-13.7 | 18 | -0.2 |

## Plots

Chemical selection:

![Chemical selection](plots/apogee_halo_chemical_selection.png)

Bulk rotation versus Galactocentric radius:

![Rotation versus radius](plots/apogee_halo_rotation_vs_radius.png)

Rotation split by height from the Galactic plane. Solid curves use `[Fe/H] < -1.0`; dashed curves use `[Fe/H] < -1.3`.

![Rotation z slices](plots/apogee_halo_rotation_z_slices.png)

Spatial distribution of the measured in-situ and accreted samples:

![R-z and z distributions](plots/apogee_halo_rz_z_distribution.png)

## References

- Belokurov, V. & Kravtsov, A. 2022, MNRAS, 514, 689.
- Vasiliev, E. & Baumgardt, H. 2021, MNRAS, 505, 5978.

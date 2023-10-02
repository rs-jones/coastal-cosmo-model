# Coastal-cosmo-model
Cosmogenic nuclide model for coastal erosion and sea level change scenarios

Written by Richard Selwyn Jones (richard.s.jones@monash.edu), using modifications of code from Swirad et al. (2020, Nat. Comms.; and refs. therein) for calculating time-dependent topographic and water shielding, and cliff retreat, and CRONUScalc (Marrero et al., 2016, Quat. Geochronology; and refs. therein) for calculating cosmogenic nuclide production. The model uses a MATLAB optimised solver to find the best-fit variables for a given case.

The model consists of various sub-models for different purposes, and to test hypotheses of increasing complexity:
- Cosmogenic inheritance model, which finds the best-fit surface erosion rate for a measured sample assuming steady-state production.
- Zero erosion platform model (no cliff retreat or down-wearing), which finds the best-fit total exposure time from a given sea-level history, with or without surface cover (e.g. beach).
- Down-wearing platform model (no cliff retreat), which finds the best-fit total exposure time and down-wearing rate from a given sea-level history, with or without surface cover (e.g. beach).
- Cliff retreat platform model, which finds the best-fit total exposure time and cliff retreat rate (accelerating, decelerating or constant), from a given sea-level history, with or without surface cover (e.g. beach). The down-wearing rate can either be specified using a present-day rate, or be a free parameter.

A misfit is calculated between predicted and measured concentrations to determine the best-fit. The measurements used for this misfit calculation can be modified (using all, or the min or max) depending on the hypothesis being tested.

Nuclide concentrations can also be calculated and plotted for specific parameters (total exposure time, cliff retreat rate, down-wearing rate, surface cover).

Inputs required to run the model (also see example input files):
- Sample data (.xlsx), with columns: 
    1. Sample name
    2. Latitude (decimal degrees)
    3. Longitude (decimal degrees)
    4. Elevation (m asl)
    5. Pressure (hPa) (zero if not known)
    6. Distance along transect from cliff (m)
    7. Sample thickness (cm)
    8. Bulk density (g/cm^3)
    9. Shielding factor for terrain (unitless)
    10. Sample 10-Be concentration, mean (atoms of 10-Be/g)
    11. Sample 10-Be concentration, 1 sigma uncertainty (atoms of 10-Be/g)
    12. Year the sample was collected (calendar year)
- Topographic shielding across the platform (.txt), with columns:
    1. Distance from cliff (m)
    2. Shielding factor for terrain (unitless)
- Platform elevation profile (.txt), with columns:
    1. Distance from cliff (m)
    2. Elevation (m AOD)
- Relative sea-level history (.txt), with columns:
    1. Time (years before present)
    2. Sea level (m relative to present AOD)
- Tides (.txt), with columns:
    1. Tidal fequency 
    2. Tidal duration (%)
    3. Central elevation of 10 cm vertical bins (m AOD)
- Tidal benchmarks (highest astronomical tide, mean high water level of spring tides, mean high water level of neap tides, mean low water level of spring tides, mean low water level of neap tides.
- Nuclide concentration of inheritance sample (mean and 1 sigma uncertainty), if exists.

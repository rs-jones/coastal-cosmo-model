% COSMOGENIC NUCLIDE MODEL FOR COASTAL EROSION AND SEA LEVEL SCENARIOS
%
% Written by Richard Selwyn Jones (richard.s.jones@monash.edu), using
% modifications of code from Swirad et al. (2020, Nat. Comms.; and refs. 
% therein) for calculating time-dependent topographic and water shielding, 
% and cliff retreat, and CRONUScalc (Marrero et al., 2016, Quat. 
% Geochronology; and refs. therein) for calculating cosmogenic nuclide 
% production. The model uses a MATLAB optimised solver to find the best-fit
% variables for a given case.
% 
% The model consists of various sub-models for different purposes, and to 
% test hypotheses of increasing complexity:
% - Cosmogenic inheritance model, which finds the best-fit surface erosion 
%     rate for a measured sample assuming steady-state production.
% - Zero platform erosion model (no cliff retreat or down-wearing), which 
%     finds the best-fit total exposure time from a given sea-level 
%     history, with or without surface cover (e.g. beach).
% - Down-wearing platform model (no cliff retreat), which finds the 
%     best-fit total exposure time and down-wearing rate from a given 
%     sea-level history, with or without surface cover (e.g. beach).
% - Cliff retreat platform model, which finds the best-fit total exposure 
%     time and cliff retreat rate (accelerating, decelerating or constant),
%     from a given sea-level history, with or without surface cover (e.g. 
%     beach). The down-wearing rate can either be specified using a
%     present-day rate, or be a free parameter.
% 
% A misfit is calculated between predicted and measured concentrations to
% determine the best-fit. The measurements used for this misfit calculation
% can be modified (using all, or the min or max) depending on the
% hypothesis being tested.
% 
% Nuclide concentrations can also be calculated and plotted for specific
% parameters (total exposure time, cliff retreat rate, down-wearing rate,
% surface cover).
% 
% Inputs required to run the model:
% - Sample data (.xlsx), with columns: 
%     1. Sample name
%     2. Latitude (decimal degrees)
%     3. Longitude (decimal degrees)
%     4. Elevation (m asl)
%     5. Pressure (hPa) (zero if not known)
%     6. Distance along transect from cliff (m)
%     7. Sample thickness (cm)
%     8. Bulk density (g/cm^3)
%     9. Shielding factor for terrain (unitless)
%     10. Sample 10-Be concentration, mean (atoms of 10-Be/g)
%     11. Sample 10-Be concentration, 1 sigma uncertainty (atoms of 10-Be/g)
%     12. Year the sample was collected (calendar year)
% - Topographic shielding across the platform (.txt), with columns:
%     1. Distance from cliff (m)
%     2. Shielding factor for terrain (unitless)
% - Platform elevation profile (.txt.), with columns:
%     1. Distance from cliff (m)
%     2. Elevation (m AOD)
% - Relative sea-level history (.txt), with columns:
%     1. Time (years before present)
%     2. Sea level (m relative to present AOD)
% - Tides (.txt), with columns:
%     1. Tidal fequency 
%     2. Tidal duration (%)
%     3. Central elevation of 10 cm vertical bins (m AOD)
% - Tidal benchmarks (highest astronomical tide, mean high water level of 
%     spring tides, mean high water level of neap tides, mean low water 
%     level of spring tides, mean low water level of neap tides.
% - Nuclide concentration of inheritance sample (mean and 1 sigma
%     uncertainty), if exists.
%
%
%%

clear; close all; % Start fresh
addpath(genpath(pwd))


% SPECIFY INPUTS

inh_meas.conc_mean = 1304; % Nuclide concentration of inheritance sample (mean, atoms/g)
inh_meas.conc_uncert = 268; % Nuclide concentration of inheritance sample (1 sigma, atoms/g)
sample_data = get_data('SampleData_Swirad.xlsx',inh_meas); % Sample information
inputs.scaling_model = 'LSDn'; % Nuclide scaling model - e.g. 'ST','LM','LSD','LSDn'
inputs.topo = load('toposhield.txt'); % Topographic shielding across platform
inputs.profile = load('profile.txt'); % Platform elevation profile
inputs.sealevel = load('sealevel.txt'); % Relative sea level through time
inputs.tides = load('tides.txt'); % Tidal data
inputs.HAT = 3.2; % Highest astronomical tide
inputs.MHWS = 2.59; % Mean high water level of spring tides
inputs.MHWN = 1.5; % Mean high water level of neap tides
inputs.MLWN = -0.75; % Mean low water level of neap tides
inputs.MLWS = -2.01; % Mean low water level of spring tides


%% RUN MODEL TO EVALUATE PLATFORM INHERITANCE

% Specify inheritance sample values
inh_meas.depth = 100; % Depth of inheritance sample from cliff top (m)

% Specify model values
erorate_initial = 0.01; % Approx. cliff surface erosion rates (mm yr; used to initiate the optimisation)


% Plot measured concentration of inheritance sample
close all; fig_h = plot_concs_depth(inh_meas);

% Find best-fit scenario from optimisation solver
bestfit = run_InheritanceSS(inputs,erorate_initial,sample_data,inh_meas,fig_h);


%% RUN MODEL FOR ZERO EROSION

% Specify model values
totaltime_initial = 7000; % Approx. total time (years; used to initiate the optimisation)
misfit_meas = 'max'; % Concentration measurements to use for calculating the misfit ('all','min','max')


% Plot platform and measured concentrations
close all; fig_h = plot_platform_concs(inputs,sample_data);

% Find best-fit scenario from optimisation solver
bestfit = run_ZeroErosion(inputs,totaltime_initial,[],sample_data,misfit_meas,fig_h);


%% RUN MODEL FOR ZERO EROSION + SURFACE COVER

% Specify model values
cover.density = 1.6; % Density of surface cover when platform above RSL+HAT (g/cm^3)
cover.depth_initial = 1; % Approx. max depth of beach berm (metres; used to initiate the optimisation)
totaltime_initial = 7000; % Approx. total time (years; used to initiate the optimisation)
misfit_meas = 'max'; % Concentration measurements to use for calculating the misfit ('all','min','max','minmax')


% Plot platform and measured concentrations
close all; fig_h = plot_platform_concs(inputs,sample_data);

% Find best-fit scenario from optimisation solver
bestfit = run_ZeroErosion(inputs,totaltime_initial,cover,sample_data,misfit_meas,fig_h);


%% RUN MODEL FOR DOWN-WEARING (+ SURFACE COVER)

% Specify model values
rateDW.present = 0.02; % Present-day down-wearing rate (mm yr)
rateDW.change_initial = 0.5; % Past rate was faster (>1), slower (0>1) or the same/constant (1) (multiplier relative to present; used to initiate the optimisation)
totaltime_initial = 7000; % Approx. total time (years; used to initiate the optimisation)
cover.density = 1.6; % Density of surface cover when platform above RSL+HAT (g/cm^3)
cover.depth_initial = []; % Approx. max depth of beach berm (metres; used to initiate the optimisation) - leave empty for no cover
misfit_meas = 'minmax'; % Concentration measurements to use for calculating the misfit ('all','min','max','minmax')


% Plot platform and measured concentrations
close all; fig_h = plot_platform_concs(inputs,sample_data);

% Find best-fit scenario from optimisation solver
bestfit = run_ErosionDW(inputs,totaltime_initial,rateDW,cover,sample_data,misfit_meas,fig_h);


%% RUN MODEL FOR CLIFF RETREAT (USING DOWN-WEARING AT PRESENT RATE)

% Specify model values
rateCR.present = 0.1; % Present-day cliff retreat rate (m yr)
rateCR.change_initial = 1; % Past rate was faster (>1), slower (0>1) or the same/constant (1) (multiplier relative to present; used to initiate the optimisation)
totaltime_initial = 7000; % Approx. total time (years; used to initiate the optimisation)
rateDW.present = 0.02; % Present-day down-wearing rate (mm yr)
misfit_meas = 'minmax'; % Concentration measurements to use for calculating the misfit ('all','min','max','minmax')


% Plot platform and measured concentrations
close all; fig_h = plot_platform_concs(inputs,sample_data);

% Find best-fit scenario from optimisation solver
bestfit = run_ErosionCRDW(inputs,totaltime_initial,rateCR,rateDW,[],sample_data,misfit_meas,fig_h);


%% RUN MODEL FOR CLIFF RETREAT (+ DOWN-WEARING + SURFACE COVER)

% Specify model values
rateCR.present = 0.025; % Present-day cliff retreat rate (m yr)
rateCR.change_initial = 1.5; % Past rate was faster (>1), slower (0>1) or the same/constant (1) (multiplier relative to present; used to initiate the optimisation)
rateDW.present = 0.02; % Present-day down-wearing rate (mm yr)
rateDW.change_initial = 0; % Past rate was faster (>1), slower (0>1) or the same/constant (1) (multiplier relative to present; used to initiate the optimisation) - leave empty to only use constant present-day down-wearing rate
totaltime_initial = 7000; % Approx. total time (years; used to initiate the optimisation)
cover.density = 1.6; % Density of surface cover when platform above RSL+HAT (g/cm^3)
cover.depth_initial = 3; % Approx. max depth of beach berm (metres; used to initiate the optimisation)
misfit_meas = 'minmax'; % Concentration measurements to use for calculating the misfit ('all','min','max','minmax')


% Plot platform and measured concentrations
close all; fig_h = plot_platform_concs(inputs,sample_data);

% Find best-fit scenario from optimisation solver
bestfit = run_ErosionCRDW(inputs,totaltime_initial,rateCR,rateDW,cover,sample_data,misfit_meas,fig_h);


%% CALCULATE CONCENTRATIONS FOR A SPECIFIC SCENARIO

% Specify scenario values
totaltime = 6800; % Approx. total time (years; used to initiate the optimisation)
rateCR.present = 0.025; % Present-day cliff retreat rate (m yr)
rateCR.change = 2; % Past rate was faster (>1), slower (0>1) or the same/constant (1) (multiplier relative to present; used to initiate the optimisation)
rateDW.present = 0.02; % Present-day down-wearing rate (mm yr)
rateDW.change = 4; % Past rate was faster (>1), slower (0>1) or the same/constant (1) (multiplier relative to present; used to initiate the optimisation) - leave empty to only use constant present-day down-wearing rate
cover.density = 1.6; % Density of surface cover when platform above RSL+HAT (g/cm^3)
cover.depth = 2; % Approx. max depth of beach berm (metres; used to initiate the optimisation)


% Plot platform and measured concentrations
close all;
fig_h = plot_platform_concs(inputs,sample_data);

% Calculate and plot scenario concentrations
calc_Scenario(inputs,totaltime,rateCR,rateDW,cover,sample_data,fig_h);


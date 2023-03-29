
function misfit = fit_ErosionCRDW(X,inputs,sample_data,misfit_meas,plot_fig)

% UNPACK
cliff_ret_multiplier = X(1);
total_time = X(2);
if total_time < 0
    misfit = 1000;
    return
end
time = 1:total_time;
if numel(X)>2
    ero_multiplier = X(3);
else
    ero_multiplier = 1;
end
if numel(X)>3
    berm_h = X(4);
    if berm_h < 0
        berm_h = 0;
    end
end


% GET PLATFORM ELEVATIONS AND RELATIVE SEA LEVEL
platform_dist = inputs.profile(:,1)';
elev_pres = inputs.profile(:,2)';
rsl = flipud(inputs.sealevel);
rsl = interp1(rsl(:,1),rsl(:,2),time'); % Interpolate for model time
rsl(isnan(rsl)) = rsl(find(~isnan(rsl),1,'last')); % Use earliest value if input shorter than total time
rsl_plusHAT = rsl+inputs.HAT; % Relative elevation of highest astronomical tide


% CALCULATE CLIFF RETREAT RATE AND EXPOSURE TIMES

% Calculate retreat rate through time
starting_rate = inputs.rateCR_pres_myr * cliff_ret_multiplier;
if starting_rate<0
    starting_rate = 0;
end
cliff_retreat = linspace(inputs.rateCR_pres_myr,starting_rate,total_time);
platform_aboveSL = min(elev_pres) > rsl_plusHAT; % Find when base of platform is above sea level (HAT)
cliff_retreat(platform_aboveSL) = 0; % Give zero retreat

% Create a table of relative cliff positions through time
cliff_pos = zeros(1,length(time));
for n = 1:length(time)-1
    cliff_pos(n+1) = cliff_pos(n) + cliff_retreat(n);
end
cliff_pos = round(cliff_pos,4);

% Create a table of across-shore exposure time
expo = zeros(1,length(inputs.profile(:,1)));
for w = 1:length(inputs.profile(:,1)) % For each point along platform
    % If starting cliff position is inland of this platform position,
    % give maximum exposure
    if cliff_pos(end)<(inputs.profile(w,1))
        expo(w) = max(time);
    
    % Otherwise, use time when this cliff position is greater/equal to
    % this platform position
    else
        this_pos_idx = find(cliff_pos>=(inputs.profile(w,1)),1);
        this_pos_t = time(this_pos_idx);
        expo(w) = this_pos_t;
    end    
end


% CALCULATE TOPOGRAPHIC SHIELDING
sTopo = shielding_topo(inputs.topo,inputs.profile,cliff_pos,expo);


% CALCULATE WATER SHIELDING
sWater = zeros(1,length(platform_dist));
sw_t = zeros(length(time),length(platform_dist));
elev_belowHAT = zeros(length(time),length(platform_dist));

% Calculate through time
for n = 1:length(time)
    elev_relRSL = elev_pres-rsl(n); % Elevation relative to sea level at time n
    % Calculate for each elevation across platform
    for w = 1:length(platform_dist)
        if expo(w)>=n % Do only if this point along platform is younger than the exposure time
            elev_belowHAT(n,w) = elev_relRSL(w) < inputs.HAT;
            % Do only if this point along the platform is below HAT
            if elev_belowHAT(n,w)
                sw_t(n,w) = shielding_water(elev_relRSL(w),inputs.tides,inputs.HAT);
            else
                sw_t(n,w) = 1;
            end
        else
            sw_t(n,w) = NaN;
        end
    end
end
for w = 1:length(platform_dist)
    sWater(w) = sum(sw_t(:,w),'omitnan')/expo(w); % Average over exposure time
end


% CALCULATE SURFACE COVER SHIELDING
if numel(X)>3
    
    % Calculate shielding through time and across the platform
    [sc_t,cover_depth] = shielding_cover(inputs,sample_data,time,platform_dist,elev_pres,rsl,berm_h);
    
    sCover = zeros(1,length(platform_dist));
    for w = 1:length(platform_dist)
        sCover(w) = (sum(sc_t(:,w)))/total_time; % Average over time
    end
    sCover(sCover>1)=1; % Limit shielding to 1
end


% DETERMINE DOWN-WEARING RATE
    
% Generate rate through time
starting_rate = inputs.rateDW_pres_gcm2yr * ero_multiplier;
ero_rate_t = linspace(inputs.rateDW_pres_gcm2yr,starting_rate,total_time);

% Calculate through time
ero_rate_t_platform = zeros(length(time),length(platform_dist));
ero_rate = zeros(1,length(platform_dist));
for n = 1:length(time)
    % Calculate across platform
    for w = 1:length(platform_dist)
        if expo(w)>=n % Do only if this point along platform is younger than the exposure time
            ero_rate_t_platform(n,w) = ero_rate_t(n);
            ero_rate_t_platform(~elev_belowHAT(n,w),w) = 0; % Remove down-wearing when platform is above RSL+HAT
            if numel(X)>3
                ero_rate_t_platform(sc_t(n,w)<1,w) = 0; % Remove down-wearing when platform is has surface cover
            end       
        else
            ero_rate_t_platform(n,w) = NaN;
        end
    end
end
for w = 1:length(platform_dist)
    ero_rate(w) = sum(ero_rate_t_platform(:,w),'omitnan')/expo(w);
end
ero_rate(ero_rate<0)=0; % Remove any negative erosion values


% PLOT SCENARIO VARIABLES
if 1
    if numel(X)<4
        sCover = [];
    end
   plot_variables(time,rsl,inputs,cliff_retreat,expo,(ero_rate./sample_data.mean_rho).*10,sTopo,sWater,sCover,elev_belowHAT);
end


% CALCULATE CONCENTRATIONS ALONG PLATFORM

% Get parameters
pars.pp = sample_data.pp;
pars.sf = sample_data.sf{1};
pars.cp = sample_data.cp{1};
pars.l = sample_data.pp.lambda10Be;
pars.top_z_gcm2 = sample_data.top_z_gcm2;
pars.bottom_z_gcm2 = sample_data.bottom_z_gcm2;


% Calculate
if numel(X)<3
    sCover = ones(size(sTopo)); % Assign no cover
end
pred_concs = calc_conc(pars,expo,ero_rate,sTopo,sWater,sCover);
concs_nonan = pred_concs; concs_nonan(isnan(concs_nonan)) = 0;


% CALCULATE MISFIT

if strcmpi(misfit_meas,'all')
    for meas = 1:length(sample_data.measured_inh(:,1))
        this_measured_idx = round(sample_data.measured_inh(meas,1))==inputs.profile(:,1)'; % Find sample along platform
        sample_misfits(meas) = ((concs_nonan(this_measured_idx) - sample_data.measured_inh(meas,2))./sample_data.measured_inh(meas,3)) .^2;
    end
    misfit = mean(sample_misfits,2);
elseif strcmpi(misfit_meas,'min')
    min_log = min(sample_data.measured_inh(:,2)) == sample_data.measured_inh(:,2);
    min_measured_idx = round(sample_data.measured_inh(min_log,1))==inputs.profile(:,1)'; % Find sample along platform
    misfit = ((concs_nonan(min_measured_idx) - sample_data.measured_inh(min_log,2))./sample_data.measured_inh(min_log,3)) .^2;
elseif strcmpi(misfit_meas,'max')
    max_log = max(sample_data.measured_inh(:,2)) == sample_data.measured_inh(:,2);
    max_measured_idx = round(sample_data.measured_inh(max_log,1))==inputs.profile(:,1)'; % Find sample along platform
    misfit = ((concs_nonan(max_measured_idx) - sample_data.measured_inh(max_log,2))./sample_data.measured_inh(max_log,3)) .^2;
elseif strcmpi(misfit_meas,'gradient') || strcmpi(misfit_meas,'minmax')
    for meas = 1:length(sample_data.measured_inh(:,1))
        this_measured_idx = round(sample_data.measured_inh(meas,1))==inputs.profile(:,1)'; % Find sample along platform
        sample_misfits(meas) = ((concs_nonan(this_measured_idx) - sample_data.measured_inh(meas,2))./sample_data.measured_inh(meas,3)) .^2;
    end
    misfit_mean = mean(sample_misfits,2);
    min_log = min(sample_data.measured_inh(:,2)) == sample_data.measured_inh(:,2);
    misfit_min = sample_misfits(min_log);
    max_log = max(sample_data.measured_inh(:,2)) == sample_data.measured_inh(:,2);
    misfit_max = sample_misfits(max_log);
    misfit_minmax = mean([misfit_min,misfit_max]);
    misfit = mean([misfit_mean,misfit_minmax]);
end


% PLOT CONCENTRATIONS FOR SCENARIO
if ~isempty(plot_fig)
    cliff_pos_idx = find(cliff_retreat==0,1)-1;
    if isempty(cliff_pos_idx) % Cliff position before recent period of retreat
        inputs.cliff_pos = round(cliff_pos(end));
    else
        inputs.cliff_pos = round(cliff_pos(cliff_pos_idx));
    end
    if ~isempty(inputs.rateDW_pres_gcm2yr) || (inputs.rateDW_pres_gcm2yr~=0 && ero_multiplier~=1)
        ero_rate_myr = (ero_rate ./ sample_data.mean_rho) ./ 100;
        inputs.profile_ero_expo = inputs.profile(:,2)' + (ero_rate_myr.*expo); % Platform profile before down-wearing
    end
    if numel(X)>3 && berm_h>0
        inputs.elev_beach = inputs.profile(:,2)' + mean(cover_depth); % Surface profile of (mean) surface cover
    end
    
    figure(plot_fig)
    plot_platform_concs(inputs,sample_data,pred_concs)
    
    drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%% Plot scenario variables %%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_variables(time,rsl,inputs,cliff_retreat,expo,ero_rate,sTopo,sWater,sCover,elev_belowHAT)
        
        figure(10);
        cols = lines(5); % colour-scale for plotting
        
        subplot(2,2,1)
        yyaxis left
        plot(time/1000,rsl,'-b','Linewidth',1);
        title('Sea level & cliff retreat')
        ylabel('Sea level (m AOD)')
        xlabel('Time (kyr BP)')
        yyaxis right
        plot(time/1000,cliff_retreat,'-k','Linewidth',1);
        ylabel('Cliff retreat (m/yr)')
        ax = gca;
        ax.YAxis(1).Color = 'b';
        ax.YAxis(2).Color = 'k';
        
        if isempty(sCover)
            subplot(2,2,2)
            plot(inputs.profile(:,1),sum(elev_belowHAT),'-k','Linewidth',1);
            title('Platform submergence')
            xlabel('Distance from the cliff (m)')
            ylabel('Years below RSL+HAT')
        else
            subplot(2,2,2)
            plot(inputs.profile(:,1),sCover,'color',cols(3,:),'Linewidth',1);
            title('Cumulative cover shielding')
            xlabel('Distance from the cliff (m)')
            ylabel('Cover shielding, sCover')
        end
        
        subplot(2,2,3)
        yyaxis left
        plot(inputs.profile(:,1),sTopo,'-','color',cols(4,:),'Linewidth',1);
        title('Cumulative shielding')
        xlabel('Distance from the cliff (m)')
        ylabel('Topographic shielding, sTopo')
        yyaxis right
        plot(inputs.profile(:,1),sWater,'color',cols(1,:),'Linewidth',1);
        ylabel('Water shielding, sWater')
        ax = gca;
        ax.YAxis(1).Color = cols(4,:);
        ax.YAxis(2).Color = cols(1,:);
        
        subplot(2,2,4)
        ero_rate(1)=NaN;
        yyaxis left
        plot(inputs.profile(:,1),expo/1000,'-','color',cols(5,:),'Linewidth',1);
        ylabel('Exposure time (kyr BP)')
        yyaxis right
        plot(inputs.profile(:,1),ero_rate,'-','Color',cols(3,:),'Linewidth',1);
        xlabel('Distance from the cliff (m)')
        ylabel('Down-wearing (mm yr)')
        title('Exposure time & erosion rates')
        ax = gca;
        ax.YAxis(1).Color = cols(5,:);
        ax.YAxis(2).Color = cols(3,:);
         
        drawnow;
    end

end

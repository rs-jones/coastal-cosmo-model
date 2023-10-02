
function misfit = fit_ErosionDW(X,inputs,sample_data,misfit_meas,plot_fig)

% UNPACK
ero_multiplier = X(1);
total_time = X(2);
if total_time < 0
    misfit = 1000;
    return
end
time = 1:total_time;
if numel(X)>2
    berm_h = X(3);
    if berm_h < 0
        berm_h = 0;
    elseif berm_h > 5
        misfit = 1000;
        return
    end
end


% GET PLATFORM ELEVATIONS AND RELATIVE SEA LEVEL
platform_dist = inputs.profile(:,1)';
elev_pres = inputs.profile(:,2)';
rsl = flipud(inputs.sealevel);
rsl = interp1(rsl(:,1),rsl(:,2),time'); % Interpolate for model time
rsl(isnan(rsl)) = rsl(find(~isnan(rsl),1,'last')); % Use earliest value if input shorter than total time


% CALCULATE EXPOSURE TIMES
expo = zeros(1,length(platform_dist)) + total_time;


% CALCULATE TOPOGRAPHIC SHIELDING
sTopo = inputs.topo(:,2)';


% CALCULATE WATER SHIELDING
sWater = zeros(1,length(platform_dist));
sw_t = zeros(length(time),length(platform_dist));
elev_belowHAT = zeros(length(time),length(platform_dist));
elev_belowMLWN = zeros(length(time),length(platform_dist));

% Calculate through time
for n = 1:length(time)
    elev_relRSL = elev_pres-rsl(n); % Elevation relative to sea level at time n
    % Calculate for each elevation across platform
    for w = 1:length(platform_dist)
        elev_belowHAT(n,w) = elev_relRSL(w) < inputs.HAT;
        % Do only if this point along the platform is below HAT
        if logical(elev_belowHAT(n,w))
            sw_t(n,w) = shielding_water(elev_relRSL(w),inputs.tides,inputs.HAT);
        else
            sw_t(n,w) = 1;
        end
        elev_belowMLWN(n,w) = elev_relRSL(w) < inputs.MLWN; % For use later
    end
end
for w = 1:length(platform_dist)
    sWater(w) = (sum(sw_t(:,w)))/expo(w);
end


% CALCULATE SURFACE COVER SHIELDING
if numel(X)>2
    
    % Calculate shielding through time and across the platform
    [sc_t,cover_depth] = shielding_cover(inputs,sample_data,time,platform_dist,elev_pres,rsl,berm_h);
    
    sCover = zeros(1,length(platform_dist));
    for w = 1:length(platform_dist)
        sCover(w) = (sum(sc_t(:,w)))/expo(w);
    end
    sCover(sCover>1)=1; % Limit shielding to 1
end


% DETERMINE DOWN-WEARING RATE

% Calculate rate through time
if ero_multiplier < 0
    starting_rate = inputs.rateDW_pres_gcm2yr * abs(ero_multiplier);
elseif ero_multiplier > 0
    starting_rate = inputs.rateDW_pres_gcm2yr * 1+ero_multiplier;
else
    starting_rate = inputs.rateDW_pres_gcm2yr * 1;
end
starting_rate = round(starting_rate,3);
ero_rate_t = linspace(inputs.rateDW_pres_gcm2yr,starting_rate,total_time);

% Calculate through time
ero_rate_t_platform = zeros(length(time),length(platform_dist));
ero_rate = zeros(1,length(platform_dist));
for n = 1:length(time)
    % Calculate across platform
    for w = 1:length(platform_dist)
        this_ero_rate = ero_rate_t(n);
        if logical(elev_belowMLWN(n,w))
            this_ero_rate = 0; % Remove down-wearing when platform is below RSL+MLWN
        end
        if numel(X)>2 && berm_h>0 && sc_t(n,w)<1
            this_ero_rate = 0; % Remove down-wearing when platform is has surface cover
        end
        ero_rate_t_platform(n,w) = this_ero_rate;
    end
end
for w = 1:length(platform_dist)
    ero_rate(w) = sum(ero_rate_t_platform(:,w))/expo(w);
end
ero_rate(ero_rate<0)=0; % Remove any negative erosion values


% PLOT SCENARIO VARIABLES
if 1
    if numel(X)<3
        sCover = [];
    end
   plot_variables(time,rsl,inputs,(ero_rate./sample_data.mean_rho).*10,sTopo,sWater,sCover,elev_belowHAT);
end


% CALCULATE CONCENTRATIONS ALONG PLATFORM

% Get parameters
pars.pp = sample_data.pp;
pars.sf10 = sample_data.sf{1};
pars.cp10 = sample_data.cp{1};
pars.l10 = sample_data.pp.lambda10Be;
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
    ero_rate_myr = (ero_rate ./ sample_data.mean_rho) ./ 100;
    inputs.profile_ero_expo = inputs.profile(:,2)' + (ero_rate_myr.*expo); % Platform profile before down-wearing
    if numel(X)>2 && berm_h>0
        inputs.elev_beach = inputs.profile(:,2)' + mean(cover_depth); % Surface profile of (mean) surface cover
    end
    
    figure(plot_fig)
    plot_platform_concs(inputs,sample_data,pred_concs)
    
    drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%% Plot scenario variables %%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_variables(time,rsl,inputs,ero_rate,sTopo,sWater,sCover,elev_belowHAT)
        
        figure(10);
        cols = lines(6); % Colour scale for plotting
        
        subplot(2,2,1)
        plot(time/1000,rsl,'-b','Linewidth',1);
        title('Relative sea level')
        ylabel('Sea level (m AOD)')
        xlabel('Time (kyr BP)')
        
        subplot(2,2,2)
        plot(inputs.profile(:,1),sum(elev_belowHAT),'-','Color',cols(6,:),'Linewidth',1);
        title('Platform submergence')
        xlabel('Distance from the cliff (m)')
        ylabel('Years below RSL+HAT')
        ax = gca; ax.YAxis(1).Color = cols(6,:);
        
        subplot(2,2,3)
        plot(inputs.profile(:,1),sTopo,'-','color',cols(4,:),'Linewidth',1); hold on;
        plot(inputs.profile(:,1),sWater,'-','color',cols(1,:),'Linewidth',1);
        if ~isempty(sCover)
            plot(inputs.profile(:,1),sCover,'-','color',cols(3,:),'Linewidth',1);
        end
        hold off;
        title('Cumulative shielding')
        xlabel('Distance from the cliff (m)')
        ylabel('Shielding factor')
        if isempty(sCover)
            legend('Topographic','Water','Location','Best');
        else
            legend('Topographic','Water','Cover','Location','Best');
        end
        
        subplot(2,2,4)
        plot(inputs.profile(:,1),ero_rate,'-','Color',cols(2,:),'Linewidth',1);
        title('Cumulative down-wearing')
        xlabel('Distance from the cliff (m)')
        ylabel('Rate (mm yr)')
        ax = gca; ax.YAxis(1).Color = cols(2,:);
                
        drawnow;
    end

end


function misfit = fit_ZeroErosion(X,inputs,sample_data,misfit_meas,plot_fig)

% UNPACK
total_time = X(1);
time = 1:total_time;
if numel(X)>1
    berm_h = X(2);
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


% CALCULATE EXPOSURE TIMES
expo = zeros(1,length(platform_dist)) + total_time;


% CALCULATE TOPOGRAPHIC SHIELDING
sTopo = inputs.topo(:,2)';


% CALCULATE WATER SHIELDING
sWater = zeros(1,length(platform_dist));
sw_t = zeros(length(time),length(platform_dist));
elev_belowHAT = zeros(length(time),length(platform_dist));

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
    end
end
for w = 1:length(platform_dist)
    sWater(w) = (sum(sw_t(:,w)))/expo(w);
end


% CALCULATE SURFACE COVER SHIELDING
if numel(X)>1
    
    % Calculate shielding through time and across the platform
    [sc_t,cover_depth] = shielding_cover(inputs,sample_data,time,platform_dist,elev_pres,rsl,berm_h);
    
    sCover = zeros(1,length(platform_dist));
    for w = 1:length(platform_dist)
        sCover(w) = (sum(sc_t(:,w)))/expo(w);
    end
    sCover(sCover>1)=1; % Limit shielding to 1
end


% PLOT SCENARIO VARIABLES
if 1
    if numel(X)==1
        sCover = [];
    end
   plot_variables(time,rsl,inputs,sTopo,sWater,sCover,elev_belowHAT);
end


% CALCULATE CONCENTRATIONS ALONG PLATFORM

% Get parameters
[~,max_idx] = max(sample_data.CC(:,3)); % Use the highest elevation sample (so not below sea level)
pars.pp = sample_data.pp;
pars.sf10 = sample_data.sf1026{max_idx};
pars.cp10 = sample_data.cp1026{max_idx};
pars.l10 = sample_data.pp.lambda10Be;
pars.top_z_gcm2 = mean(sample_data.CC(:,14)) * mean(sample_data.CC(:,6)); % Use mean of samples (convert to g cm2)
pars.bottom_z_gcm2 = mean(sample_data.CC(:,5)) * mean(sample_data.CC(:,6)); % Use mean of samples (convert to g cm2)

% Calculate
if numel(X)==1
    sCover = ones(size(sTopo)); % Assign no cover
end
ero_gcm2yr = zeros(size(sTopo)); % Assign no erosion
pred_concs = calc_conc(pars,expo,ero_gcm2yr,sTopo,sWater,sCover);
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
    if numel(X)>1 && berm_h>0
        inputs.elev_beach = inputs.profile(:,2)' + mean(cover_depth); % Surface profile of (mean) surface cover
    end
    
    figure(plot_fig)
    plot_platform_concs(inputs,sample_data,pred_concs)
    
    drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%% Plot scenario variables %%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_variables(time,rsl,inputs,sTopo,sWater,sCover,elev_belowHAT)
        
        figure(10);
        cols = lines(6); % Colour scale for plotting
        
        ax1 = subplot(2,2,1);
        plot(time/1000,rsl,'-b','Linewidth',1);
        title('Relative sea level')
        ylabel('Sea level (m AOD)')
        xlabel('Time (kyr BP)')
        
        ax3 = subplot(2,2,3);
        plot(inputs.profile(:,1),sTopo,'-','color',cols(4,:),'Linewidth',1);
        title('Cumulative topographic shielding')
        xlabel('Distance from the cliff (m)')
        ylabel('Shielding, sTopo')
        
        if isempty(sCover)
            
            ax2 = subplot(2,2,2);
            plot(inputs.profile(:,1),sum(elev_belowHAT),'-','Color',cols(6,:),'Linewidth',1);
            title('Platform submergence')
            xlabel('Distance from the cliff (m)')
            ylabel('Years below RSL+HAT')
            
            ax4 = subplot(2,2,4);
            plot(inputs.profile(:,1),sWater,'color',cols(1,:),'Linewidth',1);
            title('Cumulative water shielding')
            xlabel('Distance from the cliff (m)')
            ylabel('Shielding, sWater')
        
        else
            
            ax2 = subplot(2,2,2);
            yyaxis left
            plot(inputs.profile(:,1),sum(elev_belowHAT),'-','Color',cols(6,:),'Linewidth',1);
            ylabel('Years below RSL+HAT')
            yyaxis right
            plot(inputs.profile(:,1),sWater,'color',cols(1,:),'Linewidth',1);
            ylabel('Shielding, sWater')
            xlabel('Distance from the cliff (m)')
            title('Submergence and shielding')
            ax2.YAxis(1).Color = cols(6,:);
            ax2.YAxis(2).Color = cols(1,:);
            
            ax4 = subplot(2,2,4);
            plot(inputs.profile(:,1),sCover,'color',cols(3,:),'Linewidth',1);
            title('Cumulative cover shielding')
            xlabel('Distance from the cliff (m)')
            ylabel('Shielding, sCover')
            
        end
                
        drawnow;
    end

end

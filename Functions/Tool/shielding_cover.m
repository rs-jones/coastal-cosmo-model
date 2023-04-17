%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cover shielding %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sc,cover_depth] = shielding_cover(inputs,sample_data,time,platform_dist,elev_pres,rsl,berm_h)

attenuation = mean(sample_data.CC(:,13));
%A = 0.12;  shape = 2/3; % Bruun profile parameters (after Hurst et al., 2017)

% Determine approx. present beach slope based on tidal range
[x_high,~] = intersections(platform_dist,elev_pres,platform_dist,zeros(size(elev_pres))+inputs.HAT);
[x_low,~] = intersections(platform_dist,elev_pres,platform_dist,zeros(size(elev_pres))+inputs.MHWN);
beach_slope_dist = round(x_low(1))-round(x_high(1)); % Horizontal distance between HAT and MHWN at present

% Assign through time and across platform
sc = zeros(length(time),length(platform_dist));
elev_beach = zeros(length(time),length(platform_dist));
cover_depth = zeros(length(time),length(platform_dist));
for n = 1:length(time)
    elev_relRSL = elev_pres-rsl(n); % Elevation relative to sea level at time n
    elev_beach(n,:) = elev_relRSL + berm_h; % Berm elevation
    
    % Calculate berm slope if highest part of platform is above MHWN
    % and lowest part is below HAT
    if max(elev_relRSL) > rsl(n)+inputs.MHWN && min(elev_relRSL) < rsl(n)+inputs.HAT
        [x_high,~] = intersections(platform_dist,elev_relRSL,platform_dist,zeros(size(elev_relRSL))+rsl(n)+inputs.HAT);
        if isempty(x_high)
            x_high = 0; % Start slope from cliff if HAT above platform
        end
        elev_berm = elev_relRSL(round(x_high(1))+1) + berm_h;
        beach_slope_x = round(x_high(1)):round(x_high(1))+beach_slope_dist; % Slope from RSL+HAT
        slope_on_plat = length(platform_dist) >= beach_slope_x+1;
        slope_idx = beach_slope_x(slope_on_plat)+1; % Indices of beach slope on the platform
        slope_profile = linspace(elev_berm,elev_relRSL(round(x_low(1))+1),numel(beach_slope_x)); % Linear berm slope profile
        elev_beach(n,slope_idx) = slope_profile(slope_on_plat);
        %elev_beach(n,slope_idx) = elev_berm - A .* beach_slope_x .^shape; % Berm slope profile after Bruun (1954)
    end
    
    % Remove beach below MHWN
    beach_belowMHWN = elev_relRSL < rsl(n)+inputs.MHWN;
    if elev_relRSL(end) < rsl(n)+inputs.MHWN % Only if end of platform is under water
        elev_beach(n,beach_belowMHWN) = elev_relRSL(beach_belowMHWN);
    end
    
    % Plot cover through time
    if 0
        figure(12);
        plot(platform_dist,elev_beach(n,:),'-r'); hold on;
        plot(platform_dist,elev_relRSL,'-k'); hold off;
        title(strcat('yearsBP=',num2str(n)))
        xlabel('Distance from the cliff (m)')
        ylabel('Elevation (m above RSL)')
    end
    
    % Get cover shielding across platform
    for w = 1:length(platform_dist)
        cover_depth(n,w) = elev_beach(n,w)-elev_relRSL(w);
        if cover_depth(n,w) > 0
            sc(n,w) = get_cov_shield('manual',cover_depth(n,w)*100,attenuation,inputs.cover_density);
        else
            sc(n,w) = 1;
        end
    end
end

end

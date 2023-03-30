function calc_Scenario(inputs,total_time,rateCR,rateDW,cover,sample_data,plot_fig)
  
  % GET SCENARIO VARIABLES AND PARAMETERS
  fprintf('\nCalculating nuclide parameters and scenario concentrations...\n');
  time = 1:total_time;
  cliff_ret_multiplier = rateCR.change;
  ero_multiplier = rateDW.change;
  berm_h = cover.depth;
  
  % Create sample for parameters
  [~,highest_idx] = max(sample_data.CC(:,3)); % Use the highest elevation sample (so not below sea level)
  par_sample_data.CC = sample_data.CC(highest_idx,:);

  % Find maximum depth based on highest concentration sample
  maxconc_idx = find(max(sample_data.measured_inh(:,2)+sample_data.measured_inh(:,3))==sample_data.measured_inh(:,2)+sample_data.measured_inh(:,3));
  rho = mean(sample_data.CC(:,6));
  maxconc_z_gcm2 = sample_data.CC(maxconc_idx,5) * rho;  
  maxage=8160; % 6* half-life = 8160 ka
  if isempty(rateDW.present)
      sample_data.maxdepth = maxconc_z_gcm2*2 +10000; % sample depth*2 + safety factor
  else
      erosion = rateDW.present/10 *100; % cm/yr *100
      sample_data.maxdepth = maxconc_z_gcm2 + (maxage*1000)*erosion*rho; % sample depth + maxage*erosion(cm/yr)*density
  end

  
  % GET PLATFORM ELEVATIONS AND RELATIVE SEA LEVEL
  platform_dist = inputs.profile(:,1)';
  elev_pres = inputs.profile(:,2)';
  rsl = flipud(inputs.sealevel);
  rsl = interp1(rsl(:,1),rsl(:,2),time'); % Interpolate for model time
  rsl_plusHAT = rsl+inputs.HAT; % relative elevation of highest astronomical tide
  
  
  % CALCULATE CLIFF RETREAT RATE AND EXPOSURE TIMES
  if ~(isempty(rateCR.present) || isempty(rateCR.change))
      
      rateCR_pres_myr = rateCR.present;
      
      % Calculate retreat rate through time
      starting_rate = rateCR_pres_myr * cliff_ret_multiplier;
      if starting_rate<0
          starting_rate = 0;
      end
      cliff_retreat = linspace(rateCR_pres_myr,starting_rate,total_time);
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
              this_pos_idx = find(cliff_pos>=(inputs.profile(w,1)),1); %(w-1),1);
              this_pos_t = time(this_pos_idx);
              expo(w) = this_pos_t;
          end
      end
      
  else
      
      expo = zeros(1,length(platform_dist)) + total_time;
      warning('No cliff retreat.')
  end
  
  
  % CALCULATE TOPOGRAPHIC SHIELDING
  if ~(isempty(rateCR.present) || isempty(rateCR.change))
      sTopo = shielding_topo(inputs.topo,inputs.profile,cliff_pos,expo);
  else
      sTopo = inputs.topo(:,2)';
  end
  
  
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
  if ~isempty(berm_h) && berm_h > 0
      
      % Calculate shielding through time and across the platform
      inputs.cover_density = cover.density;
      [sc_t,cover_depth] = shielding_cover(inputs,sample_data,time,platform_dist,elev_pres,rsl,berm_h);
      
      sCover = zeros(1,length(platform_dist));
      for w = 1:length(platform_dist)
          sCover(w) = (sum(sc_t(:,w)))/total_time; % Average over time
      end
      sCover(sCover>1)=1; % Limit shielding to 1
  else
      sCover = ones(1,length(platform_dist));
      warning('No surface cover.')
  end
  
  
  % DETERMINE DOWN-WEARING RATE
  if ~(isempty(rateDW.present) || isempty(rateDW.change))
      
      rateDW_pres_gcm2yr = (rateDW.present ./10) .* rho;
      
      % Generate rate through time
      starting_rate = rateDW_pres_gcm2yr * ero_multiplier;
      ero_rate_t = linspace(rateDW_pres_gcm2yr,starting_rate,total_time);
      
      % Calculate through time
      ero_rate_t_platform = zeros(length(time),length(platform_dist));
      ero_rate = zeros(1,length(platform_dist));
      for n = 1:length(time)
          % Calculate across platform
          for w = 1:length(platform_dist)
              if expo(w)>=n % Do only if this point along platform is younger than the exposure time
                  ero_rate_t_platform(n,w) = ero_rate_t(n);
                  ero_rate_t_platform(~elev_belowHAT(n,w),w) = 0; % Remove down-wearing when platform is above RSL+HAT
                  if ~isempty(berm_h) && berm_h > 0
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
      
  else
      
      ero_rate = zeros(1,length(platform_dist));
      warning('No down-wearing.')
  end
  
  
  % CALCULATE CONCENTRATIONS ALONG PLATFORM
  
  % Get nuclide parameters
  par_sample_data = get_pars(par_sample_data,inputs.scaling_model);
  pars.pp = par_sample_data.pp;
  pars.sf = par_sample_data.sf1026{1};
  pars.cp = par_sample_data.cp1026{1};
  pars.l = par_sample_data.pp.lambda10Be;
  
  % Get sample thickness
  pars.top_z_gcm2 = mean(sample_data.CC(:,14)) * rho; % Use mean of samples (convert to g cm2)
  pars.bottom_z_gcm2 = mean(sample_data.CC(:,5)) * rho; % Use mean of samples (convert to g cm2)
    
  % Calculate
  pred_concs = calc_conc(pars,expo,ero_rate,sTopo,sWater,sCover);
  
  
  % PLOT CONCENTRATIONS FOR SCENARIO
  if ~isempty(plot_fig)
      
      if ~(isempty(rateCR.present) || isempty(rateCR.change))
          cliff_pos_idx = find(cliff_retreat==0,1)-1;
          if isempty(cliff_pos_idx) % Cliff position before recent period of retreat
              inputs.cliff_pos = round(cliff_pos(end));
          else
              inputs.cliff_pos = round(cliff_pos(cliff_pos_idx));
          end
      end
      if ~(isempty(rateDW.present) || isempty(rateDW.change))
          ero_rate_myr = (ero_rate ./ rho) ./ 100;
          inputs.profile_ero_expo = inputs.profile(:,2)' + (ero_rate_myr.*expo); % Platform profile before down-wearing
      end
      if ~isempty(berm_h) && berm_h > 0
          inputs.elev_beach = inputs.profile(:,2)' + mean(cover_depth); % Surface profile of (mean) surface cover
      end
      
      figure(plot_fig)
      plot_platform_concs(inputs,sample_data,pred_concs)
      
      drawnow;
  end
  
  fprintf('Done.');
  
end

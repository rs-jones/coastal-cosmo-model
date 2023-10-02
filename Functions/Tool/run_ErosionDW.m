function bestfit = run_ErosionDW(inputs,totaltime_initial,rateDW,cover,sample_data,misfit_meas,plot_fig)
  
  % Create sample for parameters
  fprintf('\nDown-wearing \nCalculating nuclide parameters...');
  [~,highest_idx] = max(sample_data.CC(:,3)); % Use the highest elevation sample (so not below sea level)
  par_sample_data.CC = sample_data.CC(highest_idx,:);

  % Find maximum depth based on highest concentration sample
  maxconc_idx = find(max(sample_data.measured_inh(:,2)+sample_data.measured_inh(:,3))==sample_data.measured_inh(:,2)+sample_data.measured_inh(:,3));
  rho = mean(sample_data.CC(:,6));
  maxconc_z_gcm2 = sample_data.CC(maxconc_idx,5) * rho;  
  maxage = 500; % 500 ka (6* half-life = 8160 ka)
  if ~isfield(rateDW,'change_initial') || rateDW.change_initial < 0
      max_rateDW_change = 1;
  else
      max_rateDW_change = 1+rateDW.change_initial;
  end
  max_DWrate = rateDW.present * max_rateDW_change; % Use starting down-wearing rate
  erosion = max_DWrate/10; % cm/yr
  sample_data.maxdepth = maxconc_z_gcm2 + (maxage*1000)*erosion*rho; % sample depth + maxage*erosion(cm/yr)*density

  % Get nuclide parameters, add to sample data
  par_sample_data = get_pars(par_sample_data,inputs.scaling_model);
  sample_data.pp = par_sample_data.pp;
  sample_data.sf = par_sample_data.sf1026;
  sample_data.cp = par_sample_data.cp1026;
  sample_data.mean_rho = rho;
  
  % Get sample thickness
  sample_data.top_z_gcm2 = mean(sample_data.CC(:,14)) * rho; % Use mean of samples (convert to g cm2)
  sample_data.bottom_z_gcm2 = mean(sample_data.CC(:,5)) * rho; % Use mean of samples (convert to g cm2)
  
  
  % Set optimisation rules
  opts = optimset('fminsearch');
  opts = optimset(opts,'TolFun',1e-1,'TolX',1e-1);
  
  
  % Find bestfit scenario and export result
  fprintf('\nFinding best-fit scenario...\n');
  
  inputs.rateDW_pres_gcm2yr = (rateDW.present ./10) .* rho;
  
  if isempty(cover.depth_initial)
      [optX,fmin] = fminsearch(@(X) fit_ErosionDW(X,inputs,sample_data,misfit_meas,plot_fig),[rateDW.change_initial,totaltime_initial],opts);
      
      % Export result
      if (optX(1) < 0.001) && (optX(1) > -0.001)
          optX(1) = 0;
      end
      if optX(2) <0
          optX(2) = 0;
      end
      if optX(1) > 0
          fprintf('\nBestfit down-wearing rate: past rate was %0.1f times faster (decelerating) \nBestfit total time: %.f years \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              optX(1),round(optX(2)),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      elseif optX(1) < 0
          fprintf('\nBestfit down-wearing rate: past rate was %0.1f times slower (accelerating) \nBestfit total time: %.f years \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              abs(optX(1)),round(optX(2)),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      else
          fprintf('\nBestfit down-wearing rate: same as present (%0.3f mm yr) \nBestfit total time: %.f years \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              rateDW.present,round(optX(2)),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      end
  else
      inputs.cover_density = cover.density;
      
      [optX,fmin] = fminsearch(@(X) fit_ErosionDW(X,inputs,sample_data,misfit_meas,plot_fig),[rateDW.change_initial,totaltime_initial,cover.depth_initial],opts);
      
      % Export result      
      if (optX(1) < 0.001) && (optX(1) > -0.001)
          optX(1) = 0;
      end
      if optX(2) <0
          optX(2) = 0;
      end
      if optX(1) > 0
          fprintf('\nBestfit down-wearing rate: past rate was %0.1f times faster (decelerating) \nBestfit total time: %.f years \nBestfit cover depth: %.f m \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              optX(1),round(optX(2)),optX(3),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      elseif optX(1) < 0
          fprintf('\nBestfit down-wearing rate: past rate was %0.1f times slower (accelerating) \nBestfit total time: %.f years \nBestfit cover depth: %.f m \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              abs(optX(1)),round(optX(2)),optX(3),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      else
          fprintf('\nBestfit down-wearing rate: same as present (%0.3f mm yr) \nBestfit total time: %.f years \nBestfit cover depth: %.f m \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              rateDW.present,round(optX(2)),optX(3),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      end
  end
    
  bestfit = optX;
  
end

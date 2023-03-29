function bestfit = run_ErosionCRDW(inputs,totaltime_initial,rateCR,rateDW,cover,sample_data,misfit_meas,plot_fig)
  
  % Create sample for parameters
  fprintf('\nCliff retreat (and down-wearing) \nCalculating nuclide parameters...');
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
  opts = optimset(opts,'TolFun',2,'TolX',2);
  
  
  % Find bestfit scenario and export result
  
  fprintf('\nFinding best-fit scenario...\n');
  
  inputs.rateCR_pres_myr = rateCR.present;
  inputs.rateDW_pres_gcm2yr = (rateDW.present ./10) .* rho;
  if ~isfield(rateDW,'change_initial')
      rateDW.change_initial = [];
  end
  
  % Run model to find bestfit cliff retreat rate and total time
  if isempty(rateDW.change_initial) || isempty(cover)
      [optX,fmin] = fminsearch(@(X) fit_ErosionCRDW(X,inputs,sample_data,misfit_meas,plot_fig),[rateCR.change_initial,totaltime_initial],opts);
      
      % Export result
      if optX(1) <0.005
          optX(1) = 0;
      end
      if optX(2) <0
          optX(2) = 0;
      end
      if optX(1) > 0
          fprintf('\nBestfit cliff retreat rate: past rate was %0.1f times faster (decelerating) \nBestfit total time: %.f years \nusing present-day down-wearing rate of %0.3f mm yr and no surface cover \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              optX(1),round(optX(2)),rateDW.present,fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      elseif optX(1) < 0
          fprintf('\nBestfit cliff retreat rate: past rate was %0.1f times slower (accelerating) \nBestfit total time: %.f years \nusing present-day down-wearing rate of %0.3f mm yr and no surface cover \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              optX(1)+1,round(optX(2)),rateDW.present,fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      else
          fprintf('\nBestfit cliff retreat rate: same as present (%0.3f mm yr) \nBestfit total time: %.f years \nusing present-day down-wearing rate of %0.3f mm yr and no surface cover \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
              rateCR.present,round(optX(2)),rateDW.present,fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      end
  
  % Otherwise, run model to find bestfit cliff retreat rate, total time,
  % down-wearing rate and cover depth
  else
      if isempty(rateDW.change_initial)
          rateDW.change_initial = 0;
      end
      if isempty(cover.depth_initial)
          cover.depth_initial = 0;
      end
      inputs.cover_density = cover.density;
      
      [optX,fmin] = fminsearch(@(X) fit_ErosionCRDW(X,inputs,sample_data,misfit_meas,plot_fig),[rateCR.change_initial,totaltime_initial,rateDW.change_initial,cover.depth_initial],opts);
      
      % Export result
      if optX(1) <0.005
          optX(1) = 0;
      end
      if optX(2) <0
          optX(2) = 0;
      end
      if optX(3) <0.005
          optX(3) = 0;
      end
      opt1_out = string(round(optX(1),2)+1);
      opt2_out = string(round(optX(2)));
      opt3_out = string(round(optX(3),2));
      opt4_out = string(round(optX(4),2));
      time_txt = strcat('Bestfit total time: ',{' '},opt2_out,{' '},'years');
      cover_txt = strcat('Bestfit cover depth:: ',{' '},opt4_out,{' '},'m');
      if optX(1) > 0
          cliff_txt = strcat('Bestfit cliff retreat rate: past rate was',{' '},opt1_out{1},' times faster (decelerating)');
      elseif optX(1) < 0
          cliff_txt = strcat('Bestfit cliff retreat rate: past rate was',{' '},opt1_out{1},' times slower (accelerating)');
      else
          cliff_txt = strcat('Bestfit cliff retreat rate: same as present (',string(rateCR.present),{' '},'mm yr)');
      end
      if optX(3) > 0
          dw_txt = strcat('Bestfit down-wearing rate: past rate was',{' '},opt3_out{1},' times faster (decelerating)');
      elseif optX(3) < 0
          dw_txt = strcat('Bestfit down-wearing rate: past rate was',{' '},opt3_out{1},' times slower (accelerating)');
      else
          dw_txt = strcat('Bestfit down-wearing rate: same as present (',string(rateDW.present),{' '},'mm yr)');
      end
      disp(cliff_txt{:})
      disp(time_txt{:})
      disp(dw_txt{:})
      disp(cover_txt{:})
      fprintf('(reduced chi-squared of %0.2f for %.f DOF)\n',...
              fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
  end
    
  bestfit = optX;
  
end

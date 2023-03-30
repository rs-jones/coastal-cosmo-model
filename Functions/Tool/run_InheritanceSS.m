function bestfit = run_InheritanceSS(inputs,erorate_initial,sample_data,inh_meas,plot_fig)
  
  % Get nuclide parameters based on inheritance sample
  fprintf('\nInheritance \nCalculating nuclide parameters...');
  inh_meas.sample_idx = find(inh_meas.conc_mean == sample_data.measured_raw(:,2));
  inh_meas.conc = [inh_meas.conc_mean,inh_meas.conc_uncert];
  inh_meas.rho = sample_data.CC(inh_meas.sample_idx,6);
  inh_meas.z_gcm2 = inh_meas.depth * 100 * inh_meas.rho;
  
  maxage=8160; % 6* half-life = 8160 ka
  erosion = erorate_initial/10 *10; % cm/yr *10
  if erorate_initial == 0
      sample_data.maxdepth = max(inh_meas.z_gcm2)*2 +10000; % sample depth*2 + safety factor
  else
      erosion = erorate_initial/10 *10; % cm/yr *10
      sample_data.maxdepth = max(inh_meas.z_gcm2) + (maxage*1000)*erosion*inh_meas.rho;% +10000; % sample depth + maxage*erosion(cm/yr)*density + safety factor
  end
  sample_data = get_pars(sample_data,inputs.scaling_model);
  
  
  % Set optimisation rules
  opts = optimset('fminsearch');
  if erorate_initial == 0
      opts = optimset(opts,'TolFun',1e-2,'TolX',1e-2);
  else
      opts = optimset(opts,'TolFun',1e-5,'TolX',1e-5);
  end
  
  
  % Find bestfit scenario and export result
  fprintf('\nFinding best-fit scenario...');
  
  [optX,fmin] = fminsearch(@(X) fit_Inheritance_steadystate(X,sample_data,inh_meas,plot_fig),erorate_initial,opts);
  
  if optX(1) <0
      optX(1) = 0;
  end
  
  if optX(1) == 0
      fprintf('\nBestfit surface erosion: zero mm yr \n(chi-squared of %0.2f)\n',...
          fmin);
  elseif optX(1) < 0.001
      fprintf('\nBestfit surface erosion: <0.001 mm yr \n(chi-squared of %0.2f)\n',...
          fmin);
  else
      fprintf('\nBestfit surface erosion: %0.3f mm yr \n(chi-squared of %0.2f)\n',...
          round(optX(1)),fmin);
  end
  
  bestfit = optX;
  
end

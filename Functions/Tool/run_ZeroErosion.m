function bestfit = run_ZeroErosion(inputs,totaltime_initial,cover,sample_data,misfit_meas,plot_fig)
  
  % Get nuclide parameters based on samples
  fprintf('\nZero erosion (no cliff retreat or down-wearing) \nCalculating nuclide parameters...');
  sample_data = get_pars(sample_data,inputs.scaling_model);
  
  % Set optimisation rules
  opts = optimset('fminsearch');
  opts = optimset(opts,'TolFun',1e-1,'TolX',1e-1);
  
  
  % Find bestfit scenario and export result
  fprintf('\nFinding best-fit scenario...\n');
  
  if isempty(cover)
      [optX,fmin] = fminsearch(@(X) fit_ZeroErosion(X,inputs,sample_data,misfit_meas,plot_fig),totaltime_initial,opts);
      
      fprintf('\nBestfit total time: %.f years \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
          round(optX),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
      
  else
      inputs.cover_density = cover.density;
      [optX,fmin] = fminsearch(@(X) fit_ZeroErosion(X,inputs,sample_data,misfit_meas,plot_fig),[totaltime_initial,cover.depth_initial],opts);
      
      if optX(2) <0
          optX(2) = 0;
      end
      fprintf('\nBestfit total time: %.f years \nBestfit cover depth: %.f m \n(reduced chi-squared of %0.2f for %.f DOF)\n',...
          round(optX(1)),round(optX(2)),fmin./(length(sample_data.s)-3),length(sample_data.s)-3);
  end
    
  bestfit = optX;
  
end

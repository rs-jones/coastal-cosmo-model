
function misfit = fit_Inheritance_steadystate(X,sample_data,inh_meas,plot_fig)

% UNPACK
ero_gcm2yr = (X(1) ./10) .* inh_meas.rho;
if ero_gcm2yr < 0
    ero_gcm2yr = 0;
end


% GET PARAMETERS
pars.scaling_model = sample_data.scaling_model;
pars.pp = sample_data.pp;
pars.sf = sample_data.sf1026{inh_meas.sample_idx};
pars.sp = sample_data.sp1026{inh_meas.sample_idx};
pars.cp = sample_data.cp1026{inh_meas.sample_idx};
pars.l = sample_data.pp.lambda10Be;


% CALCULATE CONCENTRATIONS ALONG PLATFORM

% Calculate depth profile
depth_cm = 0:10:(inh_meas.depth*2*100); % Make profile twice the depth of the sample
z_gcm2 = depth_cm .* inh_meas.rho;
[Ptotal,Psp,Pmu] = PR_Z(z_gcm2,pars.pp,pars.sf,pars.cp,10);
A = inh_meas.rho/pars.cp.Lambdafe;

pred_concs.sp_z = Psp ./ (pars.l + ero_gcm2yr ./ A);
pred_concs.mu_z = Pmu ./ (pars.l + ero_gcm2yr ./ A);
pred_concs.total_z = pred_concs.sp_z + pred_concs.mu_z;

% Get concentration at sample depth
bottom_z_gcm2 = inh_meas.z_gcm2 + sample_data.CC(inh_meas.sample_idx,5)*inh_meas.rho;
top_z_conc = interp1(z_gcm2,pred_concs.total_z,inh_meas.z_gcm2);
bottom_z_conc = interp1(z_gcm2,pred_concs.total_z,bottom_z_gcm2);
pred_concs.sample = mean([top_z_conc,bottom_z_conc]);


% CALCULATE MISFIT FOR INHERITANCE SAMPLE
misfit = ((pred_concs.sample - inh_meas.conc(:,1))./inh_meas.conc(:,2)) .^2;


% PLOT CONCENTRATIONS FOR SCENARIO
if ~isempty(plot_fig)
    pred_concs.depth_m = depth_cm/100;
    
    figure(plot_fig)
    plot_concs_depth(inh_meas,pred_concs)
    
    drawnow;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%% Calculate concentrations %%%%%%%%%%%%%%%%%%%%%%%%
function concs = calc_conc(pars,expo_t,ero_gcm2yr,sTopo,sWater,sCover)

% Compute depths at the end of the exposure period
top_z = pars.top_z_gcm2;
bottom_z = pars.bottom_z_gcm2;
concs = zeros(1,length(expo_t(1,:)));

% Integrate over time and depth for each point along platform
for i = 1:length(expo_t(1,:))
    
    if ~isnan(expo_t(i)) && (expo_t(i)~=0)
        % Integrate each sample from top depth to bottom depth and from
        % zero to the exposure time in this scenario.
        % Divide integral by sample thickness.
        pr_func = @(z,t) (PR_Z((z + ero_gcm2yr(i).*t),pars.pp,pars.sf,pars.cp,10)...
            .* sTopo(i) .* sWater(i) .* sCover(i)) .* exp(-pars.l.*t);
        
        concs(i) = integral2(pr_func,top_z,bottom_z,0,expo_t(i),'RelTol',1e-3,'AbsTol',1e-3)...
            ./ (bottom_z - top_z);
    else
        concs(i) = NaN;
    end
end

end

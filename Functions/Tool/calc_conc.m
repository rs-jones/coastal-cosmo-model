%%%%%%%%%%%%%%%%%%%%%%%%% Calculate concentrations %%%%%%%%%%%%%%%%%%%%%%%%
function concs = calc_conc(pars,expo_t,ero_gcm2yr,sTopo,sWater,sCover)

% Compute depths at the end of the exposure period
top_z = pars.top_z_gcm2;
bottom_z = pars.bottom_z_gcm2;

if isfield(pars,'sf10') && isfield(pars,'sf14')
    concs.Be10 = zeros(1,length(expo_t(1,:)));
    concs.C14 = zeros(1,length(expo_t(1,:)));
else
    concs = zeros(1,length(expo_t(1,:)));
end

% Integrate over time and depth for each point along platform
for i = 1:length(expo_t(1,:))
    
    if ~isnan(expo_t(i)) && (expo_t(i)~=0)
            
            % Integrate each sample from top depth to bottom depth and from
            % zero to the exposure time in this scenario.
            % Divide integral by sample thickness.
            if isfield(pars,'sf10')
                pr_func = @(z,t) (PR_Z((z + ero_gcm2yr(i).*t),pars.pp,pars.sf10,pars.cp10,10)...
                    .* sTopo(i) .* sWater(i) .* sCover(i)) .* exp(-pars.l10.*t);
                
                concs_Be10(i) = integral2(pr_func,top_z,bottom_z,0,expo_t(i),'RelTol',1e-3,'AbsTol',1e-3)...
                    ./ (bottom_z - top_z);
            end
            if isfield(pars,'sf14')
                pr_func = @(z,t) (PR_Z((z + ero_gcm2yr(i).*t),pars.pp,pars.sf14,pars.cp14,14)...
                    .* sTopo(i) .* sWater(i) .* sCover(i)) .* exp(-pars.l14.*t);
                
                concs_C14(i) = integral2(pr_func,top_z,bottom_z,0,expo_t(i),'RelTol',1e-3,'AbsTol',1e-3)...
                    ./ (bottom_z - top_z);
            end
    else
        if isfield(pars,'sf10')
            concs_Be10(i) = NaN;
        end
        if isfield(pars,'sf14')
            concs_C14(i) = NaN;
        end
    end
end

% Export result
if isfield(pars,'sf10') && isfield(pars,'sf14')
    if isfield(pars,'sf10')
        concs.Be10 = concs_Be10;
    end
    if isfield(pars,'sf14')
        concs.C14 = concs_C14;
    end
else
    if isfield(pars,'sf10')
        concs = concs_Be10;
    end
    if isfield(pars,'sf14')
        concs = concs_C14;
    end
end

end

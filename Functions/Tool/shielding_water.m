%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Water shielding %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sw = shielding_water(elev_relRSL,tides,HAT)

rho = 1.024; % Water density (g/cm3)
lambda = 160; % Attenuation length (g/cm2)

% If platform elevation is higher than highest astronomical tide
if elev_relRSL > HAT
    sw = 1; return % Give it no water shielding
end


% Evaluate for each tidal elevation (cumulatively adding shielding)
this_sw = 0;
for t = 1:length(tides(:,1))
    
    % If elev (rel. to rsl) is greater than this tidal elev (+5cm),
    % add tidal duration to the shielding value
    if elev_relRSL >= tides(t,3)+0.05
        this_sw = this_sw + tides(t,2)/100;
        
    % If elev (rel. to rsl) is within this tidal elev (+/- 5cm),
    % calculate and add shielding from water depth
    elseif elev_relRSL < tides(t,3)+0.05  &&  elev_relRSL >= tides(t,3)-0.05
        this_tide_depth = (tides(t,3)+0.05-elev_relRSL) *100; % Tide elev. (+5cm) minus platform elev. (cm)
        this_tide_sw = tides(t,2) * ((tides(t,3)+0.05-elev_relRSL)/0.1)/100 + ...
            tides(t,2) * (1-(tides(t,3)+0.05-elev_relRSL)/0.1)/100 * exp(-rho*this_tide_depth/lambda); % Water attenuation for this tidal duration
        this_sw = this_sw + this_tide_sw;
        
    % Otherwise (if elev is less than this tidal elev),
    % calculate and add shielding from water depth
    else
        this_tide_depth = (tides(t,3)-0.05-elev_relRSL) *100; % Tide elev. (-5cm) minus platform elev. (cm)
        this_tide_sw = tides(t,2)/100 * exp(-rho*this_tide_depth/lambda); % Water attenuation for this tidal duration
        this_sw = this_sw + this_tide_sw;
    end
end

sw = this_sw;

end

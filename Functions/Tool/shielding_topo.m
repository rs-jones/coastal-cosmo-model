%%%%%%%%%%%%%%%%%%%%%%%%%%% Topographic shielding %%%%%%%%%%%%%%%%%%%%%%%%%
function [st,dist] = shielding_topo(topo,profile,cliff,expo)

dist_cliff = profile(:,1); % Distance from the cliff (m)

dist = zeros(length(cliff),length(profile(:,1)));
s_topo = zeros(length(cliff),length(profile(:,1)));


% Calculate for each point along the platform
for w = 1:length(profile(:,1))
    
    % Calculate through time
    for n = 1:length(cliff)
        
        % Do only if exposure time of point on platform is greater than this time
        if expo(w)>=n
            dist(n,w) = dist_cliff(w)-cliff(n); % Determine distance to cliff, for this point and time
            
            if dist(n,w)<=0 % If distance is less than zero, use shielding next to present cliff
                s_topo(n,w) = topo(1,2);
                
            else % Otherwise, calculate shielding from distance to cliff
                for w2=2:length(profile(:,1))
                    % Relative distance of this platform point to the topo
                    % points, multiplied by the change in shielding, plus
                    % the shielding of the nearest inland point
                    if dist(n,w)>topo(w2-1,1) && dist(n,w)<=topo(w2,1)
                        this_s_topo = (dist(n,w)-topo(w2-1,1)) * (topo(w2,2)-topo(w2-1,2))...
                            + topo(w2-1,2);
                        s_topo(n,w) = this_s_topo;
                    end
                end
            end
        else
            s_topo(n,w) = NaN;
        end
    end
end

st = sum(s_topo,1,'omitnan')./expo; % Time-averaged shielding along the platform

end

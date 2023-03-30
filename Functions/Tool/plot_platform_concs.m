
function fig_h = plot_platform_concs(inputs,sample_data,pred_concs)

dist_cliff = inputs.profile(:,1);
real = inputs.profile(:,2);

if nargin == 2
    
    fig_h = figure;
    
    blue = [49,130,189]/255;
    
    subplot(2,1,1)
    plot(dist_cliff,real,'-k','LineWidth',1.5);
    hold on; yline(inputs.HAT,':','Color',blue); % HAT
    hold on; yline(inputs.MHWS,'-.','Color',blue); % MHWS
    hold on; yline(inputs.MHWN,'--','Color',blue); % MHWN
    hold on; yline(inputs.MLWN,'--','Color',blue); % MLWN
    hold on; yline(inputs.MLWS,'-.','Color',blue); % MLWS
    hold off;
    text(dist_cliff(end-10),inputs.HAT,'HAT','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MHWS,'MHWS','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MHWN,'MHWN','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MLWN,'MLWN','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MLWS,'MLWS','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    ylabel('Elevation (m AOD)');
    
    subplot(2,1,2)
    errorbar(sample_data.measured_inh(:,1),sample_data.measured_inh(:,2)/1000,sample_data.measured_inh(:,3)/1000,'color','k','LineStyle','none');
    hold on;
    scatter(sample_data.measured_inh(:,1),sample_data.measured_inh(:,2)/1000,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
    hold off;
    xlabel('Distance from the cliff (m)')
    ylabel('^{10}Be concentration (atoms/g)')
    
else
    blue = [49,130,189]/255;
    cols = lines(5);
    if isfield(inputs,'elev_beach')
        max_y = max([inputs.HAT,max(inputs.elev_beach)])+0.5;
    else
        max_y = max([inputs.HAT,max(real)])+0.5;
    end
    if isfield(inputs,'profile_ero_expo')
        min_y = min([inputs.MLWS,min(inputs.profile_ero_expo)])-0.5;
    else
        min_y = min([inputs.MLWS,min(real)])-0.5;
    end
    
    ax1 = subplot(2,1,1);
    real_h = plot(dist_cliff,real,'-k','LineWidth',1.5);
    hold on; yline(inputs.HAT,':','Color',blue); % HAT
    hold on; yline(inputs.MHWS,'-.','Color',blue); % MHWS
    hold on; yline(inputs.MHWN,'--','Color',blue); % MHWN
    hold on; yline(inputs.MLWN,'--','Color',blue); % MLWN
    hold on; yline(inputs.MLWS,'-.','Color',blue); % MLWS
    text(dist_cliff(end-10),inputs.HAT,'HAT','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MHWS,'MHWS','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MHWN,'MHWN','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MLWN,'MLWN','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    text(dist_cliff(end-10),inputs.MLWS,'MLWS','Color',blue,'HorizontalAlignment','right','VerticalAlignment','baseline');
    ylabel('Elevation (m AOD)');
    if isfield(inputs,'cliff_pos')
        cliff_pos_idx = find(dist_cliff==inputs.cliff_pos);
        if isfield(inputs,'profile_ero_expo')
            cliff_pos_base = inputs.profile_ero_expo(cliff_pos_idx);
        else
            cliff_pos_base = real(cliff_pos_idx);
        end
        cliff_h = plot([inputs.cliff_pos,inputs.cliff_pos],[cliff_pos_base,max_y],'-','Color',[.7,.7,.7],'LineWidth',3);
        if inputs.cliff_pos < dist_cliff(end)
            text(inputs.cliff_pos-5,max_y-0.5,'Cliff','Color',[.7,.7,.7],'HorizontalAlignment','right','VerticalAlignment','baseline','Rotation',90);
        end
    end
    if isfield(inputs,'profile_ero_expo')
        ero_h = plot(dist_cliff,inputs.profile_ero_expo,'--','Color',[.7,.7,.7],'LineWidth',1.5);
        if isfield(inputs,'cliff_pos')
            ero_h = plot(dist_cliff(cliff_pos_idx:end),inputs.profile_ero_expo(cliff_pos_idx:end),'-','Color',[.7,.7,.7],'LineWidth',1.5);
        end
    end
    if isfield(inputs,'elev_beach')
        beach_h = plot(dist_cliff,inputs.elev_beach,'-','Color',cols(3,:),'LineWidth',1.5);
    end
    plot(dist_cliff,real,'-k','LineWidth',1.5);
    hold off;
    ylabel('Elevation (m AOD)');
    ylim([min_y,max_y]);
    if isfield(inputs,'profile_ero_expo') && isfield(inputs,'elev_beach')
        legend([real_h,ero_h,beach_h],'Present platform','Platform before down-wearing','Surface cover (mean, above present)','Location','SouthWest');
    elseif isfield(inputs,'profile_ero_expo') && ~isfield(inputs,'elev_beach')
        legend([real_h,ero_h],'Present platform','Platform before down-wearing','Location','SouthWest');
    elseif ~isfield(inputs,'profile_ero_expo') && isfield(inputs,'elev_beach')
        legend([real_h,beach_h],'Present platform','Surface cover (mean, above present)','Location','SouthWest');
    else
        legend(real_h,'Present platform','Location','Best');
    end
    
    ax2 = subplot(2,1,2);
    errorbar(sample_data.measured_inh(:,1),sample_data.measured_inh(:,2)/1000,sample_data.measured_inh(:,3)/1000,'color','k','LineStyle','none');
    hold on;
    scatter(sample_data.measured_inh(:,1),sample_data.measured_inh(:,2)/1000,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
    hold on;
    pred_h = plot(dist_cliff,pred_concs/1000,'-g','LineWidth',1.5);
    hold off;
    xlabel('Distance from the cliff (m)')
    ylabel('^{10}Be concentration (k atoms/g)')
    legend(pred_h,'Prediction','Location','SouthEast');
    
    ax1.XLim = ax2.XLim;
    
end

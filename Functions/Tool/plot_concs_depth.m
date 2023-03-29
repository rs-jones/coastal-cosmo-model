
function fig_h = plot_concs_depth(inh_meas,pred_concs)

if nargin == 1
    
    inh_meas.conc = [inh_meas.conc_mean,inh_meas.conc_uncert];
    
    fig_h = figure;
    
    subplot(2,2,[1,3]);
    errorbar(inh_meas.conc(:,1)/1000,inh_meas.depth,inh_meas.conc(:,2)/1000,'horizontal','color','k','LineStyle','none','LineWidth',1);
    hold on;
    h_smc = scatter(inh_meas.conc(:,1)/1000,inh_meas.depth,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1);
    h_smc.SizeData = 50;
    hold off;
    ax = gca; set(ax,'YDir','reverse'); 
    xlabel('Concentration (k atoms g^{-1})'); 
    ylabel('Depth below surface (m)');
    ax.XAxisLocation = 'top';
    ylim([0,inh_meas.depth*2])
    
else 
    subplot(2,2,[1,3]);

    h_tot = plot(pred_concs.total_z/1000,pred_concs.depth_m,'-g','LineWidth',1);
    hold on;
    h_sp = plot(pred_concs.sp_z/1000,pred_concs.depth_m,'--g','LineWidth',1);
    h_mu = plot(pred_concs.mu_z/1000,pred_concs.depth_m,'-.g','LineWidth',1);
    
    errorbar(inh_meas.conc(:,1)/1000,inh_meas.depth,inh_meas.conc(:,2)/1000,'horizontal','color','k','LineStyle','none','LineWidth',1);
    h_smc = scatter(inh_meas.conc(:,1)/1000,inh_meas.depth,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1);
    h_smc.SizeData = 100;
    
    h_spc = scatter(pred_concs.sample/1000,inh_meas.depth,'MarkerEdgeColor','g','MarkerFaceColor','none','LineWidth',1); hold on;
    h_spc.SizeData = 100;
    hold off;
    
    leg = legend([h_tot,h_sp,h_mu],'Total','spallogenic','muogenic','Location','southeast');
    title(leg,'Production');
    
    if inh_meas.depth*2 <= 10
       top_z_toPlot = 2; % 2 m
    elseif inh_meas.depth*2 <= 20
       top_z_toPlot = 3; % 3 m
    else
       top_z_toPlot = 5; % 5 m
    end
    x_depth_idx = find(pred_concs.depth_m == top_z_toPlot); % Find depth
    x_conc = pred_concs.total_z(x_depth_idx); % Get concentration
    xlim([0,x_conc/1000]);
    
    ax = gca; set(ax,'YDir','reverse'); 
    xlabel('Concentration (k atoms g^{-1})'); 
    ylabel('Depth below surface (m)');
    ax.XAxisLocation = 'top';
    
end

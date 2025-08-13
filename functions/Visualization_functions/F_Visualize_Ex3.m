function F_Visualize_Ex3(Y_ens,S_obs,S_true,S_ens,xPlus_3d,Sll_2d,Cxx_3d,Cxpxp_3d)

clr_true = 'k'; marker_true = '.'; line_true = 'k-';  msize_true = 15;
clr_obs  = 'k'; marker_obs  = '^'; line_obs  = 'k-';  msize_obs = 4;
clr_mod  = '#d10000'; marker_mod  = 'o'; line_mod  = '--';  clr_mod2 = '#d10000';
clr_upd  = '#527dea'; marker_upd  = 'o'; line_upd  = ':'; clr_upd2 = 'w';
clr_ens  = [0.8 0.8 0.8];          line_ens  = '-';
lw = 2;
msize = 6;
fsize = 12;
opt_save = 1;    
opt_saveData = 0;
Ne = size(Y_ens,2);
sim = [1; 24]; % simulation phase of model


%% Visualizations

% truth and observations
figure; hold on;
set(gca,'FontSize',fsize);
for nn=1:Ne
    h1 = plot(Y_ens(:,nn),line_ens,'Color',clr_ens,'LineWidth',lw);
end
h2 = plot(S_obs,line_obs,'Marker',marker_obs,'MarkerEdgeColor',clr_obs,'MarkerFaceColor',clr_obs,'LineWidth',lw,'MarkerSize',msize);
h3 = plot(S_true,line_true,'Marker',marker_true,'MarkerEdgeColor',clr_true,'MarkerFaceColor',clr_true,'LineWidth',lw,'MarkerSize',msize_true);
legend([h3 h2 h1], 'Truth', 'Observation', 'Ensemble (Obs)')
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('S')
title('Ground truth and observation ensemble');
box on
if(opt_save==1) print(gcf,'figures/Ex3_Obs','-djpeg'); end

% model simulation (open loop)
figure; hold on;
set(gca,'FontSize',fsize)
for nn=1:Ne
    h1 = plot(S_ens(:,nn),line_ens,'Color',clr_ens,'LineWidth',lw);
end
h2 = plot(mean(S_ens,2),line_mod,'Marker',marker_mod,'MarkerEdgeColor',clr_mod,'MarkerFaceColor',clr_mod2,'LineWidth',lw,'Color',clr_mod,'MarkerSize',msize);
legend([h2 h1], 'Model mean', 'Ensemble (Model)')
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('S')
title('Water storage estimates - OL');
box on
if(opt_save==1) print(gcf,'figures/Ex3_OL','-djpeg'); end

% DA results
id = logical(ones(1,sim(2)-sim(1)+1));
xPlus_ensMean = zeros(1,sim(2)-sim(1)+1); xPlus_ensMean(id) = mean(xPlus_3d,2);
figure; hold on;
set(gca,'FontSize',fsize)
for nn=1:Ne
    h1 = plot(squeeze(xPlus_3d(1,nn,:)),line_ens,'Color',clr_ens,'LineWidth',lw);
end
h2 = plot(xPlus_ensMean(1,:),line_upd,'Marker',marker_upd,'MarkerEdgeColor',clr_upd,'MarkerFaceColor',clr_upd,'LineWidth',lw,'Color',clr_upd,'MarkerSize',msize);
legend([h2 h1], 'Model mean', 'Ensemble (Model)')
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('S')
title('Water storage estimates - DA');
box on
if(opt_save==1) print(gcf,'figures/Ex3_DA','-djpeg'); end

% ensemble prediction
figure; hold on;
set(gca,'FontSize',fsize)
plot(mean(S_ens,2),line_mod,'Marker',marker_mod,'MarkerEdgeColor',clr_mod,'MarkerFaceColor',clr_mod2,'LineWidth',lw,'Color',clr_mod,'MarkerSize',msize);
plot(xPlus_ensMean(1,:),line_upd,'Marker',marker_upd,'MarkerEdgeColor',clr_upd,'MarkerFaceColor',clr_upd2,'LineWidth',lw,'Color',clr_upd,'MarkerSize',msize);
plot(S_obs,line_obs,'Marker',marker_obs,'MarkerEdgeColor',clr_obs,'MarkerFaceColor',clr_obs,'LineWidth',lw,'MarkerSize',msize);
plot(S_true,line_true,'Marker',marker_true,'MarkerEdgeColor',clr_true,'MarkerFaceColor',clr_true,'LineWidth',lw,'MarkerSize',msize_true);
legend('OL','DA','Observation','Truth')
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('S')
title('Water storage estimates (ensemble average)');
box on
if(opt_save==1) print(gcf,'figures/Ex3_EnsAv','-djpeg'); end

% empirical standard deviations
% observations - model prediction - model update
figure; hold on;
set(gca,'FontSize',fsize)
plot(sqrt(Sll_2d),line_obs,'Marker',marker_obs,'MarkerEdgeColor',clr_obs,'MarkerFaceColor',clr_obs,'LineWidth',lw,'MarkerSize',msize);
IN = logical(ones(24,1)); test1 = zeros(24,1);
test1(IN)=Cxx_3d(1,1,:);
plot(sqrt(test1),line_mod,'Marker',marker_mod,'MarkerEdgeColor',clr_mod,'MarkerFaceColor',clr_mod2,'LineWidth',lw,'Color',clr_mod,'MarkerSize',msize);
test2 = zeros(24,1); test2(IN)=Cxpxp_3d(1,1,:);
plot(sqrt(test2),line_upd,'Marker',marker_upd,'MarkerEdgeColor',clr_upd,'MarkerFaceColor',clr_upd2,'LineWidth',lw,'Color',clr_upd,'MarkerSize',msize);
xlabel('time'); xlim([0-0.5 (sim(2)-sim(1)+1)+0.5]);
ylabel('S sigma')
title('Uncertainty (based on ensemble spread)');
box on
legend('Observations','Model Prediction','Model Update')
if(opt_save==1) print(gcf,'figures/Ex3_Uncert','-djpeg'); end

end


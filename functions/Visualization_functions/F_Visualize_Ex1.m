function F_Visualize_Ex1(P,S,S0,Q,K)

lw = 2;
msize = 15;
fsize = 20;
opt_save = 1;

% precipitation
figure
set(gca,'FontSize',fsize)
bar(1:length(P),P,'k');
xlabel('time'); xlim([0-0.5 length(P)+0.5]);
set(gca, 'XTick', (0:10:50));
ylabel('P-E')
title('Input data: precipitation - evapotranspiration');
box on
if(opt_save==1) print(gcf,'figures/Ex1_P','-djpeg'); end

% storage
figure
set(gca,'FontSize',fsize)
plot(0:length(P),[S0; S],'.-k','LineWidth',lw,'MarkerSize',msize);
xlabel('time'); xlim([0-0.5 length(P)+0.5]);
set(gca, 'XTick', (0:10:50));
ylabel('S'); ylim([0 max([S0; S])+2]);
title('Water storage estimates - single model simulation');
box on
if(opt_save==1) print(gcf,'figures/Ex1_S','-djpeg'); end

% runoff/discharge
figure
set(gca,'FontSize',fsize)
plot(1:length(P),Q,'.-k','LineWidth',lw,'MarkerSize',msize);
xlabel('time'); xlim([0-0.5 length(P)+0.5]);
set(gca, 'XTick', (0:10:50));
ylabel('R'); ylim([0 max([S0; S])+2]);
title('River discharge estimates - single model simulation');
box on
if(opt_save==1) print(gcf,'figures/Ex1_R','-djpeg'); end

% relation between S and Q
S_sim = (0:0.1:25)';
Q_sim = K*S_sim;
figure
set(gca,'FontSize',fsize)
plot(S_sim,Q_sim,'-k','LineWidth',lw,'MarkerSize',msize);
xlabel('S'); ylabel('R');
title('Relation between storage and discharge - single model simulation');
box on
axis equal
xlim([0 25]); ylim([0 8]);
if(opt_save==1) print(gcf,'figures/Ex1_StoR','-djpeg'); end

end


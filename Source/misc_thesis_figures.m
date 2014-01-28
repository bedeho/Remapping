
%% Reference frame

figure;
hold on

E=60;
R=90;

plot([-R/2 R/2], 35+[-E/2 E/2],'k','LineWidth',3);
plot([-R/2 R/2], -10+[-E/2 E/2],'k','LineWidth',3);

xlim([-R/2 R/2]);
ylim([-E/2 E/2]);
pbaspect([E R 1]);
box on

hYLabel = ylabel('Saccade (deg)');
hXLabel = xlabel('Retinal Location (deg)');

set([hYLabel hXLabel], 'FontSize', 20);
set([gca], 'FontSize', 18);
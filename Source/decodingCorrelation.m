
R_preferences = -45:1:45;
S_preferences = -30:1:30;

[X Y] = meshgrid(R_preferences, S_preferences);
    
C_R = X(:);
C_S = Y(:);

% Load finished product
x = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/analysis-basic-CLayerProbe.mat');

% plot
figure;
prewiredLocation = C_R;
decodedLocation = x.CLabeProbe_Neurons_V;

plot(decodedLocation, prewiredLocation,'o');
hold on
plot([min(decodedLocation) max(decodedLocation)],[min(prewiredLocation) max(prewiredLocation)] , 'r', 'LineWidth',2);
axis square
axis tight
hXLabel = xlabel('Decoded Visual Preference (deg)');
hYLabel = ylabel('Prewired Visual Preference (deg)');

set([hYLabel hXLabel], 'FontSize', 14);
set(gca, 'FontSize', 12);

corrcoef(prewiredLocation, decodedLocation)


% plot
figure;
prewiredLocation = C_S;
decodedLocation = x.CLabeProbe_Neurons_S;

plot(decodedLocation, prewiredLocation,'o');
hold on
plot([min(decodedLocation) max(decodedLocation)],[min(prewiredLocation) max(prewiredLocation)] , 'r', 'LineWidth',2);
axis square
axis tight
hXLabel = xlabel('Decoded Saccade Preference (deg)');
hYLabel = ylabel('Prewired Saccade Preference (deg)');

set([hYLabel hXLabel], 'FontSize', 14);
set(gca, 'FontSize', 12);

corrcoef(prewiredLocation, decodedLocation)
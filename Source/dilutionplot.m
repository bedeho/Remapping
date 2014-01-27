
values = 0.1:0.1:1;

remapping_indexes_untrained = zeros(size(values));
remapping_indexes_trained   = zeros(size(values));

stim_I = zeros(size(values));
sacc_I = zeros(size(values));
trace_I = zeros(size(values));

for i=1:length(values),
    
    % Get v value
    v = values(i)
        
    % open trained STIM
    stim_file = load(['/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/dilution-' num2str(v) '/baseline/TrainedNetwork/activity-basic-StimuliControl.mat']);
    %stim_I(i) = max(stim_file.C_Input_flat(:, 46 + 26));
    
    % open trained SACC
    sacc_file = load(['/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/dilution-' num2str(v) '/baseline/TrainedNetwork/activity-basic-SaccadeControl.mat']);
    %sacc_I(i) = max(sacc_file.C_Input_flat(:, 31 + 26));
    
    % open untrained duhamel trace
    trace_file = load(['/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/dilution-' num2str(v) '/baseline/BlankNetwork/activity-basic-DuhamelRemappingTrace.mat']);
    %trace_I(i) = max(trace_file.C_Input_flat(:, 1));
    
    % trace analysis
    trace_analysis_file_untrained = load(['/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/dilution-' num2str(v) '/baseline/BlankNetwork/analysis-basic-DuhamelRemappingTrace.mat']);
    trace_analysis_file_trained   = load(['/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/dilution-' num2str(v) '/baseline/TrainedNetwork/analysis-basic-DuhamelRemappingTrace.mat']);
    
    remapping_indexes_untrained(i) = trace_analysis_file_untrained.DuhamelRemappingTrace_Result.remapping_index;
    remapping_indexes_trained(i) = trace_analysis_file_trained.DuhamelRemappingTrace_Result.remapping_index;    
end

% Remapping index
%{
figure;
plot(values,remapping_indexes_untrained, 'r');
hold;
plot(values,remapping_indexes_trained, 'b');

hYLabel = ylabel('Remapping Index');
hXLabel = xlabel('Connectivity rate');
hLegend = legend('Untrained','Trained');

set(gca,'XTick', values);
set([hYLabel hXLabel], 'FontSize', 20);
set([gca hLegend], 'FontSize', 18);
legend('boxoff');

ylim([0 1]);
xlim([0.1 1]);
%}

% Plot

figure;
hold on;
plot(values,stim_I,'r');
plot(values,sacc_I,'b');
plot(values,trace_I,'g');

hYLabel = ylabel('Presynaptic Excitation');
hXLabel = xlabel('Connectivity rate');
hLegend = legend('Stimulus Control (trained)','Saccade Control (trained)', 'Single Step (untrained)');

set(gca,'XTick', values);
set([hYLabel hXLabel], 'FontSize', 20);
set([gca hLegend], 'FontSize', 18);
legend('boxoff');

ylim([0 1]);
xlim([0.1 1]);

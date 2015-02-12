
%
%  remappingPlots.m
%  Remapping
%
%  Created by Bedeho Mender on 27/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [remLatFig, remScatFig, indexFig] = remappingPlots(remapping_results, FaceColors, Legends)

    AxisFontSize = 12;
    LabelFontSize  = 14;

    % Iterate data sets
    lat_lower_limit = inf;
    lat_upper_limit = -inf;
    foundNaN = false;
    
    for i=1:length(remapping_results),

        remapping_result = remapping_results{i};

        % Index data
        X_Indx{i} = [remapping_result(:).sacc_index];
        Y_Indx{i} = [remapping_result(:).stim_index];

        % Latency data: turn into ms for plotting
        X_Lat{i} = 1000*[remapping_result(:).stimLatency];
        Y_Lat{i} = 1000*[remapping_result(:).remappingLatency];

        % Fix latency plot limits
        minStimLat = min(X_Lat{i});
        maxStimLat = max(X_Lat{i});

        minRemLat = min(Y_Lat{i});
        maxRemLat = max(Y_Lat{i});

        % Limits for this data set
        minLat = min(minStimLat, minRemLat);
        maxLat = max(maxStimLat, maxRemLat);

        % Limits for all data sets so far
        lat_lower_limit = min(lat_lower_limit, minLat);
        lat_upper_limit = max(lat_upper_limit, maxLat);
        
        % Did we have NaN in latency
        foundNaN = foundNaN ||any(isnan(X_Lat{i})) || any(isnan(Y_Lat{i}));

    end
    
    % 1. scatter remap latency vs. stim control latency
    Lim = [(lat_lower_limit - 5) (lat_upper_limit+5)]; %1.1*,  ticks on 10ms
    Lim = roundn(Lim,1);
    
    if(true), % ~foundNaN
    
        if(length(remapping_results) > 1),
            [remLatFig, yProjectionAxis, scatterAxis, xProjectionAxis, XLim, YLim] = scatterPlotWithMarginalHistograms(X_Lat, Y_Lat, 'XTitle', 'Stimulus Control Latency (ms)', 'YTitle', 'Remapping Latency (ms)', 'FaceColors', FaceColors, 'XLim', Lim, 'YLim', Lim, 'Legends', Legends, 'AxisFontSize', AxisFontSize, 'LabelFontSize', LabelFontSize);
        else
            [remLatFig, yProjectionAxis, scatterAxis, xProjectionAxis, XLim, YLim] = scatterPlotWithMarginalHistograms(X_Lat, Y_Lat, 'XTitle', 'Stimulus Control Latency (ms)', 'YTitle', 'Remapping Latency (ms)', 'FaceColors', FaceColors, 'XLim', Lim, 'YLim', Lim, 'AxisFontSize', AxisFontSize, 'LabelFontSize', LabelFontSize);
        end
        
        axes(scatterAxis);     
        hold on;
        plot(Lim,Lim,'--k'); % Add x=y diagonal
        ticks = Lim(1):20:Lim(2);
        set(scatterAxis, 'YTick', ticks, 'XTick', ticks);
    
    else
        remLatFig = figure;
        warning('missing latency');
        imagesc;
    end
    
    % 2. scatter stim index. vs sacc index.    
    Lim = [-1 1];
    
    if(length(remapping_results) > 1),
        [remScatFig, yProjectionAxis, scatterAxis, xProjectionAxis, XLim, YLim] = scatterPlotWithMarginalHistograms(X_Indx, Y_Indx, 'XTitle', 'Saccade Index', 'YTitle', 'Stimulus Index', 'FaceColors', FaceColors, 'XLim', Lim, 'YLim', Lim, 'Legends', Legends, 'AxisFontSize', AxisFontSize, 'LabelFontSize', LabelFontSize);
    else
        [remScatFig, yProjectionAxis, scatterAxis, xProjectionAxis, XLim, YLim] = scatterPlotWithMarginalHistograms(X_Indx, Y_Indx, 'XTitle', 'Saccade Index', 'YTitle', 'Stimulus Index', 'FaceColors', FaceColors, 'XLim', Lim, 'YLim', Lim, 'AxisFontSize', AxisFontSize, 'LabelFontSize', LabelFontSize);
    end
    
    axes(scatterAxis);     
    hold on;
    plot(Lim,Lim,'--k');     % Add x=y diagonal
    plot([0 0],[-1 1],'--g'); % x=0 bar
    plot([-1 1],[0 0],'--g'); % y=0 bar
    
    % 3. remapping index distibution
    %indexFig = figure('Units','pixels','position', [1000 1000 620 300]);
    indexFig = figure('Units','pixels','position', [1000 1000 300 300]);
    hold on;
    
    for i=1:length(remapping_results),
        
        remapping_result = remapping_results{i};
        index = fliplr(sort([remapping_result(:).remapping_index]));
        plot(index, 'Color' , FaceColors{i});

    end
    
    if(length(remapping_results) > 1),
        legend(Legends);
        legend('boxoff');
    end
    
    if(length(index) > 1);
        xlim([1 length(index)]);
    end
    
    ylim([-0.1 sqrt(2)]);

    hXLabel = xlabel('Neuron Rank');
    hYLabel = ylabel('Remapping Index');

    set(gca, 'TickDir', 'out', 'FontSize', 10); % 14 før
    set([hYLabel hXLabel], 'FontSize', 10);
    %pbaspect([1 0.5 1]);
    box on;
        
end
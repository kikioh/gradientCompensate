clearvars; close all; clc;

% specify sample, sequence and nuclei
sample = 'glucose';
seqName = 'HSQC';
nuclei = '1H';

% Load the saved .fig files
fileName = append(sample,'-quadet-HSQC-',nuclei);
fig1 = openfig(append(fileName,'-NM.fig'), 'invisible');
fig2 = openfig(append(fileName,'-GC.fig'), 'invisible');
fig3 = openfig(append(fileName,'-SL.fig'), 'invisible');

% Get handles to the axes and lines in each figure
axes1 = findobj(fig1, 'type', 'axes');
axes2 = findobj(fig2, 'type', 'axes');
axes3 = findobj(fig3, 'type', 'axes');

% Create a new figure
mergedFig = figure();
mergedFig.Position = [100,100,350,170];

% Offset for x and y axes for each sub-plot
yOffsets = [0, 0, 0];
switch sample
    case 'glucose'
        yOffsets = [0, 2, 4]*1e4;
        switch nuclei
            case '1H'; xOffsets = [0, 0.07, 0.14]; 
            case '13C'; xOffsets = [0, 1, 2]; 
        end
    case 'glycine'
        switch nuclei
            case '1H'; xOffsets = [0, 0.2, 0.4]; 
            case '13C'; xOffsets = [0, 2, 4];
        end
end

% Copy the lines from each figure to the new figure with original color and apply offsets
copyobj(findobj(axes3, 'type', 'line'), gca);
copyobj(findobj(axes2, 'type', 'line'), gca);
copyobj(findobj(axes1, 'type', 'line'), gca);

% Manually set the color and line width for specific lines
lineHandles = findall(gca, 'Type', 'Line');
for i = 1:numel(lineHandles)

    % set shift
    xdata = get(lineHandles(i), 'XData') + xOffsets(i);
    ydata = get(lineHandles(i), 'YData') + yOffsets(i);
    set(lineHandles(i), 'XData', xdata, 'YData', ydata);

    % set color
    colors = [0.3 0.3 0.3;0.8500 0.3250 0.0980;0 0.4470 0.7410];
    linewidths = [0.1, 0.1, 0.1];
    set(lineHandles(i), 'Color', colors(i,:), 'LineWidth', linewidths(i));
end

% set x axis
% axis off
set(gca,'XDir','reverse');
switch sample
    case 'glucose'
        switch nuclei
            case '1H'; xlim([3,6]);
            case '13C'; xlim([60,110]);
        end
    case 'glycine'
        switch nuclei
            case '1H'; xlim([2,6]);
            case '13C'; xlim([-60,-20]);
        end
end

set(gca, 'YColor', 'none');
xlabel(append(nuclei,' Chemical shift, ppm'));
set(gca, 'Fontsize', 8.5, 'FontName','Airal');
% save figure
% mergedFigName = append(sample,'-',seqName,'-',nuclei,'-merged-offset.fig');
% saveas(mergedFig,mergedFigName);


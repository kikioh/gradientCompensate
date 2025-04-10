% pick up the peak values for strong VS weak coupling
%
% Note: the first is value under standard hsqc, the second row is hqsc with
% spin lock
%
% mengjia.he@kit.edu, 2025.02.10

clearvars; close all; clc;

fig = figure();
fig.Position = [100,100,800,600];

% plot values
peaks = load('peaks_glucose_SW.mat');

% 1H peaks of alpha-glucose
subplot(2,2,1)
plot(peaks.H_alpha_weak,'LineWidth',1); hold on;
plot(peaks.H_alpha_strong,'LineWidth',1);
xticks(floor(min(xlim)):ceil(max(xlim))); 
xlabel('peak number');
ylabel('relative intensity');
legend('weak','strong');
set(gca,'FontSize',11,'FontName','Arial','linewidth',0.5);
title('1H projection, alpha-glucose');
ylim([0.9 1.0]);

% 13C peaks of alpha-glucose
subplot(2,2,2)
plot(peaks.C_alpha_weak,'LineWidth',1); hold on;
plot(peaks.C_alpha_strong,'LineWidth',1);
xticks(floor(min(xlim)):ceil(max(xlim))); 
xlabel('peak number');
ylabel('relative intensity');
legend('weak','strong');
set(gca,'FontSize',11,'FontName','Arial','linewidth',0.5);
title('13C projection, alpha-glucose');
ylim([0.9 1.0]);


% 1H peaks of beta-glucose
subplot(2,2,3)
plot(peaks.H_beta_weak,'LineWidth',1); hold on;
plot(peaks.H_beta_strong,'LineWidth',1);
xticks(floor(min(xlim)):ceil(max(xlim))); 
xlabel('peak number');
ylabel('relative intensity');
legend('weak','strong');
set(gca,'FontSize',11,'FontName','Arial','linewidth',0.5);
title('1H projection, beta-glucose');
ylim([0.9 1.0]);

% 13C peaks of beta-glucose
subplot(2,2,4)
plot(peaks.C_beta_weak,'LineWidth',1); hold on;
plot(peaks.C_beta_strong,'LineWidth',1);
xticks(floor(min(xlim)):ceil(max(xlim))); 
xlabel('peak number');
ylabel('relative intensity');
legend('weak','strong');
set(gca,'FontSize',11,'FontName','Arial','linewidth',0.5);
title('13C projection, beta-glucose');
ylim([0.9 1.0]);

%% read values
% % specify sample, sequence and nuclei
% sample = 'glucopyranose_beta';
% seqName = 'HSQC';
% nuclei = '13C';
% 
% % Load the saved .fig files
% fileName = append(sample,'-weakCouple-HSQC-',nuclei);
% fig1 = openfig(append(fileName,'-NM.fig'), 'invisible');
% % fig2 = openfig(append(fileName,'-GC.fig'), 'invisible');
% fig3 = openfig(append(fileName,'-SL.fig'), 'invisible');
% 
% % Get handles to the axes and lines in each figure
% axes1 = findobj(fig1, 'type', 'axes');
% % axes2 = findobj(fig2, 'type', 'axes');
% axes3 = findobj(fig3, 'type', 'axes');
% 
% % Create a new figure
% % mergedFig = figure();
% % mergedFig.Position = [100,100,800,300];
% 
% % Offset for x and y axes for each sub-plot
% xOffsets = [0, 0, 0]; 
% yOffsets = [0, 0, 0]; 
% switch sample
%     case {'glucose','glucopyranose_alpha','glucopyranose_beta','glucose-beta'}
% %         yOffsets = [0, 2, 4]*1e4;
%         switch nuclei
%             case '1H'; xOffsets = [0, 0.04, 0.14]; 
%             case '13C'; xOffsets = [0, 1, 2]; 
%         end
%     case 'glycine'
%         switch nuclei
%             case '1H'; xOffsets = [0, 0.2, 0.4]; 
%             case '13C'; xOffsets = [0, 2, 4];
%         end
% end
% 
% % Copy the lines from each figure to the new figure with original color and apply offsets
% copyobj(findobj(axes3, 'type', 'line'), gca);
% % copyobj(findobj(axes2, 'type', 'line'), gca);
% copyobj(findobj(axes1, 'type', 'line'), gca);
% 
% % Manually set the color and line width for specific lines
% lineHandles = findall(gca, 'Type', 'Line');
% peaks1 = get(lineHandles(1), 'YData');
% peaks1 = findpeaks(peaks1); 
% peaks1 = peaks1(peaks1 >= 10000);
% 
% peaks2 = get(lineHandles(2), 'YData');
% peaks2 = findpeaks(peaks2); 
% peaks2 = peaks2(peaks2 >= 10000);
% 
% % SL VS standard
% C_beta_weak = peaks2 ./peaks1;
% save('peaks_glucose_SW.mat', 'C_beta_weak', '-append');

% weak coupling value for alpha-glucose
% H_weak_alpha = [];
% C_weak_alpha = [];
% 
% % strong coupling value for alpha-glucose
% H_strong_alpha = [];
% C_strong_alpha = [];
% 
% % weak coupling value for beta-glucose
% H_weak_beta = [];
% C_weak_beta = [];
% 
% 
% % strong coupling value for beta-glucose
% H_srong_beta = [];
% C_srong_beta = [];

% %%
% H_beta_strong = load('peaks_glucose_SW.mat').H_beta_strong;
% H_beta_strong = [H_beta_strong(1),ans,H_beta_strong(2:end)];
% save('peaks_glucose_SW.mat', 'H_beta_strong', '-append');


%%
% C_beta_weak = load('peaks_glucose_SW.mat').C_beta_weak;
% C_beta_weak = [C_beta_weak(1:3),ans,C_beta_weak(4:end)];
% save('peaks_glucose_SW.mat', 'C_beta_weak', '-append');

%%
% C_beta_strong = load('peaks_glucose_SW.mat').C_beta_strong;
% C_beta_strong = [C_beta_strong(1:3),ans,C_beta_strong(4:end)];
% save('peaks_glucose_SW.mat', 'C_beta_strong', '-append');


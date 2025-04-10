% read real spectrum from topspin exported txt
%
% Mengjia He, 2024.10.25

% topspin examdata/20241206, glycine at detector 1, glucose at detector 2
% normal                        spec-80               spec-14
% GC                            spec-82               spec-16
% SL                            spec-84               spec-63

close all; clearvars; clc;

% specify sample, sequence and nuclei
sample = 'glycine';
seqName = 'HSQC';

% load FID
fileName = append('spec-',sample,'-',seqName,'.txt');
[spectrum,nRow,nCol]  = readTopSpinTxt('spec-84.txt');
spectrum = spectrum/1e8;

% define frequency axis
SW_F2 = 10;
SW_F1 = 150;

switch sample
    case 'glucose' 
         
        spectrum = flip(spectrum, 2);
        F2_ref = 1;
        F1_ref = -4;
        parameters.ylimts = [50,110];
        parameters.xlimts = [3,6];
        parameters.clim = [0.1,4];
        parameters.cbTicks = [1 2,3,4]; 
    case 'glycine'
        spectrum = flip(spectrum, 1);
        F1_ref = 65;
        F2_ref = 1.4;
        parameters.ylimts = [40,50]; 
        parameters.xlimts = [3,4];
        parameters.clim = [0.5,12];
        parameters.cbTicks = [2, 4, 6, 8, 10, 12];  
end
axis_F2 = linspace(0,SW_F2,nCol)-F2_ref;
axis_F1 = linspace(0,SW_F1,nRow)-F1_ref;


%% visulize data
% parameters.signs = 'positive';
% parameters.m = 6;
% parameters.ncont = 10;
% parameters.smin = 0.5;
parameters.delta = [0.08 1.0 0.08 1.0];
close all;
fig1 = figure('Name',fileName);
fig1.Position = [100 100 340 340*0.75];
% scale_figure([1.5 2.0]);
plot_2dSpec(axis_F2,axis_F1,spectrum,parameters);
xlabel('1H Chemical shift, ppm');
ylabel('13C Chemical shift, ppm');
set(gca,'FontSize',11,'FontName','Arial','linewidth',0.7);
% clim([-1,12]); 
% cb.Ticks = [-1,0, 2, 4, 6, 8, 10, 12];  
% axis tight;
% xlim(xlimts);
% ylim([-0.1,5]);
% yticks([]);
% set(gca, 'XDir','reverse');
% xlabel(append(nuclei,' Chemical shift, ppm'));
% set(gca,'FontSize',11,'FontName','Arial','linewidth',1);
% box off;
% set(get(gca,'YAxis'),'visible','off');
% lgd = legend('Normal', 'GC', 'SL');
% set(lgd, 'Box', 'off');
% set(lgd, 'Color', 'none');

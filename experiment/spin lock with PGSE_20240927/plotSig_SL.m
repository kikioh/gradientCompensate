% polt the PGSE experiment results, with spin locking for gradient coupling,
% read the 2D spectra, and plot the signal intensity VS gradient coupling
% strength
%
% Mengjia He, 2024.09.012

clearvars; clc; close all;

% list of gradient amplitude
grad = linspace(0,0.95,20);
numGrad = numel(grad);

% Load and process the first spectrum 
spec_gc = readTopSpinTxt('spec-20240927-18-gc.txt',numGrad)/1e8;
sigPeak_gc = max(spec_gc, [], 2);
% sigPeak_gc = sum(spec_gc(:,3200:3800),2);
% sigPeak_gc = zeros(numGrad,1);
% for m = 1:numGrad
%     sigPeak_gc(m) = max(calArea(spec_gc(m,:),0.1,0.02));
%     
% end
sigPeak_gc = sigPeak_gc  / sigPeak_gc(1);  % Normalize by the value of zero gradient

% Load and process the second spectrum 
spec_sl = readTopSpinTxt('spec-20240927-19-sl.txt',numGrad)/1e8;
sigPeak_sl = max(spec_sl, [], 2);
% sigPeak_sl = sum(spec_sl(:,3200:3800), 2);
% sigPeak_sl = zeros(numGrad,1);
% for m = 1:numGrad
%     sigPeak_sl(m) = max(calArea(spec_sl(m,:),0.1,0.02));
% end
sigPeak_sl = sigPeak_sl  / sigPeak_sl(1);  % Normalize by the value of zero gradient

% plot fitted curve
fig1 = figure('Name','signal intensity of PGSE');
fig1.Position = ([100,100,500,400]);
plot(grad, sigPeak_gc(1:numGrad), 'o','LineWidth',1.5,'color',lineColor('lightorange')); hold on;
plot(grad, sigPeak_sl(1:numGrad), '*','LineWidth',1.5,'color',lineColor('lightblue'));
xlabel('coupled gradient amplitude, a.u.');
ylabel('signal intensity, a.u.');
legend('GC', 'CL');
set(gca,'FontSize',13,'FontName','Arial');
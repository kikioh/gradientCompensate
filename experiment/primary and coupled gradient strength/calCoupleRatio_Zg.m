% calculate the gradient coupling ratio with single pulsed gradient
% (coupling) and acq experimrnt
%
% Ref: the signal intetsity comparing to non-gradient case is R =
% sin(kx)/kx, where k = ls/2 * Gmax * gamma * grad_inter * couplingRatio,
% x is the relative gradient amplitude applied in experiment
% ls - smaple size;
% Gmax      - maximum gradient 
% gamma     - gyomagnetic ratio of observe nuclei
% grad_int  - time domain intergal of gradient pulse, i.e., 2*tau/pi for sin
% shape, 0.9*tau for rectangular shape with 10% ring up and 10% down duration 
% 
% Mengjia He, 2024.09.06


clearvars; clc; close all;
%% experiment data fitting

% experiment parameters
ls = 6.5e-3;                    % sample length, m
Gmax = 103e-4*1e2;              % gradient strength at maximum current 3 A, T/m
tau = 1e-3;                     % gradient duration, s
grad_inter = 1.4*2*tau/pi;      % time integral of rectangular gradient shape
gamma = spin('1H');             % observe nuclei
grad = linspace(0,0.95,20);     % list of gradient amplitude

% Remove the zero gradient
numGrad = numel(grad);
grad = grad(2:numGrad);

% Load and process spectrum data
[spec,numRow,numCol] = readTopSpinTxt('spec-20240909-38-Cp-rec.txt');

% Trim the spectrum to match the gradient size
spec_cut = spec(1:numGrad,:)/1e8;

% take peak value
sigPeak = max(spec_cut,[],2);
sigPeak = sigPeak/sigPeak(1);
sigPeak = sigPeak(2:end);

% curve fitting
fitFunc = fittype('sin(k*x)/(k*x)', 'independent', 'x', 'coefficients', {'k'});

% set initial value to k=1
fittedModel = fit(grad', sigPeak, fitFunc, 'StartPoint', 1);
k = coeffvalues(fittedModel);

% calculate coupling coefficient
CR_exp = k/(ls/2*gamma*Gmax*grad_inter);

% plot fitted curve
figure;
plot(grad, sigPeak, 'o','LineWidth',1.5); hold on;
plot(grad, fittedModel(grad), '-','LineWidth',1.5); % fitting curve
title(append('The spillover ratio is ',num2str(round(CR_exp, 4))));
xlabel('gradient amplitude, a.u.');
ylabel('signal intensity, a.u.');
legend('experiment','fittig curve');
set(gca,'FontSize',13,'FontName','Arial');

%% formula model

% % experiment parameters
% ls = 6.5e-3;                    % sample length, m
% Gmax = 103e-4*1e2;              % gradient strength at maximum current 3 A, T/m
% tau = 1e-3;                     % gradient duration, s
% grad_inter = 1.4*2*tau/pi;      % time integral of rectangular gradient shape
% gamma = spin('1H');             % observe nuclei
% coupleRatio =  1.92e-03;        % gradient coupling ratio, prefixed value

% % set gradient amplitude
% gradCoeff = linspace(0.05,0.95,19);
% 
% % sequence gradient ratio
% grad =  gradCoeff *  coupleRatio;
% 
% % sum of phase
% phase = Gmax * gamma * grad_inter * grad;
% 
% % wavelength of helix of the phase of the coherence
% lambda = 2*pi./phase;
% 
% % signal intensity ratio
% ratio = lambda./(pi*ls) .* sin(pi*ls ./ lambda);

% % plot ratio
% fig = figure('Name','Intensity of PGSE signal');
% fig.Position = [200 100 450 350];
% plot(gradCoeff,abs(ratio),'LineWidth',1.5);
% set(gca,'FontName','Arial','FontSize',14);
% xlabel('gradient amplitude');
% ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex');


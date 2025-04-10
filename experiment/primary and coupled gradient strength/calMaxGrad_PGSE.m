% use PGSE to test gradient strength, use voxalytic parallel probe, the sample is
% 10% H2O/90% D2O, only use detector 1, the Bruker GREAT gradient system was used, 
% the gradient Y was connected to detector 1 and the gradient Z was
% connected to detector 2
%
% ref: 200 and more book, Exp. 11.5
%
% Mengjia He, 2024.09.10

clearvars; close all; clc;

% set parameters
gamma = 2.675*1e8;                      % 1H gyomagnetic ratio
tau = 1e-3;                             % gradient pulse duration
Delta = 1e-3+5e-3+14.8e-6+1e-3;         % distence between two gradient
D = 1.9e-9;                             % diffusion constant of 10% H2O/90% D2O
grad_square = (linspace(0,0.95,20)).^2; % list of gradient amplitude

% define number of sweep
numGrad = numel(grad_square);

% load and process spectrum
spec = readTopSpinTxt('spec-20240909-34-Pr-sin.txt');

% Trim the spectrum to match the gradient size
spec_cut = spec(1:numGrad,:)/1e8;

% take peak value
sigPeak = max(spec_cut,[],2);
% sigPeak = zeros(1,numGrad);
% for m = 1:numGrad
%     sigPeak(m) = calArea(spec_cut(m,:),0.5,0.001);
% end
sigPeak_ln = log(sigPeak/sigPeak(1));

% curv fitting
coefficients = polyfit(grad_square, sigPeak_ln, 1);
y_fit = polyval(coefficients, grad_square);

% calculate gradient amplitude coefficient
G_coeff = sqrt(-coefficients(1) / (gamma^2*tau^2*(Delta-tau/3)*D));

% plot curve
plot(grad_square, sigPeak_ln, 'o','LineWidth',1.5); hold on;
plot(grad_square, y_fit, '-','LineWidth',1.5); 
xlabel('square of gradient amplitude, a.u.');
ylabel('signal intensity, a.u.');
legend('experiment','fitting curve');
title(append('Gmax=',num2str(round(G_coeff*100, 2)),' G/cm'));
set(gca,'FontSize',13,'FontName','Arial');

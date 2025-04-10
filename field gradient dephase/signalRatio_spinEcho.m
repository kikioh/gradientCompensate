% calculate the signal senstivity ratio for pulsed gradient spin
% echo (PGSE) experiment under gradient coupling
%
% Ref: 200 and more book, page 467, Exp. 11.5
%
% Mengjia He, 2024.09.06

%% formula model
clearvars; clc; close all;

% experiment parameters
ls = 6.5e-3;                    % sample length, m
Gmax = 95e-4*1e2;               % gradient strength at maximum current 3 A, T/m
tau = 1e-3;                     % gradient duration, s
grad_inter = 2*tau/pi;
gamma = spin('1H');             % observe nuclei
coupleRatio = 0.5*0.67*1e-2;    % gradient coupling ratio

% quantum number of coherence
quantumNum = [1,-1];

% set gradient amplitude
gradCoeff = linspace(0.05,0.95,19);
grad1 = [0,0];                  % primary gradient
grad2 = [1,-1];                 % coupled gradient

% calculate signal dephase under coupled gradient
ratio = zeros(1,numel(gradCoeff));
for m = 1:numel(gradCoeff)
    
    % sequence gradient ratio
    grad =  gradCoeff(m) *  (grad1 - grad2*coupleRatio);
    
    % sum of phase
    phase = Gmax * gamma * grad_inter * dot(quantumNum,grad);
    
    % wavelength of helix of the phase of the coherence
    lambda = 2*pi/phase;

    % signal intensity ratio
    ratio(m) = lambda/(pi*ls) * sin(pi*ls/lambda);

end
% plot ratio
fig = figure('Name','Intensity of PGSE signal');
fig.Position = [200 100 450 350];
plot(gradCoeff,abs(ratio),'LineWidth',1.5);
set(gca,'FontName','Arial','FontSize',14);
xlabel('gradient amplitude');
ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex');



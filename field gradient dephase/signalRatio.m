% calculate the signal senstivity ratio for pulsed gradient spin
% echo experiment under gradient coupling
%
% Mengjia He, 2024.09.06

%% formula model
clearvars; close all; clc;


z = 6.5e-3;   % sample length
Gmax = 95e-4*1e2; % maximum gradient strength at 3 A   
tau = 1e-3;  % gradient duration
p = 1;  % quantum number
gama = spin('1H');

% quantum number of coherence
quantumNum = [1,-1];

% set gradient amplitude
gradY = [0.0001,0.001,0.005,linspace(0.01,0.1,10)];
gradZ = linspace(0.1,0.9,19);

% gradient induced wavelength of the helix
lamba = logspace(-2,2,1000)*z;

% calculate signal dephase under primary gradient
% lamba = 2*pi./(Gmax * gradY * (2*tau/pi)  * spin('1H') * p);

% calculate signal dephase under coupled gradient
lamba = 2*pi./(2 * (0.67*1e-2) * Gmax * gradZ * (2*tau/pi) * spin('1H') * p);

% signal senstivity ratio
% ratio = lamba/(1i*2*pi*z) .* (exp(1i*2*pi*z./lamba) - 1);
ratio = lamba/(pi*z) .* sin(pi*z./lamba);

% plot ratio
%% semilogx(lamba/z,abs(ratio),'LineWidth',1.5);
plot(gradZ,abs(ratio),'LineWidth',1.5);
set(gca,'FontName','Arial','FontSize',14);
% xlabel('$\lambda/l$','Interpreter','latex');
xlabel('gradient amplitude');
ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex');



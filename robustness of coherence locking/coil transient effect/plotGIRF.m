% plot the gradient impulse response function
%
% mengjia.he@kit.edu, 2024.04.19

clearvars; close all; clc;

% define phase delay 
phi = [0,5,10]*pi/180;
T0 = 1e-3;
t =  linspace(0,1.5*T0,100);

for m = 1:numel(phi)

    [~,gradCur] = GIRF(phi(m),12,100);

    plot(t*1e3,gradCur,'LineWidth',1); hold on
end

legend('\phi=0^\circ','\phi=5^\circ','\phi=10^\circ');
xlabel('time, ms');
ylabel('current, A');
set(gca,'FontName','Airal','FontSize',14);

% test the output of pulseResponse function under different pulse shape
%
% Mengjia He, 2023.08.04

clearvars;
close all;
clc;

% % coil parameters for helmholtz coil
% coil.R = 0.00083776;	
% coil.L = 0.000000038219;
% coil.M = 0.0000000014762;

% coil parameters for hiscore coil
coil.R = 5.1103e-4;
coil.L = 1.3318e-8;
coil.M = 1.5578e-11;

% pulse parameters
tSeq = linspace(0,1e-3,100+1);
pulse.duration = 1e-3;
pulse.pulseType = 'sin';
pulse.numSlice = 100;
pulse.amplitude = 3;
pulse.sd = 1.5e-4;
pulse.ST = 2.5e-4;
% pulseResponse('trape',pulsePara,tseq); hold on;
[i2,i1]=pulseResponse(pulse,coil);

% plot excite current
figure;
plot(tSeq*1e3, i1,'LineWidth',1.5,'Color',[212,147,151]/255);
xlabel('time, ms');
ylabel('current, A');
grid on;
set(gca,'Fontsize',13,'Fontname', 'Airal');

% plot induced current
figure;
plot(tSeq*1e3, i2*1e3,'LineWidth',1.5,'Color',[23,82,143]/255);
xlabel('time, ms');
ylabel('current, mA');
grid on;
set(gca,'Fontsize',13,'Fontname', 'Airal');
% read complex fid from topspin exported JCAMP
%
% Mengjia He, 2023.12.20

close all; clearvars; clc;

% data with cp pulse
[fid_ref,AQ] = readJCAMP2D('fid-20240923-34.dx'); 
 
% obtain acqisition parameters
TD = size(fid_ref);
Fs = TD./AQ;

%% fouier transform
numFFT1 = TD(1); numFFT2 = TD(2); 
freqRef = [125e6,500e6];
f1 = linspace(0,Fs(1),numFFT1)/freqRef(1)*1e6;
f2 = linspace(0,Fs(2),numFFT2)/freqRef(2)*1e6;
spectrum_ref = abs(fftshift(fft2(fid_ref)))/1e10;

% plotting
close all;
fig1 = figure('Name','2D spectrum');
fig1.Position = [200 100 1200 350];
parameters.delta = [0.088 1 0.05 1.0]; 
parameters.m = 6; 
plot_2dSpec(f2,f1,spectrum_ref,parameters);
% % plot 1H spectrum
% spectrum_F2 = sum(spectrum, 1);
% fig2 = figure('Name','1H spectrum');
% plot(f2,transpose(abs(spectrum_F2)));
% set(gca,'XDir','reverse');
% % plot 13C spectrum
% spectrum_F1 = sum(spectrum, 2);
% fig3 = figure('Name','13C spectrum');
% plot(f1,transpose(abs(spectrum_F1)));
% set(gca,'XDir','reverse');

%% calculate time domain sum

% fid_ref = fid_ref/1e10;
% fidSum_ref = sum(sum(abs(fid_ref)));
% fid_cp = fid_cp/1e10;
% fidSum_cp = sum(sum(abs(fid_cp)));
% fid_sl = fid_sl/1e10;
% fidSum_sl = sum(sum(abs(fid_sl)));
% 
% disp(fidSum_cp/fidSum_ref);
% 
% disp(fidSum_sl/fidSum_ref);

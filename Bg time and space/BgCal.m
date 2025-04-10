% calculate the gradient field for the parallel detector
%
% Mengjia He, 2024.02.20

close all; clearvars; clc;
%% define sample and coil

% define channel and nuclei 
num = 2;                                    % number of sample 

% create sample for field import  
options.nz = 40;
sample=createSample(num,[],'gradient',options);

% define gradient coil array;
coil=createCoil(num,[],'gradient');

% comsol gradient solution
comsol = comsolSolInf(num,[],'gradient');

%% calculate pulse gradient field

% gSeq = ones(num,1);
% options.shape = 'trapezoid';
% gPulse=createGpulse(num,gSeq,options);
% 
% [Bg,Bg_self,Bg_couple] = BgDirectArray(coil,sample,BgCoeff,gPulse);

%% calculate sequenced gradient field, for HSQC simulation

% HSQC sequence
gSeq = [2,2,-1;2,2,1]/2*0.85;
% gSeq = [0,20,-5.03;0,90,22.64]/100;

% import gradient simulation
model_GD=mphopen('C2_Hiscore_GDandRF');
BgCoeff = BgCoeffArray(coil,sample,model_GD,comsol.GD);

options.shape = 'sin';
options.numSlice = 500;
gPulse=createGpulse(num,gSeq,options);
[p500.BgSeq,p500.Bg_selfSeq,p500.Bg_coupleSeq] = BgDirectArray(coil,sample,BgCoeff,gPulse);

options.numSlice = 1000;
gPulse=createGpulse(num,gSeq,options);
[p1000.BgSeq,p1000.Bg_selfSeq,p1000.Bg_coupleSeq] = BgDirectArray(coil,sample,BgCoeff,gPulse);

%% calculate sequenced gradient field, for HMQC simulation

% HSQC sequence
% gSeq = [2,2,-1;2,2,1]/2*0.85;
% gSeq = [0,20,-5.03;0,90,22.64]/100;

% % import gradient simulation
% model_GD=mphopen('C2_Hiscore_GDandRF_editedCouple');
% BgCoeff = BgCoeffArray(coil,sample,model_GD,comsol.GD);
% 
% options.shape = 'sin';
% options.numSlice = 500;
% gPulse=createGpulse(num,gSeq,options);
% [p500.BgSeq,p500.Bg_selfSeq,p500.Bg_coupleSeq] = BgDirectArray(coil,sample,BgCoeff,gPulse);
% 
% options.numSlice = 1000;
% gPulse=createGpulse(num,gSeq,options);
% [p1000.BgSeq,p1000.Bg_selfSeq,p1000.Bg_coupleSeq] = BgDirectArray(coil,sample,BgCoeff,gPulse);

%% calculate multiple pulse gradient field, for optimization

% % set coeff to 1.22 for 0.5G
% gSeq = 0.73*ones(num,1);
% 
% % import gradient simulation
% model_GD=mphopen('C2_Hiscore_GDandRF');
% BgCoeff = BgCoeffArray(coil,sample,model_GD,comsol.GD);
% 
% B_couple = cell(1,10);
% numSlice = linspace(200,2000,10);
% options.duration = 1e-3;
% options.shape = 'sin';
% for m = 1:numel(numSlice)
%     options.numSlice = numSlice(m);
% 
%     gPulse=createGpulse(num,gSeq,options);
% 
%     [~,~,B_couple{m}] = BgDirectArray(coil,sample,BgCoeff,gPulse);
% end
% 
% % package Bg field
% Bg = struct('numP', [], 'numVox', [],'field', []);
% for m = 1:numel(numSlice)
%     B_temp = B_couple{m}(1,:);
%     Bg(m).numVox = numel(B_temp);
%     Bg(m).numP = numSlice(m);
%     Bg(m).field = B_temp;
% end

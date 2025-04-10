% simulate the gradient field as a function of time and space
%
% Mengjia He, 2024.03.06

clearvars; close all; clc;

%% define parameters for optimization
gpulse = struct; 
numGpulse = 2;
amp_gpulse = [3 0];
numSlice = 100;
numVox = 20;
l_sample = 8;

% define sample array
num = 2;
sample.num = num;                                   % number of sample
sample.nr = 1; sample.np = 1; sample.nz = numVox;   % division number
sample.axis = 'z';                                  % axis of cylinder sample

% define gradient coil array;
coilName = 'HiscoreGradient.mat';
coil.num = sample.num*2;                        % number of coils
coil.freq = 500e6; coil.Z0 = 50;
coil.R=load('HiscoreGradient.mat').R;
coil.L=load('HiscoreGradient.mat').L;
coil.M=load('HiscoreGradient.mat').M;

% import COMSOL module
comsol_G.objname = ["pt1","pt2"];   % object name of mark point in COMSOL
comsol_G.dset = 'dset2';
model_G=mphopen('C2_Hiscore_GDandRF');

% define gradient pulse
for n = 1:numGpulse
    gpulse(n).numPulse = numGpulse;
    gpulse(n).pulseType = 'sin';
    gpulse(n).duration = 1e-3;
    gpulse(n).amplitude = amp_gpulse(n);
    gpulse(n).numSlice = numSlice;
    gpulse(n).time_grad = gpulse(n).duration/gpulse(n).numSlice...
        * ones(1,gpulse(n).numSlice);
end
% calculate Bg
[Bg,Bg_induced,indCurr] = BgArray(coil,sample,model_G,comsol_G,gpulse);

%% plot direct and induced coupled gradient field
figure;
Bg_ts = zeros(numVox,numSlice);
Bg_induced_ts = zeros(numVox,numSlice);
parfor m = 1:numVox
    Bg_ts(m,:) = Bg{2,m};
    Bg_induced_ts(m,:) = Bg_induced{2,m};
end
t_axis = linspace(0,1,numSlice);
z_axis = linspace(-l_sample/2,l_sample/2,numVox);
[T,Z] = meshgrid(t_axis,z_axis);
surf(T,Z,1e6*abs(Bg_ts),'FaceAlpha',0.3); hold on
shading interp;
colormap(hmj_cmap('RdBu'));
h = surf(T,Z,1e6*abs(Bg_induced_ts)); hold on
shading interp;
xlabel('time, ms');
ylabel('z, mm');
zlabel('Bg amplitude, uT');
colorbar;
set(gca,'Fontsize',13,'Fontname', 'Airal');
xlim([-4,4]);
ylim([0,1]);

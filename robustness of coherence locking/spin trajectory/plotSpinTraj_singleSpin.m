% plot the trajectory of spin coherence, suppose the nuclei is on
% resonance, and the B1 is uniform
%
% Mengjia He, 2024.02.29

clearvars; close all; clc;

% load spin locking pulse
% data = load('SL-sin-1H-1ms-RF6kHz-BW7kHz-05G.mat');
data = load('SL-sin-13C-1ms-RF4kHzdev17-BW6kHz-05G.mat');
pulse_dt = data.pulse_dt;
para = data.para;
numP = para.numP;
numVox = para.grad.numVox;
Bg = para.grad.Bg;

% Construct Cartesian pulse waveform
[Cx,Cy]=Ap2Xy(2*pi*data.para.RF.A,data.pulse_shape);

%% Spinach parameters
% 1T magnet
sys.magnet = 1;                                 
% sys.enable={'greedy'};

% Isotopes
sys.isotopes = {'1H','13C'};

% Zeeman Interactions
inter.zeeman.scalar={0,0};                          % chemical shift, ppm

% Spin Interactions
inter.couple.scalar{1,2}=145;                       % J coupling, Hz
inter.couple.scalar{2,2}=0;                         % J coupling, Hz

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho=state(spin_system,'L+',para.nuclei);
rho=rho/norm(full(rho),'fro');

% Get the control operators
Lp=operator(spin_system,'L+',para.nuclei); 
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman'));      % Zeeman part
Hc = hamiltonian(assume(spin_system,'nmr','coupling'));         % coupling part

%% spin trajectory
rho_comp = cell(numVox,1);  rho_uncomp = cell(numVox,1); 
rhoTraj_comp = cell(numVox,1);  rhoTraj_uncomp = cell(numVox,1);
for k = 1:numVox  % voxel loop

    % uncompensated
    [rho_uncomp{k},rhoTraj_uncomp{k}]=shaped_pulse_xy(spin_system,Hc,{Hz},Bg(k),pulse_dt,rho,'expv-pwc');

    % compensated
    [rho_comp{k},rhoTraj_comp{k}]=shaped_pulse_xy(spin_system,Hc,{Lx,Ly,Hz},{Cx{1},Cy{1},Bg{k}},pulse_dt,rho,'expv-pwc');
end

% calculate fidelity
fid_comp = zeros(numVox,numP+1); fid_uncomp = zeros(numVox,numP+1);
for m = 1:numP+1
    parfor k = 1:numVox
        fid_comp(k,m)  = rho'*rhoTraj_comp{k}{m};
        fid_uncomp(k,m)  =  rho'*rhoTraj_uncomp{k}{m};
    end
end

%% plot trajectory
close all;
x_axis = linspace(0,sum(pulse_dt),numP+1);
fig = figure('Name','fidelity of J=145 Hz, ignore v0 and v1 offsets');
fig.Position = [200 100 250 350];
subplot(2,1,1)
fid_comp_min = min(real(fid_comp),[],1);
fid_comp_max = max(real(fid_comp),[],1);  
idx = fid_comp_max >= fid_comp_min;
x_fill = [x_axis(idx), fliplr(x_axis(idx))];
y_fill = [fid_comp_max(idx), fliplr(fid_comp_min(idx))];
fill(x_fill*1e3, y_fill, lineColor('matlab1'), 'EdgeColor', lineColor('matlab1'),...
    'LineWidth', 0.1,'FaceAlpha', 0.3);
xlabel('time, ms'); ylabel('<\rho|{\itS^+}>, SL');
xlim([0,1]); ylim([-1,1]);
box on;
set(gca,'FontSize',10,'FontName','Arial','Linewidth',0.4);

subplot(2,1,2)
fid_uncomp_min = min(real(fid_uncomp),[],1); 
fid_uncomp_max = max(real(fid_uncomp),[],1);  
idx = fid_uncomp_max >= fid_uncomp_min;
x_fill = [x_axis(idx), fliplr(x_axis(idx))];
y_fill = [fid_uncomp_max(idx), fliplr(fid_uncomp_min(idx))];
fill(x_fill*1e3, y_fill, lineColor('matlab2'), 'EdgeColor', lineColor('matlab2'),...
    'LineWidth', 0.1,'FaceAlpha', 0.3);
xlabel('time, ms'); ylabel('<\rho|{\itS^+}>, no SL');
xlim([0,1]); ylim([-1,1]);
box on;
set(gca,'FontSize',10,'FontName','Arial','Linewidth',0.4);


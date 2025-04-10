% spin locking robustness to coil transient effect
%
% Mengjia He, 2024.05.15

clearvars; close all; clc;

% load single spin pulse
% data = load('SL-sin-13C-1ms-RF10kHz-BW25kHz-1G.mat');
data = load('SL-sin-1H-1ms-RF6kHz-BW7kHz-05G.mat');
pulse_shape = data.pulse_shape;
pulse_dt = data.pulse_dt;
para = data.para;
numP = para.numP;
numVox = para.grad.numVox;

% set Bg field
numPhase = 11;
phi = linspace(0,40/180*pi,numPhase);

% Construct Cartesian pulse waveform
[Cx,Cy]=Ap2Xy(2*pi*para.RF.A,pulse_shape);

%% Spinach parameters
sys.magnet = 1;                                  
sys.isotopes = {'1H'};
inter.zeeman.scalar={0};      

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho=state(spin_system,'L-','all');
rho=rho/norm(full(rho),'fro');

% Get the control operators
Lp=operator(spin_system,'L+','all'); 
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman')); 
H = sparse(zeros(size(Hz)));

%% fidelity with Bg and phase 

eta = zeros(numPhase,numVox); 
for k = 1:numPhase                                      % phase offset
    Bg = GIRF(phi(k),numVox,numP);
    parfor s = 1:numVox                                 % Bg

        eta(k,s)  = rho'*shaped_pulse_xy(spin_system,H,{Lx,Ly,Hz},...
            {Cx{1},Cy{1},Bg(s,:)},pulse_dt,rho);
    end
end

%% plot fidelity
fig = figure('Name','fidelity with Bg and phase');
fig.Position = [200 100 400 300];
parameters.xlabel='\phi, degree';
parameters.ylabel='B_0 drift, uT';
parameters.ylim=[-25,25];
parameters.clim = [0.95,1];
Bg = GIRF(0,numVox,numP);
plot_pd(phi*180/pi,Bg(:,numP/2)*1e6,real(transpose(eta)),parameters);

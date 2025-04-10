% plot the trajectory of spin coherence
%
% Mengjia He, 2024.02.29

clearvars; close all; clc;

% experiment name
addpath('Q:\OneDrive\Digital NMR\Matlab\gradient compensation\pulse shape');
expName = 'SL-sin-SpAndLzSp-1ms-H30kHz-C20kHz-05G';

% load spin locking pulse
data = load(append(expName,'.mat'));
pulse_shape = data.pulse_shape;
pulse_dt = data.pulse_dt;

% load experiments parameters
para = data.para;
numP = para.numP;

% sweep B1 field
RFA = 2*pi*[para.RF.AH;para.RF.AC];
numRF = para.RF.numRF;
RFAH = 2*pi*para.RF.AH*linspace(1-para.RF.devH,1+para.RF.devH,para.RF.numRF);
RFAC = 2*pi*para.RF.AC*linspace(1-para.RF.devC,1+para.RF.devC,para.RF.numRF);

% sweep J coupling
JHC = linspace(0,250,11); numJ = numel(JHC);

% sweep resonance offset
numH = para.spin.numH*2;
BWH = linspace(-para.spin.BWH/2,para.spin.BWH/2,numH);
numC = para.spin.numC*2;
BWC = linspace(-para.spin.BWC/2,para.spin.BWC/2,numC);

% sweep Bg field
numVox = para.grad.numVox;
Bg = para.grad.Bg;

% Construct Cartesian pulse waveform
[Cx,Cy]=Ap2Xy(RFA,pulse_shape);

%% Spinach parameters
% 500 MHz magnet
sys.magnet = 1;                                     % magnet for Hamiltonian
sys.enable={'greedy'};

% % Isotopes
% sys.isotopes = {'1H','1H','13C','13C'};
% nuclei = unique(sys.isotopes,'stable');
% numIsotope = numel(nuclei);
% % chemical shift, ppm
% inter.zeeman.scalar={1,2,40,50};
% 
% % J coupling, Hz
% inter.coupling.scalar= cell(4);
% inter.coupling.scalar{1,2} = 13;
% inter.coupling.scalar{1,3} = 140;
% inter.coupling.scalar{2,3} = 156;
% inter.coupling.scalar{1,4} = 155;
% inter.coupling.scalar{2,4} = 145;
% inter.coupling.scalar{3,4} = 50;

% Spin system
sys.isotopes={'15N','1H','13C','13C','13C','15N'};

% Chemical shifts, ppm
inter.zeeman.scalar={0.0, 0.0, 55.0, 30.0, 175.0, 0.0};

% Scalar couplings, Hz (literature values)
inter.coupling.scalar=cell(6);
inter.coupling.scalar{1,3}=-400; 
inter.coupling.scalar{2,3}=140;
inter.coupling.scalar{3,4}=35;
inter.coupling.scalar{3,5}=55;
inter.coupling.scalar{3,6}=7;
inter.coupling.scalar{5,6}=-15;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=4;

JCN = linspace(-5000,5000,11);
fidelity = zeros(1,numel(JCN));
for m = 1:numel(JCN)
    inter.coupling.scalar{1,3}=JCN(m);


% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
Lm=state(spin_system,'L-','1H');
Lm=Lm/norm(full(Lm),'fro');

LmSz=state(spin_system,{'L-','Lz'},{2,3})+state(spin_system,{'L-','Lz'},{2,4})...
    +state(spin_system,{'L-','Lz'},{2,5});
LmSz=LmSz/norm(full(LmSz),'fro');

LzSm=state(spin_system,{'Lz','L-'},{2,3})+state(spin_system,{'Lz','L-'},{2,4})...
    +state(spin_system,{'Lz','L-'},{2,5});
LzSm=LzSm/norm(full(LzSm),'fro');

LzSp=state(spin_system,{'Lz','L+'},{2,3})+state(spin_system,{'Lz','L+'},{2,4})...
    +state(spin_system,{'Lz','L+'},{2,5});
LzSp=LzSp/norm(full(LzSp),'fro');

% Get the control operators
LpH=operator(spin_system,'L+','1H');
LxH=(LpH+LpH')/2; LyH=(LpH-LpH')/2i;
LzH = operator(spin_system,'Lz','1H');

LpC=operator(spin_system,'L+','13C');
LxC=(LpC+LpC')/2; LyC=(LpC-LpC')/2i;
LzC = operator(spin_system,'Lz','13C');

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman'));      % Zeeman part
Hc = hamiltonian(assume(spin_system,'nmr','coupling'));         % couple part

% tagart spin state
rho_init = LzSp;

% spin trajectory
rho=shaped_pulse_xy(spin_system,Hc,{LxH,LyH,LxC,LyC,Hz},...
    {Cx{1},Cy{1},Cx{2},Cy{2},Bg{2}},pulse_dt,rho_init);
% traj = cell2mat(traj);
% time_axis = linspace(0,sum(pulse_dt),numP+1);
% 
fidelity(m) = real(rho_init'*rho);
end
% figure
% trajan(spin_system,traj,'coherence_order',time_axis);
% trajan(spin_system,traj,'correlation_order',time_axis);
% trajan(spin_system,traj,'local_each_spin',time_axis);
% trajan(spin_system,traj,'total_each_spin',time_axis);
% disp(fidelity);
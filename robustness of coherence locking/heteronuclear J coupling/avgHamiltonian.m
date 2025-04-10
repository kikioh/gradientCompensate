% average Hamiltonian theory, vislalize the averaged J coupling
%
% Mengjia He, 2024.05.13
clearvars; close all; clc;

% experiment name
addpath('Q:\OneDrive\Digital NMR\Matlab\gradient compensation\pulse shape');
expName = 'SL-sin-SmAndLzSm-1ms-H30kHz-C20kHz-05G';

% load spin locking pulse
data = load(append(expName,'.mat'));
pulse_shape = data.pulse_shape;
pulse_dt = data.pulse_dt;
para = data.para;
numP = para.numP;
numVox = para.grad.numVox;
Bg = para.grad.Bg;

% sweep B1 field
RFA = 2*pi*[para.RF.AH;para.RF.AC];
numRF = para.RF.numRF;
RFAH = 2*pi*para.RF.AH*linspace(1-para.RF.devH,1+para.RF.devH,para.RF.numRF);
RFAC = 2*pi*para.RF.AC*linspace(1-para.RF.devC,1+para.RF.devC,para.RF.numRF);

% sweep J coupling
JHC = linspace(0,250,5); numJ = numel(JHC);

% Construct Cartesian pulse waveform
[Cx,Cy]=Ap2Xy(RFA,pulse_shape);

%% generate new spin system

% 500 MHz magnet
sys.magnet = 1;                

% isotopes
sys.isotopes = {'1H','13C','15N'};

% Interactions
inter.zeeman.scalar={0,0,0};    

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
LzSm = state(spin_system,{'Lz','L-'},{1,2});
LmSz = state(spin_system,{'L-','Lz'},{1,2});
rho_init = LzSm; rho_init=rho_init/norm(full(rho_init),'fro');

% Get the control operators
LpH=operator(spin_system,'L+','1H'); 
LxH=(LpH+LpH')/2; LyH=(LpH-LpH')/2i;
LzH = operator(spin_system,'Lz','1H');

LpC=operator(spin_system,'L+','13C'); 
LxC=(LpC+LpC')/2; LyC=(LpC-LpC')/2i;
LzC = operator(spin_system,'Lz','13C');
LzHC = operator(spin_system,{'Lz','Lz'},{1,2});

LpN=operator(spin_system,'L+','15N'); 
LxN=(LpN+LpN')/2; LyN=(LpN-LpN')/2i;
LzN = operator(spin_system,'Lz','15N');

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman'));    % Zeeman part
Hc = hamiltonian(assume(spin_system,'labframe','couplings')); % coupling part
H = Hz+Hc;

% calculate propagator
propH = cell(1,numP);
propC = cell(1,numP);
propN = cell(1,numP);
parfor q = 1:numP

    [~,~,propH{q}] = shaped_pulse_xy(spin_system,H,{LxH,LyH,LzH},...
        {Cx{1}(1:q),Cy{1}(1:q),Bg{1}(1:q)*spin('1H')},pulse_dt(1:q),rho_init);

    [~,~,propC{q}] = shaped_pulse_xy(spin_system,H,{LxC,LyC,LzC},...
        {Cx{2}(1:q),Cy{2}(1:q),Bg{1}(1:q)*spin('13C')},pulse_dt(1:q),rho_init);

    [~,~,propN{q}] = shaped_pulse_xy(spin_system,H,{LzN},...
        {Bg{1}(1:q)*spin('15N')},pulse_dt(1:q),rho_init);

end

% calculate rotation cofficient
alphaH = cell(3,3); alphaC = cell(3,3); alphaN = cell(3,3); 
TensorHC = cell(3,3);   TensorCN = cell(3,3);
for m = 1:3
    for n = 1:3
        alphaH{m,n} = zeros(1,numP);
        alphaC{m,n} = zeros(1,numP);
        alphaN{m,n} = zeros(1,numP);
        TensorHC{m,n} = zeros(1,numP);
        TensorCN{m,n} = zeros(1,numP);
    end
end

for q = 1:numP
    
    % rotation cofficient for 1H
    alphaH{1,1}(q) = trace(LxH'*(propH{q}'*LxH*propH{q}));
    alphaH{1,2}(q) = trace(LyH'*(propH{q}'*LxH*propH{q}));
    alphaH{1,3}(q) = trace(LzH'*(propH{q}'*LxH*propH{q}));

    alphaH{2,1}(q) = trace(LxH'*(propH{q}'*LyH*propH{q}));
    alphaH{2,2}(q) = trace(LyH'*(propH{q}'*LyH*propH{q}));
    alphaH{2,3}(q) = trace(LzH'*(propH{q}'*LyH*propH{q}));

    alphaH{3,1}(q) = trace(LxH'*(propH{q}'*LzH*propH{q}));
    alphaH{3,2}(q) = trace(LyH'*(propH{q}'*LzH*propH{q}));
    alphaH{3,3}(q) = trace(LzH'*(propH{q}'*LzH*propH{q}));
    
    % rotation cofficient for 13C
    alphaC{1,1}(q) = trace(LxC'*(propC{q}'*LxC*propC{q}));
    alphaC{1,2}(q) = trace(LyC'*(propC{q}'*LxC*propC{q}));
    alphaC{1,3}(q) = trace(LzC'*(propC{q}'*LxC*propC{q}));

    alphaC{2,1}(q) = trace(LxC'*(propC{q}'*LyC*propC{q}));
    alphaC{2,2}(q) = trace(LyC'*(propC{q}'*LyC*propC{q}));
    alphaC{2,3}(q) = trace(LzC'*(propC{q}'*LyC*propC{q}));

    alphaC{3,1}(q) = trace(LxC'*(propC{q}'*LzC*propC{q}));
    alphaC{3,2}(q) = trace(LyC'*(propC{q}'*LzC*propC{q}));
    alphaC{3,3}(q) = trace(LzC'*(propC{q}'*LzC*propC{q}));
    

    alphaN{1,1}(q) = trace(LxN'*(propN{q}'*LxN*propN{q}));
    alphaN{1,2}(q) = trace(LyN'*(propN{q}'*LxN*propN{q}));
    alphaN{1,3}(q) = trace(LzN'*(propN{q}'*LxN*propN{q}));

    alphaN{2,1}(q) = trace(LxN'*(propN{q}'*LyN*propN{q}));
    alphaN{2,2}(q) = trace(LyN'*(propN{q}'*LyN*propN{q}));
    alphaN{2,3}(q) = trace(LzN'*(propN{q}'*LyN*propN{q}));

    alphaN{3,1}(q) = trace(LxN'*(propN{q}'*LzN*propN{q}));
    alphaN{3,2}(q) = trace(LyN'*(propN{q}'*LzN*propN{q}));
    alphaN{3,3}(q) = trace(LzN'*(propN{q}'*LzN*propN{q}));
end

% calculate the coupling tensor
CoupleAvgHC = zeros(3,3);   CoupleAvgCN = zeros(3,3);
for m = 1:3
    for n = 1:3
        % time-dependent coupling tensor
        TensorHC{m,n} = alphaH{1,m}.*alphaC{1,n} + alphaH{2,m}.*alphaC{2,n} +...
            alphaH{3,m}.*alphaC{3,n};
        
        TensorCN{m,n} = alphaC{1,m}.*alphaN{1,n} + alphaC{2,m}.*alphaN{2,n} +...
            alphaC{3,m}.*alphaN{3,n};

        % average coupling tensor
        CoupleAvgHC(m,n) = mean(TensorHC{m,n});
        CoupleAvgCN(m,n) = mean(TensorCN{m,n});
    end
end

% apply pusle and test J-coupling impact
rho= cell(1,numVox);  rho_J= cell(1,numVox);
proj_J2real = zeros(numJ,numVox);       % projection of final to only J coupled final state
proj = zeros(numJ,numVox);              % projection of final to initial state

for k=1:numVox % voxel loop
    for q = 1:numJ
        Hc=2*pi*JHC(q)*LzHC;
        
        % case 1: RF+J+CS +gradient
        rho{k}=shaped_pulse_xy(spin_system,Hc,{LxH,LyH,LxC,LyC,Hz},...
            {Cx{1},Cy{1},Cx{2},Cy{2},Bg{k}},pulse_dt,rho_init);
        
        % case 2: only J
        rho_J{k} = evolution(spin_system,Hc,[],rho_init,numP,pulse_dt(1),'final');     
        
        proj_J2real(q,k) = trace(rho_J{k}'*rho{k});
        proj(q,k) = trace(rho_init'*rho{k});

    end
end

% visualize
figure();
hold on;
for n = 1:numJ
    plot(real(proj_J2real(n,:))); 
end
ylabel('<\rho_J|\rho_{real}>');
xlabel('voxel number');
legend('J=0','J=62.5','J=125','J=187.5','J=250');

% visualize
figure();
hold on;
for n = 1:numJ
    plot(real(proj(n,:))); 
end
xlabel('voxel number');
ylabel('<\rho_0|\rho_{real}>');
legend('J=0','J=62.5','J=125','J=187.5','J=250');

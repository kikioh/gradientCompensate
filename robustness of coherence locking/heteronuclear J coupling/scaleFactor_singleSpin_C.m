% calculate the J coupling scale factor as a function of v0 and v1
%
% Mengjia He, 2024.05.15

clearvars; close all; clc;

% load spin locking pulse
data = load('SL-sin-13C-1ms-RF4kHzdev17-BW6kHz-05G.mat');
pulse_dt = data.pulse_dt;
pulse_shape = data.pulse_shape;
para = data.para;
numP = para.numP;
numRF = 5;
numBW = 10;
para.RF.dev = 0.15;
RFA = 2*pi*para.RF.A*linspace(1-para.RF.dev,1+para.RF.dev,numRF);
omega = 2*pi*linspace(-para.spin.BW/2,para.spin.BW/2,numBW);

% generate spin system
sys.magnet = 1;
sys.isotopes = {'1H','13C'};
inter.zeeman.scalar={0,0};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init = state(spin_system,{'Lz','L-'},{1,2});
rho_init=rho_init/norm(full(rho_init),'fro');

% Get the control operators
LpH=operator(spin_system,'L+','1H');
LxH=(LpH+LpH')/2; LyH=(LpH-LpH')/2i;
LzH = operator(spin_system,'Lz','1H');

LpC=operator(spin_system,'L+','13C');
LxC=(LpC+LpC')/2; LyC=(LpC-LpC')/2i;
LzC = operator(spin_system,'Lz','13C');

% Calculate the Hamiltonian
H = hamiltonian(assume(spin_system,'nmr'));

% calculate scale factor
scaleMatrix = cell(numBW,numRF);
scaleFactor = zeros(numBW,numRF);
for s =1:numBW
    for k=1:numRF
        
        [Cx,Cy]=Ap2Xy(RFA(k),pulse_shape);

        % calculate propagator
        propH = cell(1,numP); propC = cell(1,numP);
        parfor q = 1:numP

             propH{q} = eye(size(H,1));

             [~,~,propC{q}] = shaped_pulse_xy(spin_system,H,{LxC,LyC,LzC},...
                {Cx{1}(1:q),Cy{1}(1:q),omega(s)*ones(1,q)},pulse_dt(1:q),rho_init);
        end

        % calculate rotation cofficient
        alphaH = cell(3,3); alphaC = cell(3,3);
        for m = 1:3
            for n = 1:3
                alphaH{m,n} = zeros(1,numP);
                alphaC{m,n} = zeros(1,numP);
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

        end

        % calculate the coupling tensor
        scaleMatrix{s,k} = zeros(3,3);
        for m = 1:3
            for n = 1:3
                 scaleMatrix{s,k}(m,n) = mean(alphaH{3,m}.*alphaC{3,n});
            end
        end
        scaleFactor(s,k) = norm(scaleMatrix{s,k},'fro');
    end
end

%% plot scale factor
close all;
fig = figure('name',para.name);
fig.Position = [200 100 320 300];
imagesc(omega/(2*pi*1e3),RFA/(2*pi*1e3),transpose(scaleFactor));
xlabel('13C \nu_0, kHz');
ylabel('13C \nu_1, kHz');
yticks([3.4 3.7 4.0 4.3 4.6]);
xticks([-3.0 -1.5 0 1.5 3.0]);
colormap('hot'); 
colorbar;
clim([0,0.1]);
set(gca,'YDir','normal');
set(gca,'Fontsize',11,'Fontname', 'Airal');
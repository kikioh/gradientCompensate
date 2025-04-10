% spin locking robustness to v1 and v0
%
% Mengjia He, 2024.05.15

clearvars; close all; clc;
% load pulse
% data = load('SL-sin-SmAndLzSm-1ms-H30kHz-C15kHz-05G.mat');
% pulse_shape = data.pulse_shape;
% pulse_dt = data.pulse_dt;
% para = data.para;
% numP = para.numP;

% load single spin pulse
dataH = load('SL-sin-1H-1ms-RF7kHz-BW8kHz-05G.mat');
dataC = load('SL-sin-13C-1ms-RF10kHz-BW25kHz-1G.mat');
pulse_shape = [dataH.pulse_shape;dataC.pulse_shape];
pulse_dt = dataH.pulse_dt;
numP = dataH.para.numP;
paraH = dataH.para;
paraC = dataC.para;



% set Bg field
numVox = paraH.grad.numVox;
Bg = paraH.grad.Bg;

% set v1 amplitude
RFA = 2*pi*[paraH.RF.A;paraC.RF.A];
numRF = paraH.RF.numRF;
RFAH = 2*pi*paraH.RF.A*linspace(1-paraH.RF.dev,1+paraH.RF.dev,paraH.RF.numRF);    
RFAC = 2*pi*paraC.RF.A*linspace(1-paraC.RF.dev,1+paraC.RF.dev,paraC.RF.numRF);

% sweep J coupling
JHC = linspace(0,250,11); numJ = numel(JHC);

% % sweep resonance offset
numH = paraH.spin.num;
BWH = linspace(-paraH.spin.BW/2,paraH.spin.BW/2,numH);
% numC = para.spin.numC*2;
% BWC = linspace(-para.spin.BWC/2,para.spin.BWC/2,numC);

% % sweep Bg field
% numVox = para.grad.numVox;
% Bg = para.grad.Bg;

% Construct Cartesian pulse waveform
[Cx,Cy]=Ap2Xy(RFA,pulse_shape);

%% Spinach parameters
sys.magnet = 1;                                  
sys.isotopes = {'1H','13C'};
inter.zeeman.scalar={0,0};                      

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,{'L-'},{1});
rho_init=rho_init/norm(full(rho_init),'fro');

% Get the control operators
LxH=operator(spin_system,'Lx','1H'); 
LyH=operator(spin_system,'Ly','1H'); 
LzH = operator(spin_system,'Lz','1H');

LxC=operator(spin_system,'Lx','13C'); 
LyC=operator(spin_system,'Ly','13C'); 
LzC = operator(spin_system,'Lz','13C');

LzHC = operator(spin_system,{'Lz','Lz'},{1,2});

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman')); 

%% fidelity with J and Bg
% fid_comp = zeros(numVox,numJ); fid_uncomp = zeros(numVox,numJ);
% for k = 1:numVox                                          % voxel loop
%     parfor s = 1:numJ                                      % J coupling loop
% 
%         H = 2*pi*JHC(s)*LzHC;
% 
%         % uncompensated
%         fid_uncomp(k,s) = rho_init'*shaped_pulse_xy(spin_system,H,{Hz},{Bg{k}},pulse_dt,rho_init);
% 
%         % compensated
%         fid_comp(k,s)  = rho_init'*shaped_pulse_xy(spin_system,H,{LxH,LyH,...
%             LxC,LyC,Hz},{Cx{1},Cy{1},Cx{2},Cy{2},Bg{k}},pulse_dt,rho_init);
%     end
% end
% 
% % plot fidelity
% fig = figure('Name','fidelity with J shifts, ignore v0 and v1 offsets');
% fig.Position = [200 100 1000 300];
% subplot(1,2,1)
% parameters.xlabel='J, Hz';
% parameters.ylabel='max Bg, uT';
% parameters.clim = [0.95,1];
% BGmat = cell2mat(transpose(Bg));
% plot_pd(JHC,BGmat(:,500)*1e6,real(fid_comp),parameters);
% 
% subplot(1,2,2)
% parameters.clim = [-1,1];
% plot_pd(JHC,BGmat(:,500)*1e6,real(fid_uncomp),parameters);

%% fidelity with v0 offsets
% fid_comp = zeros(numH,numC); fid_uncomp = zeros(numH,numC);
% Bg_temp = Bg{3}; 
% for k = 1:numH                                      % 1H offset
%     parfor s = 1:numC                               % 13C offset
% 
%         H = 2*pi*BWH(k)*LzH + 2*pi*BWC(s)*LzC;
% 
%         % uncompensate
%         fid_uncomp(k,s)  =  rho_init'*shaped_pulse_xy(spin_system,H,{Hz},{Bg_temp},pulse_dt,rho_init);
% 
%         % compensate
%         fid_comp(k,s)  = rho_init'*shaped_pulse_xy(spin_system,H,{LxH,LyH,...
%             LxC,LyC,Hz},{Cx{1},Cy{1},Cx{2},Cy{2},Bg_temp},pulse_dt,rho_init);
%     end
% end
% 
% % plot fidelity
% fig = figure('Name','fidelity with v0 offsets');
% fig.Position = [200 100 1000 300];
% subplot(1,2,1)
% parameters.xlabel='1H \nu_0, kHz';
% parameters.ylabel='13C \nu_0, kHz';
% % parameters.title='real(<\rho|I^+S^->), SL';
% parameters.clim = [0.95,1];
% plot_pd(BWH/1e3,BWC/1e3,real(transpose(fid_comp)),parameters);
% 
% subplot(1,2,2)
% parameters.clim = [-1,1];
% % parameters.title='real(<\rho|I^+S^->), no SL';
% plot_pd(BWH/1e3,BWC/1e3,real(transpose(fid_uncomp)),parameters);

%% fidelity within 1H channel, v0 and v1 offsets on 1H 
eta = zeros(numH,numRF);
unit_opt = eye(size(Hz,1));
Bg_temp = Bg{3}; J_temp = 143;         		% fix Bg and J       
for k = 1:numH                            	% 1H offset
    parfor s = 1:numRF                    	% RF offset

        coeffH = RFAH(s) / RFA(1);
        coeffC = RFAC(s) / RFA(2);
        H = 2*pi*J_temp*LzHC + 2*pi*BWH(k)*LzH;
        [~,~,prop_temp]  = shaped_pulse_xy(spin_system,H,...
            {LxH,LyH,LxC,LyC,Hz},{Cx{1}*coeffH,Cy{1}*coeffH,Cx{2}*coeffC,...
            Cy{2}*coeffC,Bg_temp},pulse_dt,rho_init,'expv-pwc');

        eta(k,s) = trace(unit_opt'*prop_temp)/size(unit_opt,1);
    end
end

%% plot fidelity
fig = figure('Name','fidelity with v0v1 offsets on 1H');
fig.Position = [200 100 350 350];
parameters.xlabel='1H \nu_0, kHz';
parameters.ylabel='1H \nu_1, kHz';
parameters.title=append('\eta=',num2str(mean(mean(real(eta)))));
parameters.clim = [0.95,1];
plot_pd(BWH/1e3,RFAH/(2*pi*1e3),real(transpose(eta)),parameters);

%% fidelity within 13C channel, v0 and v1 offsets on 13C
% fid_comp = zeros(numC,numRF); fid_uncomp = zeros(numC,numRF);
% Bg_temp = Bg{ceil(numVox/2)}; J_temp = JHC(ceil(numJ/2));         % fix Bg and J
% for k = 1:numC                            % 13C offset
%     parfor s = 1:numRF                    % RF offset                                                                      % RF offset
% 
%         coeffH = RFAH(s) / (2*pi*para.RF.AH);
%         coeffC = RFAC(s) / (2*pi*para.RF.AC);
%         H = 2*pi*J_temp*LzHC + 2*pi*BWC(k)*LzC;
% 
%         % uncompensated
%         fid_uncomp(k,s)  =  rho_init'*shaped_pulse_xy(spin_system,H,{Hz},{Bg_temp},pulse_dt,rho_init);
% 
%         % compensated
%         fid_comp(k,s)  = rho_init'*shaped_pulse_xy(spin_system,H,{LxH,LyH,LxC,LyC,Hz},...
%             {Cx{1}*coeffH,Cy{1}*coeffH,Cx{2}*coeffC,Cy{2}*coeffC,Bg_temp},pulse_dt,rho_init);
%     end
% end
% 
% % plot fidelity
% fig = figure('Name','fidelity with v0v1 offsets on 13C');
% fig.Position = [200 100 1000 300];
% subplot(1,2,1)
% parameters.xlabel='13C \nu_0, kHz';
% parameters.ylabel='13C \nu_1, kHz';
% % parameters.title='real(<\rho|I^+S^->), SL';
% parameters.clim = [0.95,1];
% plot_pd(BWC/1e3,RFAC/(2*pi*1e3),real(transpose(fid_comp)),parameters);
% 
% subplot(1,2,2)
% parameters.clim = [-1,1];
% plot_pd(BWC/1e3,RFAC/(2*pi*1e3),real(transpose(fid_uncomp)),parameters);


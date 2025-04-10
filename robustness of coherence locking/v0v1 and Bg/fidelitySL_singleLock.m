% spin locking robustness to v1 and v0data
%
% Mengjia He, 2024.05.15

clearvars; close all; clc;

% load single spin pulse
% data = load('SL-sin-1H-1ms-RF6kHz-BW7kHz-05G.mat');
data = load('SL-sin-13C-1ms-RF4kHzdev17-BW6kHz-05G.mat');

% read pulse parameters
pulse_shape = data.pulse_shape;
pulse_dt = data.pulse_dt;
para = data.para;

% v0 offset
numBW = 50;
BW = linspace(-para.spin.BW/2,para.spin.BW/2,numBW);

% set v1 offset
numRF = 10;
RFA = 2*pi*para.RF.A*linspace(1-para.RF.dev,1+para.RF.dev,numRF);  

% set Bg field
Bg =load('Bg-sin-05G-v21.mat').Bg(5).field;
numVox = numel(Bg);

% Construct Cartesian pulse waveform
[Cx,Cy]=Ap2Xy(2*pi*para.RF.A,pulse_shape);

%% Spinach parameters
sys.magnet = 1;                                  
sys.isotopes = {para.nuclei};
inter.zeeman.scalar={0};      

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,'L+','all');
rho_init=rho_init/norm(full(rho_init),'fro');

% Get the control operators
Lx=operator(spin_system,'Lx','all'); 
Ly=operator(spin_system,'Ly','all'); 
Lz = operator(spin_system,'Lz','all');

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman')); 

%% fidelity with v0 and v1 offsets
unit_opt = eye(size(Hz,1));
eta = zeros(numBW,numRF); Bg_temp = Bg{21};    
for k = 1:numBW                             % v0 offset
    H = 2*pi*BW(k)*Lz;
    parfor s = 1:numRF                      % v1 offset

        coeff = RFA(s) / (2*pi*para.RF.A);
    
        % compensated
        [~,~,prop_temp] = shaped_pulse_xy(spin_system,H,{Lx,Ly,Hz},...
            {Cx{1}*coeff,Cy{1}*coeff,Bg_temp},pulse_dt,rho_init,'expv-pwc');
        
         eta(k,s) = trace(unit_opt'*prop_temp)/size(unit_opt,1);

    end
end

%% plot fidelity
close all;
fig = figure('Name','Fidelity of identity propagator with v0v1 offsets');
fig.Position = [200 100 350 350];
parameters.xlabel='\nu_0, kHz';
parameters.ylabel='\nu_1, kHz';
parameters.clim = [0.95,1];
plot_pd(BW/1e3,RFA/(2*pi*1e3),real(transpose(eta)),parameters);
yticks(para.RF.A*linspace(1-para.RF.dev*2,1+para.RF.dev*2,5)/1e3);
% title(append('\eta=',num2str(mean(mean(eta)))));

%% fidelity with v0 offsets and Bg
% unit_opt = eye(size(Hz,1)); eta = zeros(numBW,numVox); 
% for k = 1:numBW                                      
%     H = 2*pi*BW(k)*Lz;
%     parfor s = 1:numVox           
% 
%         % compensated
%         [~,~,prop_temp] = shaped_pulse_xy(spin_system,H,...
%             {Lx,Ly,Hz},{Cx{1},Cy{1},Bg{s}},pulse_dt,rho_init,'expv-pwc');
% 
%          eta(k,s) = trace(unit_opt'*prop_temp)/size(unit_opt,1);
%     end
% end
% 
% %% plot fidelity
% close all;
% fig = figure('Name','Fidelity of identity propagator with v0 and Bg');
% fig.Position = [200 100 350 350];
% parameters.xlabel='\nu_0, kHz';
% parameters.ylabel='B_0 drift, uT';
% parameters.clim = [0.95,1];
% BGmat = cell2mat(transpose(Bg));
% plot_pd(BW/1e3,BGmat(:,500)*1e6,real(transpose(eta)),parameters);
% 
%  
%% fidelity with v1 offsets and Bg
% unit_opt = eye(size(Hz,1)); eta = zeros(numRF,numVox); H = sparse(zeros(size(Hz)));
% for k = 1:numRF                                      
%         coeff = RFA(k) / (2*pi*para.RF.A);
%     parfor s = 1:numVox           
% 
%         % compensated
%         [~,~,prop_temp] = shaped_pulse_xy(spin_system,H,...
%         {Lx,Ly,Hz},{Cx{1}* coeff,Cy{1}* coeff,Bg{s}},pulse_dt,rho_init,'expv-pwc');
% 
%          eta(k,s) = trace(unit_opt'*prop_temp)/size(unit_opt,1);
%     end
% end
% 
% % plot fidelity
% % close all;
% fig = figure('Name','Fidelity of identity propagator with v1 and Bg');
% fig.Position = [200 100 350 350];
% parameters.xlabel='\nu_1, kHz';
% parameters.ylabel='B_0 drift, uT';
% parameters.clim = [0.95,1];
% BGmat = cell2mat(transpose(Bg));
% plot_pd(RFA/(2*pi*1e3),BGmat(:,500)*1e6,real(transpose(eta)),parameters);
% xticks(para.RF.A*linspace(1-para.RF.dev*2,1+para.RF.dev*2,5)/1e3);


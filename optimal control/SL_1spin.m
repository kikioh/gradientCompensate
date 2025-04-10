% Sub task of spin locking optimization, and return task index if executed successfully
%
% Optimization for single spin locking
%
% Mengjia He, 2024.05.16

function result = SL_1spin(taskIndex,para,control)

% default setting for saving files
if ~isfield(para, 'save'), para.save = 'yes'; end

% report the task details 
disp(['taskIndex: ' num2str(taskIndex)]);
disp(['file name: ' para.name]);

% spin system
sys.magnet = 1;                                     	
sys.isotopes = {para.nuclei};
inter.zeeman.scalar={0};                        

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the locking state
numState = numel(para.coherence);
stateArray = cell(1,numState);
for m =1:numState

	stateArray{m}=state(spin_system,para.coherence{m},para.nuclei);
	stateArray{m}=stateArray{m}/norm(full(stateArray{m}),'fro');
end

% Get the control operators
Lx = operator(spin_system,'Lx',para.nuclei); 
Ly = operator(spin_system,'Ly',para.nuclei); 
Lz = operator(spin_system,'Lz',para.nuclei);

% Calculate the Hamiltonian
Hz = hamiltonian(assume(spin_system,'labframe','zeeman'));     
    
% drift Hamiltonian using Bg
H = cell(para.grad.numVox,para.numP);
H_stack = cell(1,para.numP);
for t = 1: para.numP   
	for k = 1:para.grad.numVox                                      
        
		% drift Hamiltonian
        H{k,t} = Hz*para.grad.Bg{k}(t);
		
		% stack drift Hamiltonian
		H_stack{t} = blkdiag(H_stack{t},H{k,t});
		
    end
	H_stack{t} = sparse(H_stack{t});
end

% stack spin vector
Lx_stack = []; Ly_stack = []; Lz_stack = [];
rho_stack = cell(1,numState);
for k = 1:para.grad.numVox                 
   
	% stack spin vector
	for m =1:numState
		rho_stack{m} = vertcat(rho_stack{m},stateArray{m});
	end

	% stack spin operators
	Lx_stack = blkdiag(Lx_stack,Lx);
	Ly_stack = blkdiag(Ly_stack,Ly);
	Lz_stack = blkdiag(Lz_stack,Lz);
end

% normilize spin state
for m =1:numState
	rho_stack{m} = rho_stack{m}/norm(full(rho_stack{m}),'fro');
end

% Waveform guess
% guess = pi*rand(1,para.numP)+pi;
guess = load(append(para.name,'.mat')).pulse_shape;
		
% control parameters
control.drifts = {{sparse(zeros(size(Hz)))}};									% drift
control.operators = {Lx,Ly};													% control operators				
control.off_ops = {Lz};															% offset opeators
control.rho_init = stateArray;													% initial state
control.rho_targ = stateArray;													% targat state
control.pulse_dt = para.tau/para.numP * ones(1,para.numP);            			% pulse duration
control.method = 'lbfgs';                                               		% optimisation method
control.parallel = 'ensemble';      											% parallelisation mode
control.penalties={'SNS'};                    									% penalty types
control.p_weights=100;                        									% penalty weights

%% first turn, sparse bandwidth and RF
control.offsets = {linspace(-para.spin.BW/2,para.spin.BW/2,31)};				% resonance offset
control.pwr_levels = 2*pi*para.RF.A*linspace(1-para.RF.dev,1+para.RF.dev,9);  	% power level

% Spinach housekeeping
% control = rmfield(control,'inherit');
spin_system = optimcon(spin_system,control);

% Run the optimization
[pulse_shape,~,para.fidelity] = fminnewton(spin_system,@grape_phase,guess);
pulse_dt = control.pulse_dt;
disp('optimization turn 1 completed');

% save experiment parameters and pulse shape
if strcmp(para.save,'yes')
	save(para.name,'para','pulse_shape','pulse_dt');
	fprintf('Saved pulse shape to %s\n', para.name);
end

%% second turn, switch on gradient, smooth shape 
control.drifts = {H_stack};                                       	% drifts
control.operators = {Lx_stack,Ly_stack};                          	% controls
control.off_ops = {Lz_stack};                                   	% offset opeators
control.rho_init = rho_stack;                                     	% starting state
control.rho_targ = rho_stack;                                     	% destination state
control.penalties={'SNS','DNS'};  									% penalty types
control.p_weights=[100,1];                        					% penalty weights

% Waveform guess
guess = load(append(para.name,'.mat')).pulse_shape;
guess = mod(guess, 2*pi);

% Spinach housekeeping
spin_system = optimcon(spin_system,control);

% Run the optimization
[pulse_shape,~,para.fidelity] = fminnewton(spin_system,@grape_phase,guess);
disp('optimization turn 2 completed');

% save experiment parameters and pulse shape
if strcmp(para.save,'yes')
	save(para.name,'para','pulse_shape','-append');
	fprintf('Saved pulse shape to %s\n', para.name);
end

%% third turn, intense bandwidth and RF, smooth shape 
control.max_iter=100;                                      									% maximum iterations
control.offsets = {linspace(-para.spin.BW/2,para.spin.BW/2,para.spin.num)};	 				% resonance offset
control.pwr_levels = 2*pi*para.RF.A*linspace(1-para.RF.dev,1+para.RF.dev,para.RF.numRF); 	% power level

% Waveform guess
guess = load(append(para.name,'.mat')).pulse_shape;

% Spinach housekeeping
spin_system = optimcon(spin_system,control);

% Run the optimization
[pulse_shape,~,para.fidelity] = fminnewton(spin_system,@grape_phase,guess);
disp('optimization turn 3 completed');

% save experiment parameters and pulse shape
if strcmp(para.save,'yes')
	save(para.name,'para','pulse_shape','-append');
	fprintf('Saved pulse shape to %s\n', para.name);
end

% return
result = taskIndex;

end



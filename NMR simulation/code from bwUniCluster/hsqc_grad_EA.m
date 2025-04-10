% Phase-sensitive HSQC pulse sequence, apply gradient field to select coherence path way.
% Syntax:
%
%              fid=hsqc_grad_EA(spin_system,parameters,H,R,K,Hz,gradient,pulse_dt)
%
% Parameters:
%
%    parameters.sweep          [F1 F2] sweep widths in each
%                              frequency direction, Hz
%
%    parameters.npoints        [F1 F2] numbers of points in
%                              each time direction
%
%    parameters.spins          {F1 F2} nuclei, e.g. {'15N','1H'}
%
%    parameters.decouple_f2    nuclei to decouple in F2, 
%                              e.g. {'15N','13C'}
%
%    parameters.decouple_f1    nuclei to decouple in F1, 
%                              e.g. {'1H','13C'}
%
%    parameters.J              primary scalar coupling, Hz
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
%
% Outputs:
%
%     fid.pos,fid.neg -  two components of the States quadrature 
%                        signal.
%
% Note: natural abundance experiments should make use of the iso-
%       tope dilution functionality. See dilute.m function.
%
% Mengjia He, 2023.11.14
%
% ref: <http://spindynamics.org/wiki/index.php?title=hsqc.m>

function fid=hsqc_grad_EA(spin_system,parameters,H,R,K,Hz,gradient,pulse_dt)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Evolution timesteps
timestep=1./parameters.sweep;

% J-coupling evolution time
delta=abs(1/(2*parameters.J));

% Initial state
rho=state(spin_system,'Lz',parameters.spins{2},'cheap');

% Detection state
coil=state(spin_system,'L+',parameters.spins{2},'cheap');

% Pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1}); Lx_F1=(Lp+Lp')/2; Ly_F1=(Lp-Lp')/2i;
Lp=operator(spin_system,'L+',parameters.spins{2}); Lx_F2=(Lp+Lp')/2; Ly_F2=(Lp-Lp')/2i;

% Pulse on F2
rho=step(spin_system,Lx_F2,rho,pi/2);

% Delta evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Inversion pulses
rho=step(spin_system,Lx_F2+Lx_F1,rho,pi);

% Delta evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Pulse on F2 
rho=step(spin_system,Ly_F2,rho,pi/2);

% Pulse on F1 with coherence selection
rho=step(spin_system,Lx_F1,rho,-pi/2) - step(spin_system,Lx_F1,rho,pi/2);

% First half of F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,parameters.npoints(1)-1,'trajectory');

% F1 dimension decoupling pulse
for n=1:numel(parameters.decouple_f1)
    Lp=operator(spin_system,'L+',parameters.decouple_f1{n});
    rho_stack=step(spin_system,(Lp+Lp')/2,rho_stack,pi);
end

% Second half of F1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,parameters.npoints(1)-1,'refocus');

% first gradient pulse
rho_stack=shaped_pulse_xy(spin_system,L,{Hz},gradient(1),pulse_dt,rho_stack);

% Inversion pulses on F1
rho_stack=step(spin_system,Lx_F1,rho_stack,pi);

% Delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,sum(pulse_dt),1,'final');

% Pulse on F1 and F2 
rho_stack=step(spin_system,Lx_F2+Lx_F1,rho_stack,pi/2);

% Delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,delta/2,1,'final');

% Inversion pulses on F1 and F2
rho_stack=step(spin_system,Lx_F2+Lx_F1,rho_stack,pi);

% second gradient pulse
rho_stack_pos=shaped_pulse_xy(spin_system,L,{Hz},gradient(2),pulse_dt,rho_stack);
rho_stack_neg=shaped_pulse_xy(spin_system,L,{Hz},{-1*gradient{2}},pulse_dt,rho_stack);

% Delta evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta/2-sum(pulse_dt),1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta/2-sum(pulse_dt),1,'final');

% F2 dimension decoupling
[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2);
[L,rho_stack_neg]=decouple(spin_system,L,rho_stack_neg,parameters.decouple_f2);

% Detection on F2
fid.pos=evolution(spin_system,L,coil,rho_stack_pos,timestep(2),parameters.npoints(2)-1,'observable');
fid.neg=evolution(spin_system,L,coil,rho_stack_neg,timestep(2),parameters.npoints(2)-1,'observable');
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
end
end

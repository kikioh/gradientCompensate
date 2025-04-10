% Main task of single spin locking optimization, setup optimizaiotn parameters.
%
% Mengjia He, 2024.05.14

% specify parallel processors and path
function SL_1spin_exam()

% start parallel pool and set path
parpool('local',40);
createPath();

% specify locking coherence
para.name = 'SL-sin-1H-1ms-RF6kHz-BW7kHz-05G';
para.coherence = {'L+','L-','Lz'};
para.nuclei = '1H';

% pulse length
para.tau = 1e-3;   				                    		% duratoin of pulse

% import gradient field
Bg =load('Bg-sin-05G-v12.mat').Bg(5);						% choose index=m for 200*m time bins
para.grad.numVox = Bg.numVox;                       		% number of voxels
para.grad.Bg = Bg.field;                            		% gradient field value, T
para.numP = Bg.numP;					            		% number of times points

% RF parameters
para.RF.A = 6e3;               
para.RF.dev = 0.2;                          
para.RF.numRF = 21;						    

% spin parameters
para.spin.BW = 7e3;                     
para.spin.num = 151; 	

% setup optimal control parameters
control.max_iter=500;                                      	% maximum iterations
control.tol_g = 1e-5;                               		% minimum gradient
control.tol_f = 0.999;                              		% minimum fidelity
control.amplitudes=ones(1,para.numP);               		% amplitude profile
control.freeze=zeros(1,para.numP);                  		% freeze mask

% execute subTasks
result = SL_1spin(1,para,control);

% close parallel pool
delete(gcp('nocreate'));

% remove path
removePath();

end
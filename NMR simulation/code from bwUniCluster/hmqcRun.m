% main function for liquid state NMR simulation
% edit the current path in line 21 and experiment function in line 28 for your need
%
% Calculation time: hours
%
% Mengjia He, 2024.04.08

% specify parallel processors and path
function hmqcRun(taskID,sample)

% start parallel computing pool
parpool('local',40);
createPath();

% parameters depending on taskID
switch taskID
	case 1, para.pulseq = 'NM'; disp('Normal pulse sequence');
	case 2, para.pulseq = 'GC'; disp('Gradient coupled pulse sequence');
	case 3, para.pulseq = 'SL'; disp('Gradient coupled pulse sequence with spin locking');
end	

% gradient pulse
para.gradFile = 'BgSeq-sin-HMQC-p50G-c04G-v400_COMSOL.mat';

% choose spin locking pulse according to locking method
para.lockMethod = 'nucleiLock';
switch para.lockMethod
	case 'nucleiLock'
		
		para.specFile = 'nucleiLock';
		para.pulseFile.H1 = 'SL-sin-1H-1ms-RF6kHz-BW5kHz-05G.mat';
		para.pulseFile.C13 = 'SL-sin-13C-1ms-RF4kHzdev17-BW6kHz-05G.mat';

	case 'coherenceLock'
		
		para.specFile = 'coherenceLock';
		para.pulseFile.LpSp = 'SL-sin-LpSp-1ms-H10BW10-C9BW15-1G.mat';
		para.pulseFile.LpSm = 'SL-sin-LpSm-1ms-H10BW10-C9BW15-1G.mat';
		para.pulseFile.Lm = 'SL-sin-LmAndLmSz-1ms-H10BW10-C9BW15-1G.mat';
		para.pulseFile.LmSp = 'SL-sin-LmSp-1ms-H10BW10-C9BW15-1G.mat';
		para.pulseFile.LmSm = 'SL-sin-LmSm-1ms-H10BW10-C9BW15-1G.mat';
		
	otherwise
		error('specify the lock as nucleiLock or coherenceLock');

end

% universal acquisition parameters
para.sample = sample;
para.spins={'13C','1H'};
para.decouple_f1={'1H'};
para.decouple_f2={'13C'};
para.axis_units='ppm';
para.npoints=[256 256];
para.zerofill=[512 512];
para.lockSwitch = 'on';
para.J=145;
		
% parameters depending on sample
switch sample

	case 'glycine'
		para.chnIndex = 1; 
		para.specFile = 'glycine-nucleiLock'; 
		para.sweep=[5000 2000];
		para.offset=[5000 2000];

	case 'glucose'
		
		para.chnIndex = 2; 
		para.specFile = 'glucose-nucleiLock';
		para.sweep=[6250 2000];
		para.offset=[10625 2000];

	otherwise 
		error('unsupported sample');
	
end

% call hmqc function
hmqc_exam(para);

% close parallel pool
delete(gcp('nocreate'));
removePath();

end

% 'C13pulseBB' means use more broadband (37.5kHz) spin locking pulse on 13C 
% 'smallRF' means use small RF amplitude spin locking pulse on 13C (9 kHz) and 1H (10 kHz)
% 'completeSS' means complete spin system
% 'smallGC' means edited gradient coupling strength
% 'lockH' means only lock the 1H spin
% 'lockC' means only lock the 13C spin
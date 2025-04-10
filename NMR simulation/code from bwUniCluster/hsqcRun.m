% main function for liquid state NMR simulation
%
% Calculation time: hours
%
% Mengjia He, 2024.11.26

% specify parallel processors and path
function hsqcRun(taskID,sample)

% start parallel computing pool
parpool('local',64);
createPath();

% parameters depending on taskID
switch taskID
	case 1, para.pulseq = 'NM'; disp('Normal pulse sequence');
	case 2, para.pulseq = 'GC'; disp('Gradient coupled pulse sequence');
	case 3, para.pulseq = 'SL'; disp('Gradient coupled pulse sequence with spin locking');
end	

% specify gradient pulse shape
para.gradFile = 'BgSeq-sin-HSQC-p50G-c04G-v400.mat';

% choose spin locking pulse according to lock method
para.pulseFile.H1 = 'SL-sin-1H-1ms-RF6kHz-BW5kHz-05G.mat';
para.pulseFile.C13 = 'SL-sin-13C-1ms-RF4kHzdev17-BW6kHz-05G.mat';

% universal acquisition parameters
para.spins={'13C','1H'};
para.decouple_f1={'1H'};
para.decouple_f2={'13C'};
para.axis_units='ppm';
para.npoints=[128 128];
para.zerofill=[512 512];
para.lockSwitch = 'on';
para.J=145;

% execute simulation sample
para.sample = sample;
switch sample
	case 'glycine'
		para.chnIndex = 1; % channel 1 for glycine and channel 2 for glucose
		para.specFile = 'glycine-quadet'; % update spectrum file name
		para.sweep=[5000 2000];
		para.offset=[5000 2000];

	case 'glucose'
		para.chnIndex = 2; % channel 1 for glycine and channel 2 for glucose
		para.specFile = 'glucose-quadet'; % update spectrum file name
		para.sweep=[6250 2000];
		para.offset=[10625 2000];

	otherwise 
		error('specify the sample as glycine or glucose');
	
end

hsqc_EA_exam(para);

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
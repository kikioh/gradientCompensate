% HMQC simulation with field gradient coupling
% Calculation time: minutes
%
% mengjia.he@kit.edu, 2024.03.28

% four gradient cases
% case NM: isolated channel, use hsqc_grad and Bg_selfSeq(1,k,:) in the parfor loop
% case GC: couple and apply gradient simultaneously, use hsqc_grad and Bg_Seq(1,k,:) in the parfor loop
% case SL: couple and apply spin locking, use hsqc_grad_SL in the parfor loop 
% case SL-off: couple and apply gradient alternately, use hsqc_grad_SL in the parfor loop and specify pulse.switch = 'off'

function hmqc_exam(para)

% specify channel 
chnIndex = para.chnIndex;

% load gradient field
data = load(para.gradFile).p1000;
BgSeq = squeeze(data.BgSeq(chnIndex,:,:));   			   % total gradient field
Bg_selfSeq = squeeze(data.Bg_selfSeq(chnIndex,:,:));       % field from self pulse
Bg_coupleSeq = squeeze(data.Bg_coupleSeq(chnIndex,:,:));   % field from couple pulse
numVox = size(BgSeq,1);    								   % number of voxels

% select spin locking pulse
SL = cell(1,3);
switch para.lockMethod
	case 'nucleiLock'
		disp('The nuclei locking pulse is applied');
		pulse_H1 = load(para.pulseFile.H1);
		pulse_C13 = load(para.pulseFile.C13);
		pulse_dt = pulse_H1.pulse_dt;
		RFdevH = pulse_H1.para.RF.dev;
		RFdevC = pulse_C13.para.RF.dev;
		
		SL{3} = {2*pi*pulse_H1.para.RF.A*cos(pulse_H1.pulse_shape),...
				 2*pi*pulse_H1.para.RF.A*sin(pulse_H1.pulse_shape),...
				 zeros(size(pulse_dt)),...
				 zeros(size(pulse_dt))};
				 
		SL{1} = {zeros(size(pulse_dt)),...
				 zeros(size(pulse_dt)),...
				 2*pi*pulse_C13.para.RF.A*cos(pulse_C13.pulse_shape),...
				 2*pi*pulse_C13.para.RF.A*sin(pulse_C13.pulse_shape)};
		SL{2} = SL{1};
			
	case 'coherenceLock'
		disp('The coherence locking pulse is applied');
	
		Lm = load(para.pulseFile.Lm);
		pulse_dt = Lm.pulse_dt;
		RFdevH = Lm.para.RF.devH;
		RFdevC = Lm.para.RF.devC;
		
		SL{3} = {2*pi*Lm.para.RF.AH*cos(Lm.pulse_shape(1,:)),...
				2*pi*Lm.para.RF.AH*sin(Lm.pulse_shape(1,:)),...
				2*pi*Lm.para.RF.AC*cos(Lm.pulse_shape(2,:)),...
				2*pi*Lm.para.RF.AC*sin(Lm.pulse_shape(2,:))};

		switch para.chnIndex 
			case 1
				pulse.coherence = 'coherenceLock-LpSp-LmSp-Lm';
						
				LpSp = load(para.pulseFile.LpSp);
				SL{1} = {2*pi*LpSp.para.RF.AH*cos(LpSp.pulse_shape(1,:)),...
						 2*pi*LpSp.para.RF.AH*sin(LpSp.pulse_shape(1,:)),...
						 2*pi*LpSp.para.RF.AC*cos(LpSp.pulse_shape(2,:)),...
						 2*pi*LpSp.para.RF.AC*sin(LpSp.pulse_shape(2,:))};
				
				LmSp = load(para.pulseFile.LmSp);						
				SL{2} = {2*pi*LmSp.para.RF.AH*cos(LmSp.pulse_shape(1,:)),...
						 2*pi*LmSp.para.RF.AH*sin(LmSp.pulse_shape(1,:)),...
						 2*pi*LmSp.para.RF.AC*cos(LmSp.pulse_shape(2,:)),...
						 2*pi*LmSp.para.RF.AC*sin(LmSp.pulse_shape(2,:))};

			case 2
				pulse.coherence = 'coherenceLock-LpSm-LmSm-Lm';
				LpSm = load(para.pulseFile.LpSm);
				SL{1} = {2*pi*LpSm.para.RF.AH*cos(LpSm.pulse_shape(1,:)),...
						 2*pi*LpSm.para.RF.AH*sin(LpSm.pulse_shape(1,:)),...
						 2*pi*LpSm.para.RF.AC*cos(LpSm.pulse_shape(2,:)),...
						 2*pi*LpSm.para.RF.AC*sin(LpSm.pulse_shape(2,:))};
				
				LmSm = load(para.pulseFile.LmSm);						
				SL{2} = {2*pi*LmSm.para.RF.AH*cos(LmSm.pulse_shape(1,:)),...
						 2*pi*LmSm.para.RF.AH*sin(LmSm.pulse_shape(1,:)),...
						 2*pi*LmSm.para.RF.AC*cos(LmSm.pulse_shape(2,:)),...
						 2*pi*LmSm.para.RF.AC*sin(LmSm.pulse_shape(2,:))};

			otherwise
				error('only support 2 channels');
		end
	
otherwise
    error('please specify locking method as nucleiLock or coherenceLock');

end
RFdevH =linspace(1-RFdevH,1+RFdevH,numVox);
RFdevC =linspace(1-RFdevC,1+RFdevC,numVox);

%% spin system
funcHandle = str2func(para.sample);
[sys,inter] = funcHandle();

% Magnet field
sys.magnet=11.74;

% Basis set
bas.formalism='sphten-liouv';
% bas.approximation='IK-2';
% bas.connectivity='scalar_couplings';
% bas.space_level=1;

bas.approximation='IK-0';
bas.level=2;

% house keeping
spin_system=create(sys,inter);

% generate isotopomers
subsystems=dilute(spin_system,'13C');

% preallocate the answer
spectrum=zeros(para.zerofill(2),para.zerofill(1),'like',1i);

% loop over isotopomers
for n=1:numel(subsystems)
    
    % build the basis
    subsystem=basis(subsystems{n},bas);

    % hamiltonian
    H = hamiltonian(assume(subsystem,'nmr'));
    H = frqoffset(subsystem,H,para);

    % Zeeman hamiltonian under 1 Telsa
    Hz = hamiltonian(assume(subsystem,'labframe','zeeman'))/sys.magnet;

    % relaxtion and kinetics
    R = sparse(zeros(size(H)));  K = sparse(zeros(size(H)));

    % preallocate the answer
    fid=zeros(para.npoints(2),para.npoints(1),'like',1i);

	switch para.pulseq
		
		case 'NM'
			parfor k=1:numVox 
				fid= fid + hmqc_grad(subsystem,para,H,R,K,Hz,Bg_selfSeq(k,:),pulse_dt);
			end	
			
		case 'GC'
			parfor k=1:numVox 
				fid= fid + hmqc_grad(subsystem,para,H,R,K,Hz,BgSeq(k,:),pulse_dt);
			end
			
		case 'SL'
			parfor k=1:numVox
				RFcoeff = {RFdevH(k),RFdevH(k),RFdevC(k),RFdevC(k)};
				fid= fid + hmqc_grad_SL(subsystem,para,H,R,K,Hz,Bg_selfSeq(k,:),Bg_coupleSeq(k,:),pulse_dt,SL,RFcoeff);
			end
		
    end
    
    % apodization
    fid=apodization(fid,'cosbell-2d');
    
    % fourier transform
    spectrum=spectrum+fftshift(fft2(fid,para.zerofill(2),para.zerofill(1)));

end

% destreaking and voxels normalized
spectrum=destreak(spectrum)/numVox;

% shift frequency for S+ coherence
if para.chnIndex == 1
    para.offset(1) = -para.offset(1);
end

%% plot 2d spectrum
fig1 = figure('Name',append(para.pulseq,'-2D spectrum'));
plot_2d(spin_system,abs(spectrum),para,20,[0.05 1.0 0.05 1.0],2,256,6,'positive');
saveas(fig1, append(para.specFile,'-HMQC-2D-',para.pulseq,'.fig'));

% plot 1H spectrum
spectrum_F2 = sum(spectrum,2);
para_temp = para;
para_temp.spins=para.spins(2);
para_temp.offset=para.offset(2);
para_temp.sweep=para.sweep(2);
para_temp.npoints=para.npoints(2);
para_temp.zerofill=para.zerofill(2);
fig2 = figure('Name',append(para.pulseq,'-1H spectrum'));
plot_1d(spin_system,abs(spectrum_F2),para_temp);
saveas(fig2, append(para.specFile,'-HMQC-1H-',para.pulseq,'.fig'));

% plot 13C spectrum
spectrum_F1 = sum(spectrum,1);
para_temp.spins=para.spins(1);
para_temp.offset=para.offset(1);
para_temp.sweep=para.sweep(1);
para_temp.npoints=para.npoints(1);
para_temp.zerofill=para.zerofill(1);
fig3 = figure('Name',append(para.pulseq,'-13C spectrum'));
plot_1d(spin_system,transpose(abs(spectrum_F1)),para_temp);
saveas(fig3, append(para.specFile,'-HMQC-13C-',para.pulseq,'.fig'));

end

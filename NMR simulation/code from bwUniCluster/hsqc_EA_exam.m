% HSQC simulation with field gradient coupling
% Calculation time: minutes
% reduced system system to acclerate calculation, data from: 
% https://bmrb.io/metabolomics/mol_summary/show_data.php?id=bmse000089
%
% mengjia.he@kit.edu, 2024.03.28

% four gradient cases
% case NM: isolated channel, use hsqc_grad and Bg_selfSeq(1,k,:) in the parfor loop
% case GC: couple and apply gradient simultaneously, use hsqc_grad and Bg_Seq(1,k,:) in the parfor loop
% case SL: couple and apply spin locking, use hsqc_grad_SL in the parfor loop 
% case SL-off: couple and apply gradient alternately, use hsqc_grad_SL in the parfor loop and specify pulse.switch = 'off'


function hsqc_EA_exam(para)

% load gradient field
data = load(para.gradFile).p1000;
BgSeq = squeeze(data.BgSeq(para.chnIndex,:,:));   			   	% total gradient field
Bg_selfSeq = squeeze(data.Bg_selfSeq(para.chnIndex,:,:));       % field from self pulse
Bg_coupleSeq = squeeze(data.Bg_coupleSeq(para.chnIndex,:,:));   % field from couple pulse
numVox = size(BgSeq,1);    										% number of voxels

% spin locking pulse for 1H
pulse_H1 = load(para.pulseFile.H1);
pulse_dt = pulse_H1.pulse_dt;
SLpulse.H1 = {2*pi*pulse_H1.para.RF.A*cos(pulse_H1.pulse_shape),2*pi*pulse_H1.para.RF.A*sin(pulse_H1.pulse_shape)};

% spin locking pulse for 13C
pulse_C13 = load(para.pulseFile.C13);
SLpulse.C13 = {2*pi*pulse_C13.para.RF.A*cos(pulse_C13.pulse_shape),2*pi*pulse_C13.para.RF.A*sin(pulse_C13.pulse_shape)};

%% spin system
funcHandle = str2func(para.sample);
[sys,inter] = funcHandle();

% Magnet field
sys.magnet=11.74;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

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

    % Zeeman hamiltonian
    Hz = hamiltonian(assume(subsystem,'labframe','zeeman'))/sys.magnet;

    % relaxtion and kinetics
    R = sparse(zeros(size(H)));  K = sparse(zeros(size(H)));

    % preallocate the answer
    temp_pos=zeros(numVox, para.npoints(2),para.npoints(1),'like',1i);
    temp_neg=zeros(numVox, para.npoints(2),para.npoints(1),'like',1i);

    switch para.pulseq

        case 'NM'

            parfor k=1:numVox
                fidvox = hsqc_grad_EA(subsystem,para,H,R,K,Hz,Bg_selfSeq(k,:),pulse_dt);
                temp_pos(k,:,:) = fidvox.pos;
                temp_neg(k,:,:) = fidvox.neg;
            end

        case 'GC'

            parfor k=1:numVox
                fidvox = hsqc_grad_EA(subsystem,para,H,R,K,Hz,BgSeq(k,:),pulse_dt);
                temp_pos(k,:,:) = fidvox.pos;
                temp_neg(k,:,:) = fidvox.neg;
            end
        case 'SL'

            parfor k=1:numVox
                fidvox = hsqc_grad_EA_SL(subsystem,para,H,R,K,Hz,Bg_selfSeq(k,:),Bg_coupleSeq(k,:),pulse_dt,SLpulse);
                temp_pos(k,:,:) = fidvox.pos;
                temp_neg(k,:,:) = fidvox.neg;
            end
    end

    % Sum up the voxels
    fid_pos = squeeze(sum(temp_pos,1));
    fid_neg = squeeze(sum(temp_neg,1));

    % Apodization
    fid_pos=apodization(fid_pos,'sqcosbell-2d');
    fid_neg=apodization(fid_neg,'sqcosbell-2d');

    % F2 Fourier transform
    f1_pos=fftshift(fft(fid_pos,para.zerofill(2),1),1);
    f1_neg=fftshift(fft(fid_neg,para.zerofill(2),1),1);

    % Form States signal
    fid=f1_pos+conj(f1_neg);

    % F1 Fourier transform
    spectrum=spectrum+fftshift(fft(fid,para.zerofill(1),2),2);

end

% destreaking and voxels normalized
spectrum=spectrum/numVox;

% shift frequency for S+ coherence
if para.chnIndex == 1
    para.offset(1) = -para.offset(1);
end

% plot 2d spectrum
fig1 = figure('Name',append(para.pulseq,'-2D spectrum'));
plot_2d(spin_system,real(spectrum),para,20,[0.05 1.0 0.05 1.0],2,256,6,'positive');
saveas(fig1, append(para.specFile,'-HSQC-2D-',para.pulseq,'.fig'));

% plot 1H spectrum
spectrum_F2 = sum(spectrum,2);
para_temp = para;
para_temp.spins=para.spins(2);
para_temp.offset=para.offset(2);
para_temp.sweep=para.sweep(2);
para_temp.npoints=para.npoints(2);
para_temp.zerofill=para.zerofill(2);
fig2 = figure('Name',append(para.pulseq,'-1H spectrum'));
plot_1d(spin_system,real(spectrum_F2),para_temp);
saveas(fig2, append(para.specFile,'-HSQC-1H-',para.pulseq,'.fig'));

% plot 13C spectrum
spectrum_F1 = sum(spectrum,1);
para_temp.spins=para.spins(1);
para_temp.offset=para.offset(1);
para_temp.sweep=para.sweep(1);
para_temp.npoints=para.npoints(1);
para_temp.zerofill=para.zerofill(1);
fig3 = figure('Name',append(para.pulseq,'-13C spectrum'));
plot_1d(spin_system,transpose(real(spectrum_F1)),para_temp);
saveas(fig3, append(para.specFile,'-HSQC-13C-',para.pulseq,'.fig'));

end

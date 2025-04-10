
close all; clearvars; clc;
%% define sample and coil

% define channel and nuclei 
num = 2;                                    % number of sample 
couple = -0.008;

% divided voxel  
nz = 400;
voxFactor = linspace(-0.5,0.5,nz);

% sin shape pulse
numP = 1000;
sinSeq = sin(linspace(pi/(2*numP),pi-pi/(2*numP),numP));

% gradient strength
maxBg = 50e-4;

% gradient sequeence for HSQC
gSeq = [1,spin('13C')/spin('1H');1,-spin('13C')/spin('1H')];

% % gradient sequeence for HMQC
% gSeq = [1,1,2*spin('13C')/spin('1H');1,1,-2*spin('13C')/spin('1H')];


% prelocate answer
numGrad = size(gSeq,2);
Bg_selfSeq = cell(num,nz,numGrad);
Bg_coupleSeq = cell(num,nz,numGrad);
BgSeq = cell(num,nz,numGrad);

for m = 1:num
    index_couple = 1 * (m==2) + 2 * (m==1);
    for n = 1:numGrad
        parfor k = 1:nz
            Bg_selfSeq{m,k,n} = maxBg  * voxFactor(k) * gSeq(m,n)* sinSeq;
            Bg_coupleSeq{index_couple,k,n} = couple * maxBg * voxFactor(k) * gSeq(m,n) * sinSeq;
        end
    end
end

for m = 1:num
    for n = 1:numGrad
        parfor k = 1:nz
            BgSeq{m,k,n} = Bg_selfSeq{m,k,n}  + Bg_coupleSeq{m,k,n};
        end
    end
end

%%
p1000.BgSeq = BgSeq;
p1000.Bg_selfSeq = Bg_selfSeq;
p1000.Bg_coupleSeq = Bg_coupleSeq;
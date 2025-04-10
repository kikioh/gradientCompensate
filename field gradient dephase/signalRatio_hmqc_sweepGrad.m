% calculate the signal senstivity ratio compare to perfect gradient sequence
%
% sweep the applied gradient strength
% Mengjia He, 2024.09.06

%% test example with HMQC sequence
% gradient pulse: max value is 75 Gauss/cm, duration is 1 ms shape is sin
% gradient ratio in detector 1 is 2:2:1, to select LpSp - LmSp - Lm
% gradient ratio in detector 2 is 2:2:-1, to select LpSm - LmSm - Lm

clearvars; close all; clc;

% calculate time integral of gradient pulse
tau = 1e-3;
grad_inter = 2*tau/pi;

% maximum gradient strength
G_max = 1.03; % Telsa/m
gamma_C = spin('13C');

% sample length
ls = 6.5e-3;

% gradient ratio
gradCoeff = linspace(0.05,1,20);
% grad1 = [2,2,1]/2;
% grad2 = [2,2,-1]/2;
grad1 = [1,-1,2]/2;
grad2 = [-1,1,-2]/2;

% quantum number of coherence
quantumNum = [5,-3,-4;3,-5,-4];

% gradient coupling ratio
coupleRatio = 0.19*1e-2;

% prelocate results
ratio = zeros(2,numel(gradCoeff));
grad = zeros(2,3);
for m = 1:numel(gradCoeff)

    grad(1,:) = gradCoeff(m) *  (grad1 - grad2*coupleRatio);
    grad(2,:) = gradCoeff(m) *  (grad2 - grad1*coupleRatio);

    for n =1:2

        % sum of phase
        phase = G_max * gamma_C * grad_inter * dot(quantumNum(n,:),grad(n,:));

        % wavelength of helix of the phase of the coherence
        lambda = 2*pi/phase;

        % signal senstivity ratio
        ratio(n,m) = lambda/(pi*ls) .* sin(pi*ls./lambda);

    end
end

fig = figure('Name','Intensity of parallel HMQC signal');
fig.Position = [200 100 450 350];
plot(gradCoeff,abs(ratio(1,:)),'LineWidth',1,'Color',[0 0.4470 0.7410 0.5]); 
hold on;
plot(gradCoeff,abs(ratio(2,:)),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
% set(gca, 'XDir','reverse');
xlabel('gradient amplitude');
ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex');
% legend('channel 1','channel 2');
set(gca,'FontName','Arial','FontSize',14);
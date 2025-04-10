%% test example with HSQC sequence
% gradient pulse: max value is 75 Gauss/cm, duration is 1 ms shape is sin
% gradient ratio in detector 1 is 2:2:-1, to select Sp - Sp - Lp - Lm
% gradient ratio in detector 2 is 2:2:1, to select Sm - Sm - Lp - Lm
% coupling ratio is 0.1% to 1%

clearvars; close all; clc;

% calculate time integral of gradient pulse
tau = 1e-3;
grad_inter = 2*tau/pi;

% maximum gradient strength
G_max = 0.75; % Telsa/m
gamma_C = spin('13C');

% sample length
ls = 8e-3;

% gradient ratio
grad1 = [2,2,1]/2;
grad2 = [2,2,-1]/2;

% quantum number of coherence
quantumNum = [-1,-1,4;1,1,4];

% gradient coupling ratio
coupleRatio = linspace(1e-3,1e-2,50);

% prelocate results
signalRatio = zeros(2,numel(coupleRatio));
grad = zeros(2,3);

for m = 1:size(signalRatio,2)

    grad(1,:) = grad1 - grad2*coupleRatio(m);
    grad(2,:) = grad2 - grad1*coupleRatio(m);

    for n =1

        % sum of phase
        phase = G_max * gamma_C * grad_inter * dot(quantumNum(n,:),grad(n,:));

        % wavelength of helix of the phase of the coherence
        lambda = abs(2*pi/phase);

        % signal senstivity ratio
        signalRatio(n,m) = lambda/(pi*ls) .* sin(pi*ls./lambda);

        % ratio(n,m) = phase;

    end
end

% plot ratio
semilogx(coupleRatio,abs(signalRatio(1,:)),'LineWidth',1.5); hold on;
semilogx(coupleRatio,ones(size(coupleRatio)),'LineWidth',1.5);
% semilogx(coupleRatio,abs(ratio(2,:)),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
% set(gca, 'XDir','reverse');
ylim([0,1.1]);
set(gca,'FontName','Arial','FontSize',14);
xlabel('spliiover ratio, %');
ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex','FontSize',17);
legend('ratio: 2:2:1','ratio: 2:2:-1');
% calculate the signal senstivity ratio compare to perfect gradient sequence
%
% Mengjia He, 2024.04.03

%% test example with HMQC sequence
% gradient pulse: max value is 75 Gauss/cm, duration is 1 ms shape is sin
% gradient ratio in detector 1 is 2:2:1, to select LpSp - LmSp - Lm
% gradient ratio in detector 2 is 2:2:-1, to select LpSm - LmSm - Lm
% coupling ratio is 0.1% to 1%

clearvars; close all; clc;

% calculate time integral of gradient pulse
tau = 1e-3;
grad_inter = 2*tau/pi;

% maximum gradient strength
G_max =1.03; % Telsa/m
gamma_C = spin('13C');

% sample length
ls = 6.5e-3;

% gradient ratio
grad1 = [2,2,1]/2;
grad2 = [2,2,-1]/2;

% quantum number of coherence
quantumNum = [5,-3,-4;3,-5,-4];

% gradient coupling ratio
coupleRatio = linspace(1e-3,1e-2,100);

% prelocate results
ratio = zeros(2,numel(coupleRatio));
grad = zeros(2,3);

for m = 1:size(ratio,2)

    grad(1,:) = grad1 - grad2*coupleRatio(m);
    grad(2,:) = grad2 - grad1*coupleRatio(m);

    for n =1:2

        % sum of phase
        phase = G_max * gamma_C * grad_inter * dot(quantumNum(n,:),grad(n,:));

        % wavelength of helix of the phase of the coherence
        lambda = 2*pi/phase;

        % signal senstivity ratio
        ratio(n,m) = lambda/(pi*ls) .* sin(pi*ls./lambda);

    end
end

% plot ratio
semilogx(coupleRatio,abs(ratio(1,:)),'LineWidth',3,'Color',[0 0.4470 0.7410 0.5]); hold on;
semilogx(coupleRatio,abs(ratio(2,:)),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
% set(gca, 'XDir','reverse');
xlabel('spillover ratio');
ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex');
legend('channel 1','channel 2');
set(gca,'FontName','Arial','FontSize',14);


%% test example with HSQC sequence
% gradient pulse: max value is 75 Gauss/cm, duration is 1 ms shape is sin
% gradient ratio in detector 1 is 2:2:-1, to select Sp - Sp - Lp - Lm
% gradient ratio in detector 2 is 2:2:1, to select Sm - Sm - Lp - Lm
% coupling ratio is 0.1% to 1%

clearvars; close all; clc;

% sweep gradient duration
tau  = [0.2,0.5,1,1.5]*1e-3;
numTau = numel(tau);

for s = 1:numTau
    % calculate time integral of gradient pulse
    tau_temp = tau(s);
    grad_inter = 2*tau_temp/pi; % for sin shape
%     grad_inter = 0.9*tau_temp; % for trapezoid shape

    % maximum gradient strength
    G_max = 0.95; % Telsa/m
    gamma_C = spin('13C');

    % sample length
    ls = 8e-3;

    % gradient ratio
    grad1 = [2,2,-1]/2;
    grad2 = [2,2,1]/2;

    % quantum number of coherence
    quantumNum = [1,1,4;-1,-1,4];

    % gradient coupling ratio
    coupleRatio = linspace(1e-3,1e-2,101);

    % prelocate results
    ratio = zeros(2,numel(coupleRatio));
    grad = zeros(2,3);

    for m = 1:size(ratio,2)

        grad(1,:) = grad1 - grad2*coupleRatio(m);
        grad(2,:) = grad2 - grad1*coupleRatio(m);

        for n =1:2

            % sum of phase
            phase = G_max * gamma_C * grad_inter * dot(quantumNum(n,:),grad(n,:));

            % wavelength of helix of the phase of the coherence
            lambda = 2*pi/phase;

            % signal senstivity ratio
            ratio(n,m) = lambda/(pi*ls) .* sin(pi*ls./lambda);

        end
    end

    % plot ratio
%     semilogx(coupleRatio,abs(ratio(1,:)),'LineWidth',3,'Color',[0 0.4470 0.7410 0.5]); hold on;
    % plot ratio
    semilogx(100*coupleRatio,abs(ratio(1,:)),'LineWidth',1); hold on;
end
% semilogx(coupleRatio,abs(ratio(2,:)),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);

% set(gca, 'XDir','reverse');
set(gca,'FontName','Arial','FontSize',14);
xlabel('spillover ratio, %');
ylabel('$\mathit{s/s_{\rm{max}}}$','Interpreter','latex','FontSize',17);
% legend('channel 1','channel 2');
legend('\tau=0.2 ms','\tau=0.5 ms','\tau=1 ms','\tau=1.5 ms');

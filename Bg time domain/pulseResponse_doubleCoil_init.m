% calculate the pulse response of double coil, when exite only one
%
% important hypothesis: homogeneous current density in the wire;
%
% coil: Helmholtz coil
%
% Mengjia He, 2023.02.28

clearvars;
close all;
clc;

% define physics constant
mu0 = 4*pi*1e-7;

% define coil
coilRadius = 0.4;
coilHeight = 0.4;
coilNum = 4;
coilCurrent = 100*[1,-1,0,0];
dx_array = 1;

% circuit parameters
R = 0.041888;
L = 0.0000018106;
M = 0.00000019772;

% Define pulse parameters
pulse_width = 1e-3;
pulse_amplitude = 1;
start_time = -1e-3;
end_time = 2e-3;

%% rectangular pulse on coil 1

% plot ideal response
t = linspace(start_time, end_time, 600);
pulse_signal = pulse_amplitude * rectpuls(t, pulse_width);
pulse_signal = circshift(pulse_signal, [0, length(t)/6]);
plot(t*1e3, pulse_signal, 'LineStyle','--','LineWidth',1.5); hold on;

% calculate the real pulse
syms x
f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3,...
    1-exp(-R/L*x), x >= 1e-3, exp(-R/L*(x-1e-3)));
x_values = linspace(start_time, end_time, 1000);
y_values = double(subs(f, x, x_values));

% plot real pulse
plot(x_values*1e3, y_values);
xlim([-1 2]);
ylim([0 1.5]);
xlabel('Time, ms');
ylabel('Amplitude, a.u.');
legend('desired pulse','real pulse');

% mutual response on coil 2
% Define the piecewise function
T0 = 1e-3;
coeff = (1+exp(R/L*T0))*T0;
f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3,...
    M*R/(L*L)*exp(-R/L*x)*x, x >= 1e-3,...
    M*R/(L*L)*exp(-R/L*x)*(coeff-exp(R/L*T0)*x));
x_values = linspace(start_time, end_time, 1000);
y_values = double(subs(f, x, x_values));

% Plot the function
figure
plot(x_values*1e3, y_values);
xlim([-1 2]);
% ylim([0 1.5]);
xlabel('Time, ms');
ylabel('Amplitude, a.u.');
% legend('desired pulse','real pulse');

%% sinc pulse on coil 1
% plot actual response
syms x
omega = pi*1e3;
% Define the piecewise function
f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3, (sin(omega*x)), x >= 1e-3, 0);
x_values = linspace(start_time, end_time, 1000);
y_values = double(subs(f, x, x_values));

% Plot the function
figure
plot(x_values*1e3, y_values,'LineStyle','--','LineWidth',1.5); hold on;
plot(x_values*1e3, y_values);
xlim([-1 2]);
ylim([0 1.5]);
xlabel('Time, ms');
ylabel('Amplitude, a.u.');
% title('sinc square response');
legend('desired pulse','real pulse');

% mutual response on coil 2
coeff = M * omega/(R^2+(omega*L)^2);
T0 = 1e-3;
% Define the piecewise function
f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3,...
    coeff*(omega*L*sin(omega*x)+R*cos(omega*x))-coeff*R*exp(-R/L*x),...
    x >= 1e-3, coeff*R*(-1-exp(-R/L*T0))*exp(-R/L*(x-T0)));
x_values = linspace(start_time, end_time, 1000);
y_values = double(subs(f, x, x_values));

% Plot the function
figure
plot(x_values*1e3, y_values);
xlim([-1 2]);
xlabel('Time, ms');
ylabel('Amplitude, a.u.');


%% sin square pulse on coil 1
% % plot actual response
% syms x
% omega = pi*1e3;
% % Define the piecewise function
% f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3, (sin(omega*x))^2, x >= 1e-3, 0);
% x_values = linspace(start_time, end_time, 1000);
% y_values = double(subs(f, x, x_values));
% 
% % Plot the function
% figure
% plot(x_values*1e3, y_values,'LineStyle','--','LineWidth',1.5); hold on;
% plot(x_values*1e3, y_values);
% xlim([-1 2]);
% ylim([0 1.5]);
% xlabel('Time, ms');
% ylabel('Amplitude, a.u.');
% % title('sinc square response');
% legend('desired pulse','real pulse');
% 
% % mutual response on coil 2
% coeff = M * omega/(R^2+(2*omega*L)^2);
% T0 = 1e-3;
% % Define the piecewise function
% f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3,...
%     coeff*(R*sin(2*omega*x)-2*omega*L*cos(2*omega*x))+coeff*2*omega*L*exp(-R/L*x),...
%     x >= 1e-3, coeff*2*omega*L*(exp(-R/L*T0)-1)*exp(-R/L*(x-T0)));
% x_values = linspace(start_time, end_time, 1000);
% y_values = double(subs(f, x, x_values));
% 
% % Plot the function
% figure
% plot(x_values*1e3, y_values);
% xlim([-1 2]);
% xlabel('Time, ms');
% ylabel('Amplitude, a.u.');
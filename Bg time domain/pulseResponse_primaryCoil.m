% calculate the gradient field of double gradient coil, using Biot-Savart
% law
%
% important hypothesis: homogeneous current density in the wire;
%
% coil: Helmholtz coil
% current: 100 A
%
% Mengjia He, 2023.02.07

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

%% rectangular pulse
% plot ideal response
% Define the parameters of the pulse
pulse_width = 1e-3;
pulse_amplitude = 1;
start_time = -1e-3;
end_time = 2e-3;

% Create a time vector
t = linspace(start_time, end_time, 600);

% Create the pulse signal
pulse_signal = pulse_amplitude * rectpuls(t, pulse_width);

% % Shift the pulse signal to set the rise at t=0
pulse_signal = circshift(pulse_signal, [0, length(t)/6]);

% Plot the pulse signal
plot(t*1e3, pulse_signal, 'LineStyle','--'); hold on;

% plot actual response
syms x

% Define the piecewise function
f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3, 1-exp(-R/L*x), x >= 1e-3, exp(-R/L*(x-1e-3)));

% Create a vector of x values
x_values = linspace(start_time, end_time, 1000);

% Evaluate the function at each x value
y_values = double(subs(f, x, x_values));

% Plot the function
plot(x_values*1e3, y_values);

xlim([-1 2]);
ylim([0 1.5]);
xlabel('Time, ms');
ylabel('Amplitude, a.u.');
title('rect pulse response');
legend('desired pulse','real pulse');

%% sinc pulse

% plot actual response
syms x

% Define the piecewise function
f(x) = piecewise(x < 0, 0, 0 <= x & x < 1e-3, (sin(pi*1e3*x))^2, x >= 1e-3, 0);

% Create a vector of x values
x_values = linspace(start_time, end_time, 1000);

% Evaluate the function at each x value
y_values = double(subs(f, x, x_values));


% Plot the function
figure
plot(x_values*1e3, y_values,'LineStyle','--','LineWidth',1.5); hold on;
plot(x_values*1e3, y_values);
xlim([-1 2]);
ylim([0 1.5]);
xlabel('Time, ms');
ylabel('Amplitude, a.u.');
title('sinc square response');
legend('desired pulse','real pulse');
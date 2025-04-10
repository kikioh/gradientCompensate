% calculate the gradient coil response
%
% mengjia.he@kit.edu, 2024.04.19
%
% reference:
% 高等数学同济大学第五版, page 681

function [Bg,current] = GIRF(phi,numVox,numSlice,Gmax)

% default number of voxels
if ~exist('numVox','var') numVox = 12; end

% default number of time slices
if ~exist('numSlice','var') numSlice = 1000; end

% default value of maximum B0 drift, T
if ~exist('Gmax','var') Gmax = 0.5e-4; end

% extract pulse parameters
T0 = 1e-3;
omega = pi/T0;
deta_t = T0/numSlice;
t_values = linspace(deta_t/2,T0-deta_t/2,numSlice);
% t_values = linspace(0,1.5*T0,numSlice);

% Q factor of the gradient coil
Q = tan(phi);

% output of gradient pulse
syms x;
if Q > 0
    f(x) = piecewise(x < 0, 0,...
        0 <= x & x < T0, Q/sqrt(Q^2+1)*exp(-omega/Q*x)+ sin(omega*x-phi),...
        x >= T0, Q/sqrt(Q^2+1)*exp(-omega/Q*x)+ sin(omega*T0-phi)*exp(-omega/Q*(x-T0)));
elseif Q == 0
    f(x) = piecewise(x < 0,0, 0 <= x & x<T0,sin(omega*x), x>=T0, 0);
end
current = double(subs(f, x, t_values));

% calculate gradient field
GScale = linspace(-Gmax/2,Gmax/2,numVox);
Bg = transpose(GScale) * current;

end

% plot two channel gradient coil and sample
%
% Mengjia He, 2023.08.08

% Define torus parameters for helmholtz gradient coil
R = 3;      			% Major radius (distance from center of tube to center of torus)
r = 0.1;      			% Minor radius (radius of the tube)
dx_array = 2.25*R;

% Generate parameter values for theta and phi
theta = linspace(0, 2*pi, 100);   % Azimuthal angle
phi = linspace(0, 2*pi, 100);     % Polar angle

% Create a meshgrid for theta and phi
[Theta, Phi] = meshgrid(theta, phi);

% Calculate torus coordinates using parametric equations
x1 = (R + r*cos(Phi)) .* cos(Theta)-dx_array/2;
x2 = (R + r*cos(Phi)) .* cos(Theta)+dx_array/2;
y1 = (R + r*cos(Phi)) .* sin(Theta)-dx_array/2;
z1 =  r * sin(Phi)+ R*0.75;
z2 = r * sin(Phi) - R*0.75;

% Plot the torus
figure;
surf(x1, y1, z1, 'EdgeColor', 'none', 'FaceColor',[.7 .7 .7],'FaceAlpha',0.5);
hold on;
surf(x1, y1, z2, 'EdgeColor', 'none', 'FaceColor',[.7 .7 .7],'FaceAlpha',0.5);
hold on;
surf(x2, y1, z1, 'EdgeColor', 'none', 'FaceColor',[.7 .7 .7],'FaceAlpha',0.5);
hold on;
surf(x2, y1, z2, 'EdgeColor', 'none', 'FaceColor',[.7 .7 .7],'FaceAlpha',0.5);
hold on;

% Plot a circle on the torus, indicate thre current line
circle_radius = 3;
circle_x1 = circle_radius * cos(theta)-dx_array/2;
circle_y1 = circle_radius * sin(theta)-dx_array/2;
circle_x2 = circle_radius * cos(theta)+dx_array/2;
circle_z1 = zeros(size(theta))+R*0.75;
circle_z2 = zeros(size(theta))-R*0.75;
plot3(circle_x1, circle_y1, circle_z1, 'Color',[0.4155 0.0037 0.1234],'LineWidth', 2);
hold on;
plot3(circle_x1, circle_y1, circle_z2, 'Color',[0.4155 0.0037 0.1234], 'LineWidth', 2);
hold on;
plot3(circle_x2, circle_y1, circle_z1, 'Color',[0.0197 0.1882 0.3804],'LineWidth', 2);
hold on;
plot3(circle_x2, circle_y1, circle_z2, 'Color',[0.0197 0.1882 0.3804], 'LineWidth', 2);
hold on;
axis off;  % Turn off axes
axis equal;   % Set aspect ratio to be equal
view(3);  % Set view to 3D

% Plot the cylinder along the z-axis, indicate the sample
cylinder_radius = 0.3;
cylinder_height = R*2;
cylinder_z = linspace(-cylinder_height/2, cylinder_height/2, 100);
cylinder_theta = linspace(0, 2*pi, 100);
[CylinderZ, CylinderTheta] = meshgrid(cylinder_z, cylinder_theta);
CylinderX1 = cylinder_radius * cos(CylinderTheta)-dx_array/2;
CylinderY1 = cylinder_radius * sin(CylinderTheta)-dx_array/2;
CylinderX2 = cylinder_radius * cos(CylinderTheta)+dx_array/2;

% sample 1
% surf(CylinderX1, CylinderY1, CylinderZ, 'EdgeColor', 'none', 'FaceColor',...
%     [0.7059 0.7804 0.9059],'FaceAlpha',0.5);
% % hold on;

% sample 2
surf(CylinderX2, CylinderY1, CylinderZ, 'EdgeColor', 'none', 'FaceColor',...
    [0.7059 0.7804 0.9059],'FaceAlpha',0.5);

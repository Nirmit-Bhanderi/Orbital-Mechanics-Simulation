% Example B: convert orbital elements to a state vector and back.
clear; clc; close all;
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));
mu = 398600.4418;
deg = pi/180;
coe0 = [10000, 0.2, 30*deg, 40*deg, 60*deg, 10*deg]; % [a,e,i,RAAN,argp,nu]
[R,V] = elements_to_state(coe0,mu);
coe1 = state_to_elements(R,V,mu);
fprintf('Position [km]:   [% .6f % .6f % .6f]\n',R);
fprintf('Velocity [km/s]: [% .6f % .6f % .6f]\n',V);
fprintf('\nRecovered elements [a e i RAAN argp nu]:\n');
fprintf('% .8g ',coe1); fprintf('\n');
fprintf('Angular elements in degrees: i=%.6f, RAAN=%.6f, argp=%.6f, nu=%.6f\n',coe1(3:6)/deg);

% Visualize the orbit represented by the recovered elements.
nuGrid = linspace(0,2*pi,500);
Rorbit = zeros(numel(nuGrid),3);
for k = 1:numel(nuGrid)
    coePlot = coe1;
    coePlot(6) = nuGrid(k);
    [Rorbit(k,:),~] = elements_to_state(coePlot,mu);
end

figure('Name','Orbit from Classical Elements','Color','w');

plot3(Rorbit(:,1), Rorbit(:,2), Rorbit(:,3), ...
    'LineWidth', 1.8);
hold on;

plot3(R(1), R(2), R(3), 'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'auto');

grid on;
grid minor;
box on;
axis equal;

xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');

title('Orbit from Classical Orbital Elements');

legend('Orbital Trajectory', 'Initial Position', ...
    'Location', 'best');

set(gca, 'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1);

view(3);

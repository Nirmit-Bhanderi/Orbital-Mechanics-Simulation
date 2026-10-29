
% Example A: Two-body orbit propagation and numerical diagnostics.
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

%% Initial conditions
mu = 398600.4418;       % Earth, km^3/s^2
R0 = [7000 0 0];        % Initial position, km
V0 = [0 7.5 1.0];       % Initial velocity, km/s
tspan = [0 2*3600];     % Propagation interval, s

%% Propagate the orbit
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, Y] = propagate_two_body(R0, V0, tspan, mu, opts);

R = Y(:, 1:3);          % Position history, km
V = Y(:, 4:6);          % Velocity history, km/s
rmag = sqrt(sum(R.^2, 2));

fprintf('Final position [km]:   [%g, %g, %g]\n', Y(end, 1:3));
fprintf('Final velocity [km/s]: [%g, %g, %g]\n', Y(end, 4:6));

%% Calculate conserved quantities
energy = 0.5*sum(V.^2, 2) - mu./rmag;  % km^2/s^2

hvec = cross(R, V, 2);
hmag = sqrt(sum(hvec.^2, 2));          % km^2/s

% Relative errors compared with initial values
energyError = abs((energy - energy(1))/energy(1));
hError = abs((hmag - hmag(1))/hmag(1));

fprintf('\nNumerical diagnostics:\n');
fprintf('Initial specific energy: %.10f km^2/s^2\n', energy(1));
fprintf('Maximum relative energy error: %.3e\n', max(energyError));
fprintf('Initial angular momentum: %.10f km^2/s\n', hmag(1));
fprintf('Maximum relative angular momentum error: %.3e\n', max(hError));

%% Figure 1: Three-dimensional orbit
fig1 = figure('Name', 'Two-Body Orbit', 'Color', 'w', ...
    'Position', [100 100 900 700]);

% Draw Earth
Re = 6378.137;  % Earth equatorial radius, km
[xe, ye, ze] = sphere(100);

surf(Re*xe, Re*ye, Re*ze, ...
    'FaceColor', [0.35 0.60 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.85);
hold on;

% Orbital trajectory
plot3(R(:,1), R(:,2), R(:,3), ...
    'Color', [0.85 0.30 0.15], 'LineWidth', 2);

% Initial and final positions
plot3(R(1,1), R(1,2), R(1,3), 'o', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', [0.20 0.65 0.30], ...
    'MarkerEdgeColor', 'k');

plot3(R(end,1), R(end,2), R(end,3), 's', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.95 0.75 0.20], ...
    'MarkerEdgeColor', 'k');

axis equal;
grid on;
box on;
view(35, 25);

xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');
title('Two-Body Orbit Propagation');

legend('Earth', 'Orbital Trajectory', ...
    'Initial Position', 'Final Position', ...
    'Location', 'northeastoutside');

set(gca, 'FontName', 'Arial', 'FontSize', 11, ...
    'LineWidth', 1, 'GridAlpha', 0.15);

camlight headlight;
lighting gouraud;
material dull;

%% Figure 2: Conservation of energy and angular momentum
fig2 = figure('Name', 'Conservation Diagnostics', 'Color', 'w', ...
    'Position', [150 100 900 650]);

tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Specific orbital energy
nexttile;
plot(t/3600, energy, 'LineWidth', 1.8, ...
    'Color', [0.10 0.40 0.75]);
grid on;
box on;
xlabel('Time (hours)');
ylabel('Energy (km^2/s^2)');
title('Specific Orbital Energy');
xlim([t(1) t(end)]/3600);

% Specific angular momentum
nexttile;
plot(t/3600, hmag, 'LineWidth', 1.8, ...
    'Color', [0.85 0.35 0.15]);
grid on;
box on;
xlabel('Time (hours)');
ylabel('Angular Momentum (km^2/s)');
title('Specific Angular Momentum Magnitude');
xlim([t(1) t(end)]/3600);

set(findall(fig2, 'Type', 'axes'), ...
    'FontName', 'Arial', 'FontSize', 11, ...
    'LineWidth', 1, 'GridAlpha', 0.15);

%% Figure 3: Relative conservation errors
fig3 = figure('Name', 'Relative Numerical Errors', 'Color', 'w', ...
    'Position', [200 150 900 550]);

semilogy(t/3600, max(energyError, eps), '-', ...
    'LineWidth', 1.9, 'Color', [0.10 0.40 0.75]);
hold on;

semilogy(t/3600, max(hError, eps), '-', ...
    'LineWidth', 1.9, 'Color', [0.85 0.35 0.15]);

grid on;
box on;

xlabel('Time (hours)');
ylabel('Relative Error');
title('Numerical Conservation Errors');

legend('Specific Energy', 'Angular Momentum', ...
    'Location', 'best');

xlim([t(1) t(end)]/3600);

set(gca, 'FontName', 'Arial', 'FontSize', 11, ...
    'LineWidth', 1, 'GridAlpha', 0.15);


% Validate Cartesian <-> classical orbital element conversions.
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

mu = 398600.4418;   % Earth, km^3/s^2
deg = pi/180;

% Each row: [a, e, i, RAAN, argPeri, trueAnomaly]
% Semimajor axis in km; angles in radians.
cases = [
    10000, 0.20, 30*deg, 40*deg, 60*deg, 10*deg;
    12000, 0.60, 10*deg, 80*deg, 25*deg, 150*deg;
    9000,  0.10, 75*deg, 20*deg, 45*deg, 270*deg;
    15000, 0.40, 60*deg, 210*deg, 130*deg, 45*deg
];

nCases = size(cases,1);
positionError = zeros(nCases,1);
velocityError = zeros(nCases,1);

fprintf('Orbital element conversion validation\n\n');

for k = 1:nCases
    coeOriginal = cases(k,:);

    % Elements -> Cartesian state
    [R0,V0] = elements_to_state(coeOriginal,mu);

    % Cartesian state -> elements
    coeRecovered = state_to_elements(R0,V0,mu);

    % Recovered elements -> Cartesian state
    [R1,V1] = elements_to_state(coeRecovered,mu);

    % Reconstruction errors
    positionError(k) = norm(R1-R0);
    velocityError(k) = norm(V1-V0);

    fprintf('Case %d: e = %.2f, inclination = %.1f deg\n', ...
        k,coeOriginal(2),coeOriginal(3)/deg);

    fprintf('  Position reconstruction error: %.3e km\n', ...
        positionError(k));

    fprintf('  Velocity reconstruction error: %.3e km/s\n\n', ...
        velocityError(k));
end

fprintf('Maximum position error: %.3e km\n',max(positionError));
fprintf('Maximum velocity error: %.3e km/s\n',max(velocityError));

%% Plot reconstruction errors
fig = figure('Name','Orbital Element Validation', ...
    'Color','w','Position',[100 100 950 650]);

% Position reconstruction error
subplot(2,1,1);

semilogy(1:nCases,max(positionError,eps),'o-', ...
    'LineWidth',1.8,'MarkerSize',7);

grid on;
grid minor;
box on;

xlabel('Test Case');
ylabel('Position Error (km)');
title('Cartesian Position Reconstruction Error');

xticks(1:nCases);
xlim([0.7 nCases+0.3]);

% Velocity reconstruction error
subplot(2,1,2);

semilogy(1:nCases,max(velocityError,eps),'s-', ...
    'LineWidth',1.8,'MarkerSize',7);

grid on;
grid minor;
box on;

xlabel('Test Case');
ylabel('Velocity Error (km/s)');
title('Cartesian Velocity Reconstruction Error');

xticks(1:nCases);
xlim([0.7 nCases+0.3]);

set(findall(fig,'Type','axes'), ...
    'FontName','Arial','FontSize',11,'LineWidth',1, ...
    'GridAlpha',0.15);

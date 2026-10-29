
% Test orbital element conversions for circular and equatorial cases.
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

mu = 398600.4418;
deg = pi/180;

% [a (km), e, i, RAAN, argPeri, nu]
cases = [
    7000, 0.0, 30*deg, 40*deg, 60*deg, 10*deg;  % Circular inclined
    8000, 0.2,  0,     0,      45*deg, 70*deg;  % Equatorial elliptical
    9000, 0.0,  0,     0,      0,      120*deg % Circular equatorial
    -12000, 1.5, 35*deg, 50*deg, 70*deg, 20*deg % Hyperbolic
    ];

labels = {
    'Circular inclined'
    'Equatorial elliptical'
    'Circular equatorial'
    'Hyperbolic'
    };

fprintf('Orbital element edge-case validation\n\n');

for k = 1:size(cases,1)
    coe0 = cases(k,:);

    [R0,V0] = elements_to_state(coe0,mu);
    coe1 = state_to_elements(R0,V0,mu);
    [R1,V1] = elements_to_state(coe1,mu);

    positionError = norm(R1-R0);
    velocityError = norm(V1-V0);

    fprintf('%s:\n',labels{k});
    fprintf('  Position error: %.3e km\n',positionError);
    fprintf('  Velocity error: %.3e km/s\n',velocityError);
    fprintf('  Recovered elements: ');
    fprintf('% .8g ',coe1);
    fprintf('\n\n');
end

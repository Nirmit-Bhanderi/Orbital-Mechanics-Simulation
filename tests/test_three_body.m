
% Automated tests for the planar three-body simulation.
clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

%% Simulation setup
G = 6.67259e-20;
m = [1e29, 1e29, 1e29];

r0 = [0 0; 300000 0; 600000 0]; % km
v0 = [0 0; 250 250; 0 0];       % km/s

y0 = [reshape(r0.', [], 1); reshape(v0.', [], 1)];

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

[t, Y] = ode45( ...
    @(t,y) three_body_rhs(t, y, m, G), ...
    [0 67000], y0, opts);

%% Calculate conserved quantities
numTimes = size(Y, 1);

totalEnergy = zeros(numTimes, 1);
totalMomentum = zeros(numTimes, 2);
totalAngularMomentum = zeros(numTimes, 1);

for k = 1:numTimes
    pos = reshape(Y(k, 1:6), 2, 3);
    vel = reshape(Y(k, 7:12), 2, 3);

    kineticEnergy = 0.5 * sum(m .* sum(vel.^2, 1));

    potentialEnergy = 0;
    for i = 1:3
        for j = i+1:3
            separation = norm(pos(:,j) - pos(:,i));
            potentialEnergy = potentialEnergy ...
                - G*m(i)*m(j)/separation;
        end
    end

    totalEnergy(k) = kineticEnergy + potentialEnergy;

    totalMomentum(k,:) = sum(vel .* m, 2).';

    totalAngularMomentum(k) = sum( ...
        m .* (pos(1,:).*vel(2,:) - pos(2,:).*vel(1,:)));
end

%% Relative conservation errors
energyError = abs( ...
    (totalEnergy - totalEnergy(1)) / totalEnergy(1));

momentumScale = norm(totalMomentum(1,:));
momentumDifference = totalMomentum - totalMomentum(1,:);
momentumError = sqrt(sum(momentumDifference.^2, 2)) ...
    / max(momentumScale, eps);

angularMomentumError = abs( ...
    (totalAngularMomentum - totalAngularMomentum(1)) ...
    / max(abs(totalAngularMomentum(1)), eps));

maxEnergyError = max(energyError);
maxMomentumError = max(momentumError);
maxAngularMomentumError = max(angularMomentumError);

%% Assertions
assert(maxEnergyError < 1e-6, ...
    'Energy conservation test failed.');

assert(maxMomentumError < 1e-10, ...
    'Linear momentum conservation test failed.');

assert(maxAngularMomentumError < 1e-7, ...
    'Angular momentum conservation test failed.');

%% Results
fprintf('Three-body conservation tests passed.\n');
fprintf('Maximum relative energy error: %.3e\n', maxEnergyError);
fprintf('Maximum relative momentum error: %.3e\n', maxMomentumError);
fprintf('Maximum relative angular momentum error: %.3e\n', ...
    maxAngularMomentumError);
fprintf('Integrated %d time points.\n', numel(t));

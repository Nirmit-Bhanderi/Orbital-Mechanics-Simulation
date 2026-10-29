
% Compare ode45 accuracy at different tolerances.
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

mu = 398600.4418;       % Earth, km^3/s^2
R0 = [7000 0 0];         % Initial position, km
V0 = [0 7.5 1.0];        % Initial velocity, km/s
tspan = [0 2*3600];       % Two hours, s

% Tolerances to compare
tolerances = [1e-3 1e-6 1e-9 1e-12];

energyError = zeros(size(tolerances));
momentumError = zeros(size(tolerances));
numPoints = zeros(size(tolerances));

for k = 1:numel(tolerances)
    tol = tolerances(k);
    opts = odeset('RelTol',tol,'AbsTol',tol*1e-3);

    [t,Y] = propagate_two_body(R0,V0,tspan,mu,opts);

    R = Y(:,1:3);
    V = Y(:,4:6);
    rmag = sqrt(sum(R.^2,2));

    energy = 0.5*sum(V.^2,2) - mu./rmag;
    h = cross(R,V,2);
    hmag = sqrt(sum(h.^2,2));

    energyError(k) = max(abs((energy-energy(1))/energy(1)));
    momentumError(k) = max(abs((hmag-hmag(1))/hmag(1)));
    numPoints(k) = numel(t);
end

% Print results
fprintf('\nTolerance study results\n');
fprintf('----------------------------------------------------------\n');
fprintf('%-12s %-18s %-18s %-12s\n', ...
    'RelTol','Energy Error','Momentum Error','Time Points');

for k = 1:numel(tolerances)
    fprintf('%-12.0e %-18.3e %-18.3e %-12d\n', ...
        tolerances(k),energyError(k),momentumError(k),numPoints(k));
end

results = table(tolerances(:), energyError(:), ...
    momentumError(:), numPoints(:), ...
    'VariableNames', {'RelTol','EnergyError', ...
    'AngularMomentumError','TimePoints'});

writetable(results, fullfile(repoRoot, ...
    'tolerance_study_results.csv'));

%% Plot tolerance study results

fig = figure('Name','Tolerance Study', ...
    'Color','w','Position',[100 100 1200 500]);

% Panel 1: Conservation errors
subplot(1,2,1);

loglog(tolerances,energyError,'o-', ...
    'LineWidth',1.8,'MarkerSize',7);
hold on;

loglog(tolerances,momentumError,'s-', ...
    'LineWidth',1.8,'MarkerSize',7);

set(gca,'XDir','reverse');
grid on;
grid minor;
box on;

xlabel('Relative Tolerance (RelTol)');
ylabel('Maximum Relative Error');
title('Conservation Accuracy');

legend('Energy Error','Angular Momentum Error', ...
    'Location','best');

set(gca,'FontName','Arial','FontSize',11,'LineWidth',1);

% Panel 2: Computational cost
subplot(1,2,2);

semilogx(tolerances,numPoints,'o-', ...
    'LineWidth',1.8,'MarkerSize',7);

set(gca,'XDir','reverse');
grid on;
grid minor;
box on;

xlabel('Relative Tolerance (RelTol)');
ylabel('Number of Solver Time Points');
title('Computational Cost');

set(gca,'FontName','Arial','FontSize',11,'LineWidth',1);

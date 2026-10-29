
% Example D: planar three-body gravitational simulation.
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

G = 6.67259e-20;                    % km^3/(kg s^2)
m = [1e29 1e29 1e29];               % kg
r0 = [0 0; 300000 0; 600000 0];     % km
v0 = [0 0; 250 250; 0 0];           % km/s

% State order: all positions, then all velocities.
y0 = [reshape(r0.',[],1); reshape(v0.',[],1)];

[t,Y] = ode45(@(t,y) three_body_rhs(t,y,m,G), ...
    [0 67000],y0,odeset('RelTol',1e-9,'AbsTol',1e-10));

% Extract positions and velocities.
X = [Y(:,1),Y(:,3),Y(:,5)];
Z = [Y(:,2),Y(:,4),Y(:,6)];

%% Figure 1: Three-body trajectories
fig1 = figure('Name','Three-Body Motion', ...
    'Color','w','Position',[100 100 1000 750]);

hold on;

plot(X(:,1),Z(:,1),'LineWidth',1.8, ...
    'DisplayName','Body 1');
plot(X(:,2),Z(:,2),'LineWidth',1.8, ...
    'DisplayName','Body 2');
plot(X(:,3),Z(:,3),'LineWidth',1.8, ...
    'DisplayName','Body 3');

% Initial positions.
plot(r0(:,1),r0(:,2),'o','MarkerSize',8, ...
    'MarkerFaceColor',[0.20 0.65 0.30], ...
    'MarkerEdgeColor','k','HandleVisibility','off');

% Final positions.
plot(X(end,:),Z(end,:),'s','MarkerSize',8, ...
    'MarkerFaceColor',[0.95 0.75 0.20], ...
    'MarkerEdgeColor','k','HandleVisibility','off');

grid on;
box on;
axis equal;

xlabel('X Position (km)');
ylabel('Y Position (km)');
title('Planar Three-Body Gravitational Motion');

legend('Location','best');

set(gca,'FontName','Arial','FontSize',11, ...
    'LineWidth',1,'GridAlpha',0.15);

%% Figure 2: Center-of-mass trajectory
comX = (X*m(:))/sum(m);
comY = (Z*m(:))/sum(m);

fig2 = figure('Name','Center of Mass', ...
    'Color','w','Position',[100 100 850 650]);

plot(comX,comY,'k--','LineWidth',1.8);
hold on;

plot(comX(1),comY(1),'o','MarkerSize',8, ...
    'MarkerFaceColor',[0.20 0.65 0.30], ...
    'MarkerEdgeColor','k');

plot(comX(end),comY(end),'s','MarkerSize',8, ...
    'MarkerFaceColor',[0.95 0.75 0.20], ...
    'MarkerEdgeColor','k');

grid on;
box on;
axis equal;

xlabel('X Position (km)');
ylabel('Y Position (km)');
title('Center-of-Mass Trajectory');

legend('Trajectory','Initial Position','Final Position', ...
    'Location','best');

set(gca,'FontName','Arial','FontSize',11, ...
    'LineWidth',1,'GridAlpha',0.15);

%% Conservation diagnostics
numTimes = size(Y,1);

totalEnergy = zeros(numTimes,1);
totalMomentum = zeros(numTimes,2);
totalAngularMomentum = zeros(numTimes,1);

for k = 1:numTimes
    pos = reshape(Y(k,1:6),2,3);
    vel = reshape(Y(k,7:12),2,3);

    % Total kinetic energy, kg*km^2/s^2.
    kineticEnergy = 0.5*sum(m .* sum(vel.^2,1));

    % Total gravitational potential energy.
    potentialEnergy = 0;

    for i = 1:3
        for j = i+1:3
            separation = norm(pos(:,j)-pos(:,i));

            potentialEnergy = potentialEnergy ...
                - G*m(i)*m(j)/separation;
        end
    end

    totalEnergy(k) = kineticEnergy + potentialEnergy;

    % Total linear momentum, kg*km/s.
    totalMomentum(k,:) = sum(vel .* m,2).';

    % Total angular momentum, kg*km^2/s.
    totalAngularMomentum(k) = sum( ...
        m .* (pos(1,:).*vel(2,:) - pos(2,:).*vel(1,:)));
end

% Relative deviations from initial values.
energyError = abs((totalEnergy-totalEnergy(1))/totalEnergy(1));

momentumScale = norm(totalMomentum(1,:));

momentumError = vecnorm( ...
    totalMomentum-totalMomentum(1,:),2,2) ...
    / max(momentumScale,eps);

angularMomentumError = abs( ...
    (totalAngularMomentum-totalAngularMomentum(1)) ...
    / max(abs(totalAngularMomentum(1)),eps));

fprintf('\nThree-body conservation diagnostics:\n');
fprintf('Maximum relative energy error: %.3e\n',max(energyError));
fprintf('Maximum relative momentum error: %.3e\n',max(momentumError));
fprintf('Maximum relative angular momentum error: %.3e\n', ...
    max(angularMomentumError));

%% Figure 3: Conservation diagnostics
fig3 = figure('Name','Three-Body Conservation Errors', ...
    'Color','w','Position',[100 100 1000 550]);

semilogy(t/3600,max(energyError,eps),'LineWidth',1.8);
hold on;

semilogy(t/3600,max(momentumError,eps),'LineWidth',1.8);
semilogy(t/3600,max(angularMomentumError,eps),'LineWidth',1.8);

grid on;
grid minor;
box on;

xlabel('Time (hours)');
ylabel('Relative Deviation');
title('Three-Body Conservation Diagnostics');

legend('Total Energy','Linear Momentum','Angular Momentum', ...
    'Location','best');

set(gca,'FontName','Arial','FontSize',11, ...
    'LineWidth',1,'GridAlpha',0.15);

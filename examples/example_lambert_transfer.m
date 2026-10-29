% Example C: solve Lambert's problem and plot a propagated transfer arc.
% This uses the worked endpoint vectors from Curtis's Example 5.2.
clear; clc; close all;
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));
mu = 398600;                       % Earth, km^3/s^2
R1 = [5000 10000 2100];             % km
R2 = [-14600 2500 7000];            % km
dt = 3600;                          % s
[V1,V2] = lambert_universal(R1,R2,dt,mu,'pro');
fprintf('Lambert departure velocity [km/s]: [% .6f % .6f % .6f]\n',V1);
fprintf('Lambert arrival velocity   [km/s]: [% .6f % .6f % .6f]\n',V2);
% Numerically propagate from the departure state to visualize the arc.
opts = odeset('RelTol',1e-10,'AbsTol',1e-11);
[t,Y] = propagate_two_body(R1,V1,[0 dt],mu,opts);



%% Figure: Lambert transfer trajectory

fig = figure('Name','Lambert Transfer','Color','w', ...
    'Position',[100 100 1000 750]);

% Draw Earth for spatial context
Re = 6378.137; % Earth equatorial radius, km
[xe,ye,ze] = sphere(100);

surf(Re*xe,Re*ye,Re*ze, ...
    'FaceColor',[0.35 0.60 0.85], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.85);
hold on;

% Transfer trajectory
plot3(Y(:,1),Y(:,2),Y(:,3), ...
    'Color',[0.85 0.30 0.15], ...
    'LineWidth',2);

% Departure position
plot3(R1(1),R1(2),R1(3),'o', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',[0.20 0.65 0.30], ...
    'MarkerEdgeColor','k');

% Arrival position
plot3(R2(1),R2(2),R2(3),'s', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',[0.95 0.75 0.20], ...
    'MarkerEdgeColor','k');

% Formatting
axis equal;
grid on;
box on;
view(35,25);

xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');

title('Lambert Transfer Trajectory');

legend('Earth','Transfer Trajectory', ...
    'Departure Position','Arrival Position', ...
    'Location','northeastoutside');

set(gca,'FontName','Arial', ...
    'FontSize',11, ...
    'LineWidth',1, ...
    'GridAlpha',0.15);

camlight headlight;
lighting gouraud;
material dull;

% Report endpoint accuracy
endpointError = norm(Y(end,1:3)-R2);
fprintf('Endpoint position error: %.6g km\n',endpointError);

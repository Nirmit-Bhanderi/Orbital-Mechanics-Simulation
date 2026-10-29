
% Study how time of flight affects a Lambert transfer.
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

mu = 398600.4418;               % Earth, km^3/s^2
R1 = [5000 10000 2100];         % Departure position, km
R2 = [-14600 2500 7000];        % Arrival position, km

% Test different times of flight, in seconds.
flightTimes = [2400 3000 3600 4200 4800];

n = numel(flightTimes);
departureSpeed = zeros(n,1);
arrivalSpeed = zeros(n,1);
endpointError = zeros(n,1);

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

%% Calculate Lambert transfers
fig1 = figure('Name','Lambert Transfer Time Study', ...
    'Color','w','Position',[100 100 1100 750]);

hold on;

for k = 1:n
    dt = flightTimes(k);

    [V1,V2] = lambert_universal(R1,R2,dt,mu,'pro');

    departureSpeed(k) = norm(V1);
    arrivalSpeed(k) = norm(V2);

    [t,Y] = propagate_two_body(R1,V1,[0 dt],mu,opts);
    endpointError(k) = norm(Y(end,1:3)-R2);

    plot3(Y(:,1),Y(:,2),Y(:,3), ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('%.0f min',dt/60));
end

%% Plot Earth and transfer endpoints
Re = 6378.137;
[xe,ye,ze] = sphere(100);

surf(Re*xe,Re*ye,Re*ze, ...
    'FaceColor',[0.35 0.60 0.85], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.75, ...
    'DisplayName','Earth');

plot3(R1(1),R1(2),R1(3),'o', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',[0.20 0.65 0.30], ...
    'MarkerEdgeColor','k', ...
    'DisplayName','Departure');

plot3(R2(1),R2(2),R2(3),'s', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',[0.95 0.75 0.20], ...
    'MarkerEdgeColor','k', ...
    'DisplayName','Arrival');

axis equal;
grid on;
box on;
view(35,25);

xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');
title('Lambert Transfer Trajectories for Different Flight Times');

legend('Location','eastoutside');

set(gca,'FontName','Arial','FontSize',11, ...
    'LineWidth',1,'GridAlpha',0.15);

camlight headlight;
lighting gouraud;
material dull;

%% Display and save numerical results
results = table(flightTimes(:)/60,departureSpeed, ...
    arrivalSpeed,endpointError, ...
    'VariableNames',{'FlightTime_min','DepartureSpeed_km_s', ...
    'ArrivalSpeed_km_s','EndpointError_km'});

fprintf('\nLambert time-of-flight study\n');
disp(results);

writetable(results,fullfile(repoRoot, ...
    'lambert_time_study_results.csv'));

%% Figure 2: Transfer performance
fig2 = figure('Name','Lambert Transfer Performance', ...
    'Color','w','Position',[100 100 1100 480]);

timeMin = flightTimes(:)/60;

% Departure and arrival speeds
subplot(1,2,1);

plot(timeMin,departureSpeed,'o-', ...
    'LineWidth',1.8,'MarkerSize',7);
hold on;

plot(timeMin,arrivalSpeed,'s-', ...
    'LineWidth',1.8,'MarkerSize',7);

grid on;
grid minor;
box on;

xlabel('Time of Flight (min)');
ylabel('Speed (km/s)');
title('Endpoint Speeds');

legend('Departure Speed','Arrival Speed', ...
    'Location','best');

set(gca,'FontName','Arial','FontSize',11,'LineWidth',1);

% Endpoint propagation error
subplot(1,2,2);

semilogy(timeMin,max(endpointError,eps),'o-', ...
    'LineWidth',1.8,'MarkerSize',7);

grid on;
grid minor;
box on;

xlabel('Time of Flight (min)');
ylabel('Endpoint Position Error (km)');
title('Lambert Solver Validation');

set(gca,'FontName','Arial','FontSize',11,'LineWidth',1);

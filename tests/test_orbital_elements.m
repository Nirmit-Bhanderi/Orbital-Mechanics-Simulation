% Automated validation tests for orbital mechanics functions.
clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

mu = 398600.4418;
deg = pi/180;

%% Test 1: Standard elliptical orbit
coe = [10000, 0.2, 30*deg, 40*deg, 60*deg, 10*deg];
[R0,V0] = elements_to_state(coe,mu);
coe2 = state_to_elements(R0,V0,mu);
[R1,V1] = elements_to_state(coe2,mu);

assert(norm(R1-R0) < 1e-6, 'Elliptical position round-trip failed.');
assert(norm(V1-V0) < 1e-9, 'Elliptical velocity round-trip failed.');

%% Test 2: Circular inclined orbit and other edge cases
cases = [
    7000, 0.0, 30*deg, 40*deg, 60*deg, 10*deg;
    8000, 0.2, 0,      0,      45*deg, 70*deg;
    9000, 0.0, 0,      0,      0,      120*deg;
   -12000,1.5, 35*deg, 50*deg, 70*deg, 20*deg
];

for k = 1:size(cases,1)
    [R0,V0] = elements_to_state(cases(k,:),mu);
    coe2 = state_to_elements(R0,V0,mu);
    [R1,V1] = elements_to_state(coe2,mu);

    assert(norm(R1-R0) < 1e-6, ...
        'Position round-trip failed for edge case %d.',k);
    assert(norm(V1-V0) < 1e-9, ...
        'Velocity round-trip failed for edge case %d.',k);
end

%% Test 3: Circular retrograde equatorial orbit
coe = [7000, 0, 180*deg, 0, 0, 60*deg];

[R0,V0] = elements_to_state(coe,mu);
coe2 = state_to_elements(R0,V0,mu);
[R1,V1] = elements_to_state(coe2,mu);

assert(norm(R1-R0) < 1e-6, ...
    'Retrograde circular position round-trip failed.');

assert(norm(V1-V0) < 1e-9, ...
    'Retrograde circular velocity round-trip failed.');

%% Test 4: Lambert solver — Curtis Example 5.2
muLambert = 398600;
R1 = [5000 10000 2100];
R2 = [-14600 2500 7000];
dt = 3600;

[V1,V2] = lambert_universal(R1,R2,dt,muLambert,'pro');

V1ref = [-5.99249, 1.92536, 3.24564];
V2ref = [-3.31246,-4.19662,-0.38529];

assert(norm(V1-V1ref) < 1e-4, ...
    'Lambert departure velocity differs from reference.');
assert(norm(V2-V2ref) < 1e-4, ...
    'Lambert arrival velocity differs from reference.');

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,Y] = propagate_two_body(R1,V1,[0 dt],muLambert,opts);
endpointError = norm(Y(end,1:3)-R2);

assert(endpointError < 1e-3, ...
    'Lambert propagation endpoint error is too large.');

fprintf('Orbital-element conversion tests passed.\n');
fprintf('Retrograde circular orbit test passed.\n');
fprintf('Lambert reference test passed. Endpoint error: %.3e km\n', ...
    endpointError);

fprintf('All orbital mechanics tests passed.\n');

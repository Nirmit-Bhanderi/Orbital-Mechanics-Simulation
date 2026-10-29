function [t, Y] = propagate_two_body(R0, V0, tspan, mu, options)
%PROPAGATE_TWO_BODY Propagate a two-body orbit using ode45.
%   [t,Y] = PROPAGATE_TWO_BODY(R0,V0,TSPAN,MU) integrates the point-mass
%   equations. R0 and V0 are 1x3 vectors in km and km/s, TSPAN is [t0 tf]
%   in seconds, and MU is in km^3/s^2. Y columns are [x y z vx vy vz].
%   Optional OPTIONS is an odeset options structure.

    arguments
        R0 (1,3) double
        V0 (1,3) double
        tspan (1,2) double
        mu (1,1) double {mustBePositive}
        options = []
    end
    if norm(R0) == 0
        error('Initial position must be nonzero.');
    end
    y0 = [R0(:); V0(:)];
    if isempty(options)
        [t, Y] = ode45(@(t,y) two_body_rhs(t,y,mu), tspan, y0);
    else
        [t, Y] = ode45(@(t,y) two_body_rhs(t,y,mu), tspan, y0, options);
    end
end

function dydt = two_body_rhs(~, y, mu)
    R = y(1:3);
    V = y(4:6);
    r = norm(R);
    if r < eps
        error('Propagation reached a near-zero radius.');
    end
    A = -mu * R / r^3;
    dydt = [V; A];
end

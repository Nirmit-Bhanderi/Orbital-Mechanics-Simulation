function [R, V] = elements_to_state(coe, mu)
%ELEMENTS_TO_STATE Convert classical orbital elements to Cartesian state.
%   [R,V] = ELEMENTS_TO_STATE(COE,MU), where COE is [a,e,i,RAAN,argPeri,nu]
%   and angles are radians. R is returned in km and V in km/s.
%   For elliptic orbits a is positive; for hyperbolic orbits a is negative.

    coe = coe(:).';
    if numel(coe) ~= 6 || mu <= 0
        error('COE must have six elements and mu must be positive.');
    end
    a = coe(1); e = coe(2); inc = coe(3);
    raan = coe(4); argp = coe(5); nu = coe(6);
    p = a*(1-e^2);
    if p <= 0 || e < 0 || (e < 1 && a <= 0) || (e > 1 && a >= 0)
        error('Inconsistent a/e values: p must be positive and conic signs consistent.');
    end
    denom = 1 + e*cos(nu);
    if denom <= 0
        error('True anomaly is outside the valid branch for this conic.');
    end
    rp = (p/denom)*[cos(nu); sin(nu); 0];
    vp = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0];
    R3_W = [cos(raan) -sin(raan) 0; sin(raan) cos(raan) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
    R3_w = [cos(argp) -sin(argp) 0; sin(argp) cos(argp) 0; 0 0 1];
    Q = R3_W*R1_i*R3_w;
    R = (Q*rp).';
    V = (Q*vp).';
end

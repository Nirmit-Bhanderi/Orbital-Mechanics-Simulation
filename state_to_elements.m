function coe = state_to_elements(R, V, mu)
%STATE_TO_ELEMENTS Convert a Cartesian state to classical orbital elements.
%   COE = STATE_TO_ELEMENTS(R,V,MU) returns [a,e,i,RAAN,argPeri,nu].
%   R is in km, V in km/s, MU in km^3/s^2; angular outputs are radians.
%   For circular/equatorial special cases, undefined angles are set to zero
%   using a documented convention rather than attempting to infer them.

    R = R(:); V = V(:);
    if numel(R) ~= 3 || numel(V) ~= 3 || norm(R) == 0 || mu <= 0
        error('R and V must be 3-vectors, R nonzero, and mu positive.');
    end
    r = norm(R); v = norm(V);
    hvec = cross(R,V); h = norm(hvec);
    if h < eps
        error('Radial trajectories do not have well-defined classical orbital elements.');
    end
    nvec = cross([0;0;1],hvec); n = norm(nvec);
    evec = ((v^2 - mu/r)*R - dot(R,V)*V)/mu;
    e = norm(evec);
    energy = v^2/2 - mu/r;
    if abs(energy) < 1e-12
        a = Inf;
    else
        a = -mu/(2*energy);
    end
    inc = acos(max(-1,min(1,hvec(3)/h)));
    tol = 1e-10;
    if n > tol
        raan = mod(atan2(nvec(2),nvec(1)),2*pi);
    else
        raan = 0;
    end

    if e > tol
        if n > tol
            % Inclined eccentric orbit
            argp = acos(max(-1,min(1,dot(nvec,evec)/(n*e))));
            if evec(3) < 0
                argp = 2*pi-argp;
            end
        else
            % Equatorial eccentric orbit: RAAN is undefined.
            % Set RAAN = 0 and preserve the periapsis direction.
            if hvec(3) >= 0
                argp = atan2(evec(2),evec(1));
            else
                argp = atan2(-evec(2),evec(1));
            end
            argp = mod(argp,2*pi);
        end
    else
        % Circular orbit: argument of periapsis is undefined.
        argp = 0;
    end

    if e > tol
        nu = acos(max(-1,min(1,dot(evec,R)/(e*r))));
        if dot(R,V) < 0, nu = 2*pi-nu; end
    elseif n > tol
        nu = acos(max(-1,min(1,dot(nvec,R)/(n*r))));
        if R(3) < 0, nu = 2*pi-nu; end
    else
        nu = mod(atan2(R(2),R(1)),2*pi);
    end
    coe = [a,e,inc,raan,argp,nu];
end

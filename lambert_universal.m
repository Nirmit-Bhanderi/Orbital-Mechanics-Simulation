function [V1, V2] = lambert_universal(R1, R2, dt, mu, direction)
%LAMBERT_UNIVERSAL Solve a single-revolution Lambert transfer (universal variables).
%   [V1,V2] = LAMBERT_UNIVERSAL(R1,R2,DT,MU,DIRECTION) finds endpoint
%   velocities for a transfer between position vectors R1 and R2 in DT seconds.
%   Distances are km, velocities km/s, MU km^3/s^2. DIRECTION is 'pro' or
%   'retro'. This compact educational solver supports the single-revolution
%   branch and should not be used as a mission-design-grade Lambert solver.

    R1 = R1(:); R2 = R2(:);
    if numel(R1) ~= 3 || numel(R2) ~= 3 || norm(R1)==0 || norm(R2)==0 || dt<=0 || mu<=0
        error('Require nonzero 3-vectors, positive time of flight, and positive mu.');
    end
    if nargin < 5, direction = 'pro'; end
    r1 = norm(R1); r2 = norm(R2);
    ctheta = max(-1,min(1,dot(R1,R2)/(r1*r2)));
    theta = acos(ctheta);
    cross12 = cross(R1,R2);
    if strcmpi(direction,'pro')
        if cross12(3) < 0, theta = 2*pi-theta; end
    elseif strcmpi(direction,'retro')
        if cross12(3) >= 0, theta = 2*pi-theta; end
    else
        error('DIRECTION must be ''pro'' or ''retro''.');
    end
    denom = 1-cos(theta);
    if abs(denom) < 1e-12
        error('Lambert transfer is singular for nearly collinear endpoint vectors.');
    end
    A = sin(theta)*sqrt(r1*r2/denom);
    if abs(A) < 1e-12
        error('Lambert geometry is degenerate for these endpoints.');
    end

    % Scan for a valid sign-changing bracket for the universal variable z.
    zgrid = linspace(-4*pi^2+1e-6, 4*pi^2-1e-6, 4001);
    Fgrid = nan(size(zgrid));
    for k = 1:numel(zgrid)
        Fgrid(k) = time_equation(zgrid(k),r1,r2,A,dt,mu);
    end
    idx = find(isfinite(Fgrid(1:end-1)) & isfinite(Fgrid(2:end)) & ...
               Fgrid(1:end-1).*Fgrid(2:end) <= 0,1,'first');
    if isempty(idx)
        error('Could not bracket a single-revolution solution. Try different geometry/time of flight or a specialized Lambert solver.');
    end
    lo = zgrid(idx); hi = zgrid(idx+1);
    flo = Fgrid(idx);
    for iter = 1:200
        mid = (lo+hi)/2;
        fm = time_equation(mid,r1,r2,A,dt,mu);
        if ~isfinite(fm)
            lo = mid; flo = time_equation(lo,r1,r2,A,dt,mu); %#ok<NASGU>
            continue
        end
        if abs(fm) < 1e-8 || abs(hi-lo) < 1e-10
            break
        end
        if flo*fm <= 0
            hi = mid;
        else
            lo = mid; flo = fm;
        end
    end
    z = mid;
    C = stumpC_local(z); S = stumpS_local(z);
    y = r1+r2+A*(z*S-1)/sqrt(C);
    if y <= 0, error('No physically valid positive-y solution found.'); end
    f = 1-y/r1;
    g = A*sqrt(y/mu);
    gdot = 1-y/r2;
    if abs(g) < 1e-12, error('Transfer has a near-zero Lagrange g coefficient.'); end
    V1 = ((R2-f*R1)/g).';
    V2 = ((gdot*R2-R1)/g).';
end

function F = time_equation(z,r1,r2,A,dt,mu)
    C = stumpC_local(z); S = stumpS_local(z);
    if C <= 0 || ~isfinite(C), F = NaN; return; end
    y = r1+r2+A*(z*S-1)/sqrt(C);
    if y < 0, F = NaN; return; end
    F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*dt;
end

function C = stumpC_local(z)
    if z > 1e-8
        C = (1-cos(sqrt(z)))/z;
    elseif z < -1e-8
        C = (cosh(sqrt(-z))-1)/(-z);
    else
        C = 1/2-z/24+z^2/720-z^3/40320;
    end
end

function S = stumpS_local(z)
    if z > 1e-8
        S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z)^3);
    elseif z < -1e-8
        S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z)^3);
    else
        S = 1/6-z/120+z^2/5040-z^3/362880;
    end
end

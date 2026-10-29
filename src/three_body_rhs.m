function dydt = three_body_rhs(~, y, masses, G)
%THREE_BODY_RHS Planar equations of motion for three mutually gravitating bodies.
%   y = [x1;y1;x2;y2;x3;y3;vx1;vy1;vx2;vy2;vx3;vy3].
%   Positions are km, velocities km/s, masses kg, G km^3/(kg*s^2).

    if numel(y) ~= 12 || numel(masses) ~= 3
        error('Expected a 12-element state and three masses.');
    end
    pos = reshape(y(1:6),2,3);
    vel = reshape(y(7:12),2,3);
    acc = zeros(2,3);
    for i = 1:3
        for j = 1:3
            if i == j, continue; end
            dr = pos(:,j)-pos(:,i);
            d = norm(dr);
            if d < 1e-9
                error('Bodies approached near-zero separation; integration stopped.');
            end
            acc(:,i) = acc(:,i) + G*masses(j)*dr/d^3;
        end
    end
    dydt = [vel(:); acc(:)];
end

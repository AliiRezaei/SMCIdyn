function [desPos,desVel,desAcc] = CircleTraj(t)

    % center of circle :
    x0 = 2;
    y0 = 1;
    z0 = 4;
    
    % radius of circle :
    r = 1;
    
    % desired trajectory :
    desPos = [r * cos(t) + x0, r * sin(t) + y0, z0]';
    desVel = [- r * sin(t), r * cos(t), 0]';
    desAcc = [- r * cos(t), - r * sin(t), 0]';

end
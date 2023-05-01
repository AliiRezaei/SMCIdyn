                       %% Rotation 4*4 Matrix

function output=R(theta,axis,varargin)
    switch axis
        case 'x'
            output=[1       0            0         0
                    0   cos(theta)   -sin(theta)   0
                    0   sin(theta)    cos(theta)   0
                    0       0             0        1];
        case 'y'
            output=[cos(theta)   0   sin(theta)   0
                    0            1       0        0
                    -sin(theta)  0   cos(theta)   0
                    0            0       0        1];
        case 'z'
            output=[cos(theta)   -sin(theta)    0    0
                    sin(theta)   cos(theta)     0    0
                    0                 0         1    0
                    0                 0         0    1];
        case 'k'
            kx=varargin{1}(1);
            ky=varargin{1}(2);
            kz=varargin{1}(3);
            
            alpha=atan2(kx/sqrt(kx^2+ky^2),ky/sqrt(kx^2+ky^2));
            beta=atan2(kz,sqrt(kx^2+ky^2));
            
            output=R(alpha,'z')*R(beta,'y')*R(theta,'z')*R(-beta,'y')*R(-alpha,'z');
            
        case 'rpy'
            
            Phi=theta(1);
            Theta=theta(2);
            Psi=theta(3);
            
            output=R(Phi,'z')*R(Theta,'y')*R(Psi,'x');
            
        case 'euler'
        
            Phi=theta(1);
            Theta=theta(2);
            Psi=theta(3);
            
            output=R(Phi,'z')*R(Theta,'y')*R(Psi,'z');

            
    end
end
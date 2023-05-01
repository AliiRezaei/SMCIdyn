function [U, s, e, de] = Controller(t, state, Dyn)

    % extract states
    q = state(1:3);
    dq = state(4:6);
    
    M = Dyn.M; % Mass Matrix
    C = Dyn.C; % Coriolis Acceleration
    G = Dyn.G; % Gravity Vector
    RobotPar = Dyn.RobotPar; % Robot Dynamical Parameters
    
    Length = RobotPar.l; % length of links
     
    DH = [q(1),Length(1),0,pi/2  % denavit-hartenberg
          q(2),0,Length(2),0
          q(3),0,Length(3),0];
    
  
    J1 = [ -sin(q(1))*(Length(3)*cos(q(2) + q(3)) + Length(2)*cos(q(2))), -cos(q(1))*(Length(3)*sin(q(2) + q(3)) + Length(2)*sin(q(2))), -Length(3)*sin(q(2) + q(3))*cos(q(1))];
    J2 = [  cos(q(1))*(Length(3)*cos(q(2) + q(3)) + Length(2)*cos(q(2))), -sin(q(1))*(Length(3)*sin(q(2) + q(3)) + Length(2)*sin(q(2))), -Length(3)*sin(q(2) + q(3))*sin(q(1))];
    J3 = [                                                 0,            Length(3)*cos(q(2) + q(3)) + Length(2)*cos(q(2)),          Length(3)*cos(q(2) + q(3))];
    J = [J1; J2; J3];

    dJ11 = - cos(q(1))*(dq(1)) *( Length(2) * cos(q(2)) + Length(3) * cos(q(2) + q(3))) + sin(q(1))*(Length(2) * sin(q(2))*dq(2) + Length(3) *sin(q(2) + q(3))*(dq(2) + dq(3)));
    dJ12 = sin(q(1))*(dq(1)) * (Length(2) * sin(q(2)) + Length(3) * sin(q(2) + q(3))) - cos(q(1))*(Length(2) * cos(q(2))*dq(2) + Length(3) * cos(q(2) + q(3))*(dq(2) + dq(3)));
    dJ13 = Length(3) * sin(q(1)) * dq(1) * sin(q(2)+q(3)) - Length(3) * cos(q(1))* cos(q(2)+q(3)) * (dq(2) + dq(3));

    dJ21 = - sin(q(1))*(dq(1)) *(Length(2) * cos(q(2)) + Length(3) * cos(q(2) + q(3))) - cos(q(1))*(Length(2) * sin(q(2))*dq(2) + Length(3) *sin(q(2) + q(3))*(dq(2) + dq(3)));
    dJ22 = - cos(q(1))*(dq(1)) * (Length(2) * sin(q(2)) + Length(3) * sin(q(2) + q(3))) - sin(q(1))*(Length(2) * cos(q(2))*dq(2) + Length(3) * cos(q(2) + q(3))*(dq(2) + dq(3)));
    dJ23 = -Length(3) * cos(q(1)) * dq(1) * sin(q(2)+q(3)) - Length(3) * sin(q(1))* cos(q(2)+q(3)) * (dq(2) + dq(3));

    dJ31 = 0;
    dJ32 = -Length(2) * sin(q(2))*dq(2) - Length(3) * sin( q(2) + q(3))*(dq(2) + dq(3));
    dJ33 = -Length(3) * sin( q(2) + q(3))*(dq(2) + dq(3));

    dJ = [dJ11,dJ12,dJ13
          dJ21,dJ22,dJ23
          dJ31,dJ32,dJ33]; % derivative of jacobian

        
    [X_des,dX_des,ddX_des] = CircleTraj(t); % desired trajectory
%     X_des = [1, 2, 3]';
%     dX_des = zeros(3,1);
%     ddX_des = zeros(3,1);
    
    X = ForwardKinematics(DH);
    X = X(1:3,end); % end-effector position
    
    dX = J * dq; % end-effector velocity
    
    e = X_des - X; % error
    de = dX_des - dX; % derivative of error
    
    % controller parameters
    lambda = 9 * eye(3);
    k = 5 * eye(3);
    eta = 10 * [1.2, 0, 0; 0, 1.75, 0; 0, 0, 1.8];
    delta = 0.07;

    s = de + lambda * e;
    RelayFunc = s / (norm(s) + delta);
    
    % disturbance upper and lower bound
    d_up = 1;
    d_low = - 1;
    d1 = (d_up - d_low) / 2;
    d2 = (d_up + d_low) / 2;

    d_c = d2 - d1 * tanh(s / delta);

    % control law
    U = M * J^-1 * (ddX_des + lambda * de + eta * RelayFunc + k * s ...
         - dJ * dq) - d_c + C * dq + G;
     
%      MaxTorque = 50;
%      
%      for m = 1:3
%          if abs(U(m)) > MaxTorque
%             U(m) = MaxTorque * sign(U(m)); 
%          end
%      end

end
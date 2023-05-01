function [H,A,Rot] = ForwardKinematics(dh)

    theta = dh(:,1);
    d = dh(:,2);
    a = dh(:,3);
    alpha = dh(:,4);
    
    n = size(dh,1);
    A = cell(1,n);

    H = eye(4);
    Rot = cell(1,n);
    for i = 1:n
       
        A{i} = R(theta(i),'z') * T(d(i),'z') * T(a(i),'x') * R(alpha(i),'x');
        H = H * A{i};
        Rot{i} = H(1:3,1:3);
        
    end
    
end
function Jc = GetJacobian(dh)

    syms q1 q2 q3 real

        [H,~,RotMat] = ForwardKinematics(dh);
        X = H(1:3,end);
        Jvc = jacobian(X,[q1,q2,q3]);
        NumLinks = size(dh,1);
        z0 = [0,0,1]';
        
        
        switch(NumLinks)
            case 1
                z1 = [0,0,0]';
                z2 = [0,0,0]';
            case 2
                z1 = RotMat{1} * [0,0,1]';
                z2 = [0,0,0]';
            case 3
                
                z1 = RotMat{1} * [0,0,1]';
                z2 = RotMat{2} * [0,0,1]';
        end
      Jwc = [z0,z1,z2];
    
    Jc = [Jvc;Jwc];

end
                     %% Transmition 4*4 Matrix
function output=T(length,axis)
    output=eye(4);
    switch axis
        case 'x'
            output=[1   0   0 length
                             0   1   0   0
                             0   0   1   0
                             0   0   0   1];
        case 'y'
            output=[1   0   0     0
                             0   1   0   length
                             0   0   1     0
                             0   0   0     1];
        case 'z'
            output=[1   0   0     0
                             0   1   0     0
                             0   0   1   length
                             0   0   0    1];
    end
end
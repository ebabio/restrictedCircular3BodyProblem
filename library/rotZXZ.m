function DCM = rotZXZ(zRot1, xRot2, zRot3)
% RotZXZ implements a 3d rotation using the XZX Euler angle rotations using
% the inputs from left to right

longitudeAscNodeRot = [ cos(zRot1)   -sin(zRot1)  0;
                        sin(zRot1)   cos(zRot1)   0;
                        0            0            1];
                    
inclinationRot = [  1   0             0;
                    0   cos(xRot2)    -sin(xRot2);
                    0   sin(xRot2)    cos(xRot2)];
                
longitudePeriapsisRot = [   cos(zRot3)  -sin(zRot3)  0;
                            sin(zRot3)  cos(zRot3)   0;
                            0           0            1];
                        
DCM = longitudeAscNodeRot * inclinationRot * longitudePeriapsisRot;

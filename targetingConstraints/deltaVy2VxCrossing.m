function [updateMap, g, gGrad] = deltaVy2VxCrossing(motion)

updateMap = [zeros(5,2)];
updateMap(4,1) = 1;

if(nargout == 2)
    motion.integrateX0([0, inf]);
    
elseif(nargout == 3)
    STMvelocity = motion.STM([0, inf], 1, [4 5]);
    gGrad = STMvelocity(2:3,4:5);
end
g = motion.xf(2:3);

end
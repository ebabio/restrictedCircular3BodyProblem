function [updateMap, g, gGrad] = lyapunovTargeting(motion)

updateMap = zeros(6,2);
updateMap(5,1) = 1;

if(nargout == 2)
    motion.integrateX0([0, inf]);
    
elseif(nargout == 3)
    STMvelocity = motion.STM([0, inf], 1, [5 7]);
    gGrad = STMvelocity([2 4], [5 7]);
end
g = motion.xf([2 4]);

end
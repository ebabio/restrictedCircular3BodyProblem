function [updateMap, g, gGrad] = haloTargeting(motion)

updateMap = zeros(6,3);
updateMap(1,1) = 1;
updateMap(5,2) = 1;

if(nargout == 2)
    motion.integrateX0([0, inf]);
    
elseif(nargout == 3)
    STMvelocity = motion.STM([0, inf], 1, [1 5 7]);
    gGrad = STMvelocity([2 4 6], [1 5 7]);
end
g = motion.xf([2 4 6]);

end
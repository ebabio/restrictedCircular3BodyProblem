function [updateMap, g, gGrad] = constantJacobiTargeting2d(motion, Cgoal)
updateMap = zeros(5,3);
updateMap(1,1) = 1;
updateMap(4,2) = 1;

% force xDot to 0 since we are picking a random point
motion.x0(3) = 0;
if(nargout == 2)
    motion.integrateX0([0, inf]);
    
elseif(nargout == 3)
    STM = motion.STM([0, inf], 1, [1 4 5]);
    fx = motion.System.force([motion.x0(1:2); 0]);
    Cgrad = [2*fx(1) -2*motion.x0(4) 0];
    gGrad = [STM(2:3,[1 4 5]); Cgrad];
end
C = motion.getJacobi(motion.x0);
g = [motion.xf(2:3); C-Cgoal];

end
classdef PeriodicOrbitFamily < OrbitFamily
    
    % Most of the functionality is implemented in Orbit Family, we just
    % modify and add functionality here
    properties
        monodromy
        periodEig
    end
    
    methods
        function obj = PeriodicOrbitFamily(initialConditionGen, orbitConstraint)
            obj@OrbitFamily(initialConditionGen, orbitConstraint);
        end
        
        function plot(obj)
            n = numel(obj.continuationParameter);
            for i=1:n
                orbitTemp = copy(obj.orbit(i));     % Integrate the other half of the orbit
                orbitTemp.clearStoppingCondition();
                orbitTemp.integrateX0([0 2*orbitTemp.tf]);
                h = orbitTemp.plot();
            end
        end
        
        function [lambda, monodromy] = getStabilityProperties(obj)
            nDim = numel(obj.orbit(1).x0)/2;
            
            % Monodromy auxiliary matrices
            if(nDim == 3)
                diagonal = (-1).^((1:6)-1);%[1 -1 -1 -1 1 1]; 
            else
                diagonal = [1 -1 -1 1];
            end
            G = diag(diagonal);
            Omega = zeros(nDim);
            Omega(2,1) = -1;
            Omega(1,2) = 1;
            H = [zeros(nDim)  -eye(nDim); eye(nDim) -2*Omega];
            
            % From half orbit
            obj.monodromy = zeros(2*nDim, 2*nDim, obj.nOrbits);
            obj.periodEig = zeros(obj.nOrbits, 2*nDim);
            for i = 1:obj.nOrbits
                STMhalf = obj.orbit(i).STM();
                obj.monodromy(:,:, i)= (G * H * STMhalf' *inv(H) * G) * STMhalf; %#ok<MINV>
                obj.periodEig(i,:) = eig(obj.monodromy(:,:, i));
            end
            lambda = obj.periodEig;
            monodromy = obj.monodromy;
        end
        
        function [parameter, x0WithTf] = parameterDependence(obj)
            nDim = numel(obj.orbit(1).x0)/2;
            
            x0WithTf = zeros(2*nDim+1, obj.nOrbits);
            for i=1:obj.nOrbits
                x0WithTf(1:2*nDim,i) = obj.orbit(i).x0;
                x0WithTf(end,i) = 2*obj.orbit(i).tf;
            end
            parameter = obj.continuationParameter;
        end
            
    end
end

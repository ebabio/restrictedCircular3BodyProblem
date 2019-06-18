classdef PeriodicOrbitFamily < OrbitFamily
    
    % Most of the functionality is implemented in Orbit Family, we just
    % modify and add functionality here
    properties
        monodromy
        periodicEig
        stabilityInd
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
            obj.periodicEig = zeros(2*nDim, obj.nOrbits);
            for i = 1:obj.nOrbits
                STMhalf = obj.orbit(i).STM();
                obj.monodromy(:,:, i)= (G * H * STMhalf' *inv(H) * G) * STMhalf; %#ok<MINV>
                eigs = eig(obj.monodromy(:,:, i));
                obj.periodicEig(:,i) = obj.orderConjugateEigenvalues(eigs);
            end
            lambda = obj.periodicEig;
            monodromy = obj.monodromy;
        end
        
        function stabilityIndices = getStabilityIndices(obj)
            obj.getStabilityProperties;
            for i=1:size(obj.periodicEig,1)/2
                index0 = 2*i-1;
                obj.stabilityInd(i,:) = sum(obj.periodicEig(index0:index0+1,:),1) / 2;
            end
            stabilityIndices = obj.stabilityInd;
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
    
    methods(Static)
        function orderedEigs = orderConjugateEigenvalues(eigs)
            %  order eigenvalues so that they are in lambda(i)*lambda(i+1)=1 pairs
            order = zeros(1,length(eigs));
            [~, indices] = sort(abs(abs(eigs)-1)); %order based on diff with unit modulus
            eigs = eigs(indices);
            
            index = 1;
            while(any(order==0))
                i = find(order==0, 1);    % assign i to the maximum remaining
                order(i)=index;
                index = index+1;
                
                for j=find(order==0)
                    areConjugate = false;
                    if(imag(eigs(i))~=0)
                        areConjugate = (eigs(i)==eigs(j)');
                    end
                    areConjugate = areConjugate || (abs(eigs(i) * eigs(j) - 1) < 5e-2); % very loose conjugate criterion
                    
                    if(areConjugate)
                        order(j)=index;
                        index = index+1;
                        break;
                    end
                end
            end
            orderedEigs = eigs(order);
        end
    end
end

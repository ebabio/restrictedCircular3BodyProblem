classdef OrbitFamily < matlab.mixin.Copyable
    
    properties
        continuationRange           % We want to find orbits in this range
        continuationNumber          % Minimum number of orbits to be found
        nOrbits                     % actual number of orbits found
        
        initialConditionGen         % a valid initial orbit guess generator
        orbitConstraint             % a valid targeting constraint for RC3BPTargeter.target
        
        initialGuessStrategy        % TBD % defines the strategy for continuation. Start from previous, Taylor, STM..
        
        orbit                       % orbits in the family
        continuationParameter       % corresponding continuation parameter for the orbits
        
        verbosity                   % options for info
    end
    
    methods
        function obj = OrbitFamily(initialConditionGen, orbitConstraint)
            obj.initialConditionGen = initialConditionGen;	% continuation parameter dependent
            obj.orbitConstraint = orbitConstraint;      % continuation parameter independnet
        end
        
        function [finalOrbit, count] = findNewOrbit(obj, orbitGuess)
            % Set up targeter
            targeter = RC3BPTargeter();
            targeter.epsilon = 1e-12;   % def 1e-12 % increasing it won't improve performance beyond integrator
            
            % Find new orbit
            [finalOrbit, count] = targeter.target(obj.orbitConstraint, orbitGuess);
        end
        
        function numberOfOrbits = findFamily(obj, orbitGuess)
            
            % initialize number of empty orbits
            auxOrbit(obj.continuationNumber) = copy(orbitGuess);
            obj.orbit =  auxOrbit(obj.continuationNumber);
            
            % determine initial orbit
            orbitGuess = obj.initialConditionGen(orbitGuess, obj.continuationRange(1));   % force
            [obj.orbit(1), iters] = obj.findNewOrbit(orbitGuess);
            if(iters == -1)
                error('first orbit is not able to converge')
            else
                count = 1;
                attemptCount = 1;
                obj.continuationParameter(count) = obj.continuationRange(1);
            end
            orbitGuess = obj.initialConditionGen(obj.orbit(1), obj.continuationRange(1)+sign(diff(obj.continuationRange))*1e-5);
            [obj.orbit(2), iters] = obj.findNewOrbit(orbitGuess);
            if(iters == -1)
                error('second orbit is not able to converge')
            end
            % iterate
            while(true)
                orbitGuess = obj.nextGuess(count, attemptCount);
                [obj.orbit(count+1), iters] = obj.findNewOrbit(orbitGuess);
                if(iters == -1)
                    attemptCount = attemptCount + 1;
                    if(attemptCount>10)
                        error('number continuation parameter decreases too high')
                    end
                else
                    count = count+1;
                    attemptCount = 1;
                    if(obj.continuationParameter(count) == obj.continuationRange(2))
                        break
                    end
                end
            end
            % it will only get here after reaching the break statement
            obj.nOrbits = count;
            numberOfOrbits = obj.nOrbits;
        end
        
        function orbitGuess = nextGuess(obj, count, attemptCount)
            obj.continuationParameter(count+1) = obj.continuationParameter(count) + diff(obj.continuationRange)/(obj.continuationNumber-1) / (2^(attemptCount - 1));
            if(diff(obj.continuationRange)>0)
                obj.continuationParameter(count+1) = min([obj.continuationParameter(count+1), obj.continuationRange(2)]);
            else
                obj.continuationParameter(count+1) = max([obj.continuationParameter(count+1), obj.continuationRange(2)]);
            end
            if(attemptCount == 1)
                orbitGuess = obj.initialConditionGen(obj.orbit, obj.continuationParameter(count+1));
            else
                orbitGuess = obj.initialConditionGen(obj.orbit(1:end-1), obj.continuationParameter(count+1));
            end            
        end
        
        function plot(obj)
            n = numel(obj.continuationParameter);
            disp('Continuation parameter:')
            disp(obj.continuationParameter);
            
            for i=1:n
                obj.orbit(i).plot();
            end
        end
        
        function decimate(obj, factor)
            obj.orbit = obj.orbit(1:factor:obj.nOrbits);
            obj.continuationParameter = obj.continuationParameter(1:factor:obj.nOrbits);
            obj.nOrbits = numel(obj.orbit);
        end
    end
end

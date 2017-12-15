classdef RC3BPTargeter < matlab.mixin.Copyable
    
    properties
        maxIter
        epsilon
        reuseGradient
        saveHistory
        
        orbit
        
        orbitHistory
        gConstraintH
        
        count
    end
    
    methods
        function obj = RC3BPTargeter()
            obj.maxIter = 60;
            obj.epsilon = 1e-12;
            obj.reuseGradient = 1;  %1 means do not reuse
            obj.saveHistory = 0;
        end
        
        function [orbit, count] = target(obj, gEqualityConstraint, initialOrbit)
            iterOrbit = copy(initialOrbit);
            
            if(obj.saveHistory)
                history(1+obj.maxIter) = copy(initialOrbit);
                obj.gConstraintH = gEqualityConstraint;
            else
                obj.orbitHistory =[];
            end
            
            % First iteration
            error = 1;
            obj.count = 0;
            % Iterate
            while(error > obj.epsilon)
                obj.count = obj.count+1;
                if(mod(obj.count-1,obj.reuseGradient) == 0) %updating gradient
                    [updateMap, f, fGrad] = gEqualityConstraint(iterOrbit);
                else    % not updating gradient
                    [updateMap, f] = gEqualityConstraint(iterOrbit);
                end
                
                if(obj.saveHistory)
                    history(obj.count) = copy(iterOrbit);
                end
                
                dx0 = - pinv(fGrad)*f; %lstsqr
                iterOrbit.addDeltaX0(updateMap*dx0);
                error = norm(f);
                
                if(obj.count > obj.maxIter)
                    obj.count = -1;
                    break
                end
            end
            
            if(obj.saveHistory)
                obj.orbitHistory =history(1:obj.count);
            end
            
            obj.orbit = copy(iterOrbit);
            orbit = copy(obj.orbit);
            count = obj.count;
        end
        
        function plotHistory(obj)
            doHold = 'off';
            if(ishold)
                doHold = 'on';
            end
            hold on
            
            for i=1:obj.count
                obj.orbitHistory(i).plot();
            end
            
            hold(gca, doHold);
        end
        
        function analyze(obj, type)            
            display('ANALYSIS OF TARGETING')
            
            [~, err] = obj.gConstraintH(obj.orbitHistory(1));
            dx0 = obj.orbitHistory(obj.count).x0 - obj.orbitHistory(1).x0;
            dt = obj.orbitHistory(obj.count).tf - obj.orbitHistory(1).tf;
            [drDim, dvDim, dtDim] = obj.orbit.System.dimensionalize(dx0(1:2), dx0(3:4), dt);
            [timeVal, units] = reasonableTime(dtDim);
            display(['number of iterations to convergence: ' num2str(obj.count)])
            display(['Initial error:' num2str(norm(err))])
            disp(err)
            
            display('Change in the initial position: ')
            display(['nondimensional: ' num2str(norm(dx0(1:2)))])
            disp(dx0(1:2))
            display(['dimensional: [km] ' num2str(norm(drDim/1e3))])
            disp(drDim/1e3)
            
            display('Change in the initial velocity: ')
            display(['nondimensional: ' num2str(norm(dx0(3:4)))])
            disp(dx0(3:4))
            display(['dimensional: [km/s] ' num2str(norm(dvDim/1e3))])
            disp(dvDim/1e3)
            
            display('Change in the final time: ')
            display(['nondimensional: ' num2str(dt)])
            display(['dimensional: ' num2str(timeVal) ' '  units])
            
            if(nargin >= 2 && type >= 2)
            % First iteration
            
            display(['iteration ' num2str(1) ':'])      
            display(['error: ' num2str(norm(err))])
            display(' ')
            
            for i = 2:obj.count
                display(['iteration ' num2str(i) ':'])                
                
                dx0 = [obj.orbitHistory(i).x0 - obj.orbitHistory(i-1).x0; obj.orbitHistory(i).tf - obj.orbitHistory(i-1).tf];
                [~, err] = obj.gConstraintH(obj.orbitHistory(i));
                [~, errPrev] = obj.gConstraintH(obj.orbitHistory(i-1));
                dErr = err - errPrev;
                
                display(['targeting change: ' num2str(norm(dErr))])
                disp(dErr)
                display(['update: ' num2str(norm(dx0))])
                disp(dx0)
                display(['new error: ' num2str(norm(err))])
                disp(err)
                
            end
            end
        end
        
    end
end

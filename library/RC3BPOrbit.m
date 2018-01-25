classdef RC3BPOrbit < matlab.mixin.Copyable
    % Important note on handle classes: handle objects behave as handles/pointers
    % To create a different object with the same values use COPY method
    
    properties
        System
        
        x0
        
        x
        t
        
        xf
        tf
        phi     % STM Matrix: not assigned from STM method to account for periodic orbits
        
        C
        dC
        
        odeSol
        odeOptions
    end
    
    methods
        function obj = RC3BPOrbit(x0, System)
            if(nargin ==2)
                obj.System = System;
                obj.x0 = x0;
                obj.x = obj.x0;
                obj.t = 0;
                
                obj.C = obj.getJacobi(x0);
                
                obj.odeOptions = odeset('Reltol', 2.22045e-14, 'AbsTol',1e-18);
            elseif(nargin ==0)
                %empty constructor
            else
                error('number of inputs not recognized');
            end
        end
        
        function [rDim, vDim, tDim] = dimensionalize(obj)
            if(size(obj.x,1) == 6)
                rDim = obj.x(1:3,:) * obj.System.l;
                vDim = obj.x(4:6,:) * obj.System.l / obj.System.t;
            else
                rDim = obj.x(1:2,:) * obj.System.l;
                vDim = obj.x(3:4,:) * obj.System.l / obj.System.t;
            end
            tDim = obj.t * obj.System.t;
        end
        
        function addDeltaX0(obj, deltaX0)
            if(size(deltaX0,1)>numel(obj.x0))   %update vector may contain final time
                obj.x0 = obj.x0 + deltaX0(1:numel(obj.x0));
                obj.tf = obj.tf + deltaX0(numel(obj.x0)+1);
            else
                obj.x0 = obj.x0 + deltaX0;
            end
        end
        function obj = setStoppingCondition(obj, stoppingCondition, direction, isterminal)
            % stoppingCondition is a boolean logical condition
            % direction takes values depending on the transition that stops
            % 1: true to false, -1: false to true, 0: any change
            if(nargin<4)
                isterminal = 1;
            end
            integrationStopEvent = @ (t,y) obj.integrationEvents(t,y, stoppingCondition, direction, isterminal);
            
            % Default integration for tolerance stringent algorithms
            obj.odeOptions = odeset('Reltol', 2.22045e-14, 'AbsTol',1e-18, 'Events',integrationStopEvent);
        end
        
        function obj = clearStoppingCondition(obj)
            obj.odeOptions = odeset('Reltol', 2.22045e-14, 'AbsTol',1e-18);
        end
        
        function obj = integrateX0(obj, tf)
            if(nargin < 2)
                tf = obj.tf;
            end
            
            obj.odeSol = ode113(@obj.odeRestricted3BP, [0 tf(end)], obj.x0, obj.odeOptions);
            
            obj.t = obj.odeSol.x;
            obj.x = obj.odeSol.y;
            
            obj.tf = obj.t(end);
            obj.xf = obj.x(:,end);
        end
        
        function phi = STM(obj, tspan, method, dimensions)
            % methods of finding the STM
            % 1: by direct integration of the state matrix
            % 2: by variation of the initial conditions
            
            if(numel(obj.x0)==4)    % check if problem is planar
                n = 4;
            elseif(numel(obj.x0)==6)
                n = 6;
            end
            
            if(nargin < 2)
                tspan = [0 obj.tf];
            end
            if(nargin < 3)
                method = 1;
            end
            if(nargin < 4)
                activeDimensions = ones(1,n);
                withTime = 0;
            else
                % transform from indexes to an indicator vector
                withTime = find(dimensions==n+1);
                activeDimensions(dimensions) = ones(1,numel(dimensions));
                activeDimensions = activeDimensions(1:n);
            end
            
            if(method == 2)
                phi0 = eye(n);
                y0 = [obj.x0; phi0(:)];
                
                obj.odeSol = ode113(@obj.odeSTM, tspan, y0, obj.odeOptions);
                
                obj.t = obj.odeSol.x;
                obj.x = obj.odeSol.y(1:n,:);
                
                obj.tf = obj.t(end);
                obj.xf = obj.x(:,end);
                
                phi = zeros(size(phi0));
                phi(:) = obj.odeSol.y(n+1:end,end);
                
                if(withTime)
                    phi(:,end+1) = odeRestricted3BP(obj, [], obj.xf);
                end
                
            elseif(method == 1)
                phi = eye(n);
                delta = 1e-8;
                
                % Reference solution
                obj.odeSol = ode113(@obj.odeRestricted3BP, tspan, obj.x0, obj.odeOptions);
                obj.t = obj.odeSol.x;
                obj.x = obj.odeSol.y(1:n,:);
                
                obj.tf = obj.t(end);
                obj.xf = obj.x(:,end);
                
                % Clear integration stopping conditions
                odeOptionsSaved = obj.odeOptions;
                obj.clearStoppingCondition();
                
                % Compute STM
                for i=find(activeDimensions==1)
                    perturbation = delta * circshift([1; zeros(n-1,1)], i-1);
                    xDelta0 = obj.x0 + perturbation;
                    [~,xDeltaf] = ode113(@obj.odeRestricted3BP, [0 obj.tf], xDelta0, obj.odeOptions);
                    phi(:,i) = (xDeltaf(end,:)' - obj.xf)/delta;
                end
                for i=find(activeDimensions==0)
                    phi(:,i) = zeros(n,1);
                end
                
                if(withTime)
                    phi(:,end+1) = odeRestricted3BP(obj, [], obj.xf);
                end
                
                % Restore ode settings
                obj.odeOptions = odeOptionsSaved;
            end
            
        end
        
        function A = stateMatrix(obj, x)
            if(numel(x)==6)
                r = x(1:3);
            elseif(numel(x)==4)
                r = [x(1:2); 0];
            end
            Urr = obj.System.forceJacobian(r);
            A = [zeros(3), eye(3); Urr, [0 2 0;-2 0 0; 0 0 0]];
            
            if(numel(x)==4)
                A(:,[3 6]) = [];
                A([3 6],:) = [];
            end
        end
        
        function xDot = odeRestricted3BP(obj, ~, x)
            % This function can be easily wrapped up for integration in
            % ODE solvers
            % Works on the ND Problem
            
            if(numel(x) == 6)
                xDot = zeros(6,1);
                xDot(1:3) = x(4:6);
                xDot(4:6) = [2*x(5); -2*x(4); 0] + obj.System.force(x(1:3));
            else
                xDot = zeros(4,1);
                xDot(1:2) = x(3:4);
                f = obj.System.force([x(1:2); 0]);
                xDot(3:4) = [2*x(4); -2*x(3)] + f(1:2);
            end
        end
        
        function C = getJacobi(obj, x)
            C = zeros(size(x,2),1);
            for i = 1:numel(C)
                if(numel(x(:,i))==4)
                    v2 = x(3:4,i)'*x(3:4,i);
                    r = [x(1:2,i); 0];
                else
                    v2 = x(4:6,i)'*x(4:6,i);
                    r = x(1:3,i);
                end
                U = obj.System.potential(r);
                C(i) = -v2 + 2*U;
            end
        end
        
        function dC = jacobiAnalysis(obj, plotErr)
            obj.C = obj.getJacobi(obj.x);
            JacobiDelta = obj.C - obj.C(1);
            obj.dC = max(abs(JacobiDelta));
            if(nargin>1 && plotErr)
                title('Jacobi Constant variation with the integration')
                xlabel('t [non-dimensional units]')
                ylabel('Jacobi Constant delta [non-dimensional units]')
                plot(obj.t, JacobiDelta)
            end
            dC = obj.dC;
        end
        
        function xInertial = inertialMotion(obj, t, x)
            % Motion in inertial frame where the bigger primary has been
            % moved to the origin
            transpose = 0;
            if((size(x,2) == 6) || (size(x,2) == 4))
                transpose =1;
                x = x';
            end
            xInertial = zeros(size(x));
            
            theta = t;
            
            for i = 1:numel(theta)
                % add velocity
                if(size(x,1) == 4)
                    x(:,i) = x(:,i) + [obj.System.mu; 0; -x(2,i); x(1,i)];
                else
                    x(:,i) = x(:,i) + [obj.System.mu; 0; 0; cross([0;0;1] + x(1:3,i))];
                end
                % rotate position/velocity
                DCM = rotZXZ(theta(i), 0, 0);
                DCMrv = [DCM zeros(3); zeros(3) DCM];
                if(size(x,1) == 4) % if problem is planar
                    DCMrv([3 6],:) = [];
                    DCMrv(:,[3 6]) = [];
                end
                xInertial(:,i) = DCMrv *x(:,i); % all transposed
                
                
            end
            
            if(transpose == 1)
                xInertial = xInertial';
            end
        end
        
        function h = plot(obj, textZVC)
            % Record hold status
            doHold = 'off';
            if(ishold)
                doHold = 'on';
            end
            hold on
            
            
            % Appropiate time units
            [timeVal, units] = reasonableTime(obj.tf*obj.System.t);
            
            nDim = numel(obj.x(:,1))/2;
            if(nDim == 2)
                
                [h(1), h(2)] = obj.System.plotSystem2d();
                h(3) = plot(obj.x(1,:), obj.x(2,:));
                
                if(nargin > 1 && textZVC)
                    h(4) = obj.System.ZVS(-1.5:0.01:1.5, -1.5:0.01:1., obj.getJacobi(obj.x0), textZVC);
                    % default interval: -1.5:0.01:1.5, -1.5:0.01:1.5
                end
                
                title(['Motion from the ICs for ' num2str(timeVal) units])
                xlabel('x [non-dimensional units]')
                ylabel('y [non-dimensional units]')
            else
                
                [h(1), h(2)] = obj.System.plotSystem3d();
                h(3) = plot3(obj.x(1,:), obj.x(2,:), obj.x(3,:));
                
                title(['Motion from the ICs for ' num2str(timeVal) units])
                xlabel('x [non-dimensional units]')
                ylabel('y [non-dimensional units]')
                zlabel('z [non-dimensional units]')
                
            end
            axis equal
            
            % Reset previous hold state
            hold(gca, doHold);
        end
    end
    
    methods (Access = private)
        
        function yDot = odeSTM(obj, t, y)
            % This function can be easily wrapped up for integration in
            % ODE solvers
            % Works on the ND Problem
            
            if(numel(y) == 42) %3d problem
                xState = y(1:6);
                phi = zeros(6);
                phi(:) = y(7:end);
            elseif (numel(y) == 20) %2d problem
                xState = y(1:4);
                phi = zeros(4);
                phi(:) = y(5:end);
            else
                error('unrecognized number of inputs')
            end
            
            xDot = odeRestricted3BP(obj, t, xState);
            phiDot = obj.stateMatrix(xState)*phi;
            yDot = [xDot; phiDot(:)];
        end
        
        function [value, isterminal, direction] = integrationEvents(~, t, y, stoppingCondition, direction, isterminal)
            value = -1+2*(~stoppingCondition(t,y));
        end
    end
end
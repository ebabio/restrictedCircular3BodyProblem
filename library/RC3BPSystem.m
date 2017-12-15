classdef RC3BPSystem
    
    properties
        % dimensional values
        l
        m
        t
        mu
    end
    
    methods
        function obj = RC3BPSystem(larger, smaller)
            global G
            obj.l = smaller.a;
            obj.m = (larger.gravP + smaller.gravP) / G;
            obj.mu = smaller.gravP / (larger.gravP + smaller.gravP);
            obj.t = sqrt(obj.l^3 / (larger.gravP + smaller.gravP));
        end
        
        function [r, v, t] = nondimensionalize(obj, r, v, t)
            r = r / obj.l;
            v = v * obj.t / obj.l;
            if(nargin >=4)
                t = t / obj.t;
            end
        end
        
        function [r, v, t] = dimensionalize(obj, r, v, t)
            r = r * obj.l;
            v = v * obj.l / obj.t;
            if(nargin >=4)
                t = t * obj.t;
            end
        end
        
        function f = force(obj, rVec)
            %if nonDimensional=1 works the nonDimensional problem, 0 the
            %dimensional problem
            % this function neglects the Coriolis force
            
            d = norm(rVec + [obj.mu; 0; 0]);
            r = norm(rVec - [1-obj.mu; 0; 0]);
            
            f = zeros(3,1);
            f(1) = -(1-obj.mu)/d^3 * (rVec(1)+obj.mu) - obj.mu/r^3 * (rVec(1)-(1-obj.mu)) + rVec(1);
            f(2) = -(1-obj.mu)/d^3 * rVec(2) - obj.mu/r^3 * rVec(2) + rVec(2);
            f(3) = -(1-obj.mu)/d^3 * rVec(3) - obj.mu/r^3 * rVec(3);
        end
        
        function U = potential(obj, rVec)
            d = norm(rVec + [obj.mu; 0; 0]);
            r = norm(rVec - [1-obj.mu; 0; 0]);
            U = (1-obj.mu)/d + obj.mu/r + 1/2*norm(rVec(1:2))^2;
        end
        
        function  Urr = forceJacobian(obj, rVec)
            d = norm(rVec + [obj.mu; 0; 0]);
            r = norm(rVec - [1-obj.mu; 0; 0]);
            
            Urr= zeros(3);
            Urr(1,1) = 1 - (1-obj.mu)/d^3 - obj.mu/r^3 + ...
                3*(1-obj.mu)*(rVec(1)+obj.mu)^2/d^5 + 3*obj.mu*(rVec(1)-1+obj.mu)^2/r^5;
            Urr(2,2) = 1 - (1-obj.mu)/d^3 - obj.mu/r^3 + ...
                3*(1-obj.mu)*rVec(2)^2/d^5 + 3*obj.mu*rVec(2)^2/r^5;
            Urr(3,3) = 1 - (1-obj.mu)/d^3 - obj.mu/r^3 + ...
                3*(1-obj.mu)*rVec(3)^2/d^5 + 3*obj.mu*rVec(3)^2/r^5;
            
            Urr(1,2) = 3*(1-obj.mu)*(rVec(1)+obj.mu)*rVec(2)/d^5 ...
                + 3*obj.mu*(rVec(1)-1+obj.mu)*rVec(2)/r^5;
            Urr(1,3) = 3*(1-obj.mu)*(rVec(1)+obj.mu)*rVec(3)/d^5 ...
                + 3*obj.mu*(rVec(1)-1+obj.mu)*rVec(3)/r^5;
            Urr(2,3) = 3*(1-obj.mu)*rVec(2)*rVec(3)/d^5 ...
                + 3*obj.mu*rVec(2)*rVec(3)/r^5;
            
            Urr(2,1) = Urr(1,2);
            Urr(3,1) = Urr(1,3);
            Urr(3,2) = Urr(2,3);
        end
        
        function [rVec, gamma] = lagrangePoint(obj,Li)
            % returns the adimensionalized Lagrange points
            % L1
            if(Li == 1)
                gamma = obj.mu;
                e=1;
                % Newton-Raphson to find the L2 (see 632 HW PB3)
                while(abs(e)>1e-12)
                    f = obj.force([1-obj.mu-gamma; 0; 0]);
                    H = obj.forceJacobian([1-obj.mu-gamma; 0; 0]);
                    e = -f(1)/H(1,1);
                    gamma = gamma - e;
                end
                rVec = [1 - obj.mu - gamma; 0; 0];
                % L2
            elseif(Li == 2)
                gamma = obj.mu;
                e=1;
                % Newton-Raphson to find the L2 (see 632 HW PB3)
                while(abs(e)>1e-12)
                    f = obj.force([1-obj.mu+gamma; 0; 0]);
                    H = obj.forceJacobian([1-obj.mu+gamma; 0; 0]);
                    e = f(1)/H(1,1);
                    gamma = gamma - e;
                end
                rVec = [1 - obj.mu + gamma; 0; 0];
                % L3
            elseif(Li == 3)
                gamma = 1-obj.mu;
                e=1;
                % Newton-Raphson to find the L2 (see 632 HW PB3)
                while(abs(e)>1e-12)
                    f = obj.force([-obj.mu-gamma; 0; 0]);
                    H = obj.forceJacobian([-obj.mu-gamma; 0; 0]);
                    e = -f(1)/H(1,1);
                    gamma = gamma - e;
                end
                rVec = [-obj.mu-gamma; 0; 0];
            elseif(Li==4)
                rVec = [-obj.mu+1/2; tand(60)/2; 0];
            elseif(Li==5)
                rVec = [-obj.mu+1/2; -tand(60)/2; 0];
            end
        end
        
        function [h1, h2] = plotSystem2d(obj)
            
            % Primaries
            x1 = [-obj.mu; 0; 0];
            x2 = [1-obj.mu; 0; 0];
            
            % Lagrange points
            for i=1:5
                L(:,i) = obj.lagrangePoint(i);
%                 C(i,:) = obj.potential(L(:,i));
            end
            
            doHold = 'off';
            if(ishold)
                doHold = 'on';
            end
            hold on
            h1 = scatter([x1(1),x2(1)],[x1(2),x2(2)], '*');
            h2 = scatter(L(1,:),L(2,:), 'x');
            hold(gca, doHold);
        end
        
        function [h1, h2] = plotSystem3d(obj)
            
            % Primaries
            x1 = [-obj.mu; 0; 0];
            x2 = [1-obj.mu; 0; 0];
            
            % Lagrange points
            for i=1:5
                L(:,i) = obj.lagrangePoint(i);
            end
            
            doHold = 'off';
            if(ishold)
                doHold = 'on';
            end
            hold on
            h1 = scatter3([x1(1),x2(1)], [x1(2),x2(2)], [x1(3),x2(3)],'*');
            h2 = scatter3(L(1,:), L(2,:), L(3,:), 'x');
            hold(gca, doHold);
        end
        
        
        function h = ZVS(obj, xLarge, y, levels, ShowText, fast)
            if(nargin <=3 || isempty(levels))
                levels =  linspace(3.3,3,5);
            end
            iterate = 1;
            if(nargin <6 || fast ==0)
                iterate =0;
            end
            
            if(iterate)
                theta = linspace(0,pi,50);
                
                % points close to the larger primary
                r0 = (1-obj.mu)/3;
                xLarge = zeros(3, length(theta));
                for i = 1:numel(theta)
                    U = @(r) (2*obj.potential(r*[cos(theta(i)); sin(theta(i)); 0]-obj.mu*[1;0;0]) - levels);
                    dUdr = @(r) (2* [cos(theta(i)); sin(theta(i)); 0]'*obj.field(r*[cos(theta(i)); sin(theta(i)); 0]-obj.mu*[1;0;0]) );
                    r = newtonRaphson(U, dUdr, r0);
                    if(~isnan(r) && r>0)
                        xLarge(:,i) = r*[cos(theta(i)); sin(theta(i)); 0]-obj.mu*[1;0;0];
                    end
                end
                
                % points close to the smaller primary
                r0 = obj.mu;
                xSmall = zeros(3, length(theta));
                for i = 1:numel(theta)
                    U = @(r) (2*obj.potential(r*[cos(theta(i)); sin(theta(i)); 0]+(1-obj.mu)*[1;0;0]) - levels);
                    dUdr = @(r) (2* [cos(theta(i)); sin(theta(i)); 0]'*obj.field(r*[cos(theta(i)); sin(theta(i)); 0]+(1-obj.mu)*[1;0;0]) );
                    r = newtonRaphson(U, dUdr, r0);
                    if(~isnan(r) && r>0)
                        xSmall(:,i) = r*[cos(theta(i)); sin(theta(i)); 0]+(1-obj.mu)*[1;0;0];
                    end
                end
                
                % points in the outer disk
                r0 = 2;
                xOut = zeros(3, length(theta));
                for i = 1:numel(theta)
                    U = @(r) (2*obj.potential(r*[cos(theta(i)); sin(theta(i)); 0]) - levels);
                    dUdr = @(r) (2* [cos(theta(i)); sin(theta(i)); 0]'*obj.field(r*[cos(theta(i)); sin(theta(i)); 0]) );
                    r = newtonRaphson(U, dUdr, r0);
                    if(~isnan(r) && r>0)
                        xOut(:,i) = r*[cos(theta(i)); sin(theta(i)); 0];
                    end
                end
                
                % points around L4
                theta = logspace(log10(0.001),log10(pi/2),20);
                theta = [theta -theta]+atan(1.7)+pi/2;
                theta = [theta theta+pi];
                r0 = 0.01;
                L4Vec = obj.findLagrangePoints(4);
                xL4 = zeros(3, length(theta));
                for i = 1:numel(theta)
                    U = @(r) (2*obj.potential(r*[cos(theta(i)); sin(theta(i)); 0] + L4Vec) - levels);
                    dUdr = @(r) (2* [cos(theta(i)); sin(theta(i)); 0]'*obj.field(r*[cos(theta(i)); sin(theta(i)); 0] + L4Vec) );
                    r = newtonRaphson(U, dUdr, r0);
                    if(~isnan(r) && r>0)
                        xL4(:,i) = r*[cos(theta(i)); sin(theta(i)); 0] + L4Vec;
                    end
                end
                
                % Eliminate weird points
                xLarge(:,xLarge(1,:)==0) = [];
                xSmall(:,xSmall(1,:)==0) = [];
                xOut(:,xOut(1,:)==0) = [];
                xL4(:,xL4(1,:)==0) = [];
                
                x = [[xLarge(1,:), NaN, xSmall(1,:), NaN, xOut(1,:), NaN, xL4(1,:)] NaN, [xLarge(1,:), NaN, xSmall(1,:), NaN, xOut(1,:), NaN, xL4(1,:)]];
                y = [[xLarge(2,:), NaN, xSmall(2,:), NaN, xOut(2,:), NaN, xL4(2,:)] NaN, -[xLarge(2,:), NaN, xSmall(2,:), NaN, xOut(2,:), NaN, xL4(2,:)]];
                scatter(x, y, '.')
                axis equal
            else
                [X,Y] = meshgrid(xLarge,y);
                Z = zeros(size(X));
                for i=1:numel(X)
                    Z(i) = obj.potential([X(i); Y(i); 0]);
                end
                Z = 2*Z;
                
                
                [~, h] = contour(X, Y, Z);
                
                h.LevelList = levels;
                if(nargin >=5 && ShowText~=0)
                    h.ShowText = 'on';
                end
                axis equal
            end
        end
        
    end
end

%% Test 4:
% Obtaining the stable and unstable manifolds for a Lyapunov L1 Orbit

%% Setup workspace
clear
% Add path
path1 = [pwd,'/library'];
path2 = [pwd,'/targetingConstraints'];
addpath(path1, path2);

Constants

format shortG
format compact

clc

%% Problem parameters
EarthMoon = RC3BPSystem(Earth, Moon);

deltaX0 = +[.01; 0; 0];

%% Miscellaneous options

stableColor = [0.3010    0.7450    0.9330];
unstableColor = [0.8500    0.3250    0.0980];

%% Find Periodic Orbit

% 1. Find linear orbit
% Information from the Lagrange Point of interest
L1Vec = EarthMoon.lagrangePoint(1);

L1C = 2*EarthMoon.potential(L1Vec);
L1Urr = EarthMoon.forceJacobian(L1Vec);

% Linearized Analysis (in the plane)
A = [zeros(2), eye(2); L1Urr(1:2,1:2) [0 2;-2 0]];
[V,D] = eig(A);

% Sort eigenvalues and rearrange eigendecomposition for Collinear points
lambda = diag(D);
P = zeros(4);
for i=1:4
    if(isreal(lambda(i)))
        if(isempty(find(P(1,:), 1)))
            P(1,i) = 1;
        else
            P(2,i) = 1;
        end
    else
        if(isempty(find(P(3,:), 1)))
            P(3,i) = 1;
        else
            P(4,i) = 1;
        end
    end
end

D = P'*D*P;
V = V*P;

% Cancelling real modes
x2m = inv(V);   %transformation from x to modes
deltaX0(3:4) = -real(x2m(1:2,3:4) \ x2m(1:2,1:2) * deltaX0(1:2));

L13d = [L1Vec; zeros(3,1)];
delta3dX0 = [deltaX0(1:2); 0; deltaX0(3:4); 0.0];
orbit = RC3BPOrbit(L13d + delta3dX0, EarthMoon);
orbit.setStoppingCondition(@(t,x) (x(2)<0), 0);
orbit.integrateX0([0, inf]);

% 2. Find periodic orbit
% Targeter settings
targeterPeriodic = RC3BPTargeter();
targeterPeriodic.saveHistory = 1;
targeterPeriodic.epsilon = 1e-12;   % def 1e-12 % increasing it won't improve performance beyond integrator
targetingCondition = @lyapunovTargeting;
[halfLyapunov, count] = targeterPeriodic.target(targetingCondition, orbit);

% Display results
display('Necessary initial conditions for Periodic Orbit:')
disp([halfLyapunov.x0; halfLyapunov.tf])

%3. Propagate Lyapunov orbit
lyapunovOrbit = copy(halfLyapunov);
lyapunovOrbit.setStoppingCondition(@(t,x) (x(2)<0), -1);
lyapunovOrbit.integrateX0([0, inf]);
dcPeriodic = lyapunovOrbit.jacobiAnalysis();

periodicError = lyapunovOrbit.x0 -lyapunovOrbit.xf;
[periodicErrorR, periodicErrorV]  = EarthMoon.dimensionalize(periodicError(1:2), periodicError(3:4));

% Display results
figure(1)
clf reset
lyapunovOrbit.plot();
% TODO: (DEBUG) remove this for final release
daspect([1 1 1]);
axis([.7 1 -.1 .1])
%% Compute monodromy matrix at initial fixed point
% Monodromy auxiliary matrices
nDim = 3;
diagonal = (-1).^((1:6)-1);
G = diag(diagonal);
Omega = zeros(nDim);
Omega(2,1) = -1;
Omega(1,2) = 1;
H = [zeros(nDim)  -eye(nDim); eye(nDim) -2*Omega];

% From half orbit
STMhalf = halfLyapunov.STM();
monodromy = (G * H * STMhalf' *inv(H) * G) * STMhalf; %#ok<MINV>
lyapunovOrbit.phi = monodromy;

% From whole orbit
monodromy2 = lyapunovOrbit.STM();

% Mondromy matrix analysis
monodromyErr = norm(monodromy-monodromy2);
[Vmonodromy, Dmonodromy] = eig(monodromy);
lambda =diag(Dmonodromy);

% Display results
display('monodromy matrix:')
disp(monodromy)
display('monodromy error wrt to full integration')
disp(monodromyErr)
display('eigenvalues of the monodromy matrix')
disp(['   ' num2str(lambda(1))])
disp(['   ' num2str(lambda(2))])
disp(lambda(3:6))


%% Find the STM for another 20 points along the orbit

% Determine the 20 points along the orbit
nFixedPoints = 20;
deltaT = linspace(-.5, .5, nFixedPoints+1) * lyapunovOrbit.tf;
% deltaT = linspace(0, 1, nFixedPoints+1) * lyapunovOrbit.tf;

fixedPointOrbit(nFixedPoints) = RC3BPOrbit;
for i=1:nFixedPoints
    fixedPointOrbit(i) = copy(lyapunovOrbit);
    if(deltaT(i) ~=0)
        fixedPointOrbit(i).phi = fixedPointOrbit(i).STM([0 deltaT(i)]);
    else
        fixedPointOrbit(i).phi = eye(numel(fixedPointOrbit(i).x0));
    end
end

% Transform STMs
fixedPointEigv = zeros([size(monodromy), nFixedPoints]);
for i=1:nFixedPoints
    if(deltaT(i) <0)
        fixedPointEigv(:,:,i) = fixedPointOrbit(i).phi * Vmonodromy;
    else
        fixedPointEigv(:,:,i) = fixedPointOrbit(i).phi * Vmonodromy;
    end
end

% Prepare data for plotting
eigvElems = 1:3;
pFixedPointEigv = zeros(nFixedPoints, numel(eigvElems), numel(fixedPointOrbit(i).x0));
xf = zeros(numel(eigvElems), nFixedPoints);
for i=1:nFixedPoints
    pFixedPointEigv(i,:,:) = normc( real( fixedPointEigv(eigvElems, :, i) ) );
    xf(:,i) = fixedPointOrbit(i).xf(1:numel(eigvElems));
end

% Display results
quiverRatio = .1;

figure(2)
clf reset
hold on
lyapunovOrbit.plot();
for j=1:2
    if(j==1)
        color = unstableColor;
    elseif(j==2)
        color = stableColor;
    else
        color = stableColor;
    end
    h = quiver3(xf(1,:), xf(2,:), xf(3,:), ...
        pFixedPointEigv(:,1,j)', pFixedPointEigv(:,2,j)', pFixedPointEigv(:,3,j)',...
        quiverRatio);
    h.Color = color;
    h = quiver3(xf(1,:), xf(2,:), xf(3,:), ...
        -pFixedPointEigv(:,1,j)', -pFixedPointEigv(:,2,j)', -pFixedPointEigv(:,3,j)',...
        quiverRatio);
    savedColor = h.Color;
    h.Color = color;
end
% TODO: (DEBUG) remove this for final release
daspect([1 1 1]);
axis([.8 .87 -.05 .05])

%% Manifolds

charDistance = EarthMoon.nondimensionalize([40, 100, 300, 1]*1e3, []);

figure(3)
clf reset

for manifold = 1:2; % 1 unstable, 2 stable manifolds
    tf = 1.4*pi;
    if(manifold == 2)
        tf = tf * -1;
    end
    
    nHalfManifolds = 2;
    periodicManifoldOrbit{manifold}(nHalfManifolds * nFixedPoints) = RC3BPOrbit;
    
    for i=1:nFixedPoints
        eigenvector = fixedPointEigv(:, manifold, i);
        unitaryDistanceV = eigenvector / norm(eigenvector(eigvElems));
        
        for j=1:nHalfManifolds
            if(j==1)
                halfManifold = 1;
            else
                halfManifold = -1;
            end
            perturbation = halfManifold * charDistance(4) * unitaryDistanceV;
            
            index = (i-1)*nHalfManifolds + j;
            periodicManifoldOrbit{manifold}(index) = RC3BPOrbit(fixedPointOrbit(i).xf + perturbation, EarthMoon);
            periodicManifoldOrbit{manifold}(index).integrateX0([0 tf]);
            JacobiErr = periodicManifoldOrbit{manifold}(index).getJacobi(periodicManifoldOrbit{manifold}(index).xf) - halfLyapunov.getJacobi(halfLyapunov.x0);
            %         disp(num2str(JacobiErr));
        end
    end
    
    % Closest flyby by the moon of the point opposite in the x axis
    if(manifold==1)
        index = 1;
    else
        index = 2;
    end
    deltaR = periodicManifoldOrbit{manifold}(index).x(1:3,:) - [1-EarthMoon.mu; 0; 0]*ones(1,size(periodicManifoldOrbit{manifold}(index).x(1:3,:),2));
    rNorm = sqrt(sum(deltaR.^2,1));
    rMin = min(rNorm);
    % Display results
    quiverRatio = .1;
    
    disp(['Closest approach to the Moon is: ' num2str(rMin) ' (' num2str(rMin*EarthMoon.l/1e3) ' km)'])
    hold on
    lyapunovOrbit.plot();
    if(manifold==1)
        color = unstableColor;
    elseif(manifold==2)
        color = stableColor;
    end
    h = quiver3(xf(1,:), xf(2,:), xf(3,:), ...
        pFixedPointEigv(:,1,j)', pFixedPointEigv(:,2,j)', pFixedPointEigv(:,3,j)',...
        quiverRatio);
    h.Color = color;
    h = quiver3(xf(1,:), xf(2,:), xf(3,:), ...
        -pFixedPointEigv(:,1,j)', -pFixedPointEigv(:,2,j)', -pFixedPointEigv(:,3,j)',...
        quiverRatio);
    savedColor = h.Color;
    h.Color = color;
    for i=1:nFixedPoints*nHalfManifolds
        h = periodicManifoldOrbit{manifold}(i).plot();
        h(3).Color = color;
    end
end

%% Test 1:
% Obtaining Lyapunov family of orbits

%% Setup workspace
clear
% Add path
path1 = [pwd,'/library'];
path2 = [pwd,'/targetingConstraints'];
addpath(path1, path2);

Constants

clc

%% Define systems
EarthMoon = RC3BPSystem(Earth, Moon);
deltaX0 = -[.001; 0; 0];

%% Find First guess
L2Vec = EarthMoon.lagrangePoint(2);
L2 = [L2Vec(1:2); zeros(2,1)];
orbit = RC3BPOrbit(L2, EarthMoon);

% In plane linear
A = orbit.stateMatrix(orbit.x0);
[V,D] = eig(A);

% Sort eigenvalues and rearrange eigendecomposition for Collinear points
stabilityIndices = diag(D);
P = zeros(4);
for i=1:4
    if(isreal(stabilityIndices(i)))
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

L23d = [L2Vec; zeros(3,1)];
delta3dX0 = [deltaX0(1:2); 0; deltaX0(3:4); 0.0];
orbit = RC3BPOrbit(L23d - delta3dX0, EarthMoon);
orbit.setStoppingCondition(@(t,x) (x(2)>0), 0);

orbit.integrateX0([0, inf]);
dC = orbit.jacobiAnalysis();

%% Find periodic orbit

% Targeter settings
targeterPeriodic = RC3BPTargeter();
targeterPeriodic.saveHistory = 1;
targeterPeriodic.epsilon = 1e-12;   % def 1e-12 % increasing it won't improve performance beyond integrator
[halfPeriodicOrbit, count] = targeterPeriodic.target(@lyapunovTargeting, orbit);

% Display results
display('Necessary initial conditions for Periodic Orbit:')
disp([halfPeriodicOrbit.x0; halfPeriodicOrbit.tf])

%% Propagate Lyapunov orbit
periodicOrbit = copy(halfPeriodicOrbit);
periodicOrbit.setStoppingCondition(@(t,x) (x(2)>0), -1);
periodicOrbit.integrateX0([0, inf]);
dcPeriodic = periodicOrbit.jacobiAnalysis();

periodicError = periodicOrbit.x0 -periodicOrbit.xf;
[periodicErrorR, periodicErrorV]  = EarthMoon.dimensionalize(periodicError(1:2), periodicError(3:4));

% Display results
figure(1)
clf reset
periodicOrbit.plot();

%% Compute monodromy matrix
% Monodromy auxiliary matrices
nDim = 3;
diagonal = (-1).^((1:6)-1);
G = diag(diagonal);
Omega = zeros(nDim);
Omega(2,1) = -1;
Omega(1,2) = 1;
H = [zeros(nDim)  -eye(nDim); eye(nDim) -2*Omega];

% From half orbit
STMhalf = halfPeriodicOrbit.STM();
monodromy = (G * H * STMhalf' *inv(H) * G) * STMhalf; %#ok<MINV>

% From whole orbit
monodromy2 = periodicOrbit.STM();

% Mondromy matrix analysis
monodromyErr = norm(monodromy-monodromy2);
[V,D] = eig(monodromy);
stabilityIndices =diag(D);

% Display results
display('monodromy matrix:')
disp(monodromy)
display('monodromy error wrt to full integration')
disp(monodromyErr)
display('eigenvalues of the monodromy matrix')
disp(stabilityIndices(1:6))

%% Find lyapunov orbit family
deltaX0Range(1) = halfPeriodicOrbit.x0(1);
deltaX0Range(2) = deltaX0Range(1)+.05;  % continuation parameter

% Find family of interest
familyLyapunov = PeriodicOrbitFamily(@lyapunovUpdate, @lyapunovTargeting);
familyLyapunov.continuationRange = deltaX0Range(1:2);
familyLyapunov.continuationNumber = round(abs(deltaX0Range(2)-deltaX0Range(1))/0.001)+1;
nOfOrbits = familyLyapunov.findFamily(orbit)
familyLyapunov.decimate(ceil(familyLyapunov.nOrbits/10))
% STM for the family
lyapStabilityIndices = familyLyapunov.getStabilityIndices();
% Orbit parameters
[~, x0WithT] = familyLyapunov.parameterDependence();
[timeVal, unit] = reasonableTime(x0WithT(end,1)*EarthMoon.t);
x0WithT(end,:) = (timeVal/x0WithT(end,1)) * x0WithT(end,:);
% Jacobi
for i=1:familyLyapunov.nOrbits
    C(i) = familyLyapunov.orbit(i).getJacobi(  familyLyapunov.orbit(i).x0 );
end

% Display results
display('Lypaunov family around L2 eigenvalues')
disp(lyapStabilityIndices)
figure(2)
clf reset
familyLyapunov.plot()
title('Lyapunov family of orbits')

figure(3)
clf reset
plot(x0WithT(1,:), x0WithT(5,:),'-x')
title('V_{y_0} as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel('V_{y_0} [non-dimensional units]')

figure(4)
clf reset
plot(x0WithT(1,:), x0WithT(7,:),'-x')
title('Period as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel(['Period [' unit ']'])

figure(5)
clf reset
plot(x0WithT(1,:), C ,'-x')
title('Jacobi Constant as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel(['C [non-dimensional units]'])

%% Find Halo Orbit Family

% select closes orbit to bifurcation as initial guess
[~, index] = min(vecnorm(lyapStabilityIndices(1:2,:)-1));
initialGuessHalo = familyLyapunov.orbit(index);

% Continuation range
Z0Range(1) = 1e-3;    % continuation parameter
Z0Range(2) = 0.2;        % for .2 and without decimation some stable Halo orbit are obtained
% the L2 halo family has a maximum in z, after this we should continue on x

% Find halo orbit family
familyHalo = PeriodicOrbitFamily(@haloUpdate, @haloTargeting);
familyHalo.continuationRange = Z0Range;
familyHalo.continuationNumber = round(abs(Z0Range(1)-Z0Range(2))/0.002)+1;
nOfOrbits = familyHalo.findFamily(initialGuessHalo)
familyHalo.decimate(ceil(familyHalo.nOrbits/10))
% STM for the family
haloStabilityIndices = familyHalo.getStabilityIndices();
% Orbit parameters
[~, x0WithT] = familyHalo.parameterDependence();
[timeVal, unit] = reasonableTime(x0WithT(end,1)*EarthMoon.t);
x0WithT(end,:) = (timeVal/x0WithT(end,1)) * x0WithT(end,:);
% Jacobi
for i=1:familyHalo.nOrbits
    C(i) = familyHalo.orbit(i).getJacobi(  familyHalo.orbit(i).x0 );
end

% Display results
display('Halo family around L2 eigenvalues')
disp(haloStabilityIndices)
figure(12)
clf reset
familyHalo.plot()
title('Halo family of orbits')

figure(13)
clf reset
plot(x0WithT(1,:), x0WithT(5,:),'-x')
title('V_{y_0} as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel('V_{y_0} [non-dimensional units]')

figure(14)
clf reset
plot(x0WithT(1,:), x0WithT(7,:),'-x')
title('Period as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel(['Period [' unit ']'])

figure(15)
clf reset
plot(x0WithT(1,:), C ,'-x')
title('Jacobi Constant as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel(['C [non-dimensional units]'])

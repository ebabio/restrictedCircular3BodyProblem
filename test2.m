%% Test 2:
% Obtaining Halo family of orbits given the initial conditions
% The family of halo orbits is propagated into the XY plane until reaching
% the bifurcation with the Lyapunov family

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

%% Starting conditions

x01 = [ .82575887090385;    0;                  0.08; ...
    0;                  .19369724986446;    0];

x02 = [ .82356490862838;    0;                  0.04; ...
    0;                  .14924319723734;    0];

%% Obtain Halo orbit

stoppingConditionPeriod = @(t,x) (x(2)>0);

% 1. X01
% a. Integrate given initial states
orbit1 = RC3BPOrbit(x01, EarthMoon);
orbit1.setStoppingCondition(stoppingConditionPeriod, -1);
orbit1.integrateX0(inf);

% b. Correct for periodic orbit
% Initial orbit guess
halo1 = RC3BPOrbit(x01, EarthMoon);
halo1.setStoppingCondition(stoppingConditionPeriod, 0);
% Setup targetting
targeterPeriodic = RC3BPTargeter();
targeterPeriodic.saveHistory = 1;
targeterPeriodic.epsilon = 1e-12;
[halo1half, count] = targeterPeriodic.target(@haloTargeting, halo1);
% Get final orbit
halo1.x0 = halo1half.x0;
halo1.tf = 2*halo1half.tf;
halo1.setStoppingCondition(stoppingConditionPeriod, -1);
halo1.integrateX0();

% Display
disp('First initial state')
disp('Error between initial and final states for given ICs: [nd]');
disp(orbit1.xf-orbit1.x0);
disp('Error between initial and final states for corrected ICs: [nd]');
disp(halo1.xf-halo1.x0);
figure(1)
clf reset
orbit1.plot();
figure(2)
clf reset
halo1.plot();

% 1. X02
% a. Integrate given initial states
orbit2 = RC3BPOrbit(x02, EarthMoon);
orbit2.setStoppingCondition(stoppingConditionPeriod, -1);
orbit2.integrateX0(inf);

% b. Correct for periodic orbit
% Initial orbit guess
halo2 = RC3BPOrbit(x02, EarthMoon);
halo2.setStoppingCondition(stoppingConditionPeriod, 0);
% Setup targetting
targeterPeriodic = RC3BPTargeter();
targeterPeriodic.saveHistory = 1;
targeterPeriodic.epsilon = 1e-12;
[halo2half, count] = targeterPeriodic.target(@haloTargeting, halo2);
% Get final orbit
halo2.x0 = halo2half.x0;
halo2.tf = 2*halo2half.tf;
halo2.setStoppingCondition(stoppingConditionPeriod, -1);
halo2.integrateX0();

% Display
disp('Second initial state')
disp('Error between initial and final states for given ICs: [nd]');
disp(orbit2.xf-orbit2.x0);
disp('Error between initial and final states for corrected ICs: [nd]');
disp(halo2.xf-halo2.x0);
figure(3)
clf reset
orbit2.plot();
figure(4)
clf reset
halo2.plot();

%% Obtain family of Halo Orbits

% Continuation range
Z0Range(1) = x01(3);    % continuation parameter
Z0Range(2) = 0;        % for .2 and without decimation some stable Halo orbit are obtained

% Find halo orbit family
familyHalo = PeriodicOrbitFamily(@haloUpdate, @haloTargeting);
familyHalo.continuationRange = Z0Range;
familyHalo.continuationNumber = round(abs(Z0Range(1)-Z0Range(2))/0.003)+1;
nOfOrbits = familyHalo.findFamily(halo1half)
familyHalo.decimate(10);
% STM for the family
lyapStabilityIndices = familyHalo.getStabilityIndices();
% Orbit parameters
[~, x0WithT] = familyHalo.parameterDependence();
[timeVal, unit] = reasonableTime(x0WithT(end,1)*EarthMoon.t);
x0WithT(end,:) = (timeVal/x0WithT(end,1)) * x0WithT(end,:);
for i=1:familyHalo.nOrbits
    C(i) = familyHalo.orbit(i).getJacobi(  familyHalo.orbit(i).x0 );
end

% Display results
display('Halo family converging to the xy-plane eigenvalues')
disp(lyapStabilityIndices)
figure(11)
clf reset
familyHalo.plot()
title('Halo family of orbits')

figure(12)
clf reset
plot(x0WithT(1,:), x0WithT(5,:),'-x')
title('V_{y_0} as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel('V_{y_0} [non-dimensional units]')

figure(13)
clf reset
plot(x0WithT(1,:), C ,'-x')
title('Jacobi Constant as a function of x_0')
xlabel('x_0 [non-dimensional units]')
ylabel(['C [non-dimensional units]'])
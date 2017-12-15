%% Test 3:
% Obtaining the surface of section for selected systems for the y=0 yDot>0
% map for a specific Jacobi constant. A stable orbit is located according to 
% the map and the initial conditions are found using a fixed Jacobi constant

%% Setup workspace
clear
% Add path
path1 = [pwd,'/library'];
path2 = [pwd,'/targetingConstraints'];
addpath(path1, path2);

Constants

clc

% Values for configuration:
% 1: Henon section of the copenhaguen problem
% 2: L1 open section of the copenhaguen problem
% 3: Henon section of the Earth Moon section
% 3: Henon section of the Earth Moon section with Lyapunov Orbits

configuration = 3;
%% Define problem
Copenhagen = RC3BPSystem(Earth, Moon);
xDotRange = 0;
if(configuration == 1)
    Copenhagen.mu = .5;
    Copenhagen.l = 1;
    Copenhagen.m = 1;
    Copenhagen.t = 1;
    
    sectionC = 4.5;
    tf = 50*pi;
    xRange = -.9:.05:-.1;
    xPeriodic = -.3;
elseif(configuration == 2)
    Copenhagen.mu = .5;
    Copenhagen.l = 1;
    Copenhagen.m = 1;
    Copenhagen.t = 1;
    
    sectionC = Copenhagen.potential(Copenhagen.lagrangePoint(1))+Copenhagen.potential(Copenhagen.lagrangePoint(2));
    tf = 30*pi;
    xRange = -1:.05:1;
    xPeriodic = -.70;
elseif(configuration == 3)
    sectionC = Copenhagen.potential(Copenhagen.lagrangePoint(1))+Copenhagen.potential(Copenhagen.lagrangePoint(2));
    tf = 50*pi;
    xRange = -1:.05:1.2;
    xPeriodic = .5;
elseif(configuration == 4)
    sectionC = 3.185176408375114;   % Lyapunov orbit at x = .83;
    tf = 20*pi;
    xRange = .75:0.01:.9;
    xDotRange = -.02:.01:0.02;
    xPeriodic = .84;
end

%% Create Phase Map for y=0,yDot>0 section


% Setup Map Grid
xRange(abs(xRange+Copenhagen.mu)<0.04) = [];  % remove small quantities
xRange(abs(xRange+Copenhagen.mu-1)<0.04) = [];  % remove small quantities

[xGrid, xDotGrid] = meshgrid(xRange, xDotRange);

% Setup Orbit
eventConditionPeriod = @(t,x) (x(2)>0);
templateOrbit = RC3BPOrbit(zeros(4,1), Copenhagen);
templateOrbit.setStoppingCondition(eventConditionPeriod, -1, 0);
templateOrbit.tf = tf;

% Create Surface of Section
tic
orbit(numel(xGrid)) = copy(templateOrbit);
n = numel(xGrid);
xSectionLoop = cell(1, n);
parfor i=1:n
    CExcedent = 2*Copenhagen.potential([xGrid(i); 0; 0]) - xDotGrid(i).^2 - sectionC;
    if(CExcedent>0)
        if(mod(i,ceil(n/10)) == 0)
            disp([num2str(round(i/n*100)) '% completed'])
        end
        yDot = sqrt(CExcedent);
        orbit(i) = copy(templateOrbit);
        orbit(i).x0 = [xGrid(i); 0+eps; xDotGrid(i); yDot];
        orbit(i).integrateX0();
        [~, msgid] = lastwarn;
        if(strcmp(msgid,'MATLAB:str2func:invalidFunctionName'))
            warning(['ICs [' num2str(xGrid(i)) ', ' num2str(xDotGrid(i)) '] do not satisfy tolerances'])
            lastwarn('')
            orbit(i).phi=0;
        else
            xSectionLoop{i} = orbit(i).odeSol.ye;
            orbit(i).phi=1;
        end
    end
end
% Remove NaNs
xSection = cat(2, xSectionLoop{:});
toc

figure(1)
clf reset
scatter(xSection(1,:), xSection(3,:),'.')
axis([-1 1.2 -2 2])

% figure(3)
% clf reset
% for i=1:n
%     if(orbit(i).phi == 1)
%         orbit(i).plot;
%     end
% end

%% Find periodic orbit for these conditions
% Initial guess come from inspection of Map
i = find(abs(xSection(1,:)-xPeriodic)<0.02 & abs(xSection(3,:))<0.02);
x0 = xSection(:,i(1));

stoppingConditionPeriod = @(t,x) (x(2)>0);
periodicHalf0 = RC3BPOrbit(x0, Copenhagen);
periodicHalf0.setStoppingCondition(stoppingConditionPeriod, -1);
periodicHalf0.integrateX0(inf);

% Setup targetting
periodicHalf0.setStoppingCondition(stoppingConditionPeriod, 0);
targeterPeriodic = RC3BPTargeter();
targeterPeriodic.saveHistory = 1;
targeterPeriodic.epsilon = 1e-12;
targetPeriodic = @(orbit) constantJacobiTargeting2d(orbit, sectionC);
[periodicHalf, count] = targeterPeriodic.target(targetPeriodic, periodicHalf0);

% Get final orbit
periodicOrbit = copy(periodicHalf);
periodicOrbit.tf = 2*periodicHalf.tf;
periodicOrbit.setStoppingCondition(stoppingConditionPeriod, -1);
periodicOrbit.integrateX0();

% Display
disp('Initial conditions for periodic orbit:')
disp(periodicOrbit.x0)

figure(2)
clf reset
periodicOrbit.plot(1)

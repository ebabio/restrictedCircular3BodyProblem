%% AAE 532 Table of Constants

%% Global constants

% as returned by google search
global G
G = 6.67259e-11;
global AU2meters
AU2meters = 149597870700;
J2000 = juliandate(datetime('2000-01-01 12:00:00'));

%% Orbit characteristics

% as listed in https://ssd.jpl.nasa.gov/

Sun.gravP   = 1.32712440018e20;
Sun.a       = 0;
Sun.radius = 695700e3;

Mercury.gravP   = G*0.330104e24;
Mercury.a   = 0.38709927*AU2meters;
Mercury.e = 0.20563593;

Venus.gravP     = G*4.86732e24;
Venus.a     = 0.72333566*AU2meters;
Venus.e = 0.00677672;
Venus.lonAscNode = deg2rad(76.68069);
Venus.inclination = deg2rad(3.39471);
Venus.lonPeriapsis = deg2rad(131.53298);
Venus.meanLon = deg2rad(181.97973);

Earth.gravP     = G*5.97219e24;
Earth.a     = 1.00000261*AU2meters; %Earth-Moon barycenter
Earth.e = 0.01671123;
Earth.radius = 6371e3;
Earth.lonAscNode = deg2rad(-11.26064);
Earth.inclination = deg2rad(0.00005);
Earth.lonPeriapsis = deg2rad(102.94719);
Earth.meanLon = deg2rad(100.46435);

Mars.gravP  = G*0.641693e24;
Mars.a      = 1.52371034*AU2meters;
Mars.e = 0.09339410;
Mars.radius = 3389.5e3;
Mars.lonAscNode = deg2rad(49.57854);
Mars.inclination = deg2rad(1.85061);
Mars.lonPeriapsis = deg2rad(336.04084);
Mars.meanLon = deg2rad(355.45332);

Jupiter.gravP	= G*1898.13e24;
Jupiter.a	= 5.20288700*AU2meters;
Jupiter.e = 0.04838624;
Jupiter.lonAscNode = deg2rad(100.55615);
Jupiter.inclination = deg2rad(1.30530);
Jupiter.lonPeriapsis = deg2rad(14.75385);
Jupiter.meanLon = deg2rad(34.40438);

Saturn.gravP	= G*568.319e24;
Saturn.a	= 9.53667594*AU2meters;
Saturn.e = 0.05386179;

Uranus.gravP	= G*86.8103e24;
Uranus.a	= 19.18916464*AU2meters;
Uranus.e = 0.04725744;

Neptune.gravP	= G*102.410e24;
Neptune.a	= 30.06992276*AU2meters;
Neptune.e = 0.00859048;

Pluto.gravP     = G*0.01309e24;
Pluto.a     = 39.48211675*AU2meters;
Pluto.e = 0.24882730;

% Satellite orbital elements from: https://ssd.jpl.nasa.gov/?sat_elem
Moon.gravP      = 4902.801e9;
Moon.radius = 1737.5e3;
Moon.a = 384400e3;
Moon.e = .0554;

Phobos.gravP    = .0007112e9;
Phobos.radius = 11.1e3;
Phobos.a = 9376e3;
Phobos.e = .0151;

Deimos.gravP = .0000985e9;
Deimos.radius = 6.2e3;
Deimos.a = 23458e3;
Deimos.e = 2e-4;

Ganymede.gravP = 9887.834e9;
Ganymede.radius = 2631.2e3;
Ganymede.a = 1070400e3;
Ganymede.e = 0.0013;

Titan.gravP = 8978.1382e9;
Titan.radius = 2574.73e3;
Titan.a = 1221865e3;
Titan.e = 0.0288;

Titania.gravP = 228.2e9;
Titania.radius = 788.9e3;
Titania.a = 436300e3;
Titania.e = 1.1e-3;

Triton.gravP = 1427.6e9;
Triton.radius = 1353.4e3;
Triton.a = 354759e3;
Triton.e = 0;


Charon.gravP = 102.3e9;
Charon.radius = 603.6e3;
Charon.a = 19591e3;
Charon.e = 0.0002;
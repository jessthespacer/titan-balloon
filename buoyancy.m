close all; clear all;

% Objective: Minimize total system mass

% Physical constants
sigm = 5.670373E-8;				% [W/m^2 K^4]
Ru = 8314.46261815324;			% [J/kmol K]

% Payload requirements
payloadMass = 10E+3;			% [kg]

% Powerplant parameters
specPowThermal = 28.87;			% [Wth/kg]
specPowElec = 6.7;				% [Wth/kg]

% Envelope parameters
emissivity = 0.1;				% [dimensionless]
arealDensity = 0.1;				% [kg/m^2]
heatTransCoeff = 1;				% [W/K m^2]

% Operational parameters
opCeiling = 10000;				% [m]

% Planetary parameters
g = 1.4;						% [m/s^2]

% Atmospheric parameters
Matm = 27.3;					% [kg/kmol]
Tatm = 98.29;					% [K]
P0 = 150000;					% [kPa]

Ratm = Ru / Matm;				% [J/kg K]
scaleH = (Ratm * Tatm / g);		% [m]

P = @(z) P0 * exp(-z / scaleH);	% [kPa]
rho = @(z, T) P(z) ./ (Ratm * T);	% [kg/m^3]

% Variables for parametric study
deltaT = 0.1:0.1:10;			% [K]

% Heat loss
qConvect = @(dT) dT .* heatTransCoeff;	% [W/m^2]
qRadiate = @(dT) emissivity .* sigm .* ((Tatm + dT).^4 - Tatm.^4);	% [W/,^2]
qLoss = @(dT) qConvect(dT) + qRadiate(dT);		%[W/m^2]

% Calculate balloon volume
envArea = @(r) 4 * pi * r.^2;		% [m^2]
envVol = @(r) 4 / 3 * pi * r.^3;	% [m^3]
sysMass = @(r, dT) envArea(r) .* (qLoss(dT) ./ specPowThermal + arealDensity) + payloadMass;	% [kg]

liftWeight = @(dT) rho(opCeiling, Tatm) - rho(opCeiling, Tatm + dT);	% [kgf/m^3]
envLift = @(r, dT) liftWeight(dT) .* envVol(r);		% [kgf]

liftDiff = @(r, dT) abs(envLift(r, dT) - sysMass(r, dT));		% [kgf]

N = size(deltaT, 2);
rad = zeros(1, N);
mass = zeros(1, N);

% Calculate required balloon volumes
for ii = 1:N
	deltT = deltaT(ii);
	fun = @(r) liftDiff(r, deltT);
	rad(ii) = fminunc(fun, 1000);
	mass(ii) = sysMass(rad(ii), deltT);
end

[minmass, ix] = min(mass);
disp(["Optimal deltaT: " num2str(deltaT(ix)) " K"]);
disp(["Minimum system mass: " num2str(mass(ix)) " kg"]);
disp(["Optimal balloon radius: " num2str(rad(ix)) " m"]);
disp(["For altitude: " num2str(opCeiling) " m"]);

% figure(1);
% plot(deltaT, rad);
% xlabel("deltaT [K]");
% ylabel("radius [m]");

% figure(2);
% plot(deltaT, mass);
% xlabel("deltaT [K]");
% ylabel("mass [kg]");
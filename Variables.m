% Fluid properties
visc = 0.01;        %*0.0001;
rho = 1;            %g/cm^3;
f=1    % Difference in contractile tension

% Geometry
Ro0 = 721/2*0.0001;
H0 = 62*0.0001;
L =0.3046;

H= H0;
Ro= Ro0;
Ri = Ro - H;

lz = 1.15;


% Fluid properties
visc = 0.01;        %*0.0001;

% Valve resistance parameters (Bertram et al 2013)
% RVn = 6*10^4;          % (dynes/cm^2)/(ml/sec) 
% RVx = 12*10^7;      % (dynes/cm^2)/(ml/sec)
RVn = 6*10^4;          % (dynes/cm^2)/(ml/sec) 
RVx =12*10^8;
Po =0;
so=0.2;

% Passive Ligation &active
b =1.5273e+04;
b11=2.6458e+04;
b21 =4.8884;
b12 =4.6117e+03;
b22 =7.1045;
b13 = 8.8856e+03;
b23 =15.4786;
alpha1 =1.5708;
alpha2 =0;
alpha3 =0.7603;
alpha4 =-0.7603;
lq_max =1.15;
lq_o =0.6119;
lq_high = 2*lq_max-lq_o;
slope0 = 13.0981;
int = -4574.4;
T_tonic=(1.015294697190801e+04)*f;
T_phasic=2.2208e+05*f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%frequency%%%%%%%%%%%%%%%%%%%
ww=2;   %contraction time
tr=28; % 28 %2.8 referactory time
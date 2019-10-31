
function [p_mod, fz, lq, lr, tq, tq_pas, tq_act, tz, h] = cmodel(ri,Tact,pex)
% Fluid properties
visc = 0.01;        %*0.0001;
rho = 1;            %g/cm^3;


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
RVx = 12*10^8;
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
T_tonic=1.015294697190801e+04;
T_phasic=2.322377944095554e+05-T_tonic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%frequency%%%%%%%%%%%%%%%%%%%
ww=2;   %contraction time
tr=28; % 28 %2.8 referactory time
% lz = 1.15;
% samp = 1; %Number of vessels in the chain
% 
% % Applied loads
% % pa = 2*980;  %2.1916*980;    % dynes/cm^2 (1 cmH2O = 980 dynes/cm^2);
% 
% % Fluid properties
% visc = 0.01;        %*0.0001;
% rho = 1;            %g/cm^3;
% 
% % Valve resistance parameters (Bertram et al 2013)
% RVn = 600*100;          % (dynes/cm^2)/(ml/sec) 
% RVx = 12*10^7;      % (dynes/cm^2)/(ml/sec)
% Po =0;
% so=0.2; % non-dimensioanl
% 
% 
%     
%    lz = 1.15;
% samp = 1; %Number of vessels in the chain
% 
% 
% 
% % Fluid properties
% visc = 0.01;        %
% 
% 
% 
% 
%     
% % Active constitutive model parameters (Alex Caulk  2015)
% 
% lq_max =1.15;
% lq_o =0.6119;
% lq_high = 2*lq_max-lq_o;
% slope0 = 13.0981;
% int = -4574.4;
% 
% % Vessel geometry
% Ro0 = 804/2*0.0001;
% H0 = 58*0.0001*1.0;
% L = 0.3046;
% H= H0;
% Ro= Ro0;
% Ri = Ro - H;
% % Passive
% b =1.9235e+04;
% b11 =1.9083e+03;
% b21 =27.2123;
% b12 = 4.7515e+04;
% b22 = 4.8137;
% b13 = 2.7769e+03;
% b23 =34.9830;
% alpha1 =1.5708;
% alpha2 = 0;
% alpha3 =  0.8492;
% alpha4 =  -0.8492;
% 
% 
% %%%%T act Values%%%%%%%%
% T_tonic=30000/2;
% % T_phasic=80000;
% % TT=486348.7903/5;
% % b=0.05
% % T_tonic=38000+0*TT*b/(1+b);
%  T_phasic=225000;
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ww=2; % contraction time
% tr=6;% %referactory time (11.3;%2.62 )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i = 1:length(ri)
    %Calculate loaded thickness from incompressibility

    coef_c = -1*H(i)*(2*Ri(i)+H(i))/lz;
    coef_b = 2*ri(i);
    coef_a = 1;
    D1 = [coef_a, coef_b, coef_c];
    E1 = roots(D1);
    F1 = E1(2);
    h(i) = F1';

    % Kinematics
    lq(i) = (ri(i)+h(i)/2)/(Ri(i)+H(i)/2);
    lr(i) = 1/(lq(i)*lz);
    Fqq = lq(i);
    Fzz = lz;
    Frr = lr(i);
    Crr = Frr^2;
    Cqq = Fqq^2;
    Czz = Fzz^2;
    Err = (Crr-1)/2;
    Eqq = (Cqq-1)/2;
    Ezz = (Czz-1)/2;

    lambda2_k1 = Cqq*(sin(alpha1))^2 + Czz*(cos(alpha1))^2;
    lambda2_k2 = Cqq*(sin(alpha2))^2 + Czz*(cos(alpha2))^2;
    lambda2_k3 = Cqq*(sin(alpha3))^2 + Czz*(cos(alpha3))^2;
    lambda2_k4 = Cqq*(sin(alpha4))^2 + Czz*(cos(alpha4))^2;

    %Strain energy function proposed by Holzapfel for four fiber families
    tr_hat = 2*Crr*b;
    tq_hat = 2*Cqq*(b + b11/2*exp(b21*(lambda2_k1-1)^2)*(Cqq*(sin(alpha1))^4 + Czz*(sin(alpha1))^2*(cos(alpha1))^2 - (sin(alpha1))^2) + ...
        b12/2*exp(b22*(lambda2_k2-1)^2)*(Cqq*(sin(alpha2))^4 + Czz*(sin(alpha2))^2*(cos(alpha2))^2 - (sin(alpha2))^2) + ...
        b13*exp(b23*(lambda2_k3-1)^2)*(Cqq*(sin(alpha3))^4 + Czz*(sin(alpha3))^2*(cos(alpha3))^2 - (sin(alpha3))^2));
    tz_hat = 2*Czz*(b + b11/2*exp(b21*(lambda2_k1-1)^2)*(Cqq*(sin(alpha1))^2*(cos(alpha1))^2 + Czz*(cos(alpha1))^4 - (cos(alpha1))^2) + ...
        b12/2*exp(b22*(lambda2_k2-1)^2)*(Cqq*(sin(alpha2))^2*(cos(alpha2))^2 + Czz*(cos(alpha2))^4 - (cos(alpha2))^2) + ...
        b13*exp(b23*(lambda2_k3-1)^2)*(Cqq*(sin(alpha3))^2*(cos(alpha3))^2 + Czz*(cos(alpha3))^4 - (cos(alpha3))^2));

    % Calculate active stress
    tq_pas(i) = tq_hat - tr_hat;
    tz(i) = tz_hat - tr_hat;

    % Calculate active stress
    if lq(i) <= lq_high && lq(i) >= lq_o;
        tq_act(i) = Tact(i)*lq(i)*(1-((lq_max-lq(i))/(lq_max-lq_o))^2);
    else
        tq_act(i) = 0;
    end
    


    tq(i) = tq_act(i) + tq_pas(i);

    p_mod(i) = h(i)/ri(i)*(tq(i) + pex) + pex;
    fz(i) = pi*h(i)*(2*ri(i)+h(i))*tz(i);

end


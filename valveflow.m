function [Q1,Q2] = valveflow(pa,pb,p1,p2,nlymph)
Variables;
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

if nlymph==1
   Q1(1)=(pa-p1(1))./(RVn + RVx.*(1./(1+exp(so.*((pa-p1(1))-Po)))));
   Q2(1)=(p2(1)-pb)./(RVn + RVx.*(1./(1+exp(so.*((p2(1)-pb)-Po)))));
else
for k=1:nlymph
       
    if k==1
    Q1(k)=(pa-p1(k))./(RVn + RVx.*(1./(1+exp(so.*((pa-p1(k))-Po)))));
    Q2(k)=(p2(k)-p1(k+1))./(RVn + RVx.*(1./(1+exp(so.*((p2(k)-p1(k+1))-Po)))));    
   
    elseif k==nlymph 
Q1(k)=(p2(k-1)-p1(k))./(RVn + RVx.*(1./(1+exp(so.*((p2(k-1)-p1(k))-Po)))));
Q2(k)=(p2(k)-pb)./(RVn + RVx.*(1./(1+exp(so.*((p2(k)-pb)-Po)))));    
    
    else
Q1(k)=(p2(k-1)-p1(k))./(RVn + RVx.*(1./(1+exp(so.*((p2(k-1)-p1(k))-Po)))));
Q2(k)=(p2(k)-p1(k+1))./(RVn + RVx.*(1./(1+exp(so.*((p2(k)-p1(k+1))-Po)))));    
    end
end
end







end




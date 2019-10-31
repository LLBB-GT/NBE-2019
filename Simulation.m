%**************Lymphangion G&R Code******************
%***************** V3*** 11/20/2016********************
for l=1:1
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
% lz = 1.15;
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
i=1;
PEX=0;
nlymph=1
while i<80
    
pa=6+0*(i-1)*0.5;%inlet pressure
pb=7+(i-1)*1;%outlet pressure
inlet_pressure(i)=pa;
outlet_pressure(i)=pb;
factor=1%0*(1+0.2*(i-1));
[mean_Q2,Qout,mean_tau_w ,mean_sigma_w,mean_diameter,Q1,Q2,r,Q_1,Q_2] = Main(pa,pb,nlymph,factor,PEX)
k=1:nlymph;
mWSS(i,k)=mean_tau_w(k); %Mean Wall Shear Stress ()
mHoop_Stress(i,k)=mean_sigma_w(k); %Mean Hoop Stress ()
% mTransmural_Pressure(i)=mean_P;  % Mean Transmural Pressure [mmHg]
% mFlow(i)=mean_Q; % Mean Flow Rate (mL/min)

% mFlow_v1(i)=mean_Q1; % Mean Flow Rate 2nd Valve(mL/min)
mFlow_v2(i)=mean_Q2 ;% Mean Flow Rate 1st Valve(mL/min)
Qoutflow(:,i)=Qout
if mFlow_v2(i)<-0.00001%mFlow_v2(i)<0.00001
    break
end
mD(i,k)=mean_diameter(k);
% mSD(i)=systolic_diameter;
% mDD(i)=diastolic_diameter;
i=i+1
end
%  i=1:40
%  mFlow_v2(i)=mean(Qoutflow(:,i))

colors = {'k','b','k','g','y','c','m'};
types = {'o','*','+','s','x','.','d'};
linetypes={'-','--',':','-','--'};
 m=1
msize=3
figure (8)
hold on
plot((outlet_pressure-inlet_pressure),mWSS,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
    'FontName','Times New Roman');
ylabel('Mean WSS (dyne/cm^2)','FontWeight','bold','FontSize',20,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
% % legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')
% 
% figure (9)
% hold on
% plot((outlet_pressure-inlet_pressure),mFlow,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
% xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% ylabel('Mean Flow (mL/min)','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
% 'LineWidth',1.5)
% figure (10)
% hold on
% plot((outlet_pressure-inlet_pressure),mTransmural_Pressure,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
% xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% ylabel('Mean Transmural Presssure [mmHg]','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
% 'LineWidth',1.5)
figure (11)
hold on
plot((outlet_pressure-inlet_pressure),mHoop_Stress,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
    'FontName','Times New Roman');
ylabel('Mean Total Hoop Stress (dyne/cm^2)','FontWeight','bold','FontSize',20,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
% figure (12)
% hold on
% plot((outlet_pressure-inlet_pressure),mD,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
% xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% ylabel('Mean Diameter (\mum)','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
% 'LineWidth',1.5)
% 
% figure (13)
% hold on
% plot((outlet_pressure-inlet_pressure),mSD,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
% % xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
% %     'FontName','Times New Roman');
% % ylabel('Mean Systolic Diameter (\mum)','FontWeight','bold','FontSize',20,...
% %     'FontName','Times New Roman');
% % set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
% % 'LineWidth',1.5)
% % 
% 
% 
% %figure (14)
% hold on
% plot((outlet_pressure-inlet_pressure),mDD,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
% xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% ylabel('Mean Diastolic& Systolic Diameters (\mum)','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
% 'LineWidth',1.5)
% 
% figure (15)
% hold on
% plot((outlet_pressure-inlet_pressure),mFlow_v1,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
% xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% ylabel('1st Valve Mean Flow (mL/min)','FontWeight','bold','FontSize',20,...
%     'FontName','Times New Roman');
% set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
% 'LineWidth',1.5)
hold;
figure (16)
hold on
plot((outlet_pressure-inlet_pressure),mFlow_v2,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('(Outflow Pressure)-(Inflow Pressure) [mmHg]','FontWeight','bold','FontSize',15,...
    'FontName','Times New Roman');
ylabel('2nd Valve Mean Flow (mL/min)','FontWeight','bold','FontSize',15,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)

if l==1
 save('n=1H=100pe00')
    elseif l==2
 save('n=1H=100pe05')
    elseif l==3
 save('n=1H=100pe10')
     elseif l==4
 save('n=1H=100pe15')
     elseif l==5
         save('n=1H=100pe20')
 end
end
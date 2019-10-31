function [mean_Q2,Qout,mean_tau_w ,mean_sigma_w,mean_diameter,Q1,Q2,r,Q_1,Q_2,Volumeo,p_cmH2O] = Main_Ligation(pa_mmhg,pb_mmhg,nlymph,~,PEX)
% 

Variables;
% lz = 1.15;
% samp = 1; %Number of vessels in the chain

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
% tr=11.3;% %referactory time (11.3;%2.62 )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0(1:nlymph)=.02 ;%Initial Guess for Diameter(cm)
pa=pa_mmhg*1333.2239;
pb=pb_mmhg*1333.2239; %inflow&outflow pressure [mmHg]
%N=201% number of time steps
dt=0.01;%2*ww/(N-1); %time step
NN=4%number of cycle
N=NN*(ww+tr)/dt+1;
l=1:ww/dt+1;
m=1:nlymph
phi=0.1
Tp(l,m)=( 1- cos(phi*(m-1)+2*pi*1/ww*(dt*(l-1))))*T_phasic/2+ T_tonic
Tt(1:tr/dt,1:nlymph)=T_tonic
Ta=[Tp' Tt' Tp' Tt' Tp' Tt' Tp' Tt']
Ta=Ta'
 
for i=1:N
    if i==1
        k=1:nlymph
p_mod=(pa+pb)/2;
p_mod(k)=pa+(pb-pa)*k/(nlymph*2) ;%%%%%%%%%%%%%%%%%%%%%%%this parameter dictates initial diameter-pressure
Tact(k)=T_tonic;
pex(k)=PEX;
f = @(R)Diameter_Initial(R,nlymph,Tact,p_mod,pex);
[R,fval] = fsolve(f,x0);
for k=1:nlymph
pm(k)=cmodel(R(k),Tact(k),pex(k));
Rf(k)=8*visc*L/(pi*R(k)^4);
end
r(i,:)=R;

for k=1:nlymph
    if nlymph==1
    xo(1)=(pa+pm(1))/2;    
    xo(2)=(pb+pm(1))/2;
    else
    if (k==1)
        xo(1)=(pa+pm(1))/2;    
        xo(2)=(pb+pm(nlymph))/2;
    elseif (k==nlymph)
        xo(2*k-1)=(pm(k)+pm(k-1))/2;
        xo(2*k)=(pm(k)+pb)/2;
    else 
        xo(2*k-1)=(pm(k)+pm(k-1))/2;
        xo(2*k)=(pm(k)+pm(k+1))/2;
    end
    end
end
xo;
Rf;
% x0 = [(pa+pm)/2,(pb+pm)/2];
f = @(p)fandv(p,pm,pa,pb,Rf,nlymph);
% x0 = [0,0];
% x = fsolve(fun,x0)fandv(p,pm,pa,pb,Rf)
[x,fval]=fsolve(f,xo);
% [y,ffval]=lsqnonlin(f,x0);

% g = @(p)fandvb(p,pm,pa,pb,Rf);
% x0 = [0,0];
% % x = fsolve(fun,x0)fandv(p,pm,pa,pb,Rf)
% [z,fval] = fsolve(g,x0);
% [w,ffval]=lsqnonlin(g,x0);
% p1=z(1);
% p2=z(2);
for k=1:nlymph
p1(i,k)=x(2*k-1);
p2(i,k)=x(2*k);
end
[Q1,Q2] = valveflow(pa,pb,p1,p2,nlymph);
% p1=x(1);
% p2=x(2);
% [Q11,Q22] = valveflow(pa,pb,p1,p2);
time(i)=0;

for k=1:nlymph
p_t(i,k)=cmodel(R(k),Tact(k),pex(k));
T(i,k)=Tact(k);
end

    else
         i
         
        ddt=dt*(i-1);
        time(i)=ddt;
%     ( 1-  cos(2*pi*ww*(ddt)))
phi=0.1
 for k=1:nlymph;
     Tact(k)=Ta(i,k)
%       Tact(k)=( 1- cos(phi*(k-1)+2*pi*1/ww*(ddt)))*T_phasic/2+ T_tonic;
      T(i,k)=Tact(k);
        pex(k)=0*( 1-  cos(phi*(k-1)+2*pi*1/ww*(ddt)))*PEX+PEX;
        Pex(i,k)=pex(k);
 end
      
     for  k=1:nlymph
       R(i,k)=r(i-1,k) ;
       Ro(k)= R(i,k);
       x1(3*k-1) = p1(i-1,k);  
       x1(3*k-2) = p2(i-1,k);  
       x1(3*k) = R(i,k);    
     end
      f = @(p)diameterp(p,pa,pb,Tact,Ro,dt,nlymph,pex);
      [x,fval] = fsolve(f,x1);
       k=1:nlymph     
      r(i,k)= x(3*k);
      p1(i,k)=x(3*k-2);
      p2(i,k)=x(3*k-1);
      p11(k)=p1(i,k);
      p22(k)= p2(i,k);
      R_t(i,k)=r(i,k);
      rr=r(i,k);
      for k=1:nlymph
   %   p_t(i,k)=cmodel(r(i,k),T(i,k),pex(k));
      [p_t(i,k),fz(i,k), lq(i,k), lr(i,k), tq(i,k), tq_pas(i,k), tq_act(i,k), tz(i,k), h(i,k)]=cmodel(r(i,k),T(i,k),pex(k));

      end
%     [p_t(i),fz(i), lq(i), lr(i), tq(i), tq_pas(i), tq_act(i), tz(i), h(i)]=cmodel(R_t,Tact);

      [Q1,Q2] = valveflow(pa,pb,p11,p22,nlymph);
%       qq(i)=Q1
%       qqq(i)=Q2
       k=1:nlymph
       Q_1(i,k)=Q1(1,k);
       Q_2(i,k)=Q2(1,k);
       Q_m(i,k)=(Q_1(i,k)+ Q_2(i,k))*0.5;
    
% figure (30)
% plot(Q_1)
% 
% figure (31)
% plot(Q_2)

end
% % figure (29)
% % plot(r)
% 
% figure (32)
% plot(p_t)

end


k=1:nlymph
wss(:,k) = 4*visc*Q_m(:,k)./(pi*r(:,k).^3);
mean_tau_w(k) = mean(wss(201:N,k)); %
mean_sigma_w(k) = mean(tq(201:N,k)); %
Volume=pi*r(201:N,k).^2*L*lz;
Volumeo=pi*r(:,k).^2*L*lz;
mean_Q2 = mean(60*Q_2(201:N,nlymph)); % Mean Flow Rate (mL/min)
Qout=60*Q_2(201:N,nlymph);

colors = {'k','b','k','g','y','c','m'};
types = {'o','*','+','s','x','.','d'};
linetypes={'-','--',':','-','--'};
m=1;
msize=0.3;
hold on
figure (1)
r_micron=10000*r;
k=1:nlymph
mean_diameter(k)=2*mean(r_micron(201:N,k));
plot(time,2*r_micron,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('Time (s)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('Diameter (\mum)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5);
legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')
% legend()
% % legend()

hold on
figure (2)
plot(time,T/10000,'LineWidth',1.5,'Color',colors{m},'LineStyle','-');
xlabel('Time (s)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('Activation Parameter(kdyne/cm^2)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)



hold on
figure (3)
p_cmH2O=p_t/980.665;
plot(time,p_cmH2O,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('Time (s)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('Pressure(cmH_2O)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
legend('Q_{avg}(Collapsible)')
legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')
% legend('3 point Model (Bertram 2014)')
% legend('Poiseuille')



hold on
figure (4)
%Q_1 (mL/s)
plot(time,Q_1*60,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('Time (s)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('1st Valve Flow Rate(mL/min)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')
% legend('3 point Model (Bertram 2014)')
% legend('Poiseuille')


hold on
figure (5)
plot(time,Q_2*60,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('Time (s)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('2nd Valve Flow Rate(mL/min)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')
% legend('3 point Model (Bertram 2014)')
% legend('Poiseuille')


hold on
figure (6)
hold on
plot(time,Q_m*60,'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('Time (s)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('Average Flow Rate(mL/min)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')

hold on

figure (7)
hold on
plot(Volumeo(401:size(Volumeo,1)),p_cmH2O(401:size(Volumeo,1)),'LineWidth',1.5,'Color',colors{m},'LineStyle',linetypes{m},'MarkerSize',msize,'Marker',types{m});
xlabel('Volume(mL)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
ylabel('Pressure(cmH_2O)Volume(mL)','FontWeight','bold','FontSize',25,...
    'FontName','Times New Roman');
set(gca,'FontName','Helvica','FontSize',13,'FontWeight','bold',...
'LineWidth',1.5)
legend('Collapsible','3 point Model (Bertram 2014)','Poiseuille')

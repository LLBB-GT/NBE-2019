
function F = diameterp(p,pa,pb,Tact,R0,dt,nLymph,pex)
Variables;


for i=1:nLymph
if nLymph == 1
    F(1)=cmodel((p(3)/1),Tact(1),pex(1))-(p(1)+p(2))/2;
    F(2)=p(2)-p(1)+(8*visc*L/(pi*(p(3)/1)^4))*(((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-pb)/(RVn + RVx*(1/(1+exp(so*((p(2)-pb)-Po)))))))/2;
    F(3)=(L*pi*((p(3)/1)^2-R0(1)^2)/dt+(-1*((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-pb)/(RVn + RVx*(1/(1+exp(so*((p(2)-pb)-Po))))))));

elseif nLymph==2
    F(1)=cmodel((p(3)/1),Tact(1),pex(1))-(p(1)+p(2))/2;
    F(2)=p(2)-p(1)+(8*visc*L/(pi*(p(3)/1)^4))*(((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-p(4))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(4))-Po)))))))/2;
    F(3)=L*pi*((p(3)/1)^2-R0(1)^2)+dt*(-1*((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-p(4))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(4))-Po)))))));
    F(4)=cmodel((p(6)/1),Tact(2),pex(2))-(p(4)+p(5))/2;
    F(5)=p(5)-p(4)+(8*visc*L/(pi*(p(6)/1)^4))*(((p(2)-p(4))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(4))-Po))))))+((p(5)-pb)/(RVn + RVx*(1/(1+exp(so*((p(5)-pb)-Po)))))))/2;
    F(6)=L*pi*((p(6)/1)^2-R0(2)^2)+dt*(-1*((p(2)-p(4))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(4))-Po))))))+((p(5)-pb)/(RVn + RVx*(1/(1+exp(so*((p(5)-pb)-Po)))))));

else
    
        fac=1;
     F(1)=fac*(cmodel((p(3)/1),Tact(1),pex(1))-(p(1)+p(2))/2);
    F(2)=fac*(p(2)-p(1)+(8*visc*L/(pi*(p(3)/1)^4))*(((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-p(4))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(4))-Po)))))))/2);
    F(3)=fac*((L*pi*((p(3)/1)^2-R0(1)^2)+dt*(-1*((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-p(4))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(4))-Po)))))))));
    
for i=2:nLymph-1
 
    F(3*i-2)=fac*(cmodel((p(3*i)/1),Tact(i),pex(i))-(p(3*i-2)+p(3*i-1))/2);
    F(3*i-1)=fac*(p(3*i-1)-p(3*i-2)+(8*visc*L/(pi*(p(3*i)/1)^4))*(((p(3*i-4)-p(3*i-2))/(RVn + RVx*(1/(1+exp(so*((p(3*i-4)-p(3*i-2))-Po))))))+((p(3*i-1)-p(3*i+1))/(RVn + RVx*(1/(1+exp(so*((p(3*i-1)-p(3*i+1))-Po)))))))/2);
    F(3*i)=fac*((L*pi*((p(3*i)/1)^2-R0(i)^2)+dt*(-1*((p(3*i-4)-p(3*i-2))/(RVn + RVx*(1/(1+exp(so*((p(3*i-4)-p(3*i-2))-Po))))))+((p(3*i-1)-p(3*i+1))/(RVn + RVx*(1/(1+exp(so*((p(3*i-1)-p(3*i+1))-Po)))))))));
 
end
  F(3*nLymph-2)=fac*(cmodel((p(3*nLymph)/1),Tact(nLymph),pex(nLymph))-(p(3*nLymph-2)+p(3*nLymph-1))/2);
    F(3*nLymph-1)=fac*(p(3*nLymph-1)-p(3*nLymph-2)+(8*visc*L/(pi*(p(3*nLymph)/1)^4))*(((p(3*nLymph-4)-p(3*nLymph-2))/(RVn + RVx*(1/(1+exp(so*((p(3*nLymph-4)-p(3*nLymph-2))-Po))))))+((p(3*nLymph-1)-pb)/(RVn + RVx*(1/(1+exp(so*((p(3*nLymph-1)-pb)-Po)))))))/2);
    F(3*nLymph)=fac*((L*pi*((p(3*nLymph)/1)^2-R0(nLymph)^2)+dt*(-1*((p(3*nLymph-4)-p(3*nLymph-2))/(RVn + RVx*(1/(1+exp(so*((p(3*nLymph-4)-p(3*nLymph-2))-Po))))))+((p(3*nLymph-1)-pb)/(RVn + RVx*(1/(1+exp(so*((p(3*nLymph-1)-pb)-Po)))))))));

    
end
end
end


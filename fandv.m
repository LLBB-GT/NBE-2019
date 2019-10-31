function F = fandv(p,pm,pa,pb,Rf,nlymph)
Variables;

if nlymph==1
     F(1)=pm(1)-(p(1)+p(2))/2;
     F(2)=p(2)-p(1)+Rf(1)*(((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-pb)/(RVn + RVx*(1/(1+exp(so*((p(2)-pb)-Po)))))))/2;
else
for k=1:nlymph
       
    if k==1
        F(1)=pm(1)-(p(1)+p(2))/2;
        F(2)=p(2)-p(1)+Rf(1)*(((pa-p(1))/(RVn + RVx*(1/(1+exp(so*((pa-p(1))-Po))))))+((p(2)-p(3))/(RVn + RVx*(1/(1+exp(so*((p(2)-p(3))-Po)))))))/2;
    elseif k==nlymph 
    F(2*k-1)=pm(k)-(p(2*k-1)+p(2*k))/2;
    F(2)=p(2*k)-p(2*k-1)+Rf(k)*(((p(2*k-2)-p(2*k-1))/(RVn + RVx*(1/(1+exp(so*((p(2*k-2)-p(2*k-1))-Po))))))+((p(2*k)-pb)/(RVn + RVx*(1/(1+exp(so*((p(2*k)-pb)-Po)))))))/2;
    
    else
    F(2*k-1)=pm(k)-(p(2*k-1)+p(2*k))/2;
    F(2*k)=p(2*k)-p(2*k-1)+Rf(k)*(((p(2*k-2)-p(2*k-1))/(RVn + RVx*(1/(1+exp(so*((p(2*k-2)-p(2*k-1))-Po))))))+((p(2*k)-p(2*k+1))/(RVn + RVx*(1/(1+exp(so*((p(2*k)-p(2*k+1))-Po)))))))/2;
    
    end
end
end

    
    
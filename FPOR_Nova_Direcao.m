function [dxn,dzn,detan,dlambdan,alfa_primal_n,alfa_dual_n] = FPOR_Nova_Direcao(w1,w2,z,lambda,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda,sigmaprimal,sigmadual);


dxn = (w1*dxp + w2*dx);

dzn = (w1*dzp + w2*dz); 

detan = (w1*detap + w2*deta);

dlambdan = (w1*dlambdap + w2*dlambda);

cand_z_n = 1;
cand_lamb_n = 1;


n1 = 1; n2 = 1;

n = length(z);

for i=1:n
    
    if z(i)>0 & dzn(i)<0
        cand_z_n(n1) = min(-z(i)/dzn(i),1);
        n1 = n1+1;
    end
    
    if lambda(i)>0 & dlambdan(i)<0
        cand_lamb_n(n2) = min(-lambda(i)/dlambdan(i),1);
        n2 = n2+1;
    end
    
end

alfa_primal_n = (min(cand_z_n))*sigmaprimal;
alfa_dual_n = (min(cand_lamb_n))*sigmadual;



end
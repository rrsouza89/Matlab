function [alfa_primal_prev,alfa_primal_cor,alfa_dual_prev,alfa_dual_cor,sigmaprimal,sigmadual] = FPOR_Calculo_Passos(z,lambda,dzp,dlambdap,dz,dlambda,G,H,mi);

tol1 = 0;
tol = 0;

cand_z_prev = 1;
cand_z_cor = 1;

cand_lamb_prev = 1;
cand_lamb_cor = 1;

n1 = 1; n2 = 1; n3 = 1; n4 = 1;

n = length(z);

for i=1:n
    
    if z(i)>tol1 & dzp(i)< tol
        cand_z_prev(n1) = min(-z(i)/dzp(i),1);
        n1 = n1+1;
    end
    
    if z(i)>tol1 & dz(i)< tol
        cand_z_cor(n2) = min(-z(i)/dz(i),1);
        n2 = n2+1;
    end
    
    if lambda(i)>tol1 & dlambdap(i)< tol
        cand_lamb_prev(n3) = min(-lambda(i)/dlambdap(i),1);
        n3 = n3+1;
    end
    
    if lambda(i)>tol1 & dlambda(i)< tol
        cand_lamb_cor(n4) = min(-lambda(i)/dlambda(i),1);
        n4 = n4+1;
    end
             
end


rg = length(G);                         
rh = length(H);                         % São desconsideradas para o cálculo do sigma as restrições de desigualdade canalizadas referentes as variáveis do problema.

% % Wright 1995

% sigmaprimal = 1-(1/(9*sqrt(rg+rh))); 
% sigmadual = 1-(1/(9*sqrt(rg+rh)));


sigmaprimal = 0.995;
sigmadual = sigmaprimal;

alfa_primal_prev = (min(cand_z_prev))*sigmaprimal;
alfa_primal_cor = (min(cand_z_cor))*sigmaprimal;

alfa_dual_prev = (min(cand_lamb_prev))*sigmadual;
alfa_dual_cor = (min(cand_lamb_cor))*sigmadual;



end


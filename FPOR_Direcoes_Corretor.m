function [deta,dx,dz,dlambda] = FPOR_Direcoes_Corretor(JG,JGT,JH,theta,mk,sk,tk,uk,pk,Zinv,L1,A)



b = [mk - pk; tk];

b = sparse(b);


aux = A\b;

n1 = length(aux);

n = length(mk);

dx = aux(1:n);

deta = aux(n+1:n1);


dz = uk - JH*dx;

dlambda = Zinv*(sk-L1*dz);

end
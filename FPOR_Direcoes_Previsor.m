function [detap,dxp,dzp,dlambdap,Dz,A] = FPOR_Direcoes_Previsor(JG,JH,JGT,theta,mk,tk,uk,skprev,pkprev,Zinv,L1)

format long

% A = (JG*thetak*JGT);
% 
% A = sparse(A);
% 
% b = (JG*thetak*(mk-pkprev)-tk); 
% 
% b = sparse(b);
% 
% detap = A\b;


A = [theta  JGT;
           JG     zeros(size(JG,1), size(JG,1))];

A = sparse(A);

b = [mk - pkprev; tk];

b = sparse(b);


aux = A\b;

n1 = length(aux);

n = length(mk);

dxp = aux(1:n);

detap = aux(n+1:n1);

%dxp = thetak*(mk-pkprev-JGT*detap);

dzp = uk - JH*dxp;

dlambdap = Zinv*(skprev-L1*dzp);

Dz = diag(dzp);

end
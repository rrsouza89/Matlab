function [x0,z,eta,lambda] = FPOR_Novo_Ponto(x0,z,eta,lambda,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda);
   
disp('ESTRATÉGIA 1');

% x0p = x0 + alfa_primal_prev*dxp;
% 
% zp = z + alfa_primal_prev*dzp;
% 
% etap = eta + detap;
% 
% lambdap = lambda + alfa_dual_prev*dlambdap;


%% ATUALIZAÇÃO PELO CORRETOR

x0 = x0 + alfa_primal_cor*dx;

z = z + alfa_primal_cor*dz;

eta = eta + deta;

lambda = lambda + alfa_dual_cor*dlambda;

end
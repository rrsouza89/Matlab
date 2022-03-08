function [x0,z,eta,lambda,PC] = FPOR_Novo_Ponto_Nova_Direcao_4(x0,z,eta,lambda,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda,Xi,sigmaprimal,sigmadual,mi,PC,w1,w2);

%disp('ESTRATÉGIA 4');

%% PONTO PROVISÓRIO PREVISOR

x0p = x0 + alfa_primal_prev*dxp;

zp = z + alfa_primal_prev*dzp;

etap = eta + detap;

lambdap = lambda + alfa_dual_prev*dlambdap;

%% PONTO PROVISÓRIO CORRETOR

x0c = x0 + alfa_primal_cor*dx;

zc = z + alfa_primal_cor*dz;

etac = eta + deta;

lambdac = lambda + alfa_dual_cor*dlambda;

%% ATUALIZAÇÃO DO PONTO 


if ((zp'+mi)*lambdap) < Xi*((zc'+mi)*lambdac)
    %disp('Atualização nova direção(Previsor)');
    
    w1; 
    w2; 
    
    [dxn,dzn,detan,dlambdan,alfa_primal_n,alfa_dual_n] = FPOR_Nova_Direcao(w1,w2,z,lambda,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda,sigmaprimal,sigmadual);
    
    x0 = x0 + alfa_primal_n*(dxn);
    
    z = z + alfa_primal_n*(dzn);
    
    eta = eta + (detan);
    
    lambda = lambda + alfa_dual_n*(dlambdan);
    
    PC = [PC;0];    
        
else
    %disp('Atualização nova direção(Corretor)');
    wi = w1;
    w1 = w2;
    w2 = wi;
    
    [dxn,dzn,detan,dlambdan,alfa_primal_n,alfa_dual_n] = FPOR_Nova_Direcao(w1,w2,z,lambda,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda,sigmaprimal,sigmadual);
    
    x0 = x0 + alfa_primal_n*(dxn);
    
    z = z + alfa_primal_n*(dzn);
    
    eta = eta + (detan);
    
    lambda = lambda + alfa_dual_n*(dlambdan);
    
    PC = [PC;1];  
        
end
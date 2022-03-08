function [theta,Beta,FLBM,TQ] = FPOR_Matriz_Theta(x0,JG,JGT,JH,JHT,K,Zinv,L1,Beta,alfa,E_sing,I,nb,Verificacao,FLBM,FLBM2,it,dam1,dx)

    %% MATRIZ THETA

    theta = K + JHT*Zinv*L1*JH; 
    
    theta = sparse(theta);

    %theta_anterior = theta;                                     % Salva Theta da itera��o anterior

    
    %% M�TODO QUASI-NEWTON

%     a = rcond(theta)        %singularidade
% 
%     if a <= E_sing
%     
%         theta = theta_anterior;
%     
%     end    

    %% TESTE QUADR�TICO
    
    % TESTE DE CONVERG�NCIA GLOBAL DIRE��O DX-1
    
    if it == 1
        dx_ant = 0.0005*x0;
    else
        dx_ant = dx;
    end
    

    if Verificacao == 1
        
        [theta,Beta,TQ] = FPOR_Teste_Quadratico(I,Beta,theta,dam1,dx_ant);
        
    
    elseif Verificacao == 2
    
        %% DECOMPOSI��O CHOLESKY
    
        Chol_it = [];
        contp = 0;
        [theta,contp] = FPOR_Teste_LDL(K,JH,JHT,Zinv,L1,I,Beta,contp,nb,theta,dam1); 
        Chol_it = [Chol_it, contp];
    
        elseif Verificacao == 3
    
        %% CORRE��O DE IN�RCIA
    
        [theta] = FPOR_Teste_Inercia(K,JH,JHT,JG,JGT,Zinv,L1,I,Beta,x0,theta);
    
    end   

    %thetak = inv(theta);
    
    
end
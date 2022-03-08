function [theta,Beta,TQ] = FPOR_Teste_Quadratico(I,Beta,theta,dam1,dx_ant);
%% TESTE QUADR�TICO

    uK = dx_ant/(norm(dx_ant));
    
    %theta = theta + Beta*I;
    
    TQ = uK'*theta*uK
    
    contit = 0;

    while TQ <= 0.01
        
        %% DETERMINA��O DO BETA E ATUALIZA��O DA FLBM
        
        Beta = Beta*dam1;
        
        theta = theta + Beta*I;
                          
        contit = contit+1;
        fprintf('Teste quadr�tico %d \n',contit)
        
        TQ = uK'*theta*uK;
        
    end
    
    %fprintf('Teste: %0.4f \n',TQ);
    %fprintf('Itera��es do procedimento: %d \n',contit);
    %fprintf('Beta: %0.4f \n\n',Beta);

    
end
function [theta,contp] = FPOR_Teste_LDL(K,JH,JHT,Zinv,L1,I,Beta,contp,nb,theta,dam1);

[A,p] = chol(theta);

% if nb == 118
%     
%     cont = 0;
%     while p~=0                      %% Condição utilizada para o 118 barras
%         if p~=0
%             Beta = Beta*dam1;
%             theta = K + JHT*Zinv*L1*JH + Beta*I;
%             cont = cont + 1;
%             [A,p] = chol(theta);
%             p;
%             contp = contp + 1;
%         end    
%     end
%     
% else
    
    while p~=0                        %% Teste utilizado para os demais sistemas                      
        Beta = Beta*dam1;
        theta = K + JHT*Zinv*L1*JH + Beta*I;             % MATRIZ THETA
        [A,p] = chol(theta);
        p;
        contp = contp + 1;
      
    end
    
    contp
    
% end

  
end
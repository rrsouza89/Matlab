%
%
function [theta] = FPOR_Teste_Inercia(K,JH,JHT,JG,JGT,Zinv,L1,I,Beta,x0,theta)

   disp('Inércia');
    
    aux = [theta  JGT;
           JG     zeros(size(JG,1), size(JG,1))];
           
    autovalores = eig(aux);
    np = sum(autovalores>0)
    nn = sum(autovalores<0)
    
    while ((np~=size(x0,1)) || (nn~=size(JG,1)))
        disp('Atualização com Beta')
        theta = K + JHT*Zinv*L1*JH + Beta*I;             % MATRIZ THETA

        aux = [theta  JGT;
               JG     zeros(size(JG,1), size(JG,1))];

        autovalores = eig(aux);
        np = sum(autovalores>0);
        nn = sum(autovalores<0); 
        Beta = Beta*1.5;
    end
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                        FUN��O LAGRANGIANA                              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FLBM,FLBM2] = FPOR_Funcao_Lagrangiana(z,delta,mi,eta,lambda,f,fe,fed,G,H,tolz,FLBM,it);

nzh = length(z);

somaln = 0; 
somalambda = 0;

for i=1:nzh
    if z(i)>tolz
        somaln = somaln + delta(i)*log(z(i));         %Bar Log
    else
        somaln = somaln + delta(i)*log(1+z(i)/mi);        % Bar Log Mod
    end
      
    somalambda = somalambda + lambda(i)*(H(i)+z(i));
end

Fbarmod = -mi*somaln

somaeta = eta'*G;

L = somaeta + somalambda



if it == 1                                              % Apenas para a primeira itera��o

    FLBM = f + fe + fed + Fbarmod + L;

    %fprintf('Fun��o Lagrangiana no ponto: %f \n\n',FLBM);
    
    FLBM2 = 0;
    
else
    
    FLBM2 = f + fe + fed + Fbarmod + L;

    fprintf('Fun��o Lagrangiana no ponto: %f \n\n',(f + fe + fed)*10000 + Fbarmod + L);
    
end

end
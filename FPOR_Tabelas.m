function [result] = FPOR_Tabelas_TQ(it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,TQ_it);

format short e;

result = [it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,TQ_it];


fprintf('\n\n')

fprintf('=================================================================================\n')
fprintf('                                RESULTADOS \n')
fprintf('---------------------------------------------------------------------------------\n')
fprintf('it\t\t F. Obj.\t F.L.\t\t Erro\t\t  mi\t\t Beta\t\t T.Quadrático\t\t     \n')
fprintf('-----\t----------\t----------  ----------  ----------  ----------  -------------\n')
fprintf('%4.0f\t %2.4f\t %2.4f\t %1.2e\t %1.2e\t %1.2e\t %1.2e \n',result');
fprintf('---------------------------------------------------------------------------------\n')

end
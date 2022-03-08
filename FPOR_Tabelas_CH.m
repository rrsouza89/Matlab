function [result] = FPOR_Tabelas_CH(it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,Chol_it,Estrategia,PC,sing_it,fp_it)

%format short e;

fprintf('\n\n')

if (Estrategia == 2)|(Estrategia == 4)|(Estrategia == 5)
    
    result = [it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,PC];

    fprintf('=============================================================================\n')
    fprintf('                          RESULTADOS \n')
    fprintf('-----------------------------------------------------------------------------\n')
    fprintf('it\t\t F. Obj.\t F.L.\t\t Erro\t\t  mi\t\t Beta\t\t Corretor\t\t     \n')
    fprintf('-----\t----------\t----------  ----------  ----------  ---------- -------\n')
    fprintf('$%4.0f$ &\t $%2.2f$ &\t $%2.2f$ &\t $%1.2e$ &\t $%1.2e$ &\t $%1.2e$ &\t  $%4.0f$ \\ \t \n',result');
    fprintf('---------------------------------------------------------------------------------\n')

    
else
    
    result = [it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,sing_it,fp_it];

    fprintf('=========================================================================================================================\n')
    fprintf('                                         RESULTADOS \n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------\n')
    fprintf('it\t\t      F. Obj.\t        F.L.\t\t    Erro\t\t      mi\t\t     Beta\t\t   Sing\t\t   Perdas\t\t     \n')
    fprintf('--------\t ----------\t     ----------\t     ----------\t     ----------\t     ------------\t ------------\t ----------- \n')
    fprintf('$%4.0f$ &\t $%2.2f$ &\t $%2.2f$ &\t $%1.2e$ &\t $%1.2e$ &\t $%1.2e$ &\t $%1.2e$ &\t $%1.2f$\\ \t \n',result');
    fprintf('-------------------------------------------------------------------------------------------------------------------------\n')    

end
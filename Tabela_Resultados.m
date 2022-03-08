function [retorno_tabelas] = Tabela_Resultados(nb,nr,ng,ntap,ngt,E,mi,tal,Beta,alfa,Vmin,Vmax,tapmin,tapmax,Estrategia,Verificacao,vM,desv,kE,c,vi,vr,v0,d,Kp,Kr,it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,TQ_it,PC,sing_it,fp_it,b,r,g,Ger,Cd_it,Cp_it,Cr_it,pg_eolica_it,Wr,fWeibull,x0,conj_ramo_tap,f,fe,fed,fp,it,tempo,GF,nodal,omegaR,omegaP)

%% TABELAS DE RESULTADOS

fprintf('============================================================================================================== \n');
fprintf('SISTEMA IEEE: %i barras \n',nb);             %% Mostra na tela
fprintf('-------------------------------------------------------------------------------------------------------------- \n');
fprintf('Número de Barras do sistema: %i \n\n',nb);     %% Mostra na tela
fprintf('Número de Linhas do sistema: %i \n\n',nr);     %% Mostra na tela
fprintf('Número de Geradores do sistema: %i \n\n',ng);  %% Mostra na tela
fprintf('Número de taps: %i\n\n', ntap);

fprintf('============================================================================================================== \n');
fprintf('PARÂMETROS INICIAIS \n');
fprintf('-------------------------------------------------------------------------------------------------------------- \n');
fprintf('E = %4.4f;     mi = %4.4f;     tal = %4.2f;       Beta = %4.4f;    alfa = %4.2f;  \n\n',E,mi,tal,Beta,alfa);
fprintf('%4.2f < V < %4.2f;     %4.2f < tap < %4.2f;      \n\n',Vmin,Vmax,tapmin,tapmax);
fprintf('Estratégia: %d;        Verificação: %d;\n\n',Estrategia,Verificacao);

fprintf('============================================================================================================== \n');
fprintf('PARÂMETROS INICIAIS EÓLICA \n');
fprintf('-------------------------------------------------------------------------------------------------------------- \n');
fprintf('Velocidade Média: %4.2f;   Desvio Padrão: %4.2f;\n\n',vM,desv);
fprintf('Fator de forma k: %4.2f;      Fator de escala c: %4.2f;\n\n',kE,c);
fprintf('v. inicial = %d;   v. nominal = %d;    v. corte = %d;\n\n',vi,vr,v0);
fprintf('d = %d;    Kp = %d;    Kr = %d;',d,Kp,Kr);

if Verificacao == 1
    [result] = FPOR_Tabelas_TQ(it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,TQ_it,Estrategia,PC);

elseif Verificacao == 2
    [result] = FPOR_Tabelas_CH(it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,Chol_it,Estrategia,PC,sing_it,fp_it);
    
end

% VARIÁVEIS

% fprintf('============================================================================================================== \n');
% fprintf('VARIÁVEIS DO SISTEMA \n');
% fprintf('-------------------------------------------------------------------------------------------------------------- \n');
% fprintf('          Magnitude       Ângulo\n');
% for i=1:nb
%    fprintf('Barra %d  %f       %f\n',i,x0(i),x0(i+nb));    
% end
% 
% if ntap ~= 0 
%     fprintf('               Tap\n')
%     for i=1:ntap
%         k = r(conj_ramo_tap(i)).NOI;
%         m = r(conj_ramo_tap(i)).NOF;
%         fprintf('Linha %d - %d    %f\n',k,m,x0(i+2*nb));
%         
%     end    
% end

fprintf('       TGER   Pot       Custo Quad       Custo PV                  Custo Total             $/MWh\n\n')

pot_total = 0;
CQuad = 0;
CPV = 0;

for i=1:ng
    pot_total = pot_total + x0(i+2*nb+ntap);
    if b(Ger(i)).TGER == 1
        fprintf('Barra %d   %d   %4.2f      & %4.2f         &  %4.2f                    & %4.2f                 & %4.2f  \n',Ger(i),b(Ger(i)).TGER,x0(i+2*nb+ntap)*100,g(i).Cter*10000,g(i).Cpv*10000,g(i).Ct*10000,((g(i).Ct)*10000)/(x0(i+2*nb+ntap)*100))
      
        CQuad = CQuad + g(i).Cter*10000;
        CPV = CPV + g(i).Cpv*10000;
        
    else
        fprintf('Barra %d   %d   %4.2f      & %4.2f         &  %4.2f                    & %4.2f                 & %4.2f  \n',Ger(i),b(Ger(i)).TGER,x0(i+2*nb+ntap)*100,g(i).Cter*10000,g(i).Cpv*10000,g(i).Cd*10000+((g(i).Cr*10000)/omegaR)+((g(i).Cp*10000)/omegaP),(g(i).Cd*10000+((g(i).Cr*10000)/omegaR)+((g(i).Cp*10000)/omegaP))/(x0(i+2*nb+ntap)*100))
        
    end
end
fprintf('    \n');
fprintf('Total:       %0.4f     & %0.4f        & %0.4f      & %0.4f          & %0.4f \n',pot_total*100,CQuad,CPV,CQuad+CPV,f*10000+fe*10000+fed*10000);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

fprintf('Perdas do sistema: %f\n\n',fp*100);

%% ENERGIA EÓLICA

if ng ~= ngt
    [result_E] = FPOR_Tabelas_Eolica(it_it,Cd_it,Cp_it,Cr_it,pg_eolica_it,Wr);

    %[retorno] = Grafico_Eolica(kE,c,g,ng,b,Ger,fWeibull);
end

soma_pot = 0;
soma_Cd = 0;
soma_Cr = 0;
soma_Cp = 0;
soma_Ct = 0;


for i=1:ng
    
    if b(Ger(i)).TGER == 2
                %k               c      Pot        Cd         Cr         Cp          Ceolica       Ctotal       Demanda       Perdas       it        tempo             k  c Pot              Cd            Cr            Cp            Ceolica       Ctotal           Demanda                                                     Perdas it       tempo
        fprintf('$%0.2f$  & $%0.2f$ & $%0.3f$  & $%0.3f$  & $%0.3f$  & $%0.3f$   & $%0.3f$     & $%0.3f$    & $%0.3f$     & $%0.3f$    & $%d$  & $%0.2f$   \n',kE,c,b(Ger(i)).PG*100,g(i).Cd*10000,(g(i).Cr*10000)/omegaR,(g(i).Cp*10000)/omegaP,g(i).Cd*10000+((g(i).Cr*10000)/omegaR)+((g(i).Cp*10000)/omegaP),(f+fed+fe)*10000,pot_total*100,fp*100,it,tempo);
    
        soma_pot = soma_pot + b(Ger(i)).PG*100;
        soma_Cd = soma_Cd + g(i).Cd*10000;
        soma_Cr = soma_Cr + (g(i).Cr*10000)/omegaR;
        soma_Cp = soma_Cp + (g(i).Cp*10000)/omegaP;
        soma_Ct = soma_Ct + g(i).Cd*10000+((g(i).Cr*10000)/omegaR)+((g(i).Cp*10000)/omegaP);
    
    end
       
end

                    %k       c           Pot        Cd         Cr         Cp          Ceolica       Ctotal       Demanda       Perdas       it        tempo          k  c Pot     Cd       Cr     Cp      Ceolica       Ctotal           Demanda       Perdas it tempo
        fprintf('$%0.2f$ & $%0.2f$ & $%0.3f$  & $%0.3f$  & $%0.3f$  & $%0.3f$   & $%0.3f$     & $%0.3f$    & $%0.3f$     & $%0.3f$    & $%d$  & $%0.2f$   \n',kE,c,soma_pot,soma_Cd,soma_Cr,soma_Cp,soma_Ct,(f+fed+fe)*10000,pot_total*100,fp*100,it,tempo);
 
        fprintf('tempo:   %0.2f \n\n ',tempo)



        
        
        
%% CUSTO INCREMENTAL

Pot = [];
Quad = [];
PV = [];
Lin = [];
Res = [];
Pen = [];
incremental = [];


for i=1:ng
   
    Pot = [Pot; (b(Ger(i)).PG)*100];
    Quad = [Quad; (g(i).Cter)*10000];
    PV = [PV; (g(i).Cpv)*10000];
    Lin = [Lin; (g(i).Cd)*10000];
    Res = [Res; ((g(i).Cr)*10000)/omegaR];
    Pen = [Pen; ((g(i).Cp)*10000)/omegaP];
    incremental = [incremental; GF(2*nb+ntap+i)*100];     
    
end        
        
        
teste = [Ger'   Pot      Quad  PV     Lin    Res    Pen     incremental     nodal];       
        
fprintf('%d  &  %0.2f  & %0.2f  & %0.2f &  %0.2f & %0.2f  & %0.2f   & %0.4f  &   %0.4f  \n', teste')




if nb == 39

            %p30   p31   p32   p33   p34   p35   p36   p37   p38   p39   Demanda Quad   PV     Term    Linear  Reserva  Penalidade   Eólica  Total    It Tempo       %p30               p31                p32                p33                p34                p35                p36                p37                p38                p39                           Demanda        Quad    PV   Term       Linear     Reserva  Penalidade   Eólica         Total    It Tempo o
    fprintf('%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f   %0.2f  %0.2f  %0.2f   %0.2f   %0.2f    %0.2f        %0.2f    %0.2f   %d  %0.2f \n', (b(Ger(1)).PG)*100, (b(Ger(2)).PG)*100, (b(Ger(3)).PG)*100, (b(Ger(4)).PG)*100, (b(Ger(5)).PG)*100, (b(Ger(6)).PG)*100, (b(Ger(7)).PG)*100, (b(Ger(8)).PG)*100, (b(Ger(9)).PG)*100, (b(Ger(10)).PG)*100, pot_total*100, CQuad,  CPV, CQuad+CPV, soma_Cd,   soma_Cr,  soma_Cp,  soma_Ct, (f+fed+fe)*10000, it, tempo)



end
% 
% 
% % % RESTRIÇÕES
% % fprintf('============================================================================================================== \n');
% % fprintf('RESTRIÇÕES \n');
% % fprintf('-------------------------------------------------------------------------------------------------------------- \n');
% % 
% % fprintf('Restrições de Igualdade: \n')
% % [G]
% % 
% % fprintf('Restrições de Desigualdade: \n')
% % 
% % [U1 H U2]
% % 
% % fprintf('Variáveis canalizadas: \n')
% % 
% % [U3 x0 U4]
% % 
% % % MULTIPLICADORES DE LAGRANGE
% % 
% % fprintf('==============================================================================================================\n')
% % fprintf('MULTIPLICADORES DE LAGRANGE \n');
% % fprintf('--------------------------------------------------------------------------------------------------------------\n')
% % 
% % fprintf('Multiplicador Rest. Igualdade: \n')
% % lamb0
% % 
% % fprintf('Multiplicador Rest. Desigualdade: \n')
% % [lamb1 lamb2]
% % 
% % fprintf('Multiplicador Rest. Desigualdade Variáveis: \n')
% % [lamb3 lamb4]
% %     
%     
% GRÁFICOS

% V = zeros(nb,1);
% T = zeros(nb,1);
% TAP = zeros(ntap,1);
% Gr = zeros(ng,1);
 %Qg = zeros(ng,1);
% 
% 
% for i=1:nb
%    V(i,1) = x0(i);
% end
% 
% for i=1:nb
%    T(i,1) = x0(i+nb);
% end
% 
% for i=1:ntap
%     TAP(i,1) = x0(i+2*nb);
% end
% 
% for i=1:ng
%     %Gr(i,1) = x0(i+2*nb+ntap);
%     Qg(i,1) = b(Ger(i)).QG;
% end
% 
% 
% 
% 
% liminf = Vmin*ones(size(V));
% limsup = Vmax*ones(size(V));
% figure(1)
% hold on
% title('V');
% plot(V.','*r');
% plot(liminf.','-');
% plot(limsup.','-');
% hold off
% 
% liminf = -0.8*ones(size(T));
% limsup = 0.8*ones(size(T));
% figure(2)
% hold on
% title('Ângulo');
% plot(T.','*r');
% plot(liminf.','-');
% plot(limsup.','-');
% hold off
%     
% liminf = 0.8*ones(size(TAP));
% limsup = 1.2*ones(size(TAP));
% figure(3)
% hold on
% title('Tap');
% plot(TAP.','*r');
% plot(liminf.','-');
% plot(limsup.','-');
% hold off
% 
% liminf = 0*ones(size(Gr));
% limsup = 8*ones(size(Gr));
% figure(4)
% hold on
% title('Potência Ativa Gerada');
% plot(Gr.','*r');
% plot(liminf.','-');
% plot(limsup.','-');
% hold off
% 
% liminf = -3*ones(size(Qg));
% limsup = 3*ones(size(Qg));
% figure(5)
% hold on
% title('Potência Reativa Injetada');
% plot(Qg.','*r');
% plot(liminf.','-');
% plot(limsup.','-');
% hold off
% 
% Qg;

retorno_tabelas = 1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                        FPO - INICIALIZAÇÃO                  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format long
clear all
close all
clc;

fprintf('FLUXO DE POTÊNCIA ÓTIMO ATIVO/REATIVO\n\n');

%% SISTEMAS ELÉTRICOS

%% PARÂMETROS INICIAIS                                                             
                
%% CASE 30 BUS

% caso = case30_E_Mishra;
% 
% E = 10^-2;                      % Precisão Epsilon para o critério de parada                        
% 
% mi = 0.05;                       % Parâmetro de barreira                           
% 
% tal = 0.15;                     % Parâmetro de atualização do mi                   
% 
% Beta = 0.1;                    % Parâmetro de atualização da matriz Theta          
% 
% alfa = tal;                    % Atualização de Beta                                                   
% 
% Xi = 0.95;                      % Parâmetro de escolha dos passos previsor e corretor               
% 
% Vmin = 0.9;     Vmax = 1.1;   % Limites mínimo e máximo das magnitudes de tensão                    
% 
% tapmin = 0.95;  tapmax = 1.05;  % Limites mínimo e máximo dos taps dos transformadores                   
% 
% Estrategia = 4;                 % Escolha da estratégia: 1 - Corretor // 2 - Diego // 3 - ColomboGondzio // 4 - Convexa // 5 - Linear
% 
% w1 = 0.95;   w2 = 0.05;
% 
% Verificacao = 1;                % Escolha do procedimento de verificação da matriz theta:   1 - Teste Quadrático // 2 - LDL // 3 - Inércia
% 
% tolz = 10^20;                   % 10^20 = Barreira Log Mod  //  10^-20 = Barreira Log
% 
% E_sing = 10^-20;                % Parâmetro para o teste de singularidade da matriz theta
% 
% %% PARÂMETROS PARTE EÓLICA 
% 
% %% Parâmetros manuais     
% 
% reg = 5;
%                                                          %        Dez-Fev         Mar-Mai          Jun-Ago              Set-Nov          Anual
% kE = 2;                         % Fator de forma               2.39            1.88              2.81                3.07             2.24
% 
% c = 10;                          % Fator de escala              5.27            3.93              6.30                7.48             5.78
% 
% fWeibull = 2;                   % Escolha da função de Weibull:        1 - Hetzer // 2 - Mishra  // 3 - (k/c)  //
% 
% omegaR = 1;                     %Ponderação do Custo de Reserva
% 
% omegaP = 1;                     %Ponderação do Custo de Penalidade
% 
% L = -10;                            %Precisão apenas para a parte contínua
% 
% %L = 10^-3;                          % Precisão da parte contínua e discreta da função de Weibull
% 
% vM = 0.00;                          % Velocidade Média
% 
% desv = 0.00;                        % Desvio Padrão
% 
% veloc = 1;                          % vetor das velocidades do vento da região
% 
% if reg == 1
%     kE = 2.39; c = 5.27;
% elseif reg == 2
%     kE = 1.88; c = 3.93;
% elseif reg == 3
%     kE = 2.81; c = 6.30;
% elseif reg == 4
%     kE = 3.07; c = 7.48;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CASE 57 BUS

% caso = case57_E_Mishra;
% 
% E = 10^-2;                      % Precisão Epsilon para o critério de parada                        
% 
% mi = 0.05;                       % Parâmetro de barreira                           
% 
% tal = 0.15;                     % Parâmetro de atualização do mi                   
% 
% Beta = 0.1;                    % Parâmetro de atualização da matriz Theta          
% 
% alfa = tal;                    % Atualização de Beta                                                   
% 
% Xi = 0.95;                      % Parâmetro de escolha dos passos previsor e corretor               
% 
% Vmin = 0.9;     Vmax = 1.1;   % Limites mínimo e máximo das magnitudes de tensão                    
% 
% tapmin = 0.95;  tapmax = 1.05;  % Limites mínimo e máximo dos taps dos transformadores                   
% 
% Estrategia = 4;                 % Escolha da estratégia: 1 - Corretor // 2 - Diego // 3 - ColomboGondzio // 4 - Convexa // 5 - Linear
% 
% w1 = 0.95;   w2 = 0.05;
% 
% Verificacao = 1;                % Escolha do procedimento de verificação da matriz theta:   1 - Teste Quadrático // 2 - LDL // 3 - Inércia
% 
% tolz = 10^20;                   % 10^20 = Barreira Log Mod  //  10^-20 = Barreira Log
% 
% E_sing = 10^-20;                % Parâmetro para o teste de singularidade da matriz theta
% 
% %% PARÂMETROS PARTE EÓLICA 
% 
% %% Parâmetros manuais     
% 
% reg = 5;
%                                                          %        Dez-Fev         Mar-Mai          Jun-Ago              Set-Nov          Anual
% kE = 2;                         % Fator de forma               2.39            1.88              2.81                3.07             2.24
% 
% c = 10;                          % Fator de escala              5.27            3.93              6.30                7.48             5.78
% 
% fWeibull = 2;                   % Escolha da função de Weibull:        1 - Hetzer // 2 - Mishra  // 3 - (k/c)  //
% 
% omegaR = 1;                     %Ponderação do Custo de Reserva
% 
% omegaP = 1;                     %Ponderação do Custo de Penalidade
% 
% L = -10;                            %Precisão apenas para a parte contínua
% 
% %L = 10^-3;                          % Precisão da parte contínua e discreta da função de Weibull
% 
% vM = 0.00;                          % Velocidade Média
% 
% desv = 0.00;                        % Desvio Padrão
% 
% veloc = 1;                          % vetor das velocidades do vento da região
% 
% if reg == 1
%     kE = 2.39; c = 5.27;
% elseif reg == 2
%     kE = 1.88; c = 3.93;
% elseif reg == 3
%     kE = 2.81; c = 6.30;
% elseif reg == 4
%     kE = 3.07; c = 7.48;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CASE 118 BUS

% caso = case118_E_Ellatar;
% 
% E = 10^-2;                      % Precisão Epsilon para o critério de parada                        
% 
% mi = 0.5;                       % Parâmetro de barreira                           
% 
% tal = 0.25;                     % Parâmetro de atualização do mi                   
% 
% Beta = 0.1;                    % Parâmetro de atualização da matriz Theta          
% 
% alfa = tal;                    % Atualização de Beta                                                   
% 
% Xi = 0.95;                      % Parâmetro de escolha dos passos previsor e corretor               
% 
% Vmin = 0.95;     Vmax = 1.05;   % Limites mínimo e máximo das magnitudes de tensão                    
% 
% tapmin = 0.9;  tapmax = 1.1;  % Limites mínimo e máximo dos taps dos transformadores                   
% 
% Estrategia = 4;                 % Escolha da estratégia: 1 - Corretor // 2 - Diego // 3 - ColomboGondzio // 4 - Convexa // 5 - Linear
% 
% w1 = 0.95;   w2 = 0.05;
% 
% Verificacao = 1;                % Escolha do procedimento de verificação da matriz theta:   1 - Teste Quadrático // 2 - LDL // 3 - Inércia
% 
% tolz = 10^20;                   % 10^20 = Barreira Log Mod  //  10^-20 = Barreira Log
% 
% E_sing = 10^-20;                % Parâmetro para o teste de singularidade da matriz theta
% 
% %% PARÂMETROS PARTE EÓLICA 
% 
% %% Parâmetros manuais     
% 
% reg = 4;
%                                                          %        Dez-Fev         Mar-Mai          Jun-Ago              Set-Nov          Anual
% kE = 2;                         % Fator de forma               2.39            1.88              2.81                3.07             2.24
% 
% c = 10;                          % Fator de escala              5.27            3.93              6.30                7.48             5.78
% 
% fWeibull = 2;                   % Escolha da função de Weibull:        1 - Hetzer // 2 - Mishra  // 3 - (k/c)  //
% 
% omegaR = 10;                     %Ponderação do Custo de Reserva
% 
% omegaP = 10;                     %Ponderação do Custo de Penalidade
% 
% L = -10;                            %Precisão apenas para a parte contínua
% 
% %L = 10^-3;                          % Precisão da parte contínua e discreta da função de Weibull
% 
% vM = 0.00;                          % Velocidade Média
% 
% desv = 0.00;                        % Desvio Padrão
% 
% veloc = 1;                          % vetor das velocidades do vento da região
% 
% if reg == 1
%     kE = 2.39; c = 5.27;
% elseif reg == 2
%     kE = 1.88; c = 3.93;
% elseif reg == 3
%     kE = 2.81; c = 6.30;
% elseif reg == 4
%     kE = 3.07; c = 7.48;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CASE 300 BUS

caso = case300_E_pinheiro;

E = 0.25;                      % Precisão Epsilon para o critério de parada                        

mi = 0.5;                       % Parâmetro de barreira                           

tal = 0.25;                     % Parâmetro de atualização do mi                   

Beta = 0.1;                    % Parâmetro de atualização da matriz Theta          

alfa = tal;                    % Atualização de Beta                                                   

Xi = 0.95;                      % Parâmetro de escolha dos passos previsor e corretor               

Vmin = 0.9;     Vmax = 1.1;   % Limites mínimo e máximo das magnitudes de tensão                    

tapmin = 0.9;  tapmax = 1.1;  % Limites mínimo e máximo dos taps dos transformadores                   

Estrategia = 4;                 % Escolha da estratégia: 1 - Corretor // 2 - Diego // 3 - ColomboGondzio // 4 - Convexa // 5 - Linear

w1 = 0.95;   w2 = 0.05;

Verificacao = 1;                % Escolha do procedimento de verificação da matriz theta:   1 - Teste Quadrático // 2 - LDL // 3 - Inércia

tolz = 10^20;                   % 10^20 = Barreira Log Mod  //  10^-20 = Barreira Log

E_sing = 10^-20;                % Parâmetro para o teste de singularidade da matriz theta

%% PARÂMETROS PARTE EÓLICA 

%% Parâmetros manuais     

reg = 5;
                                                         %        Dez-Fev         Mar-Mai          Jun-Ago              Set-Nov          Anual
kE = 2;                         % Fator de forma               2.39            1.88              2.81                3.07             2.24

c = 10;                          % Fator de escala              5.27            3.93              6.30                7.48             5.78

fWeibull = 2;                   % Escolha da função de Weibull:        1 - Hetzer // 2 - Mishra  // 3 - (k/c)  //

omegaR = 10;                     %Ponderação do Custo de Reserva

omegaP = 10;                     %Ponderação do Custo de Penalidade

L = -10;                            %Precisão apenas para a parte contínua

%L = 10^-3;                          % Precisão da parte contínua e discreta da função de Weibull

vM = 0.00;                          % Velocidade Média

desv = 0.00;                        % Desvio Padrão

veloc = 1;                          % vetor das velocidades do vento da região

if reg == 1
    kE = 2.39; c = 5.27;
elseif reg == 2
    kE = 1.88; c = 3.93;
elseif reg == 3
    kE = 2.81; c = 6.30;
elseif reg == 4
    kE = 3.07; c = 7.48;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                                      
%% ESTRUTURAS DE BARRA E RAMO    

% Estrutura para armazenamento dos dados de cada barra
barra = struct('NUM', 0, 'TIPO', 0, 'TGER', 0, 'PG', 0, 'QG', 0, 'Fluxo', 0, 'PD', 0, 'QD', 0, 'BS', 0, 'v', 1.0, 'teta', 0.0, 'VIZ', [], 'TVIZ', [], 'RAMO', [], 'QMIN', 0, 'QMAX', 0, 'm', 0);      

% Estrutura para armazenamento dos ramos das barras
ramo = struct('NOI', 0, 'NOF', 0, 'R', 0, 'X', 0, 'B', 0, 'TAP', 0, 'GKM', 0, 'BKM', 0);     

% Estrutura para armazenamento dos dados dos geradores
geradores = struct('GER', 0, 'A', 0, 'B', 0, 'C', 0, 'E', 0, 'F', 0, 'Vi', 0, 'Vr', 0, 'V0', 0, 'Wr', 0, 'd', 0, 'Kp', 0, 'Kr', 0, 'Cd', 0, 'Cp', 0, 'Cr', 0, 'Ct', 0, 'Cter', 0, 'Cpv', 0, 'CIL', 0, 'CIP', 0, 'CIR', 0);

nb = length(caso.bus(:,1));                             %% Calcula o número de barras do sistema
nr = length(caso.branch(:,1));                          %% Calcula o número de ramos do sistema
ng = length(caso.gen(:,1));                             %% Calcula o número de geradores
fprintf('Sistema IEEE: %i barras \n\n',nb);             %% Mostra na tela
fprintf('Número de Barras do sistema: %i \n\n',nb);     %% Mostra na tela
fprintf('Número de Linhas do sistema: %i \n\n',nr);     %% Mostra na tela
fprintf('Número de Geradores do sistema: %i \n\n',ng);  %% Mostra na tela
b(nb) = barra;                                          %% Cria a estrutura para todas as barras do sistema
r(nr) = ramo;                                           %% Cria a estrutura para todos os ramos do sistema
g(ng) = geradores;                                      %% Cria a estrutura para todos os geradores do sistema

%% DADOS DE BARRA
for i=1:nb                                              
    b(i).NUM  = caso.bus(i,1);                          %% Número "real" da barra
    b(i).TIPO = caso.bus(i,2);                          %% Tipo da barra // Sl // PV // PQ
    b(i).TGER = 0;                                      %% Tipo de Gerador // 0 - Não possui // 1 - Gerador Térmico // 2 - Gerador Eólico
    b(i).PG   = 0;                                      %% Potência Ativa Gerada
    b(i).QG   = 0;                                      %% Potência Reativa Gerada
    b(i).Fluxo = 0;                                     %% Fluxo Pkm da barra
    b(i).PD   = caso.bus(i,3)/100;                      %% Demanda potência ativa da barra (em pu.) (Pc)
    b(i).QD   = caso.bus(i,4)/100;                      %% Demanda potência reativa da barra (em pu.) (Qc)
    b(i).BS   = caso.bus(i,6);                          %% B shunt de barra
    b(i).v    = caso.bus(i,8);                          %% Magnitude de tensão (em pu)
    b(i).teta = (caso.bus(i,9)*pi)/180;                 %% Ângulo da tensão (em radianos)    
    b(i).QMIN = 0;                                      %% Limite mínimo de potência reativa
    b(i).QMAX = 0;                                      %% Limite máximo de potência reativa
    b(i).m    = 0;                                      %% Variável para os geradores termelétricos relacionada com a restrição da parte senoidal do ponto de válvula
end

%% RELAÇÃO ENTRE NÚMEROS DE BARRAS REAIS E SEQUÊNCIAIS

for i=1:nb
    a = b(i).NUM;
    NIN(a)= i;
end
sparse(NIN);

%% DADOS DE RAMO

for i=1:nr
    r(i).NOI = NIN(caso.branch(i,1));                   %% Nó inicial "real" do ramo
    r(i).NOF = NIN(caso.branch(i,2));                   %% Nó final "real" do ramo
    r(i).R = caso.branch(i,3);                          %% Resistência série rkm
    r(i).X = caso.branch(i,4);                          %% Reatância série xkm
    r(i).B = caso.branch(i,5)/2;                        %% Susceptância shunt bkmsh (bkmsh de linha)
    r(i).TAP = caso.branch(i,9);                        %% Tap do transformador em-fase
end

%% CÁLCULO DA CONDUTÂNCIA (GKM) E SUSCEPTÂNCIA (BKM)

for i=1:nr
    k = r(i).NOI;
    m = r(i).NOF;
    rkm = r(i).R;
    xkm = r(i).X;
    
    r(i).GKM = rkm/(rkm^2 + xkm^2);
    r(i).BKM = -xkm/(rkm^2 + xkm^2);
end


%% DADOS PG, QG, QMIN E QMAX NA ESTRUTURA DE BARRAS

 for i=1:ng
     b(NIN(caso.gen(i,1))).PG = caso.gen(i,2)/100;
     b(NIN(caso.gen(i,1))).QG = caso.gen(i,3)/100;
 end

%% BARRAS VIZINHAS

for i=1:nr
    k = r(i).NOI;
    m = r(i).NOF;
    b(k).RAMO = [b(k).RAMO, i];                         
    b(m).RAMO = [b(m).RAMO, i];
    b(k).VIZ = [b(k).VIZ, m];
    b(k).TVIZ = [b(k).TVIZ, 0];                         %% Atribui valor 0 para a barra que é nó final em relação a barra escolhida
    b(m).VIZ = [b(m).VIZ, k];
    b(m).TVIZ = [b(m).TVIZ, 1];                         %% Atribui valor 1 para a barra que é nó inicial em relação a barra escolhida
end

%% CONJUNTOS DE BARRAS

Sl = [];                        %% 3 - Barra Slack
PV = [];                        %% 2 - Barra de Geração PV
PQ = [];                        %% 1 - Barra de Carga PQ

for i=1:nb
    if caso.bus(i,2) == 3                                  
        Sl = [Sl, NIN(caso.bus(i,1))];        
    else
        if caso.bus(i,2) == 2                              
            PV = [PV, NIN(caso.bus(i,1))];            
        else
            PQ = [PQ, NIN(caso.bus(i,1))];                      
        end
    end
end

nSl = length(Sl);
nPV = length(PV);
nPQ = length(PQ);

%% CONJUNTO DOS GERADORES E DADOS DOS GERADORES

Ger = [];

for i=1:ng
    Ger = [Ger, NIN(caso.gen(i,1))];
    g(i).GER = caso.gen(i,1);
    g(i).A = caso.gencost(i,5)*10000;
    g(i).B = caso.gencost(i,6)*100;
    g(i).C = caso.gencost(i,7);
    
    if nb > 300
        g(i).E = 0;
        g(i).F = 0;
    else
        g(i).E = caso.gencost(i,8);
        g(i).F = caso.gencost(i,9)*100;
    end
    
    g(i).Cd = 0;
    g(i).Cp = 0;
    g(i).Cr = 0;
    g(i).Cter = 0;
    g(i).Cpv = 0;
end

% Caso com Eólica
for i=1:ng
    if nb > 300
        b(Ger(i)).TGER = 1;
        g(i).Vi = 0;                           % Dados dos geradores eólicos
        g(i).Vr = 0;
        g(i).V0 = 0;
        g(i).Wr = 0;
        g(i).d = 0;
        g(i).Kp = 0;
        g(i).Kr = 0;
    else
        b(Ger(i)).TGER = caso.gen(i,22);                    % Tipo de gerador
        g(i).Vi = caso.gen(i,23);                           % Dados dos geradores eólicos
        g(i).Vr = caso.gen(i,24);
        g(i).V0 = caso.gen(i,25);
        g(i).Wr = caso.gen(i,26)/100;
        g(i).d = caso.gen(i,27);
        g(i).Kp = caso.gen(i,28);
        g(i).Kr = caso.gen(i,29);
    end
end

GerT = [];                                      % Conjunto dos geradores Termelétricos

for i=1:ng
   
    if b(Ger(i)).TGER == 1
        GerT = [GerT, Ger(i)];
    end
    
end

ngt = length(GerT);                              % Número de geradores Termelétricos



%% NÚMERO DE TAPS E RAMOS COM OS TAPS

ntap = 0;
conj_ramo_tap = [];

for i=1:nr
    if caso.branch(i,9) ~= 0
        ntap = ntap+1;
        conj_ramo_tap = [conj_ramo_tap,i];                  %% Fornece os ramos que possuem taps 
     else
         r(i).TAP = 1;
    end
end

fprintf('Número de taps: %i\n\n', ntap);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% VETORES LIMITANTES

for i=1:ng
    U1(i,1) = caso.gen(i,5)/100;                % U1 - Limite mínimo de injeção líquida de potência reativa
    U2(i,1) = caso.gen(i,4)/100;                % U2 - Limite máximo de injeção líquida de potência reativa
end

                                                                                                                                %PONTO DE VÁLVULA
for i=1:ngt
    U1(ng+i,1) = 0;                % U1 - Limite mínimo da restrição ponto de válvula
    U2(ng+i,1) = 0;                % U2 - Limite máximo da restrição ponto de válvula
end

for i=1:nb
   U3(i,1) = Vmin;                               % U3 - Limite mínimo de x
   U3(i+nb,1) = -pi;
   U4(i,1) = Vmax;                               % U4 - Limite máximo de x
   U4(i+nb,1) = pi;    
end

for i=1:ntap
    U3(2*nb+i,1) = tapmin;
    U4(2*nb+i,1) = tapmax;
end

for i=1:ng
    U3(2*nb+ntap+i,1) = caso.gen(i,10)/100;     % Limite mínimo de injeção líquida de potência ativa
    U4(2*nb+ntap+i,1) = caso.gen(i,9)/100;      % Limite máximo de injeção líquida de potência ativa
end
                                                                                                                        %PONTO DE VÁLVULA
for i=1:ngt
    U3(2*nb+ntap+ng+i,1) = 0;     % U3 - Limite mínimo da restrição ponto de válvula
    U4(2*nb+ntap+ng+i,1) = 0;      % U4 - Limite máximo da restrição ponto de válvula
end

U3(nb+Sl,1) = -0.1;               % restringe o ângulo da slack a 0
U4(nb+Sl,1) = 0.1;

if nb == 3
    U3(3,1) = 1;
    U4(3,1) = 1;
    U4(2,1) = 1.2;
end

PkmMax = 0.61;
PkmMin = - PkmMax;

%% PONTO INICIAL

x0 = zeros(2*nb+ntap+ng+ngt,1);                                       %PONTO DE VÁLVULA
%x0 = zeros(2*nb+ntap+ng,1);

xv = zeros(nb,1);
xt = zeros(nb,1);
xtap = zeros(ntap,1);
xg = zeros(ng,1);
xpv = zeros(ngt,1);                                                    %PONTO DE VÁLVULA

%V
for i=1:nb
   xv(i,1) = b(i).v;
   %xv(i,1) = 1;
end
%teta
for i=1:nb
   xt(i,1) = b(i).teta;
   %xt(i,1) = 0;
end
%tap
for i=1:ntap
    xtap(i,1) = r(conj_ramo_tap(i)).TAP;
end
%Geradores
for i=1:ng
    xg(i,1) = caso.gen(i,2)/100;
end
                                                                            %PONTO DE VÁLVULA
%Variável PV
for i=1:ngt
    xpv(i,1) = b(GerT(i)).m;
end



%% Inicialização Potência ativa 300 barras e parâmetros E e F 39 e 300 barras

if nb == 57
    for i=1:ng
        PGMIN = U3(2*nb+ntap+i,1);
        PGMAX = U4(2*nb+ntap+i,1);
        A = g(i).A;
        B = g(i).B;
        C = g(i).C;
        g(i).E = 0.05*(((A*PGMAX^2 + B*PGMAX + C)+(A*PGMIN^2 + B*PGMIN + C))/2);
        g(i).F = (4*pi)/(PGMAX - PGMIN);
    end
end

if nb >= 300
    for i=1:ng
        PGMIN = U3(2*nb+ntap+i,1);
        PGMAX = U4(2*nb+ntap+i,1);
        xg(i,1) = 0.45*(PGMAX - PGMIN) + PGMIN;              % 0.45
        A = g(i).A;
        B = g(i).B;
        C = g(i).C;
        g(i).E = 0.05*(((A*PGMAX^2 + B*PGMAX + C)+(A*PGMIN^2 + B*PGMIN + C))/2);
        g(i).F = (4*pi)/(PGMAX - PGMIN);
    end
end

x0 = [xv; xt; xtap; xg; xpv];                                               %PONTO DE VÁLVULA
%x0 = [xv; xt; xtap; xg];

%% DIMENSÃO DO PROBLEMA

RI  = nSl + nPV + 2*(nPQ);                                      % Restrições de Igualdade
RDC = nSl + nPV + nb + ntap + ng;                               % Restrições de Desigualdade Canalizadas
NV  = 2*nb + ntap + ng;                                         % Número de Variáveis contínuas
%fprintf('Restrições de Igualdade: %i \n',RI);
%fprintf('Restrições de Desigualdade Canalizadas: %i \n',RDC);
%fprintf('Número de Variáveis contínuas: %i \n',NV);


%% BANCO DE DADOS

%[dbarra,dlinha] = FPOR_Banco(b,r);
                                                                                                        %PONTO DE VÁLVULA
I = eye(2*nb+ntap+ng+ngt);             % Matriz Identidade
%I = eye(2*nb+ntap+ng); 

E;
mi;
Beta;
PC = [];
TQ_it = []; 
z = 0;
eta = 0;
lambda = 0;
delta = 0;
FLBM = 0;
dx = 0;

%% INICIALIZAÇÃO

it = 0;
Erro = 1;

it_it = [];     Beta_it = [];       mi_it = [];     sing_it = [];       Erro_it = [];   alfa_primal_cor_it = [];    alfa_dual_cor_it = [];

f_it = [];     fp_it = [];      FL_it = [];       pg_eolica_it = [];      Cd_it = [];     Cp_it = [];     Cr_it = [];

x0_it = [];       eta_it = [];     lambda_it = [];       z_it = [];     TQ_it = [];

dam1 = ((1+sqrt(((sqrt(5)-1)^2)*(alfa^2)+1))/2);            % Valor utilizado para atualizar Beta qdo Dif > 0.75 

to = cputime;                   % Marca o tempo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PROCESSO ITERATIVO DO MÉTODO

while (Erro>=E) & (it<100)
    
    %disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    it = it+1;      it_it = [it_it; it];
    
    fprintf('ITERAÇÃO: %i\n\n',it)
    
    %% GRADIENTES E JACOBIANAS
    
    [f,fe,fed,fp,GF,G,H,JG,JH,JGT,JHT,seq,seqh,seqht,eta,lambda,z,Z,Zinv,L1,delta,g,vi,vr,v0,d,Cd,Cp,Cr,Kp,Kr,pg_eolica,Wr] = FPOR_Gradiente_Jacobiana(mi,b,r,g,Ger,NIN,nb,nr,ntap,ng,ngt,conj_ramo_tap,PQ,PV,Sl,x0,z,eta,lambda,delta,U1,U2,U3,U4,tolz,c,kE,L,fWeibull,it,omegaR,omegaP,PkmMin,PkmMax);

    f_it = [f_it; (f+fed+fe)*10000];     fp_it = [fp_it; fp*100];       pg_eolica_it = [pg_eolica_it; pg_eolica*100];      Cd_it = [Cd_it; Cd*10000];     Cp_it = [Cp_it; (Cp*10000)/omegaP];     Cr_it = [Cr_it; (Cr*10000)/omegaR];

    %% HESSIANAS 
    nodal = [];
    [HessG, HessH, HessF, K, nodal] = FPOR_Hessianas(seq,seqh,seqht,b,r,g,NIN,nb,nr,ntap,ng,ngt,conj_ramo_tap,PQ,PV,Sl,eta,lambda,x0,c,Ger,L,kE,Kp,Kr,fWeibull,U3,nodal,omegaR,omegaP);

    %% FUNÇÃO LAGRANGIANA BARREIRA LOG MODIFICADA

    [FLBM,FLBM2] = FPOR_Funcao_Lagrangiana(z,delta,mi,eta,lambda,f,fe,fed,G,H,tolz,FLBM,it);

    FL_it = [FL_it; FLBM2*10000];
        
    %% Atualização beta
    
    [Beta,FLBM,Dif] = FPOR_Beta(FLBM,FLBM2,Beta,it,dam1); 
    
    %% RESÍDUOS PREVISOR

    [mk,skprev,tk,uk,pkprev] = FPOR_Residuos_Previsor(G,H,GF,JGT,JHT,eta,lambda,Z,Zinv,L1,mi,delta,z);
      
    %% MATRIZ HESSIANA THETA
        
    [theta,Beta,FLBM,TQ] = FPOR_Matriz_Theta(x0,JG,JGT,JH,JHT,K,Zinv,L1,Beta,alfa,E_sing,I,nb,Verificacao,FLBM,FLBM2,it,dam1,dx);        Beta_it = [Beta_it;Beta];      TQ_it = [TQ_it;TQ];

    %% DIREÇÕES PREVISOR

    [detap,dxp,dzp,dlambdap,Dz,A] = FPOR_Direcoes_Previsor(JG,JH,JGT,theta,mk,tk,uk,skprev,pkprev,Zinv,L1);

    %% RESÍDUOS CORRETOR

    [sk,pk] = FPOR_Residuos_Corretor(Z,Zinv,lambda,mi,delta,Dz,dlambdap,JHT,L1,uk);

    %% DIREÇÕES CORRETOR

    [deta,dx,dz,dlambda] = FPOR_Direcoes_Corretor(JG,JGT,JH,theta,mk,sk,tk,uk,pk,Zinv,L1,A);

    %% CÁLCULO DOS PASSOS PRIMAL E DUAL

    [alfa_primal_prev,alfa_primal_cor,alfa_dual_prev,alfa_dual_cor,sigmaprimal,sigmadual] = FPOR_Calculo_Passos(z,lambda,dzp,dlambdap,dz,dlambda,G,H,mi);

    alfa_primal_cor_it = [alfa_primal_cor_it; alfa_primal_cor];     alfa_dual_cor_it = [alfa_dual_cor_it; alfa_dual_cor];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% ESCOLHA DA ESTRATÉGIA DE DETERMINAÇÃO DO NOVO PONTO

    if Estrategia == 1
        %% ESTRATÉGIA 1 - NOVO PONTO CORRETOR

        [x0,z,eta,lambda] = FPOR_Novo_Ponto(x0,z,eta,lambda,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda);

        elseif Estrategia == 2
            %% ESTRATÉGIA 2 - ATUALIZAÇÃO DO DIEGO

            [x0,z1,z2,z3,z4,lamb0,lamb1,lamb2,lamb3,lamb4,zp,zc,lambp,lambc,PC] = FPOR_Novo_Ponto_Diego(x0,z1,z2,z3,z4,lamb0,lamb1,lamb2,lamb3,lamb4,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dz1p,dz2p,dz3p,dz4p,dlamb0p,dlamb1p,dlamb2p,dlamb3p,dlamb4p,dxc,dz1c,dz2c,dz3c,dz4c,dlamb0c,dlamb1c,dlamb2c,dlamb3c,dlamb4c,Xi,mi,PC);

            elseif Estrategia == 3
                %% ESTRATÉGIA 3 - COLOMBO-GONDZIO
    
                [x0,z1,z2,z3,z4,lamb0,lamb1,lamb2,lamb3,lamb4] = FPOR_Novo_Ponto_CG(x0,z1,z2,z3,z4,lamb0,lamb1,lamb2,lamb3,lamb4,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dz1p,dz2p,dz3p,dz4p,dlamb0p,dlamb1p,dlamb2p,dlamb3p,dlamb4p,dxc,dz1c,dz2c,dz3c,dz4c,dlamb0c,dlamb1c,dlamb2c,dlamb3c,dlamb4c,Xi,sigmaprimal,sigmadual,mi);

                elseif Estrategia == 4    
                    %% ESTRATÉGIA 4 - NOVA DIREÇÃO COMB. CONVEXA

                    [x0,z,eta,lambda,PC] = FPOR_Novo_Ponto_Nova_Direcao_4(x0,z,eta,lambda,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dzp,detap,dlambdap,dx,dz,deta,dlambda,Xi,sigmaprimal,sigmadual,mi,PC,w1,w2);

    else    
                    %% ESTRATÉGIA 5 - NOVA DIREÇÃO COMB. LINEAR
    
                    [x0,z1,z2,z3,z4,lamb0,lamb1,lamb2,lamb3,lamb4,PC] = FPOR_Novo_Ponto_Nova_Direcao_5(x0,z1,z2,z3,z4,lamb0,lamb1,lamb2,lamb3,lamb4,alfa_primal_prev,alfa_dual_prev,alfa_primal_cor,alfa_dual_cor,dxp,dz1p,dz2p,dz3p,dz4p,dlamb0p,dlamb1p,dlamb2p,dlamb3p,dlamb4p,dxc,dz1c,dz2c,dz3c,dz4c,dlamb0c,dlamb1c,dlamb2c,dlamb3c,dlamb4c,Xi,sigmaprimal,sigmadual,mi,PC);

    end    

    x0_it = [x0_it, x0];       eta_it = [eta_it; eta];     lambda_it = [lambda_it; lambda];       z_it = [z_it; z];

    %% ATUALIZAÇÃO DOS DADOS DE BARRA E RAMO 

    [b,r] = FPOR_Atualizacao(b,r,x0,nb,ntap,ng,ngt,conj_ramo_tap,NIN,caso);

    %% ATUALIZAÇÃO DO MI E DELTA
    
    [mi,delta] = FPOR_Atualizacao_MI(mi,tal,z,lambda,H,G,Beta,delta,x0,tolz);                   mi_it = [mi_it;mi];
    

    %% CÁLCULO DO ERRO

    [Erro] = FPOR_Erro(mk,skprev,tk,uk,sk,PC,it);                                                     Erro_it = [Erro_it;Erro];
    Erro

end

tempo = cputime-to;

%% RESULTADOS - TABELAS E GRÁFICOS

[retorno_tabelas] = Tabela_Resultados(nb,nr,ng,ntap,ngt,E,mi,tal,Beta,alfa,Vmin,Vmax,tapmin,tapmax,Estrategia,Verificacao,vM,desv,kE,c,vi,vr,v0,d,Kp,Kr,it_it,f_it,FL_it,Erro_it,mi_it,Beta_it,TQ_it,PC,sing_it,fp_it,b,r,g,Ger,Cd_it,Cp_it,Cr_it,pg_eolica_it,Wr,fWeibull,x0,conj_ramo_tap,f,fe,fed,fp,it,tempo,GF,nodal,omegaR,omegaP);



















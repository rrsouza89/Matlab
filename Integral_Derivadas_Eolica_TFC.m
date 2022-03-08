
clc
clear all
format long
syms w vr vi v0 c wr d pg Kp Kr kE




%% Coeficiente do parâmetro de forma de distribuicao de Weibull

disp('Coeficiente do parâmetro de forma de distribuicao de Weibull: ')


%% FUNÇÃO DE DISTRIBUIÇÃO DE PROBABILIDADE DE WEIBULL

%% CASO W = 0

fW0 = 1 - exp(-(vi/c).^kE) + exp(-(v0/c).^kE);

%% CASO W

fWw = ((kE*((vr-vi)/vi)*vi)/c)*((((1+(w/wr)*((vr-vi)/vi))*vi)/c).^(kE-1))*exp(-(((1+(w/wr)*((vr-vi)/vi))*vi)/c).^kE);        % Hetzer

fWw2 = ((kE*(vr-vi))/(c*wr))*(((vi+(w/wr)*(vr-vi))/(c)).^(kE-1))*exp(-((vi+(w/wr)*(vr-vi))/(c)).^kE);                        % Elattar

fWw3 = (kE/c)*((((1+(w/wr)*((vr-vi)/vi))*vi)/c).^(kE-1))*exp(-(((1+(w/wr)*((vr-vi)/vi))*vi)/c).^kE);                         % (k/c)                                                                                                     % (k/c)

fWw4 = ((kE*(vr-vi))/(c*wr))*((vi*wr+w*(vr-vi))/(c*wr)).^(kE-1)*exp(-((vi*wr+w*(vr-vi))/(c*wr)).^kE);                       % Mishra

%% CASO W = Wr

fWr = exp(-(vr/c).^kE) - exp(-(v0/c).^kE);

%% INTEGRAL

disp('Função de Distribuição de Probabilidade de Weibull: ')
fW = fWr

aux = int(fW)                                     


aux_wr = subs(aux,w,wr);
aux_pg = subs(aux,w,pg);                                                    % substituindo wi na resolução da integral de w*fW(w)
aux_0 = subs(aux,w,0);                 % substituindo  0 na resolução da integral de w*fW(w)

fWpg = subs(fW,w,pg);


%% CUSTO MANUTENÇÃO

disp('Custo de Manutenção: ')
Cd = d*pg


%% CUSTO DE PENALIDADE E RESERVA

disp('Custo de Penalidade e Reserva: ')
Cw = (pg - w)*fW


%% DERIVADA PRIMEIRA DE Fw EM Pw

Cr = (Kr*(aux_pg - aux_0));

Cp = (Kp*(aux_pg - aux_wr));

Derivada_Primeira = d + Cr + Cp




%% DERIVADA SEGUNDA DE Fw EM Pw

Derivada_Segunda = (Kr*fWpg) + (Kp*fWpg)









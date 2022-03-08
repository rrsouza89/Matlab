
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

fWw = ((kE*((vr-vi)/vi)*vi)/c)*((((1+(w/wr)*((vr-vi)/vi))*vi)/c).^(kE-1))*exp(-(((1+(w/wr)*((vr-vi)/vi))*vi)/c).^kE)        % Hetzer

fWw2 = ((kE*(vr-vi))/(c*wr))*(((vi+(w/wr)*(vr-vi))/(c)).^(kE-1))*exp(-((vi+(w/wr)*(vr-vi))/(c)).^kE)                        % Elattar

fWw3 = (kE/c)*((((1+(w/wr)*((vr-vi)/vi))*vi)/c).^(kE-1))*exp(-(((1+(w/wr)*((vr-vi)/vi))*vi)/c).^kE)                         % (k/c)                                                                                                     % (k/c)

fWw4 = ((kE*(vr-vi))/(c*wr))*((vi*wr+w*(vr-vi))/(c*wr)).^(kE-1)*exp(-((vi*wr+w*(vr-vi))/(c*wr)).^kE)

%% CASO W = Wr

fWr = exp(-(vr/c).^kE) - exp(-(v0/c).^kE);

%% INTEGRAL

disp('Função de Distribuição de Probabilidade de Weibull: ')
fW = fWw4

fApprox = taylor(w*fW, w, 'Order', 2);    % ORDEM                                     

aux1 = int(fApprox,w) 
aux2 = int(fW);                     % resolução da integral indefinida de fW(w)

aux11 = subs(aux1,w,wr);
aux12 = subs(aux1,w,pg);                                                    % substituindo wi na resolução da integral de w*fW(w)
aux13 = subs(aux1,w,0);                 % substituindo  0 na resolução da integral de w*fW(w)
aux21 = subs(aux2,w,wr);  
aux22 = subs(aux2,w,pg);                                                    % substituindo wi na resolução da integral de fW(w)
aux23 = subs(aux2,w,0);


%% CUSTO MANUTENÇÃO

disp('Custo de Manutenção: ')
Cw = d*pg


%% CUSTO DE PENALIDADE

disp('Custo de Penalidade: ')
Cp = Kp*(aux11 - aux12 - pg*(aux21 - aux22))
%Cp = aux11 - aux12                            %Pfw
%Cp = aux21 - aux22                            %fw

%% CUSTO DE RESERVA

disp('Custo de Reserva: ')
Cr = Kr*(pg*(aux22 - aux23) - (aux12 - aux13))
%Cr = aux12 - aux13                            %Pfw
%Cr = aux22 - aux23                             %fw

%% FUNÇÃO CUSTO DA ENERGIA EÓLICA

disp('Função Custo da Energia Eólica: ')
Fw = Cw + Cp + Cr 

%% DERIVADA PRIMEIRA DE Fw EM Pw

disp('Derivada Primeira: ')
dFw = diff(Fw,'pg')


%% DERIVADA SEGUNDA DE Fw EM Pw

disp('Derivada Segunda: ')
d2Fw = diff(dFw,'pg')


%% COEFICIENTES Kp e Kr

% Pp = (pg*(wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - wr*exp(-(vr/c)^kE)) - (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) + (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) + (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2) - (kE*pg^3*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2*wr));
% 
% dPp = diff(Pp,'pg');
% 
% Pr = ((kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) - pg*(wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - wr*exp(-(vi/c)^kE)) + (kE*pg^3*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2*wr));
%         
% dPr = diff(Pr,'pg');










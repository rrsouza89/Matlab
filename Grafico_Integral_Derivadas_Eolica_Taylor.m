
clc
clear all
format long

%% PARÂMETROS INICIAIS

kE = 2;      c = 10;                                                        % k 2.39            1.88              2.81                3.07             2.24
                                                                            % c 5.27            3.93              6.30                7.48             5.78
vi = 3;     vr = 10.28;        v0 = 25;   

d = 0;      Kp = 2;     Kr = 5;

wr = 1;     w = 0:0.01:1;       pg = w;
%wr = 0.5;     w = 0:0.025:0.5;       pg = w;


fprintf('ESTUDO DA FUNÇÃO DE DENSIDADE DE PROBABILIDADE DE WEIBULL \n\n')

%% FUNÇÃO DE DENSIDADE DE PROBABILIDADE DE WEIBULL EM FUNÇÃO DA VELOCIDADE DO VENTO 

figure(1)
hold on

v = 0:1:100;

fv = (kE/c).*((v/c).^(kE-1)).*(exp(-(v/c).^kE));

plot(v,fv)

title('Função de Densidade de Probabilidade de Weibull')
xlabel('Vento')
ylabel('Probabilidade')

axis([0 30 0 1])


hold off

vento = [v' fv'];

fprintf('Vento        Probabilidade de Ocorrência\n\n');
fprintf('$%0.0f$    &\t     $%0.7f$ \\ \t \n',vento');
fprintf('Soma das probabilidades: %0.2f\n\n',sum(fv));


%% FUNÇÃO DE DENSIDADE DE PROBABILIDADE DE WEIBULL EM FUNÇÃO DA POTÊNCIA EÓLICA

figure(2)
hold on


fw1 = (kE./c).*((((1+(w./wr).*((vr-vi)./vi)).*vi)./c).^(kE-1)).*exp(-(((1+(w./wr).*((vr-vi)./vi)).*vi)./c).^kE);     %k/c

fw2 = ((kE.*((vr-vi)./vi).*vi)./c).*((((1+(w./wr).*((vr-vi)./vi)).*vi)./c).^(kE-1)).*exp(-(((1+(w./wr).*((vr-vi)./vi)).*vi)./c).^kE);    %Hetzer

fw3 = ((kE.*(vr-vi))./(c.*wr)).*(((vi+(w./wr).*(vr-vi))./(c)).^(kE-1)).*exp(-((vi+(w./wr).*(vr-vi))./(c)).^kE);    %Ellatar

fw4 = ((kE.*(vr-vi))./(c.*wr)).*((vi.*wr+w.*(vr-vi))./(c.*wr)).^(kE-1).*exp(-((vi.*wr+w.*(vr-vi))./(c.*wr)).^kE);    %Mishra

plot(w,fw1)
plot(w,fw2)
plot(w,fw3)
plot(w,fw4)

title('Função de Densidade de Probabilidade de Weibull')
xlabel('Potência Eólica')
ylabel('Probabilidade')

axis([0 1 0 5])

legend('f1','f2','f3','f4')

hold off

result = [w' fw1' fw2' fw3' fw4'];

fprintf('Potência            fW1                   fW2                      fW3            fW4  \n\n');
fprintf('$%0.3f$    &\t     $%0.7f$   &\t     $%0.7f$     &\t     $%0.7f$     &\t     $%0.7f$ \\ \t \n',result');
fprintf('Soma das probabilidades fw1: %0.2f\n',sum(fw1));
fprintf('Soma das probabilidades fw2: %0.2f\n',sum(fw2));
fprintf('Soma das probabilidades fw3: %0.2f\n\n',sum(fw3));
fprintf('Soma das probabilidades fw4: %0.2f\n\n',sum(fw4));



%% CÁLCULO DOS CUSTOS LINEAR, PENALIDADE E RESERVA

%% PENALIDADE

fprintf('CUSTO DE PENALIDADE\n\n');

    %% fW1 - k/c

%     fw_p = (wr.*exp(-(vr./c).^kE))./(vi - vr) - (wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE))./(vi - vr);                
% 
%     Pfw_p = (kE.*wr.^2.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*vi) - (kE.*pg.^2.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*vi);
    
    %% fW2 - Hetzer
    
%     fw_p = wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE) - wr.*exp(-(vr./c).^kE);
% 
%     Pfw_p = (kE.*pg.^2.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./(2.*vi) - (kE.*wr.^2.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./(2.*vi);
    
    %% fW3 - Ellatar
    
    fw_p = exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE) - exp(-(vr./c).^kE);
                
    Pfw_p = (kE.*pg.^2.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*wr) - (kE.*wr.*exp(-(vi./c).^kE).*(vi./c).^kE)./2 + (kE.*vr.*wr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*vi) - (kE.*pg.^2.*vr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*vi.*wr);
                

    
    pgfw_p = pg.*fw_p;

    propP = abs(Pfw_p - pgfw_p);

    Penalidade = [pg' fw_p' Pfw_p' pgfw_p' propP'];
    
    fprintf('Potência            fW                     PfW                  pg*fW                Cp     \n\n');
    fprintf('$%0.3f$    &\t     $%0.7f$     &\t     $%0.7f$     &\t $%0.7f$     &\t     $%0.7f$ \\ \t \n',Penalidade');
    
    figure(3)
    hold on
    plot(pg,Pfw_p)
    plot(pg,pgfw_p)
    plot(pg,propP)
    
    title('Custo de Penalidade')
    xlabel('Potência Eólica')
    ylabel('Probabilidade')

    axis([0 1 0 1])
    %axis([0 0.5 0 0.5])

    legend('PfW','pgfW','Cp')
    
    hold off
    
    
%% RESERVA

fprintf('CUSTO DE RESERVA\n\n');

    %% fW1 - k/c

%     fw_r = (wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE))./(vi - vr) - (wr.*exp(-(vi./c).^kE))./(vi - vr);
%                 
%     Pfw_r = (kE.*pg.^2.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*vi);
    
    %% fW2 - Hetzer
    
%     fw_r = wr.*exp(-(vi./c).^kE) - wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE);
% 
%     Pfw_r = -(kE.*pg.^2.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./(2.*vi);
    
    %% fW3 - Ellatar
    
    fw_r = exp(-(vi./c).^kE) - exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE);
                
    Pfw_r = (kE.*pg.^2.*vr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*vi.*wr) - (kE.*pg.^2.*exp(-(vi./c).^kE).*(vi./c).^kE)./(2.*wr);
                 

    
    pgfw_r = pg.*fw_r;

    propR = abs(pgfw_r - Pfw_r);

    Reserva = [pg' fw_r' pgfw_r' Pfw_r' propR'];
    
    fprintf('Potência            fW                     pg*fW                PfW                  Cr     \n\n');
    fprintf('$%0.3f$    &\t     $%0.7f$     &\t     $%0.7f$     &\t $%0.7f$     &\t     $%0.7f$ \\ \t \n',Reserva');
    
    figure(4)
    hold on
    plot(pg,pgfw_r)
    plot(pg,Pfw_r)
    plot(pg,propR)
    
    title('Custo de Reserva')
    xlabel('Potência Eólica')
    ylabel('Probabilidade')

    axis([0 1 0 1])
    %axis([0 0.5 0 0.5])

    legend('pgfW','PfW','Cr')
    
    hold off
    
    
    %% CUSTO EÓLICA
    
    Ct = pg + propP + propR;
    
    Total = [pg' pg' propP' propR' Ct']
    
    fprintf('Potência              Cd                     Cp                    Cr                  Ctotal     \n\n');
    fprintf('$%0.3f$    &\t     $%0.7f$     &\t     $%0.7f$     &\t $%0.7f$     &\t     $%0.7f$ \\ \t \n',Total');
    
    figure(5)
    hold on
    plot(pg,pg)
    plot(pg,propP)
    plot(pg,propR)
    plot(pg,Ct)
    
    title('Custo de Eólica')
    xlabel('Potência Eólica')
    ylabel('Custo')

    axis([0 1 0 1.5])
    %axis([0 0.5 0 0.7])

    legend('Cd','Cp','Cr','Ct')
    
    hold off
    
    
    %% GRADIENTE E HESSIANA
    
        %% fW1 - k/c
        
%         G = d - Kr.*((wr.*exp(-(vi./c).^kE))./(vi - vr) - (wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE))./(vi - vr) - (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c + (kE.*pg.*exp(-(vi/c).^kE).*(vi./c).^kE)./vi) - Kp.*((wr.*exp(-(vr./c).^kE))./(vi - vr) - (wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE))./(vi - vr) - (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c + (kE.*pg.*exp(-(vi./c).^kE).*(vi./c).^kE)./vi);
%     
%         H = - Kp.*((kE.*exp(-(vi./c).^kE).*(vi./c).^kE)./vi - (2.*kE.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c - (kE.^2.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(2.*kE - 2))./(c.^2.*wr) + (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*(kE - 1).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 2))./(c.^2.*wr)) - Kr.*((kE.*exp(-(vi./c).^kE).*(vi./c).^kE)./vi - (2.*kE.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c - (kE.^2.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(2.*kE - 2))./(c.^2.*wr) + (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*(kE - 1).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 2))./(c.^2.*wr));

        %% fW2 - Hetzer
        
%         G = d - Kr.*(wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE) - wr.*exp(-(vi./c).^kE) - (kE.*pg.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./vi + (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c) - Kp.*(wr.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE) - wr.*exp(-(vr./c).^kE) - (kE.*pg.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./vi + (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c);
%         
%         H = - Kp.*((2.*kE.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c - (kE.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./vi + (kE.^2.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(2.*kE - 2))./(c.^2.*wr) - (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*(kE - 1).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 2))./(c.^2.*wr)) - Kr.*((2.*kE.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./c - (kE.*exp(-(vi./c).^kE).*(vi./c).^kE.*(vi - vr))./vi + (kE.^2.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(2.*kE - 2))./(c.^2.*wr) - (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*(kE - 1).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 2))./(c.^2.*wr))
        
        
        %% fW3 - Ellatar
        
        G = d - Kr.*(exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE) - exp(-(vi./c).^kE) - (kE.*pg.*exp(-(vi./c).^kE).*(vi./c).^kE)./wr + (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./(c.*wr) + (kE.*pg.*vr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(vi.*wr)) - Kp.*(exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE) - exp(-(vr./c).^kE) - (kE.*pg.*exp(-(vi./c).^kE).*(vi./c).^kE)./wr + (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./(c.*wr) + (kE.*pg.*vr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(vi.*wr));
        
        H = - Kp.*((2.*kE.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./(c.*wr) - (kE.*exp(-(vi./c).^kE).*(vi./c).^kE)./wr + (kE.*vr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(vi.*wr) + (kE.^2.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(2.*kE - 2))./(c.^2.*wr.^2) - (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*(kE - 1).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 2))./(c.^2.*wr.^2)) - Kr.*((2.*kE.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 1))./(c.*wr) - (kE.*exp(-(vi./c).^kE).*(vi./c).^kE)./wr + (kE.*vr.*exp(-(vi./c).^kE).*(vi./c).^kE)./(vi.*wr) + (kE.^2.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(2.*kE - 2))./(c.^2.*wr.^2) - (kE.*pg.*exp(-((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^kE).*(vi - vr).^2.*(kE - 1).*((pg.*vr - pg.*vi + vi.*wr)./(c.*wr)).^(kE - 2))./(c.^2.*wr.^2));


        
        Derivada = [pg' G' H']
        
        fprintf('Potência            Derivada Primeira                 Derivada Segunda  \n\n');
        fprintf('$%0.3f$    &\t     $%0.7f$     &\t         $%0.7f$ \\ \t \n',Derivada');
        
        figure(6)
        hold on
        plot(pg,G)
           
        title('Derivada Primeira')
        xlabel('Potência Eólica')
        ylabel('Valor')

        %axis([0 1 -2.5 1.2])
        axis([0 0.5 -2.5 1.2])

        legend('G')
    
        hold off
        
        
        figure(7)
        hold on
        plot(pg,H)
           
        title('Derivada Segunda')
        xlabel('Potência Eólica')
        ylabel('Valor')

        axis([0 0.5 -12 12])
        %axis([0 0.5 0 0.7])

        legend('H')
    
        hold off
    


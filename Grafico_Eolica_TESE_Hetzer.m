        clc
        clear all 
        close all
    
%         kE = 2;     c = 5;
%         kE2 = 2;     c2 = 7;
%         kE3 = 2;    c3 = 10;
%         kE4 = 2;    c4 = 15;
        
        kE = 2;     c = 10;
        %kE2 = 5;     c2 = 10;
        %kE3 = 5;    c3 = 15;
        %kE4 = 2;    c4 = 10;
        
%         kE = 2.39;   c = 5.27;
%         kE2 = 1.88;  c2 = 3.93;
%         kE3 = 2.81;  c3 = 6.30;
%         kE4 = 3.07;  c4 = 7.48;



        
        vi = 3;
        vi2 = 3;
        vi3 = 3;
                
        vr = 15;
        vr2 = 15;
        vr3 = 15;
        
        v0 = 30;
        wr = 1;
        
        %% FUNÇÃO DE DISTRIBUIÇÃO DE WEIBULL
        
        w = 0:0.001:wr;
        
        %w1
        w0 = (1 - exp(-(vi/c).^kE) + exp(-(v0/c).^kE)); 
        fw = ((kE.*(vr-vi))/(c)).*((((1+(w/wr).*((vr-vi)/(vi))).*vi)/(c)).^(kE-1)).*exp(-((((1+(w/wr).*((vr-vi)/(vi))).*vi)/(c))).^kE);
        wr1 = (exp(-(vr/c).^kE) - exp(-(v0/c).^kE)); 
        
        %w2
%         w02 = 1 - exp(-(vi2/c2).^kE2) + exp(-(v0/c2).^kE2);
%         fw2 = (-(kE2.*exp(-((vi2.*wr - w.*(vi2 - vr2))/(c2.*wr)).^kE2).*(vi2 - vr2).*((vi2.*wr - w.*(vi2 - vr2))/(c2.*wr)).^(kE2 - 1))/(c2.*wr))/10;
%         wr2 = exp(-(vr/c2).^kE2) - exp(-(v0/c2).^kE2); 
        
        %w3         
%         w03 = 1 - exp(-(vi3/c3).^kE3) + exp(-(v0/c3).^kE3);
%         fw3 = (-(kE3.*exp(-((vi3.*wr - w.*(vi3 - vr3))/(c3.*wr)).^kE3).*(vi3 - vr3).*((vi3.*wr - w.*(vi3 - vr3))/(c3.*wr)).^(kE3 - 1))/(c3.*wr))/10;
%         wr3 = exp(-(vr/c3).^kE3) - exp(-(v0/c3).^kE3); 
        
        %w4
        %w04 = 1 - exp(-(vi4/c4).^kE4) + exp(-(v0/c4).^kE4);
        %fw4 = (-(kE4.*exp(-((vi.*wr - w.*(vi - vr))/(c4.*wr)).^kE4).*(vi - vr).*((vi.*wr - w.*(vi - vr))/(c4.*wr)).^(kE4 - 1))/(c4.*wr))/10;
                  
            
        figure(1)
        hold on
        
        plot(w,fw,'b')
        %plot(w,fw2,'r')
        %plot(w,fw3,'k')
        %plot(w,fw4,'c')
        
        plot(0,w0,'bo')
        plot(wr,wr1,'bo')
        
%         plot(0,w02,'r+')
%         plot(wr,wr2,'r+')
        
%         plot(0,w03,'k*')
%         plot(wr,wr3,'k*')
        
        
        %plot(w,fw4,'k')
        
        axis([0 wr 0 1])
        xlabel('Wind power (pu)')
        ylabel('Probability')
        
        %legend('k = 2.39 c = 5.27','k = 1.88 c = 3.93','k = 2.81 c = 6.30','k = 3.07 c = 7.48')
        %legend('k = 2 c = 10','+ k = 5 c = 10','* k = 5 c = 15')
        legend('k = 2 c = 10')
        
        %legend('k = 2 c = 5','k = 2 c = 7','k = 2 c = 10','k = 2 c = 15')
        %legend('k = 1 c = 10','k = 2 c = 10','k = 3 c = 10','k = 5 c = 10')
        %legend('vr = 10','vr = 15','vr = 20')
        
        %title('FDPW conforme a potência eólica')
        %title('Variação do parâmetro de forma k')
        %title('Variação da velocidade nominal da turbina eólica')
        %title('Custos da geração eólica')
        %title('FDPW conforme P^G_{49}')
        
        hold off
                
          
        
        
        %% FUNÇÃO DE DENSIDADE DE WEIBULL
        
%         v = 0:0.01:100;
% 
%         fv = (kE/c).*((v/c).^(kE-1)).*(exp(-(v/c).^kE));
%         fv2 = (kE2/c2).*((v/c2).^(kE2-1)).*(exp(-(v/c2).^kE2));
%         fv3 = (kE3/c3).*((v/c3).^(kE3-1)).*(exp(-(v/c3).^kE3));
%         %fv4 = (kE4/c4).*((v/c4).^(kE4-1)).*(exp(-(v/c4).^kE4));
%              
%         
%         figure(2)
%         hold on
%         plot(v,fv,'b')
%         plot(v,fv2,'r')
%         plot(v,fv3,'k')
%         %plot(v,fv4,'-')
%         axis([0 30 0 0.4])
%         xlabel('Velocidade do vento (m/s)')
%         ylabel('Probabilidade')
%         %legend('k = 2.39 c = 5.27','k = 1.88 c = 3.93','k = 2.81 c = 6.30','k = 3.07 c = 7.48')
%         legend('k = 5 c = 5','k = 5 c = 10','k = 5 c = 15')
%         title('FDPW conforme a velocidade do vento')
%         hold off
%         
%         
%         %% FUNÇÃO DE DISTRIBUIÇÃO ACUMULADA
%         
%                
%         fda = 1 - exp(-(v/c).^kE);
%         fda2 = 1 - exp(-(v/c2).^kE2);
%         fda3 = 1 - exp(-(v/c3).^kE3);
%         
%         figure(3)
%         hold on
%         plot(v,fda,'b')
%         plot(v,fda2,'r')
%         plot(v,fda3,'k')
%         axis([0 50 0 1.2])
%         xlabel('Velocidade do vento (m/s)')
%         ylabel('Probabilidade acumulada')
%         legend('k = 5 c = 5','k = 5 c = 10','k = 5 c = 15')
%         title('Função de distribuição acumulada')
%         hold off
%     
        
       
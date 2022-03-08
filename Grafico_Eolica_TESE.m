        clc
        clear all 
        close all
    
        kE = 2;     c = 5;
        kE2 = 2;    c2 = 10;
        kE3 = 2;    c3 = 15;
        
%         kE = 2.39;   c = 5.27;
%         kE2 = 1.88;  c2 = 3.93;
%         kE3 = 2.81;  c3 = 6.30;
%         kE4 = 3.07;  c4 = 7.48;



        
        vi = 3
        vr = 15
        v0 = 30
        wr = 1
        
        %% FUNÇÃO DE DISTRIBUIÇÃO DE WEIBULL
        
        w = 0:0.01:wr;
        
        %w1
        w0 = 1 - exp(-(vi/c).^kE) + exp(-(v0/c).^kE); 
        fw = (kE.*exp(-(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c).^kE).*(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c).^(kE - 1))/c;
        wr1 = exp(-(vr/c).^kE) - exp(-(v0/c).^kE); 
        
        %w2
        w02 = 1 - exp(-(vi/c2).^kE2) + exp(-(v0/c2).^kE2);
        fw2 = (kE2.*exp(-(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c2).^kE2).*(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c2).^(kE2 - 1))/c2;
        wr2 = exp(-(vr/c2).^kE2) - exp(-(v0/c2).^kE2); 
        
        %w3                 k/c
        w03 = 1 - exp(-(vi/c3).^kE3) + exp(-(v0/c3).^kE3);
        fw3 = (kE3.*exp(-(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c3).^kE3).*(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c3).^(kE3 - 1))/c3;
        wr3 = exp(-(vr/c3).^kE3) - exp(-(v0/c3).^kE3); 
        
        
        %fw4 = (kE4.*exp(-((vi - (w.*(vi - vr))/wr)/c4).^kE4).*((vi - (w.*(vi - vr))/wr)/c4).^(kE4 - 1))/c4;
                  
            
        figure(1)
        hold on
        
        plot(w,fw,'b')
        plot(w,fw2,'r')
        plot(w,fw3,'k')
        
        plot(0,w0,'bo')
        plot(wr,wr1,'bo')
        
        plot(0,w02,'r+')
        plot(wr,wr2,'r+')
        
        plot(0,w03,'k*')
        plot(wr,wr3,'k*')
        
        
        %plot(w,fw4,'k')
        
        axis([0 wr 0 0.4])
        xlabel('Potência Eólica (pu)')
        ylabel('Probabilidade')
        %legend('k = 2.39 c = 5.27','k = 1.88 c = 3.93','k = 2.81 c = 6.30','k = 3.07 c = 7.48')
        legend('o k = 2 c = 5','+ k = 2 c = 10','* k = 2 c = 15')
        title('FDPW conforme a potência eólica')
        hold off
                
          
        
        
        %% FUNÇÃO DE DENSIDADE DE WEIBULL
        
        v = 0:0.01:100;

        fv = (kE/c).*((v/c).^(kE-1)).*(exp(-(v/c).^kE));
        fv2 = (kE2/c2).*((v/c2).^(kE2-1)).*(exp(-(v/c2).^kE2));
        fv3 = (kE3/c3).*((v/c3).^(kE3-1)).*(exp(-(v/c3).^kE3));
        %fv4 = (kE4/c4).*((v/c4).^(kE4-1)).*(exp(-(v/c4).^kE4));
             
        
        figure(2)
        hold on
        plot(v,fv,'k-')
        plot(v,fv2,'k--')
        plot(v,fv3,'k-.')
        %plot(v,fv4,'-')
        axis([0 30 0 0.4])
        xlabel('Velocidade do vento (m/s)')
        ylabel('Probabilidade')
        %legend('k = 2.39 c = 5.27','k = 1.88 c = 3.93','k = 2.81 c = 6.30','k = 3.07 c = 7.48')
        legend('k = 3 c = 5','k = 3 c = 10','k = 3 c = 15')
        title('FDPW conforme a velocidade do vento')
        hold off
        
        
        %% FUNÇÃO DE DISTRIBUIÇÃO ACUMULADA
        
               
        fda = 1 - exp(-(v/c).^kE);
        fda2 = 1 - exp(-(v/c2).^kE2);
        fda3 = 1 - exp(-(v/c3).^kE3);
        
        figure(3)
        hold on
        plot(v,fda,'k-')
        plot(v,fda2,'k--')
        plot(v,fda3,'k-.')
        axis([0 50 0 1.2])
        xlabel('Velocidade do vento (m/s)')
        ylabel('Probabilidade acumulada')
        legend('k = 3 c = 5','k = 3 c = 10','k = 3 c = 15')
        title('Função de distribuição acumulada')
        hold off
    
        
        %% CUSTOS DA GERAÇÂO
        
%         fW = fw2;
%         
%         
%         
%         figure(4)
%         hold on
%         plot(w,fW,'k--')
%         
%         w=0.5;        
%         fWw = (kE2.*exp(-((vi - (w.*(vi - vr))/wr)/c2).^kE2).*((vi - (w.*(vi - vr))/wr)/c2).^(kE2 - 1))/c2;
%         
%         plot(w,fWw,'kx')
%         axis([0 wr 0 0.4])
%         
%         xlabel('Potência Eólica (pu)')
%         ylabel('Probabilidade')
%         legend('k = 2 c = 10')
%         title('Custos da geração eólica')
%         
%         hold off
        
        
        %% CUSTOS LINEAR, RESERVA E PENALIDADE
        
%         d = 1; Kr = 1; Kp = 1;  kE = 2; c = 10; wr = 1;
%         
%         w = 0:0.01:wr;
%         
%         pg = 0:0.01:wr;
%         
%         Cd = d.*w;
%         
%         
%         w1 = linspace(wr,0,101);  
%         fw1e = (kE.*exp(-(-(vi.*((w1.*(vi - vr))./(vi.*wr) - 1))./c).^kE).*(pg - w1).*(-(vi.*((w1.*(vi - vr))./(vi.*wr) - 1))./c).^(kE - 1))./c;
%         Ip = trapz(w1,fw1e)
%         Cp = (Kp*Ip)
%         
%         w = linspace(0,wr,101);                
%         fwe = (kE.*exp(-(-(vi.*((w.*(vi - vr))./(vi.*wr) - 1))./c).^kE).*(pg - w).*(-(vi.*((w.*(vi - vr))./(vi.*wr) - 1))./c).^(kE - 1))./c;
%         Ir = trapz(w,fwe)
%         Cr = (Kr*Ir)
%         
%         Ct = Cd + Cp + Cr
%         
%         figure(5)
%         hold on
%         plot(w,Cd,'k.')
%         plot(w,fwe,'k-')
%         plot(w,fw1e,'k--')
%         plot(w,Ct,'k--.')
%         
%         
%         hold off
        
        
        
    

    



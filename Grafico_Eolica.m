function [retorno] = Grafico_Eolica(kE,c,g,ng,b,Ger,fWeibull);
    
    %kE2 = 1.88; c2 = 3.93;
    %kE3 = 2.81; c3 = 6.30;
    %kE4 = 3.07; c4 = 7.48;


for i=1:ng
    
    if b(Ger(i)).TGER == 2
        
        vi = g(i).Vi;
        vr = g(i).Vr;
        v0 = g(i).V0;
        wr = g(i).Wr;
        
        %% FUNÇÃO DE DISTRIBUIÇÃO DE WEIBULL
        
        w = 0:0.01:wr;
        
                       
            if fWeibull == 1            % Hetzer
    
                fw = -(kE.*exp(-(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c).^kE).*(vi - vr).*(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c).^(kE - 1))/c;  
            
            elseif fWeibull == 2        % Mishra
                
                fw = (-(kE.*exp(-((vi.*wr - w.*(vi - vr))/(c.*wr)).^kE).*(vi - vr).*((vi.*wr - w.*(vi - vr))/(c.*wr)).^(kE - 1))/(c.*wr));
                
            elseif fWeibull == 3        % k/c
                
                fw = (kE.*exp(-(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c).^kE).*(-(vi.*((w.*(vi - vr))/(vi.*wr) - 1))/c).^(kE - 1))/c;
                
                
            end
            
            
        figure(1)
        hold on
        plot(w,fw,'-')
        %plot(w,fw2,'-')
        %plot(w,fw3,'-')
        %plot(w,fw4,'-')
        axis([0 wr 0 2])
        xlabel('Potência Eólica (pu)')
        ylabel('Probabilidade')
        %legend('k = 2.39 c = 5.27','k = 1.88 c = 3.93','k = 2.81 c = 6.30','k = 3.07 c = 7.48')
        title('FDPW em função da potência eólica')
        hold off
                
          
        
        
        %% FUNÇÃO DE DENSIDADE DE WEIBULL
        
%         v = 0:0.001:100;
% 
%         fv = (kE/c).*((v/c).^(kE-1)).*(exp(-(v/c).^kE));
%         %fv2 = (kE2/c2).*((v/c2).^(kE2-1)).*(exp(-(v/c2).^kE2));
%         %fv3 = (kE3/c3).*((v/c3).^(kE3-1)).*(exp(-(v/c3).^kE3));
%         %fv4 = (kE4/c4).*((v/c4).^(kE4-1)).*(exp(-(v/c4).^kE4));
%              
%         
%         figure(2)
%         hold on
%         plot(v,fv,'-')
%         %plot(v,fv2,'-')
%         %plot(v,fv3,'-')
%         %plot(v,fv4,'-')
%         axis([0 30 0 0.3])
%         xlabel('Velocidade do vento (m/s)')
%         ylabel('Probabilidade')
%         %legend('k = 2.39 c = 5.27','k = 1.88 c = 3.93','k = 2.81 c = 6.30','k = 3.07 c = 7.48')
%         title('FDPW em função do vento')
%         hold off
    
    end

    

end

retorno = 1;


end
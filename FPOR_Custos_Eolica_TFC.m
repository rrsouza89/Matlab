function [Cd,Cp,Cr,Kp,Kr,GF,vi,vr,v0,d,g] = FPOR_Custos_Eolica_TFC(g,c,kE,pg,i,fWeibull,L,GF,nb,ntap,omegaR,omegaP);

        pg = pg;
                       
        vi = g(i).Vi;
        vr = g(i).Vr;
        v0 = g(i).V0;
        wr = g(i).Wr;
        d = g(i).d*100;     
        Kp = g(i).Kp*100;
        Kr = g(i).Kr*100;        
        
        
        %% CUSTO LINEAR
        
        Cd = (d*pg)/10000;
        
        
        if pg < L
            
            
            %% MISHRA
        
            %% CUSTO DE PENALIDADE
        
                w1 = linspace(wr,pg);  
           
                fw1e = (pg - w1).*(exp(-(v0./c).^kE) - exp(-(vi./c).^kE) + 1);
                
                Ip = trapz(w1,fw1e);
      
                Cp = (omegaP*(Kp*Ip))/10000;            
                   
            %% CUSTO DE RESERVA  
        
                w = linspace(0,pg);
                
                fwe = (pg - w).*(exp(-(v0./c).^kE) - exp(-(vi./c).^kE) + 1);
                
                Ir = trapz(w,fwe);
                         
                Cr = (omegaR*(Kr*Ir))/10000;            
                
            %% GRADIENTE DERIVADA PRIMEIRA 
        
                GF(2*nb+ntap+i,1) = d/10000;
                
                
        elseif (wr - pg) < L
            
            
            
            %% MISHRA
        
            %% CUSTO DE PENALIDADE
        
                w1 = linspace(wr,pg);  
           
                fw1e = -(exp(-(v0./c).^kE) - exp(-(vr./c).^kE)).*(pg - w1);
                
                Ip = trapz(w1,fw1e);
      
                Cp = (omegaP*(Kp*Ip))/10000;            
                   
            %% CUSTO DE RESERVA  
        
                w = linspace(0,pg);
                
                fwe = -(exp(-(v0./c).^kE) - exp(-(vr./c).^kE)).*(pg - w);
                
                Ir = trapz(w,fwe);
                         
                Cr = (omegaR*(Kr*Ir))/10000;            
                
            %% GRADIENTE DERIVADA PRIMEIRA 
        
                GF(2*nb+ntap+i,1) = d/10000;
      
        else
        
        %% MISHRA
        
        %% CUSTO DE PENALIDADE
        
            w1 = linspace(wr,pg);  
           
            fw1e = -(kE.*exp(-((vi.*wr - w1.*(vi - vr))./(c.*wr)).^kE).*(pg - w1).*(vi - vr).*((vi.*wr - w1.*(vi - vr))./(c.*wr)).^(kE - 1))./(c.*wr);
                
            Ip = trapz(w1,fw1e);
      
            Cp = (omegaP*(Kp*Ip))/10000;            
                   
        %% CUSTO DE RESERVA  
        
            w = linspace(0,pg);
                
            fwe = -(kE.*exp(-((vi.*wr - w.*(vi - vr))./(c.*wr)).^kE).*(pg - w).*(vi - vr).*((vi.*wr - w.*(vi - vr))./(c.*wr)).^(kE - 1))./(c.*wr);
                
            Ir = trapz(w,fwe);
                         
            Cr = (omegaR*(Kr*Ir))/10000;            
                
        %% GRADIENTE DERIVADA PRIMEIRA 
        
            GF(2*nb+ntap+i,1) = d/10000 + (omegaR*(Kr*(exp(-(vi/c)^kE) - exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))))/10000 + (omegaP*(Kp*(exp(-(vr/c)^kE) - exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))))/10000;
            
            g(i).CIL = d/10000;
            g(i).CIP = (omegaP*(Kp*(exp(-(vr/c)^kE) - exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))))/10000;
            g(i).CIR = (omegaR*(Kr*(exp(-(vi/c)^kE) - exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))))/10000;
            
         end

end


function [Cd,Cp,Cr,Kp,Kr,GF,vi,vr,v0,d] = FPOR_Custos_Eolica(g,c,kE,pg,i,fWeibull,L,GF,nb,ntap);

        pg = pg;
                       
        vi = g(i).Vi;
        vr = g(i).Vr;
        v0 = g(i).V0;
        wr = g(i).Wr;
        d = g(i).d;     
        Kp = g(i).Kp;
        Kr = g(i).Kr;
        
        % Custo Linear
        Cd = (d*pg)/100
        
%         if pg < L
%             
%             %Weibull f ordem 3 fWeibull = 1 
%             %% CUSTO PENALIDADE
%             
%                 fw_p = 0
%                 
%                 Pfw_p = ((pg^2*exp(-(vi/c)^kE))/2 - (pg^2*exp(-(v0/c)^kE))/2 + (wr^2*exp(-(v0/c)^kE))/2 - (wr^2*exp(-(vi/c)^kE))/2 - pg^2/2 + wr^2/2)/100
%                 
%                 pgfw_p = pg*fw_p
%                 
%                 propP = Pfw_p - pgfw_p
%                 
%                 Cp = Kp*propP
%             
%             
%             %% CUSTO RESERVA
%                 
%                 fw_r = 0
%                 
%                 Pfw_r = ((pg^2*exp(-(v0/c)^kE))/2 - (pg^2*exp(-(vi/c)^kE))/2 + pg^2/2)/100
%                 
%                 pgfw_r = pg*fw_r
%                 
%                 %propR = pgfw_r - Pfw_r 
%                 propR = Pfw_r 
%                 
%                 Cr = Kr*propR
%             
%             %% GRADIENTE
%             GF(2*nb+ntap+i,1) = (d - Kp*(pg + pg*exp(-(v0/c)^kE) - pg*exp(-(vi/c)^kE)) - Kr*(pg + pg*exp(-(v0/c)^kE) - pg*exp(-(vi/c)^kE)))/10000;
%             
%         elseif (wr - pg) < L
%             
%             %Weibull f ordem 3 fWeibull = 1
%             %% CUSTO PENALIDADE
%             
%                 fw_p = 0
%                 
%                 Pfw_p = (exp(-(v0/c)^kE)/2 - exp(-(vr/c)^kE)/2)*pg^2 + (exp(-(vr/c)^kE)/2 - exp(-(v0/c)^kE)/2)*wr^2
%                 
%                 pgfw_p = pg*fw_p
%                 
%                 propP = Pfw_p - pgfw_p
%                 
%                 Cp = Kp*propP
%             
%             %% CUSTO RESERVA
%             
%                 fw_r = 0
%                 
%                 Pfw_r = -pg^2*(exp(-(v0/c)^kE)/2 - exp(-(vr/c)^kE)/2)
%                 
%                 pgfw_r = pg*fw_r
%                 
%                 %propR = pgfw_r - Pfw_r 
%                 propR = Pfw_r 
%                 
%                 Cr = Kr*propR
%             
%             %% GRADIENTE
%             GF(2*nb+ntap+i,1) = (d + 2*Kp*pg*(exp(-(v0/c)^kE)/2 - exp(-(vr/c)^kE)/2) + 2*Kr*pg*(exp(-(v0/c)^kE)/2 - exp(-(vr/c)^kE)/2))/10000;
%             
%         else
            
            if fWeibull == 1
        
            %% CUSTO DE PENALIDADE
            
                % Weibull f ordem 2
                fw_p = wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - wr*exp(-(vr/c)^kE)
                
                Pfw_p = (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr))/(2*vi) - (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr))/(2*vi)
                                 
            
                %Weibull f ordem 3
%                 fw_p = (wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - wr*exp(-(vr/c)^kE))
%                 
%                 Pfw_p = ((kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) - (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) - (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2) + (kE*pg^3*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2*wr))

    
                
                pgfw_p = pg*fw_p
                                
                propP = abs(Pfw_p - pgfw_p)
                
                Cp = (Kp*propP)/100
                
                
                %Weibull MÉDIA WEIBULL
                %Cp = Kp*((kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) - (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) - (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2) + (kE*pg^3*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2*wr));
                
            
            %% CUSTO DE RESERVA  
            
            
                %Weibull f ordem 2
                fw_r = wr*exp(-(vi/c)^kE) - wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)
                
                Pfw_r = -(kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr))/(2*vi)
            
                %Weibull f odem 3
%                 fw_r = (wr*exp(-(vi/c)^kE) - wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))
%                 
%                 Pfw_r = (- (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) - (kE*pg^3*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2*wr))
                
                         
                
                pgfw_r = pg*fw_r
                
                propR = abs(pgfw_r - Pfw_r)
                                         
                Cr = (Kr*propR)/100
                            
                %Weibull MÉDIA WEIBULL
                %Cr = -Kr*((kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/(2*c) + (kE*pg^3*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(3*vi^2*wr));
                
            %% GRADIENTE DERIVADA PRIMEIRA 
                        
                %Weibull f ordem 2                
                GF(2*nb+ntap+i,1) = (d - Kr*(wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - wr*exp(-(vi/c)^kE) - (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr))/vi + (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/c) - Kp*(wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - wr*exp(-(vr/c)^kE) - (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr))/vi + (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/c))/100;
                
                %Weibull f ordem 3
                %GF(2*nb+ntap+i,1) = (d + Kr*(wr*exp(-(vi/c)^kE) - wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) + (kE*pg*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/c - (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/c + (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(vi^2*wr)) + Kp*(wr*exp(-(vr/c)^kE) - wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) + (kE*pg*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/c - (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/c + (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(vi^2*wr)))/100;
                 
                                
                %Weibull MÉDIA WEIBULL
                %GF(2*nb+ntap+i,1) = (d + Kp*((kE*pg*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/c + (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(vi^2*wr)) - Kr*((kE*pg*exp(-(vi/c)^kE)*(vi/c)^(kE - 1)*(vi - vr))/c + (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE*(vi - vr)^2*(kE*(vi/c)^kE - kE + 1))/(vi^2*wr)))/10000;
            
            elseif fWeibull == 2
            
                %% CUSTO DE PENALIDADE
                
                % Weibull f ordem 2
                fw_p = exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - exp(-(vr/c)^kE)
                
                Pfw_p = (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*wr) - (kE*wr*exp(-(vi/c)^kE)*(vi/c)^kE)/2 + (kE*vr*wr*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi) - (kE*pg^2*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi*wr)
                
                
                pgfw_p = pg*fw_p
                                
                %propP = abs(Pfw_p - pgfw_p)
                propP = Pfw_p - pgfw_p
                
                Cp = (Kp*propP)/10
               
                
                %% CUSTO DE RESERVA
                
                %Weibull f ordem 2
                fw_r = exp(-(vi/c)^kE) - exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)
                
                Pfw_r = (kE*pg^2*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi*wr) - (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*wr)
                
                
                pgfw_r = pg*fw_r
                
                propR = abs(pgfw_r - Pfw_r)
                %propR = pgfw_r - Pfw_r
                                         
                Cr = (Kr*propR)/100     
                
                
                %% GRADIENTE DERIVADA PRIMEIRA
            
                
                GF(2*nb+ntap+i,1) = (d - Kr*(exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - exp(-(vi/c)^kE) - (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE)/wr + (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/(c*wr) + (kE*pg*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(vi*wr)) - Kp*(exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - exp(-(vr/c)^kE) - (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE)/wr + (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/(c*wr) + (kE*pg*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(vi*wr)))/100;
                 
            elseif fWeibull == 3
            
                %% CUSTO DE PENALIDADE
                
                % Weibull f ordem 2
                fw_p = (wr*exp(-(vr/c)^kE))/(vi - vr) - (wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))/(vi - vr)
                
                Pfw_p = (kE*wr^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi) - (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi)
                
                
                pgfw_p = pg*fw_p
                                
                propP = abs(Pfw_p - pgfw_p)
                
                Cp = (Kp*propP)
            
                
                %% CUSTO DE RESERVA
                
                %Weibull f ordem 2
                fw_r = (wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))/(vi - vr) - (wr*exp(-(vi/c)^kE))/(vi - vr)
                
                Pfw_r = (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi)
                
                
                pgfw_r = pg*fw_r
                
                propR = abs(pgfw_r - Pfw_r)
                                         
                Cr = (Kr*propR)/10               
                
            
                %% GRADIENTE DERIVADA PRIMEIRA
            
                
                GF(2*nb+ntap+i,1) = (d - Kr*((wr*exp(-(vi/c)^kE))/(vi - vr) - (wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))/(vi - vr) - (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/c + (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE)/vi) - Kp*((wr*exp(-(vr/c)^kE))/(vi - vr) - (wr*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE))/(vi - vr) - (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/c + (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE)/vi))/10;
 
                
            else fWeibull == 4
                
                %% CUSTO DE PENALIDADE
                
                % Weibull f ordem 2
                fw_p = exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - exp(-(vr/c)^kE)
                
                Pfw_p = (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*wr) - (kE*wr*exp(-(vi/c)^kE)*(vi/c)^kE)/2 + (kE*vr*wr*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi) - (kE*pg^2*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi*wr)
                
                
                pgfw_p = pg*fw_p
                                
                %propP = abs(Pfw_p - pgfw_p)
                propP = Pfw_p - pgfw_p
                
                Cp = (Kp*propP)/10
            
                
                %% CUSTO DE RESERVA
                
                %Weibull f ordem 2
                fw_r = exp(-(vi/c)^kE) - exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)
                
                Pfw_r = (kE*pg^2*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*vi*wr) - (kE*pg^2*exp(-(vi/c)^kE)*(vi/c)^kE)/(2*wr)
                
                
                pgfw_r = pg*fw_r
                
                %propR = abs(pgfw_r - Pfw_r)
                propR = pgfw_r - Pfw_r
                
                Cr = (Kr*propR)/100               
                
            
                %% GRADIENTE DERIVADA PRIMEIRA
            
                
                GF(2*nb+ntap+i,1) = (d - Kr*(exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - exp(-(vi/c)^kE) - (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE)/wr + (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/(c*wr) + (kE*pg*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(vi*wr)) - Kp*(exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE) - exp(-(vr/c)^kE) - (kE*pg*exp(-(vi/c)^kE)*(vi/c)^kE)/wr + (kE*pg*exp(-((pg*vr - pg*vi + vi*wr)/(c*wr))^kE)*(vi - vr)*((pg*vr - pg*vi + vi*wr)/(c*wr))^(kE - 1))/(c*wr) + (kE*pg*vr*exp(-(vi/c)^kE)*(vi/c)^kE)/(vi*wr)))/100;
 
                
                
                
                
            end
        
        
        %end   

end


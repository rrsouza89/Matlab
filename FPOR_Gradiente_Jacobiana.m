%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                        GRADIENTES E JACOBIANAS                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,fe,fed,fp,GF,G,H,JG,JH,JGT,JHT,seq,seqh,seqht,eta,lambda,z,Z,Zinv,L1,delta,g,vi,vr,v0,d,Cd,Cp,Cr,Kp,Kr,pg_eolica,Wr] = FPOR_Gradiente_Jacobiana(mi,b,r,g,Ger,NIN,nb,nr,ntap,ng,ngt,conj_ramo_tap,PQ,PV,Sl,x0,z,eta,lambda,delta,U1,U2,U3,U4,tolz,c,kE,L,fWeibull,it,omegaR,omegaP,PkmMin,PkmMax)

nSl = length(Sl);
nPV = length(PV);
nPQ = length(PQ);

p = (nSl+nPV);
m = 2*nPQ+p;
m1 = nSl+nPQ+nPV;                   %m1 = m"um"
                                                                                                                                        %PONTO DE VÁLVULA
G = zeros(m, 1);                    % Gradiente das restrições de igualdade
H = zeros(2*p+2*ngt+2*(2*nb+ntap+ng), 1);                    % Gradiente das restrições de desigualdade
GF = zeros(2*nb+ntap+ng+ngt,1);         % Gradiente da função objetivo

JG = zeros(m,2*nb+ntap+ng+ngt);            %Jacobiana das restrições de igualdade g(x)
JH = zeros(2*p+2*ngt+2*(2*nb+ntap+ng),2*nb+ntap+ng+ngt);            %Jacobiana das restrições de desigualdade h(x)


% G = zeros(m, 1);                    % Gradiente das restrições de igualdade
% H = zeros(2*p+2*(2*nb+ntap+ng), 1);                    % Gradiente das restrições de desigualdade
% GF = zeros(2*nb+ntap+ng,1);         % Gradiente da função objetivo
% 
% JG = zeros(m,2*nb+ntap+ng);            %Jacobiana das restrições de igualdade g(x)
% JH = zeros(2*p+2*(2*nb+ntap+ng),2*nb+ntap+ng);            %Jacobiana das restrições de desigualdade h(x)



cont_DP = 0;
cont_DQ = 0;
cont_H = 0;

seq = [];
seqh = [];
seqht = [];                                     % sequência dos geradores termelétricos das restrições de desigualdade

%% MATRIZES JACOBIANAS

for i=1:nb
    k = NIN(b(i).NUM);
    nviz = length(b(k).TVIZ);
    
    %% BARRA DE CARGA (PQ)
    if (b(k).TIPO==1)   
        
        seq = [seq,k];
        
        cont_DP = cont_DP + 1;
        cont_DQ = cont_DQ + 1;
        
        %Derivadas Primeira em relação a restrição de igualdade
        dpdvk = 0;
        dpdvm = 0;
        dpdtk = 0;
        dpdtm = 0;
            
        dqdvk = 0;
        dqdvm = 0;
        dqdtk = 0;
        dqdtm = 0;
                
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
    
        dqdvk = -2*vk*bksh;
            
        for j=1:nviz
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
                dpdvk = dpdvk + (2*gkm*vk)/akm^2 - (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaQ
                dqdvk = dqdvk + (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm - 2*vk*(bkmsh + bkm/akm^2);
                     
            else
                
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
                dpdvk = dpdvk + 2*gkm*vk - (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaQ
                dqdvk = dqdvk + (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm - 2*vk*(bkm + bkmsh);
                                   
            end
            
            % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
            dpdvm = -(vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtk = dpdtk -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtm = (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
                        
            % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaQ
            dqdvm = (vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dqdtk = dqdtk -(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dqdtm = (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            
            % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaP)
            JG(cont_DP,m) = JG(cont_DP,m) + dpdvm;
            JG(cont_DP,m+nb) = JG(cont_DP,m+nb) + dpdtm;
            
            % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaQ)
            JG(cont_DQ+m1,m) = JG(cont_DQ+m1,m) + dqdvm;
            JG(cont_DQ+m1,m+nb) = JG(cont_DQ+m1,m+nb) + dqdtm;
                        
        end %for vizinhas
        
        % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaP)
        JG(cont_DP,k) = JG(cont_DP,k) + dpdvk;
        JG(cont_DP,k+nb) = JG(cont_DP,k+nb) + dpdtk;
        
        % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaQ)
        JG(cont_DQ+m1,k) = JG(cont_DQ+m1,k) + dqdvk;
        JG(cont_DQ+m1,k+nb) = JG(cont_DQ+m1,k+nb) + dqdtk;
                       
    end %if BARRA PQ
    
    %% BARRA DE GERAÇÃO (PV)
    if (b(k).TIPO==2)  
        
        seq = [seq,k];
        seqh = [seqh,k];
        
        cont_DP = cont_DP + 1;
        cont_H = cont_H + 1;
                
        %Derivadas Primeira em relação a restrição de igualdade
        dpdvk = 0;
        dpdvm = 0;
        dpdtk = 0;
        dpdtm = 0;
                
        %Derivadas Primeira em relação a restrição de desigualdade
        dhdvk = 0;
        dhdvm = 0;
        dhdtk = 0;
        dhdtm = 0;
                              
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
    
        dhdvk = -2*vk*bksh;
            
        for j=1:nviz
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
                dpdvk = dpdvk + (2*gkm*vk)/akm^2 - (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                
                % Cálculo da Matriz Jacobiana da restrição de desigualdade
                dhdvk = dhdvk + (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm - 2*vk*(bkmsh + bkm/akm^2);
                               
            else
                
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
                dpdvk = dpdvk + 2*gkm*vk - (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                
                % Cálculo da Matriz Jacobiana da restrição de desigualdade
                dhdvk = dhdvk + (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm - 2*vk*(bkm + bkmsh);
                                   
            end
            
            % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
            dpdvm = -(vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtk = dpdtk -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtm = (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
                        
            % Cálculo da Matriz Jacobiana da restrição de desigualdade
            dhdvm = (vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dhdtk = dhdtk -(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtm = (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                        
            % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaP)
            JG(cont_DP,m) = JG(cont_DP,m) + dpdvm;
            JG(cont_DP,m+nb) = JG(cont_DP,m+nb) + dpdtm;
            
            % Armazenamento das variáveis na Matriz Jacobiana JH 
            % Qmin - Qg < 0
            JH(cont_H,m) = JH(cont_H,m) + dhdvm*(-1);            
            JH(cont_H,m+nb) = JH(cont_H,m+nb) + dhdtm*(-1);
            % Qg - Qmax < 0
            JH(cont_H+p,m) = JH(cont_H+p,m) + dhdvm;            
            JH(cont_H+p,m+nb) = JH(cont_H+p,m+nb) + dhdtm;
                                                
        end %for vizinhas
        
        % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaP)
        JG(cont_DP,k) = JG(cont_DP,k) + dpdvk;
        JG(cont_DP,k+nb) = JG(cont_DP,k+nb) + dpdtk;
        
        % Armazenamento das variáveis na Matriz Jacobiana JH 
        % Qmin - Qg < 0
        JH(cont_H,k) = JH(cont_H,k) + dhdvk*(-1);
        JH(cont_H,k+nb) = JH(cont_H,k+nb) + dhdtk*(-1);
        % Qg - Qmax < 0
        JH(cont_H+p,k) = JH(cont_H+p,k) + dhdvk;
        JH(cont_H+p,k+nb) = JH(cont_H+p,k+nb) + dhdtk;
                       
    end %if BARRA PV
    
    %% BARRA SLACK (Sl)
    if (b(k).TIPO==3)   
        
        seq = [seq,k];
        seqh = [seqh,k];
        
        cont_DP = cont_DP + 1;
        cont_H = cont_H + 1;
        
        %Derivadas Primeira em relação a restrição de igualdade
        dpdvk = 0;
        dpdvm = 0;
        dpdtk = 0;
        dpdtm = 0;
                                
        %Derivadas Primeira em relação a restrição de desigualdade
        dhdvk = 0;
        dhdvm = 0;
        dhdtk = 0;
        dhdtm = 0;
        dhdakm = 0;
                      
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
    
        dhdvk = -2*vk*bksh;
            
        for j=1:nviz
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
                dpdvk = dpdvk + (2*gkm*vk)/akm^2 - (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                
                % Cálculo da Matriz Jacobiana da restrição de desigualdade
                dhdvk = dhdvk + (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm - 2*vk*(bkmsh + bkm/akm^2);
                               
            else
                
                % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
                dpdvk = dpdvk + 2*gkm*vk - (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
                
                % Cálculo da Matriz Jacobiana da restrição de desigualdade
                dhdvk = dhdvk + (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm - 2*vk*(bkm + bkmsh);
                                    
            end
            
            % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
            dpdvm = -(vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtk = dpdtk -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtm = (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            
            % Cálculo da Matriz Jacobiana da restrição de desigualdade
            dhdvm = (vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dhdtk = dhdtk -(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtm = (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            
            % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaP)
            JG(cont_DP,m) = JG(cont_DP,m) + dpdvm;
            JG(cont_DP,m+nb) = JG(cont_DP,m+nb) + dpdtm;
                                               
            % Armazenamento das variáveis na Matriz Jacobiana JH 
            % Qmin - Qg < 0
            JH(cont_H,m) = JH(cont_H,m) + dhdvm*(-1);            
            JH(cont_H,m+nb) = JH(cont_H,m+nb) + dhdtm*(-1);
            % Qg - Qmax < 0
            JH(cont_H+p,m) = JH(cont_H+p,m) + dhdvm;            
            JH(cont_H+p,m+nb) = JH(cont_H+p,m+nb) + dhdtm;
            
            
                                                
        end %for vizinhas
        
        % Armazenamento das variáveis na Matriz Jacobiana JG (DeltaP)
        JG(cont_DP,k) = JG(cont_DP,k) + dpdvk;
        JG(cont_DP,k+nb) = JG(cont_DP,k+nb) + dpdtk;
        
        % Armazenamento das variáveis na Matriz Jacobiana JH
        % Qmin - Qg < 0
        JH(cont_H,k) = JH(cont_H,k) + dhdvk*(-1);
        JH(cont_H,k+nb) = JH(cont_H,k+nb) + dhdtk*(-1);
        % Qg - Qmax < 0
        JH(cont_H+p,k) = JH(cont_H+p,k) + dhdvk;
        JH(cont_H+p,k+nb) = JH(cont_H+p,k+nb) + dhdtk;
                                      
    end %if BARRA Slack
      
end %for nº da barra

%% RESTRIÇÃO PONTO DE VÁLVULA DESIGUALDADE

    pgmint = [];                                % Limites mínimos de potência ativa dos geradores termelétricos
    vt = [];
    
    for i=1:ng
   
        if b(seqh(i)).TGER == 1
            seqht = [seqht, seqh(i)];
            pgmint = [pgmint, U3(2*nb+ntap+i,1)];
            vt = [vt, i];
        end    
    end
    
    
%% GRADIENTE DAS RESTRIÇÕES DE DESIGUALDADE

for i=1:ngt
    
    pg = b(seqht(i)).PG;                     % Ponto de Válvula
    PGMIN = pgmint(i);
    F = g(vt(i)).F;
    
    if it == 1
        
        b(seqht(i)).m = -4*sin(F*(pg - PGMIN));        % Estimativa para a inicialização da variável auxiliar mpv relacionada ao ponto de válvula.
        
    end
    
    mpv = b(seqht(i)).m;
    
    H(2*ng+2*(2*nb+ntap+ng)+i,1) = sin(F*(PGMIN - pg)) - mpv;         % -fpv - mpv < 0 
    
    H(2*ng+2*(2*nb+ntap+ng)+ngt+i,1) = - mpv - sin(F*(PGMIN - pg));   %  fpv - mpv < 0 
            
end

%% JACOBIANA DAS RESTRIÇÕES DE DESIGUALDADE

nseqht = length(seqht);

    for i=1:ng
                  
        F = g(i).F;       
        
        for a=1:nseqht
                            
            if (Ger(i)==seqht(a))
                
                pg = b(seqht(a)).PG;                     % Ponto de Válvula
                PGMIN = pgmint(a);
                
                if (b(Ger(i)).TGER ==1)
                
                    JH(2*ng+2*(2*nb+ntap+ng)+a,2*nb+ntap+i) = -F*cos(F*(PGMIN - pg));               % Derivada em relação a PG  (-h - mpv < 0)
                    JH(2*ng+2*(2*nb+ntap+ng)+ngt+a,2*nb+ntap+i) = F*cos(F*(PGMIN - pg));          % (h - mpv < 0)
                    
                    JH(2*ng+2*(2*nb+ntap+ng)+a,2*nb+ntap+ng+a) = -1;                                 % Derivada em relação a m (-m - fpv)
                    JH(2*ng+2*(2*nb+ntap+ng)+ngt+a,2*nb+ntap+ng+a) = -1;                             % -m + fpv
                    
                end
               
            end     
        end
    end



%% FUNÇÃO OBJETIVO

% MINIMIZAÇÃO DAS PERDAS
fp = 0;

for i=1:nr
    k = r(i).NOI;
    m = r(i).NOF;
    vk = b(k).v;
    vm = b(m).v;
    tk = b(k).teta;  
    tm = b(m).teta;
    akm = r(i).TAP;
    gkm = r(i).GKM;
    
    fp = fp + gkm*(vm^2 + vk^2/akm^2 - (2*vk*vm*cos(tk - tm))/akm);
          
end

%fprintf('Perdas do sistema: %f\n\n',fp*100);


% CUSTOS DE GERAÇÃO
f = 0;
fe = 0;
fed = 0;
conta = 0;

% 30 barras       Kr = 5    Kp = 2
%xg = [1.7064; 0.4713; 0.2082; 0.1727; 0.1046; 0.1200; 0.1319];                      % PDPIEBLM
%xg = [1.557057; 0.435298; 0.196225; 0.1000; 0.1000; 0.1200; 0.40]                %    CS'

% 30 barras       Kr = 10   Kp = 2
%xg = [1.638355; 0.455225; 0.204511; 0.133914; 0.100028; 0.120000; 0.259991]                %    CS'


% 57 barras
% Kp = 2 Kr = 5
%xg = [2.184315; 1.496676; 1.410729; 1.199792; 2.400339; 1.199230; 2.428222; 0.399909]

% Kp = 2 Kr = 10
%xg = [2.174526; 1.499464; 1.412100; 1.195485; 2.427586; 1.197737; 2.436695; 0.364582]

for i=1:ng
    
    pg = b(Ger(i)).PG;
    mpv = b(Ger(i)).m;                                                   %PONTO DE VÁLVULA  
    
    %pg = xg(i);                                                         %TESTE CS
        
    if b(Ger(i)).TGER == 1
        
        conta = conta + 1;
                        
        A = g(i).A;
        B = g(i).B;
        C = g(i).C;
        E = g(i).E;
        
        g(i).Cter = (A*pg^2)/10000 + (B*pg)/10000 + C/10000;
        
        g(i).Cpv = (E*(abs(sin(F*(PGMIN - pg)))))/10000;
    
        g(i).Ct = g(i).Cter + g(i).Cpv; 
            
        f = f + g(i).Ct;                                          %PONTO DE VÁLVULA       
                
        %f = f + g(i).Cter;
    
        GF(2*nb+ntap+i,1) = B/10000 + (pg*A)/5000;                                  %PONTO DE VÁLVULA 
        GF(2*nb+ntap+ng+conta,1) = E/10000;
        
    elseif b(Ger(i)).TGER == 2
        
        [Cd,Cp,Cr,Kp,Kr,GF,vi,vr,v0,d,g] = FPOR_Custos_Eolica_TFC(g,c,kE,pg,i,fWeibull,L,GF,nb,ntap,omegaR,omegaP);
        
        pg_eolica = pg;   
        Wr = g(i).Wr;
        
        g(i).Ct = (Cd + Cp + Cr);
                          
        fed = fed + g(i).Cd;
        fe = fe + (g(i).Cp)/omegaP + (g(i).Cr)/omegaR; 
         
        g(i).Cd = Cd;
        g(i).Cp = Cp;
        g(i).Cr = Cr;

       

    end
        
    
end 

if ng == ngt
    vi = 0; vr = 0; v0 = 0; d = 0; Cd = 0; Cp = 0; Cr = 0; Kp = 0; Kr = 0; pg_eolica = 0; Wr = 0;
end

%fprintf('Custo da energia eólica no sistema: %f\n\n',fe*10000);

fprintf('Função Objetivo no ponto: %f\n\n',(f+fed+fe)*10000);



%% CÁLCULOS EM RELAÇÃO AOS TAPS

for i=1:ntap
    j = conj_ramo_tap(i);
    
    k = (r(j).NOI);
    vk = b(k).v;
    tk = b(k).teta;
    
    m = (r(j).NOF);
    vm = b(m).v;
    tm = b(m).teta;
    
    akm = r(j).TAP;
    gkm = r(j).GKM;
    bkm = r(j).BKM;
    
    % Derivadas Primeira (K - M) 
    dpdakm = 0;
    dqdakm = 0;
    dhdakm = 0;
    
    % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
    dpdakm = (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2 - (2*gkm*vk^2)/akm^3;
                    
    % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaQ
    dqdakm = (2*bkm*vk^2)/akm^3 - (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2;
    
    % Cálculo da Matriz Jacobiana da restrição de desigualdade
    dhdakm = (2*bkm*vk^2)/akm^3 - (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2;
           
    % Derivadas Primeira (M - K)
    dpdamk = 0;
    dqdamk = 0;
    dhdamk = 0;
    
    % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaP
    dpdamk = (vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2;
                    
    % Cálculo da Matriz Jacobiana da restrição de igualdade DeltaQ
    dqdamk = -(vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2;                
    
    % Cálculo da Matriz Jacobiana da restrição de desigualdade
    dhdamk = -(vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2;
    
    % Armazenamento na Matriz JG (deltaP)
    nseq = length(seq);
    for a=1:nseq
        if (k==seq(a))
            JG(a,2*nb+i) = JG(a,2*nb+i) + dpdakm;
        end
        if (m==seq(a))
            JG(a,2*nb+i) = JG(a,2*nb+i) + dpdamk;
        end
    end
    
    % Armazenamento na Matriz JG (deltaQ)
    for a=1:nPQ
        if (k==PQ(a))
            JG(a+nseq,2*nb+i) = JG(a+nseq,2*nb+i) + dqdakm;
        end
        if (m==PQ(a))
            JG(a+nseq,2*nb+i) = JG(a+nseq,2*nb+i) + dqdamk;
        end
    end
    
    % Armazenamento na Matriz JH
    nseqh = length(seqh);
    for a=1:nseqh
        if (k==seqh(a))
            JH(a,2*nb+i) = JH(a,2*nb+i) + dhdakm*(-1);
            JH(a+p,2*nb+i) = JH(a,2*nb+i) + dhdakm;
        end
        if (m==seqh(a))
            JH(a,2*nb+i) = JH(a,2*nb+i) + dhdamk*(-1);
            JH(a+p,2*nb+i) = JH(a,2*nb+i) + dhdamk;
        end
    end
    
    % Armazenamento no gradiente da função objetivo
    %GF(2*nb+i,1) = GF(2*nb+i,1) -gkm*((2*vk^2)/akm^3 - (2*vk*vm*cos(tk - tm))/akm^2);
                    
end %for TAP

%% CÁLCULO EM RELAÇÃO AOS GERADORES

%% JACOBIANA DAS RESTRIÇÕES DE IGUALDADE

    nseq = length(seq);
    
    for i=1:ng
    
        for a=1:nseq
            if (Ger(i)==seq(a))
                if (b(Ger(i)).TGER == 1)
                    JG(a,2*nb+ntap+i) = -1;
                 else
                    JG(a,2*nb+ntap+i) = -1; %% Caso seja um gerador eólico
                end
            end
        end
    end
    
 
%% CÁLCULO DAS RESTRIÇÕES

% DELTAP (G)
dp = 0;
nseq = length(seq);

for i=1:nseq
    nviz = length(b(seq(i)).TVIZ);
    vk = b(seq(i)).v;
    tk = b(seq(i)).teta;
        
    for j=1:nviz
        vm = b(b(seq(i)).VIZ(j)).v;
        tm = b(b(seq(i)).VIZ(j)).teta;
        bkm = r(b(seq(i)).RAMO(j)).BKM;
        gkm = r(b(seq(i)).RAMO(j)).GKM;
        akm = r(b(seq(i)).RAMO(j)).TAP;
                
        if b(seq(i)).TVIZ(j) == 0
            dp = dp + gkm*(vk/akm)^2-((vk*vm)/akm)*(gkm*cos(tk-tm)+bkm*sin(tk-tm));            
        else
            dp = dp + gkm*(vk^2)-((vk*vm)/akm)*(gkm*cos(tk-tm)+bkm*sin(tk-tm));            
        end
    end
    
    G(i,1) = dp - b(seq(i)).PG + b(seq(i)).PD;
    
    %%%%%%% INCLUIR RESTRIÇÃO DO FLUXO
    % -Pkm + PkmMin <= 0
    %  Pkm - PkmMax <= 0
    
    %H(2*ng+2*ngt+2*(2*nb+ntap+ng) + i, 1) = -dp + PkmMin;
    %H(2*ng+2*ngt+2*(2*nb+ntap+ng)+ nb + i, 1) = dp - PkmMax;
          
    dp = 0;
    
end

% DELTAQ (G)
dq = 0;

for i=1:nPQ
    nviz = length(b(PQ(i)).TVIZ);
    vk = b(PQ(i)).v;
    tk = b(PQ(i)).teta;
    bksh = b(PQ(i)).BS;
    
    for j=1:nviz
        vm = b(b(PQ(i)).VIZ(j)).v;
        tm = b(b(PQ(i)).VIZ(j)).teta;
        bkm = r(b(PQ(i)).RAMO(j)).BKM;
        gkm = r(b(PQ(i)).RAMO(j)).GKM;
        akm = r(b(PQ(i)).RAMO(j)).TAP;
        bkmsh = r(b(PQ(i)).RAMO(j)).B;
        
        if b(PQ(i)).TVIZ(j) == 0
            dq = dq -((bkm/(akm^2))+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));            
        else
            dq = dq -(bkm+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));            
        end
    end
    
    G(i+m1,1) = dq - b(PQ(i)).QG + b(PQ(i)).QD - ((vk^2)*bksh);
    dq = 0;
end

% (H)

nseqh = length(seqh);
somah = 0;

for i=1:nseqh
    nviz = length(b(seqh(i)).TVIZ);
    vk = b(seqh(i)).v;
    tk = b(seqh(i)).teta;
    bksh = b(seqh(i)).BS;
          
    for j=1:nviz
        vm = b(b(seqh(i)).VIZ(j)).v;
        tm = b(b(seqh(i)).VIZ(j)).teta;
        bkm = r(b(seqh(i)).RAMO(j)).BKM;
        gkm = r(b(seqh(i)).RAMO(j)).GKM;
        akm = r(b(seqh(i)).RAMO(j)).TAP;
        bkmsh = r(b(seqh(i)).RAMO(j)).B;        
        
        if b(seqh(i)).TVIZ(j) == 0
            somah = somah -((bkm/akm^2)+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));            
        else
            somah = somah -(bkm+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));            
        end
    end
    
    H(i,1) = somah + b(seqh(i)).QD - ((vk^2)*bksh);        
    somah = 0;
    
    H(i,1) = U1(i) - H(i,1);                            % Qmin - Qg < 0 
    H(i+ng,1) = H(i,1) - U2(i);                         % Qg - Qmax < 0
        
end

%% RESTRIÇÕES DE DESIGUALDADE VARIÁVEIS CANALIZADAS (V TETA TAP PG)

    for i=1:(2*nb+ntap+ng)
        H(2*ng+i,1) = -x0(i) + U3(i);
        H(2*ng+2*nb+ntap+ng+i,1) = x0(i) - U4(i);    
    end
    
%% JACOBIANA DAS RESTRIÇÕES DE DESIGUALDADE DAS VARIÁVEIS CANALIZADAS


    for i=1:(2*nb+ntap+ng)                          % V teta Tap Pg
        
        JH(2*ng+i,i) = -1;                          % - x + Umin < 0
        
        JH(2*ng+2*nb+ntap+ng+i,i) = 1;              % x - Umax < 0
                
    end
 


%% DADOS DE SAÍDA

nz = length(H);

if it == 1                                 % Apenas primeira iteração

    for i=1:nz
        z(i) = max(-H(i), 0.1);  
    end

    z = z';

    if z(1)>tolz
        delta = ones(nz,1);
    else    
        delta = 0.1*ones(nz,1);
    end
    
end

% if nb == 300
%    delta = ones(nz,1); 
%     
% end


for i=1:(nz)
    if z(i)>tolz
        Z(i,i) = z(i);
        Zinv(i,i) = 1/Z(i,i);
    else
        Z(i,i) = z(i)+mi;
        Zinv(i,i) = 1/Z(i,i);
    end    
end

JGT = JG';
JHT = JH';


if it == 1                                   % Apenas primeira iteração

    lambda = mi*Zinv*delta;  

    eta = -(JG*JGT)\(JG*(GF+JHT*lambda));
    
end

L1 = diag(lambda);

end %function

    
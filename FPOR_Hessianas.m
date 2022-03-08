%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                                 HESSIANAS                              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HessG, HessH, HessF, K, nodal] = FPOR_Hessianas(seq,seqh,seqht,b,r,g,NIN,nb,nr,ntap,ng,ngt,conj_ramo_tap,PQ,PV,Sl,eta,lambda,x0,c,Ger,L,kE,Kp,Kr,fWeibull,U3,nodal,omegaR,omegaP);

format long;

nSl = length(Sl);
nPV = length(PV);
nPQ = length(PQ);

m = 2*nPQ+nPV;
m1 = nSl+nPQ+nPV;
p = nSl+nPV;
                                                                                                                    % PONTO DE VÁLVULA
HessDP = zeros(2*nb+ntap+ng+ngt,2*nb+ntap+ng+ngt);    %Matriz Hessiana do DeltaP das restrições de igualdade
HessDQ = zeros(2*nb+ntap+ng+ngt,2*nb+ntap+ng+ngt);    %Matriz Hessiana do DeltaQ das restrições de igualdade

HessH = zeros(2*nb+ntap+ng+ngt,2*nb+ntap+ng+ngt);     %Matriz Hessiana das restrições de desigualdade

HessF = zeros(2*nb+ntap+ng+ngt,2*nb+ntap+ng+ngt);     %Matriz Hessiana da Função Objetivo


% HessDP = zeros(2*nb+ntap+ng,2*nb+ntap+ng);    %Matriz Hessiana do DeltaP das restrições de igualdade
% HessDQ = zeros(2*nb+ntap+ng,2*nb+ntap+ng);    %Matriz Hessiana do DeltaQ das restrições de igualdade
% 
% HessH = zeros(2*nb+ntap+ng,2*nb+ntap+ng);     %Matriz Hessiana das restrições de desigualdade
% 
% HessF = zeros(2*nb+ntap+ng,2*nb+ntap+ng);

nseq = length(seq);

%% HESSIANA RESTRIÇÕES DE IGUALDADE

npos = 0;

for i=1:nseq
    k = NIN(b(seq(i)).NUM);
    nviz = length(b(k).TVIZ);
    
    %% BARRA DE CARGA (PQ)
    if (b(k).TIPO==1) 
        
        npos = npos+1;
        
        %Derivadas Segunda em relação a restrições de igualdade
        %DeltaP
        dpdvkvk = 0;
        dpdvkvm = 0;
        
        dpdtkvk = 0;
        dpdtmvk = 0;
        dpdtkvm = 0;
        dpdtmvm = 0;
    
        dpdtktk = 0;
        dpdtktm = 0;
        dpdtmtm = 0;
    
        %DeltaQ
        dqdvkvk = 0;
        dqdvkvm = 0;
        
        dqdtkvk = 0;
        dqdtmvk = 0;
        dqdtkvm = 0;
        dqdtmvm = 0;
    
        dqdtktk = 0;
        dqdtktm = 0;
        dqdtmtm = 0;
        
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
        
        dqdvkvk = -2*bksh;
    
        for j=1:nviz
            
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz Hessiana DeltaP 
                dpdvkvk = dpdvkvk + (2*gkm)/akm^2;
            
                % Cálculo da Matriz Hessiana DeltaQ
                dqdvkvk = dqdvkvk - 2*bkmsh - (2*bkm)/akm^2;
                     
            else
                
                % Cálculo da Matriz Hessiana DeltaP 
                dpdvkvk = dpdvkvk + (2*gkm);
            
                % Cálculo da Matriz Hessiana DeltaQ
                dqdvkvk = dqdvkvk - 2*bkm - 2*bkmsh;
                                   
            end
            
            %Cálculo da Matriz Hessiana DeltaP
            dpdvkvm = -(gkm*cos(tk - tm) + bkm*sin(tk - tm))/akm;
            
            dpdtkvk = dpdtkvk -(vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtmvk = (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtkvm = -(vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtmvm = (vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            
            dpdtktk = dpdtktk + (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtktm = -(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtmtm = (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;   
                       
            %Cálculo da Matriz Hessiana DeltaQ
            dqdvkvm = (bkm*cos(tk - tm) - gkm*sin(tk - tm))/akm;
            
            dqdtkvk = dqdtkvk -(vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dqdtmvk = (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dqdtkvm = -(vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dqdtmvm = (vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            
            dqdtktk = dqdtktk -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dqdtktm = (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dqdtmtm = -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
                       
            %Armazenamento das variáveis na Matriz HessDP
            HessDP(k,m)= HessDP(k,m) + dpdvkvm*eta(i);
            HessDP(m,k)= HessDP(m,k) + dpdvkvm*eta(i);

            HessDP(k,m+nb)= HessDP(k,m+nb) + dpdtmvk*eta(i);
            HessDP(m,k+nb)= HessDP(m,k+nb) + dpdtkvm*eta(i);
            HessDP(m,m+nb)= HessDP(m,m+nb) + dpdtmvm*eta(i);

            HessDP(k+nb,m)= HessDP(k+nb,m) + dpdtkvm*eta(i);
            HessDP(m+nb,k)= HessDP(m+nb,k) + dpdtmvk*eta(i);
            HessDP(m+nb,m)= HessDP(m+nb,m) + dpdtmvm*eta(i);
            
            HessDP(k+nb,m+nb)= HessDP(k+nb,m+nb) + dpdtktm*eta(i);
            HessDP(m+nb,k+nb)= HessDP(m+nb,k+nb) + dpdtktm*eta(i);
            HessDP(m+nb,m+nb)= HessDP(m+nb,m+nb) + dpdtmtm*eta(i);

            %Armazenamento das variáveis na Matriz HessDQ
            HessDQ(k,m)= HessDQ(k,m) + dqdvkvm*eta(npos+m1);
            HessDQ(m,k)= HessDQ(m,k) + dqdvkvm*eta(npos+m1);

            HessDQ(k,m+nb)= HessDQ(k,m+nb) + dqdtmvk*eta(npos+m1);
            HessDQ(m,k+nb)= HessDQ(m,k+nb) + dqdtkvm*eta(npos+m1);
            HessDQ(m,m+nb)= HessDQ(m,m+nb) + dqdtmvm*eta(npos+m1);

            HessDQ(k+nb,m)= HessDQ(k+nb,m) + dqdtkvm*eta(npos+m1);
            HessDQ(m+nb,k)= HessDQ(m+nb,k) + dqdtmvk*eta(npos+m1);
            HessDQ(m+nb,m)= HessDQ(m+nb,m) + dqdtmvm*eta(npos+m1);

            HessDQ(k+nb,m+nb)= HessDQ(k+nb,m+nb) + dqdtktm*eta(npos+m1);
            HessDQ(m+nb,k+nb)= HessDQ(m+nb,k+nb) + dqdtktm*eta(npos+m1);
            HessDQ(m+nb,m+nb)= HessDQ(m+nb,m+nb) + dqdtmtm*eta(npos+m1);
                        
        end %for vizinhas
        
        %Armazenamento das variáveis na Matriz HessDP
        HessDP(k,k)= HessDP(k,k) + dpdvkvk*eta(i);
        HessDP(k,k+nb)= HessDP(k,k+nb) + dpdtkvk*eta(i);
        HessDP(k+nb,k)= HessDP(k+nb,k) + dpdtkvk*eta(i);
        HessDP(k+nb,k+nb)= HessDP(k+nb,k+nb) + dpdtktk*eta(i);
        
        %Armazenamento das variáveis na Matriz HessDQ
        HessDQ(k,k)= HessDQ(k,k) + dqdvkvk*eta(npos+m1);
        HessDQ(k,k+nb)= HessDQ(k,k+nb) + dqdtkvk*eta(npos+m1);
        HessDQ(k+nb,k)= HessDQ(k+nb,k) + dqdtkvk*eta(npos+m1);
        HessDQ(k+nb,k+nb)= HessDQ(k+nb,k+nb) + dqdtktk*eta(npos+m1);  
        
                       
    end %if BARRA PQ
    
    %% BARRA DE GERAÇÃO (PV) E SLACK
    if (b(k).TIPO==2)||(b(k).TIPO==3)  
        
        %Derivadas Segunda em relação a restrições de igualdade
        dpdvkvk = 0;
        dpdvkvm = 0;

        dpdtkvk = 0;
        dpdtmvk = 0;
        dpdtkvm = 0;
        dpdtmvm = 0;

        dpdtktk = 0;
        dpdtktm = 0;
        dpdtmtm = 0;
                      
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
        
        for j=1:nviz
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz Hessiana DeltaP 
                dpdvkvk = dpdvkvk + (2*gkm)/akm^2;
                                                 
            else
                
                % Cálculo da Matriz Hessiana DeltaP 
                dpdvkvk = dpdvkvk + (2*gkm);
                                              
            end
            
            %Cálculo da Matriz Hessiana DeltaP
            dpdtmvk = (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtmvm = (vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtktm = -(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtmtm = (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdtkvk = dpdtkvk -(vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dpdtktk = dpdtktk + (vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dpdvkvm = -(gkm*cos(tk - tm) + bkm*sin(tk - tm))/akm;
            dpdtkvm = -(vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
                                    
            %Armazenamento das variáveis na Matriz HessDP
            HessDP(k,m)= HessDP(k,m) + dpdvkvm*eta(i);
            HessDP(m,k)= HessDP(m,k) + dpdvkvm*eta(i);

            HessDP(k,m+nb)= HessDP(k,m+nb) + dpdtmvk*eta(i);
            HessDP(m,k+nb)= HessDP(m,k+nb) + dpdtkvm*eta(i);
            HessDP(m,m+nb)= HessDP(m,m+nb) + dpdtmvm*eta(i);

            HessDP(k+nb,m)= HessDP(k+nb,m) + dpdtkvm*eta(i);
            HessDP(m+nb,k)= HessDP(m+nb,k) + dpdtmvk*eta(i);
            HessDP(m+nb,m)= HessDP(m+nb,m) + dpdtmvm*eta(i);

            HessDP(k+nb,m+nb)= HessDP(k+nb,m+nb) + dpdtktm*eta(i);
            HessDP(m+nb,k+nb)= HessDP(m+nb,k+nb) + dpdtktm*eta(i);
            HessDP(m+nb,m+nb)= HessDP(m+nb,m+nb) + dpdtmtm*eta(i);  
                                                
        end %for vizinhas
        
        %Armazenamento das variáveis na Matriz HessDP
        HessDP(k,k)= HessDP(k,k) + dpdvkvk*eta(i);
        HessDP(k,k+nb)= HessDP(k,k+nb) + dpdtkvk*eta(i);
        HessDP(k+nb,k)= HessDP(k+nb,k) + dpdtkvk*eta(i);
        HessDP(k+nb,k+nb)= HessDP(k+nb,k+nb) + dpdtktk*eta(i);
        
        nodal = [nodal; eta(i)*100];
                       
    end %if BARRA PV
    
end %for nº da barra

%% HESSIANA RESTRIÇÕES DE DESIGUALDADE

nseqh = length(seqh);

for i=1:nseqh
    k = NIN(b(seqh(i)).NUM);
    nviz = length(b(k).TVIZ);

    %% BARRA SLACK (Sl)
    if (b(k).TIPO==3)   
                                
        %Derivadas Segunda em relação a restrição de desigualdade
        dhdvkvk = 0;
        dhdvkvm = 0;

        dhdtkvk = 0;
        dhdtmvk = 0;
        dhdtkvm = 0;
        dhdtmvm = 0;

        dhdtktk = 0;
        dhdtktm = 0;
        dhdtmtm = 0;
                      
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
        
        dhdvkvk = -2*bksh;
    
        for j=1:nviz
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz HessH
                dhdvkvk = dhdvkvk - 2*bkmsh - (2*bkm)/akm^2;
                               
            else
                
                % Cálculo da Matriz HessH
                dhdvkvk = dhdvkvk - 2*bkm - 2*bkmsh;
                                    
            end
            
            %Cálculo da Matriz HessH
            dhdtkvk = dhdtkvk -(vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtktk = dhdtktk -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dhdvkvm = (bkm*cos(tk - tm) - gkm*sin(tk - tm))/akm;
            dhdtmvk = (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtkvm = -(vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtmvm = (vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtktm = (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dhdtmtm = -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            
                   
            %Armazenamento das variáveis na Matriz HessH                            
            HessH(k,m)= HessH(k,m) + dhdvkvm*(lambda(ng+i)-lambda(i));
            HessH(m,k)= HessH(m,k) + dhdvkvm*(lambda(ng+i)-lambda(i));
        
            HessH(k,m+nb)= HessH(k,m+nb) + dhdtmvk*(lambda(ng+i)-lambda(i));
            HessH(m,k+nb)= HessH(m,k+nb) + dhdtkvm*(lambda(ng+i)-lambda(i));
            HessH(m,m+nb)= HessH(m,m+nb) + dhdtmvm*(lambda(ng+i)-lambda(i));
        
            HessH(k+nb,m)= HessH(k+nb,m) + dhdtkvm*(lambda(ng+i)-lambda(i));
            HessH(m+nb,k)= HessH(m+nb,k) + dhdtmvk*(lambda(ng+i)-lambda(i));
            HessH(m+nb,m)= HessH(m+nb,m) + dhdtmvm*(lambda(ng+i)-lambda(i));
        
            HessH(k+nb,m+nb)= HessH(k+nb,m+nb) + dhdtktm*(lambda(ng+i)-lambda(i));
            HessH(m+nb,k+nb)= HessH(m+nb,k+nb) + dhdtktm*(lambda(ng+i)-lambda(i));
            HessH(m+nb,m+nb)= HessH(m+nb,m+nb) + dhdtmtm*(lambda(ng+i)-lambda(i));
                                                
        end %for vizinhas
        
        %Armazenamento das variáveis na Matriz HessH
        HessH(k,k)= HessH(k,k) + dhdvkvk*(lambda(ng+i)-lambda(i));
        HessH(k,k+nb)= HessH(k,k+nb) + dhdtkvk*(lambda(ng+i)-lambda(i));
        HessH(k+nb,k)= HessH(k+nb,k) + dhdtkvk*(lambda(ng+i)-lambda(i));
        HessH(k+nb,k+nb)= HessH(k+nb,k+nb) + dhdtktk*(lambda(ng+i)-lambda(i));
                       
    end %if BARRA Slack
    
    %% BARRA DE GERAÇÃO (PV)
    if (b(k).TIPO==2)  
        
        %Derivadas Segunda em relação a restrição de desigualdade
        dhdvkvk = 0;
        dhdvkvm = 0;

        dhdtkvk = 0;
        dhdtmvk = 0;
        dhdtkvm = 0;
        dhdtmvm = 0;

        dhdtktk = 0;
        dhdtktm = 0;
        dhdtmtm = 0;
                      
        vk = b(k).v;
        tk = b(k).teta;
        bksh = b(k).BS;
        
        dhdvkvk = -2*bksh;
    
        for j=1:nviz
            m = b(k).VIZ(j);
            vm = b(m).v;
            tm = b(m).teta;
        
            akm = r(b(k).RAMO(j)).TAP;
            gkm = r(b(k).RAMO(j)).GKM;
            bkm = r(b(k).RAMO(j)).BKM;
            bkmsh = r(b(k).RAMO(j)).B;
    
            if (b(k).TVIZ(j)==0) % Barra k é nó inicial
                
                % Cálculo da Matriz HessH
                dhdvkvk = dhdvkvk - 2*bkmsh - (2*bkm)/akm^2;
                               
            else
                
                % Cálculo da Matriz HessH
                dhdvkvk = dhdvkvk - 2*bkm - 2*bkmsh;
                                    
            end
            
            %Cálculo da Matriz HessH
            dhdtkvk = dhdtkvk -(vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtktk = dhdtktk -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dhdvkvm = (bkm*cos(tk - tm) - gkm*sin(tk - tm))/akm;
            dhdtmvk = (vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtkvm = -(vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtmvm = (vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm;
            dhdtktm = (vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            dhdtmtm = -(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm;
            
                        
            %Armazenamento das variáveis na Matriz HessH
            HessH(k,m)= HessH(k,m) + dhdvkvm*(lambda(ng+i)-lambda(i));
            HessH(m,k)= HessH(m,k) + dhdvkvm*(lambda(ng+i)-lambda(i));

            HessH(k,m+nb)= HessH(k,m+nb) + dhdtmvk*(lambda(ng+i)-lambda(i));
            HessH(m,k+nb)= HessH(m,k+nb) + dhdtkvm*(lambda(ng+i)-lambda(i));
            HessH(m,m+nb)= HessH(m,m+nb) + dhdtmvm*(lambda(ng+i)-lambda(i));

            HessH(k+nb,m)= HessH(k+nb,m) + dhdtkvm*(lambda(ng+i)-lambda(i));
            HessH(m+nb,k)= HessH(m+nb,k) + dhdtmvk*(lambda(ng+i)-lambda(i));
            HessH(m+nb,m)= HessH(m+nb,m) + dhdtmvm*(lambda(ng+i)-lambda(i));

            HessH(k+nb,m+nb)= HessH(k+nb,m+nb) + dhdtktm*(lambda(ng+i)-lambda(i));
            HessH(m+nb,k+nb)= HessH(m+nb,k+nb) + dhdtktm*(lambda(ng+i)-lambda(i));
            HessH(m+nb,m+nb)= HessH(m+nb,m+nb) + dhdtmtm*(lambda(ng+i)-lambda(i));
                                                
        end %for vizinhas
        
        %Armazenamento das variáveis na Matriz HessH
        HessH(k,k)= HessH(k,k) + dhdvkvk*(lambda(ng+i)-lambda(i));
        HessH(k,k+nb)= HessH(k,k+nb) + dhdtkvk*(lambda(ng+i)-lambda(i));
        HessH(k+nb,k)= HessH(k+nb,k) + dhdtkvk*(lambda(ng+i)-lambda(i));
        HessH(k+nb,k+nb)= HessH(k+nb,k+nb) + dhdtktk*(lambda(ng+i)-lambda(i));
                       
    end %if BARRA PV
    
end % for nseqh

%% HESSIANA FUNÇÃO OBJETIVO

% MINIMIZAÇÃO DAS PERDAS
% for i=1:nr
%     
%     k = r(i).NOI;
%     m = r(i).NOF;
%     vk = b(k).v;
%     vm = b(m).v;
%     tk = b(k).teta;  
%     tm = b(m).teta;
%     akm = r(i).TAP;
%     gkm = r(i).GKM;
%     
%    % MATRIZ HESSIANA
%    
%    HessF(k,k) = HessF(k,k) + (2*gkm)/akm^2;
%    HessF(k,m) = HessF(k,m) - (2*gkm*cos(tk - tm))/akm;
%    HessF(m,k) = HessF(m,k) - (2*gkm*cos(tk - tm))/akm;
%    HessF(m,m) = HessF(m,m) + 2*gkm;
%    
%    HessF(k,k+nb) = HessF(k,k+nb) + (2*gkm*vm*sin(tk - tm))/akm;
%    HessF(m,k+nb) = HessF(m,k+nb) + (2*gkm*vk*sin(tk - tm))/akm;
%    HessF(m,m+nb) = HessF(m,m+nb) - (2*gkm*vk*sin(tk - tm))/akm;
%    HessF(k,m+nb) = HessF(k,m+nb) - (2*gkm*vm*sin(tk - tm))/akm;
%         
%    HessF(k+nb,k) = HessF(k+nb,k) + (2*gkm*vm*sin(tk - tm))/akm;
%    HessF(k+nb,m) = HessF(k+nb,m) + (2*gkm*vk*sin(tk - tm))/akm;
%    HessF(m+nb,m) = HessF(m+nb,m) - (2*gkm*vk*sin(tk - tm))/akm;
%    HessF(m+nb,k) = HessF(m+nb,k) - (2*gkm*vm*sin(tk - tm))/akm;
%         
%    HessF(k+nb,k+nb) = HessF(k+nb,k+nb) + (2*gkm*vk*vm*cos(tk - tm))/akm;
%    HessF(k+nb,m+nb) = HessF(k+nb,m+nb) - (2*gkm*vk*vm*cos(tk - tm))/akm;
%    HessF(m+nb,k+nb) = HessF(m+nb,k+nb) - (2*gkm*vk*vm*cos(tk - tm))/akm;
%    HessF(m+nb,m+nb) = HessF(m+nb,m+nb) + (2*gkm*vk*vm*cos(tk - tm))/akm;
%    
% end

% CUSTOS DE GERAÇÃO

for i=1:ng
    
    if b(Ger(i)).TGER == 1
    
        A = g(i).A;
        HessF(2*nb+ntap+i,2*nb+ntap+i) = A/5000;
        
    elseif b(Ger(i)).TGER == 2
        
        pg = b(Ger(i)).PG;
        
        vi = g(i).Vi;
        vr = g(i).Vr;
        v0 = g(i).V0;
        wr = g(i).Wr;
                                
                   
        %% MISHRA
        
        if pg < L
            
            HessF(2*nb+ntap+i,2*nb+ntap+i) = (omegaP*(Kp*(exp(-(v0/c)^kE) - exp(-(vi/c)^kE) + 1)))/10000 + (omegaR*(Kr*(exp(-(v0/c)^kE) - exp(-(vi/c)^kE) + 1)))/10000;
            
        elseif (wr-pg) < L
        
            HessF(2*nb+ntap+i,2*nb+ntap+i) = (omegaP*(- Kp*(exp(-(v0/c)^kE) - exp(-(vr/c)^kE))))/10000 + (omegaR*(- Kr*(exp(-(v0/c)^kE) - exp(-(vr/c)^kE))))/10000;
            
        else
        
            HessF(2*nb+ntap+i,2*nb+ntap+i) = (omegaP*(-(Kp*kE*exp(-((vi*wr - pg*(vi - vr))/(c*wr))^kE)*(vi - vr)*((vi*wr - pg*(vi - vr))/(c*wr))^(kE - 1))/(c*wr)))/10000 + (omegaR*(-(Kr*kE*exp(-((vi*wr - pg*(vi - vr))/(c*wr))^kE)*(vi - vr)*((vi*wr - pg*(vi - vr))/(c*wr))^(kE - 1))/(c*wr)))/10000;
                
         end   
              
       
    end
    
end
    




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
    bkmsh = r(j).B;
    
    
    %% HESSIANA RESTRIÇÕES DE IGUALDADE
    nseq = length(seq);
    
    for np=1:nseq        % DeltaP
        if (k==seq(np))  % k inicial
            
            HessDP(k,2*nb+i) = HessDP(k,2*nb+i) + ((vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2 - (4*gkm*vk)/akm^3)*eta(np);
            HessDP(m,2*nb+i) = HessDP(m,2*nb+i) + ((vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*eta(np);
            HessDP(2*nb+i,k) = HessDP(2*nb+i,k) + ((vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2 - (4*gkm*vk)/akm^3)*eta(np);
            HessDP(2*nb+i,m) = HessDP(2*nb+i,m) + ((vk*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*eta(np);

            HessDP(k+nb,2*nb+i) = HessDP(k+nb,2*nb+i) + ((vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(np);
            HessDP(m+nb,2*nb+i) = HessDP(m+nb,2*nb+i) + (-(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(np);
            HessDP(2*nb+i,k+nb) = HessDP(2*nb+i,k+nb) + ((vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(np);
            HessDP(2*nb+i,m+nb) = HessDP(2*nb+i,m+nb) + (-(vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(np);

            HessDP(2*nb+i,2*nb+i) = HessDP(2*nb+i,2*nb+i) + ((6*gkm*vk^2)/akm^4 - (2*vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^3)*eta(np);
           
        end
        
        if (m==seq(np))  % k final
            
            HessDP(k,2*nb+i) = HessDP(k,2*nb+i) + ((vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(np);
            HessDP(m,2*nb+i) = HessDP(m,2*nb+i) + ((vk*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(np);
            HessDP(2*nb+i,k) = HessDP(2*nb+i,k) + ((vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(np);
            HessDP(2*nb+i,m) = HessDP(2*nb+i,m) + ((vk*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(np);

            HessDP(k+nb,2*nb+i) = HessDP(k+nb,2*nb+i) + (-(vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(np);
            HessDP(m+nb,2*nb+i) = HessDP(m+nb,2*nb+i) + ((vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(np);
            HessDP(2*nb+i,k+nb) = HessDP(2*nb+i,k+nb) + (-(vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(np);
            HessDP(2*nb+i,m+nb) = HessDP(2*nb+i,m+nb) + ((vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(np);

            HessDP(2*nb+i,2*nb+i) = HessDP(2*nb+i,2*nb+i) + (-(2*vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^3)*eta(np);
           
        end
    end     % for nseq
       
    for q=1:nPQ        % DeltaQ
        if (k==PQ(q))  % k inicial
            
            HessDQ(k,2*nb+i) = HessDQ(k,2*nb+i) + ((4*bkm*vk)/akm^3 - (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(m1+q);
            HessDQ(m,2*nb+i) = HessDQ(m,2*nb+i) + (-(vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,k) = HessDQ(2*nb+i,k) + ((4*bkm*vk)/akm^3 - (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,m) = HessDQ(2*nb+i,m) + (-(vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*eta(m1+q);

            HessDQ(k+nb,2*nb+i) = HessDQ(k+nb,2*nb+i) + ((vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*eta(m1+q);
            HessDQ(m+nb,2*nb+i) = HessDQ(m+nb,2*nb+i) + (-(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,k+nb) = HessDQ(2*nb+i,k+nb) + ((vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,m+nb) = HessDQ(2*nb+i,m+nb) + (-(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*eta(m1+q);

            HessDQ(2*nb+i,2*nb+i) = HessDQ(2*nb+i,2*nb+i) + ((2*vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^3 - (6*bkm*vk^2)/akm^4)*eta(m1+q);
           
        end
        
        if (m==PQ(q))  % k final
            
            HessDQ(k,2*nb+i) = HessDQ(k,2*nb+i) + (-(vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(m1+q);
            HessDQ(m,2*nb+i) = HessDQ(m,2*nb+i) + (-(vk*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,k) = HessDQ(2*nb+i,k) + (-(vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,m) = HessDQ(2*nb+i,m) + (-(vk*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*eta(m1+q);

            HessDQ(k+nb,2*nb+i) = HessDQ(k+nb,2*nb+i) + (-(vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(m1+q);
            HessDQ(m+nb,2*nb+i) = HessDQ(m+nb,2*nb+i) + ((vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,k+nb) = HessDQ(2*nb+i,k+nb) + (-(vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(m1+q);
            HessDQ(2*nb+i,m+nb) = HessDQ(2*nb+i,m+nb) + ((vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*eta(m1+q);

            HessDQ(2*nb+i,2*nb+i) = HessDQ(2*nb+i,2*nb+i) + ((2*vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^3)*eta(m1+q);
           
        end
    end     % for nPQ
    
    %% HESSIANA RESTRIÇÕES DE DESIGUALDADE
    nseqh = length(seqh);
    
    for nh=1:nseqh        
        if (k==seqh(nh))  % k inicial
            
            HessH(k,2*nb+i) = HessH(k,2*nb+i) + ((4*bkm*vk)/akm^3 - (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(m,2*nb+i) = HessH(m,2*nb+i) + (-(vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,k) = HessH(2*nb+i,k) + ((4*bkm*vk)/akm^3 - (vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,m) = HessH(2*nb+i,m) + (-(vk*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));

            HessH(k+nb,2*nb+i) = HessH(k+nb,2*nb+i) + ((vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(m+nb,2*nb+i) = HessH(m+nb,2*nb+i) + (-(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,k+nb) = HessH(2*nb+i,k+nb) + ((vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,m+nb) = HessH(2*nb+i,m+nb) + (-(vk*vm*(gkm*cos(tk - tm) + bkm*sin(tk - tm)))/akm^2)*(lambda(ng+nh)-lambda(nh));

            HessH(2*nb+i,2*nb+i) = HessH(2*nb+i,2*nb+i) + ((2*vk*vm*(bkm*cos(tk - tm) - gkm*sin(tk - tm)))/akm^3 - (6*bkm*vk^2)/akm^4)*(lambda(ng+nh)-lambda(nh));
           
        end
        
        if (m==seqh(nh))  % k final
            
            HessH(k,2*nb+i) = HessH(k,2*nb+i) + (-(vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(m,2*nb+i) = HessH(m,2*nb+i) + (-(vk*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,k) = HessH(2*nb+i,k) + (-(vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,m) = HessH(2*nb+i,m) + (-(vk*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));

            HessH(k+nb,2*nb+i) = HessH(k+nb,2*nb+i) + (-(vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(m+nb,2*nb+i) = HessH(m+nb,2*nb+i) + ((vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,k+nb) = HessH(2*nb+i,k+nb) + (-(vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));
            HessH(2*nb+i,m+nb) = HessH(2*nb+i,m+nb) + ((vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2)*(lambda(ng+nh)-lambda(nh));

            HessH(2*nb+i,2*nb+i) = HessH(2*nb+i,2*nb+i) + ((2*vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^3)*(lambda(ng+nh)-lambda(nh));
           
        end
    end     % for nseqh
    
    %% HESSIANA FUNÇÃO OBJETIVO
    
%     HessF(k,2*nb+i) = HessF(k,2*nb+i)-gkm*((4*vk)/akm^3 - (2*vm*cos(tk - tm))/akm^2);
%     HessF(2*nb+i,k) = HessF(2*nb+i,k)-gkm*((4*vk)/akm^3 - (2*vm*cos(tk - tm))/akm^2);
%     HessF(m,2*nb+i) = HessF(m,2*nb+i)+(2*gkm*vk*cos(tk - tm))/akm^2;
%     HessF(2*nb+i,m) = HessF(2*nb+i,m)+(2*gkm*vk*cos(tk - tm))/akm^2;
%         
%     HessF(k+nb,2*nb+i) = HessF(k+nb,2*nb+i)-(2*gkm*vk*vm*sin(tk - tm))/akm^2;
%     HessF(2*nb+i,k+nb) = HessF(2*nb+i,k+nb)-(2*gkm*vk*vm*sin(tk - tm))/akm^2;
%     HessF(m+nb,2*nb+i) = HessF(m+nb,2*nb+i)+(2*gkm*vk*vm*sin(tk - tm))/akm^2;
%     HessF(2*nb+i,m+nb) = HessF(2*nb+i,m+nb)+(2*gkm*vk*vm*sin(tk - tm))/akm^2;
%    
%     HessF(2*nb+i,2*nb+i) = HessF(2*nb+i,2*nb+i)+gkm*((6*vk^2)/akm^4 - (4*vk*vm*cos(tk - tm))/akm^3);   
    
                        
end %for TAP

%% CÁLCULO EM RELAÇÃO AO PONTO DE VÁLVULA

%% HESSIANA RESTRIÇÃO DE DESIGUALDADE

% nseqht = length(seqht);
% 
% for i=1:ng
%     E = g(i).E;
%     F = g(i).F;
%     PGMIN = U3(2*nb+ntap+i,1);
%     
%    for a=1:nseqht
%       
%       if (Ger(i)==seqht(a))
%           pg = b(seqht(a)).PG;
%           HessH(2*nb+ntap+i,2*nb+ntap+i) = (E*F^2*sin(F*(pg - PGMIN)))*(lambda(2*ng+2*(2*nb+ntap+ng)+ngt+i)-lambda(2*ng+2*(2*nb+ntap+ng)+i));       % Derivada segunda em relação a pg restrição de desigualdade
%           HessH(2*nb+ntap+ng+a,2*nb+ntap+ng+a) = 0;                                                     % Derivada segunda em relação a m restrição de desigualdade
%           HessF(2*nb+ntap+ng+a,2*nb+ntap+ng+a) = 0;                                                     % Derivada segunda em realão a m função objetivo
%           
%       end
%    end
% end




%% DADOS DE SAÍDA

HessG = HessDP + HessDQ;

K = HessF + HessG + HessH;

end %function


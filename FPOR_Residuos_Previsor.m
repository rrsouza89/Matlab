function [mk,skprev,tk,uk,pkprev] = FPOR_Residuos_Previsor(G,H,GF,JGT,JHT,eta,lambda,Z,Zinv,L1,mi,delta,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                        RESÍDUOS PREVISOR                    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;

mk = -GF-JGT*eta-JHT*lambda;

skprev = -Z*lambda+mi*delta;

tk = -G;

uk = -H-z;

pkprev = JHT*Zinv*(skprev-L1*uk);

end
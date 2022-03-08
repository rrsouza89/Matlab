function [sk,pk] = FPOR_Residuos_Corretor(Z,Zinv,lambda,mi,delta,Dz,dlambdap,JHT,L1,uk)



sk = -Z*lambda + mi*delta - Dz*dlambdap;

pk = JHT*Zinv*(sk-L1*uk);


end

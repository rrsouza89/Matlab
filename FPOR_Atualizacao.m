function [b,r] = FPOR_Atualizacao(b,r,x0,nb,ntap,ng,ngt,conj_ramo_tap,NIN,caso);

for i=1:nb
    b(i).v = x0(i,1);
end

for i=1:nb
    b(i).teta = x0(nb+i,1);
end


for i=1:ntap
    r(conj_ramo_tap(i)).TAP = x0(2*nb+i,1);
end


for i=1:ng
    b(NIN(caso.gen(i,1))).PG = x0(2*nb+ntap+i,1);
end


cont = 0;
for i=1:ng
    if b(NIN(caso.gen(i,1))).TGER == 1 
        cont = cont + 1;
        b(NIN(caso.gen(i,1))).m = x0(2*nb+ntap+ng+cont,1);
    end
end



x0;

end
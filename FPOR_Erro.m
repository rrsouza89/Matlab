function [Erro] = FPOR_Erro(mk,skprev,tk,uk,sk,PC,it)



if PC(it) == 0
    
    Erro = norm([mk;skprev;tk;uk],inf);
    
else

    Erro = norm([mk;sk;tk;uk],inf);

end

end
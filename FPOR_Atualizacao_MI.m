function [mi,delta] = FPOR_Atualizacao_MI(mi,tal,z,lambda,H,G,Beta,delta,x0,tolz);

%% ATUALIZAÇÃO MI

if mi > 10^-11

    mi = tal*mi;

end
nz = length(z);

for i=1:nz
    if z(i)<-mi
        disp('Região relaxada');
        mi = -(1+tal)*min(z)
        
        break
    end
    
end


% r = 0.5;
% 
% if ((zp'+mi)*lambp) < 0.9*((zc'+mi)*lambc)
%     
%     mi = r*(((zp'+mi)*lambp)/(2*nb+ntap))
%     disp('mi_prev')
% else
%     mi = r*(((zc'+mi)*lambc)/(2*nb+ntap))
%     disp('mi_cor')
% end
 


%% ESTIMADORES DOS MULTIPLICADORES DE LAGRANGE

delta = lambda;

%%%%% Pinheiro

% if z(1)>tolz
%     
%     delta = delta;
%  
% else
%     
%     delta = lambda;
%     
% end

% delta1 = delta1;          
% delta2 = delta2;
% delta3 = delta3;
% delta4 = delta4;

% %%% polyak
% nh = length(lamb1);
% nx = length(lamb3);
% 
% for i=1:nh
%     delta1(i) = (mi*delta1(i))/(mi-(-z1(i)));
%     delta2(i) = (mi*delta2(i))/(mi-(-z2(i)));
% end
% 
% for i=1:nx
%     delta3(i) = (mi*delta3(i))/(mi-(-z3(i)));
%     delta4(i) = (mi*delta4(i))/(mi-(-z4(i)));
% end


end
function [Beta,FLBM,Dif] = FPOR_Beta(FLBM,FLBM2,Beta,it,dam1);



if it == 1 
    
    Beta = Beta;
    FLBM = FLBM;
    Dif = 0;
    
else

    Dif = FLBM-FLBM2;

    if Dif<.25
        Beta=Beta/dam1;
        %fprintf('Dif < 0.25 => Beta = %f\n\n',Beta)    
    end

    if Dif>.75
        Beta=Beta*dam1;
        %fprintf('Dif > 0.75 => Beta = %f\n\n',Beta)
    end

    FLBM = FLBM2;
    
end


% Levenberg-Marquardt

% if Dif<.25
%     Beta = 4*Beta;
% end
% 
% if Dif>.75
%     Beta = Beta/2;
% end

end
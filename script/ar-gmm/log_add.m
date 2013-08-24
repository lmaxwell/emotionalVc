function [ logc ] = log_add( loga ,logb )
%LOG_ADD Summary of this function goes here
%   Detailed explanation goes here
    if(loga<logb)
        swap=loga;
        loga=logb;
        logb=swap;
    end
    diff=logb-loga;
    minLogExp=-log(1.0E10);
    if(diff<minLogExp)
        if(loga<-0.5E10)
            logc=-1.0E10;
        else
            logc=loga;
        end
    else
        logc=loga+log(1.0+exp(diff));
    end
end


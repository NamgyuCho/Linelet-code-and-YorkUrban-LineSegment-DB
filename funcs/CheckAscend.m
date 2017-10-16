function [ ret ] = CheckAscend( llSet, ll_type )
    %CHECKASCEND Summary of this function goes here
    %   Detailed explanation goes here
    
    % llSet = curllSet; 
    
    if size(llSet,1) == 1
        ret = true;
    else
        [iMaj, iMin] = GetMajorAxis( ll_type );
        
        [~, iSrt] = sort(llSet(:,iMaj(1)));
        
        if llSet(iSrt(1), iMin(1)) - llSet(iSrt(end), iMin(1)) < 0
            ret = true;
        else
            ret = false;
        end
    end
end


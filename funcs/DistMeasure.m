function [ d ] = DistMeasure( h1, h2, strType )
    %DISTMEASURE Summary of this function goes here
    %   Detailed explanation goes here
    % h1 = h_inst(1,:); h2 = h_inst(2,:);
    
    if strcmp(strType, 'chi2')
        idxValid = find(h1 ~=0 | h2 ~= 0);
        
        nom = (h1(idxValid) - h2(idxValid)).^2 + eps;
        denom = (h1(idxValid) + h2(idxValid)) + eps;
        d = sum( nom ./ denom ) / 2;
    end    
end
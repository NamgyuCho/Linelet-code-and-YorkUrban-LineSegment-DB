function [ ppot ] = PairPot( LL_i, LL_j, im_grad, im_dir )
    %PAIRPOT Summary of this function goes here
    %   Detailed explanation goes here
    
    ppot = 1;
    
    % 1. length difference
    if abs(LL_i(end) - LL_j(end)) > 2,  ppot = 0;    end
        
    % 2. minor position difference
    if abs(LL_i(2) - LL_j(2)) ~= 1,  ppot = 0;    end
    
    % 3. slope difference
    if abs(LL_i(end) - LL_j(end)) > 2,  ppot = 0;    end
    
end


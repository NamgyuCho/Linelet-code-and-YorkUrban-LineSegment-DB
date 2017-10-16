function [ perpDist ] = GetPerpDist( pt0, ptSet, angle )
    %GETPERPDIST Summary of this function goes here
    %   Detailed explanation goes here
    
    dirVec = [-sin(angle) cos(angle)];
    perpDist = abs( dirVec*(ptSet - repmat( pt0, size( ptSet, 1 ), 1 ))' );
end


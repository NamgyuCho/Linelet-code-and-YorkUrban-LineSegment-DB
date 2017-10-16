function [ ret ] = LineletRearrange( llSet, ll_type )
    %LINELETREARRANGE Summary of this function goes here
    %   Detailed explanation goes here
        
    if strcmp( ll_type, 'hor' )
        major_idx1 = 1;
        major_idx2 = 3;
        minor_idx1 = 2;
        minor_idx2 = 2;
    elseif strcmp( ll_type, 'ver' )
        major_idx1 = 2;
        major_idx2 = 4;
        minor_idx1 = 1;
        minor_idx2 = 1;
    else % diag
        major_idx1 = 1;
        major_idx2 = 3;
        minor_idx1 = 2;
        minor_idx2 = 4;
    end
    
    ptSet = [];
    for i = 1:size(llSet)
        ptSet = [ptSet; ConvertLinelet2Pts(llSet(i,:), ll_type)];
    end
    
    ret = ConvertPts2Linelet( ptSet, ll_type );    
end


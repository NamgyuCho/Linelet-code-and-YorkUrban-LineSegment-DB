function [ ptSet ] = ConvertLinelet2Pts( ll_inst, ll_type )
    %CONVERTLINELET2PTS Summary of this function goes here
    %   Detailed explanation goes here
    
    ptSet = [];
    
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
    
    for i = 1:size(ll_inst, 1)
        if strcmp( ll_type, 'hor' )
            xx1 = ll_inst(i,major_idx1):ll_inst(i,major_idx2);
            yy1 = repmat(ll_inst(i,minor_idx1), 1, ll_inst(i,5));
            ptSet = [ptSet; xx1' yy1'];
        else
            xx1 = repmat(ll_inst(i,minor_idx1), 1, ll_inst(i,5));
            yy1 = ll_inst(i,major_idx1):ll_inst(i,major_idx2);
            ptSet = [ptSet; xx1' yy1'];
        end
    end
end


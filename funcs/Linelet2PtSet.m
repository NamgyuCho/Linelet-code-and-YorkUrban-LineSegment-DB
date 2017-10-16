function [ ptSet ] = Linelet2PtSet( llSet, ll_type, strMode )
    %LINELET2PTSET Summary of this function goes here
    %   Detailed explanation goes here
    
    ptSet = [];
    for i = 1:size( llSet,1 )
        if strcmp(ll_type, 'hor')%llSet(i,6) == 1
            xx = [llSet(i,1):llSet(i,3)]';
            yy = ones( llSet(i,5), 1 ) * llSet(i,2);
        elseif strcmp(ll_type, 'ver')%llSet(i,6) == 3
            xx = ones( llSet(i,5), 1 ) * llSet(i,1);
            yy = [llSet(i,2):llSet(i,4)]';
        else
            x1 = llSet(i,1);
            x2 = llSet(i,3);
            y1 = llSet(i,2);
            y2 = llSet(i,4);
            
            if x1 < x2, xx = x1:x2;
            else        xx = x1:-1:x2;
            end
            xx = xx';
            
            if y1 < y2, yy = y1:y2;
            else        yy = y1:-1:y2;
            end
            yy = yy';
        end
        ptSet = [ptSet; xx, yy];
    end
    
end


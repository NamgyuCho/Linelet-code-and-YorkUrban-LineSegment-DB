function [ ll_inst ] = ConvertPts2Linelet( ptSet, ll_type )
    %CONVERTPTS2LINELET Summary of this function goes here
    %   Detailed explanation goes here
    ll_inst = [];
    xx = ptSet(:,1);
    yy = ptSet(:,2);
    if strcmp(ll_type, 'ver')
        vUnq = unique(xx);
        for i = 1:length(vUnq)
            ymin = min(yy(xx == vUnq(i)));
            ymax = max(yy(xx == vUnq(i)));
            ll_inst = [ll_inst; vUnq(i), ymin, vUnq(i), ymax, ymax-ymin+1, 0];
        end
    else
        vUnq = unique(yy);
        for i = 1:length(vUnq)
            xmin = min(xx(yy == vUnq(i)));
            xmax = max(xx(yy == vUnq(i)));
            ll_inst = [ll_inst; xmin, vUnq(i), xmax, vUnq(i), xmax-xmin+1, 0];
        end
    end    
end


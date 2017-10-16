function [ ll_end ] = GetEndLinelet( ll_instCand )
    %GETENDLINELET Summary of this function goes here
    %   Detailed explanation goes here
    
    % ll_instCand = ll_inst(ll_idxSet{lsCandSet(k)},:);

    if size(ll_instCand,1) > 1
        [~, srt_idx] = sort(ll_instCand(:,1));
        ll_instCand = ll_instCand(srt_idx,:);
        ll_end = [ll_instCand(1,:); ll_instCand(end,:)];
    else
        ll_end = ll_instCand;
    end
end


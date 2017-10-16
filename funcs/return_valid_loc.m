function [ pts, valid_idx ] = return_valid_loc( pts, size_range, offset )
%RETURN_VALID_LOC Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(offset)
    pt_moved = pts + repmat(offset, size(pts,1), 1);
    valid_idx = pts(:,1) > 0 & pts(:,1) <= size_range(2) &...
                pts(:,2) > 0 & pts(:,2) <= size_range(1) &...
                pt_moved(:,1) > 0 & pt_moved(:,1) <= size_range(2) &...
                pt_moved(:,2) > 0 & pt_moved(:,2) <= size_range(1);
else
    valid_idx = pts(:,1) > 0 & pts(:,1) <= size_range(2) &...
                pts(:,2) > 0 & pts(:,2) <= size_range(1);
end
pts = pts(valid_idx,:);

end


function [ val, ind, sub ] = get_matrix_extreme_val( A, str_mode )
%GET_MATRIX_EXTREME_VAL Summary of this function goes here
%   Detailed explanation goes here

size_mat = size(A);

if strcmp(str_mode, 'max')
    [val, ind] = max(A(:));
else
    [val, ind] = min(A(:));
end
[rr, cc] = ind2sub(size_mat, ind);
sub = [rr cc];

end


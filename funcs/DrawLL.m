function [ im_tmp ] = DrawLL( ll_cand, size_im, ll_type )
%DRAWLL Summary of this function goes here
%   Detailed explanation goes here

if isempty(ll_type)
    ll_type = 'hor';
end

%ll_cand = ll_cand_new;
im_tmp = zeros(size_im);

for i = 1:size(ll_cand,1)
   
    if strcmp(ll_type, 'hor')
        idx = sub2ind(size_im, repmat(ll_cand(i,2), 1, ll_cand(i,5)), ll_cand(i,1):ll_cand(i,3));
    elseif strcmp(ll_type, 'ver')
        idx = sub2ind(size_im, ll_cand(i,2):ll_cand(i,4), repmat(ll_cand(i,1), 1, ll_cand(i,5)));
    else
    end
    im_tmp(idx) = im_tmp(idx) + 1;
end



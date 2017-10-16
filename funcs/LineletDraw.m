function [ im_tmp ] = LineletDraw( ll_cand, LL_type, size_im )
    %LINELETDRAW Summary of this function goes here
    %   Detailed explanation goes here

    im_tmp = zeros(size_im);
    
    for i = 1:size(ll_cand,1)
        if ll_cand(i,2) > ll_cand(i,4)
            yint = -1;
        else
            yint = 1;
        end
        
        if ll_cand(i,1) > ll_cand(i,3)
            xint = -1;
        else
            xint = 1;
        end
        
        if strcmp(LL_type, 'hor')
            idx = sub2ind(size_im, repmat(ll_cand(i,2), 1, ll_cand(i,5)), ll_cand(i,1):xint:ll_cand(i,3));
        elseif strcmp(LL_type, 'ver')
            idx = sub2ind(size_im, ll_cand(i,2):yint:ll_cand(i,4), repmat(ll_cand(i,1), 1, ll_cand(i,5)));
        else
            idx = sub2ind(size_im, ll_cand(i,2):yint:ll_cand(i,4), ll_cand(i,1):xint:ll_cand(i,3));
        end
        im_tmp(idx) = im_tmp(idx) + 1;
    end
end

% function [ im_tmp ] = DrawLL( ll_cand, size_im, ll_type )
% %DRAWLL Summary of this function goes here
% %   Detailed explanation goes here
% 
% if isempty(ll_type)
%     ll_type = 'hor';
% end
% 
% %ll_cand = ll_cand_new;
% im_tmp = zeros(size_im);
% 
% for i = 1:size(ll_cand,1)
%    
%     if strcmp(ll_type, 'hor')
%         idx = sub2ind(size_im, repmat(ll_cand(i,3), 1, ll_cand(i,4)), ll_cand(i,1):ll_cand(i,2));
%     elseif strcmp(ll_type, 'ver')
%         idx = sub2ind(size_im, repmat(ll_cand(i,3), 1, ll_cand(i,4)), ll_cand(i,1):ll_cand(i,2));
%     else
%     end
%     im_tmp(idx) = im_tmp(idx) + 1;
% end
% 
% 

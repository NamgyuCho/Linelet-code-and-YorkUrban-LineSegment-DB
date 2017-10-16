function [ valid_score ] = LineletValidation( ll_inst, im_gray, im_grad, im_dir, ll_type, func_type )
%LINELETVALIDATION Summary of this function goes here
%   Detailed explanation goes here
% ll_inst = ll_Hor; ll_type = 'hor'; func_type = 'NFA';

num_inst = size(ll_inst, 1);
valid_score = zeros(num_inst, 1);
size_im = size(im_grad);

for i = 1:num_inst
    [ ohm, ~, ~ ] = LineletOmegaGamma( ll_inst(i,:), ll_type, size_im );
    %pts = ind2sub(size_im, ohm);
    
    if strcmp(func_type, 'NFA')
%         valid_score(i) = NFA_linelet( [0 0 median(im_dird(ohm))*pi/180], ohm, im_dird);
        valid_score(i) = NFA_linelet2( ohm, im_dir );
    elseif strcmp(func_type, 'MixExpert')
        valid_score(i) = LineletForegroundProbability( ll_inst, im_gray, im_grad, im_dir );
    end
end

end


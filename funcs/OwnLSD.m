function [ line_own ] = OwnLSD( im_gray, param, BN_Model, theta_BN )
    %OWNLSD Summary of this function goes here
    %   Detailed explanation goes here
    %%
    [im_grad, im_dir, im_dird, im_gx, im_gy] = get_gradient(im_gray, param.grad.filter, param.grad.sigma_size);
    % figure; imshow(im); figure; imshow(uint8(im_grad/max(im_grad(:))*255)); colormap hot; 
    % figure; imagesc(im_grad); axis image; figure; imagesc(im_dird); axis image
    %im_grad = im_grad / 4;
    
    % Remove less significant pixels wrt the gradient magnitue
    ind = im_grad <= param.thres_1st_grad_mag;
    im_grad(ind) = 0;   im_dir(ind) = 0; im_dird(ind) = 0;  im_gx(ind) = 0;   im_gy(ind) = 0;   
    
    % Get line segment candidate map
    [~, lsc_map, ~, ~, ~] = Get_Line_Candidate(im_gray, im_grad, im_grad, im_dir, im_dird, param, false, false);
    % figure(10); imagesc(lsc_map); axis image;
    
    
    im_dird(im_dird < 0) = im_dird(im_dird < 0) + 180;
    im_dir(im_dir < 0) = im_dir(im_dir < 0) + pi;
    
    
    %% Get line segments
    if strcmp(param.est_method_main, 'voting1')
        lsc_hor = linelet_voting( im_gray, param, BN_Model, theta_BN );
        lsc_ver = linelet_voting( im_gray', param, BN_Model, theta_BN );
        lsc_ver = [lsc_ver(:,2) lsc_ver(:,1) -lsc_ver(:,3)+pi/2 lsc_ver(:,4)]; % Rotate to the target direction, vertical
        line_own = [lsc_hor; lsc_ver];

    elseif strcmp(param.est_method_main, 'voting2')
        % version 2
        lsc_hor = linelet_voting2( im_grad, im_dir, lsc_map, param, BN_Model, theta_BN );
        lsc_ver = linelet_voting2( im_grad', im_dir', lsc_map', param, BN_Model, theta_BN );
        lsc_ver = [lsc_ver(:,2) lsc_ver(:,1) -lsc_ver(:,3)+pi/2 lsc_ver(:,4)]; % Rotate to the target direction, vertical
        line_own = [lsc_hor; lsc_ver];
    elseif strcmp(param.est_method_main, 'voting3')
        
        % version 3
        line_own = linelet_voting3( im_grad, im_dir, lsc_map, param, im_dird );
    elseif strcmp(param.est_method_main, 'voting4')
        line_own = linelet_voting4( im_grad, im_dir, lsc_map, param, im_dird );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

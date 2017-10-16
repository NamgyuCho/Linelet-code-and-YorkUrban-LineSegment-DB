function [ pAng, ls_est, pFg, bValidHelmholtz ] = LineletEstimation( ...
    ll_idxSet, ll_inst, upot1, upot2, bAscend, ll_type, ll_Valid, size_im, level_t, theta_space, im_gray, im_grad, im_dir, param, bDraw )
  
    
    num_set = size(ll_idxSet, 1);
    size_theta = length(theta_space);
    pAng = zeros(num_set, size_theta);
    ls_est = zeros(num_set, 6);
    pFg = zeros(num_set, 2);
    bValidHelmholtz = false(num_set, 1);
    ppot1 = zeros(num_set, 1);
    ppot2 = zeros(num_set, size_theta);
    
    if strcmp(ll_type, 'ver')
        minorAxisIdx = 1;
        majorAxisIdx = 2;
    else
        minorAxisIdx = 2;
        majorAxisIdx = 1;
    end
    
    for i = 1:num_set
        try
            % Get target set index, instance, and potential
            curCandIdx0 = ll_idxSet{i};
            if isempty(curCandIdx0), continue; end
            
            curllSet = ll_inst(curCandIdx0, :);
            curUpot1 = upot1(curCandIdx0, :);
            curUpot2 = upot2(curCandIdx0,:);
            bAscend(i) = CheckAscend( curllSet, ll_type );                      
            
            [pAng(i, :), ls_est(i,:), ~] = LineletInference( curllSet, curUpot1, curUpot2, ll_type, bAscend(i), theta_space );
            
            % Validate -- proposed method
            pFg(i,:) = FgBgProbability( ls_est(i,:), curllSet, im_gray, im_grad, im_dir, ll_type, param );
            
            % Validate -- Helmholtz 
            ptSet = ConvertLinelet2Pts( curllSet, ll_type );
            ptIdx = sub2ind(size_im, ptSet(:,2), ptSet(:,1));
            ptOri = im_dir(ptIdx);
            if NFA_linelet(ptOri, im_dir) > param.thres_log_eps
                bValidHelmholtz(i) = true;
            else
                bValidHelmholtz(i) = false;
            end
                                   
            if bDraw
                % curllSet
                im_inst_ds = LineletDraw(curllSet, ll_type, size_im);
%                 imw = ones(size_im(1), size_im(2), 3) * 255;
                figure(3); clf; imagesc(im_inst_ds); axis image; hold on; colormap gray; title([num2str(i) ' result']);
                x1 = ls_est(i,1:2) + ls_est(i,4)/2*[cos(ls_est(i,3)) sin(ls_est(i,3))];
                x2 = ls_est(i,1:2) - ls_est(i,4)/2*[cos(ls_est(i,3)) sin(ls_est(i,3))];
                plot([x1(1) x2(1)], [x1(2) x2(2)]+1, '--', 'color', [1 0 0], 'linewidth', 1); % updated
%                 text(cent_est(1), cent_est(2), num2str(i));
                waitforbuttonpress;
            end            
        catch err
            fprintf('error at %d.\n', i);
            rethrow(err);
        end
    end
end
function [BN_cand, Pot] = LineletPot(im_grad, ll_cand, BN_Model)
% Energy = independent pot + pairwise (markov chain) pot + gradient pot
% i) independent pot -- simply indicator function
% ii) pairwise pot
% iii) gradient pot


size_im = size(im_grad);
Pot = 0;

% sort linelets wrt the length
len_max = max(ll_cand(:,4));
ll_cand(len_max - ll_cand(:,4) > 2,:) = [];
num_ll = size(ll_cand, 1);

try
    % ------------------------------------------------------------------------------------------------------------------
    % i) independent pot -- simply indicator function
    % ------------------------------------------------------------------------------------------------------------------
    BN_cand = max(BN_Model.Numbers,[], 2) - max(ll_cand(:,4)) >= 0;
    BN_cand_min = min(BN_Model.Numbers,[], 2) - min(ll_cand(:,4)) >= 0;
    idx = find(BN_cand ~= 0);
    idx1 = find(BN_cand ~= 0 & BN_cand_min ~= 0);
    pot_indep = ones(length(idx), 1);
    BN_cand = idx;
    % ------------------------------------------------------------------------------------------------------------------
    
    
    % ------------------------------------------------------------------------------------------------------------------
    % ii) pairwise pot
    % ------------------------------------------------------------------------------------------------------------------
    pot_pair = ones(length(idx), 1);
    
    if num_ll == 1
        %pot_pair(:) = 1;
    else
        for i = 1:length(idx)
            model_state = BN_Model.State{idx(i)};
            trans_matri = BN_Model.Transition{idx(i)};
            
            if numel(model_state) == 1
                observ_state = ones(num_ll,1);
            else
                [~, observ_state] = min(repmat(model_state', size(ll_cand,1),1) - repmat(ll_cand(:,4), 1, 2), [], 2);
            end
            for j = 1:num_ll - 1
                pot_pair(i) = pot_pair(i) * trans_matri(observ_state(j), observ_state(j+1));
            end
        end
    end
    % ------------------------------------------------------------------------------------------------------------------
    
    
    % ------------------------------------------------------------------------------------------------------------------
    % iii) gradient pot
    % ------------------------------------------------------------------------------------------------------------------
    pot_grad = ones(length(idx), 1);
    
    ll_coef = [];
    for i = 1:num_ll
        idx_linelet = sub2ind(size_im, repmat(ll_cand(i,3), 1, ll_cand(i,4)), ll_cand(i,1):ll_cand(i,2));
        val_grad = im_grad(idx_linelet) / max(im_grad(idx_linelet));
        ll_coef = [ll_coef; mean(abs(diff(val_grad)))];
    end
    pot_grad = 1 - abs(BN_Model.Coef(idx) - mean(ll_coef));
    % ------------------------------------------------------------------------------------------------------------------
    
    Pot = pot_indep + pot_pair + pot_grad;
catch exception
    fprintf('Err in LineletPot()\n');
    rethrow(exception);
end
% % check connection state and gradient maginitude distribution along and near the current ll
% len_ll = size(ll_cand,1);
%
% % Retrieve gradient magnitudes algon linelet
% if len_ll == 1
%     xx = ll_cand(1):ll_cand(2);
%     yy = repmat(ll_cand(3), 1, length(xx));
%
%     grad_main = im_grad(yy, xx);
%     grad_up = im_grad(yy-1, xx);
%     grad_down = im_grad(yy+1, xx);
% elseif len_ll == 2
%     slope = diff(ll_cand(:,3));
%     len_diff = diff(ll_cand(:,4)); %
%     xx = [ll_cand(1,1):ll_cand(1,2) ll_cand(2,1):ll_cand(2,2)];
%     yy = [repmat(ll_cand(1,3), 1, length(ll_cand(1,1):ll_cand(1,2))) repmat(ll_cand(2,3), 1, length(ll_cand(2,1):ll_cand(2,2)))];
% else
%     slope = diff(ll_cand(:,3));
%     xx = []; yy = [];
%     for j=1:len_ll
%         xx = [xx ll_cand(j,1):ll_cand(j,2)];
%         yy = [yy repmat(ll_cand(j,3), 1, length(ll_cand(j,1):ll_cand(j,2)))];
%     end
% end


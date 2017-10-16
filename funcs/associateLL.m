function [ ll_Hor_new, ll_Ver_new, im_Hor, im_Ver ] = associateLL(ll_Hor, ll_Ver, ll_Diag, im_grad, im_dird, size_im, bDrawIntermediateResult)
%ASSOCIATELL Summary of this function goes here
%   Detailed explanation goes here

%% Horizontally

[~, idx] = sort(ll_Hor(:,end), 'descend');
ll_Hor_old = ll_Hor(idx,:);
ll_Hor_old(ll_Hor_old(:,end) < 2, :) = [];
num_ll = size(ll_Hor_old, 1);

ll_Hor_new = [];
num_Hor_new = 1;
ll_usage = false(num_ll,1);

im_Hor = [];    im_Ver = [];

im_Hor = zeros(size_im); 

for ii = 1:num_ll
    if ll_usage(ii), continue; end
    
    stack_ll = ii;
    ll_cand = [];
    
    while ~isempty(stack_ll)
        cur_ll = stack_ll(1);   stack_ll(1) = []; % pop
        
        % Search possible connection with other linelets. At this step, 
        % also check possible connection over one empty neighbor. This is 
        % related to the minor axis update, but it is not certain yet
        idx_left = (ll_Hor_old(cur_ll,1) - ll_Hor_old(:,3) == 1) &...
            (abs(ll_Hor_old(:,2) - ll_Hor_old(cur_ll,2)) == 1); % Leftward distance
        idx_right = (ll_Hor_old(:,1) - ll_Hor_old(cur_ll,3) == 1) &...
            (abs(ll_Hor_old(:,2) - ll_Hor_old(cur_ll,2)) == 1); % Rightward distance
%         idx_ang = bAngleAligned( ll_Hor_old(cur_ll,5), ll_Hor_old(:,5), param.thres_angle_diff ); % Gradient angle
        idx_len = abs(ll_Hor_old(:,5) - ll_Hor_old(cur_ll,5)) <= 3;
        
        idx_nb = find((idx_right | idx_left) & idx_len & ~ll_usage);
        idx_nb(idx_nb == cur_ll) = [];
        
        stack_ll = [stack_ll; idx_nb];
        
        ll_usage(idx_nb) = true;
        ll_usage(cur_ll) = true;
        
        ll_cand = [ll_cand; ll_Hor_old(cur_ll,:)];
    end
    
    if isempty(ll_cand)
        continue; 
    else
        % Attch short linelet at both sides
        [~, idx1] = sort(ll_cand(:,1));
        
        % Left one
        idx_left = (ll_cand(idx1(1),1) - ll_Hor_old(:,3) == 1) &...
            (abs(ll_Hor_old(:,2) - ll_Hor_old(idx1(1),2)) == 1); % Leftward distance
        % Right one
        idx_right = (ll_Hor_old(:,1) - ll_Hor_old(idx1(end),3) == 1) &...
            (abs(ll_Hor_old(:,2) - ll_Hor_old(idx1(end),2)) == 1); % Rightward distance
        idx_nb = find((idx_right | idx_left) & ~ll_usage);
        ll_usage(idx_nb) = true;
        
        ll_cand = [ll_cand; ll_Hor_old(idx_nb,:)];
        ll_Hor_new{num_Hor_new,1} = ll_cand;
        num_Hor_new = num_Hor_new + 1;
        
        im_Hor = im_Hor + DrawLL(ll_cand, size_im, 'hor') * num_Hor_new;
        
        continue;
    end
end

if bDrawIntermediateResult
    figure; imagesc(im_Hor);
end

%% Vertically

[~, idx] = sort(ll_Ver(:,end), 'descend');
ll_Ver_old = ll_Ver(idx,:);
ll_Ver_old(ll_Ver_old(:,end) < 2, :) = [];
num_ll = size(ll_Ver_old, 1);

ll_Ver_new = [];
num_Ver_new = 1;
ll_usage = false(num_ll,1);

im_Ver = zeros(size_im(2:-1:1)); 

for ii = 1:num_ll
    if ll_usage(ii), continue; end
    
    stack_ll = ii;
    ll_cand = [];
    
    while ~isempty(stack_ll)
        cur_ll = stack_ll(1);   stack_ll(1) = []; % pop
        
        % Search possible connection with other linelets. At this step, 
        % also check possible connection over one empty neighbor. This is 
        % related to the minor axis update, but it is not certain yet
        idx_left = (ll_Ver_old(cur_ll,2) - ll_Ver_old(:,4) == 1) &...
            (abs(ll_Ver_old(:,1) - ll_Ver_old(cur_ll,1)) == 1); % Leftward distance
        idx_right = (ll_Ver_old(:,2) - ll_Ver_old(cur_ll,4) == 1) &...
            (abs(ll_Ver_old(:,1) - ll_Ver_old(cur_ll,1)) == 1); % Rightward distance
%         idx_ang = bAngleAligned( ll_Ver_old(cur_ll,5), ll_Ver_old(:,5), param.thres_angle_diff ); % Gradient angle
        idx_len = abs(ll_Ver_old(:,5) - ll_Ver_old(cur_ll,5)) <= 3;
        
        idx_nb = find((idx_right | idx_left) & idx_len & ~ll_usage);
        idx_nb(idx_nb == cur_ll) = [];
        
        stack_ll = [stack_ll; idx_nb];
        
        ll_usage(idx_nb) = true;
        ll_usage(cur_ll) = true;
        
        ll_cand = [ll_cand; ll_Ver_old(cur_ll,:)];
    end
    
    if isempty(ll_cand), 
        continue; 
    else
        % Attch short linelet at both sides
        [~, idx1] = sort(ll_cand(:,2));
        
        idx_left = (ll_Ver_old(idx1(1),2) - ll_Ver_old(:,4) == 1) &...
            (abs(ll_Ver_old(:,1) - ll_Ver_old(idx1(1),1)) == 1); % Leftward distance
        idx_right = (ll_Ver_old(:,2) - ll_Ver_old(idx1(end),4) == 1) &...
            (abs(ll_Ver_old(:,1) - ll_Ver_old(idx1(end),1)) == 1); % Rightward distance
        idx_nb = find((idx_right | idx_left) & ~ll_usage);
        ll_usage(idx_nb) = true;
        
        ll_cand = [ll_cand; ll_Ver_old(idx_nb,:)];
        
        ll_Ver_new{num_Ver_new,1} = ll_cand;
        num_Ver_new = num_Ver_new + 1;
        
        im_Ver = im_Ver + DrawLL(ll_cand, size_im, 'ver') * num_Ver_new;
        continue;
    end
end

if bDrawIntermediateResult
    figure; imagesc((im_Ver)');
end

end


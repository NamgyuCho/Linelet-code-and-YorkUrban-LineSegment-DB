function [BN_Model] = generate_BN_LUT(theta_BN)

% Generate a set of pattern of bresenham lines
%%
% tmp_mat = ones(max(size_im));
% [tr_r, tr_c] = find(tril(tmp_mat));
% 
% dx = (1:max(size_im)).^-1;
% dy = 1:max(size_im);
% 
% grad_tot = dx' * dy;
% unq_grad = unique(grad_tot);

%theta = theta_BN;
%ind_theta_mid = ceil(length(theta_BN)/2);
%AA_Coef=[];

max_len = 10000;
minor_incre_table = floor(tan(theta_BN') * [1:max_len]);
rasidual = tan(theta_BN') * [1:max_len] - minor_incre_table;

BN_Model = [];
BN_State = [];
BN_Transition = [];
BN_Freq = [];
BN_Coef = [];
BN_Max_Num = 0;
for i = 2:size(minor_incre_table,1)
    val_diff = diff(minor_incre_table(i,:));
    idx = find(val_diff == 1);
    major_no = diff(idx);
    major_no((major_no - max(major_no)) > 2) = [];
    [val_unq, ~] = count_unique(major_no);
    
    
    % Learn the transition matrix
    num_state = length(val_unq);
    [val_state, idx_state] = sort(val_unq);
    
    transition_mat = zeros(length(val_unq)) + eps;
    freq_mat = zeros(length(val_unq),1);
    
    if length(val_unq) == 1
        transition_mat = 1;
        freq_mat = 1;
    else
        for j = 1:length(major_no)-1
            if j <= length(major_no) - 1
                id_prev = find(val_state == major_no(j));
                id_next = find(val_state == major_no(j+1));
                transition_mat(id_prev, id_next) = transition_mat(id_prev, id_next) + 1;
            end
            id_curr = find(val_state == major_no(j));
            freq_mat(id_curr) = freq_mat(id_curr) + 1;
        end
        transition_mat = transition_mat / sum(transition_mat(:));
    end 
    
    BN_State{i,1} = val_state;
    BN_Transition{i,1} = transition_mat;
    BN_Freq{i,1} = freq_mat / sum(freq_mat);
    BN_Max_Num = max(BN_Max_Num, numel(val_state));
    
    % Get anti-aliasing coefficents
    rasi_diff = diff(rasidual(i,:));
    BN_Coef(i,1) = mean(rasi_diff(rasi_diff > 0));
end
BN_State{1,1} = max_len;
BN_Transition{1,1} = 1;
BN_Freq{1,1} = 1;
BN_Coef(1,1) = 0;
BN_Coef(end) = 1;

BN_Model.State = BN_State;
BN_Model.Transition = BN_Transition;
BN_Model.Freq = BN_Freq;
BN_Model.Coef = BN_Coef;


BN_Single_Num = zeros(size(minor_incre_table,1), BN_Max_Num);
for i = 1:size(minor_incre_table,1)
    num_el_diff = BN_Max_Num - numel(BN_Model.State{i});
    if num_el_diff > 0,     BN_Single_Num(i,:) = [BN_Model.State{i}' BN_Model.State{i}'];
    else                    BN_Single_Num(i,:) = BN_Model.State{i}';
    end
end
BN_Model.Numbers = BN_Single_Num;

% Divide model group wrt linelet size 
vec_unq = unique(BN_Single_Num, 'rows');
acc_num = zeros(size(vec_unq,1), 1);
num_total = size(vec_unq,1);
for i = 1:size(vec_unq,1)
    tf = (BN_Single_Num ==  repmat(vec_unq(i,:), size(BN_Single_Num,1), 1));
    idx = find(sum(tf,2) == 2);
    acc_num(num_total-i+1,1) = length(idx);
end


% if strcmp(strMode, 'fit_image')
% %     for i = 2:size_im(1)
% %         for j = 1:size_im(2)            
% %             [x,y] = bresenham(1, 1, 100, 20);
% %             [~, tmp_pat] = count_unique(y);
% %             tmp_pat_diff = diff(tmp_pat);
% %             idx = find(tmp_pat_diff == 1 | tmp_pat_diff == -1 | tmp_pat_diff == 0 | tmp_pat_diff);
% %         
% %             if numel(idx) == 0, continue; end
% %             bi_chain = tmp_pat(idx(1):idx(1)+1)';
% %             
% %             for k = 2:length(idx)
% %                 if ~ismember(tmp_pat(idx(k):idx(k)+1)', bi_chain, 'rows');
% %                     bi_chain = [bi_chain; tmp_pat(idx(k):idx(k)+1)'];
% %                 end
% %             end
% %         end
% %     end
% else
%     ang_rng = theta_BN;
%     % ang_rng = 0:thetaSampleFrequency:pi/4; %(0:44)*pi/180;
%     arrDy = sin(ang_rng);
%     size_table = size(arrDy, 2);
%     size_im_max = max(size_im);
%     
%     BN_LUT = zeros(size_table, size_im_max);
%     BN_LUT_Pattern = cell(size_table, 1);
%     %BN_LUT_Pattern{1} = [0 0];
%     
%     BN_LUT_Single_Num = cell(size_table, 1);
%     %BN_LUT_Single_Num{1} = 1;
%     %BN_LUT_Single_Num{end} = 1;
%     BN_LUT_Single_Num_mat_col = 0;
%     
%     for i = 2:size_table-1
%         BN_LUT(i,:) = floor(arrDy(i) * (1:size_im_max));
%         [~, tmp_pat] = count_unique(BN_LUT(i,:));
%         tmp_pat_diff = diff(tmp_pat);
%         idx = find(tmp_pat_diff == 1 | tmp_pat_diff == -1 | tmp_pat_diff == 0 | tmp_pat_diff);
%         
%         if numel(idx) == 0
%             bi_chain = unique(tmp_pat);
%         else
%             bi_chain = tmp_pat(idx(1):idx(1)+1)';            
%             for j = 2:length(idx)
%                 if ~ismember(tmp_pat(idx(j):idx(j)+1)', bi_chain, 'rows')
%                     bi_chain = [bi_chain; tmp_pat(idx(j):idx(j)+1)'];   
%                 end
%             end
%         end
%         BN_LUT_Pattern{i} = bi_chain;
%         BN_LUT_Single_Num{i} = unique(bi_chain(:));
%         BN_LUT_Single_Num_mat_col = max(BN_LUT_Single_Num_mat_col, numel(unique(bi_chain(:))));
%     end
%     
%     BN_LUT_Single_Num_Mat = zeros(size_table, BN_LUT_Single_Num_mat_col);
%     for i = 1:size_table
%         num_el_diff = BN_LUT_Single_Num_mat_col - numel(BN_LUT_Single_Num{i});
%         if num_el_diff > 0,     BN_LUT_Single_Num_Mat(i,:) = [BN_LUT_Single_Num{i}' zeros(1,num_el_diff)];
%         else                    BN_LUT_Single_Num_Mat(i,:) = BN_LUT_Single_Num{i}';
%         end
%     end
%     
%     BN_LUT_Pattern{1} = size_im_max;
%     %BN_LUT_Pattern{ind_theta_mid} = size_im_max;
%     BN_LUT_Pattern{end} = [1 1; 1 1; 1 1];
%     
%     BN_LUT_Single_Num_Mat(1,:) = [size_im_max 0 0];
%     %BN_LUT_Single_Num_Mat(ind_theta_mid,:) = size_im_max;
%     BN_LUT_Single_Num_Mat(end,:) = [1 0 0];
%     
%     BN_LUT_Single_Num{1} = size_im_max;
%     %BN_LUT_Single_Num{ind_theta_mid} = size_im_max;
%     BN_LUT_Single_Num{end} = 1;
%     
%     
%     BN_LUT_Max = max(BN_LUT_Single_Num_Mat, [], 2);
%     
%     BN_LUT_Bin = unique(BN_LUT_Single_Num_Mat);
% end
% 
% BN.LUT = BN_LUT;
% BN.LUT_Pattern = BN_LUT_Pattern;
% BN.LUT_Single_Num = BN_LUT_Single_Num;
% BN.LUT_Single_Num_Mat = BN_LUT_Single_Num_Mat;
% BN.LUT_Max = BN_LUT_Max;
% BN.LUT_Bin = BN_LUT_Bin;
% BN.AA_Coef = AA_Coef;

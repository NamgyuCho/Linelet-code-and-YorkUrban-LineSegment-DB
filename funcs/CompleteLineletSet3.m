function [ ll_cand_new, conn_map_new ] = CompleteLineletSet3( ll_cand, param )
%COMPLETELINELETSET2 Summary of this function goes here
%   Detailed explanation goes here

% Procedure
% Improve linelet -> Connection&Consistency Map -> Check Disconnectionable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a linelet set for a TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 1>
% ll_cand_tmp = [23 26 30 4 1]; 
% for i = 1:3, ll_cand_tmp = [ll_cand_tmp; ll_cand_tmp(end,:) + [4 4 -1 0 0]]; end
% for i = 1:3, ll_cand_tmp = [ll_cand_tmp; ll_cand_tmp(end,:) + [4 4 1 0 0]]; end
% ll_cand_tmp = [ll_cand_tmp; [ll_cand_tmp(4,:) + [4 4 -1 0 0;]]];
% ll_cand_tmp = [ll_cand_tmp; [ll_cand_tmp(4,:) + [8 8 -2 0 0;]]];
% 
% % ll_cand_tmp = [ll_cand_tmp; 22 22 31 1 1];
% % ll_cand_tmp = [ll_cand_tmp; 21 21 32 1 1];
% % ll_cand_tmp = [ll_cand_tmp; 20 20 33 1 1];
% % ll_cand_tmp = [ll_cand_tmp; 19 19 34 1 1];
% % ll_cand_tmp = [ll_cand_tmp; 18 18 35 1 1];
% % ll_cand_tmp = [ll_cand_tmp; 17 17 36 1 1];
% ll_cand = ll_cand_tmp;
% figure; im_tmp = DrawLL(ll_cand_tmp, size_im); imagesc(1-im_tmp); colormap gray
% case 2>
% ll_cand_tmp = [23 44 30 22 1; 45 46 31 2 1; 47 47 32 1 1; 48 48 33 1 1; 49 49 34 1 1]; 
% figure; im_tmp = DrawLL(ll_cand_tmp, size_im); imagesc(im_tmp);

% case 2>
% ll_cand_tmp = [57.0000   70.0000  288.0000   14.0000   -1.6709; 
%    72.0000   80.0000  289.0000    9.0000   -1.6113;
%    82.0000   93.0000  289.0000   12.0000   -1.5931;
%    95.0000  104.0000  290.0000   10.0000   -1.5990];
% figure; im_tmp = DrawLL(ll_cand_tmp, size_im); imagesc(im_tmp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure; im_tmp = DrawLL(ll_cand, size_im); imagesc(im_tmp);
num_ll = size(ll_cand, 1);
[~, idx] = sort(ll_cand(:,1));
ll_cand = ll_cand(idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i) Improve linelet: connect two disjoint linelets that share short (compared to current ones) linelet by removing it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bVisited = false(num_ll,1);
MergeFlag = zeros(num_ll,1);

% Check mergeable linelets
for iCur = 1:num_ll    
    if bVisited(iCur) || MergeFlag(iCur) == -1
        continue; 
    end
    bVisited(iCur) = true;
    
    % find possible merge candidates: one lies within a considerable range (x-value) with the same y-value    
    
%     thres_merge_dist = ceil(max([ll_cand(iCur,4)*.75,1]));
    thres_merge_dist = ceil(max([ll_cand(:,4)*.5]));
    
    % Left direction
    idxNbLeft = find(ll_cand(:,3) == ll_cand(iCur,3) & ll_cand(iCur,1) - ll_cand(:,2) <= thres_merge_dist & ll_cand(iCur,1) - ll_cand(:,2) > 0 & MergeFlag == 0);
    idxNbLeft(idxNbLeft == iCur) = [];    
    idxRemoveLeft = [];
    
    % number of merging candidate should be one
    if length(idxNbLeft) > 1
        [~, idx] = min(ll_cand(iCur,1) - ll_cand(idxNbLeft,2));
        idxNbLeft = idxNbLeft(idx);
    end
    if length(idxNbLeft) == 1
        % Also find candidates, located between two merging ones, which should be removed from the procress above
        idxRemoveLeft = find(ll_cand(:,2) < ll_cand(iCur, 1) & ll_cand(:,1) > ll_cand(idxNbLeft, 2) & MergeFlag == 0);
    end
    
    % Right direction
    idxNbRight = find(ll_cand(:,3) == ll_cand(iCur,3) & ll_cand(:,1) - ll_cand(iCur,2) <= thres_merge_dist & ll_cand(:,1) - ll_cand(iCur,2) > 0 & MergeFlag == 0);
    idxNbRight(idxNbRight == iCur) = [];    
    idxRemoveRight = [];
    
    % number of merging candidate should be one
    if length(idxNbRight) > 1
        [~, idx] = min(ll_cand(idxNbRight,1) - ll_cand(iCur,2));
        idxNbRight = idxNbRight(idx);
    end
    if length(idxNbRight) == 1
        % Also find candidates, located between two merging ones, which should be removed from the procress above
        idxRemoveRight = find(ll_cand(:,2) < ll_cand(idxNbRight, 1) & ll_cand(:,1) > ll_cand(iCur, 2) & MergeFlag == 0);
    end      
    
    % Mark mergining index
    if MergeFlag(iCur) == 0
        MergeFlag([iCur, idxNbLeft, idxNbRight]) = iCur;
        MergeFlag([idxRemoveLeft, idxRemoveRight]) = -1;
    else
        MergeFlag([idxNbLeft, idxNbRight]) = MergeFlag(iCur);
        MergeFlag([idxRemoveLeft, idxRemoveRight]) = -1;
    end
end

%%
% Merge linelet
ll_cand_new = [];
idxMerge = unique(MergeFlag);
if num_ll ~= length(idxMerge)
    for iCur = 1:length(idxMerge)
        if idxMerge(iCur) == -1, continue; end
        
        idx = find(MergeFlag == idxMerge(iCur));
        
        if length(idx) > 1
            ll_cand_new = [ll_cand_new; min(ll_cand(idx,1)), max(ll_cand(idx,2)), ll_cand(idx(1),3), max(ll_cand(idx,2))-min(ll_cand(idx,1))+1, mean(ll_cand(idx,5))];
        else
            ll_cand_new = [ll_cand_new; ll_cand(idx,:)];
        end
    end
else
    ll_cand_new = ll_cand;
end
% figure; im_tmp = DrawLL(ll_cand_new, size_im); imagesc(im_tmp);

% Stitch -- assume that all linelets are aligned according to each linelet's left position in ascending order
ll_cand_new(1:end-1,2) = ll_cand_new(2:end,1) - 1;
ll_cand_new(:,4) = ll_cand_new(:,2) - ll_cand_new(:,1) + 1;
% figure; im_tmp = DrawLL(ll_cand_new, size_im); imagesc(im_tmp);

% num_ll = size(ll_cand_new,1);
% for iCur = 1:num_ll-1
%     % Maybe, I should consider whether stitching, stretching one of the end of a linelet to the other's end. For this
%     % moment, neverthless, only consider left -> right stitching
%     if ll_cand_new(iCur+1,1) - ll_cand_new(iCur, 2) > 1
%         ll_cand_new(iCur,2) = ll_cand_new(iCur_1,1) - 1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ii) Get connection and consistency Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build linelet connection map: non-zero value of conn_map(i,j) means that linelet i and j are connected such that
%                               positive/negative for i is below/above j, and 1/2 for i is left/right to j
[conn_map, ll_cand_new] = Get_ConnMap2(ll_cand_new);

% Build linelet consistency map: 
cons_map = Get_ConsistencyMap(ll_cand_new, conn_map);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iii) Check Disconnectionable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether current linelet set, in particularly the major linelet, has  
% 1) length consistency
% 1-1) a neighbor with different length up to a tolerence, which should be modeled wrt the current linelet
%
% 2) direction consistency
% 2-1) two neighbors at one side,
% 2-2) neighbors at both ends with different direction,
% If so, divide the current set into two (or) three set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from the longest to the shorted ll
ll_cand = ll_cand_new;

% [val_srt, idx_srt] = sort(ll_cand(:,4), 'descend');
% ll_cand = ll_cand(idx_srt,:);
% conn_map = Get_ConnMap(ll_cand);
num_ll = size(ll_cand, 1);
ll_usage = false(num_ll, 1);

for i = 1:num_ll
    if ll_usage(i), continue; end
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!! THIS PART SHOULD BE MODIFED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %thes_len =10;%max(2, ll_cand(i,4)*.25);%2; % this should be modeled properely, e.g., wrt the current linelet
    
    
    % Find connected neighbors
    idx_nb = find(conn_map(i,:) ~= 0);
    if length(idx_nb) < 2, continue; end
    
    % 2-2) a neighbor with different length up to a tolerence, which should be modeled wrt the current linelet
    %     idx_nb_1 = intersect(find( abs(ll_cand(:,4) - ll_cand(i,4)) >= thes_len ), idx_nb);
    %     for j = 1:length(idx_nb_1)
    %         if length(find(conn_map(idx_nb_1(j),:) ~= 0)) > 1
    %             conn_map = DisconnectLinelet(conn_map, i, idx_nb_1(j));
    %         end
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % case 1) direction consistency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1-1) two neighbors at one side or 1-2) neighbors at both ends with different direction,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(idx_nb) - 1                    
        for k = j+1:length(idx_nb)
            % Disconnect one with less proper consistencies: size and trend
            if abs(conn_map(i,idx_nb(j))) == abs(conn_map(i,idx_nb(k)))     % 2-1)
                %if cons_map(idx_nb(j), % Consider only trend at this moment
                len_diff = [ll_cand(idx_nb(j),4) ll_cand(idx_nb(k),4)] - ll_cand(i,4);
                [~, idx_cut] = min(len_diff);
                if idx_cut == 1, idx_cut = idx_nb(j); else idx_cut = idx_nb(k); end
                conn_map = DisconnectLinelet(conn_map, i, idx_cut);
%                 '2-1', j,k
            elseif conn_map(i,idx_nb(j)) * conn_map(i,idx_nb(k)) > 0        % 2-2)
                len_diff = [ll_cand(idx_nb(j),4) ll_cand(idx_nb(k),4)] - ll_cand(i,4);
                [~, idx_cut] = min(len_diff);
                if idx_cut == 1, idx_cut = idx_nb(j); else idx_cut = idx_nb(k); end
                conn_map = DisconnectLinelet(conn_map, i, idx_cut);
%                 '2-2', j,k
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Seperate immediately
[ll_cand_new0, conn_map_new0] = DecomposeConnMap(ll_cand, conn_map);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 2) length consistency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-1) when the sum of the point-to-line distance exceeds the threshold or there is a change of the sign
% of a point-to-line distance
n_cand_set = 1;
ll_cand_new = [];
conn_map_new = [];
for iGroup = 1:length(ll_cand_new0)
    ll_cand = ll_cand_new0{iGroup};
    conn_map = conn_map_new0{iGroup};
    
    %figure; im_tmp = DrawLL(ll_cand, size_im); imagesc(im_tmp);
    
    idx_split = get_line_splitting_location(ll_cand, estimate_lineseg(ll_cand, param.est_line_seg), param.thres_split_dist);
    
    if ~isempty(idx_split)
        for k = 1:length(idx_split)
            if idx_split(k) < size(ll_cand, 1)
                conn_map = DisconnectLinelet(conn_map, idx_split(k), idx_split(k)+1);
            end
        end
        [ll_cand_new1, conn_map_new1] = DecomposeConnMap(ll_cand, conn_map);
        
        for i=1:length(ll_cand_new1)
            ll_cand_new{n_cand_set,1} = ll_cand_new1{i};
            conn_map_new{n_cand_set,1} = conn_map_new1{i};
            n_cand_set = n_cand_set + 1;
        end
    else
        ll_cand_new{n_cand_set,1} = ll_cand;
        conn_map_new{n_cand_set,1} = conn_map;
        n_cand_set = n_cand_set + 1;
    end
end
    
    
    
    
    
    %     % 2-2) a neighbor with different length up to a tolerence, which should be modeled wrt the current linelet
    %     idx_nb_1 = intersect(find( abs(ll_cand(:,4) - ll_cand(i,4)) >= thes_len ), idx_nb);
    %     for j = 1:length(idx_nb_1)
    %         if length(find(conn_map(idx_nb_1(j),:) ~= 0)) > 1
    %             conn_map = DisconnectLinelet(conn_map, i, idx_nb_1(j));
    %         end
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-organize conn_map so as to separate disconnected linelets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ll_cand_new, conn_map_new] = DecomposeConnMap(ll_cand, conn_map);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

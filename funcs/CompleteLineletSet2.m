function [ ll_cand_new, conn_map_new ] = CompleteLineletSet2( ll_cand, param )
%COMPLETELINELETSET2 Summary of this function goes here
%   Detailed explanation goes here

% To do 
% build the connection map -> validate -> check possible merging (improve linelet) -> check seprability
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% PROCESS AT 2-1 and 2-2 SHOULD BE CONSIDERED MORE
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a linelet set for a TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 1>
% ll_cand = [23 26 30 4 1]; 
% for i = 1:3, ll_cand = [ll_cand; ll_cand(end,:) + [4 4 -1 0 0]]; end
% for i = 1:2, ll_cand = [ll_cand; ll_cand(end,:) + [4 4 1 0 0]]; end
% ll_cand = [ll_cand; [ll_cand(4,:) + [4 4 -1 0 0;]]];
% case 2>
% ll_cand = [23 44 30 22 1; 45 46 31 2 1; 47 47 32 1 1; 48 48 33 1 1; 49 49 34 1 1]; 
% figure; im_tmp = DrawLL(ll_cand, size_im); figure; imagesc(im_tmp);

% case 2>
% ll_cand = [57.0000   70.0000  288.0000   14.0000   -1.6709; 
%    72.0000   80.0000  289.0000    9.0000   -1.6113;
%    82.0000   93.0000  289.0000   12.0000   -1.5931;
%    95.0000  104.0000  290.0000   10.0000   -1.5990];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_ll = size(ll_cand, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i) Build the connection map and conistency map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build linelet connection map: non-zero value of conn_map(i,j) means that linelet i and j are connected such that
%                               positive/negative for i is below/above j, and 1/2 for i is left/right to j
conn_map = Get_ConnMap(ll_cand);
cons_map = Get_ConsistencyMap(ll_cand, conn_map);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ii) Validate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iii) Check possible merging (improve linelet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_ll = size(ll_cand, 1);
% ll_extended = false(size(ll_cand,1),1);
% try
%     % When the current linelet has no connected neighbors, maybe has few pixels away.
%     for i = 1:num_ll
%         % Look for possible connection for leftward
%         idx_le = find( (ll_cand(i,1) - ll_cand(:,2) <= param.thres_linelet_merge_dist) & (ll_cand(i,1) - ll_cand(:,2) > 0) &...
%                        abs(ll_cand(:,3) - ll_cand(i,3)) <= 1 & (conn_map(i,:) == 0)' & ~ll_extended);
%                    
%         if ~isempty(idx_le)
%             if ll_cand(idx_le,4) <= ll_cand(i,4)
%                 ll_cand(idx_le, 2) = ll_cand(i,1) - 1;
%                 ll_cand(idx_le, 4) = ll_cand(idx_le, 2) - ll_cand(idx_le, 1) + 1;
%             else
%                 ll_cand(i, 1) = ll_cand(idx_le,2) + 1;
%                 ll_cand(i, 4) = ll_cand(i, 2) - ll_cand(i, 1) + 1;
%             end
%             ll_extended(idx_le) = 1;
%         end
%         
%         % Look for possible connection for rightward
%         idx_re = find( (ll_cand(:,1) - ll_cand(i,2) <= param.thres_linelet_merge_dist) & (ll_cand(:,1) - ll_cand(i,2) > 0) &...
%                        abs(ll_cand(:,3) - ll_cand(i,3)) <= 1 & (conn_map(i,:) == 0)' & ~ll_extended);
%                    
%         if ~isempty(idx_re)
%             if ll_cand(idx_re,4) <= ll_cand(i,4)
%                 ll_cand(idx_re, 1) = ll_cand(i,2) + 1;
%                 ll_cand(idx_re, 4) = ll_cand(idx_re, 2) - ll_cand(idx_re, 1) + 1;
%             else
%                 ll_cand(i, 2) = ll_cand(idx_re, 1) - 1;
%                 ll_cand(i, 4) = ll_cand(i, 2) - ll_cand(i, 1) + 1;
%             end
%             ll_extended(idx_re) = 1;
%         end
%         ll_extended(i) = 1;
%     end
% catch exception
%     fprintf('Err in the first loop of CompleteLineletSet2() with i = %d.\n', i);
%     rethrow(exception);    
% end
% conn_map = Get_ConnMap(ll_cand);
% cons_map = Get_ConsistencyMap(ll_cand, conn_map);

% When the current linelet has length 1 with two neighbors at both sides   
ll_cand_new = [];
ll_merged = false(size(ll_cand,1),1);
for i = 1:num_ll
    idx_nb = find(conn_map(i,:) ~= 0 & ll_cand(i,4) == 1); 
    if length(idx_nb) == 2 && prod(conn_map(i, idx_nb)) > 0 && ~ll_merged(i)
        % Merge neighbors
        idx_l = find(abs(conn_map(i,:)) == 1);
        idx_r = find(abs(conn_map(i,:)) == 2);
        
        ll_cand_new = [ll_cand_new; ll_cand(idx_l, 1), ll_cand(idx_r, 2), ll_cand(idx_l, 3), ll_cand(idx_r, 2) - ll_cand(idx_l, 1)+1, (ll_cand(idx_l,5)+ll_cand(idx_r,5))/2];
        
        ll_merged(idx_nb) = true;                
        %ll_merged(idx_l) = true;
        %ll_merged(idx_r) = true;
    elseif ~ll_merged(i)
        ll_cand_new = [ll_cand_new; ll_cand(i,:)];
    end    
    ll_merged(i) = true;
end


% a model to handle a case when ll becomes [0 0 0 0] should be added, maybe due to merging process?
ll_cand = ll_cand_new;
conn_map = Get_ConnMap(ll_cand);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iv) Check separability
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
[val_srt, idx_srt] = sort(ll_cand(:,4), 'descend');
ll_cand = ll_cand(idx_srt,:);
conn_map = Get_ConnMap(ll_cand);
num_ll = size(ll_cand, 1);
ll_usage = false(num_ll, 1);

for i = 1:num_ll
    if ll_usage(i), continue; end
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!! THIS PART SHOULD BE MODIFED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    thes_len = 4;%max(1, ll_cand(i,4)*.25);%2; % this should be modeled properely, e.g., wrt the current linelet
    
    
    % Find connected neighbors
    idx_nb = find(conn_map(i,:) ~= 0);
    %if length(idx_nb) < 2, continue; end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % case 1) length consistency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1-1) a neighbor with different length up to a tolerence, which should be modeled wrt the current linelet
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx_nb_1 = intersect(find( abs(ll_cand(:,4) - ll_cand(i,4)) > thes_len ), idx_nb);
    for j = 1:length(idx_nb_1)
        if length(find(conn_map(idx_nb_1(j),:) ~= 0)) > 1
            conn_map = DisconnectLinelet(conn_map, i, idx_nb_1(j)); 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % case 2) direction consistency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2-1) two neighbors at one side or 2-2) neighbors at both ends with different direction,
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
            elseif conn_map(i,idx_nb(j)) * conn_map(i,idx_nb(k)) > 0        % 2-2)
                len_diff = [ll_cand(idx_nb(j),4) ll_cand(idx_nb(k),4)] - ll_cand(i,4);
                [~, idx_cut] = min(len_diff);
                if idx_cut == 1, idx_cut = idx_nb(j); else idx_cut = idx_nb(k); end
                conn_map = DisconnectLinelet(conn_map, i, idx_cut);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-organize conn_map so as to separate disconnected linelets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ll_cand_new, conn_map_new] = DecomposeConnMap(ll_cand, conn_map);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

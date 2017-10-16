function [LL_Set_ret, valid_idx_map_LL_Set, LL_ret, valid_idx_map_LL] = PruneShortLinletSet(LL_in, thres_linelet_len, size_im)
%LL_Set_ver1, lbMap_LL_Set_ver1, LL_ver1, lbMap_LL_ver1

% LL_in = LL_8_ver2;
% thres_linelet_len = param.thres_linelet_len;

% Prune short linelet set
compo_len = cellfun(@length, LL_in.PixelIdxList); 
ind = find(compo_len > thres_linelet_len);
LL_Set_ret = cell(length(ind),1); 
[LL_Set_ret{:}] = LL_in.PixelIdxList{ind};
ind_ver = cell2mat(LL_Set_ret);

% Polish linelet set 
valid_idx_map_LL_Set = zeros(size_im);    
valid_idx_map_LL_Set(ind_ver) = 1;
LL_Set_ret = bwconncomp(valid_idx_map_LL_Set, 8);
valid_idx_map_LL_Set = bwlabel(valid_idx_map_LL_Set, 8);
%figure; imagesc(valid_idx_map)

LL_ret = bwconncomp(valid_idx_map_LL_Set>0, 4);
valid_idx_map_LL = bwlabel(valid_idx_map_LL_Set>0, 4);
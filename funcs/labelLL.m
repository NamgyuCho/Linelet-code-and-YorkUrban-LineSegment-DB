function [ ll_Hor, ll_Ver, ll_Diag, map_ll_Hor, map_ll_Ver, map_ll_Diag ] = labelLL( memX, memY, lsc_map, bDisplay )
%LABELLL Summary of this function goes here
%   Detailed explanation goes here
% memX = memX0, memY = memY0;

size_im = size(memX);
pixel_used = zeros(size_im);

%% Horizontal linelet
bw = uint8(memY == 1);

ll_Hor = []; 
map_ll_Hor = zeros(size_im);
num_ll_Hor = 1;

for i = 1:size_im(1)
    for j = 1:size_im(2)
        if bw(i,j) ~= 0 && map_ll_Hor(i,j) == 0
            if j > 1 && map_ll_Hor(i,j-1) ~= 0 
                map_ll_Hor(i,j) = map_ll_Hor(i,j-1);
            else
                map_ll_Hor(i,j) = num_ll_Hor;
                num_ll_Hor = num_ll_Hor + 1;
            end
        end
    end
    vec_tmp = map_ll_Hor(i,:);
    [vUnq, nUnq] = count_unique(vec_tmp);
    for k = 1:length(vUnq)
%         if nUnq(k) <= 1
%             vec_tmp(vec_tmp == vUnq(k)) = 0;
%         else
            if vUnq(k) ~= 0
                idx = find(vec_tmp == vUnq(k));
                lsc_id = lsc_map( sub2ind(size_im, ones(length(idx),1)*i, idx(:)) );
                %ll_Hor = [ll_Hor; idx(1) i idx(end) i length(idx) lsc_map(i, idx(1))];
                ll_Hor = [ll_Hor; idx(1) i idx(end) i length(idx) median(lsc_id)];
            end
%         end
    end
    map_ll_Hor(i,:) = vec_tmp(:);
end

if ~isempty(ll_Hor)
    [~, idx_srt] = sort(ll_Hor(:,5), 'descend');
    ll_Hor = ll_Hor(idx_srt,:);
end
    
    
%if bDisplay, figure; imagesc(map_ll_Hor); axis image; axis off; end

%% Vertical linelet        
bw = uint8(memX == 1);
ll_Ver = [];
map_ll_Ver = zeros(size_im);
num_ll_Ver = 1;

for j = 1:size_im(2)
    for i = 1:size_im(1)
        if bw(i,j) ~= 0 && map_ll_Ver(i,j) == 0
            if i > 1 && map_ll_Ver(i-1,j) ~= 0
                map_ll_Ver(i,j) = map_ll_Ver(i-1,j);
            else
                map_ll_Ver(i,j) = num_ll_Ver;
                num_ll_Ver = num_ll_Ver + 1;
            end
        end
    end
    vec_tmp = map_ll_Ver(:,j);
    [vUnq, nUnq] = count_unique(vec_tmp);
    for k = 1:length(vUnq)
%         if nUnq(k) <= 1
%             vec_tmp(vec_tmp == vUnq(k)) = 0;
%         else
            if vUnq(k) ~= 0
                idx = find(vec_tmp == vUnq(k));
%                 ll_Ver = [ll_Ver; j idx(1) j idx(end) length(idx) lsc_map(idx(1), j)];
                lsc_id = lsc_map( sub2ind(size_im, idx(:), ones(length(idx),1)*j) );
                ll_Ver = [ll_Ver; j idx(1) j idx(end) length(idx) median(lsc_id)];
            end
%         end
    end
    map_ll_Ver(:,j) = vec_tmp(:);
end

if ~isempty(ll_Ver)
    [~, idx_srt] = sort(ll_Ver(:,5), 'descend');
    ll_Ver = ll_Ver(idx_srt,:);
end

%if bDisplay, figure; imagesc(map_ll_Ver); axis image; axis off; end

ll_Diag = [];
map_ll_Diag = zeros(size_im);
end


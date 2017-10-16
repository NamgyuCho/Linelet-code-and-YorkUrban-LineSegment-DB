function [ PtSet ] = transformLinelet2PtSet( ll_cand, bHorMajor, size_im )
    %TRANSFORMLINELET Summary of this function goes here
    %   Detailed explanation goes here    
    PtSet = [];
    
    if iscell(ll_cand)
        for i_c = 1:size(ll_cand, 1)
            tmp_val = ll_cand{i_c};
            idx = [];
            
            for i_t = 1:size(tmp_val,1)
                if bHorMajor
%                     idx = sub2ind(size_im, repmat(ll_cand(i,3), 1, ll_cand(i,4)), ll_cand(i,1):ll_cand(i,2));
                    idx_tmp_val = sub2ind(size_im, repmat(tmp_val(i_t,3), 1, tmp_val(i_t,4)), tmp_val(i_t,1):tmp_val(i_t,2));                    
                else
%                     if any(tmp_val(:,1) > size_im(1)) || any(tmp_val(:,2) > size_im(1)) || any(tmp_val(:,3) > size_im(2))
%                         disp('hello');
%                     end
                    idx_tmp_val = sub2ind(size_im, tmp_val(i_t,1):tmp_val(i_t,2), repmat(tmp_val(i_t,3), 1, tmp_val(i_t,4)));                    
                end
                idx = [idx, idx_tmp_val];
            end
            PtSet{i_c,1} = idx;
        end
    else
%         if bHorMajor
%             tmp_val = [tmp_val(:,3) tmp_val(:,1) tmp_val(:,3) tmp_val(:,2)];
%         else
%             tmp_val = [tmp_val(:,1) tmp_val(:,3) tmp_val(:,2) tmp_val(:,3)];
%         end
%         PtSet = tmp_val;
        idx = [];
        tmp_val = ll_cand;
        for i_t = 1:size(tmp_val,1)
            if bHorMajor
                idx_tmp_val = sub2ind(size_im, repmat(tmp_val(i_t,3), 1, tmp_val(i_t,4)), tmp_val(i_t,1):tmp_val(i_t,2));
            else
                idx_tmp_val = sub2ind(size_im, tmp_val(i_t,1):tmp_val(i_t,2), repmat(tmp_val(i_t,3), 1, tmp_val(i_t,4)));                
            end
            idx = [idx, idx_tmp_val];
        end
        PtSet = idx;
    end        
end


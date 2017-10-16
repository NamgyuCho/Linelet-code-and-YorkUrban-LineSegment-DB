function [ gt_covered, idx_valid, pd_covered ] = line_area_intersection( l_src, l_tar )
    %LINE_AREA_INTERSECTION Calculating overlapping area (non-overlapping between targets) between two line segments wrt l_src (source)
    
    % l_src = line_gnd(i_gnd,:);
    % l_tar = line_est(idx_cand, :);
        
    idx_valid = false(size(l_tar,1),1);
    
    
    gt_covered = zeros(size(l_tar,1),1);
    pd_covered = zeros(size(l_tar,1),1);
    
    
    % Project the source to the source coordinate
    vec_base = l_src(3:4) - l_src(1:2); 
    vec_base = vec_base / norm(vec_base);
    
    vec_src = [l_src(1:2) - l_src(1:2); l_src(3:4) - l_src(1:2)] * vec_base';
    
    if vec_src(1) > vec_src(2), [vec_src(1), vec_src(2)] = swap_vec(vec_src(1), vec_src(2)); end
    
    vec_tar = [(l_tar(:,1:2) - repmat(l_src(1:2), size(l_tar,1),1)) * vec_base'...
        (l_tar(:,3:4) - repmat(l_src(1:2), size(l_tar,1),1)) * vec_base'];
    
    for i = 1:size(l_tar, 1)
        if vec_tar(i,1) > vec_tar(i,2)
            [vec_tar(i,1), vec_tar(i,2)] = swap_vec(vec_tar(i,1), vec_tar(i,2));
        end
    end
    
    % clip left area
    idx = find(vec_tar < 0); vec_tar(idx) = 0;
    
    % clip right area
    idx = find(vec_tar > max(vec_src)); vec_tar(idx) = max(vec_src);
    
    
    for k = 1:size(l_tar, 1)
        % Project the target to the source's coordinate
%         vec_tar = [l_tar(k,1:2) - l_src(1:2); l_tar(k,3:4) - l_src(1:2)] * vec_base';
%         
%         if vec_tar(1) > vec_tar(2)
%             [vec_tar(1), vec_tar(2)] = swap_vec(vec_tar(1), vec_tar(2));
%         end
        
        bValid = true;
        
        if vec_tar(k,1) >= vec_src(2) || vec_tar(k,2) <= vec_src(1)
            % case 1
            % tar:                  *----------*
            % src: *---------*
            % case 2
            % tar: *----------*
            % src:                *---------*
            gt_covered(k) = 0;
            idx_valid(k) = false;
            pd_covered(k) = 0;
            bValid = false;
        elseif vec_tar(k,1) <= vec_src(1) && vec_tar(k,2) >= vec_src(1) && vec_tar(k,2) <= vec_src(2)
            % case 3
            % tar: *-------*
            % src:    *---------*
            vec_tar(k,1) = vec_src(1);
        elseif vec_tar(k,2) >= vec_src(2) && vec_tar(k,1) >= vec_src(1) && vec_tar(k,1) <= vec_src(2)
            % case 4
            % tar:      *----------*
            % src: *---------*
            vec_tar(k,2) = vec_src(2);
        elseif vec_tar(k,1) <= vec_src(1) && vec_tar(k,2) >= vec_src(2)
            % case 5
            % tar: *-------------------*
            % src:     *---------*
            vec_tar(k,:) = vec_src;
        end
        
        %         idx = abs(vec_tar(1) - vec_tar(2)) ./ l_tar(7) > .5;
        if bValid
            idx_valid(k) = true;
            pd_covered(k) = abs(vec_tar(k,1) - vec_tar(k,2)); % area covered by target
        end
    end
    
    % get purely overlapping area -- covered ground truth area
    [~, idx] = sort(vec_tar(:,1));
    idx_valid = idx_valid(idx);
    vec_tar = vec_tar(idx,:);
    for i1 = 1:size(vec_tar,1)-1
        if ~idx_valid(i1), continue; end
        
        for i2 = i1+1:size(vec_tar,1)
            if ~idx_valid(i2), continue; end
            
            % tar: *------*
            % src:   *----------*
            if vec_tar(i2,1) <= vec_tar(i1,1) && vec_tar(i2,2) <= vec_tar(i1,2) && vec_tar(i2,2) >= vec_tar(i1,1)
                vec_tar(i2,2) = vec_tar(i1,1);
            end
            
            % tar:     *--------*
            % src: *----------*
            if vec_tar(i2,1) >= vec_tar(i1,1) && vec_tar(i2,2) >= vec_tar(i1,2) && vec_tar(i2,1) <= vec_tar(i1,2)
                vec_tar(i2,1) = vec_tar(i1,2);
            end
            
            % tar:     *----*
            % src: *----------*
            if vec_tar(i2,1) >= vec_tar(i1,1) && vec_tar(i2,2) <= vec_tar(i1,2)
                vec_tar(i2,:) = 0;
            end
        end        
    end
    
    for k = 1:size(vec_tar, 1)
        gt_covered(k) = abs(vec_tar(k,1) - vec_tar(k,2));
    end    
    
    if 0
    if size(l_tar, 1) == 1
        
        % Project the target to the source's coordinate
        vec_tar = [l_tar(1:2) - l_src(1:2); l_tar(3:4) - l_src(1:2)] * vec_base';
        
        if vec_tar(1) > vec_tar(2)
            [vec_tar(1), vec_tar(2)] = swap_vec(vec_tar(1), vec_tar(2));
        end
        
%         idx = vec_tar < 0; vec_tar(idx) = 0;
        
        if vec_tar(1) >= vec_src(2) || vec_tar(2) <= vec_src(1)
            % case 1
            % tar:                  *----------*
            % src: *---------*
            % case 2
            % tar: *----------*
            % src:                *---------*
            area = 0;
            idx_valid = false;
            area_redundant = 0; 
            return;
        elseif vec_tar(1) <= vec_src(1) && vec_tar(2) >= vec_src(1) && vec_tar(2) <= vec_src(2)
            % case 3
            % tar: *-------*
            % src:    *---------*
            vec_tar(1) = vec_src(1);
        elseif vec_tar(2) >= vec_src(2) && vec_tar(1) >= vec_src(1) && vec_tar(1) <= vec_src(2)
            % case 4
            % tar:      *----------*
            % src: *---------*
            vec_tar(2) = vec_src(2);
        elseif vec_tar(1) <= vec_src(1) && vec_tar(2) >= vec_src(2)
            % case 5
            % tar: *-------------------*
            % src:     *---------*
            vec_tar = vec_src;
        end     
        
%         idx = abs(vec_tar(1) - vec_tar(2)) ./ l_tar(7) > .5;
        idx_valid(idx) = true;

        area = abs(vec_tar(1) - vec_tar(2));
    else
        vec_tar = [(l_tar(:,1:2) - repmat(l_src(1:2), size(l_tar,1),1)) * vec_base'...
                   (l_tar(:,3:4) - repmat(l_src(1:2), size(l_tar,1),1)) * vec_base'];
               
       for i = 1:size(l_tar, 1)
           if vec_tar(i,1) > vec_tar(i,2)
               [vec_tar(i,1), vec_tar(i,2)] = swap_vec(vec_tar(i,1), vec_tar(i,2));
           end
       end
        
        % clip left area
%         idx = find(vec_tar < 0); vec_tar(idx) = 0;
        
        % clip right area
%         idx = find(vec_tar > max(vec_src)); vec_tar(idx) = max(vec_src);
        
        % only consider targets with the coverage area more than a half 
%         idx = abs(vec_tar(:,1) - vec_tar(:,2)) ./ l_tar(:,7) > .5;
        idx_valid(idx) = true;
        %vec_tar(~idx,:) = [];
        
        % align results
        for i1 = 1:size(vec_tar,1)
            if ~idx(i1), continue; end
            
            if vec_tar(i1,1) > vec_tar(i1,2)
                tmp=vec_tar(i1,2);
                vec_tar(i1,2) = vec_tar(i1,1);
                vec_tar(i1,1) = tmp;
            end
        end
        
        % Area wrt redundancy
        area_redundant(idx) = abs(vec_tar(idx,1) - vec_tar(idx,2));
        
        
        % get purely overlapping area -- remove overlap between targets
        for i1 = 1:size(vec_tar,1)
            if ~idx(i1), continue; end
            
            for i2 = i1+1:size(vec_tar,1)
                if ~idx(i2), continue; end
                
                if vec_tar(i2,2) >= vec_tar(i1,1) && vec_tar(i2,2) <= vec_tar(i1,2)
                    vec_tar(i2,2) = vec_tar(i1,1);
                end
                
                if vec_tar(i2,1) >= vec_tar(i1,1) && vec_tar(i2,1) <= vec_tar(i1,2)
                    vec_tar(i2,1) = vec_tar(i1,2);
                end
                
                if vec_tar(i2,1) >= vec_tar(i1,1) && vec_tar(i2,2) <= vec_tar(i1,2)
                    vec_tar(i2,:) = 0;
                end
            end
        end
        area(idx) = abs(vec_tar(idx,1) - vec_tar(idx,2));
    end
    %vac_src = [l_src(1:2) - l_src(1:2); l_src(3:4) - l_src(1:2)] * vec_base';
    end
end


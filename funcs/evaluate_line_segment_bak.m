function [ precision, recall, accuracy, bUsed_est ] = evaluate_line_segment( line_est, line_gnd, param )
%EVALUATE_LSC Summary of this function goes here
%   Detailed explanation goes here
% line_est = line_ht;
% --------------------------------------------------------------------------------------------------
% Initialization
% --------------------------------------------------------------------------------------------------
%%
b_plot = false;
precision = 0;
recall = 0;
accuracy = 0;

% Initialize retrieval numbers -- 1st row: pixelwise, 2nd row: line segment wise
tp_area_est = 0; tp_area_gnd = 0; 
tn_area_est = 0; tn_area_gnd = 0;

tp_inst_est = 0; tp_inst_gnd = 0;
tn_inst_est = 0; tn_inst_gnd = 0;

fn_area_est = 0; fn_area_gnd = 0;
fp_area_est = 0; fp_area_gnd = 0;

fn_inst_est = 0; fn_inst_gnd = 0;
fp_inst_est = 0; fp_inst_gnd = 0;
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% Convert a line segment to a set of indices -- gnd (x1, y1, x2, y2, center_x, center_y, length, angle)
num_gnd = size(line_gnd, 1);
idx_gnd = cell(num_gnd, 1);
num_total_pixel_gnd = 0;
bSteep_gnd = zeros(num_gnd,1);
for i = 1:num_gnd
    [tx, ty, bSteep] = bresenham(line_gnd(i,1), line_gnd(i,2), line_gnd(i,3), line_gnd(i,4));
    idx_gnd{i} = [tx ty];
    bSteep_gnd(i) = bSteep;
    num_total_pixel_gnd = num_total_pixel_gnd + length(tx);
end

% Convert a line segment to a set of indices -- est (x1, y1, x2, y2, center_x, center_y, length, angle)
num_est = size(line_est, 1);
idx_est = cell(num_est, 1);
num_total_pixel_est = 0;
bSteep_est = zeros(num_est,1);
bUsed_est = zeros(num_est, 1);
for i = 1:num_est
    [tx, ty, bSteep] = bresenham(line_est(i,1), line_est(i,2), line_est(i,3), line_est(i,4));
    
    if bSteep
        idx_est{i} = [tx ty];%; tx+1 ty; tx+2 ty; tx-1 ty; tx-2 ty];
    else
        idx_est{i} = [tx ty];%; tx ty+1; tx ty+2; tx ty-1; tx ty-2];
    end
    num_total_pixel_est = num_total_pixel_est + length(tx);
end
% --------------------------------------------------------------------------------------------------

% Pairwise center distance between  of gnd and est
pdist_gnd_est = pdist2(line_gnd(:,5:6), line_est(:,5:6));
valid_inst_est = false(1, num_est);
%num_retrieved_pixel_ret = 0;

%i_gnd = 107; %salient one
if b_plot,  figure(11); clf; end
for i_gnd = 1:num_gnd
%     if line_gnd(i_gnd,7) < 100, continue; end
%     i_gnd
    if b_plot
        imshow(im); hold on;
        plot(line_gnd(i_gnd, 1:2:3), line_gnd(i_gnd, 2:2:4), 'b-', 'linewidth', 3);
        
        
%         figure(12); clf; imshow(im); hold on;
%         idx = find(line_est(:,7) > 100);        
%         for i=1:length(idx)
%             idx(i)
%             plot(line_est(idx(i), 1:2:3), line_est(idx(i), 2:2:4), 'r-', 'linewidth', 2);
%            waitforbuttonpress
%         end
    end
    
    try
        % If the perpendicular distance and angle difference is less than
        % threshold, take it as true positive
        d = line_gnd(i_gnd, 3:4) - line_gnd(i_gnd, 1:2);
        d = d/norm(d); % direction vector of the line
        n = [-d(2),d(1)]; % unit normal vector of the line
        
        % Perpendicular distance -- line normal * distance vector between one of original line point and estimated one's
        perpendicular_dist_gnd2est = n*( (line_est(:,1:2)+line_est(:,3:4))/2 - repmat(line_gnd(i_gnd, 1:2), num_est,1))';
        idx_cand = find(abs(perpendicular_dist_gnd2est) <= param.thres_dist);
%         for k = 1:length(idx_cand), plot(line_est(idx_cand(k), 1:2:3), line_est(idx_cand(k), 2:2:4), 'r-', 'linewidth', 2); end
                
        % Angle difference
%         idx_ang1 = bAngleAligned( line_gnd(i_gnd, 8), line_est(:, 8), param.thres_ang )'; % Gradient angle            
%         idx_ang2 = bAngleAligned( line_est(:, 8), line_gnd(i_gnd, 8), param.thres_ang )'; % Gradient angle            
        %ang_diff_est = angle_diff(line_est(:, 8), repmat(line_gnd(i_gnd, 8), num_est, 1))';
        idx_ang = bAngleAligned( line_gnd(i_gnd, 8), line_est(idx_cand, 8), param.thres_ang )'; % Gradient angle            
        idx_cand = idx_cand(idx_ang);
%         for k = 1:length(idx_cand), plot(line_est(idx_cand(k), 1:2:3), line_est(idx_cand(k), 2:2:4), 'g-', 'linewidth', 2); end
        
        % Overlapping instances
%         if bSteep_gnd(i_gnd)
%             idx_overlap1 = (line_est(:,4) > line_gnd(i_gnd,2) & line_est(:,2) < line_gnd(i_gnd, 4))';
%             idx_overlap2 = (line_est(:,4) < line_gnd(i_gnd,2) & line_est(:,2) > line_gnd(i_gnd, 4))';
%         else
%             idx_overlap1 = (line_est(:,3) > line_gnd(i_gnd,1) & line_est(:,1) < line_gnd(i_gnd, 3))';
%             idx_overlap2 = (line_est(:,3) < line_gnd(i_gnd,1) & line_est(:,1) > line_gnd(i_gnd, 3))';
%         end
                
        % final candidates
%         idx_cand_est = find(abs(perpendicular_dist_gnd2est) <= param.thres_dist & (idx_ang1) & (idx_overlap1 | idx_overlap2) & ~bUsed_est');
        %idx_cand = find(abs(perpendicular_dist_gnd2est) <= param.thres_dist & (idx_ang1) & ~bUsed_est');
                        
        if isempty(idx_cand)
            % False negative
            fn_area_gnd = fn_area_gnd + line_gnd(i_gnd, 7);
            fn_inst_gnd = fn_inst_gnd + 1;
        else
            % True positive
            [len_overlap, idx_valid] = line_area_intersection(line_gnd(i_gnd,:), line_est(idx_cand,:));        
            
            if ~sum(idx_valid)
                fn_area_gnd = fn_area_gnd + line_gnd(i_gnd, 7);
                fn_inst_gnd = fn_inst_gnd + 1;
                continue; 
            end            
                    
            tp_area_est = tp_area_est + sum(len_overlap(idx_valid));
            tp_area_gnd = tp_area_gnd + line_gnd(i_gnd, 7);
            
            if sum(len_overlap(idx_valid)) / line_gnd(i_gnd, 7) >= param.thres_length_ratio
                tp_inst_est = tp_inst_est + length(idx_cand);
                tp_inst_gnd = tp_inst_gnd + 1;
                %fprintf('valid instance\n');
            end
            
            bUsed_est(idx_cand(idx_valid)) = true;
            
            %%
%             num_retrieved_tmp = 0;
%             
%             % Calculate overlaping area between the ground truth and true positive estimates
%             len_overlap = zeros(1, length(idx_cand_est));
%             for i_cand = 1:length(idx_cand_est)                
%                 len_overlap(i_cand) = line_area_intersection(line_gnd(i_gnd,1:4), line_est(idx_cand_est(i_cand),1:4));
%                 
%                 if (len_overlap(i_cand) / line_est(idx_cand_est(i_cand),7)) >= .5
%                     bUsed_est(idx_cand_est(i_cand)) = true;
%                     
%                 end
%                 
% %                 if ~bSteep_gnd(i_gnd) % when the ground truth is not a steep one (more horizontal)
% %                     if     line_est(idx_cand_est(i_cand),1) < line_gnd(i_gnd, 1) && line_est(idx_cand_est(i_cand),3) < line_gnd(i_gnd, 3)
% %                         % case 1 -- overlap one the leftside
% %                         %      ================== (gnd)
% %                         %   ----------            (est)
% %                         len_overlap(i_cand) = line_est(idx_cand_est(i_cand),3) - line_gnd(i_gnd, 1) + 1;
% %                     elseif line_est(idx_cand_est(i_cand),1) >= line_gnd(i_gnd, 1) && line_est(idx_cand_est(i_cand),3) <= line_gnd(i_gnd, 3)
% %                         % case 2 -- overlap withing the ground truth
% %                         %      ==================
% %                         %          ----------
% %                         len_overlap(i_cand) = abs( line_est(idx_cand_est(i_cand),1) - line_est(idx_cand_est(i_cand),3)) + 1;
% %                     elseif line_est(idx_cand_est(i_cand),1) >= line_gnd(i_gnd, 1) && line_est(idx_cand_est(i_cand),3) >= line_gnd(i_gnd, 3)
% %                         % case 3 -- overlap one the rightside
% %                         %      ==================
% %                         %                  ----------
% %                         len_overlap(i_cand) = line_gnd(i_gnd, 3) - line_est(idx_cand_est(i_cand),1) + 1;
% %                     elseif line_est(idx_cand_est(i_cand),1) < line_gnd(i_gnd, 1) && line_est(idx_cand_est(i_cand),3) > line_gnd(i_gnd, 3)
% %                         % case 4 -- overlap within the estimate
% %                         %      ==================
% %                         %   ------------------------
% %                         len_overlap(i_cand) = abs(line_gnd(i_gnd, 3) - line_gnd(i_gnd, 1)) + 1;
% %                     end
% %                 else % when the ground truth is a steep one (more vertical)
% %                     if     line_est(idx_cand_est(i_cand),2) < line_gnd(i_gnd, 2) && line_est(idx_cand_est(i_cand),4) < line_gnd(i_gnd, 4)
% %                         len_overlap(i_cand) = line_est(idx_cand_est(i_cand),4) - line_gnd(i_gnd, 2) + 1;
% %                     elseif line_est(idx_cand_est(i_cand),2) >= line_gnd(i_gnd, 2) && line_est(idx_cand_est(i_cand),4) <= line_gnd(i_gnd, 4)
% %                         len_overlap(i_cand) = abs( line_est(idx_cand_est(i_cand),2) - line_est(idx_cand_est(i_cand),2)) + 1;
% %                     elseif line_est(idx_cand_est(i_cand),2) >= line_gnd(i_gnd, 2) && line_est(idx_cand_est(i_cand),4) >= line_gnd(i_gnd, 4)
% %                         len_overlap(i_cand) = line_gnd(i_gnd, 4) - line_est(idx_cand_est(i_cand),2) + 1;
% %                     elseif line_est(idx_cand_est(i_cand),2) < line_gnd(i_gnd, 2) && line_est(idx_cand_est(i_cand),4) > line_gnd(i_gnd, 4)
% %                         len_overlap(i_cand) = abs(line_gnd(i_gnd, 4) - line_gnd(i_gnd, 2)) + 1;
% %                     end                    
% %                 end
%                 
% %                 [~, idx_int, ~] = intersect(idx_gnd{i_gnd}, idx_est{idx_cand_est(i_cand)}, 'rows');
% %                 idx_gnd{i_gnd}(idx_int,:) = [];
% %                 num_retrieved_tmp = num_retrieved_tmp + length(idx_int);                
%             end
%             
% %             % Calculate overlapping area between candidates - not useful anymore
% %             len_cand_overlap = [];
% %             for i1 = 1:length(idx_cand_est)
% %                 for i2 = i1:length(idx_cand_est)
% %                     if i1 ~= i2
% %                         len_cand_overlap = [len_cand_overlap; line_area_intersection(line_est(idx_cand_est(i1),1:4), line_est(idx_cand_est(i2),1:4))];
% %                     end
% %                 end
% %             end
%             
%             %num_retrieved_pixel_ret = num_retrieved_pixel_ret + num_retrieved_tmp;
%             tp_area_ratio = tp_area_ratio + num_retrieved_tmp;
%             
%             if num_retrieved_tmp / length(idx_gnd{i_gnd}) > param.thres_length_area_ratio
%                 tp_inst_est = tp_inst_est + 1;
%             end
%                         
%             % This should be changed
%             valid_inst_est(idx_cand_est) = true;    
            %%
            if b_plot, 
                for k = 1:length(idx_valid)
                    if idx_valid(k) 
                        plot(line_est(idx_cand(k), 1:2:3), line_est(idx_cand(k), 2:2:4), 'y-', 'linewidth', 2); 
                    end
                end
                waitforbuttonpress
            end
        end        
        %if b_plot, waitforbuttonpress; end
    catch err
        % Give more information of the error
        rethrow(err);
    end
end

fp_area_est = sum(line_gnd(:,7)) - tp_area_est;
fp_inst_est = num_est - tp_inst_est;


precision_area_ratio = tp_area_est / (tp_area_est + fp_area_est);
recall_area_ratio = tp_area_est / (tp_area_est + fn_area_est);
accuracy_area_ratio = (tp_area_est + tn_area_est) ./ (tp_area_est + tn_area_est + fp_area_est + fn_area_est);

precision_inst_est = tp_inst_est / num_est;%(tp_inst_est + fp_inst_est);
%recall_inst_est = tp_inst_est / num_gnd;%(tp_inst_est + fn_inst_est);
recall_inst_gnd = tp_inst_gnd / num_gnd;
accuracy_inst_est = (tp_inst_est + tn_inst_est) ./ (tp_inst_est + tn_inst_est + fp_inst_est + fn_inst_est);

precision = [precision_area_ratio precision_inst_est];
recall = [recall_area_ratio recall_inst_est];
accuracy = [accuracy_area_ratio accuracy_inst_est];

end


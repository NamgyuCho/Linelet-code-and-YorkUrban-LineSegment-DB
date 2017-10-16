function [ precision, recall, iou ] = evaluate_line_segment( line_est, line_gnd, eval_param)
%EVALUATE_LSC Summary of this function goes here
%   Line segment instance should be in a form (x1, y1, x2, y2, center_x, center_y, length, angle)

%%
b_plot = false;
precision = 0;
recall = 0;

% Initialize retrieval numbers -- 1st row: pixelwise, 2nd row: line segment wise
tp_area_est = 0; tp_area_gnd = 0; 
tn_area_est = 0; tn_area_gnd = 0;

tp_inst_est = 0; tp_inst_gnd = 0;
tn_inst_est = 0; tn_inst_gnd = 0;

fn_area_est = 0; fn_area_gnd = 0;
fp_area_est = 0; fp_area_gnd = 0;

fn_inst_est = 0; fn_inst_gnd = 0;
fp_inst_est = 0; fp_inst_gnd = 0;

tp_iou = 0;
fp_iou = 0;
fn_iou = 0;
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% Convert a line segment to a set of indices -- gnd (x1, y1, x2, y2, center_x, center_y, length, angle)
num_gnd = size(line_gnd, 1);

% Convert a line segment to a set of indices -- est (x1, y1, x2, y2, center_x, center_y, length, angle)
num_est = size(line_est, 1);
idx_est = cell(num_est, 1);
num_total_pixel_est = 0;
bSteep_est = zeros(num_est,1);
% --------------------------------------------------------------------------------------------------


%%
for i_gnd = 1:num_gnd
    try
        % If the perpendicular distance and angle difference is less than
        % threshold, take it as true positive
        d = line_gnd(i_gnd, 3:4) - line_gnd(i_gnd, 1:2);
        d = d/norm(d); % direction vector of the line
        n = [-d(2),d(1)]; % unit normal vector of the line
        
        % line structure: (x1, y1, x2, y2, center_x, center_y, length, angle)
        idx_perpd = GetPerpDist(line_gnd(i_gnd, 5:6), line_est(:, 5:6), line_gnd(i_gnd, 8))' <= eval_param.thres_dist;
        idx_ang = bAngleAligned(line_gnd(i_gnd, 8), line_est(:, 8), eval_param.thres_ang);
        idx_cand = find(idx_perpd & idx_ang);
        
        if isempty(idx_cand)
            % False negative
            fn_area_gnd = fn_area_gnd + line_gnd(i_gnd, 7);
            fn_inst_gnd = fn_inst_gnd + 1;
            fn_iou = fn_iou + line_gnd(i_gnd, 7);
        else
            % True positive
            [gt_covered, idx_valid, pd_covered] = line_area_intersection(line_gnd(i_gnd,:), line_est(idx_cand,:));
            
            if ~sum(idx_valid)
                fn_area_gnd = fn_area_gnd + line_gnd(i_gnd, 7);
                fn_inst_gnd = fn_inst_gnd + 1;
                fn_iou = fn_iou + line_gnd(i_gnd, 7);
                continue; 
            end            
                            
            if (sum(gt_covered(idx_valid)) / line_gnd(i_gnd, 7)) >= eval_param.thres_length_ratio
                tp_area_est = tp_area_est + sum(pd_covered(idx_valid));
                tp_area_gnd = tp_area_gnd + sum(gt_covered(idx_valid));
                tp_inst_est = tp_inst_est + sum(idx_valid);
                tp_inst_gnd = tp_inst_gnd + 1;
                
                tp_iou = tp_iou + sum(pd_covered(idx_valid));
                fp_iou = fp_iou + sum(line_est(idx_cand(idx_valid),7)) - sum(pd_covered(idx_valid));                
            end            
        end    
    catch err
        % Give more information of the error
        fprintf('error at evaluate_line_segment(), i_gnd: %d.\n', i_gnd);
        rethrow(err);
    end    
end

%%
fp_area_est = sum(line_est(:,7)) - tp_area_est;
fp_inst_est = num_est - tp_inst_est;

precision_area_est = tp_area_est / sum(line_est(:,7));
recall_area_gnd = tp_area_gnd / sum(line_gnd(:,7));

precision = precision_area_est;
recall = recall_area_gnd;
iou = tp_iou / (tp_iou + fp_iou + fn_iou);

end
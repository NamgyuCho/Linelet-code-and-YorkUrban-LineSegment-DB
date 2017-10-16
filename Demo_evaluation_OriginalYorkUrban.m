clear; close all

dir_db = './DB/YorkUrbanDB';
addpath(genpath('./toolbox/'));
load([dir_db '/our_annotation/Image_ID_List.mat']); % We get Image_ID_List
num_im = size(Image_ID_List, 1);

name_method = {'\fontname{times}HT', '\fontname{times}PPHT', '\fontname{times}LSD', '\fontname{times}LSDi',...
    '\fontname{times}EDLine', '\fontname{times}Linelet_H', '\fontname{times}Linelet'};
method2test = {'ht', 'htp', 'lsd', 'lsdi', 'ed', 'linelet', 'linelet'};
num_method = length(method2test);

dir_method = {'ht', 'ppht', 'lsd', 'lsdi', 'edline', 'proposed_h', 'proposed'};
var_method = {'line_ht', 'line_ht', 'line_lsd', 'line_lsd', 'line_ed', 'line_own', 'line_own'};

stats(num_method,1) = struct('prec', [], 'rec', [], 'iou', []);

bSaveResult = false;

% True positive conditions
eval_param.thres_dist = 1;
eval_param.thres_ang = pi*5/180;
eval_param.thres_length_ratio = .75; % .75 (Con1) or .5 (Con2)

for i_im = 1 : num_im
    try
        str_im = sprintf('%s/%s/%s.jpg', dir_db, Image_ID_List(i_im).name, Image_ID_List(i_im).name);
        im = imread(str_im);
        imw = im; imw(:) = 255;
        im_gray = rgb2gray(im);
        size_im = size(im_gray);
        
        % load line gnd
        str_gnd = sprintf('%s/%s/%sLinesAndVP.mat', dir_db, Image_ID_List(i_im).name, Image_ID_List(i_im).name);
        
        if ~exist(str_gnd, 'file'), continue; end
        load(str_gnd); % we 'get lines'
        line_gnd = [lines(1:2:end, 1) lines(1:2:end, 2) lines(2:2:end, 1) lines(2:2:end, 2)];
        
        % Rearrange line segment so that elements become (x1, y1, x2, y2, center_x, center_y, length, angle)
        cp = [line_gnd(:,1) + line_gnd(:,3) line_gnd(:,2) + line_gnd(:,4)]/2;
        dx = line_gnd(:,3) - line_gnd(:,1); dy = line_gnd(:,4) - line_gnd(:,2);
        line_gnd = [line_gnd, cp, sqrt(dx.^2 + dy.^2), atan2(dy, dx)];
        
        for k = 1:num_method
            % Load estimation results
            str_est = sprintf('result/%s/%s.mat', dir_method{k}, Image_ID_List(i_im).name);
            if ~exist(str_est, 'file'), continue; end
            load(str_est);
            
            if strcmp(dir_method{k}, 'proposed')
                idxH = ll_Valid_Hor1(:,1) >= ll_Valid_Hor1(:,2);
                idxV = ll_Valid_Ver1(:,1) >= ll_Valid_Ver1(:,2);
                line_own = [ls_est_Hor1(idxH,:); ls_est_Ver1(idxV,:)];
            elseif strcmp(dir_method{k}, 'proposed_h')
                idxH = bValidHelmholtz_Hor == 1;
                idxV = bValidHelmholtz_Ver == 1;
                line_own = [ls_est_Hor_Helm(idxH,:); ls_est_Ver_Helm(idxV,:)];
            end
            
            eval( sprintf('line_est = %s;', var_method{k}) );
            
            if strcmp(dir_method{k}, 'proposed') || strcmp(dir_method{k}, 'proposed_h')
                dir_vec = repmat(line_est(:,4), 1, 2)/2.*[cos(line_est(:,3)) sin(line_est(:,3))];
                x1 = line_est(:,1:2) + dir_vec;
                x2 = line_est(:,1:2) - dir_vec;
                line_est = [x1 x2 line_est(:,1:2) line_est(:,4) line_est(:,3)];
            else
                cp = [line_est(:,1) + line_est(:,3) line_est(:,2) + line_est(:,4)]/2;
                dx = line_est(:,3) - line_est(:,1); dy = line_est(:,4) - line_est(:,2);
                line_est = [line_est(:,1:4), cp, sqrt(dx.^2 + dy.^2), atan2(dy, dx)];
            end
            
            % Evaluate
            if size(line_est,1) > 0
                [pr, re, iou] = evaluate_line_segment(line_est, line_gnd, eval_param);
                
                stats(k).prec(i_im,:) = pr;
                stats(k).rec(i_im,:) = re;
                stats(k).iou(i_im,:) = iou;
            else
                stats(k).prec(i_im,:) = 0;
                stats(k).rec(i_im,:) = 0;
                stats(k).iou(i_im,:) = 0;
            end
        end
        
    catch err
        fprintf('error at i_im: %d.\n', i_im);
        rethrow(err);
    end
end

%% Display scores
AP = zeros(num_method, 1);
AR = zeros(num_method, 1);
IOU = zeros(num_method, 1);

for k = 1:num_method
    AP(k) = mean(stats(k).prec,1);
    AR(k) = mean(stats(k).rec,1);
    IOU(k) = mean(stats(k).iou,1);
end
F_sc = 2 * (AP .* AR) ./ (AP + AR);

fig = figure(1); clf;
axes1 = axes('Parent',fig,'Layer','top','FontWeight','bold','FontSize',12,...
    'FontName','Times New Roman', 'XTick', 1:num_method, 'XTickLabel',name_method);
box(axes1,'on');    hold(axes1,'on');
title('Average precision and recall')
bar([AP, AR, IOU, F_sc])
hleg = legend('AP', 'AR', 'IOU', 'F-score', 'Location', 'nw');
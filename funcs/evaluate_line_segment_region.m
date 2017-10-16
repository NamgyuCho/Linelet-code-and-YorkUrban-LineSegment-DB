function [ prec_score, rec_score, acc_score, f_score, iou_score ] = evaluate_line_segment_region( line_est, line_gnd, size_im )
    %EVALULATE_LINE_SEGMENT_REGION Summary of this function goes here
    %   Detailed explanation goes here
    % line_est = line_own_PA; 
    % --------------------------------------------------------------------------------------------------
    % Initialization
    % --------------------------------------------------------------------------------------------------
    % Convert a line segment to a set of indices -- gnd (x1, y1, x2, y2, center_x, center_y, length, angle)
    num_gnd = size(line_gnd, 1);
    num_est = size(line_est, 1);
    
    %%
    % Draw line images    
    im_gnd = zeros(size_im);
    im_est = zeros(size_im);
    
    % Ground truth image
    for i_gnd = 1:num_gnd
        [xx, yy] = bresenham(line_gnd(i_gnd,1), line_gnd(i_gnd,2), line_gnd(i_gnd,3), line_gnd(i_gnd,4)); 
        idx = xx <= 0 | xx >= size_im(2) | yy <= 0 | yy >= size_im(1);
        xx(idx) = [];   yy(idx) = [];
        
        idx = sub2ind(size_im, uint16(yy), uint16(xx));
        im_gnd(idx) = 1;
    end
    
    % Estimation image
    for i_est = 1:num_est
        [xx, yy] = bresenham(line_est(i_est,1), line_est(i_est,2), line_est(i_est,3), line_est(i_est,4)); 
        idx = xx <= 0 | xx >= size_im(2) | yy <= 0 | yy >= size_im(1);
        xx(idx) = [];   yy(idx) = [];
        idx = sub2ind(size_im, uint16(yy), uint16(xx));
        im_est(idx) = 1;
    end
    
    num_dil_step = 2;
    
    im_g_d = im_gnd;
    im_e_d = im_est;
    se = strel('disk',1);
    iou_score = zeros(num_dil_step, 1);
    prec_score = zeros(num_dil_step, 1);
    rec_score = zeros(num_dil_step, 1);
    f_score = zeros(num_dil_step, 1);
    acc_score = zeros(num_dil_step, 1);
    
    for i=1:num_dil_step
        if i > 1
            %im_g_d = imdilate(im_g_d>0, se);
            im_e_d = imdilate(im_e_d>0, se);
        end
        
        idx = im_g_d > 0 & im_e_d > 0;
%         figure; imagesc(idx);
        
        setT = numel(im_g_d);
        setA = sum(im_g_d(:));
        setB = sum(im_e_d(:));
        
        interAB = sum(idx(:));
        uniAB = setA + setB - interAB;
        
        tp = sum(idx(:));
        fp = setB - tp;
        tn = setT - tp;
        fn = setA - tp;
        
        pr = tp / (tp + fp);
        re = tp / (tp + fn);
        tnr = tn / (tn + fp);
        acc = (tp+tn) / (tp + tn + fp + fn);
        f2 = 2 * (pr * re) / (pr + re);
        
        prec_score(i) = interAB / setB;
        rec_score(i) = interAB / setA;
        acc_score(i) = acc;
        f_score(i) = 2 * (prec_score(i)*rec_score(i)) / (prec_score(i)+rec_score(i));
        iou_score(i) = interAB / uniAB;
    end
    
%     figure; imshowpair(im_g_d, im_e_d, 'montage');    
%     figure; imshowpair(im_gnd>0, (imdilate(im_gnd>0, se)), 'montage');   
end


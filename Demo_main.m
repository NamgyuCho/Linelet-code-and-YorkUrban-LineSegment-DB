%%
clear; close all; clc;

dir_db = './DB/YorkUrbanDB';
addpath(genpath('./toolbox/'));
addpath('./funcs/');
load([dir_db '/our_annotation/Image_ID_List.mat']); % We get "Image_ID_List"
num_im = size(Image_ID_List, 1);

mkdir('./visualization')
mkdir('./result/proposed')
mkdir('./result/proposed_h')

% ---------------------------------------------------------------------------
% Parameter initialization
% ---------------------------------------------------------------------------
% For preprocess 
param.thres_grad_mag = 10;              
param.thres_tau_len = 3;
param.thres_angle_diff = pi/8; % 22.5 degree
param.thres_log_eps = 0.0; % For Helmholtz validation

% For aggregation ----------------------------------------------------
param.est_aggregation = 'Kurtosis'; 
% ---------------------------------------------------------------------------

bFigureVisible = 'off';
bDebugDisplay = false;
bSaveResult = true;
bDisplayComparison = true;
bHelmholtzValidation = true;

%% 
% to run over the whole dataset, remove "22%."
for i_im = 22% 1:num_im 
    % Load image 
    tic
    str_im = sprintf('%s/%s/%s.jpg', dir_db, Image_ID_List(i_im).name, Image_ID_List(i_im).name);
    im = imread(str_im);  
    im_gray = rgb2gray(im);
    size_im = size(im_gray);
    
    % Initialization 
    theta_space = (0:.1:45) * pi/180; % line segment slope range
    BN_Model = generate_BN_LUT(theta_space);
    LineletNum = BN_Model.Numbers;
    level_t = 1;
    
    % Get gradient magnitude and direction
    myfilter = fspecial('gaussian',[3 3], .6);
    im_gray = imfilter(double(im_gray), myfilter, 'replicate');
    I = double(im_gray);
    [im_grad, im_dird] = imgradient( im_gray, 'sobel');
    im_grad = im_grad / 4;
    im_dir = im_dird * pi/180;  
    
    [~, lsc_map, ~, ~, ~] = Get_Line_Candidate(im_gray, im_grad, im_grad, im_dir, im_dird, param, false, false);
    ind = im_grad <= param.thres_grad_mag;
    im_grad(ind) = 0;
    
    % Non-maximal suppression in horizontal and vertical directionss
    [ nmsX, nmsY, memX, memY] = partial_nonmaxsup( im_grad, im_grad, false,1);
    memY(ind) = 0; memX(ind) = 0;
    
    % ------------------------------------------------
    % Linelet detection
    % ------------------------------------------------
    tic
    [ ll_Hor, ll_Ver, ll_Diag, map_Hor, map_Ver, map_Diag ] = labelLL( memX == 1, memY == 1, lsc_map, false );

    % Remove redundant linelets    
    [unq_val, unq_num] = count_unique(lsc_map(:));
    [unq_num, srt_idx] = sort(unq_num, 'descend');
    unq_val = unq_val(srt_idx);
    
    num_lsc = length(unq_val);
    bUsed_H = false(size(ll_Hor,1),1);
    bUsed_V = false(size(ll_Ver,1),1);
    
    idx_dropH = [];
    idx_dropV = [];
    map_H = zeros(size_im);
    map_V = zeros(size_im);
    
    ll_H = [];
    ll_V = [];
    
    for i = 1:max(lsc_map(:))
        if unq_val(i) == 0, continue; end
        idx_lsc = lsc_map == unq_val(i);
        map_H_tmp = double(idx_lsc) .* map_Hor;
        map_V_tmp = double(idx_lsc) .* map_Ver;
        
        % Get activated (corresponding) linelets
        idx_H = find(map_H_tmp);
        idx_tarH = map_H_tmp(idx_H);        
        idx_V = find(map_V_tmp);
        idx_tarV = map_V_tmp(idx_V);
        
        if length(idx_H) >= length(idx_V)
            map_H(idx_H) = 1;
        else
            map_V(idx_V) = 1;
        end
    end
    
    [ ll_Hor, ~, ~, map_Hor, ~, ~ ] = labelLL( map_H, map_H, lsc_map, false );
    [ ~, ll_Ver, ~, ~, map_Ver, ~ ] = labelLL( map_V, map_V, lsc_map, false );    
    fprintf('Elapsed time for linelet detection: %.2f seconds.\n', toc);
    % ------------------------------------------------
    
    % -------------------------------------------------------------------------
    % Calculate unary potential
    % -------------------------------------------------------------------------
    tic
    [upot1_Hor, upot2_Hor, ~, ll_dir_Hor] = UnaryPot( ll_Hor, 'hor', LineletNum, im_gray, im_grad, im_dir, theta_space );
    [upot1_Ver, upot2_Ver, ~, ll_dir_Ver] = UnaryPot( ll_Ver, 'ver', LineletNum, im_gray, im_grad, im_dir, theta_space );
    fprintf('Elapsed time for calculating unary potentials: %.2f seconds.\n', toc);
    % -------------------------------------------------------------------------
    
    % -------------------------------------------------------------------------
    % Grouping linelets
    % -------------------------------------------------------------------------
    tic
    [ ll_idxSet_Hor, bAscend_Hor, im_Hor_as, im_Hor_ds, ll_usage_Hor, ll_Assigned_Hor, ll_Valid_Hor ] = LineletGrouping(ll_Hor, upot1_Hor, 'hor', level_t, im_gray, im_grad, im_dir, false, param );
    [ ll_idxSet_Ver, bAscend_Ver, im_Ver_as, im_Ver_ds, ll_usage_Ver, ll_Assigned_Ver, ll_Valid_Ver ] = LineletGrouping(ll_Ver, upot1_Ver, 'ver', level_t, im_gray, im_grad, im_dir, false, param);
    fprintf('Elapsed time for grouping linelets: %.2f seconds.\n', toc);
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------
    % Estimation and validation 
    % -------------------------------------------------------------------------
    tic
    [pAng_Hor, ls_est_Hor, pFg_Hor, ~] = LineletEstimation( ll_idxSet_Hor, ll_Hor, upot1_Hor, upot2_Hor, bAscend_Hor, 'hor', ll_Valid_Hor, size_im, level_t, theta_space, im_gray, im_grad, im_dir, param, false );
    [pAng_Ver, ls_est_Ver, pFg_Ver, ~] = LineletEstimation( ll_idxSet_Ver, ll_Ver, upot1_Ver, upot2_Ver, bAscend_Ver, 'ver', ll_Valid_Ver, size_im, level_t, theta_space, im_gray, im_grad, im_dir, param, false );
    fprintf('Elapsed time for estimation and validation: %.2f seconds.\n', toc);      
    % -------------------------------------------------------------------------
    
    %% -------------------------------------------------------------------------
    % Aggregation
    % -------------------------------------------------------------------------   
    tic
    level_t = 1;
    ll_idxSet_Hor1 = ll_idxSet_Hor;
    pAng_Hor1 = pAng_Hor;
    ls_est_Hor1 = ls_est_Hor;
    bAscend_Hor1 = bAscend_Hor;
    bMerged_Hor1 = [];
    bLlMerged_Hor1 = false( size(ll_Hor, 1), 1 );
    ll_Assigned_Hor1 = ll_Assigned_Hor;
    ll_Valid_Hor1 = pFg_Hor;
    
    ll_idxSet_Ver1 = ll_idxSet_Ver;
    pAng_Ver1 = pAng_Ver;
    ls_est_Ver1 = ls_est_Ver;
    bAscend_Ver1 = bAscend_Ver;
    bMerged_Ver1 = [];
    bLlMerged_Ver1 = false( size(ll_Ver, 1), 1 );
    ll_Assigned_Ver1 = ll_Assigned_Ver;
    ll_Valid_Ver1 = pFg_Ver;
    
    bDraw = false;
    fprintf('level_t = %d', level_t);
    while 1
        [ll_idxSet_Hor1, pAng_Hor1, ls_est_Hor1, bAscend_Hor1, bMerged_Hor1, bLlMerged_Hor1, ll_Assigned_Hor1, ll_Valid_Hor1, bUpdated_Hor1] =...
            LineletAggregationHi( ll_idxSet_Hor1, pAng_Hor1, ls_est_Hor1, ll_Hor, upot1_Hor, upot2_Hor, bAscend_Hor1, 'hor',...
            level_t, theta_space, im_gray, im_grad, im_dir, bMerged_Hor1, bLlMerged_Hor1, ll_Assigned_Hor1, ll_Valid_Hor1, LineletNum, param, bDraw );
        [ll_idxSet_Ver1, pAng_Ver1, ls_est_Ver1, bAscend_Ver1, bMerged_Ver1, bLlMerged_Ver1, ll_Assigned_Ver1, ll_Valid_Ver1, bUpdated_Ver1] =...
            LineletAggregationHi( ll_idxSet_Ver1, pAng_Ver1, ls_est_Ver1, ll_Ver, upot1_Ver, upot2_Ver, bAscend_Ver1, 'ver',...
            level_t, theta_space, im_gray, im_grad, im_dir, bMerged_Ver1, bLlMerged_Ver1, ll_Assigned_Ver1, ll_Valid_Ver1, LineletNum, param, bDraw );
      
        level_t = level_t + 1;
        fprintf('...%d', level_t);
        
        if (~bUpdated_Hor1 && ~bUpdated_Ver1)
            if bDebugDisplay
                fig1 = figure(1); clf;
                set(fig1, 'Visible', bFigureVisible)
                imshow(uint8(im_gray)); hold on;
                for l = 1:2
                    if     l == 1,  line_own = ls_est_Hor1; lcol = [1 0 0]; bMerged_cur = bMerged_Hor1; ws = 1; pValid = ll_Valid_Hor1;
                    elseif l == 2,  line_own = ls_est_Ver1; lcol = [0 1 0]; bMerged_cur = bMerged_Ver1; ws = 1; pValid = ll_Valid_Ver1;
                    end
                    
                    for k = 1:size(line_own,1)
                        if pValid(k,1) < pValid(k,2)
                            continue;
                        else
                            lcoltar = lcol;
                        end
                        x1 = line_own(k,1:2) + line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
                        x2 = line_own(k,1:2) - line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
                        plot( [x1(1) x2(1)], [x1(2) x2(2)], '-', 'linewidth', ws, 'color', lcoltar );
                    end
                end
            end
            break;
        end
    end
    fprintf('\nElapsed time for aggregation: %.2f seconds.\n\n', toc);      
    
    if bSaveResult
        save(['./result/proposed/' Image_ID_List(i_im).name  '.mat'], 'ls_est_Hor1', 'ls_est_Ver1', 'll_Valid_Hor1', 'll_Valid_Ver1', 'param');
    end
    
    % Validate with the Helmholtz principle -- simply, validate aggregated (final) line segments
    if bHelmholtzValidation
        [~, ls_est_Hor_Helm, ~, bValidHelmholtz_Hor] = LineletEstimation( ll_idxSet_Hor1, ll_Hor, upot1_Hor, upot2_Hor, bAscend_Hor1, 'hor', ll_Valid_Hor1, size_im, level_t, theta_space, im_gray, im_grad, im_dir, param, false );
        [~, ls_est_Ver_Helm, ~, bValidHelmholtz_Ver] = LineletEstimation( ll_idxSet_Ver1, ll_Ver, upot1_Ver, upot2_Ver, bAscend_Ver1, 'ver', ll_Valid_Ver1, size_im, level_t, theta_space, im_gray, im_grad, im_dir, param, false );
        
        if bDebugDisplay
            fig2 = figure(2); clf;
            set(fig2, 'Visible', bFigureVisible)
            imshow(uint8(im_gray)); hold on;
            for l = 1:2
                if     l == 1,  line_own = ls_est_Hor_Helm; lcol = [1 0 0]; ws = 1; pValid = bValidHelmholtz_Hor;
                elseif l == 2,  line_own = ls_est_Ver_Helm; lcol = [0 1 0]; ws = 1; pValid = bValidHelmholtz_Ver;
                end
                
                for k = 1:size(line_own,1)
                    if ~pValid(k)
                        continue;
                    else
                        lcoltar = lcol;
                    end
                    x1 = line_own(k,1:2) + line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
                    x2 = line_own(k,1:2) - line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
                    plot( [x1(1) x2(1)], [x1(2) x2(2)], '-', 'linewidth', ws, 'color', lcoltar );
                end                
            end
        end
        
        if bSaveResult
            str_save_name = sprintf('%s', Image_ID_List(i_im).name);
            save(['./result/proposed_h/' str_save_name '.mat'], 'ls_est_Hor_Helm', 'ls_est_Ver_Helm', 'bValidHelmholtz_Hor', 'bValidHelmholtz_Ver', 'param');
        end
    end    
    % -------------------------------------------------------------------------    
    
    %% Performance comparison and additional illustration 
    if bDisplayComparison        
        fig = figure(3);clf;
        screensize = get( 0, 'Screensize' );
        set(fig, 'Position', [100, 100, 1000, 900]); 
        
        % Draw proposed
        imshow(ones(size(im_gray))*255); hold on;
        line_own = [ls_est_Hor1; ls_est_Ver1]; 
        pValid = [ll_Valid_Hor1; ll_Valid_Ver1];            
        for k = 1:size(line_own,1)
            if pValid(k,1) < pValid(k,2)
                continue;
            end
            x1 = line_own(k,1:2) + line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
            x2 = line_own(k,1:2) - line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'k-' );
        end
    end
end
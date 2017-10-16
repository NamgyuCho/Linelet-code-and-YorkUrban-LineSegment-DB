function [ p_of_agg ] = AggregationProbability( ls_estSet, ll_instSet, im_gray, im_grad, im_dird, ll_type )
    %AGGREGATIONPROBABILITY Summary of this function goes here
    %   Detailed explanation goes here
    size_im = size(im_gray);
    
    sig_fg_int = .2;    sig_bg_int = .5;
    sig_fg_gmag = .2;   sig_bg_gmag = .5;
    sig_fg_gdir = .2;   sig_bg_gdir = .5;
    
    
%     prior_fg = 0.5;     prior_bg = 0.5;
    % Geometric prior wrt angle
    tau_ang = pi/16;
    diff_ang = GetAngleDiff(ls_estSet(1,3), ls_estSet(2,3));
    p_ang = normpdf(diff_ang, 0, tau_ang*3);
    
%     x = -pi:.05:pi;      y = sigmf(x, [1 .5]);      figure; plot(x, y, 'r-')
    
    % Geometric prior wrt central distance
    dcentLs = ls_estSet(1,1:2) - ls_estSet(2,1:2);
    dcentLs = sqrt( sum( dcentLs.^2 ) );
    dbtLs = abs(dcentLs - ls_estSet(1,4)/2 - ls_estSet(2,4)/2);
    p_dist = normpdf(dbtLs, 0, .5);
    
    prior_fg = .5 + p_ang + p_dist;
    prior_bg = .5 + 2 - p_ang - p_dist;
    tmp = prior_fg + prior_bg;
    prior_fg = prior_fg / tmp;
    prior_bg = prior_bg / tmp;
    
    
    cue_weight = [.3 .3 .4]; % intensity, grad mag, grad dir1, grad dir2
    
    
    nbsize1 = 256;
    nbsize2 = 180;
    
    res1 = 8;
    res2 = 6;
    
    nbin1 = nbsize1 / res1;
    nbin2 = nbsize2 / res2;
    
    h_inst_ome = zeros(2, nbin1);   h_inst_gamU = zeros(2, nbin1);   h_inst_gamL = zeros(2, nbin1);
    h_gmag_ome = zeros(2, nbin1);   h_gmag_gamU = zeros(2, nbin1);   h_gmag_gamL = zeros(2, nbin1);
    h_gdir_ome = zeros(2, nbin2);   h_gdir_gamU = zeros(2, nbin2);   h_gdir_gamL = zeros(2, nbin2);
    
    for k = 1:size(ll_instSet, 1)
        ll_inst = ll_instSet{k};
        [~, idx_srt] = sort(ll_inst(:,1));
        ll_inst = ll_inst(idx_srt,:);
    
        for i = 1:size(ll_inst,1)
            [ome, gamU, gamL] = LineletOmegaGamma(ll_inst(i,:), ll_type, size_im);
            
            % intensity
            tmp = uint16(im_gray(ome) / res1);
            tmp( tmp == nbin1 ) = nbin1 - 1;
            h_inst_ome(k, tmp+1) = h_inst_ome(k, tmp+1) + 1;
            
            tmp = uint16(im_gray(gamU) / res1);
            tmp( tmp == nbin1 ) = nbin1 - 1;
            h_inst_gamU(k, tmp+1) = h_inst_gamU(k, tmp+1) + 1;
            
            tmp = uint16(im_gray(gamL) / res1);
            tmp( tmp == nbin1 ) = nbin1 - 1;
            h_inst_gamL(k, tmp+1) = h_inst_gamL(k, tmp+1) + 1;
            
            % gradient magnitude            
            tmp = uint16(im_grad(ome) / res1);
            tmp( tmp == nbin1 ) = nbin1 - 1;
            h_gmag_ome(k, tmp+1) = h_gmag_ome(k, tmp+1) + 1;
            
            tmp = uint16(im_grad(gamU) / res1);
            tmp( tmp == nbin1 ) = nbin1 - 1;
            h_gmag_gamU(k, tmp+1) = h_gmag_gamU(k, tmp+1) + 1;
            
            tmp = uint16(im_grad(gamL) / res1);
            tmp( tmp == nbin1 ) = nbin1 - 1;
            h_gmag_gamL(k, tmp+1) = h_gmag_gamL(k, tmp+1) + 1;
            
            % gradient direction
            tmp = int16(im_dird(ome));
            tmp(tmp <= 0) = tmp(tmp <= 0) + 180;
            tmp(tmp == 0) = tmp(tmp == 0) + 180;
            tmp = uint16(tmp / res2);
            tmp( tmp == nbin2 ) = nbin2 - 1;
            h_gdir_ome(k, tmp+1) = h_gdir_ome(k, tmp+1) + 1;
            
            tmp = int16(im_dird(gamU));
            tmp(tmp <= 0) = tmp(tmp <= 0) + 180;
            tmp(tmp == 0) = tmp(tmp == 0) + 180;
            tmp = uint16(tmp / res2);
            tmp( tmp == nbin2 ) = nbin2 - 1;
            h_gdir_gamU(k, tmp+1) = h_gdir_gamU(k, tmp+1) + 1;
            
            tmp = int16(im_dird(gamL));
            tmp(tmp <= 0) = tmp(tmp <= 0) + 180;
            tmp(tmp == 0) = tmp(tmp == 0) + 180;
            tmp = uint16(tmp / res2);
            tmp( tmp == nbin2 ) = nbin2 - 1;
            h_gdir_gamL(k, tmp+1) = h_gdir_gamL(k, tmp+1) + 1;
        end
    end
    % intensity            
    tmp = sum(h_inst_ome, 2);
    h_inst_ome = [h_inst_ome(1,:)/tmp(1); h_inst_ome(2,:)/tmp(2)];
    dVal = DistMeasure(h_inst_ome(1,:), h_inst_ome(2,:), 'chi2');
    likeH_int = [normpdf(dVal, 0, sig_fg_int), normpdf(dVal, 0, sig_bg_int)]; % [positive, negative]
    
    tmp = sum(h_inst_gamU, 2);
    h_inst_gamU = [h_inst_gamU(1,:)/tmp(1); h_inst_gamU(2,:)/tmp(2)];
    dVal = DistMeasure(h_inst_gamU(1,:), h_inst_gamU(2,:), 'chi2');
    likeH_int = [likeH_int; normpdf(dVal, 0, sig_fg_int), normpdf(dVal, 0, sig_bg_int)]; % [positive, negative]
    
    tmp = sum(h_inst_gamL, 2);
    h_inst_gamL = [h_inst_gamL(1,:)/tmp(1); h_inst_gamL(2,:)/tmp(2)];
    dVal = DistMeasure(h_inst_gamL(1,:), h_inst_gamL(2,:), 'chi2');
    likeH_int = [likeH_int; normpdf(dVal, 0, sig_fg_int), normpdf(dVal, 0, sig_bg_int)]; % [positive, negative]
    
    % gradient magnitude
    tmp = sum(h_gmag_ome, 2);
    h_gmag_ome = [h_gmag_ome(1,:)/tmp(1); h_gmag_ome(2,:)/tmp(2)];
    dVal = DistMeasure(h_gmag_ome(1,:), h_gmag_ome(2,:), 'chi2');
    likeH_gmag = [normpdf(dVal, 0, sig_fg_gmag), normpdf(dVal, 0, sig_bg_gmag)];
    
    tmp = sum(h_gmag_gamU, 2);
    h_gmag_gamU = [h_gmag_gamU(1,:)/tmp(1); h_gmag_gamU(2,:)/tmp(2)];
    dVal = DistMeasure(h_gmag_gamU(1,:), h_gmag_gamU(2,:), 'chi2');
    likeH_gmag = [likeH_gmag; normpdf(dVal, 0, sig_fg_gmag), normpdf(dVal, 0, sig_bg_gmag)];
    
    tmp = sum(h_gmag_gamL, 2);
    h_gmag_gamL = [h_gmag_gamL(1,:)/tmp(1); h_gmag_gamL(2,:)/tmp(2)];
    dVal = DistMeasure(h_gmag_gamL(1,:), h_gmag_gamL(2,:), 'chi2');
    likeH_gmag = [likeH_gmag; normpdf(dVal, 0, sig_fg_gmag), normpdf(dVal, 0, sig_bg_gmag)];
    
    % gradient direction
    tmp = sum(h_gdir_ome, 2);
    h_gdir_ome = [h_gdir_ome(1,:)/tmp(1); h_gdir_ome(2,:)/tmp(2)];
    dVal = DistMeasure(h_gdir_ome(1,:), h_gdir_ome(2,:), 'chi2');
    likeH_gdir = [normpdf(dVal, 0, sig_fg_gdir), normpdf(dVal, 0, sig_bg_gdir)];
    
    tmp = sum(h_gdir_gamU, 2);
    h_gdir_gamU = [h_gdir_gamU(1,:)/tmp(1); h_gdir_gamU(2,:)/tmp(2)];
    dVal = DistMeasure(h_gdir_gamU(1,:), h_gdir_gamU(2,:), 'chi2');
    likeH_gdir = [likeH_gdir; normpdf(dVal, 0, sig_fg_gdir), normpdf(dVal, 0, sig_bg_gdir)];
    
    tmp = sum(h_gdir_gamL, 2);
    h_gdir_gamL = [h_gdir_gamL(1,:)/tmp(1); h_gdir_gamL(2,:)/tmp(2)];
    dVal = DistMeasure(h_gdir_gamL(1,:), h_gdir_gamL(2,:), 'chi2');
    likeH_gdir = [likeH_gdir; normpdf(dVal, 0, sig_fg_gdir), normpdf(dVal, 0, sig_bg_gdir)];
    
    p_of_agg = [];
    for i =1:size(likeH_gdir,1)
        likelihood = [likeH_int(i,:); likeH_gmag(i,:); likeH_gdir(i,:)];
        term1_denom = likelihood(:,1)*prior_fg + likelihood(:,2)*prior_bg;
        term1 = likelihood(:,1)*prior_fg ./ term1_denom;
        p_of_agg(i,1) = cue_weight * term1;
    end
    
    p_of_agg = [sum(p_of_agg)/3 1-sum(p_of_agg)/3];
    
%     likelihood = [likeH_int; likeH_gmag; likeH_gdir];
%     term1_denom = likelihood(:,1)*prior_fg + likelihood(:,2)*prior_bg;
%     term1 = likelihood(:,1)*prior_fg ./ term1_denom;
%     term2 = likelihood(:,2)*prior_bg ./ term1_denom;
%     
%     p_of_agg = [p_of_agg; cue_weight * term1, cue_weight * term2];
end


function [ p_of_fg ] = FgBgProbability( ls_inst, ll_inst, im_gray, im_grad, im_dir, ll_type, param )
    %FGBGPROBABILITY Summary of this function goes here
    %   Detailed explanation goes here
    
    size_im = size(im_gray);
    
    [~, idx_srt] = sort(ll_inst(:,1));
    ll_inst = ll_inst(idx_srt,:);
    
    denom_eps = .01;
    
    if strcmp(ll_type, 'hor')
        base_factor = abs(tan(ls_inst(3))) + denom_eps;
    else
        if ls_inst(3) >= 0
            base_factor = abs(tan(ls_inst(3)-pi/2)) + denom_eps;
        else
            base_factor = abs(tan(ls_inst(3)+pi/2)) + denom_eps;
        end
    end
    
    mu_fg_gmag = base_factor;
    sig_fg_gmag = 0.3675;
    
    mu_bg_gmag = 0.9423;
    sig_bg_gmag = 0.3675;
    
    mu_fg_gmag2 = 0;
    sig_fg_gmag2 = 1;%mu_fg_gmag*2;%base_factor +.1;
    
    mu_bg_gmag2 = 0.25;
    sig_bg_gmag2 = 1;
    
    mu_fg_gdir = base_factor;
    sig_fg_gdir = 1;%mu_fg_gdir*2;%base_factor +.1; 
    
    mu_bg_gdir = 1.0854;
    sig_bg_gdir = mu_bg_gdir/3;
    
    mu_fg_gdir2 = base_factor;
    sig_fg_gdir2 = 1;%mu_fg_gdir*2;%base_factor +.1; 
    
    mu_bg_gdir2 = 1.0854;
    sig_bg_gdir2 = mu_bg_gdir2/3;
     
    total_l_len = sum(ll_inst(:,5));
    
    prior_fg = sigmf(total_l_len, [.2 20]); % prev [.2 24.63]
    prior_bg = 1-prior_fg;
    cue_weight = [.2 .2 .2 .2 .2]; % grad mag, grad dir
      
    dist_int = []; dist_gmag = []; dist_gdir = [];
    dist_gmag2 = [];
    dist_gdir2 = [];
    likeH_int = [];
    likeH_gmag = [];
    likeH_gdir = [];
    for i = 1:size(ll_inst,1)
        [ome, gamU, gamL] = LineletOmegaGamma(ll_inst(i,:), ll_type, size_im);
        
        tmp = double(im_gray(ome));
        dist_int = [dist_int, tmp / max(tmp)];
        
        tmp = [im_grad(gamU); im_grad(ome); im_grad(gamL)];
        tmp = sum(tmp,1);
        dist_gmag = [dist_gmag, tmp / max(tmp)];
        
        tmp = im_grad(ome);
        dist_gmag2 = [dist_gmag2, tmp / max(tmp)];
        
        tmp = im_dir(ome);
        dist_gdir = [dist_gdir, tmp/ max(tmp)];
        dist_gdir2 = [dist_gdir2, tmp];
    end
    dfv = diff(double(dist_int), 1, 2);
    valx = mean(abs(dfv));
    likeH_int = [normpdf(valx, mu_fg_gmag, sig_fg_gmag), normpdf(valx, mu_bg_gmag, sig_bg_gmag)];
    
    dfv = diff(double(dist_gmag), 1, 2);
    valx = mean(abs(dfv));
    likeH_gmag = [normpdf(valx, mu_fg_gmag, sig_fg_gmag), normpdf(valx, mu_bg_gmag, sig_bg_gmag)];
    
    valx = mean(sum(dist_gmag2,1));
    likeH_gmag2 = [normpdf(valx, mu_fg_gmag2, sig_fg_gmag2), normpdf(valx, mu_bg_gmag2, sig_bg_gmag2)];
    
    dfv = GetAngleDiff(dist_gdir(1:end-1), dist_gdir(2:end));
    valx = mean(abs(dfv));
    likeH_gdir = [normpdf(valx, mu_fg_gdir, sig_fg_gdir), normpdf(valx, mu_bg_gdir, sig_bg_gdir)];

    tmp = bAngleAligned(dist_gdir2(1:end-1), dist_gdir2(2:end), pi/8);
    tmp = sum(tmp) / length(tmp);
    likeH_gdir2 = [tmp 1-tmp];
    
    likelihood = [likeH_int; likeH_gmag; likeH_gmag2; likeH_gdir; likeH_gdir2];
    term1_denom = likelihood(:,1)*prior_fg + likelihood(:,2)*prior_bg;
    term1 = likelihood(:,1)*prior_fg ./ term1_denom;
    p_of_f = cue_weight * term1;
    
    p_of_fg = [p_of_f 1-p_of_f];
end
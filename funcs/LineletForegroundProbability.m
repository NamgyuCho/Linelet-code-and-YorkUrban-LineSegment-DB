function [ p_of_fg ] = LineletForegroundProbability( ll_inst, im_gray, im_grad, im_dir, ll_type )
    %LINELETFOREGROUNDPROBABILITY Summary of this function goes here
    %   Detailed explanation goes here
    size_im = size(im_gray);
    
    [~, idx_srt] = sort(ll_inst(:,1));
    ll_inst = ll_inst(idx_srt,:);
    
    sig_fg_int = 8;    sig_bg_int = 32;
    sig_fg_gmag = 8;   sig_bg_gmag = 32;
    sig_fg_gdir = pi/32; sig_bg_gdir = pi/2;
    
    
% sig_fg_gmag = 8;   sig_bg_gmag = 32;
% sig_fg_gdir = pi/32; sig_bg_gdir = pi/2;
% sig_fg_int = 8;    sig_bg_int = 32;
% prior_fg = 0.5;     prior_bg = 0.5;
    
    total_l_len = sum(ll_inst(:,5));
    
%     total_l_len = -10:.5:20;
%     prior_fg = 1 ./ (1 + exp(-.3 * total_l_len - 15));
%     figure; plot(total_l_len, prior_fg);
    
    prior_fg = sigmf(total_l_len, [.1 10]);
    prior_bg = 1-prior_fg;
%     x = 0:.5:50; y = sigmf(x, [.2 10]); figure; plot(x, y);
    
    
    cue_weight = [.2 .2 .2 .4]; % intensity, grad mag, grad dir1, grad dir2
    
    dist_int = []; dist_gmag = []; dist_gdir = [];
    dist_gdir2 = [];
    for i = 1:size(ll_inst,1)
        [ome, gamU, gamL] = LineletOmegaGamma(ll_inst(i,:), ll_type, size_im);
        
%         tmp = [im_gray(gamU); im_gray(ome); im_gray(gamL)];
        tmp = im_gray(ome);
        dist_int = [dist_int, tmp];
        
        tmp = [im_grad(gamU); im_grad(ome); im_grad(gamL)];
        tmp = im_grad(ome);
        dist_gmag = [dist_gmag, tmp];
        
%         tmp = [im_dir(gamU); im_dir(ome); im_dir(gamL)];
        tmp = im_dir(ome);
        dist_gdir = [dist_gdir, tmp];
        dist_gdir2 = [dist_gdir2, tmp];
    end
    size_obs = size(dist_int, 1);
    
    dfv = diff(double(dist_int), 1, 2);
    for i = 1:size_obs
        likeH_int(i,:) = [normpdf(std(dfv(i,:)), 0, sig_fg_int), normpdf(std(dfv(i,:)), 0, sig_bg_int)];
    end
    
    dist_gmag = sum(dist_gmag, 1);
    dfv = diff(double(dist_gmag), 1, 2);
    for i = 1:size_obs
        likeH_gmag(i,:) = [normpdf(std(dfv(i,:)), 0, sig_fg_gmag), normpdf(std(dfv(i,:)), 0, sig_bg_gmag)];
    end
    
    dfv = GetAngleDiff(dist_gdir(i, 1:end-1), dist_gdir(i, 2:end));
%     dfv = diff(double(dist_gdir), 1, 2);
    for i = 1:size_obs
        likeH_gdir(i,:) = [normpdf(std(dfv(i,:)), 0, sig_fg_gdir), normpdf(std(dfv(i,:)), 0, sig_bg_gdir)];
    end
    
    for i = 1:size_obs
        tmp = bAngleAligned(dist_gdir2(i, 1:end-1), dist_gdir2(i, 2:end), pi/8);
        tmp = sum(tmp) / length(tmp);
        likeH_gdir2(i,:) = [tmp 1-tmp];
    end
    
    p_of_f = [];
    for i =1:size_obs
        likelihood = [likeH_int(i,:); likeH_gmag(i,:); likeH_gdir(i,:); likeH_gdir2(i,:)];
        term1_denom = likelihood(:,1)*prior_fg + likelihood(:,2)*prior_bg;
        term1 = likelihood(:,1)*prior_fg ./ term1_denom;
        p_of_f(i,1) = cue_weight * term1;
    end
    
    p_of_fg = [sum(p_of_f)/size_obs 1-sum(p_of_f)/size_obs];
 %%   
%     p_of_fg = 0;
%     dist_int = [];
%     dist_gmag = [];
%     dist_gdir = [];
%     dist_gdir2 = [];
%     
%     for i = 1:size(ll_inst,1)
%         [ome, gamU, gamL] = LineletOmegaGamma(ll_inst(i,:), ll_type, size_im);
%         
%         tmp = [im_gray(gamU); im_gray(ome); im_gray(gamL)];
%         dist_int = [dist_int, tmp];
%         
%         tmp = [im_grad(gamU); im_grad(ome); im_grad(gamL)];
%         dist_gmag = [dist_gmag, tmp];
%         
%         tmp = [im_dir(gamU); im_dir(ome); im_dir(gamL)];
%         dist_gdir = [dist_gdir, tmp];
%         dist_gdir2 = [dist_gdir2, im_dir(ome)];
%     end
%     %     prob_1 = [normpdf(std(sum(dist_gmag, 1)), mean(sum(dist_gmag,1)), 1) normpdf(std(sum(dist_gmag, 1)), mean(sum(dist_gmag,1)), 128)];
%     dfv = double(diff(dist_int(2,:)));
%     %     likeH_int = [likeH_int; normpdf(std(dfv), 0, sig_fg_int), normpdf(std(dfv), 0, sig_bg_int)];
%     likeH_int = [normpdf(std(dfv), 0, sig_fg_int), normpdf(std(dfv), 0, sig_bg_int)]; % [positive, negative]
%     
%     dfv = diff(sum(dist_gmag,1));
%     %     likeH_gmag = [likeH_gmag; normpdf(std(dfv), 0, sig_fg_gmag), normpdf(std(dfv), 0, sig_bg_gmag)];
%     likeH_gmag = [normpdf(std(dfv), 0, sig_fg_gmag), normpdf(std(dfv), 0, sig_bg_gmag)];
%     
%     dfv = diff(dist_gdir(2,:));
%     %     likeH_gdir = [likeH_gdir; normpdf(std(dfv), 0, sig_fg_gdir), normpdf(std(dfv), 0, sig_bg_gdir)];
%     likeH_gdir = [normpdf(std(dfv), 0, sig_fg_gdir), normpdf(std(dfv), 0, sig_bg_gdir)];
%     
%     tmp = bAngleAligned(dist_gdir2(1:end-1), dist_gdir2(2:end), pi/8);
%     tmp = sum(tmp) / length(tmp);
%     %     likeH_gdir2 = [likeH_gdir2; tmp];
%     likeH_gdir2 = [tmp 1-tmp];
%     
%     %     sig_bg = [sig_bg; std(double(dist_int(:,2))) std(sum(dist_gmag,1)) std(dist_gdir(:,2))];
%     %     likeH_prob = [likeH_prob; max(pAng_Hor(idx_fp(k),:))];
%     % calculate p of f- given observation O
%     
%     likelihood = [likeH_int; likeH_gmag; likeH_gdir; likeH_gdir2];
%     term1_denom = likelihood(:,1)*prior_fg + likelihood(:,2)*prior_bg;
%     term1 = likelihood(:,1)*prior_fg ./ term1_denom;
%     term2 = likelihood(:,2)*prior_bg ./ term1_denom;
    
%     p_of_fg = [cue_weight * term1, cue_weight * term2];    
end
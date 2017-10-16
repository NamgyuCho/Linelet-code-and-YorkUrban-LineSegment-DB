function [ upot1, upot2, LSM_cand, ll_dir ] = UnaryPot( ll_inst, ll_type, LineletNum, im_gray, im_grad, im_dir, theta_space )
    bUseSloptChar = false;
    
    num_inst = size(ll_inst, 1);
    size_im = size(im_grad);
    size_theta = length(theta_space);
        
    denom_eps = 1; % epsilon value of division to avoid inifinity result
    
    upot1 = zeros(num_inst, 3);
    upot2 = zeros(num_inst, size_theta);    
    
    LSM_cand = cell(num_inst, 1);
    
    ll_dir = zeros(num_inst, 7);
    
    alpha = 1;
    sig = 1;
    
    for i = 1:num_inst
        try            
            % Get gamma and omega pixel sets corresponding to ith linelet
            [ome, gamU, gamL] = LineletOmegaGamma(ll_inst(i,:), ll_type, size_im);
            [yy, xx] = ind2sub(size_im, ome);
            
            
            % ----------------------------------------
            % Intrinsici properties
            dist_intrinsic = [im_gray(ome); im_grad(ome); im_dir(ome)];
            dist_grad = im_grad(ome);
            dist_int = im_gray(ome);
            dist_dir = im_dir(ome);
            % ----------------------------------------
            
            
            % ----------------------------------------
            % 1. Gradient magnitude
            % 1-1. position
            mu_x = [ll_inst(i,1) + ll_inst(i,3), ll_inst(i,2) + ll_inst(i,4)]/2;
            wx = [xx; yy] * im_grad(ome)' / sum(im_grad(ome));
            tmp = sum( abs(mu_x - wx') );
            upot_mag1 = exp( -alpha * tmp / sig );
            
            % 1-2. magnitude consistency
            if isempty(gamU) % actually, not certain 100%, exceptional cases don't happen
                nom1 = std(im_grad(ome) + im_grad(gamL));
                fprintf('exception of gamma_U at %d.\n', i);
                ll_inst(i,:)
            elseif isempty(gamL)
                nom1 = std(im_grad(ome) + im_grad(gamU));
                fprintf('exception of gamma_L at %d.\n', i);
                ll_inst(i,:)
            else
                nom1 = std(im_grad(ome) + im_grad(gamU) + im_grad(gamL));
            end
%             upot_mag2 = exp( -alpha * nom1 / sig );
            val_ideal = tan(theta_space);
            
            if ll_inst(i,5) == 1
                dv = 1;
            else
                dv = mean(abs(diff(dist_grad/max(dist_grad))));
            end
            
            upot_mag2 = (dv - val_ideal).^2 + denom_eps;
            upot_mag2 = 1./upot_mag2;
            % ----------------------------------------
            
            % ----------------------------------------
            % 2. Gradient angle
            % 2-1. direction consistency
            if isempty(gamU)
                nom1 = std(im_dir(ome)) + std(im_dir(gamL));
                nom1 = nom1/2;
            elseif isempty(gamL)
                nom1 = std(im_dir(gamU)) + std(im_dir(ome));
                nom1 = nom1/2;
            else
                nom1 = std(im_dir(gamU)) + std(im_dir(ome)) + std(im_dir(gamL));
                nom1 = nom1/3;
            end
%             upot_dir1 = exp( -alpha * nom1 / sig );

            val_ideal = tan(theta_space);
            
            if ll_inst(i,5) == 1
                dv = 1;
            else
                dv = mean(abs(diff(dist_dir/max(dist_grad))));
            end
            
            upot_dir1 = (dv - val_ideal).^2 + denom_eps;
            upot_dir1 = 1./upot_dir1; %figure; bar(upot_dir1);            
            % ----------------------------------------
            
            % ----------------------------------------
            % 3. Slope character: 1) flat, 2) ascending, 3) descending, 4) random noise
            if bUseSloptChar
            true_flat = ones(1, length(ome)) / length(ome);
            true_asce = 1:length(ome); true_asce = true_asce/sum(true_asce);
            true_desc = length(ome):-1:1; true_desc = true_desc/sum(true_desc);
            
            dist_ome  = im_grad(ome);       dist_ome  = dist_ome/sum(dist_ome);
            dist_gamU = im_grad(gamU);      dist_gamU = dist_gamU/sum(dist_gamU);
            dist_gamL = im_grad(gamL);      dist_gamL = dist_gamL/sum(dist_gamL);
            
            [~,iMin1] = min( [std(im_dir(ome)) abs(std(im_dir(ome)) - stdDir_rand)] );
            [~,iMin2] = min( [std(im_grad(ome)) abs(std(im_grad(ome)) - stdMag_rand)] );
            [~,iMin3] = min( [std(single(im_gray(ome))) abs(std(single(im_gray(ome))) - stdMag_rand)] );
            
            ll_dir(i,:) = [KL_Div( true_flat, dist_gamU ) KL_Div( true_flat, dist_gamL )...
                           KL_Div( true_asce, dist_gamU ), KL_Div( true_asce, dist_gamL )...
                           KL_Div( true_desc, dist_gamU ), KL_Div( true_desc, dist_gamL )...
                           iMin1 == 1 || iMin2 == 1 || iMin3 == 1]; % or iMin1==1 && iMin2==1];
            end
            % ----------------------------------------
                       
            % ----------------------------------------
            % 4. Difference wrt Linelet Number (model)            
            % Following procedure turns out that it causes errors, thus stop using this
            upot_LN = (max(LineletNum, [], 2) - ll_inst(i,5)).^2 + denom_eps;
            upot_LN = 1./upot_LN; %figure; bar(upot_LN);
            
            % Final potential
            upot2(i,:) = upot_mag2 + upot_dir1 + upot_LN';
            
            LSM_cand{i} = find( (min(LineletNum,[], 2) - ll_inst(i,5) >= 0) ~= 0 );
            % ----------------------------------------
        catch err
            fprintf('UnaryPot(): Error occured at %d.\n', i);
            ll_inst(i,:)
            rethrow(err);
        end
    end
end

function [ome, gamU, gamL] = LL_GammaOmega(inst, ll_type, size_im)
    %--- Variable set for testing
    %  inst = ll_inst(i,:); 
    ome = []; gamU = []; gamL = []; 
    
    if strcmp(ll_type, 'hor')
        ome = sub2ind(size_im, repmat(inst(2), [1 inst(5)]), inst(1):inst(3));
        if inst(2) == 1
            gamL = sub2ind(size_im, repmat(inst(2)+1, [1 inst(5)]), inst(1):inst(3));
        elseif inst(2) == size_im(1)
            gamU = sub2ind(size_im, repmat(inst(2)-1, [1 inst(5)]), inst(1):inst(3));
        else
            gamL = sub2ind(size_im, repmat(inst(2)-1, [1 inst(5)]), inst(1):inst(3));
            gamU = sub2ind(size_im, repmat(inst(2)+1, [1 inst(5)]), inst(1):inst(3));
        end
    elseif strcmp(ll_type, 'ver')
        ome = sub2ind(size_im, inst(2):inst(4), repmat(inst(1), [1 inst(5)]));
        if inst(1) == 1
            gamL = sub2ind(size_im, inst(2):inst(4), repmat(inst(1)+1, [1 inst(5)]));
        elseif inst(1) == size_im(2)
            gamU = sub2ind(size_im, inst(2):inst(4), repmat(inst(1)-1, [1 inst(5)]));
        else
            gamL = sub2ind(size_im, inst(2):inst(4), repmat(inst(1)-1, [1 inst(5)]));
            gamU = sub2ind(size_im, inst(2):inst(4), repmat(inst(1)+1, [1 inst(5)]));
        end
    elseif strcmp(ll_type, 'diag')
        xint = 1; yint = 1;
        if inst(1) > inst(3), xint = -1; end
        if inst(2) > inst(4), yint = -1; end
        
        ome = sub2ind(size_im, inst(2):yint:inst(4), inst(1):xint:inst(3));
        if inst(1) == 1 
            gamL = sub2ind(size_im, inst(2):yint:inst(4), inst(1)+1:xint:inst(3)+1);
        elseif inst(1) == size_im(2)
            gamU = sub2ind(size_im, inst(2):yint:inst(4), inst(1)-1:xint:inst(3)-1);
        else
            gamL = sub2ind(size_im, inst(2):yint:inst(4), inst(1)+1:xint:inst(3)+1);
            gamU = sub2ind(size_im, inst(2):yint:inst(4), inst(1)-1:xint:inst(3)-1);
        end
    end
end
function [region_candi_map, im_lc_map, lc_st] = Get_Line_Candidate_Validate_Param(im_grad, im_dir, param, bDrawRet)
%GET_LINE_CANDIDATE_VALIDATE_PARAM Summary of this function goes here
%   Detailed explanation goes here


function [region_candi_map, im_lc_map, lc_st] = Get_Line_Candidate(im_grad, im_dir, param, bDrawRet)
size_im = size(im_grad);



num_sample = [100 500 1000];
std_sampling_length = [1 3 5];
std_sampling_theta = [3 5 10 20]/180*pi;
alpha = [10 100 1000];
beta = [10 100 1000];
gamma = [10 100 1000];


if bDrawRet
    %     figure(111); imagesc(im_grad); colormap gray; hold on;
end

for i1 = 1:length(num_sample)
    for i2 = 1:length(std_sampling_length)
        for i3 = 1:length(std_sampling_theta)
            for i4 = 1:length(alpha)
                for i5 = 1:length(beta)
                    for i6 = 1:length(gamma)
                        
                        figure(111); clf; imagesc(im_grad); colormap gray; hold on;
                        [pixel_grad_val, ind] = sort(im_grad(:), 'descend');
                        [ind_r, ind_c] = ind2sub(size_im, ind);
                        angle_val = im_dir(ind);
                        im_lc_map = zeros(size_im);
                        
                        tmp_ind = find(pixel_grad_val > 2);
                        ind = ind(tmp_ind);
                        pixel_grad_val = pixel_grad_val(tmp_ind);
                        ind_r = ind_r(tmp_ind);
                        ind_c = ind_c(tmp_ind);
                        angle_val = angle_val(tmp_ind);
                        
                        lc_st = [];
                        lc_state = [];
                        pixel_used = zeros(size_im);
                        usage_count = 1;
                        
                        for i = 1:length(ind_r)
                            if pixel_used(ind_r(i), ind_c(i))
                                continue;
                            end
                            
                            ind_lc = 1;
                            lc = [ind_c(i), ind_r(i)];
                            weight = im_grad(ind_r(i), ind_c(i));
                            reg_angle = angle_val(i);
                            reg_angle1 = angle_val(i);
                            sumdx = cos( reg_angle );
                            sumdy = sin( reg_angle );
                            
                            ii = 1;
                            
                            while ii <= ind_lc
                                for xx = lc(ii, 1)-1:lc(ii, 1)+1
                                    for yy = lc(ii, 2)-1:lc(ii, 2)+1
                                        if xx <= size_im(2) && xx >= 1 && yy <= size_im(1) && yy >= 1 &&...
                                                pixel_used(yy,xx) == 0 && angle_alignment(reg_angle, im_dir(yy,xx), param.merge_angle_diff)
                                            
                                            pixel_used(yy, xx) = usage_count;
                                            lc = [lc; xx, yy];
                                            im_lc_map(yy, xx) = usage_count;
                                            weight = [weight; im_grad(yy, xx)];
                                            ind_lc = ind_lc + 1;
                                            sumdx = sumdx + cos(im_dir(yy, xx));
                                            sumdy = sumdy + sin(im_dir(yy, xx));
                                            reg_angle = atan2(sumdy, sumdx);
                                            reg_angle1 = reg_angle1 * .7 + im_dir(yy, xx) * .3;
                                        end
                                    end
                                end
                                ii = ii + 1;
                            end
                            if size(lc,1) < 10
                                continue;
                            end
                            
                            lc_st(usage_count,1).pos = lc;
                            
                            % --------------------------------------------------------
                            % Parameterize each candidate
                            center = [lc(:,1)' * weight, lc(:,2)' * weight] / sum(weight);
                            [coef, scr] = princomp(lc);
                            lengths = max(abs(scr));
                            
                            % Get main pc's angle
                            pc_angle = atan2(-coef(1,2), coef(1,1));
                            if ~angle_alignment(pc_angle, reg_angle, 45/180*pi)
                                continue;
                            end
                            
                            % Fit a line to current set of pixels
                            
                            
                            
                            x_sample = (rand(num_sample(i1), 1) * lengths(2) * 2 - lengths(2)) * cos(pc_angle);
                            y_sample = (rand(num_sample(i1), 1) * lengths(2) * 2 - lengths(2)) * sin(pc_angle);
                            xys = [x_sample y_sample] + repmat(center, num_sample(i1), 1);
                            thetas = rand(num_sample(i1), 1) * 2 * std_sampling_theta(i3) - std_sampling_theta(i3) + pc_angle;
                            lens = rand(num_sample(i1), 1) * 2 * std_sampling_length(i2) - std_sampling_length(i2) + lengths(1);
                            
                            [grad_mag_diff, grad_mag_std, grad_dir_std, overlap_ratio] = Retrieve_Grad_Mag_Diff(im_grad, im_dir, lc, xys, thetas, lens);
                            %                 cost = -grad_mag_diff + param.energy_alpha*grad_mag_std + param.energy_beta*grad_dir_std...
                            %                     -param.energy_gamma * overlap_ratio;
                            cost = -grad_mag_diff + alpha(i4)*grad_mag_std + beta(i5)*grad_dir_std...
                                -gamma(i6) * overlap_ratio;
                            [min_val, min_ind] = min(cost);
                            
                            usage_count = usage_count + 1;
                            
                            %     lc_state(usage_count,:) = [center, reg_angle, pc_angle, lengths(1), lengths(2)]; % [length, width, angle,]
                            lc_state(usage_count,:) = [xys(min_ind,:), reg_angle, thetas(min_ind), lens(min_ind), lengths(2)]; % [length, width, angle,]
                            
                            x1 = center + lengths(1)*[cos(reg_angle) sin(reg_angle)];
                            x2 = center - lengths(1)*[cos(reg_angle) sin(reg_angle)];
                            
                            x3 = center + lengths(1)*[cos(pc_angle) sin(pc_angle)];
                            x4 = center - lengths(1)*[cos(pc_angle) sin(pc_angle)];
                            
                            x5 = center + lengths(1)*[cos(reg_angle1) sin(reg_angle1)];
                            x6 = center - lengths(1)*[cos(reg_angle1) sin(reg_angle1)];
                            
                            x7 = xys(min_ind,:) + lens(min_ind)*[cos(thetas(min_ind)) sin(thetas(min_ind))];
                            x8 = xys(min_ind,:) - lens(min_ind)*[cos(thetas(min_ind)) sin(thetas(min_ind))];
                            
                            if bDrawRet
                                %         plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-'); % updated
                                plot([x3(1) x4(1)], [x3(2) x4(2)], 'g-', 'linewidth', 1); % pc
                                plot([x7(1) x8(1)], [x7(2) x8(2)], 'r-', 'linewidth', 1); % optimized
                                
                                
                                %         for k = 1:num_sample
                                %             figure(100); clf; imagesc(im_grad); colormap gray; hold on;
                                %             plot([x3(1) x4(1)], [x3(2) x4(2)], 'g-', 'linewidth', 1); % pc
                                %             x7 = xys(k,:) + lens(k)*[cos(thetas(k)) sin(thetas(k))];
                                %             x8 = xys(k,:) - lens(k)*[cos(thetas(k)) sin(thetas(k))];
                                %             plot([x7(1) x8(1)], [x7(2) x8(2)], 'r--', 'linewidth', 1); % optimized
                                %             pause(.01);
                                %         end
                            end
                        end
                        saveas(gcf, sprintf('./line_cand__num_sam(%d)__len(%d)__angle(%.4f)__alpha(%d)__beta(%d)__gamma(%d).png',...
                            num_sample(i1), std_sampling_length(i2), std_sampling_theta(i3), alpha(i4), beta(i5), gamma(i6)),...
                            'png');
                        hold off;
                    end
                end
            end
        end
    end
end
if bDrawRet,    hold off;   end

region_candi_map = pixel_used;

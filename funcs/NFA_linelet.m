function NFA_score = NFA_linelet(pt_ori, im_dir)
    TABSIZE = 100000;
    num_pts = numel(pt_ori);
    num_aligned = 0;
    ang_pivot = pt_ori(1);
    
    % 1st version
    %num_aligned = sum(bAngleAligned( ang_pivot, pt_ori, pi/8 ) > 0);
    
    % 2nd version
    num_aligned = sum(bAngleAligned( pt_ori(1:end-1), pt_ori(2:end), pi/8 ) > 0);
    
    
    inv = zeros(TABSIZE,1);  % /* table to keep computed inverse values */
    tolerance = 0.1; %      /* an error of 10% in the result is accepted */
    logNT = 4.0 * ( log10( size(im_dir, 2) ) + log10( size(im_dir, 1) ) ) / 2.0;
    M_LN10 = 2.30258509299404568402;
    p = 45 / 180;
    
    %   /* trivial cases */
    if num_pts == 0 || num_aligned == 0
        NFA_score = -logNT;
        return;
    end
    
    if num_pts == num_aligned
        NFA_score = -logNT - double(num_pts) * log10(p);
        return;
    end
    
    %   /* probability term */
    p_term = p / (1.0-p);
    log1term = log_gamma( num_pts + 1.0 ) - log_gamma( num_aligned + 1.0 )...
        - log_gamma( (num_pts-num_aligned) + 1.0 )...
        + num_aligned * log(p) + (num_pts-num_aligned) * log(1.0-p);
    term = exp(log1term);
    
    %   /* in some cases no more computations are needed */
    if double_equal(term,0.0)              %/* the first term is almost zero */
        if num_aligned > num_pts*p     %/* at begin or end of the tail?  */
            NFA_score = -log1term / M_LN10 - logNT;
            return;  %/* end: use just the first term  */
        else
            NFA_score = -logNT;
            return;                      %/* begin: the tail is roughly 1  */
        end
    end
    
    %   /* compute more terms if needed */
    bin_tail = term;
    for i = num_aligned+1:num_pts
        
        if i<TABSIZE
            if inv(i) == 0.0
                inv(i) = 1.0 / double(i);
            end
            tmp = inv(i);
        else
            tmp = 1.0 / double(i);
        end
        bin_term = (num_aligned-i+1) * tmp;
        
        mult_term = bin_term * p_term;
        term = term * mult_term;
        bin_tail = bin_tail + term;
        
        if bin_term < 1.0
            err = term * ( ( 1.0 - mult_term^(num_pts-i+1) ) / (1.0-mult_term) - 1.0);
            if err < tolerance * abs(-log10(bin_tail)-logNT) * bin_tail
                break;
            end
        end
    end
    
    NFA_score = -log10(bin_tail) - logNT;
end


function ret = log_gamma(x)
    if x > 15.0
        ret = log_gamma_windschitl(x);
    else
        ret = log_gamma_lanczos(x);
    end
end

function ret = log_gamma_lanczos(x)
    q = [ 75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511 ];
    a = (x+0.5) * log(x+5.5) - (x+5.5);
    b = 0.0;
    
    for n=1:7
        a = a - log( x + n );
        b = b + q(n) * x^n;
        
        ret = a + log(b);
    end
end

function ret = log_gamma_windschitl(x)
    ret = 0.918938533204673 + (x-0.5)*log(x) - x + 0.5*x*log( x*sinh(1/x) + 1/(810.0*x^6));
end

function ret = double_equal(a, b)
    if a == b, ret = 1; end
    
    abs_diff = abs(a-b);
    aa = abs(a);
    bb = abs(b);
    
    if aa > bb
        abs_max = aa;
    else
        abs_max = bb;
    end
    
    abs_max = max(abs_max, realmin('double'));
    
    ret = (abs_diff / abs_max) <= (100.0 * eps);
end
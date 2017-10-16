function brs_candi = Generate_Bresenham_Numbers(ang, ang_var, bX_major_oct)

% bin size
num_bin = ceil((ang_var - ang) * 180 / pi) * 2 + 1;

ang_rng = ang-ang_var:pi/180:ang+ang_var;

% arrGrad = tan(ang_rng);

arrDy = sin(ang_rng);



function offsets = get_move_offset(ang_val)

if     ang_val >= -pi/8   && ang_val <  pi/8,   offsets = [ 1  1;  1  0;  1 -1];
elseif ang_val >=  pi/8   && ang_val <  pi*3/8, offsets = [ 0  1;  1  1;  1  0];
elseif ang_val >=  pi*3/8 && ang_val <  pi*5/8, offsets = [-1  1;  0  1;  1  1];
elseif ang_val >=  pi*5/8 && ang_val <  pi*7/8, offsets = [-1  0; -1  1;  0  1];
elseif ang_val >=  pi*7/8,                      offsets = [-1 -1; -1  0; -1  1];
elseif ang_val >= -pi*3/8 && ang_val < -pi/8,   offsets = [ 1  0;  1 -1;  0 -1];
elseif ang_val >= -pi*5/8 && ang_val < -pi*3/8, offsets = [-1 -1;  0 -1;  1 -1];
elseif ang_val >= -pi*7/8 && ang_val < -pi*5/8, offsets = [-1  0; -1 -1;  0 -1];
elseif ang_val <  -pi*7/8,                      offsets = [-1 -1; -1  0; -1  1];
end

% if ang_val     ==  0,                              offsets = [ 1  0];
% elseif ang_val  >  0      && ang_val <  pi/4,   offsets = [ 1  0;  1  1];
% elseif ang_val ==  pi/4,                           offsets = [ 1  1];
% elseif ang_val  >  pi/4   && ang_val <  pi/2,   offsets = [ 1  1;  0  1];
% elseif ang_val ==  pi/2,                           offsets = [ 0  1];
% elseif ang_val  >  pi/2   && ang_val <  pi*3/4, offsets = [ 0  1; -1  1];
% elseif ang_val ==  pi*3/4,                         offsets = [-1  1];
% elseif ang_val  >  pi*3/4 && ang_val <  pi,     offsets = [-1  1; -1  0];
% elseif ang_val ==  pi,                             offsets = [-1  0];
% elseif ang_val  <  0      && ang_val > -pi/4,   offsets = [ 1  0;  1 -1];
% elseif ang_val == -pi/4,                           offsets = [ 1 -1];
% elseif ang_val  < -pi/4   && ang_val > -pi/2,   offsets = [ 1 -1;  0 -1];
% elseif ang_val == -pi/2,                           offsets = [ 0 -1];
% elseif ang_val  < -pi/2   && ang_val > -pi*3/4, offsets = [ 0 -1; -1 -1];
% elseif ang_val == -pi*3/4,                         offsets = [-1 -1];
% elseif ang_val  < -pi*3/4 && ang_val > -pi,     offsets = [-1 -1; -1  0];
% elseif ang_val == -pi,                             offsets = [-1  0];
% end
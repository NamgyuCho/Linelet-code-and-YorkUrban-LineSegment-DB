function [ ome, gamU, gamL ] = LineletOmegaGamma( inst, ll_type, size_im )
%LINELETOMEGAGAMMA Summary of this function goes here
%   Detailed explanation goes here

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
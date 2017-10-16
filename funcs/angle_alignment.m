function ret = angle_alignment(a, b, threshold)

c = b - a;

if c < 0
    c = -c;
end

if c > pi*3/2
    c = c - 2*pi;
    if c < 0
        c = -c;
    end
end

ret = c < threshold;
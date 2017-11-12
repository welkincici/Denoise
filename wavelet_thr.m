function [ clear_sig ] = wavelet_thr( raw )
%   –°≤®”≤„–÷µ»•‘Î
[d, a] = wavedec(raw, 10, 'db4');
plot(d);
d(d<0.2) = 0;
clear_sig = waverec(d,a, 'db4');
end


function [ clear_sig ] = wavelet_thr( raw )
%   С��Ӳ��ֵȥ��
[d, a] = wavedec(raw, 10, 'db4');
d(d<0.2) = 0;
clear_sig = waverec(d,a, 'db4');
end


function [ noisy_signal ] = add_noise( signal, power )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
N = length(signal);
noisy_signal = signal + wgn(N, 1, power);
end
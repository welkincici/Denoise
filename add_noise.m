function [ noisy_signal ] = add_noise( signal, power )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
N = length(signal);
noisy_signal = signal + wgn(N, 1, power);
end
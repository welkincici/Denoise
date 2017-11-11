function [ clear_sig ] = spectral_sub( raw, alpha, beta, gamma, filter )
%Spectral Subtraction
%   此处显示详细说明
mode = 1;
if nargin < 5
    if nargin == 1
        alpha = 4.5;
        beta = 0.02;
        gamma = 1;
    end
    mode = 0;
    filter=zeros(length(raw),1);
end
% WINDOW

N = min(length(raw),length(filter));

lamda =1024;
shift = 920;
steps = floor((N-lamda)/shift);
N=steps*shift+lamda;
raw = raw(1:N);
filter=filter(1:N);
win = hamming(lamda);
segments=zeros(lamda,steps);
clr_pow = zeros(lamda,steps);

% Spectral Subtractor
j = sqrt(-1);

clear_sig = zeros(length(raw), 1);

% Cut up
for i=1:steps
    start = (i-1)*shift;
    segments(:,i)=raw(start+1:start+lamda).*win;
end

seg_fft = fft(segments);
seg_lat = abs(seg_fft);
seg_ang = angle(seg_fft);
seg_pow = seg_lat.^gamma;

if mode == 0
%     noise_pow = mean(seg_pow(:,1:20), 2);
%     disp(mean(noise_pow));
%     plot(noise_pow);
%     hold on
    noise_pow = mean(seg_pow(:,1:20), 2);
    for i = 1:lamda
        clr = seg_pow(i,:) - alpha*noise_pow(i);
        clr(clr<beta*noise_pow(i)) = beta*noise_pow(i);
        clr_pow(i,:) = clr;
    end
else
    start_index=1;
    off=1;
    flag = 0;
    for i = 1:N
        if filter(i)==flag
            off=off+1;
        else
            if flag==0
                %计算noise_pow
                end_index=floor((off-lamda)/shift);
                noise_pow = mean(seg_pow(:,start_index:end_index), 2);
                
%                 noise_pow = mean(seg_pow(:,1:20), 2);
%                 disp(mean(noise_pow));
%                 plot(noise_pow);
%                 clr_pow(:,start_index:end_index)=zeros(lamda,end_index-start_index+1);
            else
                %计算谱减
                end_index=ceil((off-lamda)/shift);
%                 disp(start_index);
%                 disp(end_index);
%                 clr = seg_pow(:,start_index:end_index);
%                 for k = 1:lamda
%                     clr(k,:) = clr(k,:) - alpha*noise_pow(k);
%                     clr(clr<beta*noise_pow(k)) = beta*noise_pow(k);
%                 end
%                 clr_pow(:,start_index:end_index) = clr;
%                 figure(start_index)
%                 plot(noise_pow);
                for k = 1:lamda
%                     disp(start_index);
%                     disp(end_index);
                    clr = seg_pow(k,start_index:end_index) - alpha*noise_pow(k);
                    clr(clr<beta*noise_pow(k)) = beta*noise_pow(k);
                    clr_pow(k,start_index:end_index) = clr;
                end
            end
            flag=1-flag;
            start_index=end_index+1;
        end
    end
end
% Invert
clr_fft = (clr_pow.^(1/gamma)).*exp(j*seg_ang);
clr_seg = real(ifft(clr_fft));
wins=clear_sig;
win2=win.*win;
for i=1:steps
    start = (i-1)*shift;
%     clear_sig(start+1:start+lamda) = clear_sig(start+1:start+lamda) + clr_seg(:,i).*win;
%     wins(start+1:start+lamda)=wins(start+1:start+lamda)+win;
    wins(start+1:start+lamda)=wins(start+1:start+lamda)+win2;
    clear_sig(start+1:start+lamda) = clear_sig(start+1:start+lamda) + clr_seg(:,i).*win;
end

clear_sig=clear_sig./wins;
clear_sig(isnan(clear_sig))=0;
end


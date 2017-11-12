function [ clear_sig ] = spectral_sub( raw, alpha, beta, filter )
%Spectral Subtraction
mode = 1;
if nargin < 4
    if nargin == 1
        alpha = 4.5;
        beta = 0.02;
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
clear_sig = zeros(length(raw), 1);
j = sqrt(-1);

% Cut up
for i=1:steps
    start = (i-1)*shift;
    segments(:,i)=raw(start+1:start+lamda).*win;
end

seg_fft = fft(segments);
seg_lat = abs(seg_fft);
seg_ang = angle(seg_fft);
seg_pow = seg_lat.^2;

yitas=0.975;
yitad=0.96;
b=1;
miu=0.2;
betaalpha=0.5;
if mode == 0
    noise_pow = mean(seg_pow(:,1:20), 2);
    for i = 1:lamda
        clr = b*seg_pow(i,:) - alpha*noise_pow(i);
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
                %¼ÆËãnoise_pow
                noise_pow = seg_pow(:,start_index);
                end_index=floor((off-lamda)/shift);
                for index=start_index:end_index
                    noise_pow = (1-yitad)*seg_pow(:,index)+yitad*noise_pow;
                end
            else
                %¼ÆËãÆ×¼õ
                end_index=ceil((off-lamda)/shift);
                yip=(seg_pow(:,start_index)-noise_pow)./noise_pow;
                for index=start_index:end_index
                    for k=1:lamda
                        yip(k)=yitas*yip(k)+(1-yitas)*max(seg_pow(k,index)-noise_pow(k),0)/noise_pow(k);
                        clr_pow(k,index)=(yip(k)^2)/(yip(k)^2+betaalpha)*(seg_pow(k,index)-noise_pow(k));
                        miuY=0.25*(miu*miu*seg_pow(k,index)+clr_pow(k,index-1)+2*miu*sqrt(seg_pow(k,index)*clr_pow(k,index-1)));
                        if clr_pow(k,index)<miuY
                            clr_pow(k,index)=miuY;
                        end
                        yip(k)=clr_pow(k,index)/noise_pow(k);
                    end
                end
%                 for k = 1:lamda
%                     clr = b*seg_pow(k,start_index:end_index) - alpha*noise_pow(k);
%                     clr(clr<beta*noise_pow(k)) = beta*noise_pow(k);
%                     clr_pow(k,start_index:end_index) = clr;
%                 end
            end
            flag=1-flag;
            start_index=end_index+1;

        end
    end
end

% Invert
clr_fft = (clr_pow.^(0.5)).*exp(j*seg_ang);
clr_seg = real(ifft(clr_fft));
wins=clear_sig;
win2=win.*win;
for i=1:steps
    start = (i-1)*shift;
    wins(start+1:start+lamda)=wins(start+1:start+lamda)+win2;
    clear_sig(start+1:start+lamda) = clear_sig(start+1:start+lamda) + clr_seg(:,i).*win;
end

clear_sig=clear_sig./wins;
clear_sig(isnan(clear_sig))=0;
end


function [ clear_sig ] = gabor_thr( raw, est_sigma, lamda, shift )
%Gabor Transfomation & Blackman Window & Threshold
%   此处显示详细说明
if nargin == 2
    lamda =1024;
    shift = 65;
end

N = length(raw);
med_record = 20;
win = blackman(lamda);
steps = floor((N-lamda)/shift);


count = 0;
meds = zeros(med_record, 1);
clear_sig = zeros(steps*shift+lamda,1);

%MAIN LOOP
while count<=steps
    offset = count*shift;
    segment = raw(offset+1:offset+lamda);
    
%   Garbo Transform
    gar_sig = fft(segment.*win);
    
%   Estimate standard deviation
    med = median(imag(gar_sig(floor(lamda*3/4):lamda)));
    meds = [meds(2:20); med];
    if count>=med_record
        est_sigma = mean(meds)/(0.6745*0.55*lamda/4);
    end
    
%   Hard Thresholding
    T = 0.55*est_sigma*sqrt(lamda*log(lamda));
    gar_sig(abs(gar_sig)<=T) = 0;
    
%   Invert Gabor Transform
    clear_seg = ifft(gar_sig).*win;
    clear_sig(offset+1:offset+lamda) = clear_sig(offset+1:offset+lamda) + clear_seg;
    
    count = count+1;
end

end


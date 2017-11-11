% CONST
[raw, fsn] = audioread('«ø‘Î…˘“Ù∆µ1.wma');
raw = raw(:,1);
N = length(raw);
clear_sig = zeros(N, 1);
j = sqrt(-1);

lamda =1024;
shift = lamda/2;
steps = ceil((N-lamda)/shift);
raw = [raw; zeros(steps*shift+lamda-N,1)];
win = hamming(lamda);
segments=zeros(lamda,steps);
count = 0;

% Main Loop
while count<=steps
    offset = count*shift;
    
%    Hamming Window & FFT
    seg_fft = fft(raw(offset+1:offset+lamda).*win);
    seg_ang = angle(seg_fft);
    seg_power = abs(seg_fft).^2;
    
%     if count<=20
%         noise_power = (noise_power*count + mean(seg_power))/(count+1);
%     end
    
%     Subtract
    clear_power = seg_power - noise_power;
    clear_power(clear_power<0) = 0;
    clear_seg = real(ifft(sqrt(clear_power).*exp(j*seg_ang)));
    
%     Invert
    clear_sig(offset+1:offset+lamda) = clear_sig(offset+1:offset+lamda) + clear_seg.*win;
    
    count = count+1;
    
end

subplot(2,1,1);
plot(clear_sig);
subplot(2,1,2);
plot(raw);




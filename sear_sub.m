% CONST
[raw, fsn] = audioread('«ø‘Î…˘“Ù∆µ1.wma');
raw = raw(:,1);
N = length(raw);

% WINDOW
lamda =1024;
shift = lamda/2;
steps = ceil((N-lamda)/shift);
raw = [raw; zeros(steps*shift+lamda-N,1)];
win = hamming(lamda);
segments=zeros(lamda,steps);

% Spectral Subtractor
j = sqrt(-1);
alpha = 2;
beta = 0.2;

% Noise Estimate
HIGH = 21;
LOW = 14;

clear_sig = zeros(length(raw), 1);
clr_fft = zeros(lamda,steps);

% Cut up
for i=1:steps
    start = (i-1)*shift;
    segments(:,i)=raw(start+1:start+lamda);
end

seg_fft = fft(segments);
seg_lat = abs(seg_fft);
seg_pow = seg_lat.^2;

low_mean = zeros(steps,1);
for i=1:steps
    low_mean(i) = mean(seg_lat(4:15,i));
end

index = 1;
end_index = 1;
while index<=steps
    % Search voice & Estimate noise
    if(low_mean(index)>HIGH)
        start_index = index;
        while low_mean(start_index)>LOW
            start_index = start_index-1;
        end
        
        noise_power = mean(seg_pow(:,end_index:start_index), 2);
        
        end_index=index;
        while low_mean(end_index)>LOW
            end_index = end_index+1;
        end
        
        % Spectrum subtract
        for i= start_index:end_index
            for k = 1:lamda
                clr_pow = seg_pow(k,i) - alpha*noise_power(k);
                if clr_pow < beta*noise_power(k)
                    clr_pow = beta*noise_power(k);
                end
                clr_fft(k,i) = sqrt(clr_pow)*exp(j*angle(seg_fft(k,i)));
            end
        end
        
        index=end_index;
    end
    index = index+1;
end

% Invert
clr_seg = real(ifft(clr_fft));
for i=1:steps
    start = (i-1)*shift;
    clear_sig(start+1:start+lamda) = clear_sig(start+1:start+lamda) + clr_seg(:,i).*win;
end

plot(clear_sig);

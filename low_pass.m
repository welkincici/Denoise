[raw, fsn] = audioread('my_sample5.wav');
raw = raw(:,1);
N = length(raw);
raw_fft = fft(raw, 1024);
for i=1:1024
    if i<4 || i> 15
        raw_fft(i)=0;
    end
end

sig = real(ifft(raw));
plot(sig);
[noise1, fsn1] = audioread('noise1.mp3');
[noise2, fsn2] = audioread('noise2.mp3');
[s, fss] = audioread('»ı‘Î…˘“Ù∆µ1.wma');
[raw, fsr] = audioread('«ø‘Î…˘“Ù∆µ1.wma');
s = s(:,1);
raw = raw(:,1);
s_fft=fft(s);
raw_fft=fft(raw);
noise1 = noise1(:,1);
noise2 = noise2(:,1);
noise=[noise1;noise2];
disp(mean(abs(noise).^2));
disp(mean(abs(raw).^2));
noise_fft=fft(noise);
% disp('Noise FFT Average');
% disp(mean(abs(noise_fft)));
% disp(mean(real(noise_fft)));
% disp(mean(imag(noise_fft)));
% disp('Raw FFT Average');
% disp(mean(real(raw_fft)));
% disp(mean(imag(raw_fft)));
% disp('S FFT Average');
% disp(mean(real(s_fft)));
% disp(mean(imag(s_fft)));
% disp('Variance');
% disp(mean(abs(noise_fft).^2));

% disp(sqrt(mean(real(noise_fft).^2)));
% disp(sqrt(mean(imag(noise_fft).^2)));
% 
% N = length(noise);
% w = blackman(N);
% win_noise = w.*noise(:,1);
% disp('Windowed Noise Average');
% disp(mean(win_noise));
% disp('Windowed Noise Variance/¶“');
% disp(sqrt(mean(win_noise.^2))/sqrt(mean(noise(:,1).^2)));
% disp(sqrt(mean(noise.^2)))


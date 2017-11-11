[noise1, fsn1] = audioread('noise1.mp3');
[noise2, fsn2] = audioread('noise2.mp3');
noise1 = noise1(:,1);
noise2 = noise2(:,2);

% disp('Noise FFT Average');
% disp(mean(abs(noise_fft)));
% disp(mean(real(noise_fft)));
% disp(mean(imag(noise_fft)));
% disp('Variance');
% disp(sqrt(mean(abs(noise_fft).^2)));
% disp(sqrt(mean(real(noise_fft).^2)));
% disp(sqrt(mean(imag(noise_fft).^2)));
% 
% N = length(noise);
% w = blackman(N);
% win_noise = w.*noise(:,1);
% disp('Windowed Noise Average');
% disp(mean(win_noise));
% disp('Windowed Noise Variance/¦Ò');
% disp(sqrt(mean(win_noise.^2))/sqrt(mean(noise(:,1).^2)));
% disp(sqrt(mean(noise.^2)))


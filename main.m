[s, fss] = audioread('»ı‘Î…˘“Ù∆µ1.wma');
[raw, fsr] = audioread('«ø‘Î…˘“Ù∆µ1.wma');
% [noise, fsn] = audioread('noise2.mp3');
s = s(:,1);
raw = raw(:,1);
noisy_s=add_noise(s,-30);
% noise = noise(:,1);
clear_sig = spectral_sub(raw, 3, 0, 2);
filter = wavelet_scale_filer(clear_sig);
res=spectral_sub(raw, 3, 0.02, 2, filter);

% subplot(3,1,1);

% subplot(3,1,2);
subplot(4,1,1);
plot(noisy_s);
subplot(4,1,2);
plot(spectral_sub(noisy_s, 1, 0, 2));
subplot(4,1,3);
plot(spectral_sub(noisy_s, 3, 0, 2));
subplot(4,1,4);
plot(spectral_sub(noisy_s, 3, 0.02, 2));
% hold on
% plot(res);

% subplot(3,1,3);
% 
% plot(filter);






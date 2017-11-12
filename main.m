[s, fss] = audioread('ÈõÔëÉùÒôÆµ1.wma');
[raw, fsr] = audioread('Ç¿ÔëÉùÒôÆµ2.wma');
[noise, fsn] = audioread('noise2.mp3');
s = s(:,1);
raw = raw(:,1);
noisy_s=add_noise(s,-30);
noise = noise(:,1);
clear_sig = spectral_sub(raw, 3, 0.02);
filter = wavelet_scale_filer(raw);

res=spectral_sub(raw, 3, 0.02, filter);
% subplot(2,1,1);
% plot(raw);
% subplot(2,1,2);
% plot(res);

% hold on 
% plot(filter);
% sound(clear_sig,fsr);







[s, fss] = audioread('��������Ƶ1.wma');
[raw, fsr] = audioread('ǿ������Ƶ1.wma');
[noise, fsn] = audioread('noise2.mp3');
s = s(:,1);
raw = raw(:,1);
noisy_s=add_noise(s,-30);
noise = noise(:,1);
clear_sig = spectral_sub(raw, 3, 0.02, 2);
filter = wavelet_scale_filer(raw);

res=spectral_sub(raw, 3, 0.02, 2, filter);
plot(res);
hold on 
plot(filter);
% sound(clear_sig,fsr);







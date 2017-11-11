function [ res ] = wavelet_scale_filer( signal, wave )
%   小波尺度空间相关滤波器
if nargin < 2
    wave = 'db4';
end
start_layer = 3;
end_layer = 6;
layers = end_layer - start_layer + 1;
base = 2^(end_layer + 1); 
N = floor(length(signal)/base)*base;
signal = signal(1:N);

[A, D] = swt(signal, end_layer + 1, wave);

corr = zeros(layers, N);

for i = 1:layers
    corr(i, :)=D(i + start_layer - 1, :).*D(i + start_layer, :);
end

ratio = sqrt(sum(D(start_layer:end_layer, :).^2, 2)./sum(corr.^2, 2));

% disp(corr.^2);
% disp(sum(D(start_layer:end_layer, :).^2, 2));
% disp(sum(corr.^2, 2));
% disp(ratio);
for i = 1:layers
    corr(i, :)=abs(corr(i, :)*ratio(i));
    corr(i,:)=corr(i,:)./mean(corr(i,:));
%     for k = 1:N
%         if abs(corr(i, k))>abs(D(i + start_layer - 1, k))
%             my_filter(i, k)=1;
%         end
%     end
end

lamda = 2048;
shift = 1024;
steps = ceil((N-lamda)/shift);
detect=zeros(layers,steps);
for i=1:layers
    for k=1:steps
        start = (k-1)*shift;
        detect(i,k)=mean(corr(i,start+1:start+lamda)); %什么是最佳的判定方案？
    end
end
filter = zeros(layers,N);
for i=1:layers
    for k=1:steps
        if detect(i,k)>2
            start=(k-1)*shift;
            if start>lamda
                on=start-lamda+1;
            else
                on=1;
            end
            if start>N-2*lamda
                off=N;
            else
                off=start+2*lamda;
            end
            filter(i,on:off)=ones(1,off-on+1);
        end
    end
end
% subplot(2,1,1);
% plot(corr(4, :));
% hold on
% plot(filter(4, :));
% subplot(2,1,2);
% plot(detect(4,:));

% subplot(4,1,1);
% plot(corr(1, :));
% subplot(4,1,2);
% plot(corr(2, :));
% subplot(4,1,3);
% plot(corr(3, :));
% subplot(4,1,4);
% plot(corr(4, :));
% plot(sum(abs(corr)));
res = sum(filter);
res(res>=1)=1;
% plot(signal);
% hold on
% plot(res);
end


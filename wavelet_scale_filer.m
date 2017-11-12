function [ res ] = wavelet_scale_filer( signal, yita, wave )
%   小波尺度空间相关滤波器
if nargin < 3
    wave = 'db4';
    if nargin==1
        yita = 0.9;
    end
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
for i = 1:layers
    corr(i, :)=abs(corr(i, :).*ratio(i)./D(i + start_layer - 1, :));
%     subplot(layers,1,i);
%     plot(corr(i,:));
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
    subplot(layers,1,i);
    x=plot(detect(i, :));
    x.LineWidth=3;
end
filter = zeros(layers,N);
expand=5;
for i=1:layers
    k=expand+1;
    while k<=steps-expand
        if detect(i,k)>yita && detect(i,k)==max(detect(i,k-expand:k+expand))
            for ahead=1:min(30,k-expand-1)
                start_index=k-ahead;
                if detect(i,start_index)==min(detect(i,start_index-expand:start_index+expand))
                    break;
                end
            end
            for back=1:min(30,steps-k-expand)
                end_index=k+back;
                if detect(i,end_index)==min(detect(i,end_index-expand:end_index+expand))
                    break;
                end
            end
            k=end_index;
            on=(start_index-1)*shift+lamda;
            off=(end_index-1)*shift;
            filter(i,on:off)=ones(1,off-on+1);
        end
        k=k+1;
    end
end
% for i=1:layers
%     subplot(layers,1,i);
%     plot(corr(i, :));
%     hold on
%     x=plot(filter(i, :).*4);
%     x.LineWidth=3;
%     hold off
% end
res = sum(filter);
res(res>=1)=1;
% plot(signal);
% hold on
% plot(res);
end


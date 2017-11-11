% layers = size(my_VAR_1, 1);
% N = size(my_VAR_1, 2);
% corr = zeros(layers, N);
% my_filter = zeros(layers, N);
% for i = 1:layers-1
%     corr(i, :)=my_VAR_1(i,:).*my_VAR_1(i+1,:);
% end
% corr(layers, :)=corr(layers-1, :);
% 
% ratio = sqrt(sum(my_VAR_1.^2, 2)./sum(corr.^2, 2));
% for i = 1:layers
%     corr(i, :)=corr(i, :)*ratio(i);
%     for k = 1:N
%         if abs(corr(i,k))>abs(my_VAR_1(i,k))
%             my_filter(i,k)=1;
%         end
%     end
% end
% figure(5);
% subplot(4,1,1);
% plot(my_filter(1,:));
% subplot(4,1,2);
% plot(my_filter(2,:));
% subplot(4,1,3);
% plot(my_filter(3,:));
% subplot(4,1,4);
% plot(my_filter(4,:));
% figure(6);
% subplot(4,1,1);
% plot(my_filter(5,:));
% subplot(4,1,2);
% plot(my_filter(6,:));
% subplot(4,1,3);
% plot(my_filter(7,:));
% subplot(4,1,4);
% plot(my_filter(8,:));
% figure(7);
% subplot(4,1,1);
% plot(my_filter(9,:));
% subplot(4,1,2);
% plot(my_filter(10,:));
% subplot(4,1,3);
% plot(my_filter(11,:));
% subplot(4,1,4);
% plot(my_filter(12,:));
plot(wavelet_scale_filer(clear_sig, 9));


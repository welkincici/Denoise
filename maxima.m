end_layer = 7;
base = 2^(end_layer + 1); 
N = floor(length(raw)/base)*base;
raw = raw(1:N);

[A, D] = swt(raw, end_layer, 'db4');

maxi = zeros(size(D));

for i=1:end_layer
    pre=D(i,1);
    now=D(i,2);
    next=D(i,3);
    for k=2:N-2
        if abs(now)>abs(pre) && abs(now)>abs(next)
            maxi(i,k)=now;
        end
        pre=now;
        now=next;
        next=D(i,k+2);
    end
end
for i=1:end_layer
    figure(i);
    plot(maxi(i,:));
end
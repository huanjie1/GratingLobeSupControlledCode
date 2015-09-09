% 验证总辐射功率随阵元间距的变化
d=0:0.0125:30;
aim=90;%deg
pall=zeros(1,length(d));
N=16;
% if mod(N,2)
    
for ii=1:length(d)
    func1=@(x)2*pi*(abs(1+exp(1i*pi*d(ii)*(sin(x)-sin(aim/180*pi)))+exp(-1i*pi*d(ii)*(sin(x)-sin(aim/180*pi)))).^2).*cos(x);
    pall(ii)=integral(func1,-pi/2,pi/2);
end
plot(d,pall);hold on;
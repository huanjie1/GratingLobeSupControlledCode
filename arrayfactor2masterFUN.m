function [fomcmpue]=arrayfactor2masterFUN( N, aimdegree0, ratio, centerfreq, fbw)
% arrayfactor2masterFUN.m
% arrayfactor2masterFUN( 16, 60, 1.5, 10e9, 0.4)
% 绘制二维阵列因子，arrayfactor2master.m的函数封装形式，用于扫描阵列特性（而非波束形成网络）
% 被AFmaster.m调用
% 调用arrayfactorangFORgeneral.m, arrayfactorangFORgeneralRR.m
% clear
% N=16;
c=299792458;
diffindex=0:(N-2);
% centerfreq=10e9;
centerlambda=c/centerfreq;

spacing0=ratio*centerlambda;

devimode='squ';
% devimode='squ';*7/2

if strcmp(devimode,'squ')
    spacingdia=1*0.0000*diffindex.^2;%deviation from even spacing
    spacings=spacingdia-(min(spacingdia)+max(spacingdia))/2+spacing0;
else
    spacings=randn(1,N-1)*spacing0*0.2+spacing0;
end

% spacings(1)=(spacings(1)+spacings(2))/(1+pi/2);
% spacings(2)=(spacings(1)+spacings(2))/(1+pi/2)*pi/2;

if sum(spacings<0)>0
    error('wrong spacingdia');
end
xposition0=[0 cumsum(spacings)];
xposition=xposition0-(min(xposition0)+max(xposition0))/2;
% figure;stem(xposition,max(spacings)*ones(1,length(xposition)));hold on
% plot(linspace(-max(xposition),max(xposition),length(spacings)),spacings)

% aimdegree0=60;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN=8001;
freqaxis=linspace(-40e9,40e9,NN);
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

% dgraxis=linspace(-90,90,37); %for waterfall
dgraxis=linspace(-90,90,721);

% sigt1=cos( 2*pi*10e9*t + pi*10e9/(ts*NN)*t.^2);
% sigt1=cos( 2*pi*10e9*t).* exp(-(t/0.1e-9).^2);
sigt1=[zeros(1,2000) sigeneratorfor2d( t(2001:6001),  'lfm', fbw*centerfreq, centerfreq ) zeros(1,2000)];
% figure;plot(t,sigt1);
% title('sig waveform');

dgrsection=0;
freqsection=35.1e9;

% arrayfactorangFORgeneral( xposition, freqaxis, dgraxis, t, aimdegree0, dgrsection, freqsection, sigt1, 1 );
% cmpideal=arrayfactorangFORgeneralIDEAL( xposition, freqaxis, dgraxis, t, aimdegree0, centerfreq, dgrsection, freqsection, sigt1, 0, 1 );
cmpue=arrayfactorangFORgeneralRR( xposition, freqaxis, dgraxis, t, aimdegree0, centerfreq, dgrsection, freqsection, sigt1, 0, 0.2 );

fomcmpue=cmpFoMcalcV2( dgraxis, aimdegree0, cmpue, N, spacing0/centerlambda, 0 );
end


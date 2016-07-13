% arrayfactor2master.m
% 绘制二维阵列因子
% 调用arrayfactorangFORgeneral.m，
clear
N=16;
c=299792458;
diffindex=0:(N-2);
centerfreq=10e9;
centerlambda=c/centerfreq;
d0=1*centerlambda;

devimode='squ';
% devimode='squ';

if strcmp(devimode,'squ')
    spacingdia=1*0.00020*diffindex.^2;%deviation from even spacing
    spacings=spacingdia-(min(spacingdia)+max(spacingdia))/2+d0;
else
    spacings=randn(1,N)*d0*0.2+d0;
end


if sum(spacings<0)>0
    error('wrong spacingdia');
end
xposition0=[0 cumsum(spacings)];
xposition=xposition0-(min(xposition0)+max(xposition0))/2;
figure;stem(xposition,max(spacings)*ones(1,length(xposition)));hold on
plot(linspace(-max(xposition),max(xposition),length(spacings)),spacings)

aimdegree0=45;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN=8001;
freqaxis=linspace(-40e9,40e9,NN);
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

dgraxis=linspace(-90,90,721);

% sigt1=cos( 2*pi*10e9*t + pi*10e9/(ts*NN)*t.^2);
% sigt1=cos( 2*pi*10e9*t).* exp(-(t/0.1e-9).^2);
sigt1=sigeneratorfor2d( t,  'gaud0', 10e9, 10e9 );
figure;plot(t,sigt1);
title('sig waveform');

dgrsection=45;
freqsection=10e9;

arrayfactorangFORgeneral( xposition, freqaxis, dgraxis, t, aimdegree0, dgrsection, freqsection, sigt1, 3 );


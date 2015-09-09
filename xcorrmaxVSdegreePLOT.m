function [ gratingsup, xcmax  ] = xcorrmaxVSdegreePLOT( N0, d, aimdgr, bw, mode, fc )
% xcorrmaxVSdegreePLOT.m
% 真时延相关接收的方向图计算，栅瓣情况同时画出
% 栅瓣抑制2_16阵元_4G带宽_10G载波_相关.vsdx

if nargin<6, fc=10e9;end
if nargin<5, mode='lfm';end
if nargin<4, bw=4e9;end
if nargin<3, aimdgr=0;end
if nargin<2, d=0.06;end
if nargin<1, N0=16;end

% global r1;%用于延时和增益非理想性的仿真
% r1=randn(1,N0)*0+1;
% figure;plot(log10(r1)*20);

xcmax=[];
c=3e8;
N=N0;
drgsee0=-60:0.5:60;
for drgsee=drgsee0
    drgsee
    xcmax=[xcmax xcorrTTDarrayr1( aimdgr, drgsee, 0,N,d,bw, mode, fc)];
end

if  aimdgr<0
    deltadrg=(asin(c/fc/d+sin(aimdgr/180*pi))-aimdgr/180*pi)/2*180/pi;
else
    deltadrg=(asin(c/fc/d-sin(aimdgr/180*pi))+aimdgr/180*pi)/2*180/pi;
end

h1=xcmax;
h1(abs(drgsee0-aimdgr)<deltadrg)=0;

drg1=-60:0.25:60;
h0=sin( (N/2)*d/c*(sin(drg1/180*pi+1e-14)-sin(aimdgr/180*pi))*2*pi*fc )...
    ./sin( 0.5*d/c*(sin(drg1/180*pi+1e-14)-sin(aimdgr/180*pi))*2*pi*fc )/N;

% figure;plot(drgsee0,xcmax,drg1,max(real(xcmax))*abs(h0));
% figure;plot(drgsee0,xcmax.^2,drg1,max(real(xcmax)).^2*abs(h0));
% figure;plot(drgsee0,20*log10(xcmax+1e-12),drg1,20*log10(max(real(xcmax))*abs(h0)+1e-12));
figure;plot(drgsee0,20*log10(xcmax./max(real(xcmax))+1e-12),...
            drg1,20*log10(abs(h0)+1e-12));

gratingsup=max(xcmax)/max(h1);

end


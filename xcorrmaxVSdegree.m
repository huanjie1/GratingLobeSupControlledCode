function [ gratingsup, xcmax  ] = xcorrmaxVSdegree( N0, d, aimdgr, bw, mode, fc )
% xcorrmaxVSdegree.m
% 真时延相关接收的方向图计算，同时分析栅瓣抑制比
% 被supsweep.m调用
if nargin<6, fc=10e9;end
if nargin<5, mode='phco';end
if nargin<4, bw=5e9;end
if nargin<3, aimdgr=-60;end
if nargin<2, d=0.06;end
if nargin<1, N0=16;end

xcmax=[];
c=3e8;
N=N0;
drgsee0=-90:1:90;
for drgsee=drgsee0
    drgsee;
    xcmax=[xcmax xcorrTTDarrayr1( aimdgr, drgsee, 0,N,d,bw, mode, fc)];
end

if  aimdgr<0
    deltadrg=(asin(c/fc/d+sin(aimdgr/180*pi))-aimdgr/180*pi)/2*180/pi;
else
    deltadrg=(asin(c/fc/d-sin(aimdgr/180*pi))+aimdgr/180*pi)/2*180/pi;
end

h1=xcmax;
h1(abs(drgsee0-aimdgr)<deltadrg)=0;
plot(drgsee0,10*log10(abs(xcmax)/max(abs(xcmax))));hold on

gratingsup=max(xcmax)/max(h1)

end


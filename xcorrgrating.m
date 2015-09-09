% xcorrgrating.m
% 各阶栅瓣点的强度变化
clear
N0=16;
lambda=0.03;
d=0.16;
aimdrg=-45;
xcorrmaxVSdegreePLOT( N0, d, aimdrg);
N=N0;
xcorrt0array=[];
M0=1:8;
for m=M0
% m=5;
m
seedrg=asin(m*lambda/d+sin(aimdrg/180*pi))*180/pi;
[~,~,~,xcorrt0]=xcorrTTDarrayr1( aimdrg, seedrg, 0.4,N, d)
xcorrt0array=[xcorrt0array  xcorrt0];
end
[~,~,~,expected]=xcorrTTDarrayr1( aimdrg, aimdrg, 0.4,N, d);
figure;plot(M0,xcorrt0array,M0,expected/N0*ones(1,length(M0))); %ylim([1.25e4 1.65e4]);
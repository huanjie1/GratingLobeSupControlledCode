% arrayfactor2master.m
% 绘制二维阵列因子
% 调用arrayfactorangFORgeneral.m, arrayfactorangFORgeneralRR.m
clear
N=8;
c=299792458;
diffindex=0:(N-2);
centerfreq=10e9;
centerlambda=c/centerfreq;
% spacing0=1*centerlambda;
spacing0=0.042;

devimode='squ';
% devimode='squ';

if strcmp(devimode,'squ')
    spacingdia=1*0.0000*diffindex.^2;%deviation from even spacing
    spacings=spacingdia-(min(spacingdia)+max(spacingdia))/2+spacing0;
else
    spacings=randn(1,N-1)*spacing0*0.2+spacing0;
end


if sum(spacings<0)>0
    error('wrong spacingdia');
end
xposition0=[0 cumsum(spacings)];
xposition=xposition0-(min(xposition0)+max(xposition0))/2;
% figure;stem(xposition,max(spacings)*ones(1,length(xposition)));hold on
% plot(linspace(-max(xposition),max(xposition),length(spacings)),spacings)

aimdegree0=-60:0.001:60;
% aimdegree0=29.95:0.00001:30.05;
delaystep=1.2e-3/c;

aimdelay0=xposition.'*sin(aimdegree0/180*pi)/c;
aimdelay=aimdelay0-ones(N,1)*min(aimdelay0);
aimdelayreal=delaystep*round(aimdelay/delaystep);

% figure;imagesc(aimdegree0,1:N,aimdelay);
% xlabel('aimdegree');ylabel('antenna index');title('aim delay(s)');colorbar
% figure;imagesc(aimdegree0,1:N,aimdelay*c*1000);
% xlabel('aimdegree');ylabel('antenna index');title('aim delay(mm)');colorbar


aimdelaydiff=diff(aimdelayreal.').';
aimdegree1=aimdegree0(1:end-1);
figure;imagesc(aimdegree1,1:N,aimdelaydiff);
xlabel('aimdegree');ylabel('antenna index');title('aim delay(mm)');colorbar
aimdelaysetchange=sum(aimdelaydiff,1);
aimdelaysetchange(abs(aimdelaysetchange)>0.5*delaystep)=1;
% figure;stem(aimdegree1,aimdelaysetchange);hold on
figure;plot(aimdegree1(aimdelaysetchange>0.5),[1 diff(aimdegree1(aimdelaysetchange>0.5))])








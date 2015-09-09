function [delt, fomu, numb, ratiof, ration ] = gratinglobe( rbw, N, m )
% 相对带宽  阵元数（偶数） 第几阶栅瓣
% beta=1e-5:0.01:1; [aa,bb,cc,dd,ee]=gratinglobe( beta, 10, 3 );figure;plot(beta,aa,beta,bb,beta,cc,beta,dd,beta,ee)
% beta=1e-9:0.02:1.9; [~,bb,~,dd,~]=gratinglobe( beta, 16, 2 );figure;plot(beta,bb,beta,dd)
% gratinglobe.m
% pic003
% 能量积分：栅瓣点强度与相对带宽的关系(理论公式与数值积分) 

k=round(rbw/2*m);

% main=N*N*rbw/2;
main=N*N;

fomu= 2*(2*k*N/pi/m*sinint(N*pi)...
            +N/pi/m*sinint(2*N*m*pi*(rbw/2-k/m))...
            -1/pi/pi./(m*m*(rbw/2-k/m)).*sin(m*N*pi*(rbw/2-k/m)).^2)./rbw;%功率谱密度为1/rbw 双边需乘2

func=@(x)(sin(N*m*pi*x)./sin(m*pi*x)).^2;
numb=zeros(1,length(rbw));
for ii=1:length(rbw)
    numb(ii)=integral(func,eps,rbw(ii)/2) / rbw(ii)*2;
end

delt=fomu-numb;

ratiof=main./fomu;
ration=main./numb;

figure;plot(rbw,10*log10(fomu/max(fomu)),rbw,10*log10(numb/max(numb)));

end


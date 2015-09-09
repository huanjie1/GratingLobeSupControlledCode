% arrayfactorangFORicocn2015.m
% 能量积分：频率/时域-角度二维图，奇数阵元，ICOCN2015论文用

N=4;
c=3e8;
centerfreq=10e9;
centerlambda=c/centerfreq;
d=1*centerlambda;

degree=linspace(-90,90,360);
theta=degree/180*pi;

degree0=30;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=degree0/180*pi;

f=linspace(0,40e9,8001);
w=2*pi*f;

[ww,thth]=meshgrid(w,theta);
h0=sin( (N+0.5)*d/c*(sin(thth)-sin(theta0)).*ww )...
    ./sin( 0.5*d/c*(sin(thth)-sin(theta0)).*ww )  ...
    .*(abs(thth-theta0)>0.5/180*pi)...
    +(2*N+1)*(abs(thth-theta0)<=0.5/180*pi);
h0(h0>2*N+1)=2*N+1;
h0(h0<-2*N-1)=-2*N-1;
h=abs(h0)/(2*N+1);
% h=h0;
figure;h1=imagesc(w/2/pi/1e9,theta/pi*180,h);xlabel('Frequency/GHz');
% imagesc(w/2/pi/c*d,theta/pi*180,h);xlabel('d/{\lambda}');
ylabel('\it\theta /{\circ}');
colorbar;
set(gcf,'outerposition',get(0,'screensize'));
axt=axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           ...'YAxisLocation','none',...
           'Color','none',...
           'XColor','k');
xlim([0,max(w)/2/pi/c*d]);xlabel('\itd / {\it\lambda}','fontname','Times New Roman');
set(gca,'ytick',[])

figure;
finstr=13e9;%切片位置
[dfmin,findex]=min(abs(f-finstr));
h1=plot(theta/pi*180,h(:,findex),'linewidth',2);ylim([0,1.05]);
ylabel('Normalized AF');xlabel('{\it\theta} / {\circ}');
% xlim([-90,90]);
set(gca,'XTick',-90:45:90);
set(gca,'YTick',0:0.5:1);
grid on; set(gcf,'Position',[400,400,300,300]);

% figure;
% fstart=5e9;%宽带起止点
% fend=15e9;
% [dfsmin,fsindex]=min(abs(f-fstart));
% [dfemin,feindex]=min(abs(f-fend));
% energypatten=sum(h(:,fsindex:feindex).^2,2);
% plot(theta/pi*180,energypatten);


figure;
thetainstr=30/180*pi;%切片位置
[dthmin,thindex]=min(abs(theta-thetainstr));
hintrs=h(thindex,:);
plot(w/2/pi/1e9,hintrs,'linewidth',3);ylim([0,1.05]);
% set(gca,'XTick',-90:45:90);
set(gca,'YTick',0:0.5:1);
xlabel('Frequency / GHz');ylabel('Normalized AF');
grid on; set (gcf,'Position',[400,400,300,200]);
% plot(w/2/pi,abs(sin( (N+0.5)*d/c*(sin(thetainstr)-sin(theta0)).*w )...
%     ./sin( 0.5*d/c*(sin(thetainstr)-sin(theta0)).*w)  ));%两种方法速度差不多





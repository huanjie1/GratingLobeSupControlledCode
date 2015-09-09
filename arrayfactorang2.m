% arrayfactorang2.m
% pic001,pic002
% 能量积分：频率/时域-角度二维图，偶数阵元


clear
N=16;
c=3e8;
centerfreq=10e9;
centerlambda=c/centerfreq;
d=0.03;

degree=linspace(-90,90,721);
theta=degree/180*pi;

degree0=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=degree0/180*pi;

NN=8001;
f=linspace(-40e9,40e9,NN);
w=2*pi*f;
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;
% NN=8000;% for verify the (-1)^n of ifftshift
% f=linspace(-40e9,40e9,NN);
% w=2*pi*f;
% ts=1/2/40e9;
% t=linspace(-(NN)/2*ts,(NN)/2*ts,NN);


[ww,thth]=meshgrid(w,theta);

h0=sin( (N/2)*d/c*(sin(thth)-sin(theta0)).*ww )...
    ./sin( 0.5*d/c*(sin(thth)-sin(theta0)).*ww );
h0(h0>N)=N;
h0(h0<-N)=-N;
h0(isnan(h0))=N;

% h0=sin( (N/2)*d/c*(sin(thth)-sin(theta0)).*ww )...
%     ./(sin( 0.5*d/c*(sin(thth)-sin(theta0)).*ww )+eps);%not so good

h=abs(h0);
figure;imagesc(w((NN-1)/2:NN)/2/pi,theta/pi*180,h(:,(NN-1)/2:NN));xlabel('Frequency/GHz');
% imagesc(w/2/pi/c*d,theta/pi*180,h);xlabel('d/{\lambda}');
ylabel('\theta');
colorbar;
% set(gcf,'outerposition',get(0,'screensize'));

thetainstr=30/180*pi;%切片位置
[dthmin,thindex]=min(abs(theta-thetainstr));
hintrs=h(thindex,:);
figure;plot(w/2/pi,hintrs);
title('theta section');
% plot(w/2/pi,abs(sin( (N+0.5)*d/c*(sin(thetainstr)-sin(theta0)).*w )...
%     ./sin( 0.5*d/c*(sin(thetainstr)-sin(theta0)).*w)  ));%两种方法速度差不多

finstr=10.0e9;%切片位置
[dfmin,findex]=min(abs(f-finstr));
figure;plot(theta/pi*180,h(:,findex));
title('freq section');
hfint=abs(h(:,findex)).^2;

% sigt1=cos( 2*pi*10e9*t + pi*10e9/(ts*NN)*t.^2);
% sigt1=cos( 2*pi*10e9*t).* exp(-(t/0.1e-9).^2);
sigt1=sigeneratorfor2d( t,  'gaud0', 10e9, 10e9 );
figure;plot(t,sigt1);
title('sig waveform');
sigf1=fft_plot( sigt1, ts, NN, 2 );
sigf1out=(ones(length(degree),1)*sigf1).*h0;
% sigf1out(isnan(sigf1out))=0;
figure;imagesc(w/2/pi,theta/pi*180,abs(sigf1out));
title('sig spetrum');
% sigt1out=ifftshift(ifft(sigf1out,NN,2),2);%%%% wrong!!!!!
%%%% pay attention to the order of ifftshift and ifft
sigt1out=ifft(ifftshift(sigf1out,2),NN,2);
% energypatten2=sum(abs(sigt1out).^2,2);
% figure;plot(theta/pi*180,10*log10(energypatten2),theta/pi*180,10*log10(hfint/max(hfint)*max(energypatten2)));
sigt1outenv=abs(hilbert(sigt1out.').');
figure;imagesc(t,theta/pi*180,sigt1outenv);
title('sig envlope');
% 
% % figure;
% % fstart=8e9;%宽带起止点
% % fend=12e9;
% % [dfsmin,fsindex]=min(abs(f-fstart));
% % [dfemin,feindex]=min(abs(f-fend));
% % energypatten=sum(h(:,fsindex:feindex).^2,2);
% % plot(theta/pi*180,energypatten);





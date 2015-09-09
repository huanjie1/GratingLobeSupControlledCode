% supVSbw.m
% 能量积分与带宽

clear
N=256;
c=3e8;
fc=10e9;
d=0.065;

degree0=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=degree0/180*pi;
% deltadrg=(asin(c/fc/d+sin(0/180*pi))-0/180*pi)/2*180/pi;
deltadrg=90;
degree=linspace(-deltadrg+degree0,deltadrg+degree0,2881);
theta=degree/180*pi;

glposition=asin(sin(theta0)+1*c/fc/d);



NN=8001;
f=linspace(-40e9,40e9,NN);
w=2*pi*f;
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

[ww,thth]=meshgrid(w,theta);
h0=sin( (N/2)*d/c*(sin(thth)-sin(theta0)).*ww )...
    ./sin( 0.5*d/c*(sin(thth)-sin(theta0)).*ww );
h0(h0>N)=N;
h0(h0<-N)=-N;
h0(isnan(h0))=N;
h=abs(h0);
figure;imagesc(w/2/pi,theta/pi*180,h);xlabel('Frequency/GHz');

i0=1:21;
rbw=(i0-1)*0.4/(max(i0)-1);
mainwidth=zeros(1,length(i0));
supratio=zeros(1,length(i0));
figure;
for ii=i0
    ii
    sigt1=sigeneratorfor2d( t,  'lfm', rbw(ii)*fc, fc );
    sigf1=fft_plot( sigt1, ts, NN, 2 );
    sigf1out=(ones(length(degree),1)*sigf1).*h0;
    epn=(sum(abs(sigf1out).^2,2)/NN) / (sum(abs(sigf1).^2)/NN*N^2);
    logpattern=10*log10(epn+eps);
    plot(theta/pi*180,logpattern);hold on;
%     plot(theta/pi*180,epn);hold on;
%     drgmainlobe=theta(epn>0.2*max(epn))*180/pi;
%     mainwidth(ii)=max(drgmainlobe)-min(drgmainlobe);
%     min(drgmainlobe)
%     max(drgmainlobe)
    
    [dthmin,thindex]=min(abs(theta-glposition));
    
    supratio(ii)=max(logpattern)-logpattern(thindex);
end

figure;plot(rbw,supratio);xlabel('相对带宽');ylabel('栅瓣抑制比/dB')
    



% supVSnumINenergy.m
% 能量积分与阵元数

clear
c=3e8;
fc=10e9;
d=6.5e-2;

degree0=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=degree0/180*pi;
% deltadrg=(asin(c/fc/d+sin(0/180*pi))-0/180*pi)/2*180/pi;
deltadrg=90;
degree=linspace(-deltadrg+degree0,deltadrg+degree0,2881);
theta=degree/180*pi;

glposition=asin(sin(theta0)+1*0.03/0.065);

NN=8001;
f=linspace(-40e9,40e9,NN);
w=2*pi*f;
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

[ww,thth]=meshgrid(w,theta);


i0=1:21;
N=round(linspace(10,256,length(i0))/2)*2;
rbw=0.4;
mainwidth=zeros(1,length(i0));
supratio=zeros(1,length(i0));

figure;
for ii=i0
    N(ii)    
    h0=sin( (N(ii)/2)*d/c*(sin(thth)-sin(theta0)).*ww )...
        ./sin( 0.5*d/c*(sin(thth)-sin(theta0)).*ww );
    h0(h0>N(ii))=N(ii);
    h0(h0<-N(ii))=-N(ii);
    h0(isnan(h0))=N(ii);
    h=abs(h0);
    
    sigt1=sigeneratorfor2d( t,  'lfm', rbw*fc, fc );
    sigf1=fft_plot( sigt1, ts, NN, 2 );
    sigf1out=(ones(length(degree),1)*sigf1).*h0;
    epn=(sum(abs(sigf1out).^2,2)/NN) / (sum(abs(sigf1).^2)/NN*N(ii)^2);
    logpattern=10*log10(epn+eps);
    plot(theta/pi*180,logpattern);hold on;
%     plot(theta/pi*180,epn);hold on;
%     drgmainlobe=theta(epn>0.2*max(epn))*180/pi;
%     mainwidth(ii)=max(drgmainlobe)-min(drgmainlobe);
%     min(drgmainlobe)
%     max(drgmainlobe)
    [~,mainwidth(ii)]=threcomp( theta/pi*180, epn,[ 0.5*max(epn), inf], 'on',[-inf inf] );
    [dthmin,thindex]=min(abs(theta-glposition));
    
    supratio(ii)=max(logpattern)-logpattern(thindex);
end
figure;plot(N,mainwidth);
figure;plot(N,10.^(supratio/10));xlabel('Num of elements');ylabel('Suppression Ratio/dB')
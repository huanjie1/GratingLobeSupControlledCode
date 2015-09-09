% GratingVSbw.m
% pic003, pic004
% 能量积分：栅瓣点强度与相对带宽的关系(仿真) 

N=16;
c=3e8;
fc=10e9;
d=9e-2;
m=2;

aimdrg=-0;
deltadrg=(asin(c/fc/d+sin(0/180*pi))-0/180*pi)/2*180/pi;
seedrg=asin(m*c/fc/d+sin(aimdrg/180*pi))*180/pi;
if deltadrg+seedrg>90 || -deltadrg+seedrg<-90
    error('error m');
end
degree=linspace(-deltadrg+seedrg,deltadrg+seedrg,721);
theta=degree/180*pi;

degree0=aimdrg;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=degree0/180*pi;

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
% figure;imagesc(w/2/pi,theta/pi*180,h);xlabel('Frequency/GHz');

i0=1:20;
rbw=(i0-1)*1.9/max(i0);
glpoint=zeros(1,length(i0));
% figure;!!!!!!!!!!!!
for ii=i0
    ii
    sigt1=sigeneratorfor2d( t,  'lfm', rbw(ii)*fc, fc );
    sigf1=fft_plot( sigt1, ts, NN, 2 );
    sigf1out=(ones(length(degree),1)*sigf1).*h0;
    epn=(sum(abs(sigf1out).^2,2)/NN) / (sum(abs(sigf1).^2)/NN*N^2);
%     plot(theta/pi*180,10*log10((epn+1e-9)));hold on;!!!!!!!!!!!!!
%     plot(theta/pi*180,epn);hold on;
%     drgmainlobe=theta(epn>0.2*max(epn))*180/pi;
%     mainwidth(ii)=max(drgmainlobe)-min(drgmainlobe);
%     min(drgmainlobe)
%     max(drgmainlobe)
    glpoint(ii)=interp1(theta/pi*180,epn,seedrg);

end
% figure;plot(rbw,10*log10(glpoint));!!!!!!!!!!!
plot(rbw,10*log10(glpoint));hold on
    



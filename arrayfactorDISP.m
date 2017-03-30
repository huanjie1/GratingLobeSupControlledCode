% arrayfactorDISP.m
% 绘制扫频激光型色散OBFN的二维阵列因子，用于验证解析推导结果

clear

c=299792458;

spacing=0.039;%m

aimdeg=00;

N=4;

degvector=linspace(-90,90,721);
freqvector=linspace(0,40e9,4001);

lambda0=1550e-9;
d0=0;%ps/nm
dstep=100;%ps/nm
b20=-lambda0^2/2/pi/c*d0/1e12*1e9      /2;%1/s/s
b2step=-lambda0^2/2/pi/c*dstep/1e12*1e9    /2;%1/s/s

[wm,theta]=meshgrid(freqvector,degvector);

wm=wm*2*pi;
theta=theta/180*pi;

phi=spacing/c*(sin(theta)-sin(aimdeg/180*pi)).*wm;
H2d=exp(1i*( (b20+(N-1)/2*b2step) * wm.^2 + (N-1)/2*phi )) .*...
    sin( N/2 * (b2step*wm.^2+phi) )./...
    sin( 1/2 * (b2step*wm.^2+phi) )+...
    exp(-1i*( (b20+(N-1)/2*b2step) * wm.^2 - (N-1)/2*phi )) .*...
    sin( N/2 * (b2step*wm.^2-phi) )./...
    sin( 1/2 * (b2step*wm.^2-phi) );

figure;imagesc(freqvector,degvector,abs(H2d));colorbar

% H2dsec=cos( (b20+(N-1)/2*b2step) * wm.^2 + (N-1)/2*phi*0) .*...
%     sin( N/2 * (b2step*wm.^2+phi*0) )./...ccc
%     sin( 1/2 * (b2step*wm.^2+phi*0) );
% figure;plot(freqvector,abs(H2dsec))
%      
    
    



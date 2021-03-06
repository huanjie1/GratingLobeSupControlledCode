function [ xcorrmax, t, xcorrt, xcorrAT0 ] = xcorrTTDarrayr1( aimdgr, seedgr, fig, arraynum, d, bw, mode, fc, fs)
%                                   xcorrTTDarrayr1( 0,       0,      0,     10, 0.06, 10e9, 'lfm', 10e9,  200e9)
%   xcorrTTDarrayr1( 0,0,  0, 10, 0.06, 10e9, 'lfm', 10e9,  200e9)
% xcorrTTDarrayr1.m
% 真时延相关接收的核心计算函数，计算某一角度下的接收时域波形，并分析其最大值
% 被xcorrmaxVSdegree等调用

if nargin<9, fs=200e9;end
if nargin<8, fc=10e9;end
if nargin<7, mode='lfm';end
if nargin<6, bw=10e9;end
if nargin<5, d=0.065;end
if nargin<4, arraynum=16;end
if nargin<3, fig=1;end
if nargin<2, seedgr=0;end
if nargin<1, aimdgr=60;end
c=3e8;

ts=1/fs;
sampnum=65536*2;% be 2^k!!
seeangle=seedgr/180*pi;
aimangle=aimdgr/180*pi;

t=linspace(-ts*sampnum/2,ts*sampnum/2,sampnum);
sig0t=zeros(1,sampnum);
sig0t((sampnum/4+1):(sampnum/4*3))=sigenerator1(ts,sampnum/2,mode,bw,fc);
if 1==fig
    figure;plot(t,sig0t);
    [ sig0f, freq] = fft_plot( sig0t, ts, sampnum, 1);
else
    [ sig0f, freq] = fft_plot( sig0t, ts, sampnum, 2);
end


% delaysteparray=d*sin(aimangle)/c;
% delaystepspace=d*sin(seeangle)/c;
% delaytime=(delaystepspace-delaysteparray)*((-arraynum/2:arraynum/2-1)+0.5);
% [sigaddt,sigaddf]=timedelayftadd( delaytime, sig0f, freq);

hfspace=sin( (arraynum/2)*d/c*(sin(seeangle)-sin(aimangle+1e-14)).*2*pi*freq )...
    ./sin( 0.5*d/c*(sin(seeangle)-sin(aimangle+1e-14)).*2*pi*freq);%+1e-9 to avoid 0/0
hfelepatt=elementpattern(freq,aimangle);
hf=hfspace.*hfelepatt;
if 1==fig
   figure;plot(freq,hf);title('filter-space');
end
sigaddf=sig0f.*hf;
sigaddf(isnan(sigaddf))=0;% 0/0 @ 0
sigaddt=ifft(ifftshift(sigaddf));

if 0.7<fig
    figure;plot(t,sigaddt);
end
xcorrf=sigaddf.*conj(sig0f);
xcorrt=ifftshift(ifft(ifftshift(xcorrf)));
xcorrmax=max(abs(xcorrt));
xcorrAT0=min(abs(xcorrt(abs(t)<14e-12)));


if 0.5<fig
    figure;plot(t,abs(xcorrt));xlim([-2e-9,2e-9]);ylim([0,7e5]);title('xcorr')
end



if 1==fig
    fft_plot( xcorrt, ts, sampnum, 1);
end

end


function [ ringserpresp ] = optiresRINGcalc( ringnum, aimdelay, aimbw, ocenf, fsweep )
% optiresRINGcalc.m
% 级联微环的光谱复响应计算
% 微环数目，目标延时量，目标带宽，光中心频率（非光载波频率），光扫频范围
% 被optiresRINGsys.m调用
% 调用级联微环计算函数optiresRINGserial.m

% a=ringnum+aimdelay+aimbw+ocf+fsweep;
% 
% ringserpresp=exp(-1i*2*pi*aimdelay*fsweep);

if nargin<1
    ringnum=1;
    aimdelay=50e-12;
    aimbw=10e9;
    ocenf=193.4e12+10e9;
    fsweep=193.4e12+(-40e9:0.1e9:40e9);
end

mainfbd=ocenf+(-aimbw/2:fsweep(2)-fsweep(1):aimbw/2);

paramat0=[0.89 ocenf-2e9;...
        0.95  ocenf+2e9];
[ delayerrormean, ringserpresp, lossmean, rvec, paramatok ] = optiresRINGserial( paramat0, mainfbd, aimdelay,1);


end
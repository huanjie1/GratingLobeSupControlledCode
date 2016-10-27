function [ ringserpresp ] = optiresRINGcalc( ringnum, aimdelay, aimbw, ocenf, fsweep )
% optiresRINGcalc.m
% 级联微环的光谱复响应计算
% 微环数目，目标延时量，目标带宽，光中心频率（非光载波频率），光扫频范围
% 被optiresRINGsys.m调用
% 调用级联微环计算函数optiresRINGserial.m

% a=ringnum+aimdelay+aimbw+ocenf+fsweep;
% 
% ringserpresp=exp(-1i*2*pi*aimdelay*fsweep);

if nargin<1
    ringnum=8;
    aimdelay=800e-12;
    aimbw=5e9;
    ocenf=193.4e12+10e9;
    fsweep=193.4e12+(-40e9:0.1e9:40e9);
end

mainfbd=ocenf+(-aimbw/2:fsweep(2)-fsweep(1):aimbw/2);

fc0=linspace(-aimbw/2,aimbw/2,ringnum+2);

% paramat0=[0.89 ocenf-2e9;...
%         0.95  ocenf-0.67e9];
% %         0.9 ocenf+0.67e9;...
% %         0.89 ocenf+2e9];

% only all-pass rings are supported at the preliminary stage

paramat0=[0.5*ones(ceil(ringnum/2),1) ocenf+fc0(2:ceil(ringnum/2)+1).'];

[ delayerrormean, ~, ~, ~, paramatok ] = optiresRINGserial( paramat0, ringnum, mainfbd, aimdelay,1,ocenf);

[ ~, ringserpresp, ~, ~, ~ ] = optiresRINGserial( paramatok, ringnum, fsweep, aimdelay,0);

note1=['number of rings: ' num2str(ringnum) ';   aim delay: ' num2str(aimdelay) '; delay error: ' num2str(delayerrormean)]

end


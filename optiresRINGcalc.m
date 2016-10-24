function [ ringserpresp ] = optiresRINGcalc( ringnum, aimdelay, aimbw, ocf, fsweep )
% optiresRINGcalc.m
% 级联微环的光谱复响应计算，含参数的优化设计
% 微环数目，目标延时量，目标带宽，光中心频率，光扫频范围
% 被optiresRINGsys.m调用
% 调用级联微环计算函数optiresRINGserial.m

a=ringnum+aimdelay+aimbw+ocf+fsweep;

ringserpresp=exp(-1i*2*pi*aimdelay*fsweep);

end
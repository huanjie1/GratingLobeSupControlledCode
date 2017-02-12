function [ outputmat, yrangeout ] = histogramPlot( inputmat, xdata, ydata, yrangein, traceindex, picsw )
% histogramPlot.m
% 用于显示数字余辉——以颜色区分出现的频次（所有曲线的横轴数据须一致）
% 输入：用于显示的矩阵画布，横轴数据，纵轴数据，纵轴范围[ymin,ymax ]，第几条曲线，是否作图
% 输出：矩阵画布
% example:
clear
fs=80e9;
ts=1/fs;
t=0:ts:ts*1000;
sig=sin(2*pi*2.3e9*t);
% plot(t,sig)
plotmat=zeros(600,800);
[ outputmat, yrangeout ]=histogramPlot( plotmat, t, sig, yrange, traceindex, picsw )


if nargin<6 picsw=0; end
if nargin<5 traceindex=1; end

inputmat=zeros()



end


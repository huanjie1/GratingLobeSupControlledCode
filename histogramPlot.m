function [ outputmat, yrangeout ] = histogramPlot( inputmat, xdata, ydata, yrangein, traceindex, picsw )
% histogramPlot.m
% ������ʾ������ԡ�������ɫ���ֳ��ֵ�Ƶ�Σ��������ߵĺ���������һ�£�
% ���룺������ʾ�ľ��󻭲����������ݣ��������ݣ����᷶Χ[ymin,ymax ]���ڼ������ߣ��Ƿ���ͼ
% ��������󻭲�
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


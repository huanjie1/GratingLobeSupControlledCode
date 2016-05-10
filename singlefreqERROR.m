% singlefreqERROR.m
% 用于模拟少数阵元下实测方向图（单频）误差

clear

deg=-90:0.2:90;
theta=deg/180*pi;
theta0=-6/180*pi;%-6
d=4.2/0.75;

% amparray=[1 1 1 1];
% phaarray=-[0 1 2 3]*d*2*pi*sin(theta0);
amparray=[1 0.4];
phaarray=-[0 1 ]*d*2*pi*sin(theta0);


arrayc=amparray.*exp(1i*phaarray);

phaspace=([0 1]*2*pi*d).'*sin(theta);
% phaspace=([0 1 2 3]*2*pi*d).'*sin(theta);

arrayfactor=arrayc*exp(1i*phaspace);
arrayfactornol=arrayfactor/max(abs(arrayfactor));

% figure;plot(deg,20*log10(abs(arrayfactornol)))

elementpattern=exp(-(abs(theta)*5.5))+0.02;
elementpatternnol=elementpattern/max(abs(elementpattern));
% figure;plot(deg,20*log10(abs(elementpatternnol)));%ylim([-30 0])


allp=elementpattern.*arrayfactor;
allpn=allp/max(abs(allp));
plot(deg,20*log10(abs(allpn)));ylim([-50 0]);hold on




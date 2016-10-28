function [ delaytime ] = nearrange( seeangle, aimangle, d, arraynum, c, roomrange, offcenter )
% nearrange.m
% pic011
% 非远场距离计算，偶数阵元；位置关系参考 pic011a
% nearrange( seeangle, aimangle, d, arraynum, c);
if nargin<7, offcenter=0;end
if nargin<6, roomrange=30000;end
if nargin<5, c=3e8;end

% 接收点坐标
xrece=ones(1,arraynum)*0;
yrece=ones(1,arraynum)*roomrange;

delayarray = d*sin(aimangle)/c*((-arraynum/2:arraynum/2-1)+0.5);

rotmat=[cos(seeangle) -sin(seeangle);
        sin(seeangle) cos(seeangle)];
xarray0=((-arraynum/2:arraynum/2-1)+0.5)*d;
yarray0=ones(1,arraynum)*0+offcenter;
coorarray=rotmat*[xarray0;yarray0];
xarray=coorarray(1,:);
yarray=coorarray(2,:);

rangespace=sqrt( (xrece-xarray).^2+(yrece-yarray).^2 );
delayspace=(rangespace-mean(rangespace))/c;

delaytime=delayspace+delayarray;


end


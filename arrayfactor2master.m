% arrayfactor2master.m
% 绘制二维阵列因子
% 调用arrayfactorangFORgeneral.m, arrayfactorangFORgeneralRR.m
clear
N=16;
c=299792458;
diffindex=0:(N-2);
centerfreq=10e9;
centerlambda=c/centerfreq;
spacing0=1.5*centerlambda;

devimode='squ';
% devimode='squ';*7/2

if strcmp(devimode,'squ')
    spacingdia=1*0.0000*diffindex.^2;%deviation from even spacing
    spacings=spacingdia-(min(spacingdia)+max(spacingdia))/2+spacing0;
else
    spacings=randn(1,N-1)*spacing0*0.2+spacing0;
end

% spacings(1)=(spacings(1)+spacings(2))/(1+pi/2);
% spacings(2)=(spacings(1)+spacings(2))/(1+pi/2)*pi/2;

if sum(spacings<0)>0
    error('wrong spacingdia');
end
xposition0=[0 cumsum(spacings)];
xposition=xposition0-(min(xposition0)+max(xposition0))/2;
% figure;stem(xposition,max(spacings)*ones(1,length(xposition)));hold on
% plot(linspace(-max(xposition),max(xposition),length(spacings)),spacings)

aimdegree0=15.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN=8001;
freqaxis=linspace(-40e9,40e9,NN);
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

% dgraxis=linspace(-90,90,37); %for waterfall
dgraxis=linspace(-90,90,721);

% sigt1=cos( 2*pi*10e9*t + pi*10e9/(ts*NN)*t.^2);
% sigt1=cos( 2*pi*10e9*t).* exp(-(t/0.1e-9).^2);
sigt1=[zeros(1,2000) sigeneratorfor2d( t(2001:6001),  'lfm', 4e9, centerfreq ) zeros(1,2000)];
% figure;plot(t,sigt1);
% title('sig waveform');

dgrsection=0;
freqsection=35.1e9;

% arrayfactorangFORgeneral( xposition, freqaxis, dgraxis, t, aimdegree0, dgrsection, freqsection, sigt1, 1 );
% cmpideal=arrayfactorangFORgeneralIDEAL( xposition, freqaxis, dgraxis, t, aimdegree0, centerfreq, dgrsection, freqsection, sigt1, 0, 1 );
cmpue=arrayfactorangFORgeneralRR( xposition, freqaxis, dgraxis, t, aimdegree0, centerfreq, dgrsection, freqsection, sigt1, 9, 1 );
% fom1=cmpFoMcalc( dgraxis, aimdegree0, cmpideal, cmpue, N, 2 )
% fom1=cmpFoMcalc( dgraxis, aimdegree0, cmpideal, cmpue, N, spacing0/centerlambda, 2 )
% fomideal=cmpFoMcalcSIMP( dgraxis, aimdegree0, cmpue, N, spacing0/centerlambda, 2 )
% fomideal=cmpFoMcalcV2( dgraxis, aimdegree0, cmpideal, N, spacing0/centerlambda, 0 )
fomcmpue=cmpFoMcalcV2( dgraxis, aimdegree0, cmpue, N, spacing0/centerlambda, 1 )
% hiap=figure(991);
% lhiap=findall(hiap,'type','line');
% yalliap=get(lhiap,'ydata');
% xiap=get(lhiap,'xdata');
% 
% outputmiap=xiap{1};
% for ii2=1:length(yalliap)
%     outputmiap=[outputmiap; yalliap{ii2}];
% end
% 
% % save('outputmiapHAMMING.mat','outputmiap');
% save('outputmiap.mat','outputmiap');
% 
% hcmp=figure(992);
% lhcmp=findall(hcmp,'type','line');
% yallcmp=get(lhcmp,'ydata');
% xcmp=get(lhcmp,'xdata');
% 
% outputmcmp=xcmp{1};
% for ii2=1:length(yallcmp)
%     outputmcmp=[outputmcmp; yallcmp{ii2}];
% end
% 
% % save('outputmcmpHAMMING.mat','outputmcmp');
% save('outputmcmp.mat','outputmcmp');

% hglenv=figure(993);
% lhglenv=findall(hglenv,'type','line');
% yallhglenv=get(lhglenv,'ydata');
% xhglenv=get(lhglenv,'xdata');
% 
% outputmglenv=xhglenv{1};
% for ii2=1:length(yallhglenv)
%     outputmglenv=[outputmglenv; yallhglenv{ii2}];
% end
% 
% save('outputmglenv.mat','outputmglenv');

% hglenv=figure(999999);
% lhglenv=findall(hglenv,'type','line');
% yallhglenv=get(lhglenv,'ydata');
% xhglenv=get(lhglenv,'xdata');
% plot(xhglenv{1},yallhglenv{1}.*yallhglenv{2}.*yallhglenv{3}.*yallhglenv{4})

% for ii=1:6
%     plot(xhglenv{1},yallhglenv{ii}./max(yallhglenv{ii}),':');hold on
% end

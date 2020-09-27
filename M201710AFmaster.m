% M201710AFmaster.m
% 绘制二维阵列因子，计算方向图
% 调用M201710AFcalc.m, M201710ChanRespCalc.m, 
clear


%% 基础参数
N=50;
c=299792458;
diffindex=0:(N-2);
centerfreq=4e9;
centerlambda=c/centerfreq;
spacing0=0.5*centerlambda;

aimdegree0=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nulldegree0=20;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 一维线阵布局
devimode='squ';
% devimode='squ';*7/2

if strcmp(devimode,'squ')
    spacingdia=1*0.0000*diffindex.^2;%deviation from even spacing
    spacings=spacingdia-(min(spacingdia)+max(spacingdia))/2+spacing0;
else
    spacings=randn(1,N-1)*spacing0*0.4+spacing0;
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


%% 计算参数
NN=8001;
freqaxis=linspace(-40e9,40e9,NN);
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

% dgraxis=linspace(-90,90,37); %for waterfall
dgraxis=linspace(-90,90,721);
% dgraxis=asin(linspace(-1,1,721))*180/pi;

% sigt1=cos( 2*pi*10e9*t + pi*10e9/(ts*NN)*t.^2);
% sigt1=cos( 2*pi*10e9*t).* exp(-(t/0.1e-9).^2);
sigt1=[zeros(1,2000) sigeneratorfor2d( t(2001:6001),  'lfm', 4e9, centerfreq ) zeros(1,2000)];

sigf1=fft_plot( sigt1, ts, NN, 2 );

% figure;plot(t,sigt1);
% title('sig waveform');

dgrsection=aimdegree0;
freqsection=centerfreq;


%% 通道响应

netresponse = M201710ChanRespCalc( xposition, freqaxis, aimdegree0, centerfreq, 1 );

% % draw netresponse
% figure;imagesc(freqaxis((NN-1)/2+1:NN),xposition,angle(netresponse(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('xposition');
% figure;
% for inda=1:N
%     plot(freqaxis((NN-1)/2+1:NN),phase(netresponse(inda,(NN-1)/2+1:NN)));hold on
% end
% figure;imagesc(freqaxis((NN-1)/2+1:NN),xposition,abs(netresponse(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('xposition');
% figure;
% for inda=1:N
%     plot(freqaxis((NN-1)/2+1:NN)/1e9,abs(netresponse(inda,(NN-1)/2+1:NN)));hold on
% end

% % band pass
% ampmask=ones(N,1)*exp(-((abs(freqaxis)-10e9)/5.5e9).^18);
% plot(freqaxis,ampmask)
% netresponse=netresponse.*ampmask;

% % spatial-temporal 2D response
% st2Dresponse=fftshift(fft(conj(netresponse),length(dgraxis/180*pi),1),1); 
% figure;imagesc(freqaxis((NN-1)/2:NN)/1e9,linspace(-pi,pi,length(dgraxis/180*pi)),abs(st2Dresponse(:,(NN-1)/2:NN)));xlabel('Frequency/GHz');






%% 天线响应

antresponsearray  = M201710ANTcalc( xposition, freqaxis, dgraxis, centerfreq, 0);

% % draw antresponsearray
% antresponsearraymmax=max(max(antresponsearray));
% figure;imagesc(freqaxis((NN-1)/2+1:NN),dgraxis,angle(antresponsearray(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('\theta');
% figure;imagesc(freqaxis((NN-1)/2+1:NN)/1e9,dgraxis,abs(antresponsearray(:,(NN-1)/2+1:NN)/antresponsearraymmax));xlabel('Frequency/GHz');ylabel('\theta'); ax=gca;ax.FontSize=18;
% figure;imagesc(freqaxis((NN-1)/2+1:NN)/1e9,sin(dgraxis/180*pi),abs(antresponsearray(:,(NN-1)/2+1:NN)/antresponsearraymmax));xlabel('时间频率/GHz');ylabel('归一化空间频率/m^{-1}'); ax=gca;ax.FontSize=18;





%% 总响应
allresponse=M201710AFcalc( netresponse, antresponsearray, xposition, freqaxis, dgraxis, dgrsection, freqsection, 3 );

% allresponsemmax=max(max(allresponse));
% figure;imagesc(freqaxis((NN-1)/2+1:NN)/1e9,dgraxis,abs(allresponse(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('\theta'); ax=gca;ax.FontSize=18;
% figure;imagesc(freqaxis((NN-1)/2+1:NN)/1e9,sin(dgraxis/180*pi),abs(allresponse(:,(NN-1)/2+1:NN)/allresponsemmax));xlabel('时间频率/GHz');ylabel('归一化空间频率/m^{-1}'); ax=gca;ax.FontSize=18;


%% 相关峰值方向图及品质因数

M201710PATcalcIMG( allresponse, dgraxis, t,0, 0.5, sigt1, sigf1, 1 );

[engue, cmpue]=M201710PATcalc( allresponse, dgraxis, t, sigt1, sigf1, 1 );

figure(126);polarplot(dgraxis/180*pi,20*log10(cmpue),'linewidth',2);hold on
ax=gca;
% ax.RTick=-pmin2+[-60 -40 -20 0];
% ax.RTickLabel={'-60'; '-40'; '-20'; '0'};
ax.ThetaTickLabel={'0'; '30'; '60';'90'; '120'; '150';'±180'; '-150'; '-120';'-90'; '-60'; '-30'};
ax.ThetaZeroLocation='top';
ax.FontSize=13;
ax.LineWidth=1.7;
ax.RLim=[-70 0];
ax.RAxisLocation=-80;
ax.ThetaDir = 'clockwise';

xwidthdraw = profileAnalyzor1d(dgraxis, cmpue,'avg norm','measured',1/sqrt(2),1,1,30);
xwidthdraw2 = profileAnalyzor1d(dgraxis, sqrt(engue),'avg norm','measured',1/sqrt(2),1,1,30);





fomcmpue=cmpFoMcalcV2( dgraxis, aimdegree0, cmpue, N, spacing0/centerlambda, 1 )



%% NULL
netresponse20 = M201710ChanRespCalc( xposition, freqaxis, nulldegree0, centerfreq, 1 );
mask2=M201710AFcalc( netresponse, antresponsearray, xposition, freqaxis, nulldegree0, dgrsection, freqsection, 0 );
nullfull=M201710AFcalc( netresponse20, antresponsearray, xposition, freqaxis, nulldegree0, dgrsection, freqsection, 0 );
hfix=mask2./nullfull;

netresponse2=(ones(N,1)*hfix).*netresponse20; 
netresponsewn=netresponse-netresponse2;
allresponseng = M201710AFcalc( netresponse2, [], xposition, freqaxis, dgraxis, dgrsection, freqsection, 1 );
allresponse2 = M201710AFcalc( netresponsewn, [], xposition, freqaxis, dgraxis, dgrsection, freqsection, 1 );

% % draw netresponse
% figure;imagesc(freqaxis((NN-1)/2+1:NN),xposition,angle(netresponse2(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('xposition');
% figure;
% for inda=1:N
%     plot(freqaxis((NN-1)/2+1:NN),angle(netresponse2(inda,(NN-1)/2+1:NN)));hold on
% end
% figure;imagesc(freqaxis((NN-1)/2+1:NN),xposition,abs(netresponse2(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('xposition');
% figure;
% for inda=1:N
%     plot(freqaxis((NN-1)/2+1:NN),abs(netresponse2(inda,(NN-1)/2+1:NN)));hold on
% end

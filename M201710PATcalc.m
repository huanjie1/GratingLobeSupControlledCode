function [ energypatten2nol, xcorrpatternmaxnol] = M201710PATcalc( allresponse, dgraxis, taxis, sigt1, sigf1, fignum )
% M201710PATcalc
% 相关峰和积分方向图计算
% 系统全响应，角度轴，时间轴，目标角度，角度截面位置，频率截面位置，信号时域波形，信号频谱，作图参数
% 被M201710AFmaster.m调用
% 

NN=length(sigf1);
ts=taxis(2)-taxis(1);

sigf1out=(ones(length(dgraxis),1)*sigf1).*allresponse;
sigt1out=ifft(ifftshift(sigf1out,2),NN,2);


energypatten2=sum(abs(sigt1out).^2,2);
energypatten2nol=energypatten2/max(energypatten2);



%partial xcorr
tlmax=13.5e-9;
maxl=round(tlmax/ts);
tcor=ts*(-maxl:maxl);
xcorr0=zeros(length(dgraxis),2*maxl+1);
for thind=1:length(dgraxis)      
    xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1,maxl);
end

% %full xcorr
% xcorr0=ones(length(theta),NN*2-1);
% tcor=linspace(-taxis(end)+taxis(1),taxis(end)-taxis(1),NN*2-1);
% for thind=1:length(theta)      
%     xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1);
% end

xcorrpattern=abs(hilbert(xcorr0.').');

xcorrpatternmax=max(abs(xcorrpattern).');
xcorrpatternmaxnol=xcorrpatternmax/max(xcorrpatternmax);




if fignum>0  
    figure(991);hold on
    plot(dgraxis,10*log10(energypatten2nol),dgraxis,20*log10(xcorrpatternmaxnol));
    ylim([-40,5]);title('energy pattern');
    hold off   
end

if fignum>1
    sigt1outenv=abs(hilbert(sigt1out.').');
    figure;imagesc(taxis,dgraxis,sigt1outenv);
    title('sig envlope');
    
    figure;imagesc(tcor,dgraxis,abs(xcorrpattern));

end

if fignum>2
    thetainstr2=-17/180*pi;%切片位置
    [dthmin,thindex2]=min(abs(theta-thetainstr2));
    glenv=sigt1outenv(thindex2,:)/16;
    figure(992);hold on
    plot(taxis,glenv);title('envelop at gl');
    hold off

end

  


if fignum>5 
    
    figure;
    [tcorm,angm]=meshgrid(tcor,dgraxis);
%     surf(tcorm,angm,abs(xcorrpattern),'EdgeColor','none');
%     mesh(tcorm,angm,abs(xcorrpattern))
    waterfall(tcorm,angm,abs(xcorrpattern)) % number of angle<=37
    %waterfall(angm.',tcorm.',abs(xcorrpattern.'))  %---another direction~
    title('xcorr pattern');
end

end





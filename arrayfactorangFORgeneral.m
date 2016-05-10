function [ allresponse] = arrayfactorangFORgeneral( xposition, freqaxis, dgraxis, taxis, aimdgr, dgrsection, freqsection, sigt1, fignum )
% arrayfactorangFORgeneral.m
% 频率/角度二维图，一维布阵，任意排布
% 被arrayfactor2master.m调用

c=299792458;

spacings=diff(xposition);
if sum(spacings<0)>0
    error('wrong spacingdia');
end

aimtheta0=aimdgr/180*pi;


w=2*pi*freqaxis;
ts=taxis(2)-taxis(1);
theta=dgraxis/180*pi;
NN=length(w);

arrayresponse0=exp(1i*(-xposition*sin(aimtheta0)/c).'*w);
arrayresponse=arrayresponse0;

allresponse=ones(length(theta),length(w));
for thind=1:length(theta)
    spaceresponse=exp(1i*(xposition*sin(theta(thind))/c).'*w);
    allresponse(thind,:)=sum(arrayresponse.*spaceresponse);
end

if fignum>0
    h=abs(allresponse);
    figure;imagesc(w((NN-1)/2:NN)/2/pi,theta/pi*180,h(:,(NN-1)/2:NN));xlabel('Frequency/GHz');
    % imagesc(w/2/pi/c*d,theta/pi*180,h);xlabel('d/{\lambda}');
    ylabel('\theta');
    colorbar;
    % set(gcf,'outerposition',get(0,'screensize'));
end

if fignum>1
    thetainstr=dgrsection/180*pi;%切片位置
    [dthmin,thindex]=min(abs(theta-thetainstr));
    htheta=allresponse(thindex,:);
    figure(99999);plot(w((NN-1)/2:NN)/2/pi,abs(htheta((NN-1)/2:NN)));hold on
    xlim([0,40e9]);
    title(['theta section @ \theta = ' num2str(dgrsection) 'degree']);
end

if fignum>2
    winstr=freqsection*2*pi;%切片位置
    [dfmin,findex]=min(abs(w-winstr));
    pfreq=allresponse(:,findex);
    figure(999999);plot(theta/pi*180,abs(pfreq));hold on
    title(['freq section @ freq = ' num2str(freqsection/1e9) 'GHz']);
end

if fignum>3
    sigf1=fft_plot( sigt1, ts, NN, 2 );
    sigf1out=(ones(length(dgraxis),1)*sigf1).*allresponse;
    % sigf1out(isnan(sigf1out))=0;
    figure;imagesc(w/2/pi,theta/pi*180,abs(sigf1out));
    title('sig spetrum');
    % sigt1out=ifftshift(ifft(sigf1out,NN,2),2);%%%% wrong!!!!!
    %%%% pay attention to the order of ifftshift and ifft
    sigt1out=ifft(ifftshift(sigf1out,2),NN,2);
    % energypatten2=sum(abs(sigt1out).^2,2);
    % figure;plot(theta/pi*180,10*log10(energypatten2),theta/pi*180,10*log10(hfint/max(hfint)*max(energypatten2)));
    sigt1outenv=abs(hilbert(sigt1out.').');
    figure;imagesc(taxis,theta/pi*180,sigt1outenv);
    title('sig envlope');
end

end





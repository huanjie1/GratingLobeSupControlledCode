function [ xcorrpatternmaxnol, allresponse] = arrayfactorangFORgeneralSW( xposition, freqaxis, dgraxis, taxis, aimdgr, centerfreq, dgrsection, freqsection, sigt1, strunum, fignum, loss, exr )
% arrayfactorangFORgeneralSW.m
% 频率/角度二维图，一维布阵，任意排布，利用离散开关切换网络（含损耗）响应计算
% 一维线阵的坐标向量，频率轴，角度轴，时间轴，目标角度，角度截面位置，频率截面位置，时域信号波形，波束形成网络结构选择，作图参数，损耗，消光比
% 被arrayfactor2master.m调用


c=299792458;

xcorrpatternmaxnol=0;

spacings=diff(xposition);
if sum(spacings<0)>0
    error('wrong spacingdia');
end

aimtheta0=aimdgr/180*pi;


w=2*pi*freqaxis;
ts=taxis(2)-taxis(1);
theta=dgraxis/180*pi;
NN=length(w);
antennanum=length(xposition);

% window1=hamming(antennanum)*ones(1,NN);
window1=rectwin(antennanum)*ones(1,NN);
% window1=[1 1 1 1 0 1 1 1].'*ones(1,NN);





if 9==strunum % switch-based TTD with finite extinction ratio and loss
    dl0=-xposition*sin(aimtheta0)/c;
    dl0=dl0-min(dl0);
    delaybase=0.4*1/centerfreq; %%%
%     exr=25;%dB
    exrl=10^(-exr/20);
    losspermm=loss;%dB/mm
    neff=1.5;
    losspersl=10^(-c/neff*1e-12*1000*losspermm/20);%dB/ps
    swstate=fliplr(dec2bin(round(dl0/delaybase)));
    
    arrayresponse=ones(antennanum,NN);
    
    for anrind=1:size(swstate,1)
        for swstepind=1:size(swstate,2)
            dlstepnow=delaybase*2^(swstepind-1);
            dlpathres=(losspersl^(dlstepnow/1e-12))*exp(1i*dlstepnow*w);
            if strcmp(swstate(anrind,swstepind),'1')                
                arrayresponse(anrind,:)=arrayresponse(anrind,:).*...
                    ((1-exrl)*dlpathres+exrl);
            else
                arrayresponse(anrind,:)=arrayresponse(anrind,:).*...
                    ((1-exrl)+exrl*dlpathres);
            end
        end
    end
%     figure;imagesc(angle(arrayresponse));
%     figure;imagesc(abs(arrayresponse));

end




% figure;imagesc(w((NN-1)/2+1:NN)/2/pi,xposition,angle(arrayresponse(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('xposition');
% figure;
% for inda=1:antennanum
%     plot(w((NN-1)/2+1:NN)/2/pi,phase(arrayresponse(inda,(NN-1)/2+1:NN)));hold on
% end

%band pass
ampmask=ones(antennanum,1)*exp(-((abs(w/2/pi)-10e9)/5.5e9).^18);
% plot(w/2/pi,ampmask)
% arrayresponse=arrayresponse.*ampmask;

% st2Dresponse=fftshift(fft(conj(arrayresponse),length(theta),1),1); 
% figure;imagesc(w((NN-1)/2:NN)/2/pi/1e9,linspace(-pi,pi,length(theta)),abs(st2Dresponse(:,(NN-1)/2:NN)));xlabel('Frequency/GHz');


% % real antenna
% load('E:\课件\MWP\0PROJECTS\201410TTDnew2DTTD\TRANSonAP\col\test3\singleELE\singleELE_amp_ALL.mat');
% load('E:\课件\MWP\0PROJECTS\201410TTDnew2DTTD\TRANSonAP\col\test3\singleELE\singleELE_phase_ALL.mat');
% responseele=ampallele(:,1:10:end).*exp(1i*phaseallele(:,1:10:end));
% responseelefull0=[zeros(size(responseele,1),2700) conj(responseele(:,end:-1:1)) ...
%                 zeros(size(responseele,1),1399) ...
%                 responseele zeros(size(responseele,1),2700)];
% 
% [wm181,tm181]=meshgrid(w,-90:1:90);    
% [wm721,tm721]=meshgrid(w,-90:0.25:90);  
% responseelefull = interp2(wm181,tm181,responseelefull0,wm721,tm721,'nearest'); 
% % figure;imagesc(w/2/pi/1e9,linspace(-pi,pi,length(theta)),abs(responseelefull));


allresponse=ones(length(theta),NN);
for thind=1:length(theta)
    spaceresponse=exp(1i*(xposition*sin(theta(thind))/c).'*w);
%     spaceresponse=exp(1i*(xposition*sin(theta(thind))/c).'*w).*(ones(length(xposition),1)*responseelefull(thind,:));
    
    
%     if abs(theta(thind)-aimtheta0)/pi*180<0.2
%         figure;
%         allr01=arrayresponse.*window1.*spaceresponse;
%         for inda2=2:antennanum
%             plot(w((NN-1)/2+1:NN)/2/pi,phase(allr01(inda2,(NN-1)/2+1:NN)./allr01(inda2-1,(NN-1)/2+1:NN)));hold on
%         end
%     end
    allresponse(thind,:)=sum(arrayresponse.*window1.*spaceresponse);
end

% % alternative method to calculate allresponse,the more loop runs, the more time consumed
% window2=rectwin(antennanum)*ones(1,length(theta));
% allresponse=ones(length(theta),NN);
% for wind=1:length(w)
%     spaceresponse=exp(1i*(xposition.'*w(wind)*sin(theta)/c));
% 
%     allresponse(:,wind)=sum(...
%         arrayresponse(:,wind)*ones(1,length(theta))...
%         .*window2...
%         .*spaceresponse).';
% end


if fignum>0.5
    h=abs(allresponse);
    figure;imagesc(w((NN-1)/2:NN)/2/pi/1E9,theta/pi*180,h(:,(NN-1)/2:NN));xlabel('Frequency/GHz');
%     figure;imagesc(w/2/pi,theta/pi*180,h);xlabel('Frequency/GHz');
    % imagesc(w/2/pi/c*d,theta/pi*180,h);xlabel('d/{\lambda}');
    ylabel('\theta');
    colorbar;
    % set(gcf,'outerposition',get(0,'screensize'));
    
%     sdx1=[w((NN-1)/2) w((NN-1)/2) w(NN) w(NN)]/2/pi;
%     sdy1=[theta(1) theta(end) theta(end) theta(1)]/pi*180;
%     sdxc=[[0.8 0.8]*centerfreq,  [w(NN) w(NN)]/2/pi];
%     sdyc=[0.8 1.2 1.2 0.8]*aimtheta0/pi*180;
%     sdadx=[[1.2 1.2]*centerfreq,  [w(NN) w(NN)]/2/pi];
%     sdady=sdyc;
%     [xsd, ysd] = polybool('subtraction', sdx1, sdy1, sdxc, sdyc); % 边界必须连续    
%     aa1=patch(xsd,ysd,'black','EdgeColor','none','facealpha',0.3);
%     aa2=patch(sdadx,sdady,'black','EdgeColor','none','facealpha',0.3);
    
%     sdx1=[[w((NN-1)/2) w((NN-1)/2)]/2/pi,  [0.8 0.8]*centerfreq];
%     sdy1=[theta(1) theta(end) theta(end) theta(1)]/pi*180;
%     sdx2=[[1.2 1.2]*centerfreq,  [w(NN) w(NN)]/2/pi];
%     sdy2=sdy1;
%     aa1=patch(sdx1,sdy1,'black','EdgeColor','none','facealpha',0.3);
%     aa2=patch(sdx2,sdy2,'black','EdgeColor','none','facealpha',0.3);
    
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
    hfint=abs(allresponse(:,findex)).^2;
end

if fignum>3
    sigf1=fft_plot( sigt1, ts, NN, 2 );
    sigf1out=(ones(length(dgraxis),1)*sigf1).*allresponse;
    % sigf1out(isnan(sigf1out))=0;
%     figure;imagesc(w/2/pi,theta/pi*180,abs(sigf1out));
%     title('sig spetrum');
    % sigt1out=ifftshift(ifft(sigf1out,NN,2),2);%%%% wrong!!!!!
    %%%% pay attention to the order of ifftshift and ifft
    sigt1out=ifft(ifftshift(sigf1out,2),NN,2);
    
    
    energypatten2=sum(abs(sigt1out).^2,2);
%     energypatten2nol=energypatten2/sum(energypatten2)*length(energypatten2);
    energypatten2nol=energypatten2/max(energypatten2);
    figure(991);hold on
%     plot(theta/pi*180,10*log10(energypatten2nol),theta/pi*180,10*log10(hfint/max(hfint)*max(energypatten2nol)));
    plot(theta/pi*180,10*log10(energypatten2nol));
    ylim([-40,5]);title('energy pattern');
    hold off
    
    
    sigt1outenv=abs(hilbert(sigt1out.').');
    figure;imagesc(taxis,theta/pi*180,sigt1outenv);
    title('sig envlope');
    
    
    thetainstr2=-17/180*pi;%切片位置
    [dthmin,thindex2]=min(abs(theta-thetainstr2));
    glenv=sigt1outenv(thindex2,:)/16;
    figure(993);hold on
    plot(taxis,glenv);title('envelop at gl');
    hold off
end

if fignum>0    
%     %full xcorr
%     xcorr0=ones(length(theta),NN*2-1);
%     tcor=linspace(-taxis(end)+taxis(1),taxis(end)-taxis(1),NN*2-1);
%     for thind=1:length(theta)      
%         xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1);
%     end
    
    
    
    if fignum<=3
        sigf1=fft_plot( sigt1, ts, NN, 2 );
        sigf1out=(ones(length(dgraxis),1)*sigf1).*allresponse;
        sigt1out=ifft(ifftshift(sigf1out,2),NN,2);
    end

    %partial xcorr
    tlmax=13.5e-9;
    maxl=round(tlmax/ts);
    tcor=ts*(-maxl:maxl);
    xcorr0=zeros(length(theta),2*maxl+1);
    for thind=1:length(theta)      
        xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1,maxl);
    end
    
    xcorrpattern=abs(hilbert(xcorr0.').');
%     figure;imagesc(tcor,theta/pi*180,abs(xcorrpattern));
    
    xcorrpatternmax=max(abs(xcorrpattern).');
%     xcorrpatternmaxnol=xcorrpatternmax/sum(xcorrpatternmax)*length(xcorrpatternmax);
    xcorrpatternmaxnol=xcorrpatternmax/max(xcorrpatternmax);
%     figure(992);hold on
%     plot(theta/pi*180,20*log10(xcorrpatternmaxnol),...
%         theta/pi*180,10*log10(hfint/max(hfint)*max(xcorrpatternmaxnol)^2));%hfint is energy, 10log; but xcorrpatternmaxnol is volt,^2
%     plot(theta/pi*180,20*log10(xcorrpatternmaxnol));
%     plot(theta/pi*180,(xcorrpatternmaxnol));
%     ylim([0,1.1]);title('xcorr pattern');
%     hold off
    
    
%     figure;
%     [tcorm,angm]=meshgrid(tcor,theta/pi*180);
% %     surf(tcorm,angm,abs(xcorrpattern),'EdgeColor','none');
% %     mesh(tcorm,angm,abs(xcorrpattern))
%     waterfall(tcorm,angm,abs(xcorrpattern)) % number of angle<=37
%     %waterfall(angm.',tcorm.',abs(xcorrpattern.'))  %---another direction~
%     title('xcorr pattern');
end

end





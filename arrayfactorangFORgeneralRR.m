function [ allresponse] = arrayfactorangFORgeneralRR( xposition, freqaxis, dgraxis, taxis, aimdgr, centerfreq, dgrsection, freqsection, sigt1, strunum, fignum )
% arrayfactorangFORgeneralRR.m
% 频率/角度二维图，一维布阵，任意排布，利用各种网络响应计算
% 一维线阵的坐标向量，频率轴，角度轴，时间轴，目标角度，角度截面位置，频率截面位置，时域信号波形，波束形成网络结构选择，作图参数
% 被arrayfactor2master.m调用
% 调用光频响应函数optiresDISP.m, optiresRINGsys.m

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
antennanum=length(xposition);

% window1=hamming(antennanum)*ones(1,NN);
window1=rectwin(antennanum)*ones(1,NN);

if 0==strunum % ideal TTD
    dl0=-xposition*sin(aimtheta0)/c;
    switchmode=0;
    
    if 0==switchmode % continue
        dl=dl0;
    else if 1==switchmode % uniform delay unit %与二叉树式兼容？？？
            delaybase=10e-12;   
            dl=delaybase*round(dl0/delaybase);
        else % delay unit with different steps for different element 
            switchbitnum=4;
            delaybase=(xposition-mean(xposition))*2/c/2^switchbitnum;
            dl=delaybase.*round(dl0./delaybase);
        end            
    end
    
    arrayresponse=exp(1i*(dl).'*w);
    
end

if 0.5<strunum && strunum<3.5 %DISPERSION-BASED
    
    lambda0=1550e-9;
    w0=2*pi*c/lambda0;
    
    if 1==strunum  % multi-wavelength
        d0=300;%ps/nm
        b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
        b1=0;
        b3=0;

        w0=2*pi*c/lambda0;
        wc=(w0+xposition*sin(aimtheta0)/c/b2).';
    end

    if 2==strunum   % multi-DISP, master:wavelength 
        lc=1550.5e-9;
        wc=2*pi*c/lc;
        % bias of d0 would lower the some of the peaks of theta section
        d0=(000+xposition*sin(aimtheta0)/c/(wc-w0)/(-lambda0^2)*2*pi*c*1e12/1e9).';%ps/nm
        b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
        b1=0;
        b3=0;

    end

    if 3==strunum   % multi-DISP, master:DISP
        dstart=-100;
        dend=100;
        d0=((xposition-xposition(1))/(xposition(end)-xposition(1))*(dend-dstart)+dstart).';%ps/nm
        b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
        b1=0;
        b3=0;
        wc=w0 + (xposition(2)-xposition(1)) * sin(aimtheta0) / c / (b2(2)-b2(1));
        lc=2*pi*c/wc
    end
    
    
    
%     % closed form based (begin)    Elapsed time is 0.010651 seconds.
%     modindex=0.1;
%     AA=-1i*besselj(1,modindex);  % sideband #-1
%     BB= 1 *besselj(0,modindex);  % sideband # 0
%     CC= 1i*besselj(1,modindex);  % sideband #+1
%     
%     T= b1 +  b2*(wc-w0) + 1/2*b3*(wc-w0).^2;
%     D=   1/2*b2         + 1/2*b3*(wc-w0)  ;
%     G=                    1/6*b3          ;   
% 
%     wr=w((NN-1)/2+1:NN);
%     beat1cenv=conj(AA)*BB* exp(1i* (-T*wr + D*wr.^2 - G*ones(antennanum,1)*wr.^3) )/(-1i);%complex envelop
%     beat2cenv=conj(BB)*CC* exp(1i* (-T*wr - D*wr.^2 - G*ones(antennanum,1)*wr.^3) )/(-1i);%complex envelop    
%     % closed form based (end)  #####
    
    
    % general solution （begin）    Elapsed time is 0.014591 seconds.
%     wr=w((NN-1)/2+1:NN);    
    modindex=0.1;
    coffn1=-1i*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #-1
    coffc0= 1 *besselj(0,modindex)*ones(antennanum,(NN+1)/2);  % sideband # 0
    coffp1= 1i*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #+1
    
%     Elapsed time is 0.011766 seconds.
    respnall=optiresDISP(b1,b2,b3,w0,wc,w);
    respn1=respnall(:,(NN+1)/2:-1:1); % sideband #-1
    respc0=respnall(:,(NN+1)/2)*ones(1,(NN+1)/2); % sideband # 0
    respp1=respnall(:,(NN+1)/2:NN); % sideband #+1   

% %     Elapsed time is 0.014664 seconds.
%     respn1=optiresDISP(b1,b2,b3,w0,wc,-wr); % sideband #-1
%     respc0=optiresDISP(b1,b2,b3,w0,wc,zeros(1,length(wr))); % sideband # 0
%     respp1=optiresDISP(b1,b2,b3,w0,wc,wr); % sideband #+1
    
    beat1cenv=conj(coffn1.*respn1).*(coffc0.*respc0)/(-1i);%complex envelop
    beat2cenv=conj(coffc0.*respc0).*(coffp1.*respp1)/(-1i);%complex envelop    
    % general solution  (end)  #####

    %%%%%%%%%%%% DSB  SSB
    arrayresponser=1*beat1cenv+1*beat2cenv;
    arrayresponse=[conj(arrayresponser(:,end:-1:2)) arrayresponser];
    
end

if 4==strunum   %ideal phase shifter
    wr=w((NN-1)/2+2:NN);
    ps0=(-xposition*sin(aimtheta0)/c).'*(2*pi*centerfreq*ones(1,length(wr)));
    
    psbitnum=1;
    psbase=2*pi/2^psbitnum;
    ps=psbase*round(ps0/psbase);
    
    arrayresponser=exp(1i*ps);
    arrayresponse=[conj(arrayresponser(:,end:-1:1)) zeros(antennanum,1) arrayresponser];
end

if 5==strunum   %real phase shifter using limited TTD
   delayset=mod(-xposition*sin(aimtheta0)/c*2*pi*centerfreq,2*pi) / (2*pi*centerfreq);
   arrayresponse=exp(1i*(delayset).'*w);
end

if 6==strunum %multi-section phase shifters to approximate TTD
%     wr0=w((NN-1)/2+2:NN);
%     freqbase=1.6;%%%%%%%%%%%%%%%%%%%%%%%%%
%     wrn=freqbase.^round(log(wr0)/log(freqbase));
%     wn=[-wrn(end:-1:1) 0 wrn];
    wr0=w((NN-1)/2+2:NN)/(2*pi*centerfreq); % ensure the ideal phase @ centerfreq     
    relativebw=0.7;
    freqoctave=(2+relativebw)/(2-relativebw);
    suppnum=round(log(10)/log(freqoctave)); %section num for freq below centerfreq    
    freqlogint=round(log(wr0)/log(freqoctave));
    freqlogint(freqlogint<-suppnum)=-suppnum;
    wrn=freqoctave.^(freqlogint);
    wn=[-wrn(end:-1:1) 0 wrn]*(2*pi*centerfreq);
%     plot(w,w,w,wn);
    arrayresponse=exp(1i*(-xposition*sin(aimtheta0)/c).'*wn);
end

if 7==strunum  %with subarray; inter-subarray: ideal TTD; inner-subarray: ideal PS
    subarrayelenum=1;%%%%%%%%%%%%%%%%%%%%%%%
    
    xpositiondl=xposition;
    xpositionps=xposition;
    index7left=1;
    index7leftn=index7left+subarrayelenum-1;
    index7right=length(xposition);    
    index7rightn=index7right-subarrayelenum+1;
    
    while( xposition(index7leftn)<0 && xposition(index7rightn)>0 )
        xpositiondl(index7left:index7leftn)=mean(xposition(index7left:index7leftn));
        xpositionps(index7left:index7leftn)=xposition(index7left:index7leftn)-xpositiondl(index7left:index7leftn);
        xpositiondl(index7rightn:index7right)=mean(xposition(index7rightn:index7right));
        xpositionps(index7rightn:index7right)=xposition(index7rightn:index7right)-xpositiondl(index7rightn:index7right);
        index7left=index7left+subarrayelenum;
        index7leftn=index7leftn+subarrayelenum;
        index7right=index7right-subarrayelenum;        
        index7rightn=index7rightn-subarrayelenum;
    end
    xpositiondl(index7left:index7right)=mean(xposition(index7left:index7right));
    
%     plot(1:antennanum,xposition,1:antennanum,xpositiondl,1:antennanum,xpositionps);xlim([1,16])
    
    arrayresponsedl=exp(1i*(-xpositiondl*sin(aimtheta0)/c).'*w);
    
    wr=w((NN-1)/2+2:NN);
    arrayresponserps=exp(1i*(-xpositionps*sin(aimtheta0)/c).'*(2*pi*centerfreq*ones(1,length(wr))));
    arrayresponseps=[conj(arrayresponserps(:,end:-1:1)) zeros(antennanum,1) arrayresponserps];
    
    arrayresponse=arrayresponsedl.*arrayresponseps;
    
end

if 8==strunum % microring
    
    lc=1550.5e-9;
    wc=2*pi*c/lc;
    
    serialnum=[8 4 2 1];
    aimbw=5e9;
    
    % general solution （begin）    
%     wr=w((NN-1)/2+1:NN);    
    modindex=0.1;
    coffn1=-1i*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #-1
    coffc0= 1 *besselj(0,modindex)*ones(antennanum,(NN+1)/2);  % sideband # 0
    coffp1= 1i*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #+1
    
    respnall=optiresRINGsys(pi/3,aimtheta0,xposition,serialnum,aimbw,centerfreq,wc/2/pi,w/2/pi);
    respn1=respnall(:,(NN+1)/2:-1:1); % sideband #-1
    respc0=respnall(:,(NN+1)/2)*ones(1,(NN+1)/2); % sideband # 0
%     respc0=ones(size(respn1)); % sideband # 0
    respp1=respnall(:,(NN+1)/2:NN); % sideband #+1     
    
    beat1cenv=conj(coffn1.*respn1).*(coffc0.*respc0)/(-1i);%complex envelop
    beat2cenv=conj(coffc0.*respc0).*(coffp1.*respp1)/(-1i);%complex envelop    
    % general solution  (end)  #####

    %%%%%%%%%%%% DSB  SSB
    arrayresponser=0*beat1cenv+1*beat2cenv;
    arrayresponse=[conj(arrayresponser(:,end:-1:2)) arrayresponser];
end

% figure;imagesc(w((NN-1)/2+1:NN)/2/pi,xposition,angle(arrayresponse(:,(NN-1)/2+1:NN)));xlabel('Frequency/GHz');ylabel('xposition');
% figure;
% for inda=1:antennanum
%     plot(w((NN-1)/2+1:NN)/2/pi,phase(arrayresponse(inda,(NN-1)/2+1:NN)));hold on
% end

allresponse=ones(length(theta),NN);
for thind=1:length(theta)
    spaceresponse=exp(1i*(xposition*sin(theta(thind))/c).'*w);
%     if abs(theta(thind)-aimtheta0)/pi*180<0.2
%         figure;
%         allr01=arrayresponse.*window1.*spaceresponse;
%         for inda2=2:antennanum
%             plot(w((NN-1)/2+1:NN)/2/pi,phase(allr01(inda2,(NN-1)/2+1:NN)./allr01(inda2-1,(NN-1)/2+1:NN)));hold on
%         end
%     end
    allresponse(thind,:)=sum(arrayresponse.*window1.*spaceresponse);
end

if fignum>0
    h=abs(allresponse);
    figure;imagesc(w((NN-1)/2:NN)/2/pi,theta/pi*180,h(:,(NN-1)/2:NN));xlabel('Frequency/Hz');
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

if fignum>4    
%     %full xcorr
%     xcorr0=ones(length(theta),NN*2-1);
%     tcor=linspace(-taxis(end)+taxis(1),taxis(end)-taxis(1),NN*2-1);
%     for thind=1:length(theta)      
%         xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1);
%     end
    
    %partial xcorr
    tlmax=3.5e-9;
    maxl=round(tlmax/ts);
    tcor=ts*(-maxl:maxl);
    for thind=1:length(theta)      
        xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1,maxl);
    end
    
    xcorrpattern=abs(hilbert(xcorr0.').');
%     figure;imagesc(tcor,theta/pi*180,abs(xcorrpattern));
    
    xcorrpatternmax=max(abs(xcorrpattern).');
%     xcorrpatternmaxnol=xcorrpatternmax/sum(xcorrpatternmax)*length(xcorrpatternmax);
    xcorrpatternmaxnol=xcorrpatternmax/max(xcorrpatternmax);
    figure(992);hold on
%     plot(theta/pi*180,20*log10(xcorrpatternmaxnol),...
%         theta/pi*180,10*log10(hfint/max(hfint)*max(xcorrpatternmaxnol)^2));%hfint is energy, 10log; but xcorrpatternmaxnol is volt,^2
    plot(theta/pi*180,20*log10(xcorrpatternmaxnol));
    ylim([-40,5]);title('xcorr pattern');
    hold off
    
    
%     figure;
%     [tcorm,angm]=meshgrid(tcor,theta/pi*180);
% %     surf(tcorm,angm,abs(xcorrpattern),'EdgeColor','none');
% %     mesh(tcorm,angm,abs(xcorrpattern))
%     waterfall(tcorm,angm,abs(xcorrpattern)) % number of angle<=37
%     %waterfall(angm.',tcorm.',abs(xcorrpattern.'))  %---another direction~
%     title('xcorr pattern');
end

end





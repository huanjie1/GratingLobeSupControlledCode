function [ allresponse] = arrayfactorangFORgeneralRR( xposition, freqaxis, dgraxis, taxis, aimdgr, centerfreq, dgrsection, freqsection, sigt1, strunum, fignum )
% arrayfactorangFORgeneral.m
% Ƶ��/�Ƕȶ�άͼ��һά���������Ų�������ʵ��OBFN��Ӧ����
% һά���������������Ƶ���ᣬ�Ƕ��ᣬʱ���ᣬĿ��Ƕȣ��ǶȽ���λ�ã�Ƶ�ʽ���λ�ã�ʱ���źŲ��Σ������γ�����ṹѡ����ͼ����
% ��arrayfactor2master.m����

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
    switchmode=2;
    
    if 0==switchmode % continue
        dl=dl0;
    else if 1==switchmode % uniform delay unit
            delaybase=30e-12;   
            dl=delaybase*round(dl0/delaybase);
        else % delay unit with different steps for different element 
            switchbitnum=2;
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
        dstart=-600;
        dend=600;
        d0=((xposition-xposition(1))/(xposition(end)-xposition(1))*(dend-dstart)+dstart).';%ps/nm
        b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
        b1=0;
        b3=0;
        wc=w0 + (xposition(2)-xposition(1)) * sin(aimtheta0) / c / (b2(2)-b2(1));
        lc=2*pi*c/wc
    end
    
    T= b1 +  b2*(wc-w0) + 1/2*b3*(wc-w0).^2;
    D=   1/2*b2         + 1/2*b3*(wc-w0)  ;
    G=                    1/6*b3          ;
    
    modindex=0.1;
    AA=-1i*besselj(1,modindex);
    BB= 1 *besselj(0,modindex);
    CC= 1i*besselj(1,modindex);

    wr=w((NN-1)/2+1:NN);
    beat1cenv=conj(AA)*BB* exp(1i* (-T*wr + D*wr.^2 - G*ones(antennanum,1)*wr.^3) )/(-1i);%complex envelop
    beat2cenv=conj(BB)*CC* exp(1i* (-T*wr - D*wr.^2 - G*ones(antennanum,1)*wr.^3) )/(-1i);%complex envelop

    %%%%%%%%%%%% DSB  SSB
    arrayresponser=1*beat1cenv+1*beat2cenv;
    arrayresponse=[conj(arrayresponser(:,end:-1:2)) arrayresponser];
end

if 4==strunum   %ideal phase shifter
    wr=w((NN-1)/2+2:NN);
    ps0=(-xposition*sin(aimtheta0)/c).'*(2*pi*centerfreq*ones(1,length(wr)));
    
    psbitnum=2;
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
    wr0=w((NN-1)/2+2:NN);
    freqbase=1.1;%%%%%%%%%%%%%%%%%%%%%%%%%
    wrn=freqbase.^round(log(wr0)/log(freqbase));
    wn=[-wrn(end:-1:1) 0 wrn];
%     plot(w,w,w,wn);
    arrayresponse=exp(1i*(-xposition*sin(aimtheta0)/c).'*wn);
end

if 7==strunum  %with subarray; inter-subarray: ideal TTD; inner-subarray: ideal PS
    subarrayelenum=16;%%%%%%%%%%%%%%%%%%%%%%%
    
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

% figure;imagesc(w/2/pi,xposition,angle(arrayresponse));xlabel('Frequency/GHz');ylabel('xposition');

allresponse=ones(length(theta),NN);
for thind=1:length(theta)
    spaceresponse=exp(1i*(xposition*sin(theta(thind))/c).'*w);
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
end

if fignum>1
    thetainstr=dgrsection/180*pi;%��Ƭλ��
    [dthmin,thindex]=min(abs(theta-thetainstr));
    htheta=allresponse(thindex,:);
%     figure(99999);plot(w((NN-1)/2:NN)/2/pi,abs(htheta((NN-1)/2:NN)));hold on
%     xlim([0,40e9]);
%     title(['theta section @ \theta = ' num2str(dgrsection) 'degree']);
end

if fignum>2
    winstr=freqsection*2*pi;%��Ƭλ��
    [dfmin,findex]=min(abs(w-winstr));
    pfreq=allresponse(:,findex);
%     figure(999999);plot(theta/pi*180,abs(pfreq));hold on
%     title(['freq section @ freq = ' num2str(freqsection/1e9) 'GHz']);
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
    
    
    thetainstr2=-17/180*pi;%��Ƭλ��
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





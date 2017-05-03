function [ xcorrpatternmaxnol, allresponse] = arrayfactorangFORgeneralDISP( xposition, freqaxis, dgraxis, taxis, aimdgr, centerfreq, dgrsection, freqsection, sigt1, strunum, fignum, aux )
% arrayfactorangFORgeneralDISP.m
% Ƶ��/�Ƕȶ�άͼ��һά���������Ų�������ɫɢ��Ӧ���㣨����ɫɢ������
% һά���������������Ƶ���ᣬ�Ƕ��ᣬʱ���ᣬĿ��Ƕȣ��ǶȽ���λ�ã�Ƶ�ʽ���λ�ã�ʱ���źŲ��Σ������γ�����ṹѡ����ͼ���������Ӳ���
% ��arrayfactor2masterDISP.m����
% 

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
%         dstart=000;
%         dend=90;
%         d0=((xposition-xposition(1))/(xposition(end)-xposition(1))*(dend-dstart)+dstart).';%ps/nm
        d0=(0:length(xposition)-1).'*aux;
        b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
        b1=0;
        b3=0;
        %FIX!!!        
        wc=w0 + (xposition(2)-xposition(1)) * ...
            ( sin(aimtheta0) - (b2(2)-b2(1))*c*1*pi*centerfreq / (xposition(2)-xposition(1)) )...
            / c / (b2(2)-b2(1));
%         %no fox
%         wc=w0 + (xposition(2)-xposition(1)) * sin(aimtheta0) / c / (b2(2)-b2(1));
        lc=2*pi*c/wc
    end
    
    
    
   
    modindex=0.1;
    alphapolm=117   /180*pi;
    coffn1=-1 *besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #-1
    coffc0= 1 *besselj(0,modindex)*ones(antennanum,(NN+1)/2);  % sideband # 0
    coffp1=-1 *besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #+1
    
%     coffn1=( -cos(alphapolm) -1j*sin(alphapolm) )*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #-1
%     coffc0=(  cos(alphapolm) -1j*sin(alphapolm) )*besselj(0,modindex)*ones(antennanum,(NN+1)/2);  % sideband # 0
%     coffp1=(  cos(alphapolm) +1j*sin(alphapolm) )*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #+1
    
%     Elapsed time is 0.011766 seconds.
    respnall=optiresDISP(b1,b2,b3,w0,wc,w);
    respn1=respnall(:,(NN+1)/2:-1:1); % sideband #-1
    respc0=respnall(:,(NN+1)/2)*ones(1,(NN+1)/2); % sideband # 0
    respp1=respnall(:,(NN+1)/2:NN); % sideband #+1   

% %     Elapsed time is 0.014664 seconds.
%     respn1=optiresDISP(b1,b2,b3,w0,wc,-wr); % sideband #-1
%     respc0=optiresDISP(b1,b2,b3,w0,wc,zeros(1,length(wr))); % sideband # 0
%     respp1=optiresDISP(b1,b2,b3,w0,wc,wr); % sideband #+1
    
    beat1cenv=conj(coffn1.*respn1).*(coffc0.*respc0);%complex envelop
    beat2cenv=conj(coffc0.*respc0).*(coffp1.*respp1);%complex envelop  
    % general solution  (end)  #####

    %%%%%%%%%%%% DSB  SSB
    arrayresponser=0*beat1cenv+1*beat2cenv;
    arrayresponse=[conj(arrayresponser(:,end:-1:2)) abs(arrayresponser(:,1)) arrayresponser(:,2:end)];
    
end



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
end

if fignum>1
    thetainstr=dgrsection/180*pi;%��Ƭλ��
    [dthmin,thindex]=min(abs(theta-thetainstr));
    htheta=allresponse(thindex,:);
    figure(99999);plot(w((NN-1)/2:NN)/2/pi,abs(htheta((NN-1)/2:NN)));hold on
    xlim([0,40e9]);
    title(['theta section @ \theta = ' num2str(dgrsection) 'degree']);
end

if fignum>2
    winstr=freqsection*2*pi;%��Ƭλ��
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
    
    
    thetainstr2=-17/180*pi;%��Ƭλ��
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
    for thind=1:length(theta)      
        xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1,maxl);
    end
    
    xcorrpattern=abs(hilbert(xcorr0.').');
%     figure;imagesc(tcor,theta/pi*180,abs(xcorrpattern));
    
    xcorrpatternmax=max(abs(xcorrpattern).');
%     xcorrpatternmaxnol=xcorrpatternmax/sum(xcorrpatternmax)*length(xcorrpatternmax);
    xcorrpatternmaxnol=xcorrpatternmax/max(xcorrpatternmax);

end

end





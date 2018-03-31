function [ netresponse] = M201710ChanRespCalc( xposition, freqaxis, aimdgr, centerfreq, strunum )
% M201710ChanRespCalc.m
% 计算利用各种网络响应
% 一维线阵的坐标向量，频率轴，角度轴，时间轴，目标角度，信号中心频率，波束形成网络结构选择
% 被M201710AFmaster.m调用
% 调用光频响应函数optiresDISP.m, optiresRINGsys.m

c=299792458;

aimtheta0=aimdgr/180*pi;

w=2*pi*freqaxis;
NN=length(w);
antennanum=length(xposition);


%% back2back
if 0==strunum
    netresponse=ones(antennanum,NN);
end
    


%% ideal TTD
if 1==strunum 
    dl0=-xposition*sin(aimtheta0)/c;
    dl0=dl0-mean(dl0);
    switchmode=0;
    
    if 0==switchmode % continue
        dl=dl0;
    else
        if 1==switchmode % uniform delay unit %与二叉树式兼容？？？
            delaybase=40e-12;   
            dl=delaybase*round(dl0/delaybase);
        else % delay unit with different steps for different element 
            switchbitnum=4;
            delaybase=(xposition-mean(xposition))*2/c/2^switchbitnum;
            dl=delaybase.*round(dl0./delaybase);
        end            
    end
    
%     arrayresponse=exp(1i*(dl).'*w);


    arrayresponse0=exp(1i*(dl).'*w);
    
    ctrlnum=61;
    errlength=(NN-1)/2;
    ctrposi=round(linspace(1,errlength,ctrlnum));
    maxerrdeg=0;% degree
    
    errphaser=zeros(antennanum,errlength);
    for ind11=1:antennanum
        errphaser(ind11,:)=(interp1(ctrposi,rand(1,ctrlnum),1:errlength,'PCHIP')-0.5)*maxerrdeg*2;
    end
%     imagesc(1:errlength,1:antennanum,errphaser);        
%     plot(errphaser.');
    errresponser=exp(1i*errphaser/180*pi);
    errresponse=[conj(errresponser(:,end:-1:1)) zeros(antennanum,1) errresponser];
%     imagesc(w,1:antennanum,angle(errresponse));
    netresponse=arrayresponse0.*errresponse;
%     imagesc(w,1:antennanum,angle(arrayresponse0));
%     figure;imagesc(w,1:antennanum,angle(arrayresponse));
end


%% ideal TTD with null
if 1.5==strunum 
    
    nulldeg=-45;
    
    dl0=-xposition*sin(aimtheta0)/c;
    dln=-xposition*sin(nulldeg/180*pi)/c;
    
    
    netresponse=exp(1i*(dl0).'*w)-1/antennanum*(exp(1i*(dln).'*w)).*...
                        (ones(antennanum,1)*sum(exp(1i*(dln-dl0).'*w)));
    
    
end

%% ideal phase shifter
if 2==strunum   
    wr=w((NN-1)/2+2:NN);
    ps0=(-xposition*sin(aimtheta0)/c).'*(2*pi*centerfreq*ones(1,length(wr)));
    
    psbitnum=100;
    psbase=2*pi/2^psbitnum;
    ps=psbase*round(ps0/psbase);
    
    arrayresponser=exp(1i*ps);
    netresponse=[conj(arrayresponser(:,end:-1:1)) ones(antennanum,1) arrayresponser];
end




%% real phase shifter using limited TTD
if 3==strunum   
   delayset=mod(-xposition*sin(aimtheta0)/c*2*pi*centerfreq,2*pi) / (2*pi*centerfreq);
   netresponse=exp(1i*(delayset).'*w);
end




%% multi-section phase shifters to approximate TTD
if 4==strunum 
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
    netresponse=exp(1i*(-xposition*sin(aimtheta0)/c).'*wn);
end





%% with subarray; inter-subarray: ideal TTD; inner-subarray: ideal PS
if 5==strunum  
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
    
    netresponse=arrayresponsedl.*arrayresponseps;
    
end


%% switch-based TTD with finite extinction ratio and loss
if 6==strunum 
    dl0=-xposition*sin(aimtheta0)/c;
    dl0=dl0-min(dl0);
    delaybase=0.4*1/centerfreq; %%%
    exr=10;%dB
    exrl=10^(-exr/20);
    losspermm=0.1;%dB/mm
    neff=1.5;
    losspersl=10^(-c/neff*1e-12*1000*losspermm/20);%dB/ps
    swstate=fliplr(dec2bin(round(dl0/delaybase)));
    
    netresponse=ones(antennanum,NN);
    
    for anrind=1:size(swstate,1)
        for swstepind=1:size(swstate,2)
            dlstepnow=delaybase*2^(swstepind-1);
            dlpathres=(losspersl^(dlstepnow/1e-12))*exp(1i*dlstepnow*w);
            if strcmp(swstate(anrind,swstepind),'1')                
                netresponse(anrind,:)=netresponse(anrind,:).*...
                    ((1-exrl)*dlpathres+exrl);
            else
                netresponse(anrind,:)=netresponse(anrind,:).*...
                    ((1-exrl)+exrl*dlpathres);
            end
        end
    end
%     figure;imagesc(angle(arrayresponse));
%     figure;imagesc(abs(arrayresponse));

end





%% MWP link
if strunum > 10
    
    lambda0=1550e-9;
    w0=2*pi*c/lambda0;
    
    modindex=0.1;
    alphapolm=-10   /180*pi;
    coffn1= 1 *besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #-1
    coffc0= 1 *besselj(0,modindex)*ones(antennanum,(NN+1)/2);  % sideband # 0
    coffp1= 1 *besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #+1
    
%     coffn1=( -cos(alphapolm) -1j*sin(alphapolm) )*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #-1
%     coffc0=(  cos(alphapolm) -1j*sin(alphapolm) )*besselj(0,modindex)*ones(antennanum,(NN+1)/2);  % sideband # 0
%     coffp1=(  cos(alphapolm) +1j*sin(alphapolm) )*besselj(1,modindex)*ones(antennanum,(NN+1)/2);  % sideband #+1
%     
    
    %% DISPERSION-BASED
    if strunum < 15 
        
        if 11==strunum  % multi-wavelength
            d0=600;%ps/nm
            b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
            b1=0;
            b3=0;

            w0=2*pi*c/lambda0;
            wc=(w0+xposition*sin(aimtheta0)/c/b2).';
        end

        if 12==strunum   % multi-DISP, master:wavelength 
            lc=1550.5e-9;
            wc=2*pi*c/lc;
            % bias of d0 would lower the some of the peaks of theta section
            d0=(000+xposition*sin(aimtheta0)/c/(wc-w0)/(-lambda0^2)*2*pi*c*1e12/1e9).';%ps/nm
            b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
            b1=0;
            b3=0;

        end

        if 13==strunum   % multi-DISP, master:DISP
    %         dstart=000;
    %         dend=90;
    %         d0=((xposition-xposition(1))/(xposition(end)-xposition(1))*(dend-dstart)+dstart).';%ps/nm
            d0=(0:length(xposition)-1).'*30;
            b2=-lambda0^2/2/pi/c*d0/1e12*1e9;%1/s/s
            b1=0;
            b3=0;
    %         %FIX!!!        
    %         wc=w0 + (xposition(2)-xposition(1)) * ...
    %             ( sin(aimtheta0) - (b2(2)-b2(1))*c*1*pi*centerfreq / (xposition(2)-xposition(1)) )...
    %             / c / (b2(2)-b2(1));
            %no fox
            wc=w0 + (xposition(2)-xposition(1)) * sin(aimtheta0) / c / (b2(2)-b2(1));
            lc=2*pi*c/wc
        end
        
        
        respnall=optiresDISP(b1,b2,b3,w0,wc,w);
        
    end
    
    
    %% microring
    if 15==strunum 
        
        serialnum=[8 4 2 1];
        aimbw=5e9;
        
        respnall=optiresRINGsys(pi/3,aimtheta0,xposition,serialnum,aimbw,centerfreq,w0/2/pi,w/2/pi);
    end
    
    
    
    
    
    %%
    respn1=respnall(:,(NN+1)/2:-1:1); % sideband #-1
    respc0=respnall(:,(NN+1)/2)*ones(1,(NN+1)/2); % sideband # 0
    respp1=respnall(:,(NN+1)/2:NN); % sideband #+1   

    
    beat1cenv=conj(coffn1.*respn1).*(coffc0.*respc0);%complex envelop
    beat2cenv=conj(coffc0.*respc0).*(coffp1.*respp1);%complex envelop  

    %%%%%%%%%%%% DSB  SSB
    arrayresponser=1*beat1cenv+1*beat2cenv;
    netresponse=[conj(arrayresponser(:,end:-1:2)) abs(arrayresponser(:,1)) arrayresponser(:,2:end)];%dc must be real!!
    
    
    
end



end





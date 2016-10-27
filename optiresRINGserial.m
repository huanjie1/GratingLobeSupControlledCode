function [ delayerrormean, ringserpresp, lossmean, rvec, paramatok ] = optiresRINGserial( paramat, ringnum, fsweep, aimdelay0, optimflag, fcen)
% optiresRINGserial.m
% 级联微环的响应计算，含参数的优化设计(利用对称性，优化时输入的参数矩阵为一半微环的参数；)
% 参数矩阵（格式见下），微环的实际数目，关心的频率（行向量），目标延时（数值或向量），是否使用优化
%     参数矩阵格式（详细解释见optiresRING.m）：
%     [
%         环1的 trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%         环2的 trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%         环3的 trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%         环4的 trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%             。。。
%     ]
% 级联微环的响应，延迟的均方误差，均方损耗，实际的（等效）环半径
% 被.m调用


% figureon=1;%###########

figureon=1;

if 1==length(aimdelay0)
    aimdelay=ones(1,length(fsweep))*aimdelay0;
else if length(aimdelay0)==length(fsweep)
      aimdelay=aimdelay0;
    else
        error('bad aimdelay0');
    end
end

r0dft=100e-6;
trs2dft=0.8;
ringmodenumdft=1;
tuodft=0.96;
neffdft=1.9375;
yitadft=0.99;

if 1==optimflag
    if nargin<6
        error('bad fcen');
    end
    
    if size(paramat,1)~=ceil(ringnum/2)
        error('bad ringnum');
    end
    
    rownum=size(paramat,1);
else
    if size(paramat,1)~=ringnum
        error('bad ringnum');
    end
    rownum=ringnum;
end

if size(paramat,2)<3
    paramatADD=ones(ringnum,1)*trs2dft;
end

if size(paramat,2)<4
    paramatADD=[paramatADD ones(ringnum,1)*r0dft];
end

if size(paramat,2)<5
    paramatADD=[paramatADD ones(ringnum,1)*[ringmodenumdft tuodft neffdft yitadft]];
end

if 1==optimflag
    bd=paramat(:,2).';
    if 1==length(bd)
        if 1==ringnum
            bd=[bd bd bd];
        else
            bd=[min(fsweep) bd fcen];
        end
    else
        bd=[2*bd(1)-bd(2) bd fcen];
    end
    
    lb=[ones(rownum,1)*0.03, (bd(1:end-2).'-fcen)/1e10];%
    ub=[ones(rownum,1)*tuodft, (bd(3:end).'-fcen)/1e10];%
    
    if size(paramat,2)>2
        lb=[lb ones(rownum,1)*0.1];
        ub=[ub ones(rownum,1)*tuodft];
    end
    
    if size(paramat,2)>3
        lb=[lb ones(rownum,1)*10e-6];
        ub=[ub ones(rownum,1)*1000e-6];
    end
    
    paramat(:,2)=(paramat(:,2)-fcen)/1e10; % move & scaling
    
    er=[];
    foc=[];
    trs11=[];
    opts = optimoptions(@fmincon,'Display','none');
    [paramatok0,delayerrormean,exitflag]=fmincon(@RINGserialnested,paramat,[],[],[],[],lb,ub,[],opts);
    delayerrormean=delayerrormean*1e-12;
    
    if 1==mod(ringnum,2)
        paramatok=[paramatok0; paramatok0(end-1:-1:1,:)];
    else
        paramatok=[paramatok0; paramatok0(end:-1:1,:)];
    end

%     for ind11=1:floor(ringnum/2)
%         paramatok(end+1-ind11,2)=2*fcen-paramatok(end+1-ind11,2);
%     end
    
    paramatok((floor(ringnum/2)+1):end,2)=0-paramatok((floor(ringnum/2)+1):end,2);
    
    paramatok(:,2)=paramatok(:,2)*1e10+fcen; % move & scaling
    
    ringserpresp=0;
    lossmean=0;
    rvec=0;
    
    if 1==figureon
        figure;plot(er)
    end
    
else
    paramat=[paramat paramatADD];

    ringpresp=zeros(size(paramat,1),length(fsweep));
    rvec=zeros(1,size(paramat,1));

    for index1=1:size(paramat,1)
        [ ringpresp(index1,:), rvec(index1) ] = optiresRING( ...
            paramat(index1,1), paramat(index1,2), paramat(index1,3), paramat(index1,4), ...
            paramat(index1,5), paramat(index1,6), paramat(index1,7), fsweep, paramat(index1,8) );
        if 1==figureon
            figure(123);plot(fsweep,abs(ringpresp(index1,:)));hold on
        end
    end

    ringserpresp=prod(ringpresp,1);
    delayres0=-diff(phase(ringserpresp))/(fsweep(2)-fsweep(1))/2/pi;
    delayres=[delayres0 delayres0(end)];

    if 1==figureon
        figure(234);plot(fsweep,delayres);hold on
        figure(123);plot(fsweep,abs(ringserpresp));hold on
    end

    delayerrormean=sqrt(sum((aimdelay-delayres).^2)/length(fsweep));

    lossmean=-20*log10(sqrt(sum(abs(ringserpresp).^2)/length(fsweep)));
    
    paramatok=paramat;
end





% -------------------------------------------------------------------------------------------------------------
    function delayerrormeanNE = RINGserialnested(paramatNE0)
        
        figureonNE=0;
        
        if 1==mod(ringnum,2)
            paramatNE1=[paramatNE0; paramatNE0(end-1:-1:1,:)];
        else
            paramatNE1=[paramatNE0; paramatNE0(end:-1:1,:)];
        end
        
%         for ind11NE=1:floor(ringnum/2)
%             paramatNE1(end+1-ind11NE,2)=2*fcen-paramatNE1(end+1-ind11NE,2);
%         end
        
        paramatNE1((floor(ringnum/2)+1):end,2)=2*0-paramatNE1((floor(ringnum/2)+1):end,2);
        
        paramatNE=[paramatNE1 paramatADD];
        
        paramatNE(:,2)=paramatNE(:,2)*1e10+fcen; % move & scaling

        ringprespNE=zeros(size(paramatNE,1),length(fsweep));
        rvecNE=zeros(1,size(paramat,1));
        
        for index1NE=1:size(paramatNE,1)
            [ ringprespNE(index1NE,:), rvecNE(index1NE) ] = optiresRING( ...
                paramatNE(index1NE,1), paramatNE(index1NE,2), paramatNE(index1NE,3), paramatNE(index1NE,4), ...
                paramatNE(index1NE,5), paramatNE(index1NE,6), paramatNE(index1NE,7), fsweep, paramatNE(index1NE,8) );
            if 1==figureonNE
                figure(123);plot(fsweep,abs(ringprespNE(index1NE,:)));hold on
            end
        end

        ringserprespNE=prod(ringprespNE,1);
        delayres0NE=-diff(phase(ringserprespNE))/(fsweep(2)-fsweep(1))/2/pi;
        delayresNE=[delayres0NE delayres0NE(end)];

        if 1==figureonNE
            figure(234);plot(fsweep,delayresNE);
        end

        delayerrormeanNE=sqrt(sum((aimdelay-delayresNE).^2)/length(fsweep))*1e12;
        
        er=[er delayerrormeanNE];
        foc=[foc paramatNE(:,2)-fcen];
        trs11=[trs11 paramatNE(:,1)];
    end
% -------------------------------------------------------------------------------------------------------------


end


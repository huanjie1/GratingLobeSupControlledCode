function [ ringserpresp, delayerrormean, lossmean, rvec ] = optiresRINGserial( paramat, fsweep, aimdelay0)
% optiresRINGserial.m
% 级联微环的响应计算
% 参数矩阵（格式见下），目标谐振频率，关心的频率（行向量），目标延时（数值或向量）
%     参数矩阵格式（详细解释见optiresRING.m）：
%     [
%         环1的 trs1, faim, r0, modenum, tuo, neff, trs2, yita ;
%         环1的 trs1, faim, r0, modenum, tuo, neff, trs2, yita ;
%         环1的 trs1, faim, r0, modenum, tuo, neff, trs2, yita ;
%         环1的 trs1, faim, r0, modenum, tuo, neff, trs2, yita ;
%             。。。
%     ]
% 级联微环的响应，延迟的均方误差，均方损耗，实际的（等效）环半径
% 被.m调用

% % codes for test (begin)
% clear
% 
% aimdelay0=200e-12;
% fsweep=193.4e12+[-40:0.2:40]*1e9; 
% 
% tobeoptimum=[ ...
%     0.6, 193.4e12-25e9,...
%     0.5, 193.4e12,...
%     0.6, 193.4e12+25e9    ];    
% paramatmain=reshape(tobeoptimum,2,length(tobeoptimum)/2).';    
%     
% r0=100e-6;
% ringmodenum=1;
% tuo=0.96;
% neff=1.9375;
% trs2=0.8;
% yita=0.99;
% paramat0=ones(length(tobeoptimum)/2,1)*[r0, ringmodenum, tuo, neff, trs2, yita ];
%  
% paramat=[paramatmain,paramat0];
% 
% [ ringserpresp, delayerrormean, lossmean, rvec ] = optiresRINGserial( paramat, fsweep, aimdelay0);
% figure(123);plot(fsweep,abs(ringserpresp))
% % codes for test (end)

figureon=1;%###########

% figureon=0;

if 1==length(aimdelay0)
    aimdelay=ones(1,length(fsweep))*aimdelay0;
else if length(aimdelay0)==length(fsweep)
      aimdelay=aimdelay0;
    else
        error('bad aimdelay0');
    end
end

ringpresp=zeros(size(paramat,1),length(fsweep));
rvec=zeros(1,size(paramat,1));

for index1=1:size(paramat,1)
    [ ringpresp(index1,:), rvec(index1) ] = optiresRING( ...
        paramat(index1,1), paramat(index1,2), paramat(index1,3), paramat(index1,4), ...
        paramat(index1,5), paramat(index1,6), fsweep, paramat(index1,7), paramat(index1,8) );
    if 1==figureon
        figure(123);plot(fsweep,abs(ringpresp(index1,:)));hold on
    end
end

ringserpresp=prod(ringpresp,1);
delayres0=-diff(phase(ringserpresp))/(fsweep(2)-fsweep(1))/2/pi;
delayres=[delayres0 delayres0(end)];

if 1==figureon
    figure(234);plot(fsweep,delayres);
end

delayerrormean=sqrt(sum((aimdelay-delayres).^2)/length(fsweep));

lossmean=-20*log10(sqrt(sum(abs(ringserpresp).^2)/length(fsweep)));

end


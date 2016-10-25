function [ delayerrormean, ringserpresp, lossmean, rvec ] = optiresRINGserial( paramat, fsweep, aimdelay0)
% optiresRINGserial.m
% 级联微环的响应计算
% 参数矩阵（格式见下），目标谐振频率，关心的频率（行向量），目标延时（数值或向量）
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

r0dft=100e-6;
trs2dft=0.8;
ringmodenumdft=1;
tuodft=0.96;
neffdft=1.9375;
yitadft=0.99;

if size(paramat,2)<3
    paramat=[paramat ones(size(paramat,1),1)*trs2dft];
end

if size(paramat,2)<4
    paramat=[paramat ones(size(paramat,1),1)*r0dft];
end

if size(paramat,2)<5
    paramat=[paramat ones(size(paramat,1),1)*[ringmodenumdft tuodft neffdft yitadft]];
end


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
    figure(234);plot(fsweep,delayres);
end

delayerrormean=sqrt(sum((aimdelay-delayres).^2)/length(fsweep));

lossmean=-20*log10(sqrt(sum(abs(ringserpresp).^2)/length(fsweep)));

end


function [ delayerrormean, ringserpresp, lossmean, rvec ] = optiresRINGserial( paramat, fsweep, aimdelay0)
% optiresRINGserial.m
% ����΢������Ӧ����
% �������󣨸�ʽ���£���Ŀ��г��Ƶ�ʣ����ĵ�Ƶ�ʣ�����������Ŀ����ʱ����ֵ��������
%     ���������ʽ����ϸ���ͼ�optiresRING.m����
%     [
%         ��1�� trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%         ��2�� trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%         ��3�� trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%         ��4�� trs1, faim, trs2, r0, modenum, tuo, neff, yita ;
%             ������
%     ]
% ����΢������Ӧ���ӳٵľ�����������ģ�ʵ�ʵģ���Ч�����뾶
% ��.m����


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


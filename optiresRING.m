function [ ringpresp ] = optiresRING( tuo, ka, r, neff, wcarr, wrf )
% optiresRING.m
% ����΢����Ӧ�ļ���
% �뾶(m)����Ч�����ʣ����ϵ�������ز���Ƶ�ʣ���Ƶ��Ƶ�ʣ�row vector��
% ��arrayfactorangFORgeneralRR.m����

c=299792458;

if nargin<6
    ka=0.707;
    tuo=0.98;
    r=66e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
    neff=1.9735;    
    wcarr=2*pi*c/1551e-9;
    wrf=2*pi*(-400:0.01:400)*1e9;
end

trs=sqrt(1-ka^2);
delay=2*pi*r/(c/neff);
p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift


ringpresp=(trs-tuo*p.^2)./(1-trs*tuo*p.^2);
% ringpresp=1i*tuo*ka*p./(1-trs*tuo*p.^2);

% plot((wcarr+wrf)/2/pi,abs(ringpresp));hold on
% plot((wcarr+wrf)/2/pi,phase(ringpresp));hold on
plot((wcarr+wrf)/2/pi,[0 -diff(phase(ringpresp))/(wrf(2)-wrf(1))]);



end


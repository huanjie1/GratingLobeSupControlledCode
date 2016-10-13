function [ ringpresp ] = optiresRING( modenum, tuo, trs1, r, neff, wcarr, wrf, trs2, yita )
% optiresRING.m
% ����΢����Ӧ�ļ��㣨���ڲ���ģ�ͣ�
% ģʽѡ��1��ȫͨ��2��������,3������ʽ������ģ�����˴���ϵ�����뾶(m)����Ч�����ʣ����ز���Ƶ�ʣ���Ƶ��Ƶ�ʣ�row vector��������˴���ϵ��
% ��arrayfactorangFORgeneralRR.m����

c=299792458;

if nargin<1
    modenum=2;
end

switch modenum
    case 1  %ȫͨ
        if nargin<8
            trs2=0;
            yita=0;
        end
        if nargin<1
            trs1=0.707;
            tuo=0.7;
            r=66e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
            neff=1.9735;    
            wcarr=2*pi*c/1551e-9;
            wrf=2*pi*(-400:0.01:400)*1e9;
        end
        ka=sqrt(1-trs1^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift

        ringpresp=(trs1-tuo*p.^2)./(1-trs1*tuo*p.^2);
        % ringprespR=1i*tuo*ka*p./(1-trs1*tuo*p.^2);
        
        
    case 2  % ������
        if nargin<9
            yita=0;
        end
        if nargin<1
            trs1=0.9;
            trs2=0.9;
            tuo=1;
            r=108.9e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
            neff=1.9735;    
            wcarr=2*pi*c/1551e-9;
            wrf=2*pi*(-400:0.01:400)*1e9;
        end

        ka1=sqrt(1-trs1^2);
        ka2=sqrt(1-trs1^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift

    %     ringprespT=(trs1-tuo*trs2*p.^2)./(1-tuo*trs1*trs2*p.^2);
        ringpresp=(-sqrt(tuo)*ka1*ka2*p)./(1-tuo*trs1*trs2*p.^2);
    %     ringprespR=1i*tuo*ka1*trs2*p./(1-tuo*trs1*trs2*p.^2);
    
    
    case 3 % ����ʽ
        if nargin<1
            yita=0.99;
            trs1=0.95;
            trs2=0.99;
            tuo=1;
            r=100e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
            neff=1.9735;    
            wcarr=2*pi*c/1551e-9;
            wrf=2*pi*(-400:0.01:400)*1e9;
        end

        ka1=sqrt(1-trs1^2);
        ka2=sqrt(1-trs1^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift
        
        ringpresp=(-sqrt(tuo*yita)*ka1*ka2*p).*(trs2-tuo*trs1*p.^2)...
            ./(1-tuo*trs1*trs2*p.^2).^2;
        
    otherwise
        error('wrong modenum')  ;  
end


% if 1==modenum
%     if nargin<8
%         trs2=0;
%         yita=0;
%     end
%     if nargin<1
%         trs1=0.707;
%         tuo=0.7;
%         r=66e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
%         neff=1.9735;    
%         wcarr=2*pi*c/1551e-9;
%         wrf=2*pi*(-400:0.01:400)*1e9;
%     end
% 
%     ka=sqrt(1-trs1^2);
%     delay=2*pi*r/(c/neff);
%     p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift
% 
%     ringpresp=(trs1-tuo*p.^2)./(1-trs1*tuo*p.^2);
%     % ringprespR=1i*tuo*ka*p./(1-trs1*tuo*p.^2);
% else if 2==modenum
%     if nargin<9
%         yita=0;
%     end
%     if nargin<1
%         trs1=0.9;
%         trs2=0.9;
%         tuo=1;
%         r=108.9e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
%         neff=1.9735;    
%         wcarr=2*pi*c/1551e-9;
%         wrf=2*pi*(-400:0.01:400)*1e9;
%     end
% 
%     ka1=sqrt(1-trs1^2);
%     ka2=sqrt(1-trs1^2);
%     delay=2*pi*r/(c/neff);
%     p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift
% 
% %     ringprespT=(trs1-tuo*trs2*p.^2)./(1-tuo*trs1*trs2*p.^2);
%     ringpresp=(-sqrt(tuo)*ka1*ka2*p)./(1-tuo*trs1*trs2*p.^2);
% %     ringprespR=1i*tuo*ka1*trs2*p./(1-tuo*trs1*trs2*p.^2);
% else if 3==modenum
%     if nargin<1
%         trs1=0.9;
%         trs2=0.9;
%         tuo=1;
%         r=108.9e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
%         neff=1.9735;    
%         wcarr=2*pi*c/1551e-9;
%         wrf=2*pi*(-400:0.01:400)*1e9;
%     end
% 
%     ka1=sqrt(1-trs1^2);
%     ka2=sqrt(1-trs1^2);
%     delay=2*pi*r/(c/neff);
%     p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift
% 
% %     ringprespT=(trs1-tuo*trs2*p.^2)./(1-tuo*trs1*trs2*p.^2);
%     ringpresp=(-sqrt(tuo)*ka1*ka2*p)./(1-tuo*trs1*trs2*p.^2);
% %     ringprespR=1i*tuo*ka1*trs2*p./(1-tuo*trs1*trs2*p.^2);
%     else
%     end
% end




figure(22222);hold on
subplot(3,1,1);plot((wcarr+wrf)/2/pi,abs(ringpresp));title('amp');hold on
subplot(3,1,2);plot((wcarr+wrf)/2/pi,phase(ringpresp));title('phase');hold on
subplot(3,1,3);plot((wcarr+wrf)/2/pi,[0 -diff(phase(ringpresp))/(wrf(2)-wrf(1))]);title('delay')



end


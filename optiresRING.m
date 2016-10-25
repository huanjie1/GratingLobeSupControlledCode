function [ ringpresp, r ] = optiresRING( trs1, faim, trs2, r0, modenum, tuo, neff, fsweep, yita )
% optiresRING.m
% ����΢����Ӧ�ļ��㣨���ڲ���ģ�ͣ�
% ģʽѡ��1��ȫͨ��2��������,3������ʽ������ģ�����˴���ϵ��������˴���ϵ��������ģʽ2��������Ч���뾶(m)��
%            ��Ч�����ʣ�Ŀ��г���Ƶ�ʣ�С��0��Ϊ�޴�Լ������ɨƵƵ�Σ�row vector����������ģ�����ģʽ3��
% ����fsweep�µ���Ӧ������������г�����ȷ����ʵ�ʰ뾶
% ��arrayfactorangFORgeneralRR.m��optiresRINGserial.m����

c=299792458;
figureenable=0;

if nargin<1
    modenum=3;
    figureenable=1;
end

switch modenum
    case 1  %ȫͨ
        if nargin<1
            trs1=0.8;
            tuo=0.96;
            r0=200e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
            neff=1.9735;
            faim=193e12;            
            fsweep=193e12+(-400:0.01:400)*1e9;
        end
        
        if faim<0
            r=r0;
        else
            r= 1/(faim / round( faim / (1/(2*pi*r0/(c/neff))) ) )*(c/neff)/2/pi; % ensure resonate @ faim
        end
        ka=sqrt(1-trs1^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*2*pi*(fsweep)*delay); %half phase-shift

        ringpresp=(trs1-tuo*p.^2)./(1-trs1*tuo*p.^2);
        % ringprespR=1i*tuo*ka*p./(1-trs1*tuo*p.^2);
        
        
    case 2  % ������
        if nargin<1
            trs1=0.9;
            trs2=0.8;
            tuo=0.96;
            r0=120e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
            neff=1.9735; 
            faim=c/1550e-9;                        
            fsweep=193e12+(-400:0.01:400)*1e9;
        end
        if faim<0
            r=r0;
        else
            r= 1/(faim / round( faim / (1/(2*pi*r0/(c/neff))) ) )*(c/neff)/2/pi; % ensure resonate @ faim
        end
        ka1=sqrt(1-trs1^2);
        ka2=sqrt(1-trs2^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*2*pi*(fsweep)*delay); %half phase-shift

    %     ringprespT=(trs1-tuo*trs2*p.^2)./(1-tuo*trs1*trs2*p.^2);
        ringpresp=(-sqrt(tuo)*ka1*ka2*p)./(1-tuo*trs1*trs2*p.^2);
    %     ringprespR=1i*tuo*ka1*trs2*p./(1-tuo*trs1*trs2*p.^2);
    
    
    case 3 % ����ʽ
        if nargin<1
            yita=0.99;
            trs1=0.65;
            trs2=0.95;
            tuo=0.96;
            r0=151e-6; % ��r��һ�������Լ���仯��������Լ��ǡ����neff���ʹ��ĳһƵ�㴦��p����
            neff=1.9735;    
            faim=c/1550e-9;
            fsweep=193e12+(-400:0.01:400)*1e9;
        end
        if faim<0
            r=r0;
        else
            r= 1/(faim / round( faim / (1/(2*pi*r0/(c/neff))) ) )*(c/neff)/2/pi; % ensure resonate @ faim
        end
        ka1=sqrt(1-trs1^2);
        ka2=sqrt(1-trs2^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*2*pi*(fsweep)*delay); %half phase-shift
        
        ringpresp=(-sqrt(tuo*yita)*ka1*ka2*p).*(trs2-tuo*trs1*p.^2)...
            ./(1-tuo*trs1*trs2*p.^2).^2;
        
    otherwise
        error('wrong modenum')  ;  
end



if 1==figureenable
    figure(22222);hold on
    subplot(3,1,1);plot((fsweep),abs(ringpresp));title('amp');hold on
    subplot(3,1,2);plot((fsweep),phase(ringpresp));title('phase');hold on
    subplot(3,1,3);plot((fsweep),[0 -diff(phase(ringpresp))/(fsweep(2)-fsweep(1))/2/pi]);title('delay');
end




end


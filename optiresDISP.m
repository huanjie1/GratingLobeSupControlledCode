function [ dispresp ] = optiresDISP( beta1, beta2, beta3, wref, wcarr, wrf )
% optiresDISP.m
% ������ɫɢ��Ӧ�ļ���
% ̩��չ��ϵ��(column vector suppoorted)���ο����Ƶ�ʣ����ز��Ľ�Ƶ��(column vector suppoorted)����Ƶ��Ƶ��(row vector suppoorted)
% ��arrayfactorangFORgeneralRR.m����

b0=0;

if length(beta2)>length(wcarr)
    antnum=length(beta2);
    wcarra0=wcarr*ones(antnum,1);
    beta2a0=beta2;
else
    antnum=length(wcarr);
    wcarra0=wcarr;
    beta2a0=beta2*ones(antnum,1);
end
beta1a0=beta1*ones(antnum,1);
beta3a0=beta3*ones(antnum,1);

wcarra=wcarra0*ones(1,length(wrf));
beta1a=beta1a0*ones(1,length(wrf));
beta2a=beta2a0*ones(1,length(wrf));
beta3a=beta3a0*ones(1,length(wrf));
wdiav = wcarra+ones(antnum,1)*wrf-wref;

dispresp=exp(-1i*...%refer to the '-' in ɫɢ����ʱ
        (   b0+...
            beta1a.*wdiav.^1+...
        1/2*beta2a.*wdiav.^2+...
        1/6*beta3a.*wdiav.^3)...
        );




end


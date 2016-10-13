function [ ringpresp ] = optiresRING( modenum, tuo, trs1, r, neff, wcarr, wrf, trs2, yita )
% optiresRING.m
% 理想微环响应的计算（基于参量模型）
% 模式选择（1：全通；2：上下载,3：反射式），损耗，输入端传输系数，半径(m)，有效折射率，光载波角频率，射频角频率（row vector），输出端传输系数
% 被arrayfactorangFORgeneralRR.m调用

c=299792458;

if nargin<1
    modenum=2;
end

switch modenum
    case 1  %全通
        if nargin<8
            trs2=0;
            yita=0;
        end
        if nargin<1
            trs1=0.707;
            tuo=0.7;
            r=66e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
            neff=1.9735;    
            wcarr=2*pi*c/1551e-9;
            wrf=2*pi*(-400:0.01:400)*1e9;
        end
        ka=sqrt(1-trs1^2);
        delay=2*pi*r/(c/neff);
        p=exp(-0.5i*(wcarr+wrf)*delay); %half phase-shift

        ringpresp=(trs1-tuo*p.^2)./(1-trs1*tuo*p.^2);
        % ringprespR=1i*tuo*ka*p./(1-trs1*tuo*p.^2);
        
        
    case 2  % 上下载
        if nargin<9
            yita=0;
        end
        if nargin<1
            trs1=0.9;
            trs2=0.9;
            tuo=1;
            r=108.9e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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
    
    
    case 3 % 反射式
        if nargin<1
            yita=0.99;
            trs1=0.95;
            trs2=0.99;
            tuo=1;
            r=100e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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
%         r=66e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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
%         r=108.9e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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
%         r=108.9e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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


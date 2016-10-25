function [ ringpresp, r ] = optiresRING( trs1, faim, trs2, r0, modenum, tuo, neff, fsweep, yita )
% optiresRING.m
% 理想微环响应的计算（基于参量模型）
% 模式选择（1：全通；2：上下载,3：反射式），损耗，输入端传输系数，输出端传输系数（用于模式2），（等效）半径(m)，
%            有效折射率，目标谐振光频率（小于0则为无此约束），扫频频段（row vector），反射损耗（用于模式3）
% 环在fsweep下的响应，由所需中心谐振点所确定的实际半径
% 被arrayfactorangFORgeneralRR.m，optiresRINGserial.m调用

c=299792458;
figureenable=0;

if nargin<1
    modenum=3;
    figureenable=1;
end

switch modenum
    case 1  %全通
        if nargin<1
            trs1=0.8;
            tuo=0.96;
            r0=200e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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
        
        
    case 2  % 上下载
        if nargin<1
            trs1=0.9;
            trs2=0.8;
            tuo=0.96;
            r0=120e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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
    
    
    case 3 % 反射式
        if nargin<1
            yita=0.99;
            trs1=0.65;
            trs2=0.95;
            tuo=0.96;
            r0=151e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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


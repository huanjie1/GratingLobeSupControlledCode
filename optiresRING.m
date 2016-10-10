function [ ringpresp ] = optiresRING( tuo, ka, r, neff, wcarr, wrf )
% optiresRING.m
% 理想微环响应的计算
% 半径(m)，有效折射率，耦合系数，光载波角频率，射频角频率（row vector）
% 被arrayfactorangFORgeneralRR.m调用

c=299792458;

if nargin<6
    ka=0.707;
    tuo=0.98;
    r=66e-6; % 若r以一定的最大公约数变化，则该最大公约数恰能与neff配合使在某一频点处的p不变
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


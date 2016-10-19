function [ ringserpresp ] = optiresRINGserial( fcaim, bwaim, ringnum, fsweep, ringmode )


figureenable=0;

if nargin<1
    fcaim=193.4e12;
    bwaim=25e9;
    fsweep=193.4e12+(-400:0.01:400)*1e9;
    ringnum=3;
    ringmode=1;
    figureenable=1;
end

[ring1,r1]=optiresRING( ringmode, 0.9, 0.6, 100e-6, 1.9735, fcaim-bwaim/2, fsweep, 0, 0 );
[ring2,r2]=optiresRING( ringmode, 0.9, 0.5, 100e-6, 1.9735, fcaim, fsweep, 0, 0 );
[ring3,r3]=optiresRING( ringmode, 0.9, 0.6, 100e-6, 1.9735, fcaim+bwaim/2, fsweep, 0, 0 );
r1
r2
r3
ringserpresp=ring1.*ring2.*ring3;

if 1==figureenable
    figure(33333);hold on
    subplot(3,1,1);plot((fsweep),abs(ringserpresp));title('amp');hold on
    subplot(3,1,2);plot((fsweep),phase(ringserpresp));title('phase');hold on
    subplot(3,1,3);plot((fsweep),[0 -diff(phase(ringserpresp))/(fsweep(2)-fsweep(1))/2/pi]);title('delay');hold on
    figure(33333);hold on
    subplot(3,1,1);plot((fsweep),abs(ring1));title('amp');hold on
    subplot(3,1,2);plot((fsweep),phase(ring1));title('phase');hold on
    subplot(3,1,3);plot((fsweep),[0 -diff(phase(ring1))/(fsweep(2)-fsweep(1))/2/pi]);title('delay');hold on
    figure(33333);hold on
    subplot(3,1,1);plot((fsweep),abs(ring2));title('amp');hold on
    subplot(3,1,2);plot((fsweep),phase(ring2));title('phase');hold on
    subplot(3,1,3);plot((fsweep),[0 -diff(phase(ring2))/(fsweep(2)-fsweep(1))/2/pi]);title('delay');hold on
    figure(33333);hold on
    subplot(3,1,1);plot((fsweep),abs(ring3));title('amp');hold on
    subplot(3,1,2);plot((fsweep),phase(ring3));title('phase');hold on
    subplot(3,1,3);plot((fsweep),[0 -diff(phase(ring3))/(fsweep(2)-fsweep(1))/2/pi]);title('delay');hold on
end

end


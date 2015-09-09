function [ signalt] = sigenerator1( ts, NN, mode, bw, fc )
% sigenerator1.m
% pic010 及多个与“相关”有关的图
% 信号波形产生  
%   
t=linspace(-ts*NN/2,ts*NN/2,NN);%don't forget '/2'
% T0=1e-10;
% signalt=-2*(t).^1.*exp(-((t)/T0).^2);

if strcmp('phco',mode)
    %Phase coding
    % bits=[1 0 1 0 1 1 0 0];
    bits=[1 0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1  0 1 0 ];
    % bits=[1 1 1 1 0 0 0 0 1 0 1 0 0 1 1 0];
    bitrate=bw/2;
    bitwidth=1/bitrate;
    phase=1.5*pi*bits(mod(round(t/bitwidth),length(bits))+1);
    signalt=cos(2*pi*fc*t+phase);
else if strcmp('lfm',mode)
    %LFM
    tw=ts*NN;
    kk=bw/tw;
    signalt=cos(2*pi*fc*t+pi*kk*t.^2);
    else
        signalt=0;
    end
end



end
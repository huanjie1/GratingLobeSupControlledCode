function [ signalt] = sigeneratorfor2d( t,  mode, bw, fc )
%   sigeneratorfor2d.m
%   信号波形产生

switch mode
    case 'phco'
        %Phase coding
        % bits=[1 0 1 0 1 1 0 0];
        bits=[1 0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1  0 1 0 ];
        % bits=[1 1 1 1 0 0 0 0 1 0 1 0 0 1 1 0];
        bitrate=bw/2;
        bitwidth=1/bitrate;
        phase=pi*bits(mod(round(t/bitwidth),length(bits))+1);
        signalt=cos(2*pi*fc*t+phase);
    
    case 'lfm'
        %LFM
        tw=max(t)-min(t);
        kk=bw/tw;
        signalt=cos(2*pi*fc*t+pi*kk*t.^2);
        
    case 'gaud0'
        %GAUSS
        signalt=cos(2*pi*fc*t).* exp(-(t*bw).^2);
        
    case 'gaud1'
        %GAUSS d1 (derivation first, then up-convertion)
%         signalt=-2*pi*fc*sin(2*pi*fc*t).* exp(-(t/(1/bw)).^2)...
%             -2*bw*bw*t.*cos(2*pi*fc*t).* exp(-(t/(1/bw)).^2);
        signalt= -2*bw*bw*t.*cos(2*pi*fc*t).* exp(-(t*bw).^2);
        
        
end




% if strcmp('phco',mode)
%     %Phase coding
%     % bits=[1 0 1 0 1 1 0 0];
%     bits=[1 0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1  0 1 0 ];
%     % bits=[1 1 1 1 0 0 0 0 1 0 1 0 0 1 1 0];
%     bitrate=bw/2;
%     bitwidth=1/bitrate;
%     phase=pi*bits(mod(round(t/bitwidth),length(bits))+1);
%     signalt=cos(2*pi*fc*t+phase);
% else if strcmp('lfm',mode)
%     %LFM
%     tw=max(t)-min(t);
%     kk=bw/tw;
%     signalt=cos(2*pi*fc*t+pi*kk*t.^2);
% 
% else if strcmp('gaud0',mode)
%     %GAUSS
%     signalt=cos( 2*pi*fc*t).* exp(-(t/(1/bw)).^2);
%     end
% else if strcmp('gaud0',mode)
%     %GAUSS
%     signalt=cos( 2*pi*fc*t).* exp(-(t/(1/bw)).^2);
% else if strcmp('gaud0',mode)
%     %GAUSS
%     signalt=cos( 2*pi*fc*t).* exp(-(t/(1/bw)).^2);
%     else
%         signalt=0;
%     end
% end



end
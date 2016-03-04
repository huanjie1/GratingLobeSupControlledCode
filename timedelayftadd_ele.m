function [ sigoutt, sigoutf ] = timedelayftadd_ele( delay, siginf, freq, seeangle)
% 延时叠加函数,计入天线阵元的影响
% global r1;
% r1=randn(1,length(delay))*0.25+1;
r1=ones(1,length(delay));
% r1=[ 0.7063    1.1841    1.1117    0.6808    1.1343    1.1174    0.7039    0.9821   ...
%     1.4005    0.8766   0.6631    1.0397    0.8600    0.9103    1.1025    0.9863    0.9046];
% r2=randn(1,length(delay))*50e-12;
r2=0;
% tic
load('antenna_element.mat');
antennafrposi=interp2(freqdata,angledata,freqanglelinnor,freq(length(freq)/2+1:end),seeangle,'linear',0);
antennafr=[antennafrposi(end:-1:1) antennafrposi];
% toc
% sigoutf=siginf .* (ones(1,length(delay))* exp(-1i*2*pi * (delay.') * freq ));
sigoutf=antennafr .* siginf .* (r1 * exp(-1i*2*pi * ((delay+r2).') * freq ));
sigoutt=ifft(ifftshift(sigoutf));


end


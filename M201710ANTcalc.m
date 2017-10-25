function [ antresponsearray ] = M201710ANTcalc( xposition, freqaxis, dgraxis, centerfreq, antnum)
% M201710ANTcalc.m
% 单元频率相关方向图的计算
% 一维线阵的坐标向量，频率轴，角度轴，天线选择
% 被M201710AFmaster.m调用


c=299792458;



w=2*pi*freqaxis;
theta=dgraxis/180*pi;
NN=length(w);
antennanum=length(xposition);


%% ideal element
if 0==antnum
    antresponsearray=[];
end

%% uniform half-dipole (centerfreq) array -.- -.- -.- .... -.- -.- -.-
if 1==antnum
    [wm,dgrm]=meshgrid(w,theta);
    antresponsearray=( cos(wm/4/centerfreq.*sin(dgrm))- cos(wm/4/centerfreq) ) ./ (cos(dgrm)+1e-9);
end


%% uniform, one measured 
if 2==antnum
    % real antenna
%     load('E:\课件\MWP\0PROJECTS\201410TTDnew2DTTD\TRANSonAP\col\test3\singleELE\singleELE_amp_ALL.mat');
%     load('E:\课件\MWP\0PROJECTS\201410TTDnew2DTTD\TRANSonAP\col\test3\singleELE\singleELE_phase_ALL.mat');
%     responseele=ampallele(:,1:10:end).*exp(1i*phaseallele(:,1:10:end))...
%         .*exp(1i*2*pi*ones(size(ampallele,1),1)*linspace(7e9,13e9,601)*4.865/c);%delay compensation 5m
%     responseelefull0=[zeros(size(responseele,1),2700) conj(responseele(:,end:-1:1)) ...
%                     zeros(size(responseele,1),1399) ...
%                     responseele zeros(size(responseele,1),2700)];
%     save('../antannaResponse38.mat','responseelefull0');

    load('../antannaResponse38.mat');
    [wm0,dgrm0]=meshgrid(w,-90:1:90); 
    [wm,dgrm]=meshgrid(w,dgraxis); 
    antresponsearray = interp2(wm0,dgrm0,responseelefull0,wm,dgrm,'nearest');
%     figure;imagesc(w/2/pi/1e9,linspace(-pi,pi,length(theta)),abs(responseelefull));
    
end




%% uniform, all measured 
if 9==antnum
    antresponsearray=[];
end




end





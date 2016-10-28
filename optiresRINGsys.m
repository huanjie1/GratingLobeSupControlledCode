function [ ringarrayresponse ] = optiresRINGsys( maxscan, aimtheta, xposition, serialnum, aimbw, rfcf, ocf, rfaxis )
% optiresRINGsys.m
% ������ʽ��ͨ��΢�����еĹ��׸���Ӧ�����������������ʱҪ�󵼳�������ʱ��
% �趨���ɨ��ǣ�rad����Ŀ��ǣ�rad������Ԫ���꣨����������������΢�����������������Ӹ����Ͷ��壩��...
%     Ŀ����ʱ������Ƶ����Ƶ�ʣ����ز�Ƶ�ʣ���ƵƵ����
% ��arrayfactorangFORgeneralRR.m����
% ����΢�����㺯��optiresRINGcalc.m

% �˺���Ŀǰֻ�����ڵȼ������
xpref=linspace(min(xposition),max(xposition),length(xposition));
xpd=xposition-xpref;
if sum(xpd.^2)>0.0001*max(abs(xposition)).^2
    error('bad xposition');
end

delaystep=10000e-12; % switch step
c=299792458;

aimdelay0=xposition.'*sin(aimtheta)/c;
delayoffset=-(xposition-min(xposition)).'*sin(abs(maxscan))/c;
aimdelay=aimdelay0-aimdelay0(1)-delayoffset;
if 0==prod(aimdelay>=0)
    error('bad delayoffset');
end

antanum=length(aimdelay0);
levelnum=ceil(log(antanum)/log(2));

if length(serialnum)<levelnum
    serialnum=[serialnum ones(1,levelnum-length(serialnum))];
end

leveldelay=aimdelay(fix(antanum/2)+1)-aimdelay(fix(antanum/2));

levelringmat=zeros(levelnum,length(rfaxis));
for ind1=1:levelnum
    levelringmat(ind1,:)=optiresRINGcalc(...
        serialnum(end-ind1+1), mod(leveldelay*2^(ind1-1),delaystep*2^(ind1-1)), ...
        aimbw, ocf+rfcf, ocf+rfaxis).*...
        exp(-1i*2*pi*...
        floor((leveldelay*2^(ind1-1))/(delaystep*2^(ind1-1)))*delaystep*2^(ind1-1)*(ocf+rfaxis));
end

binbase=2.^((1:levelnum)-1);
elemat=mod(floor((((1:antanum)-1).'*ones(1,levelnum))./(ones(antanum,1)*binbase)),2);
% ringarrayresponse0=elemat*levelringmat;  % WRONG not ADD 
% using log would be a waste of time? NO!!

% % Elapsed time is 0.003-4 seconds.
ringarrayresponse0=exp( elemat * log(levelringmat) );

ringarrayresponse0=ringarrayresponse0./abs(ringarrayresponse0);


% % Elapsed time is 0.003-4 seconds.
% ringarrayresponse0=ones(antanum,length(rfaxis));
% for ind2=1:antanum
%     for ind3=1:levelnum
%         if 1==elemat(ind2,ind3)
%             ringarrayresponse0(ind2,:)=ringarrayresponse0(ind2,:).*levelringmat(ind3,:);
%         end
%     end
% end

ringarrayresponse=ringarrayresponse0.*exp(-1i*2*pi*delayoffset*(ocf+rfaxis));

% figure;imagesc(rfaxis,1:levelnum,abs(levelringmat));title('mag.');
% leveldelay=[zeros(levelnum,1) -diff(angle(levelringmat),1,2)/((rfaxis(2)-rfaxis(1)))/2/pi];
% leveldelay(abs(leveldelay)>1e-8)=0;
% figure;imagesc(rfaxis,1:levelnum,leveldelay);title('delay');
% 
% figure;imagesc(rfaxis,1:antanum,abs(ringarrayresponse0));title('mag.');
% ringarrayresponse0delay=[zeros(antanum,1) -diff(angle(ringarrayresponse0),1,2)/((rfaxis(2)-rfaxis(1)))/2/pi];
% ringarrayresponse0delay(abs(ringarrayresponse0delay)>1e-8)=0;
% figure;imagesc(rfaxis,1:antanum,ringarrayresponse0delay);title('delay');


end



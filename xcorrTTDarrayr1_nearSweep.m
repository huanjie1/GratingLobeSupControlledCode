% xcorrTTDarrayr1_nearSweep.m
% pic011
% ��ؽ��ɨ�裬��Զ����ż����Ԫ

clear

aimdgr=-45;
N=16;
d=0.03;
drgsee0=-90:1:90;
% engp=[];
% for drgsee=drgsee0
%     drgsee
%     [~,~,~,~,eng]=xcorrTTDarrayr1_near( aimdgr, drgsee, 0,N,d);
%     engp=[engp eng];
% end
% 
% plot(drgsee0,engp,'linewidth',2,'color','r');

xcormax=[];
for drgsee=drgsee0
    drgsee
    xcormaxn=xcorrTTDarrayr1_near( aimdgr, drgsee, 0,N,d);
    xcormax=[xcormax xcormaxn];
end
xcormaxnor=xcormax./max(xcormax);
plot(drgsee0,20*log10(xcormaxnor),'linewidth',2,'color','b');
hold on;
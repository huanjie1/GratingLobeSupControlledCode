% engiTTDarrayr1_nearSweep.m
% 
% �������ֽ��ɨ�裬��Զ����ż����Ԫ
% �޸���xcorrTTDarrayr1_nearSweep.m

clear

aimdgr=45.5;
N=16;
d=0.065;
drgsee0=-90:1:90;
% engp=[];
% for drgsee=drgsee0
%     drgsee
%     [~,~,~,~,eng]=xcorrTTDarrayr1_near( aimdgr, drgsee, 0,N,d);
%     engp=[engp eng];
% end
% 
% plot(drgsee0,engp,'linewidth',2,'color','r');

engi=[];
for drgsee=drgsee0
    drgsee
    [~,~,~,~, engin]=xcorrTTDarrayr1_near( aimdgr, drgsee, 0,N,d);
    engi=[engi engin];
end

plot(drgsee0,engi);
hold on;
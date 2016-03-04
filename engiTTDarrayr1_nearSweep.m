% engiTTDarrayr1_nearSweep.m
% 
% 能量积分结果扫描，非远场，偶数阵元
% 修改自xcorrTTDarrayr1_nearSweep.m

clear

aimdgr=45;
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

engi=[];
for drgsee=drgsee0
    drgsee
    [~,~,~,~, engin]=xcorrTTDarrayr1_near( aimdgr, drgsee, 0,N,d);
    engi=[engi engin];
end
enginor=engi/max(engi);
plot(drgsee0,10*log10(enginor));
hold on;
% xcorrTTDarrayr1_nearSweep.m
% pic011
% 相关结果扫描，非远场，偶数阵元

clear

aimdgr=-60;
N=16;
d=0.06;
drgsee0=-90:0.5:90;
bw=4e9;
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
    xcormaxn=xcorrTTDarrayr1_near( aimdgr, drgsee, 0,N,d,bw);
    xcormax=[xcormax xcormaxn];
end
max(xcormax)
xcormaxnor=xcormax./max(xcormax);
plot(drgsee0,20*log10(xcormaxnor),'linewidth',2,'color','b');
hold on;

mout=[drgsee0.' 20*log10(xcormaxnor.')];

fid2=fopen(['xcorr' num2str(aimdgr) '_' num2str(d) '_' num2str(bw/1e9) '.csv'],'w');
fprintf(fid2,'%.6f,%.6f\n',mout.');
fclose(fid2);
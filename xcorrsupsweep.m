% xcorrsupsweep.m
% 从不同维度扫描相关后的栅瓣抑制比变化


% fc=10e9;
% mode='lfm';
% aimdgr=0;
% d=0.065;
% sup=[];
% bw=4e9;
% N0=round(linspace(10,256,21)/2)*2;
% for ii=N0
%     ii
%     sup=[sup xcorrmaxVSdegree( ii, d, aimdgr, bw, mode, fc )];
% end
% figure;plot(N0,sup)



% fc=10e9;
% mode='lfm';
% aimdgr=0;
% d=0.03:0.0025:0.1;
% sup=[];
% N0=16;
% bw=10e9;
% for ii=d
%     ii
%     sup=[sup xcorrmaxVSdegree( N0, ii, aimdgr, bw, mode, fc )];
% end
% figure;plot(d,sup)



% 
% fc=10e9;
% mode='lfm';
% aimdgr=-60:5:60;
% d=0.065;
% sup=[];
% N0=16;
% bw=10e9;
% for ii=aimdgr
%     ii
%     sup=[sup xcorrmaxVSdegree( N0, d, ii, bw, mode, fc )];
% end
% figure;plot(aimdgr,sup)


fc=10e9;
mode='lfm';
aimdgr=0;
d=0.065;
sup=[];
N0=256;
bw=linspace(0,4e9,21);
for ii=bw
    ii
    sup=[sup xcorrmaxVSdegree( N0, d, aimdgr, ii, mode, fc )];
end
figure;plot(bw,sup)


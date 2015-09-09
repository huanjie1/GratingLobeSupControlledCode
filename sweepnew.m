% sweepnew.m
% pic007,,pic008
% 能量积分：主瓣方向性与相对带宽和阵元间距在恒孔径、恒阵元数下的关系（扫描控制）

kk0=1:40;
% d=0.12/max(kk0)*kk0;
d=0.02+0.08/max(kk0)*kk0;
N=round(0.8./d);
epnmaxlogmat=[];
typpatmat=[];

for kk=kk0
%     [epnmaxlog, rbw, typpattern, degree]= dirVSbwfunc( 16,d(kk),'lfm' );
    [epnmaxlog, rbw, typpattern, degree]= dirVSbwfunc( N(kk),d(kk),'lfm' );
    epnmaxlogmat=[epnmaxlogmat;epnmaxlog];
    typpatmat=[typpatmat;typpattern];
end
    
figure;imagesc(rbw,d,epnmaxlogmat);
figure;imagesc(degree,d,typpatmat);
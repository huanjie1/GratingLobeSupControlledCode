% sweepnew.m
% pic007,,pic008
% �������֣����귽��������Դ������Ԫ����ں�׾�������Ԫ���µĹ�ϵ��ɨ����ƣ�

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
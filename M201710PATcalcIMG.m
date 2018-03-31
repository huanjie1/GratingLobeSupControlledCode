function [ xcorrpattern2D, tcor ] = M201710PATcalcIMG( allresponse, dgraxis, taxis, t0ref, rscanmax, sigt1, sigf1, fignum )
% M201710PATcalcIMG
% 面向二维成像的角度-距离方向图（二维滤波器的响应，即点扩散函数）
% 系统全响应，角度轴，信号时间轴，参考时间零点（用于补偿波束形成网络的基础时延），扫描距离最大偏移量，信号时域波形，信号频谱，作图参数
% 被M201710AFmaster.m调用
% 

c=299792458;

NN=length(sigf1);
ts=taxis(2)-taxis(1);

sigf1out=(ones(length(dgraxis),1)*sigf1).*allresponse;
sigt1out=ifft(ifftshift(sigf1out,2),NN,2);





%partial xcorr
tlmax=abs(t0ref)+rscanmax/c;
maxl=round(tlmax/ts);
tcor0=ts*(-maxl:maxl);
xcorr0=zeros(length(dgraxis),2*maxl+1);
for thind=1:length(dgraxis)      
    xcorr0(thind,:)=xcorr(sigt1out(thind,:),sigt1,maxl);
end

[m1,idx1]=min(abs(t0ref-rscanmax/c-tcor0));
[m2,idx2]=min(abs(t0ref+rscanmax/c-tcor0));
tcor=tcor0(idx1:idx2);

xcorrpattern2D=abs(hilbert(xcorr0(:,idx1:idx2).'));

if fignum>0.5
    figure;
%     imagesc(dgraxis, tcor*c, xcorrpattern2D);
    [dgd,rd]=meshgrid(dgraxis,tcor*c);
    surf(dgd,rd,xcorrpattern2D,'EdgeColor','none')
    xlabel('degree');ylabel('Range(m)')
    
end



end





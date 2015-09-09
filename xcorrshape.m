% xcorrshape.m
% 不同观察角度下的相关后波形的动态变化

aimdgr=0;
d=0.065;
filename='1.gif';
N=16;
% colour1=['r','b','k','y','g','c','m'];
% ii=0;
% for seedgr=-6:1:6
%     ii=ii+1;
%     [ ~, t, xcorrt ] = xcorrTTDarrayr1( aimdgr, seedgr, 0, N, d);
%     plot(t,xcorrt,colour1(mod(ii,length(colour1))+1),'LineWidth',2);xlim([-2e-9,2e-9]);ylim([0,7e5]);title('xcorr');
%     hold on;
% end
ii=0;
for seedgr=-60:1:60
    ii=ii+1;
    [ ~, t, xcorrt ] = xcorrTTDarrayr1( aimdgr, seedgr, 0, N, d);
    figure(ii);plot(t,xcorrt,'LineWidth',2);xlim([-2e-9,2e-9]);ylim([0,7e5]);title('xcorr');
    string0={['AIM:' num2str(aimdgr) ' degree'] ;[ ' @ ' num2str(seedgr) ' degree']};
    text(1e-9, 5e5, string0);
    
    
    frame=getframe(ii);
    im=frame2im(frame);%制作gif文件，图像必须是index索引图像
    [I,map]=rgb2ind(im,256);
    if ii==1;
        imwrite(I,map,filename,'gif','Loopcount',inf,...
            'DelayTime',0.1);%loopcount只是在i==1的时候才有用
    else
        imwrite(I,map,filename,'gif','WriteMode','append',...
            'DelayTime',0.1);%layTime用于设置gif文件的播放快慢
    end
    
    if 0==mod(ii,10)
        close all
    end
end
close all

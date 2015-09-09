% directivityVSbw_odd_numFORicocn2015.m
% 能量积分：方向性与带宽的关系，奇数阵元，ICOCN2015论文用
N=4;
c=3e8;
fc=20e9;
d=3e-2;

ND=720;%when 721,notch at grating lobes??
degree=linspace(-90,90,ND);
theta=degree/180*pi;

degree0=30;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=degree0/180*pi;

NN=8001;
f=linspace(-40e9,40e9,NN);
w=2*pi*f;
ts=1/2/40e9;
t=-(NN-1)/2*ts : ts : ts*(NN-1)/2;

[ww,thth]=meshgrid(w,theta);
h0=sin( (N+0.5)*d/c*(sin(thth)-sin(theta0)).*ww )...
    ./sin( 0.5*d/c*(sin(thth)-sin(theta0)).*ww );
h0(h0>2*N+1)=2*N+1;
h0(h0<-2*N-1)=-2*N-1;
h0(isnan(h0))=2*N+1;
h=abs(h0);
figure;imagesc(w/2/pi,theta/pi*180,h);xlabel('Frequency/GHz');

i0=1:10;
rbw=(i0-1)*2/max(i0);
epnmax=zeros(1,length(i0));
eall=zeros(1,length(i0));
figure;
gainout=theta.'/pi*180;
supratio=[];
for ii=i0
    ii
    sigt1=sigeneratorfor2d( t,  'lfm', rbw(ii)*fc, fc );
    eone=sum(sigt1.^2)*4*pi; % 3D energy
    sigf1=fft_plot( sigt1, ts, NN, 2 );
    sigf1out=(ones(length(degree),1)*sigf1).*h0;
%     sigt1out=ifft(ifftshift(sigf1out,2),NN,2);  
%     ep=sum(abs(sigt1out).^2,2);
    ep=sum(abs(sigf1out).^2,2)/NN; % Parseval Th.
%     epn=ep/sum(ep);
    eall(ii)=2*pi*cos(theta)*ep*pi/ND/eone;% /ND: numerical intergal
    epn=4*pi*ep/eall(ii)/eone;
    
%     %%%%%%%%%%%GIF
%     a1=figure(ii+99999);%ii为循环变量
%     set(a1,'color','white');
%     plot(theta/pi*180,10*log10(epn+1e-9));
%     ylim([-20 20]);%选取绘图区的并集
%     text(-90,18,...
%         ['rbw = ',  num2str(rbw(ii))],...
%         'VerticalAlignment','bottom',...
%         'HorizontalAlignment','left');%扫描参数的说明，位置自行设定
%     frame=getframe(ii+99999);
%     im=frame2im(frame);%制作gif文件，图像必须是index索引图像
%     close (ii+99999)
%     [I,map]=rgb2ind(im,256);
%     if ii==1;
%         imwrite(I,map,'filename.gif','gif','Loopcount',inf,...
%             'DelayTime',0.1);
%     else
%         imwrite(I,map,'filename.gif','gif','WriteMode','append',...
%             'DelayTime',0.15);%layTime用于设置gif文件的播放快慢
%     end
%     %%%%%%%GIF
    gaindb=10*log10(epn+1e-9);
    supratio=[supratio max(epn) / nearesty(theta/pi*180,epn,0)];
    gainout=[gainout gaindb];
    plot(theta/pi*180,gaindb);hold on;
    epnmax(ii)=max(epn);
%     dlt=sum(ep)-N*sum(abs(sigf1).^2)
end
fid1=fopen(['E:\课件\MWP\0PROJECTS\201410TTDnew2DTTD\TRANS ON AP\icocn2015\', 'gainout','.csv'],'w');
fprintf(fid1,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,\n',gainout.');
fclose(fid1);

figure;plot(rbw,10*log10(epnmax));
figure;plot(rbw,eall);  
figure;plot(rbw,supratio);

supratio=[rbw.' supratio.'];
fid2=fopen(['E:\课件\MWP\0PROJECTS\201410TTDnew2DTTD\TRANS ON AP\icocn2015\', 'supratio','.csv'],'w');
fprintf(fid2,'%.4f,%.4f\n',supratio.');
fclose(fid2);


% 计算延时误差下的方向图
% 部分代码来自arrayfactorangFORgeneralRR.m
clear
c=299792458;

varstd=2e-12;

dgraxis=-90:0.1:90;
aimdgr=0;
freqaxis=10e9;
d=0.015;
antennanum=8;
xposition=0:d:d*(antennanum-1);
expnum=1000;
rmat=randn(expnum,antennanum)*varstd;
figure;h1=histogram(rmat);
h1.BinLimits= [-30e-12 30e-12];
title(['\sigma = ' num2str(varstd) 'ps'])

ax1=gca;
ax1.FontSize=6;
af1=gcf;
af1.PaperUnits = 'centimeters';
af1.PaperPosition = [1 0 4 3];
af1.PaperPositionMode = 'manual';
print([num2str(varstd) 'err.png'],'-dpng')

aimtheta0=aimdgr/180*pi;
w=2*pi*freqaxis;
theta=dgraxis/180*pi;
NN=length(w);
window1=rectwin(antennanum)*ones(1,length(theta));

dl0=-xposition*sin(aimtheta0)/c;


figure;
peakrem=zeros(1,expnum);
for rind=1:expnum
    
    dl=dl0+rmat(rind,:);

    arrayresponse=exp(1i*(dl).'*w);


    allresponse=ones(length(theta),NN);
    for wind=1:length(w)
        spaceresponse=exp(1i*(xposition.'*w(wind)*sin(theta)/c));

        allresponse(:,wind)=sum(...
            arrayresponse(:,wind)*ones(1,length(theta))...
            .*window1...
            .*spaceresponse).';
    end

    plot(dgraxis,20*log10(abs(allresponse)/antennanum));hold on
    ylim([-60,5]);xlim([-90,90]);
    grid on
    [~,maxind]=max(abs(allresponse));
    peakrem(rind)=dgraxis(maxind);
    
end

af2=gcf;
af2.PaperUnits = 'centimeters';
af2.PaperPosition = [1 0 10 8];
af2.PaperPositionMode = 'manual';
print([num2str(varstd) 'pattern.png'],'-dpng')

varstdpatternpeak=sqrt(sum((peakrem-mean(peakrem)).^2)/expnum);
figure;h2=histogram(peakrem);
h2.BinLimits= [-3 3];
title(['\sigma = ' num2str(varstdpatternpeak) 'deg'])
ax3=gca;
ax3.FontSize=8;
af3=gcf;
af3.PaperUnits = 'centimeters';
af3.PaperPosition = [1 0 7 5];
af3.PaperPositionMode = 'manual';
print([num2str(varstd) 'patternerr.png'],'-dpng')







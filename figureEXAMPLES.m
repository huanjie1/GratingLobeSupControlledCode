% latex:
% labels = {'$f_1(x)$', '$f_2(x)$', '$f_3(x)$', '$f_4(x)$'};
% legend(labels, 'interpreter', 'latex', 'FontSize', 12)





% % plot basis
% 
% fid1=fopen('K2.5data.csv');
% data1=textscan(fid1,'%f%f','headerlines',0,'delimiter', ',');
% fclose(fid1);
% N=data1{1}.';
% BW2P5=data1{2}.';
% 
% lear
% load('outputmiap.mat');
% load('outputmcmp.mat');
% drg=outputmiap(1,:);
% iap=outputmiap(end:-1:2,:);%end:-1:2
% cmp=outputmcmp(end:-1:2,:);
% 
% [dmin,dindex]=min(abs(drg-(-17)));
% gliapat=outputmiap(end:-1:2,dindex);
% glcmpat=outputmcmp(end:-1:2,dindex);
% 
% freqbw=0:0.5:4;
% 
% plot(freqbw,gliapat,'--*','LineWidth',2);hold on;
% plot(freqbw,glcmpat,'--o','LineWidth',2);hold off
% grid on;xlim([0 4.1]);ylim([-17 0])
% 
% figure;imagesc(freq(1001:5001)/1e9,angleaxis,abs(responsearrcalcfig));caxis([0,1]);
% figurewidth=9;
% figureheight=7;
% af=gcf;
% af.Units='centimeters';
% af.Position(3)=figurewidth;
% af.Position(4)=figureheight;
% 
% ax=gca;
% ax.FontSize=18;
% ax.XTick=[8 9 10 11 12];
% % ax.YTick=[ -15 -10 -5 0];
% % ax.LineWidth=1.5;
% % ax.OuterPosition=[0 0.07 1 0.93];
% xlabel({'Frequency (GHz)';'(a)'},'FontSize',18,'FontWeight','bold','FontName','Arial');
% ylabel({'{\theta} (degree)'},'FontSize',18,'FontWeight','bold','FontName','Arial');% 1y
% 
% l1=legend('Direct','Indirect','Location','northwest');
% set(l1,'FontSize',13)
% l1.Position(1)=l1.Position(1)+0.09;
% l1.Position(2)=l1.Position(2);
% 
% saveas(gca,['4F.fig']);
% af.PaperUnits = 'centimeters';
% af.PaperPosition = [1 0 figurewidth figureheight];
% af.PaperPositionMode = 'manual';
% print(['4F.png'],'-dpng','-r400')




% % plotyy
% 
% figure;
% [hAx,hLine1,hLine2] = plotyy([freq',freq',freq',freq'],...
%     [amplitude1.',amplitude2.',amplitude3.',amplitude4.'],...
%     [freq',freq',freq',freq'],...
%     [dpha1',dpha2',dpha3',dpha4']);
% 
% hLine2(1).LineStyle = '--';
% hLine2(2).LineStyle = '--';
% hLine2(3).LineStyle = '--';
% hLine2(4).LineStyle = '--';
% hLine1(1).LineWidth = 2;
% hLine1(2).LineWidth = 2;
% hLine1(3).LineWidth = 2;
% hLine1(4).LineWidth = 2;
% hLine2(1).Color=hLine1(1).Color;
% hLine2(2).Color=hLine1(2).Color;
% hLine2(3).Color=hLine1(3).Color;
% hLine2(4).Color=hLine1(4).Color;
% 
% xlabel('freq/Hz')
% set(hAx(1),'YTick',linspace(-47,-27,5),'YLim',[-47 -27],'FontSize',10);
% ylabel(hAx(1),{'Amplitude / dB'});
% set(hAx(2),'YTick',linspace(-90,30,5),'YLim',[-90 30],'FontSize',10);
% ylabel(hAx(2),{'Phase Deviation / deg'});
% title('Response of the Channels')
% grid on
% annotation('ellipse',[.45 .19 .03 .23],'LineStyle','-');
% annotation('arrow',[.465 .4] ,[.42 .42],'LineStyle','-');
% annotation('ellipse',[.6 .63 .03 .2],'LineStyle','-');
% annotation('arrow',[.615 .68] ,[.63 .63],'LineStyle','-');
% 
% l1=legend('channel1','channel2','channel3','channel4','Location','southeast');
% set(l1,'FontSize',15,'box','off')
% rectx=min(hAx(1).XTick)+(max(hAx(1).XTick)-min(hAx(1).XTick))*...
%     ([
%         l1.Position(1) ...
%         l1.Position(1)+l1.Position(3) ...
%         l1.Position(1)+l1.Position(3) ...
%         l1.Position(1)...
%     ]-1*hAx(1).Position(1))/hAx(1).Position(3);
% recty=min(hAx(1).YTick)+(max(hAx(1).YTick)-min(hAx(1).YTick))*...
%     ([
%         l1.Position(2) ...
%         l1.Position(2) ...
%         l1.Position(2)+l1.Position(4) ...
%         l1.Position(2)+l1.Position(4)...
%     ]-1*hAx(1).Position(2))/hAx(1).Position(4);
% delete(l1)
% patch(rectx,recty,'white','facealpha',0.6,'LineStyle','-');%can only cover plots on hAx(1)
% l1=legend('channel1','channel2','channel3','channel4','Location','southeast');
% set(l1,'FontSize',15,'box','off')









% % 3D figure
% 
% for ii1=1:length(freqbw)
%     henv=plot3(t*1e9,freqbw(ii1)*ones(1,length(t)),sigenv(ii1,:),'linewidth',2);hold on
% end
% grid on
% hold off
% 
% view(-23,60) 
% zlim([-0.02,1.02]);xlim([-50.1 50.1]);
% 
% 
% [angm,bwm]=meshgrid(drg,freqbw);
% figure;hiap=waterfall(angm,bwm,iap);view(-25,33)
% set(hiap,'linewidth',1.5)
% zlim([-30,0.2]);xlim([-89.9 89.9]);colormap('jet');caxis([-50,0])
% 
% figurewidth=9;
% figureheight=7;
% af=gcf;
% af.Units='centimeters';
% af.Position(3)=figurewidth;
% af.Position(4)=figureheight;
% 
% ax=gca;
% ax.FontSize=12;
% ax.LineWidth=1.5;
% ax.XTick=[-60 0 60];
% ax.YTick=[0 1 2 3 4];
% ax.OuterPosition=[0 0.07 1 0.93];
% xl=xlabel({'Observation Angle (degree)'},'FontSize',12,'FontWeight','bold','FontName','Arial');
% % set(xl,'Position',[0 -1.5 -30],'Rotation',10)
% yl=ylabel({'Bandwidth (GHz)'},'FontSize',12,'FontWeight','bold','FontName','Arial');
% % set(yl,'Position',[-123 -0.5 -30],'Rotation',-40)
% zlabel({'Level (dB)'},'FontSize',12,'FontWeight','bold','FontName','Arial');
% tx1=title(['Integrated Antenna Pattern'],'FontSize',12,'FontWeight','bold','FontName','Arial');
% 
% set(xl,'Position',[-25 -1.5 -30],'Rotation',9)
% set(yl,'Position',[-130 -1.2 -30],'Rotation',-40)
% saveas(gca,['7A.fig']);
% 
% af.PaperUnits = 'centimeters';
% af.PaperPosition = [1 0 figurewidth figureheight];
% af.PaperPositionMode = 'manual';
% print(['7A.png'],'-dpng','-r400')
% 

clear

fn='xcorr45_0.045';

fid1=fopen([fn '_0' '.csv']);
data1=textscan(fid1,'%f%f','headerlines',0,'delimiter', ',');
fclose(fid1);
seeangle=data1{1}.';
xcorr0=data1{2}.';

fid1=fopen([fn '_1' '.csv']);
data1=textscan(fid1,'%f%f','headerlines',0,'delimiter', ',');
fclose(fid1);
xcorr1=data1{2}.';

fid1=fopen([fn '_2' '.csv']);
data1=textscan(fid1,'%f%f','headerlines',0,'delimiter', ',');
fclose(fid1);
xcorr2=data1{2}.';

fid1=fopen([fn '_4' '.csv']);
data1=textscan(fid1,'%f%f','headerlines',0,'delimiter', ',');
fclose(fid1);
xcorr4=data1{2}.';

plot(seeangle,xcorr0,seeangle,xcorr1,seeangle,xcorr2,seeangle,xcorr4,'LineWidth',1.2);
xlim([-90,90]);
ylim([-35,2]);
grid on

figurewidth=9;
figureheight=7;
af=gcf;
af.Units='centimeters';
af.Position(3)=figurewidth;
af.Position(4)=figureheight;

ax=gca;
ax.FontSize=12;
ax.XTick=[-90 -60 -30 0 30 60 90];
% ax.YTick=[ -15 -10 -5 0];
ax.LineWidth=1.2;
% ax.OuterPosition=[0 0.07 1 0.93];
xlabel({'观察角度{\theta}(°)';'(c)'},'FontSize',14,'FontName','Microsoft YaHei');
ylabel({' 归一化相关峰值方向图(dB)'},'FontSize',13,'FontName','Microsoft YaHei');% 1y

l1=legend('单音','1 GHz','2 GHz','4 GHz','Location','northeast');
set(l1,'FontSize',6)
l1.Position(1)=l1.Position(1)+0.13;
l1.Position(2)=l1.Position(2);

% saveas(gca,['4F.fig']);
af.PaperUnits = 'centimeters';
af.PaperPosition = [1 0 figurewidth figureheight];
af.PaperPositionMode = 'manual';
print([fn '.png'],'-dpng','-r400')









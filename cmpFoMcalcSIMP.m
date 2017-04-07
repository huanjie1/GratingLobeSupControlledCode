function [ cmpFoM ] = cmpFoMcalcSIMP( drgaxis, mainlobedir, cmpue, antnum, dlratio, fignum )
%   cmpFoMcalcSIMP.m
%   计算相关峰值方向图的优值（等距布阵，单音激励适用）  CMP under evaluation

% spurth=10^(-13/20); %peak sidelobe of rectangular window
% spurth=10^(-13.262/20); %peak sidelobe of rectangular window  0.217233628
% spurth=0.2472;
[~,spurth00] = fminsearch(@(x)sin(antnum*x/2)/sin(x/2)/antnum, 2*pi/antnum*1.5);
spurth0=abs(spurth00);
spurth=spurth0*ones(1,length(drgaxis));
sindsin=sin(antnum*pi*dlratio*(sin(mainlobedir/180*pi)-sin(drgaxis/180*pi)))...
    ./sin(pi*dlratio*(sin(mainlobedir/180*pi)-sin(drgaxis/180*pi)))/antnum;
sindsinenv=1./abs(sin(pi*dlratio*(sin(mainlobedir/180*pi)-sin(drgaxis/180*pi))))/antnum;
sindsinenv(isnan(sindsinenv))=1;
spurth=((spurth<=sindsinenv).*spurth+(spurth>sindsinenv).*sindsinenv)*1.003;%>1 for a float calculation error toleration

beamwidthth=1/sqrt(2);

[~,mlindex]=min(abs(drgaxis-mainlobedir));
% idealcmpnm=idealcmp./idealcmp(mlindex);
idealcmpnm=abs(sindsin);
cmpuenm=cmpue/cmpue(mlindex);

% % main lobe deviation------------------------------------------------------------------------------------------
idealdia=max(idealcmpnm);
cmpuedia=max(cmpuenm);


% % main lobe width------------------------------------------------------------------------------------------
for idealmlleftindex=mlindex:-1:2
    if idealcmpnm(idealmlleftindex)>beamwidthth && idealcmpnm(idealmlleftindex-1)<=beamwidthth
        break;
    end
end
for idealmlrightindex=mlindex:1:(length(drgaxis)-1)
    if idealcmpnm(idealmlrightindex)>beamwidthth && idealcmpnm(idealmlrightindex+1)<=beamwidthth
        break;
    end
end
idealmlw=drgaxis(idealmlrightindex)-drgaxis(idealmlleftindex);

for cmpuemlleftindex=mlindex:-1:2
    if cmpuenm(cmpuemlleftindex)>beamwidthth && cmpuenm(cmpuemlleftindex-1)<=beamwidthth
        break;
    end
end
for cmpuemlrightindex=mlindex:1:(length(drgaxis)-1)
    if cmpuenm(cmpuemlrightindex)>beamwidthth && cmpuenm(cmpuemlrightindex+1)<=beamwidthth
        break;
    end
end
cmpuemlw=drgaxis(cmpuemlrightindex)-drgaxis(cmpuemlleftindex);


% % total spur------------------------------------------------------------------------------------------
idealcmpnmth=idealcmpnm;
idealcmpnmth(idealcmpnmth<spurth)=0;
idealcmpnmth(idealmlleftindex:idealmlrightindex)=0;
idealslgl=sum(idealcmpnmth(idealcmpnmth>spurth)-spurth(idealcmpnmth>spurth))*(drgaxis(2)-drgaxis(1));
idealcmpbaddrg=drgaxis(idealcmpnmth>spurth);

cmpuenmth=cmpuenm;
cmpuenmth(cmpuenmth<spurth)=0;
cmpuenmth(cmpuemlleftindex:cmpuemlrightindex)=0;
cmpueslgl=sum(cmpuenmth(cmpuenmth>spurth)-spurth(cmpuenmth>spurth))*(drgaxis(2)-drgaxis(1));
cmpuebaddrg=drgaxis(cmpuenmth>spurth);


%---------------
cmpFoM=(idealdia/cmpuedia)*(idealmlw/cmpuemlw)*(idealslgl/cmpueslgl);
%---------------

if fignum>0
    figure;
    plot(drgaxis,idealcmpnm);hold on 
%     plot(drgaxis,sindsin,drgaxis,sindsinenv);ylim([-1,1])
%     quiver(beamwidthth,beamwidthth,drgaxis(mlindex),drgaxis(idealmlleftindex),...
%         'MaxHeadSize',0.9,'Color','black' );
%     quiver(beamwidthth,beamwidthth,drgaxis(mlindex),drgaxis(idealmlrightindex),...
%         'MaxHeadSize',0.9,'Color','black' );
    line([drgaxis(idealmlleftindex),drgaxis(idealmlrightindex)],[beamwidthth beamwidthth],...
        'linestyle','--','color','black')
    text(drgaxis(mlindex),beamwidthth*1.1,['beamwidth=' num2str(idealmlw) '°'],...
        'FontSize',14,'HorizontalAlignment','center');
%     area(drgaxis,idealcmpnmth,'FaceAlpha',0.3);
    patch([idealcmpbaddrg, fliplr(idealcmpbaddrg)], ...
        [max(idealcmpnmth(idealcmpnmth>spurth),spurth(idealcmpnmth>spurth)),...
        fliplr(min(idealcmpnmth(idealcmpnmth>spurth),spurth(idealcmpnmth>spurth)))], 'r','FaceAlpha',0.3)
%     line([-90 90],[spurth0 spurth0],'linestyle','--','color','black')
    plot(drgaxis,spurth,'linestyle','--','color','black');
    
    figure;
    plot(drgaxis,cmpuenm);hold on
    line([drgaxis(cmpuemlleftindex),drgaxis(cmpuemlrightindex)],[beamwidthth beamwidthth],...
        'linestyle','--','color','black')
    text(drgaxis(mlindex),beamwidthth*1.1,['beamwidth=' num2str(cmpuemlw) '°'],...
        'FontSize',14,'HorizontalAlignment','center');
%     area(drgaxis,cmpuenmth,'FaceAlpha',0.3);
    patch([cmpuebaddrg, fliplr(cmpuebaddrg)], ...
        [max(cmpuenmth(cmpuenmth>spurth),spurth(cmpuenmth>spurth)),...
        fliplr(min(cmpuenmth(cmpuenmth>spurth),spurth(cmpuenmth>spurth)))], 'r','FaceAlpha',0.3)
%     line([-90 90],[spurth0 spurth0],'linestyle','--','color','black')
    plot(drgaxis,spurth,'linestyle','--','color','black');
end
    

end


function [ cmpFoM ] = cmpFoMcalcV2( drgaxis, mainlobedir, cmpue, antnum, dlratio, fignum )
%   cmpFoMcalcV2.m
%   计算相关峰值方向图的优值  CMP under evaluation


sindsin=sin(antnum*pi*dlratio*(sin(mainlobedir/180*pi)-sin(drgaxis/180*pi)))...
    ./sin(pi*dlratio*(sin(mainlobedir/180*pi)-sin(drgaxis/180*pi)))/antnum;
sindsin(isnan(sindsin))=1;

idealcmpnm1=abs(sindsin)/(sum(abs(sindsin))*(drgaxis(2)-drgaxis(1)));
cmpuenm1=cmpue/(sum(abs(cmpue))*(drgaxis(2)-drgaxis(1)));



% % total spur (var)------------------------------------------------------------------------------------------
idealcmpvar=sum((drgaxis-mainlobedir).^2.*idealcmpnm1)*(drgaxis(2)-drgaxis(1));
cmpuevar=sum((drgaxis-mainlobedir).^2.*cmpuenm1)*(drgaxis(2)-drgaxis(1));

% % main lobe deviation------------------------------------------------------------------------------------------
[~,mlindex]=min(abs(drgaxis-mainlobedir));
idealdia=idealcmpnm1(mlindex);
cmpuedia=cmpuenm1(mlindex);

% % maxspur------------------------------------------------------------------------------------------
% ideallmax=[false diff(sign(diff(idealcmpnm1)))==-2 false];
ideallmin=[false diff(sign(diff(idealcmpnm1)))==2 false];
ideallminindex=find(ideallmin);
[~,ideallmindd1]=min(abs(ideallminindex-mlindex));
if ideallminindex(ideallmindd1)>mlindex
    if 1==ideallmindd1
        regionleft=1;
    else
        regionleft=ideallminindex(ideallmindd1-1);
    end
    regionright=ideallminindex(ideallmindd1);
else if ideallminindex(ideallmindd1)<mlindex
        regionleft=ideallminindex(ideallmindd1);
        if length(ideallminindex)==ideallmindd1
            regionright=length(idealcmpnm1);
        else
            regionright=ideallminindex(ideallmindd1+1);
        end
    else
        regionleft=ideallminindex(ideallmindd1-1);
        regionright=ideallminindex(ideallmindd1+1);
    end
end
% ideallmax(regionleft:regionright)=false;
% idealmaxspur=max(idealcmpnm1(ideallmax));
idealmaxspur=max(idealcmpnm1([1:regionleft regionright:end]));

cmpuelmax=[false diff(sign(diff(cmpuenm1)))== -2 false];
% cmpuelmin=[true diff(sign(diff(cmpuenm1)))== 2 true];
% cmpuelminindex=find(cmpuelmin);
% [~,cmpuelmindd1]=min(abs(cmpuelminindex-mlindex));
% if cmpuelminindex(cmpuelmindd1)>mlindex
%     cmpuelmax(cmpuelminindex(cmpuelmindd1-1):cmpuelminindex(cmpuelmindd1))=false;
% else if cmpuelminindex(cmpuelmindd1)<mlindex
%         cmpuelmax(cmpuelminindex(cmpuelmindd1):cmpuelminindex(cmpuelmindd1+1))=false;
%     else
%         cmpuelmax(cmpuelminindex(cmpuelmindd1-1):cmpuelminindex(cmpuelmindd1+1))=false;
%     end
% end
cmpuelmax(regionleft:regionright)=false;
% cmpuemaxspur=max(cmpuenm1([1:regionleft regionright:end]));
cmpuemaxspur=max(cmpuenm1(cmpuelmax));

if fignum>0
    figure;plot(drgaxis,cmpuenm1,drgaxis,idealcmpnm1,drgaxis,0.01*cmpuelmax,':');
    text(-90,0.02,{['cmpuevar = ' num2str(cmpuevar)];...
                   ['idealcmpvar = ' num2str(idealcmpvar)];...
                   ['cmpuedia = ' num2str(cmpuedia)];...
                   ['idealdia = ' num2str(idealdia)];...
                   ['cmpuemaxspur = ' num2str(cmpuemaxspur)]; ...
                   ['idealmaxspur = ' num2str(idealmaxspur)];})
end
    


cmpFoM = (cmpuevar/idealcmpvar)^(-1)*(cmpuedia/idealdia)*(cmpuemaxspur/idealmaxspur)^(-1);

    

end


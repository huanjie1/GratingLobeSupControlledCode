% AFmaster.m
% arrayfactor2masterFUN.m�ĵ���ɨ�躯��������ɨ���������ԣ����ǲ����γ����磩



arraynum=16;
ratio=1.5;
fc=10e9;
fbw=0.4;

filenamepre=['an' num2str(arraynum) ...
          '_r' num2str(ratio) ... 
          '_fc' num2str(fc/1e9)...
          '_fbw' num2str(fbw)...
          '_step0.1new'];

mainlobedrghalf=linspace(-60,0,121);
cmpfomarrayhalf=zeros(1,length(mainlobedrghalf));

for ind1=1:length(mainlobedrghalf)
   mainlobe= mainlobedrghalf(ind1)
   tic
   cmpfomarrayhalf(ind1)=arrayfactor2masterFUN( arraynum, mainlobedrghalf(ind1), ratio, fc, fbw);
   toc
end

mainlobedrg=[mainlobedrghalf -mainlobedrghalf(end-1:-1:1)];
cmpfomarray=[cmpfomarrayhalf cmpfomarrayhalf(end-1:-1:1)];

figure;plot(mainlobedrg,cmpfomarray);

mout=[mainlobedrg;cmpfomarray];
fid1=fopen([filenamepre '.csv'],'w');
fprintf(fid1,'%.9e,%.9e\n',mout);
fclose(fid1);


    
function [ ringserpresp ] = optiresRINGcalc( ringnum, aimdelay, aimbw, ocf, fsweep )
% optiresRINGcalc.m
% ����΢���Ĺ��׸���Ӧ���㣬���������Ż����
% ΢����Ŀ��Ŀ����ʱ����Ŀ�����������Ƶ�ʣ���ɨƵ��Χ
% ��optiresRINGsys.m����
% ���ü���΢�����㺯��optiresRINGserial.m

a=ringnum+aimdelay+aimbw+ocf+fsweep;

ringserpresp=exp(-1i*2*pi*aimdelay*fsweep);

end
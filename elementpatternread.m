clear

freqangle=[];
fid1=fopen('8-12-8000.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-8500.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-9000.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-9500.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-10000.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-10500.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-11000.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-11500.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

fid1=fopen('8-12-12000.txt');
data1=textscan(fid1,'%f %f','headerlines',9);
fclose(fid1);
drgsee=data1{1}.';amplitude1=data1{2}.';
% amplitude1n=amplitude1/max(amplitude1);
% plot(drgsee(1:end-2),amplitude1(1:end-2));hold on;
freqangle=[freqangle amplitude1.'];

freqanglelin=10.^(freqangle(1:end-1,:)/20);
freqanglelinnor=freqanglelin./max(max(freqanglelin));
imagesc(freqanglelinnor)

freqdata=8e9:0.5e9:12e9;
angledata=(drgsee(1:end-1)/180*pi).';

% example for Interpolation
freqneed=7e9:0.1e9:13e9;
angleneed=(-pi:0.05:pi).';
% 1 is assigned toall queries that lie outside the domain of the sample points
vq=interp2(freqdata,angledata,freqanglelinnor,freqneed,angleneed,'linear',1);
figure;imagesc(vq)
save('antenna_element','freqdata','angledata','freqanglelinnor');

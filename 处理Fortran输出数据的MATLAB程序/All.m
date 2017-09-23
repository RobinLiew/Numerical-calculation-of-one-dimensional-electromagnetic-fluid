clc
clear
%global d; %focus
%cd('C:\Users\lenovo\Documents\Visual Studio 2010\Projects\deira3\deira3\data\sphere\N3')

Time=load('Time.txt');%read the time when printing
N=length(Time);%N is the number of Output
d.Time=Time;
d.Time=10.*d.Time;
d.N=N;
%file=load('step1=1.txt');
PBR=load('PBR.txt')
d.PBR=PBR
for I=1:1:N
    rab=strcat('step1=',num2str(I),'.txt');
    file=load(rab);
    d.J=file(:,1);%cell number
    d.R(I,:)=file(:,2);%radius
    d.U(I,:)=file(:,3);%velocity
    d.rho(I,:)=file(:,4);%density
    d.Te(I,:)=file(:,5);%electric temperature
    d.Ti(I,:)=file(:,6);%ion temperature
    d.Tr(I,:)=file(:,7);%radiation temperature
    d.Bz(I,:)=file(:,8);%axis magnetic field
    d.P(I,:)=file(:,9);%pressure  
end
d.nn=2.408*10^23.*d.rho;


for I=1:1:N
    rab=strcat('step2=',num2str(I),'.txt');
    file=load(rab);
    d.Ioniz(I,:)=file(:,2);%ionization
    d.Cv_e(I,:)=file(:,3);%capability
    d.Cv_i(I,:)=file(:,4);%
    d.Cv_r(I,:)=file(:,5);%
    d.EALF(I,:)=file(:,6);%
    d.HIAL(I,:)=file(:,7);%
    d.XD(I,:)=file(:,8);%
    d.XT(I,:)=file(:,9);%
    d.Q(I,:)=file(:,10);% 
end
d.layer1=load('layer1.txt');
d.layer2=load('layer2.txt');
d.layer3=load('layer3.txt');
d.layer4=load('layer4.txt');
d.layer5=load('layer5.txt');
d.LayerAll=load('LayerAll.txt');
%d.All=LayerAll;

%跟踪的点数，及序数
d.Nplot=5;
d.Nsub(1)=20;
d.Nsub(2)=50;
d.Nsub(3)=29;
d.Nsub(4)=36;
d.Nsub(5)=43;
d


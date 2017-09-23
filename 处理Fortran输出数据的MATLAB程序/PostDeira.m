function varargout = PostDeira(varargin)
% POSTDEIRA MATLAB code for PostDeira.fig
%      POSTDEIRA, by itself, creates a new POSTDEIRA or raises the existing
%      singleton*.
%
%      H = POSTDEIRA returns the handle to a new POSTDEIRA or the handle to
%      the existing singleton*.
%
%      POSTDEIRA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTDEIRA.M with the given input arguments.
%
%      POSTDEIRA('Property','Value',...) creates a new POSTDEIRA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PostDeira_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PostDeira_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PostDeira

% Last Modified by GUIDE v2.5 12-Jun-2014 15:35:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PostDeira_OpeningFcn, ...
                   'gui_OutputFcn',  @PostDeira_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PostDeira is made visible.
function PostDeira_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PostDeira (see VARARGIN)

% Choose default command line output for PostDeira
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PostDeira wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PostDeira_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --- Executes on button press in LoadButt.
function LoadButt_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
clear
global d; %focus
%cd('C:\Users\lenovo\Documents\Visual Studio 2010\Projects\deira3\deira3\data\sphere\N3')

Time=load('Time.txt');%read the time when printing
N=length(Time);%N is the number of Output
d.Time=Time;
d.Time=10.*d.Time;
d.N=N;
%file=load('step1=1.txt');

for I=1:1:d.N
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


for I=1:1:d.N
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
d.Nsub(1)=10;
d.Nsub(2)=20;
d.Nsub(3)=29;
d.Nsub(4)=36;
d.Nsub(5)=43;
d



    






function InterPoint_Callback(hObject, eventdata, handles)
% hObject    handle to InterPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterPoint as text
%        str2double(get(hObject,'String')) returns contents of InterPoint as a double


% --- Executes during object creation, after setting all properties.
function InterPoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RButt.
function RButt_Callback(hObject, eventdata, handles)
% hObject    handle to RButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotR(d)



% --- Executes on button press in UButt.
function UButt_Callback(hObject, eventdata, handles)
% hObject    handle to UButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotU(d);


% --- Executes on button press in XDButt.
function XDButt_Callback(hObject, eventdata, handles)
% hObject    handle to XDButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotXD(d)



% --- Executes on button press in XTButt.
function XTButt_Callback(hObject, eventdata, handles)
% hObject    handle to XTButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotXT(d);


% --- Executes on button press in XHeButt.
function XHeButt_Callback(hObject, eventdata, handles)
% hObject    handle to XHeButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotXHe(d)


% --- Executes on button press in XHButt.
function XHButt_Callback(hObject, eventdata, handles)
% hObject    handle to XHButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in XBButt.
function XBButt_Callback(hObject, eventdata, handles)
% hObject    handle to XBButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in TeButt.
function TeButt_Callback(hObject, eventdata, handles)
% hObject    handle to TeButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotTe(d)


% --- Executes on button press in TiButt.
function TiButt_Callback(hObject, eventdata, handles)
% hObject    handle to TiButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotTi(d);


% --- Executes on button press in TrButt.
function TrButt_Callback(hObject, eventdata, handles)
% hObject    handle to TrButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotTr(d)


% --- Executes on button press in EALButt.
function EALButt_Callback(hObject, eventdata, handles)
% hObject    handle to EALButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotEAL(d);


% --- Executes on button press in Ep3Butt.
function Ep3Butt_Callback(hObject, eventdata, handles)
% hObject    handle to Ep3Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotEp3(d);


% --- Executes on button press in Ep14Butt.
function Ep14Butt_Callback(hObject, eventdata, handles)
% hObject    handle to Ep14Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotEp14(d);



% --- Executes on button press in rhoButt.
function rhoButt_Callback(hObject, eventdata, handles)
% hObject    handle to rhoButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotrho(d);

% --- Executes on button press in Rho1Butt.
function Rho1Butt_Callback(hObject, eventdata, handles)
% hObject    handle to Rho1Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotrho1(d);

% --- Executes on button press in BzButt.
function BzButt_Callback(hObject, eventdata, handles)
% hObject    handle to BzButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotB(d);

% --- Executes on button press in PButt.
function PButt_Callback(hObject, eventdata, handles)
% hObject    handle to PButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotP(d);

% --- Executes on button press in GButt.
function GButt_Callback(hObject, eventdata, handles)
% hObject    handle to GButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotG(d);

% --- Executes on button press in RhoR1Butt.
function RhoR1Butt_Callback(hObject, eventdata, handles)
% hObject    handle to RhoR1Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotRhoR1(d);

% --- Executes on button press in RhoR2Butt.
function RhoR2Butt_Callback(hObject, eventdata, handles)
% hObject    handle to RhoR2Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotRhoR2(d);

% --- Executes on button press in n1Butt.
function n1Butt_Callback(hObject, eventdata, handles)
% hObject    handle to n1Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotn1(d);

% --- Executes on button press in n2Butt.
function n2Butt_Callback(hObject, eventdata, handles)
% hObject    handle to n2Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
plotn2(d);



% *******************函数库BEGIN*******************
function plotR(d)
plot(d.Time,d.R(:,d.Nsub(1)),d.Time,d.R(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
%legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),num2str(d.Njz(1,3)),num2str(d.Njz(1,4)),num2str(d.Njz(1,5)),'Location','NorthWest');
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
title('R-t');
xlabel('t(ns)');
ylabel('R(mm)');


function plotU(d)
plot(d.Time,d.U(:,d.Nsub(1)),d.Time,d.U(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
title('U-t');
xlabel('t(ns)');
ylabel('U(10**7cm/s)');



function plotTe(d)
plot(d.Time,d.Te(:,d.Nsub(1)),d.Time,d.Te(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
title('Te-t');
xlabel('t(ns)');
ylabel('T(keV)');


function plotTi(d)
plot(d.Time,d.Ti(:,d.Nsub(1)),d.Time,d.Ti(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
title('Ti-t');
xlabel('t(ns)');
ylabel('T(keV)');


function plotTr(d)
plot(d.Time,d.Tr(:,d.Nsub(1)),d.Time,d.Tr(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
title('Tr-t');
xlabel('t(ns)');
ylabel('T(keV)');


function plotEAL(d)
plot(d.Time,d.EAL(:,d.Nsub(1)),d.Time,d.EAL(:,d.Nsub(2)),d.Time,d.EAL(:,d.Nsub(3)),d.Time,d.EAL(:,d.Nsub(4)),d.Time,d.EAL(:,d.Nsub(5)))
legend('1','2','3','4','5','Location','NorthWest')
title('EAL-t');
xlabel('t(ns)');
ylabel('E(10**14 ergs cm^3)');


function plotEp3(d)
plot(d.Time,d.Ep3)
legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),num2str(d.Njz(1,3)),num2str(d.Njz(1,4)),num2str(d.Njz(1,5)),'Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
title('Ep3-t');
xlabel('t(ns)');
ylabel('E=10**14 ergs cm^3');


function plotEp14(d)
plot(d.Time,d.EAL(:,d.Nsub(1)),d.Time,d.EAL(:,d.Nsub(2)),d.Time,d.EAL(:,d.Nsub(3)),d.Time,d.EAL(:,d.Nsub(4)),d.Time,d.EAL(:,d.Nsub(5)))
legend('1','2','3','4','5','Location','NorthWest')
title('Ep14-t');
xlabel('t(ns)');
ylabel('E=10**14 ergs cm^3');


function plotXD(d)
plot(d.Time,d.XD(:,d.Nsub(1)),d.Time,d.XD(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
set(gca,'ylim',[-0.1,1.1]);
title('XD-t');
xlabel('t(ns)');
ylabel('X');
% function plotXD(d)
% plot(d.Time,d.XD(:,1),d.Time,d.XD(:,2))
% legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),'Location','NorthWest')
% %axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
% set(gca,'ylim',[-0.1,1.1]);
% title('XD-t');
% xlabel('t(ns)');
% ylabel('X');

function plotXT(d)
plot(d.Time,d.XT(:,d.Nsub(1)),d.Time,d.XT(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
set(gca,'ylim',[-0.1,1.1]);
title('XT-t');
xlabel('t(ns)');
ylabel('X');

% function plotXT(d)
% plot(d.Time,d.XT(:,1),d.Time,d.XT(:,2))
% legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),'Location','NorthWest')
% %axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
% set(gca,'ylim',[-0.1,1.1]);
% title('XT-t');
% xlabel('t(ns)');
% ylabel('X');


% function plotXHe(d)
% plot(d.Time,d.XHe)
% legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),num2str(d.Njz(1,3)),num2str(d.Njz(1,4)),num2str(d.Njz(1,5)),'Location','NorthWest')
% %axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
% title('XHe-t');
% xlabel('t(ns)');
% ylabel('X');
function plotXHe(d)
plot(d.Time,d.XHe(:,1),d.Time,d.XHe(:,2))
legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),'Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
title('XHe-t');
xlabel('t(ns)');
ylabel('X');




function plotrho(d)
plot(d.Time,d.rho(:,d.Nsub(1)),d.Time,d.rho(:,d.Nsub(2)))
legend('1','2','Location','NorthWest')
title('rho-t');
xlabel('t(ns)');
ylabel('rho(g cm^{-3})');

function plotB(d)
plot(d.Time,d.Bz(:,d.Nsub(1)),d.Time,d.Bz(:,d.Nsub(2)))
legend('1','2','3','Location','NorthWest')
title('Bz-t');
xlabel('t(ns)');
ylabel('Bz(10^{7}G)');

function plotP(d)
%plot(d.Time,log(d.P(:,d.Nsub(1)))./log(10))
plot(d.Time,d.P(:,d.Nsub(1)),d.Time,d.P(:,d.Nsub(2)))
legend('1','2','3','Location','NorthWest')
title('P-t');
xlabel('t(ns)');
ylabel('P(100Mbar)'); %1bar=10^5Pa,100Mbar=10^13Pa


function plotG(d)
plot(d.LayerAll(:,46))
%legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),num2str(d.Njz(1,3)),num2str(d.Njz(1,4)),num2str(d.Njz(1,5)),'Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
title('G-t');
xlabel('t(ns)');
ylabel('G');

function plotRhoR1(d)
plot(d.Time,d.layer1(:,5));
%legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),num2str(d.Njz(1,3)),num2str(d.Njz(1,4)),num2str(d.Njz(1,5)),'Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
legend('DT gas','Location','NorthWest')
title('\rhoR-t');
xlabel('t(ns)');
ylabel('\rhoR');

function plotRhoR2(d)
plot(d.Time,d.layer2(:,5));
%legend(num2str(d.Njz(1,1)),num2str(d.Njz(1,2)),num2str(d.Njz(1,3)),num2str(d.Njz(1,4)),num2str(d.Njz(1,5)),'Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
legend('DT ice','Location','NorthWest')
title('\rhoR-t');
xlabel('t(ns)');
ylabel('\rhoR');

function plotn1(d)
plot(d.Time,2.408*10^23.*d.rho(:,1))
legend('DT gas','Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
title('n-t');
xlabel('t(ns)');
ylabel('n(cm^{-3})');

function plotn2(d)
plot(d.Time,2.408*10^23.*d.rho(:,2))
legend('DT ice','Location','NorthWest')
%axis([min(d.CellN1),max(d.CellN1),min(min(d.R)),max(max(d.R))])
title('n-t');
xlabel('t(ns)');
ylabel('n(cm^{-3})');

% *******************函数库END*******************

function varargout = SimithChart(varargin)
% SIMITHCHART MATLAB code for SimithChart.fig
%      SIMITHCHART, by itself, creates a new SIMITHCHART or raises the existing
%      singleton*.
%
%      H = SIMITHCHART returns the handle to a new SIMITHCHART or the handle to
%      the existing singleton*.
%
%      SIMITHCHART('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMITHCHART.M with the given input arguments.
%
%      SIMITHCHART('Property','Value',...) creates a new SIMITHCHART or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SimithChart_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SimithChart_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SimithChart

% Last Modified by GUIDE v2.5 02-Jun-2016 18:57:46

% Begin initialization code - DO NOT EDIT




gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SimithChart_OpeningFcn, ...
                   'gui_OutputFcn',  @SimithChart_OutputFcn, ...
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


% --- Executes just before SimithChart is made visible.
function SimithChart_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SimithChart (see VARARGIN)

% Choose default command line output for SimithChart
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SimithChart wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SimithChart_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Z0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Z0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z0Edit as text
%        str2double(get(hObject,'String')) returns contents of Z0Edit as a double
set(handles.ConfirmButton,'value',1)


% --- Executes during object creation, after setting all properties.
function Z0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZLEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ZLEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZLEdit as text
%        str2double(get(hObject,'String')) returns contents of ZLEdit as a double


% --- Executes during object creation, after setting all properties.
function ZLEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZLEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddElementButton.
function AddElementButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddElementButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowAddElement();

function ShowAddElement()
%a=SimithChartDialog()

[selection,isok] = listdlg('PromptString','Please add an element expected:',...%提示文件
                'SelectionMode','single',...%多选
                'ListString',{'Transmission Line','Resistance','Capacitance','Inductance','Normalized impedance','Normalized reactance'},...%列表选项
                'ListSize',[200 100],...%列表框大小默认 [160 300]
                'InitialValue',3,...%初始选项，默认为第一个1
                'name','Add an element');
if(isok)
    switch(selection)
        case 1
            prompt = {'Transmission  Line Length(factor of wavelength):'};
            def = {'0'};
        case 2
            prompt = {'Resistance(\Omega)'};
            def = {'0'};
    end
    title = 'Input Value';
    lines = 1;
    answer = inputdlg(prompt, title, lines, def);
    assignin('base', 'imfile', answer{1});
    assignin('base', 'cmap', answer{2});
end      

function [ L,D ] = tuneWithParallelSingleStub( ZL,Z0,ifoc )
%returns 1*2 array of each variable indicating two solution of the tuning
%two solutions are same if RL=Z0
%all elements are coeffcients of lambda
RL = real(ZL);
XL = imag(ZL);

if RL == Z0
    t1 = -XL/(2*Z0);
    T = [t1];
else
    t1 = (XL+sqrt(RL*((Z0-RL)^2+XL^2)/Z0))/(RL-Z0);
    t2 = (XL-sqrt(RL*((Z0-RL)^2+XL^2)/Z0))/(RL-Z0);
    T = [t1 t2];
end
D = mod(atan(T)/(2*pi)+0.5,0.5);
B = (RL^2*T-(Z0-XL*T).*(XL+Z0*T))./(Z0*(RL^2+power(XL+Z0*T,2)));
if ifoc
    L = mod(-1/(2*pi)*atan(B*Z0)+0.5,0.5);
else
    L = mod(1/(2*pi)*atan(1./(B*Z0))+0.5,0.5);
end

function [ L1,L2 ] = tuneWithDoubleStubs( ZL,Z0,dl1,d12,ifoc )
%returns 1*2 array of each variable indicating two solution of the tuning
%two solutions are same as the textbook indicates
%all elements are coeffcients of lambda
ZL = Z0*(ZL+j*Z0*tan(2*pi*dl1))/(Z0+j*ZL*tan(2*pi*dl1));
Y0 = 1/Z0;
YL = 1/ZL;
GL = real(YL);
BL = imag(YL);
t = tan(2*pi*d12);

temp=(Y0*(1+t^2)/t^2);
if GL<0||(GL>temp&&abs(GL-temp)>10e-10)
    L1 = [0,0];
    L2 = [0,0];
elseif GL>temp&&abs(GL-temp)<10e-10
    GL=temp; 
    B1 = [-BL+(Y0+sqrt((1+t^2)*GL*Y0-GL^2*t^2))/t -BL+(Y0-sqrt((1+t^2)*GL*Y0-GL^2*t^2))/t];
    B2 = [(Y0*sqrt(Y0*GL*(1+t^2)-GL^2*t^2)+GL*Y0)/(GL*t) (-Y0*sqrt(Y0*GL*(1+t^2)-GL^2*t^2)+GL*Y0)/(GL*t)];
    if ifoc
        L1 = mod(1/(2*pi)*atan(B1./Y0)+0.5,0.5);
        L2 = mod(1/(2*pi)*atan(B2./Y0)+0.5,0.5);
    else
        L1 = mod(-1/(2*pi)*atan(Y0/B1)+0.5,0.5);
        L2 = mod(-1/(2*pi)*atan(Y0/B2)+0.5,0.5);
    end
else
    B1 = [-BL+(Y0+sqrt((1+t^2)*GL*Y0-GL^2*t^2))/t -BL+(Y0-sqrt((1+t^2)*GL*Y0-GL^2*t^2))/t];
    B2 = [(Y0*sqrt(Y0*GL*(1+t^2)-GL^2*t^2)+GL*Y0)/(GL*t) (-Y0*sqrt(Y0*GL*(1+t^2)-GL^2*t^2)+GL*Y0)/(GL*t)];
    if ifoc
        L1 = mod(1/(2*pi)*atan(B1./Y0)+0.5,0.5);
        L2 = mod(1/(2*pi)*atan(B2./Y0)+0.5,0.5);
    else
        L1 = mod(-1/(2*pi)*atan(Y0./B1)+0.5,0.5);
        L2 = mod(-1/(2*pi)*atan(Y0./B2)+0.5,0.5);
    end

end

function  Axes_Initial(hObject)

set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');

%set(handles.axes10,'xtick',[],'ytick',[]);
hold on;
    
range=[0 0.3 1 2.5 10];
%RL=str2num(ZL)/str2num(Z0);
for RL=range
    [Fr,Fi] = getResistanceCircle(RL);
    plot(Fr,Fi,'color',[.5 .5 .5]);
    xlim([-1 1]);
    ylim([-1 1]);
end

range=[-2 -1 -0.3 0 0.3 1 2];
%RL=str2num(ZL)/str2num(Z0);
for XL=range

    [Fr, Fi] = getReactanceCircle(XL);

    plot(Fr,Fi,'color',[.5 .5 .5]);

    xlim([-1 1]);
    ylim([-1 1]);
end

function plotSections(handles,x,y,color,w)
Object=handles.MainAxes;
var=get(handles.DirectButton,'Value');
global rate;
if ~var
    if(isempty(rate))
        rate=1;
    end
    %计算帧数
    max=40*2;
    %帧率1/40
    dt=1/40*rate;
    %块长
    sl=length(x)/max;
    %分块画图
    for k=1:max
        init=(k-1)*sl+1;
        if k>=2
            plot(Object,x(init-1:init+sl-1),y(init-1:init+sl-1),'color',color,'LineWidth',w);
        else
            plot(Object,x(init:init+sl-1),y(init:init+sl-1),'color',color,'LineWidth',w);
        end
        pause(dt);
    end
    %获得图像所有元素
    h1=get(Object,'Children');
    %分块部分删除
    delete(h1(1:max));
end
%画一个总块
plot(Object,x,y,'color',color,'LineWidth',w);


function [Fr, Fi] = getResistanceCircle(RL)
    theta=linspace(0,2*pi,800);
    r=1/(1+RL);
    Fr=r*cos(theta)+RL/(1+RL);
    Fi=r*sin(theta); 
    
function [Fr,Fi] =getResistancePart(z1,z2)
    RL=real(z1);
    xc=RL/(1+RL);
    [m1,n1]=getLocation(z1);
    [m2,n2]=getLocation(z2);
    theta1=atan2(n1,m1-xc);
    theta2=atan2(n2,m2-xc);
    if((theta2-theta1)>=pi)
        theta2=theta2-2*pi;
    elseif((theta2-theta1)<=-pi)
        theta2=theta2+2*pi;
    end
    theta=linspace(theta1,theta2,800);
    r=1/(1+RL);
    Fr=r*cos(theta)+xc;
    Fi=r*sin(theta);
    
function [Fr,Fi] =getConductanceCircle(RL)
    [Fr,Fi]=getResistanceCircle(RL);
    Fr=-Fr;
    
function [Fr,Fi] =getConductancePart(z1,z2)
    [Fr1,Fi1]=getLocation(z1);
    [Fr2,Fi2]=getLocation(z2);
    z1=getImpedance(-Fr1,Fi1);
    z2=getImpedance(-Fr2,Fi2);
    [Fr,Fi] =getResistancePart(z1,z2);
    Fr=-Fr;

function [Fr, Fi] = getReactancePart(z1,z2)
m=imag(z1);
start=real(z1);
stop=real(z2);
k=linspace(start,stop,800);
[Fr,Fi]=getLocation(k+1j*m);

function [Fr, Fi] = getSusceptancePart(z1,z2)
m=imag(1/z1);
start=real(1/z1);
stop=real(1/z2);
k=linspace(start,stop,800);
k=1./(k+1j*m);
[Fr,Fi]=getLocation(k);

    
function [Fr, Fi] = getReactanceCircle(XL)
    if(XL==0)
       Fr=linspace(-1,1,800);
       Fi=zeros(1,800);
       return;
    end
    [Fr,Fi]=getLocation(1j*XL);
    if Fi>0
        t=atan2(Fi-1/XL,Fr-1);
        if t<0
            theta=linspace(-pi/2,t,800);
        else
            theta=linspace(-pi/2,t-2*pi,800);
        end
    else
        t=atan2(Fi-1/XL,Fr-1);
        if t<0
            theta=linspace(pi/2,t+2*pi,800);
        else
            theta=linspace(pi/2,t,800);
        end
    end

    r=abs(1/(XL));

    Fr=r*cos(theta)+1;
    Fi=r*sin(theta)+1/XL;
    %Fr(Fr.^2+Fi.^2>1)=NaN;
    
function [Fr, Fi] = getVSWRCircle(ZL)
    [m,n]= getLocation(ZL);
    offset=atan2(n,m);
    theta=linspace(0,2*pi,800)+offset;
    theta=fliplr(theta);

    r=sqrt(m^2+n^2);
    Fr=r*cos(theta);
    Fi=r*sin(theta);
function [Fr, Fi] = getVSWRPart(z1,z2)
    [m1,n1]=getLocation(z1);
    [m2,n2]=getLocation(z2);
    theta1=atan2(n1,m1);
    theta2=atan2(n2,m2);
    if(theta2>theta1)
        theta2=theta2-2*pi;
    end
    theta=linspace(theta1,theta2,800);
    r=sqrt(m1^2+n1^2);
    Fr=r*cos(theta);
    Fi=r*sin(theta);
    
function z= getImpedance(Fr,Fi)
F=Fr+Fi*1j;
z=(1+F)/(1-F);

function [m,n]= getLocation(z)
F=(z-1)./(z+1);
m=real(F);
n=imag(F);
% Hint: get(hObject,'Value') returns toggle state of Parallel


% --- Executes during object creation, after setting all properties.
function MainAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MainAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
Axes_Initial(hObject)
% Hint: place code in OpeningFcn to populate MainAxes


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ForwardButton.
function ForwardButton_Callback(hObject, eventdata, handles)
% hObject    handle to ForwardButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global history;
local=history.local;
history.length;

DisableButtonState(handles);

if CanGoForward()    
    m=[1 0;0 1];
    

    
    %如果下一个是回到原点
    if cell2mat(history.mode(local+1))==0
        temp=cell2mat(history.infor(1));
        [Fr, Fi]=getLocation(temp(1)/temp(2));
        ScatterHere(handles,Fr,Fi,cell2mat(history.text(local+1)));
        GoForward();
        PlotElement(handles,local+1,[]);
        UpdateHistoryState(handles)
        return;
    end
    %其他情况，先算起点 
    for k=local:-1:1
        if(cell2mat(history.mode(k))==0)   
            temp=cell2mat(history.infor(1));
            m=m*temp;
            break;
        else
            temp=GetMatrix(cell2mat(history.infor(k)));
            m=m*temp;
        end
    end
    PlotElement(handles,local+1,cell2mat(history.infor(local+1)));
    z1=m(1)/m(2);
    temp=GetMatrix(cell2mat(history.infor(local+1)));
    m=temp*m;
    z2=m(1)/m(2);
    mode=cell2mat(history.mode(local+1));
    GoForward();
    
	[Fr,Fi]=GetCurve(z1,z2,mode);
    plotSections(handles,Fr,Fi,rand(1,3),3);
    [Fr,Fi]=getLocation(z2);
    ScatterHere(handles,Fr,Fi,cell2mat(history.text(local+1)));
    

end
UpdateHistoryState(handles)
function [x,y]=PlotElement(handles,index,infor)
if isempty(infor)
       x=[-0.5,0,0,-0.15,-0.15,0.15,0.15,-0.15,-0.15,0,0,-0.5;...
         -0.5,0,0,-0.15,-0.15,0.15,0.15,-0.15,-0.15,0,0,-0.5];
       y=[0.75,0.75,0.65,0.65,0.35,0.35,0.65,0.65,0.35,0.35,0.25,0.25;...
        0.75,0.75,0.65,0.65,0.35,0.35,0.65,0.65,0.35,0.35,0.25,0.25;];
else
line=infor.line;
serial=infor.serial;
    if line
        x=[-0.5,-0.4 -0.4 0.4 0.4 -0.4,-0.4 0.4 0.4, 0.5];
        y=[0.75,0.75,0.8,0.8,0.7,0.7,0.8,0.8,0.75,0.75];

        x(2,1:10)=x(1,1:10);
        y(2,1:10)=y(1,1:10)-0.5;
    elseif serial
        x=[-0.5,-0.15 -0.15 0.15 0.15 -0.15,-0.15 0.15 0.15, 0.5];
        y=[0.75,0.75,0.9,0.9,0.6,0.6,0.9,0.9,0.75,0.75];

        x(2,1:10)=linspace(-0.5,0.5,10);
        y(2,1:10)=ones(1,10)*0.25;
    else
        x=[-0.5,0,0,-0.15,-0.15,0,0,-0.5;...
            0.5,0,0,0.15,0.15,0,0,0.5];
        y=[0.75,0.75,0.65,0.65,0.35,0.35,0.25,0.25;...
            0.75,0.75,0.65,0.65,0.35,0.35,0.25,0.25];
    end
end
x=x-index+1;

plot(handles.PlotAxes,x(1,1:end),y(1,1:end),x(2,1:end),y(2,1:end),'Color','b');


function UpdateHistoryState(handles)
if CanGoForward()
    set(handles.ForwardButton,'Enable','On');
else
    set(handles.ForwardButton,'Enable','Off');
end
if CanGoBack()
    set(handles.BackButton,'Enable','On');
else
    set(handles.BackButton,'Enable','Off');
end
global history;
local=history.local;
set(handles.HistoryListBox,'Value',local+2);

function DisableButtonState(handles)
set(handles.ForwardButton,'Enable','Off')
set(handles.BackButton,'Enable','Off')


% --- Executes on button press in BackButton.
function BackButton_Callback(hObject, eventdata, handles)
% hObject    handle to BackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global history;
local=history.local;

DisableButtonState(handles);

if CanGoBack()
    h=get(handles.PlotAxes,'Children');
    delete(h(1:2));
    if cell2mat(history.mode(local))==0
        h=get(handles.MainAxes,'Children');
        delete(h(1:2));
    else
        h=get(handles.MainAxes,'Children');
        delete(h(1:3));
    end
    GoBack();
end

UpdateHistoryState(handles);


% --- Executes on button press in InitializationButton.
function InitializationButton_Callback(hObject, eventdata, handles)
% hObject    handle to InitializationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

HasConfirm(handles,'off')

axes(handles.MainAxes);
cla reset
Axes_Initial(handles.MainAxes);
PopHistoryAll(handles);
UpdateHistoryState(handles);

h=get(handles.PlotAxes,'Children');
delete(h);


% --- Executes on key press with focus on Z0Edit and none of its controls.
function Z0Edit_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Z0Edit (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ConfirmButton.
function ConfirmButton_Callback(hObject, eventdata, handles)
% hObject    handle to ConfirmButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
HasConfirm(handles,'on');

Z0=str2double(get(handles.Z0Edit,'String'));
ZL=str2double(get(handles.ZLEdit,'String'));
z=ZL/Z0;
[Fr,Fi]=getLocation(z);
axes(handles.MainAxes);
hold on;

ScatterHere(handles,Fr,Fi,'Z_L');
global history;

if length(history)==0;
    t.local=1;
    t.length=1;
    t.mode={0};
    t.infor={[z;1]};
    t.text={'Z_L'};
    history=t;
    AddHistotyList(handles,'SetImpedance');
else
    history.local=1;
    history.infor(1)={[z;1]};
end
PlotElement(handles,1,[]);

%     history.local=1;
%     history.infor={[z;1]};

UpdateHistoryState(handles);


function AddHistory(handles,m,mode,des,showtext)
global history;
history.local=history.local;
history.length=history.length+1;
history.mode=[history.mode mode];
history.infor=[history.infor m];
history.text=[history.text des];

AddHistotyList(handles,showtext)

function value=PopHistory(handles)
global history;
history.length
if(history.length==1)
    value=0;
    return;
end

history.length=history.length-1;

if history.local>history.length
    history.local=history.length;
end
history.mode=history.mode(1:end-1);
history.text=history.text(1:end-1);
history.infor=history.infor(1:end-1);

value=1;
RemoveHistotyList(handles);


function value=RemoveHistoryOnce(handles,index)
global history;
if(history.length==1)
    value=0;
    return;
end

history.length=history.length-1;
history.local=index-1;
history.mode=[history.mode(1:index-1) history.mode(index+1:end-1)];
history.text=[history.text(1:index-1) history.text(index+1:end-1)];
history.infor=[history.infor(1:index-1) history.infor(index+1:end-1)];

value=1;

RemoveHistotyList(handles);

function PopHistoryAll(handles)
while(PopHistory(handles))
end


function GoBack()
global history;
if history.local>1
    history.local=history.local-1;
end

function GoForward()
global history;
if history.local<history.length
    history.local=history.local+1;
end
function v=CanGoBack()
global history;
v=(history.local>1);

function v=CanGoForward()
global history;
v=(history.local<history.length);



% --- Executes on key press with focus on ZLEdit and none of its controls.
function ZLEdit_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to ZLEdit (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in ListBox.
function ListBox_Callback(hObject, eventdata, handles)
% hObject    handle to ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles)

function SetState(handles)
index=get(handles.ListBox,'value');
switch index
    case {1,3,6,7}
        set(handles.FreqEdit,'Enable','Off');
    otherwise
        set(handles.FreqEdit,'Enable','On');
end
if index>2
    set(handles.SCButton,'Enable','Off');
    set(handles.OCButton,'Enable','Off');
    set(handles.LineButton,'Enable','Off');
    set(handles.LoadButton,'Enable','Off');
    set(handles.SerialButton,'Enable','On');
    set(handles.ParallelButton,'Enable','On');
else
    set(handles.LineButton,'Enable','On');
    set(handles.LoadButton,'Enable','On');
    if get(handles.LoadButton,'Value')
        set(handles.SerialButton,'Enable','On');
        set(handles.ParallelButton,'Enable','On');
        pa=get(handles.ParallelButton,'value');
        if(pa)
            set(handles.SCButton,'Enable','on');
            set(handles.OCButton,'Enable','on');
        else
            set(handles.SCButton,'Enable','Off');
            set(handles.OCButton,'Enable','Off');
        end
    else
        set(handles.SCButton,'Enable','Off');
        set(handles.OCButton,'Enable','Off');
        set(handles.SerialButton,'Enable','Off');
        set(handles.ParallelButton,'Enable','Off');
    end

end

% Hints: contents = cellstr(get(hObject,'String')) returns ListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListBox


% --- Executes during object creation, after setting all properties.
function ListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ValueEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ValueEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ValueEdit as text
%        str2double(get(hObject,'String')) returns contents of ValueEdit as a double


% --- Executes during object creation, after setting all properties.
function ValueEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ValueEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FreqEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqEdit as text
%        str2double(get(hObject,'String')) returns contents of FreqEdit as a double


% --- Executes during object creation, after setting all properties.
function FreqEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PreviewButton.
function PreviewButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AddButton.
function AddButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value=str2num(get(handles.ValueEdit,'String'));
freq=str2num(get(handles.FreqEdit,'String'));
index=get(handles.ListBox,'Value');
serial=get(handles.SerialButton,'Value');
line=get(handles.LineButton,'Value');
sc=get(handles.SCButton,'Value');
z0=str2num(get(handles.Z0Edit,'String'));
infor=GetInfor(z0,value,freq,index,serial,line,sc);

[matrix,modes,des]=GetMatrix(infor);


AddHistory(handles,infor,modes,'Z_t_e_m_p',des);
UpdateHistoryState(handles);

function ScatterHere(handles,x,y,des)
des={des};
scatter(handles.MainAxes,x,y,75,'Filled');
if get(handles.ZShowBox,'Value');
    z=getImpedance(x,y);
    des=[des,num2str(z,4)];
end
if get(handles.YShowBox,'Value');
    z=getImpedance(x,y);
    des=[des,num2str(1/z,4)];
end
text(x,y,des);

function result=GetInfor(z0,value,freq,index,serial,line,sc)
result.z0=z0;
result.value=value;
result.freq=freq;
result.index=index;
result.serial=serial;
result.line=line;
result.sc=sc;
if index>2
    result.line=0;
end

%Mode=1 vswr
%Mode=2 电导圈Conductance电纳变
%mode=3 电阻圈resistance电抗变
%mode=4 点纳圈susceptance电导变
%Mode=5 电抗圈reactance电阻变
%Mode=6 直接

function [matrix,modes,des]=GetMatrix(infor)
z0=infor.z0;
value=infor.value;
freq=infor.freq;
index=infor.index;
serial=infor.serial;
line=infor.line;
sc=infor.sc;

switch index
	case 1
        theta=2*pi*value;
        if line
         	A=cos(theta);
            B=1j*sin(theta);
            C=1j*sin(theta);
            D=cos(theta);
            matrix=[A B;C D];
            modes=1;%VSWR
            des='Add TL as Line';
        else
             if sc
                    Z=1j*tan(theta);
             else
                    Z=-1j*cot(theta);
             end
             if serial
                matrix=[1 Z; 0 1];
                modes=3;
                des='Add Serial TL as Load';
             else
                matrix=[1 0; 1/Z 1];
                modes=2;
                des='Add Parallel TL as Load';
            end
        end
	case 2
         lambda=3e8/freq;
         theta=2*pi/lambda*value;
         if line
                A=cos(theta);
                B=1j*sin(theta);
                C=1j*sin(theta);
                D=cos(theta);
                matrix=[A B;C D];
                modes=1;
                des='Add TL as Line';
         else
            if sc
                Z=1j*tan(theta);
            else
                Z=-1j*cot(theta);
            end
            if serial
                matrix=[1 Z; 0 1];
                modes=3;
                des='Add Serial TL as Load';
            else
                matrix=[1 0; 1/Z 1];
                modes=2;
                des='Add Parallel TL as Load';
            end
         end
    case 3
        Z=value/z0;
        if serial
            matrix=[1 Z;0 1];
            modes=5;
            des='Add Serial Resistance';
        else
            matrix=[1 0; 1/Z 1];
            modes=4;
            des='Add Parallel Resistance';
        end
    case 4
        Z=1/(2j*pi*freq*value);
        if serial
            matrix=[1 Z;0 1];
            modes=3;
            des='Add Serial Capacitance';
        else
            matrix=[1 0; 1/Z 1];
            modes=2;
            des='Add Parallel Capacitance';
        end
    case 5
        Z=2j*pi*freq*value;
        if serial
            matrix=[1 Z;0 1];
            modes=3;
            des='Add Serial Inductance';
        else
            matrix=[1 0; 1/Z 1];
            modes=2;
            des='Add Parallel Inductance';
        end
    case 6
        Z=value;
        if serial
            matrix=[1 Z;0 1];
            modes=6;
            des='Add Serial Impedance(Normalized)';
        else
            matrix=[1 0; 1/Z 1];
            modes=6;
            des='Add Parallel Impedance(Normalized)';
        end
    case 7
        Z=1/value;
        if serial
            matrix=[1 Z;0 1];
            modes=6;
            des='Add Serial Adimittance(Normalized)';
        else
            matrix=[1 0; 1/Z 1];
            modes=6;
            des='Add Parallel Adimittance(Normalized)';
        end
end

%Mode=1 vswr
%Mode=2 电导圈Conductance电纳变
%mode=3 电阻圈resistance电抗变
%mode=4 点纳圈susceptance电导变
%Mode=5 电抗圈reactance电阻变
%Mode=6 直接

function [Fr,Fi]=GetCurve(z1,z2,modes)
modes
if modes==1
    [Fr,Fi]=getVSWRPart(z1,z2);
elseif modes==2
    [Fr,Fi]=getConductancePart(z1,z2);
elseif modes==3
    [Fr,Fi]=getResistancePart(z1,z2);
elseif modes==4
    [Fr,Fi]=getSusceptancePart(z1,z2);
elseif modes==5
    [Fr,Fi]=getReactancePart(z1,z2);
else
    r1=real(z1);
    i2=imag(z2);
    z3=r1+1j*i2;
    [Fr1,Fi1]=getResistancePart(z1,z3);
    [Fr2,Fi2]=getReactancePart(z3,z2);
    Fr = [Fr1 Fr2];
    Fi = [Fi1 Fi2];
   
end





% --- Executes on button press in LineButton.
function LineButton_Callback(hObject, eventdata, handles)
% hObject    handle to LineButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles);

% Hint: get(hObject,'Value') returns toggle state of LineButton


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles)
% Hint: get(hObject,'Value') returns toggle state of LoadButton


% --- Executes on button press in SCButton.
function SCButton_Callback(hObject, eventdata, handles)
% hObject    handle to SCButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles)
% Hint: get(hObject,'Value') returns toggle state of SCButton


% --- Executes on button press in OCButton.
function OCButton_Callback(hObject, eventdata, handles)
% hObject    handle to OCButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles)
% Hint: get(hObject,'Value') returns toggle state of OCButton


% --- Executes on button press in SerialButton.
function SerialButton_Callback(hObject, eventdata, handles)
% hObject    handle to SerialButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles)
% Hint: get(hObject,'Value') returns toggle state of SerialButton


% --- Executes on button press in ParallelButton.
function ParallelButton_Callback(hObject, eventdata, handles)
% hObject    handle to ParallelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetState(handles)
% Hint: get(hObject,'Value') returns toggle state of ParallelButton


% --- Executes on selection change in HistoryListBox.
function HistoryListBox_Callback(hObject, eventdata, handles)
% hObject    handle to HistoryListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HistoryListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HistoryListBox

UpdateHistoryState(handles);



% --- Executes during object creation, after setting all properties.
function HistoryListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HistoryListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in SSSButton.
function AddHistotyList(handles,str)
value=get(handles.HistoryListBox,'String');
len=length(value);
str=sprintf('[%d]%s',len-2+1,str);
set(handles.HistoryListBox,'String',[value;str]);

function RemoveHistotyList(handles)
str=get(handles.HistoryListBox,'String');

set(handles.HistoryListBox,'String',str(1:end-1));
UpdateHistoryState(handles);

function HasConfirm(handles,state)
set(handles.SSSButton,'Enable',state)
set(handles.DSSButton,'Enable',state)
set(handles.AddButton,'Enable',state)




function SSSButton_Callback(hObject, eventdata, handles)
% hObject    handle to SSSButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[Selection,ok] = listdlg('liststring',{'OC','SC'}, 'promptstring','Mode','name','Select','listsize',[180 40],'selectionmode','Single');

if ~ok
    return;
end

PopHistoryAll(handles);

Z0=get(handles.Z0Edit, 'String');
ZL=get(handles.ZLEdit, 'String');
Z0=str2num(Z0);
ZL=str2num(ZL);
z=ZL/Z0;
%R=1 reactance circle

ifoc=Selection==1;
[L,D]=tuneWithParallelSingleStub(ZL,Z0,ifoc );

global history;
local=history.local;

set(handles.ResultEdit,'String',{'L=',num2str(L),'D=',num2str(D)}')

[Fr, Fi] = getResistanceCircle(1);
plotSections(handles,-Fr,Fi,[1 0 0],1);

for k=1:length(D)

if k~=1
    AddHistory(handles,[NaN; NaN],0,'Z_L','Set Impedance');
end
    %m=GetMatrix(
    m=GetInfor(Z0,D(k),0,1,0,1,0);
    AddHistory(handles,m,1,'Z_1','Add a transmission line');
   
    %GetMatrix(value,freq,index,serial,line,sc)
    m=GetInfor(Z0,L(k),0,1,0,0,~ifoc);
    AddHistory(handles,m,2,'Z_2','Add a parallel TL as Load');

end
history.local=local;
UpdateHistoryState(handles);



% --- Executes on button press in DSSButton.
function DSSButton_Callback(hObject, eventdata, handles)
% hObject    handle to DSSButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

result=inputdlg({'D_STUB1','D_STUB2','1 For OC; 0 For SC'},'Input parameters',1,{'1/8','1/8','1'});  
if isempty(result)
    return;
end

PopHistoryAll(handles);

global history;
local=history.local;

Z0=get(handles.Z0Edit, 'String');
ZL=get(handles.ZLEdit, 'String');
Z0=str2num(Z0);
ZL=str2num(ZL);

z=ZL/Z0;

D1=str2num(cell2mat(result(1)));
D2=str2num(cell2mat(result(2)));
ifoc=str2num(cell2mat(result(3)));

[L1, L2]= tuneWithDoubleStubs( ZL,Z0,D1,D2,ifoc );
set(handles.ResultEdit,'String',{'L1=',num2str(L1),'L2=',num2str(L2)}')


%R=1 resistance circle
[Fr, Fi] = getConductanceCircle(1);
plotSections(handles,Fr,Fi,[1 0 0],1);

% circle
F=Fr+Fi*1j;
F=F*exp(1j*D2/0.5*2*pi);
Fr=real(F);
Fi=imag(F);
plotSections(handles,Fr,Fi,[0 1 0],1);



%get LL




for k=1:length(L1)

if k~=1
    AddHistory(handles,[NaN;NaN],0,'Z_L','Set Impedance');
    %AddHistotyList(handles,'Set Impedance')
end

    %m=GetMatrix(
    m=GetInfor(Z0,D1,0,1,0,1,0);
    AddHistory(handles,m,1,'Z_1','Add a transmission line');
    %AddHistotyList(handles,)
   
    %GetMatrix(value,freq,index,serial,line,sc)
    %m=GetMatrix(
    m=GetInfor(Z0,L1(k),0,1,0,0,~ifoc);
    AddHistory(handles,m,2,'Z_2','Add a parallel TL as Load');
    %AddHistotyList(handles,'Add a parallel TL as Load')

    %m=GetMatrix(
    m=GetInfor(Z0,D2,0,1,0,1,0);
    AddHistory(handles,m,1,'Z_1','Add a transmission line');
    %AddHistotyList(handles,'Add a transmission line')
   
    %GetMatrix(value,freq,index,serial,line,sc)
    %m=GetMatrix(
    m=GetInfor(Z0,L2(k),0,1,0,0,~ifoc);
    AddHistory(handles,m,2,'Z_2','Add a transmission line');
    %AddHistotyList(handles,'Add a parallel TL as Load')

    
end
history.local=local;
UpdateHistoryState(handles);



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
movegui(gcf,'center');


function ResultEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ResultEdit as text
%        str2double(get(hObject,'String')) returns contents of ResultEdit as a double


% --- Executes during object creation, after setting all properties.
function ResultEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function SSSButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SSSButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
clear global history;
clear global rate;
delete(hObject);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global rate;
rate=get(hObject,'Value');


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in YShowBox.
function YShowBox_Callback(hObject, eventdata, handles)
% hObject    handle to YShowBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of YShowBox


% --- Executes on button press in ZShowBox.
function ZShowBox_Callback(hObject, eventdata, handles)
% hObject    handle to ZShowBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ZShowBox


% --------------------------------------------------------------------
function HelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AboutMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AboutMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str={'This MALAB application is created by','FANG Ge','and','HAN Fangzhou','From South University of Science and Technology of China'};
msgbox(str,'Simth Chart','none') 


% --- Executes on button press in ResetButton.
function ResetButton_Callback(hObject, eventdata, handles)
% hObject    handle to ResetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over HistoryListBox.
function HistoryListBox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to HistoryListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OperateMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OperateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DeleteMenu_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

result=inputdlg({'[end] means last one'},'Input index',1,{'end'});  
if isempty(result)
    return;
end
if  strcmp(result,'end')
    PopHistory(handles);
else
    d=str2int(result);
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PopHistory(handles);



% --- Executes during object creation, after setting all properties.
function PlotAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');
xlim([-8.5,0.5])
ylim([0,1])
hold on
% Hint: place code in OpeningFcn to populate PlotAxes

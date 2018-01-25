function varargout = control_gui(varargin)
% CONTROL_GUI MATLAB code for control_gui.fig
%      CONTROL_GUI, by itself, creates a new CONTROL_GUI or raises the existing
%      singleton*.
%
%      H = CONTROL_GUI returns the handle to a new CONTROL_GUI or the handle to
%      the existing singleton*.
%
%      CONTROL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTROL_GUI.M with the given input arguments.
%
%      CONTROL_GUI('Property','Value',...) creates a new CONTROL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before control_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to control_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help control_gui

% Last Modified by GUIDE v2.5 28-Jun-2015 23:10:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @control_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @control_gui_OutputFcn, ...
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


% --- Executes just before control_gui is made visible.
function control_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to control_gui (see VARARGIN)

%% Initialise Session

handles.devices = daq.getDevices;				% scan DAQ devices
handles.session = daq.createSession('ni');		% create session
handles.session.Rate = 5000;
handles.digi_2_value = 0;
handles.gearRatio = 1296/1;
handles.motorEnabled = 0;

%% Create Channels

% Digital Output:	P0.0 - J5 1 Digital In 1 - rpmField
digi_1 = addDigitalChannel(handles.session, 'dev1', 'Port0/Line0', 'OutputOnly');
% Digital Output:	P0.1 - J5 2 Digital In 2 - 
digi_2 = addDigitalChannel(handles.session, 'dev1', 'Port0/Line1', 'OutputOnly');
% Digital I/O:		P0.2 - J5 3 Digital I/O 1
digi_3 = addDigitalChannel(handles.session, 'dev1', 'Port0/Line2', 'InputOnly');
% Digital I/O:		P0.3 - J5 4 Digital I/O 2
digi_4 = addDigitalChannel(handles.session, 'dev1', 'Port0/Line3', 'InputOnly');

% Analogue Output:	AO.0 - J6 1 Analogue In 1+
analOut_1 = addAnalogOutputChannel(handles.session, 'dev1', 'ao0', 'Voltage');
% Analogue Output:	AO.1 - J6 3 Analogue In 2+
analOut_2 = addAnalogOutputChannel(handles.session, 'dev1', 'ao1', 'Voltage');
% Analogue Input:	AI.0 - J6 5 Analogue Out 1
analIn_1 = addAnalogInputChannel(handles.session, 'dev1', 'ai0', 'Voltage');
% Analogue Input:	AI.1 - J6 6 Analogue Out 2
analIn_2 = addAnalogInputChannel(handles.session, 'dev1', 'ai1', 'Voltage');

% Choose default command line output for control_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes control_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = control_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in enableButton.
function enableButton_Callback(hObject, eventdata, handles)
motorEnabledInput = get(handles.enableButton,'Value');
handles.rpm = str2double(get(handles.rpmField, 'String'));
handles.currentlimit = str2double(get(handles.currentlimitField, 'String'));
guidata(hObject, handles);
control_gui_DataSyncHost(hObject, 'motorEnabled', motorEnabledInput, handles)

function rpmField_Callback(hObject, eventdata, handles)
rpmInput = str2double(get(handles.rpmField, 'String'));
control_gui_DataSyncHost(hObject, 'rpm', rpmInput, handles)

% --- Executes on slider movement.
function rpmSlider_Callback(hObject, eventdata, handles)
rpmInput = round(get(handles.rpmSlider, 'Value'), -2, 'Decimal');
control_gui_DataSyncHost(hObject, 'rpm', rpmInput, handles)

function currentlimitField_Callback(hObject, eventdata, handles)
currentlimitInput = str2double(get(handles.currentlimitField, 'String'));
control_gui_DataSyncHost(hObject, 'currentlimit', currentlimitInput, handles)

% --- Executes on slider movement.
function currentlimitSlider_Callback(hObject, eventdata, handles)
currentlimitInput = round(get(handles.currentlimitSlider, 'Value'), 1, 'Decimal');
control_gui_DataSyncHost(hObject, 'currentlimit', currentlimitInput, handles)

function degField_Callback(hObject, eventdata, handles)
degInput = str2double(get(handles.degField,'String'));
control_gui_DataSyncHost(hObject, 'deg', degInput, handles)

function radField_Callback(hObject, eventdata, handles)
radInput = str2double(get(handles.radField,'String'));
control_gui_DataSyncHost(hObject, 'rad', radInput, handles)

function revField_Callback(hObject, eventdata, handles)
revInput = str2double(get(handles.revField,'String'));
control_gui_DataSyncHost(hObject, 'rev', revInput, handles)

% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate_motor_const(handles)

function control_gui_DataSyncHost(hObject, structSpecifier, valueInput, handles)
switch structSpecifier
	case 'motorEnabled'
		handles.motorEnabled = valueInput;
		refresh_motor(handles)
	case 'rpm'
		if (-40000 <= valueInput) && (valueInput <= 40000)
			set(handles.rpmSlider,'Value', valueInput)
			set(handles.rpmField, 'String', num2str(valueInput))
			handles.rpm = valueInput;
			display('New RPM Accepted')
			refresh_motor(handles)
		else
			display('Enter Valid RPM')
		end
		handles.rpm = valueInput;
	case 'currentlimit'
		if (0 <= valueInput) && (valueInput <= 12)		
			set(handles.currentlimitSlider,'Value', valueInput)
			set(handles.currentlimitField, 'String', num2str(valueInput))
			handles.currentlimit = valueInput;
			display('New Current Limit Accepted')
			refresh_motor(handles)
		else
			display('Enter Valid Current Limit')
		end
		handles.currentlimit = valueInput;
	case 'deg'
		handles.angle = valueInput/360*(2*pi);
		set(handles.degField, 'String', num2str(handles.angle/(2*pi)*360))
		set(handles.radField, 'String', num2str(handles.angle))
		set(handles.revField, 'String', num2str(handles.angle/(2*pi)))
	case 'rad'
		handles.angle = valueInput;
		set(handles.degField, 'String', num2str(handles.angle/(2*pi)*360))
		set(handles.radField, 'String', num2str(handles.angle))
		set(handles.revField, 'String', num2str(handles.angle/(2*pi)))
	case 'rev'
		handles.angle = valueInput*(2*pi);
		set(handles.degField, 'String', num2str(handles.angle/(2*pi)*360))
		set(handles.radField, 'String', num2str(handles.angle))
		set(handles.revField, 'String', num2str(handles.angle/(2*pi)))
guidata(hObject, handles);
end

